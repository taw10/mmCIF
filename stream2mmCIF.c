/*
 * stream2mmCIF.c
 *
 * Convert a stream to mmCIF with unmerged intensities
 *
 * Copyright © 2022 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2022 Thomas White <taw@physics.org>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include <stdio.h>
#include <stdarg.h>
#include <assert.h>

#include <unicode/ustring.h>
#include <cif.h>

#include <crystfel/stream.h>
#include <crystfel/utils.h>
#include <crystfel/image.h>


static void cif_call(int rval, const char *complaint)
{
	if ( rval != CIF_OK ) {
		fprintf(stderr, "%s (%i)\n", complaint, rval);
		abort();
	}
}


static int add_string_val(cif_container_tp *container,
                          const char *name, const char *val)
{
	UChar *name_u;
	UChar *val_u;
	cif_value_tp *val_cif = NULL;

	cif_call(cif_cstr_to_ustr(name, -1, &name_u), "Create string");
	cif_call(cif_cstr_to_ustr(val, -1, &val_u), "Create string");

	cif_call(cif_value_create(CIF_CHAR_KIND, &val_cif),
	         "Initialise value");
	cif_call(cif_value_init_char(val_cif, val_u),
	         "Failed to set value");
	cif_call(cif_container_set_value(container, name_u, val_cif),
	         "Failed to set value");
	return 0;
}


static cif_loop_tp *create_loop(cif_container_tp *container,
                                const char *category,
                                ...)
{
	va_list ap;
	int n_names = 0;
	const char *name;
	cif_loop_tp *loop;
	UChar *names[32];
	UChar *cat_u;

	va_start(ap, category);

	name = va_arg(ap, char *);
	while ( name != NULL ) {

		assert(n_names < 32);
		cif_call(cif_cstr_to_ustr(name, -1, &names[n_names]),
		         "Failed to create name string");

		n_names++;
		name = va_arg(ap, char *);
	}

	va_end(ap);

	names[n_names] = NULL;

	cif_call(cif_cstr_to_ustr(category, -1, &cat_u), "Create string");
	cif_call(cif_container_create_loop(container, cat_u, names, &loop),
	         "Failed to create loop");

	return loop;
}


static void cif_packet_set_int(cif_packet_tp *pkt, const char *name, int val)
{
	UChar *name_u;
	cif_value_tp *val_cif;

	cif_call(cif_cstr_to_ustr(name, -1, &name_u),
	         "Create name string");
	cif_call(cif_value_create(CIF_NUMB_KIND, &val_cif),
	         "Initialise value");
	cif_call(cif_value_init_numb(val_cif, val, 0, 0, 0),
	         "Failed to set value");
	cif_call(cif_packet_set_item(pkt, name_u, val_cif),
	         "Failed to set packet value");
}


static void cif_packet_set_float(cif_packet_tp *pkt, const char *name,
                                 double val)
{
	UChar *name_u;
	cif_value_tp *val_cif;

	cif_call(cif_cstr_to_ustr(name, -1, &name_u),
	         "Create name string");
	cif_call(cif_value_create(CIF_NUMB_KIND, &val_cif),
	         "Initialise value");
	cif_call(cif_value_init_numb(val_cif, val, 0, 3, 1),
	         "Failed to set value");
	cif_call(cif_packet_set_item(pkt, name_u, val_cif),
	         "Failed to set packet value");
}


static void add_crystal_info(cif_loop_tp *crystal_loop, int crystal_id,
                             int image_id)
{
	cif_packet_tp *packet;

	cif_call(cif_packet_create(&packet, NULL),
	         "Failed to create loop packet");

	cif_packet_set_int(packet, "_diffrn_batch.diffrn_id", crystal_id);
	cif_packet_set_int(packet, "_diffrn_batch.pdbx_image_id", image_id);
	cif_packet_set_int(packet, "_diffrn_batch.cell_id", crystal_id);

	cif_call(cif_loop_add_packet(crystal_loop, packet),
	         "Failed to add packet to loop");

	cif_packet_free(packet);
}


static void add_cell(cif_loop_tp *loop, int crystal_id, UnitCell *cell)
{
	cif_packet_tp *packet;
	double a, b, c, al, be, ga;

	cif_call(cif_packet_create(&packet, NULL),
	         "Failed to create loop packet");

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);

	cif_packet_set_int(packet, "_diffrn_cell.id", crystal_id);
	cif_packet_set_float(packet, "_diffrn_cell.length_a", a*1e10);
	cif_packet_set_float(packet, "_diffrn_cell.length_b", b*1e10);
	cif_packet_set_float(packet, "_diffrn_cell.length_c", c*1e10);
	cif_packet_set_float(packet, "_diffrn_cell.angle_alpha", rad2deg(al));
	cif_packet_set_float(packet, "_diffrn_cell.angle_beta", rad2deg(be));
	cif_packet_set_float(packet, "_diffrn_cell.angle_gamma", rad2deg(ga));

	cif_call(cif_loop_add_packet(loop, packet),
	         "Failed to add packet to loop");

	cif_packet_free(packet);
}


int main(int argc, char *argv[])
{
	Stream *st;
	struct image *image;
	struct cif_write_opts_s *write_opts;
	cif_tp *cif = NULL;
	FILE *fh;
	cif_loop_tp *refl_loop;
	cif_loop_tp *crystal_loop;
	cif_loop_tp *cell_loop;
	long int refln_id;
	long int image_id;
	long int crystal_id;

	if ( argc < 3 ) {
		fprintf(stderr, "Syntax: %s <input.stream> <output.cif>\n",
		argv[0]);
		return 1;
	}

	cif_call(cif_create(&cif),
	         "Failed to create CIF structure\n");

	st = stream_open_for_read(argv[1]);
	if ( st == NULL ) {
		fprintf(stderr, "Failed to open '%s'\n", argv[1]);
		return 1;
	}

	/* Start data block */
	cif_block_tp *block = NULL;
	U_STRING_DECL(u_data_name, "fromstream", 10);
	cif_call(cif_create_block(cif, u_data_name, &block),
	         "Failed to create CIF data block\n");

	/* Audit info */
	char tmp[128];
	snprintf(tmp, 128, "Converted from %s using stream2mmCIF", argv[1]);
	add_string_val(block, "_audit_creation_method", tmp);

	refl_loop = create_loop(block, "_diffrn_refln",
	                        "_diffrn_refln.id",
	                        "_diffrn_refln.diffrn_id",
	                        "_diffrn_refln.index_h",
	                        "_diffrn_refln.index_k",
	                        "_diffrn_refln.index_l",
	                        "_diffrn_refln.pdbx_image_id",
	                        "_diffrn_refln.intensity_net",
	                        "_diffrn_refln.intensity_sigma",
	                        "_diffrn_refln.detector_x",
	                        "_diffrn_refln.detector_y",
	                        NULL);

	crystal_loop = create_loop(block, "_diffrn_batch",
	                           "_diffrn_batch.diffrn_id",
	                           "_diffrn_batch.pdbx_image_id",
	                           "_diffrn_batch.cell_id",
	                           NULL);

	cell_loop = create_loop(block, "_diffrn_cell",
	                        "_diffrn_cell.id",
	                        "_diffrn_cell.length_a",
	                        "_diffrn_cell.length_b",
	                        "_diffrn_cell.length_c",
	                        "_diffrn_cell.angle_alpha",
	                        "_diffrn_cell.angle_beta",
	                        "_diffrn_cell.angle_gamma",
	                        NULL);

	image = NULL;
	image_id = 1;
	crystal_id = 1;
	refln_id = 1;
	do {

		int i;

		image = stream_read_chunk(st, STREAM_REFLECTIONS);
		if ( image == NULL ) continue;

		if ( image->n_crystals == 0 ) continue;

		for ( i=0; i<image->n_crystals; i++ ) {

			Reflection *refl;
			RefListIterator *iter;
			Crystal *cr = image->crystals[i];
			add_crystal_info(crystal_loop, crystal_id, image_id);
			add_cell(cell_loop, crystal_id, crystal_get_cell(cr));

			for ( refl = first_refl(crystal_get_reflections(cr), &iter);
			      refl != NULL;
			      refl = next_refl(refl, iter) )
			{
				signed int h, k, l;
				cif_packet_tp *packet;

				get_indices(refl, &h, &k,  &l);

				cif_call(cif_packet_create(&packet, NULL),
				         "Failed to create reflection packet");

				cif_packet_set_int(packet,
				                   "_diffrn_refln.id", refln_id);
				cif_packet_set_int(packet,
				                   "_diffrn_refln.diffrn_id",
				                   crystal_id);
				cif_packet_set_int(packet,
				                   "_diffrn_refln.index_h", h);
				cif_packet_set_int(packet,
				                   "_diffrn_refln.index_k", k);
				cif_packet_set_int(packet,
				                   "_diffrn_refln.index_l", l);
				cif_packet_set_int(packet,
				                   "_diffrn_refln.pdbx_image_id",
				                   image_id);
				cif_packet_set_float(packet,
				                   "_diffrn_refln.intensity_net",
				                   get_intensity(refl));
				cif_packet_set_float(packet,
				                   "_diffrn_refln.intensity_sigma",
				                   get_esd_intensity(refl));

				cif_call(cif_loop_add_packet(refl_loop, packet),
				         "Failed to add reflection packet");

				cif_packet_free(packet);
				refln_id++;

			}

			crystal_id++;
		}

		image_id++;

	} while ( image != NULL );

	//cif_loop_free(refl_loop);
	cif_loop_free(crystal_loop);

	stream_close(st);

	fh = fopen(argv[2], "w");
	if ( fh == NULL ) {
		fprintf(stderr, "Failed to open '%s'\n", argv[2]);
		return 1;
	}
	cif_call(cif_write_options_create(&write_opts),
	         "Failed to create write options");
	cif_call(cif_write(fh, write_opts, cif),
	         "Failed to write CIF");
	fclose(fh);

	cif_call(cif_destroy(cif), "Cleanup failed");

	return 0;
}
