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

	cif_value_free(val_cif);
	free(name_u);
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
	int i;

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

	for ( i=0; i<n_names; i++ ) {
		free(names[i]);
	}
	free(cat_u);

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
	cif_value_free(val_cif);
	free(name_u);
}


static void cif_packet_set_float(cif_packet_tp *pkt, const char *name,
                                 double val, int dp)
{
	UChar *name_u;
	cif_value_tp *val_cif;

	cif_call(cif_cstr_to_ustr(name, -1, &name_u),
	         "Create name string");
	cif_call(cif_value_create(CIF_NUMB_KIND, &val_cif),
	         "Initialise value");
	cif_call(cif_value_init_numb(val_cif, val, 0, dp, 1),
	         "Failed to set value");
	cif_call(cif_packet_set_item(pkt, name_u, val_cif),
	         "Failed to set packet value");
	cif_value_free(val_cif);
	free(name_u);
}


static void cif_packet_set_string(cif_packet_tp *pkt, const char *name,
                                  const char *val)
{
	UChar *name_u;
	UChar *val_u;
	cif_value_tp *val_cif;

	cif_call(cif_cstr_to_ustr(name, -1, &name_u),
	         "Create name string");
	cif_call(cif_value_create(CIF_CHAR_KIND, &val_cif),
	         "Initialise value");
	cif_call(cif_cstr_to_ustr(val, -1, &val_u),
	         "Create value string");
	cif_call(cif_value_init_char(val_cif, val_u),
	         "Failed to set value");
	cif_call(cif_packet_set_item(pkt, name_u, val_cif),
	         "Failed to set packet value");
	cif_value_free(val_cif);
	free(name_u);
}


static void add_batch_info(cif_loop_tp *batch_loop,
                           int crystal_id,
                           int image_id)
{
	cif_packet_tp *packet;

	cif_call(cif_packet_create(&packet, NULL),
	         "Failed to create loop packet");

	cif_packet_set_int(packet, "_diffrn_batch.id", crystal_id);
	cif_packet_set_int(packet, "_diffrn_batch.diffrn_id", 1);
	cif_packet_set_int(packet, "_diffrn_batch.cell_id", crystal_id);
	cif_packet_set_int(packet, "_diffrn_batch.space_group_id", 1);
	cif_packet_set_int(packet, "_diffrn_batch.orient_matrix_id", crystal_id);
	cif_packet_set_int(packet, "_diffrn_batch.wavelength_id", image_id);
	cif_packet_set_int(packet, "_diffrn_batch.pdbx_image_id", image_id);

	cif_call(cif_loop_add_packet(batch_loop, packet),
	         "Failed to add batch packet to loop");

	cif_packet_free(packet);
}


static void add_image_info(cif_loop_tp *wavelength_loop,
                           int image_id, double wl)
{
	cif_packet_tp *packet;

	cif_call(cif_packet_create(&packet, NULL),
	         "Failed to create loop packet");

	cif_packet_set_int(packet, "_diffrn_radiation_wavelength.id", image_id);
	cif_packet_set_float(packet, "_diffrn_radiation_wavelength.wavelength",
	                     wl*1e10, 5);  /* In Angstroms */

	cif_call(cif_loop_add_packet(wavelength_loop, packet),
	         "Failed to add image packet to loop");

	cif_packet_free(packet);
}


static void add_spg(cif_loop_tp *spg_loop,
                    const char *HM_alt, int IT_num,
                    const char *Hall, const char *crystal_system)
{
	cif_packet_tp *packet;

	cif_call(cif_packet_create(&packet, NULL),
	         "Failed to create loop packet");

	cif_packet_set_int(packet, "_space_group.id", 1);
	cif_packet_set_string(packet, "_space_group.name_H-M_alt", HM_alt);
	cif_packet_set_int(packet, "_space_group.IT_number", IT_num);
	cif_packet_set_string(packet, "_space_group.name_Hall", Hall);
	cif_packet_set_string(packet, "_space_group.crystal_system", crystal_system);

	cif_call(cif_loop_add_packet(spg_loop, packet),
	         "Failed to add SPG packet to loop");

	cif_packet_free(packet);
}


static void add_cell(cif_loop_tp *cell_loop, int crystal_id, UnitCell *cell)
{
	cif_packet_tp *packet;
	double a, b, c, al, be, ga;

	cif_call(cif_packet_create(&packet, NULL),
	         "Failed to create loop packet");

	cell_get_parameters(cell, &a, &b, &c, &al, &be, &ga);

	cif_packet_set_int(packet, "_diffrn_cell.id", crystal_id);
	cif_packet_set_float(packet, "_diffrn_cell.length_a", a*1e10, 7);
	cif_packet_set_float(packet, "_diffrn_cell.length_b", b*1e10, 7);
	cif_packet_set_float(packet, "_diffrn_cell.length_c", c*1e10, 7);
	cif_packet_set_float(packet, "_diffrn_cell.angle_alpha", rad2deg(al), 7);
	cif_packet_set_float(packet, "_diffrn_cell.angle_beta", rad2deg(be), 7);
	cif_packet_set_float(packet, "_diffrn_cell.angle_gamma", rad2deg(ga), 7);

	cif_call(cif_loop_add_packet(cell_loop, packet),
	         "Failed to add cell packet to loop");

	cif_packet_free(packet);
}


static void add_orient(cif_loop_tp *orient_loop, int crystal_id, UnitCell *cell)
{
	cif_packet_tp *packet;
	double ax, ay, az, bx, by, bz, cx, cy, cz;

	cif_call(cif_packet_create(&packet, NULL),
	         "Failed to create loop packet");

	cell_get_reciprocal(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

	cif_packet_set_int(packet, "_diffrn_orient_matrix.id", crystal_id);
	cif_packet_set_int(packet, "_diffrn_orient_matrix.diffrn_id", 1);
	cif_packet_set_float(packet, "_diffrn_orient_matrix.matrix[1][1]", ax/1e10, 8);
	cif_packet_set_float(packet, "_diffrn_orient_matrix.matrix[1][2]", ay/1e10, 8);
	cif_packet_set_float(packet, "_diffrn_orient_matrix.matrix[1][3]", az/1e10, 8);
	cif_packet_set_float(packet, "_diffrn_orient_matrix.matrix[2][1]", bx/1e10, 8);
	cif_packet_set_float(packet, "_diffrn_orient_matrix.matrix[2][2]", by/1e10, 8);
	cif_packet_set_float(packet, "_diffrn_orient_matrix.matrix[2][3]", bz/1e10, 8);
	cif_packet_set_float(packet, "_diffrn_orient_matrix.matrix[3][1]", cx/1e10, 8);
	cif_packet_set_float(packet, "_diffrn_orient_matrix.matrix[3][2]", cy/1e10, 8);
	cif_packet_set_float(packet, "_diffrn_orient_matrix.matrix[3][3]", cz/1e10, 8);
	cif_packet_set_string(packet, "_diffrn_orient_matrix.type", "UB matrix");

	cif_call(cif_loop_add_packet(orient_loop, packet),
	         "Failed to add orientation packet to loop");

	cif_packet_free(packet);
}


static void add_reflection(cif_loop_tp *refl_loop, int refln_id, int crystal_id,
                           signed int h, signed int k, signed int l,
                           double fs, double ss,
                           double intensity, double esd_intensity,
                           const char *panel_name)
{
	cif_packet_tp *packet;

	cif_call(cif_packet_create(&packet, NULL),
	         "Failed to create reflection packet");

	cif_packet_set_int(packet,
	                   "_diffrn_refln.id", refln_id);
	cif_packet_set_int(packet,
	                   "_diffrn_refln.batch_id",
	                   crystal_id);
	cif_packet_set_int(packet,
	                   "_diffrn_refln.diffrn_id",
	                   1);
	cif_packet_set_int(packet,
	                   "_diffrn_refln.index_h", h);
	cif_packet_set_int(packet,
	                   "_diffrn_refln.index_k", k);
	cif_packet_set_int(packet,
	                   "_diffrn_refln.index_l", l);
	cif_packet_set_float(packet,
	                     "_diffrn_refln.intensity_net",
	                     intensity, 2);
	cif_packet_set_float(packet,
	                     "_diffrn_refln.intensity_sigma",
	                     esd_intensity, 2);
	cif_packet_set_float(packet,
	                     "_diffrn_refln.pdbx_detector_calc_fast",
	                     fs, 2);
	cif_packet_set_float(packet,
	                     "_diffrn_refln.pdbx_detector_calc_slow",
	                     ss, 2);
	cif_packet_set_string(packet,
	                      "_diffrn_refln.pdbx_detector_array_id",
	                      panel_name);

	cif_call(cif_loop_add_packet(refl_loop, packet),
	         "Failed to add reflection packet");

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
	cif_loop_tp *batch_loop;
	cif_loop_tp *cell_loop;
	cif_loop_tp *orient_loop;
	cif_loop_tp *spg_loop;
	cif_loop_tp *wavelength_loop;
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
	                        "_diffrn_refln.batch_id",
	                        "_diffrn_refln.diffrn_id",
	                        "_diffrn_refln.index_h",
	                        "_diffrn_refln.index_k",
	                        "_diffrn_refln.index_l",
	                        "_diffrn_refln.intensity_net",
	                        "_diffrn_refln.intensity_sigma",
	                        "_diffrn_refln.pdbx_detector_calc_fast",
	                        "_diffrn_refln.pdbx_detector_calc_slow",
	                        "_diffrn_refln.pdbx_detector_array_id",
	                        NULL);

	batch_loop = create_loop(block, "_diffrn_batch",
	                         "_diffrn_batch.id",
	                         "_diffrn_batch.diffrn_id",
	                         "_diffrn_batch.cell_id",
	                         "_diffrn_batch.space_group_id",
	                         "_diffrn_batch.orient_matrix_id",
	                         "_diffrn_batch.wavelength_id",
	                         "_diffrn_batch.pdbx_image_id",
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

	wavelength_loop = create_loop(block, "_diffrn_radiation_wavelength",
	                              "_diffrn_radiation_wavelength.id",
	                              "_diffrn_radiation_wavelength.wavelength",
	                              NULL);

	spg_loop = create_loop(block, "_space_group",
	                       "_space_group.id",
	                       "_space_group.name_H-M_alt",
	                       "_space_group.IT_number",
	                       "_space_group.name_Hall",
	                       "_space_group.crystal_system",
	                       NULL);

	add_spg(spg_loop, "C 2 2 21", 20, "C 2c 2", "Orthorhombic");

	orient_loop = create_loop(block, "_diffrn_orient_matrix",
	                          "_diffrn_orient_matrix.id",
	                          "_diffrn_orient_matrix.diffrn_id",
	                          "_diffrn_orient_matrix.matrix[1][1]",
	                          "_diffrn_orient_matrix.matrix[1][2]",
	                          "_diffrn_orient_matrix.matrix[1][3]",
	                          "_diffrn_orient_matrix.matrix[2][1]",
	                          "_diffrn_orient_matrix.matrix[2][2]",
	                          "_diffrn_orient_matrix.matrix[2][3]",
	                          "_diffrn_orient_matrix.matrix[3][1]",
	                          "_diffrn_orient_matrix.matrix[3][2]",
	                          "_diffrn_orient_matrix.matrix[3][3]",
	                          "_diffrn_orient_matrix.type",
	                          NULL);

	image = NULL;
	image_id = 1;
	crystal_id = 1;
	refln_id = 1;
	do {

		int i;

		image = stream_read_chunk(st, STREAM_REFLECTIONS | STREAM_DATA_DETGEOM);
		if ( image == NULL ) continue;

		if ( image->n_crystals == 0 ) {
			image_free(image);
			continue;
		}

		add_image_info(wavelength_loop, image_id, image->lambda);

		for ( i=0; i<image->n_crystals; i++ ) {

			Reflection *refl;
			RefListIterator *iter;
			Crystal *cr = image->crystals[i];

			add_batch_info(batch_loop, crystal_id, image_id);
			add_cell(cell_loop, crystal_id, crystal_get_cell(cr));
			add_orient(orient_loop, crystal_id, crystal_get_cell(cr));

			for ( refl = first_refl(crystal_get_reflections(cr), &iter);
			      refl != NULL;
			      refl = next_refl(refl, iter) )
			{
				signed int h, k, l;
				double fs, ss;
				int panel;

				get_indices(refl, &h, &k,  &l);
				get_detector_pos(refl, &fs, &ss);
				panel = get_panel_number(refl);

				add_reflection(refl_loop, refln_id,
				               crystal_id, h, k,l,
				               fs, ss,
				               get_intensity(refl),
				               get_esd_intensity(refl),
				               image->detgeom->panels[panel].name);

				refln_id++;

			}

			crystal_id++;
		}

		image_id++;

		image_free(image);

	} while ( image != NULL );

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

	free(write_opts);
	cif_loop_free(refl_loop);
	cif_loop_free(batch_loop);
	cif_loop_free(spg_loop);
	cif_loop_free(wavelength_loop);
	cif_loop_free(orient_loop);
	cif_loop_free(cell_loop);
	cif_call(cif_container_destroy(block),
	         "Failed to clean up block");
	cif_call(cif_destroy(cif), "Cleanup failed");

	return 0;
}
