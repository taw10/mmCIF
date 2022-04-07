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


#include <cif.h>
#include <crystfel/stream.h>

int main(int argc, char *argv[])
{
	Stream *st;
	struct image *image;
	struct cif_write_opts_s *write_opts;
	cif_tp *cif = NULL;
	FILE *fh;

	if ( argc < 3 ) {
		fprintf(stderr, "Syntax: %s <input.stream> <output.cif>\n",
		argv[0]);
		return 1;
	}

	st = stream_open_for_read(argv[1]);
	if ( st == NULL ) {
		fprintf(stderr, "Failed to open '%s'\n", argv[1]);
		return 1;
	}

	image = NULL;
	do {

		image = stream_read_chunk(st, STREAM_REFLECTIONS);
		if ( image == NULL ) continue;

	} while ( image != NULL );

	stream_close(st);

	if ( cif_create(&cif) != CIF_OK ) {
		fprintf(stderr, "Failed to create CIF structure\n");
		return 1;
	}

	fh = fopen(argv[2], "w");
	if ( fh == NULL ) {
		fprintf(stderr, "Failed to open '%s'\n", argv[2]);
		return 1;
	}
	if ( cif_write_options_create(&write_opts) != CIF_OK ) {
		fprintf(stderr, "Failed to create write options\n");
		return 1;
	}
	if ( cif_write(fh, write_opts, cif) != CIF_OK ) {
		fprintf(stderr, "Failed to write CIF\n");
		return 1;
	}
	fclose(fh);

	if ( cif_destroy(cif) != CIF_OK ) {
		fprintf(stderr, "Cleanup failed\n");
	}

	return 0;
}
