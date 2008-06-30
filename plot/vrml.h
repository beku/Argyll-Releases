
#ifndef VRML_H

/*
 * Simple diagnostic VRML function library for debugging
 *
 * Copyright 2005 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

struct _vrml {

/* Private: */
	FILE *fp;

	double lcent;		/* Gamut L center, usually 50 */

	/* Expandable point array */
	int npoints;
	int paloc;
	struct { double pp[3]; double cc[3]; } *pary;

/* Public: */

	/* Methods */

	/* Finish writing the file and free ourselves */
	void (*del)(struct _vrml *s);

	/* Add a spherical marker point to the plot. col == NULL for natural color  */
	void (*add_marker)(struct _vrml *s, double pos[3], double col[3], double rad);

	/* Start building up verticies that will be converted to lines */
	void (*start_line_set)(struct _vrml *s);

	/* Add a verticy (color automatic from Lab position) */
	void (*add_vertex)(struct _vrml *s, double pos[3]);

	/* Add a verticy with color */
	void (*add_col_vertex)(struct _vrml *s, double pos[3], double col[3]);

	/* Convert the verticies to lines, ppset verticies per line */
	void (*make_lines)(struct _vrml *s, int ppset);

	/* Helper :- convert a Lab value to RGB */
	void (*Lab2RGB)(struct _vrml *s, double *out, double *in);

#ifdef GAMUT_H		/* If gamut.h is #included ahead of us */
	/* Create a gamut surface from the given gamut */
	/* trans is trasparency, cc is surface color, cc[0] < 0.0 for natural */
	void (*make_gamut_surface)(struct _vrml *s, gamut *g, double trans, double cc[3]);
#endif /* GAMUT_H */

}; typedef struct _vrml vrml;

/* Constructor */
vrml *new_vrml(char *name, int doaxes);

#define VRML_H
#endif /* VRML_H */
