
#ifndef VRML_H

/*
 * Simple diagnostic VRML function library for debugging
 *
 * Copyright 2005 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

struct _vrml {

/* Private: */
	FILE *fp;

	double lcent;		/* Gamut L center, usually 50 */

	/* Expandable point arrays */
	struct {
		int npoints;
		int paloc;
		struct {
			double pp[3];			/* Vertex position */
			double cc[3];			/* Vertex color */
			int last;				/* Last vertex of line flag */
		} *pary;
	} set[10];		/* Ten sets */

	/* Expandable triangle vertex index */
	int ntris;
	int taloc;
	struct { int set; int ix[3]; } *tary;

/* Public: */

	/* Methods */

	/* Finish writing the file and free ourselves */
	void (*del)(struct _vrml *s);

	/* Add a spherical marker point to the plot. col == NULL for natural color  */
	/* (Need to do this before or after start_line_set()/dd_vertex()/make_lines() !) */
	void (*add_marker)(struct _vrml *s, double pos[3], double col[3], double rad);

	/* Add a cone marker to the plot. col == NULL for natural color  */
	/* (Need to do this before or after start_line_set()/dd_vertex()/make_lines() !) */
	void (*add_cone)(struct _vrml *s, double p0[3], double p1[3], double col[3], double rad);

	/* Add a text marker to the plot. col == NULL for natural color  */
	/* (Need to do this before or after start_line_set()/dd_vertex()/make_lines() !) */
	void (*add_text)(struct _vrml *s, char *text, double p[3], double col[3], double size);


	/* Start building up verticies that will be converted to lines */
	/* Set can be from 0 - 9 */
	void (*start_line_set)(struct _vrml *s, int set);

	/* Add a verticy (color automatic from Lab position) */
	void (*add_vertex)(struct _vrml *s, int set, double pos[3]);

	/* Add a verticy with color */
	void (*add_col_vertex)(struct _vrml *s, int set, double pos[3], double col[3]);

	/* Turn the last added vertex into the last vertex of the line */
	void (*make_last_vertex)(struct _vrml *s, int set);

	/* Convert the verticies to lines, ppset verticies per line (or using .last flag) */
	/* Use large ppset for just .last flag */
	void (*make_lines)(struct _vrml *s, int set, int ppset);

	/* Add a triangles vertexes defined by vertex index */
	void (*add_triangle)(struct _vrml *s, int set, int ix[3]);

	/* Convert the triangle vertexes to triangles with overall color */
	void (*make_triangles)(struct _vrml *s, int set, double trans, double col[3]);

	/* Convert the triangle vertexes to triangles using vertex colors */
	void (*make_triangles_vc)(struct _vrml *s, int set, double trans);

	/* Clear verticies and triangles */
	void (*clear)(struct _vrml *s);
	
	/* Helper :- convert a Lab value to RGB */
	void (*Lab2RGB)(struct _vrml *s, double *out, double *in);

#ifdef GAMUT_H		/* If gamut.h is #included ahead of us */
	/* Create a solid gamut surface from the given gamut */
	/* trans is trasparency, cc is surface color, cc[0] < 0.0 for natural */
	void (*make_gamut_surface)(struct _vrml *s, gamut *g, double trans, double cc[3]);

	/* Create a solid or wireframe gamut surface from the given gamut */
	/* trans is trasparency, cc is surface color, cc[0] < 0.0 for natural */
	void (*make_gamut_surface_2)(struct _vrml *s, gamut *g, double trans, int wire, double cc[3]);

	/* Add cusp markers from the gamut surface */
	void (*add_cusps)(struct _vrml *s, gamut *g, double trans, double cc[3]);
#endif /* GAMUT_H */

}; typedef struct _vrml vrml;

/* Constructor */
vrml *new_vrml(char *name, int doaxes);

#define VRML_H
#endif /* VRML_H */
