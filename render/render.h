
#ifndef RENDER2D_H
#define RENDER2D_H

/* 
 * render2d
 *
 * Simple 2D raster rendering support.
 *
 * Author:  Graeme W. Gill
 * Date:    28/12/2005
 *
 * Copyright 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* Mathematical coordinate in mm are used for primitives, ie. the origin is */
/* the bottom left corner. */
/* Device color values range from 0.0 to 1.0 */

#define MXCH2D 4		/* Maximum color channels */

/* Color type */
/* Shouldn't this be an xcolorants mask ? */
typedef enum {
	w_2d,			/* Video style grey */
	k_2d,			/* Printing style grey */
	rgb_2d,			/* RGB */
	cmyk_2d			/* CMYK */
} colort2d;

/* Pixel depth */
typedef enum {
	bpc8_2d,		/* 8 bits per component */
	bpc16_2d		/* 16 bits per component */
} depth2d;

typedef double color2d[MXCH2D]; 

/* Primitive type */
//typedef enum {
//	rect		/* Solid color rectangle */
//} primt2d;

/* ------------------------------------ */
#define PRIM_STRUCT							\
/*	primt2d tag; */			/* Type of primitive */	\
	struct _prim2d *next;	/* Linked list to next primitive */ \
	double x0, y0, x1, y1;	/* Extent, top & left inclusive, bot & right non-inclusive */ \
	void (*del)(struct _prim2d *s);		/* Delete the object */ \
							/* Render the object at location. Return nz if in primitive */ \
	int (*rend)(struct _prim2d *s, color2d rv, double x, double y);

struct _prim2d {
	PRIM_STRUCT
}; typedef struct _prim2d prim2d;

/* ------------------------------------ */
/* Solid rectange primitive */
struct _rect2d {
	PRIM_STRUCT
	color2d c;			/* Color of rectangle */
}; typedef struct _rect2d rect2d;

prim2d *new_rect2d(double x, double y, double w, double h, color2d c);

/* ------------------------------------ */
/* Vertex shaded rectange */
struct _rectvs2d {
	PRIM_STRUCT
	color2d c[4];	/* Bot left, bot right, top left, top right */
	int x_blend;	/* Blending rule flags, 0 = linear, 1 = spline, 2 = sine */
	int y_blend;
	int y_sine;
}; typedef struct _rectvs2d rectvs2d;

prim2d *new_rectvs2d(double x, double y, double w, double h, color2d c[4]);

/* ------------------------------------ */
/* Vertex shaded triangle */
struct _trivs2d {
	PRIM_STRUCT
	double be[3][3];	/* baricentric equations */
	color2d c[3];		/* Color of each vertex */
}; typedef struct _trivs2d trivs2d;

prim2d *new_trivs2d(double v[3][2], color2d c[3]);

/* ------------------------------------ */
/* Render object */

struct _render2d {

/* Private: */
	double w, h;		/* Page size in mm */
	double hres, vres;	/* Page pixel resolution in pixels/mm */
	int pw, ph;			/* Page size in pixels */
	colort2d csp;		/* Color space */
	int      ncc;		/* Number of color components */
	depth2d  dpth;		/* Depth of the components */

	color2d defc;		/* Default color value */

	prim2d *head;		/* Start of list of primitives in rendering order */

/* Public: */
	/* Methods */
	void (*del)(struct _render2d *s);					/* Free ourselves and all primitives */

	void (*set_defc)(struct _render2d *s, color2d c);	/* Set the default color */

	void (*add)(struct _render2d *s, prim2d *p);		/* Add a primitive */

	int (*write)(struct _render2d *s, char *filename);	/* Render and write to a TIFF file */
}; typedef struct _render2d render2d;

/* Constructor */
/* Sizes are in mm, resolutions are in pixels/mm */
render2d *new_render2d(double w, double h, double hres, double vres, colort2d csp, depth2d dpth);

#endif /* RENDER2D_H */

