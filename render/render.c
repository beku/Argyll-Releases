
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

#undef DEBUG

#define verbo stdout

#include <stdio.h>
#include <math.h>
#include "copyright.h"
#include "config.h"
#include "numsup.h"
#include "tiffio.h"
#include "render.h"

/* ------------------------------------------------------------- */
/* Main class implementation */

/* Free ourselves and all primitives */
static void render2d_del(render2d *s) {
	prim2d *th, *nx;

	/* Delete all the primitives */
	for (th = s->head; th != NULL; th = nx) {
		nx = th->next;
		th->del(th);
	}

	/* And then ourselves */
	free(s);
}

/* Add a primitive */
static void render2d_add(render2d *s, prim2d *p) {
	if (p == NULL)
		error("render2d: Adding NULL primitive");

	p->next = s->head;
	s->head = p;
}

/* Set the default color */
static void render2d_set_defc(render2d *s, color2d c) {
	int j;

	for (j = 0; j < MXCH2D; j++)
		s->defc[j] = c[j];
}

/* Compute the length of a double nul terminated string, including */
/* the nuls. */
static int zzstrlen(char *s) {
	int i;
	for (i = 0;; i++) {
		if (s[i] == '\000' && s[i+1] == '\000')
			return i+2;
	}
	return 0;
}

/* Render and write to a TIFF file */
/* Return NZ on error */
static int render2d_write(render2d *s, char *filename) {
	TIFF *wh = NULL;
	uint16 samplesperpixel = 0, bitspersample = 0;
	uint16 photometric = 0;
	uint16 inkset = 0xffff;
	char *inknames = NULL;
	tdata_t *outbuf;
	double rx, ry;			/* Read x & y */
	int x, y;				/* Pixel x & y */

	switch (s->csp) {
		case w_2d:			/* Video style grey */
			samplesperpixel = 1;
			photometric = PHOTOMETRIC_MINISBLACK;
			break;
		case k_2d:			/* Printing style grey */
			samplesperpixel = 1;
			photometric = PHOTOMETRIC_MINISWHITE;
			break;
		case rgb_2d:		/* RGB */
			samplesperpixel = 3;
			photometric = PHOTOMETRIC_RGB;
			break;
		case cmyk_2d:		/* CMYK */
			samplesperpixel = 4;
			photometric = PHOTOMETRIC_SEPARATED;
			inkset = INKSET_CMYK;
			inknames = "cyan\000magenta\000yellow\000\000";
			break;
		default:
			error("Illegal colorspace in file '%s'",filename);
	}
	s->ncc = samplesperpixel;

	switch (s->dpth) {
		case bpc8_2d:		/* 8 bits per component */
			bitspersample = 8;
			break;
		case bpc16_2d:		/* 16 bits per component */
			bitspersample = 16;
			break;
		default:
			error("Illegal bits per component in file '%s'",filename);
	}

	if ((wh = TIFFOpen(filename, "w")) == NULL)
		error("Can\'t create TIFF file '%s'!",filename);
	
	TIFFSetField(wh, TIFFTAG_IMAGEWIDTH,  s->pw);
	TIFFSetField(wh, TIFFTAG_IMAGELENGTH, s->ph);
	TIFFSetField(wh, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(wh, TIFFTAG_SAMPLESPERPIXEL, samplesperpixel);
	TIFFSetField(wh, TIFFTAG_BITSPERSAMPLE, bitspersample);
	TIFFSetField(wh, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(wh, TIFFTAG_PHOTOMETRIC, photometric);

	if (inknames != NULL) {
		int inlen = zzstrlen(inknames);
		TIFFSetField(wh, TIFFTAG_INKSET, inkset);
		TIFFSetField(wh, TIFFTAG_INKNAMES, inlen, inknames);
	}
	TIFFSetField(wh, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(wh, TIFFTAG_RESOLUTIONUNIT, RESUNIT_CENTIMETER);
	TIFFSetField(wh, TIFFTAG_XRESOLUTION, 10.0 * s->hres);
	TIFFSetField(wh, TIFFTAG_YRESOLUTION, 10.0 * s->vres);
	TIFFSetField(wh, TIFFTAG_IMAGEDESCRIPTION, "Test chart created with Argyll");

	/* Allocate one line buffer */
	outbuf = _TIFFmalloc(TIFFScanlineSize(wh));

	/* Render each line and write it */
	/* (This could be speeded up by using a depth bufer, */
	/* tagging each primitive with its order depth to, */
	/* allow arbitrary primitive testing, and then */
	/* using a Y edge list to reduce the number */
	/* of primitives tested at each line.) */
	for (y = 0; y < s->ph; y++) {
		ry = ((s->ph-1) - y) / s->vres;

		for (x = 0; x < s->pw; x++) {
			int j;
			color2d rv;
			prim2d *th;

			rx = x / s->hres;

			for (th = s->head; th != NULL; th = th->next) {
				if (ry < th->y0 || ry > th->y1
				 || rx < th->x0 || rx > th->x1)
					continue;
				if (th->rend(th, rv, rx, ry))
					break;
			}
			if (th == NULL) {	/* Apply default color */
				for (j = 0; j < s->ncc; j++)
					rv[j] = s->defc[j];
			}
			if (s->dpth == bpc8_2d) {
				unsigned char *p = ((unsigned char *)outbuf) + x * s->ncc;
				for (j = 0; j < s->ncc; j++)
					p[j] = (int)(255.0 * rv[j] + 0.5);
			} else {
				unsigned short *p = ((unsigned short *)outbuf) + x * s->ncc;
				for (j = 0; j < s->ncc; j++)
					p[j] = (int)(65525.0 * rv[j] + 0.5);
			}
		}

		if (TIFFWriteScanline(wh, outbuf, y, 0) < 0)
			error ("Failed to write TIFF file '%s' line %d",filename,y);
	}

	_TIFFfree(outbuf);
	TIFFClose(wh);		/* Close Output file */

	return 0;
}

/* Constructor */
render2d *new_render2d(
double w,
double h,
double hres,
double vres,
colort2d csp,
depth2d dpth
) {
	render2d *s;

	if ((s = (render2d *)calloc(1, sizeof(render2d))) == NULL) {
		return NULL;
	}

	s->w = w;
	s->h = h;
	s->hres = hres;
	s->vres = vres;
	s->csp = csp;
	s->dpth = dpth;

	s->del = render2d_del;
	s->set_defc = render2d_set_defc;
	s->add = render2d_add;
	s->write = render2d_write;

	/* Figure the raster size and actuall size */
	s->pw = (int)(s->hres * w + 0.5);
	s->ph = (int)(s->vres * h + 0.5);
	s->w = s->pw * s->hres;
	s->h = s->ph * s->vres;

	return s;
}

/* ------------------------------------------------------------- */
/* Primitive implementations */

/* Primitive destructor for no sub-allocation */
static void prim2d_del(prim2d *s) {
	free(s);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Flat shaded rectangle */

/* Render the rectangle object at location. Return nz if in primitive */
static int rect2d_rend(prim2d *ss, color2d rv, double x, double y) {
	rect2d *s = (rect2d *)ss;
	int j;

	if (y < s->y0 || y > s->y1
	 || x < s->x0 || x > s->x1)
		return 0;

	for (j = 0; j < MXCH2D; j++)
		rv[j] = s->c[j];

	return 1;
}

prim2d *new_rect2d(
double x,
double y,
double w,
double h,
color2d c
) {
	int j;
	rect2d *s;

	if ((s = (rect2d *)calloc(1, sizeof(rect2d))) == NULL) {
		return NULL;
	}

	s->del = prim2d_del; 
	s->rend = rect2d_rend; 

	s->x0 = x;
	s->y0 = y;
	s->x1 = x + w;
	s->y1 = y + h;

	for (j = 0; j < MXCH2D; j++)
		s->c[j] = c[j];

	return (prim2d *)s;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Vertex shaded rectangle */

/* Render the rectangle object at location. Return nz if in primitive */
static int rectvs2d_rend(prim2d *ss, color2d rv, double x, double y) {
	rectvs2d *s = (rectvs2d *)ss;
	double bx, by, b[4];
	int i, j;

	if (y < s->y0 || y > s->y1
	 || x < s->x0 || x > s->x1)
		return 0;

	/* Compute linear blend */
	bx = (x - s->x0)/(s->x1 - s->x0);
	by = (y - s->y0)/(s->y1 - s->y0);

	if (s->x_blend == 1) {
		bx = bx * bx * (3.0 - 2.0 * bx);	/* Cubic spline */
	} else if (s->x_blend == 2) {
		bx = bx - 0.5;						/* Sine */
		bx *= 3.141592654;
		bx = sin(bx);
		bx = 0.5 + (0.5 * bx);
	}
	if (s->y_blend == 1) {
		by = by * by * (3.0 - 2.0 * by);	/* Spline */
	} else if (s->y_blend == 2) {
		double ty;
		ty = by * by * (3.0 - 2.0 * by);	/* spline at y == 1, linear at y == 0 */
		by = by * ty + (1.0 - by) * by;
	} else if (s->y_blend == 3) {
		double ty;
		ty = by * by * (3.0 - 2.0 * by);	/* linear at y == 1, spline at y == 0 */
		by = (1.0 - by) * ty + by * by;
	}

	/* Compute 2d blend */
	b[0] = (1.0 - by) * (1.0 - bx);
	b[1] = (1.0 - by) * bx;
	b[2] = by * (1.0 - bx);
	b[3] = by * bx;

	/* Compute the resulting color */
	for (j = 0; j < MXCH2D; j++) {
		rv[j] = 0.0;
		for (i = 0; i < 4; i++)
			rv[j] += b[i] * s->c[i][j];
		rv[j] = rv[j];
	}

	return 1;
}

prim2d *new_rectvs2d(
double x,
double y,
double w,
double h,
color2d c[4]
) {
	int i, j;
	rectvs2d *s;

	if ((s = (rectvs2d *)calloc(1, sizeof(rectvs2d))) == NULL) {
		return NULL;
	}

	s->del = prim2d_del; 
	s->rend = rectvs2d_rend; 

	s->x0 = x;
	s->y0 = y;
	s->x1 = x + w;
	s->y1 = y + h;

	for (i = 0; i < 4; i++)
		for (j = 0; j < MXCH2D; j++)
			s->c[i][j] = c[i][j];

	return (prim2d *)s;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Vertex shaded triangle */

/* Render the triangle object at location. Return nz if in primitive */
static int trivs2d_rend(prim2d *ss, color2d rv, double x, double y) {
	trivs2d *s = (trivs2d *)ss;
	double b[3];
	int i, j;

	/* Compute the baricentric values for the input point, */
	/* and reject this point if it is outside the triangle */
	for (i = 0; i < 3; i++) {
		b[i] = s->be[i][0] * x + s->be[i][1] * y + s->be[i][2]; 
		if (b[i] < 0.0 || b[i] > 1.0)
			return 0;
	}

#ifdef NEVER		/* Experiment with smoothing hue blending */
	{
		double ss = b[1] + b[2];
		if (ss > 1e-6) {
			b[1] /= ss;
			b[2] /= ss;
//			b[1] = b[1] * b[1] * (3.0 - 2.0 * b[1]);

			b[1] = b[1] - 0.5;
			b[1] *= 3.141592654;
			b[1] = sin(b[1]);
			b[1] = 0.5 + (0.5 * b[1]);

			b[2] = 1.0 - b[1];
			b[1] *= ss;
			b[2] *= ss;
		}
	}
#endif

	/* Compute the resulting color */
	for (j = 0; j < MXCH2D; j++) {
		rv[j] = 0.0;
		for (i = 0; i < 3; i++)
			rv[j] += b[i] * s->c[i][j];
	}

	return 1;
}

static int inverse3x3(double out[3][3], double in[3][3]);

prim2d *new_trivs2d(
double v[3][2],			/* Vertex locations */
color2d c[3]			/* Corresponding colors */
) {
	int i, j;
	trivs2d *s;
	double tt[3][3];

	if ((s = (trivs2d *)calloc(1, sizeof(trivs2d))) == NULL) {
		return NULL;
	}

	s->del = prim2d_del; 
	s->rend = trivs2d_rend; 

	/* Set the bouding box values */
	s->x0 = s->y0 = 1e38;	
	s->x1 = s->y1 = -1e38;	
	for (i = 0; i < 3; i++) {
		if (v[i][0] < s->x0)
			s->x0 = v[i][0];
		if (v[i][1] < s->y0)
			s->y0 = v[i][1];
		if (v[i][0] > s->x1)
			s->x1 = v[i][0];
		if (v[i][1] > s->y1)
			s->y1 = v[i][1];
	}

	/* Compute baricentric equations */
	for (i = 0; i < 3; i++) {
		tt[0][i] = v[i][0];
		tt[1][i] = v[i][1];
		tt[2][i] = 1.0;
	}
	if (inverse3x3(s->be, tt))
		error("trivs2d: Matrix inversion failed");

	/* Copy vertex colors */
	for (i = 0; i < 3; i++) {
		for (j = 0; j < MXCH2D; j++)
			s->c[i][j] = c[i][j];
	}

	return (prim2d *)s;
}


/* ==================================================== */
/* Misc. support functions.                      */

/* 
	Matrix Inversion
	by Richard Carling
	from "Graphics Gems", Academic Press, 1990
*/

/* 
 *   adjoint( original_matrix, inverse_matrix )
 * 
 *     calculate the adjoint of a 3x3 matrix
 *
 *      Let  a   denote the minor determinant of matrix A obtained by
 *           ij
 *
 *      deleting the ith row and jth column from A.
 *
 *                    i+j
 *     Let  b   = (-1)    a
 *          ij            ji
 *
 *    The matrix B = (b  ) is the adjoint of A
 *                     ij
 */

#define det2x2(a, b, c, d) (a * d - b * c)

static void adjoint(
double out[3][3],
double in[3][3]
) {
    double a1, a2, a3, b1, b2, b3, c1, c2, c3;

    /* assign to individual variable names to aid  */
    /* selecting correct values  */

	a1 = in[0][0]; b1 = in[0][1]; c1 = in[0][2];
	a2 = in[1][0]; b2 = in[1][1]; c2 = in[1][2];
	a3 = in[2][0]; b3 = in[2][1]; c3 = in[2][2];

    /* row column labeling reversed since we transpose rows & columns */

    out[0][0]  =   det2x2(b2, b3, c2, c3);
    out[1][0]  = - det2x2(a2, a3, c2, c3);
    out[2][0]  =   det2x2(a2, a3, b2, b3);
        
    out[0][1]  = - det2x2(b1, b3, c1, c3);
    out[1][1]  =   det2x2(a1, a3, c1, c3);
    out[2][1]  = - det2x2(a1, a3, b1, b3);
        
    out[0][2]  =   det2x2(b1, b2, c1, c2);
    out[1][2]  = - det2x2(a1, a2, c1, c2);
    out[2][2]  =   det2x2(a1, a2, b1, b2);
}

/*
 * double = det3x3(  a1, a2, a3, b1, b2, b3, c1, c2, c3 )
 * 
 * calculate the determinant of a 3x3 matrix
 * in the form
 *
 *     | a1,  b1,  c1 |
 *     | a2,  b2,  c2 |
 *     | a3,  b3,  c3 |
 */

static double det3x3(double in[3][3]) {
    double a1, a2, a3, b1, b2, b3, c1, c2, c3;
    double ans;

	a1 = in[0][0]; b1 = in[0][1]; c1 = in[0][2];
	a2 = in[1][0]; b2 = in[1][1]; c2 = in[1][2];
	a3 = in[2][0]; b3 = in[2][1]; c3 = in[2][2];

    ans = a1 * det2x2(b2, b3, c2, c3)
        - b1 * det2x2(a2, a3, c2, c3)
        + c1 * det2x2(a2, a3, b2, b3);
    return ans;
}

#define SMALL_NUMBER	1.e-8
/* 
 *   inverse( original_matrix, inverse_matrix )
 * 
 *    calculate the inverse of a 4x4 matrix
 *
 *     -1     
 *     A  = ___1__ adjoint A
 *         det A
 */

/* Return non-zero if not invertable */
static int inverse3x3(
double out[3][3],
double in[3][3]
) {
    int i, j;
    double det;

    /*  calculate the 3x3 determinant
     *  if the determinant is zero, 
     *  then the inverse matrix is not unique.
     */
    det = det3x3(in);

    if ( fabs(det) < SMALL_NUMBER)
        return 1;

    /* calculate the adjoint matrix */
    adjoint(out, in);

    /* scale the adjoint matrix to get the inverse */
    for (i = 0; i < 3; i++)
        for(j = 0; j < 3; j++)
		    out[i][j] /= det;
	return 0;
}

