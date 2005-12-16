#ifndef NEARSMTH_H
#define NEARSMTH_H

/* 
 * nearsmth
 *
 * Gamut mapping support routine that creates a list of
 * nearest mapped points from a source to destination
 * gamut, smoothed to retain reasonably even spacing.
 *
 * Author:  Graeme W. Gill
 * Date:    17/1/2002
 * Version: 1.00
 *
 * Copyright 2002 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */


/* Returned point values. */
/* Note that source == input point, destination == output point */
struct _nearsmth {
	double sdv[3];		/* Smoothed destination value (output) */
	double drv[3];		/* Destination radial value */
	double dr;			/* Destination radial value radius */
	double sv[3];		/* Source value (input) */
	double sr;			/* Source radius */
	struct _nearsmth *n[MAXGAMN+1];	/* List of pointers to neigbors, NULL terminated */
	double nv[MAXGAMN+1][3];		/* source vector to neighbors) */
	int debug;
}; typedef struct _nearsmth nearsmth;


/* Return a list of points. Free list after use */
/* Return NULL on error */
nearsmth *near_smooth(
	int verb,			/* Verbose flag */
	int *npp,			/* Return the number of points returned */
	gamut *sc_gam,		/* Source colorspace gamut */
	gamut *s_gam,		/* Source image gamut (== sc_gam if none) */
	gamut *d_gam,		/* Destination colorspace gamut */
	double a_weight,	/* Absolute error overall weight */
	double a_lweight,	/* Absolute luminance error weight */
	double a_cweight0,	/* Absolute chroma error weight for L = 100 */
	double a_cweight1,	/* Absolute chroma error weight for L = 0 */
	double a_hweight,	/* Absolute hue error weight */
	double r_weight,	/* Relative error overall weight */
	double r_lweight,	/* Relative luminance error weight */
	double r_cweight,	/* Relative chroma error weight */
	double r_hweight,	/* Relative hue error weight */
	double d_weight,	/* Radial error overall weight */
	double d_lweight,	/* Radial luminance error weight */
	double d_cweight,	/* Radial chroma error weight */
	double d_hweight,	/* Radial hue error weight */
	double overshoot,	/* vector overshoot */
	int   usecomp,		/* Flag indicating whether smoothed compressed value will be used */
	int   useexp		/* Flag indicating whether smoothed expanded value will be used */
);

	
#endif /* NEARSMTH_H */

