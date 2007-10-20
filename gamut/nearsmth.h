#ifndef NEARSMTH_H
#define NEARSMTH_H

/* 
 * nearsmth
 *
 * Gamut mapping support routine that creates a list of
 * guide vectors from the source to destination
 * gamut, smoothed to retain reasonably even spacing.
 *
 * Author:  Graeme W. Gill
 * Date:    17/1/2002
 * Version: 1.00
 *
 * Copyright 2002 - 2006 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */


#define NNB 6				/* Maxmimum number of near points */

/* Returned point values. */
/* Note that source == input point, destination == output point */
struct _nearsmth {

  /* Public: */
	double sv[3];		/* Source value (input) */
	double sr;			/* Source radius */
	double drv[3];		/* Destination radial value (starting point & reference) */
	double dr;			/* Destination radial value radius */
	double sdv[3];		/* Smoothed destination value (output) */

 /* Private to nearsmth: */
	double _sv[3];		/* Original (non elevated) Source value (input) */
	double _sr;			/* Original radius */
	double _sdv[3];		/* Raw optimised destination value (not fixed for swap offset) */
	double sdve;		/* Smoothed destination error */
	gamut *sgam;		/* Source gamut */
	gamut *dgam;		/* Destination gamut */
	double nsv[3];		/* Normalised to rad. 1.0 source value (used for neighbor finding) */
	struct _nearsmth *n[NNB+1];	/* List of pointers to neigbors, NULL terminated */
	int swap;			/* src & dst are swapped */
	int debug;
}; typedef struct _nearsmth nearsmth;


/* Structures to hold weightings */

/* Color hextant */
typedef enum {
	gmm_end     = 0,		/* Mark the end of the list */
	gmm_ignore  = 1,		/* Ignore this entry */
	gmm_red     = 2,
	gmm_yellow  = 3,
	gmm_green   = 4,
	gmm_cyan    = 5,
	gmm_blue    = 6,
	gmm_magenta = 7,
	gmm_default = 8			/* Default parameter set */
} gmm_chex;

/* Single weight */
typedef struct {
	double l,c,h;
} iweight;

/* Group of complete weights */
/* (Remeber to alter near_wblend() !!) */
typedef struct {

	gmm_chex ch;			/* Color hextant this applies to */

	/* Absolute weighting */
	struct {
		double o;			/* Overall weight */
		iweight w;			/* Component weights */
	} a;

	/* Relative weighting */
	struct {
		double o;			/* Overall weight */
		iweight w;			/* Component weights */
	} r;

	/* radial weighting */
	struct {
		double o;			/* Overall weight */
		iweight w;			/* Component weights */
	} l;

	/* Cusp alignment control */
	double cw;				/* Cusp alignment weighting, 0 = none, 1 = full */
	double cr;				/* Cusp alignment radius minumum */

	/* Additional chrominance elevation */
	double e;

} gammapweights;


/* Blend a two groups of weights into one, given two weightings */
void near_wblend(gammapweights *dst,
	             gammapweights *src1, double wgt1, gammapweights *src2, double wgt2);

/* Return a list of points. Free list after use */
/* Return NULL on error */
nearsmth *near_smooth(
	int verb,			/* Verbose flag */
	int *npp,			/* Return the number of points returned */
	gamut *sc_gam,		/* Source colorspace gamut */
	gamut *s_gam,		/* Source image gamut (== sc_gam if none) */
	gamut *d_gam,		/* Destination colorspace gamut */
	double d_bp[3],		/* Destination target black point - may be NULL */
	gammapweights *wh,  /* Structure holding weights */
	int   usecomp,		/* Flag indicating whether smoothed compressed value will be used */
	int   useexp		/* Flag indicating whether smoothed expanded value will be used */
);

/* Expand the compact form of weights into the explicit form. */
/* The explicit form is red, yellow, green, cyan, blue, magenta, default */
void expand_weights(gammapweights out[7], gammapweights *in);
	
/* Blend a two expanded groups of individual weights into one */
void near_xwblend(
gammapweights *dst,
gammapweights *src1, double wgt1,
gammapweights *src2, double wgt2
);

#endif /* NEARSMTH_H */

