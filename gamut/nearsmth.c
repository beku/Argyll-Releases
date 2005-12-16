
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

/*

	Algorithm description:

	We want to return a list of points lying on the destination gamut
	surface, that are as near as possible to the grid points defining
	the source gamut surface, with the constraint that they are not
	too close together compared to the distance apart they would
	be if they were simply the radial mapping of the source point
	onto the destination gamut.

	In this way we hope to get a "smoothed" version of the
	nearest points. It is hoped that this will prevent
	"dead" areas in the final gamut mapping caused by many source
	colors being mapped onto a single nearest destination color,
	while still ensuring that most saturated source colors are
	mapped to most saturated destination colors.

	Initialy we start with destination mapping points that
	are simply the radially mapped source points. We then
	run an optimisation process on each point in turn,
	moving it to balance its closeness to its ideal "closest target point",
	against what it does to change the distance to its neighbor points.

	In this way we hope to allow the (initially) radial mapped
	points to be "squeezed" together somewhat, to move them
	closer to the nearest mapped target points, without squeezing
	them so much that they effectively become a many to one
	mapping.

 */

/*
 * TTBD:
 *       Make Abs L error weight inv. prop. to Abs C ??
 *       Improve error handling.
 *       Merge with gammap ?
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "gamut.h"
#include "numlib.h"
#include "nearsmth.h"

static double wde(double in1[3], double in2[3], double lweight, double cweight, double hweight);

struct _smthopt {
	nearsmth *p;			/* Point being optimised */
	gamut *d_gam;			/* Destination gamut surface */
	double a_lweight;		/* Absolute luminance error weight */
	double a_cweight0;		/* Absolute chroma error weight at L = 100 */
	double a_cweight1;		/* Absolute chroma error weight at L = 0 */
	double a_hweight;		/* Absolute hue error weight */
	double r_lweight;		/* Relative luminance error weight */
	double r_cweight;		/* Relative chroma error weight */
	double r_hweight;		/* Relative hue error weight */
	double d_lweight;		/* Radial luminance error weight */
	double d_cweight;		/* Radial chroma error weight */
	double d_hweight;		/* Radial hue error weight */
	int debug;				/* debug flag */
}; typedef struct _smthopt smthopt;

/* Powell optimisation function */
static double optfunc(
void *fdata,
double *_tp
) {
	smthopt *s = (smthopt *)fdata;
	nearsmth *p = s->p;
	double rv;
	int i, j;
	double dtp[3];		/* Point in question mapped to dst surface */
	double anv[3];		/* Average neighborhood target point */
	double a_cw;

	s->d_gam->radial(s->d_gam, dtp, _tp);	/* Map to dst surface to check current location */

	/* Compute absolute delta E squared to the source point */
#define LMIN 50.0
	a_cw = p->sv[0] > LMIN ? (p->sv[0] - LMIN)/(100.0 - LMIN) : 0.0;
	a_cw = (1.0 - a_cw) * s->a_cweight0 + a_cw * s->a_cweight1,
//printf("~1 L = %f, a_cw = %f\n",p->sv[0],a_cw);
	rv = wde(dtp, p->sv, s->a_lweight, a_cw, s->a_hweight);

	/* Compute delta E squared to relative target point */
	anv[0] = anv[1] = anv[2] =  0.0;
	for (i = 0; i < MAXGAMN; i++) {
		nearsmth *np = p->n[i];		/* Pointer to neighbor */

		if (np == NULL)
			break;

		anv[0] += np->sdv[0] - p->nv[i][0];
		anv[1] += np->sdv[1] - p->nv[i][1];
		anv[2] += np->sdv[2] - p->nv[i][2];
	}
	anv[0] /= (double)i;
	anv[1] /= (double)i;
	anv[2] /= (double)i;
	rv += wde(dtp, anv, s->r_lweight, s->r_cweight, s->r_hweight);

	/* Compute delta E squared to radial mapping of the source point */
	rv += wde(dtp, p->drv, s->d_lweight, s->d_cweight, s->d_hweight);

	if (s->debug)
		printf("debug: rv = %f from %f %f %f\n",rv, dtp[0], dtp[1], dtp[2]);

	return rv;
}

/* Compute the weighted delta E squared */
static double wde(
double in1[3],
double in2[3],
double lweight,
double cweight,
double hweight
) {
	double desq, dhsq;
	double dlsq, dcsq;

	/* Compute delta L squared and delta E squared */
	{
		double dl, da, db;
		dl = in1[0] - in2[0];
		dlsq = dl * dl;		/* dl squared */
		da = in1[1] - in2[1];
		db = in1[2] - in2[2];

		desq = dlsq + da * da + db * db;
	}

	/* compute delta chromanance squared */
	{
		double c1, c2, dc;

		/* Compute chromanance for the two colors */
		c1 = sqrt(in1[1] * in1[1] + in1[2] * in1[2]);
		c2 = sqrt(in2[1] * in2[1] + in2[2] * in2[2]);

		dc = c2 - c1;
		dcsq = dc * dc;

		/* ~~ Doing dcsq = sqrt(dcsq); here seemed */
		/* to improve the saturation result. How did this work ?? */
	}

	/* Compute delta hue squared */
	if ((dhsq = desq - dlsq - dcsq) < 0.0)
		dhsq = 0.0;


	return lweight * dlsq + cweight * dcsq + hweight * dhsq;
}

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
double a_cweight0,	/* Absolute chroma error weight at L = 100  */
double a_cweight1,	/* Absolute chroma error weight at L = 0  */
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
					/* Return sdv == sv if flag is 0 */
) {
	int nspts;		/* Number of source gamut points */
	nearsmth *smp;	/* Absolute delta E weighting */
	int ix, i, j, k;
	smthopt opts;	/* optimisation structure */
	int pass;
	double mxmv;						/* Maximum a point gets moved */
	
	nspts = sc_gam->nverts(sc_gam);		/* Source surface points */

//printf("~1 No src points = %d\n",nspts);
	if ((smp = (nearsmth *)malloc(nspts * sizeof(nearsmth))) == NULL) { 
		fprintf(stderr,"gamut map: Malloc of near smooth points failed\n");
		*npp = 0;
		return NULL;
	}

	opts.d_gam = d_gam;
	opts.a_lweight = a_weight * a_lweight;
	opts.a_cweight0 = a_weight * a_cweight0;
	opts.a_cweight1 = a_weight * a_cweight1;
	opts.a_hweight = a_weight * a_hweight;
	opts.r_lweight = r_weight * r_lweight;
	opts.r_cweight = r_weight * r_cweight;
	opts.r_hweight = r_weight * r_hweight;
	opts.d_lweight = d_weight * d_lweight;
	opts.d_cweight = d_weight * d_cweight;
	opts.d_hweight = d_weight * d_hweight;
	opts.debug = 0;

	/* First pass fills in the original point values */
	for (ix = i = 0; i < nspts; i++) {
		double csr, csv[3];		/* Colorspace gamut source point */
		double imr, imv[3];		/* Image gamut source point */
		int nix[MAXGAMN+1];		/* Neighbor indexes */

		ix = sc_gam->getvertn(sc_gam, nix, &csr, csv, ix); /* Get the vertex value */

		/* Lookup equivalent point on image gamut */
		if (s_gam != sc_gam)
			imr = s_gam->radial(s_gam, imv, csv);

		if (s_gam == sc_gam || imr > csr) {	/* No image gamut or strange image gamut */
			imr = csr;
			imv[0] = csv[0];
			imv[1] = csv[1];
			imv[2] = csv[2];
		}

		/* Lookup radialy equivalent point on destination gamut */
		smp[i].dr = d_gam->radial(d_gam, smp[i].drv, csv);

		/* Default setup a non-mapping of destination point to destination point */
		smp[i].sr = smp[i].dr;
		smp[i].sdv[0] = smp[i].sv[0] = smp[i].drv[0];
		smp[i].sdv[1] = smp[i].sv[1] = smp[i].drv[1];
		smp[i].sdv[2] = smp[i].sv[2] = smp[i].drv[2];

		/* Decide what we will use as our source point (sv) we use to relate */
		/* our output (sdv) to, in terms of absolute and relative error */
		if (usecomp && imr > smp[i].dr) {	/* Image is outside destination */
			smp[i].sr = imr;			/* We need to compress image to destination */
			smp[i].sv[0] = imv[0];
			smp[i].sv[1] = imv[1];
			smp[i].sv[2] = imv[2];
		} else if (useexp && smp[i].dr > csr) { /* Dest is outside Image & src colorspace */
			smp[i].sr = csr;			/* We need to expand colorspace to destination */
			smp[i].sv[0] = csv[0];
			smp[i].sv[1] = csv[1];
			smp[i].sv[2] = csv[2];
		}

		/* Copy the pointers to the neighbors and save neighbor vectors */
		for (j = 0; j < MAXGAMN; j++) {
			if (nix[j] < 0)
				break;
			smp[i].n[j] = &smp[nix[j]];
		}
		smp[i].n[j] = NULL;
	}

	/* If we can skip any optimisation */
	if (usecomp == 0 && useexp == 0) {
		*npp = nspts;
		return smp;
	}

	/* We need to generate the smoothed destination values */

	/* Second part saves the original (relative weighted) vectors to neighbours */
	for (i = 0; i < nspts; i++) {
		
		for (j = 0; j < MAXGAMN; j++) {
			double dd;

			if (smp[i].n[j] == NULL)
				break;
			
			for (k = 0; k < 3; k++) 
				smp[i].nv[j][k] = smp[i].n[j]->sv[k] - smp[i].sv[k];
		}
	}

	/* Now we want to optimise the location of the smoothed points */
	mxmv = 1e6;
	for (pass = 0; pass < 500 && mxmv > 0.01; pass++) {	/* Until we have converged */
		double s[3] = { 1.0, 1.0, 1.0 };	/* search area */
		double pv[3];						/* Previous value */

		mxmv = 0.0;
		for (i = 0; i < nspts; i++) {		/* Move all the points */
			double err, mv;

//printf("~1 moving point %f %f %f\n",smp[i].sdv[0], smp[i].sdv[1], smp[i].sdv[2]);

			for (j = 0; j < 3; j++)
				pv[j] = smp[i].sdv[j];

			opts.p = &smp[i];	/* Point to optimise */

			/* Optimise the point */
			err = powell(3, smp[i].sdv, s, 0.1, 100, optfunc, (void *)(&opts));

			if (err == -1.0) {
fprintf(stderr,"~1 powell failed in nearsmth()\n");
				opts.debug = 1;
				/* Optimise the point */
				err = powell(3, smp[i].sdv, s, 0.1, 100, optfunc, (void *)(&opts));

				free(smp);
				*npp = 0;
				return NULL;
			}

			/* Remap it to the destinaton gamut surface */
			d_gam->radial(d_gam, smp[i].sdv, smp[i].sdv);

			/* See how much it moved */
			mv = 0.0;
			for (j = 0; j < 3; j++) {
				double tt = pv[j] - smp[i].sdv[j];
				mv += tt * tt;
			}
			if (mv > mxmv)
				mxmv = mv;

//printf("~1     to %f %f %f\n",smp[i].sdv[0], smp[i].sdv[1], smp[i].sdv[2]);
		}
		mxmv = sqrt(mxmv);
		if (verb) {
			printf("."); fflush(stdout);
		}

//printf("~1 mxmv = %f\n",mxmv);
	}

	/* Apply the given overshoot to mapping vectors */
	if (overshoot != 1.0) {

		for (ix = i = 0; i < nspts; i++) {
			double rat = overshoot - 1.0;
			for (j = 0; j < 3; j++) {
				smp[i].sdv[j] += rat * (smp[i].sdv[j] - smp[i].sv[j]);
			}
		}
	}
	if (verb)
		printf("\n");

	*npp = nspts;
	return smp;
}

	


















