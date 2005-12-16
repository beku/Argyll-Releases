
/* 
 * Argyll Color Correction System
 *
 * Perceptual space random test point class
 *
 * Author: Graeme W. Gill
 * Date:   12/9/2004
 *
 * Copyright 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */


/* TTBD:

 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#ifdef DEBUG
#include "plot.h"
#endif
#include "numlib.h"
#include "sort.h"
#include "plot.h"
#include "icc.h"
#include "xcolorants.h"
#include "targen.h"
#include "prand.h"

/* ----------------------------------------------------- */

/* Default convert the nodes device coordinates into approximate perceptual coordinates */
/* (usually overriden by caller supplied function) */
static void
default_prand(void *od, double *p, double *d) {
	prand *s = (prand *)od;
	int e;
	double tt;

	/* Default Do nothing - copy device to perceptual. */
	for (e = 0; e < s->di; e++) {
		p[e] = d[e] * 100.0;
	}
}

/* Return the largest distance of the point outside the device gamut. */
/* This will be 0 if inside the gamut, and > 0 if outside.  */
static double
prand_in_dev_gamut(prand *s, double *d) {
	int e;
	int di = s->di;
	double tt, dd = 0.0;
	double ss = 0.0;

	for (e = 0; e < di; e++) {
		ss += d[e];

		tt = 0.0 - d[e];
		if (tt > 0.0) {
			if (tt > dd)
				dd = tt;
		}
		tt = d[e] - 1.0; 
		if (tt > 0.0) {
			if (tt > dd)
				dd = tt;
		}
	}
	tt = ss - s->ilimit;
	if (tt > 0.0) {
		if (tt > dd)
			dd = tt;
	}
	return dd;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Reverse lookup function :- perceptual to device coordinates */

/* Structure to hold data for optimization function */
struct _edatas {
	prand *s;			/* prand structure */
	double *ptp;		/* Perceptual target point */
	}; typedef struct _edatas edatas;

/* Definition of the optimization functions handed to powell() */

/* This one returns error from perceptual target point, and */
/* an error >= 50000 on being out of device gamut */
static double efunc(void *edata, double p[]) {
	edatas *ed = (edatas *)edata;
	prand *s = ed->s;
	int e, di = s->di;
	double rv, pp[MXTD];
	if ((rv = (prand_in_dev_gamut(s, p))) > 0.0) {
		rv = rv * 5000.0 + 100000.0;		/* Discourage being out of gamut */
	} else {
		s->percept(s->od, pp, p);
		for (rv = 0.0, e = 0; e < di; e++) {
			double tt = pp[e] - ed->ptp[e];
			rv += tt * tt;
		}
	}
//printf("rv = %f from %f %f\n",rv,p[0],p[1]);
	return rv;
}

/* Given a point in perceptual space, an approximate point */
/* in device space, return the device value corresponding to */
/* the perceptual value, plus the clipped perceptual value. */
/* Return 1 if the point is out of gamut. */
static int
prand_from_percept(
prand *s,
double *d,			/* return device position */
double *p			/* given perceptual value, return clipped. */
) {
	int e, di = s->di;
	edatas ed;
	double pp[MXTD];	/* Clipped perceptual */
	double sr[MXTD];	/* Search radius */
	double tt;
	double drad = 50.0;	/* Search radius */
	double ptol = 0.0001;	/* Tolerance */
	ed.s = s;
	ed.ptp = p;			/* Set target perceptual point */

	for (e = 0; e < di; e++) {
		sr[e] = drad;			/* Device space search radius */
	}
	if ((tt = powell(di, d, sr,  ptol, 500, efunc, (void *)&ed)) < 0.0 || tt >= 50000.0) {
		error("prand: powell failed, tt = %f\n",tt);
	}
	s->percept(s->od, pp, d);	/* Lookup clipped perceptual */
	
	/* Compute error between target and possibly clipped perceptual */
	tt = 0.0;
	for (e = 0; e < di; e++) {
		double t = p[e] - pp[e];
		p[e] = pp[e];
		tt += t * t;
	}
	tt = sqrt(tt);
//printf("~1 perc %f %f -> %f %f dev %f %f, err = %f\n",ed.ptp[0],ed.ptp[1],p[0],p[1],d[0],d[1],tt);
	if (tt > 1.0)
		return 1;		/* More than 1 delta E */
	return 0;
}

/* --------------------------------------------------- */
/* Seed the object with the initial fixed points */

static void
prand_add_fixed(
prand *s,
fxpos *fxlist,			/* List of existing fixed points */
int fxno				/* Number in fixed list */
) {
	int e, di = s->di;
	int i, j;

	/* Add fixed points if there are any */
	if (fxno > 0) {

		for (i = 0; (i < fxno) && (i < s->tinp); i++) {
			prnode *p = &s->n[i];	/* Destination for point */

			for (e = 0; e < di; e++)
				p->p[e] = fxlist[i].p[e];

			p->fx = 1;			/* is a fixed point */
			s->percept(s->od, p->v, p->p);
			s->np = s->fnp = i+1;
		}
	}
}

/* Seed the object with the movable incremental farthesr points. */
static void
prand_seed(prand *s) {
	int e, di = s->di;

	printf("\n");

	/* Seed the non-fixed points */
	for (; s->np < s->tinp;) {
		prnode *p = &s->n[s->np];		/* Next node */

		for (e = 0; e < di; e++) {
			if (e == 1 || e == 2)
				p->v[e] = d_rand(-128.0, 128.0);
			else
				p->v[e] = d_rand(0.0, 100.0);
		}
		if (prand_from_percept(s, p->p, p->v) == 0) {
			s->np++;
			printf("\rAdded %d/%d",s->np,s->tinp); fflush(stdout);
		}
	}
	printf("\n");
}

/* --------------------------------------------------- */
/* Support accessing the list of generated sample points */

/* Rest the read index */
static void
prand_reset(prand *s) {
	s->rix = 0;
}

/* Read the next non-fixed point value */
/* Return nz if no more */
static int
prand_read(prand *s, double *p, double *f) {
	int e, di = s->di;

	/* Advance to next non-fixed point */
	while(s->rix < s->np && s->n[s->rix].fx)
		s->rix++;
	
	if (s->rix >= s->np)
		return 1;

	/* Return point info to caller */
	for (e = 0; e < s->di; e++) {
		if (p != NULL)
			p[e] = s->n[s->rix].p[e];
		if (f != NULL)
			f[e] = s->n[s->rix].v[e];
	}
	s->rix++;

	return 0;
}

/* --------------------------------------------------- */
/* Main object creation/destruction */

/* Destroy ourselves */
static void
prand_del(prand *s) {
	int i;

	free(s->n);

	free (s);
}

/* Creator */
prand *new_prand(
int di,					/* Dimensionality of device space */
double ilimit,			/* Ink limit (sum of device coords max) */
int tinp,				/* Total number of points to generate, including fixed */
fxpos *fxlist,			/* List of existing fixed points (may be NULL) */
int fxno,				/* Number of existing fixes points */
void (*percept)(void *od, double *out, double *in),		/* Perceptual lookup func. */
void *od				/* context for Perceptual function */
) {
	int i;
	prand *s;

	if ((s = (prand *)calloc(sizeof(prand), 1)) == NULL)
		error ("prand: malloc failed");

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	s->di = di;

	if (tinp < fxno)	/* Make sure we return at least the fixed points */
		tinp = fxno;

	s->tinp = tinp;		/* Target total number of points */
	s->ilimit = ilimit;

	/* Init method pointers */
	s->reset = prand_reset;
	s->read  = prand_read;
	s->del   = prand_del;

	/* If no perceptual function given, use default */
	if (percept == NULL) {
		s->percept = default_prand;
		s->od = s;
	} else {
		s->percept = percept;
		s->od = od;
	}

	/* Allocate the space for the target number of points */
	if ((s->n = (prnode *)calloc(sizeof(prnode), s->tinp)) == NULL)
		error ("prand: malloc failed on sample nodes");
	s->np = s->fnp = 0;

	/* Setup the fixed points */
	prand_add_fixed(s, fxlist, fxno);

	if (tinp > fxno) { /* Create the perceptual space random points */
		prand_seed(s);
	}

	prand_reset(s);		/* Reset read index */

	return s;
}

/* =================================================== */
















