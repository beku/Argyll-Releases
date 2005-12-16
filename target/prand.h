
#ifndef PRAND_H

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

/* A sample point node */
struct _prnode {
	int fx;				/* nz if point is fixed */
	double p[MXTD];		/* Device coordinate position */
	double v[MXTD];		/* Subjective value (Labk) */
}; typedef struct _prnode prnode;


/* Main simplex latice object */
struct _prand {
/* private: */
	int di;			/* Point dimensionality */
	double ilimit;	/* Ink limit - limit on sum of p[] */
	int fnp;		/* Number of existing fixed points in list */
	int tinp;		/* target number of total points in list, including fixed points */

	int np;			/* Number of points currently in list */
	prnode *n;		/* tinp list of points */

	/* Perceptual function handed in */
	void (*percept)(void *od, double *out, double *in);
	void *od;		/* Opaque data for perceptual point */
	
	/* Other info */
	int rix;			/* Next read index */

/* public: */
	/* Initialise, ready to read out all the points */
	void (*reset)(struct _prand *s);

	/* Read the next set of non-fixed points values */
	/* return non-zero when no more points */
	int (*read)(struct _prand *s, double *d, double *p);

	/* Destroy ourselves */
	void (*del)(struct _prand *s);

	}; typedef struct _prand prand;

/* Constructor */
extern prand *new_prand(int di, double ilimit, int npoints,
	fxpos *fxlist, int fxno,
	void (*percept)(void *od, double *out, double *in), void *od);

#define PRAND_H
#endif /* PRAND_H */