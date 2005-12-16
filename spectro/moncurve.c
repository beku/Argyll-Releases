
/* 
 * Argyll Color Correction System
 * Monotonic curve class for display calibration.
 *
 * Author: Graeme W. Gill
 * Date:   30/10/2005
 *
 * Copyright 2005 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the LICENCE.TXT file for licencing details.
 *
 * This is based on the monotonic curve equations used elsewhere,
 * but currently intended to support the display calibration process.
 * moncurve is not currently general, missing:
 *
 *  input scaling
 *  output scaling
 */

#undef DIAG

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#if defined(__IBMC__) && defined(_M_IX86)
#include <float.h>
#endif
#include "numlib.h"
#include "moncurve.h"

#define POWTOL 1e-4
#define MAXITS 10000

/* SHAPE_BASE doesn't seem extremely critical.  It is centered in +/- 1 magnitude */
/* 10 x more filters out noise reasonably heaviliy, 10 x less gives noticable */
/* overshoot Range 0.00001 .. 0.001 */
/* SHAPE_HBASE is more critical. */
/* Range 0.00005 .. 0.001 */
#define SHAPE_BASE  0.00001		/* 0 & 1 harmonic weight */
#define SHAPE_HBASE 0.00001		/* 2nd and higher additional weight */

static void mcv_del(mcv *p);
static void mcv_fit(mcv *p, int verb, int order, mcvco *d, int ndp, double smooth);
static void mcv_force_scale(mcv *p, double target);
static double mcv_interp(struct _mcv *p, double in);
double mcv_dinterp(mcv *p, double *dv, double vv);
static double mcv_inv_interp(struct _mcv *p, double in);

/* Allocate a new, uninitialised mcv */
/* Note thate black and white points aren't allocated */
mcv *new_mcv(void) {
	mcv *p;

	if ((p = (mcv *)calloc(1, sizeof(mcv))) == NULL)
		return NULL;

	/* Init method pointers */
	p->del         = mcv_del;
	p->fit         = mcv_fit;
	p->force_scale = mcv_force_scale;
	p->interp      = mcv_interp;
	p->inv_interp  = mcv_inv_interp;

	p->luord = 0;
	p->pms = NULL;
	p->upms = NULL;
	return p;
}

/* Delete an mcv */
static void mcv_del(mcv *p) {
	if (p->pms != NULL)
		free(p->pms);
	free(p);
}

/* Shaper+Matrix optimisation function handed to powell() */
static double mcv_opt_func(void *edata, double *v) {
	mcv *p = (mcv *)edata;
	double rv = 0.0, smv;
	double out;
	int i;

	p->upms = v;

	/* For all our data points */
	for (i = 0; i < p->ndp; i++) {
		double ev;
		int j;

		/* Apply our function */
		out = p->interp(p, p->d[i].p);

		ev = out - p->d[i].v;

		rv += ev * ev;
	}
	p->upms = p->pms;

	/* Normalise error to be an average delta E squared */
	rv /= (double)p->ndp;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	smv = 0.0;
	for (i = 2; i < p->luord; i++) {
		double w, tt = v[i];
		tt = v[i];
		tt *= tt;

		/* Weigh to supress ripples */
		if (i <= 3) {	/* First or second curves */
			w = SHAPE_BASE;
		} else {
			w = SHAPE_BASE + (i-3) * SHAPE_HBASE * p->smooth;
		}
		smv += w * tt;
	}
	rv += smv;

//printf("~1 rv = %f (%f)\n",rv,smv);
	return rv;
}

/* Shaper+Matrix optimisation function handed to conjgrad() */
static double mcv_dopt_func(void *edata, double *dv, double *v) {
	mcv *p = (mcv *)edata;
	double rv = 0.0, smv;
	double out;
	int i, j;

	/* Zero the dv's */
	for (j = 0; j < p->luord; j++)
		dv[j] = 0.0;

	p->upms = v;

	/* For all our data points */
	for (i = 0; i < p->ndp; i++) {
		double ev;

		/* Apply our function with dv's */
		out = mcv_dinterp(p, p->dv, p->d[i].p);

		ev = out - p->d[i].v;
		rv += ev * ev;

		/* Sum the dv's */
		for (j = 0; j < p->luord; j++)
			dv[j] += 2.0 * ev * p->dv[j];

	}
	p->upms = p->pms;

	/* Normalise error to be an average delta E squared */
	rv /= (double)p->ndp;
	for (j = 0; j < p->luord; j++)
		dv[j] /= (double)p->ndp;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	smv = 0.0;
	for (i = 2; i < p->luord; i++) {
		double w, tt = v[i];
		tt = v[i];

		/* Weigh to supress ripples */
		if (i <= 3) {	/* First or second curves */
			w = SHAPE_BASE;
		} else {
			w = SHAPE_BASE + (i-3) * SHAPE_HBASE * p->smooth;
		}
		dv[i] += 2.0 * w * tt;
		tt *= tt;
		smv += w * tt;
	}
	rv += smv;

//printf("~1 rv = %f (%f)\n",rv,smv);
	return rv;
}

#ifdef NEVER
/* Check partial derivative function */

static double mcv_opt_func(void *edata, double *v) {
	mcv *p = (mcv *)edata;
	int i;
	double dv[200];
	double rv, drv;
	double trv;
	
	rv = mcv_opt_func_(edata, v);
	drv = mcv_dopt_func(edata, dv, v);

	if (fabs(rv - drv) > 1e-6)
		printf("######## RV MISMATCH is %f should be %f ########\n",rv,drv);

	/* Check each parameter delta */
	for (i = 0; i < p->luord; i++) {
		double del;

		v[i] += 1e-7;
		trv = mcv_opt_func_(edata, v);
		v[i] -= 1e-7;
		
		/* Check that del is correct */
		del = (trv - rv)/1e-7;
		if (fabs(dv[i] - del) > 0.04) {
//printf("~1 del = %f from (trv %f - rv %f)/0.1\n",del,trv,rv);
			printf("######## EXCESSIVE at v[%d] is %f should be %f ########\n",i,dv[i],del);
		}
	}
	return rv;
}
#endif

/* Fit the curve to the given points */
static void mcv_fit(mcv *p,
	int verb,		/* Vebosity level, 0 = none */
	int order,		/* Number of curve orders, 1..MCV_MAXORDER */
	mcvco *d,		/* Array holding scattered initialisation data */
	int ndp,		/* Number of data points */
	double smooth	/* Degree of smoothing, 1.0 = normal */			
) {
	int i;
	double *sa;		/* Search area */
	double *pms;	/* Parameters to optimise */

	p->verb = verb;
	p->smooth = smooth;
	p->luord = order+2;		/* Add two for offset and scale */

	if (p->pms != NULL)
		free(p->pms);
	if ((p->pms = (double *)calloc(p->luord, sizeof(double))) == NULL)
		error ("Malloc failed");
	if ((pms = (double *)calloc(p->luord, sizeof(double))) == NULL)
		error ("Malloc failed");
	if ((sa = (double *)calloc(p->luord, sizeof(double))) == NULL)
		error ("Malloc failed");
	if ((p->dv = (double *)calloc(p->luord, sizeof(double))) == NULL)
		error ("Malloc failed");
	p->upms = p->pms;

	/* Set offset and scale to reasonable values */
	p->pms[0] = 1e38;
	p->pms[1] = -1e38;
	for (i = 0; i < ndp; i++) {
		if (d[i].v < p->pms[0])
			p->pms[0] = d[i].v;
		if (d[i].v > p->pms[1])
			p->pms[1] = d[i].v;
	}
	p->pms[1] -= p->pms[0];

	/* Use powell to minimise the sum of the squares of the */
	/* input points to the curvem, plus a parameter damping factor. */
	p->d = d;
	p->ndp = ndp;

	for (i = 0; i < p->luord; i++)
		sa[i] = 0.2;

#ifdef NEVER
	if (powell(p->luord, p->pms, sa, POWTOL, MAXITS, mcv_opt_func, (void *)p) < 0.0)
		error ("Powell failed");
#else
	if (conjgrad(p->luord, p->pms, sa, POWTOL, MAXITS, mcv_opt_func, mcv_dopt_func, (void *)p) < 0.0)
		error ("Conjgrad failed");
#endif

	free(p->dv);
	p->dv = NULL;
	free(sa);
	free(pms);
}

/* Scale the the output so that the value for input 1.0, */
/* is the given target value. */
static void mcv_force_scale(mcv *p,
	double target	/* Target output value */
) {
	if (p->luord > 1) {
		target -= p->upms[0];	/* Offset */
		p->upms[1] = target;	/* Scale */
	}
}

/* Translate a value through the curve */
static double mcv_interp(struct _mcv *p,
	double vv	/* Input value */
) {
	double g;
	int ord;

	/* Process all the shaper orders from low to high. */
	/* [These shapers were inspired by a Graphics Gem idea */
	/* (Gems IV, VI.3, "Fast Alternatives to Perlin's Bias and */
	/*  Gain Functions, pp 401). */
	/*  They have the nice properties that they are smooth, and */
	/*  are monotonic. The control parameter has been */
	/*  altered to have a range from -oo to +oo rather than 0.0 to 1.0 */
	/*  so that the search space is less non-linear. */
	for (ord = 2; ord < p->luord; ord++) {
		int nsec;			/* Number of sections */
		double sec;			/* Section */

		g = p->upms[ord];	/* Parameter */

		nsec = ord-1;		/* Increase sections for each order */

		vv *= (double)nsec;

		sec = floor(vv);
		if (((int)sec) & 1)
			g = -g;			/* Alternate action in each section */
		vv -= sec;
		if (g >= 0.0) {
			vv = vv/(g - g * vv + 1.0);
		} else {
			vv = (vv - g * vv)/(1.0 - g * vv);
		}
		vv += sec;
		vv /= (double)nsec;
	}
	/* Do order 0 & 1 */
	if (p->luord > 1)
		vv *= p->upms[1];	/* Scale */

	if (p->luord > 0)
		vv += p->upms[0];	/* Offset */

	return vv;
}

/* Transfer function with partial derivative */
/* with respect to the parameters. */
double mcv_dinterp(mcv *p,
double *dv,			/* Return derivative wrt each parameter */
double vv			/* Source of value */
) {
	double g;
	int i, ord;

	/* Process all the shaper orders from low to high. */
	for (ord = 2; ord < p->luord; ord++) {
		double dsv;		/* del for del in g */
		double ddv;		/* del for del in vv */
		int nsec;		/* Number of sections */
		double sec;		/* Section */

		g = p->upms[ord];			/* Parameter */

		nsec = ord-1;	/* Increase sections for each order */

		vv *= (double)nsec;

		sec = floor(vv);
		if (((int)sec) & 1) {
			g = -g;				/* Alternate action in each section */
		}
		vv -= sec;
		if (g >= 0.0) {
			double tt = g - g * vv + 1.0;
			dsv = (vv * vv - vv)/(tt * tt);
			ddv = (g + 1.0)/(tt * tt);
			vv = vv/tt;
		} else {
			double tt = 1.0 - g * vv;
			dsv = (vv * vv - vv)/(tt * tt);
			ddv = (1.0 - g)/(tt * tt);
			vv = (vv - g * vv)/tt;
		}

		vv += sec;
		vv /= (double)nsec;
		dsv /= (double)nsec;
		if (((int)sec) & 1)
			dsv = -dsv;

		dv[ord] = dsv;
		for (i = ord - 1; i >= 2; i--)
			dv[i] *= ddv;
	}
	/* Do order 0, the scale */
	if (p->luord > 1) {
		dv[1] = vv;
		vv *= p->upms[1];
	}
	if (p->luord > 0) {
		dv[0] = 1.0;
		vv += p->upms[0];	/* Offset */
	}

	return vv;
}

/* Translate a value through backwards the curve */
static double mcv_inv_interp(struct _mcv *p,
	double vv	/* Input value */
) {
	double g;
	int ord;

	/* Process everything in reverse order to mcv_interp */

	/* Do order 0 & 1, the offset and scale */
	if (p->luord > 0)
		vv -= p->upms[0];

	if (p->luord > 1)
		vv /= p->upms[1];

	for (ord = p->luord-1; ord > 1; ord--) {
		int nsec;			/* Number of sections */
		double sec;			/* Section */

		g = -p->upms[ord];	/* Inverse parameter */

		nsec = ord-1;		/* Increase sections for each order */

		vv *= (double)nsec;

		sec = floor(vv);
		if (((int)sec) & 1)
			g = -g;			/* Alternate action in each section */
		vv -= sec;
		if (g >= 0.0) {
			vv = vv/(g - g * vv + 1.0);
		} else {
			vv = (vv - g * vv)/(1.0 - g * vv);
		}
		vv += sec;
		vv /= (double)nsec;
	}

	return vv;
}







