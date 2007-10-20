
/* 
 * Argyll Color Correction System
 * Monotonic curve class for display calibration.
 *
 * Author: Graeme W. Gill
 * Date:   30/10/2005
 *
 * Copyright 2005 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * This is based on the monotonic curve equations used elsewhere,
 * but currently intended to support the display calibration process.
 * moncurve is not currently general, missing:
 *
 *  input scaling
 *  output scaling
 */

#undef DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#if defined(__IBMC__) && defined(_M_IX86)
#include <float.h>
#endif
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "moncurve.h"

#define POWTOL 1e-4			/* Powell optimiser tollerance (was 1e-5 ?) */
#define MAXITS 10000

#undef TEST_PDE				/* Ckeck partial derivative calcs */

/* SHAPE_BASE doesn't seem extremely critical.  It is centered in +/- 1 magnitude */
/* 10 x more filters out noise reasonably heaviliy, 10 x less gives noticable */
/* overshoot Range 0.00001 .. 0.001 */
/* SHAPE_HBASE is more critical. */
/* Range 0.00005 .. 0.001 */
#define SHAPE_BASE  0.00001		/* 0 & 1 harmonic weight */
#define SHAPE_HBASE 0.00005		/* 2nd and higher additional weight */

static void mcv_del(mcv *p);
static void mcv_fit(mcv *p, int verb, int order, mcvco *d, int ndp, double smooth);
static void mcv_force_0(mcv *p, double target);
static void mcv_force_1(mcv *p, double target);
static void mcv_force_scale(mcv *p, double target);
static int mcv_get_params(mcv *p, double **rp);
static double mcv_interp(struct _mcv *p, double in);
static double mcv_inv_interp(struct _mcv *p, double in);

static double mcv_interp_p(struct _mcv *p, double *pms, double in);
static double mcv_shweight_p(mcv *p, double *v, double smooth);
double mcv_dinterp_p(mcv *p, double *pms, double *dv, double vv);
static double mcv_dshweight_p(mcv *p, double *v, double *dv, double smooth);

/* Allocate a new, uninitialised mcv */
/* Note thate black and white points aren't allocated */
mcv *new_mcv(void) {
	mcv *p;

	if ((p = (mcv *)calloc(1, sizeof(mcv))) == NULL)
		return NULL;

	/* Init method pointers */
	p->del         = mcv_del;
	p->fit         = mcv_fit;
	p->force_0     = mcv_force_0;
	p->force_1     = mcv_force_1;
	p->force_scale = mcv_force_scale;
	p->get_params  = mcv_get_params;
	p->interp      = mcv_interp;
	p->inv_interp  = mcv_inv_interp;
	p->interp_p    = mcv_interp_p;
	p->shweight_p  = mcv_shweight_p;
	p->dinterp_p   = mcv_dinterp_p;
	p->dshweight_p = mcv_dshweight_p;

	p->luord = 0;
	p->pms = NULL;
	return p;
}

/* Create a new mcv initiated with the given parameters */
mcv *new_mcv_p(double *pp, int np) {
	int i;
	mcv *p;

	if ((p = new_mcv()) == NULL)
		return p;

	p->luord = np;
	if ((p->pms = (double *)calloc(p->luord, sizeof(double))) == NULL)
		error ("Malloc failed");

	for (i = 0; i < np; i++)
		p->pms[i] = *pp++;

	return p;
}

/* Delete an mcv */
static void mcv_del(mcv *p) {
	if (p->pms != NULL)
		free(p->pms);
	free(p);
}

#ifdef TEST_PDE
#define mcv_opt_func mcv_opt_func_
#endif

/* Shaper+Matrix optimisation function handed to powell() */
static double mcv_opt_func(void *edata, double *v) {
	mcv *p = (mcv *)edata;
	double totw = 0.0;
	double rv = 0.0, smv;
	double out;
	int i;

#ifdef NEVER
	printf("params =");
	for (i = 0; i < p->luord; i++)
		printf(" %f",v[i]);
	printf("\n");
#endif

	/* For all our data points */
	for (i = 0; i < p->ndp; i++) {
		double ev;

		/* Apply our function */
		out = p->interp_p(p, v, p->d[i].p);

		ev = out - p->d[i].v;

		rv += p->d[i].w * ev * ev;
		totw += p->d[i].w;
	}

	/* Normalise error to be an average delta E squared */
	rv /= totw;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	rv += smv = mcv_shweight_p(p, v, p->smooth);

#ifdef NEVER
	printf("rv = %f (%f)\n",rv,smv);
#endif
	return rv;
}

/* Shaper+Matrix optimisation function handed to conjgrad() */
static double mcv_dopt_func(void *edata, double *dv, double *v) {
	mcv *p = (mcv *)edata;
	double totw = 0.0;
	double rv = 0.0, smv;
	double out;
	int i, j;

#ifdef NEVER
	printf("params =");
	for (i = 0; i < p->luord; i++)
		printf(" %f",v[i]);
	printf("\n");
#endif

	/* Zero the dv's */
	for (j = 0; j < p->luord; j++)
		dv[j] = 0.0;

	/* For all our data points */
	for (i = 0; i < p->ndp; i++) {
		double ev;

		/* Apply our function with dv's */
		out = p->dinterp_p(p, v, p->dv, p->d[i].p);

		ev = out - p->d[i].v;
		rv += p->d[i].w * ev * ev;

		/* Sum the dv's */
		for (j = 0; j < p->luord; j++)
			dv[j] += p->d[i].w * 2.0 * ev * p->dv[j];

		totw += p->d[i].w;
	}

	/* Normalise error to be an average delta E squared */
	rv /= totw;
	for (j = 0; j < p->luord; j++)
		dv[j] /= totw; 

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles", */
	/* with partial derivatives */
	rv += smv = mcv_dshweight_p(p, v, dv, p->smooth);

#ifdef NEVER
	printf("drv = %f (%f)\n",rv,smv);
#endif
	return rv;
}

#ifdef TEST_PDE
/* Check partial derivative function */

#undef mcv_opt_func

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
#endif	/* TEST_PDE */

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

	/* Set offset and scale to reasonable values */
	p->pms[0] = 1e38;			/* Locate min, and make that offset */
	p->pms[1] = -1e38;			/* Locate max */
	for (i = 0; i < ndp; i++) {
		if (d[i].v < p->pms[0])
			p->pms[0] = d[i].v;
		if (d[i].v > p->pms[1])
			p->pms[1] = d[i].v;
#ifdef DEBUG
		printf("point %d is %f %f\n",i,d[i].p,d[i].v);
#endif
	}
	p->pms[1] -= p->pms[0];		/* make max into scale */

	/* Use powell to minimise the sum of the squares of the */
	/* input points to the curvem, plus a parameter damping factor. */
	p->d = d;
	p->ndp = ndp;

	for (i = 0; i < p->luord; i++)
		sa[i] = 0.2;

#ifdef NEVER
	if (powell(NULL, p->luord, p->pms, sa, POWTOL, MAXITS, mcv_opt_func, (void *)p) != 0)
		error ("Mcv fit powell failed");
#else
	if (conjgrad(NULL, p->luord, p->pms, sa, POWTOL, MAXITS, mcv_opt_func, mcv_dopt_func, (void *)p) != 0)
		error ("Mcv fit conjgrad failed");
#endif

	free(p->dv);
	p->dv = NULL;
	free(sa);
	free(pms);
}

/* The native values from the curve parameters are 0 - 1.0, */
/* then the scale is applied, then the offset added, so the */
/* output always ranges from (offset) to (offset + scale). */

/* Offset the the output so that the value for input 0.0, */
/* is the given value. Don't change the output for 1.0 */
void mcv_force_0(
	mcv *p,
	double target	/* Target output value */
) {
	if (p->luord > 0) {
		target -= p->pms[0];	/* Change */
		if (p->luord > 1)
			p->pms[1] -= target;	/* Adjust scale to leave 1.0 output untouched */
		p->pms[0] += target;	/* Adjust offset */
	}
}

/* Scale the the output so that the value for input 1.0, */
/* is the given target value. Don't change the output for 0.0 */
static void mcv_force_1(
	mcv *p,
	double target	/* Target output value */
) {
	if (p->luord > 1) {
		target -= p->pms[0];	/* Offset */
		p->pms[1] = target;	/* Scale */
	}
}

/* Scale the the output so that the value for input 1.0, */
/* is the given target value. Scale the value for 0 in proportion. */
static void mcv_force_scale(
	mcv *p,
	double target	/* Target output value */
) {
	if (p->luord > 1) {
		p->pms[0] *= target/(p->pms[0] + p->pms[1]);	/* Offset */
		p->pms[1] = target - p->pms[0];	/* Scale */
	}
}

/* Return the number of parameters and the parameters in */
/* an allocated array. free() when done */
static int mcv_get_params(mcv *p, double **rp) {
	double *pp;
	int np, i;

	np = p->luord;

	if ((pp = (double *)malloc(np * sizeof(double))) == NULL)
		error("mcb_get_params malloc failed");
	
	*rp = pp;
	
	for (i = 0; i < np; i++)
		*pp++ = p->pms[i];
		
	return np;
}

/* Translate a value through the curve */
/* using the currently set pms */
static double mcv_interp(struct _mcv *p,
	double vv	/* Input value */
) {
	return mcv_interp_p(p, p->pms, vv);
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
		vv -= p->pms[0];

	if (p->luord > 1)
		vv /= p->pms[1];

	for (ord = p->luord-1; ord > 1; ord--) {
		int nsec;			/* Number of sections */
		double sec;			/* Section */

		g = -p->pms[ord];	/* Inverse parameter */

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

/* Translate a value through the curve */
/* using the given parameters */
static double mcv_interp_p(
	mcv *p,
	double *pms,	/* Parameters to use */
	double vv		/* Input value */
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

		g = pms[ord];	/* Parameter */

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
		vv *= pms[1];	/* Scale */

	if (p->luord > 0)
		vv += pms[0];	/* Offset */

	return vv;
}

/* Return the shaper parameters regularizing weight */
static double mcv_shweight_p(mcv *p, double *v, double smooth) {
	double smv;
	int i;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	/* Note:- we start at 2, to skip offset and scale. */
	smv = 0.0;
	for (i = 2; i < p->luord; i++) {
		double w, tt = v[i];
		tt = v[i];
		tt *= tt;

		/* Weigh to supress ripples */
		if (i <= 3) {	/* First or second curves */
			w = SHAPE_BASE;
		} else {
			w = SHAPE_BASE + (i-3) * SHAPE_HBASE * smooth;
		}
		smv += w * tt;
	}
	return smv;
}

/* Transfer function with partial derivative */
/* with respect to the given parameters. */
double mcv_dinterp_p(mcv *p,
double *pms,		/* Parameters to use */
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

		g = pms[ord];			/* Parameter */

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
		vv *= pms[1];
	}
	if (p->luord > 0) {
		dv[0] = 1.0;
		vv += pms[0];	/* Offset */
	}

	return vv;
}

/* Return the shaper parameters regularizing weight, */
/* and add in partial derivatives. */
/* Weight error and derivatrive by wht */
static double mcv_dshweight_p(mcv *p, double *v, double *dv, double smooth) {
	double smv;
	int i;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles", */
	/* with partial derivatives */
	smv = 0.0;
	for (i = 2; i < p->luord; i++) {
		double w, tt = v[i];
		tt = v[i];

		/* Weigh to supress ripples */
		if (i <= 3) {	/* First or second curves */
			w = SHAPE_BASE;
		} else {
			w = SHAPE_BASE + (i-3) * SHAPE_HBASE * smooth;
		}
		dv[i] += 2.0 * w * tt;
		tt *= tt;
		smv += w * tt;
	}

	return smv;
}




