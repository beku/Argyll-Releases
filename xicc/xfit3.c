
/*
 * Clut per channel curve fitting - algorithm 3
 *
 * Author:  Graeme W. Gill
 * Date:    27/5/2007
 * Version: 1.00
 *
 * Copyright 2000 - 2007 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * Based on the rspl.c and xlut.c code.
 */

/*
 * This module provides curve and matrix fitting functionality used to
 * create per channel input and/or output curves for clut profiles.
 *
 * ~~~99 when finished ~~~~
 * This version creates a set of directional rspl's and then uses
 * these as idealized "perfectly well known" device behaviour
 * to create the input curves, then the real resolution rspl
 * is created, and the input scattered points are then used with
 * it to create the sub-grid and output curves.
 */

#include <sys/types.h>
#include <string.h>
#include <ctype.h>
#ifdef __sun
#include <unistd.h>
#endif
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "icc.h"
#include "rspl.h"
#include "xicc.h"
#include "plot.h"


#define USE_CIE94_DE	/* Use CIE94 delta E measure when creating in/out curves */

#undef DEBUG 			/* Verbose debug information */
#undef DEBUG_PLOT 		/* Plot in & out curves */

#define POWTOL 5e-5			/* Shaper Powell optimiser tollerance */
#define MAXITS 2000			/* Shaper number of itterations before giving up */

/* Weights in shaper parameters, to minimise unconstrained "wiggles" */
#define SHAPE_WEIGHT 0.000001	/* Overal shaper weight contribution - err on side of smoothness */
#define SHAPE_BASE  0.001		/* 0 & 1 harmonic weight */
#define SHAPE_HBASE 0.01		/* 2nd and higher additional weight */

#define UDEV_WEIGHT 0.00005		/* Overal grid point uniformity weight spacing */

#define NTP 1					/* Number of test points mid grid */
#define MXSGORD 4				/* Sub-grid shaper orders to use */
#define MXLUORD 10				/* Input/output curve shaper harmonic orders to use */

/*
 * TTBD:
 *
 */

/* - - - - - - - - - - - - - - - - - */
/* Lookup a value though a sub-grid curve */
static double xfit3_subcurve(xfit3 *p, double in, int chan) {
	double rv;
	int i;

	/* Normalize */
	in = (in - p->in_min[chan])/(p->in_max[chan] - p->in_min[chan]);
	/* Quantize to grid */
	in *= (p->gres[chan]-1.0);
	i = (int)floor(in);
	if (i > (p->gres[chan]-2))
		i = (p->gres[chan]-2);
	in -= (double)i;				/* In value is now 0.0 - 1.0 between grid */

	rv = icxTransFunc(p->v + p->sub_offs[chan] + i * p->sluord[chan], p->sluord[chan], in); 

	/* Undo quantize to grid */
	rv += (double)i;
	rv /= (p->gres[chan]-1.0);
	rv = rv  * (p->in_max[chan] - p->in_min[chan]) + p->in_min[chan];

	return rv;
}

/* Inverse Lookup a value though a sub-grid curve */
static double xfit3_invsubcurve(xfit3 *p, double in, int chan) {
	double rv;
	int i;

	/* Normalize */
	in = (in - p->in_min[chan])/(p->in_max[chan] - p->in_min[chan]);
	/* Quantize to grid */
	in *= (p->gres[chan]-1.0);
	i = (int)floor(in);
	if (i > (p->gres[chan]-2))
		i = (p->gres[chan]-2);
	in -= (double)i;				/* In value is now 0.0 - 1.0 between grid */

	rv = icxInvTransFunc(p->v + p->sub_offs[chan] + i * p->sluord[chan], p->sluord[chan], in); 

	/* Undo quantize to grid */
	rv += (double)i;
	rv /= (p->gres[chan]-1.0);
	rv = rv  * (p->in_max[chan] - p->in_min[chan]) + p->in_min[chan];

	return rv;
}

/* - - - - - - - - - - - - - - - - - */
/* Overall transfer curve functions */

/* Lookup a value though the input and sub-grid curves */
static double xfit3_incurve(xfit3 *p, double in, int chan) {
	double rv;
	if (p->tcomb & oc3_i)
		rv = icxSTransFunc(p->v + p->in_offs[chan], p->iluord[chan], in,  
			                       p->in_min[chan], p->in_max[chan]);
	else
		rv = in;

	if (p->tcomb & oc3_s)
		rv = xfit3_subcurve(p, rv, chan);
	return rv;
}

/* Inverse Lookup a value though the input and sub-grid curve */
static double xfit3_invincurve(xfit3 *p, double in, int chan) {
	double rv;
	if (p->tcomb & oc3_s)
		rv = xfit3_invsubcurve(p, in, chan);
	else
		rv = in;
	if (p->tcomb & oc3_i)
		rv = icxInvSTransFunc(p->v + p->in_offs[chan], p->iluord[chan], rv,  
			                       p->in_min[chan], p->in_max[chan]);
	return rv;
}

/* Lookup a value though an output curve */
static double xfit3_outcurve(xfit3 *p, double in, int chan) {
	double rv;
	if (p->tcomb & oc3_o)
		rv = icxSTransFunc(p->v + p->out_offs[chan], p->oluord[chan], in,  
		                       p->out_min[chan], p->out_max[chan]);
	else
		rv = in;
	return rv;
}

/* Inverse Lookup a value though an output curve */
static double xfit3_invoutcurve(xfit3 *p, double in, int chan) {
	double rv;
	if (p->tcomb & oc3_o)
		rv = icxInvSTransFunc(p->v + p->out_offs[chan], p->oluord[chan], in,  
		                       p->out_min[chan], p->out_max[chan]);
	else
		rv = in;
	return rv;
}


/* Delta E squared replacement for device value output */
static double deverr(double *in1, double *in2, int di) {
	int i;
	double tt, rv = 0.0;
	for (i = 0; i < di; i++) {
		tt = 100.0 * (in1[i] - in2[i]);
		rv += tt * tt;
	}
	return rv;
}

/* - - - - - - - - - */
/* Convert all channels call */

/* callback to convert fdi values from out' to out values */ 
static void ocurv(void *cntx, double *out, double *in) {
	xfit3 *p = (xfit3 *)cntx;
	int f;

	for (f = 0; f < p->fdi; f++) {
		out[f] = icxSTransFunc(p->v + p->out_offs[f], p->oluord[f], in[f],  
			                                 p->out_min[f], p->out_max[f]);
	}
}

/* callback to convert fdi values from out to out' values */ 
static void iocurv(void *cntx, double *out, double *in) {
	xfit3 *p = (xfit3 *)cntx;
	int f;

	for (f = 0; f < p->fdi; f++) {
		out[f] = icxInvSTransFunc(p->v + p->out_offs[f], p->oluord[f], in[f],  
		                                        p->out_min[f], p->out_max[f]);
	}
}

/* - - - - - - - - - */

#define MAX_NTP 5

/* General optimisation assessment routine */
/* Given a grid cell span, and using the current input and */
/* output curves in xfit3, estimate the resulting linear */
/* interpolation error. */
static double span_err(
	xfit3 *p,
	double *pspande,		/* if not NULL return the span width delta E */
							/* as measured by delta E of the span */
	int ch,					/* Axis being evaluated */
	int ntp,				/* Number of sub-grid test points, min 1 */
	double x0, double x1	/* Span range in input space, x0 < x1 */
) {
	int di = p->di;
	int fdi = p->fdi;
	int i, j, ee, e, f, n;
	int gc[MXRI];
	co lup[MAX_NTP+2];		/* in to out sample points */
	co lupd[MAX_NTP+2];		/* in' to out' sample points */ 
	double S, Sx, Sy[MXDO], Syy[MXDO], Sxx, Sxy[MXDO];		/* Regression parameters */
	double a[MXDO], b[MXDO];	/* Regression line parameters */
	double rv = 0.0;		/* Sum of delta E's */
	int nde = 0;			/* Number of delta E's summed */
	double spande = 0.0;	/* Average span delta E */
	int nspans = 0;			/* Number of spans to average */

	if (ntp > MAX_NTP)
		error("xfit3 span_err got ntp > %d",MAX_NTP);

	ntp += 2;		/* Include end points */
	S = ntp;

	/* Go through each interpolation in the major axis direction */	
	/* for each rspl node in the orthogonal directions. */
	for (ee = 0; ee < di; ee++)
		gc[ee] = 0;	/* init coords */

	for (ee = 0; ee < di;) {		/* Until index completed */

		Sx = Sxx = 0.0;
		for (f = 0; f < p->fdi; f++)
			Sy[f] = Sxy[f] = 0.0;

		for (n = 0; n < ntp; n++) {	/* Lookup points along span */
		
			for (e = 0; e < di; e++)
				lup[n].p[e] = (gc[e]/(p->orthres-1.0) * (p->imax[e] - p->imin[e])) + p->imin[e]; 
			lup[n].p[ch] = (n/(ntp-1.0) * (x1 - x0)) + x0; 

			p->iaxs[ch]->interp(p->iaxs[ch], &lup[n]);

			/* Convert the sample point into an in' to out' mapping, */
			/* and keep linear regression sums */
			for (e = 0; e < p->di; e++) {
				lupd[n].p[e] = icxSTransFunc(p->v + p->in_offs[e], p->iluord[e], lup[n].p[e],  
			                                                      p->in_min[e], p->in_max[e]);
				if (p->tcomb & oc3_s)
					lupd[n].p[e] = xfit3_subcurve(p, lupd[n].p[e], e);
			}
			Sx += lupd[n].p[ch];
			Sxx += lupd[n].p[ch] * lupd[n].p[ch];
			for (f = 0; f < p->fdi; f++) {
				lupd[n].v[f] = icxInvSTransFunc(p->v + p->out_offs[f], p->oluord[f], lup[n].v[f],  
			                                                        p->out_min[f], p->out_max[f]);
				Sy[f] += lupd[n].v[f];
				Syy[f] += lupd[n].v[f] * lupd[n].v[f];
				Sxy[f] += lupd[n].p[ch] * lupd[n].v[f];
			}
		}

		/* compute span DE here by taking the delta E between the first and last values */
		if (pspande != NULL) {
			if (p->flags & XFIT_FM_INPUT) {	/* x as a function of y: x = a + by  */
				spande += p->to_de2(p->cntx2, lup[0].p, lup[ntp-1].p);
				nspans++;
			} else {
				spande += p->to_de2(p->cntx2, lup[0].v, lup[ntp-1].v);
				nspans++;
			}
		}

		/* Compute linear regression fit and delta E squared */
		if (p->flags & XFIT_FM_INPUT) {	/* x as a function of y: x = a + by  */
			double d, lp[MXDI];

			for (f = 0; f < p->fdi; f++) {
				d = S * Syy[f] - Sy[f] * Sy[f];
				a[f] = (Syy[f] * Sx - Sy[f] * Sxy[f])/d;
				b[f] = (S * Sxy[f] - Sy[f] * Sx)/d;
			}

			for (n = 0; n < ntp; n++) {

				for (e = 0; e < di; e++)
					lp[e] = lup[n].p[e];

				for (f = 0; f < p->fdi; f++) {
					lp[ch] = a[f] + b[f] * lupd[n].v[ch];	/* Interpolated value, convert to in */

					if (p->tcomb & oc3_s)
						lp[ch] = xfit3_invsubcurve(p, lp[ch], ch);
					lp[ch] = icxInvSTransFunc(p->v + p->in_offs[ch], p->oluord[ch], lp[ch],
			                                            p->in_min[ch], p->in_max[ch]);
					rv += p->to_de2(p->cntx2, lup[n].p, lp);
					/* Note we're accumulating fdi times as many delta E's, */
					/* but only in the ch direction. This should be roughly the same */
					/* as one delta E of fdi changed components ? */
				}
			}
		} else {		/* y as a function of x: y = a + bx */
			double d, lv[MXDO];

			d = S * Sxx - Sx * Sx;
			for (f = 0; f < p->fdi; f++) {
				a[f] = (Sxx * Sy[f] - Sx * Sxy[f])/d;
				b[f] = (S * Sxy[f] - Sx * Sy[f])/d;
			}

			for (n = 0; n < ntp; n++) {
				for (f = 0; f < p->fdi; f++) {
					lv[f] = a[f] + b[f] * lupd[n].p[ch];	/* Interpolated value, convert to out */
					lv[f] = icxSTransFunc(p->v + p->out_offs[f], p->oluord[f], lv[f],
			                                            p->out_min[f], p->out_max[f]);
				}
//printf("~1 n = %d, lv %f %f %f, tv %f %f %f\n",n,lv[0],lv[1],lv[2],lup[n].v[0],lup[n].v[1],lup[n].v[2]);
				/* We're accumulating one delta E per test point, */
				/* but the delta E includes delta's in all fdi directions at once */
				rv += p->to_de2(p->cntx2, lup[n].v, lv);
			}
		}
		nde += ntp;

		/* Increment index and start points orthogonal to ch */
		for (ee = 0; ee < di; ee++) {
			if (ee == ch) 		/* Don't go in direction xi */
				continue;
			gc[ee]++;
			if (gc[ee] < p->orthres)
				break;	/* No carry */
			gc[ee] -= p->orthres;			/* Reset coord */
		}
	}

	/* Normalize rv to be intergral over span range */
	rv *= (x1 - x0)/((p->in_max[ch] - p->in_min[ch]) * (double)nde);

	/* Return delta E of span width */
	if (pspande != NULL) {
		spande /= (double)nspans;
		*pspande = spande;
	}

	return rv;
}

/* Output/sub-grid curve cached test point creation routine. */
/* Given a grid cell span, and using the current input and */
/* output curves in xfit3, allocate and initialize all the test points */
/* for testing within that span interpolation error. */
/* Return nz on error */
static int span_setup(
	xfit3 *p,
	int ch,					/* Axis being evaluated */
	int ntp,				/* Number of sub-grid test points, min 1 */
	double x0, double x1,	/* Span range in input space, x0 < x1 */
	int sx					/* Span index, range 0..gres[ch]-2 */
) {
	int di = p->di;
	int fdi = p->fdi;
	int i, j, ee, e, f, n;
	int gc[MXRI];
	co *lup;				/* in to out sample points */
	co *lupd;				/* in' to out' sample points */ 
	int ns = 1;
	xfit3_ctp *ctp;

//printf("~1 span_setup ch %d, span index %d, %f-%f\n",ch,sx,x0,x1);
	/* Compute the number of sets */
	for (ee = 0; ee < di; ee++) {
		if (ee == ch)
			continue;
		ns *= p->orthres;
	}

	/* Comute number of points in each set including grid points */
	ntp += 2;

	/* Initialise arrays if they are new */
	if (p->ctp == NULL) {
		if ((p->ctp = (xfit3_ctp **)calloc(p->di, sizeof(xfit3_ctp *))) == NULL)
			return 2;
	}
	if (p->ctp[ch] == NULL) {
		if ((p->ctp[ch] = (xfit3_ctp *)calloc(p->gres[ch]-1, sizeof(xfit3_ctp))) == NULL)
			return 2;
	}
	ctp = &p->ctp[ch][sx];

	if (ctp->lup == NULL) {
		if ((ctp->lup = (co *)calloc(ns * ntp, sizeof(co))) == NULL)
			return 2;
		if ((ctp->lupd = (co *)calloc(ns * ntp, sizeof(co))) == NULL)
			return 2;

		ctp->ns = ns;
		ctp->ntp = ntp;
		ctp->x0 = x0;
		ctp->x1 = x1;

	} else {
		/* Check the size hasn't changed */
		if (ctp->ns != ns || ctp->ntp != ntp
		 || ctp->x0 != x0 || ctp->x1 != x1) {
			warning("xfit3, span_setup assert: ns or np has changed");
			return 1;
		}
	}

	lup = ctp->lup;
	lupd = ctp->lupd;

	/* Go through each interpolation in the major axis direction */	
	/* for each rspl node in the orthogonal directions. */
	for (ee = 0; ee < di; ee++)
		gc[ee] = 0;	/* init coords */

	for (ee = 0; ee < di;) {		/* Until index completed */
		for (n = 0; n < ntp; n++) {	/* Lookup points along span */
		
			for (e = 0; e < di; e++)
				lup[n].p[e] = (gc[e]/(p->orthres-1.0) * (p->imax[e] - p->imin[e])) + p->imin[e]; 
			lup[n].p[ch] = (n/(ntp-1.0) * (x1 - x0)) + x0; 

			p->iaxs[ch]->interp(p->iaxs[ch], &lup[n]);

			/* Convert the sample point into an in' to out' mapping using current curves */
			for (e = 0; e < p->di; e++) {
				lupd[n].p[e] = icxSTransFunc(p->v + p->in_offs[e], p->iluord[e], lup[n].p[e],  
			                                                      p->in_min[e], p->in_max[e]);
				if (p->tcomb & oc3_s)
					lupd[n].p[e] = xfit3_subcurve(p, lupd[n].p[e], e);
			}
			for (f = 0; f < p->fdi; f++)
				lupd[n].v[f] = icxInvSTransFunc(p->v + p->out_offs[f], p->oluord[f], lup[n].v[f],  
			                                                        p->out_min[f], p->out_max[f]);
//printf("~1 ch %d, span %f-%f, set %d, p = %f %f %f %f, v = %f %f %f\n",
//ch,ctp->x0,ctp->x1, (lup-ctp->lup)/ntp,
//lup[n].p[0], lup[n].p[1], lup[n].p[2], lup[n].p[3],
//lup[n].v[0], lup[n].v[1], lup[n].v[2]);
		}
		lup += ntp;
		lupd += ntp;
	
		/* Increment index and start points orthogonal to ch */
		for (ee = 0; ee < di; ee++) {
			if (ee == ch) 		/* Don't go in direction xi */
				continue;
			gc[ee]++;
			if (gc[ee] < p->orthres)
				break;	/* No carry */
			gc[ee] -= p->orthres;			/* Reset coord */
		}
	}

	return 0;
}

/* Output curve error evaluation routine using cached test points */
/* Compute on the expecations that the output curves may have changed */
static double out_span_err(
	xfit3 *p,
	int ch,					/* Axis being evaluated */
	xfit3_ctp *ctp			/* test points for this span on this axis */
) {
	int di = p->di;
	int fdi = p->fdi;
	int i, j, ee, e, f, n;
	int ntp;				/* Number of test points including grid points */
	co *lup;				/* in to out sample points */
	co *lupd;				/* in' to out' sample points */ 
	double S, Sx, Sy[MXDO], Syy[MXDO], Sxx, Sxy[MXDO];		/* Regression parameters */
	double a[MXDO], b[MXDO];	/* Regression line parameters */
	double rv = 0.0;		/* Sum of delta E's */
	int nde = 0;			/* Number of delta E's summed */

//printf("~1 out_span_err ch %d\n",ch);
	ntp = ctp->ntp;			/* Number of test points per set */
	S = ntp;

	/* Go through each interpolation in the major axis direction */	
	/* for each rspl node in the orthogonal directions. */

	lup = ctp->lup;
	lupd = ctp->lupd;

	for (i = 0; i < ctp->ns; i++) {		/* Until done all sets */

		Sx = Sxx = 0.0;
		for (f = 0; f < p->fdi; f++)
			Sy[f] = Sxy[f] = 0.0;

		/* Compute linear regression sums */
		for (n = 0; n < ntp; n++) {	
			Sx += lupd[n].p[ch];
			Sxx += lupd[n].p[ch] * lupd[n].p[ch];
			for (f = 0; f < p->fdi; f++) {
				if (f == p->och)	/* Update out mapping for the channel being optimized */
					lupd[n].v[f] = icxInvSTransFunc(p->v + p->out_offs[f], p->oluord[f],
					                                 lup[n].v[f], p->out_min[f], p->out_max[f]);
				Sy[f] += lupd[n].v[f];
				Syy[f] += lupd[n].v[f] * lupd[n].v[f];
				Sxy[f] += lupd[n].p[ch] * lupd[n].v[f];
			}
//printf("~1 ch %d, span %f-%f, set %d, p = %f %f %f %f, v = %f %f %f\n",
//ch,ctp->x0,ctp->x1,i,
//lup[n].p[0], lup[n].p[1], lup[n].p[2], lup[n].p[3],
//lup[n].v[0], lup[n].v[1], lup[n].v[2]);
		}

		/* Compute linear regression fit and delta E squared */
		if (p->flags & XFIT_FM_INPUT) {	/* x as a function of y: x = a + by  */
			double d, lp[MXDI];

			for (f = 0; f < p->fdi; f++) {
				d = S * Syy[f] - Sy[f] * Sy[f];
				a[f] = (Syy[f] * Sx - Sy[f] * Sxy[f])/d;
				b[f] = (S * Sxy[f] - Sy[f] * Sx)/d;
			}

			for (n = 0; n < ntp; n++) {

				for (e = 0; e < di; e++)
					lp[e] = lup[n].p[e];

				for (f = 0; f < p->fdi; f++) {
					lp[ch] = a[f] + b[f] * lupd[n].v[ch];	/* Interpolated value, convert to in */
					lp[ch] = icxInvSTransFunc(p->v + p->in_offs[ch], p->oluord[ch], lp[ch],
			                                            p->in_min[ch], p->in_max[ch]);
					rv += p->to_de2(p->cntx2, lup[n].p, lp);
					/* Note we're accumulating fdi times as many delta E's, */
					/* but only in the ch direction. This should be roughly the same */
					/* as one delta E of fdi changed components ? */
				}
			}
		} else {		/* y as a function of x: y = a + bx */
			double d, lv[MXDO];

			d = S * Sxx - Sx * Sx;
			for (f = 0; f < p->fdi; f++) {
				a[f] = (Sxx * Sy[f] - Sx * Sxy[f])/d;
				b[f] = (S * Sxy[f] - Sx * Sy[f])/d;
			}

			for (n = 0; n < ntp; n++) {
				for (f = 0; f < p->fdi; f++) {
					lv[f] = a[f] + b[f] * lupd[n].p[ch];	/* Interpolated value, convert to out */
					lv[f] = icxSTransFunc(p->v + p->out_offs[f], p->oluord[f], lv[f],
			                                            p->out_min[f], p->out_max[f]);
				}
//printf("~1 n = %d, lv %f %f %f, tv %f %f %f\n",n,lv[0],lv[1],lv[2],lup[n].v[0],lup[n].v[1],lup[n].v[2]);
				/* We're accumulating one delta E per test point, */
				/* but the delta E includes delta's in all fdi directions at once */
				rv += p->to_de2(p->cntx2, lup[n].v, lv);
			}
		}
		nde += ntp;
		lup += ntp;
		lupd += ntp;
	}

	/* Normalize rv to be intergral over span range */
	rv *= (ctp->x1 - ctp->x0)/((p->in_max[ch] - p->in_min[ch]) * (double)nde);

	return rv;
}


/* Output curve error evaluation routine using cached test points. */
/* Compute on the expecations that the input sub-curves may have changed */
static double out_span_err2(
	xfit3 *p,
	int ch,					/* Axis being evaluated */
	xfit3_ctp *ctp			/* test points for this span on this axis */
) {
	int di = p->di;
	int fdi = p->fdi;
	int i, j, ee, e, f, n;
	int ntp;				/* Number of test points including grid points */
	co *lup;				/* in to out sample points */
	co *lupd;				/* in' to out' sample points */ 
	double S, Sx, Sy[MXDO], Syy[MXDO], Sxx, Sxy[MXDO];		/* Regression parameters */
	double a[MXDO], b[MXDO];	/* Regression line parameters */
	double rv = 0.0;		/* Sum of delta E's */
	int nde = 0;			/* Number of delta E's summed */

//printf("~1 out_span_err ch %d\n",ch);
	ntp = ctp->ntp;			/* Number of test points per set */
	S = ntp;

	/* Go through each interpolation in the major axis direction */	
	/* for each rspl node in the orthogonal directions. */
	lup = ctp->lup;
	lupd = ctp->lupd;

	for (i = 0; i < ctp->ns; i++) {		/* Until done all sets */

		Sx = Sxx = 0.0;
		for (f = 0; f < p->fdi; f++)
			Sy[f] = Sxy[f] = 0.0;

		/* Compute linear regression sums */

		for (n = 0; n < ntp; n++) {	/* Lookup points along span */
		
			/* Recompute in' sample point values, */
			/* and keep linear regression sums */
			for (e = 0; e < p->di; e++) {
				lupd[n].p[e] = icxSTransFunc(p->v + p->in_offs[e], p->iluord[e], lup[n].p[e],  
			                                                      p->in_min[e], p->in_max[e]);
				if (p->tcomb & oc3_s)
					lupd[n].p[e] = xfit3_subcurve(p, lupd[n].p[e], e);
			}
			Sx += lupd[n].p[ch];
			Sxx += lupd[n].p[ch] * lupd[n].p[ch];
			for (f = 0; f < p->fdi; f++) {
				Sy[f] += lupd[n].v[f];
				Syy[f] += lupd[n].v[f] * lupd[n].v[f];
				Sxy[f] += lupd[n].p[ch] * lupd[n].v[f];
			}
//printf("~1 ch %d, span %f-%f, set %d, p = %f %f %f %f, v = %f %f %f\n",
//ch,ctp->x0,ctp->x1,i,
//lup[n].p[0], lup[n].p[1], lup[n].p[2], lup[n].p[3],
//lup[n].v[0], lup[n].v[1], lup[n].v[2]);
		}

		/* Compute linear regression fit and delta E squared */
		if (p->flags & XFIT_FM_INPUT) {	/* x as a function of y: x = a + by  */
			double d, lp[MXDI];

			for (f = 0; f < p->fdi; f++) {
				d = S * Syy[f] - Sy[f] * Sy[f];
				a[f] = (Syy[f] * Sx - Sy[f] * Sxy[f])/d;
				b[f] = (S * Sxy[f] - Sy[f] * Sx)/d;
			}

			for (n = 0; n < ntp; n++) {

				for (e = 0; e < di; e++)
					lp[e] = lup[n].p[e];

				for (f = 0; f < p->fdi; f++) {
					lp[ch] = a[f] + b[f] * lupd[n].v[ch];	/* Interpolated value, convert to in */
					lp[ch] = icxInvSTransFunc(p->v + p->in_offs[ch], p->oluord[ch], lp[ch],
			                                            p->in_min[ch], p->in_max[ch]);
					rv += p->to_de2(p->cntx2, lup[n].p, lp);
					/* Note we're accumulating fdi times as many delta E's, */
					/* but only in the ch direction. This should be roughly the same */
					/* as one delta E of fdi changed components ? */
				}
			}
		} else {		/* y as a function of x: y = a + bx */
			double d, lv[MXDO];

			d = S * Sxx - Sx * Sx;
			for (f = 0; f < p->fdi; f++) {
				a[f] = (Sxx * Sy[f] - Sx * Sxy[f])/d;
				b[f] = (S * Sxy[f] - Sx * Sy[f])/d;
			}

			for (n = 0; n < ntp; n++) {
				for (f = 0; f < p->fdi; f++) {
					lv[f] = a[f] + b[f] * lupd[n].p[ch];	/* Interpolated value, convert to out */
					lv[f] = icxSTransFunc(p->v + p->out_offs[f], p->oluord[f], lv[f],
			                                            p->out_min[f], p->out_max[f]);
				}
//printf("~1 n = %d, lv %f %f %f, tv %f %f %f\n",n,lv[0],lv[1],lv[2],lup[n].v[0],lup[n].v[1],lup[n].v[2]);
				/* We're accumulating one delta E per test point, */
				/* but the delta E includes delta's in all fdi directions at once */
				rv += p->to_de2(p->cntx2, lup[n].v, lv);
			}
		}
		nde += ntp;
		lup += ntp;
		lupd += ntp;
	}

	/* Normalize rv to be intergral over span range */
	rv *= (ctp->x1 - ctp->x0)/((p->in_max[ch] - p->in_min[ch]) * (double)nde);

	return rv;
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Input shaper optimisation function handed to powell() */
static double xfit3ifunc(void *edata, double *v) {
	xfit3 *p = (xfit3 *)edata;
	int i;
	double *pp, de = 0.0, rv;
	double pgx, gx;	
	double tt, w, avuerr = 0.0, avudev, smv = 0.0;

//printf("~1 chan %d\n",p->och);

	/* Copy the values into place */
	pp = &p->v[p->opt_off];
	for (i = 0; i < p->opt_cnt; i++)
		*pp++ = v[i];

	/* Compute the locations of the grid points, and */
	/* sum the expected errors from those locations */
	pgx = -1.0;
	for (i = 0; i < p->gres[p->och]; i++) { 

		gx = i/(p->gres[p->och]-1.0);				/* Raw grid node location on clut */
//printf("~1 grid point %d, raw gx = %f\n",i,gx);
		gx = icxInvTransFunc(p->v + p->in_offs[p->och], p->oluord[p->och], gx);
		gx = gx * (p->in_max[p->och] - p->in_min[p->och]) + p->in_min[p->och];
//printf("~1 gx after input curve = %f\n",gx);
		
		if (pgx >= 0.0) {
			double ser;
			/* Estimated error squared due to this span */ 
//			Discrete version with cache
//			ser = p->iaxs[p->och]->gcso_fwd_err(p->iaxs[p->och], &p->uerrv[i], pgx, gx);

			/* General version */
			ser = span_err(p, &p->uerrv[i], p->och, NTP, pgx, gx);

//printf("~1 span %f -> %f : uer = %f, rv = %f\n",pgx,gx,p->uerrv[i],ser);
			avuerr += p->uerrv[i];
			de += ser;
		}
		pgx = gx;
	}

	/* Compute the average deviation squared from the average span delta E */
	/* This encourages even grid spacing in perceptual space. */
	avuerr /= (p->gres[p->och]-1.0);
	for (avudev = 0.0, i = 1; i < p->gres[p->och]; i++) { 
		double tt;
		tt = p->uerrv[i] - avuerr;
		avudev += tt * tt;
	}
	avudev /= (p->gres[p->och]-1.0);
	avudev = sqrt(avudev) * UDEV_WEIGHT;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	for (i = 0; i < p->opt_cnt; i++) {
		tt = v[i];
		tt *= tt;
		if (i <= 1)	/* First or second curves */
			w = SHAPE_BASE;
		else
			w = SHAPE_BASE + (i-1) * SHAPE_HBASE;
		smv += w * tt;
	}
	smv *= SHAPE_WEIGHT;
	rv = de + avudev + smv;
#ifdef DEBUG
printf("~1 xfit3ifunc returning (sm %f, ue %f, de %f) %f\n",smv, avudev, de, rv);
#endif

	if (p->verb)
		printf("."), fflush(stdout);
	return rv;
}

/* Fake version of above that simulate partial de's for conjgrad */
static double dxfit3ifunc(void *edata, double *dv, double *v) {
	xfit3 *p = (xfit3 *)edata;
	int i;
	double rv;

	rv = xfit3ifunc(edata, v);

	/* Compute each parameter delta */
	for (i = 0; i < p->opt_cnt; i++) {
		double del;

		v[i] += 1e-7;
		dv[i] = (xfit3ifunc(edata, v) - rv)/1e-7;
		v[i] -= 1e-7;
	}
	return rv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Output shaper optimisation function handed to powell() */
/* Slow general one that uses span_err() */
static double xfit3ofunc(void *edata, double *v) {
	xfit3 *p = (xfit3 *)edata;
	int e, i;
	double *pp, de = 0.0, rv;
	double pgx, gx;	
	double tt, w, smv = 0.0;

//printf("~1 chan %d\n",p->och);

	/* Copy the values into place */
	pp = &p->v[p->opt_off];
	for (i = 0; i < p->opt_cnt; i++)
		*pp++ = v[i];

	/* For each input axis */
	for (e = 0; e < p->di; e++) {

		/* Compute the locations of the grid points, and */
		/* sum the expected errors from those locations */
		pgx = -1.0;
		for (i = 0; i < p->gres[e]; i++) { 
	
			gx = i/(p->gres[e]-1.0);				/* Raw grid node location on clut */
//printf("~1 grid point %d, raw gx = %f\n",i,gx);
			gx = icxInvTransFunc(p->v + p->in_offs[e], p->oluord[e], gx);
			gx = gx * (p->in_max[e] - p->in_min[e]) + p->in_min[e];
//printf("~1 gx after input curve = %f\n",gx);
		
			if (pgx >= 0.0) {
				double ser;
				ser = span_err(p, NULL, e, NTP, pgx, gx);
//printf("~1 span %f -> %f : rv = %f\n",pgx,gx,ser);
				de += ser;
			}
			pgx = gx;
		}
	}
	de /= (double)p->di;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	for (i = 0; i < p->opt_cnt; i++) {
		tt = v[i];
		tt *= tt;
		if (i <= 1)	/* First or second curves */
			w = SHAPE_BASE;
		else
			w = SHAPE_BASE + (i-1) * SHAPE_HBASE;
		smv += w * tt;
	}
	smv *= SHAPE_WEIGHT;
	rv = de + smv;
#ifdef DEBUG
printf("~1 xfit3ofunc returning (sm %f, de %f) %f\n",smv, de, rv);
#endif

	if (p->verb)
		printf("."), fflush(stdout);
	return rv;
}

/* Fake version of above that simulate partial de's for conjgrad */
static double dxfit3ofunc(void *edata, double *dv, double *v) {
	xfit3 *p = (xfit3 *)edata;
	int i;
	double rv;

	rv = xfit3ofunc(edata, v);

	/* Compute each parameter delta */
	for (i = 0; i < p->opt_cnt; i++) {
		double del;

		v[i] += 1e-7;
		dv[i] = (xfit3ofunc(edata, v) - rv)/1e-7;
		v[i] -= 1e-7;
	}
	return rv;
}


/* - - - - - - - - - - - - - - - - - */

/* Output shaper optimisation function handed to powell() */
/* Fast one that uses cache of test points p->ctp */
static double xfit3ofunc2(void *edata, double *v) {
	xfit3 *p = (xfit3 *)edata;
	int e, i, j;
	double *pp, de = 0.0, rv;
	double pgx, gx;	
	double tt, w, smv = 0.0;

//printf("~1 chan %d\n",p->och);

	/* Copy the values into place */
	pp = &p->v[p->opt_off];
	for (i = 0; i < p->opt_cnt; i++)
		*pp++ = v[i];

	/* For each input axis */
	for (e = 0; e < p->di; e++) {

		/* For each grid span */
		for (i = 0; i < (p->gres[e]-1); i++)  {
			de += out_span_err(p, e, &p->ctp[e][i]);
		}
	}
	de /= (double)p->di;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	for (i = 0; i < p->opt_cnt; i++) {
		tt = v[i];
		tt *= tt;
		if (i <= 1)	/* First or second curves */
			w = SHAPE_BASE;
		else
			w = SHAPE_BASE + (i-1) * SHAPE_HBASE;
		smv += w * tt;
	}
	smv *= SHAPE_WEIGHT;
	rv = de + smv;
#ifdef DEBUG
printf("~1 xfit3ofunc returning (sm %f, de %f) %f\n",smv, de, rv);
#endif

	if (p->verb)
		printf("."), fflush(stdout);
	return rv;
}

/* Fake version of above that simulate partial de's for conjgrad */
static double dxfit3ofunc2(void *edata, double *dv, double *v) {
	xfit3 *p = (xfit3 *)edata;
	int i;
	double rv;

	rv = xfit3ofunc2(edata, v);

	/* Compute each parameter delta */
	for (i = 0; i < p->opt_cnt; i++) {
		double del;

		v[i] += 1e-7;
		dv[i] = (xfit3ofunc2(edata, v) - rv)/1e-7;
		v[i] -= 1e-7;
	}
	return rv;
}


/* - - - - - - - - - - - - - - - - - */

// ~~999
/* Sub-grid shaper optimisation function handed to powell() */
/* Uses cache of test points p->ctp */
static double xfit3sfunc2(void *edata, double *v) {
	xfit3 *p = (xfit3 *)edata;
	int e, i, j;
	double *pp, de = 0.0, rv;
	double pgx, gx;	
	double tt, w, smv = 0.0;

//printf("~1 chan %d\n",p->och);

	/* Copy the values into place */
	pp = &p->v[p->opt_off];
	for (i = 0; i < p->opt_cnt; i++)
		*pp++ = v[i];

	de = out_span_err2(p, p->och, &p->ctp[p->och][p->osp]);
	de *= (p->gres[p->och]-1.0);		/* Scale to whole axis */

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	for (i = 0; i < p->opt_cnt; i++) {
		tt = v[i];
		tt *= tt;
		if (i <= 1)	/* First or second curves */
			w = SHAPE_BASE;
		else
			w = SHAPE_BASE + (i-1) * SHAPE_HBASE;
		smv += w * tt;
	}
	smv *= SHAPE_WEIGHT;
// ~~~999
//	rv = de + smv;
	rv = de;
#ifdef DEBUG
printf("~1 xfit3sfunc2 returning (sm %f, de %f) %e\n",smv, de, rv);
#endif

	if (p->verb)
		printf("."), fflush(stdout);
	return rv;
}

/* Fake version of above that simulate partial de's for conjgrad */
static double dxfit3sfunc2(void *edata, double *dv, double *v) {
	xfit3 *p = (xfit3 *)edata;
	int i;
	double rv;

	rv = xfit3sfunc2(edata, v);

	/* Compute each parameter delta */
	for (i = 0; i < p->opt_cnt; i++) {
		double del;

		v[i] += 1e-7;
		dv[i] = (xfit3sfunc2(edata, v) - rv)/1e-7;
		v[i] -= 1e-7;
	}
	return rv;
}


/* - - - - - - - - - - - - - - - - - */

/* Output curve symetry optimisation function handed to powell() */
/* Just the order 0 value will be adjusted */
static double symoptfunc(void *edata, double *v) {
	xfit3 *p = (xfit3 *)edata;
	double out[1], in[1] = { 0.0 };
	int ch = p->och;		/* Output channel being adjusted for symetry */
	double rv;

	/* Copy the parameter being tested back into xfit3 */
	p->v[p->out_offs[ch]] = v[0];
	*out = icxSTransFunc(p->v + p->out_offs[ch], p->oluord[ch], *in,
							   p->out_min[ch], p->out_max[ch]);

	rv = out[0] * out[0];

#ifdef DEBUG
printf("~1 symoptfunc returning %f\n",rv);
#endif
	return rv;
}

/* - - - - - - - - - */

/* Set up for an optimisation run: */
/* Figure out parameters being optimised, */
/* copy them to start values, */
/* init and scale the search radius */
static void setup_xfit3(
xfit3 *p,
double *wv,		/* Return parameters to hand to optimiser */
double *sa,		/* Return search radius to hand to optimiser */
double transrad /* Nominal transfer curve search radius, 0.0 - 3.0 */
) {
	int i;
	
	/* Input curve */
	if (p->optt == 0) {
		p->opt_off = p->in_offs[p->och],
		p->opt_cnt = p->iluord[p->och];

		for (i = 0; i < p->opt_cnt; i++) {
			*wv++ = p->v[p->opt_off + i];
			*sa++ = transrad;
		}

	/* Sub-grid input curves */
	} else if (p->optt == 1) {
		p->opt_off = p->sub_offs[p->och] + p->osp * p->sluord[p->och],
		p->opt_cnt = p->sluord[p->och];

		for (i = 0; i < p->opt_cnt; i++) {
			*wv++ = p->v[p->opt_off + i];
			*sa++ = transrad;
		}

	/* Output curve */
	} else {
		p->opt_off = p->out_offs[p->och];
		p->opt_cnt = p->oluord[p->och];

		for (i = 0; i < p->opt_cnt; i++) {
			*wv++ = p->v[p->opt_off + i];
			*sa++ = transrad;
		}
	}
#ifdef NEVER
printf("~1 optt    = %d\n",p->optt);
printf("~1 och     = %d\n",p->och);
printf("~1 osp     = %d\n",p->osp);
printf("~1 opt_off = %d\n",p->opt_off);
printf("~1 opt_cnt = %d\n\n",p->opt_cnt);
#endif /* NEVER */
}

#ifdef DEBUG
/* Diagnostic */
static void dump_xfit3(
xfit3 *p
) {
	int i, e, f;
	double *b;			/* Base of parameters for this section */
	int di, fdi;
	di   = p->di;
	fdi  = p->fdi;

	/* Input curve */
	b = p->v + p->in_off;
	for (e = 0; e < di; b += p->iluord[e], e++) {
		printf("in %d = ",e);
		for (i = 0; i < p->iluord[e]; i++)
			printf("%f ",b[i]);
		printf("\n");
	}

	/* Sub-grid input curve */
	b = p->v + p->sub_off;
	for (e = 0; e < di; b += p->sluord[e], e++) {
		printf("in %d = ",e);
		for (i = 0; i < p->sluord[e] * (p->gres[e]-1); i++)
			printf("%f ",b[i]);
		printf("\n");
	}

	/* Output curve */
	b = p->v + p->out_off;
	for (f = 0; f < fdi; b += p->oluord[f], f++) {
		printf("out %d = ",f);
		for (i = 0; i < p->oluord[f]; i++)
			printf("%f ",b[i]);
		printf("\n");
	}
}
#endif /* DEBUG */

/* Do the fitting. */
/* return nz on error */
/* 1 = malloc error */
int xfit3_fit(
	xfit3 *p,
	int flags,				/* Flag values */
	int di,					/* Input dimensions */
	int fdi,				/* Output dimensions */
	int nodp,				/* Number of data points */
	cow *ipoints,			/* Array of data points to fit */
	int gres[MXDI],			/* clut resolutions being optimised for */
	double in_min[MXDI],	/* Input value scaling minimum */
	double in_max[MXDI],	/* Input value scaling maximum */
	double out_min[MXDO],	/* Output value scaling minimum */
	double out_max[MXDO],	/* Output value scaling maximum */
	int iord[],				/* Order of input shaper curve for each dimension */
	int sord[],				/* Order of input sub-grid shaper curve for each dimension */
	int oord[],				/* Order of output shaper curve for each dimension */
	optcomb3 tcomb,			/* Target elements to fit. */
	void *cntx2,			/* Context of to_de2 callback */
	double (*to_de2)(void *cntxw, double *in1, double *in2)	
							/* callback to convert in or out values to fit metric squared */
) {
	int i, e, f;
	int ee;
	double *b;			/* Base of parameters for this section */
	int poff;
	int maxxres = 0;	/* Maximum xres used */

	p->flags  = flags;
	if (flags & XFIT_VERB)
		p->verb = 1;
	else
		p->verb = 0;
	p->di      = di;
	p->fdi     = fdi;
	p->nodp    = nodp;
	p->ipoints = ipoints;
	p->tcomb   = tcomb;
	for (e = 0; e < di; e++)
		p->gres[e] = gres[e];
	p->cntx2   = cntx2;		/* to_de2 context */
	p->to_de2  = to_de2;

	/* Sanity protect shaper orders and save scaling factors. */
	for (e = 0; e < di; e++) {
		if (iord[e] > MXLUORD)
			p->iluord[e] = MXLUORD;
		else
			p->iluord[e] = iord[e];
		if (sord[e] > MXSGORD)
			p->sluord[e] = MXSGORD;
		else
			p->sluord[e] = sord[e];
		p->in_min[e] = in_min[e];
		p->in_max[e] = in_max[e];
	}
	for (f = 0; f < fdi; f++) {
		if (oord[f] > MXLUORD)
			p->oluord[f] = MXLUORD;
		else
			p->oluord[f] = oord[f];
		p->out_min[f] = out_min[f];
		p->out_max[f] = out_max[f];
	}

	/* Compute parameter offset and count information */
	p->in_off = 0;
	for (poff = p->in_off, p->in_cnt = 0, e = 0; e < di; e++) {
		p->in_offs[e] = poff;
		p->in_cnt += p->iluord[e];
		poff += p->iluord[e];
	}

	p->sub_off = p->in_off + p->in_cnt;
	for (poff = p->sub_off, p->sub_cnt = 0, e = 0; e < di; e++) {
		p->sub_offs[e] = poff;
		p->sub_cnt += p->sluord[e] * (p->gres[e]-1);
		poff += p->sluord[e] * (p->gres[e]-1);
	}

	p->out_off = p->sub_off + p->sub_cnt;
	for (poff = p->out_off, p->out_cnt = 0, f = 0; f < fdi; f++) {
		p->out_offs[f] = poff;
		p->out_cnt += p->oluord[f];
		poff += p->oluord[f];
	}

	p->tot_cnt = p->in_cnt + p->sub_cnt + p->out_cnt;

	/* Allocate space for parameter values */
	if (p->v != NULL) {
		free(p->v);
		p->v = NULL;
	}
	if (p->wv != NULL) {
		free(p->wv);
		p->wv = NULL;
	}
	if (p->sa != NULL) {
		free(p->sa);
		p->sa = NULL;
	}
	if ((p->v = (double *)calloc(p->tot_cnt, sizeof(double))) == NULL)
		return 1;
	if ((p->wv = (double *)calloc(p->tot_cnt, sizeof(double))) == NULL)
		return 1;
	if ((p->sa = (double *)calloc(p->tot_cnt, sizeof(double))) == NULL)
		return 1;

	/* Setup input curves to be linear */
	b = p->v + p->in_off;
	for (e = 0; e < di; b += p->iluord[e], e++) {
		for (i = 0; i < p->iluord[e]; i++) {
			b[i] = 0.0;
		}
	}

	/* Setup input sub-grid curves to be linear */
	b = p->v + p->sub_off;
	for (e = 0; e < di; b += (p->sluord[e] * (p->gres[e]-1)), e++) {
		for (i = 0; i < (p->sluord[e] * (p->gres[e]-1)); i++)
			b[i] = 0.0;
	}

	/* Setup output curves to be linear */
	b = p->v + p->out_off;
		
	for (f = 0; f < fdi; b += p->oluord[f], f++) {
		for (i = 0; i < p->oluord[f]; i++) {
			b[i] = 0.0;
		}
	}

printf("~1 tcomb7 = 0x%x\n",p->tcomb);
// ~~99
//	if ((p->tcomb & oc3_iso) == 0 )
//		return 0; 		/* Nothing to do */

	/* Create a rspl that is oriented for high resolution */
	/* in each input dimension to act as detail model for */
	/* optimization of the clut grid point locations */
#ifdef NEVER
	if (di == 1)
		p->orthres = 10;
	else if (di == 2)
		p->orthres = 10;
	else if (di == 3)
		p->orthres = 8;
	else if (di == 4)
		p->orthres = 5;
	else
		p->orthres = 4;

#else		/* Faster */
	if (di == 1)
		p->orthres = 10;
	else if (di == 2)
		p->orthres = 8;
	else if (di == 3)
		p->orthres = 5;
	else if (di == 4)
		p->orthres = 4;
	else
		p->orthres = 3;
#endif

	for (ee = 0; ee < di; ee++) {
		int xres[MXDI];
		double avgdev[MXDI];

		if ((p->iaxs[ee] = new_rspl(di, fdi)) == NULL)
			return 1;

		for  (e = 0; e < di; e++) {
			xres[e] = p->orthres;
			avgdev[e] = 0.005;
		}
		/* Make axis oriented resolution higher than intended */
		/* final resolution so that the interpolation errors can */
		/* be estimated, but not so high as to be too slow to compute. */
		xres[ee] = gres[ee] * 6;
		if (xres[ee] > 128) {
			xres[ee] = 128;
			if (128/gres[ee] < 4) {
				xres[ee] = gres[ee] * 4;
				if (xres[ee] > 256)
					xres[ee] = 256;
			}
		}
		if (xres[ee] > maxxres)
			maxxres = xres[ee];

printf("~1 creating rspl for axis %d, res %d\n",ee,xres[ee]);
		p->iaxs[ee]->fit_rspl_w(p->iaxs[ee], RSPL_SYMDOMAIN, ipoints, nodp, NULL, NULL, xres,
			NULL, NULL, 0.2, avgdev);

printf("~1 done\n");
	}

	/* Make sure we know exactly the input range of the rspls */
	p->iaxs[0]->get_in_range(p->iaxs[0], p->imin, p->imax);

	/* Allocate optimization scatch space */
	if (p->uerrv != NULL) {
		free(p->uerrv);
		p->uerrv = NULL;
	}
	if ((p->uerrv = (double *)calloc(maxxres, sizeof(double))) == NULL)
		return 1;

	/* - - - - - - - - - - - */
	for (ee = 0; ee < di; ee++) {
		/* Setup for (fast, cached rspl implemented) input curve optimization */

printf("~1 setting up rspl gcso for axis %d\n",ee);
		if (p->iaxs[ee]->gcso_setup(
			p->iaxs[ee],	/* rspl */
			ee,				/* index of input axis of interest */
			0,				/* Compute error in input space rather than output space */
			(void *)p,		/* Context of callbacks */
			ocurv,			/* callback to convert fdi values from out' to out values */
			iocurv,			/* callback to convert fdi values from out to out' values */
			p->cntx2,		/* Context to to_de2 callback */
			p->to_de2 		/* callback to convert error values to delta E */
		) != 0)
			return 1;
printf("~1 done\n");
	}

	/* Optimise the input curves */
	if ((p->tcomb & oc3_i) != 0 ) {

		for (e = 0; e < di; e++) {
		
			if (p->verb)
				printf("\nAbout to optimise input curve %d\n",e);

			p->optt = 0;
			p->och = e;
			setup_xfit3(p, p->wv, p->sa, 2.0); 
#ifdef NEVER
			if (powell(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit3ifunc, (void *)p) != 0)
				warning("set_icxLuLut: Powell failed to converge");
#else
			/* fake conjgrad to see how speed is */
			if (conjgrad(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit3ifunc, dxfit3ifunc, (void *)p) != 0)
				warning("set_icxLuLut: Powell failed to converge");
#endif
			for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
				p->v[p->opt_off + i] = p->wv[i];
		}
#ifdef DEBUG
printf("\nAfter input opt:\n");
dump_xfit3(p);
#endif
	}

	/* - - - - - - - - - - - */
	/* Optimise the input sub-grid curves */
// ~~999
	if ((p->tcomb & oc3_s) != 0 ) {

printf("~1 about to setup output test points\n");
		/* Setup the test points for all the input dimensions and grid spans */
		for (e = 0; e < di; e++) {
			double pgx, gx;	

			/* Compute the locations of the grid points */
			pgx = -1.0;
			for (i = 0; i < p->gres[p->och]; i++) { 

				gx = i/(p->gres[e]-1.0);				/* Raw grid node location on clut */
				gx = icxInvTransFunc(p->v + p->in_offs[e], p->oluord[e], gx);
				gx = gx * (p->in_max[e] - p->in_min[p->och]) + p->in_min[e];
	
				if (pgx >= 0.0) {
					if (span_setup(p, e, NTP, pgx, gx, i-1) != 0) 
						return 1;
				}
				pgx = gx;
			}
		}
printf("~1 output test points have been set\n");

		/* Do sub-grid optimisation */

		/* For each input axis */
		for (e = 0; e < p->di; e++) {
			int n;

			/* For each grid span */
			for (n = 0; n < (p->gres[e]-1); n++)  {
				p->optt = 1;
				p->och = e;
				p->osp = n;

printf("~1 optimising sub grid axis %d, grid %d\n",e,n);
				setup_xfit3(p, p->wv, p->sa, 2.0); 
#ifdef NEVER
				if (powell(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit3sfunc2, (void *)p) != 0)
					warning("set_icxLuLut: Powell failed to converge");
#else
				if (conjgrad(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit3sfunc2, dxfit3sfunc2, (void *)p) != 0)
					warning("set_icxLuLut: Powell failed to converge");
#endif
				for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
					p->v[p->opt_off + i] = p->wv[i];
			}
		}
	}

	/* - - - - - - - - - - - */
	/* Optimise the output curves */
	if ((p->tcomb & oc3_o) != 0 ) {

		for (f = 0; f < fdi; f++) {

printf("~1 about to setup output test points\n");
			/* Setup the test points for all the input dimensions and grid spans */
			for (e = 0; e < di; e++) {
				double pgx, gx;	

				/* Compute the locations of the grid points */
				pgx = -1.0;
				for (i = 0; i < p->gres[p->och]; i++) { 
	
					gx = i/(p->gres[e]-1.0);				/* Raw grid node location on clut */
					gx = icxInvTransFunc(p->v + p->in_offs[e], p->oluord[e], gx);
					gx = gx * (p->in_max[e] - p->in_min[p->och]) + p->in_min[e];
		
					if (pgx >= 0.0) {
						if (span_setup(p, e, NTP, pgx, gx, i-1) != 0) 
							return 1;
					}
					pgx = gx;
				}
			}
printf("~1 output test points have been set\n");

			if (p->verb)
				printf("\nAbout to optimise output curve %d\n", f);

			p->optt = 2;
			p->och = f;
			setup_xfit3(p, p->wv, p->sa, 2.0); 

#ifdef NEVER
			if (powell(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit3ofunc2, (void *)p) != 0)
				warning("set_icxLuLut: Powell failed to converge");
#else
			if (conjgrad(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit3ofunc2, dxfit3ofunc2, (void *)p) != 0)
				warning("set_icxLuLut: Powell failed to converge");
#endif
			for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
				p->v[p->opt_off + i] = p->wv[i];
		}
#ifdef DEBUG
printf("\nAfter output opt:\n");
dump_xfit3(p);
#endif

	}
	/* - - - - - - - - - - - */

#ifdef NEVER
	/* Adjust output curve white point */
	if (p->flags & XFIT_OUT_ZERO) {

		if (p->verb)
			printf("\nAbout to adjust a and b output curves for white point\n");

		for (f = 1; f < 3; f++) {
			p->och = f;
			p->wv[0] = p->v[p->out_offs[f]];	/* Current parameter value */
			p->sa[0] = 0.1;					/* Search radius */
			if (powell(NULL, 1, p->wv, p->sa, 0.0000001, 1000, symoptfunc, (void *)p) != 0)
				error("set_icxLuLut: Powell failed to converge");
			p->v[p->out_offs[f]] = p->wv[0];	/* Copy results back */
		}
	}
#endif /* NEVER */

	if (p->verb)
		printf("\n");

#ifdef DEBUG
printf("\nFinal parameters:\n");
dump_xfit3(p);
#endif

#ifdef DEBUG_PLOT
	{
#define	XRES 100
		double xx[XRES];
		double y1[XRES];

		for (e = 0; e < p->di; e++) {
			printf("Input curve channel %d\n",e);
			for (i = 0; i < XRES; i++) {
				double x;
				x = i/(double)(XRES-1);
				xx[i] = x = x * (p->in_max[e] - p->in_min[e]) + p->in_min[e];
				y1[i] = icxSTransFunc(p->v + p->in_offs[e], p->iluord[e], x,  
			                                    p->in_min[e], p->in_max[e]);
			}
			do_plot(xx,y1,NULL,NULL,XRES);

			printf("Input sub-grid curve channel %d\n",e);
			for (i = 0; i < XRES; i++) {
				double x;
				x = i/(double)(XRES-1);
				xx[i] = x = x * (p->in_max[e] - p->in_min[e]) + p->in_min[e];
				y1[i] = xfit3_subcurve(p, x, e);
			}
			do_plot(xx,y1,NULL,NULL,XRES);
		}
	
		for (f = 0; f < p->fdi; f++) {
			printf("Output curve channel %d\n",f);
			for (i = 0; i < XRES; i++) {
				double x;
				x = i/(double)(XRES-1);
				xx[i] = x = x * (p->out_max[f] - p->out_min[f]) + p->out_min[f];
				y1[i] = p->outcurve(p, x, f);
			}
			do_plot(xx,y1,NULL,NULL,XRES);
		}
	}
#endif /* DEBUG_PLOT */

	return 0;
}

/* We're done with an xfit3 */
static void xfit3_del(xfit3 *p) {
	int e;
	if (p->v != NULL)
		free(p->v);
	if (p->wv != NULL)
		free(p->wv);
	if (p->sa != NULL)
		free(p->sa);
	for (e = 0; e < p->di; e++) {
		if (p->iaxs[e] != NULL)
			p->iaxs[e]->del(p->iaxs[e]);
	}
	if (p->uerrv != NULL)
		free(p->uerrv);
	if (p->ctp != NULL) {
		for (e = 0; e < p->di; e++) {
			if (p->ctp[e] != NULL) {
				int i;
				for (i = 0; i < (p->gres[e]-1); i++) {
					xfit3_ctp *ctp = &p->ctp[e][i];
					if (ctp->lup != NULL)
						free(ctp->lup);
					if (ctp->lupd != NULL)
						free(ctp->lupd);
				}
				free(p->ctp[e]);
			}
		}
		free(p->ctp);
	}
	free(p);
}

/* Create a transform fitting object */
/* return NULL on error */
xfit3 *new_xfit3(
) {
	xfit3 *p;

	if ((p = (xfit3 *)calloc(1, sizeof(xfit3))) == NULL) {
		return NULL;
	}

	/* Set method pointers */
	p->fit         = xfit3_fit;
	p->incurve     = xfit3_incurve;
	p->invincurve  = xfit3_invincurve;
	p->outcurve    = xfit3_outcurve;
	p->invoutcurve = xfit3_invoutcurve;
	p->del         = xfit3_del;

	return p;
}



