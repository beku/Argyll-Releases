
/* 
 * Argyll Color Correction System
 * Multi-dimensional regularized splines,
 * grid cell size optimisation support.
 *
 * Author: Graeme W. Gill
 * Date:   2007/5/23
 *
 * Copyright 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* This file contains suplimentary code for a rspl */
/* that allows it to be used as a source of optimization data */
/* for determining the near optimal lower resolution grid spacing */
/* for a given input axis. The rspl being used will normally have */
/* been set to have a high resolution along the axis of interest, */
/* and a low resolution for the other axes. */
/* The caller will be using the functions here to determine */
/* the best locations for the lower resolution grid nodes along */
/* the axis of interest. */

/* TTBD:
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>

#include "rspl_imp.h"
#include "numlib.h"
#include "counters.h"	/* Counter macros */

extern void error(char *fmt, ...), warning(char *fmt, ...);

#undef DEBUG

#undef NEVER
#define ALWAYS

/* Convention is to use:
   i to index grid points u.a
   n to index data points d.a
   e to index position dimension di
   f to index output function dimension fdi
   j misc and cube corners
   k misc
 */

/* Completely free data associated with gcso */
void free_gcso(rspl *s) {

	if (s->gcso.sy != NULL) {
		free(s->gcso.sy);
		s->gcso.sy = NULL;
	}
	if (s->gcso.sxy != NULL) {
		free(s->gcso.sxy);
		s->gcso.sxy = NULL;
	}
	if (s->gcso.de != NULL) {
		free_dmatrix(s->gcso.de, 0, s->g.res[s->gcso.xi]-1, 0, s->g.res[s->gcso.xi]-1);
		s->gcso.de = NULL;
	}
	if (s->gcso.sde != NULL) {
		free_dmatrix(s->gcso.sde, 0, s->g.res[s->gcso.xi]-1, 0, s->g.res[s->gcso.xi]-1);
		s->gcso.sde = NULL;
	}
}

/* Setup the gcso to return particular error etimates */
/* Return nz on error */
int gcso_setup(
	struct _rspl *s,		/* this */
	int xi,					/* index of input axis of interest */
	int inde,				/* Compute error in input space rather than output space */
	void *cntx,				/* Context of callbacks */
	void (*ocurv)(void *cntx, double *out, double *in),
							/* callback to convert fdi values from out' to out values */
	void (*iocurv)(void *cntx, double *out, double *in),
							/* callback to convert fdi values from out to out' values */ 
	void *cntx2,			/* Context of to_de2 callback */
	double (*to_de2)(void *cntx2, double *in1, double *in2)	
							/* callback to convert error values to delta E squared */
) {
	int di = s->di;
	int fdi = s->fdi;
	int i, j, ee, e, f, n;
	int gc[MXRI];
	float *sdv, *ssy, *ssxy;

	free_gcso(s);			/* Free any exsting data */

	s->gcso.xi = xi;
	s->gcso.xires = s->g.res[xi];
	s->gcso.inde = inde;
	s->gcso.cntx = cntx;
	s->gcso.ocurv = ocurv;
	s->gcso.iocurv = iocurv;
	s->gcso.cntx2 = cntx2;
	s->gcso.to_de2 = to_de2;

	/* Compute index coordinate increments into linear grid for each dimension */
	/* ie. 1, gres, gres^2, gres^3 */
	for (e = 0; e < di; e++) {
		s->gcso.ci[e]  = s->g.ci[e];			/* Duplicate of g.ci[] */
		s->gcso.fci[e] = s->g.ci[e] * s->fdi;	/* In floats */
	}

	/* Allocate arrays for holding pre-comuted linear regression sums */
	if ((s->gcso.sy = (float *)malloc(sizeof(float) * s->g.no * s->fdi)) == NULL) {
		warning("rspl malloc failed - gcso grid points");
		return 1;
	}
	if ((s->gcso.sxy = (float *)malloc(sizeof(float) * s->g.no * s->fdi)) == NULL) {
		warning("rspl malloc failed - gcso grid points");
		return 1;
	}

	/* Compute linear regression sums. Sums are in out' (linear interp.) space */
	/* Do all lines in direction xi */
	/* Init counter and start of line pointers */
	for (e = 0; e < di; e++)
		gc[e] = 0;	/* init coords */
	sdv = s->g.a;				/* Start pointers */
	ssy = s->gcso.sy;
	ssxy = s->gcso.sxy;

//printf("~1 sy base = 0x%x, increment = 0x%x\n",s->gcso.sy,s->gcso.fci[xi]);

	for (e = 0; e < di;) {		/* Until index completed */
		float *dv, *sy, *sxy;
		double vsy[MXDO], vsxy[MXDO];

		dv = sdv;
		sy = ssy;
		sxy = ssxy;
		for (f = 0; f < fdi; f++)
			vsy[f] = vsxy[f] = 0.0;

		for (n = 0; n < s->gcso.xires; n++) {
			double vx, odv[MXDO], idv[MXDO];
			vx = n/(s->gcso.xires-1.0);

			for (f = 0; f < fdi; f++)
				idv[f] = dv[f];
			s->gcso.iocurv(s->gcso.cntx, odv, idv);

			for (f = 0; f < fdi; f++) {
				sy[f] = vsy[f] += odv[f]; 
				sxy[f] = vsxy[f] += vx * odv[f]; 
			}
			dv  += s->g.fci[xi];
			sy  += s->gcso.fci[xi];
			sxy += s->gcso.fci[xi];
		}

		/* Increment index and start points orthogonal to xi */
		for (e = 0; e < di; e++) {
			if (e == xi) 		/* Don't go in direction xi */
				continue;
			gc[e]++;
			sdv  += s->g.fci[e];
			ssy  += s->gcso.fci[e];
			ssxy += s->gcso.fci[e];
			if (gc[e] < s->g.res[e])
				break;	/* No carry */
			gc[e] -= s->g.res[e];			/* Reset coord */
			sdv  -= s->g.res[e] * s->g.fci[e];
			ssy  -= s->g.res[e] * s->gcso.fci[e];
			ssxy -= s->g.res[e] * s->gcso.fci[e];
		}
	}

	/* Allocate and init cache of delta E's for a given quantized span. */
	/* Only the triangle j >= i is populated */
	s->gcso.de = dmatrix(0, s->gcso.xires-1, 0, s->gcso.xires-1);
	for (i = 0; i < s->gcso.xires; i++) {
		for (j = 0; j < s->gcso.xires; j++) {
			s->gcso.de[i][j] = -1.0;		/* Uninitialized value */
		}
	}
	/* Span delta E cache */
	s->gcso.sde = dmatrix(0, s->gcso.xires-1, 0, s->gcso.xires-1);

	return 0;
}

#define square( xx )  ((xx) * (xx))

/* Compute the estimated error for the given integer grid span. */
static double gcso_fwd_err_q(
	struct _rspl *s,		/* this */
	double *pspande,			/* if not NULL return the span width delta E */
							/* as measured by delta E of the span */
	int xx0, int xx1		/* Span range, 0 - gres[xi]-1, xx0 - gres[xi]-1 */
) {
	int di = s->di;
	int fdi = s->fdi;
	int e, f, n;
	int gc[MXRI];		/* Grid index */
	float *sdv, *ssy, *ssxy;
	double S, Sx, Sy, Sxx, Sxy;		/* Regression parameters */
	int xi = s->gcso.xi;
	double rv = 0.0;
	double spande = 0.0;			/* Average span delta E */
	int nspans = 0;					/* Number of spans to average */

//printf("\n~1 gcso_fwd_err_q called with %d %d\n",xx0,xx1);

	/* Do all lines in direction xi */
	/* Init counter and start of line pointers */
	for (e = 0; e < di; e++)
		gc[e] = 0;				/* init coords */
	sdv = s->g.a;				/* Start pointers */
	ssy = s->gcso.sy;
	ssxy = s->gcso.sxy;

	S = xx1 - xx0 + 1;			/* Number of samples */

	if (S <= 2)
		return 0.0;				/* Will be zero error */

	Sx = (xx1+1)/2.0 * xx1/(s->gcso.xires-1.0);
	if (xx0 > 0)
	   Sx -= (xx0)/2.0 * (xx0-1)/(s->gcso.xires-1.0);
	Sxx = ((xx1+1)/3.0 + (xx1+1)/(6.0 * xx1)) * square(xx1/(s->gcso.xires-1.0));
	if (xx0 > 1)
		Sxx -= (xx0/3.0 + xx0/(6.0 * (xx0-1.0))) * square((xx0-1)/(s->gcso.xires-1.0));
	
	for (e = 0; e < di;) {			/* Until index completed */
		int so1, so0;				/* Sum table offsets */
		float *dv;					/* Value and sum table bases */
		double a[MXDO], b[MXDO];	/* Offset + slope */

		so1 = xx1 * s->gcso.fci[xi];
		so0 = (xx0-1) * s->gcso.fci[xi];

		for (f = 0; f < fdi; f++) {
			double d;

			Sy = *(ssy + so1 + f);
			if (xx0 > 0)
				Sy -= *(ssy + so0 + f);

			Sxy = *(ssxy + so1 + f);
			if (xx0 > 0)
				Sxy -= *(ssxy + so0 + f);
//printf("~1 f = %d, Sx = %f, Sxx = %f, Sy = %f, f, Sxy = %f\n",f, Sx,Sxx,Sy,Sxy);

			/* Compute linear regression line params */
			d = S * Sxx - Sx * Sx;
			a[f] = (Sxx * Sy - Sx * Sxy)/d;
			b[f] = (S * Sxy - Sx * Sy)/d;
		}

#ifdef NEVER
// ~~~~ check
{
		double _Sx[MXDO], _Sxx[MXDO], _Sy[MXDO], _Sxy[MXDO];
		float *sy, *sxy;
	
		for (f = 0; f < fdi; f++)
			_Sx[f] = _Sxx[f] = _Sy[f] = _Sxy[f] = 0.0;

		dv = sdv + xx0 * s->g.fci[xi];
		sy = ssy + xx0 * s->gcso.fci[xi];
		sxy = ssxy + xx0 * s->gcso.fci[xi];
		for (n = xx0; n <= xx1; n++) {
			double vx;

			vx = n/(s->gcso.xires-1.0);
			for (f = 0; f < fdi; f++) {
//				printf("n = %d, vx = %f, f = %d, &sy = 0x%x, sy = %f, sxy = %f\n",n,vx,f,&sy[f],sy[f],sxy[f]);
				if (n >= xx0 && n <= xx1) {
					_Sx[f] += vx;
					_Sxx[f] += vx * vx;
					_Sy[f] += dv[f];
					_Sxy[f] += vx * dv[f];
				}
			}
			dv  += s->g.fci[xi];
			sy  += s->gcso.fci[xi];
			sxy += s->gcso.fci[xi];
		}
		for (f = 0; f < fdi; f++)
			printf("~1 f = %d, _Sx = %f, _Sxx = %f, _Sy = %f, _Sxy = %f\n",f,_Sx[f],_Sxx[f],_Sy[f],_Sxy[f]);
}
// ~~~~~~~~~~~~~~~~~~`
#endif /* NEVER */

		/* compute span DE here by taking the delta E between the xx0 and xx1 values */
		if (s->gcso.inde) {
			double iv0[MXDO], iv1[MXDO];

			/* Compute input values from grid indexes */
			for (e = 0; e < di; e++)
				iv0[e] = iv1[e] = s->g.l[e] + gc[e] * s->g.w[e];
			iv0[xi] = s->g.l[e] + xx0 * s->g.w[e];
			iv1[xi] = s->g.l[e] + xx1 * s->g.w[e];

			spande += s->gcso.to_de2(s->gcso.cntx2, iv0, iv1);
			nspans++;
		} else {
			double dv0[MXDO], dv1[MXDO];
			
			dv = sdv + xx0 * s->g.fci[xi];
			for (f = 0; f < fdi; f++)
			dv0[f] = dv[f];
			dv = sdv + xx1 * s->g.fci[xi];
			for (f = 0; f < fdi; f++)
				dv1[f] = dv[f];

			spande += s->gcso.to_de2(s->gcso.cntx2, dv0, dv1);
			nspans++;
		}

		/* Sum the errors between the linear regression lines */
		/* and the actual rspl values */
		dv = sdv + xx0 * s->g.fci[xi];
		for (n = xx0; n <= xx1; n++) {
			double vx, fdv[MXDO], ldv[MXDO];

			vx = n/(s->gcso.xires-1.0);

			if (s->gcso.inde) {
				for (f = 0; f < fdi; f++)
					ldv[f] = a[f] + b[f] * vx;			/* Linear interp out' value */
				s->gcso.ocurv(s->gcso.cntx, ldv, ldv);	/* Linear interp out value */
				
				/* ~~~99 the idea then is to search the output values for a match */
				/* to ldv[], and then take the delta of that input value to vx. */
				error("~~99 gcso not implemented rspl.gcso.inde mode yet!");

			} else {

				for (f = 0; f < fdi; f++) {
					fdv[f] = dv[f];				/* target out value */
					ldv[f] = a[f] + b[f] * vx;	/* Linear interp out' value */
				}
				s->gcso.ocurv(s->gcso.cntx, ldv, ldv);	/* Linear interp out value */

//printf("~1 vx ldv %f, %f %f %f <=> fdv %f %f %f, de %f\n",
//vx,ldv[0],ldv[1],ldv[2],fdv[0],fdv[1],fdv[2], s->gcso.to_de2(s->gcso.cntx2, ldv, fdv));
				rv += s->gcso.to_de2(s->gcso.cntx2, ldv, fdv);	/* Delta to output value */
			}
			dv += s->g.fci[xi];
		}

		/* Increment index and start points orthogonal to xi */
		for (e = 0; e < di; e++) {
			if (e == xi) 		/* Don't go in direction xi */
				continue;
			gc[e]++;
			sdv  += s->g.fci[e];
			ssy  += s->gcso.fci[e];
			ssxy += s->gcso.fci[e];
			if (gc[e] < s->g.res[e])
				break;	/* No carry */
			gc[e] -= s->g.res[e];			/* Reset coord */
			sdv  -= s->g.res[e] * s->g.fci[e];
			ssy  -= s->g.res[e] * s->gcso.fci[e];
			ssxy -= s->g.res[e] * s->gcso.fci[e];
		}
	}
	/* Normalize so that the error is the average of the orthogonal grids */
	/* integrated error in the primary direction over the span. */
	/* Various quantities cancel out to give this: */
	rv /= (double)s->g.no;
	rv *= 100.0;				/* Arbitrary scaling factor */

	spande /= (double)nspans;

//printf("~1 gcso_fwd_err_q returning (%f) %f\n",spande,rv);
	if (pspande != NULL)
		*pspande = spande;
	return rv;
}
#undef square

/* Return the estimated error for the given grid span on the axis of interest */
/* Note that this is for quick initial input curve creation */
/* for placement, and doesn't take into account any existing */
/* input curve in the linear regression fitting or interpolation */
/* use for the error estimate. */
static double gcso_fwd_err(
	struct _rspl *s,		/* this */
	double *spande,			/* if not NULL return the span width delta E */
							/* as measured by delta E of the span */
	double x0, double x1	/* Span range in input space, x0 < x1 */
) {
	double rv;
	int xx0, xx1;				/* Base indexes */
	double w0, w1;				/* Blend weights */
	double v00, v01, v10, v11;	/* Corner values */
	double sv00, sv01, sv10, sv11;	/* Corner span width DE values */

	if (x1 < x0) {				/* oops */
		double tt = x0;
		x0 = x1; x1 = tt;
	}

//printf("~1 Input %f %f\n",x0,x1);

#ifdef NEVER
	/* Interpolated result */
	/* (This seems to be worse than using quantized) */

	x0 = (x0 - s->g.l[s->gcso.xi])/s->g.w[s->gcso.xi];
	xx0 = (int)floor(x0);			/* Coordinate */
	if (xx0 > (s->gcso.xires-2))
		xx0 = (s->gcso.xires-2);
	w0 = 1.0 - x0 + (double)xx0;	/* weight */

	x1 = (x1 - s->g.l[s->gcso.xi])/s->g.w[s->gcso.xi];
	xx1 = (int)floor(x1);			/* Coordinate */
	if (xx1 > (s->gcso.xires-2))
		xx1 = (s->gcso.xires-2);
	w1 = 1.0 - x1 + (double)xx1;	/* weight */

//printf("~1 base %d weight %f, base %d weight %f\n",xx0,w0,xx1,w1);

	/* Lookup corner weights */
	if ((v00 = s->gcso.de[xx0][xx1]) < 0.0)
		v00  = s->gcso.de[xx0][xx1] = gcso_fwd_err_q(s, &s->gcso.sde[xx0][xx1], xx0, xx1);
	if ((v01 = s->gcso.de[xx0+1][xx1]) < 0.0)
		v01  = s->gcso.de[xx0+1][xx1] = gcso_fwd_err_q(s, &s->gcso.sde[xx0+1][xx1], xx0+1, xx1);
	if ((v10 = s->gcso.de[xx0][xx1+1]) < 0.0)
		v10  = s->gcso.de[xx0][xx1+1] = gcso_fwd_err_q(s, &s->gcso.sde[xx0][xx1+1], xx0, xx1+1);
	if ((v11 = s->gcso.de[xx0+1][xx1+1]) < 0.0)
		v11  = s->gcso.de[xx0+1][xx1+1] = gcso_fwd_err_q(s, &s->gcso.sde[xx0+1][xx1+1], xx0+1, xx1+1);
	sv00 = s->gcso.sde[xx0][xx1];
	sv01 = s->gcso.sde[xx0+1][xx1];
	sv10 = s->gcso.sde[xx0][xx1+1];
	sv11 = s->gcso.sde[xx0+1][xx1+1];
	
//printf("~1 Corner values %f %f %f %f\n",v00,v01,v10,v11);

	/* Compute interpolated error value */
	rv =        w0  *        w1  * v00
	   + (1.0 - w0) *        w1  * v01
	   +        w0  * (1.0 - w1) * v10
	   + (1.0 - w0) * (1.0 - w1) * v11;

	*spande =        w0  *        w1  * sv00
	        + (1.0 - w0) *        w1  * sv01
	        +        w0  * (1.0 - w1) * sv10
	        + (1.0 - w0) * (1.0 - w1) * sv11;
#else
	/* Quantized result */
	xx0 = (int)((x0 - s->g.l[s->gcso.xi])/s->g.w[s->gcso.xi] + 0.5);
	xx1 = (int)((x1 - s->g.l[s->gcso.xi])/s->g.w[s->gcso.xi] + 0.5);
	if ((rv = s->gcso.de[xx0][xx1]) < 0.0)
		rv  = s->gcso.de[xx0][xx1] = gcso_fwd_err_q(s, &s->gcso.sde[xx0][xx1], xx0, xx1);
//	s->gcso.de[xx0][xx1] = -1.0;			/* Disable caching */
	*spande = s->gcso.sde[xx0][xx1];
#endif

//printf("~1 Final value %f\n",rv);
	return rv;
}

/* Initialise the gcso data */ 
/* return nz on error */
int init_gcso(rspl *s) {

	s->gcso_setup   = gcso_setup;
	s->gcso_fwd_err = gcso_fwd_err;

	return 0;
}










