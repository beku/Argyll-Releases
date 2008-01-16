

/*
 * International Color Consortium color transform expanded support
 *
 * Author:  Graeme W. Gill
 * Date:    2/7/00
 * Version: 1.00
 *
 * Copyright 2000 - 2007 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * Based on the xlut.c code.
 */

/*
 * This module provides curve and matrix fitting functionality used to
 * create per channel input and/or output curves for clut profiles,
 * and optionally creates the rspl to fit the data as well.
 * 
 * The approach used is to initialy creating input and output shaper
 * curves that minimize the overall delta E of the test point set,
 * when a linear matrix is substituted for a clut is used.
 * (This is the same as pre-V0.70 approach. )
 *
 * The residual error between the test point set and the shaper/matrix/shaper
 * model is then computed, and mapped against each input channel, and
 * a positioning mapping curve created that aims to map rspl grid
 * locations more densly where residual errors are high, and more
 * sparsely where they are low.
 * The input shaper and positioning curves are then combined together,
 * so that the positioning curve determines which input values map
 * to rspl grid points, and the shaper curve determines the mapping
 * between the grid points. The width of the grid cells is computed
 * in the shaper mapped input space, and then fed to the rspl fitting
 * code, so that its smoothness evaluation function can take account
 * of the non-uniform spacing of the grid points.
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
#include "xfit.h"
#include "sort.h"


#undef DEBUG 			/* Verbose debug information */
#undef DEBUG_PLOT 		/* Plot in & out curves */
#undef SPECIAL_FORCE	/* Check rspl nodes against linear XYZ model */
#undef SPECIAL_FORCE_GAMMA		/* Force correct gamma shaper curves */
#define SPECIAL_TEST_GAMMA		/* Use gamma on reference model */
#undef SPECIAL_TEST_LAB	

#define USE_CIE94_DE	/* Use CIE94 delta E measure when creating in/out curves */
#undef NODDV			/* Use slow non d/dv powell */

#define POWTOL 5e-4			/* Shaper Powell optimiser tollerance */
#define MAXITS 2000			/* Shaper number of itterations before giving up */
#define PDDEL  1e-6			/* Fake partial derivative del */

/* Weights for shaper in/out curve parameters, to minimise unconstrained "wiggles" */
#define SHAPE_WEIGHT	1.0		/* Overal shaper weight contribution - err on side of smoothness */
#define SHAPE_HW01		0.2		/* 0 & 1 harmonic weights */
#define SHAPE_HBREAK    4		/* Harmonic that has HWBR */
#define SHAPE_HWBR  	0.8		/* Base weight of harmonics HBREAK up */
#define SHAPE_HWINC  	0.4		/* Increase in weight for each harmonic above HWBR */

/* Weights for the positioning curve parameters */
#define PSHAPE_MINE 0.02		/* Minum background residual error level */

/*
 * TTBD:
 *
 */

/* - - - - - - - - - - - - - - - - - */
/* Lookup a value though an input position curve */
static double xfit_poscurve(xfit *p, double in, int chan) {
	double rv = in;
	if (p->tcomb & oc_p)
		rv = icxSTransFunc(p->v + p->pos_offs[chan], p->iluord[chan], rv,  
		                       p->in_min[chan], p->in_max[chan]);
	return rv;
}

/* Lookup di values though the input position curves */
static void xfit_poscurves(xfit *p, double *out, double *in) {
	int e;
	
	for (e = 0; e < p->di; e++) {
		double val = in[e];
		if (p->tcomb & oc_i)
			val = icxSTransFunc(p->v + p->pos_offs[e], p->iluord[e], val,  
			                       p->in_min[e], p->in_max[e]);
		out[e] = val;
	}
}

/* Inverse Lookup a value though an input position curve */
static double xfit_invposcurve(xfit *p, double in, int chan) {
	double rv = in;
	if (p->tcomb & oc_i)
		rv = icxInvSTransFunc(p->v + p->pos_offs[chan], p->iluord[chan], rv,  
		                       p->in_min[chan], p->in_max[chan]);
	return rv;
}

/* Inverse Lookup di values though input position curves */
static void xfit_invposcurves(xfit *p, double *out, double *in) {
	int e;
	
	for (e = 0; e < p->di; e++) {
		double val = in[e];
		if (p->tcomb & oc_i)
			val = icxInvSTransFunc(p->v + p->pos_offs[e], p->iluord[e], val,  
			                       p->in_min[e], p->in_max[e]);
		out[e] = val;
	}
}

/* - - - - - - - - - - - - - - - - - */
/* Lookup a value though input shape curve */
static double xfit_shpcurve(xfit *p, double in, int chan) {
	double rv = in;

	if (p->tcomb & oc_i) {
#ifdef SPECIAL_FORCE_GAMMA	
		double gam;
		if (chan == 0)
			gam = 1.9;
		else if (chan == 1)
			gam = 2.0;
		else if (chan == 2)
			gam = 2.1;
		else
			gam = 1.0;
		rv = pow(in, gam);
#else
		rv = icxSTransFunc(p->v + p->shp_offs[chan], p->iluord[chan], rv,  
			                       p->in_min[chan], p->in_max[chan]);
#endif
	}
	return rv;
}

/* Lookup a value though shape curves */
static void xfit_shpcurves(xfit *p, double *out, double *in) {
	int e;
	
	for (e = 0; e < p->di; e++) {
		double val = in[e];
		if (p->tcomb & oc_i)
			val = icxSTransFunc(p->v + p->shp_offs[e], p->iluord[e], val,  
			                       p->in_min[e], p->in_max[e]);
		out[e] = val;
	}
}

/* Inverse Lookup a value though a shape curve */
static double xfit_invshpcurve(xfit *p, double in, int chan) {
	double rv = in;

	if (p->tcomb & oc_i) {
#ifdef SPECIAL_FORCE_GAMMA	
		double gam;
		if (chan == 0)
			gam = 1.9;
		else if (chan == 1)
			gam = 2.0;
		else if (chan == 2)
			gam = 2.1;
		else
			gam = 1.0;
		rv = pow(rv, 1.0/gam);
#else
		rv = icxInvSTransFunc(p->v + p->shp_offs[chan], p->iluord[chan], rv,  
		                       p->in_min[chan], p->in_max[chan]);
#endif
	}
	return rv;
}

/* Inverse Lookup a value though shape curves */
static void xfit_invshpcurves(xfit *p, double *out, double *in) {
	int e;
	
	for (e = 0; e < p->di; e++) {
		double val = in[e];
		if (p->tcomb & oc_i)
			val = icxInvSTransFunc(p->v + p->shp_offs[e], p->iluord[e], val,  
			                       p->in_min[e], p->in_max[e]);
		out[e] = val;
	}
}

/* - - - - - - - - - - - - - - - - - */

/* Lookup values through the shaper/matrix/shaper model */
static void xfit_shmatsh(xfit *p, double *out, double *in) {
	double tin[MXDI];
	int e, f;
	
	for (e = 0; e < p->di; e++)
		tin[e] = icxSTransFunc(p->v + p->shp_offs[e], p->iluord[e], in[e],  
			                       p->in_min[e], p->in_max[e]);

	icxCubeInterp(p->v + p->mat_off, p->fdi, p->di, out, tin);

	for (f = 0; f < p->fdi; f++)
		out[f] = icxSTransFunc(p->v + p->out_offs[f], p->oluord[f], out[f],
			                       p->out_min[f], p->out_max[f]);
}

/* - - - - - - - - - - - - - - - - - - - - */
/* Combined input positioning & shaper transfer curve functions */

int db = 0;

/* Lookup a value though the input positioning and shaper curves */
static double xfit_inpscurve(xfit *p, double in, int chan) {
	double rv;
	/* Just shaper curve */
	if ((p->tcomb & oc_ip) == oc_i) {
		rv = icxSTransFunc(p->v + p->shp_offs[chan], p->iluord[chan], in,  
			                       p->in_min[chan], p->in_max[chan]);

	/* shaper and positioning */
	} else if ((p->tcomb & oc_ip) == oc_ip) {
		double nin, npind;		/* normalized in, normalized position in' */
		int six;				/* Span index */
		double npind0, npind1;	/* normalized position in' span values */
		double npin0, npin1;	/* normalized position in span values */
		double nsind0, nsind1;	/* normalized shaper in' span values */
		double nsind;			/* normalised shaper in' value */

		/* Normalize */
		nin = (in - p->in_min[chan])/(p->in_max[chan] - p->in_min[chan]);

if (db) printf("\n~1 inpscurve: cha %d, input value %f, norm %f\n",chan,in,nin);

		/* Locate the span the input point will be in after positioning lookup */
		npind = icxTransFunc(p->v + p->pos_offs[chan], p->iluord[chan], nin); 

		/* Quantize position space value to grid */
		six = (int)floor(npind * (p->gres[chan]-1.0));
		if (six > (p->gres[chan]-2))
			six = (p->gres[chan]-2);

		/* Compute span position position in' values */
		npind0 = six / (p->gres[chan]-1.0);
		npind1 = (six + 1.0) / (p->gres[chan]-1.0);

if (db) printf("~1 npind %f, six %d, npind0 %f, npind1 %f\n",npind,six,npind0,npind1);

		/* Compute span in values */
		npin0 = icxInvTransFunc(p->v + p->pos_offs[chan], p->iluord[chan], npind0); 
		npin1 = icxInvTransFunc(p->v + p->pos_offs[chan], p->iluord[chan], npind1); 
	
if (db) printf("~1 npin0 %f, npin1 %f\n",npin0,npin1);

		/* Compute shaper space values of in' and spane */
#ifdef NEVER
		nsind = icxTransFunc(p->v + p->shp_offs[chan], p->iluord[chan], nin); 
		nsind0 = icxTransFunc(p->v + p->shp_offs[chan], p->iluord[chan], npin0); 
		nsind1 = icxTransFunc(p->v + p->shp_offs[chan], p->iluord[chan], npin1); 
#else
		nsind = xfit_shpcurve(p, nin, chan); 
		nsind0 = xfit_shpcurve(p, npin0, chan); 
		nsind1 = xfit_shpcurve(p, npin1, chan); 
#endif

if (db) printf("~1 nsind %f, nsind0 %f, nsind1 %f\n",nsind,nsind0,nsind1);

		/* Offset and scale shaper in' value to match position span */
		rv = (nsind - nsind0)/(nsind1 - nsind0) * (npind1 - npind0) + npind0;

if (db) printf("~1 scale offset ind %f\n",rv);

		/* de-normalize */
		rv = rv * (p->in_max[chan] - p->in_min[chan]) + p->in_min[chan];
if (db) printf("~1 returning \n",rv);
	
	/* Just positioning curve */
	} else if ((p->tcomb & oc_ip) == oc_p) {
		rv = icxSTransFunc(p->v + p->pos_offs[chan], p->iluord[chan], in,  
				                       p->in_min[chan], p->in_max[chan]);
	} else {
		rv = in;
	}
	return rv;
}

/* Lookup di values though the input positioning and shaper curves */
static void xfit_inpscurves(xfit *p, double *out, double *in) {
	int e;

	for (e = 0; e < p->di; e++)
		out[e] = xfit_inpscurve(p, in[e], e);
}

/* Inverse Lookup a value though the input positioning and shaper curves */
static double xfit_invinpscurve(xfit *p, double in, int chan) {
	double rv;
	/* Just shaper curve */
	if ((p->tcomb & oc_ip) == oc_i) {
		rv = icxInvSTransFunc(p->v + p->shp_offs[chan], p->iluord[chan], in,  
			                       p->in_min[chan], p->in_max[chan]);

	/* shaper and positioning */
	} else if ((p->tcomb & oc_ip) == oc_ip) {
		double nind, nin;		/* normalized in', normalized in */
		int six;				/* Span index */
		double npind0, npind1;	/* normalized position in' span values */
		double npin0, npin1;	/* normalized position in span values */
		double nsind0, nsind1;	/* normalized shaper in' span values */
		double nsind;			/* normalized shaper in' value */

		/* Normalize */
		nind = (in - p->in_min[chan])/(p->in_max[chan] - p->in_min[chan]);

if (db) printf("\n~1 invinpscurve: cha %d, input value %f, norm %f\n",chan,in,nind);

		/* Quantize to grid */
		six = (int)floor(nind * (p->gres[chan]-1.0));
		if (six > (p->gres[chan]-2))
			six = (p->gres[chan]-2);

		/* Compute span in' values */
		npind0 = six / (p->gres[chan]-1.0);
		npind1 = (six + 1.0) / (p->gres[chan]-1.0);

if (db) printf("~1 six %d, npind0 %f, npind1 %f\n",six,npind0,npind1);

		/* Lookup span in values through position curve */
		npin0 = icxInvTransFunc(p->v + p->pos_offs[chan], p->iluord[chan], npind0); 
		npin1 = icxInvTransFunc(p->v + p->pos_offs[chan], p->iluord[chan], npind1); 
	
if (db) printf("~1 npin0 %f, npin1 %f\n",npin0,npin1);

		/* Compute span shaper in' values */
#ifdef NEVER
		nsind0 = icxTransFunc(p->v + p->shp_offs[chan], p->iluord[chan], npin0); 
		nsind1 = icxTransFunc(p->v + p->shp_offs[chan], p->iluord[chan], npin1); 
#else
		nsind0 = xfit_shpcurve(p, npin0, chan); 
		nsind1 = xfit_shpcurve(p, npin1, chan); 
#endif

		/* Offset and scale position in' value to match shaper span */
		nsind = (nind - npind0)/(npind1 - npind0) * (nsind1 - nsind0) + nsind0;

if (db) printf("~1 nsind %f, nsind0 %f, nsind1 %f\n",nsind,nsind0,nsind1);

		/* Invert through shaper curve */
#ifdef NEVER
		nin = icxInvTransFunc(p->v + p->shp_offs[chan], p->iluord[chan], nsind); 
#else
		nin = xfit_invshpcurve(p, nsind, chan); 
#endif

		/* de-normalize */
		rv = nin * (p->in_max[chan] - p->in_min[chan]) + p->in_min[chan];

if (db) printf("\n~1 nin = %f, returning %f\n",nin,rv);
	
	/* Just positioning curve */
	} else if ((p->tcomb & oc_ip) == oc_p) {
		rv = icxInvSTransFunc(p->v + p->pos_offs[chan], p->iluord[chan], in,  
			                       p->in_min[chan], p->in_max[chan]);
	} else {
		rv = in;
	}
	return rv;
}

#ifdef NEVER
/* Check that inverse is working */
static double _xfit_inpscurve(xfit *p, double in, int chan) {
	double inv, rv, iinv;

	inv = in;
	rv = xfit_inpscurve(p, in, chan); 
	iinv = xfit_invinpscurve(p, rv, chan);

	if (fabs(in - iinv) > 1e-5)
		warning("xfit_inpscurve check, got %f, should be %f\n",iinv,in);
	
	return rv;
}
#define xfit_inpscurve _xfit_inpscurve
#endif /* NEVER */

/* Inverse Lookup di values though the input positioning and shaper curves */
static void xfit_invinpscurves(xfit *p, double *out, double *in) {
	int e;

	for (e = 0; e < p->di; e++)
		out[e] = xfit_invinpscurve(p, in[e], e);
}

/* - - - - - - - - - - - - - - - - - - - - */
/* Combined output transfer curve functions */

/* Lookup a value though an output curve */
static double xfit_outcurve(xfit *p, double in, int chan) {
	double rv;
	if (p->tcomb & oc_o)
		rv = icxSTransFunc(p->v + p->out_offs[chan], p->oluord[chan], in,  
		                       p->out_min[chan], p->out_max[chan]);
	else
		rv = in;
	return rv;
}

/* Lookup fdi values though the output curves */
static void xfit_outcurves(xfit *p, double *out, double *in) {
	int f;
	
	for (f = 0; f < p->fdi; f++) {
		double val = in[f];
		if (p->tcomb & oc_o)
			val = icxSTransFunc(p->v + p->out_offs[f], p->oluord[f], val,  
			                       p->out_min[f], p->out_max[f]);
		out[f] = val;
	}
}

/* Inverse Lookup a value though an output curve */
static double xfit_invoutcurve(xfit *p, double in, int chan) {
	double rv;
	if (p->tcomb & oc_o)
		rv = icxInvSTransFunc(p->v + p->out_offs[chan], p->oluord[chan], in,  
		                       p->out_min[chan], p->out_max[chan]);
	else
		rv = in;
	return rv;
}

/* Inverse Lookup fdi values though output curves */
static void xfit_invoutcurves(xfit *p, double *out, double *in) {
	int f;
	
	for (f = 0; f < p->fdi; f++) {
		double val = in[f];
		if (p->tcomb & oc_o)
			val = icxInvSTransFunc(p->v + p->out_offs[f], p->oluord[f], val,  
			                       p->out_min[f], p->out_max[f]);
		out[f] = val;
	}
}

/* - - - - - - - - - */

/* Convert an output value from absolute */
/* to relative using the current white point. */
static void xfit_abs_to_rel(xfit *p, double *out, double *in) {
	if (p->flags & XFIT_OUT_WP_REL) {
		if (p->flags & XFIT_OUT_LAB) {
			icmLab2XYZ(&icmD50, out, in);
			icmMulBy3x3(out, p->fromAbs, out);
			icmXYZ2Lab(&icmD50, out, out);
		} else {
			icmMulBy3x3(out, p->fromAbs, in);
		}
	}
}

/* - - - - - - - - - */

/* return a weighting for the magnitude of the in and out */
/* shaping parameters. This is to reduce unconstrained "wiggles" */
static double shapmag(
xfit  *p			/* Base of optimisation structure */
) {
	double tt, w;
	double *b;			/* Base of parameters for this section */
	int di =  p->di;
	int fdi = p->fdi;
	int e, f, k;
	double iparam = 0.0;
	double oparam = 0.0;
	double dd;

	if (p->opt_msk & oc_i) {
		dd = SHAPE_WEIGHT/(double)(di);
		b = p->v + p->shp_off;
		for (e = 0; e < di; e++) {
			for (k = 0; k < p->iluord[e]; k++) {
				if (k <= 1) {
					w = SHAPE_HW01;
				} else if (k <= SHAPE_HBREAK) {
					double bl = (k - 1.0)/(SHAPE_HBREAK - 1.0);
					w = (1.0 - bl) * SHAPE_HW01 + bl * SHAPE_HWBR;
				} else {
					w = SHAPE_HWBR + (k-SHAPE_HBREAK) * SHAPE_HWINC;
				}
				tt = *b++;
				tt *= tt;
				iparam += w * tt;
			}
		}
		iparam *= dd;
	}

	if (p->opt_msk & oc_o) {
		dd = SHAPE_WEIGHT/(double)(fdi);
		b = p->v + p->out_off;
		for (f = 0; f < fdi; f++) {
			for (k = 0; k < p->oluord[f]; k++) {
				if (k <= 1) {
					w = SHAPE_HW01;
				} else if (k <= SHAPE_HBREAK) {
					double bl = (k - 1.0)/(SHAPE_HBREAK - 1.0);
					w = (1.0 - bl) * SHAPE_HW01 + bl * SHAPE_HWBR;
				} else {
					w = SHAPE_HWBR + (k-SHAPE_HBREAK) * SHAPE_HWINC;
				}
				tt = *b++;
				tt *= tt;
				oparam += w * tt;
			}
		}
		oparam *= dd;
	}
	return iparam + oparam;
}

/* return a weighting for the magnitude of the in and out */
/* shaping parameters. This is to reduce unconstrained "wiggles" */
/* Also sum the partial derivative for the parameters involved */
static double dshapmag(
xfit  *p,			/* Base of optimisation structure */
double *dav			/* Sum del's */
) {
	double tt, w;
	double *b, *c;			/* Base of parameters for this section */
	int di =  p->di;
	int fdi = p->fdi;
	int e, f, k;
	double iparam = 0.0;
	double oparam = 0.0;
	double dd;

	if (p->opt_msk & oc_i) {
		dd = SHAPE_WEIGHT/(double)(di);
		b = p->v + p->shp_off;
		c = dav + p->shp_off;
		for (e = 0; e < di; e++) {
			for (k = 0; k < p->iluord[e]; k++) {
				if (k <= 1) {
					w = SHAPE_HW01;
				} else if (k <= SHAPE_HBREAK) {	
					double bl = (k - 1.0)/(SHAPE_HBREAK - 1.0);
					w = (1.0 - bl) * SHAPE_HW01 + bl * SHAPE_HWBR;
				} else {
					w = SHAPE_HWBR + (k-SHAPE_HBREAK) * SHAPE_HWINC;
				}
				tt = *b++;
				*c++ += 2.0 * dd * w * tt;
				tt *= tt;
				iparam += w * tt;
				
			}
		}
		iparam *= dd;
	}

	if (p->opt_msk & oc_o) {
		dd = SHAPE_WEIGHT/(double)(fdi);
		b = p->v + p->out_off;
		c = dav + p->out_off;
		for (f = 0; f < fdi; f++) {
			for (k = 0; k < p->oluord[f]; k++) {
				if (k <= 1) {
					w = SHAPE_HW01;
				} else if (k <= SHAPE_HBREAK) {
					double bl = (k - 1.0)/(SHAPE_HBREAK - 1.0);
					w = (1.0 - bl) * SHAPE_HW01 + bl * SHAPE_HWBR;
				} else {
					w = SHAPE_HWBR + (k-SHAPE_HBREAK) * SHAPE_HWINC;
				}
				tt = *b++;
				*c++ += 2.0 * dd * w * tt;
				tt *= tt;
				oparam += w * tt;
			}
		}
		oparam *= dd;
	}
	return iparam + oparam;
}


/* Shaper+Matrix optimisation function handed to powell() */
static double xfitfunc(void *edata, double *v) {
	xfit *p = (xfit *)edata;
	double tw = 0.0;				/* Total weight */
	double ev = 0.0, rv, smv;
	double tin[MXDI], out[MXDO];
	int di = p->di;
	int fdi = p->fdi;
	int i, e, f;

	/* Copy the parameters being optimised into xfit structure */
	for (i = 0; i < p->opt_cnt; i++) { 
//printf("~1 param %d = %f\n",i,v[i]);
		p->v[p->opt_off + i] = v[i];
	}

	/* For all our data points */
	for (i = 0; i < p->nodp; i++) {
		double del;

		/* Apply input shaper channel curves */
		for (e = 0; e < di; e++)
			tin[e] = icxSTransFunc(p->v + p->shp_offs[e], p->iluord[e], p->rpoints[i].p[e],  
			                       p->in_min[e], p->in_max[e]);

		/* Apply matrix cube interpolation */
		icxCubeInterp(p->v + p->mat_off, fdi, di, out, tin);

		/* Apply output channel curves */
		for (f = 0; f < fdi; f++)
			out[f] = icxSTransFunc(p->v + p->out_offs[f], p->oluord[f], out[f],  
			                       p->out_min[f], p->out_max[f]);

		/* Evaluate the error squared */
		if (p->flags & XFIT_FM_INPUT) {
			double pp[MXDI];
			for (e = 0; e < di; e++)
				pp[e] = p->rpoints[i].p[e];
			for (f = 0; f < fdi; f++) {
				double t1 = p->rpoints[i].v[f] - out[f];	/* Error in output */
	
				/* Create input point offset by equivalent delta to output point */
				/* error, in proportion to the partial derivatives for that output. */
				for (e = 0; e < di; e++)
					pp[e] += t1 * p->piv[i].ide[f][e];
			}
			del = p->to_de2(p->cntx2, pp, p->rpoints[i].p);
		} else {
			del = p->to_de2(p->cntx2, out, p->rpoints[i].v);
		}

		tw += p->rpoints[i].w;
		ev += p->rpoints[i].w * del;
	}

	/* Normalise error to be an average delta E squared */
	ev /= tw;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	smv = shapmag(p);
	rv = ev + smv;

#ifdef DEBUG
printf("~1(sm %f, ev %f)xfitfunc returning %f\n",smv,ev,rv);
#endif

#ifdef NODDV
	if (p->verb)
		printf("."), fflush(stdout);
#endif
	return rv;
}

/* Shaper+Matrix optimisation function with partial derivatives, */
/* handed to conjgrad() */
static double dxfitfunc(void *edata, double *dv, double *v) {
	xfit *p = (xfit *)edata;
	double tw = 0.0;				/* Total weight */
	double ev = 0.0, rv, smv;
	double tin[MXDI], out[MXDO];

	double dav[MXPARMS];				/* Overall del due to del param vals */

	double dtin_iv[MXDI * MXLUORD];		/* Del in itrans out due to del itrans param vals */
	double dmato_mv[1 << MXDI];			/* Del in mat out due to del in matrix param vals */
	double dmato_tin[MXDO * MXDI];		/* Del in mat out due to del in matrix input values */
	double dout_ov[MXDO * MXLUORD];		/* Del in otrans out due to del in otrans param values */
	double dout_mato[MXDO];				/* Del in otrans out due to del in otrans input values */

	double dout_de[2][MXDIDO];			/* Del in DE due to two output values */

	int di = p->di;
	int fdi = p->fdi;
	int i, jj, k, e, ee, f, ff;

	/* Copy the parameters being optimised into xfit structure */
	for (i = 0; i < p->opt_cnt; i++)
		p->v[p->opt_off + i] = v[i];

	/* Zero the accumulated partial derivatives */
	/* We compute deriv for all parameters (not just current optimised) */
	for (i = 0; i < p->tot_cnt; i++)
		dav[i] = 0.0;

	/* For all our data points */
	for (i = 0; i < p->nodp; i++) {
		double del;

		/* Apply input channel curves */
		for (e = 0; e < di; e++)
			tin[e] = icxdpSTransFunc(p->v + p->shp_offs[e], &dtin_iv[p->shp_offs[e] - p->shp_off],
		                         p->iluord[e], p->rpoints[i].p[e], p->in_min[e], p->in_max[e]);

		/* Apply matrix cube interpolation */
		icxdpdiCubeInterp(p->v + p->mat_off, dmato_mv, dmato_tin, fdi, di, out, tin);

		/* Apply output channel curves */
		for (f = 0; f < fdi; f++)
			out[f] = icxdpdiSTransFunc(p->v + p->out_offs[f],
			                           &dout_ov[p->out_offs[f] - p->out_off], &dout_mato[f],
			                           p->oluord[f], out[f], p->out_min[f], p->out_max[f]);


		/* Convert to Delta E and compute pde's into dout_de squared */
		if (p->flags & XFIT_FM_INPUT) {
			double tdout_de[2][MXDIDO];
			double pp[MXDI];
			int n,j,k;
			for (e = 0; e < di; e++)
				pp[e] = p->rpoints[i].p[e];
			for (f = 0; f < fdi; f++) {
				double t1 = p->rpoints[i].v[f] - out[f];	/* Error in output */
	
				/* Create input point offset by equivalent delta to output point */
				/* error, in proportion to the partial derivatives for that output. */
				for (e = 0; e < di; e++)
					pp[e] += t1 * p->piv[i].ide[f][e];
			}
			del = p->to_dde2(p->cntx2, tdout_de, pp, p->rpoints[i].p);

			/* Compute partial derivative (is this correct ??) */
			for (e = 0; e < di; e++) {
				dout_de[0][e] = 0.0;
				for (f = 0; f < fdi; f++) {
					dout_de[0][e] += tdout_de[0][e] * p->piv[i].ide[f][e];
				}
			}
		} else {
			del = p->to_dde2(p->cntx2, dout_de, out, p->rpoints[i].v);
		}

		/* Accumulate total weighted delta E squared */
		tw += p->rpoints[i].w;
		ev += p->rpoints[i].w * del;

		for (ff = 0; ff < fdi; ff++) {				/* Parameter output chanel */
			for (k = 0; k < p->oluord[ff]; k++) {	/* Param within channel */
				double vv = 0.0;
				jj = p->out_offs[ff] - p->out_off + k;	/* Overall output trans param */

				vv += dout_de[0][ff] * dout_ov[jj];
				dav[p->out_off + jj] += p->rpoints[i].w * vv;
			}
		}

		/* Compute and accumulate partial difference values for each parameter value */
		if (p->opt_msk & oc_i) {
			/* Input transfer parameters */
			for (ee = 0; ee < di; ee++) {				/* Parameter input chanel */
				for (k = 0; k < p->iluord[ee]; k++) {	/* Param within channel */
					double vv = 0.0;
					jj = p->shp_offs[ee] - p->shp_off + k;	/* Overall input trans param */

					for (ff = 0; ff < 3; ff++) {		/* Lab channels */
						vv += dout_de[0][ff] * dout_mato[ff]
						    * dmato_tin[ff * di + ee] * dtin_iv[jj];
					}
					dav[p->shp_off + jj] += p->rpoints[i].w * vv;
				}
			}
		}

		/* Matrix parameters */
		if (p->opt_msk & oc_m) {
			for (ff = 0; ff < fdi; ff++) {				/* Parameter output chanel */
				for (ee = 0; ee < (1 << di); ee++) {	/* Matrix input combination chanel */
					double vv = 0.0;
					jj = ff * (1 << di) + ee;			/* Matrix Parameter index */

					vv += dout_de[0][ff] * dout_mato[ff] * dmato_mv[ee];
					dav[p->mat_off + jj] += p->rpoints[i].w * vv;
				}
			}
		}

		if (p->opt_msk & oc_o) {
			/* Output transfer parameters */
			for (ff = 0; ff < fdi; ff++) {				/* Parameter output chanel */
				for (k = 0; k < p->oluord[ff]; k++) {	/* Param within channel */
					double vv = 0.0;
					jj = p->out_offs[ff] - p->out_off + k;	/* Overall output trans param */

					vv += dout_de[0][ff] * dout_ov[jj];
					dav[p->out_off + jj] += p->rpoints[i].w * vv;
				}
			}
		}
	}

	/* Normalise error to be an average delta E squared */
	ev /= tw;
	for (i = 0; i < p->tot_cnt; i++)
		dav[i] /= tw;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	smv = dshapmag(p, dav);
	rv = ev + smv;

	/* Copy the del for parameters being optimised to return array */
	for (i = 0; i < p->opt_cnt; i++)
		dv[i] = dav[p->opt_off + i];

#ifdef DEBUG
printf("~1(sm %f, ev %f)dxfitfunc returning %f\n",smv,ev,rv);
#endif

	if (p->verb)
		printf("."), fflush(stdout);
	return rv;
}

#ifdef NEVER
/* Check partial derivative function within xfitfunc() */

static double _xfitfunc(void *edata, double *v) {
	xfit *p = (xfit *)edata;
	int i;
	double dv[MXPARMS];
	double rv, drv;
	double trv;
	int verb;
	
	rv = xfitfunc(edata, v);
	verb = p->verb;
	p->verb = 0;
	drv = dxfitfunc(edata, dv, v);
	p->verb = verb;

	if (fabs(rv - drv) > 1e-6)
		printf("######## RV MISMATCH is %f should be %f ########\n",rv,drv);

	/* Check each parameter delta */
	for (i = 0; i < p->opt_cnt; i++) {
		double del;

		v[i] += 1e-7;
		trv = xfitfunc(edata, v);
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

#define xfitfunc _xfitfunc
#endif

/* - - - - - - - - - */

/* Output curve symetry optimisation function handed to powell() */
/* Just the order 0 value will be adjusted */
static double symoptfunc(void *edata, double *v) {
	xfit *p = (xfit *)edata;
	double out[1], in[1] = { 0.0 };
	int ch = p->opt_ch;		/* Output channel being adjusted for symetry */
	double rv;

	/* Copy the parameter being tested back into xfit */
	p->v[p->out_offs[ch]] = v[0];
	*out = icxSTransFunc(p->v + p->out_offs[ch], p->oluord[ch], *in,
							   p->out_min[ch], p->out_max[ch]);

	rv = out[0] * out[0];

#ifdef DEBUG
printf("~1symoptfunc returning %f\n",rv);
#endif
	return rv;
}

/* - - - - - - - - - */

/* Set up for an optimisation run: */
/* Figure out parameters being optimised, */
/* copy them to start values, */
/* init and scale the search radius */
static void setup_xfit(
xfit *p,
double *wv,		/* Return parameters to hand to optimiser */
double *sa,		/* Return search radius to hand to optimiser */
double transrad,/* Nominal transfer curve radius, 0.0 - 3.0 */
double matrad	/* Nominal matrix radius, 0.0 - 1.0 */
) {
	int i;
	p->opt_off = -1;
	p->opt_cnt = 0;
	
	if (p->opt_msk & oc_i) {
		if (p->opt_off < 0)
			p->opt_off = p->shp_off;
		p->opt_cnt += p->shp_cnt;
	
		for (i = 0; i < p->shp_cnt; i++) {
			*wv++ = p->v[p->shp_off + i];
			*sa++ = transrad;
		}
	}
	if (p->opt_msk & oc_m) {
		if (p->opt_off < 0)
			p->opt_off = p->mat_off;
		p->opt_cnt += p->mat_cnt;

		for (i = 0; i < p->mat_cnt; i++) {
			*wv++ = p->v[p->mat_off + i];
			*sa++ = matrad;
		}
	}
	if (p->opt_msk & oc_o) {
		if (p->opt_off < 0)
			p->opt_off = p->out_off;
		p->opt_cnt += p->out_cnt;

		for (i = 0; i < p->out_cnt; i++) {
			*wv++ = p->v[p->out_off + i];
			*sa++ = transrad;
		}
	}
	if (p->opt_cnt > MXPARMS)
		error("setup_xfit: asert, %d exceeded MXPARMS %d",p->opt_cnt,MXPARMS);

#ifdef NEVER
printf("~1 opt_msk = 0x%x\n",p->opt_msk);
printf("~1 opt_off = %d\n",p->opt_off);
printf("~1 opt_cnt = %d\n\n",p->opt_cnt);
#endif /* NEVER */
}

#ifdef DEBUG
/* Diagnostic */
static void dump_xfit(
xfit *p
) {
	int i, e, f;
	double *b;			/* Base of parameters for this section */
	int di, fdi;
	di   = p->di;
	fdi  = p->fdi;

	/* Input positioning curve */
	b = p->v + p->pos_off;
	for (e = 0; e < di; b += p->iluord[e], e++) {
		printf("pos %d = ",e);
		for (i = 0; i < p->iluord[e]; i++)
			printf("%f ",b[i]);
		printf("\n");
	}

	/* Input shaper curve */
	b = p->v + p->shp_off;
	for (e = 0; e < di; b += p->iluord[e], e++) {
		printf("shp %d = ",e);
		for (i = 0; i < p->iluord[e]; i++)
			printf("%f ",b[i]);
		printf("\n");
	}

	/* Matrix */
	b = p->v + p->mat_off;
	for (e = 0; e < (1 << di); e++) {
		printf("mx %d = ",e);
		for (f = 0; f < fdi; f++)
			printf("%f ",*b++);
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

/* - - - - - - - - - */

/* Setup the pseudo inverse information for each test point, */
/* using the current model. */
static void setup_piv(xfit *p) {
	int di =  p->di;
	int fdi = p->fdi;
	int i, e, f;

	/* Estimate in -> out partial derivatives */
	for (i = 0; i < p->nodp; i++) {
		double pd[MXDO][MXDI];
		double pp[MXDI];
		double vv[MXDIDO];

		/* Estimate in -> out partial derivatives */
		for (e = 0; e < di; e++)
			pp[e] = p->ipoints[i].p[e];

		/* Apply input shaper channel curves */
		for (e = 0; e < di; e++)
			vv[e] = icxSTransFunc(p->v + p->shp_offs[e], p->iluord[e], pp[e],  
			                       p->in_min[e], p->in_max[e]);

		/* Apply matrix cube interpolation */
		icxCubeInterp(p->v + p->mat_off, fdi, di, vv, vv);

		/* Apply output channel curves */
		for (f = 0; f < fdi; f++)
			vv[f] = icxSTransFunc(p->v + p->out_offs[f], p->oluord[f], vv[f],  
			                       p->out_min[f], p->out_max[f]);

		
		for (e = 0; e < di; e++) {
			double tt[MXDIDO];
	
			pp[e] += 1e-4;

			/* Apply input shaper channel curves */
			for (e = 0; e < di; e++)
				tt[e] = icxSTransFunc(p->v + p->shp_offs[e], p->iluord[e], pp[e],  
				                       p->in_min[e], p->in_max[e]);
	
			/* Apply matrix cube interpolation */
			icxCubeInterp(p->v + p->mat_off, fdi, di, tt, tt);
	
			/* Apply output channel curves */
			for (f = 0; f < fdi; f++)
				tt[f] = icxSTransFunc(p->v + p->out_offs[f], p->oluord[f], tt[f],  
				                       p->out_min[f], p->out_max[f]);

	
			for (f = 0; f < p->fdi; f++)
				pd[f][e] = (tt[f] - vv[f])/1e-4;
		
			pp[e] -= 1e-4;
		}

		/* Compute a psudo inverse matrix to map rout delta E to */
		/* in delta E in proportion to the pd magnitude. */
		for (f = 0; f < fdi; f++) {
			double ss = 0.0;

			for (e = 0; e < di; e++)		/* Sum of pd's ^4 */
				ss += pd[f][e] * pd[f][e] * pd[f][e] * pd[f][e];
			ss = sqrt(ss);
			if (ss > 1e-8) {
				for (e = 0; e < di; e++)
					p->piv[i].ide[f][e] = pd[f][e]/ss;
			} else {		/* Hmm. */
				for (e = 0; e < di; e++)
					p->piv[i].ide[f][e] = 0.0;
			}
		}
	}
}

/* - - - - - - - - - */

/* Function to pass to rspl to re-set output values, */
/* to make them relative to the white and black points */
static void
conv_rspl_out(
	void *pp,			/* relativectx structure */
	double *out,		/* output value */
	double *in			/* input value */
) {
	xfit *p    = (xfit *)pp;
	double tt[3];

	/* Convert the clut values to output values */
	xfit_outcurves(p, tt, out);

	if (p->flags & XFIT_OUT_LAB) {

		icmLab2XYZ(&icmD50, tt, tt);
		icmMulBy3x3(out, p->cmat, tt);
		icmXYZ2Lab(&icmD50, out, out);

	} else {	/* We are all in XYZ */
		icmMulBy3x3(out, p->cmat, tt);
	}

	/* And then convert them back to clut values */
	xfit_invoutcurves(p, out, out);
}

//#ifdef SPECIAL_TEST
/* - - - - - - - - - */

/* Execute the linear XYZ device model */
static void domodel(double *out, double *in) {
	double tmp[3];
	int i, j;
	double col[3][3];	/* sRGB additive colorant values in XYZ :- [out][in]  */

	col[0][0] = 0.412424;	/* X from R */
	col[0][1] = 0.357579;	/* X from G */
	col[0][2] = 0.180464;	/* X from B */
	col[1][0] = 0.212656;	/* Y from R */
	col[1][1] = 0.715158;	/* Y from G */
	col[1][2] = 0.0721856;	/* Y from B */
	col[2][0] = 0.0193324;	/* Z from R */
	col[2][1] = 0.119193;	/* Z from G */
	col[2][2] = 0.950444;	/* Z from B */

#ifdef SPECIAL_TEST_GAMMA	
	tmp[0] = pow(in[0], 1.9);
	tmp[1] = pow(in[1], 2.0);
	tmp[2] = pow(in[2], 2.1);
#else
	tmp[0] = in[0];
	tmp[1] = in[1];
	tmp[2] = in[2];
#endif

	for (j = 0; j < 3; j++) {
		out[j] = 0.0;
		for (i = 0; i < 3; i++)
			out[j] += col[j][i] * tmp[i];
	}
}

/* Function to pass to rspl to set nodes against */
/* synthetic model. */
static void
set_rspl_out1(
	void *pp,			/* relativectx structure */
	double *out,		/* output value */
	double *in			/* input value */
) {
	xfit *p  = (xfit *)pp;
	double tt[3], tout[3];

	/* Convert the input' values to input values */
	xfit_invinpscurves(p, tt, in);

	/* Synthetic linear rgb->XYZ model */
	domodel(tout, tt);

	/* Apply abs->rel white point adjustment */
	icmMulBy3x3(tout, p->cmat, tout);

#ifdef SPECIAL_TEST_LAB
	icmXYZ2Lab(&icmD50, tout, tout);
#endif
	/* And then convert them back to clut values */
	xfit_invoutcurves(p, tout, tout);

#ifdef DEBUG
printf("~1 changing %f %f %f -> %f %f %f\n", out[0], out[1], out[2], tout[0], tout[1], tout[2]);
#endif

	out[0] = tout[0];
	out[1] = tout[1];
	out[2] = tout[2];
}

//#endif /* SPECIAL_TEST */

/* - - - - - - - - - */
/* Do the fitting. */
/* return nz on error */
/* 1 = malloc or other error */
int xfit_fit(
	struct _xfit *p, 
	int flags,				/* Flag values */
	int di,					/* Input dimensions */
	int fdi,				/* Output dimensions */
	int rsplflags,			/* clut rspl creation flags */
	double *wp,				/* if flags & XFIT_OUT_WP_REL, */
							/* Initial white point, returns final wp */
	double *dw,				/* Device white value to adjust to be D50 */
	cow *ipoints,			/* Array of data points to fit - referece taken */
	int nodp,				/* Number of data points */
	double in_min[MXDI],	/* Input value scaling/domain minimum */
	double in_max[MXDI],	/* Input value scaling/domain maximum */
	int gres[MXDI],			/* clut resolutions being optimised for/returned */
	double out_min[MXDO],	/* Output value scaling/range minimum */
	double out_max[MXDO],	/* Output value scaling/range maximum */
	double smooth,			/* clut rspl smoothing factor */
	double oavgdev[MXDO],	/* Average output value deviation */
	int iord[],				/* Order of input pos/shaper curve for each dimension */
	int sord[],				/* Order of input sub-grid shaper curve (not used) */
	int oord[],				/* Order of output shaper curve for each dimension */
	optcomb tcomb,			/* Flag - target elements to fit. */
	void *cntx2,			/* Context of callbacks */
							/* Callback to convert two fit values delta E squared */
	double (*to_de2)(void *cntx, double *in1, double *in2),
							/* Same as above, with partial derivatives */
	double (*to_dde2)(void *cntx, double dout[2][MXDIDO], double *in1, double *in2)
) {
	int i, e, f;
	double *b;				/* Base of parameters for this section */
	int poff;

	if (tcomb & oc_io)		/* If we're doing anything, we need the matrix */
		tcomb |= oc_m;

	p->flags  = flags;
	if (flags & XFIT_VERB)
		p->verb = 1;
	else
		p->verb = 0;
	p->di      = di;
	p->fdi     = fdi;
	p->wp      = wp;		/* Take reference, so modified wp can be returned */
	p->dw      = dw;
	p->nodp    = nodp;
	p->ipoints = ipoints;
	p->tcomb   = tcomb;
	p->cntx2   = cntx2;
	for (e = 0; e < di; e++)
		p->gres[e] = gres[e];
	p->to_de2  = to_de2;
	p->to_dde2 = to_dde2;

#ifdef DEBUG
	printf("xfitc called with flags = 0x%x, di = %d, fdi = %d, nodp = %d, tcomb = 0x%x\n",flags,di,fdi,nodp,tcomb);
#endif

//printf("~1 out min = %f %f %f max = %f %f %f\n", out_min[0], out_min[1], out_min[2], out_max[0], out_max[1], out_max[2]);

	/* Sanity protect shaper orders and save scaling factors. */
	for (e = 0; e < di; e++) {
		if (iord[e] > MXLUORD)
			p->iluord[e] = MXLUORD;
		else
			p->iluord[e] = iord[e];
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
	p->shp_off = 0;
	for (poff = p->shp_off, p->shp_cnt = 0, e = 0; e < di; e++) {
		p->shp_offs[e] = poff;
		p->shp_cnt += p->iluord[e];
		poff      += p->iluord[e];
	}

	p->mat_off = p->shp_off + p->shp_cnt;
	for (poff = p->mat_off, p->mat_cnt = 0, f = 0; f < fdi; f++) {
		p->mat_offs[f] = poff;
		p->mat_cnt += (1 << di);
		poff       += (1 << di);
	}

	p->out_off = p->mat_off + p->mat_cnt;
	for (poff = p->out_off, p->out_cnt = 0, f = 0; f < fdi; f++) {
		p->out_offs[f] = poff;
		p->out_cnt += p->oluord[f];
		poff       += p->oluord[f];
	}

	p->pos_off = p->out_off + p->out_cnt;
	for (poff = p->pos_off, p->pos_cnt = 0, e = 0; e < di; e++) {
		p->pos_offs[e] = poff;
		p->pos_cnt += p->iluord[e];
		poff       += p->iluord[e];
	}

	p->tot_cnt = p->shp_cnt + p->mat_cnt + p->out_cnt + p->pos_cnt;
	if (p->tot_cnt > MXPARMS)
		error("xfit_fit: assert tot_cnt exceeds MXPARMS");

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

	/* Setup initial white point abs->rel conversions */
	if ((p->flags & XFIT_OUT_WP_REL) != 0) { 
		icmXYZNumber _wp;

		icmAry2XYZ(_wp, p->wp);

		/* Absolute->Aprox. Relative Adaptation matrix */
		icmChromAdaptMatrix(ICM_CAM_BRADFORD, icmD50, _wp, p->fromAbs);
	
		/* Aproximate relative to absolute conversion matrix */
		icmChromAdaptMatrix(ICM_CAM_BRADFORD, _wp, icmD50, p->toAbs);

		if (p->verb) {
			double lab[3];
			icmXYZ2Lab(&icmD50, lab, p->wp);
			printf("Initial White Point XYZ %f %f %f, Lab %f %f %f\n",
			    p->wp[0], p->wp[1], p->wp[2], lab[0], lab[1], lab[2]);
		}
	} else {
		icmSetUnity3x3(p->fromAbs);
		icmSetUnity3x3(p->toAbs);
	}

	/* Setup input position/shape curves to be linear initially */
	b = p->v + p->shp_off;
	for (e = 0; e < di; b += p->iluord[e], e++) {
		for (i = 0; i < p->iluord[e]; i++) {
			b[i] = 0.0;
		}
	}

	/* Setup matrix to be pure colorant' values initially */
	b = p->v + p->mat_off;
	for (e = 0; e < (1 << di); e++) {	/* For each colorant combination */
		int j, k, bk = 0;
		double bdif = 1e6;
		double ov[MXDO];
	
		/* Search the patch list to find the one closest to this input combination */
		for (k = 0; k < p->nodp; k++) {
			double dif = 0.0;
			for (j = 0; j < di; j++) {
				double tt;
				if (e & (1 << j))
					tt = p->in_max[j] - p->ipoints[k].p[j];
				else
					tt = p->in_min[j] - p->ipoints[k].p[j];
				dif += tt * tt;
			}
			if (dif < bdif) {		/* best so far */
				bdif = dif;
				bk = k;
				if (dif < 0.001)
					break;			/* Don't bother looking further */
			}
		}

		xfit_abs_to_rel(p, ov, p->ipoints[bk].v);

		for (f = 0; f < fdi; f++)
			b[f * (1 << di) + e] = ov[f];
	}


	/* Setup output curves to be linear initially */
	b = p->v + p->out_off;
	for (f = 0; f < fdi; b += p->oluord[f], f++) {
		for (i = 0; i < p->oluord[f]; i++) {
			b[i] = 0.0;
		}
	}

	/* Setup positioning curves to be linear initially */
	b = p->v + p->pos_off;
	for (e = 0; e < di; b += p->iluord[e], e++) {
		for (i = 0; i < p->iluord[e]; i++) {
			b[i] = 0.0;
		}
	}

	/* Create copy of input points with output converted to white relative */
	if (p->rpoints == NULL) {
		if ((p->rpoints = (cow *)malloc(p->nodp * sizeof(cow))) == NULL)
			return 1;
	}
	for (i = 0; i < p->nodp; i++) {
		p->rpoints[i].w = p->ipoints[i].w;
		for (e = 0; e < di; e++)
			p->rpoints[i].p[e] = p->ipoints[i].p[e];
		for (f = 0; f < fdi; f++)
			p->rpoints[i].v[f] = p->ipoints[i].v[f];

		/* out -> rout */
		xfit_abs_to_rel(p, p->rpoints[i].v, p->rpoints[i].v);
	}
  
	/* Allocate array of pseudo-inverse matricies */
	if ((p->flags & XFIT_FM_INPUT) != 0 && p->piv == NULL) {
		if ((p->piv = (xfit_piv *)malloc(p->nodp * sizeof(xfit_piv))) == NULL)
			return 1;
	}

	/* Allocate array of span DE's for current opt channel */
	{
		int lres = 0;
		for (e = 0; e < di; e++) {
			if (p->gres[e] > lres)
				lres = p->gres[e];
		}
		if ((p->uerrv = (double *)malloc(lres * sizeof(double))) == NULL)
			return 1;
	}

	/* Do the fitting one part at a time, then together */
	/* Shaper curves are created if poistion or shaper curves are requested */

	if ((p->tcomb & oc_ipo) != 0
	 && (p->tcomb & oc_m) == oc_m) {	/* Only bother with matrix if in and/or out */

		if (p->verb)
			printf("About to optimise temporary matrix\n");

		/* Setup pseudo-inverse if we need it */
		if (p->flags & XFIT_FM_INPUT)
			setup_piv(p);

#ifdef DEBUG
printf("\nBefore matrix opt:\n");
dump_xfit(p);
#endif
		/* Optimise matrix on its own */
		p->opt_ch = -1;
		p->opt_msk = oc_m;
		setup_xfit(p, p->wv, p->sa, 0.0, 0.5); 

#ifdef NODDV
		if (powell(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfitfunc, (void *)p) != 0)
			warning("set_icxLuLut: Powell failed to converge");
#else
		if (conjgrad(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfitfunc, dxfitfunc, (void *)p) != 0)
			warning("set_icxLuLut: Conjgrad failed to converge");
#endif
		for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
			p->v[p->opt_off + i] = p->wv[i];

#ifdef DEBUG
printf("\nAfter matrix opt:\n");
dump_xfit(p);
#endif
	}

	/* Optimise input and matrix together */
	if ((p->tcomb & oc_im) == oc_im) {

		if (p->verb)
			printf("\nAbout to optimise input curves and matrix\n");

		/* Setup pseudo-inverse if we need it */
		if (p->flags & XFIT_FM_INPUT)
			setup_piv(p);

		p->opt_ch = -1;
		p->opt_msk = oc_im;
		setup_xfit(p, p->wv, p->sa, 0.5, 0.3); 
#ifdef NODDV
		if (powell(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfitfunc, (void *)p) != 0)
			warning("set_icxLuLut: Powell failed to converge");
#else
		if (conjgrad(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfitfunc, dxfitfunc, (void *)p) != 0)
			warning("set_icxLuLut: Conjgrad failed to converge");
#endif
		for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
			p->v[p->opt_off + i] = p->wv[i];
#ifdef DEBUG
printf("\nAfter input and matrix opt:\n");
dump_xfit(p);
#endif
	}

	/* Optimise the matrix and output curves together */
	if ((p->tcomb & oc_mo) == oc_mo) {

		if (p->verb)
			printf("\nAbout to optimise output curves and matrix\n");

		/* Setup pseudo-inverse if we need it */
		if (p->flags & XFIT_FM_INPUT)
			setup_piv(p);

		p->opt_ch = -1;
		p->opt_msk = oc_mo;
		setup_xfit(p, p->wv, p->sa, 0.3, 0.3); 
#ifdef NODDV
		if (powell(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfitfunc, (void *)p) != 0)
			warning("set_icxLuLut: Powell failed to converge");
#else
		if (conjgrad(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfitfunc, dxfitfunc, (void *)p) != 0)
			warning("set_icxLuLut: Conjgrad failed to converge");
#endif
		for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
			p->v[p->opt_off + i] = p->wv[i];
#ifdef DEBUG
printf("\nAfter output opt:\n");
dump_xfit(p);
#endif

		/* Optimise input and matrix together again, after altering matrix */
		if ((p->tcomb & oc_im) == oc_im) {

			if (p->verb)
				printf("\nAbout to optimise input curves and matrix again\n");

			/* Setup pseudo-inverse if we need it */
			if (p->flags & XFIT_FM_INPUT)
				setup_piv(p);

			p->opt_ch = -1;
			p->opt_msk = oc_im;
			setup_xfit(p, p->wv, p->sa, 0.2, 0.2); 
#ifdef NODDV
			if (powell(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfitfunc, (void *)p) != 0)
				warning("set_icxLuLut: Powell failed to converge");
#else
			if (conjgrad(NULL, p->opt_cnt, p->v, p->sa, POWTOL, MAXITS, xfitfunc, dxfitfunc, (void *)p) != 0)
				warning("set_icxLuLut: Conjgrad failed to converge");
#endif
			for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
				p->v[p->opt_off + i] = p->wv[i];
		}
#ifdef DEBUG
printf("\nAfter 2nd input and matrix opt:\n");
dump_xfit(p);
#endif

#ifndef NODDV
		/* Optimise all together */
		if ((p->tcomb & oc_imo) == oc_imo) {

			if (p->verb)
				printf("\nAbout to optimise input, matrix and output together\n");

			/* Setup pseudo-inverse if we need it */
			if (p->flags & XFIT_FM_INPUT)
				setup_piv(p);

			p->opt_ch = -1;
			p->opt_msk = oc_imo;
			setup_xfit(p, p->wv, p->sa, 0.1, 0.1); 
			if (conjgrad(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfitfunc, dxfitfunc, (void *)p) != 0)
				warning("set_icxLuLut: Conjgrad failed to converge");
			for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
				p->v[p->opt_off + i] = p->wv[i];

			/* Setup final pseudo-inverse from shaper/matrix/out */
			if (p->flags & XFIT_FM_INPUT)
				setup_piv(p);
		}

#ifdef DEBUG
printf("\nAfter all together opt:\n");
dump_xfit(p);
#endif

#endif /* !NODDV */

		/* Adjust output curve white point. */
		/* This is for the benefit of the B2A table */
		if (p->flags & XFIT_OUT_ZERO) {

			if (p->verb)
				printf("\nAbout to adjust a and b output curves for white point\n");

			for (f = 1; f < 3 && f < p->fdi; f++) {
				p->opt_ch = f;
				p->wv[0] = p->v[p->out_offs[f]];	/* Current parameter value */
				p->sa[0] = 0.1;					/* Search radius */
				if (powell(NULL, 1, p->wv, p->sa, 0.0000001, 1000, symoptfunc, (void *)p) != 0)
					error("set_icxLuLut: Powell failed to converge");
				p->v[p->out_offs[f]] = p->wv[0];	/* Copy results back */
			}
		}
	}

	/* In case we don't generate position curves, */
	/* copy the input curves to the position, so that */
	/* iwidth is computed correctly */
	{
		double *bb;
		
		b = p->v + p->shp_off;
		bb = p->v + p->pos_off;
		for (e = 0; e < di; b += p->iluord[e], bb += p->iluord[e], e++) {
			for (i = 0; i < p->iluord[e]; i++) {
				bb[i] = b[i];
			}
		}
	}


	/* If we want position curves, generate them */
	if (p->tcomb & oc_p) {
		int ee;

		if (p->verb)
			printf("\nAbout to create grid position input curves\n");

		/* Allocate in->rout duplicate point set */
		if (p->rpoints == NULL) {
			if ((p->rpoints = (cow *)malloc(p->nodp * sizeof(cow))) == NULL)
				return 1;
		}

		/* Create a set of 1D rspl setup points that contains */
		/* the residual error */
		for (i = 0; i < p->nodp; i++) {
			double tv[MXDO];		/* Target output value */
			double mv[MXDO];		/* Model output value */
			double ev;

			xfit_abs_to_rel(p, tv, p->ipoints[i].v);
			xfit_shmatsh(p, mv, p->ipoints[i].p);

			/* Evaluate the error squared */
			if (p->flags & XFIT_FM_INPUT) {
				double pp[MXDI];
				for (e = 0; e < di; e++)
					pp[e] = p->ipoints[i].p[e];
				for (f = 0; f < fdi; f++) {
					double t1 = tv[f] - mv[f];	/* Error in output */
		
					/* Create input point offset by equivalent delta to output point */
					/* error, in proportion to the partial derivatives for that output. */
					for (e = 0; e < di; e++)
						pp[e] += t1 * p->piv[i].ide[f][e];
				}
				ev = p->to_de2(p->cntx2, pp, p->ipoints[i].p);
			} else {
				ev = p->to_de2(p->cntx2, mv, tv);
			}
			p->rpoints[i].v[0] = ev;
			p->rpoints[i].w    = p->ipoints[i].w;
		}

		/* Do each input axis in turn */
		for (ee = 0; ee < p->di; ee++) {
			rspl *resid;
			double imin[1],imax[1],omin[1],omax[1]; 
			int resres[1] = { 1024 };
#define NPGP 100
			mcv *posc;
			mcvco pgp[NPGP];
			double vo, vs;
			double *pms;
			
			/* Create a rspl that plots the residual error */
			/* vs the axis value */
			for (i = 0; i < p->nodp; i++)
				p->rpoints[i].p[0] = p->ipoints[i].p[ee];

			imin[0] = in_min[ee];
			imax[0] = in_max[ee];
			omin[0] = 0.0;
			omax[0] = 0.0;
			
			if ((resid = new_rspl(1, 1)) == NULL)
				return 1;

			resid->fit_rspl_w(resid, 0, p->rpoints, p->nodp, imin, imax, resres,
				omin, omax, 2.0, NULL, NULL);

#ifdef DEBUG_PLOT
			{
#define	XRES 100
				double xx[XRES];
				double y1[XRES];

				printf("Input residual error channel %d\n",ee);
				for (i = 0; i < XRES; i++) {
					co pp;
					double x;
					x = i/(double)(XRES-1);
					xx[i] = x = x * (imax[0] - imin[0]) + imin[0];
					pp.p[0] = xx[i];
					resid->interp(resid, &pp);
					y1[i] = pp.v[0];
					if (y1[i] < 0.0)
						y1[i] = 0.0;
					y1[i] = pow(y1[i], 0.5);	/* Convert from error^2 to error */
				}
				do_plot(xx,y1,NULL,NULL,XRES);
			}
#endif /* DEBUG_PLOT */

			/* Create a set of guide points that contain */
			/* the accumulated residual error vs. the axis value */
			for (i = 0; i < NPGP; i++) {
				co pp;
				double vv;

				pp.p[0] = i/(NPGP-1.0);
				resid->interp(resid, &pp);

				pgp[i].p = (pp.p[0] - in_min[ee])/(in_max[ee] - in_min[ee]);
				pgp[i].w = 1.0; 

				vv = pp.v[0];
				if (vv < 0.0)
					vv = 0.0;
				vv = pow(vv, 0.5);		/* Convert from error^2 to error */
				vv += PSHAPE_MINE;		/* In case error is near zero */
				vv = pow(vv, 1.00);		/* Agressivness of grid distribution */

				if (i == 0)
					pgp[i].v = vv;
				else
					pgp[i].v = pgp[i-1].v + vv;
			}
			resid->del(resid);

			/* Normalize the output range */
			vo = pgp[0].v;
			vs = pgp[NPGP-1].v - vo;
			
			for (i = 0; i < NPGP; i++) {
				pgp[i].v = (pgp[i].v - vo)/vs;
				
//printf("~1 guide point %d: %f -> %f\n",i,pgp[i].p,pgp[i].v);
			}
			/* Fit the non-monotonic parameters to the guide points */
			if ((posc = new_mcv_noos()) == NULL)
				return 1;

			posc->fit(posc, 0, p->iluord[ee], pgp, NPGP, 1.0);

#ifdef DEBUG_PLOT
			{
#define	XRES 100
				double xx[XRES];
				double y1[XRES];

				printf("Position curve %d\n",ee);
				for (i = 0; i < XRES; i++) {
					xx[i] = i/(double)(XRES-1);
					y1[i] = posc->interp(posc, xx[i]);
				}
				do_plot(xx,y1,NULL,NULL,XRES);
			}
#endif /* DEBUG_PLOT */

			/* Transfer parameters to xfit pos (skip offset and scale) */
			posc->get_params(posc, &pms);

			for (i = 0; i < p->iluord[ee]; i++) {
				p->v[p->pos_offs[ee] + i] = pms[i+2];
			}
//p->v[p->in_offs[ee]] = -1.5;

			free(pms);
			posc->del(posc);
		}
	}

	if (p->verb)
		printf("\n");

#ifdef DEBUG
printf("\nFinal parameters:\n");
dump_xfit(p);
#endif

#ifdef DEBUG_PLOT
	{
#define	XRES 100
		double xx[XRES];
		double y1[XRES];

		for (e = 0; e < p->di; e++) {
			printf("Input position curve channel %d\n",e);
			for (i = 0; i < XRES; i++) {
				double x;
				x = i/(double)(XRES-1);
				xx[i] = x = x * (p->in_max[e] - p->in_min[e]) + p->in_min[e];
				y1[i] = xfit_poscurve(p, x, e);
			}
			do_plot(xx,y1,NULL,NULL,XRES);
		}

		for (e = 0; e < p->di; e++) {
			printf("Input shape curve channel %d\n",e);
			for (i = 0; i < XRES; i++) {
				double x;
				x = i/(double)(XRES-1);
				xx[i] = x = x * (p->in_max[e] - p->in_min[e]) + p->in_min[e];
				y1[i] = xfit_shpcurve(p, x, e);
			}
			do_plot(xx,y1,NULL,NULL,XRES);
		}
	
		for (e = 0; e < p->di; e++) {
			printf("Combined input curve channel %d\n",e);
			for (i = 0; i < XRES; i++) {
				double x;
				x = i/(double)(XRES-1);
				xx[i] = x = x * (p->in_max[e] - p->in_min[e]) + p->in_min[e];
				y1[i] = p->incurve(p, x, e);
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

	/* Create final clut rspl using the established pos/shape/output curves */
	/* and white point */ 
	if (flags & XFIT_MAKE_CLUT) {
		double *iwidth[MXDI];

		/* Create an in' -> rout' scattered test point set */
		if (p->rpoints == NULL) {
			if ((p->rpoints = (cow *)malloc(p->nodp * sizeof(cow))) == NULL)
				return 1;
		}
		for (i = 0; i < p->nodp; i++) {
			p->rpoints[i].w = p->ipoints[i].w;
			xfit_inpscurves(p, p->rpoints[i].p, p->ipoints[i].p);
			for (f = 0; f < fdi; f++)
				p->rpoints[i].v[f] = p->ipoints[i].v[f];
			xfit_abs_to_rel(p, p->rpoints[i].v, p->rpoints[i].v);
			xfit_invoutcurves(p, p->rpoints[i].v, p->rpoints[i].v);
// ~~999
//printf("~1 point %d, w %f, %f %f %f %f -> %f %f %f\n",
//i,p->rpoints[i].w,p->rpoints[i].p[0], p->rpoints[i].p[1], p->rpoints[i].p[2], p->rpoints[i].p[3],
//p->rpoints[i].v[0], p->rpoints[i].v[1], p->rpoints[i].v[2]);
		}

		/* Create iwidth[] arrays, that hold the shaper space */
		/* grid spacing due to the poistioning curves. */
		/* This tells the rspl scattered data interpolator */
		/* the grid spacing that smoothness should be */
		/* measured against. */
#ifdef DEBUG
printf("~1 about to setup iwidth\n");
#endif
		for (e = 0; e < p->di; e++) {
			double lv, cv;
//printf("~1 e = %d\n",e);
			if ((iwidth[e] = (double *)malloc((p->gres[e]-1) * sizeof(double))) == NULL)
				return 1;
//printf("~1 about to do %d spans\n",p->gres[e]);
			for (i = 0; i < p->gres[e]; i++) {
				cv = (double)i/p->gres[e];

//printf("~1 i = %d, pos space = %f\n",i,cv);
				/* Inverse lookup grid position through positioning curve */
				/* to give device space type value */
				cv = icxInvTransFunc(p->v + p->pos_offs[e], p->iluord[e], cv); 

//printf("~1 dev space = %f\n",cv);
				/* Forward lookup device type value through the shaper curve */
				/* to give value in shape linearized space. */
				cv = icxTransFunc(p->v + p->shp_offs[e], p->iluord[e], cv); 
//printf("~1 shape space = %f\n",cv);

				if (i > 0) {
					double wi = cv - lv;
					if (wi < 1e-5)
						wi = 1e-5;		/* Hmm. */
//wi = 1.0;
					iwidth[e][i-1] = wi;
#ifdef DEBUG
printf("~1 iwidth[%d][%d] = %f\n",e,i-1,wi);
#endif
				}
				lv = cv;
			}
		}

		if (p->clut != NULL)
			p->clut->del(p->clut);
		if ((p->clut = new_rspl(di, fdi)) == NULL)
			return 1;

		if (p->verb)
			printf("Create final clut from scattered data\n");

// ~~~999
//out_max[0] = 100.390625;

		p->clut->fit_rspl_w(p->clut, rsplflags, p->rpoints, p->nodp, in_min, in_max, gres,
			out_min, out_max, smooth, oavgdev, iwidth);
		if (p->verb) printf("\n");

		for (e = 0; e < p->di; e++)
			free(iwidth[e]);

		/* Now do a final white point fixup */
		if (p->flags & XFIT_OUT_WP_REL) {
			co wcc;
			icmXYZNumber _wp;

			if (p->verb)
				printf("Doing White point fine tune:\n");
				
			/* See what current device white maps to */
			xfit_inpscurves(p, wcc.p, dw);
			p->clut->interp(p->clut, &wcc);
			xfit_outcurves(p, wcc.v, wcc.v);

//printf("~1 device white %f %f %f %f -> %f %f %f\n",
//dw[0], dw[1], dw[2], dw[3], wcc.v[0], wcc.v[1], wcc.v[2]);

			if (p->flags & XFIT_OUT_LAB)
				icmLab2XYZ(&icmD50, wcc.v, wcc.v);

//printf("~1 XYZ white %f %f %f\n", wcc.v[0], wcc.v[1], wcc.v[2]);

			/* Fixup conversion matricies so that they transform */
			/* the absolute input values to perfect relative. */
			icmAry2XYZ(_wp, wcc.v);

			/* Absolute->Relative Adaptation matrix */
			icmChromAdaptMatrix(ICM_CAM_BRADFORD | ICM_CAM_MULMATRIX, icmD50, _wp, p->fromAbs);
		
			/* Relative to absolute conversion matrix */
			icmChromAdaptMatrix(ICM_CAM_BRADFORD | ICM_CAM_MULMATRIX, _wp, icmD50, p->toAbs);
	
			/* return accurate white point to caller */
			icmXYZ2Ary(p->wp, icmD50);
			icmMulBy3x3(p->wp, p->toAbs, p->wp);

			if (p->verb) {
				double labwp[3];
				icmXYZ2Lab(&icmD50, labwp, wcc.v);
				printf("Before fine tune, WP = XYZ %f %f %f, Lab %f %f %f\n",
				wcc.p[0], wcc.p[1],wcc.p[2], labwp[0], labwp[1], labwp[2]);
			}

			/* Matrix needed to correct from approximate to exact D50 */
			icmChromAdaptMatrix(ICM_CAM_BRADFORD, icmD50, _wp, p->cmat);
	
			p->clut->re_set_rspl(
				p->clut,		/* this */
				0,				/* Combination of flags */
				(void *)p,		/* Opaque function context */
				conv_rspl_out	/* Function to set from */
			);
		}

#ifdef SPECIAL_FORCE
		/* Replace the rspl nodes with ones directly computed */
		/* from the synthetic linear RGB->XYZ model */
		{
			double twp[3];
			icmXYZNumber _wp;

			/* See what current device white maps to */
			twp[0] = twp[1] = twp[2] = 1.0;
			domodel(twp, twp);

			icmAry2XYZ(_wp, twp);

			/* Matrix needed to correct to D50 */
			icmChromAdaptMatrix(ICM_CAM_BRADFORD, icmD50, _wp, p->cmat);
	
			p->clut->re_set_rspl(
				p->clut,		/* this */
				0,				/* Combination of flags */
				(void *)p,		/* Opaque function context */
				set_rspl_out1	/* Function to set from */
			);
		}
#endif /* SPECIAL_FORCE */

	/* Evaluate the residual error now, with the rspl in place */
#ifdef DEBUG_PLOT
	{
		int ee;
		double maxe = 0.0;
		double avee = 0.0;

		/* Allocate in->rout duplicate point set */
		if (p->rpoints == NULL) {
			if ((p->rpoints = (cow *)malloc(p->nodp * sizeof(cow))) == NULL)
				return 1;
		}

		/* Create a set of 1D rspl setup points that contains */
		/* the residual error */
		for (i = 0; i < p->nodp; i++) {
			double tv[MXDO];		/* Target output value */
			double mv[MXDO];		/* Model output value */
			co pp;
			double ev;

			xfit_abs_to_rel(p, tv, p->ipoints[i].v);

			for (e = 0; e < p->di; e++)
				pp.p[e] = p->ipoints[i].p[e];

			xfit_inpscurves(p, pp.p, pp.p);
			p->clut->interp(p->clut, &pp);
			xfit_outcurves(p, pp.v, pp.v);

			for (f = 0; f < p->fdi; f++)
				mv[f] = pp.v[f];

			/* Evaluate the residual error suqared */
			if (p->flags & XFIT_FM_INPUT) {
				double pp[MXDI];
				for (e = 0; e < di; e++)
					pp[e] = p->ipoints[i].p[e];
				for (f = 0; f < fdi; f++) {
					double t1 = tv[f] - mv[f];	/* Error in output */
		
					/* Create input point offset by equivalent delta to output point */
					/* error, in proportion to the partial derivatives for that output. */
					for (e = 0; e < di; e++)
						pp[e] += t1 * p->piv[i].ide[f][e];
				}
				ev = p->to_de2(p->cntx2, pp, p->ipoints[i].p);
			} else {
				ev = p->to_de2(p->cntx2, mv, tv);
			}
//printf("~1 point %d, loc %f %f %f %f\n",i, p->rpoints[i].p[0], p->rpoints[i].p[1], p->rpoints[i].p[2], p->rpoints[i].p[3]);
//printf("        targ %f %f %f, is %f %f %f\n", tv[0], tv[1], tv[2], mv[0], mv[1], mv[2]);

			p->rpoints[i].v[0] = ev;
			p->rpoints[i].w    = p->ipoints[i].w;

			ev = sqrt(ev);
			if (ev > maxe)
				maxe = ev;
			avee += ev;
		}
		printf("Max err = %f, avg err = %f\n",maxe, avee/(double)p->nodp);

		/* Evaluate each input axis in turn */
		for (ee = 0; ee < p->di; ee++) {
			rspl *resid;
			double imin[1],imax[1],omin[1],omax[1]; 
			int resres[1] = { 1024 };
			
			/* Create a rspl that gives the residual error squared */
			/* vs. the axis value */
			for (i = 0; i < p->nodp; i++)
				p->rpoints[i].p[0] = p->ipoints[i].p[ee];

			imin[0] = in_min[ee];
			imax[0] = in_max[ee];
			omin[0] = 0.0;
			omax[0] = 0.0;
			
			if ((resid = new_rspl(1, 1)) == NULL)
				return 1;

			resid->fit_rspl_w(resid, 0, p->rpoints, p->nodp, imin, imax, resres,
				omin, omax, 2.0, NULL, NULL);

			{
#define	XRES 100
				double xx[XRES];
				double y1[XRES];

				printf("Finale residual error vs. input channel %d\n",ee);
				for (i = 0; i < XRES; i++) {
					co pp;
					double x;
					x = i/(double)(XRES-1);
					xx[i] = x = x * (imax[0] - imin[0]) + imin[0];
					pp.p[0] = xx[i];
					resid->interp(resid, &pp);
					y1[i] = sqrt(fabs(pp.v[0]));
				}
				do_plot(xx,y1,NULL,NULL,XRES);
			}
			resid->del(resid);
		}
	}
#endif /* DEBUG_PLOT */

	/* Special test code to figure out what's wrong with position curves */
#ifdef NEVER
	{
		double rgb[3], xyz[3];
		co pp;

extern int rspldb;
db = 1;
rspldb = 1;
		printf("~1 gres = %d %d %d\n", gres[0], gres[1], gres[2]);
		for (i = 0; i < 3; i++) {

			if (i == 0) {
				/* Test point on a grid point */
				printf("\n~1 ##########################################\n");
				printf("~1 testing input at first diagonal grid point\n");
	
				rgb[0] = 1.0/(gres[0]-1.0);
				rgb[1] = 1.0/(gres[0]-1.0);
				rgb[2] = 1.0/(gres[0]-1.0);

				printf("~1 target rgb' = %f %f %f\n", rgb[0], rgb[1], rgb[2]);
				xfit_invinpscurves(p, rgb, rgb);

			} else if (i == 1) {
				/* Test point half way through a grid point */
				printf("\n~1 ##########################################\n");
				printf("\n~1 testing half way through diagonal grid point\n");
	
				rgb[0] = 0.5/(gres[0]-1.0);
				rgb[1] = 0.5/(gres[0]-1.0);
				rgb[2] = 0.5/(gres[0]-1.0);

				printf("~1 target rgb' = %f %f %f\n", rgb[0], rgb[1], rgb[2]);
				xfit_invinpscurves(p, rgb, rgb);

			} else {
				printf("\n~1 ##########################################\n");
				printf("\n~1 testing worst case point\n");
	
				rgb[0] = 0.039915;
				rgb[1] = 0.053148;
				rgb[2] = 0.230610;
			}

			pp.p[0] = rgb[0];
			pp.p[1] = rgb[1];
			pp.p[2] = rgb[2];

			printf("~1 rgb = %f %f %f\n", pp.p[0], pp.p[1], pp.p[2]);

			xfit_inpscurves(p, pp.p, pp.p);
			
			printf("~1 rgb' = %f %f %f\n", pp.p[0], pp.p[1], pp.p[2]);

			p->clut->interp(p->clut, &pp);

			printf("~1 xyz' = %f %f %f\n", pp.v[0], pp.v[1], pp.v[2]);

			xfit_outcurves(p, pp.v, pp.v);

			printf("~1 xyz = %f %f %f\n", pp.v[0], pp.v[1], pp.v[2]);

			/* Synthetic linear rgb->XYZ model */
			domodel(xyz, rgb);

			/* Apply abs->rel white point adjustment */
			icmMulBy3x3(xyz, p->mat, xyz);

			printf("~1 ref = %f %f %f, de = %f\n", xyz[0], xyz[1], xyz[2],icmXYZLabDE(&icmD50,xyz,pp.v));
		}

db = 0;
rspldb = 0;
//		exit(0);
	}
#endif /* NEVER */

	}
	return 0;
}

/* We're done with an xfit */
static void xfit_del(xfit *p) {
	if (p->v != NULL)
		free(p->v);
	if (p->wv != NULL)
		free(p->wv);
	if (p->sa != NULL)
		free(p->sa);
	if (p->rpoints != NULL)
		free(p->rpoints);
	if (p->piv != NULL)
		free(p->piv);
	if (p->uerrv != NULL)
		free(p->uerrv);
	free(p);
}

/* Create a transform fitting object */
/* return NULL on error */
xfit *new_xfit(
) {
	xfit *p;

	if ((p = (xfit *)calloc(1, sizeof(xfit))) == NULL) {
		return NULL;
	}

	/* Set method pointers */
	p->fit         = xfit_fit;
	p->incurve     = xfit_inpscurve;
	p->invincurve  = xfit_invinpscurve;
	p->outcurve    = xfit_outcurve;
	p->invoutcurve = xfit_invoutcurve;
	p->del = xfit_del;

	return p;
}


