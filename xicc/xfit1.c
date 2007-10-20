
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
 * create per channel input and/or output curves for clut profiles.
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

#undef NODDV			/* Use slow non d/dv powell */

#define POWTOL 5e-4			/* Shaper Powell optimiser tollerance */
#define MAXITS 2000			/* Shaper number of itterations before giving up */

/* Weights in shaper parameters, to minimise unconstrained "wiggles" */
#define SHAPE_WEIGHT 5000.0		/* Overal shaper weight contribution - err on side of smoothness */
#define SHAPE_BASE  0.00001		/* 0 & 1 harmonic weight */
#define SHAPE_HBASE 0.0001		/* 2nd and higher additional weight */

#define MXLUORD 10		/* Shaper harmonic orders to use */

/*
 * TTBD:
 *
 */

/* - - - - - - - - - */

/* Lookup a value though an input curve */
static double xfit1_incurve(xfit1 *p, double in, int chan) {
	double rv;
	rv = icxSTransFunc(p->v + p->in_offs[chan], p->iluord[chan], in,  
			                       p->in_min[chan], p->in_max[chan]);
	return rv;
}

/* Inverse Lookup a value though an input curve */
static double xfit1_invincurve(xfit1 *p, double in, int chan) {
	double rv;
	rv = icxInvSTransFunc(p->v + p->in_offs[chan], p->iluord[chan], in,  
			                       p->in_min[chan], p->in_max[chan]);
	return rv;
}

/* Lookup a value though an output curve */
static double xfit1_outcurve(xfit1 *p, double in, int chan) {
	double rv;
	rv = icxSTransFunc(p->v + p->out_offs[chan], p->oluord[chan], in,  
		                       p->out_min[chan], p->out_max[chan]);
	return rv;
}

/* Inverse Lookup a value though an output curve */
static double xfit1_invoutcurve(xfit1 *p, double in, int chan) {
	double rv;
	rv = icxInvSTransFunc(p->v + p->out_offs[chan], p->oluord[chan], in,  
		                       p->out_min[chan], p->out_max[chan]);
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

/* Delta E squared replacement for device value output */
/* including partial derivatives. */
static double ddeverr(double dout_de[2][MXDO], double *in1, double *in2, int fdi) {
	int i;
	double tt, rv = 0.0;
	for (i = 0; i < fdi; i++) {
		tt = 100.0 * (in1[i] - in2[i]);
		dout_de[0][i] =  2.0 * 100.0 * tt;
		dout_de[1][i] = -2.0 * 100.0 * tt;
		rv += tt * tt;
	}
	rv += tt * tt;
	
	return rv;
}

/* Convert from device/unknown to delta E squared using provided function */
/* including partial derivatives of in1. */
static double dlookup(xfit1 *p, double dout_lu[2][MXDO], double *in1, double *in2, int fdi) {

	int e;
	double out, dout, rin[MXDI];

	for (e = 0; e < fdi; e++)
		rin[e] = in1[e];

	/* Reference conversion and output */
	out = p->to_de2(p->cntx2, in1, in2);

	/* Partial derivatives of in1 - the hard way */
	for (e = 0; e < fdi; e++) {
		rin[e] += 1e-5;
		dout = p->to_de2(p->cntx2, rin, in2);
		dout_lu[0][e] = (dout - out)/1e-5;
		dout_lu[1][e] = 0.0;		/* Don't compute partial de of in2 since it's fixed */
		rin[e] -= 1e-5;
	}

	return out;
}

/* - - - - - - - - - */

/* return a weighting for the magnitude of the in and out */
/* shaping parameters. This is to reduce unconstrained "wiggles" */
static double shapmag(
xfit1  *p			/* Base of optimisation structure */
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
		b = p->v + p->in_off;
		for (e = 0; e < di; e++) {
			for (k = 0; k < p->iluord[e]; k++) {
				tt = *b++;
				tt *= tt;
				if (k <= 1)	/* First or second curves */
					w = SHAPE_BASE;
				else
					w = SHAPE_BASE + (k-1) * SHAPE_HBASE;
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
				tt = *b++;
				tt *= tt;
				if (k <= 1)	/* First or second curves */
					w = SHAPE_BASE;
				else
					w = SHAPE_BASE + (k-1) * SHAPE_HBASE;
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
xfit1  *p,			/* Base of optimisation structure */
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
		b = p->v + p->in_off;
		c = dav + p->in_off;
		for (e = 0; e < di; e++) {
			for (k = 0; k < p->iluord[e]; k++) {
				if (k <= 1)	/* First or second curves */
					w = SHAPE_BASE;
				else
					w = SHAPE_BASE + (k-1) * SHAPE_HBASE;
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
				if (k <= 1)	/* First or second curves */
					w = SHAPE_BASE;
				else
					w = SHAPE_BASE + (k-1) * SHAPE_HBASE;
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
static double xfit1func(void *edata, double *v) {
	xfit1 *p = (xfit1 *)edata;
	double tw = 0.0;				/* Total weight */
	double rv = 0.0, smv;
	double tin[MXDI], out[MXDO];
	int di = p->di;
	int fdi = p->fdi;
	int i, e, f;

	/* Copy the parameters being optimised into xfit1 structure */
	for (i = 0; i < p->opt_cnt; i++)
		p->v[p->opt_off + i] = v[i];

	/* For all our data points */
	for (i = 0; i < p->nodp; i++) {
		double ev;

		/* Apply input channel curves */
		for (e = 0; e < di; e++)
			tin[e] = icxSTransFunc(p->v + p->in_offs[e], p->iluord[e], p->points[i].p[e],  
			                       p->in_min[e], p->in_max[e]);

		/* Apply matrix cube interpolation */
		icxCubeInterp(p->v + p->mat_off, fdi, di, out, tin);

		/* Apply output channel curves */
		for (f = 0; f < fdi; f++)
			out[f] = icxSTransFunc(p->v + p->out_offs[f], p->oluord[f], out[f],  
			                       p->out_min[f], p->out_max[f]);

		if (p->flags & XFIT_FM_INPUT)
			error("xfit1 can't handle XFIT_FM_INPUT flag!");

		if ((p->flags & XFIT_FM_MASK) == XFIT_FM_LU) {
			ev = p->to_de2(p->cntx2, out, p->points[i].v);

		} else if ((p->flags & XFIT_FM_MASK) == XFIT_FM_DEV) {
			ev = deverr(out, p->points[i].v, fdi);

		} else if ((p->flags & XFIT_FM_MASK) == XFIT_FM_XYZ) {
			icmLab2XYZ(&icmD50, out, out);

#ifdef USE_CIE94_DE
			ev = icmCIE94sq(out, p->points[i].v);
#else
			ev = icmLabDEsq(out, p->points[i].v);
#endif
		} else {
			
#ifdef USE_CIE94_DE
			ev = icmCIE94sq(out, p->points[i].v);
#else
			ev = icmLabDEsq(out, p->points[i].v);
#endif
		}

		tw += p->points[i].w;
		rv += p->points[i].w * ev;
	}

	/* Normalise error to be an average delta E squared */
	rv /= tw;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	smv = shapmag(p);
	rv += smv;

#ifdef DEBUG
printf("~1(%f)xfit1func returning %f\n",smv,rv);
#endif

#ifdef NODDV
	if (p->verb)
		printf("."), fflush(stdout);
#endif
	return rv;
}

/* Shaper+Matrix optimisation function with partial derivatives, */
/* handed to conjgrad() */
static double dxfit1func(void *edata, double *dv, double *v) {
	xfit1 *p = (xfit1 *)edata;
	double tw = 0.0;				/* Total weight */
	double rv = 0.0, smv;
	double tin[MXDI], out[MXDO];

	double dav[MXPARMS];				/* Overall del due to del param vals */

	double dtin_iv[MXDI * MXLUORD];		/* Del in itrans out due to del itrans param vals */
	double dmato_mv[1 << MXDI];			/* Del in mat out due to del in matrix param vals */
	double dmato_tin[MXDO * MXDI];		/* Del in mat out due to del in matrix input values */
	double dout_ov[MXDO * MXLUORD];		/* Del in otrans out due to del in otrans param values */
	double dout_mato[MXDO];				/* Del in otrans out due to del in otrans input values */

	double dout_de[2][MXDO];			/* Del in DE due to two output values */

	int di = p->di;
	int fdi = p->fdi;
	int i, jj, k, e, ee, f, ff;

	/* Copy the parameters being optimised into xfit1 structure */
	for (i = 0; i < p->opt_cnt; i++)
		p->v[p->opt_off + i] = v[i];

	/* Zero the accumulated partial derivatives */
	/* We compute deriv for all parameters (not just current optimised) */
	for (i = 0; i < p->tot_cnt; i++)
		dav[i] = 0.0;

	/* For all our data points */
	for (i = 0; i < p->nodp; i++) {
		double ev;

		/* Apply input channel curves */
		for (e = 0; e < di; e++)
			tin[e] = icxdpSTransFunc(p->v + p->in_offs[e], &dtin_iv[p->in_offs[e] - p->in_off],
		                         p->iluord[e], p->points[i].p[e], p->in_min[e], p->in_max[e]);

		/* Apply matrix cube interpolation */
		icxdpdiCubeInterp(p->v + p->mat_off, dmato_mv, dmato_tin, fdi, di, out, tin);

		/* Apply output channel curves */
		for (f = 0; f < fdi; f++)
			out[f] = icxdpdiSTransFunc(p->v + p->out_offs[f],
			                           &dout_ov[p->out_offs[f] - p->out_off], &dout_mato[f],
			                           p->oluord[f], out[f], p->out_min[f], p->out_max[f]);


		/* Convert to Delta E and compute pde's into dout_de */
		if (p->flags & XFIT_FM_INPUT)
			error("xfit1 can't handle XFIT_FM_INPUT flag!");

		if ((p->flags & XFIT_FM_MASK) == XFIT_FM_LU) {

			ev = dlookup(p, dout_de, out, p->points[i].v, fdi);

		} else if ((p->flags & XFIT_FM_MASK) == XFIT_FM_DEV) {
			ev = ddeverr(dout_de, out, p->points[i].v, fdi);

		} else if ((p->flags & XFIT_FM_MASK) == XFIT_FM_XYZ) {
			double dout_lab[3][3];
			double dout_tde[2][3];

			icxdXYZ2Lab(&icmD50, out, dout_lab, out);

#ifdef USE_CIE94_DE
			ev = icxdCIE94sq(dout_tde, out, p->points[i].v);
#else
			ev = icxdLabDEsq(dout_tde, out, p->points[i].v);
#endif
			/* Concatentate de's from conversion to delta E */
			for (ff = 0; ff < 1; ff++) {		/* New output */
				for (ee = 0; ee < 3; ee++) {	/* Input */
					dout_de[ff][ee] = 0.0;
					for (k = 0; k < 3; k++) {		/* Lab or device channels */
						dout_de[ff][ee] += dout_lab[ff][k] * dout_tde[k][ee];
					}
				}
			}
		} else {
			double dout_tde[2][3];
			
#ifdef USE_CIE94_DE
			ev = icxdCIE94sq(dout_tde, out, p->points[i].v);
#else
			ev = icxdLabDEsq(dout_tde, out, p->points[i].v);
#endif
			for (ff = 0; ff < 1; ff++) {
				for (ee = 0; ee < 3; ee++) {
					dout_de[ff][ee] = dout_tde[ff][ee];
				}
			}
		}

		/* Accumulate total weighted delta E squared */
		tw += p->points[i].w;
		rv += p->points[i].w * ev;

		for (ff = 0; ff < fdi; ff++) {				/* Parameter output chanel */
			for (k = 0; k < p->oluord[ff]; k++) {	/* Param within channel */
				double vv = 0.0;
				jj = p->out_offs[ff] - p->out_off + k;	/* Overall output trans param */

				vv += dout_de[0][ff] * dout_ov[jj];
				dav[p->out_off + jj] += p->points[i].w * vv;
			}
		}

		/* Compute and accumulate partial difference values for each parameter value */
		if (p->opt_msk & oc_i) {
			/* Input transfer parameters */
			for (ee = 0; ee < di; ee++) {				/* Parameter input chanel */
				for (k = 0; k < p->iluord[ee]; k++) {	/* Param within channel */
					double vv = 0.0;
					jj = p->in_offs[ee] - p->in_off + k;	/* Overall input trans param */

					for (ff = 0; ff < 3; ff++) {		/* Lab channels */
						vv += dout_de[0][ff] * dout_mato[ff]
						    * dmato_tin[ff * di + ee] * dtin_iv[jj];
					}
					dav[p->in_off + jj] += p->points[i].w * vv;
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
					dav[p->mat_off + jj] += p->points[i].w * vv;
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
					dav[p->out_off + jj] += p->points[i].w * vv;
				}
			}
		}
	}

	/* Normalise error to be an average delta E squared */
	rv /= tw;
	for (i = 0; i < p->tot_cnt; i++)
		dav[i] /= tw;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	rv += smv = dshapmag(p, dav);

	/* Copy the del for parameters being optimised to return array */
	for (i = 0; i < p->opt_cnt; i++)
		dv[i] = dav[p->opt_off + i];

#ifdef DEBUG
printf("~1(%f)dxfit1func returning %f\n",smv,rv);
#endif

	if (p->verb)
		printf("."), fflush(stdout);
	return rv;
}

#ifdef NEVER
/* Check partial derivative function within xfit1func() */

static double _xfit1func(void *edata, double *v) {
	xfit1 *p = (xfit1 *)edata;
	int i;
	double dv[MXPARMS];
	double rv, drv;
	double trv;
	int verb;
	
	rv = xfit1func(edata, v);
	verb = p->verb;
	p->verb = 0;
	drv = dxfit1func(edata, dv, v);
	p->verb = verb;

	if (fabs(rv - drv) > 1e-6)
		printf("######## RV MISMATCH is %f should be %f ########\n",rv,drv);

	/* Check each parameter delta */
	for (i = 0; i < p->opt_cnt; i++) {
		double del;

		v[i] += 1e-7;
		trv = xfit1func(edata, v);
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

#define xfit1func _xfit1func
#endif


/* Output curve symetry optimisation function handed to powell() */
/* Just the order 0 value will be adjusted */
static double symoptfunc(void *edata, double *v) {
	xfit1 *p = (xfit1 *)edata;
	double out[1], in[1] = { 0.0 };
	int ch = p->symch;		/* Output channel being adjusted for symetry */
	double rv;

	/* Copy the parameter being tested back into xfit1 */
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
static void setup_xfit1(
xfit1 *p,
double *wv,		/* Return parameters to hand to optimiser */
double *sa,		/* Return search radius to hand to optimiser */
double transrad,/* Nominal transfer curve radius, 0.0 - 3.0 */
double matrad	/* Nominal matrix radius, 0.0 - 1.0 */
) {
	int i;
	p->opt_off = -1;
	p->opt_cnt = 0;
	
	if (p->flags & XFIT_OUT_DEV)	/* scale search radius */
		matrad *= 100.0;

	if (p->opt_msk & oc_i) {
		p->opt_off = p->in_off;
		p->opt_cnt += p->in_cnt;

		for (i = 0; i < p->in_cnt; i++) {
			*wv++ = p->v[p->in_off + i];
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
#ifdef NEVER
printf("~1 opt_msk = 0x%x\n",p->opt_msk);
printf("~1 opt_off = %d\n",p->opt_off);
printf("~1 opt_cnt = %d\n\n",p->opt_cnt);
#endif /* NEVER */
}

#ifdef DEBUG
/* Diagnostic */
static void dump_xfit1(
xfit1 *p
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

/* Do the fitting. */
/* return nz on error */
/* 1 = malloc or other error */
int xfit1_fit(
	xfit1 *p,
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
	int dumy[],				/* Dummy array - not used */
	int oord[],				/* Order of output shaper curve for each dimension */
	optcomb tcomb,			/* Target elements to fit. */
	void *cntx2,			/* Context of callback */
	double (*to_de2)(void *cntx2, double *in1, double *in2)	
							/* callback to convert device value to fit metric */
) {
	int i, e, f;
	double *b;			/* Base of parameters for this section */
	int poff;

	if (flags & XFIT_FM_INPUT) {
		warning("xfit1 can't handle XFIT_FM_INPUT flag!");
		return 1;
	}

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
	p->cntx2   = cntx2;
	p->to_de2  = to_de2;

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
	p->in_off = 0 + p->mat_cnt;
	for (poff = p->in_off, p->in_cnt = 0, e = 0; e < di; e++) {
		p->in_offs[e] = poff;
		p->in_cnt += p->iluord[e];
		poff += p->iluord[e];
	}

	p->mat_off = p->in_off + p->in_cnt;
	p->mat_cnt = fdi * (1 << di);

	p->out_off = p->mat_off + p->mat_cnt;
	for (poff = p->out_off, p->out_cnt = 0, f = 0; f < fdi; f++) {
		p->out_offs[f] = poff;
		p->out_cnt += p->oluord[f];
		poff += p->oluord[f];
	}

	p->tot_cnt = p->mat_cnt + p->in_cnt + p->out_cnt;

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

	/* Create copy of input points with output converted to Lab */
	if (p->points != NULL) {
		free(p->points);
		p->points = NULL;
	}
	if ((p->points = (cow *)calloc(nodp, sizeof(cow))) == NULL)
		return 1;

	for (i = 0; i < p->nodp; i++) {
		p->points[i].w = p->ipoints[i].w; 
		for (e = 0; e < di; e++)
			p->points[i].p[e] = p->ipoints[i].p[e]; 
		for (f = 0; f < fdi; f++)
			p->points[i].v[f] = p->ipoints[i].v[f]; 
	
		/* Convert value to Lab if CIE output */

		if ((p->flags & XFIT_FM_MASK) == XFIT_FM_XYZ)
			icmLab2XYZ(&icmD50, p->points[i].v, p->points[i].v);
	}

	/* Setup input curves to be linear initially */
	b = p->v + p->in_off;
		
	for (e = 0; e < di; b += p->iluord[e], e++) {
		for (i = 0; i < p->iluord[e]; i++) {
			b[i] = 0.0;
		}
	}

	/* Setup matrix to be pure colorant values initially */
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
					tt = p->in_max[j] - p->points[k].p[j];
				else
					tt = p->in_min[j] - p->points[k].p[j];
				dif += tt * tt;
			}
			if (dif < bdif) {		/* best so far */
				bdif = dif;
				bk = k;
				if (dif < 0.001)
					break;			/* Don't bother looking further */
			}
		}

		for (f = 0; f < fdi; f++) {
			ov[f] = p->points[bk].v[f];
		}

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

	/* Do the fitting one part at a time, then together */
	if ((p->tcomb & oc_io) != 0
	 && (p->tcomb & oc_m) == oc_m) {	/* Only bother with matrix if in and/or out */

		if (p->verb)
			printf("About to optimise temporary matrix\n");

#ifdef DEBUG
printf("\nBefore matrix opt:\n");
dump_xfit1(p);
#endif
		/* Optimise matrix on its own */
		p->opt_msk = oc_m;
		setup_xfit1(p, p->wv, p->sa, 0.0, 0.5); 

#ifdef NODDV
		if (powell(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit1func, (void *)p) != 0)
			warning("set_icxLuLut: Powell failed to converge");
#else
		if (conjgrad(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit1func, dxfit1func, (void *)p) != 0)
			warning("set_icxLuLut: Conjgrad failed to converge");
#endif
		for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
			p->v[p->opt_off + i] = p->wv[i];

#ifdef DEBUG
printf("\nAfter matrix opt:\n");
dump_xfit1(p);
#endif
	}

	/* Optimise input and matrix together */
	if ((p->tcomb & oc_im) == oc_im) {

		if (p->verb)
			printf("\nAbout to optimise input curves and matrix\n");

		p->opt_msk = oc_im;
		setup_xfit1(p, p->wv, p->sa, 0.5, 0.3); 
#ifdef NODDV
		if (powell(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit1func, (void *)p) != 0)
			warning("set_icxLuLut: Powell failed to converge");
#else
		if (conjgrad(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit1func, dxfit1func, (void *)p) != 0)
			warning("set_icxLuLut: Conjgrad failed to converge");
#endif
		for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
			p->v[p->opt_off + i] = p->wv[i];
#ifdef DEBUG
printf("\nAfter input and matrix opt:\n");
dump_xfit1(p);
#endif
	}

	/* Optimise the matrix and output curves together */
	if ((p->tcomb & oc_mo) == oc_mo) {

		if (p->verb)
			printf("\nAbout to optimise output curves and matrix\n");

		p->opt_msk = oc_mo;
		setup_xfit1(p, p->wv, p->sa, 0.3, 0.3); 
#ifdef NODDV
		if (powell(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit1func, (void *)p) != 0)
			warning("set_icxLuLut: Powell failed to converge");
#else
		if (conjgrad(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit1func, dxfit1func, (void *)p) != 0)
			warning("set_icxLuLut: Conjgrad failed to converge");
#endif
		for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
			p->v[p->opt_off + i] = p->wv[i];
#ifdef DEBUG
printf("\nAfter output opt:\n");
dump_xfit1(p);
#endif


		/* Optimise input and matrix together again, after altering matrix */
		if ((p->tcomb & oc_im) == oc_im) {

			if (p->verb)
				printf("\nAbout to optimise input curves and matrix again\n");

			p->opt_msk = oc_im;
			setup_xfit1(p, p->wv, p->sa, 0.2, 0.2); 
#ifdef NODDV
			if (powell(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit1func, (void *)p) != 0)
				warning("set_icxLuLut: Powell failed to converge");
#else
			if (conjgrad(NULL, p->opt_cnt, p->v, p->sa, POWTOL, MAXITS, xfit1func, dxfit1func, (void *)p) != 0)
				warning("set_icxLuLut: Conjgrad failed to converge");
#endif
			for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
				p->v[p->opt_off + i] = p->wv[i];
		}
#ifdef DEBUG
printf("\nAfter 2nd input and matrix opt:\n");
dump_xfit1(p);
#endif

#ifndef NODDV
		/* Optimise all together */
		if ((p->tcomb & oc_imo) == oc_imo) {

			if (p->verb)
				printf("\nAbout to optimise input, matrix and output together\n");

			p->opt_msk = oc_imo;
			setup_xfit1(p, p->wv, p->sa, 0.1, 0.1); 
			if (conjgrad(NULL, p->opt_cnt, p->wv, p->sa, POWTOL, MAXITS, xfit1func, dxfit1func, (void *)p) != 0)
				warning("set_icxLuLut: Conjgrad failed to converge");
			for (i = 0; i < p->opt_cnt; i++)		/* Copy optimised values back */
				p->v[p->opt_off + i] = p->wv[i];
		}

#ifdef DEBUG
printf("\nAfter all together opt:\n");
dump_xfit1(p);
#endif

#endif /* !NODDV */

		/* Adjust output curve white point */
		if (p->flags & XFIT_OUT_ZERO) {

			if (p->verb)
				printf("\nAbout to adjust a and b output curves for white point\n");

			for (f = 1; f < 3; f++) {
				p->symch = f;
				p->wv[0] = p->v[p->out_offs[f]];	/* Current parameter value */
				p->sa[0] = 0.1;					/* Search radius */
				if (powell(NULL, 1, p->wv, p->sa, 0.0000001, 1000, symoptfunc, (void *)p) != 0)
					error("set_icxLuLut: Powell failed to converge");
				p->v[p->out_offs[f]] = p->wv[0];	/* Copy results back */
			}
		}
	}

	if (p->verb)
		printf("\n");

#ifdef DEBUG
printf("\nFinal parameters:\n");
dump_xfit1(p);
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

	return 0;
}

/* We're done with an xfit1 */
static void xfit1_del(xfit1 *p) {
	if (p->v != NULL)
		free(p->v);
	if (p->wv != NULL)
		free(p->wv);
	if (p->sa != NULL)
		free(p->sa);
	if (p->points != NULL)
		free(p->points);
	free(p);
}

/* Create a transform fitting object */
/* return NULL on error */
xfit1 *new_xfit1(
) {
	xfit1 *p;

	if ((p = (xfit1 *)calloc(1, sizeof(xfit1))) == NULL) {
		return NULL;
	}

	/* Set method pointers */
	p->fit         = xfit1_fit;
	p->incurve     = xfit1_incurve;
	p->invincurve  = xfit1_invincurve;
	p->outcurve    = xfit1_outcurve;
	p->invoutcurve = xfit1_invoutcurve;
	p->del = xfit1_del;

	return p;
}



