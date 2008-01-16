
/*
 * International Color Consortium color transform expanded support
 *
 * Author:  Graeme W. Gill
 * Date:    2/7/00
 * Version: 1.00
 *
 * Copyright 2000, 2001 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * This is the third major version of xlut.c (originally called xlut2.c)
 * Based on the old xlut.c code (preserved as xlut1.c)
 * This version uses xfit.c to do the curve and rspl fitting.
 */

/*
 * This module provides the expanded icclib functionality
 * for lut based profiles.
 * This file is #included in xicc.c, to keep its functions private.
 *
 * This version creates both input and output 1D luts by
 * optimising the accuracy of the profile for a linear clut.
 *
 */

#include "xfit.h"

#define USE_CIE94_DE	/* Use CIE94 delta E measure when creating in/out curves */

#undef DEBUG 			/* Verbose debug information */
#undef DEBUG_RLUT 		/* Print values being reverse lookup up */
#undef INK_LIMIT_TEST	/* Turn off input tables for ink limit testing */
#undef CHECK_ILIMIT		/* Do sanity checks on meeting ink limit */

#undef NODDV			/* Use slow non d/dv powell */

#define POWTOL 5e-4			/* Shaper Powell optimiser tollerance */
#define MAXITS 2000			/* Shaper number of itterations before giving up */

/* Weights in shaper parameters, to minimise unconstrained "wiggles" */
#define SHAPE_WEIGHT 5000.0		/* Overal shaper weight contribution - err on side of smoothness */
#define SHAPE_BASE  0.00001		/* 0 & 1 harmonic weight */
#define SHAPE_HBASE 0.0001		/* 2nd and higher additional weight */

/*
 * TTBD:
 *
 *       Reverse lookup of Lab
 *       Make NEARCLIP the default ??
 *
 *       XYZ vector clipping isn't implemented.
 *
 *       Some of the error handling is crude. Shouldn't use
 *       error(), should return status.
 */

#ifndef _CAT2
#define _CAT2(n1,n2)  n1 ## n2
#define CAT2(n1,n2) _CAT2(n1,n2)
#endif



static double icxLimitD(void *lcntx, double *in);		/* For input' */
static double icxLimit(void *lcntx, double *in);		/* For input */
static int icxLuLut_init_clut_camclip(icxLuLut *p);

/* ========================================================== */
/* xicc lookup code.                                          */
/* ========================================================== */

/* Forward and Backward Multi-Dimensional Interpolation type conversion */
/* Return 0 on success, 1 if clipping occured, 2 on other error */

/* Components of overall lookup, in order */

int icxLuLut_in_abs(icxLuLut *p, double *out, double *in) {
	int rv = 0;

	if (p->ins == icxSigJabData) {
		p->cam->cam_to_XYZ(p->cam, out, in);
		rv |= ((icmLuLut *)p->plu)->in_abs((icmLuLut *)p->plu, out, out);
	} else {
		rv |= ((icmLuLut *)p->plu)->in_abs((icmLuLut *)p->plu, out, in);
	}

	return rv;
}

/* Possible matrix lookup */
/* input->input (not distinguishing matrix altered input values) */
int icxLuLut_matrix(icxLuLut *p, double *out, double *in) {
	int rv = 0;
	rv |= ((icmLuLut *)p->plu)->matrix((icmLuLut *)p->plu, out, in);
	return rv;
}

/* Do input -> input' lookup */
int icxLuLut_input(icxLuLut *p, double *out, double *in) {
#ifdef NEVER
	return ((icmLuLut *)p->plu)->input((icmLuLut *)p->plu, out, in);
#else
	int rv = 0;
	co tc;
	int i;
	for (i = 0; i < p->inputChan; i++) {
		tc.p[0] = in[i];
		rv |= p->inputTable[i]->interp(p->inputTable[i], &tc);
		out[i] = tc.v[0];
	}
	return rv;
#endif
}

/* Do input'->output' lookup, with aux return */
int icxLuLut_clut_aux(icxLuLut *p,
double *out,	/* output' value */
double *oink,	/* If not NULL, return amount input is over the ink limit, 0 if not */
double *auxv,	/* If not NULL, return aux value used (packed) */
double *in		/* input' value */
) {
	int rv = 0;
	co tc;
	int i;

	for (i = 0; i < p->inputChan; i++)
		tc.p[i] = in[i];
	rv |= p->clutTable->interp(p->clutTable, &tc);
	for (i = 0; i < p->outputChan; i++)
		out[i] = tc.v[i];

	if (auxv != NULL) {
		int ee = 0;
		for (i = 0; i < p->clutTable->di; i++) {
			double v = in[i];
			if (p->auxm[i] != 0) {
				auxv[ee] = v;
				ee++;
			}
		}
	}

	if (oink != NULL) {
		double lim = 0.0;

		if (p->ink.tlimit >= 0.0 || p->ink.klimit >= 0.0) {
			lim = icxLimitD((void *)p, in);
			if (lim < 0.0)
				lim = 0.0;
		}
		*oink = lim;
	}

	return rv;
}

/* Do input'->output' lookup */
int icxLuLut_clut(icxLuLut *p, double *out, double *in) {
#ifdef NEVER
	return ((icmLuLut *)p->plu)->clut((icmLuLut *)p->plu, out, in);
#else
	return icxLuLut_clut_aux(p, out, NULL, NULL, in);
#endif
}

/* Do output'->output lookup */
int icxLuLut_output(icxLuLut *p, double *out, double *in) {
	int rv = 0;

	if (p->mergeclut == 0) {
#ifdef NEVER
		rv = ((icmLuLut *)p->plu)->output((icmLuLut *)p->plu, out, in);
#else
		co tc;
		int i;
		for (i = 0; i < p->outputChan; i++) {
			tc.p[0] = in[i];
			rv |= p->outputTable[i]->interp(p->outputTable[i], &tc);
			out[i] = tc.v[0];
		}
#endif
	} else {
		int i;
		for (i = 0; i < p->outputChan; i++)
			out[i] = in[i];
	}
	return rv;
}

/* Relative to absolute conversion + PCS to PCS override (Effective PCS) conversion */
int icxLuLut_out_abs(icxLuLut *p, double *out, double *in) {
	int rv = 0;
	if (p->mergeclut == 0) {
		rv |= ((icmLuLut *)p->plu)->out_abs((icmLuLut *)p->plu, out, in);

		if (p->outs == icxSigJabData) {
			p->cam->XYZ_to_cam(p->cam, out, out);
		}
	} else {
		int i;
		for (i = 0; i < p->outputChan; i++)
			out[i] = in[i];
	}

	return rv;
}

/* Overall lookup */
static int
icxLuLut_lookup (
icxLuBase *pp,		/* This */
double *out,		/* Vector of output values */
double *in			/* Vector of input values */
) {
	icxLuLut *p = (icxLuLut *)pp;
	int rv = 0;
	double temp[MAX_CHAN];

	rv |= p->in_abs  (p, temp, in);
	rv |= p->matrix  (p, temp, temp);
	rv |= p->input   (p, temp, temp);
	rv |= p->clut    (p, out,  temp);
	if (p->mergeclut == 0) {
		rv |= p->output  (p, out,  out);
		rv |= p->out_abs (p, out,  out);
	}
	return rv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Given a relative XYZ or Lab PCS value, convert in the fwd direction into */ 
/* the nominated output PCS (ie. Absolute, Jab etc.) */
/* (This is used in generating gamut compression in B2A tables) */
void icxLuLut_fwd_relpcs_outpcs(
icxLuBase *pp,
icColorSpaceSignature is,		/* Input space, XYZ or Lab */
double *out, double *in) {
	icxLuLut *p = (icxLuLut *)pp;

	if (is == icSigLabData && p->natpcs == icSigXYZData) {
		icmLab2XYZ(&icmD50, out, in);
		icxLuLut_out_abs(p, out, out);
	} else if (is == icSigXYZData && p->natpcs == icSigLabData) {
		icmXYZ2Lab(&icmD50, out, in);
		icxLuLut_out_abs(p, out, out);
	} else {
		icxLuLut_out_abs(p, out, in);
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Components of inverse lookup, in order */

#ifdef DEBUG_RLUT
#undef DBR
#define DBR(xxx) printf xxx ;
#else
#undef DBR
#define DBR(xxx) 
#endif

/* Utility function - compute the clip vector direction. */
/* return NULL if vector clip isn't used. */
double *icxClipVector(
icxClip *p,			/* Clipping setup information */
double *in,			/* Target point */
double *cdirv		/* Space for returned clip vector */
) {
	int f;
	if (p->nearclip != 0)
		return NULL;			/* Doing nearest clipping, not vector */

	/* Default is simple vector clip */
	for (f = 0; f < p->fdi; f++)
		cdirv[f] = p->ocent[f] - in[f];	/* Clip towards output gamut center */

	if (p->ocentl != 0.0) {			/* Graduated vector clip */
		double cvl, nll;

		/* Compute (negative) clip vector length */
		for (cvl = 0.0, f = 0; f < p->fdi; f++) {
			cvl += cdirv[f] * cdirv[f];
		}
		cvl = sqrt(cvl);
		if (cvl > 1e-8) {
			/* Dot product of clip vector and clip center line */
			for (nll = 0.0, f = 0; f < p->fdi; f++)
				nll -= cdirv[f] * p->ocentv[f];	/* (Fix -ve clip vector) */
			nll /= (p->ocentl * p->ocentl);	/* Normalised location along line */

			/* Limit to line */
			if (nll < 0.0)
				nll = 0.0;
			else if (nll > 1.0)
				nll = 1.0;

			if (p->LabLike) {
				/* Aim more towards center for saturated targets */
				double sat = sqrt(in[1] * in[1] + in[2] * in[2]);
				nll += sat/150.0 * (0.5 - nll);
			}

			/* Compute target clip direction */
			for (f = 0; f < p->fdi; f++)
				cdirv[f] = p->ocent[f] + nll * p->ocentv[f] - in[f];
		}
	}

	return cdirv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* PCS override (Effective PCS) to PCS conversion + absolute to relative conversion */
int icxLuLut_inv_out_abs(icxLuLut *p, double *out, double *in) {
	int rv = 0;

	DBR(("icxLuLut_inv_out_abs got PCS %f %f %f\n",in[0],in[1],in[2]))

	if (p->mergeclut == 0) {
		if (p->outs == icxSigJabData) {
			p->cam->cam_to_XYZ(p->cam, out, in);
			rv |= ((icmLuLut *)p->plu)->inv_out_abs((icmLuLut *)p->plu, out, out);
		} else {
			rv |= ((icmLuLut *)p->plu)->inv_out_abs((icmLuLut *)p->plu, out, in);
		}
	} else {
		int i;
		for (i = 0; i < p->outputChan; i++)
			out[i] = in[i];
	}
	DBR(("icxLuLut_inv_out_abs returning PCS %f %f %f\n",out[0],out[1],out[2]))
	return rv;
}

/* Do output->output' inverse lookup */
int icxLuLut_inv_output(icxLuLut *p, double *out, double *in) {
	int rv = 0;
	DBR(("icxLuLut_inv_output got PCS %f %f %f\n",in[0],in[1],in[2]))
	if (p->mergeclut == 0) {
#ifdef NEVER
		rv = ((icmLuLut *)p->plu)->inv_output((icmLuLut *)p->plu, out, in);
#else
		int i,j;
		int nsoln;				/* Number of solutions found */
		co pp[MAX_INVSOLN];		/* Room for all the solutions found */
		double cdir;

		for (i = 0; i < p->outputChan; i++) {
			pp[0].p[0] = p->outputClipc[i];
			pp[0].v[0] = in[i];
			cdir = p->outputClipc[i] - in[i];	/* Clip towards output range */

			nsoln = p->outputTable[i]->rev_interp (
				p->outputTable[i], 	/* this */
				RSPL_NEARCLIP,		/* Clip to nearest (faster than vector) */
				MAX_INVSOLN,		/* Maximum number of solutions allowed for */
				NULL, 				/* No auxiliary input targets */
				&cdir,				/* Clip vector direction and length */
				pp);				/* Input and output values */

			if (nsoln & RSPL_DIDCLIP)
				rv = 1;

			nsoln &= RSPL_NOSOLNS;		/* Get number of solutions */

			if (nsoln == 1) { /* Exactly one solution */
				j = 0;
			} else if (nsoln == 0) {	/* Zero solutions. This is unexpected. */
				error("xlut: Unexpected failure to find reverse solution for output table");
				return 2;
			} else {		/* Multiple solutions */
				/* Use a simple minded resolution - choose the one closest to the center */
				double bdist = 1e300;
				int bsoln = 0;
/* Don't expect this - 1D luts are meant to be monotonic */
printf("~~~1 got %d reverse solutions\n",nsoln);
printf("~~~1 solution 0 = %f\n",pp[0].p[0]);
printf("~~~1 solution 1 = %f\n",pp[1].p[0]);
				for (j = 0; j < nsoln; j++) {
					double tt;
					tt = pp[i].p[0] - p->outputClipc[i];
					tt *= tt;
					if (tt < bdist) {	/* Better solution */
						bdist = tt;
						bsoln = j;
					}
				}
				j = bsoln;
			}
			out[i] = pp[j].p[0];
		}

#endif /* NEVER */
	} else {
		int i;
		for (i = 0; i < p->outputChan; i++)
			out[i] = in[i];
	}
	DBR(("icxLuLut_inv_output returning PCS' %f %f %f\n",out[0],out[1],out[2]))
	return rv;
}

/* Ink limit+gamut limit calculation function for xLuLut. */
/* Returns < 0.0 if input value is within limits, */
/* > 0.0 if outside limits. Value is proportinal to distance to limits. */
/* We implement device gamut check to improve utility outside rspl, */
/* in optimisation routines. */
static double icxLimit(
void *lcntx,
double *in
) {
	icxLuLut *p = (icxLuLut *)lcntx;
	double tlim, klim;
	double ovr, val;
	int e;

//printf("~~ limit got %f %f %f %f\n", in[0], in[1], in[2], in[3]);

	if ((tlim = p->ink.tlimit) < 0.0)
		tlim = (double)p->inputChan + 1.0;

	if ((klim = p->ink.klimit) < 0.0)
		klim = 2.0;

	/* Compute amount outside total limit */
	{
		double sum;
		for (sum = 0.0, e = 0; e < p->inputChan; e++)
			sum += in[e];
		val = sum - tlim;
	}

	/* Compute amount outside black limit */
	if (p->ink.klimit >= 0.0) {
		double kval = 0.0;
		switch(p->natis) {
			case icSigCmykData:
				kval = in[3] - klim;
				break;
			default:
				error("xlut: Unknown colorspace when black limit specified");
		}
		if (kval > val)
			val = kval;
	}

	/* Compute amount outside device value limits 0.0 - 1.0 */
	for (ovr = 0.0, e = 0; e < p->inputChan; e++) {
		if (in[e] < 0.0) {
			if (-in[e] > ovr)
				ovr = -in[e];
		} else if (in[e] > 1.0) {
			if ((in[e] - 1.0) > ovr)
				ovr = in[e] - 1.0;
		}
	}
	if (ovr > val)
		val = ovr;

//printf("~~ limit returning %f\n", val);
	return val;
}

/* Same as above, but works with input' values */
static double icxLimitD(
void *lcntx,
double *ind
) {
	icxLuLut *p = (icxLuLut *)lcntx;
	double in[MAX_CHAN];
	co tc;
	int e;

	/* Convert input' to input through revinput Luts */
	for (e = 0; e < p->inputChan; e++) {
		tc.p[0] = ind[e];
		p->revinputTable[e]->interp(p->revinputTable[e], &tc);
		in[e] = tc.v[0];
	}

	return icxLimit(lcntx, in);
}



/* helper function that creates our standard K locus curve value, */
/* given the curve parameters, and the normalised L value. */
/* No filtering version */
static double icxKcurveNF(double L, icxInkCurve *c) {
	double rv;

	/* Constrain L to 0.0 .. 1.0, just in case... */
	if (L < 0.0)
		L = 0.0;
	else if (L > 1.0)
		L = 1.0;

	/* Invert sense of L, so that 0.0 = white, 1.0 = black */
	L = 1.0 - L;

	if (L <= c->Kstpo && L <= c->Kenpo) {
		/* We are at white level */
		rv = c->Kstle;
	} else if (L >= c->Kstpo && L >= c->Kenpo) {
		/* We are at black level */
		rv = c->Kenle;
	} else {
		double g;
		/* We must be on the curve from start to end levels */

		if (c->Kstpo > c->Kenpo) {
			rv = (L - c->Kenpo)/(c->Kstpo - c->Kenpo);
		} else {
			rv = (L - c->Kstpo)/(c->Kenpo - c->Kstpo);
		}

		g = c->Kshap/2.0;

		/* A value of 0.5 will be tranlated to g */
		rv = rv/((1.0/g - 2.0) * (1.0 - rv) + 1.0);

		/* Transition between start end end levels */
		rv = rv * (c->Kenle - c->Kstle) + c->Kstle;
	}

	/* To be safe */
	if (rv < 0.0)
		rv = 0.0;
	else if (rv > 1.0)
		rv = 1.0;

	return rv;
}

/* Same as above, but implement a convolution filter */
static double icxKcurve(double L, icxInkCurve *c) {
	double rv = 0.0;

	if (c->Ksmth == 0.0)
		return icxKcurveNF(L, c);

	/* Use a simple triangle filter */
	rv  =  0.20 * icxKcurveNF(L - 1.0 * c->Ksmth, c);
	rv +=  0.44 * icxKcurveNF(L - 0.7 * c->Ksmth, c);
	rv +=  0.76 * icxKcurveNF(L - 0.3 * c->Ksmth, c);

	rv +=  1.0  * icxKcurveNF(L,                  c);

	rv +=  0.76 * icxKcurveNF(L + 0.3 * c->Ksmth, c);
	rv +=  0.44 * icxKcurveNF(L + 0.7 * c->Ksmth, c);
	rv +=  0.20 * icxKcurveNF(L + 1.0 * c->Ksmth, c);

	return rv/3.8;
}

/* Do output'->input' lookup with aux details. */
/* Note than out[] will be used as the inking value if icxKrule is */
/* icxKvalue, icxKlocus, icxKl5l or icxKl5lk, and that the icxKrule value will */
/* be in the input' space. Note that the ink limit will be computed after converting */
/* input' to input, auxt will override the inking rule, and auxr reflects the */
/* available auxiliary range there was to choose from. auxv was the actual auxiliary used. */
/* Note that the aux values will reflect linearised auxiliary values. */
/* Retuns clip status. */
int icxLuLut_inv_clut_aux(
icxLuLut *p,
double *out,	/* Function return values, plus aux value or locus target input if auxt == NULL */
double *auxv,	/* If not NULL, return aux value used (packed) */
double *auxr,	/* If not NULL, return aux locus range (packed, 2 at a time) */
double *auxt,	/* If not NULL, specify the aux target for this lookup (override ink) */
double *in		/* Function input values to invert (== clut output' values) */
) {
	co pp[MAX_INVSOLN];		/* Room for all the solutions found */
	int nsoln;			/* Number of solutions found */
	double *cdir, cdirv[MXDO];	/* Clip vector direction and length */
	int e,f,i;
	int fdi = p->clutTable->fdi;
	int flags = 0;		/* reverse interp flags */
	int xflags = 0;		/* extra clip/noclip flags */
	int crv = 0;		/* Return value - set to 1 if clipped */

	if (p->nearclip != 0)
		flags |= RSPL_NEARCLIP;			/* Use nearest clipping rather than clip vector */

	DBR(("inv_clut_aux input is %f %f %f\n",in[0], in[1], in[2]))

	if (auxr != NULL) {		/* Set a default locus range */
		int ee = 0;
		for (e = 0; e < p->clutTable->di; e++) {
			if (p->auxm[e] != 0) {
				auxr[ee++] = 1e60;
				auxr[ee++] = -1e60;
			}
		}
	}

	/* Setup for reverse lookup */
	for (f = 0; f < fdi; f++)
		pp[0].v[f] = in[f];				/* Target value */

	/* Compute clip vector, if any */
	cdir = icxClipVector(&p->clip, in, cdirv);

	if (p->clutTable->di > fdi) {	/* ie. CMYK->Lab, there will be ambiguity */
		double min[MXDI], max[MXDI];	/* Auxiliary locus range */

		/* Compute auxiliary locus on the fly. This is in dev' == input' space. */
		nsoln = p->clutTable->rev_locus(
			p->clutTable,	/* rspl object */
			p->auxm,		/* Auxiliary mask */
			pp,				/* Input target and output solutions */
			min, max);		/* Returned locus of valid auxiliary values */
		
		if (nsoln == 0) {
			xflags |= RSPL_WILLCLIP;	/* No valid locus, so we expect to have to clip */
		}

		else {  /* Got a valid locus */

			if (auxr != NULL) {		/* Report the locus range */
				int ee = 0;
				for (e = 0; e < p->clutTable->di; e++) {
					if (p->auxm[e] != 0) {
						auxr[ee++] = min[e];
						auxr[ee++] = max[e];
					}
				}
			}

			if (auxt != NULL) {		/* overiding auxiliary target */
				int ee = 0;
				for (e = 0; e < p->clutTable->di; e++) {
					if (p->auxm[e] != 0) {
						double iv = auxt[ee++];
						if (iv < min[e])
							iv = min[e];
						else if (iv > max[e])
							iv = max[e];
						pp[0].p[e] = iv;
					}
				}
			} else if (p->ink.k_rule == icxKvalue) {
				/* Implement the auxiliary inking rule */
				/* Target auxiliary values are provided in out[] */
				for (e = 0; e < p->clutTable->di; e++) {
					if (p->auxm[e] != 0) {
						double iv = out[e];		/* out[] holds aux target value */
						if (iv < min[e])
							iv = min[e];
						else if (iv > max[e])
							iv = max[e];
						pp[0].p[e] = iv;
					}
				}
			} else if (p->ink.k_rule == icxKlocus) {
				/* Set target auxliary input values from values in out[] and locus */
				for (e = 0; e < p->clutTable->di; e++) {
					if (p->auxm[e] != 0) {
						double ii, iv;
						ii = out[e];							/* Input ink locus */
						iv = min[e] + ii * (max[e] - min[e]);	/* Output ink from locus */
						if (iv < min[e])
							iv = min[e];
						else if (iv > max[e])
							iv = max[e];
						pp[0].p[e] = iv;
					}
				}
			} else { /* p->ink.k_rule == icxKluma5 || icxKluma5k || icxKl5l || icxKl5lk */
				/* Auxiliaries are driven by a rule and the output values */
				double tin[MAX_CHAN], rv, L;

				/* If we've got a mergeclut, then the PCS' is the same as the */
				/* effective PCS, and we need to convert to native PCS */
				if (p->mergeclut) {
					p->mergeclut = 0;					/* Hack to be able to use inv_out_abs() */
					icxLuLut_inv_out_abs(p, tin, in);
					p->mergeclut = 1;

				} else {
					/* Convert native PCS' to native PCS values */
					p->output(p, tin, in);	
				}

				/* Figure out Luminance number */
				if (p->natos == icSigXYZData) {
					icmXYZ2Lab(&icmD50, tin, tin);

				} else if (p->natos != icSigLabData) {	/* Hmm. that's unexpected */
					error("Assert: xlut K locus, unexpected native pcs of 0x%x\n",p->natos);
				}
				L = 0.01 * tin[0];

				/* Normalise L to its possible range from min to max */
				L = (L - p->Lmin)/(p->Lmax - p->Lmin);

				/* Convert L to curve value */
				rv = icxKcurve(L, &p->ink.c);

				if (p->ink.k_rule == icxKluma5) {	/* Curve is locus value */

					/* Set target black as K fraction within locus */

					for (e = 0; e < p->clutTable->di; e++) {
						if (p->auxm[e] != 0) {
							pp[0].p[e] = min[e] + rv * (max[e] - min[e]);
						}
					}

				} else if (p->ink.k_rule == icxKluma5k) {	/* Curve is K value */

					for (e = 0; e < p->clutTable->di; e++) {
						if (p->auxm[e] != 0) {
							double iv = rv;
							if (iv < min[e])			/* Clip to locus */
								iv = min[e];
							else if (iv > max[e])
								iv = max[e];
							pp[0].p[e] = iv;
						}
					}

				} else { /* icxKl5l || icxKl5lk */
					/* Create second curve, and use input locus to */
					/* blend between */

					double rv2;		/* Upper limit */

					/* Convert L to max curve value */
					rv2 = icxKcurve(L, &p->ink.x);

					if (rv2 < rv) {		/* Ooops - better swap. */
						double tt;
						tt = rv;
						rv = rv2;
						rv2 = tt;
					}

					for (e = 0; e < p->clutTable->di; e++) {
						if (p->auxm[e] != 0) {
							if (p->ink.k_rule == icxKl5l) {
								double ii;
								ii = out[e];				/* Input K locus */
								if (ii < 0.0)
									ii = 0.0;
								else if (ii > 1.0)
									ii = 1.0;
								ii = (1.0 - ii) * rv + ii * rv2;/* Blend between locus rule curves */
								/* Out ink from output locus */
								pp[0].p[e] = min[e] + ii * (max[e] - min[e]);
							} else {
								double iv;
								iv = out[e];				/* Input K level */
								if (iv < rv)				/* Constrain to curves */
									iv = rv;
								else if (iv > rv2)
									iv = rv2;
								pp[0].p[e] = iv; 
							}
						}
					}
				}
			}

			xflags |= RSPL_EXACTAUX;	/* Since we confine aux to locus */
		}

#ifdef DEBUG_RLUT
		printf("inv_clut_aux computed aux values ");
		for (e = 0; e < p->clutTable->di; e++) {
			if (p->auxm[e] != 0)
				printf("%d: %f ",e,pp[0].p[e]);
		}
		printf("\n");
#endif /* DEBUG_RLUT */

		/* Find reverse solution with target auxiliaries */
		nsoln = p->clutTable->rev_interp(
			p->clutTable, 	/* rspl object */
			flags | xflags,	/* Combine all the flags */
			MAX_INVSOLN, 	/* Maxumum solutions to return */
			p->auxm, 		/* Auxiliary input chanel mask */
			cdir,			/* Clip vector direction and length */
			pp);			/* Input target and output solutions */
							/* returned solutions in pp[0..retval-1].p[] */

	} else {
		DBR(("inv_clut_aux needs no aux value\n"))

		/* Color spaces don't need auxiliaries to choose from solution locus */
		nsoln = p->clutTable->rev_interp(
			p->clutTable, 	/* rspl object */
			flags,			/* No extra flags */
			MAX_INVSOLN, 	/* Maxumum solutions to return */
			NULL, 			/* No auxiliary input targets */
			cdir,			/* Clip vector direction and length */
			pp);			/* Input target and output solutions */
							/* returned solutions in pp[0..retval-1].p[] */
	}
	if (nsoln & RSPL_DIDCLIP)
		crv = 1;			/* Clipped on PCS inverse lookup */

	nsoln &= RSPL_NOSOLNS;		/* Get number of solutions */

	DBR(("inv_clut_aux got %d rev_interp solutions, clipflag = %d\n",nsoln,crv))

	/* If we clipped and we should clip in CAM Jab space, compute reverse */
	/* clip solution using our additional CAM space. */
	/* (Note that we don't support vector clip in CAM space at the moment) */
	if (crv != 0 && p->camclip && p->nearclip) {
		double tin[MXDO];	/* CAM space value to be inverted */
		co cpp;				/* Alternate CAM space solution */
		double cdist;		/* CAM clip distance */
		double bf;			/* Blend factor */

		DBR(("inv_clut_aux got clip, compute CAM clip\n"))

		if (nsoln != 1) {	/* This would be unexpected */
			error("~~~1 Unexpected failure to return 1 solution on clip for input to output table");
		}

		if (p->cclutTable == NULL) {	/* we haven't created this yet, so do so */
			if (icxLuLut_init_clut_camclip(p))
				error("~~~1 Error in creating CAM rspl for camclip");
		}

		/* Setup for reverse lookup */
		DBR(("inv_clut_aux cam clip PCS in %f %f %f\n",in[0],in[1],in[2]))

		/* Convert from PCS' to (XYZ) PCS */
		((icmLuLut *)p->absxyzlu)->output((icmLuLut *)p->absxyzlu, tin, in);
		DBR(("inv_clut_aux cam clip PCS' -> PCS %f %f %f\n",tin[0],tin[1],tin[2]))

		((icmLuLut *)p->absxyzlu)->out_abs((icmLuLut *)p->absxyzlu, tin, tin);
		DBR(("inv_clut_aux cam clip abs XYZ PCS %f %f %f\n",tin[0],tin[1],tin[2]))

		p->cam->XYZ_to_cam(p->cam, tin, tin);
		DBR(("inv_clut_aux cam clip PCS after XYZtoCAM %f %f %f\n",tin[0],tin[1],tin[2]))

		for (f = 0; f < fdi; f++)
			cpp.v[f] = tin[f];

		if (p->clutTable->di > fdi) {	/* ie. CMYK->Lab, there will be ambiguity */

			nsoln = p->cclutTable->rev_interp(
				p->cclutTable, 	/* rspl object */
				flags | xflags | RSPL_WILLCLIP,	/* Combine all the flags + clip ?? */
				1, 				/* Maximum solutions to return */
				p->auxm, 		/* Auxiliary input chanel mask */
				cdir,			/* Clip vector direction and length */
				&cpp);			/* Input target and output solutions */

		} else {
			nsoln = p->cclutTable->rev_interp(
				p->cclutTable, 	/* rspl object */
				flags | RSPL_WILLCLIP,	/* Because we know it will clip ?? */
				1, 				/* Maximum solutions to return */
				NULL, 			/* No auxiliary input targets */
				cdir,			/* Clip vector direction and length */
				&cpp);			/* Input target and output solutions */
		}

		nsoln &= RSPL_NOSOLNS;		/* Get number of solutions */

		if (nsoln != 1) {	/* This would be unexpected */
			error("~~~1 Unexpected failure to return 1 solution on CAM clip for input to output table");
		}

		/* Compute the CAM clip distances */
		for (cdist = 0.0, f = 0; f < fdi; f++) {
			double tt;
			tt = cpp.v[f] - tin[f];
			cdist += tt * tt;
		}
		cdist = sqrt(cdist);

		/* Use magic number to set blend distance, and compute a blend factor. */
		/* Blend over 4 delta E */
		bf = cdist/4.0;						/* 0.0 for PCS result, 1.0 for CAM result */
		if (bf > 1.0)
			bf = 1.0;
//printf("~1 raw blend = %f\n",bf);
		bf = bf * bf * (3.0 - 2.0 * bf);	/* Convert to spline blend */
		DBR(("cdist %f, spline blend %f\n",cdist,bf))

		/* Blend between solution values for PCS and CAM clip result. */
		/* We're hoping that the solutions are close, and expect them to be */
		/* that way when we're close to the gamut surface (since the cell */
		/* vertex values should be exact, irrespective of which space they're */
		/* in), but weird things could happen away from the surface. Weird */
		/* things can happen anyway with "clip to nearest", since this is not */
		/* guaranteed to be a smooth function, depending on the gamut surface */
		/* geometry, without taking some precaution such as clipping to a */
		/* convex hull "wrapper" to create a clip vector, which we're not */
		/* current doing. */
		DBR(("Clip blend between:\n"))
		DBR(("Lab: %f %f %f %f and\n",pp[0].p[0], pp[0].p[1], pp[0].p[2], pp[0].p[3]))
		DBR(("Jab: %f %f %f %f\n",cpp.p[0], cpp.p[1], cpp.p[2], cpp.p[3]))

		for (e = 0; e < p->clutTable->di; e++) {
			out[e] = (1.0 - bf) * pp[0].p[e] + bf * cpp.p[e];
		}

	/* Not CAM clip case */
	} else {

		if (nsoln == 1) { /* Exactly one solution */
			i = 0;
		} else if (nsoln == 0) {	/* Zero solutions. This is unexpected. */
			double in_v[MXDO];
			p->output(p, in_v, pp[0].v);		/* Get ICC inverse input values */
			p->out_abs(p, in_v, in_v);
			error("~~~1 Unexpected failure to find reverse solution for input to output table for value %f %f %f (ICC input %f %f %f)",pp[0].v[0],pp[0].v[1],pp[0].v[2], in_v[0], in_v[1], in_v[2]);
			return 2;
		} else {		/* Multiple solutions */
			/* Use a simple minded resolution - choose the one closest to the center */
			double bdist = 1e300;
			int bsoln = 0;
//printf("~~~1 got multiple reverse solutions\n");
			for (i = 0; i < nsoln; i++) {
				double ss;
				for (ss = 0.0, e = 0; e < p->clutTable->di; e++) {
					double tt;
					tt = pp[i].p[e] - p->icent[e];
					tt *= tt;
					if (tt < bdist) {	/* Better solution */
						bdist = tt;
						bsoln = i;
					}
				}
			}
			i = bsoln;
		}
		for (e = 0; e < p->clutTable->di; e++) {
			out[e] = pp[i].p[e];			/* Solution */
		}
	}

	/* Sanitise auxiliary locus range and auxiliary value return */
	if (auxr != NULL || auxv != NULL) {
		int ee = 0;
		for (e = 0; e < p->clutTable->di; e++) {
			double v = out[e];			/* Solution */
			if (p->auxm[e] != 0) {
				if (auxr != NULL) {			/* Make sure locus encloses actual value */
					if (auxr[2 * ee] > v)
						auxr[2 * ee] = v;
					if (auxr[2 * ee + 1] < v)
						auxr[2 * ee + 1] = v;
				}
				if (auxv != NULL) {
					auxv[ee] = v;
				}
				ee++;
			}
		}
	}

#ifdef CHECK_ILIMIT	/* Do sanity checks on meeting ink limit */
if (p->ink.tlimit >= 0.0 || p->ink.klimit >= 0.0) {
	double sum = icxLimitD((void *)p, out);
	if (sum > 0.0)
		printf("xlut assert%s: icxLuLut_inv_clut returned outside limits by %f > tlimit %f\n",crv ? " (clip)" : "", sum, p->ink.tlimit);
}
#endif

	DBR(("inv_clut_aux returning %f %f %f %f\n",out[0],out[1],out[2],out[3]))
	return crv;
}

/* Do output'->input' lookup, simple version */
/* Note than out[] will carry inking value if icxKrule is icxKvalue of icxKlocus */
/* and that the icxKrule value will be in the input' space. */
/* Note that the ink limit will be computed after converting input' to input */
int icxLuLut_inv_clut(icxLuLut *p, double *out, double *in) {
	return icxLuLut_inv_clut_aux(p, out, NULL, NULL, NULL, in);
}

/* Given the proposed auxiliary input' values in in[di], */
/* and the target output' (ie. PCS') values in out[fdi], */
/* return the auxiliary input' values as a proportion of their possible locus */
/* in locus[di]. */
/* This is generally used on a source CMYK profile to convey the black intent */
/* to destination CMYK profile. */
int icxLuLut_clut_aux_locus(icxLuLut *p, double *locus, double *out, double *in) {
	co pp[1];		/* Room for all the solutions found */
	int nsoln;		/* Number of solutions found */
	int e,f;

	if (p->clutTable->di > p->clutTable->fdi) {	/* ie. CMYK->Lab, there will be ambiguity */
		double min[MXDI], max[MXDI];	/* Auxiliary locus range */

		/* Setup for auxiliary locus lookup */
		for (f = 0; f < p->clutTable->fdi; f++) {
			pp[0].v[f] = out[f];			/* Target value */
		}
	
		/* Compute auxiliary locus */
		nsoln = p->clutTable->rev_locus(
			p->clutTable,	/* rspl object */
			p->auxm,		/* Auxiliary mask */
			pp,				/* Input target and output solutions */
			min, max);		/* Returned locus of valid auxiliary values */
		
		if (nsoln == 0) {
			for (e = 0; e < p->clutTable->di; e++)
				locus[e] = 0.0;		/* Return some safe values */
		} else {  /* Got a valid locus */

			/* Figure out the proportion of the locus */
			for (e = 0; e < p->clutTable->di; e++) {
				if (p->auxm[e] != 0) {
					double iv = in[e];
					if (iv <= min[e])
						locus[e] = 0.0;
					else if (iv >= max[e])
						locus[e] = 1.0;
					else {
						double lpl = max[e] - min[e];	/* Locus path length */
						if (lpl > 1e-6)
							locus[e] = (iv - min[e])/lpl;
						else
							locus[e] = 0.0;
					}
				}
			}
		}
	} else {
		/* There should be no auxiliaries */
		for (e = 0; e < p->clutTable->di; e++)
			locus[e] = 0.0;		/* Return some safe values */
	}
	return 0;
}

/* Do input' -> input inverse lookup */
int icxLuLut_inv_input(icxLuLut *p, double *out, double *in) {
#ifdef NEVER
	return ((icmLuLut *)p->plu)->inv_input((icmLuLut *)p->plu, out, in);
#else
	int rv = 0;
	int i,j;
	int nsoln;				/* Number of solutions found */
	co pp[MAX_INVSOLN];		/* Room for all the solutions found */
	double cdir;

	DBR(("inv_input got DEV' %f %f %f %f\n",in[0],in[1],in[2],in[3]))

	for (i = 0; i < p->inputChan; i++) {
		pp[0].p[0] = p->inputClipc[i];
		pp[0].v[0] = in[i];
		cdir = p->inputClipc[i] - in[i];	/* Clip towards output range */

		nsoln = p->inputTable[i]->rev_interp (
			p->inputTable[i], 	/* this */
			RSPL_NEARCLIP,		/* Clip to nearest (faster than vector) */
			MAX_INVSOLN,		/* Maximum number of solutions allowed for */
			NULL, 				/* No auxiliary input targets */
			&cdir,				/* Clip vector direction and length */
			pp);				/* Input and output values */

		if (nsoln & RSPL_DIDCLIP)
			rv = 1;

		nsoln &= RSPL_NOSOLNS;		/* Get number of solutions */

		if (nsoln == 1) { /* Exactly one solution */
			j = 0;
		} else if (nsoln == 0) {	/* Zero solutions. This is unexpected. */
			error("~~~1 Unexpected failure to find reverse solution for input table");
			return 2;
		} else {		/* Multiple solutions */
			/* Use a simple minded resolution - choose the one closest to the center */
			double bdist = 1e300;
			int bsoln = 0;
/* Don't expect this - 1D luts are meant to be monotonic */
printf("~~~1 got %d reverse solutions\n",nsoln);
printf("~~~1 solution 0 = %f\n",pp[0].p[0]);
printf("~~~1 solution 1 = %f\n",pp[1].p[0]);
			for (j = 0; j < nsoln; j++) {
				double tt;
				tt = pp[i].p[0] - p->inputClipc[i];
				tt *= tt;
				if (tt < bdist) {	/* Better solution */
					bdist = tt;
					bsoln = j;
				}
			}
			j = bsoln;
		}
		out[i] = pp[j].p[0];
	}

	DBR(("inv_input returning DEV %f %f %f %f\n",out[0],out[1],out[2],out[3]))
	return rv;
#endif /* NEVER */
}

/* Possible inverse matrix lookup */
int icxLuLut_inv_matrix(icxLuLut *p, double *out, double *in) {
	int rv = 0;
	rv |= ((icmLuLut *)p->plu)->inv_matrix((icmLuLut *)p->plu, out, in);
	return rv;
}

int icxLuLut_inv_in_abs(icxLuLut *p, double *out, double *in) {
	int rv = 0;
	rv |= ((icmLuLut *)p->plu)->inv_in_abs((icmLuLut *)p->plu, out, in);

	if (p->ins == icxSigJabData) {
		p->cam->XYZ_to_cam(p->cam, out, out);
	}

	return rv;
}

/* Overall inverse lookup */
/* Note that if k_rule is icxKvalue, the inking value is in lut input space */
static int
icxLuLut_inv_lookup (
icxLuBase *pp,		/* This */
double *out,		/* Vector of output values/input auxiliary values */
double *in			/* Vector of input values */
) {
	icxLuLut *p = (icxLuLut *)pp;
	int rv = 0;
	int i;
	double temp[MAX_CHAN];
	double tout[MAX_CHAN];			/* Use as out[] to clut in case out == in */
	if (p->mergeclut == 0) {		/* Do this if it's not merger with clut */
		rv |= p->inv_out_abs (p, temp, in);
		rv |= p->inv_output  (p, temp, temp);
	} else {
		for (i = 0; i < p->outputChan; i++)
			temp[i] = in[i];
	}
	if (p->ink.k_rule == icxKvalue) {	/* aux K value is provided in out[] */
		co tc;
		/* Convert aux targets from input to input' space ready for inv_clut */
		for (i = 0; i < p->inputChan; i++) {
			tc.p[0] = out[i];
			p->inputTable[i]->interp(p->inputTable[i], &tc);
			tout[i] = tc.v[0];
		}
	} else if (p->ink.k_rule == icxKlocus) {	/* Carry aux K locus above input */
		/* ~~~888 This isn't right - we need to convert between dev and dev' */
		/* locus values. This means we need the dev locus ranges here, to */
		/* recompute the proportions ~~~888 */
		for (i = p->outputChan; i < p->inputChan; i++)
			tout[i] = out[i];
	}
	rv |= p->inv_clut    (p, tout, temp);
	rv |= p->inv_input   (p, out, tout);
	rv |= p->inv_matrix  (p, out, out);
	rv |= p->inv_in_abs  (p, out, out);
	return rv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Given a nominated output PCS (ie. Absolute, Jab etc.), convert it in the bwd */
/* direction into a relative XYZ or Lab PCS value */
/* (This is used in generating gamut compression in B2A tables) */
void icxLuLut_bwd_outpcs_relpcs(
icxLuBase *pp,
icColorSpaceSignature os,		/* Output space, XYZ or Lab */
double *out, double *in) {
	icxLuLut *p = (icxLuLut *)pp;

	icxLuLut_inv_out_abs(p, out, in);
	if (os == icSigXYZData && p->natpcs == icSigLabData) {
		icmLab2XYZ(&icmD50, out, out);
	} else if (os == icSigXYZData && p->natpcs == icSigLabData) {
		icmXYZ2Lab(&icmD50, out, out);
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Return LuLut information */
static void icxLuLut_get_info(
	icxLuLut     *p,		/* this */
	icmLut       **lutp,	/* Pointer to icc lut type return value */
	icmXYZNumber *pcswhtp,	/* Pointer to profile PCS white point return value */
	icmXYZNumber *whitep,	/* Pointer to profile absolute white point return value */
	icmXYZNumber *blackp	/* Pointer to profile absolute black point return value */
) {
	((icmLuLut *)p->plu)->get_info((icmLuLut *)p->plu, lutp, pcswhtp, whitep, blackp);
}

/* Return the underlying Lut matrix */
static void
icxLuLut_get_matrix (
	icxLuLut *p,
	double m[3][3]
) {
	((icmLuLut *)p->plu)->get_matrix((icmLuLut *)p->plu, m);
}

static void
icxLuLut_free(
icxLuBase *pp
) {
	icxLuLut *p = (icxLuLut *)pp;
	int i;

	for (i = 0; i < p->inputChan; i++) {
		if (p->inputTable[i] != NULL)
			p->inputTable[i]->del(p->inputTable[i]);
		if (p->revinputTable[i] != NULL)
			p->revinputTable[i]->del(p->revinputTable[i]);
	}

	if (p->clutTable != NULL)
		p->clutTable->del(p->clutTable);

	if (p->cclutTable != NULL)
		p->cclutTable->del(p->cclutTable);

	for (i = 0; i < p->outputChan; i++) {
		if (p->outputTable[i] != NULL)
			p->outputTable[i]->del(p->outputTable[i]);
	}

	if (p->plu != NULL)
		p->plu->del(p->plu);

	if (p->cam != NULL)
		p->cam->del(p->cam);

	if (p->absxyzlu != NULL)
		p->absxyzlu->del(p->absxyzlu);

	free(p);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

static gamut *icxLuLutGamut(icxLuBase *plu, double detail); 

/* Do the basic icxLuLut creation and initialisation */
static icxLuLut *
alloc_icxLuLut(
	xicc                  *xicp,
	icmLuBase             *plu,			/* Pointer to Lu we are expanding (ours) */
	int                   flags			/* clip, merge flags */
) {
	icxLuLut *p;						/* Object being created */
	icmLuLut *luluto = (icmLuLut *)plu;	/* Lookup Lut type object */

	if ((p = (icxLuLut *) calloc(1,sizeof(icxLuLut))) == NULL)
		return NULL;

	p->pp                = xicp;
	p->plu               = plu;
	p->del               = icxLuLut_free;
	p->lutspaces         = icxLutSpaces;
	p->spaces            = icxLuSpaces;
	p->get_native_ranges = icxLu_get_native_ranges;
	p->get_ranges        = icxLu_get_ranges;
	p->rel_wh_bk_points  = icxLuRel_wh_bk_points;
	p->get_gamut         = icxLuLutGamut;
	p->fwd_relpcs_outpcs = icxLuLut_fwd_relpcs_outpcs;
	p->bwd_outpcs_relpcs = icxLuLut_bwd_outpcs_relpcs;
	p->nearclip = 0;				/* Set flag defaults */
	p->mergeclut = 0;
	p->noisluts = 0;
	p->noipluts = 0;
	p->nooluts = 0;
	p->intsep = 0;

	p->lookup   = icxLuLut_lookup;
	p->in_abs   = icxLuLut_in_abs;
	p->matrix   = icxLuLut_matrix;
	p->input    = icxLuLut_input;
	p->clut     = icxLuLut_clut;
	p->clut_aux = icxLuLut_clut_aux;
	p->output   = icxLuLut_output;
	p->out_abs  = icxLuLut_out_abs;

	p->inv_lookup   = icxLuLut_inv_lookup;
	p->inv_in_abs   = icxLuLut_inv_in_abs;
	p->inv_matrix   = icxLuLut_inv_matrix;
	p->inv_input    = icxLuLut_inv_input;
	p->inv_clut     = icxLuLut_inv_clut;
	p->inv_clut_aux = icxLuLut_inv_clut_aux;
	p->inv_output   = icxLuLut_inv_output;
	p->inv_out_abs  = icxLuLut_inv_out_abs;

	p->clut_locus   = icxLuLut_clut_aux_locus;

	p->get_info   = icxLuLut_get_info;
	p->get_matrix = icxLuLut_get_matrix;

	/* Setup all the rspl analogs of the icc Lut */
	/* NOTE: We assume that none of this relies on the flag settings, */
	/* since they will be set on our return. */

	/* Get details of underlying, native icc color space */
	p->plu->lutspaces(p->plu, &p->natis, NULL, &p->natos, NULL, &p->natpcs);

	/* Get other details of conversion */
	p->plu->spaces(p->plu, NULL, &p->inputChan, NULL, &p->outputChan, NULL, NULL, NULL, NULL);

	/* Sanity check */
	if (p->inputChan > MXDI) {
		sprintf(p->pp->err,"xicc can only handle input channels of %d or less",MXDI);
		p->inputChan = MXDI;		/* Avoid going outside array bounds */
		p->pp->errc = 1;
		p->del((icxLuBase *)p);
		return NULL;
	}
	if (p->outputChan > MXDO) {
		sprintf(p->pp->err,"xicc can only handle output channels of %d or less",MXDO);
		p->outputChan = MXDO;		/* Avoid going outside array bounds */
		p->pp->errc = 1;
		p->del((icxLuBase *)p);
		return NULL;
	}

	/* Get pointer to icmLut */
	luluto->get_info(luluto, &p->lut, NULL, NULL, NULL);

	return p;
}

/* Initialise the clut ink limiting and black */
/* generation information. */
/* return 0 or error code */
static int
setup_ink_icxLuLut(
icxLuLut *p,			/* Object being initialised */
icxInk   *ink,			/* inking details (NULL for default) */
int       setLminmax	/* Figure the L locus for inking rule */
) {
	if (ink) {
		p->ink = *ink;	/* Copy the structure */
	} else {
		p->ink.tlimit = 3.0;			/* default ink limit of 300% */
		p->ink.klimit = -1.0;			/* default no black limit */
		p->ink.k_rule = icxKluma5;		/* default K generation rule */
		p->ink.c.Ksmth = ICXINKDEFSMTH;	/* Default smoothing */
		p->ink.c.Kstle = 0.0;		/* Min K at white end */
		p->ink.c.Kstpo = 0.0;		/* Start of transition is at white */
		p->ink.c.Kenle = 1.0;		/* Max K at black end */
		p->ink.c.Kenpo = 1.0;		/* End transition at black */
		p->ink.c.Kshap = 1.0;		/* Linear transition */
	}

	/* Normalise total and black ink limits */
    if (p->ink.tlimit <= 1e-4 || p->ink.tlimit >= (double)p->inputChan)
    	p->ink.tlimit = -1.0;		/* Turn ink limit off if not effective */
    if (p->inputChan < 4 || p->ink.klimit < 0.0 || p->ink.klimit >= 1.0)
    	p->ink.klimit = -1.0;		/* Turn black limit off if not effective */
	
	/* Set the ink limit information for any reverse interpolation. */
	/* Calling this will clear the reverse interpolaton cache. */
	p->clutTable->rev_set_limit(
		p->clutTable,		/* this */
		p->ink.tlimit >= 0.0 || p->ink.klimit >= 0.0 ? icxLimitD : NULL,
		                    /* Optional input space limit() function. */
		                	/* Function should evaluate in[0..di-1], and return */
		                	/* numbert that is not to exceed 0.0. NULL if not used. */
		(void *)p,			/* Context passed to limit() */
		0.0					/* Value that limit() is not to exceed */
	);

	/* Duplicate in the CAM clip rspl if it exists */
	if (p->cclutTable != NULL) {
		p->cclutTable->rev_set_limit(
			p->cclutTable,		/* this */
			p->ink.tlimit >= 0.0 || p->ink.klimit >= 0.0 ? icxLimitD : NULL,
			                    /* Optional input space limit() function. */
			                	/* Function should evaluate in[0..di-1], and return */
			                	/* number that is not to exceed 0.0. NULL if not used. */
			(void *)p,			/* Context passed to limit() */
			0.0					/* Value that limit() is not to exceed */
		);
	}

	/* Figure Lmin and Lmax for icxKluma5 curve basis */
	if (setLminmax
	 && p->clutTable->di > p->clutTable->fdi) {	/* If K generation makes sense */
		double wh[3], bk[3];
		int mergeclut;			/* Save/restore mergeclut value */

		/* Get white/black in effective xlu pcs space */
		p->rel_wh_bk_points((icxLuBase *)p, wh, bk);

		/* Convert from effective PCS to native relative XYZ or Lab PCS */
		mergeclut = p->mergeclut;			/* Hack to be able to use inv_out_abs() */
		p->mergeclut = 0;					/* if mergeclut is active. */
		icxLuLut_inv_out_abs(p, wh, wh);
		icxLuLut_inv_out_abs(p, bk, bk);
		p->mergeclut = mergeclut;			/* Restore */

		/* Convert to Lab PCS */
		if (p->natos == icSigXYZData) {	/* Always do K rule in L space */
			icmXYZ2Lab(&icmD50, wh, wh);
			icmXYZ2Lab(&icmD50, bk, bk);
		}
		p->Lmax = 0.01 * wh[0];
		p->Lmin = 0.01 * bk[0];
	} else {

		/* Some sane defaults */
		p->Lmax = 1.0;
		p->Lmin = 0.0;
	}

	return 0;
}

	
/* Initialise the clut clipping information, ink limiting */
/* and auxiliary parameter settings for all the rspl. */
/* return 0 or error code */
static int
setup_clip_icxLuLut(
icxLuLut *p			/* Object being initialised */
) {
	double tmin[MXDIDO], tmax[MXDIDO]; 
	int i;

	/* Setup for inversion of multi-dim clut */

	/* Set auxiliaries */
	for (i = 0; i < p->inputChan; i++)
		p->auxm[i] = 0;

	if (p->inputChan > p->outputChan) {
		switch(p->natis) {
			case icSigCmykData:
				p->auxm[3] = 1;		/* K is the auxiliary channel */
				break;
			default:
				p->pp->errc = 2;
				sprintf(p->pp->err,"Unknown colorspace %s when setting auxliaries",
				                icm2str(icmColorSpaceSignature, p->natis));
				return p->pp->errc;
				break;
		}
	}

	p->auxlinf = NULL;		/* Input space auxiliary linearisation function  - not implemented */
	p->auxlinf_ctx = NULL;	/* Opaque context for auxliliary linearisation function */

	/* Aproximate center of clut input gamut - used for */
	/* resolving multiple reverse solutions. */
	p->clutTable->get_in_range(p->clutTable, tmin, tmax);
	for (i = 0; i < p->clutTable->di; i++) {
		p->icent[i] = (tmin[i] + tmax[i])/2.0;
	}

	/* Compute clip setup information relating to clut output gamut. */
	if (p->nearclip != 0			/* Near clip requested */
	 || p->inputChan == 1) {		/* or vector clip won't work */
		p->clip.nearclip = 1;

	} else {		/* Vector clip */
		icColorSpaceSignature clutos = p->natos;

		p->clip.nearclip = 0;
		p->clip.LabLike = 0;
		p->clip.fdi = p->clutTable->fdi;

		switch(clutos) {
			case icxSigJabData:
			case icSigLabData: {
				co pp;				/* Room for all the solutions found */
				int nsoln;			/* Number of solutions found */
				double cdir[MXDO];	/* Clip vector direction and length */

				p->clip.LabLike = 1;

				/* Find high clip target */
				for (i = 0; i < p->inputChan; i++)
					pp.p[i] = 0.0;					/* Set aux values */
				pp.v[0] = 105.0; pp.v[1] = pp.v[2] = 0.0; 	/* PCS Target value */
				cdir[0] = cdir[1] = cdir[2] = 0.0;	/* Clip Target */

				p->inv_output(p, pp.v, pp.v);				/* Compensate for output curve */
				p->inv_output(p, cdir, cdir);
			
				cdir[0] -= pp.v[0];							/* Clip vector */
				cdir[1] -= pp.v[1];
				cdir[2] -= pp.v[2];

				/* PCS -> Device with clipping */
				nsoln = p->clutTable->rev_interp(
					p->clutTable, 	/* rspl object */
					0,				/* No flags - might be in gamut, might vector clip */
					1,			 	/* Maxumum solutions to return */
					p->auxm, 		/* Auxiliary input targets */
					cdir,			/* Clip vector direction and length */
					&pp);			/* Input target and output solutions */
									/* returned solutions in pp[0..retval-1].p[] */
				nsoln &= RSPL_NOSOLNS;	/* Get number of solutions */

				if (nsoln != 1) {
					p->pp->errc = 2;
					sprintf(p->pp->err,"Failed to find high clip target for Lab space");
					return p->pp->errc;
				}

				p->clip.ocent[0] = pp.v[0] - 0.001;					/* Got main target */
				p->clip.ocent[1] = pp.v[1];
				p->clip.ocent[2] = pp.v[2];

				/* Find low clip target */
				pp.v[0] = -5.0; pp.v[1] = pp.v[2] = 0.0; 	/* PCS Target value */
				cdir[0] = 100.0; cdir[1] = cdir[2] = 0.0;	/* Clip Target */

				p->inv_output(p, pp.v, pp.v);				/* Compensate for output curve */
				p->inv_output(p, cdir, cdir);
				cdir[0] -= pp.v[0];							/* Clip vector */
				cdir[1] -= pp.v[1];
				cdir[2] -= pp.v[2];

				/* PCS -> Device with clipping */
				nsoln = p->clutTable->rev_interp(
					p->clutTable, 	/* rspl object */
					RSPL_WILLCLIP,	/* Since there was no locus, we expect to have to clip */
					1,			 	/* Maxumum solutions to return */
					NULL, 			/* No auxiliary input targets */
					cdir,			/* Clip vector direction and length */
					&pp);			/* Input target and output solutions */
									/* returned solutions in pp[0..retval-1].p[] */
				nsoln &= RSPL_NOSOLNS;		/* Get number of solutions */
				if (nsoln != 1) {
					p->pp->errc = 2;
					sprintf(p->pp->err,"Failed to find low clip target for Lab space");
					return p->pp->errc;
				}

				p->clip.ocentv[0] = pp.v[0] + 0.001 - p->clip.ocent[0];	/* Raw line vector */
				p->clip.ocentv[1] = pp.v[1] - p->clip.ocent[1];
				p->clip.ocentv[2] = pp.v[2] - p->clip.ocent[2];

				/* Compute vectors length */
				for (p->clip.ocentl = 0.0, i = 0; i < 3; i++)
					p->clip.ocentl += p->clip.ocentv[i] * p->clip.ocentv[i];
				p->clip.ocentl = sqrt(p->clip.ocentl);
				if (p->clip.ocentl <= 1e-8)
					p->clip.ocentl = 0.0;

				break;
				}
			case icSigXYZData:
				// ~~~~~~1 need to add this.

			default:
				/* Do a crude approximation, that may not work. */
				p->clutTable->get_out_range(p->clutTable, tmin, tmax);
				for (i = 0; i < p->clutTable->fdi; i++) {
					p->clip.ocent[i] = (tmin[i] + tmax[i])/2.0;
				}
				p->clip.ocentl = 0.0;
				break;
		}
	}
	return 0;
}

/* Function to pass to rspl to set input/output transfer functions */
static void
icxLuLut_inout_func(
	void *pp,			/* icxLuLut */
	double *out,		/* output value */
	double *in			/* inut value */
) {
	icxLuLut *p      = (icxLuLut *)pp;			/* this */
	icmLuLut *luluto = (icmLuLut *)p->plu;		/* Get icmLuLut object */
	double tin[MAX_CHAN];
	double tout[MAX_CHAN];
	int i;

	if (p->iol_out == 0) {			/* fwd input */
#ifdef INK_LIMIT_TEST
		tout[p->iol_ch] = in[0];
#else
		for (i = 0; i < p->inputChan; i++)
			tin[i] = 0.0;
		tin[p->iol_ch] = in[0];
		luluto->input(luluto, tout, tin);
#endif
	} else if (p->iol_out == 1) {	/* fwd output */
		for (i = 0; i < p->outputChan; i++)
			tin[i] = 0.0;
		tin[p->iol_ch] = in[0];
		luluto->output(luluto, tout, tin);
	} else {						/* bwd input */
#ifdef INK_LIMIT_TEST
		tout[p->iol_ch] = in[0];
#else
		for (i = 0; i < p->inputChan; i++)
			tin[i] = 0.0;
		tin[p->iol_ch] = in[0];
		/* This won't be valid if matrix is used !!! */
		luluto->inv_input(luluto, tout, tin);
		luluto->inv_in_abs(luluto, tout, tout);
#endif
	}
	out[0] = tout[p->iol_ch];
}

/* Function to pass to rspl to set clut up, when mergeclut is set */
static void
icxLuLut_clut_merge_func(
	void *pp,			/* icxLuLut */
	double *out,		/* output value */
	double *in			/* input value */
) {
	icxLuLut *p      = (icxLuLut *)pp;			/* this */
	icmLuLut *luluto = (icmLuLut *)p->plu;		/* Get icmLuLut object */

	luluto->clut(luluto, out, in);
	luluto->output(luluto, out, out);
	luluto->out_abs(luluto, out, out);

	if (p->outs == icxSigJabData)
		p->cam->XYZ_to_cam(p->cam, out, out);
}

/* Implimenation of Lut create from icc. */
/* Note that xicc_get_luobj() will have set the pcsor & */
/* intent to consistent values if Jab and/or icxAppearance */
/* has been requested. */
/* It will also have created the underlying icm lookup object */
/* that is used to create and implement the icx one. The icm */
/* will be used to translate from native to effective PCS, unless */
/* the effective PCS is Jab, in which case the icm will be set to */
/* have an effective PCS of XYZ. Since native<->effecive PCS conversion */
/* is done at the to/from_abs() stage, none of it affects the individual */
/* conversion steps, which will all talk the native PCS (unless merged). */
static icxLuBase *
new_icxLuLut(
xicc                  *xicp,
int                   flags,		/* clip, merge flags */
icmLuBase             *plu,			/* Pointer to Lu we are expanding (ours) */
icmLookupFunc         func,			/* Functionality requested */
icRenderingIntent     intent,		/* Rendering intent */
icColorSpaceSignature pcsor,		/* PCS override (0 = def) */
icxViewCond           *vc,			/* Viewing Condition (NULL if not using CAM) */
icxInk                *ink			/* inking details (NULL for default) */
) {
	icxLuLut *p;						/* Object being created */
	icmLuLut *luluto = (icmLuLut *)plu;	/* Lookup Lut type object */

	int i;

	/* Do basic creation and initialisation */
	if ((p = alloc_icxLuLut(xicp, plu, flags)) == NULL)
		return NULL;

	/* Set LuLut "use" specific creation flags: */
	if (flags & ICX_CLIP_NEAREST) {
		p->nearclip = 1;
	}

	if (flags & ICX_MERGE_CLUT)
		p->mergeclut = 1;

	/* We're only implementing this under specific conditions. */
	if (flags & ICX_CAM_CLIP
	 && func == icmFwd
	 && !(p->mergeclut != 0 && pcsor == icxSigJabData))		/* Don't need camclip if merged Jab */
		p->camclip = 1;

	if (flags & ICX_INT_SEPARATE) {
fprintf(stderr,"~1 Internal optimised 4D separations not yet implemented!\n");
		p->intsep = 1;
	}

	/* Init the CAM model if it will be used */
	if (pcsor == icxSigJabData || p->camclip) {
		if (vc != NULL)		/* One has been provided */
			p->vc  = *vc;		/* Copy the structure */
		else
			xicc_enum_viewcond(xicp, &p->vc, -1, NULL, 0);	/* Use a default */
		p->cam = new_icxcam(cam_default);
		p->cam->set_view(p->cam, p->vc.Ev, p->vc.Wxyz, p->vc.Yb, p->vc.La, p->vc.Lv, p->vc.Yf, p->vc.Fxyz, XICC_USE_HK);
	} else {
		p->cam = NULL;
	}
	
	/* Remember the effective intent */
	p->intent = intent;

	/* Get the effective spaces of underlying icm */
	plu->spaces(plu, &p->ins, NULL, &p->outs, NULL, NULL, NULL, NULL, &p->pcs);

	/* Override with pcsor */
	/* We assume that any profile that has a cie color as a "device" color */
	/* intends it to stay that way, and not be overridden. */
	if (pcsor == icxSigJabData) {
		p->pcs = pcsor;		

		if (xicp->pp->header->deviceClass == icSigAbstractClass) {
				p->ins = pcsor;
				p->outs = pcsor;

		} else if (xicp->pp->header->deviceClass != icSigLinkClass) {
			if (func == icmBwd || func == icmGamut || func == icmPreview)
				p->ins = pcsor;
			if (func == icmFwd || func == icmPreview)
				p->outs = pcsor;
		}
	}

	/* In general the native and effective ranges of the icx will be the same as the */
	/* underlying icm lookup object. */
	p->plu->get_lutranges(p->plu, p->ninmin, p->ninmax, p->noutmin, p->noutmax);
	p->plu->get_ranges(p->plu, p->inmin,  p->inmax,  p->outmin,  p->outmax);

	/* If we have a Jab PCS override, reflect this in the effective icx range. */
	/* Note that the ab ranges are nominal. They will exceed this range */
	/* for colors representable in L*a*b* PCS */
	if (p->ins == icxSigJabData) {
		p->inmin[0] = 0.0;		p->inmax[0] = 100.0;
		p->inmin[1] = -128.0;	p->inmax[1] = 128.0;
		p->inmin[2] = -128.0;	p->inmax[2] = 128.0;
	} else if (p->outs == icxSigJabData) {
		p->outmin[0] = 0.0;		p->outmax[0] = 100.0;
		p->outmin[1] = -128.0;	p->outmax[1] = 128.0;
		p->outmin[2] = -128.0;	p->outmax[2] = 128.0;
	} 

	/* If we have a merged clut, reflect this in the icx native PCS range. */
	/* Merging merges output processing (irrespective of whether we are using */
	/* the forward or backward cluts) */
	if (p->mergeclut != 0) {
		int i;
		for (i = 0; i < p->outputChan; i++) {
			p->noutmin[i] = p->outmin[i];
			p->noutmax[i] = p->outmax[i];
		}
	}

	/* ------------------------------- */
	/* Create rspl based input lookups */
	for (i = 0; i < p->inputChan; i++) {
		if ((p->inputTable[i] = new_rspl(1, 1)) == NULL) {
			p->pp->errc = 2;
			sprintf(p->pp->err,"Creation of input table rspl failed");
			p->del((icxLuBase *)p);
			return NULL;
		}
		p->iol_out = 0;		/* Input lookup */
		p->iol_ch = i;		/* Chanel */
		p->inputTable[i]->set_rspl(p->inputTable[i], 0,
		           (void *)p, icxLuLut_inout_func,
		           &p->ninmin[i], &p->ninmax[i], (int *)&p->lut->inputEnt, &p->ninmin[i], &p->ninmax[i]);
	}

	/* Setup center clip target for inversion */
	for (i = 0; i < p->inputChan; i++) {
		p->inputClipc[i] = (p->ninmin[i] + p->ninmax[i])/2.0;
	}

	/* Create rspl based reverse input lookups used in ink limit function. */
	for (i = 0; i < p->inputChan; i++) {
		int gres = 256;
		if ((p->revinputTable[i] = new_rspl(1, 1)) == NULL) {
			p->pp->errc = 2;
			sprintf(p->pp->err,"Creation of reverse input table rspl failed");
			p->del((icxLuBase *)p);
			return NULL;
		}
		p->iol_out = 2;		/* Input lookup */
		p->iol_ch = i;		/* Chanel */
		p->revinputTable[i]->set_rspl(p->revinputTable[i], 0,
		           (void *)p, icxLuLut_inout_func,
		           &p->ninmin[i], &p->ninmax[i], &gres, &p->ninmin[i], &p->ninmax[i]);
	}

	/* ------------------------------- */
	{ 
		int gres[MXDI];

		for (i = 0; i < p->inputChan; i++)
			gres[i] = p->lut->clutPoints;

		/* Create rspl based multi-dim table */
		if ((p->clutTable = new_rspl(p->inputChan, p->outputChan)) == NULL) {
			p->pp->errc = 2;
			sprintf(p->pp->err,"Creation of clut table rspl failed");
			p->del((icxLuBase *)p);
			return NULL;
		}

		if (p->mergeclut == 0) {	/* Do this if it's not merged with clut, */
			p->clutTable->set_rspl(p->clutTable, 0, (void *)luluto,
			           (void (*)(void *, double *, double *))luluto->clut,
		               p->ninmin, p->ninmax, gres, p->noutmin, p->noutmax);

		} else {	/* If mergeclut */
			p->clutTable->set_rspl(p->clutTable, 0, (void *)p,
			           icxLuLut_clut_merge_func,
		               p->ninmin, p->ninmax, gres, p->noutmin, p->noutmax);

		}

		/* clut clipping is setup separately */
	}

	/* ------------------------------- */
	/* Create rspl based output lookups */
	for (i = 0; i < p->outputChan; i++) {
		if ((p->outputTable[i] = new_rspl(1, 1)) == NULL) {
			p->pp->errc = 2;
			sprintf(p->pp->err,"Creation of output table rspl failed");
			p->del((icxLuBase *)p);
			return NULL;
		}
		p->iol_out = 1;		/* Output lookup */
		p->iol_ch = i;		/* Chanel */
		p->outputTable[i]->set_rspl(p->outputTable[i], 0,
		           (void *)p, icxLuLut_inout_func,
		           &p->noutmin[i], &p->noutmax[i], (int *)&p->lut->outputEnt, &p->noutmin[i], &p->noutmax[i]);
	}

	/* Setup center clip target for inversion */
	for (i = 0; i < p->outputChan; i++) {
		p->outputClipc[i] = (p->noutmin[i] + p->noutmax[i])/2.0;
	}

	/* ------------------------------- */

	/* Setup all the clipping, ink limiting and auxiliary stuff, */
	/* in case a reverse call is used. Only do this if we know */
	/* the reverse stuff isn't going to fail due to channel limits. */
	if (p->clutTable->within_restrictedsize(p->clutTable)) {

		if (setup_ink_icxLuLut(p, ink, 1) != 0) {
			p->del((icxLuBase *)p);
			return NULL;
		}
	
		if (setup_clip_icxLuLut(p) != 0) {
			p->del((icxLuBase *)p);
			return NULL;
		}
	}

	return (icxLuBase *)p;
}


/* Function to pass to rspl to set clut up, when camclip is going to be used. */
/* We use the temporary icm fwd absolute xyz lookup, then convert to CAM Jab. */
static void
icxLuLut_clut_camclip_func(
	void *pp,			/* icxLuLut */
	double *out,		/* output value */
	double *in			/* inut value */
) {
	icxLuLut *p      = (icxLuLut *)pp;			/* this */
	icmLuLut *luluto = (icmLuLut *)p->absxyzlu;

	luluto->clut(luluto, out, in);
	luluto->output(luluto, out, out);
	luluto->out_abs(luluto, out, out);
	p->cam->XYZ_to_cam(p->cam, out, out);
}

/* Initialise the additional CAM space clut rspl, used to compute */
/* reverse lookup CAM clipping results when the camclip flag is set. */
/* return error code. */
static int
icxLuLut_init_clut_camclip(
icxLuLut *p) {
	int e, gres[MXDI];

	/* Setup so clut contains transform to CAM Jab */
	/* (camclip is only used in fwd or invfwd direction lookup) */
	double cmin[3], cmax[3];
	cmin[0] = 0.0;		cmax[0] = 100.0;	/* Nominal Jab output ranges */
	cmin[1] = -128.0;	cmax[1] = 128.0;
	cmin[2] = -128.0;	cmax[2] = 128.0;

	/* Get icm lookup we need for seting up and using CAM icx clut */
	if ((p->absxyzlu = p->pp->pp->get_luobj(p->pp->pp, icmFwd, icAbsoluteColorimetric,
	                                    icSigXYZData, icmLuOrdNorm)) == NULL) {
		p->pp->errc = p->pp->pp->errc;		/* Copy error to xicc */
		strcpy(p->pp->err, p->pp->pp->err);
		return p->pp->errc;
	}

	/* Create CAM rspl based multi-dim table */
	if ((p->cclutTable = new_rspl(p->inputChan, p->outputChan)) == NULL) {
		p->pp->errc = 2;
		sprintf(p->pp->err,"Creation of clut table rspl failed");
		return p->pp->errc;
	}

	for (e = 0; e < p->inputChan; e++)
		gres[e] = p->lut->clutPoints;

	/* Setup our special CAM space rspl */
	p->cclutTable->set_rspl(p->cclutTable, 0, (void *)p,
	           icxLuLut_clut_camclip_func,
               p->ninmin, p->ninmax, gres, cmin, cmax);

	/* Duplicate the ink limit information for any reverse interpolation. */
	p->cclutTable->rev_set_limit(
		p->cclutTable,		/* this */
		p->ink.tlimit >= 0.0 || p->ink.klimit >= 0.0 ? icxLimitD : NULL,
		                    /* Optional input space limit() function. */
		                	/* Function should evaluate in[0..di-1], and return */
		                	/* number that is not to exceed 0.0. NULL if not used. */
		(void *)p,			/* Context passed to limit() */
		0.0					/* Value that limit() is not to exceed */
	);
	return 0;
}

/* ========================================================== */
/* xicc creation code                                         */
/* ========================================================== */

/* Callback for computing delta E squared for two output (PCS) */
/* values, passed as a callback to xfit */
static double xfit_to_de2(void *cntx, double *in1, double *in2) {
	icxLuLut *p = (icxLuLut *)cntx;
	double rv;

	if (p->pcs == icSigLabData) {
#ifdef USE_CIE94_DE
		rv = icmCIE94sq(in1, in2);
#else
		rv = icmLabDEsq(in1, in2);
#endif
	} else {
		double lab1[3], lab2[3];
		icmXYZ2Lab(&icmD50, lab1, in1);
		icmXYZ2Lab(&icmD50, lab2, in2);
#ifdef USE_CIE94_DE
		rv = icmCIE94sq(lab1, lab2);
#else
		rv = icmLabDEsq(lab1, lab2);
#endif
	}
	return rv;
}

/* Same as above plus partial derivatives */
static double xfit_to_dde2(void *cntx, double dout[2][MXDIDO], double *in1, double *in2) {
	icxLuLut *p = (icxLuLut *)cntx;
	double rv;

	if (p->pcs == icSigLabData) {
		int i,j,k;
		double tdout[2][3];
#ifdef USE_CIE94_DE
		rv = icxdCIE94sq(tdout, in1, in2);
#else
		rv = icxdLabDEsq(tdout, in1, in2);
#endif
		for (k = 0; k < 2; k++) {
			for (j = 0; j < 3; j++)
				dout[k][j] = tdout[k][j];
		}
	} else {
		double lab1[3], lab2[3];
		double dout12[2][3][3];
		double tdout[2][3];
		int i,j,k;

		icxdXYZ2Lab(&icmD50, lab1, dout12[0], in1);
		icxdXYZ2Lab(&icmD50, lab2, dout12[1], in2);
#ifdef USE_CIE94_DE
		rv = icxdCIE94sq(tdout, lab1, lab2);
#else
		rv = icxdLabDEsq(tdout, lab1, lab2);
#endif
		/* Compute partial derivative (is this correct ??) */
		for (k = 0; k < 2; k++) {
			for (j = 0; j < 3; j++) {
				dout[k][j] = 0.0;
				for (i = 0; i < 3; i++) {
					dout[k][j] += tdout[k][i] * dout12[k][i][j];
				}
			}
		}
	}
	return rv;
}

#ifdef NEVER
/* Check partial derivative function within xfit_to_dde2() */

static double _xfit_to_dde2(void *cntx, double dout[2][MXDIDO], double *in1, double *in2) {
	icxLuLut *pp = (icxLuLut *)cntx;
	int k, i;
	double rv, drv;
	double trv;
	
	rv = xfit_to_de2(cntx, in1, in2);
	drv = xfit_to_dde2(cntx, dout, in1, in2);

	if (fabs(rv - drv) > 1e-6)
		printf("######## DDE2: RV MISMATCH is %f should be %f ########\n",rv,drv);

	/* Check each parameter delta */
	for (k = 0; k < 2; k++) {
		for (i = 0; i < 3; i++) {
			double *in;
			double del;
	
			if (k == 0)
				in = in1;
			else
				in = in2;

			in[i] += 1e-9;
			trv = xfit_to_de2(cntx, in1, in2);
			in[i] -= 1e-9;
			
			/* Check that del is correct */
			del = (trv - rv)/1e-9;
			if (fabs(dout[k][i] - del) > 0.04) {
				printf("######## DDE2: EXCESSIVE at in[%d][%d] is %f should be %f ########\n",k,i,dout[k][i],del);
			}
		}
	}
	return rv;
}

#define xfit_to_dde2 _xfit_to_dde2

#endif

/* Context for rspl setting input and output curves */
typedef struct {
	int iix;
	int oix;
	xfit *xf;		/* Optimisation structure */
} curvectx;

/* Function to pass to rspl to set input and output */
/* transfer function for xicc lookup function */
static void
set_linfunc(
	void *cc,			/* curvectx structure */
	double *out,		/* Device output value */
	double *in			/* Device input value */
) {
	curvectx *c = (curvectx *)cc;		/* this */
	xfit *p = c->xf;

	if (c->iix >= 0) {				/* Input curve */
		*out = p->incurve(p, *in, c->iix);
	} else if (c->oix >= 0) {		/* Output curve */
		*out = p->outcurve(p, *in, c->oix);
	}
}

/* Function to pass to rspl to set inverse input transfer function, */
/* used for ink limiting calculation. */
static void
icxLuLut_invinput_func(
	void *cc,			/* curvectx structure */
	double *out,		/* Device output value */
	double *in			/* Device input value */
) {
	curvectx *c = (curvectx *)cc;		/* this */
	xfit *p = c->xf;

	*out = p->invincurve(p, *in, c->iix);
}


/* Functions to pass to icc settables() to setup icc A2B Lut: */

/* Input table */
static void set_input(void *cntx, double *out, double *in) {
	icxLuLut *p = (icxLuLut *)cntx;

	if (p->noisluts != 0 && p->noipluts != 0) {	/* Input table must be linear */
		int i;
		for (i = 0; i < p->inputChan; i++)
			out[i] = in[i];
	} else {
		if (p->input(p, out, in) > 1)
			error ("%d, %s",p->pp->errc,p->pp->err);
	}
}

/* clut */
static void set_clut(void *cntx, double *out, double *in) {
	icxLuLut *p = (icxLuLut *)cntx;
	int didc = 0;		/* Clipped */
	int f;
#ifdef DEBUG
	double ucout[MXDO];
#endif

	if (p->clut(p, out, in) > 1)
		error ("%d, %s",p->pp->errc,p->pp->err);

	/* Convert from efective pcs to natural pcs */
	if (p->pcs != p->plu->icp->header->pcs) {
		if (p->pcs == icSigLabData)
			icmLab2XYZ(&icmD50, out, out);
		else
			icmXYZ2Lab(&icmD50, out, out);
	}

#ifdef DEBUG
	ucout[0] = out[0];
	ucout[1] = out[1];
	ucout[2] = out[2];
#endif

	/* Do some crude clipping of the PCS */
	/* Should yell and scream if this happens ??? */
	for (f = 0; f < 3; f++) {
		if (out[f] < p->noutmin[f]) {
			out[f] = p->noutmin[f];
			didc = 1;
		} else if (out[f] > p->noutmax[f]) {
			out[f] = p->noutmax[f];
			didc = 1;
		}
	}
#ifdef DEBUG
	if (didc) {
		printf("set_clut callback clipped %f %f %f -> %f %f %f\n",
		ucout[0], ucout[1], ucout[2], out[0], out[1], out[2]);
	}
#endif
}

/* output */
static void set_output(void *cntx, double *out, double *in) {
	icxLuLut *p = (icxLuLut *)cntx;

	if (p->nooluts != 0) {	/* Output table must be linear */
		int i;
		for (i = 0; i < p->outputChan; i++)
			out[i] = in[i];
	} else {
		if (p->output(p, out, in) > 1)
			error ("%d, %s",p->pp->errc,p->pp->err);
	}
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Routine to figure out a suitable black point for CMYK */

/* Structure to hold optimisation information */
typedef struct {
	icxLuLut *p;			/* Object being created */
	double toAbs[3][3];		/* To abs from aprox relative */
	double p1[3];			/* white pivot point in abs Lab */
	double p2[3];			/* Point on vector towards black */
} bfinds;

/* Optimise device values to minimise L, while remaining */
/* within the ink limit, and staying in line between p1 (white) and p2 (K) */
static double bfindfunc(void *adata, double pv[]) {
	bfinds *b = (bfinds *)adata;
	double rv = 0.0;
	double tt[3], Lab[3];
	co bcc;
	double lr, ta, tb, terr;	/* L ratio, target a, target b, target error */
	double ovr;

	/* See if over ink limit or outside device range */
	ovr = icxLimit((void *)b->p, pv);		/* > 0.0 if outside device gamut or ink limit */

	/* Compute the absolute Lab value: */
	b->p->input(b->p, bcc.p, pv);						/* Through input tables */
	b->p->clutTable->interp(b->p->clutTable, &bcc);		/* Through clut */
	b->p->output(b->p, bcc.v, bcc.v);					/* Through the output tables */

	if (b->p->pcs != icSigXYZData) 	/* Convert PCS to XYZ */
		icmLab2XYZ(&icmD50, bcc.v, bcc.v);

	/* Convert from relative to Absolute colorimetric */
	icmMulBy3x3(tt, b->toAbs, bcc.v);
	icmXYZ2Lab(&icmD50, Lab, tt);	/* Convert to Lab */

#ifdef DEBUG
printf("~1 p1 =  %f %f %f, p2 = %f %f %f\n",b->p1[0],b->p1[1],b->p1[2],b->p2[0],b->p2[1],b->p2[2]);
printf("~1 device value %f %f %f %f, Lab = %f %f %f\n",pv[0],pv[1],pv[2],pv[3],Lab[0],Lab[1],Lab[2]);
#endif

	/* Primary aim is to minimise L value */
	rv = Lab[0];

	/* See how out of line from p1 to p2 we are */
	lr = (Lab[0] - b->p1[0])/(b->p2[0] - b->p1[0]);		/* Distance towards p2 from p1 */
	ta = lr * (b->p2[1] - b->p1[1]) + b->p1[1];			/* Target a value */
	tb = lr * (b->p2[2] - b->p1[2]) + b->p1[2];			/* Target b value */

	terr = (ta - Lab[1]) * (ta - Lab[1])
	     + (tb - Lab[2]) * (tb - Lab[2]);
	
#ifdef DEBUG
printf("~1 target error %f\n",terr);
#endif
	rv += 100.0 * terr;

#ifdef DEBUG
printf("~1 out of range error %f\n",ovr);
#endif
	rv += 200 * ovr;

#ifdef DEBUG
printf("~1 black find tc ret %f\n",rv);
#endif
	return rv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Create icxLuLut and underlying fwd Lut from scattered data */
/* The scattered data is assumed to map Device -> native PCS */
/* NOTE:- in theory once this icxLuLut is setup, it can be */
/* called to translate color values. In practice I suspect */
/* that the icxLuLut hasn't been setup completely enough to allows this. */
/* Might be easier to close it and re-open it ? */
static icxLuBase *
set_icxLuLut(
xicc               *xicp,
icmLuBase          *plu,			/* Pointer to Lu we are expanding (ours) */	
icmLookupFunc      func,			/* Functionality requested */
icRenderingIntent  intent,			/* Rendering intent */
int                flags,			/* white/black point flags */
int                nodp,			/* Number of points */
cow                *ipoints,		/* Array of input points (Lab or XYZ normalized to 1.0) */
double             smooth,			/* RSPL smoothing factor, -ve if raw */
double             avgdev,			/* reading Average Deviation as a prop. of the input range */
icxViewCond        *vc,				/* Viewing Condition (NULL if not using CAM) */
icxInk             *ink,			/* inking details (NULL for default) */
int                quality			/* Quality metric, 0..3 */
) {
	icxLuLut *p;						/* Object being created */
	icc *icco = xicp->pp;				/* Underlying icc object */
	int luflags = 0;					/* icxLuLut alloc clip, merge flags */
	int pcsy;							/* Effective PCS L or Y chanel index */
	double pcsymax;						/* Effective PCS L or Y maximum value */
	icmHeader *h = icco->header;		/* Pointer to icc header */
	int maxchan;						/* max(inputChan, outputChan) */
	int rsplflags = 0;					/* Flags for scattered data rspl */
	int e, f, i, j;
	double dwhite[MXDI], dblack[MXDI];	/* Device white and black values */
	double wp[3];			/* Absolute White point in XYZ */
	double bp[3];			/* Absolute Black point in XYZ */
	double oavgdev[MXDO];	/* Average output value deviation */
	int gres[MXDI];			/* RSPL/CLUT resolution */
	xfit *xf = NULL;		/* Curve fitting class instance */

	if (flags & ICX_VERBOSE)
		rsplflags |= RSPL_VERBOSE;

	if (flags & ICX_EXTRA_FIT)
		rsplflags |= RSPL_EXTRAFIT;		/* Try extra hard to fit data */

	luflags = flags;		/* Transfer straight though ? */

	/* Do basic creation and initialisation */
	if ((p = alloc_icxLuLut(xicp, plu, luflags)) == NULL)
		return NULL;

	/* Set LuLut "create" specific flags: */
	if (flags & ICX_NO_IN_SHP_LUTS)
		p->noisluts = 1;

	if (flags & ICX_NO_IN_POS_LUTS)
		p->noipluts = 1;

	if (flags & ICX_NO_OUT_LUTS)
		p->nooluts = 1;

	/* Get the effective spaces of underlying icm, and set icx the same */
	plu->spaces(plu, &p->ins, NULL, &p->outs, NULL, NULL, &p->intent, NULL, &p->pcs);

	/* For set_icx the effective pcs has to be the same as the native pcs */

	if (p->pcs == icSigXYZData) {
		pcsy = 1;	/* Y of XYZ */
		pcsymax = 1.0;
	} else {
		pcsy = 0;	/* L or Lab */
		pcsymax = 100.0;
	}

	maxchan = p->inputChan > p->outputChan ? p->inputChan : p->outputChan;

	/* Translate overall average deviation into output channel deviation */
	/* (This is for rspl scattered data fitting smoothness adjustment) */
	/* (This could do with more tuning) */
	if (p->pcs == icSigXYZData) {
		oavgdev[0] = 0.60 * avgdev;
		oavgdev[1] = 1.00 * avgdev;
		oavgdev[2] = 0.60 * avgdev;
	} else if (p->pcs == icSigLabData) {
		oavgdev[0] = 1.00 * avgdev;
		oavgdev[1] = 0.60 * avgdev;
		oavgdev[2] = 0.60 * avgdev;
	} else
	{
		for (f = 0; f < p->outputChan; f++)
			oavgdev[f] = avgdev;
	}

	/* In general the native and effective ranges of the icx will be the same as the */
	/* underlying icm lookup object. */
	p->plu->get_lutranges(p->plu, p->ninmin, p->ninmax, p->noutmin, p->noutmax);
	p->plu->get_ranges(p->plu, p->inmin,  p->inmax,  p->outmin,  p->outmax);

	/* ??? Does this ever happen with set_icxLuLut() ??? */
	/* If we have a Jab PCS override, reflect this in the effective icx range. */
	/* Note that the ab ranges are nominal. They will exceed this range */
	/* for colors representable in L*a*b* PCS */
	if (p->ins == icxSigJabData) {
		p->inmin[0] = 0.0;		p->inmax[0] = 100.0;
		p->inmin[1] = -128.0;	p->inmax[1] = 128.0;
		p->inmin[2] = -128.0;	p->inmax[2] = 128.0;
	} else if (p->outs == icxSigJabData) {
		p->outmin[0] = 0.0;		p->outmax[0] = 100.0;
		p->outmin[1] = -128.0;	p->outmax[1] = 128.0;
		p->outmin[2] = -128.0;	p->outmax[2] = 128.0;
	} 

	/* ------------------------------- */

	if (flags & ICX_VERBOSE)
		printf("Estimating white point\n");

	icmXYZ2Ary(wp, icmD50);		/* Set a default value - D50 */
	icmXYZ2Ary(bp, icmBlack);	/* Set a default value - absolute black */

	if (flags & (ICX_SET_WHITE | ICX_SET_BLACK)) {

		/* Figure out as best we can the device white and black points */

		if (h->deviceClass == icSigInputClass) {
			/* We assume that the input target is well behaved, */
			/* and that it includes a white and black point patch, */
			/* and that they have the extreme L/Y values */

			wp[pcsy] = -1e60;
			bp[pcsy] =  1e60;

			/* Discover the white and black patches */
			for (i = 0; i < nodp; i++) {
				if (ipoints[i].v[pcsy] > wp[pcsy]) {
					wp[0] = ipoints[i].v[0];
					wp[1] = ipoints[i].v[1];
					wp[2] = ipoints[i].v[2];
					for (e = 0; e < p->inputChan; e++)
						dwhite[e] = ipoints[i].p[e];
				}
				if (ipoints[i].v[pcsy] < bp[pcsy]) {
					bp[0] = ipoints[i].v[0];
					bp[1] = ipoints[i].v[1];
					bp[2] = ipoints[i].v[2];
					for (e = 0; e < p->inputChan; e++)
						dblack[e] = ipoints[i].p[e];
				}
			}
			if (p->pcs != icSigXYZData) {
				icmLab2XYZ(&icmD50, wp, wp);
				icmLab2XYZ(&icmD50, bp, bp);
			}

		} else {	/* Output or Display device */
			/* We assume that the output target is well behaved, */
			/* and that it includes a white point patch. */
			int nw = 0;

			wp[0] = wp[1] = wp[2] = 0.0;

			switch (h->colorSpace) {
	
				case icSigCmykData:
					for (i = 0; i < nodp; i++) {
						if (ipoints[i].p[0] < 0.001
						 && ipoints[i].p[1] < 0.001
						 && ipoints[i].p[2] < 0.001
						 && ipoints[i].p[3] < 0.001) {
							wp[0] += ipoints[i].v[0];
							wp[1] += ipoints[i].v[1];
							wp[2] += ipoints[i].v[2];
							nw++;
						}
					}
					for (e = 0; e < p->inputChan; e++) {
						dwhite[e] = 0.0;
						dblack[e] = 1.0;
					}
					break;
				case icSigCmyData:
					for (i = 0; i < nodp; i++) {
						if (ipoints[i].p[0] < 0.001
						 && ipoints[i].p[1] < 0.001
						 && ipoints[i].p[2] < 0.001) {
							wp[0] += ipoints[i].v[0];
							wp[1] += ipoints[i].v[1];
							wp[2] += ipoints[i].v[2];
							nw++;
						}
					}
					for (e = 0; e < p->inputChan; e++) {
						dwhite[e] = 0.0;
						dblack[e] = 1.0;
					}
					break;
				case icSigRgbData:
					for (i = 0; i < nodp; i++) {
						if (ipoints[i].p[0] > 0.999
						 && ipoints[i].p[1] > 0.999
						 && ipoints[i].p[2] > 0.999) {
							wp[0] += ipoints[i].v[0];
							wp[1] += ipoints[i].v[1];
							wp[2] += ipoints[i].v[2];
							nw++;
						}
					}
					for (e = 0; e < p->inputChan; e++) {
						dwhite[e] = 1.0;
						dblack[e] = 0.0;
					}
					break;
	
				case icSigGrayData: {	/* Could be additive or subtractive */
					double minwp[3], maxwp[3];
					int nminwp = 0, nmaxwp = 0;

					minwp[0] = minwp[1] = minwp[2] = 0.0;
					maxwp[0] = maxwp[1] = maxwp[2] = 0.0;

					/* Look for both */
					for (i = 0; i < nodp; i++) {
						if (ipoints[i].p[0] < 0.001)
							minwp[0] += ipoints[i].v[0];
							minwp[1] += ipoints[i].v[1];
							minwp[2] += ipoints[i].v[2]; {
							nminwp++;
						}
						if (ipoints[i].p[0] > 0.999)
							maxwp[0] += ipoints[i].v[0];
							maxwp[1] += ipoints[i].v[1];
							maxwp[2] += ipoints[i].v[2]; {
							nmaxwp++;
						}
					}
					if (nminwp > 0) {			/* Subtractive */
						wp[0] = minwp[0];
						wp[1] = minwp[1];
						wp[2] = minwp[2];
						nw = nminwp;
						for (e = 0; e < p->inputChan; e++) {
							dwhite[e] = 0.0;
							dblack[e] = 1.0;
						}
						if (minwp[pcsy]/nminwp < (0.5 * pcsymax))
							nw = 0;					/* Looks like a mistake */
					}
					if (nmaxwp > 0				/* Additive */
					 && (nminwp == 0 || maxwp[pcsy]/nmaxwp > minwp[pcsy]/nminwp)) {
						wp[0] = maxwp[0];
						wp[1] = maxwp[1];
						wp[2] = maxwp[2];
						nw = nmaxwp;
						for (e = 0; e < p->inputChan; e++) {
							dwhite[e] = 1.0;
							dblack[e] = 0.0;
						}
						if (maxwp[pcsy]/nmaxwp < (0.5 * pcsymax))
							nw = 0;					/* Looks like a mistake */
					}
					break;
				}

				default:
					xicp->errc = 1;
					sprintf(xicp->err,"set_icxLuLut: can't handle color space %s",
					                           icm2str(icmColorSpaceSignature, h->colorSpace));
					p->del((icxLuBase *)p);
					return NULL;
					break;
			}

			if (nw == 0) {
				xicp->errc = 1;
				sprintf(xicp->err,"set_icxLuLut: can't handle test points without a white patch");
				p->del((icxLuBase *)p);
				return NULL;
			}
			wp[0] /= (double)nw;
			wp[1] /= (double)nw;
			wp[2] /= (double)nw;
			if (p->pcs != icSigXYZData) 	/* Convert white point to XYZ */
				icmLab2XYZ(&icmD50, wp, wp);
		}

		if (flags & ICX_VERBOSE) {
			double lwp[3];
			icmXYZ2Lab(&icmD50, lwp, wp);
			printf("Approximate White point XYZ = %f %f %f, Lab = %f %f %f\n",
			        wp[0],wp[1],wp[2],lwp[0],lwp[1],lwp[2]);
		}

	/* Else not ICX_SET_WHITE */
	} else {
		icmXYZ2Ary(wp, icmD50);		/* Set a default value - D50 */
	}

	if (h->colorSpace == icSigGrayData) {	/* Don't use device or PCS curves for monochrome */
		p->noisluts = p->noipluts = p->nooluts = 1;
	}

	if ((flags & ICX_VERBOSE) && (p->noisluts == 0 || p->noipluts == 0 || p->nooluts == 0))
		printf("Creating optimised per channel curves\n");

	/* Set the target CLUT grid resolution so in/out curves can be optimised for it */
	for (e = 0; e < p->inputChan; e++)
		gres[e] = p->lut->clutPoints;

	/* Setup and then create xfit object that does most of the work */
	{
		int xfflags = 0;		/* xfit flags */
		double in_min[MXDI];	/* Input value scaling minimum */
		double in_max[MXDI];	/* Input value scaling maximum */
		double out_min[MXDO];	/* Output value scaling minimum */
		double out_max[MXDO];	/* Output value scaling maximum */
		int iluord, sluord, oluord;
		int iord[MXDI];			/* Input curve orders */
		int sord[MXDI];			/* Input sub-grid curve orders */
		int oord[MXDO];			/* Output curve orders */

		optcomb tcomb = oc_ipo;	/* Create all by default */

		if ((xf = CAT2(new_, xfit)()) == NULL) {
			p->pp->errc = 2;
			sprintf(p->pp->err,"Creation of xfit object failed");
			p->del((icxLuBase *)p);
			return NULL;
		}
			
		/* Setup for optimising run */
		if (p->noisluts)
			tcomb &= ~oc_i;

		if (p->noipluts)
			tcomb &= ~oc_p;

		if (p->nooluts)
			tcomb &= ~oc_o;

		if (flags & ICX_VERBOSE)
			xfflags |= XFIT_VERB;

		if (flags & ICX_SET_WHITE) {
			xfflags |= XFIT_OUT_WP_REL;
			if (p->pcs != icSigXYZData)
				xfflags |= XFIT_OUT_LAB;
		}

		/* With current B2A code, make sure a & b curves */
		/* pass through zero. */
		if (p->pcs == icSigLabData) {
			xfflags |=XFIT_OUT_ZERO;
		}

		/* Let xfit create the clut */
		xfflags |= XFIT_MAKE_CLUT;

		/* Set the curve orders for input (device) and output (PCS) */
		if (quality >= 3) {				/* Ultra high */
			iluord = 25;			
			sluord = 4;			
			oluord = 25;			
		} else if (quality == 2) {		/* High */
			iluord = 20;			
			sluord = 3;			
			oluord = 20;			
		} else if (quality == 1) {		/* Medium */
			iluord = 17;			
			sluord = 2;			
			oluord = 17;			
		} else {						/* Low */
			iluord = 10;			
			sluord = 1;			
			oluord = 10;			
		}
		for (e = 0; e < p->inputChan; e++) {
			iord[e] = iluord;
			sord[e] = sluord;
			in_min[e] = p->inmin[e];
			in_max[e] = p->inmax[e];
		}

		for (f = 0; f < p->outputChan; f++) {
			oord[f] = oluord;
			out_min[f] = p->outmin[f];
			out_max[f] = p->outmax[f];

			/* Hack to prevent a convex L curve pushing */
			/* the clut L values above the maximum value */
			/* that can be represented, causing clipping. */
			/* Do this by making sure that the L curve pivots */
			/* through 100.0 to 100.0 */
			if (f == 0 && p->pcs == icSigLabData) {
				if (out_min[f] < 0.0001 && out_max[f] > 100.0) {
					out_max[f] = 100.0;	
				}
			}
		}

		/* Create input, sub and output per channel curves (if configured), */
		/* adjust for white point to make relative (if configured), */
		/* and create clut rspl using xfit class. */
		/* The true white point for the returned curves and rspl is returned. */
		if (xf->fit(xf, xfflags, p->inputChan, p->outputChan,
			rsplflags, wp, dwhite, 
		    ipoints, nodp, in_min, in_max, gres, out_min, out_max,
		    smooth, oavgdev, iord, sord, oord, tcomb,
		   (void *)p, xfit_to_de2, xfit_to_dde2) != 0) {
			p->pp->errc = 2;
			sprintf(p->pp->err,"xfit fitting failed");
			xf->del(xf);
			p->del((icxLuBase *)p);
			return NULL;
		}  

		/* - - - - - - - - - - - - - - - */
		/* Set the xicc input curve rspl */
		for (e = 0; e < p->inputChan; e++) {
			curvectx cx;
	
			cx.xf = xf;
			cx.oix = -1;
			cx.iix = e;

			if ((p->inputTable[e] = new_rspl(1, 1)) == NULL) {
				p->pp->errc = 2;
				sprintf(p->pp->err,"Creation of input table rspl failed");
				xf->del(xf);
				p->del((icxLuBase *)p);
				return NULL;
			}

			p->inputTable[e]->set_rspl(p->inputTable[e], 0,
			           (void *)&cx, set_linfunc,
    			       &p->ninmin[e], &p->ninmax[e],
			           (int *)&p->lut->inputEnt,
			           &p->ninmin[e], &p->ninmax[e]);
		}

		/* - - - - - - - - - - - - - - - */
		/* Set the xicc output curve rspl */
		for (f = 0; f < p->outputChan; f++) {
			curvectx cx;

			cx.xf = xf;
			cx.iix = -1;
			cx.oix = f;

			if ((p->outputTable[f] = new_rspl(1, 1)) == NULL) {
				p->pp->errc = 2;
				sprintf(p->pp->err,"Creation of output table rspl failed");
				xf->del(xf);
				p->del((icxLuBase *)p);
				return NULL;
			}

			p->outputTable[f]->set_rspl(p->outputTable[f], 0,
		                (void *)&cx, set_linfunc,
			            &p->noutmin[f], &p->noutmax[f],
			            (int *)&p->lut->outputEnt,
			            &p->noutmin[f], &p->noutmax[f]);

		}
	}

	if (flags & ICX_VERBOSE)
		printf("Creating fast inverse input lookups\n");

	/* Create rspl based reverse input lookups used in ink limit function. */
	for (e = 0; e < p->inputChan; e++) {
		int res = 256;
		curvectx cx;

		cx.xf = xf;
		cx.oix = -1;
		cx.iix = e;

		if ((p->revinputTable[e] = new_rspl(1, 1)) == NULL) {
			p->pp->errc = 2;
			sprintf(p->pp->err,"Creation of reverse input table rspl failed");
			xf->del(xf);
			p->del((icxLuBase *)p);
			return NULL;
		}
		p->revinputTable[e]->set_rspl(p->revinputTable[e], 0,
		           (void *)&cx, icxLuLut_invinput_func,
		           &p->ninmin[e], &p->ninmax[e], &res, &p->ninmin[e], &p->ninmax[e]);
	}


	if (flags & ICX_VERBOSE)
		printf("Compensate scattered data for input curves\n");

	/* ------------------------------- */
	/* Set clut lookup table from xfit */
	p->clutTable = xf->clut;
	xf->clut = NULL;

	/* Setup all the clipping, ink limiting and auxiliary stuff, */
	/* in case a reverse call is used. Need to avoid relying on inking */
	/* stuff that makes use of the white/black points, since they haven't */
	/* been set up properly yet. */
	if (setup_ink_icxLuLut(p, ink, 0) != 0) {
		xf->del(xf);
		p->del((icxLuBase *)p);
		return NULL;
	}

	/* Deal with finalizing white/black points */
	if (flags & (ICX_SET_WHITE | ICX_SET_BLACK)) {

		if ((flags & ICX_SET_WHITE) && (flags & ICX_VERBOSE)) {
			double lwp[3];
			icmXYZ2Lab(&icmD50, lwp, wp);
			printf("White point XYZ = %f %f %f, Lab = %f %f %f\n",
			        wp[0],wp[1],wp[2],lwp[0],lwp[1],lwp[2]);
		}

		/* Lookup the black point */
		if (flags & ICX_SET_BLACK) { /* Black Point Tag: */
			co bcc;

			if (flags & ICX_VERBOSE)
				printf("Find black point\n");

			/* For CMYK devices, we choose a black point that is in */
			/* the same Lab vector direction as K, with the minimum L value. */
			if (h->deviceClass != icSigInputClass
			 && h->colorSpace == icSigCmykData) {
				bfinds bfs;					/* Callback context */
				double sr[MXDO];			/* search radius */
				double tt[MXDO];			/* Temporary */
				int trial;
				double brv;

				/* Setup callback function context */
				bfs.p = p;

				/* !!! we should use an accessor funcion of xfit !!! */
				for (i = 0; i < 3; i++) {
					for (j = 0; j < 3; j++) {
						bfs.toAbs[i][j] = xf->toAbs[i][j];
					}
				}

				/* Lookup abs Lab value of white point */
				icmXYZ2Lab(&icmD50, bfs.p1, wp);

				/* Now figure abs Lab value of K only, as the direction */
				/* to use for the rich black. */
				for (e = 0; e < p->inputChan; e++)
					bcc.p[e] = 0.0;
				bcc.p[3] = 1.0;

				p->input(p, bcc.p, bcc.p);			/* Through input tables */
				p->clutTable->interp(p->clutTable, &bcc);	/* Through clut */
				p->output(p, bcc.v, bcc.v);		/* Through the output tables */

				if (p->pcs != icSigXYZData) 	/* Convert PCS to XYZ */
					icmLab2XYZ(&icmD50, bcc.v, bcc.v);

				/* Convert from relative to Absolute colorimetric */
				icmMulBy3x3(tt, xf->toAbs, bcc.v);
				icmXYZ2Lab(&icmD50, bfs.p2, tt); /* Convert K only black point to Lab */
 
				if (flags & ICX_VERBOSE)
					printf("K only value (Lab) = %f %f %f\n",bfs.p2[0], bfs.p2[1], bfs.p2[2]);

				/* Find the device black point using optimization */
				/* Do several trials to avoid local minima */
				for (j = 0; j < p->inputChan; j++) { 
					tt[j] = bcc.p[j] = 0.5;		/* Starting point */
					sr[j] = 0.1;
				}
				brv = 1e38;
				for (trial = 0; trial < 20; trial++) {
					double rv;			/* Temporary */

					if (powell(&rv, p->inputChan, tt, sr, 0.00001, 500, bfindfunc, (void *)&bfs) == 0) {
//printf("~1 trial %d, rv %f bp %f %f %f %f\n",trial,rv,tt[0],tt[1],tt[2],tt[3]);
						if (rv < brv) {
							brv = rv;
							for (j = 0; j < p->inputChan; j++)
								bcc.p[j] = tt[j];
						}
					}
					for (j = 0; j < p->inputChan; j++) {
						tt[j] = bcc.p[j] + d_rand(-0.3, 0.3);
						if (tt[j] < 0.0)
							tt[j] = 0.0;
						else if (tt[j] > 1.0)
							tt[j] = 1.0;
					}
				}
				if (brv > 1000.0)
					error ("set_icxLuLut: Black point powell failed");

				for (j = 0; j < p->inputChan; j++) { /* Make sure device values are in range */
					if (bcc.p[j] < 0.0)
						bcc.p[j] = 0.0;
					else if (bcc.p[j] > 1.0)
						bcc.p[j] = 1.0;
				}
				/* Now have device black in bcc.p[] */

			/* Else not a CMYK output device, */
			/* use the previously determined device black value. */
			} else {
				for (e = 0; e < p->inputChan; e++)
					bcc.p[e] = dblack[e];
			}

			/* Lookup the PCS for the device black: */
			p->input(p, bcc.p, bcc.p);			/* Through input tables */

			p->clutTable->interp(p->clutTable, &bcc);	/* Through clut */
			p->output(p, bcc.v, bcc.v);		/* Through the output tables */

			if (p->pcs != icSigXYZData) 	/* Convert PCS to XYZ */
				icmLab2XYZ(&icmD50, bcc.v, bcc.v);

			/* Convert from relative to Absolute colorimetric */
			icmMulBy3x3(bp, xf->toAbs, bcc.v);

			/* Got XYZ black point in bp[] */
			if (flags & ICX_VERBOSE) {
				double labbp[3];
				icmXYZ2Lab(&icmD50, labbp, bp);
				printf("Black point XYZ = %f %f %f, Lab = %f %f %f\n",
				        bp[0],bp[1],bp[2],labbp[0],labbp[1],labbp[2]);
			}
		}

		/* And write them */
		if (flags & ICX_SET_WHITE) { /* White Point Tag: */
			icmXYZArray *wo;
			if ((wo = (icmXYZArray *)icco->read_tag(
			           icco, icSigMediaWhitePointTag)) == NULL)  {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_white_black: couldn't find white tag");
				xf->del(xf);
				p->del((icxLuBase *)p);
				return NULL;
			}
			if (wo->ttype != icSigXYZArrayType) {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_white_black: white tag has wrong type");
				xf->del(xf);
				p->del((icxLuBase *)p);
				return NULL;
			}

			wo->size = 1;
			wo->allocate((icmBase *)wo);	/* Allocate space */
			wo->data[0].X = wp[0];
			wo->data[0].Y = wp[1];
			wo->data[0].Z = wp[2];
		}
		if (flags & ICX_SET_BLACK) { /* Black Point Tag: */
			icmXYZArray *wo;
			if ((wo = (icmXYZArray *)icco->read_tag(
			           icco, icSigMediaBlackPointTag)) == NULL)  {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_white_black: couldn't find black tag");
				xf->del(xf);
				p->del((icxLuBase *)p);
				return NULL;
				}
			if (wo->ttype != icSigXYZArrayType) {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_white_black: black tag has wrong type");
				xf->del(xf);
				p->del((icxLuBase *)p);
				return NULL;
			}

			wo->size = 1;
			wo->allocate((icmBase *)wo);	/* Allocate space */
			wo->data[0].X = bp[0];
			wo->data[0].Y = bp[1];
			wo->data[0].Z = bp[2];
		}
		if ((flags & ICX_SET_WHITE) || (flags & ICX_SET_BLACK)) {
			/* Make sure ICC white/black point lookup notices the new white and black points */
			p->plu->init_wh_bk(p->plu);
		}

		/* Setup the clut clipping, ink limiting and auxiliary stuff again */
		/* since re_set_rspl will have invalidated */
		if (setup_ink_icxLuLut(p, ink, 1) != 0) {
			xf->del(xf);
			p->del((icxLuBase *)p);
			return NULL;
		}
	}

	/* Done with xfit now */
	xf->del(xf);
	xf = NULL;

	if (setup_clip_icxLuLut(p) != 0) {
		p->del((icxLuBase *)p);
		return NULL;
	}

	/* ------------------------------- */

	/* Use our rspl's to set the icc Lut AtoB table values. */
	/* Use helper function to do the hard work. */
	if (p->lut->set_tables(p->lut, ICM_CLUT_SET_EXACT, (void *)p,
			h->colorSpace, 				/* Input color space */
			h->pcs,						/* Output color space */
			set_input,					/* Input transfer function, Dev->Dev' */
			NULL, NULL,					/* Use default Maximum range of Dev' values */
			set_clut,					/* Dev' -> PCS' transfer function */
			NULL, NULL,					/* Use default Maximum range of PCS' values */
			set_output					/* Linear output transform PCS'->PCS */
	) != 0)
		error("Setting 16 bit %s->%s Lut failed: %d, %s",
			     icm2str(icmColorSpaceSignature, h->colorSpace),
			     icm2str(icmColorSpaceSignature, h->pcs),
		         p->pp->pp->errc,p->pp->pp->err);

	/* Init a CAM model in case it will be used (ie. in profile with camclip flag) */
	if (vc != NULL)		/* One has been provided */
		p->vc  = *vc;		/* Copy the structure */
	else
		xicc_enum_viewcond(xicp, &p->vc, -1, NULL, 0);	/* Use a default */
	p->cam = new_icxcam(cam_default);
	p->cam->set_view(p->cam, p->vc.Ev, p->vc.Wxyz, p->vc.Yb, p->vc.La, p->vc.Lv, p->vc.Yf, p->vc.Fxyz, XICC_USE_HK);
	
	if (flags & ICX_VERBOSE)
		printf("Done A to B table creation\n");

	return (icxLuBase *)p;
}

/* ========================================================== */
/* Gamut boundary code.                                       */
/* ========================================================== */


/* Context for creating gamut boundary points fro, xicc */
typedef struct {
	gamut *g;			/* Gamut being created */
	icxLuLut *x;		/* xLut we are working from */
	icxLuBase *flu;		/* Forward xlookup */
	double in[MAX_CHAN];	/* Device input values */
} lutgamctx;

/* Function to hand to zbrent to find a clut input' value at the ink limit */
/* Returns value < 0.0 when within gamut, > 0.0 when out of gamut */
static double icxLimitFind(void *fdata, double tp) {
	int i;
	double in[MAX_CHAN];
	lutgamctx *p = (lutgamctx *)fdata;
	double tt;

	for (i = 0; i < p->x->inputChan; i++) {
		in[i] = tp * p->in[i];		/* Scale given input value */
	}
	
	tt = icxLimitD((void *)p->x, in);

	return tt;
}

/* Function to pass to rspl to create gamut boundary from */
/* forward xLut transform grid points. */
static void
lutfwdgam_func(
	void *pp,			/* lutgamctx structure */
	double *out,		/* output' value at clut grid point (ie. PCS' value) */
	double *in			/* input' value at clut grid point (ie. device' value) */
) {
	lutgamctx *p    = (lutgamctx *)pp;
	double pcso[3];	/* PCS output value */

	/* Figure if we are over the ink limit. */
	if (   (p->x->ink.tlimit >= 0.0 || p->x->ink.klimit >= 0.0)
	    && icxLimitD((void *)p->x, in) > 0.0) {
		int i;
		double sf;

		/* We are, so use the bracket search to discover a scale */
		/* for the clut input' value that will put us on the ink limit. */

		for (i = 0; i < p->x->inputChan; i++)
			p->in[i] = in[i];

		if (zbrent(&sf, 0.0, 1.0, 1e-4, icxLimitFind, pp) != 0) {
			return;		/* Give up */
		}

		/* Compute ink limit value */
		for (i = 0; i < p->x->inputChan; i++)
			p->in[i] = sf * in[i];
		
		/* Compute the clut output for this clut input */
		p->x->clut(p->x, pcso, p->in);	
		p->x->output(p->x, pcso, pcso);	
		p->x->out_abs(p->x, pcso, pcso);	
	} else {	/* No ink limiting */
		/* Convert the clut PCS' values to PCS output values */
		p->x->output(p->x, pcso, out);
		p->x->out_abs(p->x, pcso, pcso);	
	}

	/* Expand the gamut surface with this point */
	p->g->expand(p->g, pcso);

	/* Leave out[] unchanged */
}


/* Function to pass to rspl to create gamut boundary from */
/* backwards Lut transform. */
static void
lutbwdgam_func(
	void *pp,			/* lutgamctx structure */
	double *out,		/* output value */
	double *in			/* input value */
) {
	lutgamctx *p    = (lutgamctx *)pp;
	double devo[MAX_CHAN];	/* Device output value */
	double pcso[3];	/* PCS output value */

	/* Convert the clut values to device output values */
	p->x->output(p->x, devo, out);		/* (Device never uses out_abs()) */

	/* Convert from device values to PCS values */
	p->flu->lookup(p->flu, pcso, devo);

	/* Expand the gamut surface with this point */
	p->g->expand(p->g, pcso);

	/* Leave out[] unchanged */
}

/* Given an xicc lookup object, return a gamut object. */
/* Note that the PCS must be Lab or Jab */
/* An icxLuLut type must be icmFwd or icmBwd, */
/* and for icmFwd, the ink limit (if supplied) will be applied. */
/* Return NULL on error, check errc+err for reason */
static gamut *icxLuLutGamut(
icxLuBase   *plu,		/* this */
double       detail		/* gamut detail level, 0.0 = def */
) {
	xicc     *p = plu->pp;				/* parent xicc */
	icxLuLut *luluto = (icxLuLut *)plu;	/* Lookup xLut type object */
	icColorSpaceSignature ins, pcs;
	icmLookupFunc func;
	icRenderingIntent intent;
	double white[3], black[3];
	int inn;
	gamut *gam;

	/* get some details */
	plu->spaces(plu, &ins, &inn, NULL, NULL, NULL, &intent, &func, &pcs);

	if (func != icmFwd && func != icmBwd) {
		p->errc = 1;
		sprintf(p->err,"Creating Gamut surface for anything other than Device <-> PCS is not supported.");
		return NULL;
	}

	if (pcs != icSigLabData && pcs != icxSigJabData) {
		p->errc = 1;
		sprintf(p->err,"Creating Gamut surface PCS of other than Lab or Jab is not supported.");
		return NULL;
	}

	if (func == icmFwd) {
		lutgamctx cx;

		cx.g = gam = new_gamut(detail, pcs == icxSigJabData);
		cx.x = luluto;

		luluto->clutTable->scan_rspl(
			luluto->clutTable,	/* this */
			0,					/* Combination of flags */
			(void *)&cx,		/* Opaque function context */
			lutfwdgam_func		/* Function to set from */
		);

		if (detail == 0.0)
			detail = 10.0;

		/* If the gamut is more than cursary, add some more detail surface points */
		if (detail < 20.0) {
			int res;
			DCOUNT(co, MAX_CHAN, inn, 0, 0, 2);
		
			res = (int)(500.0/detail);	/* Establish an appropriate sampling density */
			if (res < 10)
				res = 10;

			/* Itterate over all the faces in the device space */
			DC_INIT(co);
			while(!DC_DONE(co)) {		/* Count through the corners of hyper cube */
				int e, m1, m2;
				double in[MAX_CHAN];
				double out[3];
		
				for (e = 0; e < inn; e++)
					in[e] = (double)co[e];		/* Base value */

   				/* Figure if we are over the ink limit. */
				if ((luluto->ink.tlimit >= 0.0 || luluto->ink.klimit >= 0.0)
			        && icxLimit((void *)luluto, in) > 0.0) {
					DC_INC(co);
					continue;		/* Skip points over limit */
				}

				/* Scan only device surface */
				for (m1 = 0; m1 < inn; m1++) {		/* Choose first coord to scan */
					if (co[m1] != 0)
						continue;					/* Not at lower corner */
					for (m2 = m1 + 1; m2 < inn; m2++) {	/* Choose second coord to scan */
						int x, y;
		
						if (co[m2] != 0)
							continue;					/* Not at lower corner */
		
						for (e = 0; e < inn; e++)
							in[e] = (double)co[e];		/* Base value */
		
						for (x = 0; x < res; x++) {				/* step over surface */
							in[m1] = x/(res - 1.0);
							for (y = 0; y < res; y++) {
								in[m2] = y/(res - 1.0);

				   				/* Figure if we are over the ink limit. */
								if (   (luluto->ink.tlimit >= 0.0 || luluto->ink.klimit >= 0.0)
							        && icxLimit((void *)luluto, in) > 0.0) {
									continue;		/* Skip points over limit */
								}
		
								luluto->lookup((icxLuBase *)luluto, out, in);
								gam->expand(gam, out);
							}
						}
					}
				}
				/* Increment index within block */
				DC_INC(co);
			}
		}

		/* Now set the cusp points by itterating through colorant 0 & 100% combinations */
		/* If we know what sort of space it is: */
		if (ins == icSigRgbData || ins == icSigCmyData || ins == icSigCmykData) {
			DCOUNT(co, 3, 3, 0, 0, 2);

			gam->setcusps(gam, 0, NULL);
			DC_INIT(co);
			while(!DC_DONE(co)) {
				int e;
				double in[MAX_CHAN];
				double out[3];
		
				if (!(co[0] == 0 && co[1] == 0 && co[2] == 0)
				 && !(co[0] == 1 && co[1] == 1 && co[2] == 1)) {	/* Skip white and black */
					for (e = 0; e < 3; e++)
						in[e] = (double)co[e];
					in[e] = 0;					/* K */
		
					/* Always use the device->PCS conversion */
					if (luluto->lookup((icxLuBase *)luluto, out, in) > 1)
						error ("%d, %s",p->errc,p->err);
					gam->setcusps(gam, 3, out);
				}

				DC_INC(co);
			}
			gam->setcusps(gam, 2, NULL);
		} else {	/* Do all ink combinations and hope we can sort it out */
			DCOUNT(co, MAX_CHAN, inn, 0, 0, 2);

			gam->setcusps(gam, 0, NULL);
			DC_INIT(co);
			while(!DC_DONE(co)) {
				int e;
				double in[MAX_CHAN];
				double out[3];
		
				for (e = 0; e < inn; e++)
					in[e] = (double)co[e];
	
	   			/* Figure if we are over the ink limit. */
				if ((luluto->ink.tlimit >= 0.0 || luluto->ink.klimit >= 0.0)
			        && icxLimit((void *)luluto, in) > 0.0) {
					DC_INC(co);
					continue;		/* Skip points over limit */
				}
	
				luluto->lookup((icxLuBase *)luluto, out, in);
				gam->setcusps(gam, 1, out);

				DC_INC(co);
			}
			gam->setcusps(gam, 2, NULL);
		}

	} else { /* Must be icmBwd */
		lutgamctx cx;
	
		/* Get an appropriate device to PCS conversion */
		switch (intent) {
			/* If it is relative */
			case icmDefaultIntent:					/* Shouldn't happen */
			case icPerceptual:
			case icRelativeColorimetric:
			case icSaturation:
				intent = icRelativeColorimetric;	/* Choose relative */
				break;
			/* If it is absolute */
			case icAbsoluteColorimetric:
			case icxAppearance:
			case icxAbsAppearance:
				break;								/* Leave unchanged */
			default:
				break;
		}
		if ((cx.flu = p->get_luobj(p, 0, icmFwd, intent, pcs, icmLuOrdNorm,
		                              &plu->vc, NULL)) == NULL) {
			return NULL;	/* oops */
		}

		cx.g = gam = new_gamut(detail, pcs == icxSigJabData);

		cx.x = luluto;

		luluto->clutTable->scan_rspl(
			luluto->clutTable,	/* this */
			0,					/* Combination of flags */
			(void *)&cx,		/* Opaque function context */
			lutbwdgam_func		/* Function to set from */
		);

		cx.flu->del(cx.flu);	/* Done with the fwd conversion */

	}

	/* Put the white and black points in the gamut */
	plu->rel_wh_bk_points(plu, white, black);
	gam->setwb(gam, white, black);

#ifdef NEVER	/* Not sure if this is a good idea ?? */
	/* Since we know this is a profile, force the space and gamut points to be the same */
	gam->getwb(gam, NULL, NULL, white, black);	/* Get the actual gamut white and black points */
	gam->setwb(gam, white, black);				/* Put it back as colorspace one */
#endif

	return gam;
}

#ifdef DEBUG
#undef DEBUG 	/* Limit extent to this file */
#endif






