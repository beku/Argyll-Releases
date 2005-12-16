
/*
 * International Color Consortium color transform expanded support
 *
 * Author:  Graeme W. Gill
 * Date:    2/7/00
 * Version: 1.00
 *
 * Copyright 2000, 2001 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the LICENCE.TXT file for licencing details.
 *
 * This is the second major version of xlut.c (originally called xlut2.c)
 * Based on the old xlut.c code (preserved as xlut1.c)
 */

/*
 * This module provides the expanded icclib functionality
 * for lut based profiles.
 * This file is #included in xicc.c, to keep its functions private.
 *
 * This version creates both input and output 1D luts by
 * optimising the accuracy of the profile for a linear clut.
 *
 * Note that the black point finding code needs improving.
 *
 */

#define USE_CIE94_DE	/* Use CIE94 delta E measure when creating in/out curves */

#undef DEBUG 			/* Verbose debug information */
#undef DEBUG_PLOT 		/* Plot in & out curves */
#undef INK_LIMIT_TEST	/* Turn off input tables for ink limit testing */
#undef CHECK_ILIMIT		/* Do sanity checks on meeting ink limit */

#undef NODDV			/* Use slow non d/dv powell */

/* Weights in shaper parameters, to minimise unconstrained "wiggles" */
#define SHAPE_WEIGHT 1.0	/* Overal shaper weight contribution */
#define SHAPE_BASE  1.0		/* 0 & 1 harmonic parameter weight */
#define SHAPE_HBASE 2.0		/* 3rd harmonic and above base parameter weight */

#undef HACK_ILORD /* 10 */			/* Override input per channel curve order */
#undef HACK_OLORD /* 6 */			/* Override output per channel curve order */

/*
 * TTBD:
 *       Reverse lookup of Lab
 *       Make NEARCLIP the default ??
 *
 *       XYZ vector clipping isn't implemented.
 *
 *       Some of the error handling is crude. Shouldn't use
 *       error(), should return status.
 */

static double icxLimit(void *lcntx, double *in);
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
			lim = icxLimit((void *)p, in);
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
int icxLuLut_inv_out_abs(icxLuLut *p, double *out, double *in) {
	int rv = 0;

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
	return rv;
}

/* Do output->output' inverse lookup */
int icxLuLut_inv_output(icxLuLut *p, double *out, double *in) {
	int rv = 0;
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
	return rv;
}

/* Ink limit+gamut limit calculation function for xLuLut. */
/* Returns < 0.0 if output value is within limits, */
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
	double ind[MAX_CHAN];
	co tc;
	int e;

//printf("~~ limit got %f %f %f %f\n", in[0], in[1], in[2], in[3]);

	if ((tlim = p->ink.tlimit) < 0.0)
		tlim = (double)p->inputChan + 1.0;

	if ((klim = p->ink.klimit) < 0.0)
		klim = 2.0;

	/* Convert input' to input through revinput Luts */
	for (e = 0; e < p->inputChan; e++) {
		tc.p[0] = in[e];
		p->revinputTable[e]->interp(p->revinputTable[e], &tc);
		ind[e] = tc.v[0];
	}

	/* Compute amount outside total limit */
	{
		double sum;
		for (sum = 0.0, e = 0; e < p->inputChan; e++)
			sum += ind[e];
		val = sum - tlim;
	}

	/* Compute amount outside black limit */
	if (p->ink.klimit >= 0.0) {
		double kval;
		switch(p->natis) {
			case icSigCmykData:
				kval = ind[3] - klim;
				break;
			default:
				error("xlut: Unknown colorspace when black limit specified");
		}
		if (kval > val)
			val = kval;
	}

	/* Compute amount outside device value limits 0.0 - 1.0 */
	for (ovr = 0.0, e = 0; e < p->inputChan; e++) {
		if (ind[e] < 0.0) {
			if (-ind[e] > ovr)
				ovr = -ind[e];
		} else if (ind[e] > 1.0) {
			if ((ind[e] - 1.0) > ovr)
				ovr = ind[e] - 1.0;
		}
	}
	if (ovr > val)
		val = ovr;

//printf("~~ limit returning %f\n", val);
	return val;
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

/* Do output'->input' lookup with aux details */
/* Note than out[] will be used as the inking value if icxKrule is */
/* icxKvalue, icxKlocus or icxKl5l, and that the icxKrule value will be in the input' space. */
/* Note that the ink limit will be computed after converting input' to input */
/* auxt will override the inking rule, and auxr reflects the available auxiliary */
/* range there was to choose from. auxv was the actual auxiliary used. */
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
	double min[MXDI], max[MXDI];	/* Auxiliary locus range */
	int sauxt = 0;		/* Smooth auxiliary target was set */

	if (p->nearclip != 0)
		flags |= RSPL_NEARCLIP;			/* Use nearest clipping rather than clip vector */

//printf("~~1 Input is %f %f %f\n",in[0], in[1], in[2]);

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

		/* Compute auxiliary locus on the fly */
		nsoln = p->clutTable->rev_locus(
			p->clutTable,	/* rspl object */
			p->auxm,		/* Auxiliary mask */
			pp,				/* Input target and output solutions */
			min, max);		/* Returned locus of valid auxiliary values */
		
		if (nsoln == 0)
			xflags |= RSPL_WILLCLIP;	/* No valid locus, so we expect to have to clip */

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
			} else { /* p->ink.k_rule == icxKluma5 || icskl5l */
				/* Auxiliaries are driven by a rule and the output values */
				double rv, L;

				/* Figure out Luminance number */
				if (p->natos == icSigLabData
				 || p->natos == icxSigJabData) {
					L = 0.01 * in[0];
				} else if (p->natos == icSigXYZData) {
					double tt[3];
					icmXYZ2Lab(&icmD50, tt, in);
					L = 0.01 * tt[0];
				} else {	/* Hmm. How do we handle this ? */
					L = 0.5;
				}

				/* Normalise L to its possible range from min to max */
				L = (L - p->Lmin)/(p->Lmax - p->Lmin);

				/* Convert L to locus curve value */
				rv = icxKcurve(L, &p->ink.c);

				if (p->ink.k_rule == icxKluma5) {

					/* Set target black as K fraction within locus */

					for (e = 0; e < p->clutTable->di; e++) {
						if (p->auxm[e] != 0) {
							pp[0].p[e] = min[e] + rv * (max[e] - min[e]);
						}
					}

				} else {
					/* Create second curve, and use input locus to */
					/* blend between */

					double rv2;		/* Upper limit */

					/* Convert L to locus curve value */
					rv2 = icxKcurve(L, &p->ink.x);

					for (e = 0; e < p->clutTable->di; e++) {
						if (p->auxm[e] != 0) {
							double ii, iv;
							ii = out[e];					/* Input ink locus */
							if (ii < 0.0)
								ii = 0.0;
							else if (ii > 1.0)
								ii = 1.0;
							ii = (1.0 - ii) * rv + ii * rv2;/* Blend between locus rule curves */
							/* Out ink from output locus */
							pp[0].p[e] = min[e] + ii * (max[e] - min[e]);
						}
					}
				}
			}

			xflags |= RSPL_EXACTAUX;	/* Since we confine aux to locus */
		}

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
//printf("~~1 No auxiliary needed\n");
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

	/* If we clipped and we should clip in CAM Jab space, compute reverse */
	/* clip solution using our additional CAM space. */
	/* (Note that we don't support vector clip in CAM space at the moment) */
	if (crv != 0 && p->camclip && p->nearclip) {
		double tin[MXDO];	/* CAM space value to be inverted */
		co cpp;				/* Alternate CAM space solution */
		double cdist;		/* CAM clip distance */
		double bf;			/* Blend factor */

		if (nsoln != 1) {	/* This would be unexpected */
			error("~~~1 Unexpected failure to return 1 solution on clip for input to output table");
		}

		if (p->cclutTable == NULL) {	/* we haven't created this yet, so do so */
			if (icxLuLut_init_clut_camclip(p))
				error("~~~1 Error in creating CAM rspl for camclip");
		}

		/* Setup for reverse lookup */
		((icmLuLut *)p->absxyzlu)->output((icmLuLut *)p->absxyzlu, tin, in);
		((icmLuLut *)p->absxyzlu)->out_abs((icmLuLut *)p->absxyzlu, tin, tin);
		p->cam->XYZ_to_cam(p->cam, tin, tin);

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
//printf("~1 cdist = %f\n",cdist);
		bf = cdist/2.0;						/* 0.0 for PCS result, 1.0 for CAM result */
		if (bf > 1.0)
			bf = 1.0;
//printf("~1 raw blend = %f\n",bf);
		bf = bf * bf * (3.0 - 2.0 * bf);	/* Convert to spline blend */
//printf("~1 spline blend = %f\n",bf);

		/* Blend between solution values for PCS and CAM clip result. */
		/* We're hoping that the solutions are close, and expect them to be */
		/* that way when we're close to the gamut surface (since the cell */
		/* vertex values should be exact, irrespective of which space they're */
		/* in), but weird things could happen away from the surface. Weird */
		/* things can happen anyway with "clip to nearest", since this is not */
		/* guaranteed to be a smooth function, depending on the gamut surface */
		/* geometry. */
// ~1
//printf("~1 Clip blend between:\n");
//printf("~1 %f %f %f %f and\n",pp[0].p[0], pp[0].p[1], pp[0].p[2], pp[0].p[3]);
//printf("~1 %f %f %f %f\n",cpp.p[0], cpp.p[1], cpp.p[2], cpp.p[3]);

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
	double sum = icxLimit((void *)p, out);
	if (sum > 0.0)
		printf("xlut assert%s: icxLuLut_inv_clut returned outside limits by %f > tlimit %f\n",crv ? " (clip)" : "", sum, p->ink.tlimit);
}
#endif

//printf("~~1 Output is %f %f %f %f\n",out[0], out[1], out[2], out[3]);
	return crv;
}

/* Do output'->input' lookup, simple version */
/* Note than out[] will carry inking value if icxKrule is icxKvalue of icxKlocus */
/* and that the icxKrule value will be in the input' space. */
/* Note that the ink limit will be computed after converting input' to input */
int icxLuLut_inv_clut(icxLuLut *p, double *out, double *in) {
	return icxLuLut_inv_clut_aux(p, out, NULL, NULL, NULL, in);
}

/* Given the input' values in in[di], and the target output' values in out[fdi], */
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
	p->noiluts = 0;
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
		p->ink.tlimit >= 0.0 || p->ink.klimit >= 0.0 ? icxLimit : NULL,
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
			p->ink.tlimit >= 0.0 || p->ink.klimit >= 0.0 ? icxLimit : NULL,
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

		/* Get white/black in effective xlu pcs space */
		p->rel_wh_bk_points((icxLuBase *)p, wh, bk);

		/* Note we are assuming that pcs == output space (usual, but not enforced though!) */
		/* Convert from the overall output space to the clut output space, */
		/* allowing for output curves and colorspace change. */
		icxLuLut_inv_out_abs(p, wh, wh);
		icxLuLut_inv_output(p, wh, wh);
		icxLuLut_inv_out_abs(p, bk, bk);
		icxLuLut_inv_output(p, bk, bk);

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
	int *auxmp = NULL;
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
	double *in			/* inut value */
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
/* intent to consistent and values if Jab and/or icxAppearance */
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
	 && !p->mergeclut
	 && func == icmFwd)
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
			xicc_enum_viewcond(xicp, &p->vc, -1, 0);	/* Use a default */
		p->cam = new_icxcam(cam_default);
		p->cam->set_view(p->cam, p->vc.Ev, p->vc.Wxyz, p->vc.Yb, p->vc.La, p->vc.Lv, p->vc.Yf, p->vc.Fxyz, 1);
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
		           &p->ninmin[i], &p->ninmax[i], &p->lut->inputEnt, &p->ninmin[i], &p->ninmax[i]);
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
		           &p->noutmin[i], &p->noutmax[i], &p->lut->outputEnt, &p->noutmin[i], &p->noutmax[i]);
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
	cmin[0] = 0.0;		cmax[0] = 100.0;	/* Nominal Jab ranges */
	cmin[1] = -128.0;	cmax[1] = 128.0;
	cmin[2] = -128.0;	cmax[2] = 128.0;

	/* Get icm lookup we need for seting up and useng CAM icx clut */
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
		p->ink.tlimit >= 0.0 || p->ink.klimit >= 0.0 ? icxLimit : NULL,
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

#define MXLUORD 10		/* Shaper harmonic orders to use */
#define PCSD 3			/* PCS dimensions */

#define MXPARMS (((1 << MXDI) * PCSD) + (MXDI + PCSD) * (MXLUORD))


/* Optimisation mask legal combinations */
typedef enum {
	oc_i   = 1,		/* Input */
	oc_m   = 2,		/* Matrix */
	oc_im  = 3,		/* Input and matrix */
	oc_o   = 4,		/* Output */
	oc_mo  = 6,		/* Matrix and output */
	oc_imo = 7		/* Input, matrix and output */
} optcomb;

/* Context for optimising input and output luts */
typedef struct {
	int verb;				/* Verbose */
	optcomb opt_msk;		/* Optimisation mask: 3 = i+m, 2 = m, 6 = m+o, 7 = i+m+o */
	int opt_off;			/* Optimisation parameters offset */
	int opt_cnt;			/* Optimisation parameters count */
	int iluord[MXDI];		/* Input Shaper order actualy used (must be <= MXLUORD) */
	int oluord[PCSD];		/* Output Shaper order actualy used (must be <= MXLUORD) */
	double in_min[MXDI];	/* Input value scaling minimum */
	double in_max[MXDI];	/* Input value scaling maximum */
	double out_min[PCSD];	/* Output value scaling minimum */
	double out_max[PCSD];	/* Output value scaling maximum */
	int in_off;				/* Input  parameters offset */
	int in_offs[MXDI];		/* Input  parameter offsets for each channel from v[0] */
	int in_cnt;				/* Input  parameters count */
	int mat_off;			/* Matrix parameters offset from v[0] */
	int mat_cnt;			/* Matrix parameters count */
	int out_off;			/* Output parameters offset from v[0] */
	int out_offs[PCSD];		/* Output  parameter offsets for each channel from v[0] */
	int out_cnt;			/* Output parameters count */
	int tot_cnt;			/* Total parameter count */
	int symch;				/* Output channel being adjusted for symetry */
	icxLuLut *x;			/* xicc LuLut lookup object being created */
	double v[MXPARMS];		/* Holder for parameters */
							/* Optimisation parameters are layed out:         */
							/*                                                */
							/* Input curves:, di groups of iluord[e] parameters */
							/*                                                */
							/* Matrix: fdi groups of 2 ^ di parameters        */
							/*                                                */
							/* Output curves:, fdi groups of oluord[f] parameters */
							/*                                                */
	co *points;				/* List of test points as dev->Lab                */
	int nodp;				/* Number of data points                          */
} luopt;


/* Init the initial parameters and search area */
static void init_luopt(
luopt *p,
icxLuLut *x,
int *iord,			/* Order of input shaper curve */
int *oord			/* Order of output shaper curve */
) {
	int i, e, f;
	double *b;			/* Base of parameters for this section */
	int di, fdi;
	int poff;

	p->x = x;
	di   = x->inputChan;
	fdi  = x->outputChan;
	p->opt_msk = oc_imo;

	/* Sanity protect shaper orders and set scaling factors. */
	for (e = 0; e < di; e++) {
		if (iord[e] > MXLUORD)
			p->iluord[e] = MXLUORD;
		else
			p->iluord[e] = iord[e];
		p->in_min[e] = p->x->ninmin[e];
		p->in_max[e] = p->x->ninmax[e];
	}
	for (f = 0; f < fdi; f++) {
		if (oord[f] > MXLUORD)
			p->oluord[f] = MXLUORD;
		else
			p->oluord[f] = oord[f];
		p->out_min[f] = p->x->noutmin[f];
		p->out_max[f] = p->x->noutmax[f];

		/* Hack to prevent a convex L curve pushing */
		/* the clut L values above the maximum value */
		/* that can be represented, causing clipping. */
		/* Do this by making sure that the L curve pivots */
		/* through 100.0 to 100.0 */
		if (f == 0 && p->x->pcs == icSigLabData) {
			if (p->out_min[f] < 0.0001 && p->out_max[f] > 100.0) {
				p->out_max[f] = 100.0;	
			}
		}
	}

	/* Set parameter offset and count information */
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
		double ov[PCSD];
	
		/* Search the patch list to find the one closest to this colorant combination */
		for (k = 0; k < p->nodp; k++) {
			double dif = 0.0;
			for (j = 0; j < di; j++) {
				double tt;
				if (e & (1 << j))
					tt = 1.0 - p->points[k].p[j];
				else
					tt = p->points[k].p[j];
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
		if (x->pcs == icSigXYZData) 		/* Convert from Lab to XYZ PCS */
			icmLab2XYZ(&icmD50, ov, ov);
			
		for (f = 0; f < fdi; f++) {
			b[f * (1 << di) + e] = ov[f];
		}
	}

	/* Setup output curves to be liniear initially */
	b = p->v + p->out_off;
		
	for (f = 0; f < fdi; b += p->oluord[f], f++) {
		for (i = 0; i < p->oluord[f]; i++) {
			b[i] = 0.0;
		}
	}

#ifdef NEVER
printf("~1 in_off = %d\n",p->in_off);
printf("~1 in_cnt = %d\n",p->in_cnt);
printf("~1 mat_off = %d\n",p->mat_off);
printf("~1 mat_cnt = %d\n",p->mat_cnt);
printf("~1 out_off = %d\n",p->out_off);
printf("~1 out_cnt = %d\n",p->out_cnt);
printf("~1 tot_cnt = %d\n\n",p->tot_cnt);
#endif /* NEVER */
}

/* Set up for an optimisation run: */
/* Figure out parameters being optimised, */
/* copy them to start values, */
/* init and scale the search radius */
static void setup_luopt(
luopt *p,
double *v,		/* Return parameters hand to optimiser */
double *sa,		/* Return search radius to hand to optimiser */
double transrad,/* Nominal transfer curve radius, 0.0 - 3.0 */
double matrad	/* Nominal matrix radius, 0.0 - 1.0 */
) {
	int i;
	p->opt_off = -1;
	p->opt_cnt = 0;
	
	if (p->x->pcs == icSigLabData) /* scale search radius */
		matrad *= 100.0;

	if (p->opt_msk & oc_i) {
		p->opt_off = p->in_off;
		p->opt_cnt += p->in_cnt;

		for (i = 0; i < p->in_cnt; i++) {
			*v++ = p->v[p->in_off + i];
			*sa++ = transrad;
		}
	}
	if (p->opt_msk & oc_m) {
		if (p->opt_off < 0)
			p->opt_off = p->mat_off;
		p->opt_cnt += p->mat_cnt;

		for (i = 0; i < p->mat_cnt; i++) {
			*v++ = p->v[p->mat_off + i];
			*sa++ = matrad;
		}
	}
	if (p->opt_msk & oc_o) {
		if (p->opt_off < 0)
			p->opt_off = p->out_off;
		p->opt_cnt += p->out_cnt;

		for (i = 0; i < p->out_cnt; i++) {
			*v++ = p->v[p->out_off + i];
			*sa++ = transrad;
		}
	}
#ifdef NEVER
printf("~1 opt_msk = 0x%x\n",p->opt_msk);
printf("~1 opt_off = %d\n",p->opt_off);
printf("~1 opt_cnt = %d\n\n",p->opt_cnt);
#endif /* NEVER */
}

/* Diagnostic */
static void dump_luopt(
luopt *p
) {
	int i, e, f;
	double *b;			/* Base of parameters for this section */
	int di, fdi;
	di   = p->x->inputChan;
	fdi  = p->x->outputChan;

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

/* return a weighting for the magnitude of the in and out */
/* shaping parameters. This is to reduce unconstrained "wiggles" */
static double shapmag(
luopt  *p			/* Base of optimisation structure */
) {
	double tt, w;
	double *b;			/* Base of parameters for this section */
	int di =  p->x->inputChan;
	int fdi = p->x->outputChan;
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
				w = (k < 2) ? SHAPE_BASE : k * SHAPE_HBASE;
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
				w = (k < 2) ? SHAPE_BASE : k * SHAPE_HBASE;
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
luopt  *p,			/* Base of optimisation structure */
double *dav			/* Sum del's */
) {
	double tt, w;
	double *b, *c;			/* Base of parameters for this section */
	int di =  p->x->inputChan;
	int fdi = p->x->outputChan;
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
				w = (k < 2) ? SHAPE_BASE : k * SHAPE_HBASE;
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
				w = (k < 2) ? SHAPE_BASE : k * SHAPE_HBASE;
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
static double luoptfunc(void *edata, double *v) {
	luopt *p = (luopt *)edata;
	double rv = 0.0, smv;
	double tin[MXDI], out[PCSD];
	int di = p->x->inputChan;
	int fdi = p->x->outputChan;
	int i, e, f;

	/* Copy the parameters being optimised into luopt structure */
	for (i = 0; i < p->opt_cnt; i++)
		p->v[p->opt_off + i] = v[i];

	/* For all our data points */
	for (i = 0; i < p->nodp; i++) {
		double ev;
		int j;

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

		if (p->x->pcs == icSigXYZData) 		/* Convert to Lab */
			icmXYZ2Lab(&icmD50, out, out);
	
		/* Accumulate total delta E squared */
#ifdef USE_CIE94_DE
		ev = icmCIE94sq(out, p->points[i].v);
#else
		ev = icmLabDEsq(out, p->points[i].v);
#endif
		rv += ev;
	}

	/* Normalise error to be an average delta E squared */
	rv /= (double)p->nodp;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	smv = shapmag(p);
	rv += smv;

#ifdef DEBUG
printf("~1(%f)luoptfunc returning %f\n",smv,rv);
#endif

//	if (p->verb)
//		printf("."), fflush(stdout);
	return rv;
}

/* Shaper+Matrix optimisation function with partial derivatives, */
/* handed to conjgrad() */
static double dluoptfunc(void *edata, double *dv, double *v) {
	luopt *p = (luopt *)edata;
	double rv = 0.0, smv;
	double tin[MXDI], out[PCSD];

	double dav[MXPARMS];				/* Overall del due to del param vals */

	double dtin_iv[MXDI * MXLUORD];		/* Del in itrans out due to del itrans param vals */
	double dmato_mv[1 << MXDI];			/* Del in mat out due to del in matrix param vals */
	double dmato_tin[PCSD * MXDI];		/* Del in mat out due to del in matrix input values */
	double dout_ov[PCSD * MXLUORD];		/* Del in otrans out due to del in otrans param values */
	double dout_mato[PCSD];				/* Del in otrans out due to del in otrans input values */
	double dout_lab[PCSD][PCSD];		/* Del in out due to possible XYZ to Lab conversion */
	double de_dout[2][PCSD];			/* Del in delta E due to input Lab values */

	int di = p->x->inputChan;
	int fdi = p->x->outputChan;
	int i, jj, k, e, ee, f, ff;
	int dolab = 0;						/* Convert XYZ to Lab */

	if (p->x->pcs == icSigXYZData) 		/* Convert to Lab */
		dolab = 1;

	/* Copy the parameters being optimised into luopt structure */
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
		if (1 || p->opt_msk & oc_i) {
			for (e = 0; e < di; e++)
				tin[e] = icxdpSTransFunc(p->v + p->in_offs[e], &dtin_iv[p->in_offs[e] - p->in_off],
			                         p->iluord[e], p->points[i].p[e], p->in_min[e], p->in_max[e]);
		} else {
			for (e = 0; e < di; e++)
				tin[e] = icxSTransFunc(p->v + p->in_offs[e], p->iluord[e], p->points[i].p[e],  
			                       p->in_min[e], p->in_max[e]);
		}

		/* Apply matrix cube interpolation */
		if (1 || (p->opt_msk & oc_i) || (p->opt_msk & oc_m))
			icxdpdiCubeInterp(p->v + p->mat_off, dmato_mv, dmato_tin, fdi, di, out, tin);
		else
			icxCubeInterp(p->v + p->mat_off, fdi, di, out, tin);

		/* Apply output channel curves */
		for (f = 0; f < fdi; f++)
			out[f] = icxdpdiSTransFunc(p->v + p->out_offs[f],
			                           &dout_ov[p->out_offs[f] - p->out_off], &dout_mato[f],
			                           p->oluord[f], out[f], p->out_min[f], p->out_max[f]);

		if (dolab)				 		/* Convert to Lab */
			icxdXYZ2Lab(&icmD50, out, dout_lab, out);
	
		/* Compute dela E */
#ifdef USE_CIE94_DE
		ev = icxdCIE94sq(de_dout, out, p->points[i].v);
#else
		ev = icxdLabDEsq(de_dout, out, p->points[i].v);
#endif

		/* Accumulate total delta E squared */
		rv += ev;

		/* Compute and accumulate partial difference values for each parameter value */

		if (p->opt_msk & oc_i) {
			/* Input transfer parameters */
			for (ee = 0; ee < di; ee++) {				/* Parameter input chanel */
				for (k = 0; k < p->iluord[ee]; k++) {	/* Param within channel */
					double vv = 0.0;
					jj = p->in_offs[ee] - p->in_off + k;	/* Overall input trans param */

					if (dolab) {
						for (ff = 0; ff < 3; ff++) {		/* Lab channels */
							for (f = 0; f < 3; f++) {		/* XYZ channels */
								vv += de_dout[0][ff] * dout_lab[ff][f]
								    * dout_mato[f] * dmato_tin[f * di + ee] * dtin_iv[jj];
							}
						}
					} else {
						for (ff = 0; ff < 3; ff++) {		/* Lab channels */
							vv += de_dout[0][ff] * dout_mato[ff]
							    * dmato_tin[ff * di + ee] * dtin_iv[jj];
						}
					}
					dav[p->in_off + jj] += vv;
				}
			}
		}

		/* Matrix parameters */
		if (p->opt_msk & oc_m) {
			for (ff = 0; ff < fdi; ff++) {				/* Parameter output chanel */
				for (ee = 0; ee < (1 << di); ee++) {	/* Matrix input combination chanel */
					double vv = 0.0;
					jj = ff * (1 << di) + ee;			/* Matrix Parameter index */

					if (dolab) {
						for (f = 0; f < 3; f++) {		/* XYZ channels */
							vv += de_dout[0][f] * dout_lab[f][ff]
							    * dout_mato[ff] * dmato_mv[ee];
						}
					} else {
						vv += de_dout[0][ff] * dout_mato[ff] * dmato_mv[ee];
					}
					dav[p->mat_off + jj] += vv;
				}
			}
		}

		if (p->opt_msk & oc_o) {
			/* Output transfer parameters */
			for (ff = 0; ff < fdi; ff++) {				/* Parameter output chanel */
				for (k = 0; k < p->oluord[ff]; k++) {	/* Param within channel */
					double vv = 0.0;
					jj = p->out_offs[ff] - p->out_off + k;	/* Overall output trans param */

					if (dolab) {
						for (f = 0; f < 3; f++) {		/* Lab channels */
							vv += de_dout[0][f] * dout_lab[f][ff] * dout_ov[jj];
						}
					} else {
						vv += de_dout[0][ff] * dout_ov[jj];
					}
					dav[p->out_off + jj] += vv;
				}
			}
		}
	}

	/* Normalise error to be an average delta E squared */
	rv /= (double)p->nodp;
	for (i = 0; i < p->tot_cnt; i++)
		dav[i] /= (double)p->nodp;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	rv += smv = dshapmag(p, dav);

	/* Copy the del for parameters being optimised to return array */
	for (i = 0; i < p->opt_cnt; i++)
		dv[i] = dav[p->opt_off + i];

#ifdef DEBUG
printf("~1(%f)dluoptfunc returning %f\n",smv,rv);
#endif

	if (p->verb)
		printf("."), fflush(stdout);
	return rv;
}

#ifdef NEVER
/* Check partial derivative function within luoptfunc() */

static double luoptfunc(void *edata, double *v) {
	luopt *p = (luopt *)edata;
	int i;
	double dv[MXPARMS];
	double rv, drv;
	double trv;
	
	rv = luoptfunc_(edata, v);
	drv = dluoptfunc(edata, dv, v);

	if (fabs(rv - drv) > 1e-6)
		printf("######## RV MISMATCH is %f should be %f ########\n",rv,drv);

	/* Check each parameter delta */
	for (i = 0; i < p->opt_cnt; i++) {
		double del;

		v[i] += 1e-7;
		trv = luoptfunc_(edata, v);
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

/* Output curve symetry optimisation function handed to powell() */
/* Just the order 0 value will be adjusted */
static double symoptfunc(void *edata, double *v) {
	luopt *p = (luopt *)edata;
	double out[1], in[1] = { 0.0 };
	double parms[MXLUORD];
	int ch = p->symch;		/* Output channel being adjusted for symetry */
	double rv;
	int i;

	/* Copy the parameter being tested back into luopt */
	p->v[p->out_offs[ch]] = v[0];
	*out = icxSTransFunc(p->v + p->out_offs[ch], p->oluord[ch], *in,
							   p->out_min[ch], p->out_max[ch]);

	rv = out[0] * out[0];

#ifdef DEBUG
printf("~1symoptfunc returning %f\n",rv);
#endif
	return rv;
}

/* -------------------------------------------- */

/* Context for figuring input and output curves */
typedef struct {
	int iix;
	int oix;
	luopt *os;		/* Optimisation structure */
} curvectx;

/* Function to pass to rspl to set input transfer function */
static void
set_linfunc(
	void *cc,			/* curvectx structure */
	double *out,		/* Device output value */
	double *in			/* Device input value */
) {
	curvectx *c = (curvectx *)cc;		/* this */
	luopt *p = c->os;

	if (c->iix >= 0) {				/* Input curve */
		*out = icxSTransFunc(p->v + p->in_offs[c->iix], p->iluord[c->iix], *in,
		                     p->in_min[c->iix], p->in_max[c->iix]);
	} else if (c->oix >= 0) {		/* Output curve */
		*out = icxSTransFunc(p->v + p->out_offs[c->oix], p->oluord[c->oix], *in,
							   p->out_min[c->oix], p->out_max[c->oix]);
	}
}

/* Function to pass to rspl to set inverse input transfer function */
static void
icxLuLut_invinput_func(
	void *pp,			/* icxLuLut */
	double *out,		/* Device output value */
	double *in			/* Device input value */
) {
	icxLuLut *p      = (icxLuLut *)pp;			/* this */
	double tin[MAX_CHAN];
	double tout[MAX_CHAN];
	int i;

	for (i = 0; i < p->inputChan; i++)
		tin[i] = 0.0;
	tin[p->iol_ch] = in[0];
	p->inv_input(p, tout, tin);				/* Use rspl inverse just setup */

	out[0] = tout[p->iol_ch];
}


/* Functions to pass to icc settables() to setup icc Lut */
/* Input table */
static void set_input(void *cntx, double *out, double *in) {
	icxLuLut *p = (icxLuLut *)cntx;

	if (p->noiluts != 0) {	/* Input table must be linear */
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
/* Context for making clut relative to white and black points */
typedef struct {
	icxLuLut *x;		/* Object being created */
	icmHeader *h;		/* Pointer to icc header */
	double mat[3][3];	/* XYZ White point aprox relative to accurate relative transform matrix */
} relativectx;

/* Function to pass to rspl to re-set output values, */
/* to make them relative to the white and black points */
static void
reset_relative_func(
	void *pp,			/* relativectx structure */
	double *out,		/* output value */
	double *in			/* input value */
) {
	int f;
	relativectx *p    = (relativectx *)pp;
	double tt[3];

	/* Convert the clut values to output values */
	p->x->output(p->x, tt, out);	

	if (p->x->pcs == icSigLabData) {	/* We are in Lab */

		icmLab2XYZ(&icmD50, tt, tt);

		/* Convert from aproximate to accurate Relative colorimetric */
		out[0] = p->mat[0][0] * tt[0] + p->mat[0][1] * tt[1] + p->mat[0][2] * tt[2];
		out[1] = p->mat[1][0] * tt[0] + p->mat[1][1] * tt[1] + p->mat[1][2] * tt[2];
		out[2] = p->mat[2][0] * tt[0] + p->mat[2][1] * tt[1] + p->mat[2][2] * tt[2];

		icmXYZ2Lab(&icmD50, out, out);

	} else {	/* We are all in XYZ */

		/* Convert from aproximate to accurate Relative colorimetric */
		out[0] = p->mat[0][0] * tt[0] + p->mat[0][1] * tt[1] + p->mat[0][2] * tt[2];
		out[1] = p->mat[1][0] * tt[0] + p->mat[1][1] * tt[1] + p->mat[1][2] * tt[2];
		out[2] = p->mat[2][0] * tt[0] + p->mat[2][1] * tt[1] + p->mat[2][2] * tt[2];
	}

	/* And then convert them back to clut values */
	p->x->inv_output(p->x, out, out);	
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Routine to figure out a suitable black point for CMYK */

/* Structure to hold optimisation information */
typedef struct {
	icxLuLut *p;			/* Object being created */
	double fromAbs[3][3];	/* From abs to aprox relative */
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
	int j;

	/* See if over ink limit or outside device range */
	ovr = icxLimit((void *)b->p, pv);		/* > 0.0 if outside device gamut or ink limit */

	/* Compute the asbolute Lab value: */
	b->p->input(b->p, bcc.p, pv);						/* Through input tables */
	b->p->clutTable->interp(b->p->clutTable, &bcc);	/* Through clut */
	b->p->output(b->p, bcc.v, bcc.v);					/* Through the output tables */

	if (b->p->pcs != icSigXYZData) 	/* Convert PCS to XYZ */
		icmLab2XYZ(&icmD50, bcc.v, bcc.v);

	/* Convert from aproximate relative to Absolute colorimetric */
	tt[0] = b->toAbs[0][0] * bcc.v[0] + b->toAbs[0][1] * bcc.v[1]
	                                  + b->toAbs[0][2] * bcc.v[2];
	tt[1] = b->toAbs[1][0] * bcc.v[0] + b->toAbs[1][1] * bcc.v[1]
	                                  + b->toAbs[1][2] * bcc.v[2];
	tt[2] = b->toAbs[2][0] * bcc.v[0] + b->toAbs[2][1] * bcc.v[1]
		                              + b->toAbs[2][2] * bcc.v[2];
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
co                 *ipoints,		/* Array of input points */
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
	co *points;				/* Working copy of ipoints */
	double apxwp[3];		/* Aproximate XYZ White point */
	double fromAbs[3][3];	/* From abs to aprox relative */
	double toAbs[3][3];		/* To abs from aprox relative */
	int rv;

#ifdef DEBUG_PLOT
	#define	XRES 100
	double xx[XRES];
	double y1[XRES];
#endif /* DEBUG */

	if (flags & ICX_VERBOSE)
		rsplflags |= RSPL_VERBOSE;

	if (flags & ICX_EXTRA_FIT)
		rsplflags |= RSPL_EXTRAFIT;		/* Try extra hard to fit data */

	luflags = flags;		/* Transfer straight though ? */

	/* Do basic creation and initialisation */
	if ((p = alloc_icxLuLut(xicp, plu, luflags)) == NULL)
		return NULL;

	/* Set LuLut "create" specific flags: */
	if (flags & ICX_NO_IN_LUTS)
		p->noiluts = 1;

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

	/* ------------------------------- */

	if (flags & ICX_VERBOSE)
		printf("Estimating white point\n");

	if (flags & ICX_SET_WHITE) {
		rspl *wpest;		/* Temporary device -> CIE interpolation */
		int gres[MXDI];
		co cc;

		/* Allocate a temporary spline interpolation structure */
		if ((wpest = new_rspl(p->inputChan, p->outputChan)) == NULL) {
			p->pp->errc = 2;
			sprintf(p->pp->err,"Creation of temporary rspl failed");
			p->del((icxLuBase *)p);
			return NULL;
		}

		for (e = 0; e < p->inputChan; e++)
			gres[e] = p->lut->clutPoints > 7 ? 7 : p->lut->clutPoints;

		/* Create the temporary interpolation */
		wpest->fit_rspl(wpest,
	               0,
		           ipoints, nodp,
		           p->ninmin, p->ninmax,
	               gres,
		           p->noutmin, p->noutmax, smooth, avgdev);


		/* Figure out the device values for white */
		if (h->deviceClass == icSigInputClass) {

			/* We assume that the input target is well behaved, */
			/* and that it includes a white and black point patch, */
			/* and that they have the extreme L/Y values */
	
			cc.v[pcsy] = -1e60;
	
			/* Discover the white points */
			for (i = 0; i < nodp; i++) {
				if (ipoints[i].v[pcsy] > cc.v[pcsy]) {
					cc.v[0] = ipoints[i].v[0];
					cc.v[1] = ipoints[i].v[1];
					cc.v[2] = ipoints[i].v[2];
					for (e = 0; e < p->inputChan; e++)
						cc.p[e] = ipoints[i].p[e];
				}
			}

		} else {	/* Output or Display device */

			switch(h->colorSpace) {
	
				case icSigCmykData:
				case icSigCmyData:
					for (e = 0; e < p->inputChan; e++)
						cc.p[e] = 0.0;
					break;
				case icSigRgbData:
					for (e = 0; e < p->inputChan; e++)
						cc.p[e] = 1.0;
					break;
	
				case icSigGrayData: {	/* Could be additive or subtractive */
					/* Use heuristic guess */
					if (ipoints[0].v[pcsy] > ipoints[nodp-1].v[pcsy]) {
						for (e = 0; e < p->inputChan; e++)
							cc.p[e] = 0.0;	/* Subractive */
					} else {
						for (e = 0; e < p->inputChan; e++)
							cc.p[e] = 1.0;		/* Additive */
					}
					break;
				}

				default:
					xicp->errc = 1;
					sprintf(xicp->err,"set_icxLuLut: can't handle color space %s",
					                           icm2str(icmColorSpaceSignature, h->colorSpace));
					wpest->del(wpest);
					p->del((icxLuBase *)p);
					return NULL;
					break;
			}
		}

		wpest->interp(wpest, &cc);

		for (f = 0; f < p->outputChan; f++)
			apxwp[f] = cc.v[f];

		if (p->pcs != icSigXYZData) 	/* Convert white point to XYZ */
			icmLab2XYZ(&icmD50, apxwp, apxwp);

		if (flags & ICX_VERBOSE) {
			double apxlwp[3];
			icmXYZ2Lab(&icmD50, apxlwp, apxwp);
			printf("Approximate White point XYZ = %f %f %f, Lab = %f %f %f\n",
			        apxwp[0],apxwp[1],apxwp[2],apxlwp[0],apxlwp[1],apxlwp[2]);
		}

		wpest->del(wpest);		/* Done with trial white point interpolation */
	} else {

		icmXYZ2Ary(apxwp, icmD50);		/* Set a default value - D50 */
	}

	/* Setup data points that are adjusted for the approximate white point */
	{
		double tt[3];
		icmXYZNumber swp;

		/* Allocate the array passed to fit_rspl() */
		if ((points = malloc(sizeof(co) * nodp)) == NULL) {
			p->pp->errc = 2;
			sprintf(p->pp->err,"Allocation of scattered coordinate array failed");
			p->del((icxLuBase *)p);
			return NULL;
		}

		icmAry2XYZ(swp, apxwp);

		/* Absolute->Aprox. Relative Adaptation matrix */
		icmChromAdaptMatrix(ICM_CAM_BRADFORD, icmD50, swp, fromAbs);

		/* Aproximate relative to absolute conversion matrix */
		icmChromAdaptMatrix(ICM_CAM_BRADFORD, swp, icmD50, toAbs);

		/* Setup transformed points for the device */
		for (i = 0; i < nodp; i++) {
			for (e = 0; e < p->inputChan; e++)
				points[i].p[e] = ipoints[i].p[e];
			for (f = 0; f < p->outputChan; f++)
				tt[f] = ipoints[i].v[f];

			if (p->pcs == icSigLabData)	/* Convert to XYZ for chromatic shift */
				icmLab2XYZ(&icmD50, tt, tt);

			/* Convert from Absolute to aproximate Relative colorimetric */
			points[i].v[0] = fromAbs[0][0] * tt[0] + fromAbs[0][1] * tt[1] + fromAbs[0][2] * tt[2];
			points[i].v[1] = fromAbs[1][0] * tt[0] + fromAbs[1][1] * tt[1] + fromAbs[1][2] * tt[2];
			points[i].v[2] = fromAbs[2][0] * tt[0] + fromAbs[2][1] * tt[1] + fromAbs[2][2] * tt[2];

			/* Convert to Lab for optimisation comparison */
			icmXYZ2Lab(&icmD50, points[i].v, points[i].v);
		}
	}

	if (h->colorSpace == icSigGrayData) {
		p->noiluts = p->nooluts = 1;	/* Don't use device or PCS curves for monochrome */
	}

	if (flags & ICX_VERBOSE)
		printf("Creating optimised input and output curves\n");

	{
		luopt os;
		double v[MXPARMS];		/* Working parameters */
		double sa[MXPARMS];		/* Search area */
		int iluord[MXDI];		/* Input curve orders */
		int oluord[MXDO];		/* Output curve orders */

		/* Setup for optimising run */
		if (flags & ICX_VERBOSE)
			os.verb = 1;
		else
			os.verb = 0;
		os.points = points;
		os.nodp   = nodp;

		{
			int luord;

			if (quality >= 3) {				/* Ultra high */
				luord = 10;			
			} else if (quality == 2) {		/* High */
				luord = 8;			
			} else if (quality == 1) {		/* Medium */
				luord = 6;			
			} else {						/* Low */
				luord = 4;			
			}
#ifdef HACK_ILORD
luord = HACK_ILORD;			
printf("HACK - luord set to %d!\n",luord);
#endif

			/* Set the curve order for input (device) */
			for (e = 0; e < p->inputChan; e++) {
				iluord[e] = luord;
			}

			if (quality >= 3) {				/* Ultra high */
				luord = 6;			
			} else if (quality == 2) {		/* High */
				luord = 5;			
			} else if (quality == 1) {		/* Medium */
				luord = 4;			
			} else {						/* Low */
				luord = 3;			
			}
#ifdef HACK_OLORD
luord = HACK_OLORD;			
printf("HACK - luord set to %d!\n",luord);
#endif

			/* Set curve order for output (PCS) */
			for (f = 0; f < p->outputChan; f++) {
				oluord[f] = luord;
			}
		}

		init_luopt(&os, p, iluord, oluord);
	
		if (p->noiluts == 0 || p->nooluts == 0) {	/* No point otherwise */

			if (flags & ICX_VERBOSE)
				printf("About to optimise temporary matrix\n");

#ifdef DEBUG
printf("\nBefore matrix opt:\n");
dump_luopt(&os);
#endif
			/* Optimise matrix on its own */
			os.opt_msk = oc_m;
			setup_luopt(&os, v, sa, 0.0, 0.5); 
#ifdef NODDV
			if (powell(os.opt_cnt, v, sa, 0.05, 1000, luoptfunc, (void *)&os) < 0.0)
				error ("set_icxLuLut: Powell failed");
#else
			if (conjgrad(os.opt_cnt, v, sa, 0.05, 1000, luoptfunc, dluoptfunc, (void *)&os) < 0.0)
				error ("set_icxLuLut: Conjgrad failed");
#endif
			for (i = 0; i < os.opt_cnt; i++)		/* Copy optimised values back */
				os.v[os.opt_off + i] = v[i];

#ifdef DEBUG
printf("\nAfter matrix opt:\n");
dump_luopt(&os);
#endif
		}

		/* Optimise input and matrix together */
		if(p->noiluts == 0) {

			if (flags & ICX_VERBOSE)
				printf("\nAbout to optimise input curves and matrix\n");

			os.opt_msk = oc_im;
			setup_luopt(&os, v, sa, 0.5, 0.3); 
#ifdef NODDV
			if (powell(os.opt_cnt, v, sa, 0.01, 1000, luoptfunc, (void *)&os) < 0.0)
				error ("set_icxLuLut: Powell failed");
#else
			if (conjgrad(os.opt_cnt, v, sa, 0.01, 1000, luoptfunc, dluoptfunc, (void *)&os) < 0.0)
				error ("set_icxLuLut: Conjgrad failed");
#endif
			for (i = 0; i < os.opt_cnt; i++)		/* Copy optimised values back */
				os.v[os.opt_off + i] = v[i];
#ifdef DEBUG
printf("\nAfter input and matrix opt:\n");
dump_luopt(&os);
#endif
		}

		/* Optimise the matrix and output curves together */
		if(p->nooluts == 0) {

			if (flags & ICX_VERBOSE)
				printf("\nAbout to optimise output curves and matrix\n");

			os.opt_msk = oc_mo;
			setup_luopt(&os, v, sa, 0.3, 0.3); 
#ifdef NODDV
			if (powell(os.opt_cnt, v, sa, 0.005, 1000, luoptfunc, (void *)&os) < 0.0)
				error ("set_icxLuLut: Powell failed");
#else
			if (conjgrad(os.opt_cnt, v, sa, 0.005, 1000, luoptfunc, dluoptfunc, (void *)&os) < 0.0)
				error ("set_icxLuLut: Conjgrad failed");
#endif
			for (i = 0; i < os.opt_cnt; i++)		/* Copy optimised values back */
				os.v[os.opt_off + i] = v[i];
//~2
#ifdef DEBUG
printf("\nAfter output opt:\n");
dump_luopt(&os);
#endif
	
	
			/* Optimise input and matrix together again, after altering matrix */
			if (p->noiluts == 0) {

				if (flags & ICX_VERBOSE)
					printf("\nAbout to optimise input curves and matrix again\n");

				os.opt_msk = oc_im;
				setup_luopt(&os, v, sa, 0.2, 0.2); 
#ifdef NODDV
				if (powell(os.opt_cnt, v, sa, 0.002, 2000, luoptfunc, (void *)&os) < 0.0)
					error ("set_icxLuLut: Powell failed");
#else
				if (conjgrad(os.opt_cnt, v, sa, 0.005, 1000, luoptfunc, dluoptfunc, (void *)&os) < 0.0)
					error ("set_icxLuLut: Conjgrad failed");
#endif
				for (i = 0; i < os.opt_cnt; i++)		/* Copy optimised values back */
					os.v[os.opt_off + i] = v[i];
			}
//~2
#ifdef DEBUG
printf("\nAfter 2nd input and matrix opt:\n");
dump_luopt(&os);
#endif

#ifndef NODDV
			/* Optimise all together */
			if (p->noiluts == 0) {

				if (flags & ICX_VERBOSE)
					printf("\nAbout to optimise input, matrix and output together\n");

				os.opt_msk = oc_imo;
				setup_luopt(&os, v, sa, 0.1, 0.1); 
				if (conjgrad(os.opt_cnt, v, sa, 0.001, 1000, luoptfunc, dluoptfunc, (void *)&os) < 0.0)
					error ("set_icxLuLut: Conjgrad failed");
				for (i = 0; i < os.opt_cnt; i++)		/* Copy optimised values back */
					os.v[os.opt_off + i] = v[i];
			}

#ifdef DEBUG
printf("\nAfter all together opt:\n");
dump_luopt(&os);
#endif

#endif /* !NODDV */

			/* Adjust output curve white point */
			if (p->pcs == icSigLabData) {

				if (flags & ICX_VERBOSE)
					printf("\nAbout to adjust a and b output curves for white point\n");

				for (f = 1; f < 3; f++) {
					os.symch = f;
					v[0] = os.v[os.out_offs[f]];	/* Current parameter value */
					sa[0] = 0.1;					/* Search radius */
					if (powell(1, v, sa, 0.0000001, 500, symoptfunc, (void *)&os) < 0.0)
						error ("set_icxLuLut: Powell failed");
					os.v[os.out_offs[f]] = v[0];	/* Copy results back */
				}
			}
		}

		if (os.verb)
			printf("\n");

		/* - - - - - - - - - - - - - - - */
		/* Set the xicc input curve rspl */
		for (e = 0; e < p->inputChan; e++) {
			curvectx cx;
	
			cx.os = &os;
			cx.oix = -1;
			cx.iix = e;

			if ((p->inputTable[e] = new_rspl(1, 1)) == NULL) {
				p->pp->errc = 2;
				sprintf(p->pp->err,"Creation of input table rspl failed");
				p->del((icxLuBase *)p);
				return NULL;
			}

			p->inputTable[e]->set_rspl(p->inputTable[e], 0,
			           (void *)&cx, set_linfunc,
    			       &p->ninmin[e], &p->ninmax[e],
			           &p->lut->inputEnt,
			           &p->ninmin[e], &p->ninmax[e]);
#ifdef DEBUG_PLOT
			/* Display the result scattered fit */
			for (i = 0; i < XRES; i++) {
				double x;
				co c;
				x = i/(double)(XRES-1);
				xx[i] = (x * p->ninmax[e] - p->ninmin[e]) + p->ninmin[e];
				c.p[0] = xx[i];
				p->inputTable[e]->interp(p->inputTable[e], &c);
				y1[i] = c.v[0];
				}
			do_plot(xx,y1,NULL,NULL,XRES);
#endif /* DEBUG_PLOT */
		}

		/* - - - - - - - - - - - - - - - */
		/* Set the xicc output curve rspl */

		/* Allow for a bigger than normal input and output range, to */
		/* give some leaway in accounting for approximate white point shifted */
		/* profile creation. */
		for (f = 0; f < p->outputChan; f++) {
			double min[1], max[1], exval;
			int entries;
			curvectx cx;

			cx.os = &os;
			cx.iix = -1;
			cx.oix = f;

			/* Expand in and out range by 1.05 */
			exval = (p->noutmax[f] - p->noutmin[f]);
			min[0] = p->noutmin[f] - exval * 0.05 * 0.5;
			max[0] = p->noutmax[f] + exval * 0.05 * 0.5;
  	      	entries = (int)(1.05 * (double)p->lut->outputEnt + 0.5);

			if ((p->outputTable[f] = new_rspl(1, 1)) == NULL) {
				p->pp->errc = 2;
				sprintf(p->pp->err,"Creation of output table rspl failed");
				p->del((icxLuBase *)p);
				return NULL;
			}

			p->outputTable[f]->set_rspl(p->outputTable[f], 0,
		           (void *)&cx, set_linfunc,
					min, max, &entries, min, max);

#ifdef DEBUG_PLOT
			/* Display the result scattered fit */
			for (i = 0; i < XRES; i++) {
				double x;
				co c;
				x = i/(double)(XRES-1);
				xx[i] = (x * p->noutmax[f] - p->noutmin[f]) + p->noutmin[f];
				c.p[0] = xx[i];
				p->outputTable[f]->interp(p->outputTable[f], &c);
				y1[i] = c.v[0];
				}
			do_plot(xx,y1,NULL,NULL,XRES);
#endif /* DEBUG_PLOT */

		}
	}

	if (flags & ICX_VERBOSE)
		printf("Creating fast inverse input lookups\n");

	/* Setup center clip target for input inversion */
	for (i = 0; i < p->inputChan; i++) {
		p->inputClipc[i] = (p->ninmin[i] + p->ninmax[i])/2.0;
	}

	/* Create rspl based reverse input lookups used in ink limit function. */
	for (e = 0; e < p->inputChan; e++) {
		int gres = 256;

		if ((p->revinputTable[e] = new_rspl(1, 1)) == NULL) {
			p->pp->errc = 2;
			sprintf(p->pp->err,"Creation of reverse input table rspl failed");
			p->del((icxLuBase *)p);
			return NULL;
		}
		p->iol_out = 2;		/* Input lookup */
		p->iol_ch = e;		/* Chanel */
		p->revinputTable[e]->set_rspl(p->revinputTable[e], 0,
		           (void *)p, icxLuLut_invinput_func,
		           &p->ninmin[e], &p->ninmax[e], &gres, &p->ninmin[e], &p->ninmax[e]);
	}


	if (flags & ICX_VERBOSE)
		printf("Compensate scattered data for input curves\n");

	/* Setup center clip target for output inversion */
	for (i = 0; i < p->outputChan; i++) {
		p->outputClipc[i] = (p->noutmin[i] + p->noutmax[i])/2.0;
	}

	if (flags & ICX_VERBOSE)
		printf("Compensate scattered data for output curve\n");

	/* Setup transformed points for the device */
	for (i = 0; i < nodp; i++) {
		double tt[3];
		co cc;
		int nsoln;
		double cdir;

		for (e = 0; e < p->inputChan; e++)
			points[i].p[e] = ipoints[i].p[e];
		for (f = 0; f < p->outputChan; f++)
			tt[f] = ipoints[i].v[f];

		if (p->pcs == icSigLabData)	/* Convert to XYZ for chromatic shift */
			icmLab2XYZ(&icmD50, tt, tt);

		/* Convert from Absolute to aproximate Relative colorimetric */
		points[i].v[0] = fromAbs[0][0] * tt[0] + fromAbs[0][1] * tt[1] + fromAbs[0][2] * tt[2];
		points[i].v[1] = fromAbs[1][0] * tt[0] + fromAbs[1][1] * tt[1] + fromAbs[1][2] * tt[2];
		points[i].v[2] = fromAbs[2][0] * tt[0] + fromAbs[2][1] * tt[1] + fromAbs[2][2] * tt[2];

		if (p->pcs == icSigLabData)	/* Convert back to Lab */
			icmXYZ2Lab(&icmD50, points[i].v, points[i].v);

		/* Input values forward though input curve */
		for (e = 0; e < p->inputChan; e++) {
			cc.p[0] = points[i].p[e];
			p->inputTable[e]->interp(p->inputTable[e], &cc);
			points[i].p[e] = cc.v[0];
		}

		/* Output values backwards though output curve */
		for (f = 0; f < p->outputChan; f++) {
			cc.v[0] = points[i].v[f];

			/* Clip towards center of output range */
			cdir = (p->noutmin[f] + p->noutmax[f]) * 0.5 - cc.v[0];

			nsoln = p->outputTable[f]->rev_interp (
				p->outputTable[f],	/* this */
				0,					/* No flags */
				1,					/* Maximum number of solutions allowed for */
				NULL, 				/* No auxiliary input targets */
				&cdir,				/* Clip vector direction and length */
				&cc);				/* Input and output values */
			nsoln &= RSPL_NOSOLNS;	/* Get number of solutions */
			if (nsoln == 0)
				error("set_icxLuLut: Unexpected failure to find reverse solution for output linearisation");

			points[i].v[f] = cc.p[0];
		}
	}

	/* ------------------------------- */
	{
		int gres[MXDI];

		if (flags & ICX_VERBOSE)
			printf("Create clut from scattered data\n");

		if ((p->clutTable = new_rspl(p->inputChan, p->outputChan)) == NULL) {
			p->pp->errc = 2;
			sprintf(p->pp->err,"Creation of clut table rspl failed");
			p->del((icxLuBase *)p);
			return NULL;
		}

		for (e = 0; e < p->inputChan; e++)
			gres[e] = p->lut->clutPoints;

		/* Initialise from scattered data */
		/* Return non-zero if result is non-monotonic */
		/* ~~~~ should warn if non-monotonic ??? ~~~~ */
		p->clutTable->fit_rspl(
			p->clutTable,	/* this */
			rsplflags,		/* Combination of flags */
			points,			/* Array holding position and function values of data points */
			nodp,			/* Number of data points */
			p->ninmin,		/* Grid low scale - will expand to enclose data, NULL = default 0.0 */
			p->ninmax,		/* Grid high scale - will expand to enclose data, NULL = default 1.0 */
			gres,			/* Spline grid resolution, ncells = gres-1 */
			p->noutmin,		/* Data value low normalize, NULL = default 0.0 */
			p->noutmax,		/* Data value high normalize - NULL = default 1.0 */
			smooth,			/* Smoothing factor, nominal = 1.0 */
		    avgdev			/* reading Average Deviation as a proportion of the input range */
		);
	}

#ifdef DEBUG
	/* Sanity check the rspl result */

#endif	// DEBUG

	/* Setup all the clipping, ink limiting and auxiliary stuff, */
	/* in case a reverse call is used. Need to avoid relying on inking */
	/* stuff that makes use of the white/black points, since they haven't */
	/* beem set up yet. */
	if (setup_ink_icxLuLut(p, ink, 0) != 0) {
		p->del((icxLuBase *)p);
		return NULL;
	}

	/* Deal with white/black points */
	if (flags & (ICX_SET_WHITE | ICX_SET_BLACK)) {
		double dwhite[MXDI], dblack[MXDI];	/* Device white and black values */
		double rwp[3];		/* Relative White point in XYZ */
		double wp[3];		/* Absolute White point in XYZ */
		double bp[3];		/* Absolute Black point in XYZ */

		if (flags & ICX_VERBOSE)
			printf("\nFind white & black points\n");

		icmXYZ2Ary(wp, icmD50);		/* Set a default value - D50 */
		icmXYZ2Ary(bp, icmBlack);	/* Set a default value - absolute black */

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
		} else {

			/* Assume Output or Monitor class */
			switch(h->colorSpace) {

				case icSigCmykData:
				case icSigCmyData:
					for (e = 0; e < p->inputChan; e++) {
						dwhite[e] = 0.0;
						dblack[e] = 1.0;
					}
					break;

				case icSigRgbData:

					for (e = 0; e < p->inputChan; e++) {
						dwhite[e] = 1.0;
						dblack[e] = 0.0;
					}

					break;

				case icSigGrayData: {	/* Could be additive or subtractive */
					/* Use heuristic guess */
					if (ipoints[0].v[pcsy] > ipoints[nodp-1].v[pcsy]) {
						for (e = 0; e < p->inputChan; e++) {
							dwhite[e] = 0.0;	/* Subractive */
							dblack[e] = 1.0;
						}
					} else {
						for (e = 0; e < p->inputChan; e++) {
							dwhite[e] = 1.0;		/* Additive */
							dblack[e] = 0.0;
						}
					}
					break;
				}
			}
		}

		if (flags & ICX_SET_WHITE	 /* White point to find */
		 || ((flags & ICX_SET_BLACK) && (h->colorSpace == icSigCmykData))) {
			co wcc;		/* White point to lookup */

			for (e = 0; e < p->inputChan; e++)
				wcc.p[e] = dwhite[e];

			/* Look this up through the input tables */
			p->input(p, wcc.p, wcc.p);

			p->clutTable->interp(p->clutTable, &wcc);

			/* Look this up through the output tables */
			p->output(p, wcc.v, wcc.v);	
					
			if (p->pcs != icSigXYZData) 	/* Convert PCS to XYZ */
				icmLab2XYZ(&icmD50, wcc.v, wcc.v);

			rwp[0] = wcc.v[0];	/* Remember the relative XYZ white point */
			rwp[1] = wcc.v[1];
			rwp[2] = wcc.v[2];

			/* Convert from aproximate relative to Absolute colorimetric */
			wp[0] = toAbs[0][0] * wcc.v[0] + toAbs[0][1] * wcc.v[1]
			                               + toAbs[0][2] * wcc.v[2];
			wp[1] = toAbs[1][0] * wcc.v[0] + toAbs[1][1] * wcc.v[1]
			                               + toAbs[1][2] * wcc.v[2];
			wp[2] = toAbs[2][0] * wcc.v[0] + toAbs[2][1] * wcc.v[1]
			                               + toAbs[2][2] * wcc.v[2];

			if (flags & ICX_VERBOSE) {
				double labwp[3];
				icmXYZ2Lab(&icmD50, labwp, wp);
				printf("White point XYZ = %f %f %f, Lab = %f %f %f\n",
				        wp[0],wp[1],wp[2],labwp[0],labwp[1],labwp[2]);
			}

		}

		if (flags & ICX_SET_BLACK) { /* Black Point Tag: */
			co bcc;

			/* For CMYK devices, we choose a black point that is in */
			/* the same Lab vector direction as K, with the minimum L value. */
			if (h->deviceClass != icSigInputClass
			 && h->colorSpace == icSigCmykData) {
				bfinds bfs;					/* Callback context */
				double sr[MPP_MXINKS];		/* search radius */
				double tt[MXDO];			/* Temporary */

				/* Setup callback function context */
				bfs.p = p;

				for (i = 0; i < 3; i++) {
					for (j = 0; j < 3; j++) {
						bfs.fromAbs[i][j] = fromAbs[i][j];
						bfs.toAbs[i][j]   = toAbs[i][j];
					}
				}

				/* Lookup abs Lab value of white point */
				icmXYZ2Lab(&icmD50, bfs.p1, wp);

				/* Now figure abs Lab value of K only */
				for (e = 0; e < p->inputChan; e++)
					bcc.p[e] = 0.0;
				bcc.p[3] = 1.0;

				p->input(p, bcc.p, bcc.p);			/* Through input tables */
				p->clutTable->interp(p->clutTable, &bcc);	/* Through clut */
				p->output(p, bcc.v, bcc.v);		/* Through the output tables */

				if (p->pcs != icSigXYZData) 	/* Convert PCS to XYZ */
					icmLab2XYZ(&icmD50, bcc.v, bcc.v);

				/* Convert from aproximate relative to Absolute colorimetric */
				tt[0] = toAbs[0][0] * bcc.v[0] + toAbs[0][1] * bcc.v[1]
				                               + toAbs[0][2] * bcc.v[2];
				tt[1] = toAbs[1][0] * bcc.v[0] + toAbs[1][1] * bcc.v[1]
				                               + toAbs[1][2] * bcc.v[2];
				tt[2] = toAbs[2][0] * bcc.v[0] + toAbs[2][1] * bcc.v[1]
				                               + toAbs[2][2] * bcc.v[2];

				icmXYZ2Lab(&icmD50, bfs.p2, tt); /* Convert K only black point to Lab */
 
				if (flags & ICX_VERBOSE)
					printf("K only value (Lab) = %f %f %f\n",bfs.p2[0], bfs.p2[1], bfs.p2[2]);

				/* Find the device black point */
				for (j = 0; j < p->inputChan; j++) { 
					bcc.p[j] = 0.5;
					sr[j] = 0.05;
				}
				if (powell(p->inputChan, bcc.p, sr, 0.0001, 1000, bfindfunc, (void *)&bfs) < 0.0)
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

			/* Convert from aproximate relative to Absolute colorimetric */
			bp[0] = toAbs[0][0] * bcc.v[0] + toAbs[0][1] * bcc.v[1]
			                               + toAbs[0][2] * bcc.v[2];
			bp[1] = toAbs[1][0] * bcc.v[0] + toAbs[1][1] * bcc.v[1]
			                               + toAbs[1][2] * bcc.v[2];
			bp[2] = toAbs[2][0] * bcc.v[0] + toAbs[2][1] * bcc.v[1]
			                               + toAbs[2][2] * bcc.v[2];

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
				p->del((icxLuBase *)p);
				return NULL;
			}
			if (wo->ttype != icSigXYZArrayType) {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_white_black: white tag has wrong type");
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
				p->del((icxLuBase *)p);
				return NULL;
				}
			if (wo->ttype != icSigXYZArrayType) {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_white_black: black tag has wrong type");
				p->del((icxLuBase *)p);
				return NULL;
			}

			wo->size = 1;
			wo->allocate((icmBase *)wo);	/* Allocate space */
			wo->data[0].X = bp[0];
			wo->data[0].Y = bp[1];
			wo->data[0].Z = bp[2];
		}

		/* Now fixup the clut to make this Lut exactly relative */
		if (flags & ICX_SET_WHITE) {
			icmXYZNumber nrwp;
			relativectx cx;

			if (flags & ICX_VERBOSE)
				printf("Fixup clut for white point\n");

			icmAry2XYZ(nrwp, rwp);

			cx.x = p; 		/* Object being created */
			cx.h = h;		/* Pointer to icc header */

			/* Create aprox relative to accurate relative correction matrix */
			icmChromAdaptMatrix(ICM_CAM_BRADFORD, icmD50, nrwp, cx.mat);

			/* Adjust the clut node values to fix white point */
			p->clutTable->re_set_rspl(
				p->clutTable,		/* this */
				0,					/* Combination of flags */
				(void *)&cx,		/* Opaque function context */
				reset_relative_func /* Function to set from */
			);

			/* Verify the white point */
			if (flags & ICX_VERBOSE) {
				double chlwp[3];
				double chwp[3];			/* Check white point */
				double ch2Abs[3][3];	/* Check to abs from relative */
				icmXYZNumber nwp;
				co wcc;					/* White point to lookup */

				icmAry2XYZ(nwp, wp);

				/* Accurate relative to absolute conversion matrix */
				icmChromAdaptMatrix(ICM_CAM_BRADFORD, nwp, icmD50, ch2Abs);
	
				for (e = 0; e < p->inputChan; e++)
					wcc.p[e] = dwhite[e];


				p->input(p, wcc.p, wcc.p);
				p->clutTable->interp(p->clutTable, &wcc);
				p->output(p, wcc.v, wcc.v);	
						
				if (p->pcs != icSigXYZData) 	/* Convert PCS to XYZ */
					icmLab2XYZ(&icmD50, wcc.v, wcc.v);
	
				/* Convert from relative to Absolute colorimetric */
				chwp[0] = ch2Abs[0][0] * wcc.v[0] + ch2Abs[0][1] * wcc.v[1]
				                               + ch2Abs[0][2] * wcc.v[2];
				chwp[1] = ch2Abs[1][0] * wcc.v[0] + ch2Abs[1][1] * wcc.v[1]
				                               + ch2Abs[1][2] * wcc.v[2];
				chwp[2] = ch2Abs[2][0] * wcc.v[0] + ch2Abs[2][1] * wcc.v[1]
				                               + ch2Abs[2][2] * wcc.v[2];

				icmXYZ2Lab(&icmD50, chlwp, chwp);

				printf("Check White point XYZ = %f %f %f, Lab = %f %f %f\n",
				        chwp[0],chwp[1],chwp[2],chlwp[0],chlwp[1],chlwp[2]);
			}
		}

		/* Setup the clut clipping, ink limiting and auxiliary stuff */
		if (setup_ink_icxLuLut(p, ink, 1) != 0) {	/* re_set_rspl will have invalidated */
			p->del((icxLuBase *)p);
			return NULL;
		}
	}

	if (setup_clip_icxLuLut(p) != 0) {
		p->del((icxLuBase *)p);
		return NULL;
	}

	/* ------------------------------- */

	/* Free the coordinate lists */
	free(points);

	/* Use our rspl's to set the icc Lut AtoB table values. */
	/* Use helper function to do the hard work. */
	if (p->lut->set_tables(p->lut, (void *)p,
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
		xicc_enum_viewcond(xicp, &p->vc, -1, 0);	/* Use a default */
	p->cam = new_icxcam(cam_default);
	p->cam->set_view(p->cam, p->vc.Ev, p->vc.Wxyz, p->vc.Yb, p->vc.La, p->vc.Lv, p->vc.Yf, p->vc.Fxyz, 1);
	
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

/* Function to hand to zbrent to find a clut input value at the ink limit */
/* Returns value < 0.0 when within gamut, > 0.0 when out of gamut */
static double icxLimitFind(void *fdata, double tp) {
	int i;
	double in[MAX_CHAN];
	lutgamctx *p = (lutgamctx *)fdata;
	double tt;

	for (i = 0; i < p->x->inputChan; i++) {
		in[i] = tp * p->in[i];		/* Scale given input value */
	}
	
	tt = icxLimit((void *)p->x, in);

	return tt;
}

/* Function to pass to rspl to create gamut boundary from */
/* forward xLut transform grid points. */
static void
lutfwdgam_func(
	void *pp,			/* lutgamctx structure */
	double *out,		/* output value at clut grid point (ie. PCS value) */
	double *in			/* input value at clut grid point (ie. device value) */
) {
	int f;
	lutgamctx *p    = (lutgamctx *)pp;
	double pcso[3];	/* PCS output value */

	/* Figure if we are over the ink limit. */
	if (   p->x->ink.tlimit >= 0.0 && p->x->ink.klimit >= 0.0
	    && icxLimit((void *)p->x, in) > 0.0) {
		int i;
		double sf;

		/* We are, so use the bracket search to discover a scale */
		/* for the clut input value that will put us on the ink limit. */

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
		/* Convert the clut values to PCS output values */
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
	int f;
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
	icColorSpaceSignature pcs;
	icmLookupFunc func;
	icRenderingIntent intent;
	double white[3], black[3];
	gamut *gam;

	/* get some details */
	plu->spaces(plu, NULL, NULL, NULL, NULL, NULL, &intent, &func, &pcs);

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

		cx.g = gam = new_gamut(detail);
		cx.x = luluto;

		luluto->clutTable->scan_rspl(
			luluto->clutTable,	/* this */
			0,					/* Combination of flags */
			(void *)&cx,		/* Opaque function context */
			lutfwdgam_func		/* Function to set from */
		);

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
		}
		if ((cx.flu = p->get_luobj(p, 0, icmFwd, intent, pcs, icmLuOrdNorm,
		                              &plu->vc, NULL)) == NULL) {
			return NULL;	/* oops */
		}

		cx.g = gam = new_gamut(detail);

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






