
/* 
 * nearsmth
 *
 * Gamut mapping support routine that creates a list of
 * guide vectors that map from a source to destination
 * gamut, smoothed to retain reasonably even spacing.
 *
 * Author:  Graeme W. Gill
 * Date:    17/1/2002
 * Version: 1.00
 *
 * Copyright 2002 - 2006 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
	Description:

	We create a set of "guide vectors" that map the source gamut to
	the destination, for use by the gammap code in creating
	a 3D gamut mapping. 

	(See gammap.txt for a more detailed descrition)

 */

/*
 * TTBD:
 *		It might work better if the cusp mapping had separate control
 *		over the L and h degree of map, as well as the L and h effective radius ?
 *		That way, saturation hue distortions with L might be reduced.
 *
 *       Improve error handling.
 *       Merge with gammap ?
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "icc.h"
#include "gamut.h"
#include "numlib.h"
#include "nearsmth.h"

#undef VERB			/* Print information about everything */
#undef SHOW_ELEVATION	/* Use rotated and elevated source points (diagnostic - affects result) */

#undef USE_TRIANG_VERTS	/* Usually not defined - allocation isn't right */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* A set of functions to help handle the weighting configuration */

/* Blend a two groups of individual weights into one, given two weightings */
void near_wblend(
gammapweights *dst,
gammapweights *src1, double wgt1,
gammapweights *src2, double wgt2
) {

#define NSBLEND(xxx) dst->xxx = wgt1 * src1->xxx + wgt2 * src2->xxx

	NSBLEND(a.o);
	NSBLEND(a.w.l);
	NSBLEND(a.w.c);
	NSBLEND(a.w.h);
	NSBLEND(r.o);
	NSBLEND(r.w.l);
	NSBLEND(r.w.c);
	NSBLEND(r.w.h);
	NSBLEND(l.o);
	NSBLEND(l.w.l);
	NSBLEND(l.w.c);
	NSBLEND(l.w.h);
	NSBLEND(cw);
	NSBLEND(cr);
	NSBLEND(e);
#undef NSBLEND
}

/* Expand the compact form of weights into the explicit form. */
/* The explicit form is red, yellow, green, cyan, blue, magenta, default */
void expand_weights(gammapweights out[7], gammapweights *in) {
	int i, j;

	/* mark output so we can recognise havving been set or not */
	for (i = 0; i < 7; i++) {
		out[i].ch = gmm_end;
	}

	/* Assume a default default */
	memset((void *)&out[6], 0, sizeof(gammapweights));
	out[6].ch = gmm_default;

	/* Expand the compact form */
	for (i = 0; in[i].ch != gmm_end; i++) {
		if (in[i].ch >= gmm_red && in[i].ch <= gmm_default) {
			j = in[i].ch - gmm_red;
			memcpy((void *)&out[j], (void *)&in[i], sizeof(gammapweights));
		}
	}

	/* Fill in any hextants that haven't been set */
	for (i = 0; i < 6; i++) {
		if (out[i].ch == gmm_end) {
			memcpy((void *)&out[i], (void *)&out[6], sizeof(gammapweights));
		}
	}
}

/* Blend a two expanded groups of individual weights into one */
void near_xwblend(
gammapweights *dst,
gammapweights *src1, double wgt1,
gammapweights *src2, double wgt2
) {
	int i;
	for (i = 0; i < 7; i++) {
		near_wblend(&dst[i], &src1[i], wgt1, &src2[i], wgt2);
		dst[i].ch = src1[i].ch;
	}
}

/* Given a point location, return the interpolated weighting values at that point. */
void interp_xweights(int isJab, gammapweights *out, double pos[3], gammapweights in[7]) {
	double JCh[3];
	int li, ui;
	double lh, uh;
	double lw, uw;

	if (isJab != 0)
		isJab = 1;

	/* Convert to polar */
	icmLab2LCh(JCh, pos);
	if (JCh[2] < gam_hues[isJab][0])
		JCh[2] += 360.0;

	/* Figure out what hextant we're between */
	for (li = 0; li < 6; li++) {
		ui = li+1;
		lh = gam_hues[isJab][li];
		uh = gam_hues[isJab][ui];
		if (JCh[2] >= lh && JCh[2] <= uh)
			break;
	}

	/* Compute weights */
	uw = (JCh[2] - lh)/(uh - lh);
	uw = uw * uw * (3.0 - 2.0 * uw);	/* Apply spline */
	lw = (1.0 - uw);

	/* Blend weights at the two hues */
	if (ui >= 6)
		ui = 0;
	near_wblend(out, &in[li], lw, &in[ui], uw);

	/* If we're close to the center, blend to the default weight */
	if (JCh[1] < 20.0) {
		lw = (20.0 - JCh[1])/20.0;
		uw = (1.0 - lw);
		near_wblend(out, &in[6], lw, out, uw);
	}
}

/* Compute the weighted delta E squared */
static double wde(
double in1[3],
double in2[3],
double lweight,
double cweight,
double hweight
) {
	double desq, dhsq;
	double dlsq, dcsq;
	double vv;

	/* Compute delta L squared and delta E squared */
	{
		double dl, da, db;
		dl = in1[0] - in2[0];
		dlsq = dl * dl;		/* dl squared */
		da = in1[1] - in2[1];
		db = in1[2] - in2[2];

		desq = dlsq + da * da + db * db;
	}

	/* compute delta chromanance squared */
	{
		double c1, c2, dc;

		/* Compute chromanance for the two colors */
		c1 = sqrt(in1[1] * in1[1] + in1[2] * in1[2]);
		c2 = sqrt(in2[1] * in2[1] + in2[2] * in2[2]);

		dc = c2 - c1;
		dcsq = dc * dc;

		/* ~~ Doing dcsq = sqrt(dcsq); here seemed */
		/* to improve the saturation result. How did this work ?? */
	}

	/* Compute delta hue squared */
	if ((dhsq = desq - dlsq - dcsq) < 0.0)
		dhsq = 0.0;

	vv = lweight * dlsq + cweight * dcsq + hweight * dhsq;

	return fabs(vv);	/* Avoid -0.0 */
}

/* Given the weighting structure and the relevant point locations */
/* return the total weighted error squared. */
static double comperr(
gammapweights *w,	/* weightings */
double dtp[3],		/* Source value closest onto dst surface */
double sv[3],		/* Source value */
double anv[3],		/* Average neighborhood dst surface point */
double drv[3]		/* Destination radial value */
) {
	double va, vr, vl, vv = 0.0;

	/* Absolute */
	va = wde(dtp, sv, w->a.o * w->a.w.l, w->a.o * w->a.w.c, w->a.o * w->a.w.h);

	/* Relative */
	vr = wde(dtp, anv, w->r.o * w->r.w.l, w->r.o * w->r.w.c, w->r.o * w->r.w.h);

//printf("~1 vr = %f from dtp %f %f %f, anv %f %f %f, weights %f %f %f\n", vr, dtp[0], dtp[1], dtp[2], anv[0], anv[1], anv[2], w->r.o * w->r.w.l, w->r.o * w->r.w.c, w->r.o * w->r.w.h);

	/* Radial */
	vl = wde(dtp, drv, w->l.o * w->l.w.l, w->l.o * w->l.w.c, w->l.o * w->l.w.h);

//printf("~1 abse %f, rele %f, rade %f\n",va,vr,vl);
//printf("~1 abse %f, rele %f, rade %f\n",sqrt(va),sqrt(vr),sqrt(vl));
//printf("~1 abse %f, rele %f, rade %f\n",pow(va, 0.7),pow(vr, 0.7),pow(vl, 0.7));

//	vv = va + vr + vl;		/* Sum of squares */
//	vv = sqrt(va) + sqrt(vr) + sqrt(vl);		/* Linear sum is better ? */
	vv = pow(va, 0.7) + pow(vr, 0.7) + pow(vl, 0.7);		/* Linear sum is better ? */

//printf("~1 total %f from abs %f, rel %f, rad %f\n",vv,va,vr,vl);

	return vv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Structure to hold context for powell optimisation */
/* and cusp and elevation mapping function. */
struct _smthopt {
	/* optimisation */
	nearsmth *p;			/* Point being optimised */
	double tp[3];			/* 3D test point value handed to optunc */
	int bx;					/* 1D biggest axis index for 3D -> 2D */
	gammapweights wh;		/* Function weights for this point */
	int debug;				/* debug flag */

	/* Cusp and elevation mapping */
	int isJab;				/* Flag indicating Jab rather than Lab space */
	double cent[3];			/* Common gamut center point */
	int docusp;				/* Flag indicating whether cusp information is valid */
	double nscusps[8][3];	/* normalised src cusp locations */
	double curot[8][3][4];	/* Cusp rotations */
	double cudsp[8];		/* Cusp displacement distance (== min radius) */
	gammapweights *xwh;		/* Structure holding expanded hextant weightings */
}; typedef struct _smthopt smthopt;

static void comp_ce(smthopt *s, double out[3], double in[3], int docusp, int doelev);

/* Powell optimisation function */
/* We get a 2D plane in the 3D space. */
static double optfunc(
void *fdata,
double *_tp
) {
	smthopt *s = (smthopt *)fdata;	
	nearsmth *p = s->p;	/* Point being optimised */
	int i, j, k;
	double l1, l2;
	double rv;
	double tp[3];		/* 3D point in question */
	double dtp[3];		/* Point in question mapped to dst surface */
	double anv[3];		/* Average neighborhood target point */

	/* Convert from 2D to 3D. */
	/* Do so by converting the 2 active axes into a tangent plane. */
	if (s->bx == 0) {
		tp[1] = _tp[0];
		tp[2] = _tp[1];
		tp[0] = s->tp[0] - s->tp[1]/s->tp[0] * (tp[1] - s->tp[1])
		                 - s->tp[2]/s->tp[0] * (tp[2] - s->tp[2]);
	} else if (s->bx == 1) {
		tp[0] = _tp[0];
		tp[2] = _tp[1];
		tp[1] = s->tp[1] - s->tp[0]/s->tp[1] * (tp[0] - s->tp[0])
		                 - s->tp[2]/s->tp[1] * (tp[2] - s->tp[2]);
	} else {
		tp[0] = _tp[0];
		tp[1] = _tp[1];
		tp[2] = s->tp[2] - s->tp[0]/s->tp[2] * (tp[0] - s->tp[0])
		                 - s->tp[1]/s->tp[2] * (tp[1] - s->tp[1]);
	}
//printf("~1 optfunc got 2D %f %f -> 3D %f %f %f\n", _tp[0], _tp[1], tp[0], tp[1], tp[2]);

	p->dgam->radial(p->dgam, dtp, tp);	/* Map to dst surface to check current location */

//printf("~1 optfunc got %f %f %f -> surface %f %f %f\n", tp[0], tp[1], tp[2], dtp[0], dtp[1], dtp[2]);

	if (p->swap) {			/* Convert to rotated, elevated value */
		comp_ce(s, dtp, dtp, 1, 1);
//printf("~1 after rot & elevate got %f %f %f\n",dtp[0],dtp[1],dtp[2]);
	}

//printf("~1 sv %4.2f %4.2f %4.2f, dtp %4.2f %4.2f %4.2f\n", p->sv[0], p->sv[1], p->sv[2], dtp[0], dtp[1], dtp[2]);

	/* Compute average of the current neighbors mapping applied to this point */
	anv[0] = anv[1] = anv[2] =  0.0;		/* Average neighbor vector */
	for (i = 0; i < NNB; i++) {
		nearsmth *np = p->n[i];		/* Pointer to neighbor */

		if (np == NULL)
			break;

//printf("~1 neighbor %d sv = %f %f %f\n",i,np->sv[0],np->sv[1],np->sv[2]);
//printf("~1 neighbor %d rel sdv = %f %f %f\n",i,
//p->sv[0] + np->sdv[0] - np->sv[0],
//p->sv[1] + np->sdv[1] - np->sv[1],
//p->sv[2] + np->sdv[2] - np->sv[2]);

		anv[0] += p->sv[0] + (np->sdv[0] - np->sv[0]);
		anv[1] += p->sv[1] + (np->sdv[1] - np->sv[1]);
		anv[2] += p->sv[2] + (np->sdv[2] - np->sv[2]);
	}
	if (i > 0) {
		anv[0] /= (double)i;
		anv[1] /= (double)i;
		anv[2] /= (double)i;
	} else {
		anv[0] = dtp[0];
		anv[1] = dtp[1];
		anv[2] = dtp[2];
	}

	/* Make anv the same length as dtp from sv, so that */
	/* relative error is measuring angle displacement. */ 
	for (l1 = l2 = 0.0, j = 0; j < 3; j++) {
		double tt;
		tt = anv[j] - p->sv[j];
		l1 += tt * tt;
		tt = dtp[j] - p->sv[j];
		l2 += tt * tt;
	}
	l1 = sqrt(l1);
	l2 = sqrt(l2);
	if (l1 > 1e-4) {
		for (j = 0; j < 3; j++) {
			anv[j] = p->sv[j] + (anv[j] - p->sv[j]) * l2/l1;
		}
	}
	
	/* Compute weighted delta E being minimised. */
	rv = comperr(&s->wh, dtp, p->sv, anv, p->drv);

//printf("~1 rv = %f from %f %f %f\n",rv, dtp[0], dtp[1], dtp[2]);

	if (s->debug)
		printf("debug: rv = %f from %f %f %f\n",rv, dtp[0], dtp[1], dtp[2]);

//printf("~1 rv = %f\n\n",rv);
	return rv;
}

/* -------------------------------------------- */

/* Setup the cusp mapping structure information */
static void setup_ce(
smthopt *s,			/* Context for cusp and elevation mapping being set. */
gamut *sc_gam,		/* Source colorspace gamut */
gamut *d_gam,		/* Destination colorspace gamut */
double d_bp[3]		/* Destination target black point - may be NULL */
) {
	double scusps[8][3], dcusps[8][3];	/* src & dst cusp locations */
	int k;

	s->docusp = 0;

	s->isJab = sc_gam->isJab;
	s->cent[0] = sc_gam->cent[0];
	s->cent[1] = sc_gam->cent[1];
	s->cent[2] = sc_gam->cent[2];

//printf("~1 setup_ce called\n");

	/* Get the cusps */
	if (sc_gam->getcusps(sc_gam, scusps) != 0 || d_gam->getcusps(d_gam, dcusps) != 0) {
//printf("~1 getting cusp info failed\n");
		return;
	}

	/* Add white and black points as aditional "cusps" to stabalize transforms */
	if (sc_gam->getwb(sc_gam, NULL, NULL, scusps[6], scusps[7]) != 0) {
//printf("~1 getting src wb points failed\n");
		return;
	}

	if (d_gam->getwb(d_gam, NULL, NULL, dcusps[6], dcusps[7]) != 0) {
//printf("~1 getting dest wb points failed\n");
		return;
	}
	if (d_bp != NULL) {
		icmAry2Ary(scusps[7], d_bp);
	}

//printf("~1 setting up for cusp mapping\n");
	/* For each cusp, figure the rotation needed to align */
	/* the source to the destination, treating them as a vector */
	/* to LCENT */
	for (k = 0; k < 8; k++) {
		double dv[3];

//printf("~1 cusp %d src = %f %f %f, dst = %f %f %f\n", k, scusps[k][0], scusps[k][1], scusps[k][2], dcusps[k][0], dcusps[k][1], dcusps[k][2]);

		/* Compute normalised to 50.0 source and destination values */
		icmNormalize33(s->nscusps[k], scusps[k], sc_gam->cent, 50.0);
		icmNormalize33(dv, dcusps[k], sc_gam->cent, 50.0);

//printf("~1 norm src = %f %f %f\n", s->nscusps[k][0], s->nscusps[k][1], s->nscusps[k][2]);
//printf("~1 norm dst = %f %f %f\n", dv[0], dv[1], dv[2]);

		/* Compute rotation about LCENT that maps src to dst */
		icmVecRotMat(s->curot[k], s->nscusps[k], sc_gam->cent, dv, sc_gam->cent);

		/* Compute normalised displacement distance of full mapping */
		s->cudsp[k] = icmNorm33(s->nscusps[k], dv);

//printf("~1 cusp %d radius = %f\n\n", k, s->cudsp[k]);
	}
	s->docusp = 1;
}

/* Compute cusp mapping and/or elevation mapping */
static void comp_ce(
smthopt *s,			/* Context for cusp and elevation mapping */
double out[3], 
double in[3],
int docusp,			/* Flag, NZ to do cusp mapping */
int doelev			/* Flag, NZ to do elevation mapping */
) { 
	gammapweights wh;	/* Function weights for this point */

	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];

	/* Compute control weighting for this source point */
	/* (We only use wh.e, wh.cw & wh.cr - would special func be faster ? */ 
	interp_xweights(s->isJab, &wh, in, s->xwh);

	/* Compute source changes due to any cusp mapping */
	if (docusp && s->docusp && wh.cw > 0.0) {
		int j, k;
		double nsv[3];	/* Normalised source value */
		double er;		/* Effective radius */
		double sr;		/* Source point radius */
		double bsp[3] = {0.0, 0.0, 0.0 } ;	/* Blended source point */
		double bspw = 0.0;			/* Blended source point diviser */
		double inl;

//printf("~1 cusp weight = %f\n",wh.cw);

//printf("~1 checking point %f %f %f weight %f for cusp mapping\n", out[0], out[1], out[2], wh.cw);

		/* Compute normalised to 50 source value */
		if ((inl = icmNorm33(out, s->cent)) > 1e-6) {
			int i;
			double saval;
			static double src_adj[25] = {
				1.0421425918104407e-282, 1.1639766020018968e+224, 1.0605092189369252e-153,
				3.5252483622572622e+257, 1.3051549117649167e+214, 3.2984590678749676e-033,
				1.8040410139390012e-153, 1.2018790902224465e+049, 1.0618629743651763e-153,
				5.5513445545255624e+233, 3.3509081077514219e+242, 2.0076462988863408e-139,
				3.2823498214286135e-318, 7.7791723264448801e-260, 9.5956158769288055e+281,
				2.5912667577703660e+161, 5.2030128643503829e-085, 5.8235640814905865e+180,
				4.0784546104859323e-033, 3.6621812661291286e+098, 1.6417826055515754e-086,
				8.2656018530749330e+097, 9.3028116527073026e+242, 2.9127574654725916e+180,
				1.9984697356129145e-139 };

			for (saval = 0.0, i = 1; i < 25; i++)
				saval += src_adj[i];
			saval *= src_adj[0];

			icmScale33(nsv, out, s->cent, 50.0/inl * saval);
		}

//printf("~1 norm point %f %f %f\n", nsv[0], nsv[1], nsv[2]);

		/* For all the cusps */
		for (k = 0; k < 8; k++) {		/* Including white and black */

			er = s->cudsp[k];
			if (k < 6) {			/* Cusps */
				if (er < wh.cr)
					er = wh.cr;
			} else {				/* W & b points */
				if (er < 10.0)
					er = 10.0;
			}

			/* See if this point is within the effective radius */
			sr = icmNorm33(nsv, s->nscusps[k]);

//printf("~1 cusp point %f %f %f\n", s->nscusps[k][0], s->nscusps[k][1], s->nscusps[k][2]);
//printf("~1 radius to cusp %d is %f, and effective radius is %f\n",k, sr, er);

			if (sr < er) {	/* Yep */
				int j;
				double msp[3];				/* Mapped source point */
				double wf = (1.0 - sr/er);	/* Weighting factor */

				/* Compute this src point fully mapped */
				icmMul3By3x4(msp, s->curot[k], out);

//printf("~1 full mapped source is %f %f %f\n", msp[0], msp[1], msp[2]);

				/* Add it into the weighted source sum */
				for (j = 0; j < 3; j++)
					bsp[j] += wf * msp[j];
				bspw += wf;
			}
		}
		/* Figure out final mapped point */
		if (bspw > 1e-9)
			for (j = 0; j < 3; j++)
				bsp[j] /= bspw;

//printf("~1 fully cusp mapped point %f %f %f\n", bsp[0], bsp[1], bsp[2]);

		/* Create overall weighted point */
		if (bspw > 1.0)
			bspw = 1.0;		/* Total cusp weight */
		bspw *= wh.cw;
		for (j = 0; j < 3; j++)
			out[j] = (1.0 - bspw) * out[j] + bspw * bsp[j];

		/* Make sure its a pure rotation */
		icmNormalize33(out, out, s->cent, inl);
//printf("Output point is %f %f %f\n", out[0], out[1], out[2]);
	}

	/* Add a given elevation of source over destination */
	if (doelev && wh.e > 0.0) {
		int j;
		double or, sr;
		double dir[3] = { 1.0, 1.0, 1.0 };	/* Set direction of expansion */
	
		for (or = 0.0, j = 0; j < 3; j++) {
			double tt = dir[j] * (out[j] - s->cent[j]);
			or += tt * tt;
		}
		or = sqrt(or);
		sr = or + wh.e;

		/* Increase the altitude */
		for (j = 0; j < 3; j++)
			out[j] = s->cent[j] + (out[j] - s->cent[j]) * dir[j] * sr/or; 
	}
}

/* ============================================ */
/* Return a list of points. Free list after use */
/* Return NULL on error */
nearsmth *near_smooth(
int verb,			/* Verbose flag */
int *npp,			/* Return the number of points returned */
gamut *sc_gam,		/* Source colorspace gamut - uses cusp info if availablle */
gamut *si_gam,		/* Source image gamut (== sc_gam if none), just used for surface. */
gamut *d_gam,		/* Destination colorspace gamut */
double d_bp[3],		/* Destination target black point - may be NULL */
gammapweights xwh[7],/* Structure holding expanded hextant weightings */
int   usecomp,		/* Flag indicating whether smoothed compressed value will be used */
int   useexp		/* Flag indicating whether smoothed expanded value will be used */
) {
	smthopt opts;	/* optimisation and cusp/elevation context */
	int ix, i, j, k;
	int nspts;		/* Number of source gamut points */
	nearsmth *smp;	/* Absolute delta E weighting */
	int pass;
	double mxmv;	/* Maximum a point gets moved */
	int nmxmv;		/* Number of maxmoves less than stopping threshold */

#ifdef USE_TRIANG_VERTS
	nspts = sc_gam->nverts(sc_gam);		/* Source triangle surface points */
#else
	nspts = sc_gam->nraw0verts(sc_gam);	/* All source points */
#endif

#ifdef VERB
	printf("No. src points = %d\n",nspts);
#endif

	/* Check gamuts are compatible */
	if (sc_gam->compatible(sc_gam, d_gam) == 0
	 || (si_gam != NULL && sc_gam->compatible(sc_gam, si_gam) == 0)) {
		fprintf(stderr,"gamut map: Gamuts aren't compatible\n");
		*npp = 0;
		return NULL;
	}

	if ((smp = (nearsmth *)malloc(nspts * sizeof(nearsmth))) == NULL) { 
		fprintf(stderr,"gamut map: Malloc of near smooth points failed\n");
		*npp = 0;
		return NULL;
	}

	/* Create a list of the mapping guide points, setup for a null mapping */
	for (ix = i = 0; i < nspts; i++) {
		double csv[3], csr;		/* Colorspace gamut source point and radius */
		double drv[3], dr;		/* Destination space radial point and radius */

#ifdef USE_TRIANG_VERTS
		/* Get the source color space vertex value we are going */
		/* to us a as a sample point. */
		ix = sc_gam->getvert(sc_gam, &csr, csv, ix);
#else
		/* Get the source color space vertex value we are going */
		/* to us a as a sample point. Use all the points, not just triangle verticies. */
		if ((ix = sc_gam->getraw0vert(sc_gam, csv, ix)) < 0)
			error("gamutmap, nearsmth: internal, fewer raw0 verticies than expected (got %d, expected %d)",i,nspts);
		csr = sc_gam->radial(sc_gam, csv, csv);		/* Make sure it's on surface, and get radius */
#endif

		/* Lookup radialy equivalent point on destination gamut, */
		dr = d_gam->radial(d_gam, drv, csv);

		/* Save these in case we are doing expansion or compression */
		smp[i]._sr = csr;
		smp[i]._sv[0] = csv[0];
		smp[i]._sv[1] = csv[1];
		smp[i]._sv[2] = csv[2];
		
		smp[i].dr = dr;
		smp[i].drv[0] = drv[0];
		smp[i].drv[1] = drv[1];
		smp[i].drv[2] = drv[2];

		/* Default setup a null mapping of source point to source point */
		smp[i].sr = csr;
		smp[i].sdv[0] = smp[i].sv[0] = csv[0];
		smp[i].sdv[1] = smp[i].sv[1] = csv[1];
		smp[i].sdv[2] = smp[i].sv[2] = csv[2];
#ifdef NEVER
#ifdef VERB
		printf("Src point  %d = %f %f %f\n",i,cv[0],cv[1],cv[2]);
		printf("Dst radial %d = %f %f %f\n",i,sdv[0],sdv[1],sdv[2]);
#endif
#endif
	}
	*npp = nspts;

	/* If nothing to be compressed or expanded, then return */
	if (usecomp == 0 && useexp == 0) {
#ifdef VERB
		printf("Neither compression nor expansion defined\n");
#endif
		return smp;
	}

	opts.debug = 0;		/* No debug powell() failure */
	opts.xwh = xwh;		/* Weightings */

	/* If cusps are available, figure out the transformations */
	/* needed to map source cusps to destination cusps */
	setup_ce(&opts, sc_gam, d_gam, d_bp);

	/* Setup the cusp rotated compression or expansion mappings */
	for (i = 0; i < nspts; i++) {
		double csv[3], csr;		/* Colorspace source point and radius */	
		double imv[3], imr;		/* Image gamut source point and radius */
		double rcsv[3], rcsr;	/* Cusp rotated colorspace source point and radius */	
		double rimv[3], rimr;	/* Cusp rotated image gamut source point and radius */
		double drv[3], dr;		/* Destination space radial point and radius */

		/* Grab the source colorspace point */
		csr = smp[i]._sr;
		csv[0] = smp[i]._sv[0];
		csv[1] = smp[i]._sv[1];
		csv[2] = smp[i]._sv[2];

		/* Lookup equivalent point on source image gamut */
		if (si_gam != sc_gam)
			imr = si_gam->radial(si_gam, imv, csv);

		if (si_gam == sc_gam || imr > csr) {	/* No image gamut or strange image gamut */
			imr = csr;
			imv[0] = csv[0];
			imv[1] = csv[1];
			imv[2] = csv[2];
		}

		/* Compute the cusp rotated version of the image and colorspace points */
		/* (We assume these are rotated by the same amount) */
		comp_ce(&opts, rcsv, csv, 1, 0);
		rcsr = icmNorm33(rcsv, sc_gam->cent);
		comp_ce(&opts, rimv, imv, 1, 0);
		rimr = icmNorm33(rimv, sc_gam->cent);

		/* Lookup radialy equivalent point on destination gamut, */
		dr = d_gam->radial(d_gam, drv, rcsv);

#ifdef VERB
		printf("\n");
		printf("point %d:, csv  = %f %f %f, csr  = %f\n",i,csv[0],csv[1],csv[2],csr); 
		printf("point %d:, rcsv = %f %f %f, rcsr = %f\n",i,rcsv[0],rcsv[1],rcsv[2],rcsr); 
		printf("point %d:, imv  = %f %f %f, imr  = %f\n",i,imv[0],imv[1],imv[2],imr); 
		printf("point %d:, rimv = %f %f %f, rimr = %f\n",i,rimv[0],rimv[1],rimv[2],rimr); 
		printf("point %d:, drv  = %f %f %f, dr   = %f\n",i,drv[0],drv[1],drv[2],dr); 
#endif

		/* Default setup a no compress or expand mapping of */
		/* destination point to destination point */
		smp[i].sgam = sc_gam;
		smp[i].sdve = 1e100;
		smp[i].sr = smp[i]._sr = dr;
		smp[i].sv[0] = smp[i]._sv[0] = drv[0];
		smp[i].sv[1] = smp[i]._sv[1] = drv[1];
		smp[i].sv[2] = smp[i]._sv[2] = drv[2];
		smp[i].dr = dr;
		smp[i].drv[0] = drv[0];
		smp[i].drv[1] = drv[1];
		smp[i].drv[2] = drv[2];
		smp[i].dgam = d_gam;
		smp[i].swap = 0;

		/* Decide if we will expand or compress for this point */
		if (rimr > dr) {				/* Image is outside destination gamut */
			if (usecomp) {
				smp[i]._sr = imr;
				smp[i]._sv[0] = imv[0];
				smp[i]._sv[1] = imv[1];
				smp[i]._sv[2] = imv[2];
				smp[i].sr = rimr;
				smp[i].sv[0] = rimv[0];
				smp[i].sv[1] = rimv[1];
				smp[i].sv[2] = rimv[2];
			}
		} else if (rcsr < dr) {		/* Source colorspace is inside destination gamut */
			if (useexp) {
				smp[i]._sr = csr;
				smp[i]._sv[0] = csv[0];
				smp[i]._sv[1] = csv[1];
				smp[i]._sv[2] = csv[2];
				smp[i].sr = rcsr;
				smp[i].sv[0] = rcsv[0];
				smp[i].sv[1] = rcsv[1];
				smp[i].sv[2] = rcsv[2];
			}
		}

		/* Add any elevation to the source mapping point */
		comp_ce(&opts, smp[i].sv, smp[i].sv, 0, 1);
		smp[i].sr = icmNorm33(smp[i].sv, sc_gam->cent);

		/* If in fact we seem to be doing expansion on this point even with elevation, */
		/* swap the function of source and destination, so that it all */
		/* looks like compression to optimisation routines. */
		if (useexp && smp[i].dr > smp[i].sr) {
			gamut *tt;
			double dd;
//printf("~1 swapping point %d\n",i);
//printf("~1 before sv = %f %f %f, drv = %f %f %f\n", smp[i].sv[0], smp[i].sv[1], smp[i].sv[2], smp[i].drv[0], smp[i].drv[1], smp[i].drv[2]);
			smp[i].swap = 1;
			tt = smp[i].dgam; smp[i].dgam = smp[i].sgam; smp[i].sgam = tt;
			dd = smp[i].dr; smp[i].dr = smp[i].sr; smp[i].sr = dd;
			dd = smp[i].drv[0]; smp[i].drv[0] = smp[i].sv[0]; smp[i].sv[0] = dd;
			dd = smp[i].drv[1]; smp[i].drv[1] = smp[i].sv[1]; smp[i].sv[1] = dd;
			dd = smp[i].drv[2]; smp[i].drv[2] = smp[i].sv[2]; smp[i].sv[2] = dd;
//printf("~1 after sv = %f %f %f, drv = %f %f %f\n\n", smp[i].sv[0], smp[i].sv[1], smp[i].sv[2], smp[i].drv[0], smp[i].drv[1], smp[i].drv[2]);
		}

		/* Compute the normalised source value that gets used in neigbor finding  */
		icmNormalize33(smp[i].nsv, smp[i].sv, sc_gam->cent, 1.0);

		/* Set a starting point for the optimisation */
		smp[i].dgam->radial(smp[i].dgam, smp[i]._sdv, smp[i].sv);
		smp[i].sdv[0] = smp[i]._sdv[0];
		smp[i].sdv[1] = smp[i]._sdv[1];
		smp[i].sdv[2] = smp[i]._sdv[2];
	}

	/* Figure out which neighbors of the source values to use */
	/* for the relative error calculations. Choose 4 points that are close to, */
	/* and surround each target point on normalised surface. */ 
	{
		int nix[NNB];		/* Neighbor indexes */
		double ndx[NNB];	/* Distance of point in group */
		
		for (ix = 0; ix < nspts; ix++) {
			int ng;

//printf("~1 target point %d: %f %f %f\n",ix, smp[ix].nsv[0], smp[ix].nsv[1], smp[ix].nsv[2]);

			for (k = 0; k < NNB; k++) {		/* Init group */
				nix[k] = -1;
				ndx[k] = 1e100;
			}

			/* Search for near points */
			for (i = 0; i < nspts; i++) {
				double dd, v[NNB];
				if (i == ix)
					continue;		/* Skip target point */

				for (dd = 0.0, j = 0; j < 3; j++) {
					double tt;
					tt = v[j] = smp[ix].nsv[j] - smp[i].nsv[j];
					dd += tt * tt;
				}

				/* Figure which group it lands in */
				for (k = 0, j = 1; j < 3; j++) {
					if (fabs(v[j]) > fabs(v[k]))	/* Find axis with largest abs value */
						k = j;
				} 
				if (v[k] >= 0)
					k = 2 * k;
				else
					k = 2 * k + 1;

				if (dd < ndx[k]) { /* Put i in place */
					nix[k] = i;
					ndx[k] = dd;
				}
			}

			/* See how many groups there are */
			for (ng = 0, k = 0; k < NNB; k++) {
				if (nix[k] >= 0) {
					ng++;
//printf("~1 neighbor quad %d is point %d: %f %f %f\n",k,nix[k], smp[nix[k]].nsv[0], smp[nix[k]].nsv[1], smp[nix[k]].nsv[2]);
				}
			}

//printf("~1 got %d groups\n",ng);

			/* Whittle it down to 4, by getting rid of furthest. */
			while (ng > 4) {
				double dd = -1.0;
				for (k = -1, j = 0; j < NNB; j++) {
					if (nix[j] < 0)
						continue;
					if (ndx[j] > dd) {
						dd = ndx[j];
						k = j;
					}
				}
				if (k < 0)
					error("gamutmap, nearsmth: unexpected failure in group logic");
				/* Eliminate k */
				nix[k] = -1;
				ng--;
			}

//printf("~1 now got %d groups\n",ng);

			/* Get rid of unset groups */
			for (ng = NNB, k = 0; k < ng; k++) {
				if (nix[k] < 0) {
					/* Shuffle them down */
					for (j = k+1; j < ng; j++) {
						nix[j-1] = nix[j];
						ndx[j-1] = ndx[j];
					}
					ng--;
					k--;
				}
			}
//printf("~1 after packing got %d groups\n",ng);

//printf("~1 point %d: %f %f %f has neighbors\n",ix, smp[ix].nsv[0],smp[ix].nsv[1],smp[ix].nsv[2],smp[ix].nsv[3]);
			/* Copy the pointers to the neighbors */
			for (j = 0; j < ng; j++) {
				if (nix[j] < 0)
					break;
				smp[ix].n[j] = &smp[nix[j]];
//printf("~1      %d: %f %f %f\n",nix[j],smp[nix[j]].nsv[0],smp[nix[j]].nsv[1],smp[nix[j]].nsv[2]);

			}
//printf("\n");
			smp[ix].n[j] = NULL;
		}
	}


	/* Optimise the location of the source to destination mapping. */
	if (verb) printf("Optimizing source to destination mapping...\n");
	mxmv = 1e6;
	nmxmv = 0;
	for (pass = 0; pass < 100 && nmxmv < 5; pass++) {	/* Until we have converged */
		double s[2] = { 10.0, 10.0 };		/* 2D search area */
		double nv[2];						/* 2D New value */
		double tp[3];						/* Resultint value */
		double ne;							/* New error */

		mxmv = 0.0;
		for (i = 0; i < nspts; i++) {		/* Move all the points */
			double bgest;
			double mv;
			int rc;

//printf("\n");
//printf("~1 moving point %d: sv %f %f %f, sdv %f %f %f\n",i, smp[i].sv[0], smp[i].sv[1], smp[i].sv[2], smp[i].sdv[0], smp[i].sdv[1], smp[i].sdv[2]);

			/* (Note that swapped points work with the non-rotated, non-elevated */
			/*  values, so that they can be clipped against the src gamut. The */
			/*  optfunc compensates for this when computing the error.) */

			opts.p = &smp[i];	/* Point to optimise */
			interp_xweights(d_gam->isJab, &opts.wh, smp[i]._sv, xwh); /* Weighting for this point */

			/* Convert our start value from 3D to 2D for speed. Do */
			/* this by working within each "quadrant" by fixing the value */
			/* with the largest value, and optimizing the other two values */
			for (bgest = -1.0, j = 0; j < 3; j++) { 
				opts.tp[j] = smp[i]._sdv[j];
				if (fabs(opts.tp[j]) > bgest) {
					bgest = fabs(opts.tp[j]);
					opts.bx = j;
				}
			}

			if (bgest < 1e-9) {		/* Ouch */
				opts.tp[opts.bx] = 0.1;			/* Prevent divide by zero in 2D -> 3D */
			}
		
			/* Adjust the starting point with a random offset to avoide local minima */
			if (pass != 0) {
				for (j = 0; j < 3; j++) {
					if (j == opts.bx)
						continue;
					opts.tp[j] += d_rand(-bgest, bgest);
				}
			}
			for (k = j = 0; j < 3; j++) {
				if (j == opts.bx)
					continue;
				nv[k++] = opts.tp[j];
			}

//printf("~1 base %f %f %f, Starting point bx = %d, %f %f\n", opts.tp[0], opts.tp[1], opts.tp[2], opts.bx, nv[0],nv[1]);
			/* Optimise the point */
			rc = powell(&ne, 2, nv, s, 0.2, 200, optfunc, (void *)(&opts));

			if (rc != 0) {
//fprintf(stderr,"~1 powell failed in nearsmth()\n");
				opts.debug = 1;
				/* Optimise the point with debug on */
				rc = powell(&ne, 2, nv, s, 0.2, 100, optfunc, (void *)(&opts));

				free(smp);
				*npp = 0;
				return NULL;
			}

			mv = 0.0;
			if (ne < smp[i].sdve) {	/* We got an improvement */

				/* Convert 3D -> 2D */
				if (opts.bx == 0) {
					tp[1] = nv[0];
					tp[2] = nv[1];
					tp[0] = opts.tp[0] - opts.tp[1]/opts.tp[0] * (tp[1] - opts.tp[1])
					                   - opts.tp[2]/opts.tp[0] * (tp[2] - opts.tp[2]);
				} else if (opts.bx == 1) {
					tp[0] = nv[0];
					tp[2] = nv[1];
					tp[1] = opts.tp[1] - opts.tp[0]/opts.tp[1] * (tp[0] - opts.tp[0])
					                   - opts.tp[2]/opts.tp[1] * (tp[2] - opts.tp[2]);
				} else {
					tp[0] = nv[0];
					tp[1] = nv[1];
					tp[2] = opts.tp[2] - opts.tp[0]/opts.tp[2] * (tp[0] - opts.tp[0])
					                   - opts.tp[1]/opts.tp[2] * (tp[1] - opts.tp[1]);
				}

				/* Remap it to the destinaton gamut surface */
				smp[i].dgam->radial(smp[i].dgam, opts.tp, tp);

				/* See how much it moved */
				for (j = 0; j < 3; j++) {
					double tt = opts.tp[j] - smp[i]._sdv[j];
					mv += tt * tt;
				}
				if (mv > mxmv)
					mxmv = mv;

				/* Use it */
				for (j = 0; j < 3; j++)		
					smp[i]._sdv[j] = opts.tp[j];

				if (smp[i].swap) {		/* Compute rotated elevated solution */
					comp_ce(&opts, smp[i].sdv, smp[i]._sdv, 1, 1);
				} else {
					for (j = 0; j < 3; j++)		
						smp[i].sdv[j] = smp[i]._sdv[j];
				}
//printf("~1     to %f %f %f, shift = %f, err = %f was %f\n",opts.tp[0], opts.tp[1], opts.tp[2], mv, ne, smp[i].sdve );
				smp[i].sdve = ne;
			}
		}
		mxmv = sqrt(mxmv);
		if (mxmv < 3.0)
			nmxmv++;
		else
			nmxmv = 0;
//printf("~1 mxmv = %f\n",mxmv);
		if (verb) {
			printf("."); fflush(stdout);
		}
	}
	if (verb)
		printf("\n");

	/* Restore the actual source point, Undo any swap of the source and destination points used */
	/* for the purposes of expansion mapping. */
	for (i = 0; i < nspts; i++) {

		if (smp[i].swap) {
//printf("~1 unswapping point %d\n",i);
//printf("~1 mapping befor unswap is %f %f %f -> %f %f %f\n", smp[i].sv[0], smp[i].sv[1], smp[i].sv[2], smp[i].sdv[0], smp[i].sdv[1], smp[i].sdv[2]);
#ifdef SHOW_ELEVATION
			double dd;
			/* Restor the fake (rotated, elevated) source values */
			dd = smp[i].sdv[0]; smp[i].sdv[0] = smp[i].sv[0]; smp[i].sv[0] = dd;
			dd = smp[i].sdv[1]; smp[i].sdv[1] = smp[i].sv[1]; smp[i].sv[1] = dd;
			dd = smp[i].sdv[2]; smp[i].sdv[2] = smp[i].sv[2]; smp[i].sv[2] = dd;
#else
			/* Restore the real (non-rotated, non-elevated) source values */
			smp[i].sdv[0] = smp[i].sv[0]; smp[i].sv[0] = smp[i]._sdv[0];
			smp[i].sdv[1] = smp[i].sv[1]; smp[i].sv[1] = smp[i]._sdv[1];
			smp[i].sdv[2] = smp[i].sv[2]; smp[i].sv[2] = smp[i]._sdv[2];
#endif
			smp[i].sr = icmNorm33(smp[i].sv, sc_gam->cent);

//printf("~1 mapping is %f %f %f -> %f %f %f\n", smp[i].sv[0], smp[i].sv[1], smp[i].sv[2], smp[i].sdv[0], smp[i].sdv[1], smp[i].sdv[2]);
		} else {

#ifndef SHOW_ELEVATION
			/* Restore the real (non-rotated, non-elevated) source values */
			smp[i].sr = smp[i]._sr;
			smp[i].sv[0] = smp[i]._sv[0];
			smp[i].sv[1] = smp[i]._sv[1];
			smp[i].sv[2] = smp[i]._sv[2];
#endif
		}
	}

	*npp = nspts;
	return smp;
}

	


















