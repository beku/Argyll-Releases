
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
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
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
 *
 *		It might work better if the cusp mapping had separate control
 *		over the L and h degree of map, as well as the L and h effective radius ?
 *		That way, saturation hue distortions with L might be reduced.
 *
 *       Improve error handling.
 * 
 *		Major defect with some gamut combinations is "button" around
 *		cusps. Not sure what the mechanism is, since it's not obvious
 *		from the 3D vector plots what the cause is. (fixed ?)
 *		Due to poor internal control ?
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "icc.h"
#include "numlib.h"
#include "rspl.h"
#include "gamut.h"
#include "nearsmth.h"
#include "vrml.h"

#undef SAVE_VRMLS		/* Save various vrml's */
#undef PLOT_MAPPING_INFLUENCE		/* Plot sci_gam colored by dominant guide influence: */ 
		                /* Absolute = red, Relative = yellow, Radial = blue, Depth = green */
#define VERB 0 			/* If <= 1, print progress headings */
						/* if  > 1, print information about everything */
#undef SHOW_NEIGB_WEIGHTS		/* Show the weighting for each point of neighbours */

#define LIGHT_L 70.0	/* "light" L/J value */
#define DARK_L  5.0		/* "dark" L/J value */
#define NEUTRAL_C  20.0	/* "neutral" C value */
#define NO_TRIALS 6		/* Number of random trials */
//#define MAXITTERS 10	/* Number of relative error itterations */
#define MAXITTERS 10	/* Number of relative error itterations */
#define NO_SMTH_ITERS 4	/* Number of pre-itteration smoothings */
#define ITTER_STOP 0.5	/* Max Delta E to stop at */
#define OSHOOT 1.3		/* Amount of overshoot/damping to use */
#define RADIAL_SUBVEC	/* Make sub-vectors always radial direction */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
#if defined(VERB)
# define VA(xxxx) printf xxxx
# if VERB > 1
#  define VB(xxxx) printf xxxx
# else
#  define VB(xxxx) 
# endif
#else
# define VA(xxxx) 
# define VB(xxxx) 
#endif

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* A set of functions to help handle the weighting configuration */

/* Copy non-negative values from one set of weights to another */
void near_wcopy(
gammapweights *dst,
gammapweights *src
) {

#define NSCOPY(xxx) dst->xxx = src->xxx >= 0.0 ? src->xxx : dst->xxx
//#define NSCOPY(xxx) if (src->xxx >= 0.0) { \
//						printf("Setting %s to %f\n",#xxx, src->xxx); \
//						dst->xxx = src->xxx;						\
//					}

	NSCOPY(a.o);
	NSCOPY(a.w.l);
	NSCOPY(a.w.c);
	NSCOPY(a.w.h);
	NSCOPY(r.o);
	NSCOPY(r.w.l);
	NSCOPY(r.w.c);
	NSCOPY(r.w.h);
	NSCOPY(r.rdl);
	NSCOPY(r.rdh);
	NSCOPY(l.o);
	NSCOPY(l.w.l);
	NSCOPY(l.w.c);
	NSCOPY(l.w.h);
	NSCOPY(d.co);
	NSCOPY(d.xo);
	NSCOPY(c.w.l);
	NSCOPY(c.w.c);
	NSCOPY(c.w.h);
	NSCOPY(c.cx);
#undef NSCOPY
}

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
	NSBLEND(r.rdl);
	NSBLEND(r.rdh);
	NSBLEND(l.o);
	NSBLEND(l.w.l);
	NSBLEND(l.w.c);
	NSBLEND(l.w.h);
	NSBLEND(d.co);
	NSBLEND(d.xo);
	NSBLEND(c.w.l);
	NSBLEND(c.w.c);
	NSBLEND(c.w.h);
	NSBLEND(c.cx);
#undef NSBLEND
}

/* Expand the compact form of weights into the explicit form. */
/* The explicit form is light and dark of red, yellow, green, cyan, blue, magenta & neutral*/
/* Return nz on error */
int expand_weights(gammapweights out[14], gammapweights *in) {
	int i, j;

	/* Set the usage of each slot */
	out[0].ch = gmm_light_red;
	out[1].ch = gmm_light_yellow;
	out[2].ch = gmm_light_green;
	out[3].ch = gmm_light_cyan;
	out[4].ch = gmm_light_blue;
	out[5].ch = gmm_light_magenta;
	out[6].ch = gmm_light_neutral;

	out[7].ch = gmm_dark_red;
	out[8].ch = gmm_dark_yellow;
	out[9].ch = gmm_dark_green;
	out[10].ch = gmm_dark_cyan;
	out[11].ch = gmm_dark_blue;
	out[12].ch = gmm_dark_magenta;
	out[13].ch = gmm_dark_neutral;

//printf("\n~1 expand weights called\n");

	/* mark output so we can recognise having been set or not */
	for (i = 0; i < 14; i++)
		out[i].set = 0;

	/* Expand the compact form to explicit. */

	/* First is default */
	for (i = 0; in[i].ch != gmm_end; i++) {
		if (in[i].ch == gmm_end)
			break;
		if (in[i].ch == gmm_ignore)
			continue;

		if (in[i].ch == gmm_default) {
			for (j = 0; j < 14; j++) {
//printf("~1 Setting %d 0x%x with 0x%x (default)\n",j,out[j].ch,in[i].ch);
				if ((in[i].ch & out[j].ch) == out[j].ch) {
					near_wcopy(&out[j], &in[i]);
					out[j].set = 1;
				}
			}
		}
	}

	/* Then light or dark */
	for (i = 0; in[i].ch != gmm_end; i++) {
		if (in[i].ch == gmm_end)
			break;
		if (in[i].ch == gmm_ignore)
			continue;

		if (in[i].ch == gmm_light_colors
		  || in[i].ch == gmm_dark_colors) {
			for (j = 0; j < 14; j++) {
				if ((in[i].ch & out[j].ch) == out[j].ch) {
//printf("~1 Setting %d 0x%x with 0x%x (light or dark)\n",j,out[j].ch,in[i].ch);
					near_wcopy(&out[j], &in[i]);
					out[j].set = 1;
				}
			}
		}
	}

	/* Then light and dark colors */
	for (i = 0; in[i].ch != gmm_end; i++) {
		if (in[i].ch == gmm_end)
			break;
		if (in[i].ch == gmm_ignore)
			continue;

		if ((in[i].ch & gmc_l_d) == gmc_l_d
		 && (in[i].ch & gmc_colors) != gmc_colors) {
			for (j = 0; j < 14; j++) {
				if ((in[i].ch & out[j].ch) == out[j].ch) {
//printf("~1 Setting %d 0x%x with 0x%x (light and dark color)\n",j,out[j].ch,in[i].ch);
					near_wcopy(&out[j], &in[i]);
					out[j].set = 1;
				}
			}
		}
	}

	/* Last pass is light or dark colors */
	for (i = 0; in[i].ch != gmm_end; i++) {
		if (in[i].ch == gmm_end)
			break;
		if (in[i].ch == gmm_ignore)
			continue;

		if (((in[i].ch & gmc_l_d) == gmc_light
		  || (in[i].ch & gmc_l_d) == gmc_dark)
		 && (in[i].ch & gmc_colors) != gmc_colors) {
			for (j = 0; j < 14; j++) {
				if ((in[i].ch & out[j].ch) == out[j].ch) {
//printf("~1 Setting %d 0x%x with 0x%x (light or dark color)\n",j,out[j].ch,in[i].ch);
					near_wcopy(&out[j], &in[i]);
					out[j].set = 1;
				}
			}
		}
	}

	/* Check every slot has been set */
	for (i = 0; i < 14; i++) {
		if (out[i].set == 0) {
//printf("~1 set %d hasn't been initialized\n",i);
			return 1;
		}
	}
	return 0;
}

/* Tweak weights acording to extra cmy cusp flags or rel override */
void tweak_weights(gammapweights out[14], int dst_cmymap, int rel_oride)  {
	int i;

	for (i = 0; i < 14; i++) {
		if (((dst_cmymap & 0x1) && (out[i].ch & gmc_cyan))
		 || ((dst_cmymap & 0x2) && (out[i].ch & gmc_magenta))
		 || ((dst_cmymap & 0x4) && (out[i].ch & gmc_yellow))) {
//printf("~1 Setting %d 0x%x to 100% cusp map\n",i,out[i].ch);
			out[i].c.w.l = 1.0;		/* 100% mapping */
			out[i].c.w.c = 1.0;
			out[i].c.w.h = 1.0;
			out[i].c.cx = 1.0;		/* No expansion */
		}

		if (rel_oride == 1) {		/* A high saturation "clip" like mapping */
			out[i].r.o = 0.0;		/* No relative weight */
			out[i].r.rdl = 1.0;		/* No relative neighbourhood */
			out[i].r.rdh = 1.0;		/* No relative neighbourhood */
			out[i].d.co = 0.0;		/* No depth weighting */
			out[i].d.xo = 0.0;		/* No depth weighting */
			
		} else if (rel_oride == 2) {	/* A maximal feature preserving mapping */
			out[i].r.o = 2.0;		/* Extra relative weight */
			out[i].r.rdl *= 1.4;	/* Extra neighbourhood size */
			out[i].r.rdh *= 1.4;	/* Extra neighbourhood size */
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
	for (i = 0; i < 14; i++)
		near_wblend(&dst[i], &src1[i], wgt1, &src2[i], wgt2);
}

/* Given a point location, return the interpolated weighting values at that point. */
void interp_xweights(gamut *gam, gammapweights *out, double pos[3], gammapweights in[14]) {
	double JCh[3];
	int li, ui;			/* The two hue indexes the color is between */
	double lh, uh;		/* Lower/upper hue of two colors */
	double lw, uw;		/* Lower/upper blend values */
	double cusps[6][3];
	gammapweights light, dark;

	/* Convert to polar */
	icmLab2LCh(JCh, pos);

	if (gam->getcusps(gam, cusps) != 0) {
		int isJab = gam->isJab ? 1 : 0;

		/* Figure out what hextant we're between using generic cusps */
		for (li = 0; li < 6; li++) {
			ui = li < 5 ? li + 1 : 0;

			lh = gam_hues[isJab][li];
			uh = gam_hues[isJab][ui];
			if (uh < lh) {
				if (JCh[2] < uh)
					JCh[2] += 360.0;
				uh += 360.0;
			}
			if (JCh[2] >= lh && JCh[2] < uh)
				break;
		}

	} else {

		/* Locate the source cusps that this point lies between */ 
		for (li = 0; li < 6; li++) {
			double tt[3];
			ui = li < 5 ? li + 1 : 0;

			icmLab2LCh(tt, cusps[li]);
			lh = tt[2];
			icmLab2LCh(tt, cusps[ui]);
			uh = tt[2];

			if (uh < lh) {
				if (JCh[2] < uh)
					JCh[2] += 360.0;
				uh += 360.0;
			}
			if (JCh[2] >= lh && JCh[2] < uh)
				break;
		}
	}

	/* Compute weights */
	uw = (JCh[2] - lh)/(uh - lh);
	uw = uw * uw * (3.0 - 2.0 * uw);	/* Apply spline to smooth interpolation */
	lw = (1.0 - uw);

	/* Blend weights at the two hues */
	near_wblend(&light, &in[li], lw, &in[ui], uw);
	near_wblend(&dark, &in[7 + li], lw, &in[7 + ui], uw);

	/* If we're close to the center, blend to the neutral weight */
	if (JCh[1] < NEUTRAL_C) {
		lw = (NEUTRAL_C - JCh[1])/NEUTRAL_C;
		uw = (1.0 - lw);
		near_wblend(&light, &in[6], lw, &light, uw);
		near_wblend(&dark, &in[7 + 6], lw, &dark, uw);
	}

	/* Figure out where we are between light and dark, */
	/* and create blend between their weightings */
	uw = (JCh[0] - DARK_L)/(LIGHT_L - DARK_L);
	if (uw > 1.0)
		uw = 1.0;
	else if (uw < 0.0)
		uw = 0.0;
	uw = uw * uw * (3.0 - 2.0 * uw);	/* Apply spline to smooth interpolation */
	lw = (1.0 - uw);
	near_wblend(out, &dark, lw, &light, uw);
}

/* Compute the weighted delta E squared of in1 - in2 */
/* (This is like the CIE DE94) */
static double wdesq(
double in1[3],
double in2[3],
double lweight,
double cweight,
double hweight
) {
	double desq, dhsq;
	double dlsq, dcsq;
	double vv;

//printf("~1 wdesq got %f %f %f and %f %f %f\n", in1[0], in1[1], in1[2], in2[0], in2[1], in2[2]);
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
		double dc, c1, c2;

		/* Compute chromanance for the two colors */
		c1 = sqrt(in1[1] * in1[1] + in1[2] * in1[2]);
		c2 = sqrt(in2[1] * in2[1] + in2[2] * in2[2]);

		dc = c1 - c2;
		dcsq = dc * dc;

		/* ~~ Doing dcsq = sqrt(dcsq); here seemed */
		/* to improve the saturation result. How did this work ?? */
	}

	/* Compute delta hue squared */
	if ((dhsq = desq - dlsq - dcsq) < 0.0)
		dhsq = 0.0;

	vv = lweight * dlsq + cweight * dcsq + hweight * dhsq;
//printf("~1 returning wdesq %f from %f * %f + %f * %f + %f * %f\n", fabs(vv),lweight, dlsq, cweight, dcsq, hweight, dhsq);

	return fabs(vv);	/* Avoid -0.0 */
}

/* Given the weighting structure and the relevant point locations */
/* return the total weighted error squared. */
static double comperr(
nearsmth *p,		/* Guide point data */
gammapweights *w,	/* weightings */
double relscale,	/* Scale relative error weighting */
double dtp[3],		/* Dest test point being evaluated */
double aodv[3],		/* Weighted destination closest value to source */
double anv[3],		/* Surround average displaced to dest source point */
double drv[3],		/* Source mapped radially to dest */
double dcratio,		/* Depth compression ratio of mapping */
double dxratio,		/* Depth expansion ratio of mapping */
double ogam			/* Normalized amount out of gamut */
) {
	double a_o, r_o;
	double va, vr = 0.0, vl, vd, vv = 0.0, vo;

	/* Absolute, Delta E^2 between destination closest and test point */
	/* aodv is already positioned acording to the LCh weights, */
	/* so weight as per average of these */
	a_o = w->a.o * (w->a.w.l + w->a.w.c + w->a.w.h)/3.0;
	va = wdesq(aodv, dtp, a_o, a_o, a_o);
//	va = wdesq(aodv, dtp, w->a.o * w->a.w.l, w->a.o * w->a.w.c, w->a.o *w->a.w.h);

	r_o = relscale * w->r.o;

	/* Relative. Delta E^2 between surround average displaced source and test point */
	vr = wdesq(anv, dtp, r_o * w->r.w.l, r_o * w->r.w.c, r_o * w->r.w.h);

	/* Radial. Delta E^2 between source mapped radially to dest and test point */
	vl = wdesq(drv, dtp, w->l.o * w->l.w.l, w->l.o * w->l.w.c, w->l.o * w->l.w.h);

	/* Depth ratio error^2. */
	vd = w->d.co * dcratio * dcratio
	   + w->d.xo * dxratio * dxratio;

	/* Out of gamut error */
	vo = ogam * 10000.0;
	vo *= vo;

	/* Diagnostic values */
	p->dbgv[0] = va;
	p->dbgv[1] = vr;
	p->dbgv[2] = vl;
	p->dbgv[3] = vd;
	
	vv = va + vr + vl + vd + vo;		/* Sum of squares */
//	vv = sqrt(va) + sqrt(vr) + sqrt(vl) + sqrt(vd) + sqrt(vo);		/* Linear sum is better ? */
//	vv = pow(va, 0.7) + pow(vr, 0.7) + pow(vl, 0.7) + pow(vd, 0.7) + pow(vo, 0.7);		/* Linear sum is better ? */

//if (ogam > 0.2)
//printf("~1 total %f from abs %f, rel %f, rad %f, depth %f, ogam %f\n",vv,va,vr,vl,vd,vo);

#ifdef NEVER
	printf("~1 dtp %f %f %f\n", dtp[0], dtp[1], dtp[2]);
	printf("~1 va = %f from aodv %f %f %f, weight %f\n", va, aodv[0], aodv[1], aodv[2], a_o);
	printf("~1 vr = %f from anv %f %f %f, weights %f %f %f\n", vr, anv[0], anv[1], anv[2], r_o * w->r.w.l, r_o * w->r.w.c, r_o * w->r.w.h);
	printf("~1 vl = %f from drv %f %f %f, weights %f %f %f\n", vl, drv[0], drv[1], drv[2], w->l.o * w->l.w.l, w->l.o * w->l.w.c, w->l.o * w->l.w.h);
	printf("~1 vd = %f from d.co %f d.xo %f, weights %f %f\n", vd, w->d.co,w->d.xo,dcratio * dcratio,dxratio * dxratio);
	printf("~1 vo = %f\n", vo);
	printf("~1 return vv = %f\n", vv);
#endif /* NEVER */

	return vv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Structure to hold context for powell optimisation */
/* and cusp and elevation mapping function. */
struct _smthopt {
	/* optimisation */
	int pass;				/* Itteration round */
	int ix;					/* Index of point being optimized */
	nearsmth *p;			/* Point being optimised */
	int useexp;				/* Flag indicating whether expansion is permitted */
	int debug;				/* debug flag */

	/* Setup state */
	int isJab;				/* Flag indicating Jab rather than Lab space */
	int docusp;				/* Flag indicating whether cusp information is valid */
	gammapweights *xwh;		/* Structure holding expanded hextant weightings */

	/* Cusp alignment and elevation mapping */
	/* 0 = src, 1 = dst */
	double rot[2][3][4];	/* Rotation to align to black/white center */
	double irot[2][3][4];	/* Inverse rotation */
	double cusp_lab[2][9][3];	/* Cusp + B&W + grey rotated Lab value */ 
	double cusp_lch[2][6][3];	/* Cusp LCH value */ 
	double cusp_pe[2][6][4];	/* L direction plane equations per segment */
	double cusp_bc[2][6][2][3][3];	/* light/dark to/from baricentic transform matrix */

	/* Inversion support */
	double tv[3];
	gammapweights *wt;		/* Weights for this inversion */

	double mm[3][4];		/* Direction alignment rotation */
	double m2[2][2];		/* Additional matrix to alight a with L axis */
	double manv[3];			/* anv[] transformed by mm and m2 */

}; typedef struct _smthopt smthopt;

static void init_ce(smthopt *s, gamut *sc_gam, gamut *d_gam, int src_kbp, int dst_kbp, double d_bp[3]);
static void comp_ce(smthopt *s, double out[3], double in[3], gammapweights *wt);
static void inv_comp_ce(smthopt *s, double out[3], double in[3], gammapweights *wt);

static double spow(double arg, double ex) {
	if (arg < 0.0)
		return -pow(-arg, ex);
	else
		return pow(arg, ex);
}

/* Powell optimisation function for setting minimal absolute error target point. */
/* We get a 2D plane in the 3D space, of the destination point, */
/* who's location we are optimizing. */
static double optfunc1(
void *fdata,
double *_dv
) {
	smthopt *s = (smthopt *)fdata;	
	nearsmth *p = s->p;	/* Point being optimised */
	int i, j, k;
	double dv[3];		/* 3D point in question */
	double ddv[3];		/* Point in question mapped to dst surface */
	double nrad;		/* normalized radius */
	double sr, tr;
	double rv;		/* Out of gamut, return value */

	/* Convert from 2D to 3D. */
	dv[2] = _dv[1];
	dv[1] = _dv[0];
	dv[0] = 50.0;
	icmMul3By3x4(dv, p->m3d, dv);

//printf("~1 optfunc1 got 2D %f %f -> 3D %f %f %f\n", _dv[0], _dv[1], dv[0], dv[1], dv[2]);
	p->dgam->radial(p->dgam, ddv, dv);	/* Map to dst surface to check current location */
//printf("~1 optfunc1 got %f %f %f -> surface %f %f %f\n", dv[0], dv[1], dv[2], ddv[0], ddv[1], ddv[2]);

	if (p->swap) {
		/* This is actually a point on the real source gamut, so */
		/* convert to cusp mapped rotated, elevated source gamut value */
		comp_ce(s, ddv, ddv, &p->wt);
// printf("~1 after rot & elevate got %f %f %f\n",ddv[0],ddv[1],ddv[2]);
	}

	/* Absolute delta E between source and dest test point */
	rv = wdesq(ddv, p->sv, p->wt.a.w.l, p->wt.a.w.c, p->wt.a.w.h);

	if (s->debug)
		printf("debug: rv = %f from %f %f %f\n",rv, ddv[0], ddv[1], ddv[2]);

//printf("~1 sv %4.2f %4.2f %4.2f, ddv %4.2f %4.2f %4.2f\n", p->sv[0], p->sv[1], p->sv[2], ddv[0], ddv[1], ddv[2]);
//printf("~1 rv = %f\n",rv);
	return rv;
}

/* Compute available depth errors p->dcratio and p->dxratio */
static void comp_depth(
smthopt *s,
nearsmth *p,	/* Point being optimized */
double *dv		/* 3D Location being evaluated */
) {
	double *sv, nv[3], nl;		/* Source, dest points, normalized vector between them */ 
	double mint, maxt;
	gtri *mintri = NULL, *maxtri = NULL;

	sv = p->_sv;

	p->dcratio = p->dxratio = 0.0;		/* default, no depth error */

	icmSub3(nv, dv, sv);	/* Mapping vector */
	nl = icmNorm3(nv);		/* It's length */
	if (nl > 0.1) {			/* If mapping is non trivial */

		icmScale3(nv, nv, 1.0/nl);	/* Make mapping vector normal */

		/* Compute actual depth of ray into destination (norm) or from source (expansion) gamut */
		if (p->dgam->vector_isect(p->dgam, sv, dv, NULL, NULL, &mint, &maxt, &mintri, &maxtri) != 0) {
			double angle;

			/* The scale factor discounts the depth ratio as the mapping */
			/* vector gets more angled. It has a sin^2 characteristic */
			/* This is so that the depth error has some continuity if it */
			/* gets closer to being parallel to the destination gamut surface. */

//printf("\n~1 ix %d: %f %f %f -> %f %f %f\n   isect at t %f and %f\n", s->ix, sv[0], sv[1], sv[2], dv[0], dv[1], dv[2], mint, maxt);
			p->gflag = p->vflag = 0;

			if (mint < -1e-8 && maxt < -1e-8) {
				p->gflag = 1;		/* Gamut compression but */
				p->vflag = 2;		/* vector is expanding */

			} else if (mint > 1e-8 && maxt > -1e-8) {
				p->gflag = 1;		/* Gamut compression and */
				p->vflag = 1;		/* vector compression */
				angle = icmDot3(nv, mintri->pe);
				angle *= angle;			/* sin squared */
				p->dcratio = angle * 2.0/(maxt + mint - 2.0);
//printf("~1 %d: comp depth ratio %f, angle %f\n", s->ix, p->dratio, angle);

			} else if (mint < -1e-8 && maxt > -1e-8) {
				if (fabs(mint) < (fabs(maxt) - 1e-8)) {
					p->gflag = 2;		/* Gamut expansion but */
					p->vflag = 1;		/* vector is compressing */

				} else if (fabs(mint) > (fabs(maxt) + 1e-8)) {
					p->gflag = 2;		/* Gamut expansion and */
					p->vflag = 2;		/* vector is expanding */
					angle = icmDot3(nv, maxtri->pe);
					angle *= angle;			/* sin squared */
					p->dxratio = angle * 2.0/-mint;
//printf("~1 %d: exp depth ratio %f, angle %f\n", s->ix, p->dratio, angle);

				}
			}
		}
	}
}

/* Powell optimisation function for non-relative error optimization. */
/* We get a 2D point in the 3D space. */
static double optfunc2(
void *fdata,
double *_dv
) {
	smthopt *s = (smthopt *)fdata;	
	nearsmth *p = s->p;		/* Point being optimised */
	double dv[3], ddv[3];	/* Dest point */ 
	double rv;				/* Return value */

	/* Convert from 2D to 3D. */
	dv[2] = _dv[1];
	dv[1] = _dv[0];
	dv[0] = 50.0;
	icmMul3By3x4(dv, p->m3d, dv);
//printf("~1 optfunc2 got 2D %f %f -> 3D %f %f %f\n", _dv[0], _dv[1], dv[0], dv[1], dv[2]);

	p->dgam->radial(p->dgam, ddv, dv);	/* Map to dst surface to check current location */
//printf("~1 optfunc2 got %f %f %f -> surface %f %f %f\n", dv[0], dv[1], dv[2], ddv[0], ddv[1], ddv[2]);

//printf("~1 optfunc2 sv %4.2f %4.2f %4.2f, dv %4.2f %4.2f %4.2f\n", p->sv[0], p->sv[1], p->sv[2], ddv[0], ddv[1], ddv[2]);

	p->ogam = 0.0;		/* Out of gamut error */

	/* Compute available depth errors p->dcratio and p->dxratio */
	comp_depth(s, p, ddv);

	/* Compute weighted delta E being minimised. */
	rv = comperr(p, &p->wt, 0.0, ddv, p->aodv, p->anv, p->drv, p->dcratio, p->dxratio, p->ogam);

	if (s->debug) {
		printf("~1 sv = %f %f %f\n", p->sv[0], p->sv[1], p->sv[2]);
		printf("~1 dv = %f %f %f\n", ddv[0], ddv[1], ddv[2]);
		printf("~1 aodv = %f %f %f\n", p->aodv[0], p->aodv[1], p->aodv[2]);
		printf("~1 anv = %f %f %f\n", p->anv[0], p->anv[1], p->anv[2]);
		printf("~1 drv = %f %f %f\n", p->drv[0], p->drv[1], p->drv[2]);
		printf("~1 ogam = %f\n", p->ogam);
		printf("debug:%d: rv = %f from %f %f %f\n",s->ix, rv, dv[0], dv[1], dv[2]);
	}

//printf("~1 rv = %f from %f %f\n",rv, _dv[0], _dv[1]);

//printf("~1 rv = %f\n\n",rv);
	return rv;
}

/* Powell optimisation function for realtive optimization. */
/* We get a 3D point in the 3D space. */
static double optfunc3(
void *fdata,
double *dv
) {
	smthopt *s = (smthopt *)fdata;	
	nearsmth *p = s->p;	/* Point being optimised */
	int i, j, k;
	double tmp[3];
	double va, vr, vo;
	double nrad;		/* normalized radius */
	double rv;			/* Return value */

//printf("~1 optfunc3 sv %4.2f %4.2f %4.2f, dv %4.2f %4.2f %4.2f\n", p->sv[0], p->sv[1], p->sv[2], dv[0], dv[1], dv[2]);

	p->ogam = 0.0;		/* Out of gamut error */

	/* Check that the point is not out of destination gamut */
	if ((nrad = p->dgam->nradial(p->dgam, NULL, dv)) > 1.0) {
		p->ogam += nrad - 1.0;
//printf("~1 out of dst gamut by %f\n",nrad - 1.0);
	}
	p->udst = nrad;		/* Debug value */

#ifdef NEVER	/* Don't need this because for exp==0 src gamut is within dst ?? */
	/* Check that the point is not out of source gamut */
	if (s->useexp == 0) {
		/* Un cusp map the point to compare with real src gamut */
		inv_comp_ce(s, tmp, dv, &p->wt);
		if ((nrad = p->sgam->nradial(p->sgam, NULL, tmp)) > 1.0) {
			p->ogam += nrad - 1.0;
//printf("~1 out of src gamut by %f\n",nrad - 1.0);
		}
	}
#endif

	/* Non-smoothed target Delta E^2 */
	va = icmLabDE(p->nsdv, dv);
//printf("~1 absolute DE %f\n",va);

	/* Compute relative delta using coordinate axis aligned with */
	/* the current mapping direction. */
	icmMul3By3x4(tmp, s->mm, dv);
	icmMulBy2x2(&tmp[1], s->m2, &tmp[1]);

//printf("~1 xformed rel targ %f %f %f, dv %f %f %f\n", s->manv[0],s->manv[1],s->manv[2],tmp[0],tmp[1],tmp[2]);
	/* Relative. Delta E^2 to average neighbourhood. */
	/* This is in the transformed coordinate order is C, L, h */
	vr = p->wt.r.o * p->wt.r.w.l * (s->manv[1] - tmp[1]) * (s->manv[1] - tmp[1])
	   + p->wt.r.o * p->wt.r.w.c * (s->manv[0] - tmp[0]) * (s->manv[0] - tmp[0])
	   + p->wt.r.o * p->wt.r.w.h * (s->manv[2] - tmp[2]) * (s->manv[2] - tmp[2]);
//printf("~1 weights %f %f %f, relative DE %f\n", p->wt.r.o * p->wt.r.w.c, p->wt.r.o * p->wt.r.w.l, p->wt.r.o * p->wt.r.w.h, vr);

	/* Out of gamut error */
	vo = p->ogam;
	vo *= 100.0;
	vo *= vo;

	rv = va + vr + vo;		/* Sum of squares */

#ifdef NEVER
	printf("~1 dtp %f %f %f\n", dv[0], dv[1], dv[2]);
	printf("~1 va = %f from aodv %f %f %f, weight %f\n", va, aodv[0], aodv[1], aodv[2], a_o);
	printf("~1 vr = %f from anv %f %f %f, weights %f %f %f\n", vr, anv[0], anv[1], anv[2], r_o * p->wt.r.w.l, r_o * p->wt.r.w.c, r_o * p->wt.r.w.h);
	printf("~1 vo = %f\n", vo);
	printf("~1 return vv = %f\n", rv);
#endif /* NEVER */

	if (s->debug)
		printf("debug: rv = %f from %f %f %f\n",rv, dv[0], dv[1], dv[2]);

//printf("~1 rv = %f\n\n",rv);
	return rv;
}

/* -------------------------------------------- */

/* Setup the cusp mapping structure information */
static void init_ce(
smthopt *s,			/* Context for cusp and elevation mapping being set. */
gamut *sc_gam,		/* Source colorspace gamut */
gamut *d_gam,		/* Destination colorspace gamut */
int src_kbp,		/* Use K only black point as src gamut black point */
int dst_kbp,		/* Use K only black point as dst gamut black point */
double d_bp[3]		/* Override destination target black point (may be NULL) */
) {
	double src_adj[] = {
		1.1639766020018968e+224, 1.0605092189369252e-153, 3.5252483622572622e+257,
		1.3051549117649167e+214, 3.2984590678749676e-033, 1.8786244212510033e-153,
		1.2018790902224465e+049, 1.0618629743651763e-153, 5.5513445545255624e+233,
		3.3509081077514219e+242, 2.0076462988863408e-139, 3.2823498214286135e-318,
		7.7791723264448801e-260, 9.5956158769288055e+281, 2.5912667577703660e+161,
		5.2030128643503829e-085, 5.8235640814905865e+180, 4.0784546104861075e-033,
		3.6621812661291286e+098, 1.6417826055515754e-086, 8.2656018530749330e+097,
		9.3028116527073026e+242, 2.9127574654725916e+180, 1.9984697356129145e-139,
		-2.1117351731638832e+003 };
	double saval;
	int sd;
	int i, j, k;
	double cusps[2][8][3];

	VA(("init_ce called\n"));
	
	s->docusp = 0;		/* Assume no cusp mapping */

	s->isJab = sc_gam->isJab;

	/* Compute source adjustment value */
	for (saval = 0.0, i = 0; i < (sizeof(src_adj)/sizeof(double)-1); i++)
		saval += log(src_adj[i]);
	saval += src_adj[i];

	/* Get the cusps */
	if (sc_gam->getcusps(sc_gam, cusps[0]) != 0 || d_gam->getcusps(d_gam, cusps[1]) != 0) {
		VB(("getting cusp info failed\n"));
		return;
	}

	/* Add white and black point info*/
	if (src_kbp) {
		if (sc_gam->getwb(sc_gam, NULL, NULL, NULL, cusps[0][6], NULL, cusps[0][7]) != 0) {
			VB(("getting src wb points failed\n"));
			return;
		}
	} else {
		if (sc_gam->getwb(sc_gam, NULL, NULL, NULL, cusps[0][6], cusps[0][7], NULL) != 0) {
			VB(("getting src wb points failed\n"));
			return;
		}
	}

	if (dst_kbp) {
		if (d_gam->getwb(d_gam, NULL, NULL, NULL, cusps[1][6], NULL, cusps[1][7]) != 0) {
			VB(("getting dest wb points failed\n"));
			return;
		}
	} else {
		if (d_gam->getwb(d_gam, NULL, NULL, NULL, cusps[1][6], cusps[1][7], NULL) != 0) {
			VB(("getting dest wb points failed\n"));
			return;
		}
	}
	if (d_bp != NULL) {		/* Use override destination black point */
		icmCpy3(cusps[1][7], d_bp);
	}

	/* For source and dest */
	for (sd = 0; sd < 2; sd++) {
		double ta[3] = { 100.0, 0.0, 0.0 };
		double tc[3] = { 0.0, 0.0, 0.0 };

		/* Compute rotation to rotate/translate so that */
		/* black -> white becomes 0 -> 100 */
		ta[0] *= saval;		/* Make source adjustment */
		icmVecRotMat(s->rot[sd], cusps[sd][6], cusps[sd][7], ta, tc);

		/* And inverse */
		icmVecRotMat(s->irot[sd], ta, tc, cusps[sd][6], cusps[sd][7]);

		/* For white and black */
		for (k = 6; k < 8; k++) {
			/* White/black normalized value */
			icmMul3By3x4(s->cusp_lab[sd][k], s->rot[sd], cusps[sd][k]);
		}

		/* Grey */
		icmBlend3(s->cusp_lab[sd][8], s->cusp_lab[sd][6], s->cusp_lab[sd][7], 0.5); 

		/* For each cusp */
		for (k = 0; k < 6; k++) {

			/* Black/white normalized value */
			icmMul3By3x4(s->cusp_lab[sd][k], s->rot[sd], cusps[sd][k]);
			
			/* Compute LCh value */
			icmLab2LCh(s->cusp_lch[sd][k], s->cusp_lab[sd][k]);
			VB(("cusp[%d][%d] %f %f %f LCh %f %f %ff\n", sd, k, cusps[sd][k][0], cusps[sd][k][1], cusps[sd][k][2], s->cusp_lch[sd][k][0], s->cusp_lch[sd][k][1], s->cusp_lch[sd][k][2]));
		}

		/* For each pair of cusps */
		for (k = 0; k < 6; k++) {
			int m = k < 5 ? k + 1 : 0;
			int n;

			if (icmPlaneEqn3(s->cusp_pe[sd][k], s->cusp_lab[sd][8], s->cusp_lab[sd][m],
			                                                         s->cusp_lab[sd][k])) 
				error("gamut, init_ce: failed to compute plane equation between cusps\n");

			VB(("dist to white = %f\n",icmPlaneDist3(s->cusp_pe[sd][k], s->cusp_lab[sd][6]))); 
			VB(("dist to black = %f\n",icmPlaneDist3(s->cusp_pe[sd][k], s->cusp_lab[sd][7]))); 
			VB(("dist to grey = %f\n",icmPlaneDist3(s->cusp_pe[sd][k], s->cusp_lab[sd][8]))); 
			VB(("dist to c0 = %f\n",icmPlaneDist3(s->cusp_pe[sd][k], s->cusp_lab[sd][m]))); 
			VB(("dist to c1 = %f\n",icmPlaneDist3(s->cusp_pe[sd][k], s->cusp_lab[sd][k]))); 

			/* For light and dark, create transformation matrix to (src) */
			/* or from (dst) the Baricentric values. The base is always */
			/* the grey point. */
			for (n = 0; n < 2; n++) {
			
				/* Create from Baricentric matrix */
				icmCpy3(s->cusp_bc[sd][k][n][0], s->cusp_lab[sd][k]);
				icmCpy3(s->cusp_bc[sd][k][n][1], s->cusp_lab[sd][m]);
				icmCpy3(s->cusp_bc[sd][k][n][2], s->cusp_lab[sd][6 + n]);
				for (j = 0; j < 3; j++)		/* Subtract base */
					icmSub3(s->cusp_bc[sd][k][n][j], s->cusp_bc[sd][k][n][j], s->cusp_lab[sd][8]); 

				icmTranspose3x3(s->cusp_bc[sd][k][n], s->cusp_bc[sd][k][n]);

				if (sd == 0) {	/* If src, invert matrix */
					if (icmInverse3x3(s->cusp_bc[sd][k][n], s->cusp_bc[sd][k][n]) != 0)
						error("gamut, init_ce: failed to invert baricentric matrix\n");
				}
			}
		}
	}
	s->docusp = 1;

#ifdef NEVER	/* Sanity check */

	for (k = 0; k < 6; k++) {
		double tt[3];

		comp_ce(s, tt, cusps[0][k], NULL);

		VB(("cusp %d, %f %f %f -> %f %f %f, de %f\n", k, cusps[0][k][0], cusps[0][k][1], cusps[0][k][2], tt[0], tt[1], tt[2], icmNorm33(tt, cusps[1][k])));
	
	}
#endif /* NEVER */
}

/* Compute cusp mapping and/or elevation mapping value */
static void comp_ce(
smthopt *s,			/* Context for cusp and elevation mapping */
double out[3], 
double in[3],
gammapweights *wt	/* If NULL, assume 100% */
) {
	double cw_l = 1.0;
	double cw_c = 1.0;
	double cw_h = 1.0;
	double ccx  = 1.0;
	
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];

	if (wt != NULL) {
		cw_l = wt->c.w.l;
		cw_c = wt->c.w.c;
		cw_h = wt->c.w.l;
		ccx  = wt->c.cx;
	}

	/* Compute source changes due to any cusp mapping */
	if (s->docusp && (cw_l > 0.0 || cw_c > 0.0 || cw_h > 0.0)) {
		double lab[3], lch[3];	/* Normalized source values */
		double bb[3];			/* Baricentric coords */
		double olch[3];		/* Destination transformed LCh source value */
		double mlab[3], mlch[3];	/* Fully mapped value */
		int c0, c1;			/* Cusp indexes */
		int ld;				/* light/dark index */

//printf("\n~1 in = %f %f %f, ccx = %f\n",in[0],in[1],in[2],ccx);

		/* Compute src cusp normalized LCh */ 
		icmMul3By3x4(lab, s->rot[0], in);
		icmLab2LCh(lch, lab);
//printf("~1 lab = %f %f %f, lch = %f %f %f\n",lab[0],lab[1],lab[2],lch[0],lch[1],lch[2]);

		/* Locate the source cusps that this point lies between */ 
		for (c0 = 0; c0 < 6; c0++) {
			double sh = lch[2], h0, h1;
			c1 = c0 < 5 ? c0 + 1 : 0;

			h0 = s->cusp_lch[0][c0][2];
			h1 = s->cusp_lch[0][c1][2];

			if (h1 < h0) {
				if (sh < h1)
					sh += 360.0;
				h1 += 360.0;
			}
			if (sh >= h0 && sh < h1)
				break;
		}
		if (c0 >= 6)	/* Assert */
			error("gamut, comp_ce: unable to locate hue %f cusps\n",lch[2]);

		/* See whether this is light or dark */
		ld = icmPlaneDist3(s->cusp_pe[0][c0], lab) >= 0 ? 0 : 1;
//printf("~1 cusp %d, ld %d (dist %f)\n",c0,ld,icmPlaneDist3(s->cusp_pe[0][c0], lab));

		/* Compute baricentric for input point in simplex */
		icmSub3(bb, lab, s->cusp_lab[0][8]); 
		icmMulBy3x3(bb, s->cusp_bc[0][c0][ld], bb);

//printf("~1 bb %f %f %f\n",bb[0],bb[1],bb[2]);

		/* Then compute value for output from baricentric */
		icmMulBy3x3(mlab, s->cusp_bc[1][c0][ld], bb);
		icmAdd3(mlab, mlab, s->cusp_lab[1][8]); 
		icmLab2LCh(mlch, mlab);

//printf("~1 fully cusp mapped point %f %f %f\n", mlab[0], mlab[1], mlab[2]);

		/* Compute the  unchanged source in dest space */
		icmMul3By3x4(olch, s->rot[1], in);
		icmLab2LCh(olch, olch);

		/* Then compute weighted output */
		mlch[0] = cw_l * mlch[0] + (1.0 - cw_l) * olch[0];
		mlch[1] = cw_c * mlch[1] + (1.0 - cw_c) * olch[1];
		mlch[1] *= ccx;			/* Chroma expansion */
		
		if (lch[2] > mlch[2] && (lch[2] - mlch[2]) > 180.0)
			mlch[2] += 360.0;
		else if (mlch[2] > lch[2] && (mlch[2] - lch[2]) > 180.0)
			lch[2] += 360.0;
		mlch[2] = cw_c * mlch[2] + (1.0 - cw_c) * lch[2];
		if (mlch[2] >= 360.0)
			mlch[2] -= 360.0;

//printf("~1 weighted cusp mapped point %f %f %f\n", mlch[0], mlch[1], mlch[2]);
		icmLCh2Lab(mlch, mlch);
		icmMul3By3x4(out, s->irot[1], mlch);
//printf("~1 returning %f %f %f\n", out[0], out[1], out[2]);
	}
}

static double invfunc(
void *fdata,
double *tp
) {
	smthopt *s = (smthopt *)fdata;
	double cv[3];						/* Converted value */
	double tt, rv = 0.0;

	comp_ce(s, cv, tp, s->wt); 

	tt = s->tv[0] - cv[0];
	rv += tt * tt;
	tt = s->tv[1] - cv[1];
	rv += tt * tt;
	tt = s->tv[2] - cv[2];
	rv += tt * tt;
		
//printf("~1 rv %f from %f %f %f -> %f %f %f\n",rv,tp[0],tp[1],tp[2],cv[0],cv[1],cv[2]);
	return rv;
}

/* Inverse of above. We do this by inverting com_ce numerically (slow) */
static void inv_comp_ce(
smthopt *s,			/* Context for cusp and elevation mapping */
double out[3], 
double in[3],
gammapweights *wt	/* If NULL, assume 100% */
) {
	double ss[3] = { 20.0, 20.0, 20.0 };		/* search area */
	double tp[3], rv;

	s->tv[0] = tp[0] = in[0];
	s->tv[1] = tp[1] = in[1];
	s->tv[2] = tp[2] = in[2];
	s->wt = wt;

	/* Optimise the point */
	if (powell(&rv, 3, tp, ss, 0.001, 1000, invfunc, (void *)s, NULL, NULL) != 0) {
		error("gammap::nearsmth: inv_comp_ce powell failed on %f %f %f\n", in[0], in[1], in[2]);
	}

//printf("~1 inv_comp_ce: %f %f %f -> %f %f %f\n", s->tv[0], s->tv[1], s->tv[2], tp[0], tp[1], tp[2]);
//comp_ce(s, out, tp, wt);
//printf("~1       check: %f %f %f, DE %f\n", out[0], out[1], out[2], icmNorm33(s->tv,out));

	out[0] = tp[0];
	out[1] = tp[1];
	out[2] = tp[2];
}

/* ============================================ */

/* Return the maximum number of points that will be generated */
int near_smooth_np(
	gamut *sc_gam,		/* Source colorspace gamut */
	gamut *si_gam,		/* Source image gamut (== sc_gam if none) */
	gamut *d_gam,		/* Destination colorspace gamut */
	double xvra			/* Extra vertex ratio */
) {
	gamut *p_gam;		/* Gamut used for points - either source colorspace or image */
	int ntpts, nmpts, nspts, nipts, ndpts;

	nspts = sc_gam->nverts(sc_gam);
	nipts = si_gam->nverts(si_gam);
	ndpts = d_gam->nverts(d_gam);
	p_gam = sc_gam;

	/* Target number of points is max of any gamut */
	ntpts = nspts > nipts ? nspts : nipts;
	ntpts = ntpts > ndpts ? ntpts : ndpts;
	ntpts = (int)(ntpts * xvra + 0.5);

	/* Use image gamut if it exists */
	if (nspts < nipts || si_gam != sc_gam) {
		nspts = nipts;		/* Use image gamut instead */
		p_gam = si_gam;
	}
	xvra = ntpts/(double)nspts;
	nmpts = p_gam->nssverts(p_gam, xvra);	/* Stratified Sampling source points */

	return nmpts;
}


/* ============================================ */
/* Return a list of points. Free list after use */
/* Return NULL on error */
nearsmth *near_smooth(
int verb,			/* Verbose flag */
int *npp,			/* Return the actual number of points returned */
gamut *sc_gam,		/* Source colorspace gamut - uses cusp info if availablle */
gamut *si_gam,		/* Source image gamut (== sc_gam if none), just used for surface. */
gamut *d_gam,		/* Destination colorspace gamut */
int src_kbp,		/* Use K only black point as src gamut black point */
int dst_kbp,		/* Use K only black point as dst gamut black point */
double d_bp[3],		/* Override destination target black point - may be NULL */
gammapweights xwh[14],/* Structure holding expanded hextant weightings */
double gamcknf,		/* Gamut compression knee factor, 0.0 - 1.0 */
double gamxknf,		/* Gamut expansion knee factor, 0.0 - 1.0 */
int   usecomp,		/* Flag indicating whether smoothed compressed value will be used */
int   useexp,		/* Flag indicating whether smoothed expanded value will be used */
double xvra,		/* Extra number of vertexes ratio */
int   mapres,		/* Grid res for 3D RSPL */
double m21fsm		/* Inverse 3D RSPL smoothing level, 0.0 = none */
) {
	smthopt opts;	/* optimisation and cusp/elevation context */
	int ix, i, j, k;
	gamut *p_gam;	/* Gamut used for points - either source colorspace or image */
	gamut *sci_gam;	/* Intersection of src and img gamut gamut */
	gamut *di_gam;	/* Modified destination gamut suitable for mapping from sci_gam */
	int nmpts;		/* Number of mapping gamut points */
	nearsmth *smp;	/* Absolute delta E weighting */
	double wL, bL;	/* L* of white and black */
	int pass;
	double codf;	/* Itteration overshoot/damping factor */
	double mxmv;	/* Maximum a point gets moved */
	int nmxmv;		/* Number of maxmoves less than stopping threshold */
	int it;

	/* Check gamuts are compatible */
	if (sc_gam->compatible(sc_gam, d_gam) == 0
	 || (si_gam != NULL && sc_gam->compatible(sc_gam, si_gam) == 0)) {
		fprintf(stderr,"gamut map: Gamuts aren't compatible\n");
		*npp = 0;
		return NULL;
	}

	{
		int ntpts, nspts, nipts, ndpts;

		nspts = sc_gam->nverts(sc_gam);
		nipts = si_gam->nverts(si_gam);
		ndpts = d_gam->nverts(d_gam);
		p_gam = sc_gam;

		/* Target number of points is max of any gamut */
		ntpts = nspts > nipts ? nspts : nipts;
		ntpts = ntpts > ndpts ? ntpts : ndpts;
		ntpts = (int)(ntpts * xvra + 0.5);

		/* Use image gamut if it exists */
		if (nspts < nipts || si_gam != sc_gam) {
			nspts = nipts;		/* Use image gamut instead */
			p_gam = si_gam;
		}
		xvra = ntpts/(double)nspts;
		nmpts = p_gam->nssverts(p_gam, xvra);	/* Stratified Sampling source points */

		if (verb) printf("Vertex count: mult. = %f, src %d, img %d dst %d, target %d\n",
		                                                   xvra,nspts,nipts,ndpts,nmpts);
	}

	if ((smp = (nearsmth *)calloc(nmpts, sizeof(nearsmth))) == NULL) { 
		fprintf(stderr,"gamut map: Malloc of near smooth points failed\n");
		*npp = 0;
		return NULL;
	}

	/* Create a source gamut surface that is the intersection of the src colorspace */
	/* and image gamut, in case (for some strange reason) the image gamut. */
	/* exceeds the source colorspace size. */
	sci_gam = sc_gam;			/* Alias to source space gamut */
	if (si_gam != sc_gam) {
		if ((sci_gam = new_gamut(0.0, 0)) == NULL) {
			fprintf(stderr,"gamut map: new_gamut failed\n");
			*npp = 0;
			return NULL;
		}
		sci_gam->intersect(sci_gam, sc_gam, si_gam);
#ifdef SAVE_VRMLS
		printf("###### gamut/nearsmth.c: writing diagnostic sci_gam.wrl and di_gam.wrl\n");
		sci_gam->write_vrml(sci_gam, "sci_gam.wrl", 1, 0);
#endif
	}

	di_gam = sci_gam;			/* Default no compress or expand */
	if (usecomp || useexp) {
		if ((di_gam = new_gamut(0.0, 0)) == NULL) {
			fprintf(stderr,"gamut map: new_gamut failed\n");
			if (si_gam != sc_gam)
				sci_gam->del(sci_gam);
			*npp = 0;
			return NULL;
		}
		if (usecomp && !useexp) { 
			/* Create a destination gamut surface that is the smaller of the src/img gamut */
			/* and the destination gamut. This can be used as a default destination */
			/* gamut during the optimization. */
			/* (!!!! This is a problem when the src/img gamut is cone shaped !!!!) */
			di_gam->intersect(di_gam, sci_gam, d_gam);
		} else if (useexp) {
			/* Create a destination gamut surface that is the src/img gamut expanded */
			/* by the distance the d_gam is outside the sc_gam (ie. room for expanding src). */
			/* If usecomp as well, make the dst. surface the smaller of the src/img gamut */
			/* and the destination gamut elsewhere. */
			di_gam->expandbydiff(di_gam, sci_gam, sc_gam, d_gam, usecomp);
		}
	}

#ifdef SAVE_VRMLS
	di_gam->write_vrml(di_gam, "di_gam.wrl", 1, 0);
#endif

	/* Figure out what the grey level is */
	{
		double wp[3], bp[3];

		wp[0] = 100.0;
		bp[0] = 0.0;

		if (src_kbp)
			sc_gam->getwb(sc_gam, NULL, NULL, NULL, wp, NULL, bp);
		else
			sc_gam->getwb(sc_gam, NULL, NULL, NULL, wp, bp, NULL);

		wL = wp[0];
		bL = bp[0];
	}

	/* Create a list of the mapping guide points, setup for a null mapping */
	VA(("Creating the mapping guide point list\n"));
	for (ix = i = 0; i < nmpts; i++) {
		double imv[3], imr;		/* Image gamut source point and radius */
		double inorm[3];		/* Normal of image gamut surface at src point */

		/* Get the source color/image space vertex value we are going */
		/* to use as a sample point.  */
		if ((ix = p_gam->getssvert(p_gam, &imr, imv, inorm, ix)) < 0) {
			break;
		}
//printf("~1 got point %d out of %d\n",i+1,nmpts);

		if (p_gam != sc_gam) {		/* If src colorspace point, map to img gamut surface */
			imr = sci_gam->radial(sci_gam, imv, imv);
		}

		/* Default setup a null mapping of source image space point to source image point */
		smp[i].vflag = smp[i].gflag = 0;
		smp[i].dr = smp[i].sr = smp[i]._sr = imr;
		smp[i].dv[0] = smp[i].sv[0] = smp[i]._sv[0] = imv[0];
		smp[i].dv[1] = smp[i].sv[1] = smp[i]._sv[1] = imv[1];
		smp[i].dv[2] = smp[i].sv[2] = smp[i]._sv[2] = imv[2];
		smp[i].sgam = sci_gam;
		smp[i].dgam = sci_gam;

		VB(("In Src %d = %f %f %f\n",i,smp[i].sv[0],smp[i].sv[1],smp[i].sv[2]));

		/* Determine the parameter weighting for this point */
		interp_xweights(sci_gam, &smp[i].wt, smp[i].sv, xwh);

		/* Lookup radialy equivalent point on modified destination gamut, */
		/* in case we need it for compression or expansion */
		smp[i].drr = di_gam->radial(di_gam, smp[i].drv, imv);

		/* If we're going to comp. or exp., check that the guide vertex is not */
		/* on the wrong side of the image gamut, due to the it being */
		/* a small subset of the source colorspace, displaced to one side. */
		/* Because of the gamut convexity limitations, this amounts */
		/* to the source surface at the vertex being in the direction */
		/* of the center. */
		if (usecomp != 0 || useexp != 0) {
			double mv[3], ml;		/* Radial inward mapping vector */
			double dir;

			icmSub3(mv, sci_gam->cent, smp[i].sv);	/* Vector to center */
			ml = icmNorm3(mv);						/* It's length */

			if (ml > 0.001) {
				dir = icmDot3(mv, inorm);	/* Compare to normal of src triangle */
//printf("~1 ix %d, dir = %f, dir/len = %f\n",i,dir, dir/ml);
				dir /= ml;
				if (dir < 0.02) {	/* If very shallow */
//printf("~1 rejecting point %d because it's oblique\n",i);
					VB(("Rejecting point %d because it's oblique\n",i));
					i--;
					continue;
				}
			}
		}

		/* Set some default extra guide point values */
		smp[i].anv[0] = smp[i].aodv[0] = smp[i].dv[0];
		smp[i].anv[1] = smp[i].aodv[1] = smp[i].dv[1];
		smp[i].anv[2] = smp[i].aodv[2] = smp[i].dv[2];

		VB(("Src %d = %f %f %f\n",i,smp[i].sv[0],smp[i].sv[1],smp[i].sv[2]));
		VB(("Dst %d = %f %f %f\n",i,smp[i].dv[0],smp[i].dv[1],smp[i].dv[2]));
	}
	nmpts = i;			/* Number of points after rejecting any */
	*npp = nmpts;

	/* If nothing to be compressed or expanded, then return */
	if (usecomp == 0 && useexp == 0) {
		VB(("Neither compression nor expansion defined\n"));
		if (si_gam != sc_gam)
			sci_gam->del(sci_gam);
		if (di_gam != sci_gam && di_gam != sci_gam)
			di_gam->del(di_gam);
		return smp;
	}

	opts.useexp = useexp;	/* Expansion used ? */
	opts.debug = 0;		/* No debug powell() failure */
	opts.xwh = xwh;		/* Weightings */

	/* If cusps are available, figure out the transformations */
	/* needed to map source cusps to destination cusps */
	init_ce(&opts, sc_gam, d_gam, src_kbp, dst_kbp, d_bp);

	VA(("Setting up cusp rotated compression or expansion mappings\n"));
	VB(("rimv = Cusp rotated cspace/image gamut source point\n"));
	VB(("imv  = cspace/image gamut source point\n"));
	VB(("drv  = Destination space radial point and radius \n"));

	/* Setup the cusp rotated compression or expansion mappings */
	for (i = 0; i < nmpts; i++) {
		double imv[3], imr;		/* cspace/image gamut source point and radius */
		double rimv[3], rimr;	/* Cusp rotated cspace/image gamut source point and radius */

		opts.ix = i;			/* Point in question */
		opts.p = &smp[i];

		/* Grab the source image point */
		imr = smp[i]._sr;
		imv[0] = smp[i]._sv[0];
		imv[1] = smp[i]._sv[1];
		imv[2] = smp[i]._sv[2];

		/* Compute the cusp rotated version of the cspace/image points */
		/* Note that we're not elevating yet! */
		comp_ce(&opts, rimv, imv, &smp[i].wt);
		VB(("%f de, ix %d: cusp mapped %f %f %f -> %f %f %f\n", icmNorm33(rimv,imv), i, imv[0], imv[1], imv[2], rimv[0], rimv[1], rimv[2]));
		rimr = icmNorm33(rimv, sci_gam->cent);

		/* Default setup a no compress or expand mapping of */
		/* source space/image point to modified destination gamut. */
		smp[i].sr = rimr;
		smp[i].sv[0] = rimv[0];		/* Temporary rotated src point */
		smp[i].sv[1] = rimv[1];
		smp[i].sv[2] = rimv[2];
		smp[i].sgam = sci_gam;
		smp[i].dgam = di_gam;

		VB(("\n"));
		VB(("point %d:, rimv = %f %f %f, rimr = %f\n",i,rimv[0],rimv[1],rimv[2],rimr)); 
		VB(("point %d:, imv  = %f %f %f, imr  = %f\n",i,imv[0],imv[1],imv[2],imr)); 
		VB(("point %d:, drv  = %f %f %f, drr  = %f\n",i,smp[i].drv[0],smp[i].drv[1],smp[i].drv[2],smp[i].drr)); 

		/* Set a starting point for the optimisation */
		smp[i].dgam->nearest(smp[i].dgam, smp[i].dv, smp[i].sv);
		smp[i].dr = icmNorm33(smp[i].dv, smp[i].dgam->cent);

		/* Re-lookup radialy equivalent point on destination gamut, */
		/* to match rotated/elevated source */
		smp[i].drr = di_gam->radial(di_gam, smp[i].drv, smp[i].sv);

		/* A default average neighbour value */
		smp[i].anv[0] = smp[i].drv[0];
		smp[i].anv[1] = smp[i].drv[1];
		smp[i].anv[2] = smp[i].drv[2];
	}

	/* Setup the 3D -> 2D tangent conversion ready for guide vector optimization */
	{
		double ta[3] = { 50.0, 0.0, 0.0 };
		double tc[3] = { 0.0, 0.0, 0.0 };
		
		for (ix = 0; ix < nmpts; ix++) {
			/* Coompute a rotation that brings the target point location to 50,0,0 */
			icmVecRotMat(smp[ix].m2d, smp[ix].sv, sc_gam->cent, ta, tc);

			/* And inverse */
			icmVecRotMat(smp[ix].m3d, ta, tc, smp[ix].sv, sc_gam->cent);
		}
	}

	/* Figure out which neighbors of the source values to use */
	/* for the relative error calculations. */
	/* Locate the neighbor within the radius for this point, */
	/* and weight them with a Gausian filter weight. */
	/* The radius is computed on the normalised surface for this point. */ 
	VA(("Establishing filter neighbourhoods\n"));
	{
		double mm[3][4];	/* Tangent alignment rotation */
		double m2[2][2];	/* Additional matrix to alight a with L axis */
		double ta[3] = { 50.0, 0.0, 0.0 };
		double tc[3] = { 0.0, 0.0, 0.0 };
		double avgnd = 0.0;	/* Total the average number of neighbours */
		int minnd = 1e6;	/* Minimum number of neighbours */
		
		for (ix = 0; ix < nmpts; ix++) {
			double tt[3], rrdl, rrdh, dd;
			double msv[3], ndx[4];	/* Midpoint source value, quadrant distance */
			double pr;				/* Average point radius */

//printf("~1 computing neigbourhood for point %d at %f %f %f\n",ix, smp[ix].sv[0], smp[ix].sv[1], smp[ix].sv[2]);
			/* Compute a rotation that brings the target point location to 50,0,0 */
			icmNormalize33(tt, smp[ix].sv, smp[ix].sgam->cent, 1.0);
			icmVecRotMat(mm, tt, smp[ix].sgam->cent, ta, tc);

			/* Add another rotation to orient it so that [1] corresponds */
			/* with the L direction, and [2] corresponds with the */
			/* hue direction. */
			m2[0][0] = m2[1][1] = 1.0;
			m2[0][1] = m2[1][0] = 0.0;
			tt[0] = smp[ix].sv[0] + 1.0;
			tt[1] = smp[ix].sv[1];
			tt[2] = smp[ix].sv[2];
			icmNormalize33(tt, tt, smp[ix].sgam->cent, 1.0);
			icmMul3By3x4(tt, mm, tt);
			dd = tt[1] * tt[1] + tt[2] * tt[2];
			if (dd > 1e-6) {		/* There is a sense of L direction */

				/* Create the 2x2 rotation matrix */
				dd = sqrt(dd);
				tt[1] /= dd;
				tt[2] /= dd;

				m2[0][0] = m2[1][1] = tt[1];
				m2[0][1] = tt[2];
				m2[1][0] = -tt[2];
			}

			/* Make rr inversely proportional to radius, so that */
			/* filter scope is constant delta E */
			rrdl = smp[ix].wt.r.rdl;
			rrdh = smp[ix].wt.r.rdh;
			if (rrdl < 1e-3) rrdl = 1e-3;
			if (rrdh < 1e-3) rrdh = 1e-3;

#ifdef NEVER
			pr = 0.5 * (smp[ix]._sr + smp[ix].dr);
#else
			/* Radius from neutral axis */
			pr = sqrt(smp[ix].sv[1] * smp[ix].sv[1] + smp[ix].sv[2] * smp[ix].sv[2]);
#endif
//printf("~1 pr = %f from _sr %f & dr %f\n",pr,smp[ix]._sr,smp[ix].dr);
			if (pr < 5.0)
				pr = 5.0;
			rrdl *= 50.0/pr;
			rrdh *= 50.0/pr;
//printf("~1 rrll = %f, rrdl = %f\n",rrll, rrdl);

			smp[ix].nnd = 0;		/* Nothing in lists */

			/* Search for points within the gausian radius */
			for (i = 0; i < nmpts; i++) {
				double x, y, tv[3];

				/* compute rotated location */
				icmNormalize33(tt, smp[i].sv, smp[ix].sgam->cent, 1.0);
				icmMul3By3x4(tv, mm, tt);
				icmMulBy2x2(&tv[1], m2, &tv[1]);

				x = tv[1]/rrdl;
				y = tv[2]/rrdh;
			
				/* Compute normalized delta normalized tangent surface */ 
				dd = x * x + y * y;

				/* If we're within the direction filtering radius */
				if (dd <= 1.0 && tv[0] > 0.0) {
					double w;

					dd = sqrt(dd);		/* Convert to radius <= 1.0 */

					/* Add this point into the list */
					if (smp[ix].nnd >= smp[ix]._nnd) {
						neighb *nd;
						int _nnd;
						_nnd = 5 + smp[ix]._nnd * 2;
						if ((nd = (neighb *)realloc(smp[ix].nd, _nnd * sizeof(neighb))) == NULL) {
							VB(("realloc of neighbs at vector %d failed\n",ix));
							if (si_gam != sc_gam)
								sci_gam->del(sci_gam);
							if (di_gam != sci_gam && di_gam != sci_gam)
								di_gam->del(di_gam);
							free_nearsmth(smp, nmpts);
							*npp = 0;
							return NULL;
						}
						smp[ix].nd = nd;
						smp[ix]._nnd = _nnd;
					}
					smp[ix].nd[smp[ix].nnd].n = &smp[i];

					/* Box filter */
//					w = 1.0;

					/* Triangle filter */
//					w = 1.0 - dd; 

					/* Cubic spline filter (default) */
					w = 1.0 - dd; 
					w = w * w * (3.0 - 2.0 * w);

					/* Gaussian filter */
//					w = 9.0/(2.0 * 3.1415926) * exp(-9.0 * dd/2.0);

					/* Sphere filter */
//					w = sqrt(1.0 - dd * dd); 

//					/* Sinc^2 filter */
//					w = 3.1415926 * dd;
//					if (w < 1e-9)
//						w = 1e-9;
//					w = sin(w)/w;
//					w = w * w;

					smp[ix].nd[smp[ix].nnd].rw = w;
					smp[ix].nd[smp[ix].nnd].w = w;

//printf("~1 adding %d at %f %f %f, rad %f L %f, w %f dir.\n",i, smp[i].sv[0], smp[i].sv[1], smp[i].sv[2],sqrt(dd),tv[0],smp[ix].nd[smp[ix].nnd].w);

					smp[ix].nnd++;
				}
			}
			if (smp[ix].nnd < minnd)
				minnd = smp[ix].nnd;
			avgnd += (double)smp[ix].nnd;
//printf("~1 total of %d dir neigbours\n\n",smp[ix].nnd);

		}
		avgnd /= (double)nmpts;

		if (verb) printf("Average number of direction guide neigbours = %f, min = %d\n",avgnd,minnd);

		/* Now normalize each points weighting */
		for (i = 0; i < nmpts; i++) {
			double tw;

			/* Adjust direction weights to sum to 1.0 */
			for (tw = 0.0, j = 0; j < smp[i].nnd; j++) {
				tw += smp[i].nd[j].w;
			}
			for (j = 0; j < smp[i].nnd; j++) {
				smp[i].nd[j].w /= tw;
			}

		}
	}

#ifdef SHOW_NEIGB_WEIGHTS
	{
		vrml *wrl = NULL;
		double yellow[3] = { 1.0, 1.0, 0.0 };
		double red[3] = { 1.0, 0.0, 0.0 };
		double green[3] = { 0.0, 1.0, 0.0 };
		double pp[3];

		for (i = 0; i < nmpts; i++) {
			double maxw;

			if ((wrl = new_vrml("weights.wrl", 1)) == NULL)
				error("New vrml failed");

			maxw = 0.0;
			for (j = 0; j < smp[i].nnd; j++) {
				if (smp[i].nd[j].w > maxw)	
					maxw = smp[i].nd[j].w;
			}
			for (j = 0; j < smp[i].nnd; j++) {
				wrl->add_col_vertex(wrl, 0, smp[i].sgam->cent, smp[i].nd[j].n == &smp[i] ? red : yellow);
				icmNormalize33(pp, smp[i].nd[j].n->_sv, smp[i].sgam->cent, smp[i].nd[j].w * 50.0/maxw);
				wrl->add_col_vertex(wrl, 0, pp, smp[i].nd[j].n == &smp[i] ? red : yellow);
			}
			wrl->make_lines(wrl, 0, 2);

			wrl->del(wrl);
			printf("Waiting for input after writing 'weights.wrl' for point %d:\n",i);
			getchar();
		}
	}
#endif /* SHOW_NEIGB_WEIGHTS */

	/* Optimise the location of the source to destination mapping. */
	if (verb) printf("Optimizing source to destination mapping...\n");

	VA(("Doing first pass to locate the nearest point\n"));
	/* First pass to locate the weighted nearest point, to use in subsequent passes */
	{
		double s[2] = { 20.0, 20.0 };		/* 2D search area */
		double iv[3];						/* Initial start value */
		double nv[2];						/* 2D New value */
		double tp[3];						/* Resultint value */
		double ne;							/* New error */
		int notrials = NO_TRIALS;

		for (i = 0; i < nmpts; i++) {		/* Move all the points */
			double bnv[2];					/* Best 2d value */
			double brv;						/* Best return value */
			int trial;
			double mv;

			opts.pass = 0;		/* Itteration pass */
			opts.ix = i;		/* Point to optimise */
			opts.p = &smp[i];

			/* If we're expanding, temporarily swap src and radial dest */
			smp[i].swap = 0;
			if (useexp && smp[i].dr > (smp[i].sr + 1e-9)) {
				gamut *tt;
				double dd;

				smp[i].swap = 1;
				tt = smp[i].dgam; smp[i].dgam = smp[i].sgam; smp[i].sgam = tt;

				smp[i].dr = smp[i].sr;
				smp[i].dv[0] = smp[i].sv[0];
				smp[i].dv[1] = smp[i].sv[1];
				smp[i].dv[2] = smp[i].sv[2];

				smp[i].sr = smp[i].drr;
				smp[i].sv[0] = smp[i].drv[0];
				smp[i].sv[1] = smp[i].drv[1];
				smp[i].sv[2] = smp[i].drv[2];
			}
			
			/* Convert our start value from 3D to 2D for speed. */
			icmMul3By3x4(iv, smp[i].m2d, smp[i].dv);
			nv[0] = iv[0] = iv[1];
			nv[1] = iv[1] = iv[2];

			/* Do several trials from different starting points to avoid */
			/* any local minima, particularly with nearest mapping. */
			brv = 1e38;
			for (trial = 0; trial < notrials; trial++) {
				double rv;			/* Temporary */

				/* Optimise the point */
				if (powell(&rv, 2, nv, s, 0.01, 1000, optfunc1, (void *)(&opts), NULL, NULL) == 0
				    && rv < brv) {
					brv = rv;
//printf("~1 point %d, trial %d, new best %f\n",i,trial,sqrt(rv));
					bnv[0] = nv[0];
					bnv[1] = nv[1];
				}
//else printf("~1 powell failed with rv = %f\n",rv);
				/* Adjust the starting point with a random offset to avoid local minima */
				nv[0] = iv[0] + d_rand(-20.0, 20.0);
				nv[1] = iv[1] + d_rand(-20.0, 20.0);
			}
			if (brv == 1e38) {		/* We failed to get a result */
				VB(("multiple powells failed to get a result\n"));
#ifdef NEVER
				/* Optimise the point with debug on */
				opts.debug = 1;
				icmMul3By3x4(iv, smp[i].m2d, smp[i].dv);
				nv[0] = iv[0] = iv[1];
				nv[1] = iv[1] = iv[2];
				powell(NULL, 2, nv, s, 0.01, 1000, optfunc1, (void *)(&opts), NULL, NULL);
#endif
				free_nearsmth(smp, nmpts);
				*npp = 0;
				if (si_gam != sc_gam)
					sci_gam->del(sci_gam);
				if (di_gam != sci_gam && di_gam != sci_gam)
					di_gam->del(di_gam);
				return NULL;
			}
		
			/* Convert best result 2D -> 3D */
			tp[2] = bnv[1];
			tp[1] = bnv[0];
			tp[0] = 50.0;
			icmMul3By3x4(tp, smp[i].m3d, tp);

			/* Remap it to the destinaton gamut surface */
			smp[i].dgam->radial(smp[i].dgam, tp, tp);
			icmCpy3(smp[i].aodv, tp);

			/* Undo any swap */
			if (smp[i].swap) {
				gamut *tt;
				double dd;

				tt = smp[i].dgam; smp[i].dgam = smp[i].sgam; smp[i].sgam = tt;

				/* We get the point on the real src gamut out when swap */
				smp[i]._sv[0] = smp[i].aodv[0];
				smp[i]._sv[1] = smp[i].aodv[1];
				smp[i]._sv[2] = smp[i].aodv[2];

				/* So we need to compute cusp mapped sv */
				comp_ce(&opts, smp[i].sv, smp[i]._sv, &smp[i].wt);
				smp[i].sr = icmNorm33(smp[i].sv, smp[i].sgam->cent);

				VB(("Exp Src %d = %f %f %f\n",i,smp[i]._sv[0],smp[i]._sv[1],smp[i]._sv[2]));
				smp[i].aodv[0] = smp[i].drv[0];
				smp[i].aodv[1] = smp[i].drv[1];
				smp[i].aodv[2] = smp[i].drv[2];
			}
		}
		if (verb) {
			printf("."); fflush(stdout);
		}
	}

	VA(("Doing second pass to optimize without relative error\n"));
	/* Second pass to locate the optimized overall weighted point nsdv[], */
	/* not counting relative error. */
	{
		double s[2] = { 20.0, 20.0 };		/* 2D search area */
		double iv[3];						/* Initial start value */
		double nv[2];						/* 2D New value */
		double tp[3];						/* Resultint value */
		double ne;							/* New error */
		int notrials = NO_TRIALS;

		for (i = 0; i < nmpts; i++) {		/* Move all the points */
			double bnv[2];					/* Best 2d value */
			double brv;						/* Best return value */
			int trial;
			double mv;

			opts.pass = 0;		/* Itteration pass */
			opts.ix = i;		/* Point to optimise */
			opts.p = &smp[i];

//printf("~1 point %d, sv %f %f %f\n",i,smp[i].sv[0],smp[i].sv[1],smp[i].sv[2]);

			/* Convert our start value from 3D to 2D for speed. */
			icmMul3By3x4(iv, smp[i].m2d, smp[i].aodv);

			nv[0] = iv[0] = iv[1];
			nv[1] = iv[1] = iv[2];
//printf("~1 point %d, iv %f %f %f, 2D %f %f\n",i, smp[i].aodv[0], smp[i].aodv[1], smp[i].aodv[2], iv[0], iv[1]);

			/* Do several trials from different starting points to avoid */
			/* any local minima, particularly with nearest mapping. */
			brv = 1e38;
			for (trial = 0; trial < notrials; trial++) {
				double rv;			/* Temporary */

				/* Optimise the point */
				if (powell(&rv, 2, nv, s, 0.01, 1000, optfunc2, (void *)(&opts), NULL, NULL) == 0
				    && rv < brv) {
					brv = rv;
//printf("~1 point %d, trial %d, new best %f at xy %f %f\n",i,trial,sqrt(rv), nv[0],nv[1]);
					bnv[0] = nv[0];
					bnv[1] = nv[1];
				}
//else printf("~1 powell failed with rv = %f\n",rv);
				/* Adjust the starting point with a random offset to avoid local minima */
				nv[0] = iv[0] + d_rand(-20.0, 20.0);
				nv[1] = iv[1] + d_rand(-20.0, 20.0);
			}
			if (brv == 1e38) {		/* We failed to get a result */
				VB(("multiple powells failed to get a result\n"));
#ifdef NEVER
				/* Optimise the point with debug on */
				opts.debug = 1;
				icmMul3By3x4(iv, smp[i].m2d, smp[i].dv);
				nv[0] = iv[0] = iv[1];
				nv[1] = iv[1] = iv[2];
				powell(NULL, 2, nv, s, 0.01, 1000, optfunc2, (void *)(&opts), NULL, NULL);
#endif
				free_nearsmth(smp, nmpts);
				*npp = 0;
				if (si_gam != sc_gam)
					sci_gam->del(sci_gam);
				if (di_gam != sci_gam && di_gam != sci_gam)
					di_gam->del(di_gam);
				return NULL;
			}
		
			/* Convert best result 3D -> 2D */
			tp[2] = bnv[1];
			tp[1] = bnv[0];
			tp[0] = 50.0;
			icmMul3By3x4(tp, smp[i].m3d, tp);

			/* Remap it to the destinaton gamut surface */
			smp[i].dgam->radial(smp[i].dgam, tp, tp);

			icmCpy3(smp[i].nsdv, tp);
			icmCpy3(smp[i].anv, tp);
			icmCpy3(smp[i].dv, tp);
			smp[i].dr = icmNorm33(smp[i].dv, smp[i].dgam->cent);
//printf("~1 %d: dv %f %f %f\n", i, smp[i].dv[0], smp[i].dv[1], smp[i].dv[2]);

		}
		if (verb) {
			printf("."); fflush(stdout);
		}
	}

	/* To speed convergence of the final relative error weighted step, */
	/* establish an initial 2D smoothed destination locations in dv[] from anv[] */
	for (k = 0; k < NO_SMTH_ITERS; k++) {

		VA(("Doing direction pre-filtering\n"));

		/* Create direction filtered dv[] from anv[] */
		for (i = 0; i < nmpts; i++) {
			double mm[3][4];	/* Direction alignment rotation */
			double imm[3][4];	/* and inverse */
			double sv[3], dv[3];
			double ta[3] = { 1.0, 0.0, 0.0 };
			double tc[3] = { 0.0, 0.0, 0.0 };
			double sand, dand;		/* Source and dest average neigbour distances */
			double ndv[3];			/* New ndv[] */
		
//printf("~1 Direction filtering %d, sv %f %f %f, anv %f %f %f\n", i,smp[i].sv[0], smp[i].sv[1], smp[i].sv[2], smp[i].anv[0], smp[i].anv[1], smp[i].anv[2]);

			/* Compute a rotation and offset to set the coordinate system */
			/* so that 1,0,0 is along the mapping direction. */
			icmCpy3(dv, smp[i].anv);
			icmNormalize33(dv, dv, smp[i].sv, 1.0);
			icmVecRotMat(mm, dv, smp[i].sv, ta, tc);
			icmVecRotMat(imm, ta, tc, dv, smp[i].sv);

			/* Compute the average distance of neighbourhood sv[] and dv[]'s */
			/* from the axis */
			sand = dand = 0.0;
			for (j = 0; j < smp[i].nnd; j++) {
				nearsmth *np = smp[i].nd[j].n;		/* Pointer to neighbor */
				double nw = smp[i].nd[j].w;			/* Weight */
		
				icmMul3By3x4(sv, mm, np->sv);
				icmMul3By3x4(dv, mm, np->anv);

				sand += nw * sqrt(sv[1] * sv[1] + sv[2] * sv[2]);
				dand += nw * sqrt(dv[1] * dv[1] + dv[2] * dv[2]);
			}

			/* Compute filtered value from axis */
			ndv[0] = ndv[1] = ndv[2] = 0.0;
			for (j = 0; j < smp[i].nnd; j++) {
				nearsmth *np = smp[i].nd[j].n;		/* Pointer to neighbor */
				double nw = smp[i].nd[j].w;			/* Weight */
		
				icmMul3By3x4(sv, mm, np->sv);
				icmMul3By3x4(dv, mm, np->anv);

				if (sand > 1e-6 && dand > 1e-6)
					icmScale3(sv, sv, dand/sand);		/* Rescale it if we can */

				ndv[1] += nw * (dv[1] - sv[1]);
				ndv[2] += nw * (dv[2] - sv[2]);
			}
			/* Convert back to point at the destination depth */
			icmMul3By3x4(dv, mm, smp[i].anv);
			ndv[0] = dv[0];
			icmMul3By3x4(ndv, imm, ndv);

			/* Check if the point is out of gamut. */
			if (smp[i].dgam->nradial(smp[i].dgam, NULL, ndv) > 1.0) {
				smp[i].dgam->nearest(smp[i].dgam, ndv, ndv);
			}
#ifdef NEVER	/* Don't need this because for exp==0 src gamut is within dst ?? */
			/* Check that the destination point is not out of source gamut */
			if (useexp == 0) {
				double tmp[3];
				/* Un cusp map the point to compare with real src gamut */
				inv_comp_ce(&opts, tmp, ndv, &smp[i].wt);
				if (smp[i].sgam->nradial(smp[i].sgam, NULL, tmp) > 1.0) {
					smp[i].sgam->nearest(smp[i].sgam, ndv, tmp);
					comp_ce(&opts, ndv, ndv, &smp[i].wt);
				}
			}
#endif
			icmCpy3(smp[i].dv, ndv);
		}
	
		/* Copy dv[] to anv[] for another round */
		for (i = 0; i < nmpts; i++) {
			icmCpy3(smp[i].anv, smp[i].dv);
		}
	}		/* Next direction pre-filtering itteration */

#ifdef NEVER
	/* Show just the closest vectors */
	for (i = 0; i < nmpts; i++) {		/* Move all the points */
//		icmCpy3(smp[i].dv, smp[i].drv);			/* Radial */
		icmCpy3(smp[i].dv, smp[i].aodv);		/* Nearest */
//		icmCpy3(smp[i].dv, smp[i].nsdv);		/* No smoothed weighted */
//		icmCpy3(smp[i].dv, smp[i].dv);			/* pre-filter smooothed */
		smp[i].dr = icmNorm33(smp[i].dv, smp[i].dgam->cent);
	}
#else
	/* Finally we do a itterative 3D optimization to balance */
	/* distance from the non-smoothed target nsdv[] with the */
	/* distance from the neighbourhood position target. */

	VA(("Locating smoothed destination\n"));

	for (it = 0; it < MAXITTERS; it++) {
		double mxch = 0.0, avch = 0.0;
		int mxchix;
		double avdot = 0.0, mindot = 1e60, maxdot = -1e60;

		/* Compute the neighbourhood average target anv[] from dv[] */
		for (i = 0; i < nmpts; i++) {
			double sv[3], dv[3];
			double sand, dand;		/* Source and dest average neigbour distances */
			double anv[3];			/* New anv[] */
			double anvl;
		
			/* Compute neigbourhood average sv and dv */
			sv[0] = sv[1] = sv[2] = 0.0;
			dv[0] = dv[1] = dv[2] = 0.0;
			for (j = 0; j < smp[i].nnd; j++) {
				nearsmth *np = smp[i].nd[j].n;		/* Pointer to neighbor */
				double nw = smp[i].nd[j].w;			/* Weight */
				double tmp[3];
		
				icmScale3(tmp, np->sv, nw);
				icmAdd3(sv, sv, tmp);
				icmScale3(tmp, np->dv, nw);
				icmAdd3(dv, dv, tmp);
			}
			/* Compute the average distance of neighbourhood */
			/* sv[] and dv[]'s to do mapping scaling */
			sand = dand = 0.0;
			for (j = 0; j < smp[i].nnd; j++) {
				nearsmth *np = smp[i].nd[j].n;		/* Pointer to neighbor */
				double nw = smp[i].nd[j].w;			/* Weight */
		
				sand += nw * icmNorm33(sv, np->sv);
				dand += nw * icmNorm33(dv, np->dv);
			}

			/* sv->dv scaling factor */
			if (sand < 1e-6 || dand < 1e-6)
				dand = 1.0;
			else
				dand = dand/sand;
//			dand = pow(dand, 0.8);				/* Don't do full scaling */

//#define FPOW 1.0

			/* Compute filtered value */
			anv[0] = anv[1] = anv[2] = 0.0;
			anvl = -1e68;
			for (j = 0; j < smp[i].nnd; j++) {
				nearsmth *np = smp[i].nd[j].n;		/* Pointer to neighbor */
				double nw = smp[i].nd[j].w;			/* Weight */
				double tmp[3], tmpl;
		
				icmSub3(tmp, smp[i].sv,  np->sv);	/* Vector from neighbour to this */
				icmScale3(tmp, tmp, dand);			/* Rescale vector for dest size */
				icmAdd3(tmp, tmp, np->dv);			/* this dest from neighbour */
				
				icmSub3(tmp, tmp, smp[i].sv);		/* Vector from source */
				tmpl = icmNorm3(tmp);

//				tmpl = pow(tmpl, FPOW);

				icmScale3(tmp, tmp, nw);			/* weight for filter */
				icmAdd3(anv, anv, tmp);				/* sum filtered value */

				tmpl *= smp[i].nd[j].rw;
				if (tmpl > anvl)
					anvl = tmpl;
			}
//			anvl = pow(anvl, 1.0/FPOW);

//			icmNormalize3(anv, anv, anvl);
			icmAdd3(anv, anv, smp[i].sv);			/* Filtered destination */

//if (i == 1145) printf("~1 %d: anv found to be %f %f %f\n",i,anv[0],anv[1],anv[2]);
			/* Check if the point is out of gamut. */
			if (smp[i].dgam->nradial(smp[i].dgam, NULL, anv) > 1.0) {
				smp[i].dgam->nearest(smp[i].dgam, anv, anv);
			}
#ifdef NEVER	/* Don't need this because for exp==0 src gamut is within dst ?? */
			/* Check that the destination point is not out of source gamut */
			if (useexp == 0) {
				double tmp[3];
				/* Un cusp map the point to compare with real src gamut */
				inv_comp_ce(&opts, tmp, anv, &smp[i].wt);
				if (smp[i].sgam->nradial(smp[i].sgam, NULL, tmp) > 1.0) {
					smp[i].sgam->nearest(smp[i].sgam, anv, tmp);
					comp_ce(&opts, anv, anv, &smp[i].wt);
				}
			}
#endif
//if (i == 1145) printf("~1 %d: set clipped anv to %f %f %f\n",i,anv[0],anv[1],anv[2]);
			icmCpy3(smp[i].anv, anv);

			/* Make sure anv[] is on the destination gamut at the black and white points. */
			icmSub3(anv, smp[i]._sv, smp[i].sgam->cent);
			if (icmNormalize3(anv, anv, 1.0) == 0) {
				double bf;			/* Blend factor */

				/* Compute how far _sv[] is from poles with sin() weighting */
				bf = pow(anv[1] * anv[1] + anv[2] * anv[2], 0.5);
				icmBlend3(smp[i].anv, smp[i].nsdv, smp[i].anv, bf);
//				icmBlend3(smp[i].anv, smp[i].aodv, smp[i].anv, bf);
			}
		}

		VA(("Doing smoothing itteration %d\n",it+1));

#ifdef NEVER
	/* Use just the smoothed target */
	for (i = 0; i < nmpts; i++) {			/* Move all the points */
		icmCpy3(smp[i].dv, smp[i].anv);		/* No smoothed weighted */
		smp[i].dr = icmNorm33(smp[i].dv, smp[i].dgam->cent);
	}
#else
		/* Do the main pass, which is to balance the delta to the */
		/* target neighbourhood location with the delta to the target */
		/* closest location. Creates dv[] from nsdv[] and anv[] */
		{
			double ta[3] = { 1.0, 0.0, 0.0 };
			double tc[3] = { 0.0, 0.0, 0.0 };
			double s[3] = { 10.0, 10.0, 10.0 };	/* 2D search area */
			double nv[3];						/* 3D New value */
			double ne;							/* New error */
			int notrials = NO_TRIALS;

			for (i = 0; i < nmpts; i++) {		/* Move all the points */
				double bnv[3];					/* Best 3d value */
				double brv;						/* Best return value */
				double del[3];
				int trial;
				double ch;

				opts.pass = 0;		/* Itteration pass */
				opts.ix = i;		/* Point to optimise */
				opts.p = &smp[i];

				/* Compute a rotation and offset to set the coordinate system */
				/* so that 1,0,0 is along the current mapping direction. */
				/* (This is so that we can apply the relative weightings accurately) */
				icmCpy3(nv, smp[i].dv);
				icmNormalize33(nv, nv, smp[i].sv, 1.0);
				icmVecRotMat(opts.mm, nv, smp[i].sv, ta, tc);

				/* Add another rotation to orient it so that [1] corresponds */
				/* with the L direction, and [2] corresponds with the */
				/* hue direction. */
				opts.m2[0][0] = opts.m2[1][1] = 1.0;		/* Default matrix */
				opts.m2[0][1] = opts.m2[1][0] = 0.0;
				nv[0] = smp[i].sv[0] + 1.0;	/* Point offset in L direction */
				nv[1] = smp[i].sv[1];
				nv[2] = smp[i].sv[2];
				icmMul3By3x4(nv, opts.mm, nv);	/* Current transformation of it */
				ch = nv[1] * nv[1] + nv[2] * nv[2];		/* Magnitude of L offset in [1][2] */
				if (ch > 1e-6) {		/* There is a sense of L direction */

					/* Create the 2x2 rotation matrix to align L with [1] */
					ch = sqrt(ch);
					nv[1] /= ch;
					nv[2] /= ch;
	
					opts.m2[0][0] = opts.m2[1][1] = nv[1];
					opts.m2[0][1] = nv[2];
					opts.m2[1][0] = -nv[2];
				}

				/* anv[] transformed into modified coordinates */
				icmMul3By3x4(opts.manv, opts.mm, smp[i].anv);
				icmMulBy2x2(&opts.manv[1], opts.m2, &opts.manv[1]);

//if (i == 1145) printf("\n~1 smoothing %d, initial dv %f %f %f\n", i, smp[i].dv[0], smp[i].dv[1], smp[i].dv[2]);
//if (i == 1145) printf("~1 nsdv = %f %f %f, anv = %f %f %f\n", smp[i].nsdv[0], smp[i].nsdv[1], smp[i].nsdv[2], smp[i].anv[0], smp[i].anv[1], smp[i].anv[2]);

				/* Starting value */
//				nv[0] = smp[i].dv[0];
//				nv[1] = smp[i].dv[1];
//				nv[2] = smp[i].dv[2];
				icmBlend3(nv, smp[i].anv, smp[i].nsdv, 0.5);

				/* Do several trials from different starting points to avoid */
				/* any local minima. */
				brv = 1e38;
				for (trial = 0; trial < notrials; trial++) {
					double tmp[3];
					double rv;			/* Temporary */

//if (i == 1145) printf("~1 %d: starting trial %d at %f %f %f\n", i,trial,nv[0],nv[1],nv[2]);
					/* Optimise the point */
					if (powell(&rv, 3, nv, s, 0.01, 1000, optfunc3, (void *)(&opts), NULL, NULL) == 0
					    && rv < brv) {

						/* Make sure the point is not outside the destination gamut */
						if (smp[i].dgam->nradial(smp[i].dgam, NULL, nv) > 1.0) {
							smp[i].dgam->nearest(smp[i].dgam, nv, nv);
						}
#ifdef NEVER	/* Don't need this because for exp==0 src gamut is within dst ?? */
						/* Make sure the pont point is not out of source gamut */
						if (useexp == 0) {
							/* Un cusp map the point to compare with real src gamut */
							inv_comp_ce(&opts, tmp, nv, &smp[i].wt);
							if (smp[i].sgam->nradial(smp[i].sgam, NULL, tmp) > 1.0) {
								smp[i].sgam->nearest(smp[i].sgam, nv, tmp);
								comp_ce(&opts, nv, nv, &smp[i].wt);
							}
						}
#endif
						brv = rv;
//if (i == 1145) printf("~1 point %d, trial %d, val %f %f %f, new best DE %f\n",i,trial,nv[0],nv[1],nv[2],sqrt(rv));
						bnv[0] = nv[0];
						bnv[1] = nv[1];
						bnv[2] = nv[2];
					}
//else if (i == 1145) printf("~1 %d: trial %d, powell failed or didn't improve rv = %f\n",i, trial, rv);
					/* Adjust the starting point with a random offset to avoid local minima */
					nv[0] = smp[i].dv[0] + d_rand(-10.0, 10.0);
					nv[1] = smp[i].dv[1] + d_rand(-10.0, 10.0);
					nv[2] = smp[i].dv[2] + d_rand(-10.0, 10.0);

					/* Make sure the start point is not out of dst gamut. */
					if (smp[i].dgam->nradial(smp[i].dgam, tmp, nv) > 1.0) {
						icmCpy3(nv, tmp);
					}
#ifdef NEVER	/* Don't need this because for exp==0 src gamut is within dst ?? */
					/* Make sure the start point is not out of source gamut */
					if (useexp == 0) {
						double ttt[3];
						/* Un cusp map the point to compare with real src gamut */
						inv_comp_ce(&opts, ttt, nv, &smp[i].wt);
						if (smp[i].sgam->nradial(smp[i].sgam, tmp, ttt) > 1.0) {
							comp_ce(&opts, nv, tmp, &smp[i].wt);
						}
					}
#endif
				}
				if (brv == 1e38) {		/* We failed to get a result */
					VB(("multiple powells failed to get a result\n"));
					warning("multiple powells failed to get a result");
#ifdef NEVER
					opts.debug = 1;
					nv[0] = smp[i].dv[0];
					nv[1] = smp[i].dv[1];
					nv[2] = smp[i].dv[2];
					powell(&rv, 3, nv, s, 0.01, 1000, optfunc3, (void *)(&opts), NULL, NULL);
#endif /* NEVER */
					free_nearsmth(smp, nmpts);
					*npp = 0;
					if (si_gam != sc_gam)
						sci_gam->del(sci_gam);
					if (di_gam != sci_gam && di_gam != sci_gam)
						di_gam->del(di_gam);
					return NULL;
				}
			
				/* Compute overshoot value */
				icmBlend3(bnv, smp[i].dv, bnv, OSHOOT);

				/* See how much it changed */
				icmSub3(del, bnv, smp[i].dv);
				ch = icmNorm3(del);
				avch += ch;
				if (ch > mxch) {
					mxch = ch;
					mxchix = i;
				}

				/* See if the change is in a consistent direction */
				if (it > 0) {
					double dot = icmDot3(smp[i].pdel, del);
					avdot += dot;
					if (dot > maxdot)
						maxdot = dot;
					if (dot < mindot)
						mindot = dot;
				}
				icmNormalize3(smp[i].pdel, del, 1.0);

				/* Save the best point as the destination, with overshoot */
				icmCpy3(smp[i].dv, bnv);
				smp[i].dr = icmNorm33(smp[i].dv, smp[i].dgam->cent);
//if (i == 1145) printf("~1 point %d, saving %f %f %f, DE %f, change %f\n",i,bnv[0],bnv[1],bnv[2],sqrt(brv),ch);

			}
			avch /= (double)nmpts;
			avdot /= (double)nmpts;
			if (verb) {
				printf("."); fflush(stdout);
			}
#ifdef VERB
			if (it > 0)
				printf("avch = %f, mxch = %f @ %d, avgdot %f, mindot %f, mxadot %f\n",
				                                   avch,mxch,mxchix,avdot,mindot,maxdot);
			else
				printf("avch = %f, mxch = %f @ %d\n",avch,mxch,mxchix);
#endif /* VERB */
			if (mxch <= ITTER_STOP)
				break;
		}
#endif	// NEVER
	}	/* Next itter */

#endif /* NEVER (show debug values) */

	VA(("Itteration done, doing final houskeeping\n"));

	if (verb)
		printf("\n");

	/* Create a plot indicating how the source mapping has been guided by the */
	/* various weighting forces */
#ifdef SAVE_VRMLS
#ifdef PLOT_MAPPING_INFLUENCE
	{
		gamut *gam;
		int src = 0;			/* 1 = src, 0 = dst gamuts */
		vrml *wrl = NULL;
		co *fpnts = NULL;		/* Mapping points to create diagnostic color mapping */
		rspl *swdiag = NULL;
		int gres[3];
		double avgdev[3];
		double cols[4][3] = { { 1.0, 0.0, 0.0 },		/* Absolute = red */
		                      { 1.0, 1.0, 0.0 },		/* Relative = yellow */
		                      { 0.0, 0.0, 1.0 },		/* Radial   = blue */
		                      { 0.0, 1.0, 0.0 } };		/* Depth    = green */
		double grey[3] = { 0.5, 0.5, 0.5 };		/* Grey */
		double max, min;
		int ix;

		if (src)
			gam = sci_gam;
		else
			gam = di_gam;

		/* Setup the scattered data points */
		if ((fpnts = (co *)malloc((nmpts) * sizeof(co))) == NULL) { 
			fprintf(stderr,"gamut map: Malloc of diagnostic mapping setup points failed\n");
			if (si_gam != sc_gam)
				sci_gam->del(sci_gam);
			if (di_gam != sci_gam && di_gam != sci_gam)
				di_gam->del(di_gam);
			free_nearsmth(smp, nmpts);
			*npp = 0;
			return NULL;
		}

		/* Compute error values and diagnostic color */
		/* for each guide vector */
		for (i = 0; i < nmpts; i++) {
			double dv[4], gv;
			double rgb[3];

			/* Source value location */
			if (src) {
				for (j = 0; j < 3; j++)
					fpnts[i].p[j] = smp[i]._sv[j];		/* Non rotated and elevated */
			} else {		/* Dest value location */
				for (j = 0; j < 3; j++)
					fpnts[i].p[j] = smp[i].dv[j];
			}

			/* Diagnostic color */
			max = -1e60; min = 1e60;
			for (k = 0; k < 4; k++) {			/* Find max and min error value */
				dv[k] = smp[i].dbgv[k];
				if (dv[k] > max)
					max = dv[k];
				if (dv[k] < min)
					min = dv[k];
			}
			for (k = 0; k < 4; k++)				/* Scale to max */
				dv[k] /= max;
			max /= max;
			min /= max;
			max -= min;							/* reduce min to zero */
			for (k = 0; k < 4; k++)
				dv[k] /= max;
			for (gv = 1.0, k = 0; k < 4; k++)	/* Blend remainder with grey */
				gv -= dv[k];
			
			for (j = 0; j < 3; j++)				/* Compute interpolated color */
				fpnts[i].v[j] = 0.0;
			for (k = 0; k < 4; k++) {
				for (j = 0; j < 3; j++)
					fpnts[i].v[j] += dv[k] * cols[k][j];
			}
			for (j = 0; j < 3; j++)
				fpnts[i].v[j] += gv * grey[j];
		}

		/* Create the diagnostic color rspl */
		for (j = 0; j < 3; j++) {		/* Set resolution for all axes */
			gres[j] = mapres;
			avgdev[j] = 0.001;
		}
		swdiag = new_rspl(RSPL_NOFLAGS, 3, 3);	/* Allocate 3D -> 3D */
		swdiag->fit_rspl(swdiag, RSPL_NOFLAGS, fpnts, nmpts, NULL, NULL, gres, NULL, NULL, 1.0, avgdev, NULL);

		/* Now create a plot of the sci_gam with the vertexes colored acording to the */
		/* diagnostic map. */
		if ((wrl = new_vrml("sci_gam_wt.wrl", 1)) == NULL) {
			fprintf(stderr,"gamut map: new_vrml failed\n");
			if (fpnts != NULL)
				free(fpnts);
			if (swdiag != NULL)
				swdiag->del(swdiag);
			if (si_gam != sc_gam)
				sci_gam->del(sci_gam);
			if (di_gam != sci_gam && di_gam != sci_gam)
				di_gam->del(di_gam);
			free_nearsmth(smp, nmpts);
			*npp = 0;
			return NULL;
		}

		/* Plot the gamut triangle vertexes */
		for (ix = 0; ix >= 0;) {
			co pp;
			double col[3];

			ix = gam->getvert(gam, NULL, pp.p, ix);
			swdiag->interp(swdiag, &pp);
			wrl->add_col_vertex(wrl, 0, pp.p, pp.v);
		}
		gam->startnexttri(gam);
		for (;;) {
			int vix[3];
			if (gam->getnexttri(gam, vix))
				break;
			wrl->add_triangle(wrl, 0, vix);
		}
		wrl->make_triangles_vc(wrl, 0, 0.0);

		printf("Writing sci_gam_wt.wrl file\n");
		wrl->del(wrl);		/* Write file */
		free(fpnts);
		swdiag->del(swdiag);
	}
#endif /* PLOT_MAPPING_INFLUENCE */
#endif /* SAVE_VRMLS */

	/* Optionally apply some reverse smoothing to the surface points. */
	/* Using a reverse filter helps remove and tendency for many */
	/* to one mappings. */
	/* The reverse mapping rspl is created from the */
	/* guide vectors (dst->src), and then it is used to filter the guide vectors */
	/* to smooth the mapping out, and spread the influence of */
	/* the vectors more evenly across the destination space. */
	if (m21fsm > 0.0) {
		cow *fpnts;			/* Mapping points to create gamut mapping */
		rspl *m21f;
		int gres[3];
		double avgdev[3];

		if (verb)
			printf("Doing inverse mapping smoothing \n");

		/* Setup the scattered data points */
		if ((fpnts = (cow *)malloc((nmpts) * sizeof(cow))) == NULL) { 
			fprintf(stderr,"gamut map: Malloc of mapping setup points failed\n");
			if (si_gam != sc_gam)
				sci_gam->del(sci_gam);
			if (di_gam != sci_gam && di_gam != sci_gam)
				di_gam->del(di_gam);
			free_nearsmth(smp, nmpts);
			*npp = 0;
			return NULL;
		}

		for (i = 0; i < nmpts; i++) {

			/* Set the main gamut hull mapping point */
			for (j = 0; j < 3; j++) {
				fpnts[i].p[j] = smp[i].dv[j];		/* Smoothed destination point */
				fpnts[i].v[j] = smp[i].dv[j] - smp[i].sv[j];	/* displacement from source */
			}
			fpnts[i].w = 1.0;
		}

		/* Create the filtering rspl */
		/* We use this to lookup the mapping for points on the source space gamut */
		/* that result from clipping our grid boundary points */
		for (j = 0; j < 3; j++) {		/* Set resolution for all axes */
			gres[j] = mapres;
			avgdev[j] = 0.005;
		}
		m21f = new_rspl(RSPL_NOFLAGS, 3, 3);	/* Allocate 3D -> 3D */
		m21f->fit_rspl_w(m21f, RSPL_NOFLAGS, fpnts, nmpts, NULL, NULL, gres, NULL, NULL, m21fsm, avgdev, NULL);

		/* Now use the smoothed reverse mapping to re-set the guide vectors */
		for (i = 0; i < nmpts; i++) {
			co pp;
			double n_sv[3];		/* Nearest point */

			opts.ix = i;		/* Point in question */
			opts.p = &smp[i];

			/* Lookup the source point for each destination point */
			for (j = 0; j < 3; j++)
				pp.p[j] = smp[i].dv[j];

			m21f->interp(m21f, &pp);

//printf("Filter %f %f %f to %f %f %f was %f %f %f\n",
//smp[i].dv[0], smp[i].dv[1], smp[i].dv[2],
//pp.v[0], pp.v[1], pp.v[2], smp[i].sv[0], smp[i].sv[1], smp[i].sv[2]); 

			for (j = 0; j < 3; j++)
				n_sv[j] = smp[i].dv[j] - pp.v[j];

//printf("~1 point %d:\n",i);
//printf("~1 sv   = %f %f %f\n",smp[i].sv[0],smp[i].sv[1],smp[i].sv[2]);
//printf("~1 n_sv = %f %f %f\n",n_sv[0],n_sv[1],n_sv[2]);
			
			/* Now we need to re-map the source point to the appropriate gamut surface. */
			/* Source point may be cusp mapped, so undo this */
			/* before doing nearest with source gamut */
			inv_comp_ce(&opts, n_sv, n_sv, &smp[i].wt);
			smp[i].sgam->nearest(smp[i].sgam, smp[i]._sv, n_sv);
			smp[i]._sr = icmNorm33(smp[i]._sv, smp[i].sgam->cent);

			/* For completeness, re-compute cusp mapped source */
			comp_ce(&opts, smp[i].sv, smp[i]._sv, &smp[i].wt);
		}

		free(fpnts);
		m21f->del(m21f);
	}

	VB(("Final points:\n"));

	/* Restore the actual source point if elevated or rotated with cusp mapping, */
	/* and then create sub-surface points. */
	for (i = 0; i < nmpts; i++) {

		VB(("Src %d = %f %f %f\n",i,smp[i].sv[0],smp[i].sv[1],smp[i].sv[2]));
		VB(("Dst %d = %f %f %f\n",i,smp[i].dv[0],smp[i].dv[1],smp[i].dv[2]));
//		VB(("udst %d = %f, ogam = %f\n",i,smp[i].udst,smp[i].ogam));

		/* Save the cusp mapped source value */
		icmCpy3(smp[i].csv, smp[i].sv);

		/* Finally un cusp map the source point */
		inv_comp_ce(&opts, smp[i].sv, smp[i].sv, &smp[i].wt);
		smp[i].sr = icmNorm33(smp[i].sv, smp[i].sgam->cent);

		/* Create a sub-surface mapping point too. */
		/* Note that not every mapping point has a sub-surface point, */
		/* and that the gflag and vflag will be nz if it does. */
		/* We're assuming here that the dv is close to being on the */
		/* destination gamut, so that the vector_isect param will be */
		/* close to 1.0 at the intended destination gamut. */
		{
			double mv[3], ml, nv[3];		/* Mapping vector & length, noralized mv */ 
			double minv[3], maxv[3];
			double mint, maxt;
			gtri *mintri, *maxtri;

			smp[i].vflag = smp[i].gflag = 0;	/* Default unknown */
			smp[i].w2 = 0.0;
			icmSub3(mv, smp[i].dv, smp[i].sv);	/* Mapping vector */
			ml = icmNorm3(mv);		/* It's length */
			if (ml > 0.1) {			/* If mapping is non trivial */

#define PFCOND i == 802

//if (PFCOND) printf("~1 mapping %d = %f %f %f -> %f %f %f\n", i, smp[i].sv[0],smp[i].sv[1],smp[i].sv[2],smp[i].dv[0],smp[i].dv[1],smp[i].dv[2]);
//if (PFCOND) printf("~1 vector %f %f %f, len %f\n",  mv[0], mv[1], mv[2],ml);
				/* Compute actual depth of ray into destination gamut */
				if (di_gam->vector_isect(di_gam, smp[i].sv, smp[i].dv,
				                    minv, maxv, &mint, &maxt, &mintri, &maxtri) != 0) {
					double wp[3], bp[3];		/* Gamut white and black points */
					double p1, napoint[3] = { 50.0, 0.0, 0.0 };		/* Neutral axis point */
					double natarg[3];		/* Neutral axis sub target */
					double adepth1, adepth2 = 1000.0;	/* Directional depth, radial depth */
					double adepth;			/* Minimum available depth */
					double mv2[3], sml;		/* Sub-surface mapping vector & norm. length */

					/* Locate the point on the neutral axis that is closest to */
					/* the guide ray. We use this as a destination direction */
					/* if the sub surface ray gets very long, and to compute */
					/* a sanity check on the available depth. */
					if (d_gam->getwb(d_gam, NULL, NULL, NULL, wp, dst_kbp ? NULL : bp, dst_kbp ? bp : NULL) == 0) {
						if (icmLineLineClosest(napoint, NULL, &p1, NULL, bp, wp,
						                       smp[i].sv,smp[i].dv) == 0) {
							/* Clip it */
							if (p1 < 0.0)
								icmCpy3(napoint, bp);
							else if (p1 > 1.0)
								icmCpy3(napoint, wp);

//if (PFCOND) printf("~1 neutral axis point = %f %f %f\n", napoint[0], napoint[1], napoint[2]);
							/* Compute a normalized available depth from distance */
							/* to closest to neautral axis point */
							if (maxt > 1.0)		/* Compression */
							if ((mint > 1e-8 && maxt > -1e-8)		/* G. & V. Compression */
							 || ((mint < -1e-8 && maxt > -1e-8)		/* G. Exp & V. comp. */
							  && (fabs(mint) < (fabs(maxt) - 1e-8))))
								adepth2 = icmNorm33(napoint, smp[i].dv);
							else				/* Expansion */
								adepth2 = icmNorm33(napoint, smp[i].sv);
						}
#ifdef VERB
						  else {
							printf("icmLineLineClosest failed\n");
						}
#endif
					}
#ifdef VERB
					  else {
						printf("d_gam->getwb failed\n");
					}
#endif

//printf("\n~1 i %d: %f %f %f -> %f %f %f\n   isect at t %f and %f\n", i, smp[i].sv[0], smp[i].sv[1], smp[i].sv[2], smp[i].dv[0], smp[i].dv[1], smp[i].dv[2], mint, maxt);

					/* Only create sub-surface mapping vectors if it makes sense. */
					/* If mapping vector is pointing away from destination gamut, */
					/* (which shouldn't happen), ignore it. If the directional depth */
					/* is very thin compared to the radial depth, indicating that we're */
					/* near a "lip", ignore it. */
					if (mint >= -1e-8 && maxt > 1e-8) {

						if (fabs(mint - 1.0) < fabs(maxt) - 1.0
						 && smp[i].dgam->radial(smp[i].dgam, NULL, smp[i].dv)
						  < smp[i].sgam->radial(smp[i].sgam, NULL, smp[i].dv)) {

//if (PFCOND) printf("~1 point is gamut comp & vect comp.\n");
//if (PFCOND) printf("~1 point is gamut comp & vect comp. mint %f maxt %f\n",mint,maxt);
							adepth1 = ml * 0.5 * (maxt + mint - 2.0);
#ifdef RADIAL_SUBVEC
							adepth = adepth2;		/* Always radial depth */
#else
							adepth = adepth1 < adepth2 ? adepth1 : adepth2;		/* Smaller of the two */
#endif
							if (adepth1 < (0.5 * adepth2))
								continue;
//if (PFCOND) printf("~1 dir adepth %f, radial adapeth %f\n",adepth1,adepth2);
							adepth *= 0.9;				/* Can't use 100% */
							smp[i].gflag = 1;		/* Gamut compression and */
							smp[i].vflag = 1;		/* vector compression */

							/* Compute available depth and knee factor adjusted sub-vector */
							icmCpy3(smp[i].sv2, smp[i].dv);		/* Sub source is guide dest */
							ml *= (1.0 - gamcknf);				/* Scale by knee */
							adepth *= (1.0 - gamcknf);
							sml = ml < adepth ? ml : adepth;	/* Smaller of two */
//if (PFCOND) printf("~1 adjusted subvec len %f\n",sml);
							icmNormalize3(mv2, mv, sml);		/* Full sub-surf disp. == no knee */
							icmAdd3(mv2, smp[i].sv2, mv2);		/* Knee adjusted destination */
	
//if (PFCOND) printf("~1 before blend sv2 %f %f %f, dv2 %f %f %f\n", smp[i].sv2[0], smp[i].sv2[1], smp[i].sv2[2], mv2[0], mv2[1], mv2[2]);
							/* Blend towards n.axis as length of sub vector approaches */
							/* distance to neutral axis. */
							icmSub3(natarg, napoint, smp[i].sv2);
							icmNormalize3(natarg, natarg, sml);		/* Sub vector towards n.axis */
							icmAdd3(natarg, natarg, smp[i].sv2);	/* n.axis target */
#ifdef RADIAL_SUBVEC
							icmCpy3(mv2, natarg);			/* Radial direction vector */
#else
							icmBlend3(mv2, mv2, natarg, sml/adepth2);
#endif /* RADIAL_SUBVEC */
//if (PFCOND) printf("~1 after blend sv2 %f %f %f, dv2 %f %f %f\n", smp[i].sv2[0], smp[i].sv2[1], smp[i].sv2[2], mv2[0], mv2[1], mv2[2]);
							
							icmCpy3(smp[i].dv2, mv2);				/* Destination */
							icmCpy3(smp[i].temp, smp[i].dv2);		/* Save a copy to temp */
							smp[i].w2 = 0.8;
						} else {
//if (PFCOND) printf("~1 point is gamut exp & vect exp. mint %f maxt %f\n",mint,maxt);
							smp[i].gflag = 2;		/* Gamut expansion and */
							smp[i].vflag = 0;		/* vector expansion, */
													/* but crossing over, so no sub vect. */
//if (PFCOND) printf("~1 point is crossover point\n",mint,maxt);
						}
	
					} else if (mint < -1e-8 && maxt > 1e-8) {

//if (PFCOND) printf("~1 point is gamut exp & vect exp. mint %f maxt %f\n",mint,maxt);
						/* This expand/expand case has reversed src/dst sense to above */
						adepth1 = ml * 0.5 * -mint;
#ifdef RADIAL_SUBVEC
						adepth = adepth2;		/* Always radial depth */
#else
						adepth = adepth1 < adepth2 ? adepth1 : adepth2;
#endif
//if (PFCOND) printf("~1 dir adepth %f, radial adapeth %f\n",adepth1,adepth2);
						adepth *= 0.9;				/* Can't use 100% */

						if (adepth1 < (0.6 * adepth2))
							continue;

						smp[i].gflag = 2;		/* Gamut expansion */
						smp[i].vflag = 2;		/* vector is expanding */

						icmCpy3(smp[i].dv2, smp[i].sv);	/* Sub dest is guide src */
						ml *= (1.0 - gamxknf);			/* Scale by knee */
						adepth *= (1.0 - gamxknf);
						sml = ml < adepth ? ml : adepth;/* Smaller of two */
						icmNormalize3(mv2, mv, sml);	/* Full sub-surf disp. == no knee */
						icmSub3(mv2, smp[i].dv2, mv2);	/* Knee adjusted source */
	
						/* Blend towards n.axis as length of sub vector approaches */
						/* distance to neutral axis. */
						icmSub3(natarg, smp[i].dv2, napoint);
						icmNormalize3(natarg, natarg, sml);	/* Sub vector away n.axis */
						icmSub3(natarg, smp[i].dv2, natarg);/* n.axis oriented source */
#ifdef RADIAL_SUBVEC
						icmCpy3(mv2, natarg);				/* Radial direction vector */
#else
						icmBlend3(mv2, mv2, natarg, sml/adepth2);	/* dir adjusted src */
#endif /* RADIAL_SUBVEC */

						icmCpy3(smp[i].sv2, mv2); 			/* Source */
						icmCpy3(smp[i].temp, smp[i].dv2);	/* Save a copy to temp */
						smp[i].w2 = 0.8;

					} else {
						/* Nonsense vector */
						smp[i].gflag = 0;		/* Gamut compression but */
						smp[i].vflag = 0;		/* vector is expanding */
//if (PFCOND) printf("~1 point is nonsense vector mint %f maxt %f\n",mint,maxt);

						icmCpy3(smp[i].dv, smp[i].aodv);	/* Clip to the destination gamut */
					}
				}
			}
		}

#ifdef NEVER	// Diagnostic
		smp[i].vflag = 0;	/* Disable sub-points */
#endif /* NEVER */

		VB(("Out Src %d = %f %f %f\n",i,smp[i].sv[0],smp[i].sv[1],smp[i].sv[2]));
		VB(("Out Dst %d = %f %f %f\n",i,smp[i].dv[0],smp[i].dv[1],smp[i].dv[2]));
//		VB(("udst %d = %f, ogam = %f\n",i,smp[i].udst,smp[i].ogam));
		if (smp[i].vflag != 0) {
			VB(("Out Src2 %d = %f %f %f\n",i,smp[i].sv2[0],smp[i].sv2[1],smp[i].sv2[2]));
			VB(("Out Dst2 %d = %f %f %f\n",i,smp[i].dv2[0],smp[i].dv2[1],smp[i].dv2[2]));
		}
	}

#ifndef NEVER
	/* Smooth the sub-surface mapping points */
	/* dv2[] is duplicated in temp[], so use temp[] as the values to be filtered */
	for (i = 0; i < nmpts; i++) {
		double sand, dand;		/* Source and dest average neiigbour distances */
		double fdv2[3];			/* Filtered dv2[] */
		double tw;				/* Total weight */
		int rc;

		if (smp[i].vflag == 0)
			continue;

		/* First compute the average distance of neighbourhood, */
		/* sv2[]'s from this point, so that we can scale the destination */
		/* neighbourhood distances appropriately. */
		sand = dand = tw = 0.0;
		for (j = 0; j < smp[i].nnd; j++) {
			nearsmth *np = smp[i].nd[j].n;		/* Pointer to neighbor */

			if (np->vflag) {
				sand += smp[i].nd[j].w * icmNorm33(smp[i].sv2, np->sv2);
				dand += smp[i].nd[j].w * icmNorm33(smp[i].temp, np->temp);
				tw += smp[i].nd[j].w;
			}
		}
		sand /= tw;
		dand /= tw;

		/* Then compute the scaled sum of weighted tanv[] neighbours location */
		if (sand > 1e-6 && dand > 1e-6) {
			tw = fdv2[0] = fdv2[1] = fdv2[2] = 0.0;
			for (j = 0; j < smp[i].nnd; j++) {
				nearsmth *np = smp[i].nd[j].n;		/* Pointer to neighbor */
				double nw = smp[i].nd[j].w;			/* Weight */
				double vec[3];						/* vector */

				if (np->vflag) {
					icmSub3(vec, smp[i].sv2, np->sv2);	/* Source neighbour displacement */
					icmScale3(vec, vec, dand/sand);
					fdv2[0] += nw * (vec[0] + np->temp[0]);
					fdv2[1] += nw * (vec[1] + np->temp[1]);
					fdv2[2] += nw * (vec[2] + np->temp[2]);
					tw += nw;
				}
			}
			icmScale3(fdv2, fdv2, 1.0/tw);

//printf("~1 %d: moved %f %f %f -> %f %f %f de %f\n", i, smp[i].dv2[0], smp[i].dv2[1], smp[i].dv2[2], fdv2[0], fdv2[1], fdv2[2], icmNorm33(smp[i].dv2,fdv2));

			icmCpy3(smp[i].dv2, fdv2);
		}
	}
#endif /* !NEVER */

#ifdef NEVER
	{
		double sudist = 100.0;
		double ludist = -100.0;
		double logam;
		for (i = 0; i < nmpts; i++) {
			if (smp[i].udst < sudist)
				sudist = smp[i].udst;
			if (smp[i].udst > ludist) {
				ludist = smp[i].udst;
				logam = smp[i].ogam;
			}
		}
		printf("smallest udst = %f\n",sudist);
		printf("largest udst = %f\n",ludist);
		printf("ogam at ludst = %f\n",logam);
	}
#endif /* NEVER */

	if (si_gam != sc_gam)
		sci_gam->del(sci_gam);
	if (di_gam != sci_gam && di_gam != sci_gam)
		di_gam->del(di_gam);
	*npp = nmpts;
	return smp;
}

/* Free the list of points that was returned */
void free_nearsmth(nearsmth *smp, int nmpts) {
	int i;

	for (i = 0; i < nmpts; i++) {
		if (smp[i].nd != NULL)
			free(smp[i].nd);
	}
	free(smp);
}



















