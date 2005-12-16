
/* 
 * Argyll Gamut Mapping Library
 *
 * Author:  Graeme W. Gill
 * Date:    1/10/00
 * Version: 2.00
 *
 * Copyright 2000 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 *
 * For a discussion of gamut mapping strategy used,
 * see gammap.txt
 */

/*
 * TTBD:
 *       Check that color "knee" stuff calculations are correct.
 *
 *       Improve error handling.
 */

#undef DEBUG		/* Plot L map */
#define VERBOSE		/* Print out interesting information */
#define XRES 100	/* Res of plot */
#undef PLOT_GAMVEC	/* Save the gamut mapping points as "gammap.wrl" */
#undef PLOT_GAMUTS	/* Save input and output gamuts as src.wrl, img.wrl, dst.wrl, gmsrc.wrl */

#define INDEPL		/* Remap L independenly to everything else (this is good) */
#undef RADIAL_HACK	/* Hack to use radial mapping rather than nearest */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "gamut.h"
#include "rspl.h"
#include "gammap.h"
#include "nearsmth.h"
#ifdef DEBUG
#include "plot.h"
#endif

/* Context for enhancing the saturation of the clut values */
typedef struct {
	gamut *dst;			/* Destination colorspace gamut */
	double wp[3], bp[3];/* Destination colorspace white and black points */
	double satenh;		/* Saturation engancement value */
} adjustsat;

/* Context for making clut relative to white and black points */
typedef struct {
	double sbp[3], dbp[3];	/* Mapping source and dest */
	double mat[3][3];
} adjustwb;

static void adjust_wb_func(void *pp, double *out, double *in);
static void adjust_sat_func(void *pp, double *out, double *in);

#if defined(PLOT_GAMVEC)

FILE *start_vrml(char *name, int doaxes);
void start_line_set(FILE *wrl);
void add_vertex(FILE *wrl, double pp[3]);
void make_lines(FILE *wrl, int ppset);
void end_vrml(FILE *wrl);

#endif

#define USE_GLUMKNF			/* Enable luminence knee function points */
#define USE_GAMKNF			/* Enable gamut boundary knee function points */

#define OVERSHOOT 1.0		/* Gammut mapping vector overshoot (experemental) */

/* The smoothed near weighting control values: */

/* ----------- */
/* PERCEPTUAL: */

/* Weighting of absolute error of destination from source */
#define P_A_WEIGHT 1.0		/* Absolute error overall weight */
#define P_A_LWEIGHT 3.0		/* Absolute luminance error weight */
#define P_A_CWEIGHT1 1.2	/* Absolute chroma error weight at L == 100 (ie. Yellow) */
#define P_A_CWEIGHT0 0.5	/* Absolute chroma error weight at L == 50 */
#define P_A_HWEIGHT 3.1		/* Absolute hue error weight */

/* Weighting of relative error of destination points to each */
/* other, compared to source points to each other. */
#define P_R_WEIGHT 1.0		/* Relative error overall weight */
#define P_R_LWEIGHT 1.0		/* Relative luminance error weight */
#define P_R_CWEIGHT 0.5		/* Relative chroma error weight */
#define P_R_HWEIGHT 1.0		/* Relative hue error weight */

/* Weighting of error between destination point and source */
/* point radially mapped to destination. (Generally not used other than for testing) */
#define P_D_WEIGHT  0.0		/* Radial error overall weight */
#define P_D_LWEIGHT 1.0		/* Radial luminance error weight */
#define P_D_CWEIGHT 0.5		/* Radial chroma error weight */
#define P_D_HWEIGHT 1.0		/* Radial hue error weight */

/* ----------- */
/* SATURATION: */
/* Weighting of absolute error of destination from source */
#define S_A_WEIGHT 1.0		/* Absolute error overall weight */
#define S_A_LWEIGHT 1.0		/* Absolute luminance error weight */
#define S_A_CWEIGHT1 1.5	/* Absolute chroma error weight at L == 100 (ie. Yellow) */
#define S_A_CWEIGHT0 1.0	/* Absolute chroma error weight at L == 50 */
#define S_A_HWEIGHT 0.7		/* Absolute hue error weight */

/* Weighting of relative error of destination points to each */
/* other, compared to source points to each other. */
#define S_R_WEIGHT 0.0		/* Relative error overall weight */
#define S_R_LWEIGHT 1.0		/* Relative luminance error weight */
#define S_R_CWEIGHT 0.5		/* Relative chroma error weight */
#define S_R_HWEIGHT 1.0		/* Relative hue error weight */

/* Weighting of error between destination point and source */
/* point radially mapped to destination. (Generally not used other than for testing) */
#define S_D_WEIGHT  0.0		/* Radial error overall weight */
#define S_D_LWEIGHT 1.0		/* Radial luminance error weight */
#define S_D_CWEIGHT 0.5		/* Radial chroma error weight */
#define S_D_HWEIGHT 1.0		/* Radial hue error weight */

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#ifdef NEVER // ###### TESTING VALUES #######

/* Weighting of absolute error of destination from source */
#define A_WEIGHT 1.0		/* Absolute error overall weight */
#define A_LWEIGHT 4.0		/* Absolute luminance error weight */
#define A_CWEIGHT1 1.5		/* Absolute chroma error weight at L == 100 */
							/* Designed to ensure yellow maps to pure yellow */
#define A_CWEIGHT0 0.3		/* Absolute chroma error weight at L == 50 */
#define A_HWEIGHT 4.0		/* Absolute hue error weight */

/* Weighting of relative error of destination points to each */
/* other, compared to source points to each other. */
#define R_WEIGHT 1.0		/* Relative error overall weight */
#define R_LWEIGHT 1.0		/* Relative luminance error weight */
#define R_CWEIGHT 0.5		/* Relative chroma error weight */
#define R_HWEIGHT 1.0		/* Relative hue error weight */

/* Weighting of error between destination point and source */
/* point radially mapped to destination. (Generally not used other than for testing) */
#define D_WEIGHT  1.0		/* Radial error overall weight - allow it to dominate */
#define D_LWEIGHT 1.0		/* Radial luminance error weight */
#define D_CWEIGHT 0.5		/* Radial chroma error weight */
#define D_HWEIGHT 1.0		/* Radial hue error weight */

#endif	// NEVER
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

/*
 * Notes:
 *       The "knee" shape produced by the rspl (thin plate spline) code
 *       is not what one would expect for expansion. It is not
 *       symetrical with compression, and is less "sharp". This
 *       is due to the rspl "smoothness" criteria being flawed.
 *       Rather than smoothness being measured as curvature, it
 *       is measured as grid difference values, and this means
 *       that the spline gets "stiffer" as it increases in slope.
 *       Possibly rspl could be improved in this respect ???
 */

static void del_gammap(gammap *s);
static void domap(gammap *s, double *out, double *in);
static void map_trans(void *cntx, double out[3], double in[3]);

/* Return a gammap to map from the input space to the output space */
/* Return NULL on error. */
gammap *new_gammap(
	int verb,			/* Verbose flag */
	gamut *sc_gam,		/* Source colorspace gamut */
	gamut *s_gam,		/* Source image gamut (NULL if none) */
	gamut *d_gam,		/* Destination colorspace gamut */
	double greymf,		/* Grey axis hue matching factor, 0.0 - 1.0 */
	double glumwcpf,	/* Grey axis White luminance compression factor, 0.0 - 1.0 */
	double glumwexf,	/* Grey axis White luminance expansion factor,   0.0 - 1.0 */
	double glumbcpf,	/* Grey axis Black luminance compression factor, 0.0 - 1.0 */
	double glumbexf,	/* Grey axis Black luminance expansion factor,   0.0 - 1.0 */
	double glumknf,		/* Grey axis luminance knee factor, 0.0 - 1.0 */
	double gamcpf,		/* Gamut compression factor, 0.0 - 1.0 */
	double gamexf,		/* Gamut expansion factor, 0.0 - 1.0 */
	double gamknf,		/* Gamut knee factor, 0.0 - 1.0 */
	double gampwf,		/* Gamut Perceptual Map weighting factor, 0.0 - 1.0 */
	double gamswf,		/* Gamut Saturation Map weighting factor, 0.0 - 1.0 */
	double satenh,		/* Saturation enhancement value, 0.0 - Inf */
	int    mapres,		/* Gamut map resolution, typically 9 - 33 */
	double *mn,			/* If not NULL, set minimum mapping input range */
	double *mx			/* for rspl grid */
) {
	gammap *s;			/* This */
	gamut *scl_gam;		/* Source colorspace gamut with L mapping applied */
	gamut *sl_gam;		/* Source image gamut with L mapping applied */
	double s_cs_wp[3];	/* Source colorspace white point */
	double s_cs_bp[3];	/* Source colorspace black point */
	double s_ga_wp[3];	/* Source (image) gamut white point */
	double s_ga_bp[3];	/* Source (image) gamut black point */
	double sl_cs_wp[3];	/* Source L mapped colorspace white point */
	double sl_cs_bp[3];	/* Source L mapped colorspace black point */
	double d_cs_wp[3];	/* Destination colorspace white point */
	double d_cs_bp[3];	/* Destination colorspace black point */
	double d_gx_wp[3];	/* Matching factor weighted destination grey axis white point */
	double d_gx_bp[3];	/* Matching factor weighted destination grey axis black point */
	double s_mt_wp[3];	/* Overall source mapping target white point (used for finetune) */
	double s_mt_bp[3];	/* Overall source mapping target black point (used for finetune) */
	double d_mt_wp[3];	/* Overall destination mapping white point (used for finetune) */
	double d_mt_bp[3];	/* Overall destination mapping black point (used for finetune) */
	cow lpnts[10];	/* Mapping points to create grey axis map */
	int ngreyp = 0;	/* Number of grey axis mapping points */
	rspl *grey;		/* Grey axis mapping rspl */
	cow *gpnts;		/* Mapping points to create gamut mapping */
	int ngamp = 0;	/* Number of gamut mapping points */
	rspl *map;		/* Whole gamut mapping */

#ifdef VERBOSE
	printf("gamut_match parameters:\n");
	printf("Grey axis hue matching factor: %f\n",greymf);
	printf("Grey axis luminance white compression factor: %f\n",glumwcpf);
	printf("Grey axis luminance white expansion factor: %f\n",glumwexf);
	printf("Grey axis luminance black compression factor: %f\n",glumbcpf);
	printf("Grey axis luminance black expansion factor: %f\n",glumbexf);
	printf("Grey axis luminance knee factor: %f\n",glumknf);
	printf("Gamut compression factor: %f\n",gamcpf);
	printf("Gamut expansion factor: %f\n",gamexf);
	printf("Gamut knee factor: %f\n",gamknf);
	printf("Gamut mapping Perceptual Mapping weighting factor: %f\n",gampwf);
	printf("Gamut mapping Saturation Mapping weighting factor: %f\n",gamswf);
	printf("Saturation enhancement factor %f\n",satenh);
	printf("Gamut map resolution: %d\n",mapres);
	if (s_gam != NULL)
		printf("Image gamut supplied\n");
#endif /* VERBOSE */

#ifdef RADIAL_HACK
	printf("!!!! RADIAL HACK IS DEFINED !!!!\n");
#endif

	/* Allocate the object */
	if ((s = (gammap *)calloc(1, sizeof(gammap))) == NULL) {
		fprintf(stderr,"gammap: calloc failed on gammap object\n");
		exit (-1);
	}

	/* Setup methods */
	s->del = del_gammap;
	s->domap = domap;

	/* Now create everything */

	if (s_gam == NULL) {
		s_gam = sc_gam;		/* Source space is source gamut */

		/* Grab the white and black points */
		/* Note that the (image) gamut white & black points should always be */
		/* within their colorspace */
		if (s_gam->getwb(s_gam, s_cs_wp, s_cs_bp, s_ga_wp, s_ga_bp)) {
			fprintf(stderr,"gamut map: Unable to read source white and black points\n");
			free(s);
			return NULL;
		}
	} else {
		if (sc_gam->getwb(sc_gam, s_cs_wp, s_cs_bp, NULL, NULL)) {
			fprintf(stderr,"gamut map: Unable to read source colorspace white and black points\n");
			free(s);
			return NULL;
		}
		if (s_gam->getwb(s_gam, NULL, NULL, s_ga_wp, s_ga_bp)) {
			fprintf(stderr,"gamut map: Unable to read source gamut white and black points\n");
			free(s);
			return NULL;
		}

		/* Guard against silliness. Image must be within colorspace */
		if (s_ga_wp[0] > s_cs_wp[0]) {
			int j;
			double t;
#ifdef VERBOSE
			printf("Fixing wayward image white point\n");
#endif
			t = (s_cs_wp[0] - s_ga_bp[0])/(s_ga_wp[0] - s_ga_bp[0]);
			for (j = 0; j < 3; j++)
				s_ga_wp[j] = s_ga_bp[j] + t * (s_ga_wp[j] - s_ga_bp[j]);

		}
		if (s_ga_bp[0] < s_cs_bp[0]) {
			int j;
			double t;
#ifdef VERBOSE
			printf("Fixing wayward image black point\n");
#endif
			t = (s_cs_bp[0] - s_ga_wp[0])/(s_ga_bp[0] - s_ga_wp[0]);
			for (j = 0; j < 3; j++)
				s_ga_bp[j] = s_ga_wp[j] + t * (s_ga_bp[j] - s_ga_wp[j]);
		}
	}

	if (d_gam->getwb(d_gam, d_cs_wp, d_cs_bp, NULL, NULL)) {
		fprintf(stderr,"gamut map: Unable to read destination white and black points\n");
		free(s);
		return NULL;
	}

#ifdef VERBOSE
	printf("Src colorspace white/black are %f %f %f, %f %f %f\n",
	s_cs_wp[0], s_cs_wp[1], s_cs_wp[2], s_cs_bp[0], s_cs_bp[1], s_cs_bp[2]);

	printf("Src gamut white/black are %f %f %f, %f %f %f\n",
	s_ga_wp[0], s_ga_wp[1], s_ga_wp[2], s_ga_bp[0], s_ga_bp[1], s_ga_bp[2]);

	printf("Dst colorspace white/black are %f %f %f, %f %f %f\n",
	d_cs_wp[0], d_cs_wp[1], d_cs_wp[2], d_cs_bp[0], d_cs_bp[1], d_cs_bp[2]);
#endif /* VERBOSE */

	/* ------------------------------------ */
	/* Figure out the destination grey axis */
	{
		int j;

		/* The grey axis mapping target is a blend between the */
		/* source colorspace (NOT gamut) and the destination */
		/* colorspace. */

		if (greymf > 0.99) {	// Don't mess about and loose accuracy
			for (j = 0; j < 3; j++) {
				d_gx_wp[j] = d_cs_wp[j];
				d_gx_bp[j] = d_cs_bp[j];
			}

		} else {				// Compute blended destination
			double sb[3], sw[3];	/* Normalised source grey axis, L = 100, 0 */
			double db[3], dw[3];	/* Normalised destination grey axis, L = 100, 0 */
			double t;

			/* Compute normalised to L 100..0 source and destination */
			/* white and black points, given them a common basis, so */
			/* that they can be blended togetherm to create an output */
			/* output axis direction. */
			t = (100.0 - s_cs_bp[0])/(s_cs_wp[0] - s_cs_bp[0]);
			for (j = 0; j < 3; j++)
				sw[j] = s_cs_bp[j] + t * (s_cs_wp[j] - s_cs_bp[j]);

			t = (0.0 - s_cs_bp[0])/(s_cs_wp[0] - s_cs_bp[0]);
			for (j = 0; j < 3; j++)
				sb[j] = s_cs_bp[j] + t * (s_cs_wp[j] - s_cs_bp[j]);
			
			t = (100.0 - d_cs_bp[0])/(d_cs_wp[0] - d_cs_bp[0]);
			for (j = 0; j < 3; j++)
				dw[j] = d_cs_bp[j] + t * (d_cs_wp[j] - d_cs_bp[j]);

			t = (0.0 - d_cs_bp[0])/(d_cs_wp[0] - d_cs_bp[0]);
			for (j = 0; j < 3; j++)
				db[j] = d_cs_bp[j] + t * (d_cs_wp[j] - d_cs_bp[j]);
			
			/* Make destination axis direction a blend between the */
			/* source colorspace and destination, acording to greymf. */
			for (j = 0; j < 3; j++) {
				d_gx_wp[j] = greymf * dw[j] + (1.0 - greymf) * sw[j];
				d_gx_bp[j] = greymf * db[j] + (1.0 - greymf) * sb[j];
			}

			/* Map normalised direction vector back to the actual destination */
			/* gamut limits. */
			if (d_gam->vector_isect(d_gam, d_gx_bp, d_gx_wp, d_gx_bp, d_gx_wp, NULL, NULL) == 0) {
				fprintf(stderr,"gamut: vector_isect failed!\n");
				exit(-1);
			}
		}
#ifdef VERBOSE
		printf("Destination grey axis target wp = %f %f %f, bp = %f %f %f\n",
				d_gx_wp[0], d_gx_wp[1], d_gx_wp[2], d_gx_bp[0], d_gx_bp[1], d_gx_bp[2]);
#endif
	}

#ifdef NEVER
//~~9	test extremes
s_cs_wp[0] = 100.0;
s_cs_bp[0] = 30.0;
d_gx_wp[0] = 80.0;
d_gx_bp[0] = 10.0;
glumknf	= 1.0;
#endif /* NEVER */

	/* Create the mapping points needed to build the 1D L mapping rspl. */
	/* If we have a gamut (ie. image) range that is smaller than the */
	/* range of the colorspace, then there is no point allowing for */
	/* points outside the gamut range in creating the L mapping. */
	{
		double swL, dwL;			/* Source and destination white point L */
		double sbL, dbL;			/* Source and destination black point L */
		int j;
		double t;

		/* Setup white point mapping */
		if (s_cs_wp[0] <= d_gx_wp[0]) {	/* Needs possible expansion */
			swL = s_cs_wp[0];
			dwL = glumwexf * d_gx_wp[0] + (1.0 - glumwexf) * s_cs_wp[0];

		} else {
			if (s_ga_wp[0] > d_gx_wp[0]) {	/* Gamut or colorspace needs compression */
				
				swL = (1.0 - glumwcpf) * d_gx_wp[0] + glumwcpf * s_ga_wp[0];
				dwL = d_gx_wp[0];

			} else {	/* Neither needed */
				swL = s_ga_wp[0];
				dwL = s_ga_wp[0];
			}
		}

		/* Setup black point mapping */
		if (s_cs_bp[0] >= d_gx_bp[0]) {	/* Needs possible expansion */
			sbL = s_cs_bp[0];
			dbL = glumbexf * d_gx_bp[0] + (1.0 - glumbexf) * s_cs_bp[0];

		} else {
			if (s_ga_bp[0] < d_gx_bp[0]) {	/* Gamut or colorspace needs compression */
				
				sbL = (1.0 - glumbcpf) * d_gx_bp[0] + glumbcpf * s_ga_bp[0];
				dbL = d_gx_bp[0];

			} else {	/* Neither needed */
				sbL = s_ga_bp[0];
				dbL = s_ga_bp[0];
			}
		}

		/* White point end */
		lpnts[ngreyp].p[0] = swL;
		lpnts[ngreyp].v[0] = dwL;
		lpnts[ngreyp++].w  = 1.0;

		/* Black point end */
		lpnts[ngreyp].p[0] = sbL;
		lpnts[ngreyp].v[0] = dbL;
		lpnts[ngreyp++].w  = 1.0;

#ifdef USE_GLUMKNF
		if (glumknf < 0.05)
#endif /* USE_GLUMKNF */
		{			/* make sure curve is firmly anchored */
			lpnts[ngreyp].p[0] = 0.3 * lpnts[ngreyp-1].p[0] + 0.7 * lpnts[ngreyp-2].p[0];
			lpnts[ngreyp].v[0] = 0.3 * lpnts[ngreyp-1].v[0] + 0.7 * lpnts[ngreyp-2].v[0];
			lpnts[ngreyp++].w  = 1.0;	

			lpnts[ngreyp].p[0] = 0.7 * lpnts[ngreyp-2].p[0] + 0.3 * lpnts[ngreyp-3].p[0];
			lpnts[ngreyp].v[0] = 0.7 * lpnts[ngreyp-2].v[0] + 0.3 * lpnts[ngreyp-3].v[0];
			lpnts[ngreyp++].w  = 1.0;	
		}
#ifdef USE_GLUMKNF
		else {		/* There is at least some weight in knee points */
			double wt1, wt2, wt;
			double bt1, bt2, bt;

			/* Center point */
			lpnts[ngreyp].p[0] = 55.0;
			lpnts[ngreyp].v[0] = 55.0;
			lpnts[ngreyp++].w  = glumknf;	
	
			wt1 = (1.0 * swL + 55.0)/2.0;
			wt2 = (1.0 * dwL + 55.0)/2.0;
			wt = wt1 < wt2 ? wt1 : wt2;

			bt1 = (1.0 * sbL + 55.0)/2.0;
			bt2 = (1.0 * dbL + 55.0)/2.0;
			bt = bt1 > bt2 ? bt1 : bt2;

			/* Middle 1:1 points to cause "knee" curve */
			lpnts[ngreyp].p[0] = wt;
			lpnts[ngreyp].v[0] = wt;
			lpnts[ngreyp++].w  = glumknf * glumknf;	
		
			lpnts[ngreyp].p[0] = bt;
			lpnts[ngreyp].v[0] = bt;
			lpnts[ngreyp++].w  = glumknf * glumknf;	
		}
#endif /* USE_GLUMKNF */

		/* Remember our source and destinatio mapping targets */
		/* so that we can use them for fine tuning later. */
		t = (swL - s_cs_bp[0])/(s_cs_wp[0] - s_cs_bp[0]);
		for (j = 0; j < 3; j++)
			s_mt_wp[j] = s_cs_bp[j] + t * (s_cs_wp[j] - s_cs_bp[j]);

		t = (sbL - s_cs_bp[0])/(s_cs_wp[0] - s_cs_bp[0]);
		for (j = 0; j < 3; j++)
			s_mt_bp[j] = s_cs_bp[j] + t * (s_cs_wp[j] - s_cs_bp[j]);

		t = (dwL - d_gx_bp[0])/(d_gx_wp[0] - d_gx_bp[0]);
		for (j = 0; j < 3; j++)
			d_mt_wp[j] = d_gx_bp[j] + t * (d_gx_wp[j] - d_gx_bp[j]);

		t = (dbL - d_gx_bp[0])/(d_gx_wp[0] - d_gx_bp[0]);
		for (j = 0; j < 3; j++)
			d_mt_bp[j] = d_gx_bp[j] + t * (d_gx_wp[j] - d_gx_bp[j]);
	}

	{
		datai il, ih;
		datao ol, oh;
		int gres = 256;

		/* Create a (possibly temporary) 1D rspl, that is used to */
		/* form the grey axis L mapping. */
		grey = new_rspl(1, 1);	/* Allocate 1D -> 1D */
	
		il[0] = -1.0;		/* Set possible input range */
		ih[0] = 101.0;
		ol[0] = 0.0;		/* Set normalisation output range */
		oh[0] = 100.0;

#ifdef NEVER		/* Dump out the L mapping points */
		{
			int i;
			for (i = 0; i < ngreyp; i++)
				printf("%d %f -> %f (w %f)\n",i,lpnts[i].p[0],lpnts[i].v[0],lpnts[i].w);
		}
#endif
		/* Create spline from the data points */
		if (grey->fit_rspl_w(grey, 0, lpnts, ngreyp, il, ih, &gres, ol, oh, 1.0, 0.005)) {
			fprintf(stderr,"Warning: Grey axis mapping is non-monotonic!\n");
		}
	}

#ifdef DEBUG
	{	/* Plot the 1D mapping */
		double xx[XRES];
		double y1[XRES];
		int i;

		for (i = 0; i < XRES; i++) {
			double x;
			co cp;		/* Conversion point */
			x = s_cs_bp[0] + (i/(double)(XRES-1)) * (s_cs_wp[0] - s_cs_bp[0]);
			xx[i] = x;
			cp.p[0] = x;
			grey->interp(grey, &cp);
			y1[i] = cp.v[0];
		}
		do_plot(xx,y1,NULL,NULL,XRES);
	}
#endif /* DEBUG */

	{
#ifdef INDEPL
		co cp;
		int i, ix;

		/* We want to remap L independently to everything else, */
		/* so transform source gamut through the L mapping */

		scl_gam = new_gamut(sc_gam->getsres(sc_gam));

		for (ix = 0;;) {
			double p[3];

			if ((ix = sc_gam->getrawvert(sc_gam, p, ix)) < 0)
				break;

			/* Remap L forwards through src to dst mapping */
			cp.p[0] = p[0];			/* L value */
			grey->interp(grey, &cp);
			p[0] = cp.v[0];
			scl_gam->expand(scl_gam, p);
		}

		if (sc_gam == s_gam)
			sl_gam = scl_gam;

		else {
			sl_gam = new_gamut(s_gam->getsres(s_gam));

			for (ix = 0;;) {
				double p[3];
	
				if ((ix = s_gam->getrawvert(s_gam, p, ix)) < 0)
					break;
	
				/* Remap L forwards through src to dst mapping */
				cp.p[0] = p[0];			/* L value */
				grey->interp(grey, &cp);
				p[0] = cp.v[0];
				sl_gam->expand(sl_gam, p);
			}
		}

		/* Setup starting point, L mapped source white points */
		for (i = 0; i < 3; i++) {
			sl_cs_wp[i] = s_cs_wp[i];
			sl_cs_bp[i] = s_cs_bp[i];
		}

		/* Create L mapped versions of src colorspace white/black points */
		cp.p[0] = sl_cs_wp[0];
		grey->interp(grey, &cp);
		sl_cs_wp[0] = cp.v[0];

		cp.p[0] = sl_cs_bp[0];
		grey->interp(grey, &cp);
		sl_cs_bp[0] = cp.v[0];
	}

#else
	scl_gam = sc_gam;		/* Use unmapped source gamut */
	sl_gam = s_gam;			/* Use unmapped image gamut */
#endif

	{
		int nspts;		/* Number of source gamut surface points */
		int ndpts;		/* Number of destination gamut surface points */
		int i, j;
		datai il, ih;
		datao ol, oh;
		nearsmth *nsm;	/* Returned list of near smooth points */
		int nnsm;		/* Number of near smoothed points */
#ifdef PLOT_GAMVEC
		FILE *wrl;
#endif

		nspts = sl_gam->nverts(sl_gam);		/* Source surface points */
		ndpts = d_gam->nverts(d_gam);		/* Destination surface points */

		if ((gpnts = (cow *)malloc((4 * mapres + 2 * (nspts + ndpts)) * sizeof(cow))) == NULL) { 
			fprintf(stderr,"gamut map: Malloc of mapping setup points failed\n");
			grey->del(grey);
			if (sl_gam != s_gam)
				sl_gam->del(sl_gam);
			free(s);
			return NULL;
		}

		/* ------------------------------------------- */
		/* Finish off the grey axis mapping by creating the */
		/* grey axis 3D->3D mapping points */
		/* Add in the grey axis mapping */
		/* We use 4 times the grid density, and create */
		/* points that span the source colorspace (this may exceed) */
		/* the source image gamut, and map to points outside the */
		/* destination gamut) */
		for (i = 0; i < (4 * mapres); i++) {
			double t;
			co cp;			/* Lookup point */

			/* Create source grey axis point */
			/* Range is nominalay to cover possible input points */
			t = i/(4.0 * mapres - 1.0);
			for (j = 0; j < 3; j++)
				gpnts[ngamp].p[j] = sl_cs_bp[j] + t * (sl_cs_wp[j] - sl_cs_bp[j]);

#ifdef INDEPL	/* L values are already mapped */
			gpnts[ngamp].v[0] = gpnts[ngamp].p[0];
#else
			/* Lookup destination L value (convert src L value to dst L value) */
			cp.p[0] = gpnts[ngamp].p[0];
			grey->interp(grey, &cp);
			gpnts[ngamp].v[0] = cp.v[0];
#endif

			/* Figure destination point on destination grey axis */
			t = (gpnts[ngamp].v[0] - d_gx_bp[0])/(d_gx_wp[0] - d_gx_bp[0]);
			for (j = 1; j < 3; j++)
				gpnts[ngamp].v[j] = d_gx_bp[j] + t * (d_gx_wp[j] - d_gx_bp[j]);
			
			gpnts[ngamp++].w  = 2.0;		/* Weighting grey axis slightly */
		}

		/* ---------------------------------------------------- */
		/* Now deal with the gamut edges. */
		/* For compression, create a mapping for each vertex of */
		/* the source gamut (image) surface towards the destination gamut */

		/* Create the near point mapping */
		nsm = near_smooth(verb, &nnsm, scl_gam, sl_gam, d_gam,
		    gampwf * P_A_WEIGHT   + gamswf * S_A_WEIGHT,
			gampwf * P_A_LWEIGHT  + gamswf * S_A_LWEIGHT,
			gampwf * P_A_CWEIGHT0 + gamswf * S_A_CWEIGHT0,
			gampwf * P_A_CWEIGHT1 + gamswf * S_A_CWEIGHT1,
			gampwf * P_A_HWEIGHT  + gamswf * S_A_HWEIGHT,

		    gampwf * P_R_WEIGHT  + gamswf * S_R_WEIGHT,
			gampwf * P_R_LWEIGHT + gamswf * S_R_LWEIGHT,
			gampwf * P_R_CWEIGHT + gamswf * S_R_CWEIGHT,
			gampwf * P_R_HWEIGHT + gamswf * S_R_HWEIGHT,

		    gampwf * P_D_WEIGHT  + gamswf * S_D_WEIGHT,
			gampwf * P_D_LWEIGHT + gamswf * S_D_LWEIGHT,
			gampwf * P_D_CWEIGHT + gamswf * S_D_CWEIGHT,
			gampwf * P_D_HWEIGHT + gamswf * S_D_HWEIGHT,

		    OVERSHOOT,
		    gamcpf > 1e-6, gamexf > 1e-6);
		if (nsm == NULL) {
			fprintf(stderr,"Creating smoothed near points failed\n");
			grey->del(grey);
			free(gpnts);
			if (sl_gam != s_gam)
				sl_gam->del(sl_gam);
			free(s);
			return NULL;
		}

#ifdef PLOT_GAMVEC
		wrl = start_vrml("gammap.wrl", 1);
		start_line_set(wrl);
#endif 

		for (i = 0; i < nnsm; i++) {
			double div[3];

			if (nsm[i].sr >= nsm[i].dr) {		/* Compression needed */

				/* Compute compression destination value */
				for (j = 0; j < 3; j++)				/* Compute compressed value */
#ifdef RADIAL_HACK	
					div[j] = gamcpf * nsm[i].drv[j] + (1.0 - gamcpf) * nsm[i].sv[j];
#else
					div[j] = gamcpf * nsm[i].sdv[j] + (1.0 - gamcpf) * nsm[i].sv[j];
#endif
//printf("Compression:\n");
//printf("Src point = %f %f %f radius %f\n",nsm[i].sv[0], nsm[i].sv[1], nsm[i].sv[2], nsm[i].sr);
//printf("Dst point = %f %f %f radius %f\n",nsm[i].sdv[0], nsm[i].sdv[1], nsm[i].sdv[2], nsm[i].dr);
//printf("Blended dst point = %f %f %f\n",div[0], div[1], div[2]);
//printf("\n");
				/* Set the mapping point */
				for (j = 0; j < 3; j++) {
					gpnts[ngamp].p[j] = nsm[i].sv[j];
					gpnts[ngamp].v[j] = div[j];
				}
				gpnts[ngamp++].w  = 1.0;

#ifdef PLOT_GAMVEC
				add_vertex(wrl, nsm[i].sv);
//				add_vertex(wrl, div);
				add_vertex(wrl, nsm[i].sdv);
#endif

#ifdef USE_GAMKNF
				/* Create a "knee" point */
				{
					double t;
					co cp;			/* Lookup point */
					double skp[3], dkp[3];

					/* Since the div is an associated point within both gamuts, */
					/* we're going to use it as a basis for the "half way" knee point. */
					/* To stop this point affecting the grey axis mapping, we need */
					/* to adjust the half way mapping sympatheticaly to the grey axis mapping */

					/* Find the div's corresponding point on the source L axis */
					t = (div[0] - sl_cs_bp[0])/(sl_cs_wp[0] - sl_cs_bp[0]);
					for (j = 0; j < 3; j++)
						skp[j] = sl_cs_bp[j] + t * (sl_cs_wp[j] - sl_cs_bp[j]);
			
#ifdef INDEPL
					dkp[0] = skp[0];	/* sl_gam is already L mapped */
#else
					/* Lookup destination L value */
					cp.p[0] = skp[0];
					grey->interp(grey, &cp);
					dkp[0] = cp.v[0];
#endif

					/* Figure destination point on destination grey axis */
					t = (dkp[0] - d_gx_bp[0])/(d_gx_wp[0] - d_gx_bp[0]);
					for (j = 1; j < 3; j++)
						dkp[j] = d_gx_bp[j] + t * (d_gx_wp[j] - d_gx_bp[j]);
					
					/* Make mapping point 1/2 way between div and corresponding grey axis point */
					for (j = 0; j < 3; j++) {
						gpnts[ngamp].p[j] = (skp[j] + div[j])/2.0;
						gpnts[ngamp].v[j] = (dkp[j] + div[j])/2.0;
					}
					gpnts[ngamp++].w = gamknf * gamknf;		/* Knee weight */
				}
#endif /* USE_GAMKNF */
			} else {	/* Expansion needed */

				/* Compute expansion destination value */
				for (j = 0; j < 3; j++)				/* Compute compressed value */
#ifdef RADIAL_HACK	
					div[j] = gamexf * nsm[i].drv[j] + (1.0 - gamexf) * nsm[i].sv[j];
#else
					div[j] = gamexf * nsm[i].sdv[j] + (1.0 - gamexf) * nsm[i].sv[j];
#endif
//printf("Expansion:\n");
//printf("Src point = %f %f %f radius %f\n",nsm[i].sv[0], nsm[i].sv[1], nsm[i].sv[2], nsm[i].sr);
//printf("Dst point = %f %f %f radius %f\n",nsm[i].sdv[0], nsm[i].sdv[1], nsm[i].sdv[2], nsm[i].dr);
//printf("Blended dst point = %f %f %f\n",div[0], div[1], div[2]);
//printf("\n");
				/* Set the mapping point */
				for (j = 0; j < 3; j++) {
					gpnts[ngamp].p[j] = nsm[i].sv[j];
					gpnts[ngamp].v[j] = div[j];
				}
				gpnts[ngamp++].w  = 1.0;
	
#ifdef PLOT_GAMVEC
				add_vertex(wrl, nsm[i].sv);
//				add_vertex(wrl, div);
				add_vertex(wrl, nsm[i].sdv);
#endif
	
#ifdef USE_GAMKNF
				/* Create a "knee" point */
				{
					double t;
					co cp;			/* Lookup point */
					double skp[3], dkp[3];
	
					/* Since the sv is an associated point within both gamuts, */
					/* we're going to use it as a basis for the "half way" knee point. */
					/* To stop this point affecting the grey axis mapping, we need */
					/* to adjust the half way mapping sympatheticaly to the grey axis mapping */

					/* Find the sv's corresponding point on the source L axis */
					t = (nsm[i].sv[0] - sl_cs_bp[0])/(sl_cs_wp[0] - sl_cs_bp[0]);
					for (j = 0; j < 3; j++)
						skp[j] = sl_cs_bp[j] + t * (sl_cs_wp[j] - sl_cs_bp[j]);
			
#ifdef INDEPL
					dkp[0] = skp[0];	/* sl is already L mapped */
#else
					/* Lookup destination L value */
					cp.p[0] = skp[0];
					grey->interp(grey, &cp);
					dkp[0] = cp.v[0];
#endif
			
					/* Figure destination point on destination grey axis */
					t = (dkp[0] - d_gx_bp[0])/(d_gx_wp[0] - d_gx_bp[0]);
					for (j = 1; j < 3; j++)
						dkp[j] = d_gx_bp[j] + t * (d_gx_wp[j] - d_gx_bp[j]);
					
					/* Make mapping point 1/2 way between sv and corresponding grey axis point */
					for (j = 0; j < 3; j++) {
						gpnts[ngamp].p[j] = (skp[j] + nsm[i].sv[j])/2.0;
						gpnts[ngamp].v[j] = (dkp[j] + nsm[i].sv[j])/2.0;
					}
					gpnts[ngamp++].w = gamknf * gamknf;		/* Knee weight */
				}
#endif /* USE_GAMKNF */
			}
		}

		free(nsm);

#ifdef PLOT_GAMVEC
		make_lines(wrl, 2);
		end_vrml(wrl);
#endif
		/* --------------------------- */
		/* Compute the bounding values */
		for (j = 0; j < 3; j++) {
			il[j] = ol[j] =  1e60;
			ih[j] = oh[j] = -1e60;
		}
		for (i = 0; i < ngamp; i++) {
			for (j = 0; j < 3; j++) {
				if (gpnts[i].p[j] < il[j])
					il[j] = gpnts[i].v[j];
				if (gpnts[i].p[j] > ih[j])
					ih[j] = gpnts[i].v[j];
				if (gpnts[i].v[j] < ol[j])
					ol[j] = gpnts[i].v[j];
				if (gpnts[i].v[j] > oh[j])
					oh[j] = gpnts[i].v[j];
			}
		}

		/* Override with argument input range 
		if (mn != NULL && mx != NULL) {
			for (j = 0; j < 3; j++) {
				if (mn[j] < il[j])
					il[j] = mn[j];
				if (mx[j] > ih[j])
					ih[j] = mx[j];
			}
		}

		/* Create the gamut mapping rspl. */
		map = new_rspl(3, 3);	/* Allocate 3D -> 3D */
	
		/* The mapping tends to get "lumpy" as we crank the resolution up, */
		/* so increase smoothing with resolution */
		{
			double smth;
			int gres[MXDI];

#ifdef NEVER
			smth = 1.0 + 7.0 * (mapres - 9.0)/(33.0 - 9.0);
			if (smth < 1.0)
				smth = 1.0;
			else if (smth > 10.0)
				smth = 10.0;
#else
			smth = 2.0;
#endif

#ifdef NEVER		/* Dump out all the mapping points */
			{
				for (i = 0; i < ngamp; i++) {
					printf("%d: %f %f %f -> %f %f %f\n",i,
						gpnts[i].p[0], gpnts[i].p[1], gpnts[i].p[2],
						gpnts[i].v[0], gpnts[i].v[1], gpnts[i].v[2]);
				}
			}
#endif
			for (j = 0; j < 3; j++)
				gres[j] = mapres;

			/* Create 3D->3D spline from the data points */
			if (map->fit_rspl_w(map, 0, gpnts, ngamp, il, ih, gres, ol, oh, smth, 0.005)) {
				fprintf(stderr,"Warning: Gamut mapping is non-monotonic!\n");
			}
		}

		/* Final "safety": */
		/* Go through all the rspl data points and clip them to the destination gamut */
		/* and make sure that surface points map outside destination ?? */
		{
			// ~~~ do we need something like this ???
		}

		/* Put the mapping rspl's into place in the gammap */
		{
#ifdef INDEPL
			s->grey = grey;
#else
			grey->del(grey);
			s->grey = NULL;
#endif
			s->map = map;
		}

		/* If requested, enhance the saturation of the output values. */
		if (satenh > 0.0) {
			adjustsat cx;		/* Adjustment context */

			/* Compute what our source white and black points actually maps to */
			s->domap(s, cx.wp, s_mt_wp);
			s->domap(s, cx.bp, s_mt_bp);

			cx.dst = d_gam; 
			cx.satenh = satenh; 

			/* Saturation enhance the output values */
			map->re_set_rspl(
				map,				/* this */
				0,					/* Combination of flags */
				(void *)&cx,		/* Opaque function context */
				adjust_sat_func /* Function to set from */
			);
		}

		/* Test the gamut white and black point mapping, and "fine tune" */
		/* the mapping, to ensure an accurate transform of the white */
		/* and black points to the destination colorspace. */
		/* This compensates for any inacuracy introduced in the */
		/* various rspl mappings. */
		{
			adjustwb cx;		/* Adjustment context */
			double a_wp[3];		/* actual white point */
			double a_bp[3];		/* actual black point */
			int j;
			double sv[3], dv[3];	/* Source and destination vectors */

			/* Check what the source white and black points actually maps to */
			s->domap(s, a_wp, s_mt_wp);
			s->domap(s, a_bp, s_mt_bp);

#ifdef VERBOSE
			printf("White is %f %f %f, should be %f %f %f\n",
			a_wp[0], a_wp[1], a_wp[2], d_mt_wp[0], d_mt_wp[1], d_mt_wp[2]);
			printf("Black is %f %f %f, should be %f %f %f\n",
			a_bp[0], a_bp[1], a_bp[2], d_mt_bp[0], d_mt_bp[1], d_mt_bp[2]);
#endif /* VERBOSE */

			/* Setup the fine tune target */
			for (j = 0; j < 3; j++) {
				cx.sbp[j] = a_bp[j];
				sv[j] = a_wp[j] - cx.sbp[j];

				cx.dbp[j] = d_mt_bp[j];
				dv[j] = d_mt_wp[j] - cx.dbp[j];
			}

			/* Compute rotation/scale relative white point matrix */
			icmRotMat(cx.mat, sv, dv);

			/* Fine tune the 3D->3D mapping */
			map->re_set_rspl(
				map,		/* this */
				0,					/* Combination of flags */
				(void *)&cx,		/* Opaque function context */
				adjust_wb_func /* Function to set from */
			);

#ifdef VERBOSE
			/* Check what the source white and black points actually maps to */
			s->domap(s, a_wp, s_mt_wp);
			s->domap(s, a_bp, s_mt_bp);

			printf("After fine tuning:\n");
			printf("White is %f %f %f, should be %f %f %f\n",
			a_wp[0], a_wp[1], a_wp[2], d_mt_wp[0], d_mt_wp[1], d_mt_wp[2]);
			printf("Black is %f %f %f, should be %f %f %f\n",
			a_bp[0], a_bp[1], a_bp[2], d_mt_bp[0], d_mt_bp[1], d_mt_bp[2]);
#endif /* VERBOSE */
		}

#ifdef PLOT_GAMUTS
		sc_gam->write_vrml(sc_gam, "src.wrl", 1);
		s_gam->write_vrml(s_gam, "img.wrl", 1);
		d_gam->write_vrml(d_gam, "dst.wrl", 1);
		sc_gam->write_trans_vrml(sc_gam, "gmsrc.wrl", 1, map_trans, s);
#endif
	}

	free(gpnts);
	if (scl_gam != sc_gam) {
		scl_gam->del(scl_gam);
		if (scl_gam != sl_gam)
			sl_gam->del(sl_gam);
	}

	return s;
}

#ifdef PLOT_GAMUTS

/* Debug */
static void map_trans(void *cntx, double out[3], double in[3]) {
	gammap *map = (gammap *)cntx;

	map->domap(map, out, in);
}

#endif

/* Object methods */
static void del_gammap(
gammap *s
) {
	if (s->grey != NULL)
		s->grey->del(s->grey);
	s->map->del(s->map);

	free(s);
}

/* Apply the gamut mapping to the given color value */
static void domap(
gammap *s,
double *out,
double *in
) {
	co cp;

	/* If there is a separate grey mapping */
	if (s->grey != NULL) {
		cp.p[0] = in[0];
		s->grey->interp(s->grey, &cp);
		cp.p[0] = cp.v[0];
		cp.p[1] = in[1];
		cp.p[2] = in[2];
		s->map->interp(s->map, &cp);
		out[0] = cp.v[0];
		out[1] = cp.v[1];
		out[2] = cp.v[2];

	/* else there is a single mapping */
	} else {
		cp.p[0] = in[0];
		cp.p[1] = in[1];
		cp.p[2] = in[2];
		s->map->interp(s->map, &cp);
		out[0] = cp.v[0];
		out[1] = cp.v[1];
		out[2] = cp.v[2];
	}
}


/* Function to pass to rspl to alter output values, */
/* to enhance the saturation. */
static void
adjust_sat_func(
	void *pp,			/* adjustsat structure */
	double *out,		/* output value to be adjusted */
	double *in			/* corresponding input value */
) {
	adjustsat *p = (adjustsat *)pp;
	double cp[3];		/* Center point */
	double rr, t1[3], p1;
	double t2[3], p2;

	/* Locate center point on the white/black axis corresponding to this color */
	cp[0] = out[0];
	rr = (out[0] - p->bp[0])/(p->wp[0] - p->bp[0]);	/* Relative location on the white/black axis */
	cp[1] = p->bp[1] + rr * (p->wp[1] - p->bp[1]);
	cp[2] = p->bp[2] + rr * (p->wp[2] - p->bp[2]);
	
	/* Locate the point on the destination gamut surface in the direction */
	/* from the center point to the point being processed. */
	if (p->dst->vector_isect(p->dst, cp, out, t2, t1, &p2, &p1) != 0) {

		if (p1 > 1.0) {		/* If this point is within gamut */
			double ep1, bf;

//printf("\n");
//printf("~1 cp %f %f %f input %f %f %f\n",cp[0],cp[1],cp[2], out[0], out[1], out[2]);
//printf("~1 min %f %f %f mint %f\n",t2[0],t2[1],t2[2],p2);
//printf("~1 max %f %f %f maxt %f\n",t1[0],t1[1],t1[2],p1);

			p1 = 1.0/p1;		/* Position of out from cp to t1 */

#ifdef NEVER
			/* Enhanced parameter value */
			ep1 = (p1 + p->satenh * p1)/(1.0 + p->satenh * p1);
			/* Make blend between linear p1 and enhanced p1, */
			/* to reduce effects on near neutrals. */
			p1 = (1.0 - p1) * p1 + p1 * ep1;
#else
			/* Compute Enhanced p1 */
			ep1 = (p1 + p->satenh * p1)/(1.0 + p->satenh * p1);

			/* Make blend factor between linear p1 and enhanced p1, */
			/* to reduce effects on near neutrals. */
			{
				double pp = 4.0;		/* Sets where the 50% transition is */
				double g = 2.0;			/* Sets rate of transition */
				double sec, vv = p1;
		
				vv = vv/(pp - pp * vv + 1.0);

				vv *= 2.0;
				sec = floor(vv);
				if (((int)sec) & 1)
					g = -g;				/* Alternate action in each section */
				vv -= sec;
				if (g >= 0.0) {
					vv = vv/(g - g * vv + 1.0);
				} else {
					vv = (vv - g * vv)/(1.0 - g * vv);
				}
				vv += sec;
				vv *= 0.5;

				bf = (vv + pp * vv)/(1.0 + pp * vv);
			}
			/* Do the blend */
			p1 = (1.0 - bf) * p1 + bf * ep1;
#endif
			/* Compute enhanced values position */
			out[0] = cp[0] + (t1[0] - cp[0]) * p1;
			out[1] = cp[1] + (t1[1] - cp[1]) * p1;
			out[2] = cp[2] + (t1[2] - cp[2]) * p1;
//printf("~1 output %f %f %f, param %f\n",out[0],out[1],out[2],p1);
		}
	}
}

/* Function to pass to rspl to re-set output values, */
/* to adjust the white and black points */
static void
adjust_wb_func(
	void *pp,			/* adjustwb structure */
	double *out,		/* output value to be adjusted */
	double *in			/* corresponding input value */
) {
	int f;
	adjustwb *p    = (adjustwb *)pp;
	double sl, dl, t1[3], t2[3];

	/* Do a linear mapping from swp -> dwp and sbp -> dbp */

	/* Translate input color to be black relative */
	for (f = 0; f < 3; f++)
		t1[f] = out[f] - p->sbp[f];

	/* Map relative white points */
	t2[0] = p->mat[0][0] * t1[0] + p->mat[0][1] * t1[1] + p->mat[0][2] * t1[2];
	t2[1] = p->mat[1][0] * t1[0] + p->mat[1][1] * t1[1] + p->mat[1][2] * t1[2];
	t2[2] = p->mat[2][0] * t1[0] + p->mat[2][1] * t1[1] + p->mat[2][2] * t1[2];

	/* Translate to the destination black point */
	for (f = 0; f < 3; f++)
		out[f] = t2[f] + p->dbp[f];
}

#if defined(PLOT_GAMVEC)
/* ------------------------------------------------ */
/* Some simple functions to do basix VRML work */

static int npoints = 0;
static int paloc = 0;
static struct { double pp[3]; } *pary;

static void Lab2RGB(double *out, double *in);

FILE *start_vrml(char *name, int doaxes) {
	FILE *wrl;
	struct {
		double x, y, z;
		double wx, wy, wz;
		double r, g, b;
	} axes[5] = {
		{ 0, 0,  50-GAMUT_LCENT,  2, 2, 100,  .7, .7, .7 },	/* L axis */
		{ 50, 0,  0-GAMUT_LCENT,  100, 2, 2,   1,  0,  0 },	/* +a (red) axis */
		{ 0, -50, 0-GAMUT_LCENT,  2, 100, 2,   0,  0,  1 },	/* -b (blue) axis */
		{ -50, 0, 0-GAMUT_LCENT,  100, 2, 2,   0,  1,  0 },	/* -a (green) axis */
		{ 0,  50, 0-GAMUT_LCENT,  2, 100, 2,   1,  1,  0 },	/* +b (yellow) axis */
	};
	int i;
	
	if ((wrl = fopen(name,"w")) == NULL)
		error("Error opening VRML file '%s'\n",name);

	npoints = 0;

	fprintf(wrl,"#VRML V2.0 utf8\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"# Created by the Argyll CMS\n");
	fprintf(wrl,"Transform {\n");
	fprintf(wrl,"children [\n");
	fprintf(wrl,"	NavigationInfo {\n");
	fprintf(wrl,"		type \"EXAMINE\"        # It's an object we examine\n");
	fprintf(wrl,"	} # We'll add our own light\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    DirectionalLight {\n");
	fprintf(wrl,"        direction 0 0 -1      # Light illuminating the scene\n");
	fprintf(wrl,"        direction 0 -1 0      # Light illuminating the scene\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    Viewpoint {\n");
	fprintf(wrl,"        position 0 0 340      # Position we view from\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");
	if (doaxes != 0) {
		fprintf(wrl,"# Lab axes as boxes:\n");
		for (i = 0; i < 5; i++) {
			fprintf(wrl,"Transform { translation %f %f %f\n", axes[i].x, axes[i].y, axes[i].z);
			fprintf(wrl,"\tchildren [\n");
			fprintf(wrl,"\t\tShape{\n");
			fprintf(wrl,"\t\t\tgeometry Box { size %f %f %f }\n",
			                  axes[i].wx, axes[i].wy, axes[i].wz);
			fprintf(wrl,"\t\t\tappearance Appearance { material Material ");
			fprintf(wrl,"{ diffuseColor %f %f %f} }\n", axes[i].r, axes[i].g, axes[i].b);
			fprintf(wrl,"\t\t}\n");
			fprintf(wrl,"\t]\n");
			fprintf(wrl,"}\n");
		}
		fprintf(wrl,"\n");
	}

	return wrl;
}

void
start_line_set(FILE *wrl) {

	fprintf(wrl,"\n");
	fprintf(wrl,"Shape {\n");
	fprintf(wrl,"  geometry IndexedLineSet { \n");
	fprintf(wrl,"    coord Coordinate { \n");
	fprintf(wrl,"	   point [\n");
}

void add_vertex(FILE *wrl, double pp[3]) {

	fprintf(wrl,"%f %f %f,\n",pp[1], pp[2], pp[0]-GAMUT_LCENT);
	
	if (paloc < (npoints+1)) {
		paloc = (paloc + 10) * 2;
		if (pary == NULL)
			pary = malloc(paloc * 3 * sizeof(double));
		else
			pary = realloc(pary, paloc * 3 * sizeof(double));

		if (pary == NULL)
			error ("Malloc failed");
	}
	pary[npoints].pp[0] = pp[0];
	pary[npoints].pp[1] = pp[1];
	pary[npoints].pp[2] = pp[2];
	npoints++;
}


void make_lines(FILE *wrl, int ppset) {
	int i, j;

	fprintf(wrl,"      ]\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"  coordIndex [\n");

	for (i = 0; i < npoints;) {
		for (j = 0; j < ppset; j++, i++) {
			fprintf(wrl,"%d, ", i);
		}
		fprintf(wrl,"-1,\n");
	}
	fprintf(wrl,"    ]\n");

	/* Color */
	fprintf(wrl,"            colorPerVertex TRUE\n");
	fprintf(wrl,"            color Color {\n");
	fprintf(wrl,"              color [			# RGB colors of each vertex\n");

	for (i = 0; i < npoints; i++) {
		double rgb[3], Lab[3];
		Lab[0] = pary[i].pp[0];
		Lab[1] = pary[i].pp[1];
		Lab[2] = pary[i].pp[2];
		Lab2RGB(rgb, Lab);
		fprintf(wrl,"                %f %f %f,\n", rgb[0], rgb[1], rgb[2]);
	}
	fprintf(wrl,"              ] \n");
	fprintf(wrl,"            }\n");
	/* End color */

	fprintf(wrl,"  }\n");
	fprintf(wrl,"} # end shape\n");

}

void end_vrml(FILE *wrl) {

	fprintf(wrl,"\n");
	fprintf(wrl,"  ] # end of children for world\n");
	fprintf(wrl,"}\n");

	if (fclose(wrl) != 0)
		error("Error closing VRML file\n");
}


/* Convert a gamut Lab value to an RGB value for display purposes */
static void
Lab2RGB(double *out, double *in) {
	double L = in[0], a = in[1], b = in[2];
	double x,y,z,fx,fy,fz;
	double R, G, B;

	/* Scale so that black is visible */
	L = L * (100 - 40.0)/100.0 + 40.0;

	/* First convert to XYZ using D50 white point */
	if (L > 8.0) {
		fy = (L + 16.0)/116.0;
		y = pow(fy,3.0);
	} else {
		y = L/903.2963058;
		fy = 7.787036979 * y + 16.0/116.0;
	}

	fx = a/500.0 + fy;
	if (fx > 24.0/116.0)
		x = pow(fx,3.0);
	else
		x = (fx - 16.0/116.0)/7.787036979;

	fz = fy - b/200.0;
	if (fz > 24.0/116.0)
		z = pow(fz,3.0);
	else
		z = (fz - 16.0/116.0)/7.787036979;

	x *= 0.9642;	/* Multiply by white point, D50 */
	y *= 1.0;
	z *= 0.8249;

	/* Now convert to sRGB values */
	R = x * 3.2410  + y * -1.5374 + z * -0.4986;
	G = x * -0.9692 + y * 1.8760  + z * 0.0416;
	B = x * 0.0556  + y * -0.2040 + z * 1.0570;

	if (R < 0.0)
		R = 0.0;
	else if (R > 1.0)
		R = 1.0;

	if (G < 0.0)
		G = 0.0;
	else if (G > 1.0)
		G = 1.0;

	if (B < 0.0)
		B = 0.0;
	else if (B > 1.0)
		B = 1.0;

	R = pow(R, 1.0/2.2);
	G = pow(G, 1.0/2.2);
	B = pow(B, 1.0/2.2);

	out[0] = R;
	out[1] = G;
	out[2] = B;
}

#endif





















