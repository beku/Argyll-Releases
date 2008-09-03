
/* 
 * Argyll Gamut Mapping Library
 *
 * Author:  Graeme W. Gill
 * Date:    1/10/00
 * Version: 2.00
 *
 * Copyright 2000 - 2006 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * For a discussion of gamut mapping strategy used,
 * see gammap.txt
 */

/*
 * TTBD:
 *       Improve error handling.
 */

#define VERBOSE			/* Print out extra interesting information when verbose */
#define XRES 100		/* Res of plot */
#undef PLOT_LMAP		/* Plot L map */
#undef PLOT_GAMVEC		/* Save the gamut mapping points as "gammap.wrl" */
#undef PLOT_GAMUTS		/* Save (part mapped) input and output gamuts as */
						/* src.wrl, img.wrl, dst.wrl, gmsrc.wrl */
#undef PLOT_3DKNEES		/* Plot the 3D compression knees */
#define CHECK_NEARMAP	/* Check how accurately near map vectors are represented by rspl */

#define USE_GLUMKNF			/* Enable luminence knee function points */
#define USE_GAMKNF			/* Enable gamut boundary knee function points */
#define USE_BOUND			/* Enable grid boundary anchor points */


/* Optional marker points for gamut mapping diagnosotic */
struct {
	int type;		/* 1 = src point (xlate), 2 = dst point (no xlate), 0 = end marker */
	double pos[3];	/* Position, (usually in Jab space) */
	double col[3];	/* RGB color */
} markers[] = {
	{ 0, },								/* End marker */
	{ 2, { 41.222695, 63.911823, 37.695310 }, { 0.0, 1.0, 0.3 } },	/* destination in green */
	{ 1, { 41.951770, 60.220284, 34.788195 }, { 1.0, 0.3, 0.3 } },	/* source in red (Red) */
	{ 2, { 41.222695, 63.911823, 37.695310 }, { 0.3, 1.3, 0.3 } },		/* Dest in green */
	{ 1, { 85.117353, -60.807580, -22.195118 }, { 0.3, 0.3, 1 } },		/* Cyan Source (Blue) */
	{ 2, { 61.661622, -38.164411, -18.090824 }, { 1.0, 0.3, 0.3 } },	/* CMYK destination (Red) */
	{ 0 }								/* End marker */
};

/* Degree to which the hue & saturation of the black point axes should be aligned: */
#define GREYBPHSMF 0.0

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "counters.h"
#include "icc.h"
#include "xicc.h"
#include "gamut.h"
#include "rspl.h"
#include "gammap.h"
#include "nearsmth.h"
#include "vrml.h"
#ifdef PLOT_LMAP
#include "plot.h"
#endif

/* Callback context for enhancing the saturation of the clut values */
typedef struct {
	gamut *dst;			/* Destination colorspace gamut */
	double wp[3], bp[3];/* Destination colorspace white and black points */
	double satenh;		/* Saturation engancement value */
} adjustsat;

/* Callback context for making clut relative to white and black points */
typedef struct {
	double mat[3][4];
} adjustwb;

static void inv_grey_func(void *pp, double *out, double *in);
static void adjust_wb_func(void *pp, double *out, double *in);
static void adjust_sat_func(void *pp, double *out, double *in);

/* The smoothed near weighting control values. */
/* These weightings setup the detailed behaviour of the */
/* gamut mapping for the fully perceptual and saturation intents. */

/* Perceptual mapping weights */
gammapweights pweights[] = {
	{
		gmm_default,	/* Non hue specific defaults */
	
		{		/* Weighting of absolute error of destination from source */
			1.0,	/* Absolute error overall weight */
			{
				1.4,	/* Absolute luminance error weight */
//				1.0,	/* (V0.7) Absolute chroma error weight */
//				1.0		/* (V0.7) Absolute hue error weight */
				0.8,	/* Absolute chroma error weight */
				1.4		/* Absolute hue error weight */
			}
		},
		{		/* Weighting of relative error of destination points to each */
				/* other, compared to source points to each other. */
			1.0,	/* Relative error overall weight */
			{
				15.0,	/* Relative luminance error weight */
				8.0,	/* Relative chroma error weight */
				5.0		/* Relative hue error weight */
			}
		},
		{		/* Weighting of error between destination point and source */
				/* point radially mapped to center of destination. */
			1.0,	/* Radial error overall weight */
			{
				0.3,	/* Radial luminance error weight */
				0.1,	/* Radial chroma error weight */
				0.2		/* Radial hue error weight */
			}
		},
	
//		0.18,			/* (V0.7) Cusp alignment weighting */
		0.10,			/* Cusp alignment weighting */
		60.0,			/* Cusp alignment minimum radius */

		5.0				/* Extra elevation */
	},
	{
		gmm_yellow,		/* Treat yellow differently, to get purer result. */
	
		{		/* Weighting of absolute error of destination from source */
			1.0,	/* Absolute error overall weight */
			{
				1.3,	/* Absolute luminance error weight */
				1.6,	/* Absolute chroma error weight */
				1.0		/* Absolute hue error weight */
			}
		},
		{		/* Weighting of relative error of destination points to each */
				/* other, compared to source points to each other. */
			1.0,	/* Relative error overall weight */
			{
				15.0,	/* Relative luminance error weight */
				8.0,	/* Relative chroma error weight */
				5.0		/* Relative hue error weight */
			}
		},
		{		/* Weighting of error between destination point and source */
				/* point radially mapped to center of destination. */
			1.0,	/* Radial error overall weight */
			{
				0.1,	/* Radial luminance error weight */
				0.1,	/* Radial chroma error weight */
				0.2		/* Radial hue error weight */
			}
		},
	
		0.7,			/* Cusp alignment weighting */
		60.0,			/* Cusp alignment minimum radius */

		15.0				/* Extra elevation */
	},
	{
		gmm_end
	}
};
double psmooth = 1.0;		/* Level of RSPL smoothing for perceptual */

/* Saturation mapping weights */
gammapweights sweights[] = {
	{
		gmm_default,	/* Non hue specific defaults */
	
		{		/* Weighting of absolute error of destination from source */
			1.0,	/* Absolute error overall weight */
			{
				1.3,	/* Absolute luminance error weight */
				1.5,	/* Absolute chroma error weight */
				0.5		/* Absolute hue error weight */
			}
		},
		{		/* Weighting of relative error of destination points to each */
				/* other, compared to source points to each other. */
			1.0,	/* Relative error overall weight */
			{
				20.0,	/* Relative luminance error weight */
				10.0,	/* Relative chroma error weight */
				7.0		/* Relative hue error weight */
			}
		},
		{		/* Weighting of error between destination point and source */
				/* point radially mapped to center of destination. */
			1.0,	/* Radial error overall weight */
			{
				0.3,	/* Radial luminance error weight */
				0.1,	/* Radial chroma error weight */
				0.2		/* Radial hue error weight */
			}
		},
	
		0.6,			/* Cusp alignment weighting */
		70.0,			/* Cusp alignment minimum radius */

		10.0				/* Extra elevation */
	},
	{
		gmm_yellow,		/* Treat yellow differently, to get purer result. */
	
		{		/* Weighting of absolute error of destination from source */
			1.0,	/* Absolute error overall weight */
			{
				1.3,	/* Absolute luminance error weight */
				3.0,	/* Absolute chroma error weight */
				0.5		/* Absolute hue error weight */
			}
		},
		{		/* Weighting of relative error of destination points to each */
				/* other, compared to source points to each other. */
			1.0,	/* Relative error overall weight */
			{
				20.0,	/* Relative luminance error weight */
				10.0,	/* Relative chroma error weight */
				7.0		/* Relative hue error weight */
			}
		},
		{		/* Weighting of error between destination point and source */
				/* point radially mapped to center of destination. */
			1.0,	/* Radial error overall weight */
			{
				0.1,	/* Radial luminance error weight */
				0.1,	/* Radial chroma error weight */
				0.2		/* Radial hue error weight */
			}
		},
	
		1.0,			/* Cusp alignment weighting */
		70.0,			/* Cusp alignment minimum radius */

		15.0				/* Extra elevation */
	},
	{
		gmm_end
	}
};
double ssmooth = 1.8;		/* Level of RSPL smoothing for saturation */

/*
 * Notes:
 *       The "knee" shape produced by the rspl (thin plate spline) code
 *       is not what one would expect for expansion. It is not
 *       symetrical with compression, and is less "sharp". This
 *       is due to the rspl "smoothness" criteria being based on
 *       grid value difference rather than smoothness being measured,
 *       as curvature. This means that the spline gets "stiffer" as
 *       it increases in slope.
 *       Possibly rspl could be improved in this respect ???
 *       (Doesn't matter for L compression now, because rspl is
 *       being inverted for expansion).
 */

static void del_gammap(gammap *s);
static void domap(gammap *s, double *out, double *in);
#ifdef PLOT_GAMUTS
static void map_trans(void *cntx, double out[3], double in[3]);
#endif

/* Return a gammap to map from the input space to the output space */
/* Return NULL on error. */
gammap *new_gammap(
	int verb,			/* Verbose flag */
	gamut *sc_gam,		/* Source colorspace gamut */
	gamut *si_gam,		/* Source image gamut (NULL if none) */
	gamut *d_gam,		/* Destination colorspace gamut */
	icxGMappingIntent *gmi,	/* Gamut mapping specification */
	int    mapres,		/* Gamut map resolution, typically 9 - 33 */
	double *mn,			/* If not NULL, set minimum mapping input range */
	double *mx,			/* for rspl grid. */
	char *diagname		/* If non-NULL, write a gamut mapping diagnostic WRL */
) {
	gmm_BPmap bph = gmm_bendBP;		/* Prefered algorithm */
//	gmm_BPmap bph = gmm_clipBP;		/* Alternatives tried */
//	gmm_BPmap bph = gmm_BPadpt;
//	gmm_BPmap bph = gmm_noBPadpt;

	gammap *s;			/* This */
	gamut *scl_gam;		/* Source colorspace gamut with L mapping applied */
	gamut *sil_gam;		/* Source image gamut with L mapping applied */

	double s_cs_wp[3];	/* Source colorspace white point */
	double s_cs_bp[3];	/* Source colorspace black point */
	double s_ga_wp[3];	/* Source (image) gamut white point */
	double s_ga_bp[3];	/* Source (image) gamut black point */
	double d_cs_wp[3];	/* Destination colorspace white point */
	double d_cs_bp[3];	/* Destination colorspace black point */

	double sr_cs_wp[3];	/* Source rotated colorspace white point */
	double sr_cs_bp[3];	/* Source rotated colorspace black point */
	double sr_ga_wp[3];	/* Source rotated (image) gamut white point */
	double sr_ga_bp[3];	/* Source rotated (image) gamut black point */
	double dr_cs_wp[3];	/* Target (gmi->greymf aligned) white point */
	double dr_cs_bp[3];	/* Target (gmi->greymf aligned) black point */
	double dr_be_bp[3];	/* Bend at end Target black point (Same as dr_cs_bp[] otherwise) */
						/* == end target destination black point */

	double sl_cs_wp[3];	/* Source rotated and L mapped colorspace white point */
	double sl_cs_bp[3];	/* Source rotated and L mapped colorspace black point */

	double s_mt_wp[3];	/* Overall source mapping target white point (used for finetune) */
	double s_mt_bp[3];	/* Overall source mapping target black point (used for finetune) */
	double d_mt_wp[3];	/* Overall destination mapping white point (used for finetune) */
	double d_mt_bp[3];	/* Overall destination mapping black point (used for finetune) */

	int defrgrid = 6;	/* mapping range surface default anchor point resolution */
	int nres = 512;		/* Neutral axis resolution */
	cow lpnts[10];		/* Mapping points to create grey axis map */
	int revrspl = 0;	/* Reverse grey axis rspl construction */
	int ngreyp = 0;		/* Number of grey axis mapping points */
	cow *gpnts;			/* Mapping points to create gamut mapping */
	int ngamp = 0;		/* Number of gamut mapping points */
	int j;

#if defined(PLOT_LMAP) || defined(PLOT_GAMVEC) || defined(PLOT_GAMUTS) || defined(PLOT_3DKNEES)
	fprintf(stderr,"##### A gammap.c PLOT is #defined ####\n");
#endif

	if (verb) {
		xicc_dump_gmi(gmi);
		printf("Gamut map resolution: %d\n",mapres);
		if (si_gam != NULL)
			printf("Image gamut supplied\n");
		switch(bph) {
			case gmm_clipBP:	printf("Neutral axis no-adapt extend and clip\n"); break;
			case gmm_BPadpt:	printf("Neutral axis fully adapt\n"); break;
			case gmm_bendBP:	printf("Neutral axis no-adapt extend and bend\n"); break;
			case gmm_noBPadpt:	printf("Neutral axis no-adapt\n"); break;
		}
	}

	/* Allocate the object */
	if ((s = (gammap *)calloc(1, sizeof(gammap))) == NULL)
		error("gammap: calloc failed on gammap object");

	/* Setup methods */
	s->del = del_gammap;
	s->domap = domap;

	/* Now create everything */

	/* Grab the white and black points */

	if (sc_gam->getwb(sc_gam, NULL, NULL, s_cs_wp, s_cs_bp)) {
		fprintf(stderr,"gamut map: Unable to read source colorspace white and black points\n");
		free(s);
		return NULL;
	}
	if (si_gam == NULL) {
		si_gam = sc_gam;		/* Source space is source gamut */
		for (j = 0; j < 3; j++) {
			s_ga_wp[j] = s_cs_wp[j];
			s_ga_bp[j] = s_cs_bp[j];
		}
	} else {

		if (si_gam->getwb(si_gam, NULL, NULL, s_ga_wp, s_ga_bp)) {
			fprintf(stderr,"gamut map: Unable to read source gamut white and black points\n");
			free(s);
			return NULL;
		}

		/* Guard against silliness. Image must be within colorspace */
		if (s_ga_wp[0] > s_cs_wp[0]) {
			int j;
			double t;
#ifdef VERBOSE
			if (verb)
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
			if (verb)
				printf("Fixing wayward image black point\n");
#endif
			t = (s_cs_bp[0] - s_ga_wp[0])/(s_ga_bp[0] - s_ga_wp[0]);
			for (j = 0; j < 3; j++)
				s_ga_bp[j] = s_ga_wp[j] + t * (s_ga_bp[j] - s_ga_wp[j]);
		}
	}

	if (d_gam->getwb(d_gam, NULL, NULL, d_cs_wp, d_cs_bp)) {
		fprintf(stderr,"gamut map: Unable to read destination white and black points\n");
		free(s);
		return NULL;
	}

#ifdef VERBOSE
	if (verb) {
		printf("Src colorspace white/black are %f %f %f, %f %f %f\n",
		s_cs_wp[0], s_cs_wp[1], s_cs_wp[2], s_cs_bp[0], s_cs_bp[1], s_cs_bp[2]);
	
		printf("Src gamut white/black are      %f %f %f, %f %f %f\n",
		s_ga_wp[0], s_ga_wp[1], s_ga_wp[2], s_ga_bp[0], s_ga_bp[1], s_ga_bp[2]);
	
		printf("Dst colorspace white/black are %f %f %f, %f %f %f\n",
		d_cs_wp[0], d_cs_wp[1], d_cs_wp[2], d_cs_bp[0], d_cs_bp[1], d_cs_bp[2]);
	}
#endif /* VERBOSE */

	/* ------------------------------------ */
	/* Figure out the destination grey axis alignment */
	/* This is all done using colorspace white & black points */
	{
		double t, svl, dvl;
		double wrot[3][3];			/* Rotation about 0,0,0 to match white points */
		double sswp[3], ssbp[3];	/* Temporary source white & black points */
		double fawp[3], fabp[3];	/* Fully adapted destination white & black */
		double hawp[3], habp[3];	/* Half adapted destination white & black */

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		/* The first task is to decide what our target destination */
		/* white and black points are going to be. */


		/* Figure out what our initial target destination white point is going to be: */

		/* Compute source white and black points with same L value as the destination */
		t = (d_cs_wp[0] - s_cs_bp[0])/(s_cs_wp[0] - s_cs_bp[0]);
		for (j = 0; j < 3; j++)
			sswp[j] = s_cs_bp[j] + t * (s_cs_wp[j] - s_cs_bp[j]);

		t = (d_cs_bp[0] - s_cs_wp[0])/(s_cs_bp[0] - s_cs_wp[0]);
		for (j = 0; j < 3; j++)
			ssbp[j] = s_cs_wp[j] + t * (s_cs_bp[j] - s_cs_wp[j]);

		/* The raw grey axis alignment target is a blend between the */
		/* source colorspace (NOT gamut) and the destination */
		/* colorspace. */

		for (j = 0; j < 3; j++) {
			dr_cs_wp[j] = gmi->greymf * d_cs_wp[j] + (1.0 - gmi->greymf) * sswp[j];
			dr_cs_bp[j] = gmi->greymf * d_cs_bp[j] + (1.0 - gmi->greymf) * ssbp[j];
		}

#ifdef VERBOSE
		if (verb) {
			printf("Target (blended) dst wp/bp   = %f %f %f, %f %f %f\n",
				dr_cs_wp[0], dr_cs_wp[1], dr_cs_wp[2], dr_cs_bp[0], dr_cs_bp[1], dr_cs_bp[2]);
		}
#endif /* VERBOSE */

		/* Compute full adaptation target destinations */
		for (j = 0; j < 3; j++) {
			fawp[j] = dr_cs_wp[j];			/* White fully adapted */
			fabp[j] = dr_cs_bp[j];			/* Black fully adapted */
		}

		/* Clip the target grey axis to the destination gamut */
		if (d_gam->vector_isect(d_gam, fabp, fawp, fabp, fawp, NULL, NULL) == 0)
			error("gamut: vector_isect failed!");

		/* To work around the problem that vector_isect() is not entirely accurate, */
		/* special case the situation where gmi->greymf == 1.0 */
		if (gmi->greymf > 0.99) {
			for (j = 0; j < 3; j++) {
				fawp[j] = d_cs_wp[j];
				fabp[j] = d_cs_bp[j];
			}
		}

		/* Compute half adapted (full white, not black) target destinations */
		for (j = 0; j < 3; j++)
			hawp[j] = dr_cs_wp[j];			/* White fully adapted */

		/* Compute the rotation matrix that maps the source white point */
		/* onto the target white point. */
		icmRotMat(wrot, sswp, dr_cs_wp);

		/* Compute the target black point as the rotated source black point */ 
		icmMulBy3x3(habp, wrot, s_cs_bp);
		
		/* Now intersect the target white and black points with the destination */
		/* colorspace gamut to arrive at the best possible in gamut values for */
		/* the target white and black points. */

		if (d_gam->vector_isect(d_gam, habp, hawp, habp, hawp, NULL, NULL) == 0)
			error("gamut: vector_isect failed!");

		/* To work around the problem that vector_isect() is not entirely accurate, */
		/* special case the situation where gmi->greymf == 1.0 */
		if (gmi->greymf > 0.99) {
			for (j = 0; j < 3; j++) {
				hawp[j] = d_cs_wp[j];
			}
		}

		/* Now decide the detail of the white and black alignment */
		if (bph == gmm_BPadpt) {	/* Adapt to destination white and black */

			/* Use the fully adapted white and black points */
			for (j = 0; j < 3; j++) {
				dr_cs_wp[j] = fawp[j];
				dr_cs_bp[j] = fabp[j];
			}

			/* Set bent black point target to be the same as our actual */
			/* black point target, so that the "bend" code does nothing. */ 
			for (j = 0; j < 3; j++)
				dr_be_bp[j] = dr_cs_bp[j];

		} else {					/* Adapt to destination white but not black */

			/* Use the half adapted white and black points */
			for (j = 0; j < 3; j++) {
				dr_cs_wp[j] = hawp[j];
				dr_cs_bp[j] = habp[j];
			}

#ifdef VERBOSE
		if (verb) {
			printf("Adapted target wp/bp         = %f %f %f, %f %f %f\n",
				dr_cs_wp[0], dr_cs_wp[1], dr_cs_wp[2], dr_cs_bp[0], dr_cs_bp[1], dr_cs_bp[2]);
		}
#endif
			if (bph == gmm_clipBP || bph == gmm_bendBP) {

				/* Extend the target black point to accomodate the */
				/* bent or clipped destination space L* range */
				if (fabp[0] < dr_cs_bp[0]) {
					t = (fabp[0] - dr_cs_wp[0])/(dr_cs_bp[0] - dr_cs_wp[0]);
					for (j = 0; j < 3; j++)
						dr_cs_bp[j] = dr_cs_wp[j] + t * (dr_cs_bp[j] - d_cs_wp[j]);
				}
			}
		
			if (bph == gmm_bendBP) {

				/* Set the destination "bent" black target value to */
				/* be the fully adapted target black point, for use */
				/* by the 3D neutral axis bend code. */
				for (j = 0; j < 3; j++)
					dr_be_bp[j] = fabp[j];

			} else {

				/* Set the bent black point target to be the same as our actual */
				/* black point target, so that the "bend" code does nothing. */ 
				for (j = 0; j < 3; j++)
					dr_be_bp[j] = dr_cs_bp[j];
			}
		}

#ifdef VERBOSE
		if (verb) {
			printf("Adapted & extended tgt wp/bp = %f %f %f, %f %f %f\n",
				dr_cs_wp[0], dr_cs_wp[1], dr_cs_wp[2], dr_cs_bp[0], dr_cs_bp[1], dr_cs_bp[2]);
			printf("Bend target                                              bp = %f %f %f\n",
		        dr_be_bp[0], dr_be_bp[1], dr_be_bp[2]);
		}
#endif /* VERBOSE */

		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
		/* Now we need to figure out what origin alignment is needed, as well as */
		/* making sure the vectors are the same length to avoid rescaling. */
		/* (Scaling is meant to be done with the L curve.) */

		/* Create temporary source white point that has the same L as the */
		/* target destination white point. */
		t = (dr_cs_wp[0] - s_cs_bp[0])/(s_cs_wp[0] - s_cs_bp[0]);
		for (j = 0; j < 3; j++)
			sswp[j] = s_cs_bp[j] + t * (s_cs_wp[j] - s_cs_bp[j]);

		/* Create temporary source black point that will form a vector to the src white */
		/* point with same length as the target destination black->white vector. */
		for (svl = dvl = 0.0, j = 0; j < 3; j++) {
			double tt;
			tt = sswp[j] - s_cs_bp[j]; 
			svl += tt * tt;
			tt = dr_cs_wp[j] - dr_cs_bp[j];
			dvl += tt * tt;
		}
		svl = sqrt(svl);
		dvl = sqrt(dvl);
		for (j = 0; j < 3; j++)
			ssbp[j] = sswp[j] + dvl/svl * (s_cs_bp[j] - sswp[j]); 

#ifdef VERBOSE
		if (verb) {
			printf("Rotate matrix src wp/bp      = %f %f %f, %f %f %f\n",
				sswp[0], sswp[1], sswp[2], ssbp[0], ssbp[1], ssbp[2]);
			printf("Rotate matrix dst wp/bp      = %f %f %f, %f %f %f\n",
				dr_cs_wp[0], dr_cs_wp[1], dr_cs_wp[2], dr_cs_bp[0], dr_cs_bp[1], dr_cs_bp[2]);
		}
#endif /* VERBOSE */

		/* Now create the general rotation and translation to map the source grey */
		/* axis to our destination grey axis. */
		icmVecRotMat(s->grot, sswp, ssbp, dr_cs_wp, dr_cs_bp);

		/* And create the inverse as well: */
		icmVecRotMat(s->igrot, dr_cs_wp, dr_cs_bp, sswp, ssbp);

		/* Create rotated versions of source colorspace & image white and */
		/* black points for use from now on, given that rotation will */
		/* be applied first to everything. */
		icmMul3By3x4(sr_cs_wp, s->grot, s_cs_wp);
		icmMul3By3x4(sr_cs_bp, s->grot, s_cs_bp);
		icmMul3By3x4(sr_ga_wp, s->grot, s_ga_wp);
		icmMul3By3x4(sr_ga_bp, s->grot, s_ga_bp);

#ifdef VERBOSE
		if (verb) {
			printf("Rotated source grey axis wp/bp %f %f %f, %f %f %f\n",
				sr_cs_wp[0], sr_cs_wp[1], sr_cs_wp[2], sr_cs_bp[0], sr_cs_bp[1], sr_cs_bp[2]);
			printf("Rotated gamut grey axis wp/bp  %f %f %f, %f %f %f\n",
				sr_ga_wp[0], sr_ga_wp[1], sr_ga_wp[2], sr_ga_bp[0], sr_ga_bp[1], sr_ga_bp[2]);
			printf("Destination axis target wp/bp  %f %f %f, %f %f %f\n",
				dr_cs_wp[0], dr_cs_wp[1], dr_cs_wp[2], dr_cs_bp[0], dr_cs_bp[1], dr_cs_bp[2]);
		}
#endif
	}

#ifdef NEVER
sr_cs_wp[0] = 100.0;
sr_cs_bp[0] = 30.0;
dr_cs_wp[0] = 80.0;
dr_cs_bp[0] = 10.0;
glumknf	= 1.0;
#endif /* NEVER */

	/* Create the mapping points needed to build the 1D L mapping rspl. */
	/* If we have a gamut (ie. image) range that is smaller than the */
	/* L range of the colorspace, then use its white and black L values */
	/* as the source to be compressed to the destination L range. */
	/* We expand only a colorspace range, not a gamut/image range. */
	{
		double swL, dwL;			/* Source and destination white point L */
		double sbL, dbL;			/* Source and destination black point L */
		int j;
		double t;

		/* Setup white point mapping */
		if (sr_cs_wp[0] <= dr_cs_wp[0]) {	/* Needs possible expansion */
			swL = sr_cs_wp[0];
			dwL = gmi->glumwexf * dr_cs_wp[0] + (1.0 - gmi->glumwexf) * sr_cs_wp[0];

		} else {
			if (sr_ga_wp[0] > dr_cs_wp[0]) {	/* Gamut or colorspace needs compression */
				
				swL = (1.0 - gmi->glumwcpf) * dr_cs_wp[0] + gmi->glumwcpf * sr_ga_wp[0];
				dwL = dr_cs_wp[0];

			} else {	/* Neither needed */
				swL = sr_ga_wp[0];
				dwL = sr_ga_wp[0];
			}
		}

		/* Setup black point mapping */
		if (sr_cs_bp[0] >= dr_cs_bp[0]) {	/* Needs possible expansion */
			sbL = sr_cs_bp[0];
			dbL = gmi->glumbexf * dr_cs_bp[0] + (1.0 - gmi->glumbexf) * sr_cs_bp[0];

		} else {
			if (sr_ga_bp[0] < dr_cs_bp[0]) {	/* Gamut or colorspace needs compression */
				
				sbL = (1.0 - gmi->glumbcpf) * dr_cs_bp[0] + gmi->glumbcpf * sr_ga_bp[0];
				dbL = dr_cs_bp[0];

			} else {	/* Neither needed */
				sbL = sr_ga_bp[0];
				dbL = sr_ga_bp[0];
			}
		}

		/* To ensure symetry between compression and expansion, always create RSPL */
		/* for compression and its inverse, and then swap grey and igrey rspl to compensate. */
		if ((dwL - dbL) > (swL - sbL))
			revrspl = 1;

		/* White point end */
		lpnts[ngreyp].p[0] = swL;
		lpnts[ngreyp].v[0] = dwL;
		lpnts[ngreyp++].w  = 10.0;		/* Must go through here */

		/* Black point end */
		lpnts[ngreyp].p[0] = sbL;
		lpnts[ngreyp].v[0] = dbL;
		lpnts[ngreyp++].w  = 10.0;		/* Must go through here */

//printf("~1 white loc %f, val %f\n",swL,dwL);
//printf("~1 black loc %f, val %f\n",sbL,dbL);

#ifdef USE_GLUMKNF
		if (gmi->glumknf < 0.05)
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
			double cppos = 0.50;		/* Center point ratio between black and white */
			double cplv;				/* Center point location and value */
			double kppos = 0.30;		/* Knee point ration between white/black & center */
			double kwl, kbl, kwv, kbv;	/* Knee point values and locations */
			double kwx, kbx;			/* Knee point extra  */
			

//printf("sbL = %f, swL = %f\n",sbL,swL);
//printf("dbL = %f, dwL = %f\n",dbL,dwL);

			/* Center point */
			cplv = cppos * (swL - sbL) + sbL;
//printf("~1 computed cplv = %f\n",cplv);

#ifdef NEVER	/* Don't use a center point */
			lpnts[ngreyp].p[0] = cplv;
			lpnts[ngreyp].v[0] = cplv;
			lpnts[ngreyp++].w  = 0.5;	
#endif
	
//printf("~1 black half diff = %f\n",dbL - sbL); 
//printf("~1 white half diff = %f\n",dwL - swL); 

			/* Knee point locations */
			kwl = kppos * (cplv - swL) + swL;
			kbl = kppos * (cplv - sbL) + sbL;
			
			/* Extra compression for white and black knees */
			kwx = 0.6 * (dbL - sbL) + 1.0 * (swL - dwL);
			kbx = 1.0 * (dbL - sbL) + 0.6 * (swL - dwL);

//kwx = 0.0;
//kbx = 0.0;
//glumknf = 0.0;

			/* Knee point values */
			kwv = (dwL + kwx - cplv) * (kwl - cplv)/(swL - cplv) + cplv;
			if (kwv > dwL)		/* Sanity check */
				kwv = dwL;

			kbv = (dbL - kbx - cplv) * (kbl - cplv)/(sbL - cplv) + cplv;
			if (kbv < dbL)		/* Sanity check */
				kbv = dbL;


//printf("~1 kbl = %f, kbv = %f\n",kbl, kbv);
//printf("~1 kwl = %f, kwv = %f\n",kwl, kwv);

			/* Emphasise points to cause "knee" curve */
			lpnts[ngreyp].p[0] = kwl;
			lpnts[ngreyp].v[0] = kwv;
			lpnts[ngreyp++].w  = gmi->glumknf * gmi->glumknf;	
		
			lpnts[ngreyp].p[0] = kbl;
			lpnts[ngreyp].v[0] = kbv;
			lpnts[ngreyp++].w  = gmi->glumknf * gmi->glumknf;	
		}
#endif /* USE_GLUMKNF */

		/* Remember our source and destinatio mapping targets */
		/* so that we can use them for fine tuning later. */

		/* We scale the source and target white and black */
		/* points to match the L values of the source and destination */
		/* L curve mapping, as this is how we have chosen the */
		/* white and black point mapping for the link. */
		/* Put them back in pre-rotated space, so that we can */
		/* check the overall transform of the white and black points. */
		t = (swL - sr_cs_bp[0])/(sr_cs_wp[0] - sr_cs_bp[0]);
		for (j = 0; j < 3; j++)
			s_mt_wp[j] = sr_cs_bp[j] + t * (sr_cs_wp[j] - sr_cs_bp[j]);
		icmMul3By3x4(s_mt_wp, s->igrot, s_mt_wp);

		t = (sbL - sr_cs_wp[0])/(sr_cs_bp[0] - sr_cs_wp[0]);
		for (j = 0; j < 3; j++)
			s_mt_bp[j] = sr_cs_wp[j] + t * (sr_cs_bp[j] - sr_cs_wp[j]);
//printf("~1 check black point rotated = %f %f %f\n",s_mt_bp[0],s_mt_bp[1],s_mt_bp[2]);
		icmMul3By3x4(s_mt_bp, s->igrot, s_mt_bp);
//printf("~1 check black point prerotated = %f %f %f\n",s_mt_bp[0],s_mt_bp[1],s_mt_bp[2]);

		t = (dwL - dr_cs_bp[0])/(dr_cs_wp[0] - dr_cs_bp[0]);
		for (j = 0; j < 3; j++)
			d_mt_wp[j] = dr_cs_bp[j] + t * (dr_cs_wp[j] - dr_cs_bp[j]);

		/* Overal black point takes into account possible bend to fully adapted bp */
		t = (dbL - dr_cs_wp[0])/(dr_be_bp[0] - dr_cs_wp[0]);
		for (j = 0; j < 3; j++)
			d_mt_bp[j] = dr_cs_wp[j] + t * (dr_be_bp[j] - dr_cs_wp[j]);
	}

	/* We now create the 1D rspl L map, that compresses or expands the luminence */
	/* range, independent of grey axis alignment, or gamut compression. */
	/* Because the rspl isn't symetrical when we swap X & Y, and we would */
	/* like a conversion from profile A to B to be the inverse of profile B to A */
	/* (as much as possible), we contrive here to always create a compression */
	/* RSPL, and create an inverse for it, and swap the two of them so that */
	/* the transform is correct and has an accurate inverse available. */
	{
		datai il, ih;
		datao ol, oh;
		double avgdev[MXDO];
		int gres = 256;

		if (revrspl) {		/* Invert creation and usage for symetry between compress and exp. */
			int i;
			for (i = 0; i < ngreyp; i++) {
				double tt = lpnts[i].p[0];			/* Swap source and dest */
				lpnts[i].p[0] = lpnts[i].v[0];
				lpnts[i].v[0] = tt;
			}
		}

		/* Create a 1D rspl, that is used to */
		/* form the overall L compression mapping. */
		s->grey = new_rspl(RSPL_NOFLAGS, 1, 1);	/* Allocate 1D -> 1D */
	
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
		/* Create spline from the data points, with appropriate smoothness. */
		avgdev[0] = 0.005;
		if (s->grey->fit_rspl_w(s->grey, 0, lpnts, ngreyp, il, ih, &gres, ol, oh, 5.0, avgdev, NULL)) {
			fprintf(stderr,"Warning: Grey axis mapping is non-monotonic - may not be very smooth !\n");
		}

		/* Create an inverse mapping too, for reverse gamut and/or expansion. */
		il[0] = -1.0;		/* Set possible input range */
		ih[0] = 101.0;
		ol[0] = 0.0;		/* Set normalisation output range */
		oh[0] = 100.0;

		s->igrey = new_rspl(RSPL_NOFLAGS, 1, 1);	/* Allocate 1D -> 1D */
		/* Create it from inverse lookups of s->grey */
		s->igrey->set_rspl(s->igrey, 0, (void *)s->grey, inv_grey_func, il, ih, &gres, ol, oh);

		if (revrspl) {		/* Swap to compensate for expansion */
			rspl *tt = s->grey;
			s->grey = s->igrey;
			s->igrey = tt;
		}
	}

#ifdef PLOT_LMAP
	{	/* Plot the 1D mapping */
		double xx[XRES];
		double y1[XRES];
		int i;

		for (i = 0; i < XRES; i++) {
			double x;
			co cp;		/* Conversion point */
			x = sr_cs_bp[0] + (i/(double)(XRES-1)) * (sr_cs_wp[0] - sr_cs_bp[0]);
			xx[i] = x;
			cp.p[0] = x;
			s->grey->interp(s->grey, &cp);
			y1[i] = cp.v[0];
		}
		do_plot(xx,y1,NULL,NULL,XRES);
	}
#endif /* PLOT_LMAP */

	{
		co cp;
		double cusps[6][3];
		double wp[3], bp[3];
		double t;
		int i, ix;

		/* We want to rotate and then map L independently of everything else, */
		/* so transform source gamut through the rotation and L mapping */

		/* Create L mapped versions of rotated src colorspace white/black points */
		cp.p[0] = sr_cs_wp[0];
		s->grey->interp(s->grey, &cp);

		t = (cp.v[0] - sr_cs_bp[0])/(sr_cs_wp[0] - sr_cs_bp[0]);
		for (j = 0; j < 3; j++)
			sl_cs_wp[j] = sr_cs_bp[j] + t * (sr_cs_wp[j] - sr_cs_bp[j]);

		cp.p[0] = sr_cs_bp[0];
		s->grey->interp(s->grey, &cp);
		t = (cp.v[0] - sr_cs_wp[0])/(sr_cs_bp[0] - sr_cs_wp[0]);
		for (j = 0; j < 3; j++)
			sl_cs_bp[j] = sr_cs_wp[j] + t * (sr_cs_bp[j] - sr_cs_wp[j]);

#ifdef VERBOSE
		if (verb) {
			printf("Mapped source grey axis wp/bp  %f %f %f, %f %f %f\n",
				sl_cs_wp[0], sl_cs_wp[1], sl_cs_wp[2], sl_cs_bp[0], sl_cs_bp[1], sl_cs_bp[2]);
		}
#endif

		scl_gam = new_gamut(sc_gam->getsres(sc_gam), sc_gam->getisjab(sc_gam));

		for (ix = 0;;) {
			double p[3];

			if ((ix = sc_gam->getrawvert(sc_gam, p, ix)) < 0)
				break;

			/* Rotate and map gamut surface values */
			icmMul3By3x4(p, s->grot, p);
			cp.p[0] = p[0];			/* L value */
			s->grey->interp(s->grey, &cp);
			p[0] = cp.v[0];
			scl_gam->expand(scl_gam, p);
		}
		/* Translate cusps */
		if (sc_gam->getcusps(sc_gam, cusps) == 0) {
			scl_gam->setcusps(scl_gam, 0, NULL);
			for (i = 0; i < 6; i++) {
				double p[3];
				icmMul3By3x4(p, s->grot, cusps[i]);
				cp.p[0] = p[0];			/* L value */
				s->grey->interp(s->grey, &cp);
				p[0] = cp.v[0];
				scl_gam->setcusps(scl_gam, 1, p);
			}
			scl_gam->setcusps(scl_gam, 2, NULL);
		}
		/* Translate white and black points */
		if (sc_gam->getwb(sc_gam, wp, bp, NULL, NULL) == 0) {
			icmMul3By3x4(wp, s->grot, wp);
			cp.p[0] = wp[0];			/* L value */
			s->grey->interp(s->grey, &cp);
			wp[0] = cp.v[0];

			icmMul3By3x4(bp, s->grot, bp);
			cp.p[0] = bp[0];			/* L value */
			s->grey->interp(s->grey, &cp);
			bp[0] = cp.v[0];

			scl_gam->setwb(scl_gam, wp, bp);
		}

		if (sc_gam == si_gam)
			sil_gam = scl_gam;

		else {
			sil_gam = new_gamut(si_gam->getsres(si_gam), si_gam->getisjab(si_gam));

			for (ix = 0;;) {
				double p[3];
	
				if ((ix = si_gam->getrawvert(si_gam, p, ix)) < 0)
					break;
	
				/* Rotate and map gamut surface values */
				icmMul3By3x4(p, s->grot, p);
				cp.p[0] = p[0];			/* L value */
				s->grey->interp(s->grey, &cp);
				p[0] = cp.v[0];
				sil_gam->expand(sil_gam, p);
			}
			/* Cusps, w & b points for image gamut aren't used by nearsmth */
		}
	}

	/* Create all the 3D->3D gamut mapping points and 3D rspl */
	{
		int nspts;		/* Number of source gamut surface points */
		int rgridpts;	/* Number of range surface grid points */
		int i, j;
		datai il, ih;
		datao ol, oh;
		int gres[MXDI];
		double avgdev[MXDO];
		nearsmth *nsm;			/* Returned list of near smooth points */
		int nnsm;				/* Number of near smoothed points */
		double brad = 0.0;		/* Black bend radius */
		gammapweights xpweights[7], xsweights[7];	/* Explicit perceptial and sat. weights */
		gammapweights xwh[7]; 	/* Structure holding blended weights */
		double smooth = 1.0;	/* Level of 3D RSPL smoothing, blend of psmooth and ssmooth */
		vrml *wrl = NULL;

#ifdef PLOT_3DKNEES
typedef struct {
	double v0[3], v1[3];
} p3dk_lpoint;
		p3dk_lpoint *p3dk_locus;
		int p3dk_ix = 0;
#endif /* PLOT_3DKNEES */

#ifdef NEVER
		nspts = scl_gam->nverts(scl_gam);	/* Source surface points */
#else
		nspts = scl_gam->nraw0verts(scl_gam);	/* Source surface points */
#endif

		rgridpts = 0;
		if (defrgrid >= 2) {
			rgridpts = defrgrid * defrgrid * defrgrid
			         - (defrgrid -2) * (defrgrid -2) * (defrgrid -2);
		}
		if ((gpnts = (cow *)malloc((nres + 2 * nspts + rgridpts) * sizeof(cow))) == NULL) { 
			fprintf(stderr,"gamut map: Malloc of mapping setup points failed\n");
			s->grey->del(s->grey);
			s->igrey->del(s->igrey);
			if (sil_gam != scl_gam)
				sil_gam->del(sil_gam);
			scl_gam->del(scl_gam);
			free(s);
			return NULL;
		}

#ifdef PLOT_3DKNEES
		if ((p3dk_locus = (p3dk_lpoint *)malloc((2 * nspts) * sizeof(p3dk_lpoint))) == NULL)
			error("gamut: Diagnostic array p3dk_locus malloc failed");
#endif /* PLOT_3DKNEES */

		/* ------------------------------------------- */
		/* Finish off the grey axis mapping by creating the */
		/* grey axis 3D->3D mapping points */
		/* We use 4 times the grid density, and create */
		/* points that span the source colorspace (this may exceed) */
		/* the source image gamut, and map to points outside the */
		/* destination gamut) */

		/* See how much to bend the black - compute the color difference */
		/* (brad will be 0 for non gmm_bendBP because dr_be_bp dr_cs_bp */
		for (brad = 0.0, i = 1; i < 3; i++) {
			double tt = dr_cs_bp[i] - dr_be_bp[i];
			brad += tt * tt;
		}
		brad = sqrt(brad);

		for (i = 0; i < nres; i++) {
			double t;
			double dv[3];	/* Straight destination value */
			double bv[3];	/* Bent destination value */

			/* Create source grey axis point */
			t = i/(nres - 1.0);

			/* Cover 0.0 to 100.0 */
			t = ((100.0 * t) - sl_cs_bp[0])/(sl_cs_wp[0] - sl_cs_bp[0]);
			for (j = 0; j < 3; j++)
				gpnts[ngamp].p[j] = sl_cs_bp[j] + t * (sl_cs_wp[j] - sl_cs_bp[j]);

			/* L values are the same, as they have been mapped prior to 3D */
			gpnts[ngamp].v[0] = gpnts[ngamp].p[0];

			/* Figure destination point on destination grey axis */
			t = (gpnts[ngamp].v[0] - dr_cs_wp[0])/(dr_cs_bp[0] - dr_cs_wp[0]);
			for (j = 0; j < 3; j++)
				dv[j] = dr_cs_wp[j] + t * (dr_cs_bp[j] - dr_cs_wp[j]);
			
			/* Figure destination point on bent grey axis */
			t = (gpnts[ngamp].v[0] - dr_cs_wp[0])/(dr_be_bp[0] - dr_cs_wp[0]);
			for (j = 0; j < 3; j++)
				bv[j] = dr_cs_wp[j] + t * (dr_be_bp[j] - dr_cs_wp[j]);
			
			/* Figure out a blend value between the straight value */
			/* and the bent value, so that it curves smoothly from */
			/* one to the other. (This seems to look worse, so not used.) */
			if (brad > 0.001) {
				double ty;
				t = ((dr_cs_bp[0] + brad)  - gpnts[ngamp].v[0])/brad;
				if (t < 0.0)
					t = 0.0;
				else if (t > 1.0)
					t = 1.0;
				/* Make it a spline ? */
				t = t * t * (3.0 - 2.0 * t);
				ty = t * t * (3.0 - 2.0 * t);	/* spline blend value */
				t = (1.0 - t) * ty + t * t;		/* spline at t == 0, linear at t == 1 */
			} else {
				t = 0.0;	/* stick to straight, it will be close anyway. */
			}

			for (j = 0; j < 3; j++)		/* full bend when t == 1 */
				gpnts[ngamp].v[j] = t * bv[j] + (1.0 - t) * dv[j];

#ifdef NEVER
			printf("Grey axis %d maps %f %f %f -> %f %f %f\n",ngamp,
			gpnts[ngamp].p[0], gpnts[ngamp].p[1], gpnts[ngamp].p[2],
			gpnts[ngamp].v[0], gpnts[ngamp].v[1], gpnts[ngamp].v[2]);
#endif
			gpnts[ngamp++].w = 1.0;		/* Weighting grey axis */
		}

		/* ---------------------------------------------------- */
		/* Now deal with the gamut edges. */
		/* For compression, create a mapping for each vertex of */
		/* the source gamut (image) surface towards the destination gamut */
		/* For expansion, do the opposite. */

		/* Convert from compact to explicit hextant weightings */
		expand_weights(xpweights, pweights);
		expand_weights(xsweights, sweights);

		/* Create weights as blend between perceptual and saturation */
		near_xwblend(xwh, xpweights, gmi->gampwf, xsweights, gmi->gamswf);
		if ((gmi->gampwf + gmi->gamswf) > 0.1)
			smooth = (gmi->gampwf * psmooth) + (gmi->gamswf * ssmooth);

		/* Create the near point mapping, which is our fundamental gamut */
		/* hull to gamut hull mapping. */
		nsm = near_smooth(verb, &nnsm, scl_gam, sil_gam, d_gam, dr_be_bp, xwh, 
		    gmi->gamcpf > 1e-6, gmi->gamexf > 1e-6);
		if (nsm == NULL) {
			fprintf(stderr,"Creating smoothed near points failed\n");
			free(gpnts);
			s->grey->del(s->grey);
			s->igrey->del(s->igrey);
			if (sil_gam != scl_gam)
				sil_gam->del(sil_gam);
			scl_gam->del(scl_gam);
			free(s);
			return NULL;
		}
		
		if (diagname != NULL)
			wrl = new_vrml(diagname, 1);
#ifdef PLOT_GAMVEC
		else
			wrl = new_vrml("gammap.wrl", 1);
#endif 

		if (wrl != NULL) {
			for (i = 0; ; i++) {	/* Add diagnostic markers */
				double pp[3];
				co cp;
				if (markers[i].type == 0)
					break;
	
				if (markers[i].type == 1) {
					/* Rotate and map marker points the same as the src gamuts */
					icmMul3By3x4(pp, s->grot, markers[i].pos);
					cp.p[0] = pp[0];			/* L value */
					s->grey->interp(s->grey, &cp);
					pp[0] = cp.v[0];
				} else {
					pp[0] = markers[i].pos[0];
					pp[1] = markers[i].pos[1];
					pp[2] = markers[i].pos[2];
				}
				wrl->add_marker(wrl, pp, markers[i].col, 1.0);
			}
			wrl->start_line_set(wrl);
		}

		/* --------------------------- */
		/* Compute the input bounding values */
		for (j = 0; j < 3; j++) {
			il[j] = ol[j] =  1e60;
			ih[j] = oh[j] = -1e60;
		}
		/* From grey axis points */
		for (i = 0; i < ngamp; i++) {
			for (j = 0; j < 3; j++) {
				if (gpnts[i].p[j] < il[j])
					il[j] = gpnts[i].p[j];
				if (gpnts[i].p[j] > ih[j])
					ih[j] = gpnts[i].p[j];
			}
		}

		/* From source gamut boundary of near point mapping, */
		/* we compute the input range here, so that we can add */
		/* mapping points on the surface of the range cube. */
		for (i = 0; i < nnsm; i++) {
			for (j = 0; j < 3; j++) {
				if (nsm[i].sv[j] < il[j])
					il[j] = nsm[i].sv[j];;
				if (nsm[i].sv[j] > ih[j])
					ih[j] = nsm[i].sv[j];
			}
		}

#ifdef NEVER
		if (verb) {
			fprintf(stderr,"Input bounding box after grey axis and source mapping points:\n");
			fprintfstderr,("%f -> %f, %f -> %f, %f -> %f\n",
			il[0], ih[0], il[1], ih[1], il[2], ih[2]);
		}
#endif
		/* Expand to make sure grid covers the entire source gamut */
		{
			double tmx[3], tmn[3];
			sc_gam->getrange(sc_gam, tmn, tmx);
			for (j = 0; j < 3; j++) {
				if (tmn[j] < il[j])
					il[j] = tmn[j];
				if (tmx[j] > ih[j])
					ih[j] = tmx[j];
			}
		}

#ifdef NEVER
		if (verb) {
			fprintf(stderr,"Input bounding box after input colorspace boundary:\n");
			fprintf(stderr,"%f -> %f, %f -> %f, %f -> %f\n",
			il[0], ih[0], il[1], ih[1], il[2], ih[2]);
		}
#endif

		/* Expand to input range given as input arguments */ 
		if (mn != NULL && mx != NULL) {

			for (j = 0; j < 3; j++) {
				if (mn[j] < il[j])
					il[j] = mn[j];
				if (mx[j] > ih[j])
					ih[j] = mx[j];
			}

#ifdef NEVER
			if (verb) {
				fprintf(stderr,"After ovverride, input bounding box for 3D gamut mapping is:\n");
				fprintf(stderr,"%f -> %f, %f -> %f, %f -> %f\n",
				il[0], ih[0], il[1], ih[1], il[2], ih[2]);
			}
#endif
		}

		/* Now expand the bounding box by aprox 5% margin, but scale grid res */
		/* to match, so that the natural or given boundary still lies on the grid. */
		{
			int xmapres;
			double scale;

			xmapres = (int) ((mapres-1) * 0.05 + 0.5);
			if (xmapres < 1)
				xmapres = 1;

			scale = (double)(mapres-1 + xmapres)/(double)(mapres-1);

			for (j = 0; j < 3; j++) {
				double low, high;
				high = ih[j];
				low = il[j];
				ih[j] = (scale * (high - low)) + low;
				il[j] = (scale * (low - high)) + high;
			}

			mapres += 2 * xmapres;
#ifdef NEVER
			if (verb) {
				fprintf(stderr,"After incresing mapres to %d, input bounding box for 3D gamut mapping is:\n",mapres);
				fprintf(stderr,"%f -> %f, %f -> %f, %f -> %f\n",
				il[0], ih[0], il[1], ih[1], il[2], ih[2]);
			}
#endif
		}

		/* --------------------------- */
		/* Now computue our 3D mapping points from the near point mapping */
		/* NOTE: the knee points aren't directly connected with the gamut shell point, */
		/* since they are positioned at a constant L from them, rather than */
		/* in the gamut mapping direction or spherically, it's just a convenient way */
		/* of locating a set of inner shell points that have a weaker 1:1 mapping */
		/* that (hopefully) makes any compression/expansion progressive, rather than linear. */
		for (i = 0; i < nnsm; i++) {
			double div[3];			/* sdv moderated by gamcpf */
			double knpos = 0.8;		/* Gamut Knee position */

			if (nsm[i].sr >= nsm[i].dr) {		/* Compression needed */

				/* Compute compression destination value which is a blend */
				/* between the source value and the fully mapped destination value. */
				for (j = 0; j < 3; j++)				/* Compute compressed value */
					div[j] = gmi->gamcpf * nsm[i].sdv[j] + (1.0 - gmi->gamcpf) * nsm[i].sv[j];

#ifdef NEVER
				printf("Compression:\n");
				printf("Src point = %f %f %f radius %f\n",nsm[i].sv[0], nsm[i].sv[1], nsm[i].sv[2], nsm[i].sr);
				printf("Dst point = %f %f %f radius %f\n",nsm[i].sdv[0], nsm[i].sdv[1], nsm[i].sdv[2], nsm[i].dr);
				printf("Blended dst point = %f %f %f\n",div[0], div[1], div[2]);
				printf("\n");
#endif	/* NEVER */
				/* Set the main gamut hull mapping point */
				for (j = 0; j < 3; j++) {
					gpnts[ngamp].p[j] = nsm[i].sv[j];
					gpnts[ngamp].v[j] = div[j];
				}
				gpnts[ngamp++].w  = 1.0;

				if (wrl != NULL) {
					wrl->add_vertex(wrl, nsm[i].sv);	/* Source value */
//					wrl->add_vertex(wrl, div);			/* Blended destination value */
					wrl->add_vertex(wrl, nsm[i].sdv);	/* Full compression dest. value */
				}

				/* Create an inner "knee" point and an outer boundary point. */
				/* (Given the large smoothness factor used for the 3D */
				/*  mapping, the "knee" is probably not being very effective here ?) */
				{
					double t;
					double sxkp[3], dxkp[3];	/* Axis intercept of knee point direction */
					double skp[3], dkp[3];		/* Inner knee point */

					/* Since the div is an associated point within both gamuts, */
					/* we're going to use it as a basis for the "half way" knee point. */
					/* To stop this point affecting the grey axis mapping, we need */
					/* to adjust the half way mapping sympatheticaly to the grey axis mapping */

					/* Find the div's corresponding point on the source L axis */
					t = (div[0] - sl_cs_bp[0])/(sl_cs_wp[0] - sl_cs_bp[0]);
					for (j = 0; j < 3; j++)
						sxkp[j] = sl_cs_bp[j] + t * (sl_cs_wp[j] - sl_cs_bp[j]);
			
					dxkp[0] = sxkp[0];	/* sil_gam is already L mapped */

					/* Figure destination point on destination grey axis */
					t = (dxkp[0] - dr_cs_bp[0])/(dr_cs_wp[0] - dr_cs_bp[0]);
					for (j = 1; j < 3; j++)
						dxkp[j] = dr_cs_bp[j] + t * (dr_cs_wp[j] - dr_cs_bp[j]);
					
#ifdef PLOT_3DKNEES
					p3dk_locus[p3dk_ix].v0[0] = sxkp[0]; 
					p3dk_locus[p3dk_ix].v0[1] = sxkp[1]; 
					p3dk_locus[p3dk_ix].v0[2] = sxkp[2]; 
					p3dk_locus[p3dk_ix].v1[0] = nsm[i].sv[0]; 
					p3dk_locus[p3dk_ix].v1[1] = nsm[i].sv[1]; 
					p3dk_locus[p3dk_ix].v1[2] = nsm[i].sv[2]; 
					p3dk_ix++;
#endif /* PLOT_3DKNEES */

					/* Compute the knee point const. L plane radially inwards from div[] */
					for (j = 0; j < 3; j++) {
						skp[j] = (1.0 - knpos) * sxkp[j] + knpos * div[j];
						dkp[j] = (1.0 - knpos) * dxkp[j] + knpos * div[j];
					}
#ifdef USE_GAMKNF
					/* Add this knee point only if there is some room for it */
					if ((skp[1] * skp[1] + skp[2] * skp[2]) > (5.0 * 5.0)) {	/* > 5 DE */
						for (j = 0; j < 3; j++) {
							gpnts[ngamp].p[j] = skp[j];
							gpnts[ngamp].v[j] = dkp[j];
						}
						gpnts[ngamp++].w = gmi->gamknf * gmi->gamknf;		/* Knee weight */
					}
#endif /* USE_GAMKNF */
				}
			} else {	/* Expansion needed */

				/* Compute expansion destination value */
				for (j = 0; j < 3; j++)				/* Compute compressed value */
					div[j] = gmi->gamexf * nsm[i].sdv[j] + (1.0 - gmi->gamexf) * nsm[i].sv[j];
#ifdef NEVER
				printf("Expansion:\n");
				printf("Src point = %f %f %f radius %f\n",nsm[i].sv[0], nsm[i].sv[1], nsm[i].sv[2], nsm[i].sr);
				printf("Dst point = %f %f %f radius %f\n",nsm[i].sdv[0], nsm[i].sdv[1], nsm[i].sdv[2], nsm[i].dr);
				printf("Blended dst point = %f %f %f\n",div[0], div[1], div[2]);
				printf("\n");
#endif	/* NEVER */
				/* Set the mapping point */
				for (j = 0; j < 3; j++) {
					gpnts[ngamp].p[j] = nsm[i].sv[j];
					gpnts[ngamp].v[j] = div[j];
				}
				gpnts[ngamp++].w  = 1.0;
	
				if (wrl != NULL) {
					wrl->add_vertex(wrl, nsm[i].sv);		/* Source value */
//					wrl->add_vertex(wrl, div);			/* Blended destination value */
					wrl->add_vertex(wrl, nsm[i].sdv);	/* Full compression dest. value */
				}
	
				/* Create a "knee" point */
				{
					double t;
					double sxkp[3], dxkp[3];	/* Axis intercept of knee point direction */
					double skp[3], dkp[3];		/* Inner knee point */
	
					/* Since the sv is an associated point within both gamuts, */
					/* we're going to use it as a basis for the "half way" knee point. */
					/* To stop this point affecting the grey axis mapping, we need */
					/* to adjust the half way mapping sympatheticaly to the grey axis mapping */

					/* Find the sv's corresponding point on the source L axis */
					t = (nsm[i].sv[0] - sl_cs_bp[0])/(sl_cs_wp[0] - sl_cs_bp[0]);
					for (j = 0; j < 3; j++)
						sxkp[j] = sl_cs_bp[j] + t * (sl_cs_wp[j] - sl_cs_bp[j]);
			
					dxkp[0] = sxkp[0];	/* sl is already L mapped */
			
					/* Figure destination point on destination grey axis */
					t = (dxkp[0] - dr_cs_bp[0])/(dr_cs_wp[0] - dr_cs_bp[0]);
					for (j = 1; j < 3; j++)
						dxkp[j] = dr_cs_bp[j] + t * (dr_cs_wp[j] - dr_cs_bp[j]);
					
#ifdef PLOT_3DKNEES
					p3dk_locus[p3dk_ix].v0[0] = sxkp[0]; 
					p3dk_locus[p3dk_ix].v0[1] = sxkp[1]; 
					p3dk_locus[p3dk_ix].v0[2] = sxkp[2]; 
					p3dk_locus[p3dk_ix].v1[0] = nsm[i].sv[0]; 
					p3dk_locus[p3dk_ix].v1[1] = nsm[i].sv[1]; 
					p3dk_locus[p3dk_ix].v1[2] = nsm[i].sv[2]; 
					p3dk_ix++;
#endif /* PLOT_3DKNEES */

					/* Compute the knee point const. L plane radially inwards from sv[] */
					for (j = 0; j < 3; j++) {
						skp[j] = (1.0 - knpos) * sxkp[j] + knpos * nsm[i].sv[j];
						dkp[j] = (1.0 - knpos) * dxkp[j] + knpos * nsm[i].sv[j];
					}
#ifdef USE_GAMKNF
					/* Add this knee point only if there is some room for it */
					if ((skp[1] * skp[1] + skp[2] * skp[2]) > (5.0 * 5.0)) {	/* > 5 DE */
						for (j = 0; j < 3; j++) {
							gpnts[ngamp].p[j] = skp[j];
							gpnts[ngamp].v[j] = dkp[j];
						}
						gpnts[ngamp++].w = gmi->gamknf * gmi->gamknf;		/* Knee weight */
					}
#endif /* USE_GAMKNF */
				}
			}
		}

		if (wrl != NULL) {
			double cc[3] = { 0.7, 0.7, 0.7 }; 

			wrl->make_lines(wrl, 2);
	
			wrl->make_gamut_surface(wrl, sil_gam, 0.6, cc);
			cc[0] = -1.0;
			wrl->make_gamut_surface(wrl, d_gam, 0.2, cc);

			wrl->del(wrl);
			wrl = NULL;
		}

		/* Create preliminary gamut mapping rspl, without grid boundary values. */
		/* We use this to lookup the mapping for points on the source space gamut */
		/* that result from clipping our grid boundary points */
		for (j = 0; j < 3; j++) {		/* Set resolution for all axes */
			gres[j] = mapres;
			avgdev[j] = 0.005;
		}
#ifdef USE_BOUND
		s->map = new_rspl(RSPL_NOFLAGS, 3, 3);	/* Allocate 3D -> 3D */
		s->map->fit_rspl_w(s->map, 0, gpnts, ngamp, il, ih, gres, ol, oh, smooth, avgdev, NULL);

		/* Add input range grid surface anchor points to improve clipping behaviour. */
		if (defrgrid >= 2) {
			DCOUNT(gc, 3, 3, 0, 0, defrgrid);
			double cent[3];

			sc_gam->getcent(d_gam, cent);

			DC_INIT(gc);
			for (;;) {
				/* If point is on the grid surface */
				if (   gc[0] == 0 || gc[0] == (defrgrid-1)
					|| gc[1] == 0 || gc[1] == (defrgrid-1)
					|| gc[2] == 0 || gc[2] == (defrgrid-1)) {
					double grid2gamut, gamut2cent, ww;
					co cp;

					/* Clip the point to the closest location on the source */
					/* colorspace gamut. */
					for (j = 0; j < 3; j++)
						gpnts[ngamp].p[j] = il[j] + gc[j]/(defrgrid-1.0) * (ih[j] - il[j]);
					sc_gam->nearest(sc_gam, cp.p, gpnts[ngamp].p);

					/* Then lookup the equivalent gamut mapped value */
					s->map->interp(s->map, &cp);

					for (j = 0; j < 3; j++)
						gpnts[ngamp].v[j] = cp.v[j];

					/* Compute the distance of the grid surface point to the to the */
					/* source colorspace gamut, as well as the distance from there */
					/* to the gamut center point. */
					for (grid2gamut = gamut2cent = 0.0, j = 0; j < 3; j++) {
						double tt;
						tt = gpnts[ngamp].p[j] - cp.p[j];
						grid2gamut += tt * tt;
						tt = cp.p[j] - cent[j];
						gamut2cent += tt * tt;
					}
					grid2gamut = sqrt(grid2gamut);
					gamut2cent = sqrt(gamut2cent);
					
					/* Make the weighting inversely related to distance, */
					/* to reduce influence on in gamut mapping shape, */
					/* while retaining some influence at the edge of the */
					/* grid. */
					ww = grid2gamut / gamut2cent;
					if (ww > 1.0)
						ww = 1.0;
					
					/* A low weight seems to be enough ? */
					/* the lower the better in terms of geting best hull mapping fidelity */
					gpnts[ngamp++].w = 0.1 * ww;
				}
				DC_INC(gc);
				if (DC_DONE(gc))
					break;
			}
		}
#endif /* USE_BOUND */

		/* --------------------------- */
		/* Compute the output bounding values, and check input range hasn't changed */
		for (i = 0; i < ngamp; i++) {
			for (j = 0; j < 3; j++) {
				if (gpnts[i].p[j] < (il[j]-1e-5) || gpnts[i].p[j] > (ih[j]+1e-5))
					warning("gammap internal: input bounds has changed! %f <> %f <> %f",il[j],gpnts[i].p[j],ih[j]);
				if (gpnts[i].v[j] < ol[j])
					ol[j] = gpnts[i].v[j];
				if (gpnts[i].v[j] > oh[j])
					oh[j] = gpnts[i].v[j];
			}
		}

		/* --------------------------- */

#ifdef NEVER		/* Dump out all the mapping points */
		{
			for (i = 0; i < ngamp; i++) {
				printf("%d: %f %f %f -> %f %f %f\n",i,
					gpnts[i].p[0], gpnts[i].p[1], gpnts[i].p[2],
					gpnts[i].v[0], gpnts[i].v[1], gpnts[i].v[2]);
			}
		}
#endif

		/* Create the final gamut mapping rspl. */
		if (s->map != NULL)
			s->map->del(s->map);
		s->map = new_rspl(RSPL_NOFLAGS, 3, 3);	/* Allocate 3D -> 3D */
		if (s->map->fit_rspl_w(s->map, 0, gpnts, ngamp, il, ih, gres, ol, oh, smooth, avgdev, NULL)) {
			if (verb)
				fprintf(stderr,"Warning: Gamut mapping is non-monotonic - may not be very smooth !\n");
		}

		/* return the min and max of the input values valid in the grid */
		s->map->get_in_range(s->map, s->imin, s->imax); 

#ifdef CHECK_NEARMAP
		/* Check how accurate gamut shell mapping is against nsm */
		if (verb) {
			double de, avgde = 0.0, maxde = 0.0;		/* DE stats */
			for (i = 0; i < nnsm; i++) {
				double div[3];			/* Ref destination nearmap value */
				co cp;

				if (nsm[i].sr >= nsm[i].dr) {		/* Compression needed */

					for (j = 0; j < 3; j++) {			/* Compute compressed value */
						cp.p[j] = nsm[i].sv[j];
						div[j] = gmi->gamcpf * nsm[i].sdv[j] + (1.0 - gmi->gamcpf) * nsm[i].sv[j];
					}

				} else {	/* Expansion needed */

					/* Compute expansion destination value */
					for (j = 0; j < 3; j++) {
						cp.p[j] = nsm[i].sv[j];
						div[j] = gmi->gamexf * nsm[i].sdv[j] + (1.0 - gmi->gamexf) * nsm[i].sv[j];
					}
				}
				s->map->interp(s->map, &cp);
				
				de = icmLabDE(div, cp.v);
				avgde += de;
				if (de > maxde)
					maxde = de;
			}
			avgde /= nnsm;
			fprintf(stderr,"Gamut hull mapping errors: = avg %f, max %f\n",avgde,maxde);
		}
#endif /* CHECK_NEARMAP */

		free(nsm);

		/* If requested, enhance the saturation of the output values. */
		if (gmi->satenh > 0.0) {
			adjustsat cx;		/* Adjustment context */

			/* Compute what our source white and black points actually maps to */
			s->domap(s, cx.wp, s_mt_wp);
			s->domap(s, cx.bp, s_mt_bp);

			cx.dst = d_gam; 
			cx.satenh = gmi->satenh; 

			/* Saturation enhance the output values */
			s->map->re_set_rspl(
				s->map,				/* this */
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

			/* Check what the source white and black points actually maps to */
			s->domap(s, a_wp, s_mt_wp);
			s->domap(s, a_bp, s_mt_bp);

#ifdef VERBOSE
			if (verb) {
				printf("White is %f %f %f, should be %f %f %f\n",
				a_wp[0], a_wp[1], a_wp[2], d_mt_wp[0], d_mt_wp[1], d_mt_wp[2]);
				if (bph == gmm_bendBP)
					printf("Black is %f %f %f, ideal is %f %f %f\n",
					a_bp[0], a_bp[1], a_bp[2], d_mt_bp[0], d_mt_bp[1], d_mt_bp[2]);
				else
					printf("Black is %f %f %f, should be %f %f %f\n",
					a_bp[0], a_bp[1], a_bp[2], d_mt_bp[0], d_mt_bp[1], d_mt_bp[2]);
			}
#endif /* VERBOSE */

			/* Setup the fine tune transform */

			/* We've decided not to fine tune the black point if we're */
			/* bending to the destination black, as the bend is not */
			/* followed perfectly (too sharp, or in conflict with */
			/* the surface mapping ?) and we don't want to shift */
			/* mid neutrals due to this. */

			/* Compute rotation/scale relative white point matrix */
			if (bph == gmm_bendBP)
				icmVecRotMat(cx.mat, a_wp, a_bp, d_mt_wp, a_bp);		/* just wp */
			else
				icmVecRotMat(cx.mat, a_wp, a_bp, d_mt_wp, d_mt_bp);		/* wp & bp */

			/* Fine tune the 3D->3D mapping */
			s->map->re_set_rspl(
				s->map,		/* this */
				0,					/* Combination of flags */
				(void *)&cx,		/* Opaque function context */
				adjust_wb_func /* Function to set from */
			);

#ifdef VERBOSE
			if (verb) {
				/* Check what the source white and black points actually maps to */
				s->domap(s, a_wp, s_mt_wp);
				s->domap(s, a_bp, s_mt_bp);
	
				printf("After fine tuning:\n");
				printf("White is %f %f %f, should be %f %f %f\n",
				a_wp[0], a_wp[1], a_wp[2], d_mt_wp[0], d_mt_wp[1], d_mt_wp[2]);
				if (bph == gmm_bendBP)
					printf("Black is %f %f %f, ideal is %f %f %f\n",
					a_bp[0], a_bp[1], a_bp[2], d_mt_bp[0], d_mt_bp[1], d_mt_bp[2]);
				else
					printf("Black is %f %f %f, should be %f %f %f\n",
					a_bp[0], a_bp[1], a_bp[2], d_mt_bp[0], d_mt_bp[1], d_mt_bp[2]);
			}
#endif /* VERBOSE */
		}

#ifdef PLOT_GAMUTS
		scl_gam->write_vrml(scl_gam, "src.wrl", 1, 0);
		sil_gam->write_vrml(sil_gam, "img.wrl", 1, 0);
		d_gam->write_vrml(d_gam, "dst.wrl", 1, 0);
		sc_gam->write_trans_vrml(sc_gam, "gmsrc.wrl", 1, 0, map_trans, s);
#endif

#ifdef PLOT_3DKNEES
		/* Plot one graph per 3D gamut boundary mapping point */
		for (j = 0; j < p3dk_ix; j++) {
			double xx[XRES];
			double yy[XRES];

			printf("Vector %f %f %f -> %f %f %f\n", p3dk_locus[j].v0[0], p3dk_locus[j].v0[1], p3dk_locus[j].v0[2], p3dk_locus[j].v1[0], p3dk_locus[j].v1[1], p3dk_locus[j].v1[2]);

			for (i = 0; i < XRES; i++) {
				double v;
				co cp;		/* Conversion point */
				v = (i/(double)(XRES-1.0));
				cp.p[0] = p3dk_locus[j].v0[0] + v * (p3dk_locus[j].v1[0] - p3dk_locus[j].v0[0]);
				cp.p[1] = p3dk_locus[j].v0[1] + v * (p3dk_locus[j].v1[1] - p3dk_locus[j].v0[1]);
				cp.p[2] = p3dk_locus[j].v0[2] + v * (p3dk_locus[j].v1[2] - p3dk_locus[j].v0[2]);
				xx[i] = sqrt(cp.p[1] * cp.p[1] + cp.p[2] * cp.p[2]);
				s->map->interp(s->map, &cp);
				yy[i] = sqrt(cp.v[1] * cp.v[1] + cp.v[2] * cp.v[2]);
			}
			do_plot(xx,yy,NULL,NULL,XRES);
		}
		free(p3dk_locus);
#endif /* PLOT_3DKNEES */
	}

#ifdef NEVER
	/* Test some ProPhoto trouble values */
	{
		int i;
		double vals[2][3] = {
			{ 37.1, -100.3, -58.0 },	/* Comes out too white */
			{ 36.9, -97.6,  -59.6 }		/* Comes out as blue */
		};
		for (i = 0; i < 2; i++) {
			double mm[3];
			s->domap(s, mm, vals[i]);
			printf("Maps %d: %f %f %f -> %f %f %f\n",i,
			   vals[i][0], vals[i][1], vals[i][2], mm[0], mm[1], mm[2]);
		}
	}
#endif

	free(gpnts);
	if (sil_gam != scl_gam)
		sil_gam->del(sil_gam);
	scl_gam->del(scl_gam);

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
	s->grey->del(s->grey);
	s->igrey->del(s->igrey);
	s->map->del(s->map);

	free(s);
}

/* Apply the gamut mapping to the given color value */
static void domap(
gammap *s,
double *out,
double *in
) {
	double rin[3];
	co cp;

	if (s->dbg) printf("domap: got input %f %f %f\n",in[0],in[1],in[2]);
	icmMul3By3x4(rin, s->grot, in);		/* Rotate */

	if (s->dbg) printf("domap: after rotate %f %f %f\n",rin[0],rin[1],rin[2]);
	cp.p[0] = rin[0];
	s->grey->interp(s->grey, &cp);		/* L map */

	if (s->dbg) printf("domap: after L map %f %f %f\n",cp.v[0],rin[1],rin[2]);

	/* Clip out of range a, b proportionately */
	if (rin[1] < s->imin[1] || rin[1] > s->imax[1]
	 || rin[2] < s->imin[2] || rin[2] > s->imax[2]) {
		double as = 1.0, bs = 1.0;
		if (rin[1] < s->imin[1])
			as = s->imin[1]/rin[1];
		else if (rin[1] > s->imax[1])
			as = s->imax[1]/rin[1];
		if (rin[2] < s->imin[2])
			bs = s->imin[2]/rin[2];
		else if (rin[2] > s->imax[2])
			bs = s->imax[2]/rin[2];
		if (bs < as)
			as = bs;
		rin[1] *= as;
		rin[2] *= as;
	}

	cp.p[0] = cp.v[0];					/* 3D map */
	cp.p[1] = rin[1];
	cp.p[2] = rin[2];
	s->map->interp(s->map, &cp);

	out[0] = cp.v[0];
	out[1] = cp.v[1];
	out[2] = cp.v[2];
	if (s->dbg) printf("domap: after 3D map %f %f %f\n\n",out[0],out[1],out[2]);
}

/* Function to pass to rspl to invert grey curve */
static void inv_grey_func(
	void *cntx,
	double *out,
	double *in
) {
	rspl *fwd = (rspl *)cntx;
	int nsoln;		/* Number of solutions found */
	co pp[2];		/* Room for all the solutions found */

	pp[0].p[0] = 
	pp[0].v[0] = in[0];

	nsoln = fwd->rev_interp(
		fwd,
		RSPL_NEARCLIP,		/* Clip to nearest (faster than vector) */
		2,					/* Maximum number of solutions allowed for */
		NULL, 				/* No auxiliary input targets */
		NULL,				/* Clip vector direction and length */
		pp);				/* Input and output values */

	nsoln &= RSPL_NOSOLNS;		/* Get number of solutions */

	if (nsoln != 1)
		error("gammap: Unexpected failure to find reverse solution for grey axis lookup");

	out[0] = pp[0].p[0];
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
	adjustwb *p = (adjustwb *)pp;

	/* Do a linear mapping from swp -> dwp and sbp -> dbp, */
	/* to compute the adjusted value. */
	icmMul3By3x4(out, p->mat, out);
}





















