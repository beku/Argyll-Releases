
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
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * For a discussion of gamut mapping strategy used,
 * see gammap.txt
 */

/*
 * TTBD:
 *       Improve error handling.
 */

#define VERBOSE			/* (defined) Print out extra interesting information when verbose is set */
#define XRES 100		/* Res of plot */

	/* What do display when user requests disgnostic .wrl */
#define PLOT_SRC_GMT	/* [Def] Plot the source surface to "gammap.wrl" as well */
#define PLOT_DST_GMT	/* [Def] Plot the dest surface to "gammap.wrl" as well */
#undef PLOT_SRC_CUSPS	/* [Und] Plot the source surface cusps to "gammap.wrl" as well */
#undef PLOT_DST_CUSPS	/* [Und] Plot the dest surface cusps to "gammap.wrl" as well */
#undef PLOT_TRANSSRC_CUSPS	/* [Und] Plot the gamut mapped source surface cusps to "gammap.wrl" */
#define PLOT_AXES		/* [Und] Plot the axes to "gammap.wrl" as well */
#undef SHOW_VECTOR_INDEXES	/* [Und] Show the mapping vector index numbers */
#define SHOW_MAP_VECTORS	/* [Def] Show the mapping vectors */
#undef SHOW_SUB_SURF	/* [Und] Show the sub-surface mapping vector */
#undef SHOW_CUSPMAP		/* [Und] Show the cusp mapped vectors rather than final vectors */
#define SHOW_ACTUAL_VECTORS		/* [Def] Show how the source vectors actually map */
#undef SHOW_ACTUAL_VEC_DIFF		/* [Und] Show how the difference between guide and actual vectors */

#undef PLOT_LMAP		/* [undef] Plot L map */
#undef PLOT_GAMUTS		/* Save (part mapped) input and output gamuts as */
						/* src.wrl, img.wrl, dst.wrl, gmsrc.wrl */
#undef PLOT_3DKNEES		/* Plot the 3D compression knees */
#define CHECK_NEARMAP	/* Check how accurately near map vectors are represented by rspl */

#define RSPLFLAGS (0)	/* Default rspl flags */
#define MAINRSPLFLAGS (0 /* | RSPL_EXTRAFIT2 */ )		/* Flags for main mapping rspl */

#define USE_GLUMKNF			/* [Define] Enable luminence knee function points */
#define USE_GAMKNF			/* [Define] Enable 3D knee function points */
#define USE_BOUND			/* [Define] Enable grid boundary anchor points */

#undef SHOW_NEIGBORS		/* Show nearsmth neigbors in gammap.wrl */

/* Optional marker points for gamut mapping diagnosotic */
struct {
	int type;		/* 1 = src point (xlate), 2 = dst point (no xlate) */
					/* 0 = end marker */
	double pos[3];	/* Position, (usually in Jab space) */
	double col[3];	/* RGB color */
} markers[] = {
	{ 0, },								/* End marker */
	{ 1, { 67.575411, -37.555250, -36.612862 }, { 1.0, 0.3, 0.3 } },	/* bad source in red (Red) */
	{ 1, { 61.003078, -44.466554, 1.922585 }, { 0.0, 1.0, 0.3 } },	/* good source in green */
	{ 2, { 49.294793, 50.749543, -51.383167 }, { 1.0, 0.0, 0.0 } },	
	{ 2, { 42.783425, 49.089363, -37.823712 }, { 0.0, 1.0, 0.0 } },	
	{ 2, { 41.222695, 63.911823, 37.695310 }, { 0.0, 1.0, 0.3 } },	/* destination in green */
	{ 1, { 41.951770, 60.220284, 34.788195 }, { 1.0, 0.3, 0.3 } },	/* source in red (Red) */
	{ 2, { 41.222695, 63.911823, 37.695310 }, { 0.3, 1.3, 0.3 } },		/* Dest in green */
	{ 1, { 85.117353, -60.807580, -22.195118 }, { 0.3, 0.3, 1 } },		/* Cyan Source (Blue) */
	{ 2, { 61.661622, -38.164411, -18.090824 }, { 1.0, 0.3, 0.3 } },	/* CMYK destination (Red) */
	{ 0 }								/* End marker */
};

/* Optional marker rings for gamut mapping diagnosotic */
struct {
	int type;		/* 1 = src ring point, 2 = ignore, */
					/* 0 = end marker */
	double ppoint[3];	/* Location of a point on the plane in source space */
	double pnorm[3];	/* Plane normal direction in source space */
	int    nverts;		/* Number of points to make ring */
	double rad;			/* Relative Radius from neutral to source surface (0.0 - 1.0) */
	double scol[3];	/* Source RGB color */
	double dcol[3];	/* Destination RGB color */
} rings[] = {
	{ 0 },								/* End marker */
	{ 1,
		{ 60.0, 0.0, 0.0 }, { 1.0, 0.8, 0.0 },		/* plane point and normal */
		100, 1.0,									/* 20 vertexes at source radius */
		{ 0.0, 1.0, 0.0 },			/* Green source */
		{ 1.0, 0.0, 0.0 }			/* Red destination */
	},
	{ 1,
		{ 60.0, 0.0, 0.0 }, { 1.0, 0.8, 0.0 },		/* plane point and normal */
		100, 0.9,									/* 20 vertexes at source radius */
		{ 0.0, 1.0, 0.0 },			/* Green source */
		{ 1.0, 0.0, 0.0 }			/* Red destination */
	},
	{ 1,
		{ 60.0, 0.0, 0.0 }, { 1.0, 0.8, 0.0 },		/* plane point and normal */
		100, 0.8,									/* 20 vertexes at source radius */
		{ 0.0, 1.0, 0.0 },			/* Green source */
		{ 1.0, 0.0, 0.0 }			/* Red destination */
	},
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
#include "numlib.h"
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

#define XVRA 3.0	/* Extra mapping vertex ratio over tri verts from gamut */

/* The smoothed near weighting control values. */
/* These weightings setup the detailed behaviour of the */
/* gamut mapping for the fully perceptual and saturation intents. */
/* They are ordered here by increasing priority. A -ve value is ignored */

/* Perceptual mapping weights, where smoothness and proportionality are important.. */
gammapweights pweights[] = {
	{
		gmm_default,	/* Non hue specific defaults */
		{			/* Weighting of absolute error of destination from source */
			1.0,	/* Absolute error overall weight */
			{
				40.0,	/* Absolute luminance error weight */
				25.0,	/* Absolute chroma error weight */
				80.0	/* Absolute hue error weight */
			}
		},
		{			/* Weighting of relative error of destination points to each */
					/* other, compared to source points to each other. */
			1.0,	/* Relative error overall weight */
			{
				30.0,	/* Relative luminance error weight */
				20.0,	/* Relative chroma error weight */
				30.0	/* Relative hue error weight */
			},
			25.0, 30.0	/* Relative Smoothing radius L* H* */
		},
		{			/* Weighting of error between destination point and source */
					/* point radially mapped towards center of destination. */
			
			0.1,	/* Radial error overall weight */
			{
				30.0,	/* Radial luminance error weight */
				10.0,	/* Radial chroma error weight */
				20.0	/* Radial hue error weight */
			}
		},
		{		/* Weighting of excessive compression error, which is */
				/* the src->dst vector length over the available dst depth. */
				/* The depth is half the distance to the intersection of the */
				/* vector to the other side of the gamut. (doesn't get triggered much ?) */
			100.0,		/* Compression depth weight */
			100.0		/* Expansion depth weight */
		},
	
		{
			{
				0.1,	/* Cusp luminance alignment weighting 0 = none, 1 = full */
				0.0,	/* Cusp chroma alignment weighting    0 = none, 1 = full */
				0.1		/* Cusp hue alignment weighting       0 = none, 1 = full */
			},
			1.03		/* Chroma expansion 1 = none */
		}
	},
	{
		gmm_light_yellow,		/* Treat yellow differently, to get purer result. */
		{			/* Weighting of absolute error of destination from source */
			-1.0,	/* Absolute error overall weight */
			{
				30.0,	/* Absolute luminance error weight */
				50.0,	/* Absolute chroma error weight */
				25.0	/* Absolute hue error weight */
			}
		},
		{			/* Weighting of relative error of destination points to each */
					/* other, compared to source points to each other. */
			-1.0,	/* Relative error overall weight */
			{
				20.0,	/* Relative luminance error weight */
				5.0,	/* Relative chroma error weight */
				20.0		/* Relative hue error weight */
			},
			15.0, 25.0	/* Relative smoothing radius */
		},
		{			/* Weighting of error between destination point and source */
					/* point radially mapped towards center of destination. */
			-1.0,	/* Radial error overall weight */
			{
				-1.0,	/* Radial luminance error weight */
				-1.0,	/* Radial chroma error weight */
				-1.0		/* Radial hue error weight */
			}
		},
		{		/* Weighting of excessive compression error, which is */
				/* the src->dst vector length over the available dst depth. */
				/* The depth is half the distance to the intersection of the */
				/* vector to the other side of the gamut. (doesn't get triggered much ?) */
			-1.0,		/* Compression depth weight */
			-1.0		/* Expansion depth weight */
		},
	
		{
			{
				0.8,	/* Cusp luminance alignment weighting 0 = none, 1 = full */
				0.5,	/* Cusp chroma alignment weighting    0 = none, 1 = full */
				0.8		/* Cusp hue alignment weighting       0 = none, 1 = full */
			},
			1.15		/* Chroma expansion 1 = none */
		}
	},
	{
		gmm_dark_colors,		/* Make dark colors have more constant L */
		{			/* Weighting of absolute error of destination from source */
			-1.0,	/* Absolute error overall weight */
			{
				100.0,	/* Absolute luminance error weight */
				20.0,	/* Absolute chroma error weight */
				80.0	/* Absolute hue error weight */
			}
		},
		{			/* Weighting of relative error of destination points to each */
					/* other, compared to source points to each other. */
			-1.0,	/* Relative error overall weight */
			{
				-1.0,	/* Relative luminance error weight */
				-1.0,	/* Relative chroma error weight */
				-1.0		/* Relative hue error weight */
			},
			10.0, 30.0	/* Relative Smoothing radius L* H* */
		},
		{			/* Weighting of error between destination point and source */
					/* point radially mapped towards center of destination. */
			-1.0,	/* Radial error overall weight */
			{
				-1.0,	/* Radial luminance error weight */
				-1.0,	/* Radial chroma error weight */
				-1.0		/* Radial hue error weight */
			}
		},
		{		/* Weighting of excessive compression error, which is */
				/* the src->dst vector length over the available dst depth. */
				/* The depth is half the distance to the intersection of the */
				/* vector to the other side of the gamut. (doesn't get triggered much ?) */
			-1.0,		/* Compression depth weight */
			-1.0		/* Expansion depth weight */
		},
	
		{
			{
				-1.0,	/* Cusp luminance alignment weighting 0 = none, 1 = full */
				-1.0,	/* Cusp chroma alignment weighting    0 = none, 1 = full */
				-1.0		/* Cusp hue alignment weighting       0 = none, 1 = full */
			},
			-1.0		/* Chroma expansion 1 = none */
		}
	},
	{
		gmm_end,
	}
};
double pm21fsm = 0.0;		/* Level of inverse RSPL smoothing for perceptual, 0 = none */
double psmooth = 7.0;		/* Level of RSPL smoothing for perceptual, 1 = nominal */

/* Saturation mapping weights, where saturation has priority over smoothness */
gammapweights sweights[] = {
	{
		gmm_default,	/* Non hue specific defaults */
		{			/* Weighting of absolute error of destination from source */
			1.0,	/* Absolute error overall weight */
			{
				45.0,	/* Absolute luminance error weight */
				25.0,	/* Absolute chroma error weight */
				60.0	/* Absolute hue error weight */
			}
		},
		{			/* Weighting of relative error of destination points to each */
					/* other, compared to source points to each other. */
			1.0,	/* Relative error overall weight */
			{
				30.0,	/* Relative luminance error weight */
				20.0,	/* Relative chroma error weight */
				30.0	/* Relative hue error weight */
			},
			25.0, 25.0	/* Relative smoothing radius, L* H* */
		},
		{			/* Weighting of error between destination point and source */
					/* point radially mapped towards center of destination. */
			0.1,	/* Radial error overall weight */
			{
				10.0,	/* Radial luminance error weight */
				10.0,	/* Radial chroma error weight */
				20.0	/* Radial hue error weight */
			}
		},
		{		/* Weighting of excessive compression error, which is */
				/* the src->dst vector length over the available dst depth. */
				/* The depth is half the distance to the intersection of the */
				/* vector to the other side of the gamut. (doesn't get triggered much ?) */
			100.0,		/* Compression depth weight */
			100.0		/* Expansion depth weight */
		},
	
		{
			{
				0.6,	/* Cusp luminance alignment weighting 0 = none, 1 = full */
				0.5,	/* Cusp chroma alignment weighting    0 = none, 1 = full */
				0.6		/* Cusp hue alignment weighting       0 = none, 1 = full */
			},
			1.05		/* Chroma expansion 1 = none */
		},
	},
	{
		gmm_light_yellow,		/* Treat yellow differently, to get purer result. */
		{			/* Weighting of absolute error of destination from source */
			-1.0,	/* Absolute error overall weight */
			{
				-1.0,	/* Absolute luminance error weight */
				50.0,	/* Absolute chroma error weight */
				-1.0	/* Absolute hue error weight */
			}
		},
		{			/* Weighting of relative error of destination points to each */
					/* other, compared to source points to each other. */
			-1.0,	/* Relative error overall weight */
			{
				-1.0,	/* Relative luminance error weight */
				10.0,	/* Relative chroma error weight */
				-1.0	/* Relative hue error weight */
			},
			10.0, 20.0	/* Relative smoothing radius */
		},
		{			/* Weighting of error between destination point and source */
					/* point radially mapped towards center of destination. */
			-1.0,	/* Radial error overall weight */
			{
				-1.0,	/* Radial luminance error weight */
				-1.0,	/* Radial chroma error weight */
				-1.0	/* Radial hue error weight */
			}
		},
		{		/* Weighting of excessive compression error, which is */
				/* the src->dst vector length over the available dst depth. */
				/* The depth is half the distance to the intersection of the */
				/* vector to the other side of the gamut. (doesn't get triggered much ?) */
			-1.0,		/* Compression depth weight */
			-1.0		/* Expansion depth weight */
		},
		{
			{
				1.0,	/* Cusp luminance alignment weighting 0 = none, 1 = full */
				1.0,	/* Cusp chroma alignment weighting    0 = none, 1 = full */
				1.0		/* Cusp hue alignment weighting       0 = none, 1 = full */
			},
			1.15		/* Chroma expansion 1 = none */
		}
	},
	{
		gmm_dark_colors,		/* Make dark colors have more constant L */
		{			/* Weighting of absolute error of destination from source */
			-1.0,	/* Absolute error overall weight */
			{
				60.0,	/* Absolute luminance error weight */
				25.0,	/* Absolute chroma error weight */
				60.0	/* Absolute hue error weight */
			}
		},
		{			/* Weighting of relative error of destination points to each */
					/* other, compared to source points to each other. */
			-1.0,	/* Relative error overall weight */
			{
				-1.0,	/* Relative luminance error weight */
				-1.0,	/* Relative chroma error weight */
				-1.0		/* Relative hue error weight */
			},
			10.0, 20.0	/* Relative Smoothing radius L* H* */
		},
		{			/* Weighting of error between destination point and source */
					/* point radially mapped towards center of destination. */
			-1.0,	/* Radial error overall weight */
			{
				-1.0,	/* Radial luminance error weight */
				-1.0,	/* Radial chroma error weight */
				-1.0		/* Radial hue error weight */
			}
		},
		{		/* Weighting of excessive compression error, which is */
				/* the src->dst vector length over the available dst depth. */
				/* The depth is half the distance to the intersection of the */
				/* vector to the other side of the gamut. (doesn't get triggered much ?) */
			-1.0,		/* Compression depth weight */
			-1.0		/* Expansion depth weight */
		},
	
		{
			{
				-1.0,	/* Cusp luminance alignment weighting 0 = none, 1 = full */
				-1.0,	/* Cusp chroma alignment weighting    0 = none, 1 = full */
				-1.0		/* Cusp hue alignment weighting       0 = none, 1 = full */
			},
			-1.0		/* Chroma expansion 1 = none */
		}
	},
	{
		gmm_end
	}
};
double sm21fsm = 0.0;		/* Level of inverse RSPL smoothing for perceptual, 0 = none */
double ssmooth = 7.0;		/* Level of RSPL smoothing for saturation */

/*
 * Notes:
 *       The "knee" shape produced by the rspl (regular spline) code
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
	int src_kbp,		/* Use K only black point as src gamut black point */
	int dst_kbp,		/* Use K only black point as dst gamut black point */
	int dst_cmymap,		/* masks C = 1, M = 2, Y = 4 to force 100% cusp map */
	int rel_oride,		/* 0 = normal, 1 = clip like, 2 = max relative */
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
	gamut *scl_gam;		/* Source colorspace gamut with rotation and L mapping applied */
	gamut *sil_gam;		/* Source image gamut with rotation and L mapping applied */

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
	double xvra = XVRA;	/* Extra ss vertex ratio to src gamut vertex count */
	int j;

#if defined(PLOT_LMAP) || defined(PLOT_GAMUTS) || defined(PLOT_3DKNEES)
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
	if (src_kbp) {
		// ~~99 Hmm. Shouldn't this be colorspace rather than gamut ????
		if (sc_gam->getwb(sc_gam, NULL, NULL, NULL, s_cs_wp, NULL, s_cs_bp)) {
//		if (sc_gam->getwb(sc_gam, s_cs_wp, NULL, s_cs_bp, NULL, NULL, NULL))
			fprintf(stderr,"gamut map: Unable to read source colorspace white and black points\n");
			free(s);
			return NULL;
		}
	} else {
		if (sc_gam->getwb(sc_gam, NULL, NULL, NULL, s_cs_wp, s_cs_bp, NULL)) {
//		if (sc_gam->getwb(sc_gam, s_cs_wp, s_cs_bp, NULL, NULL, NULL, NULL))
			fprintf(stderr,"gamut map: Unable to read source colorspace white and black points\n");
			free(s);
			return NULL;
		}
	}

	/* If source space is source gamut */
	if (si_gam == NULL) {
		si_gam = sc_gam;
		for (j = 0; j < 3; j++) {
			s_ga_wp[j] = s_cs_wp[j];
			s_ga_bp[j] = s_cs_bp[j];
		}

	/* Else have explicit sourcegamut */
	} else {

		if (src_kbp) {
			if (si_gam->getwb(si_gam, NULL, NULL, NULL, s_ga_wp, NULL, s_ga_bp)) {
				fprintf(stderr,"gamut map: Unable to read source gamut white and black points\n");
				free(s);
				return NULL;
			}
		} else {
			if (si_gam->getwb(si_gam, NULL, NULL, NULL, s_ga_wp, s_ga_bp, NULL)) {
				fprintf(stderr,"gamut map: Unable to read source gamut white and black points\n");
				free(s);
				return NULL;
			}
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

	if (dst_kbp) {
		if (d_gam->getwb(d_gam, NULL, NULL, NULL, d_cs_wp, NULL, d_cs_bp)) {
			fprintf(stderr,"gamut map: Unable to read destination white and black points\n");
			free(s);
			return NULL;
		}
	} else {
		if (d_gam->getwb(d_gam, NULL, NULL, NULL, d_cs_wp, d_cs_bp, NULL)) {
			fprintf(stderr,"gamut map: Unable to read destination white and black points\n");
			free(s);
			return NULL;
		}
	}

#ifdef VERBOSE
	if (verb) {
		if (src_kbp)
			printf("Using Src K only black point\n");

		if (dst_kbp)
			printf("Using Dst K only black point\n");

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
		double hawp[3], habp[3];	/* Half (full white, not black) adapted destination w & b */

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
		if (d_gam->vector_isect(d_gam, fabp, fawp, fabp, fawp, NULL, NULL, NULL, NULL) == 0)
			error("gamut: vector_isect failed!");

		/* To work around the problem that vector_isect() is not entirely accurate, */
		/* special case the situation where gmi->greymf == 1.0 */
		if (gmi->greymf > 0.99) {
			for (j = 0; j < 3; j++) {
				fawp[j] = d_cs_wp[j];
				fabp[j] = d_cs_bp[j];
			}
		}

		/* If dst_kbp is set, then clipping to the dest gamut doesn't do what we want, */
		/* since it extends the black to a full composite black point. */
		/* A "K only" gamut is hard to define, so do a hack: */
		/* scale fabp[] towards fawp[] so that it has the same L as */
		/* the destination K only black point. */
		if (dst_kbp && fabp[0] < d_cs_bp[0]) {
			t = (d_cs_bp[0] - fawp[0])/(fabp[0] - fawp[0]);

			for (j = 0; j < 3; j++)
				fabp[j] = fawp[j] + t * (fabp[j] - fawp[j]);
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
		if (d_gam->vector_isect(d_gam, habp, hawp, habp, hawp, NULL, NULL, NULL, NULL) == 0)
			error("gamut: vector_isect failed!");

		/* To work around the problem that vector_isect() is not entirely accurate, */
		/* special case the situation where gmi->greymf == 1.0 */
		if (gmi->greymf > 0.99) {
			for (j = 0; j < 3; j++) {
				hawp[j] = d_cs_wp[j];
			}
		}

		/* If dst_kbp is set, then clipping to the dest gamut doesn't do what we want, */
		/* since it extends the black to a full composite black point. */
		/* A "K only" gamut is hard to define, so do a hack: */
		/* scale habp[] towards hawp[] so that it has the same L as */
		/* the destination K only black point. */
		if (dst_kbp && habp[0] < d_cs_bp[0]) {
			t = (d_cs_bp[0] - hawp[0])/(habp[0] - hawp[0]);

			for (j = 0; j < 3; j++)
				habp[j] = hawp[j] + t * (habp[j] - hawp[j]);
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

			/* Use the half adapted (full white, not black) white and black points */
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
		/* be applied first to all source points. */
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
			double kppos = 0.30;		/* Knee point ratio between white/black & center */
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
		if ((s->grey = new_rspl(RSPL_NOFLAGS, 1, 1)) == NULL)	/* Allocate 1D -> 1D */
			error("gamut: grey new_rspl failed");
	
		il[0] = -1.0;		/* Set possible input range */
		ih[0] = 101.0;
		ol[0] = 0.0;		/* Set normalisation output range */
		oh[0] = 100.0;

#ifdef NEVER		/* Dump out the L mapping points */
		{
			int i;
			printf("1D rspl L mapping points:\n");
			for (i = 0; i < ngreyp; i++)
				printf("%d %f -> %f (w %f)\n",i,lpnts[i].p[0],lpnts[i].v[0],lpnts[i].w);
		}
#endif
		/* Create spline from the data points, with appropriate smoothness. */
		avgdev[0] = 0.005;
		if (s->grey->fit_rspl_w(s->grey, RSPLFLAGS, lpnts, ngreyp, il, ih, &gres, ol, oh, 5.0, avgdev, NULL)) {
			fprintf(stderr,"Warning: Grey axis mapping is non-monotonic - may not be very smooth ?\n");
		}

		/* Create an inverse mapping too, for reverse gamut and/or expansion. */
		il[0] = -1.0;		/* Set possible input range */
		ih[0] = 101.0;
		ol[0] = 0.0;		/* Set normalisation output range */
		oh[0] = 100.0;

		if ((s->igrey = new_rspl(RSPL_NOFLAGS, 1, 1)) == NULL)	/* Allocate 1D -> 1D */
			error("gamut: igrey new_rspl failed");
			
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
		double wp[3], bp[3], kp[3];
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

		scl_gam = new_gamut(sc_gam->getsres(sc_gam), sc_gam->getisjab(sc_gam), sc_gam->getisrast(sc_gam));
		scl_gam->setnofilt(scl_gam);

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
		if (sc_gam->getwb(sc_gam, wp, bp, kp, NULL, NULL, NULL) == 0) {
			icmMul3By3x4(wp, s->grot, wp);
			cp.p[0] = wp[0];			/* L value */
			s->grey->interp(s->grey, &cp);
			wp[0] = cp.v[0];

			icmMul3By3x4(bp, s->grot, bp);
			cp.p[0] = bp[0];			/* L value */
			s->grey->interp(s->grey, &cp);
			bp[0] = cp.v[0];

			icmMul3By3x4(kp, s->grot, kp);
			cp.p[0] = kp[0];			/* L value */
			s->grey->interp(s->grey, &cp);
			kp[0] = cp.v[0];

			scl_gam->setwb(scl_gam, wp, bp, kp);
		}

		if (sc_gam == si_gam)
			sil_gam = scl_gam;

		else {
			sil_gam = new_gamut(si_gam->getsres(si_gam), si_gam->getisjab(si_gam), si_gam->getisrast(si_gam));
			sil_gam->setnofilt(sil_gam);

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

	/* Create all the 3D->3D gamut mapping points and 3D rspl, */
	/* if there is any compression or expansion to do. */
	if (gmi->gamcpf > 1e-6 || gmi->gamexf > 1e-6) {
		int nspts;		/* Number of source gamut surface points */
		int rgridpts;	/* Number of range surface grid points */
		int i, j;
		datai il, ih;
		datao ol, oh;
		int gres[MXDI];
		double avgdev[MXDO];
		nearsmth *nsm = NULL;	/* Returned list of near smooth points */
		int nnsm;				/* Number of near smoothed points */
		double brad = 0.0;		/* Black bend radius */
		gammapweights xpweights[14], xsweights[14];	/* Explicit perceptial and sat. weights */
		gammapweights xwh[14]; 	/* Structure holding blended weights */
		double m21fsm = 1.0;	/* Level of inverse RSPL smoothing, blend of pm21fsm and sm21fsm */
		double smooth = 1.0;	/* Level of 3D RSPL smoothing, blend of psmooth and ssmooth */
		vrml *wrl = NULL;		/* Gamut mapping illustration (hulls + guide vectors) */
		cgats *locus = NULL;	/* Diagnostic locus to plot in wrl, NULL if none */

#ifdef PLOT_3DKNEES
typedef struct {
	double v0[3], v1[3];
} p3dk_lpoint;
		p3dk_lpoint *p3dk_locus;
		int p3dk_ix = 0;
#endif /* PLOT_3DKNEES */

		if ((gmi->gampwf + gmi->gamswf) > 0.1)
			m21fsm = (gmi->gampwf * pm21fsm) + (gmi->gamswf * sm21fsm);

		/* Get the maximum number of points that will be created */
		nspts = near_smooth_np(scl_gam, sil_gam, d_gam, xvra);
		
		rgridpts = 0;
#ifdef USE_BOUND
		if (defrgrid >= 2) {
			rgridpts = defrgrid * defrgrid * defrgrid
			         - (defrgrid -2) * (defrgrid -2) * (defrgrid -2);
		}
#endif

		if ((gpnts = (cow *)malloc((nres + 3 * nspts + rgridpts) * sizeof(cow))) == NULL) { 
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
			double dv[3];		/* Straight destination value */
			double bv[3];		/* Bent destination value */
			double wt = 1.0;	/* Default grey axis point weighting */

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

				wt *= (1.0 + t * brad);	/* Increase weigting with the bend */

			} else {
				t = 0.0;	/* stick to straight, it will be close anyway. */
			}

			for (j = 0; j < 3; j++)		/* full bend when t == 1 */
				gpnts[ngamp].v[j] = t * bv[j] + (1.0 - t) * dv[j];
			gpnts[ngamp].w = wt;

#ifdef NEVER
			printf("Grey axis %d maps %f %f %f -> %f %f %f wit %f\n",ngamp,
			gpnts[ngamp].p[0], gpnts[ngamp].p[1], gpnts[ngamp].p[2],
			gpnts[ngamp].v[0], gpnts[ngamp].v[1], gpnts[ngamp].v[2],
			gpnts[ngamp].w);
#endif
			ngamp++;
		}

		/* ---------------------------------------------------- */
		/* Deal with gamut hull guide vector creation. */

		/* For compression, create a mapping for each vertex of */
		/* the source gamut (image) surface towards the destination gamut */
		/* For expansion, do the opposite. */

		/* Convert from compact to explicit hextant weightings */
		if (expand_weights(xpweights, pweights)
		 || expand_weights(xsweights, sweights)) {
			fprintf(stderr,"gamut map: expand_weights() failed\n");
			s->grey->del(s->grey);
			s->igrey->del(s->igrey);
			if (sil_gam != scl_gam)
				sil_gam->del(sil_gam);
			scl_gam->del(scl_gam);
			free(s);
			return NULL;
		}

		/* Create weights as blend between perceptual and saturation */
		near_xwblend(xwh, xpweights, gmi->gampwf, xsweights, gmi->gamswf);
		if ((gmi->gampwf + gmi->gamswf) > 0.1)
			smooth = (gmi->gampwf * psmooth) + (gmi->gamswf * ssmooth);

		/* Tweak gamut mappings according to extra cmy cusp flags or rel override */ 
		if (dst_cmymap != 0 || rel_oride != 0) {
			tweak_weights(xwh, dst_cmymap, rel_oride); 
		}

		/* Create the near point mapping, which is our fundamental gamut */
		/* hull to gamut hull mapping. */
		nsm = near_smooth(verb, &nnsm, scl_gam, sil_gam, d_gam, src_kbp, dst_kbp,
		                  dr_be_bp, xwh, gmi->gamcknf, gmi->gamxknf,
		                  gmi->gamcpf > 1e-6, gmi->gamexf > 1e-6,
		                  xvra, mapres, m21fsm);
		if (nsm == NULL) {
			fprintf(stderr,"Creating smoothed near points failed\n");
			s->grey->del(s->grey);
			s->igrey->del(s->igrey);
			if (sil_gam != scl_gam)
				sil_gam->del(sil_gam);
			scl_gam->del(scl_gam);
			free(s);
			return NULL;
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

		/* ---------------------------------------------------- */
		/* Setup for diagnostic plot, that will have elements added */
		/* as we create the final 3D gamut mapping rspl */

		{
			int doaxes = 0;

#ifdef PLOT_AXES
			doaxes = 1;
#endif
			if (diagname != NULL)
				wrl = new_vrml(diagname, doaxes);
			else
				wrl = new_vrml("gammap.wrl", doaxes);
		}

		if (wrl != NULL) {
			/* See if there is a diagnostic locus to plot too */
			if ((locus = new_cgats()) == NULL)
				error("Failed to create cgats object");

			locus->add_other(locus, "TS");

			if (locus->read_name(locus, "locus.ts")) {
				locus->del(locus);
				locus = NULL;
			} else {
				if (verb)
					printf("!! Found diagnostic locus.ts file !!\n");
				/* locus will be added later */
			}

			/* Add diagnostic markers from markers structure */
			for (i = 0; ; i++) {
				double pp[3];
				co cp;
				if (markers[i].type == 0)
					break;
	
				if (markers[i].type == 1) {		/* Src point - do luminance mapping */
					/* Rotate and map marker points the same as the src gamuts */
					icmMul3By3x4(pp, s->grot, markers[i].pos);		/* Grey axis rotation */
					cp.p[0] = pp[0];			/* L value */
					s->grey->interp(s->grey, &cp);					/* Grey axis mapping */
					pp[0] = cp.v[0];
				} else {
					pp[0] = markers[i].pos[0];
					pp[1] = markers[i].pos[1];
					pp[2] = markers[i].pos[2];
				}
				wrl->add_marker(wrl, pp, markers[i].col, 1.0);
			}
		}

		/* --------------------------- */
		/* Now computue our 3D mapping points from the near point mapping. */
		for (i = 0; i < nnsm; i++) {
			double cpexf;			/* The effective compression or expansion factor */

			if (nsm[i].vflag == 0) {			/* Unclear whether compression or expansion */
				/* Use larger to the the two factors */
				cpexf = gmi->gamcpf > gmi->gamexf ? gmi->gamcpf : gmi->gamexf;
				
			} else if (nsm[i].vflag == 1) {		/* Compression */
				cpexf = gmi->gamcpf;

			} else if (nsm[i].vflag == 2) {		/* Expansion */
				cpexf = gmi->gamexf;

			} else {
				error("gammap: internal, unknown guide point flag");
			}

			/* Compute destination value which is a blend */
			/* between the source value and the fully mapped destination value. */
			icmBlend3(nsm[i].div, nsm[i].sv, nsm[i].dv, cpexf);

#ifdef NEVER
			printf("%s mapping:\n",nsm[i].vflag == 0 ? "Unclear" : nsm[i].vflag == 1 ? "Compression" : "Expansion");
			printf("Src point = %f %f %f radius %f\n",nsm[i].sv[0], nsm[i].sv[1], nsm[i].sv[2], nsm[i].sr);
			printf("Dst point = %f %f %f radius %f\n",nsm[i].dv[0], nsm[i].dv[1], nsm[i].dv[2], nsm[i].dr);
			printf("Blended dst point = %f %f %f\n",nsm[i].div[0], nsm[i].div[1], nsm[i].div[2]);
#endif	/* NEVER */
			/* Set the main gamut hull mapping point */
			for (j = 0; j < 3; j++) {
				gpnts[ngamp].p[j] = nsm[i].sv[j];
				gpnts[ngamp].v[j] = nsm[i].div[j];
			}
			gpnts[ngamp++].w = 1.0;		/* Main gamut surface mapping point */

#ifdef USE_GAMKNF
			/* Add sub surface mapping point if available */
			if (nsm[i].vflag != 0) {	/* Sub surface point is available */

				/* Compute destination value which is a blend */
				/* between the source value and the fully mapped destination value. */
				icmBlend3(nsm[i].div2, nsm[i].sv2, nsm[i].dv2, cpexf);

#ifdef NEVER
				printf("Src2 point = %f %f %f radius %f\n",nsm[i].sv2[0], nsm[i].sv2[1], nsm[i].sv2[2], nsm[i].sr);
				printf("Dst2 point = %f %f %f radius %f\n",nsm[i].dv2[0], nsm[i].dv2[1], nsm[i].dv2[2], nsm[i].dr);
				printf("Blended dst2 point = %f %f %f\n",nsm[i].div2[0], nsm[i].div2[1], nsm[i].div2[2]);
				printf("\n");
#endif	/* NEVER */
				/* Set the sub-surface gamut hull mapping point */
				for (j = 0; j < 3; j++) {
					gpnts[ngamp].p[j] = nsm[i].sv2[j];
					gpnts[ngamp].v[j] = nsm[i].div2[j];
				}
				gpnts[ngamp++].w = nsm[i].w2;		/* Sub-suface mapping points */
			}
#endif /* USE_GAMKNF */
		}

		/* Create preliminary gamut mapping rspl, without grid boundary values. */
		/* We use this to lookup the mapping for points on the source space gamut */
		/* that result from clipping our grid boundary points */
#ifdef USE_BOUND
		for (j = 0; j < 3; j++) {		/* Set resolution for all axes */
			gres[j] = mapres/2;
			avgdev[j] = 0.005;
		}
		s->map = new_rspl(RSPL_NOFLAGS, 3, 3);	/* Allocate 3D -> 3D */
		s->map->fit_rspl_w(s->map, RSPLFLAGS, gpnts, ngamp, il, ih, gres, ol, oh, smooth, avgdev, NULL);

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
					gpnts[ngamp++].w = 0.05 * ww;
				}
				DC_INC(gc);
				if (DC_DONE(gc))
					break;
			}
		}
#else /* !USE_BOUND */
		printf("!!!! Warning - gammap boundary points disabled !!!!\n");
#endif /* !USE_BOUND */

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
		/* [ The smoothing is not as useful as it should be, because */
		/*   if it is increased it tends to push colors out of gamut */
		/*   where they get clipped. Some cleverer scheme which makes */
		/*   sure that smoothness errs on the side of more compression */
		/*   is needed. ] */
		if (s->map != NULL)
			s->map->del(s->map);
		if (verb)
			printf("Creating rspl..\n");
		for (j = 0; j < 3; j++) {		/* Set resolution for all axes */
			gres[j] = mapres;
			avgdev[j] = 0.005;
		}
		s->map = new_rspl(RSPL_NOFLAGS, 3, 3);	/* Allocate 3D -> 3D */
		if (s->map->fit_rspl_w(s->map, MAINRSPLFLAGS, gpnts, ngamp, il, ih, gres, ol, oh, smooth, avgdev, NULL)) {
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
						div[j] = gmi->gamcpf * nsm[i].dv[j] + (1.0 - gmi->gamcpf) * nsm[i].sv[j];
					}

				} else {	/* Expansion needed */

					/* Compute expansion destination value */
					for (j = 0; j < 3; j++) {
						cp.p[j] = nsm[i].sv[j];
						div[j] = gmi->gamexf * nsm[i].dv[j] + (1.0 - gmi->gamexf) * nsm[i].sv[j];
					}
				}
				s->map->interp(s->map, &cp);
				
				de = icmLabDE(div, cp.v);
				avgde += de;
				if (de > maxde)
					maxde = de;
			}
			avgde /= nnsm;
			printf("Gamut hull fit to guides: = avg %f, max %f\n",avgde,maxde);
		}
#endif /* CHECK_NEARMAP */

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

			if (verb)
				printf("Fine tuning white and black point mapping\n");

			/* Check what the source white and black points actually maps to */
			s->domap(s, a_wp, s_mt_wp);
			s->domap(s, a_bp, s_mt_bp);

#ifdef VERBOSE
			if (verb) {
				printf("White is %f %f %f, should be %f %f %f\n",
				a_wp[0], a_wp[1], a_wp[2], d_mt_wp[0], d_mt_wp[1], d_mt_wp[2]);
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
			/* We do fine tune it if dst_kbp is set though, since */
			/* we would like perfect K only out. */

			/* Compute rotation/scale relative white point matrix */
			icmVecRotMat(cx.mat, a_wp, a_bp, d_mt_wp, d_mt_bp);		/* wp & bp */

			/* Fine tune the 3D->3D mapping */
			s->map->re_set_rspl(
				s->map,				/* this */
				0,					/* Combination of flags */
				(void *)&cx,		/* Opaque function context */
				adjust_wb_func 		/* Function to set from */
			);

#ifdef VERBOSE
			if (verb) {
				/* Check what the source white and black points actually maps to */
				s->domap(s, a_wp, s_mt_wp);
				s->domap(s, a_bp, s_mt_bp);
	
				printf("After fine tuning:\n");
				printf("White is %f %f %f, should be %f %f %f\n",
				a_wp[0], a_wp[1], a_wp[2], d_mt_wp[0], d_mt_wp[1], d_mt_wp[2]);
				printf("Black is %f %f %f, should be %f %f %f\n",
				a_bp[0], a_bp[1], a_bp[2], d_mt_bp[0], d_mt_bp[1], d_mt_bp[2]);
			}
#endif /* VERBOSE */
		}

		if (wrl != NULL) {
			int arerings = 0;
			double cc[3] = { 0.7, 0.7, 0.7 }; 
			double nc[3] = { 1.0, 0.4, 0.7 }; 	/* Pink for neighbors */
			int nix = -1;						/* Index of point to show neighbour */

#ifdef SHOW_NEIGBORS
#ifdef NEVER
			/* Show all neighbours */
			wrl->start_line_set(wrl, 0);
			for (i = 0; i < nnsm; i++) {
				for (j = 0; j < XNNB; j++) {
					nearsmth *np = nsm[i].n[j];		/* Pointer to neighbor */
			
					if (np == NULL)
						break;

					wrl->add_col_vertex(wrl, 0, nsm[i].sv, nc);	/* Source value */
					wrl->add_col_vertex(wrl, 0, np->sv, nc);	/* Neighbpor value */
				}
			}
			wrl->make_lines(wrl, 0, 2);
#else
			/* Show neighbours of points near source markers */
			for (i = 0; ; i++) {	/* Add diagnostic markers */
				double pp[3];
				co cp;
				int ix, bix;
				double bdist = 1e6;

				if (markers[i].type == 0)
					break;
	
				if (markers[i].type != 1)
					continue;
				
				/* Rotate and map marker point the same as the src gamuts */
				icmMul3By3x4(pp, s->grot, markers[i].pos);
				cp.p[0] = pp[0];			/* L value */
				s->grey->interp(s->grey, &cp);
				pp[0] = cp.v[0];
//printf("~1 looking for closest point to marker %d at %f %f %f\n",i,pp[0],pp[1],pp[2]);

				/* Locate the nearest source point */
				for (ix = 0; ix < nnsm; ix++) {
					double dist = icmNorm33(pp, nsm[ix].sv);
					if (dist < bdist) {
						bdist = dist;
						bix = ix;
					}
				}
//printf("~1 closest src point ix %d at %f %f %f\n",bix,nsm[bix].sv[0],nsm[bix].sv[1],nsm[bix].sv[2]);
//printf("~1 there are %d neighbours\n",nsm[bix].nnb);

				wrl->start_line_set(wrl, 0);
				for (j = 0; j < nsm[bix].nnb; j++) {
					nearsmth *np = nsm[bix].n[j].n;	/* Pointer to neighbor */

					wrl->add_col_vertex(wrl, 0, nsm[bix].sv, nc);	/* Source value */
					wrl->add_col_vertex(wrl, 0, np->sv, nc);		/* Neighbpor value */
				}
				wrl->make_lines(wrl, 0, 2);
			}
#endif
#endif /* SHOW_NEIGBORS */

			/* Add the source and dest gamut surfaces */
#ifdef PLOT_SRC_GMT
			wrl->make_gamut_surface_2(wrl, sil_gam, 0.6, 0, cc);
#endif /* PLOT_SRC_GMT */
#ifdef PLOT_DST_GMT
			cc[0] = -1.0;
			wrl->make_gamut_surface(wrl, d_gam, 0.2, cc);
#endif /* PLOT_DST_GMT */
#ifdef PLOT_SRC_CUSPS
			wrl->add_cusps(wrl, sil_gam, 0.6, NULL);
#endif /* PLOT_SRC_CUSPS */
#ifdef PLOT_DST_CUSPS
			wrl->add_cusps(wrl, d_gam, 0.2, NULL);
#endif /* PLOT_DST_CUSPS */

#ifdef PLOT_TRANSSRC_CUSPS
			/* Add transformed source cusp markers */
			{
				int i;
				double cusps[6][3];
				double ccolors[6][3] = {
					{ 1.0, 0.1, 0.1 },		/* Red */
					{ 1.0, 1.0, 0.1 },		/* Yellow */
					{ 0.1, 1.0, 0.1 },		/* Green */
					{ 0.1, 1.0, 1.0 },		/* Cyan */
					{ 0.1, 0.1, 1.0 },		/* Blue */
					{ 1.0, 0.1, 1.0 }		/* Magenta */
				};
			
				if (sc_gam->getcusps(sc_gam, cusps) == 0) {

					for (i = 0; i < 6; i++) {
						double val[3];

						s->domap(s, val, cusps[i]);
						wrl->add_marker(wrl, val, ccolors[i], 2.5);
					}
				}
			}
#endif

#if defined(SHOW_MAP_VECTORS) || defined(SHOW_SUB_SURF) || defined(SHOW_ACTUAL_VECTORS) || defined(SHOW_ACTUAL_VEC_DIFF)
			/* Start of guide vector plot */
			wrl->start_line_set(wrl, 0);

			for (i = 0; i < nnsm; i++) {
				double cpexf;			/* The effective compression or expansion factor */
				double yellow[3] = { 1.0, 1.0, 0.0 };
				double red[3]    = { 1.0, 0.0, 0.0 };
				double green[3]  = { 0.0, 1.0, 0.0 };
				double lgrey[3]  = { 0.8, 0.8, 0.8 };
				double purp[3]   = { 0.6, 0.0, 1.0 };
				double blue[3]   = { 0.2, 0.2, 1.0 };
				double *ccc;
				double mdst[3];

#if defined(SHOW_ACTUAL_VECTORS) || defined(SHOW_ACTUAL_VEC_DIFF)
# ifdef SHOW_ACTUAL_VECTORS
				wrl->add_col_vertex(wrl, 0, nsm[i].sv, yellow);
# else	/* SHOW_ACTUAL_VEC_DIFF */
				wrl->add_col_vertex(wrl, 0, nsm[i].div, yellow);
# endif
				s->domap(s, mdst, nsm[i].sv);
				wrl->add_col_vertex(wrl, 0, mdst, red);

#else
# ifdef SHOW_MAP_VECTORS
				ccc = yellow;

				if (nsm[i].gflag == 0)
					ccc = green;			/* Mark "no clear direction" vectors in green->red */
#  ifdef SHOW_CUSPMAP
				wrl->add_col_vertex(wrl, 0, nsm[i].csv, ccc);	/* Cusp mapped source value */
#  else
				wrl->add_col_vertex(wrl, 0, nsm[i].sv, ccc);	/* Source value */
#  endif
				wrl->add_col_vertex(wrl, 0, nsm[i].div, red);	/* Blended destination value */
# endif /* SHOW_MAP_VECTORS */

# ifdef SHOW_SUB_SURF
				if (nsm[i].vflag != 0) {	/* Sub surface point is available */

					wrl->add_col_vertex(wrl, 0, nsm[i].sv2, lgrey); /* Subs-surf Source value */
					wrl->add_col_vertex(wrl, 0, nsm[i].div2, purp); /* Blended destination value */
				}
# endif /* SHOW_SUB_SURF */
#endif	/* !SHOW_ACTUAL_VECTORS */
			}
			wrl->make_lines(wrl, 0, 2);			/* Guide vectors */
#endif	/* Show vectors */

#ifdef SHOW_VECTOR_INDEXES
			for (i = 0; i < nnsm; i++) {
				double cream[3] = { 0.7, 0.7, 0.5 };
				char buf[100];
				sprintf(buf, "%d", i);
				wrl->add_text(wrl, buf, nsm[i].sv, cream, 0.5);
			}
#endif /* SHOW_VECTOR_INDEXES */

			/* add the locus from locus.ts file */ 
			if (locus != NULL) {
				int table, npoints;
				char *fnames[3] = { "LAB_L", "LAB_A", "LAB_B" };
				int ix[3];
				double v0[3], v1[3];
				double rgb[3];

				/* Each table holds a separate locus */
				for (table = 0; table < locus->ntables; table++) {

				    if ((npoints = locus->t[table].nsets) <= 0)
				        error("No sets of data in diagnostic locus");

					for (j = 0; j < 3; j++) {
						if ((ix[j] = locus->find_field(locus, 0, fnames[j])) < 0)
							error ("Locus file doesn't contain field %s",fnames[j]);
						if (locus->t[table].ftype[ix[j]] != r_t)
							error ("Field %s is wrong type",fnames[j]);
					}

					/* Source locus */
					rgb[0] = 1.0;
					rgb[1] = 0.5;
					rgb[2] = 0.5;
					for (i = 0; i < npoints; i++) {
						co cp;

						for (j = 0; j < 3; j++)
							v1[j] = *((double *)locus->t[table].fdata[i][ix[j]]);

						/* Rotate and locus verticies the same as the src gamuts */
						icmMul3By3x4(v1, s->grot, v1);
						cp.p[0] = v1[0];			/* L value */
						s->grey->interp(s->grey, &cp);
						v1[0] = cp.v[0];

						if (i > 0 )
							wrl->add_cone(wrl, v0, v1, rgb, 0.5);
						icmAry2Ary(v0,v1);
					}

					/* Destination locus */
					rgb[0] = 1.0;
					rgb[1] = 1.0;
					rgb[2] = 1.0;
					for (i = 0; i < npoints; i++) {
						co cp;

						for (j = 0; j < 3; j++)
							v1[j] = *((double *)locus->t[table].fdata[i][ix[j]]);

						s->domap(s, v1, v1);

						if (i > 0 )
							wrl->add_cone(wrl, v0, v1, rgb, 0.5);
						icmAry2Ary(v0,v1);
					}
				}
				
				locus->del(locus);
				locus = NULL;
			}

			/* Add any ring mapping diagnostics */
			for (i = 0; ; i++) {
				if (rings[i].type == 0)
					break;

				if (rings[i].type == 2)
					continue;
	
				if (rings[i].type == 1) {
					double pconst;
					double cpoint[3];
					double mat[3][4];		/* translate to our plane */
					double imat[3][4];		/* translate from our plane */
					double s1[3], s0[3], t1[3];
					int j;
					double maxa, mina;
					double maxb, minb;

					if (arerings == 0) {
						arerings = 1;
						wrl->start_line_set(wrl, 1);	/* Source ring */
						wrl->start_line_set(wrl, 2);	/* Destination ring */
					}

					if (icmNormalize3(rings[i].pnorm, rings[i].pnorm, 1.0))
						error("Ring %d diagnostic plane normal failed",i);

					pconst = -icmDot3(rings[i].ppoint, rings[i].pnorm);

					/* Locate intersection of source neautral axis and plane */
					if (icmVecPlaneIsect(cpoint, pconst, rings[i].pnorm, s_cs_wp, s_cs_bp))
						error("Ring %d diagnostic center point intersection failed",i);

					/* Compute the rotation and translation between */
					/* a plane in ab and the plane we are using */
					s0[0] = s0[1] = s0[2] = 0.0;
					s1[0] = 1.0, s1[1] = s1[2] = 0.0;
					t1[0] = cpoint[0] + rings[i].pnorm[0];
					t1[1] = cpoint[1] + rings[i].pnorm[1];
					t1[2] = cpoint[2] + rings[i].pnorm[2];
					icmVecRotMat(mat, s1, s0, t1, cpoint);
					icmVecRotMat(imat, t1, cpoint, s1, s0);

					/* Do a min/max of a circle of vectors so as to */
					/* establish an offset to the centroid for this slice */
					maxa = maxb = -1e60;
					mina = minb = 1e60;
					for (j = 0; j < 20; j++) {
						double ang = 2 * 3.1415926 * j/(20 - 1.0);
						double vec[3], isect[3];
						double pp[3];
						co cp;
						int k;

						vec[0] = 0.0;
						vec[1] = sin(ang);
						vec[2] = cos(ang);
						icmMul3By3x4(vec, mat, vec);

						/* Intersect it with the source gamut */
						if (si_gam->vector_isect(si_gam, vec, cpoint, isect,
						                 NULL, NULL, NULL, NULL, NULL) == 0) {
							continue;
						}

						/* Translate back to plane */
						icmMul3By3x4(pp, imat, isect);

						if (pp[1] > maxa)
							maxa = pp[1];
						if (pp[1] < mina)
							mina = pp[1];
						if (pp[2] > maxb)
							maxb = pp[2];
						if (pp[2] < minb)
							minb = pp[2];
					}
					/* Move center to centroid of min/max box */
					t1[0] = 0.0;
					t1[1] = (maxa + mina) * 0.5;
					t1[2] = (maxb + minb) * 0.5;
					if (t1[1] < -200.0 || t1[1] > 200.0
					 || t1[2] < -200.0 || t1[2] > 200.0)
						error("Failed to locate centroid of slice");
					icmMul3By3x4(cpoint, mat, t1);
					
//printf("~1 ring centroid point = %f %f %f\n", cpoint[0],cpoint[1],cpoint[2]);
		
					/* Recompute the rotation and translation between */
					/* a plane in ab and the plane we are using */
					s0[0] = s0[1] = s0[2] = 0.0;
					s1[0] = 1.0, s1[1] = s1[2] = 0.0;
					t1[0] = cpoint[0] + rings[i].pnorm[0];
					t1[1] = cpoint[1] + rings[i].pnorm[1];
					t1[2] = cpoint[2] + rings[i].pnorm[2];
					icmVecRotMat(mat, s1, s0, t1, cpoint);
					icmVecRotMat(imat, t1, cpoint, s1, s0);

//printf("~1 generating %d ring verts\n",rings[i].nverts);
					/* Create a circle of vectors in the plane from the center */
					/* point, to intersect with the source gamut surface. */
					/* (Duplicate start and end vertex) */
					for (j = 0; j <= rings[i].nverts; j++) {
						double ang = 2 * 3.1415926 * j/((double) rings[i].nverts);
						double vec[3], isect[3];
						double pp[3];
						co cp;
						int k;

						vec[0] = 0.0;
						vec[1] = sin(ang);
						vec[2] = cos(ang);
						icmMul3By3x4(vec, mat, vec);

						/* Intersect it with the source gamut */
						if (si_gam->vector_isect(si_gam, vec, cpoint, isect,
						                 NULL, NULL, NULL, NULL, NULL) == 0) {
							warning("Ring %d vect %d diagnostic vector intersect failed",i,j);
							continue;
						}

//printf("~1 vec %d = %f %f %f\n",j,isect[0],isect[1],isect[2]);

						/* Scale them to the ratio */
						for (k = 0; k < 3; k++)
							vec[k] = isect[k] * rings[i].rad + (1.0 - rings[i].rad) * cpoint[k];

//printf("~1 rad vec %d = %f %f %f\n",j,vec[0],vec[1],vec[2]);

						/* Transform them into rotated and scaled destination space */
						icmMul3By3x4(vec, s->grot, vec);
						cp.p[0] = vec[0];			/* L value */
						s->grey->interp(s->grey, &cp);
						vec[0] = cp.v[0];
//printf("~1 trans vec %d = %f %f %f\n",j,vec[0],vec[1],vec[2]);

						/* Add to plot */
						wrl->add_col_vertex(wrl, 1, vec, rings[i].scol);
//printf("~1 src vec %d = %f %f %f\n",j,vec[0],vec[1],vec[2]);

						/* Gamut map and add to plot */
						s->domap(s, vec, vec);
//printf("~1 dst vec %d = %f %f %f\n",j,vec[0],vec[1],vec[2]);
						wrl->add_col_vertex(wrl, 2, vec, rings[i].dcol);
					}
					wrl->make_last_vertex(wrl, 1);		/* Source ring */
					wrl->make_last_vertex(wrl, 2);		/* Destination ring */
				}
				if (arerings) {
					wrl->make_lines(wrl, 1, 1000000);		/* Source ring */
					wrl->make_lines(wrl, 2, 1000000);		/* Destination ring */
				}
			}

			wrl->del(wrl);		/* Write and delete */
			wrl = NULL;
		}

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


		free_nearsmth(nsm, nnsm);
	}

#ifdef PLOT_GAMUTS
	scl_gam->write_vrml(scl_gam, "src.wrl", 1, 0);
	sil_gam->write_vrml(sil_gam, "img.wrl", 1, 0);
	d_gam->write_vrml(d_gam, "dst.wrl", 1, 0);
	sc_gam->write_trans_vrml(sc_gam, "gmsrc.wrl", 1, 0, map_trans, s);
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
	if (s->grey != NULL)
		s->grey->del(s->grey);
	if (s->igrey != NULL)
		s->igrey->del(s->igrey);
	if (s->map != NULL)
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

	/* If there is a 3D->3D mapping */
	if (s->map != NULL) {
		int e;

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
	
		for (e = 0; e < s->map->fdi; e++)
			out[e] = cp.v[e];
	
		if (s->dbg) printf("domap: after 3D map %s\n\n",icmPdv(s->map->fdi, out));
	} else {
		out[0] = cp.v[0];
		out[1] = rin[1];
		out[2] = rin[2];
	}
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
	if (p->dst->vector_isect(p->dst, cp, out, t2, t1, &p2, &p1, NULL, NULL) != 0) {

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





















