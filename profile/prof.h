#ifndef PROF_H
#define PROF_H
/* 
 * ICC Profile creation library.
 *
 * Author:  Graeme W. Gill
 * Date:    11/10/00
 * Version: 1.00
 *
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This library provide high level routines to create device ICC
 * profiles from argyll cgats patch test data.
 */


/* Profile algorithm type */
typedef enum {
	prof_default          = 0,		/* Default for type of device */
	prof_clutLab          = 1,		/* Lab clut. */
	prof_clutXYZ          = 2,		/* XYZ clut. */
	prof_gammat           = 3,		/* XYZ gamut + matrix */
	prof_shamat           = 4,		/* XYZ shaper + matrix */
	prof_gam1mat          = 5,		/* XYZ shared TRC gamut + matrix */
	prof_sha1mat          = 6		/* XYZ shared TRC shaper + matrix */
} prof_atype;

/* Output or Display device */
void make_output_icc(
	prof_atype ptype,		/* Profile output type */
	int mtxtoo,				/* NZ if matrix tags should be created for Display XYZ cLUT */
	icmICCVersion iccver,	/* ICC profile version to create */
	int verb,				/* Vebosity level, 0 = none */
	int iquality,			/* A2B table quality, 0..2 */
	int oquality,			/* B2A table quality, 0..2 */
	int noiluts,			/* nz to supress creation of input (Device) shaper luts */
	int noisluts,			/* nz to supress creation of input sub-grid (Device) shaper luts */
	int nooluts,			/* nz to supress creation of output (PCS) shaper luts */
	int nocied,				/* nz to supress inclusion of .ti3 data in profile */
	int noptop,				/* nz to use colorimetic source gamut to make perceptual table */ 
	int nostos,				/* nz to use colorimetic source gamut to make perceptual table */
	int gamdiag,			/* Make gamut mapping diagnostic wrl plots */
	int verify,				/* nz to print verification */
	icxInk *ink,			/* Ink limit/black generation setup */
	char *in_name,			/* input .ti3 file name */
	char *file_name,		/* output icc name */
	cgats *icg,				/* input cgats structure */
	int spec,				/* Use spectral data flag */
	icxIllumeType illum,	/* Spectral illuminant */
	xspect *cust_illum,		/* Possible custom illumination */
	icxObserverType observ,	/* Spectral observer */
	int fwacomp,			/* FWA compensation requested */
	double smooth,			/* RSPL smoothing factor, -ve if raw */
	double avgdev,			/* reading Average Deviation as a proportion of the input range */
	char *ipname,			/* input icc profile - enables gamut map, NULL if none */
	char *sgname,			/* source image gamut - NULL if none */
	char *absname,			/* abstract profile name - NULL if none */
	int sepsat,				/* Create separate Saturation B2A */
	icxViewCond *ivc_p,		/* Input Viewing Parameters for CIECAM97s */
	icxViewCond *ovc_p,		/* Output Viewing Parameters for CIECAM97s (enables CAM clip) */
	int ivc_e,				/* Input Enumerated viewing condition */
	int ovc_e,				/* Output Enumerated viewing condition */
	icxGMappingIntent *pgmi,/* Perceptual gamut mapping intent */
	icxGMappingIntent *sgmi,/* Saturation gamut mapping intent */
	profxinf *pi			/* Optional Profile creation extra data */
);

/* Input device */
void make_input_icc(
	prof_atype ptype,		/* Profile output type */
	icmICCVersion iccver,	/* ICC profile version to create */
	int verb,				/* Vebosity level, 0 = none */
	int iquality,			/* A2B table quality, 0..2 */
	int oquality,			/* B2A table quality, 0..2 */
	int noiluts,			/* nz to supress creation of input (Device) shaper luts */
	int noisluts,			/* nz to supress creation of input sub-grid (Device) shaper luts */
	int nooluts,			/* nz to supress creation of output (PCS) shaper luts */
	int nocied,				/* nz to supress inclusion of .ti3 data in profile */
	int verify,				/* nz to print verification */
	int nsabs,				/* nz for non-standard absolute output */
	double iwpscale,		/* >= 0.0 for media white point scale factor */
	int dob2a,				/* nz to create a B2A table as well */
	char *in_name,			/* input .ti3 file name */
	char *file_name,		/* output icc name */
	cgats *icg,				/* input cgats structure */
	int spec,				/* Use spectral data flag */
	icxIllumeType illum,	/* Spectral illuminant */
	xspect *cust_illum,		/* Possible custom illumination */
	icxObserverType observ,	/* Spectral observer */
	double smooth,			/* RSPL smoothing factor, -ve if raw */
	double avgdev,			/* reading Average Deviation as a proportion of the input range */
	profxinf *pi			/* Optional Profile creation extra data */
);

#endif /* PROF_H */
