
/* 
 * Argyll Color Correction System
 * Output device profile creator.
 *
 * Author: Graeme W. Gill
 * Date:   11/10/00
 *
 * Copyright 2000-2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/*
 * This program takes in the scattered test chart
 * points, and interpolates them into a gridded
 * forward ICC device profile for an output or display device,
 * using a clut based profle.
 * It also creates backward tables based on the forward grid.
 *
 * Preview profiles are not currently generated.
 * 
 * The gamut clut should be implemented with xicc/rspl
 */

/*
 * TTBD:
 *
 *  Need to make this more of a library:
 *  Fix error handling
 *  fix verbose output
 *  hand icc object back rather than writing file ?
 */

#define DEBUG
#define verbo stdout

#define IMP_MONO			/* Turn on development code */

#undef USE_CORRECT_GAMUT		/* Use correct gamut boundary code. Not working well yet */
#define NO_B2A_PCS_CURVES		/* PCS curves seem to make B2A less accurate. why ? */
#undef USE_CAM_CLIP_OPT			/* Clip out of gamut in CAM space rather than XYZ or L*a*b* */
#undef USE_EXTRA_FITTING		/* Turn on extra effort at fitting if medium or higher quality. */
								/* Reduces noise filtering ?? - introduce smoothness control instead */

#include <stdio.h>
#include "cgats.h"
#include "icc.h"
#include "xicc.h"
#include "numlib.h"
#include "rspl.h"
#include "insttypes.h"
#include "prof.h"
#include "gamut.h"
#include "gammap.h"

/*
   Basic algorithm outline:

 Printer:
   Figure out the input curves to give
   the flattest grid.

   Figure out the grid values.

   Use them to generate all the A2B tables.

   Use the inverse rspl code to compute
   all the B2A profiles.

*/

/*
	Notes:

	The shared gamma/shaper support is for silly applications
	like which can't handle display profiles that have per
	chanel gamma's, per chanel curves or clut based display profiles.

*/

/* NOTE:-
	It's interesting that the white and black points recorded in the tags,
	generally won't quite match the white and black points returned by
	looking up the profile in absolute mode.

	For a Matrix profile, in the case of the white point this is
	because we're not using the ICC 16 bit quantized value to
	create the relative transform matrix, and in the case of
	the black point, it can never be a perfect match because the black
	point returned by a profile lookup will be the quantized black
	point of the matrix, transformed by the rel->abs matrix, which
	generally won't be equal to an ICC quantized value.
	It might help the latter case if we were at least able to convert
	the profile quantized black point into the absolute black point
	via the rel->abs transform, and then quantize it.

	Of course all of this will be worse in the Lut type profile,
	due to the white and black points being stored in a different
	quantized space (XYZ vs. Lab) than the Lut grid point values!
 */


/* ---------------------------------------- */

/* structure to support output icc B2A Lut initialisation calbacks */
/* Note that we don't cope with matrix - assume it's unity. */

typedef struct {
	int verb;
	int total, count, last;	/* Progress count information */
	int noPCScurves;		/* Flag set if we don't want PCS curves */
	icColorSpaceSignature pcsspace;	/* The PCS colorspace */
	icColorSpaceSignature devspace;	/* The device colorspace */
	icxLuLut *x;			/* A2B icxLuLut we are inverting in std PCS */

	int ntables;			/* Number of tables being set. 1 = colorimetric */
							/* 2 = colorimetric + saturation, 3 = all intents */
	int ochan;				/* Number of output channels for B2A */
	gammap *pmap;			/* Perceptual CAM to CAM Gamut mapping, NULL if no mapping */
	gammap *smap;			/* Saturation CAM to CAM Gamut mapping, NULL if no mapping */
	icxLuBase *ix;			/* Source profile std PCS to CAM conversion */
	icxLuBase *ox;			/* Destination profile CAM to std PCS conversion */

	icRenderingIntent abs_intent;	/* Desired abstract profile rendering intent */
	icxLuBase *abs_luo;		/* abstract profile tranform in PCS, NULL if none */

    gamut *gam;				/* Output gamut object for setting gamut Lut */
	int wantLab;			/* 0 if want is XYZ PCS, 1 want is Lab PCS */
} out_b2a_callback;

/* Utility to handle abstract profile application to PCS */
/* PCS in creating output table is always XYZ or Lab relative colorimetric, */
/* and abstract profile is absolute or relative, and will be */
/* XYZ if absolute, and PCS if relative. */
static void do_abstract(out_b2a_callback *p, double out[3], double in[3]) {
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];

	if (p->abs_intent == icAbsoluteColorimetric) {
		if (p->pcsspace == icSigLabData) {
			icmLab2XYZ(&icmD50, out, out);
		}
		p->x->plu->XYZ_Rel2Abs(p->x->plu, out, out);
	}

	p->abs_luo->lookup(p->abs_luo, out, out);

	if (p->abs_intent == icAbsoluteColorimetric) {
		p->x->plu->XYZ_Abs2Rel(p->x->plu, out, out);
		if (p->pcsspace == icSigLabData) {
			icmXYZ2Lab(&icmD50, out, out);
		}
	}
}


/* Input table is the inverse of the AtoB output table */
/* Input PCS output PCS' */
void out_b2a_input(void *cntx, double out[3], double in[3]) {
	out_b2a_callback *p = (out_b2a_callback *)cntx;

	if (p->noPCScurves) {
		out[0] = in[0];
		out[1] = in[1];
		out[2] = in[2];
	} else {
		if (p->x->inv_output(p->x, out, in) > 1)
			error("%d, %s",p->x->pp->errc,p->x->pp->err);
	}
}


/* clut - multitable */
/* Input PCS' output Dev' */
void out_b2a_clut(void *cntx, double *out, double in[3]) {
	out_b2a_callback *p = (out_b2a_callback *)cntx;
	double inn[3], in1[3];
	int tn;

//printf("~1 out_b2a_clut got %f %f %f\n",in[0],in[1],in[2]);

	inn[0] = in1[0] = in[0];		/* in[] may be aliased with out[] */
	inn[1] = in1[1] = in[1];		/* so take two copies.  */
	inn[2] = in1[2] = in[2];

	if (p->abs_luo != NULL) {	/* Abstract profile to apply to first table */

		if (!p->noPCScurves) {	/* Convert from PCS' to PCS */
			if (p->x->output(p->x, in1, in1) > 1)
				error("%d, %s",p->x->pp->errc,p->x->pp->err);
		}
		do_abstract(p, in1, in1);		/* Abstract profile to apply to first table */
		/* We now have PCS */
	}

	if (p->noPCScurves || p->abs_luo != NULL) {	/* We were given PCS or have converted to PCS */

		/* PCS to PCS' */
		if (p->x->inv_output(p->x, in1, in1) > 1)
			error("%d, %s",p->x->pp->errc,p->x->pp->err);
	}

	/* Invert AtoB clut (PCS' to Dev') Colorimetric */
	if (p->x->inv_clut(p->x, out, in1) > 1)
		error("%d, %s",p->x->pp->errc,p->x->pp->err);

	if (p->ntables > 1) {		/* Do first part once for both intents */

		/* Starting with original input inn[] */
		if (p->noPCScurves) {	/* We have PCS because PCS Curve is combined with clut */
			in1[0] = inn[0];	/* We were actually given PCS */
			in1[1] = inn[1];
			in1[2] = inn[2];

		} else {				/* We were given PCS' */
			/* Convert from PCS' to PCS */
			if (p->x->output(p->x, in1, inn) > 1)
				error("%d, %s",p->x->pp->errc,p->x->pp->err);
		}

		/* Convert from PCS to CAM */
		p->ix->fwd_relpcs_outpcs(p->ix, p->pcsspace, in1, in1);

		/* Apply gamut mapping in CAM space for remaining tables */
		for (tn = 1; tn < p->ntables; tn++) {
			double in2[3];

			out += p->ochan;	/* next table/intent */

			if (tn == 1)
				p->pmap->domap(p->pmap, in2, in1);	/* Perceptual mapping */
			else
				p->smap->domap(p->smap, in2, in1);	/* Saturation mapping */

			/* Convert from CAM to PCS */
			p->ox->bwd_outpcs_relpcs(p->ox, p->pcsspace, in2, in2);

			if (p->abs_luo != NULL)
				do_abstract(p, in2, in2);	/* Abstract profile to other tables after gamut map */

			/* Convert from PCS to PCS' */
			if (p->x->inv_output(p->x, in2, in2) > 1)
				error("%d, %s",p->x->pp->errc,p->x->pp->err);

			/* Invert AtoB clut (PCS' to Dev') */
			if (p->x->inv_clut(p->x, out, in2) > 1)
				error("%d, %s",p->x->pp->errc,p->x->pp->err);
		}
	}

	if (p->verb) {		/* Output percent intervals */
		int pc;
		p->count++;
		pc = p->count * 100.0/p->total + 0.5;
		if (pc != p->last) {
			printf("\r%2d%%",pc); fflush(stdout);
			p->last = pc;
		}
	}
}

/* Output table is the inverse of the AtoB input table */
/* Input Dev' output Dev */
void out_b2a_output(void *cntx, double out[4], double in[4]) {
	out_b2a_callback *p = (out_b2a_callback *)cntx;

	if (p->x->inv_input(p->x, out, in) > 1)
		error("%d, %s",p->x->pp->errc,p->x->pp->err);
}

#ifdef USE_CORRECT_GAMUT
/* PCS' -> distance to gamut boundary */
static void PCSp_bdist(void *cntx, double out[1], double in[3]) {
	out_b2a_callback *p = (out_b2a_callback *)cntx;
	double pcs[3];				/* PCS value of input */
	double rad[3];				/* Radial coords of input point */
	double nrad[3];				/* Radial Lab value */
	double rrad;				/* Radius of radial point on gamut surface */
	double gdist;				/* Out of gamut distance */
	
printf("~1 bdist got PCS %f %f %f\n",in[0],in[1],in[2]);
	/* Do PCS' -> PCS */
	if (p->x->inv_output(p->x, pcs, in) > 1)
		error("%d, %s",p->x->pp->errc,p->x->pp->err);

	/* If PCS is XYZ, convert to Lab */
	if (p->wantLab == 0) {
		icmXYZ2Lab(&icmD50, pcs, pcs);
printf("~1 bdist converted to Lab %f %f %f\n",pcs[0],pcs[1],pcs[2]);
	}

	/* Compute radial distance of input point */
	gamut_rect2radial(rad, pcs);
printf("~1 bdist got radial distance %f\n",rad[0]);

	/* Do PCS -> radial intersection on the gamut */
	rrad = p->gam->radial(p->gam, nrad, pcs);
printf("~1 bdist got nearest point on gamut %f %f %f\n",nrad[0],nrad[1],nrad[2]);
	
	gdist = rad[0] - rrad;

printf("~1 bdist got gdist %f\n",gdist);

	/* Distance in PCS space will be roughly -128 -> 128 */
	/* Clip to range -30 - +30, then scale to 0.0 - 1.0 */
	if (gdist < -30.0)
		gdist = -30.0;
	else if (gdist > 30.0)
		gdist = 30.0;

	out[0] = (gdist + 30.0)/60.0;
printf("~1 bdist returning %f\n",out[0]);
}

/* The output table for a Gamut lut is usually a special function, returning */
/* a value of 0 for all inputs <= 0.5, and then outputing between */
/* 0.0 and 1.0 for the input range 0.5 to 1.0. This is so a graduated */
/* "gamut boundary distance" number from the multi-d lut can be */
/* translated into the ICC "0.0 if in gamut, > 0.0 if not" number. */
static void gamut_output(void *cntx, double out[1], double in[1]) {
	double iv, ov;
	iv = in[0];
	if (iv <= 0.5)
		ov = 0.0;
	else
		ov = (iv - 0.5) * 2.0;
	out[0] = ov;
}

#endif /* USE_CORRECT_GAMUT */

/* Make an output device profile, where a forward mapping is from */
/* CMYK to Lab space */
void
make_output_icc(
	prof_atype ptype,		/* Profile algorithm type */
	int verb,				/* Vebosity level, 0 = none */
	int iquality,			/* A2B table quality, 0..3 */
	int oquality,			/* B2A table quality, 0..3 */
	int noiluts,			/* nz to supress creation of input (Device) shaper luts */
	int nooluts,			/* nz to supress creation of output (PCS) shaper luts */
	int verify,				/* nz to print verification */
	icxInk *ink,			/* Ink limit/black generation setup (NULL if n/a) */
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
	icxViewCond *ivc_p,		/* Input Viewing Parameters for CAM */
	icxViewCond *ovc_p,		/* Output Viewing Parameters for CAM */
	int ivc_e,				/* Input Enumerated viewing condition */
	int ovc_e,				/* Output Enumerated viewing condition */
	icxGMappingIntent *pgmi,/* Perceptual gamut mapping intent */
	icxGMappingIntent *sgmi,/* Saturation gamut mapping intent */

	profxinf *xpi			/* Optional Profile creation extra data */
) {
	int isdisp;				/* nz if this is a display device, 0 if output */
	icmFile *wr_fp;
	icc *wr_icco;
	int npat;				/* Number of patches */
	co *tpat;				/* Patch input values */
	int i, j, rv = 0;
	icColorSpaceSignature devspace = icmSigDefaultData;	/* The device colorspace */
	int isAdditive = 0;		/* 0 if subtractive, 1 if additive colorspace */
	int isLab = 0;			/* 0 if input is XYZ, 1 if input is Lab */
	int wantLab = 0;		/* 0 if want is XYZ, 1 want is Lab */
	int isLut = 0;			/* 0 if shaper+ matrix, 1 if lut type */
	int isShTRC = 0;		/* 0 if separate gamma/shaper TRC, 1 if shared */
	int devchan = 0;		/* Number of device chanels */
	int a2binres = 0;		/* A2B input (device) table resolution */
	int a2bres = 0;			/* A2B clut resolution */
	int a2boutres = 0;		/* A2B output (PCS) table resolution */
	int b2ainres = 0;		/* B2A input (PCS) table resolution */
	int b2ares = 0;			/* B2A clut resolution */
	int b2aoutres = 0;		/* B2A output (device) table resolution */
	double cal[3][256];		/* Display calibration, cal[0][0] < 0.0 if invalid */

	{
		int ti;

		if ((ti = icg->find_kword(icg, 0, "DEVICE_CLASS")) < 0)
			error ("Input file doesn't contain keyword DEVICE_CLASS");

		if (strcmp(icg->t[0].kdata[ti],"DISPLAY") == 0)
			isdisp = 1;
		else
			isdisp = 0;
	}

	if (ptype == prof_clutLab) {		/* Lab lut */
		wantLab = 1;
		isLut = 1;
	} else if (ptype == prof_clutXYZ) {	/* XYZ lut */
		wantLab = 0;
		isLut = 1;
	} else {
		wantLab = 0;			/* gamma/shaper + matrix profile must be XYZ */
		isLut = 0;

		if (ptype == prof_gam1mat	
		 || ptype == prof_sha1mat) {
			isShTRC = 1;		/* Single curve */
		}
	}

	/* Figure out what sort of device colorspace it is */
	{
		int ti;

		if ((ti = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
			error("Input file doesn't contain keyword COLOR_REPS");

		if (strcmp(icg->t[0].kdata[ti],"CMYK_XYZ") == 0) {
			devspace = icSigCmykData;
			devchan = 4;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMYK_LAB") == 0) {
			devspace = icSigCmykData;
			devchan = 4;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMY_XYZ") == 0) {
			devspace = icSigCmyData;
			devchan = 3;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMY_LAB") == 0) {
			devspace = icSigCmyData;
			devchan = 3;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"RGB_XYZ") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"RGB_LAB") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 1;
			isAdditive = 1;
#ifdef IMP_MONO
		} else if (strcmp(icg->t[0].kdata[ti],"K_XYZ") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"K_LAB") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"W_XYZ") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"W_LAB") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 1;
			isAdditive = 1;
#endif /* IMP_MONO */
		} else 
			error("Output device input file has unhandled color representation '%s'",
			                                                     icg->t[0].kdata[ti]);
		/* Figure out some suitable table sizes */
		if (devchan >= 4) {
			if (iquality >= 3) {
		    	a2binres = 2048;
		    	a2bres = 23;
		    	a2boutres = 2048;
			} else if (iquality == 2) {
		    	a2binres = 2048;
		    	a2bres = 17;
		    	a2boutres = 2048;
			} else if (iquality == 1) {
		    	a2binres = 1024;
		    	a2bres = 9;
		    	a2boutres = 1024;
			} else {
		    	a2binres = 512;
		    	a2bres = 5;
		    	a2boutres = 512;
			}
		} else if (devchan >= 2) {
			if (iquality >= 3) {
		    	a2binres = 2048;
		    	a2bres = 45;
		    	a2boutres = 2048;
			} else if (iquality == 2) {
		    	a2binres = 2048;
		    	a2bres = 33;
		    	a2boutres = 2048;
			} else if (iquality == 1) {
		    	a2binres = 1024;
		    	a2bres = 17;
		    	a2boutres = 1024;
			} else {
		    	a2binres = 512;
		    	a2bres = 9;
		    	a2boutres = 512;
			}
		} else {	/* devchan == 1 */
			if (iquality >= 3) {
		    	a2binres = 256;
		    	a2bres = 2048;
		    	a2boutres = 256;
			} else if (iquality == 2) {
		    	a2binres = 256;
		    	a2bres = 1024;
		    	a2boutres = 256;
			} else if (iquality == 1) {
		    	a2binres = 256;
		    	a2bres = 512;
		    	a2boutres = 256;
			} else {
		    	a2binres = 256;
		    	a2bres = 256;
		    	a2boutres = 256;
			}
		}

		if (devchan >= 2) {
			if (oquality >= 3) {	/* Ultra High */
		    	b2ainres = 2048;
		    	b2ares = 45;
		    	b2aoutres = 2048;
			} else if (oquality == 2) {
		    	b2ainres = 2048;
		    	b2ares = 33;		/* High */
		    	b2aoutres = 2048;
			} else if (oquality == 1) {
		    	b2ainres = 1024;
		    	b2ares = 17;		/* Medium */
		    	b2aoutres = 1024;
			} else if (oquality >= 0) {
		    	b2ainres = 512;
		    	b2ares = 9;			/* Low */
		    	b2aoutres = 512;
			} else {				/* Special, Extremely low quality */
		    	b2ainres = 64;
		    	b2ares = 3;
		    	b2aoutres = 64;
			}
		} else {	/* devchan == 1 */
			if (oquality >= 3) {	/* Ultra High */
		    	b2ainres = 256;
		    	b2ares = 3;
		    	b2aoutres = 4096;
			} else if (oquality == 2) {
		    	b2ainres = 256;
		    	b2ares = 3;		/* High */
		    	b2aoutres = 2048;
			} else if (oquality == 1) {
		    	b2ainres = 256;
		    	b2ares = 3;		/* Medium */
		    	b2aoutres = 1024;
			} else if (oquality >= 0) {
		    	b2ainres = 256;
		    	b2ares = 3;			/* Low */
		    	b2aoutres = 512;
			} else {				/* Special, Extremely low quality */
		    	b2ainres = 64;
		    	b2ares = 3;
		    	b2aoutres = 64;
			}
		}
	}

	/* Read a display calibration set if we are given one */
	/* It will be in the second table with other type "CAL" */
	if (isdisp && icg->ntables >= 2 && icg->t[1].tt == tt_other && icg->t[0].oi != 1) {
		int ncal;
		int fi, ii, ri, gi, bi;
		
		if ((ncal = icg->t[1].nsets) <= 0)
			error ("No display calibration set data in .ti3");
	
		if (ncal != 256)
			error ("Expect 256 data sets in file '%s'");
	
		if ((fi = icg->find_kword(icg, 1, "DEVICE_CLASS")) < 0)
			error ("Display calibration doesn't contain keyword COLOR_REPS");
		if (strcmp(icg->t[1].kdata[fi],"DISPLAY") != 0)
			error ("Display Calibration doesn't have DEVICE_CLASS of DISPLAY");

		if ((ii = icg->find_field(icg, 1, "RGB_I")) < 0)
			error ("Display Calibration doesn't contain field RGB_I");
		if (icg->t[1].ftype[ii] != r_t)
			error ("Field RGB_R in display calibration is wrong type");
		if ((ri = icg->find_field(icg, 1, "RGB_R")) < 0)
			error ("Display Calibration doesn't contain field RGB_R");
		if (icg->t[1].ftype[ri] != r_t)
			error ("Field RGB_R in display calibration is wrong type");
		if ((gi = icg->find_field(icg, 1, "RGB_G")) < 0)
			error ("Display Calibration doesn't contain field RGB_G");
		if (icg->t[1].ftype[gi] != r_t)
			error ("Field RGB_G in display calibration is wrong type");
		if ((bi = icg->find_field(icg, 1, "RGB_B")) < 0)
			error ("Display Calibration doesn't contain field RGB_B");
		if (icg->t[1].ftype[bi] != r_t)
			error ("Field RGB_B in display calibration is wrong type");
		for (i = 0; i < ncal; i++) {
			cal[0][i] = *((double *)icg->t[1].fdata[i][ri]);
			cal[1][i] = *((double *)icg->t[1].fdata[i][gi]);
			cal[2][i] = *((double *)icg->t[1].fdata[i][bi]);
		}
	} else {
		cal[0][0] = -1.0;	/* Not used */
	}

	/* Open up the file for writing */
	if ((wr_fp = new_icmFileStd_name(file_name,"w")) == NULL)
		error("Write: Can't open file '%s'",file_name);

	if ((wr_icco = new_icc()) == NULL)
		error("Write: Creation of ICC object failed");

	/* Add all the tags required */

	/* The header: */
	{
		icmHeader *wh = wr_icco->header;

		/* Values that must be set before writing */
		if (isdisp)
			wh->deviceClass     = icSigDisplayClass;
		else
			wh->deviceClass     = icSigOutputClass;
    	wh->colorSpace      = devspace;					/* Device space it is */
		if (wantLab) 
	    	wh->pcs         = icSigLabData;
		else
    		wh->pcs         = icSigXYZData;				/* Must be XYZ for matrix based profile */
    	wh->renderingIntent = icRelativeColorimetric;	/* For want of something */

		/* Values that should be set before writing */
		if (xpi != NULL && xpi->manufacturer != 0L)
			wh->manufacturer = xpi->manufacturer;
		else
			wh->manufacturer = str2tag("????");

		if (xpi != NULL && xpi->model != 0L)
			wh->model = xpi->model;
		else
	    	wh->model        = str2tag("????");

		/* Values that may be set before writing */
		if (xpi != NULL && xpi->creator != 0L)
			wh->creator = xpi->creator;
	}
	/* Profile Description Tag: */
	{
		icmTextDescription *wo;
		char *dst, dstm[200];			/* description */

		if (xpi != NULL && xpi->profDesc != NULL)
			dst = xpi->profDesc;
		else {
			sprintf(dstm, "This is a Lut style %s - %s Output Profile",
			devspace == icSigCmykData ? "CMYK" :
			            devspace == icSigCmyData ? "CMY" : "RGB",
			wantLab ? "Lab" : "XYZ");
			dst = dstm;
		}

		if ((wo = (icmTextDescription *)wr_icco->add_tag(
		           wr_icco, icSigProfileDescriptionTag,	icSigTextDescriptionType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, dst);		/* Copy the string in */
	}
	/* Copyright Tag: */
	{
		icmText *wo;
		char *crt;

		if (xpi != NULL && xpi->copyright != NULL)
			crt = xpi->copyright;
		else
			crt = "Copyright, the creator of this profile";

		if ((wo = (icmText *)wr_icco->add_tag(
		           wr_icco, icSigCopyrightTag,	icSigTextType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(crt)+1; 	/* Allocated and used size of text, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->data, crt);		/* Copy the text in */
	}
	/* Device Manufacturers Description Tag: */
	if (xpi != NULL && xpi->deviceMfgDesc != NULL) {
		icmTextDescription *wo;
		char *dst = xpi->deviceMfgDesc;

		if ((wo = (icmTextDescription *)wr_icco->add_tag(
		           wr_icco, icSigDeviceMfgDescTag,	icSigTextDescriptionType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, dst);		/* Copy the string in */
	}
	/* Model Description Tag: */
	if (xpi != NULL && xpi->modelDesc != NULL) {
		icmTextDescription *wo;
		char *dst = xpi->modelDesc;

		if ((wo = (icmTextDescription *)wr_icco->add_tag(
		           wr_icco, icSigDeviceModelDescTag,	icSigTextDescriptionType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, dst);		/* Copy the string in */
	}
	/* White Point Tag: */
	{
		icmXYZArray *wo;
		/* Note that tag types icSigXYZType and icSigXYZArrayType are identical */
		if ((wo = (icmXYZArray *)wr_icco->add_tag(
		           wr_icco, icSigMediaWhitePointTag, icSigXYZArrayType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = 1;
		wo->allocate((icmBase *)wo);	/* Allocate space */
		wo->data[0] = icmD50;			/* Set a default value - D50 */
	}
	/* Black Point Tag: */
	{
		icmXYZArray *wo;
		if ((wo = (icmXYZArray *)wr_icco->add_tag(
		           wr_icco, icSigMediaBlackPointTag, icSigXYZArrayType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = 1;
		wo->allocate((icmBase *)wo);	/* Allocate space */
		wo->data[0] = icmBlack;			/* Set a default value - absolute black */
	}

	/* vcgt tag */
	if (cal[0][0] >= 0.0) {
		int j, i;
		icmVideoCardGamma *wo;
		wo = (icmVideoCardGamma *)wr_icco->add_tag(wr_icco,
		                           icSigVideoCardGammaTag, icSigVideoCardGammaType);
		if (wo == NULL)
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->tagType = icmVideoCardGammaTableType;
		wo->u.table.channels = 3;             /* rgb */
		wo->u.table.entryCount = 256;         /* full lut */
		wo->u.table.entrySize = 2;            /* 16 bits */
		wo->allocate((icmBase*)wo);
		for (j = 0; j < 3; j++)
			for (i = 0; i < 256; i++)
				((unsigned short*)wo->u.table.data)[256 * j + i]
				                                = (int)(cal[j][i] * 65535.0 + 0.5);
	}
	if (isLut) {		/* Lut type profile */

// ~~~~ Note. Assumption here that display Lut profile only has
// ~~~~ colorimetric AtoB0 is true for < ICC V2.4, but not 
// ~~~~ for >= ICC V2.4. How should this be dealt with ?

		/* 16 bit dev -> pcs lut: (A2B) */
		{
			icmLut *wo;
	
			if (isdisp) {	/* Only A2B0, no intent */
				if ((wo = (icmLut *)wr_icco->add_tag(
				           wr_icco, icSigAToB0Tag,	icSigLut16Type)) == NULL) 
					error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
			} else {
				/* Intent 1 = relative colorimetric */
				if ((wo = (icmLut *)wr_icco->add_tag(
				           wr_icco, icSigAToB1Tag,	icSigLut16Type)) == NULL) 
					error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
			}
		
			wo->inputChan = devchan;
			wo->outputChan = 3;
	    	wo->clutPoints = a2bres;
	    	wo->inputEnt   = a2binres;
	    	wo->outputEnt  = a2boutres;
			wo->allocate((icmBase *)wo);/* Allocate space */
	
			/* icxLuLut will set tables values */
		}

		if (!isdisp) {	/* Output device needs all the intents */

			/* 16 bit dev -> pcs lut - link intent 0 to intent 1 */
			{
				icmLut *wo;
				/* Intent 0 = perceptual */
				if ((wo = (icmLut *)wr_icco->link_tag(
				           wr_icco, icSigAToB0Tag,	icSigAToB1Tag)) == NULL) 
					error("link_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
			}
	
			/* 16 dev -> pcs bit lut - link intent 2 to intent 1 */
			{
				icmLut *wo;
				/* Intent 2 = saturation */
				if ((wo = (icmLut *)wr_icco->link_tag(
				           wr_icco, icSigAToB2Tag,	icSigAToB1Tag)) == NULL) 
					error("link_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
			}
		}

		/* 16 bit pcs -> dev lut: (B2A) */
		{
			icmLut *wo;

			if (isdisp) {	/* Only B2A0, no intent */
				if ((wo = (icmLut *)wr_icco->add_tag(
				           wr_icco, icSigBToA0Tag,	icSigLut16Type)) == NULL) 
					error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

			} else {

				/* Intent 1 = relative colorimetric */
				if ((wo = (icmLut *)wr_icco->add_tag(
				           wr_icco, icSigBToA1Tag,	icSigLut16Type)) == NULL) 
					error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
			}

			wo->inputChan = 3;
			wo->outputChan = devchan;
		    wo->clutPoints = b2ares;
		    wo->inputEnt   = b2ainres;
		    wo->outputEnt  = b2aoutres;
			wo->allocate((icmBase *)wo);/* Allocate space */
		}

		if (!isdisp) {	/* Output needs all the intents */

			if (ipname == NULL) {		/* No gamut mapping */
				icmLut *wo;

				/* link intent 0 = perceptual to intent 1 = colorimetric */
				if ((wo = (icmLut *)wr_icco->link_tag(
				           wr_icco, icSigBToA0Tag,	icSigBToA1Tag)) == NULL) 
					error("link_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

				/* Intent 2 = saturation */
				/* link intent 2 = saturation to intent 1 = colorimetric */
				if ((wo = (icmLut *)wr_icco->link_tag(
				           wr_icco, icSigBToA2Tag,	icSigBToA1Tag)) == NULL) 
					error("link_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

			} else {				/* We have gamut mapping */
				icmLut *wo;

				/* Intent 0 = perceptual */
				if ((wo = (icmLut *)wr_icco->add_tag(
				           wr_icco, icSigBToA0Tag,	icSigLut16Type)) == NULL) 
					error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

				wo->inputChan = 3;
				wo->outputChan = devchan;
			    wo->inputEnt   = b2ainres;
			    wo->clutPoints = b2ares;
			    wo->outputEnt  = b2aoutres;
				wo->allocate((icmBase *)wo);/* Allocate space */

				if (sepsat == 0) {		/* No separate gamut mapping for saturation */
					/* link intent 2 = saturation to intent 0 = perceptual */
					if ((wo = (icmLut *)wr_icco->link_tag(
					           wr_icco, icSigBToA2Tag,	icSigBToA0Tag)) == NULL) 
						error("link_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
				} else {
					/* Intent 2 = saturation */
					if ((wo = (icmLut *)wr_icco->add_tag(
					           wr_icco, icSigBToA2Tag,	icSigLut16Type)) == NULL) 
						error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

					wo->inputChan = 3;
					wo->outputChan = devchan;
				    wo->inputEnt   = b2ainres;
				    wo->clutPoints = b2ares;
				    wo->outputEnt  = b2aoutres;
					wo->allocate((icmBase *)wo);/* Allocate space */
				}
			}

			/* 16 bit pcs -> gamut lut: */
			{
				icmLut *wo;

				if ((wo = (icmLut *)wr_icco->add_tag(
				           wr_icco, icSigGamutTag,	icSigLut16Type)) == NULL) 
					error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

				wo->inputChan = 3;
				wo->outputChan = 1;
				wo->inputEnt = 256;
				wo->clutPoints = 16;
				wo->outputEnt = 256;
				wo->allocate((icmBase *)wo);/* Allocate space */

				/* ~~~ This isn't being set at the moment. Needs to be fixed. */
			}
		}

	} else {	/* shaper + matrix type */

		/* Red, Green and Blue Colorant Tags: */
		{
			icmXYZArray *wor, *wog, *wob;
			if ((wor = (icmXYZArray *)wr_icco->add_tag(
			           wr_icco, icSigRedColorantTag, icSigXYZArrayType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);
			if ((wog = (icmXYZArray *)wr_icco->add_tag(
			           wr_icco, icSigGreenColorantTag, icSigXYZArrayType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);
			if ((wob = (icmXYZArray *)wr_icco->add_tag(
			           wr_icco, icSigBlueColorantTag, icSigXYZArrayType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);

			wor->size = wog->size = wob->size = 1;
			wor->allocate((icmBase *)wor);	/* Allocate space */
			wog->allocate((icmBase *)wog);
			wob->allocate((icmBase *)wob);

			/* Setup some sane dummy values */
			/* icxMatrix will override these later */
			wor->data[0].X = 1.0; wor->data[0].Y = 0.0; wor->data[0].Z = 0.0;
			wog->data[0].X = 0.0; wog->data[0].Y = 1.0; wog->data[0].Z = 0.0;
			wob->data[0].X = 0.0; wob->data[0].Y = 0.0; wob->data[0].Z = 1.0;
		}

		/* Red, Green and Blue Tone Reproduction Curve Tags: */
		{
			icmCurve *wor, *wog, *wob;
			int i;
			if ((wor = (icmCurve *)wr_icco->add_tag(
			           wr_icco, icSigRedTRCTag, icSigCurveType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);

			if (isShTRC) {	/* Make all TRCs shared */
				if ((wog = (icmCurve *)wr_icco->link_tag(
				           wr_icco, icSigGreenTRCTag, icSigRedTRCTag)) == NULL) 
					error("link_tag failed: %d, %s",rv,wr_icco->err);
				if ((wob = (icmCurve *)wr_icco->link_tag(
				           wr_icco, icSigBlueTRCTag, icSigRedTRCTag)) == NULL) 
					error("link_tag failed: %d, %s",rv,wr_icco->err);

			} else {		/* Else individual */
				if ((wog = (icmCurve *)wr_icco->add_tag(
				           wr_icco, icSigGreenTRCTag, icSigCurveType)) == NULL) 
					error("add_tag failed: %d, %s",rv,wr_icco->err);
				if ((wob = (icmCurve *)wr_icco->add_tag(
				           wr_icco, icSigBlueTRCTag, icSigCurveType)) == NULL) 
					error("add_tag failed: %d, %s",rv,wr_icco->err);
			}
	
			if (ptype == prof_shamat || ptype == prof_sha1mat) {	/* Shaper */
				wor->flag = wog->flag = wob->flag = icmCurveSpec; 
				wor->size = wog->size = wob->size = 256;			/* Number of entries */
			} else {						/* Gamma */
				wor->flag = wog->flag = wob->flag = icmCurveGamma;
				wor->size = wog->size = wob->size = 1;				/* Must be 1 for gamma */
			}
			wor->allocate((icmBase *)wor);	/* Allocate space */
			wog->allocate((icmBase *)wog);
			wob->allocate((icmBase *)wob);

			/* icxMatrix will set curve values */
		}
	}

	if ((npat = icg->t[0].nsets) <= 0)
		error("No sets of data");

	if (verb) {
		fprintf(verbo,"No of test patches = %d\n",npat);
	}

	/* Allocate arrays to hold test patch input and output values */
	if ((tpat = (co *)malloc(sizeof(co) * npat)) == NULL)
		error("Malloc failed - tpat[]");

	/* Read in the CGATs fields */
	{
		int ti, ii, ci, mi, yi, ki;
		int Xi, Yi, Zi;

		/* Read the ink limit */
		if (ink != NULL && (ii = icg->find_kword(icg, 0, "TOTAL_INK_LIMIT")) >= 0) {
			double ilimit = -1;
			ilimit = atof(icg->t[0].kdata[ii]);
			if (ilimit > 1e-4 && ilimit <= 400.0 && ink->tlimit < 0.0) {
				ink->tlimit = ilimit/100.0;	/* Set requested ink limit */
			}
		}
		/* Should targen/.ti3 file allow for BLACK_INK_LIMIT ?? */
		
		/* A problem here is that if the .ti3 is corrupted, then */
		/* often this results in the field type being "wrong", */
		/* rather than a more inteligable message. */

		if (devspace == icSigGrayData) {

			if (isAdditive) {
				if ((ci = icg->find_field(icg, 0, "GRAY_W")) < 0)
					error("Input file doesn't contain field GRAY_W");
				if (icg->t[0].ftype[ci] != r_t)
					error("Field GRAY_W is wrong type - corrupted file ?");
			} else {
				if ((ci = icg->find_field(icg, 0, "GRAY_K")) < 0)
					error("Input file doesn't contain field GRAY_K");
				if (icg->t[0].ftype[ci] != r_t)
					error("Field GRAY_K is wrong type - corrupted file ?");
			}
			mi = yi = ki = ci;

		} else if (devspace == icSigRgbData) {

			if ((ci = icg->find_field(icg, 0, "RGB_R")) < 0)
				error("Input file doesn't contain field RGB_R");
			if (icg->t[0].ftype[ci] != r_t)
				error("Field RGB_R is wrong type - corrupted file ?");
			if ((mi = icg->find_field(icg, 0, "RGB_G")) < 0)
				error("Input file doesn't contain field RGB_G");
			if (icg->t[0].ftype[mi] != r_t)
				error("Field RGB_G is wrong type - corrupted file ?");
			if ((yi = icg->find_field(icg, 0, "RGB_B")) < 0)
				error("Input file doesn't contain field RGB_B");
			if (icg->t[0].ftype[yi] != r_t)
				error("Field RGB_B is wrong type - corrupted file ?");
			ki = yi;

		} else if (devspace == icSigCmyData) {

			if ((ci = icg->find_field(icg, 0, "CMY_C")) < 0)
				error("Input file doesn't contain field CMY_C");
			if (icg->t[0].ftype[ci] != r_t)
				error("Field CMY_C is wrong type - corrupted file ?");
			if ((mi = icg->find_field(icg, 0, "CMY_M")) < 0)
				error("Input file doesn't contain field CMY_M");
			if (icg->t[0].ftype[mi] != r_t)
				error("Field CMY_M is wrong type - corrupted file ?");
			if ((yi = icg->find_field(icg, 0, "CMY_Y")) < 0)
				error("Input file doesn't contain field CMY_Y");
			if (icg->t[0].ftype[yi] != r_t)
				error("Field CMY_Y is wrong type - corrupted file ?");
			ki = yi;

		} else {	/* Assume CMYK */

			if ((ci = icg->find_field(icg, 0, "CMYK_C")) < 0)
				error("Input file doesn't contain field CMYK_C");
			if (icg->t[0].ftype[ci] != r_t)
				error("Field CMYK_C is wrong type - corrupted file ?",icg->t[0].ftype[ci],r_t);
			if ((mi = icg->find_field(icg, 0, "CMYK_M")) < 0)
				error("Input file doesn't contain field CMYK_M");
			if (icg->t[0].ftype[mi] != r_t)
				error("Field CMYK_M is wrong type - corrupted file ?");
			if ((yi = icg->find_field(icg, 0, "CMYK_Y")) < 0)
				error("Input file doesn't contain field CMYK_Y");
			if (icg->t[0].ftype[yi] != r_t)
				error("Field CMYK_Y is wrong type - corrupted file ?");
			if ((ki = icg->find_field(icg, 0, "CMYK_K")) < 0)
				error("Input file doesn't contain field CMYK_K");
			if (icg->t[0].ftype[ki] != r_t)
				error("Field CMYK_K is wrong type - corrupted file ?");
		}

		if (spec == 0) { 		/* Using instrument tristimulous value */

			if (isLab) {		/* Expect Lab */
				if ((Xi = icg->find_field(icg, 0, "LAB_L")) < 0)
					error("Input file doesn't contain field LAB_L");
				if (icg->t[0].ftype[Xi] != r_t)
					error("Field LAB_L is wrong type - corrupted file ?");
				if ((Yi = icg->find_field(icg, 0, "LAB_A")) < 0)
					error("Input file doesn't contain field LAB_A");
				if (icg->t[0].ftype[Yi] != r_t)
					error("Field LAB_A is wrong type - corrupted file ?");
				if ((Zi = icg->find_field(icg, 0, "LAB_B")) < 0)
					error("Input file doesn't contain field LAB_B");
				if (icg->t[0].ftype[Zi] != r_t)
					error("Field LAB_B is wrong type - corrupted file ?");

			} else { 		/* Expect XYZ */
				if ((Xi = icg->find_field(icg, 0, "XYZ_X")) < 0)
					error("Input file doesn't contain field XYZ_X");
				if (icg->t[0].ftype[Xi] != r_t)
					error("Field XYZ_X is wrong type - corrupted file ?");
				if ((Yi = icg->find_field(icg, 0, "XYZ_Y")) < 0)
					error("Input file doesn't contain field XYZ_Y");
				if (icg->t[0].ftype[Yi] != r_t)
					error("Field XYZ_Y is wrong type - corrupted file ?");
				if ((Zi = icg->find_field(icg, 0, "XYZ_Z")) < 0)
					error("Input file doesn't contain field XYZ_Z");
				if (icg->t[0].ftype[Zi] != r_t)
					error("Field XYZ_Z is wrong type - corrupted file ?");
			}

			for (i = 0; i < npat; i++) {
				tpat[i].p[0] = *((double *)icg->t[0].fdata[i][ci]) / 100.0;
				tpat[i].p[1] = *((double *)icg->t[0].fdata[i][mi]) / 100.0;
				tpat[i].p[2] = *((double *)icg->t[0].fdata[i][yi]) / 100.0;
				tpat[i].p[3] = *((double *)icg->t[0].fdata[i][ki]) / 100.0;
				if (tpat[i].p[0] > 1.0
				 || tpat[i].p[1] > 1.0
				 || tpat[i].p[2] > 1.0
				 || tpat[i].p[3] > 1.0) {
					double bgst = 0.0;
					int j, bj;
					for (j = 0; j < 4; j++) {
						if (tpat[i].p[j] > bgst) {
							bgst = tpat[i].p[j];
							bj = j;
						}
					}
					error("Device value field value exceeds 100.0 (%d:%d:%f) !",i,bj,bgst * 100.0);
				}
				tpat[i].v[0] = *((double *)icg->t[0].fdata[i][Xi]);
				tpat[i].v[1] = *((double *)icg->t[0].fdata[i][Yi]);
				tpat[i].v[2] = *((double *)icg->t[0].fdata[i][Zi]);
				if (!isLab) {
					tpat[i].v[0] /= 100.0;		/* Normalise XYZ to range 0.0 - 1.0 */
					tpat[i].v[1] /= 100.0;
					tpat[i].v[2] /= 100.0;
				}
				if (!isLab && wantLab) { /* Convert test patch result XYZ to PCS (D50 Lab) */
					icmXYZ2Lab(&icmD50, tpat[i].v, tpat[i].v);
				} else if (isLab && !wantLab) {
					icmLab2XYZ(&icmD50, tpat[i].v, tpat[i].v);
				}
			}

		} else { 		/* Using spectral data */
			int ii;
			xspect sp;
			char buf[100];
			int  spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
			xsp2cie *sp2cie;	/* Spectral conversion object */

			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_BANDS")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_BANDS");
			sp.spec_n = atoi(icg->t[0].kdata[ii]);
			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_START_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_START_NM");
			sp.spec_wl_short = atof(icg->t[0].kdata[ii]);
			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_END_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_END_NM");
			sp.spec_wl_long = atof(icg->t[0].kdata[ii]);
			sp.norm = 100.0;

			/* Find the fields for spectral values */
			for (j = 0; j < sp.spec_n; j++) {
				int nm;
		
				/* Compute nearest integer wavelength */
				nm = (int)(sp.spec_wl_short + ((double)j/(sp.spec_n-1.0))
				            * (sp.spec_wl_long - sp.spec_wl_short) + 0.5);
				
				sprintf(buf,"SPEC_%03d",nm);

				if ((spi[j] = icg->find_field(icg, 0, buf)) < 0)
					error("Input file doesn't contain field %s",buf);
			}

			if (isdisp) {
				illum = icxIT_none;		/* Displays are assumed to be self luminous */
				cust_illum = NULL;
			}

			/* Create a spectral conversion object */
			if ((sp2cie = new_xsp2cie(illum, cust_illum, observ, NULL,
			                          wantLab ? icSigLabData : icSigXYZData)) == NULL)
				error("Creation of spectral conversion object failed");

			/* If Fluorescent Whitening Agent compensation is enabled */
			/* (isdisp must be 0) */
			if (fwacomp) {
				double nw = 0.0;		/* Number of media white patches */
				xspect mwsp;			/* Medium spectrum */
				instType itype;			/* Spectral instrument type */
				xspect *insp;			/* Instrument illuminant */

				mwsp = sp;		/* Struct copy */

	 			if ((ti = icg->find_kword(icg, 0, "TARGET_INSTRUMENT")) < 0)
					error ("Can't find target instrument needed for FWA compensation");

				if ((itype = inst_enum(icg->t[0].kdata[ti])) == instUnknown)
					error ("Unrecognised target instrument '%s'", icg->t[0].kdata[ti]);

				if ((insp = inst_illuminant(itype)) == NULL)
					error ("Instrument doesn't have an FWA illuminent");

				/* Find the media white spectral reflectance */
				for (j = 0; j < mwsp.spec_n; j++)
					mwsp.spec[j] = 0.0;

				/* Compute the mean of all the media white patches */
				for (i = 0; i < npat; i++) {
	
					if (*((double *)icg->t[0].fdata[i][ci]) < 1e-4
					 && *((double *)icg->t[0].fdata[i][mi]) < 1e-4
					 && *((double *)icg->t[0].fdata[i][yi]) < 1e-4
					 && *((double *)icg->t[0].fdata[i][ki]) < 1e-4) {
	
						/* Read the spectral values for this patch */
						for (j = 0; j < mwsp.spec_n; j++) {
							mwsp.spec[j] += *((double *)icg->t[0].fdata[i][spi[j]]);
						}
						nw++;
					}
				}
				if (nw == 0.0)
					error ("Can't find a media white patch to init FWA");

				for (j = 0; j < mwsp.spec_n; j++)
					mwsp.spec[j] /= nw;	/* Compute average */

				if (sp2cie->set_fwa(sp2cie, insp, &mwsp)) 
					error ("Set FWA on sp2cie failed");
			}

			for (i = 0; i < npat; i++) {

				tpat[i].p[0] = *((double *)icg->t[0].fdata[i][ci]) / 100.0;
				tpat[i].p[1] = *((double *)icg->t[0].fdata[i][mi]) / 100.0;
				tpat[i].p[2] = *((double *)icg->t[0].fdata[i][yi]) / 100.0;
				tpat[i].p[3] = *((double *)icg->t[0].fdata[i][ki]) / 100.0;

				if (tpat[i].p[0] > 1.0
				 || tpat[i].p[1] > 1.0
				 || tpat[i].p[2] > 1.0
				 || tpat[i].p[3] > 1.0) {
					error("Device value field value exceeds 100.0 !");
				}

				/* Read the spectral values for this patch */
				for (j = 0; j < sp.spec_n; j++) {
					sp.spec[j] = *((double *)icg->t[0].fdata[i][spi[j]]);
				}

				/* Convert it to CIE space */
				sp2cie->convert(sp2cie, tpat[i].v, &sp);

			}

			sp2cie->del(sp2cie);		/* Done with this */

		}
	}	/* End of reading in CGATs file */

	if (isLut) {
		xicc *wr_xicc;			/* extention object */
		icxLuBase *AtoB;		/* AtoB ixcLu */

		/* Create A2B clut */
		{
			int flags = 0;

			/* Wrap with an expanded icc */
			if ((wr_xicc = new_xicc(wr_icco)) == NULL)
				error("Creation of xicc failed");
		
			if (noiluts)
				flags |= ICX_NO_IN_LUTS;

			if (nooluts)
				flags |= ICX_NO_OUT_LUTS;

			if (verb)
				flags |= ICX_VERBOSE;

			/* Setup Device -> PCS conversion (Fwd) object from scattered data. */
			if ((AtoB = wr_xicc->set_luobj(
			               wr_xicc, icmFwd, isdisp ? icmDefaultIntent : icRelativeColorimetric,
			               icmLuOrdNorm,
#ifdef USE_EXTRA_FITTING
			               (iquality >= 1 ? ICX_EXTRA_FIT : 0) |	/* Extra if medium or higher */
#endif
			               flags | ICX_SET_WHITE | ICX_SET_BLACK, 		/* Flags */
			               npat, tpat, smooth, avgdev, NULL, ink, iquality)) == NULL)
				error("%d, %s",wr_xicc->errc, wr_xicc->err);

			AtoB->del(AtoB);		/* Done with lookup */
		}

		/* Create B2A clut */
		{
			icmFile *fp = NULL;
			icc *src_icco = NULL;
			xicc *src_xicc = NULL;	/* Source profile */
			icxViewCond ivc;	/* Input Viewing Condition for CAM */
			icxViewCond ovc;	/* Output Viewing Condition for CAM */
			icmLut *wo[3];
			icmFile *abs_fp;	/* Abstract profile transform: */
			icc *abs_icc;
			xicc *abs_xicc;

			out_b2a_callback cx;

			if (verb)
				printf("Creating B to A tables\n");

			if (ipname != NULL) {		/* There is a source profile to determine gamut mapping */

				/* Open up the profile for reading */
				if ((fp = new_icmFileStd_name(ipname,"r")) == NULL)
					error ("Can't open file '%s'",ipname);

				if ((src_icco = new_icc()) == NULL)
					error ("Creation of ICC object failed");

				if ((rv = src_icco->read(src_icco,fp,0)) != 0)
					error ("%d, %s",rv,src_icco->err);

				/* Wrap with an expanded icc */
				if ((src_xicc = new_xicc(src_icco)) == NULL)
					error ("Creation of src_xicc failed");
			}

			/* Figure out the final src & dst viewing conditions */
			for (i = 0; i < 2; i++) {
				xicc *x;
				icxViewCond *v, *vc;
				int es;
		
				if (i == 0) {			/* Input */
					v = ivc_p;			/* Override parameters */
					es = ivc_e;
					vc = &ivc;			/* Target parameters */
					if (src_xicc == NULL)
						continue;		/* Source viewing conditions won't be used */
					x = src_xicc;
				} else {				/* Output */
					v = ovc_p;			/* Override parameters */
					es = ovc_e;
					vc = &ovc;			/* Target parameters */
					x = wr_xicc;
				}
				
				/* Set the default */
				xicc_enum_viewcond(x, vc, -1, 0);
		
				/* Override the viewing conditions */
				if (es >= 0)
					if (xicc_enum_viewcond(x, vc, es, 0))
						error ("%d, %s",x->errc, x->err);
				if (v->Ev >= 0)
					vc->Ev = v->Ev;
				if (v->Wxyz[0] >= 0.0 && v->Wxyz[1] > 0.0 && v->Wxyz[2] >= 0.0) {
					/* Normalise XYZ to current media white */
					vc->Wxyz[0] = v->Wxyz[0]/v->Wxyz[1] * vc->Wxyz[1];
					vc->Wxyz[2] = v->Wxyz[2]/v->Wxyz[1] * vc->Wxyz[1];
				} 
				if (v->Wxyz[0] >= 0.0 && v->Wxyz[1] >= 0.0 && v->Wxyz[2] < 0.0) {
					/* Convert Yxy to XYZ */
					double x = v->Wxyz[0];
					double y = v->Wxyz[1];	/* If Y == 1.0, then X+Y+Z = 1/y */
					double z = 1.0 - x - y;
					vc->Wxyz[0] = x/y * vc->Wxyz[1];
					vc->Wxyz[2] = z/y * vc->Wxyz[1];
				}
				if (v->La >= 0.0)
					vc->La = v->La;
				if (v->Yb >= 0.0)
					vc->Yb = v->Yb;
				if (v->Yf >= 0.0)
					vc->Yf = vc->Yf;
				if (v->Fxyz[0] >= 0.0 && v->Fxyz[1] > 0.0 && v->Fxyz[2] >= 0.0) {
					/* Normalise XYZ to current media white */
					vc->Fxyz[0] = v->Fxyz[0]/v->Fxyz[1] * vc->Fxyz[1];
					vc->Fxyz[2] = v->Fxyz[2]/v->Fxyz[1] * vc->Fxyz[1];
				}
				if (v->Fxyz[0] >= 0.0 && v->Fxyz[1] >= 0.0 && v->Fxyz[2] < 0.0) {
					/* Convert Yxy to XYZ */
					double x = v->Fxyz[0];
					double y = v->Fxyz[1];	/* If Y == 1.0, then X+Y+Z = 1/y */
					double z = 1.0 - x - y;
					vc->Fxyz[0] = x/y * vc->Fxyz[1];
					vc->Fxyz[2] = z/y * vc->Fxyz[1];
				}
			}

			/* Get a suitable forward conversion object to invert */
			/* By creating a separate one to the one created using scattered data, */
			/* we ge the chance to set ICX_CAM_CLIP. */
			{
				int flags = 0;
	
				if (verb)
					flags |= ICX_VERBOSE;
	
				flags |= ICX_CLIP_NEAREST;		/* Not vector clip */
#ifdef USE_CAM_CLIP_OPT
				flags |= ICX_CAM_CLIP;			/* Clip in CAM Jab space rather than Lab */
#endif
	
				if ((AtoB = wr_xicc->get_luobj(wr_xicc, flags, icmFwd,
				                  isdisp ? icmDefaultIntent : icRelativeColorimetric, icSigLabData,
                                  icmLuOrdNorm, &ovc, ink)) == NULL)
					error ("%d, %s",wr_xicc->errc, wr_xicc->err);
			}

			/* setup context ready for B2A table setting */
			cx.verb = verb;
#ifdef NO_B2A_PCS_CURVES
			cx.noPCScurves = 1;
#else
			cx.noPCScurves = 0;
#endif
			cx.pcsspace = wantLab ? icSigLabData : icSigXYZData;
			cx.devspace = devspace;
			cx.x = (icxLuLut *)AtoB;		/* A2B icxLuLut created from scattered data */

			cx.ix = NULL;		/* PCS to CAM conversion */
			cx.ox = NULL;		/* CAM to PCS conversion */
			cx.pmap = NULL;		/* perceptual gamut map */
			cx.smap = NULL;		/* Saturation gamut map */

			/* Open up the abstract profile if supplied, and setup luo */
			if (absname != NULL) {
				if ((abs_fp = new_icmFileStd_name(absname,"r")) == NULL)
					error ("Can't open abstract profile file '%s'",absname);
				
				if ((abs_icc = new_icc()) == NULL)
					error ("Creation of Abstract profile ICC object failed");
		
				/* Read header etc. */
				if ((rv = abs_icc->read(abs_icc,abs_fp,0)) != 0)
					error ("%d, %s",rv,abs_icc->err);
				abs_icc->header;
		
				if (abs_icc->header->deviceClass != icSigAbstractClass)
					error("Abstract profile isn't an abstract profile");
		
				/* Take intended abstract intent from profile itself */
				if ((cx.abs_intent = abs_icc->header->renderingIntent) != icAbsoluteColorimetric)
					cx.abs_intent = icRelativeColorimetric;
		
				/* Wrap with an expanded icc */
				if ((abs_xicc = new_xicc(abs_icc)) == NULL)
					error ("Creation of abstract profile xicc failed");

				/* The abstract profile intent is assumed to determine how it gets applied. */
				/* Make abstract PCS XYZ if icAbsoluteColorimetric is needed. */
				if ((cx.abs_luo = abs_xicc->get_luobj(abs_xicc, 0, icmFwd, cx.abs_intent,
			        (cx.pcsspace == icSigLabData && cx.abs_intent == icRelativeColorimetric)
					             ? icSigLabData : icSigXYZData,
					icmLuOrdNorm, NULL, NULL)) == NULL)
					error ("%d, %s",abs_icc->errc, abs_icc->err);
			} else {
				cx.abs_luo = NULL;
			}
			
			if (isdisp) {	/* Only B2A0, no intent */
				if ((wo[0] = (icmLut *)wr_icco->read_tag(
				           wr_icco, icSigBToA0Tag)) == NULL) 
					error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
				cx.ntables = 1;

			} else {		/* Output profile, 3 intent tables */
				/* Intent 1 = relative colorimetric */
				if ((wo[0] = (icmLut *)wr_icco->read_tag(
				           wr_icco, icSigBToA1Tag)) == NULL) 
					error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
				cx.ntables = 1;

				if (src_xicc) {		/* Creating separate perceptual and Saturation tables */
					icRenderingIntent intent;
					gamut *csgam, *igam = NULL, *ogam;
					double gres;			/* Gamut surface feature resolution */
					int    mapres;			/* Mapping rspl resolution */
					double mn[3] = { 0.0, -150.0, -150.0 };	/* Set Jab range */
					double mx[3] = { 100.0, 150.0, 150.0 };
			
					if (verb)
						printf("Creating Gamut Mapping\n");
			
					if (oquality == 3) {	/* Ultra High */
			  	 		gres = 8.0;
			  	 		mapres = 43;
					} else if (oquality == 2) {	/* High */
			  	 		gres = 9.0;
			  	 		mapres = 33;
					} else if (oquality == 1) {	/* Medium */
			  	 		gres = 10.0;
			  	 		mapres = 25;
					} else if (oquality == 0) {	/* Low quality */
			  	 		gres = 12.0;
			  	 		mapres = 17;
					} else {					/* Extremely low */
			  	 		gres = 15.0;
			  	 		mapres = 9;
					}

					if (pgmi->usecas != sgmi->usecas)
						error("Can't handle precept and sat table intents with different CAM spaces");
					intent = icAbsoluteColorimetric;
					if (pgmi->usecas > 0)
						intent = icxAppearance;
					if (pgmi->usecas == 2) {
						double mxw;
						intent = icxAbsAppearance;

						/* Make absolute common white point average between the two */
						ivc.Wxyz[0] = 0.5 * (ivc.Wxyz[0] + ovc.Wxyz[0]);
						ivc.Wxyz[1] = 0.5 * (ivc.Wxyz[1] + ovc.Wxyz[1]);
						ivc.Wxyz[2] = 0.5 * (ivc.Wxyz[2] + ovc.Wxyz[2]);
		
						/* And scale it Y to be equal to 1.0 */
						mxw = 1.0/ivc.Wxyz[1];
						ivc.Wxyz[0] *= mxw;
						ivc.Wxyz[1] *= mxw;
						ivc.Wxyz[2] *= mxw;
		
						/* set output view conditions the same as the input */
						ovc = ivc;		/* Structure copy */
					}

					/* Get lookup object simply for fwd_relpcs_outpcs() */
					if ((cx.ix = src_xicc->get_luobj(src_xicc, 0, icmFwd, intent, icSigLabData,
	                                  icmLuOrdNorm, &ivc, NULL)) == NULL)
						error ("%d, %s",src_xicc->errc, src_xicc->err);

					/* Get lookup object simply for bwd_outpcs_relpcs() */
					if ((cx.ox = wr_xicc->get_luobj(wr_xicc, 0, icmFwd, intent, icSigLabData,
	                                  icmLuOrdNorm, &ovc, NULL)) == NULL)
						error ("%d, %s",wr_xicc->errc, wr_xicc->err);

					/* Create the source colorspace gamut surface */
					if (verb)
						printf(" Finding Source Colorspace Gamut\n");
			
					if ((csgam = cx.ix->get_gamut(cx.ix, gres)) == NULL)
						error ("%d, %s",src_xicc->errc, src_xicc->err);
			
					/* Read image source gamut if provided */
					if (sgname != NULL) {	/* Optional source gamut - ie. from an images */
			
						if (verb)
							printf(" Loading Image Source Gamut '%s'\n",sgname);
			
						igam = new_gamut(gres);
			
						if (igam->read_gam(igam, sgname))
							error("Reading source gamut '%s' failed",sgname);
					}

					/* Creat the destination gamut surface */
					if (verb)
						printf(" Finding Destination Gamut\n");
			
					if ((ogam = cx.ox->get_gamut(cx.ox, gres)) == NULL)
						error ("%d, %s",wr_xicc->errc, wr_xicc->err);
			
					if (verb)
						printf(" Creating Gamut match\n");

					/* setup perceptual gamut mapping */
					cx.pmap = new_gammap(verb, csgam, igam, ogam, pgmi->greymf,
		                    pgmi->glumwcpf, pgmi->glumwexf, pgmi->glumbcpf, pgmi->glumbexf,
		                    pgmi->glumknf,
		                    pgmi->gamcpf, pgmi->gamexf, pgmi->gamknf,
		                    pgmi->gampwf, pgmi->gamswf, pgmi->satenh,
		                    mapres, mn, mx);
					if (cx.pmap == NULL)
						error ("Failed to make perceptual gamut map transform");

					/* Intent 0 = perceptual */
					if ((wo[1] = (icmLut *)wr_icco->read_tag(
					           wr_icco, icSigBToA0Tag)) == NULL) 
						error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
					cx.ntables = 2;

					if (sepsat) {
						/* setup saturation gamut mapping */
						cx.smap = new_gammap(verb, csgam, igam, ogam, sgmi->greymf,
			                    sgmi->glumwcpf, sgmi->glumwexf, sgmi->glumbcpf, sgmi->glumbexf,
			                    sgmi->glumknf,
			                    sgmi->gamcpf, sgmi->gamexf, sgmi->gamknf,
			                    sgmi->gampwf, sgmi->gamswf, sgmi->satenh,
			                    mapres, mn, mx);
						if (cx.smap == NULL)
							error ("Failed to make saturation gamut map transform");

						/* Intent 2 = saturation */
						if ((wo[2] = (icmLut *)wr_icco->read_tag(
						           wr_icco, icSigBToA2Tag)) == NULL) 
							error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
						cx.ntables = 3;
					}
					csgam->del(csgam);
					csgam = NULL;
					if (igam != NULL) {
						igam->del(igam);
						igam = NULL;
					}
					ogam->del(ogam);
					ogam = NULL;
				}
			}

			cx.ochan = wo[0]->outputChan;

			/* we setup an exact inverse, colorimetric style, plus gamut mapping */
			/* for perceptual and saturation intents */
			/* Use helper function to do the hard work. */

			if (cx.verb) {
				cx.count = 0;
				cx.last = -1;
				for (cx.total = 1, i = 0; i < wo[0]->inputChan; i++, cx.total *= wo[0]->clutPoints)
					; 
				printf(" 0%%"); fflush(stdout);
			}

			if (icmSetMultiLutTables(
			        cx.ntables,
			        wo,
					&cx,						/* Context */
					cx.pcsspace,				/* Input color space */
					devspace, 					/* Output color space */
					out_b2a_input,				/* Input transform PCS->PCS' */
					NULL, NULL,					/* Use default Lab' range */
					out_b2a_clut,				/* Lab' -> Device' transfer function */
					NULL, NULL,					/* Use default Device' range */
					out_b2a_output) != 0)		/* Output transfer function, Device'->Device */
				error("Setting 16 bit PCS->Device Lut failed: %d, %s",wr_icco->errc,wr_icco->err);
			if (cx.verb) {
				printf("\n");
			}

			if (cx.abs_luo != NULL) {		/* Free up abstract transform */
				cx.abs_luo->del(cx.abs_luo);
				abs_xicc->del(abs_xicc);
				abs_icc->del(abs_icc);
				abs_fp->del(abs_fp);
			}

			if (cx.ix != NULL)
				cx.ix->del(cx.ix), cx.ix = NULL;
			if (cx.ox != NULL)
				cx.ox->del(cx.ox), cx.ox = NULL;

			if (src_xicc != NULL)
				src_xicc->del(src_xicc), src_xicc = NULL;
			if (src_icco != NULL)
				src_icco->del(src_icco), src_icco = NULL;
			if (fp != NULL)
				fp->del(fp), fp = NULL;

			if (verb)
				printf("Done B to A tables\n");
		}

#ifdef USE_CORRECT_GAMUT
		/* Create Gamut clut */
		{
			icxLuBase *xluo;		/* Forward ixcLu */
			icmLut *wo;
			out_b2a_callback cx;

			cx.verb = verb;
			cx.wantLab = wantLab;			/* Copy PCS flag over */
			cx.pcsspace = wantLab ? icSigLabData : icSigXYZData;
			cx.devspace = devspace;
			cx.x = (icxLuLut *)AtoB;		/* A2B icxLuLut */

			/* Need to switch AtoB to be override Lab PCS */
			/* Do this the dirty way (Should really create set method) */
			/* ??? does this work, or does change have to be in underlying icc ??? */
			cx.x->ins = cx.x->pcs = icSigLabData;
			/* ~~~ or fix get_gamut to independently override ?? */

			/* Creat a gamut surface */
			if ((cx.gam = AtoB->get_gamut(AtoB, 0.0)) == NULL)
				error("Get_gamut failed: %d, %s",AtoB->pp->errc,AtoB->pp->err);
			
			if ((wo = (icmLut *)wr_icco->read_tag(
			           wr_icco, icSigGamutTag)) == NULL) 
				error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

			if (cx.verb) {
				cx.count = 0;
				cx.last = -1;
				for (cx.total = 1, i = 0; i < wo->inputChan; i++, cx.total *= wo->clutPoints)
					; 
				printf(" 0%%"); fflush(stdout);
			}
			if (wo->set_tables(wo,
					&cx,				/* Context */
					cx.pcsspace,		/* Input color space */
					devspace, 			/* Output color space */
					out_b2a_input,		/* Input transform PCS->PCS' */
					NULL, NULL,			/* Use default Lab' range */
					PCSp_bdist,			/* Lab' -> Boundary distance */
					NULL, NULL,			/* Use default Device' range */
					gamut_output) != 0)	/* Boundary distance to out of gamut value */
				error("Setting 16 bit PCS->Device Lut failed: %d, %s",wr_icco->errc,wr_icco->err);
			if (cx.verb) {
				printf("\n");
			}

			cx.gam->del(cx.gam);		/* Done with gamut object */
		}
#endif /* USE_CORRECT_GAMUT */

		/* Free up xicc stuff */
		AtoB->del(AtoB);
		wr_xicc->del(wr_xicc);

	} else {		/* Gamma/Shaper + matrix profile */
		xicc *wr_xicc;			/* extention object */
		icxLuBase *xluo;		/* Forward ixcLu */
		int flags = 0;

		/* Wrap with an expanded icc */
		if ((wr_xicc = new_xicc(wr_icco)) == NULL)
			error("Creation of xicc failed");
		
		if (verb)
			flags |= ICX_VERBOSE;

		/* Setup Device -> XYZ conversion (Fwd) object from scattered data. */
		if ((xluo = wr_xicc->set_luobj(
		               wr_xicc, icmFwd, isdisp ? icmDefaultIntent : icRelativeColorimetric,
		               icmLuOrdNorm,
		               flags | ICX_SET_WHITE | ICX_SET_BLACK, 		/* Flags */
		               npat, tpat, smooth, avgdev, NULL, ink, iquality)) == NULL)
			error("%d, %s",wr_xicc->errc, wr_xicc->err);

		/* Free up xicc stuff */
		xluo->del(xluo);
		wr_xicc->del(wr_xicc);
	}

	/* Write the file (including all tags) out */
	if ((rv = wr_icco->write(wr_icco,wr_fp,0)) != 0) {
		error("Write file: %d, %s",rv,wr_icco->err);
	}

	/* Close the file */
	wr_icco->del(wr_icco);
	wr_fp->del(wr_fp);

	/* Check the forward profile accuracy against the data points */
	if (verb || verify) {
		icmFile *rd_fp;
		icc *rd_icco;
		icmLuBase *luo;
		double merr = 0.0;
		double aerr = 0.0;
		double nsamps = 0.0;

		/* Open up the file for reading */
		if ((rd_fp = new_icmFileStd_name(file_name,"r")) == NULL)
			error("Write: Can't open file '%s'",file_name);

		if ((rd_icco = new_icc()) == NULL)
			error("Write: Creation of ICC object failed");

		/* Read the header and tag list */
		if ((rv = rd_icco->read(rd_icco,rd_fp,0)) != 0)
			error("Read: %d, %s",rv,rd_icco->err);

		/* Get the Fwd table */
		if ((luo = rd_icco->get_luobj(rd_icco, icmFwd, icAbsoluteColorimetric,
		                              icSigLabData, icmLuOrdNorm)) == NULL) {
			error("%d, %s",rd_icco->errc, rd_icco->err);
		}

		for (i = 0; i < npat; i++) {
			double out[3], ref[3];
			double mxd;

			/* Lookup the profiles PCS for out test patch point */
			if (luo->lookup(luo, out, tpat[i].p) > 1)
				error("%d, %s",rd_icco->errc,rd_icco->err);
		
			/* Our tpat data might be in XYZ, so generate an Lab ref value */
			if (!wantLab) { /* Convert test patch result XYZ to PCS (D50 Lab) */
				icmXYZ2Lab(&icmD50, ref, tpat[i].v);

			} else {
				ref[0] = tpat[i].v[0];
				ref[1] = tpat[i].v[1];
				ref[2] = tpat[i].v[2];
			}

			if (verify && verb) {
				if (devspace == icSigCmykData) {
					printf("[%f] %f %f %f %f -> %f %f %f should be %f %f %f\n",
					       icmLabDE(ref, out),
					       tpat[i].p[0],tpat[i].p[1],tpat[i].p[2],tpat[i].p[3],
					       out[0],out[1],out[2],
					       ref[0],ref[1],ref[2]);
				} else {
					printf("[%f] %f %f %f -> %f %f %f should be %f %f %f\n",
					       icmLabDE(ref, out),
					       tpat[i].p[0],tpat[i].p[1],tpat[i].p[2],
					       out[0],out[1],out[2],
					       ref[0],ref[1],ref[2]);
				}
			}

			/* Check the result */
			mxd = icmLabDE(ref, out);
			if (mxd > merr)
				merr = mxd;

			aerr += mxd;
			nsamps++;
		}
		printf("profile check complete, peak err = %f, avg err = %f\n",merr,aerr/nsamps);

		/* Done with lookup object */
		luo->del(luo);

		/* Close the file */
		rd_icco->del(rd_icco);
		rd_fp->del(rd_fp);
	}
}

