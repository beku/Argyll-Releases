
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
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
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
 *	For flexibility, it might be nice to allow the user to override the
 *  default intent for all 3 output tables, allowing fully custom usage.
 *  Need to break dependence between intents CAM to do this.
 *
 *  Need to make this more of a library:
 *  Fix error handling
 *  fix verbose output
 *  hand icc object back rather than writing file ?
 */

/*
	Outline of code flow:

	profout:
		Create ICC profile and all the tags. Table tags are initialy not set.

		Read in the CGTATS data and convert spectral to PCS if needed.

		Wrap the icc in an xicc, and then create an xicc lookup object
		for the device->PCS tables by calling xicc->set_luobj().

		For a CLUT type profile, create a gamut mapping object and
		setup all the other bits and peices needed to convert a color
		from PCS to device, then use (icc) icmSetMultiLutTables() which will
		call back the out_b2a_input(), out_b2a_clut() and out_b2a_output()
		functions.
 */

#undef DEBUG				/* Print B2A processing information */
#undef DEBUG_ONE			/* Test a particular value */

#define verbo stdout

#undef IMP_MONO					/* [Undef] Turn on development code */

#define NO_B2A_PCS_CURVES		/* [Define] PCS curves seem to make B2A less accurate. Why ? */
#define USE_CAM_CLIP_OPT		/* [Define] Clip out of gamut in CAM space rather than PCS */
#define USE_LEASTSQUARES_APROX	/* [Define] Use least squares fitting approximation */
#undef USE_EXTRA_FITTING		/* [Undef] Turn on data point error compensation */
#undef USE_2PASSSMTH			/* [Undef] Turn on Gaussian smoothing */
#undef DISABLE_GAMUT_TAG		/* [Undef] To dispable gamut tag */
#undef WARN_CLUT_CLIPPING		/* [Undef] Print warning if setting clut clips */

#include <stdio.h>
#include "numlib.h"
#include "icc.h"
#include "cgats.h"
#include "xicc.h"
#include "counters.h"
#include "rspl.h"
#include "insttypes.h"
#include "prof.h"
#include "gamut.h"
#include "gammap.h"

#ifndef MAX_CAL_ENT
#define MAX_CAL_ENT 4096
#endif

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


#ifdef DEBUG
#undef DBG
#define DBG(xxx) printf xxx ;
#else
#undef DBG
#define DBG(xxx) 
#endif

/* ---------------------------------------- */

/* structure to support output icc B2A Lut initialisation calbacks. */
/* Note that we don't cope with a LUT matrix - assume it's unity. */

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
	icxLuBase *ixp;			/* Source profile perceptual PCS to CAM conversion */
	icxLuBase *ox;			/* Destination profile CAM to std PCS conversion */
							/* (This is NOT used for the colorimetric B2A table creation!) */

	icRenderingIntent abs_intent;	/* Desired abstract profile rendering intent */
	icxLuBase *abs_luo;		/* abstract profile transform in PCS, NULL if none */

	double xyzscale[2];		/* < 1.0 if XYZ is to be scaled in destination space */
							/* for perceptual [0], and saturation [1] */
	double swxyz[3];		/* Source white point in XYZ */

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

/* --------------------------------------------------------- */

/* Extra non-linearity applied to BtoA XYZ PCS */
/* This distributes the LUT indexes more evenly in */
/* perceptual space, greatly improving the B2A accuracy of XYZ LUT */
/* Since typically XYZ doesn't use the full range of 0-2.0 allowed */
/* for in the encoding, we scale the cLUT index values to use the 0-1.3 range */
static void xyzcurve(double *out, double *in) {
	int i;
	double sc = 2.0/1.3 * 65535.0/32768.0;

	/* Use an L* like curve, scaled to the maximum XYZ valu */
	out[0] = in[0]/sc;
	out[1] = in[1]/sc;
	out[2] = in[2]/sc;
	for (i = 0; i < 3; i++) {
		if (out[i] > 0.08)
			out[i] = pow((out[i] + 0.16)/1.16, 3.0);
		else
			out[i] = out[i]/9.032962896;
	}
	out[0] = out[0] * sc;
	out[1] = out[1] * sc;
	out[2] = out[2] * sc;
}

static void invxyzcurve(double *out, double *in) {
	int i;
	double sc = 2.0/1.3 * 65535.0/32768.0;

	out[0] = in[0]/sc;
	out[1] = in[1]/sc;
	out[2] = in[2]/sc;
	for (i = 0; i < 3; i++) {
		if (out[i] > 0.008856451586)
			out[i] = 1.16 * pow(out[i],1.0/3.0) - 0.16;
		else
			out[i] = 9.032962896 * out[i];
	}
	out[0] = out[0] * sc;
	out[1] = out[1] * sc;
	out[2] = out[2] * sc;
}
/* --------------------------------------------------------- */

/* sRGB device gamma encoded value to linear value 0.0 .. 1.0 */
static double gdv2dv(double iv) {
	double ov;

	if (iv < 0.04045)
		ov = iv/12.92;
	else
		ov = pow((iv + 0.055)/1.055, 2.4);
	return ov;
}


/* --------------------------------------------------------- */
/* NOTE :- the assumption that each stage of the BtoA is a mirror */
/* of the AtoB makes for inflexibility. */
/* Perhaps it would be better to remove this asumption from the */
/* out_b2a_clut processing ? */
/* To do this we then need inv_out_b2a_input(), and */
/* inv_out_b2a_output(), and we need to clearly distinguish between */
/* AtoB PCS' & DEV', and BtoA PCS' & DEV', since they are not */
/* necessarily the same... */


/* B2A Input table is the inverse of the AtoB output table */
/* Input PCS output PCS'' */
void out_b2a_input(void *cntx, double out[3], double in[3]) {
	out_b2a_callback *p = (out_b2a_callback *)cntx;

	DBG(("out_b2a_input got         PCS %f %f %f\n",in[0],in[1],in[2]))

	/* PCS to PCS' */
	if (p->noPCScurves) {
		out[0] = in[0];
		out[1] = in[1];
		out[2] = in[2];
	} else {
		if (p->x->inv_output(p->x, out, in) > 1)
			error("%d, %s",p->x->pp->errc,p->x->pp->err);
	}
	/* PCS' to PCS'' */
	if (p->pcsspace == icSigXYZData)	/* Apply XYZ non-linearity curve */
		invxyzcurve(out, out);

	DBG(("out_b2a_input returning PCS'' %f %f %f\n",out[0],out[1],out[2]))
}

/* clut - multitable */
/* Input PCS' output Dev' */
/* We're applying any abstract profile after gamut mapping, */
/* on the assumption is is primarily being used to "correct" the */
/* output device. Ideally the gamut mapping should take the change */
/* the abstract profile has on the output device into account, but */
/* currently we're not doing this.. */
void out_b2a_clut(void *cntx, double *out, double in[3]) {
	out_b2a_callback *p = (out_b2a_callback *)cntx;
	double inn[3], in1[3];
	int tn;

	DBG(("\nout_b2a_clut got      PCS'' %f %f %f\n",in[0],in[1],in[2]))

	in1[0] = in[0];		/* in[] may be aliased with out[] */
	in1[1] = in[1];		/* so take a copy.  */
	in1[2] = in[2];

	DBG(("out_b2a_clut got       PCS' %f %f %f\n",in[0],in[1],in[2]))

	if (p->pcsspace == icSigXYZData)		/* Undo effects of extra XYZ non-linearity curve */
		xyzcurve(in1, in1);

	inn[0] = in1[0];	/* Copy of PCS' for 2nd and 3rd tables */
	inn[1] = in1[1];
	inn[2] = in1[2];

	if (p->abs_luo != NULL) {	/* Abstract profile to apply to first table only. */

		if (!p->noPCScurves) {	/* Convert from PCS' to PCS */
			if (p->x->output(p->x, in1, in1) > 1)
				error("%d, %s",p->x->pp->errc,p->x->pp->err);
		}
		do_abstract(p, in1, in1);		/* Abstract profile to apply to first table */
		DBG(("through abstract prof   PCS %f %f %f\n",in1[0],in1[1],in1[2]))
		/* We now have PCS */
	}

	if (p->noPCScurves || p->abs_luo != NULL) {	/* We were given PCS or have converted to PCS */

		/* PCS to PCS' */
		if (p->x->inv_output(p->x, in1, in1) > 1)
			error("%d, %s",p->x->pp->errc,p->x->pp->err);

		DBG(("noPCScurves = %d, abs_luo = 0x%x\n",p->noPCScurves,p->abs_luo))
		DBG(("convert to PCS' got         %f %f %f\n",in1[0],in1[1],in1[2]))
	}

	/* Invert AtoB clut (PCS' to Dev') Colorimetric */
	/* to producte the colorimetric tables output. */
	/* (Note that any aux target if we were using one, would be in Dev space) */
	if (p->x->inv_clut(p->x, out, in1) > 1)
		error("%d, %s",p->x->pp->errc,p->x->pp->err);

	DBG(("convert PCS' to DEV' got    %f %f %f %f\n",out[0],out[1],out[2],out[3]))

	if (p->ntables > 1) {		/* Do first part once for both intents */

		DBG(("\n"))

		/* Starting with original input inn[] PCS' (no abstract profile applied to other tables) */
		in1[0] = inn[0];
		in1[1] = inn[1];
		in1[2] = inn[2];

		if (!p->noPCScurves) {					/* Convert from PCS' to PCS */
			if (p->x->output(p->x, in1, in1) > 1)
				error("%d, %s",p->x->pp->errc,p->x->pp->err);
		}
		DBG(("convert PCS' to PCS got %f %f %f\n",in1[0],in1[1],in1[2]))

		/* Convert from PCS to CAM/Gamut mapping space */
		p->ixp->fwd_relpcs_outpcs(p->ixp, p->pcsspace, in1, in1);

		DBG(("convert PCS to CAM got %f %f %f\n",in1[0],in1[1],in1[2]))

		/* Apply gamut mapping in CAM space for remaining tables */
		/* and create the output values */
		for (tn = 1; tn < p->ntables; tn++) {
			double in2[3];

			out += p->ochan;	/* next table/intent */
			in2[0] = in1[0];	/* Copy in1[] so it can be used for both tables */
			in2[1] = in1[1];
			in2[2] = in1[2];

			/* Do luminence scaling if requested */
			if (p->xyzscale[tn-1] < 1.0) {
				double xyz[3];

				DBG(("got xyzscale = %f\n",p->xyzscale[tn-1]))
				DBG(("PCS %f %f %f\n",in2[0], in2[1], in2[2]))
				/* Convert our CAM to XYZ */
				/* We're being bad in delving inside the xluo, but we'll fix it latter */
				p->ox->cam->cam_to_XYZ(p->ox->cam, xyz, in2);

				DBG(("XYZ %f %f %f\n",xyz[0], xyz[1], xyz[2]))
				/* Scale it */
				xyz[0] *= p->xyzscale[tn-1];
				xyz[1] *= p->xyzscale[tn-1];
				xyz[2] *= p->xyzscale[tn-1];

				DBG(("scaled XYZ %f %f %f\n",xyz[0], xyz[1], xyz[2]))
				/* Convert back to CAM */
				/* We're being bad in delving inside the xluo, but we'll fix it latter */
				p->ox->cam->XYZ_to_cam(p->ox->cam, in2, xyz);

				DBG(("scaled PCS %f %f %f\n",in2[0], in2[1], in2[2]))
			}

			if (tn == 1) {
				DBG(("percep gamut map in %f %f %f\n",in2[0],in2[1],in2[2]))
				p->pmap->domap(p->pmap, in2, in2);	/* Perceptual mapping */
				DBG(("percep gamut map out %f %f %f\n",in2[0],in2[1],in2[2]))
			} else {
				DBG(("sat gamut map in %f %f %f\n",in2[0],in2[1],in2[2]))
				p->smap->domap(p->smap, in2, in2);	/* Saturation mapping */
				DBG(("sat gamut map got %f %f %f\n",in2[0],in2[1],in2[2]))
			}

			/* Convert from Gamut maping/CAM space to PCS */
			p->ox->bwd_outpcs_relpcs(p->ox, p->pcsspace, in2, in2);

			DBG(("convert CAM to PCS got %f %f %f\n",in2[0],in2[1],in2[2]))

			if (p->abs_luo != NULL)
				do_abstract(p, in2, in2);	/* Abstract profile to other tables after gamut map */

			/* Convert from PCS to PCS' */
			if (p->x->inv_output(p->x, in2, in2) > 1)
				error("%d, %s",p->x->pp->errc,p->x->pp->err);
				DBG(("convert PCS to PCS' got %f %f %f\n",in2[0],in2[1],in2[2]))

			/* Invert AtoB clut (PCS' to Dev') */
			/* to producte the perceptual or saturation tables output. */
			/* (Note that any aux target if we were using one, would be in Dev space) */
			if (p->x->inv_clut(p->x, out, in2) > 1)
				error("%d, %s",p->x->pp->errc,p->x->pp->err);
			DBG(("convert PCS' to DEV' got %f %f %f %f\n",out[0],out[1],out[2],out[3]))
		}
	}

	DBG(("out_b2a_clut returning DEV' %f %f %f\n",out[0],out[1],out[2]))

	if (p->verb) {		/* Output percent intervals */
		int pc;
		p->count++;
		pc = (int)(p->count * 100.0/p->total + 0.5);
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

	DBG(("out_b2a_output got      DEV' %f %f %f\n",in[0],in[1],in[2]))

	if (p->x->inv_input(p->x, out, in) > 1)
		error("%d, %s",p->x->pp->errc,p->x->pp->err);

	DBG(("out_b2a_output returning DEV %f %f %f\n",out[0],out[1],out[2]))
}

/* --------------------------------------------------------- */

/* PCS' -> distance to gamut boundary */
static void PCSp_bdist(void *cntx, double out[1], double in[3]) {
	out_b2a_callback *p = (out_b2a_callback *)cntx;
	double pcs[3];				/* PCS value of input */
	double nrad;				/* normalised radial value of point */
	double np[3];				/* Nearest point */
	double gdist;				/* Out of gamut distance */
	
//printf("~1 bdist got PCS %f %f %f\n",in[0],in[1],in[2]);
	/* Do PCS' -> PCS */
	if (p->x->inv_output(p->x, pcs, in) > 1)
		error("%d, %s",p->x->pp->errc,p->x->pp->err);

	/* If PCS is XYZ, convert to Lab */
	if (p->wantLab == 0) {
		icmXYZ2Lab(&icmD50, pcs, pcs);
//printf("~1 bdist converted to Lab %f %f %f\n",pcs[0],pcs[1],pcs[2]);
	}

	/* Check if this point is in or out of gamut */
	nrad = p->gam->nradial(p->gam, np, pcs);

//printf("~1 nrad = %f\n",nrad);

	/* Radial rather than nearest distance seems the best overall. */

	gdist = icmNorm33(np, pcs);
	if (nrad <= 1.0)
		gdist = -gdist;		/* -ve delta E if within gamut */
	
//printf("~1 gdist %f\n",gdist);

	/* Distance in PCS space will be roughly -128 -> 128 */
	/* Clip to range -20 - +20, then scale to 0.0 - 1.0 */
	if (gdist < -20.0)
		gdist = -20.0;
	else if (gdist > 20.0)
		gdist = 20.0;

	out[0] = (gdist + 20.0)/40.0;
//printf("~1 bdist returning %f\n",out[0]);

	if (p->verb) {		/* Output percent intervals */
		int pc;
		p->count++;
		pc = (int)(p->count * 100.0/p->total + 0.5);
		if (pc != p->last) {
			printf("\r%2d%%",pc); fflush(stdout);
			p->last = pc;
		}
	}
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

/* -------------------------------------------------------------- */
/* powell() callback to set XYZ white scaling factor for */
/* Absolute Appearance mode with scaling intent */

static double xyzoptfunc(void *cntx, double *v) {
	out_b2a_callback *p = (out_b2a_callback *)cntx;
	double swxyz[3], jab[3], dev[MAX_CHAN];
	double rv;
	int rc;

	rv = 2.0 - v[0];	/* Make Y as large as possible */

	/* If we wanted to use this function to maximise the brightness */
	/* we would not limit the scale to 1.0 */
	if (v[0] > 1.0) {
		rv += 1000.0;
		return rv;
	}
	if (v[0] < 0.0) {
		rv += 100.0;
		return rv;
	}
	swxyz[0] = v[0] * p->swxyz[0];
	swxyz[1] = v[0] * p->swxyz[1];
	swxyz[2] = v[0] * p->swxyz[2];

//printf("~1 scaled white XYZ = %f %f %f\n", swxyz[0], swxyz[1], swxyz[2]);

	/* We're being bad in delving inside the xluo, but we'll fix it latter */
	p->ox->cam->XYZ_to_cam(p->ox->cam, jab, swxyz);

//printf("~1 scaled white Jab = %f %f %f\n", jab[0], jab[1], jab[2]);

	/* Run the target PCS backwards through the output space to see if it clips */
	rc = p->ox->inv_lookup(p->ox, dev, jab);
//printf("~1 device = %f %f %f, rc = %d\n", dev[0], dev[1], dev[2],rc);
	if (rc != 0)
		rv += 500.0;

//printf("~1 xyzoptfunc rv %f from xyzscale %f\n\n",rv,v[0]);
	return rv;
}

/* -------------------------------------------------------------- */
/* Make an output device profile, where a forward mapping is from */
/* RGB/CMYK to XYZ/Lab space */
void
make_output_icc(
	prof_atype ptype,		/* Profile algorithm type */
	int mtxtoo,				/* NZ if matrix tags should be created for Display XYZ cLUT */
	icmICCVersion iccver,	/* ICC profile version to create */
	int verb,				/* Vebosity level, 0 = none */
	int iquality,			/* A2B table quality, 0..3 */
	int oquality,			/* B2A table quality, 0..3 */
	int noisluts,			/* nz to supress creation of input (Device) shaper luts */
	int noipluts,			/* nz to supress creation of input (Device) position luts */
	int nooluts,			/* nz to supress creation of output (PCS) shaper luts */
	int nocied,				/* nz to supress inclusion of .ti3 data in profile */
	int noptop,				/* nz to use colorimetic source gamut to make perceptual table */ 
	int nostos,				/* nz to use colorimetic source gamut to make saturation table */
	int gamdiag,			/* Make gamut mapping diagnostic wrl plots */
	int verify,				/* nz to print verification */
	icxInk *oink,			/* Ink limit/black generation setup (NULL if n/a) */
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
	icxViewCond *ivc_p,		/* Input Viewing Parameters for CAM */
	icxViewCond *ovc_p,		/* Output Viewing Parameters for CAM */
	int ivc_e,				/* Input Enumerated viewing condition */
	int ovc_e,				/* Output Enumerated viewing condition */
	icxGMappingIntent *pgmi,/* Perceptual gamut mapping intent */
	icxGMappingIntent *sgmi,/* Saturation gamut mapping intent */
	profxinf *xpi			/* Optional Profile creation extra data */
) {
	int isdisp;				/* nz if this is a display device, 0 if output */
	double dispLuminance = 0.0;	/* Display luminance. 0 if not known */
	int allintents;			/* nz if all intents should possibly be created */
	icmFile *wr_fp;
	icc *wr_icco;
	int npat;				/* Number of patches */
	cow *tpat;				/* Patch input values */
	int i, j, rv = 0;
	icColorSpaceSignature devspace = icmSigDefaultData;	/* The device colorspace */
	inkmask imask = 0;		/* Device inkmask */ 
	int isAdditive = 0;		/* 0 if subtractive, 1 if additive colorspace */
	int isLab = 0;			/* 0 if input is XYZ, 1 if input is Lab */
	int wantLab = 0;		/* 0 if want XYZ, 1 want Lab */
	int isLut = 0;			/* 0 if shaper+ matrix, 1 if lut type */
	int isShTRC = 0;		/* 0 if separate gamma/shaper TRC, 1 if shared */
	int devchan = 0;		/* Number of device chanels */
	int a2binres = 0;		/* A2B input (device) table resolution */
	int a2bres = 0;			/* A2B clut resolution */
	int a2boutres = 0;		/* A2B output (PCS) table resolution */
	int b2ainres = 0;		/* B2A input (PCS) table resolution */
	int b2ares = 0;			/* B2A clut resolution */
	int b2aoutres = 0;		/* B2A output (device) table resolution */
	xcal *cal = NULL;		/* Calibration if present, NULL if none */
	icxInk iink;			/* Source profile ink limit values */

	iink.tlimit = -1.0;		/* default to unknown */
	iink.klimit = -1.0;

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

	{
		int ti;

		if ((ti = icg->find_kword(icg, 0, "DEVICE_CLASS")) < 0)
			error ("Input file doesn't contain keyword DEVICE_CLASS");

		if (strcmp(icg->t[0].kdata[ti],"DISPLAY") == 0) {
			isdisp = 1;
			if (isLut && ipname != NULL)
				allintents = 1;
			else
				allintents = 0;		/* Only the default intent */
		} else {
			isdisp = 0;
			allintents = 1;
		}
	}

	/* See if display luminance data is present */
	if (isdisp) {
		int ti;

		if ((ti = icg->find_kword(icg, 0, "LUMINANCE_XYZ_CDM2")) >= 0) {
			if (sscanf(icg->t[0].kdata[ti], " %*lf %lf %*lf ",&dispLuminance) != 1)
				dispLuminance = 0.0;
		}
	}

	/* Figure out what sort of device colorspace it is */
	{
		int ti;

		if ((ti = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
			error("Input file doesn't contain keyword COLOR_REPS");

		if (strcmp(icg->t[0].kdata[ti],"CMYK_XYZ") == 0) {
			imask = ICX_CMYK;
			devspace = icSigCmykData;
			devchan = 4;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMYK_LAB") == 0) {
			imask = ICX_CMYK;
			devspace = icSigCmykData;
			devchan = 4;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMY_XYZ") == 0) {
			imask = ICX_CMY;
			devspace = icSigCmyData;
			devchan = 3;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMY_LAB") == 0) {
			imask = ICX_CMY;
			devspace = icSigCmyData;
			devchan = 3;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"RGB_XYZ") == 0) {
			imask = ICX_RGB;
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"RGB_LAB") == 0) {
			imask = ICX_RGB;
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 1;
			isAdditive = 1;
		} else if (   strcmp(icg->t[0].kdata[ti],"iRGB_XYZ") == 0) {
			imask = ICX_IRGB;
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"iRGB_LAB") == 0) {
			imask = ICX_IRGB;
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 1;
			isAdditive = 1;
#ifdef IMP_MONO
		} else if (strcmp(icg->t[0].kdata[ti],"K_XYZ") == 0) {
			imask = ICX_K;
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"K_LAB") == 0) {
			imask = ICX_K;
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"W_XYZ") == 0) {
			imask = ICX_W;
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"W_LAB") == 0) {
			imask = ICX_W;
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 1;
			isAdditive = 1;
#endif /* IMP_MONO */
		} else 
			error("Output device input file has unhandled color representation '%s'",
			                                                     icg->t[0].kdata[ti]);
		/* Figure out some suitable table sizes */
		if (devchan >= 4) {			/* devchan == 4 or greater */
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
		} else if (devchan >= 2) {		/* devchan == 2 or 3 */
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

	/* See if there is a calibration in the .ti3, and read it if there is */
	{
		int oi, tab;

		if ((oi = icg->get_oi(icg, "CAL")) < 0)
			error("Expect CAL type to be registered");

		for (tab = 0; tab < icg->ntables; tab++) {
			if (icg->t[tab].tt == tt_other && icg->t[tab].oi == oi) {
				break;
			}
		}
		if (tab < icg->ntables) {

			if ((cal = new_xcal()) == NULL) {
				error("new_xcal failed");
			}
			if (cal->read_cgats(cal, icg, tab, in_name) != 0)  {
				error("%s",cal->err);
			}

			if (cal->devmask != imask)
				error("Calibrate colorspace %s doesn't match .ti3 %s",
				       icx_inkmask2char(cal->devmask, 1), 
				       icx_inkmask2char(imask, 1));

			if ((isdisp && cal->devclass != icSigDisplayClass)
			 || (!isdisp && cal->devclass != icSigOutputClass)) 
				error("Calibration class %s doesn't match .ti3 %s",
				       icm2str(icmProfileClassSignature,cal->devclass),
				       isdisp ? "Display" : "Output");
		}
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
    	wh->renderingIntent = icRelativeColorimetric;	/* Should let user set this!!! */

		/* Values that should be set before writing */
		if (xpi != NULL && xpi->manufacturer != 0L)
			wh->manufacturer = xpi->manufacturer;
		else
			wh->manufacturer = str2tag("????");

		if (xpi != NULL && xpi->model != 0L)
			wh->model        = xpi->model;
		else
	    	wh->model        = str2tag("????");

		/* Values that may be set before writing */
		if (xpi != NULL && xpi->creator != 0L)
			wh->creator = xpi->creator;

#ifdef NT
		wh->platform = icSigMicrosoft;
#endif
#ifdef __APPLE__
		wh->platform = icSigMacintosh;
#endif
#if defined(UNIX) && !defined(__APPLE__)
		wh->platform = icmSig_nix;
#endif
	}

	/* mtxtoo only applies to Display cLUT profiles */
	if (!isLut
		|| wr_icco->header->deviceClass != icSigDisplayClass
	    || wr_icco->header->pcs != icSigXYZData) {
		mtxtoo = 0;
	}

	/* Set the version of ICC profile we want */
	if (isdisp && allintents) {
		if (iccver < icmVersion2_4) {
			iccver = icmVersion2_4;		/* Need 2.4.0 for Display intents */
			if (verb)
				fprintf(verbo,"Bumped ICC version to 2.4.0 to accomodate multiple Display intents\n");
		}
	}
	if (wr_icco->set_version(wr_icco, iccver) != 0)
		error("set_version failed: %d, %s",wr_icco->errc,wr_icco->err);

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
	/* Display Luminance tag */
	if (isdisp && dispLuminance > 0.0) {
		{
			icmXYZArray *wo;;
	
			if ((wo = (icmXYZArray *)wr_icco->add_tag(
			           wr_icco, icSigLuminanceTag, icSigXYZArrayType)) == NULL) 
				error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

			wo->size = 1;
			wo->allocate((icmBase *)wo);	/* Allocate space */
			wo->data[0].X = 0.0;			/* Set a default value */
			wo->data[0].Y = dispLuminance;	/* Set a default value */
			wo->data[0].Z = 0.0;			/* Set a default value */
		}
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

	/* Colorant Table Tag: */
	{
		unsigned int i;
		icmColorantTable *wo;
		if ((wo = (icmColorantTable *)wr_icco->add_tag(
		           wr_icco, icSigColorantTableTag, icSigColorantTableType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

   		wo->count = icmCSSig2nchan(devspace);
   		if (wo->count != (unsigned long)icx_noofinks(imask))
			error("Interna: device colorspace and inkmask conflict!");

		wo->allocate((icmBase *)wo);	/* Allocate ColorantTable structures */

		for (i = 0; i < wo->count; i++) {
			inkmask iimask;			/* Individual ink mask */
			char *name;
			
			iimask = icx_index2ink(imask, i);
			name = icx_ink2string(iimask); 
			if (strlen(name) > 31)
				error("Internal: colorant name exceeds 31 characters");
			strcpy(wo->data[i].name, name);
		}
		/* Fill in the colorant PCS values when we've got something to lookup */
	}

	/* vcgt tag */
	if (isdisp && cal != NULL) {		/* We've been given vcgt information */
		int j, i;
		int ncal = 256;					/* This is safe with other s/w */
		icmVideoCardGamma *wo;
		wo = (icmVideoCardGamma *)wr_icco->add_tag(wr_icco,
		                           icSigVideoCardGammaTag, icSigVideoCardGammaType);
		if (wo == NULL)
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->tagType = icmVideoCardGammaTableType;
		wo->u.table.channels = 3;			/* rgb */
		wo->u.table.entryCount = ncal;		/* full lut */
		wo->u.table.entrySize = 2;			/* 16 bits */
		wo->allocate((icmBase*)wo);
		for (j = 0; j < 3; j++) {
			for (i = 0; i < ncal; i++) {
				double cc, vv = i/(ncal - 1.0);
				cc = cal->interp_ch(cal, j, vv);
				if (cc < 0.0)
					cc = 0.0;
				else if (cc > 1.0)
					cc = 1.0;
				((unsigned short*)wo->u.table.data)[ncal * j + i] = (int)(cc * 65535.0 + 0.5);
			}
		}
	}

	if (isLut) {		/* Lut type profile */

		/* Up to and including ICC Version 2.3, Display LUT profiles were assumed */
		/* to have AtoB0 with no interpretation of the intent, which */
		/* implies Relative Colorimetric. */

		/* 16 bit dev -> pcs lut: (A2B) */
		{
			icmLut *wo;
	
			if (!allintents) {	/* Only A2B0, no intent */
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

		if (allintents) {	/* All the intents may be needed */

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

			if (!allintents) {	/* Only B2A0, no intent */
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

		if (allintents) {	/* All the intents may be needed */

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
				wo->clutPoints = b2ares;
				wo->outputEnt = 256;
				wo->allocate((icmBase *)wo);/* Allocate space */
			}
		}
	}

	/* shaper + matrix type tags */
	if (!isLut
	|| (   wr_icco->header->deviceClass == icSigDisplayClass
	    && wr_icco->header->pcs == icSigXYZData)) {

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

			wor->data[0].X = 1.0; wor->data[0].Y = 0.0; wor->data[0].Z = 0.0;
			wog->data[0].X = 0.0; wog->data[0].Y = 1.0; wog->data[0].Z = 0.0;
			wob->data[0].X = 0.0; wob->data[0].Y = 0.0; wob->data[0].Z = 1.0;

			/* Setup deliberately wrong dummy values (red and green swapped). */
			/* icxMatrix may override override these later */
			wor->data[0].X = 0.385147; wor->data[0].Y = 0.716873; wor->data[0].Z = 0.097076;
			wog->data[0].X = 0.436066; wog->data[0].Y = 0.222488; wog->data[0].Z = 0.013916;
			wob->data[0].X = 0.143066; wob->data[0].Y = 0.060608; wob->data[0].Z = 0.714096;
		}

		/* Red, Green and Blue Tone Reproduction Curve Tags: */
		{
			icmCurve *wor, *wog, *wob;
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
	
			/* Shaper */
			if (ptype == prof_shamat || ptype == prof_sha1mat || ptype == prof_clutXYZ) {
				wor->flag = wog->flag = wob->flag = icmCurveSpec; 
				wor->size = wog->size = wob->size = 256;			/* Number of entries */
			} else {						/* Gamma */
				wor->flag = wog->flag = wob->flag = icmCurveGamma;
				wor->size = wog->size = wob->size = 1;				/* Must be 1 for gamma */
			}
			wor->allocate((icmBase *)wor);	/* Allocate space */
			wog->allocate((icmBase *)wog);
			wob->allocate((icmBase *)wob);

			/* Setup a default sRGB like curve. */
			/* icxMatrix will may override curve values later */
			for (i = 0; i < wor->size; i++) {
				wor->data[i] = 
				wog->data[i] = 
				wob->data[i] = gdv2dv(i/(wor->size-1.0));
			}

		}
	}
	/* .ti3 Sample data use to create profile, plus any calibration curves: */
	if (nocied == 0) {
		icmText *wo;
		char *crt;
		FILE *fp;

		if ((wo = (icmText *)wr_icco->add_tag(
		           wr_icco, icmMakeTag('t','a','r','g'), icSigTextType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

#if defined(O_BINARY) || defined(_O_BINARY)
	    if ((fp = fopen(in_name, "rb")) == NULL)
#else
	    if ((fp = fopen(in_name, "r")) == NULL)
#endif
			error("Unable to open input file '%s' for reading",in_name);

		if (fseek(fp, 0, SEEK_END))
			error("Unable to seek to end of file '%s'",in_name);
		wo->size = ftell(fp) + 1;		/* Size needed + null */
		wo->allocate((icmBase *)wo);/* Allocate space */

		if (fseek(fp, 0, SEEK_SET))
			error("Unable to seek to end of file '%s'",in_name);

		if (fread(wo->data, 1, wo->size-1, fp) != wo->size-1)
			error("Failed to read file '%s'",in_name);
		wo->data[wo->size-1] = '\000';
		fclose(fp);

		/* Duplicate for compatibility */
		if (wr_icco->link_tag(
		         wr_icco, icmMakeTag('D','e','v','D'), icmMakeTag('t','a','r','g')) == NULL) 
			error("link_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
		if (wr_icco->link_tag(
		         wr_icco, icmMakeTag('C','I','E','D'), icmMakeTag('t','a','r','g')) == NULL) 
			error("link_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
	}

	if ((npat = icg->t[0].nsets) <= 0)
		error("No sets of data");

	if (verb)
		fprintf(verbo,"No of test patches = %d\n",npat);

	/* Allocate arrays to hold test patch input and output values */
	if ((tpat = (cow *)malloc(sizeof(cow) * npat)) == NULL)
		error("Malloc failed - tpat[]");

	/* Read in the CGATs fields */
	{
		int ti, ii, ci, mi, yi, ki;
		int Xi, Yi, Zi;

		/* Read the ink limit */
		if (oink != NULL && (ii = icg->find_kword(icg, 0, "TOTAL_INK_LIMIT")) >= 0) {
			double ilimit = -1;
			ilimit = atof(icg->t[0].kdata[ii]);
			if (ilimit > 1e-4 && ilimit <= 400.0 && oink->tlimit < 0.0) {
				oink->tlimit = ilimit/100.0;	/* Set requested ink limit */
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
				tpat[i].w = 1.0;
				tpat[i].p[0] = *((double *)icg->t[0].fdata[i][ci]) / 100.0;
				tpat[i].p[1] = *((double *)icg->t[0].fdata[i][mi]) / 100.0;
				tpat[i].p[2] = *((double *)icg->t[0].fdata[i][yi]) / 100.0;
				tpat[i].p[3] = *((double *)icg->t[0].fdata[i][ki]) / 100.0;
				if (tpat[i].p[0] > 1.0
				 || tpat[i].p[1] > 1.0
				 || tpat[i].p[2] > 1.0
				 || tpat[i].p[3] > 1.0) {
					double bgst = 0.0;
					int j, bj = 0;
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
			if (!isdisp && fwacomp) {
				double nw = 0.0;		/* Number of media white patches */
				xspect mwsp;			/* Medium spectrum */
				instType itype;			/* Spectral instrument type */
				xspect insp;			/* Instrument illuminant */

				mwsp = sp;		/* Struct copy */

	 			if ((ti = icg->find_kword(icg, 0, "TARGET_INSTRUMENT")) < 0)
					error ("Can't find target instrument needed for FWA compensation");

				if ((itype = inst_enum(icg->t[0].kdata[ti])) == instUnknown)
					error ("Unrecognised target instrument '%s'", icg->t[0].kdata[ti]);

				if (inst_illuminant(&insp, itype) != 0)
					error ("Instrument doesn't have an FWA illuminent");

				/* Find the media white spectral reflectance */
				for (j = 0; j < mwsp.spec_n; j++)
					mwsp.spec[j] = 0.0;

				/* Compute the mean of all the media white patches */
				for (i = 0; i < npat; i++) {
					int use = 0;

					if (devspace == icSigGrayData) {
						if (isAdditive) {
							if (*((double *)icg->t[0].fdata[i][ci]) > (100.0 - 0.1))
								use = 1;
						} else {
							if (*((double *)icg->t[0].fdata[i][ci]) < 0.1)
								use = 1;
						}
					} else if (devspace == icSigRgbData) {
						if (*((double *)icg->t[0].fdata[i][ci]) > (100.0 - 0.1)
						 && *((double *)icg->t[0].fdata[i][mi]) > (100.0 - 0.1)
						 && *((double *)icg->t[0].fdata[i][yi]) > (100.0 - 0.1))
							use = 1;
					} else if (devspace == icSigCmyData) {
						if (*((double *)icg->t[0].fdata[i][ci]) < 0.1 
						 && *((double *)icg->t[0].fdata[i][mi]) < 0.1 
						 && *((double *)icg->t[0].fdata[i][yi]) < 0.1)
							use = 1;
					} else {	/* Assume CMYK */

						if (*((double *)icg->t[0].fdata[i][ci]) < 0.1
						 && *((double *)icg->t[0].fdata[i][mi]) < 0.1
						 && *((double *)icg->t[0].fdata[i][yi]) < 0.1
						 && *((double *)icg->t[0].fdata[i][ki]) < 0.1) {
							use = 1;
						}
					}

					if (use) {
						/* Read the spectral values for this patch */
						for (j = 0; j < mwsp.spec_n; j++) {
							mwsp.spec[j] += *((double *)icg->t[0].fdata[i][spi[j]]);
						}
						nw++;
					}
				}

				if (nw == 0.0) {
					warning("Can't find a media white patch to init FWA");

					/* Track the maximum reflectance for any band to determine white. */
					/* This might give bogus results if there is no white patch... */
					for (i = 0; i < npat; i++) {
						for (j = 0; j < mwsp.spec_n; j++) {
							double rv = *((double *)icg->t[0].fdata[i][spi[j]]);
							if (rv > mwsp.spec[j])
								mwsp.spec[j] = rv;
						}
					}
					nw++;
				}
				for (j = 0; j < mwsp.spec_n; j++) {
					mwsp.spec[j] /= nw;	/* Compute average */
				}

				/* (Note that sp and mwsp.norm is set to 100.0) */
				if (sp2cie->set_fwa(sp2cie, &insp, &mwsp)) 
					error ("Set FWA on sp2cie failed");

				if (verb) {
					double FWAc;
					sp2cie->get_fwa_info(sp2cie, &FWAc);
					fprintf(verbo,"FWA content = %f\n",FWAc);
				}
			}

			for (i = 0; i < npat; i++) {

				tpat[i].w = 1.0;
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
		
			flags |= ICX_CLIP_NEAREST;		/* This will avoid clip caused rev setup */

			if (noisluts)
				flags |= ICX_NO_IN_SHP_LUTS;

			if (noipluts)
				flags |= ICX_NO_IN_POS_LUTS;

			if (nooluts)
				flags |= ICX_NO_OUT_LUTS;

			if (verb)
				flags |= ICX_VERBOSE;

			/* Setup Device -> PCS conversion (Fwd) object from scattered data. */
			if ((AtoB = wr_xicc->set_luobj(
			               wr_xicc, icmFwd, !allintents ? icmDefaultIntent : icRelativeColorimetric,
			               icmLuOrdNorm,
#ifdef USE_EXTRA_FITTING
			               ICX_EXTRA_FIT |
#endif
#ifdef USE_2PASSSMTH
			               ICX_2PASSSMTH |
#endif
			               flags | ICX_SET_WHITE | ICX_SET_BLACK, 		/* Flags */
			               npat, tpat, dispLuminance, -1.0, smooth, avgdev,
			               NULL, oink, cal, iquality)) == NULL)
				error("%d, %s",wr_xicc->errc, wr_xicc->err);

			AtoB->del(AtoB);		/* Done with lookup */
		}

		/* Create B2A clut */
		{
			icc *src_icco = NULL;
			xicc *src_xicc = NULL;	/* Source profile */
			icxViewCond ivc;		/* Input Viewing Condition for CAM */
			icxViewCond ovc;		/* Output Viewing Condition for CAM */
			icmLut *wo[3];
			icmFile *abs_fp = NULL;	/* Abstract profile transform: */
			icc *abs_icc = NULL;
			xicc *abs_xicc = NULL;

			out_b2a_callback cx;

			if (verb)
				printf("Setting up B to A table lookup\n");

			if (ipname != NULL) {		/* There is a source profile to determine gamut mapping */

				/* Open up the profile for reading */
				if ((src_icco = read_embeded_icc(ipname)) == NULL)
					error ("Can't open file '%s'",ipname);

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
				xicc_enum_viewcond(x, vc, -1, NULL, 0, NULL);
		
				/* Override the viewing conditions */
				if (es >= 0)
					if (xicc_enum_viewcond(x, vc, es, NULL, 0, NULL) == -2)
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

			/* Get a suitable forward conversion object to invert. */
			/* By creating a separate one to the one created using scattered data, */
			/* we get the chance to set ICX_CAM_CLIP. */
			{
				int flags = 0;
	
				if (verb)
					flags |= ICX_VERBOSE;
	
				flags |= ICX_CLIP_NEAREST;		/* Not vector clip */
#ifdef USE_CAM_CLIP_OPT
				flags |= ICX_CAM_CLIP;			/* Clip in CAM Jab space rather than Lab */
#else
				warning("!!!! USE_CAM_CLIP_OPT in profout.c is off !!!!");
#endif

				if ((AtoB = wr_xicc->get_luobj(wr_xicc, flags, icmFwd,
				                  !allintents ? icmDefaultIntent : icRelativeColorimetric,
				                  wantLab ? icSigLabData : icSigXYZData,
                                  icmLuOrdNorm, &ovc, oink)) == NULL)
					error ("%d, %s",wr_xicc->errc, wr_xicc->err);
			}

			/* setup context ready for B2A table setting */
			cx.verb = verb;
			cx.pcsspace = wantLab ? icSigLabData : icSigXYZData;
#ifdef NO_B2A_PCS_CURVES
			cx.noPCScurves = 1;		/* Don't use PCS curves */
#else
			cx.noPCScurves = 0;
#endif
			cx.devspace = devspace;
			cx.x = (icxLuLut *)AtoB;		/* A2B icxLuLut created from scattered data */

			cx.ixp = NULL;		/* Perceptual PCS to CAM conversion */
			cx.ox = NULL;		/* CAM to PCS conversion */
			cx.pmap = NULL;		/* perceptual gamut map */
			cx.smap = NULL;		/* Saturation gamut map */

			cx.abs_luo = NULL;
			cx.xyzscale[0] = 1.0;
			cx.xyzscale[1] = 1.0;
			cx.gam = NULL;
			cx.wantLab = wantLab;			/* Copy PCS flag over */

			/* Open up the abstract profile if supplied, and setup luo */
			if (absname != NULL) {
				if ((abs_fp = new_icmFileStd_name(absname,"r")) == NULL)
					error ("Can't open abstract profile file '%s'",absname);
				
				if ((abs_icc = new_icc()) == NULL)
					error ("Creation of Abstract profile ICC object failed");
		
				/* Read header etc. */
				if ((rv = abs_icc->read(abs_icc,abs_fp,0)) != 0)
					error ("%d, %s",rv,abs_icc->err);
		
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
				if ((cx.abs_luo = abs_xicc->get_luobj(abs_xicc, ICX_CLIP_NEAREST, icmFwd,
				                  cx.abs_intent,
			        (cx.pcsspace == icSigLabData && cx.abs_intent == icRelativeColorimetric)
					             ? icSigLabData : icSigXYZData,
					icmLuOrdNorm, NULL, NULL)) == NULL)
					error ("%d, %s",abs_icc->errc, abs_icc->err);
			}
			
			if (!allintents) {	/* Only B2A0, no intent */
				if ((wo[0] = (icmLut *)wr_icco->read_tag(
				           wr_icco, icSigBToA0Tag)) == NULL) 
					error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
				cx.ntables = 1;

			} else {		/* All 3 intent tables */
				/* Intent 1 = relative colorimetric */
				if ((wo[0] = (icmLut *)wr_icco->read_tag(
				           wr_icco, icSigBToA1Tag)) == NULL) 
					error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
				cx.ntables = 1;

				if (src_xicc) {		/* Creating separate perceptual and Saturation tables */
					icRenderingIntent intentp;	/* Gamut mapping space perceptual selection */
					icRenderingIntent intents;	/* Gamut mapping space saturation selection */
					icRenderingIntent intento;	/* Gamut mapping space output selection */
					gamut *csgamp;			/* Incoming colorspace perceptual gamut */
					gamut *csgams;			/* Incoming colorspace saturation gamut */
					gamut *igam = NULL;		/* Incoming image gamut */
					gamut *ogam;			/* Destination colorspace gamut */
					double gres;			/* Gamut surface feature resolution */
					int    mapres;			/* Mapping rspl resolution */
			
					if (verb)
						printf("Creating Gamut Mapping\n");
			
					/* Gamut mapping will extend given grid res to encompas */
					/* source gamut by a margin. */
					if (oquality == 3) {	/* Ultra High */
			  	 		gres = 8.0;
			  	 		mapres = 41;
					} else if (oquality == 2) {	/* High */
			  	 		gres = 8.0;
			  	 		mapres = 33;
					} else if (oquality == 1) {	/* Medium */
			  	 		gres = 10.0;
			  	 		mapres = 25;
					} else if (oquality == 0) {	/* Low quality */
			  	 		gres = 12.0;
			  	 		mapres = 17;
					} else {					/* Extremely low */
			  	 		gres = 14.0;
			  	 		mapres = 9;
					}

					/* We could lift this restriction by allowing for separate */
					/* cx.ix and cx.ox for each intent, but that would be expensive!. */
					if (sepsat && (pgmi->usecas & 0xff) != (sgmi->usecas & 0xff))
						error("Can't handle percept and sat table intents with different CAM spaces");
					/* Default perceptual input gamut mapping space is absolute perceptual */
					intentp = noptop ? icAbsoluteColorimetric : icmAbsolutePerceptual;

					/* But override this for apperance space gamut mapping */
					if ((pgmi->usecas & 0xff) != 0x0) {
						intentp = noptop ? icxAppearance : icxPerceptualAppearance;
					}
					if (sepsat) {
						/* Default saturation gamut mapping space is absolute saturation */
						intents = nostos ? icAbsoluteColorimetric : icmAbsoluteSaturation;

						/* But override this for apperance space gamut mapping */
						if ((sgmi->usecas & 0xff) != 0x0) {
							intents = nostos ? icxAppearance : icxSaturationAppearance;
						}
					}
					/* Default output gamut mapping space is absolute colorimetric */
					intento = icAbsoluteColorimetric;

					/* But override this for apperance space gamut mapping */
					if ((pgmi->usecas & 0xff) != 0x0) {
						intento = icxAppearance;
					}

					if ((pgmi->usecas & 0xff) == 0x2) {
						double mxw;
						intentp = intents = intento = icxAbsAppearance;

						/* Make absolute common white point average between the two */
						ivc.Wxyz[0] = 0.5 * (ivc.Wxyz[0] + ovc.Wxyz[0]);
						ivc.Wxyz[1] = 0.5 * (ivc.Wxyz[1] + ovc.Wxyz[1]);
						ivc.Wxyz[2] = 0.5 * (ivc.Wxyz[2] + ovc.Wxyz[2]);
		
						/* And scale it's Y to be equal to 1.0 */
						mxw = 1.0/ivc.Wxyz[1];
						ivc.Wxyz[0] *= mxw;
						ivc.Wxyz[1] *= mxw;
						ivc.Wxyz[2] *= mxw;
		
						/* set output view conditions the same as the input */
						ovc = ivc;		/* Structure copy */
					}

					/* Unlike icclink, we've not provided a way for the user */
					/* to set the source profile ink limit, so estimate it */
					/* from the profile */
					icxDefaultLimits(src_xicc, &iink.tlimit, iink.tlimit,
					                           &iink.klimit, iink.klimit);

					/* Get lookup object simply for fwd_relpcs_outpcs() */
					/* and perceptual input gamut shell creation. */
					/* Note that the intent=Appearance will trigger Jab CAM, */
					/* overriding icSigLabData.. */
#ifdef NEVER
					printf("~1 input space flags = 0x%x\n",ICX_CLIP_NEAREST);
					printf("~1 input space intent = %s\n",icx2str(icmRenderingIntent,intentp));
					printf("~1 input space pcs = %s\n",icx2str(icmColorSpaceSignature,icSigLabData));
					printf("~1 input space viewing conditions =\n"); xicc_dump_viewcond(&ivc);
					printf("~1 input space inking =\n"); xicc_dump_inking(&iink);
#endif
					if ((cx.ixp = src_xicc->get_luobj(src_xicc,ICX_CLIP_NEAREST
					       , icmFwd, intentp, icSigLabData, icmLuOrdNorm, &ivc, &iink)) == NULL)
						error ("%d, %s",src_xicc->errc, src_xicc->err);

					/* Create the source colorspace gamut surface */
					if (verb)
						printf(" Finding Source Colorspace Perceptual Gamut\n");
			
					if ((csgamp = cx.ixp->get_gamut(cx.ixp, gres)) == NULL)
						error ("%d, %s",src_xicc->errc, src_xicc->err);
			
					if (sepsat) {
						icxLuBase *ixs = NULL;	/* Source profile saturation lookup for gamut */
						/* Get lookup object for saturation input gamut shell creation */
						/* Note that the intent=Appearance will trigger Jab CAM, */
						/* overriding icSigLabData.. */
						if ((ixs = src_xicc->get_luobj(src_xicc, ICX_CLIP_NEAREST
						   , icmFwd, intents, icSigLabData, icmLuOrdNorm, &ivc, &iink)) == NULL)
						error ("%d, %s",src_xicc->errc, src_xicc->err);

						if (verb)
							printf(" Finding Source Colorspace Saturation Gamut\n");

						if ((csgams = ixs->get_gamut(ixs, gres)) == NULL)
							error ("%d, %s",src_xicc->errc, src_xicc->err);
						ixs->del(ixs);
					}

					/* Read image source gamut if provided */
					if (sgname != NULL) {	/* Optional source gamut - ie. from an images */
						int isJab = 0;
			
						if ((pgmi->usecas & 0xff) != 0)
							isJab = 1;

						if (verb)
							printf(" Loading Image Source Gamut '%s'\n",sgname);
			
						igam = new_gamut(gres, 0);
			
						if (igam->read_gam(igam, sgname))
							error("Reading source gamut '%s' failed",sgname);

						if (igam->getisjab(igam) != isJab) {
							/* Should really convert to/from Jab here! */
							warning("Image gamut is wrong colorspace for gamut mapping (Lab != Jab)");
							/* This will actually error in the gamut mapping code */
							/* Note that we're not checking relative/absolute colorspace here. */
							/* At the moment it's up to the user to get this right. */
						}
					}

					/* Get lookup object for bwd_outpcs_relpcs(), */
					/* and output gamut shell creation */
					/* Note that the intent=Appearance will trigger Jab CAM, */
					/* overriding icSigLabData.. */
					if ((cx.ox = wr_xicc->get_luobj(wr_xicc, ICX_CLIP_NEAREST
						   , icmFwd, intento,
					                  icSigLabData, icmLuOrdNorm, &ovc, oink)) == NULL)
						error ("%d, %s",wr_xicc->errc, wr_xicc->err);

					/* Creat the destination gamut surface */
					if (verb)
						printf(" Finding Destination Gamut\n");
			
					if ((ogam = cx.ox->get_gamut(cx.ox, gres)) == NULL)
						error ("%d, %s",wr_xicc->errc, wr_xicc->err);
			
					if (verb)
						printf(" Creating Gamut match\n");

					/* The real range of Lab 0..100,-128..128,1-28..128 cube */
					/* when mapped to CAM is ridiculously large (ie. */
					/* 0..100, -288..265, -112..533), so we don't attempt to */
					/* set a gamut mapping grid range based on this. Instead */
					/* rely on the gamut map code to set a reasonable grid range */
					/* around the source gamut, and to cope reasonably with */
					/* values outside the grid range. */

					/* setup perceptual gamut mapping */
					cx.pmap = new_gammap(verb, csgamp, igam, ogam, pgmi, 0, 0, 0, 0, mapres,
					                     NULL, NULL, gamdiag ? "gammap_p.wrl" : NULL
					);
					if (cx.pmap == NULL)
						error ("Failed to make perceptual gamut map transform");

					/* Intent 0 = perceptual */
					if ((wo[1] = (icmLut *)wr_icco->read_tag(
					           wr_icco, icSigBToA0Tag)) == NULL) 
						error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
					cx.ntables = 2;

					if (sepsat) {
						/* setup saturation gamut mapping */
						cx.smap = new_gammap(verb, csgams, igam, ogam, sgmi, 0, 0, 0, 0, mapres,
						                     NULL, NULL, gamdiag ? "gammap_s.wrl" : NULL
						);
						if (cx.smap == NULL)
							error ("Failed to make saturation gamut map transform");

						/* Intent 2 = saturation */
						if ((wo[2] = (icmLut *)wr_icco->read_tag(
						           wr_icco, icSigBToA2Tag)) == NULL) 
							error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
						cx.ntables = 3;
					}
					csgamp->del(csgamp);
					csgamp = NULL;
					if (sepsat) {
						csgams->del(csgams);
						csgams = NULL;
					}
					if (igam != NULL) {
						igam->del(igam);
						igam = NULL;
					}
					ogam->del(ogam);
					ogam = NULL;
				}
			}
			cx.ochan = wo[0]->outputChan;

			/* If we've got a request for Absolute Appearance mode with scaling */
			/* to avoid clipping the source white point, compute the needed XYZ scaling factor. */
			if (src_xicc != NULL && allintents
			  && ((pgmi->usecas & 0x100) || (sgmi->usecas & 0x100))) {
				double swcam[3];
				double xyzscale[1], sa[1];

				/* Grab the source white point in CAM space */
				cx.ixp->efv_wh_bk_points(cx.ixp, swcam, NULL, NULL);

				/* Convert it to destination XYZ */
				/* We're being bad in delving inside the xluo, but we'll fix it latter */
				cx.ox->cam->cam_to_XYZ(cx.ox->cam, cx.swxyz, swcam);

//printf("~1 Source white Jab = %f %f %f\n", swcam[0], swcam[1], swcam[2]);
//printf("~1 Source white XYZ = %f %f %f\n", cx.swxyz[0], cx.swxyz[1], cx.swxyz[2]);

				/* Compute the bigest scale factor less than or equal to 1.0, */
				/* that doesn't clip the cx.swxyz[] on the destination gamut */
				sa[0] = 0.1;
				xyzscale[0] = 0.5;
				if (powell(NULL, 1, xyzscale, sa, 1e-6, 2000,
				                             xyzoptfunc, (void *)&cx, NULL, NULL) != 0) {
					warning("make_output_icc: XYZ scale powell failed to converge - set scale to 1.0");
				} else {
					if (pgmi->usecas & 0x100) {
						cx.xyzscale[0] = xyzscale[0];
						if (cx.verb) printf("Set Perceptual XYZ scale factor to %f\n",xyzscale[0]);
					}
					if (sgmi->usecas & 0x100) {
						cx.xyzscale[1] = xyzscale[0];
						if (cx.verb) printf("Set Saturation XYZ scale factor to %f\n",xyzscale[0]);
					}
				}
			}

// ====================================================================
#ifdef NEVER
// ~~99
			/* DEVELOPMENT CODE - not complete */
			/* Setup optimised B2A per channel curves */
			{
				xfit *xf;				/* Curve fitting class instance */
				int xfflags = 0;		/* xfit flags */
				double in_min[MXDI];	/* Input value scaling minimum */
				double in_max[MXDI];	/* Input value scaling maximum */
				double out_min[MXDO];	/* Output value scaling minimum */
				double out_max[MXDO];	/* Output value scaling maximum */
				int iluord, oluord;
				int iord[MXDI];			/* Input curve orders */
				int oord[MXDO];			/* Output curve orders */
				int nodp;				/* Number of inverse data points */
				co *points;				/* List of inverse points as PCS->dev */

				optcomb tcomb = oc_imo;	/* Create all by default */

				if ((xf = new_xfit()) == NULL) {
					error("profout: Creation of xfit object failed");
				}
					
				/* Setup for optimising run */
				if (nooluts)				/* Use option flags - swap sense for B2A */
					tcomb &= ~oc_i;

				if (noiluts)
					tcomb &= ~oc_o;

				if (verb)
					xfflags |= XFIT_VERB;

				xfflags |= XFIT_OUT_DEV;			/* Outpupt is device */
													/* (Switch to XFIT_OUT_LU latter ?) */
				if (cx.pcsspace == icSigLabData)
					xfflags |= XFIT_IN_ZERO;		/* Adjust a & b to zero */

~~~~~~~~~~~~~~~~~~~
				/* Set the curve order for input (PCS) */
				if (oquality >= 3) {				/* Ultra high */
					iluord = 25;			
					nodp = 2000;
				} else if (oquality == 2) {		/* High */
					iluord = 20;			
					nodp = 1000;
				} else if (oquality == 1) {		/* Medium */
					iluord = 17;			
					nodp = 750;
				} else {						/* Low */
					iluord = 10;			
					nodp = 500;
				}
				for (e = 0; e < p->inputChan; e++) {
					iord[e] = iluord;
					in_min[e] = p->inmin[e];
					in_max[e] = p->inmax[e];

					/* Hack to prevent a convex L curve pushing */
					/* the clut L values above the maximum value */
					/* that can be represented, causing clipping. */
					/* Do this by making sure that the L curve pivots */
					/* through 100.0 to 100.0 */
					if (e == 0 && cx.pcsspace == icSigLabData) {
						if (in_min[e] < 0.0001 && in_max[e] > 100.0) {
							in_max[e] = 100.0;	
						}
					}
				}

				/* Set curve order for output (Device) */
				if (oquality >= 3) {			/* Ultra high */
					oluord = 25;			
				} else if (oquality == 2) {		/* High */
					oluord = 20;			
				} else if (oquality == 1) {		/* Medium */
					oluord = 17;			
				} else {						/* Low */
					oluord = 10;			
				}
				for (f = 0; f < p->outputChan; f++) {
					oord[f] = oluord;
					out_min[f] = p->outmin[f];
					out_max[f] = p->outmax[f];

				}

				/* Create the sample points */
				~~~~~~~~~~~
				malloc

				for (i = 0; i < nodp;) {
					generate random pcs value

					if not within gamut
						continue;

					lookup through overall conversion
					i++;
				}

				/* Fit input and output curves to our data points */
				if (xf->fit(xf, xfflags, p->inputChan, p->outputChan, nodp, points,
				            in_min, in_max, out_min, out_max, iord, oord, tcomb, NULL) != 0) {
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
						p->del((icxLuBase *)p);
						return NULL;
					}

					p->inputTable[e]->set_rspl(p->inputTable[e], RSPL_NOFLAGS,
					           (void *)&cx, set_linfunc,
    			       &p->ninmin[e], &p->ninmax[e],
					           &p->lut->inputEnt,
					           &p->ninmin[e], &p->ninmax[e]);
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

					cx.xf = xf;
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

					p->outputTable[f]->set_rspl(p->outputTable[f], RSPL_NOFLAGS,
				           (void *)&cx, set_linfunc,
							min, max, &entries, min, max);

				}

				xf->del(xf);
			}
#endif /* NEVER (Setup optimised B2A per channel curves) */
// ====================================================================

			/* We now setup an exact inverse, colorimetric style, plus gamut mapping */
			/* for perceptual and saturation intents */
			/* Use helper function to do the hard work. */

			if (cx.verb) {
				unsigned int ui;
				int extra;
				cx.count = 0;
				cx.last = -1;
				for (cx.total = 1, ui = 0; ui < wo[0]->inputChan; ui++, cx.total *= wo[0]->clutPoints)
					; 
				/* Add in cell center points */
				for (extra = 1, ui = 0; ui < wo[0]->inputChan; ui++, extra *= (wo[0]->clutPoints-1))
					;
				cx.total += extra;
				printf("Creating B to A tables\n");
				printf(" 0%%"); fflush(stdout);
			}

#ifdef DEBUG_ONE
#define DBGNO 1		/* Up to 10 */

			/* Test a single given PCS (Rel D50 Lab) -> cmyk value */
			{
				double in[10][MAX_CHAN];
				double out[MAX_CHAN];
				in[0][0] = 29.565320;			/* sRGB blue */
				in[0][1] = 68.301300;
				in[0][2] = -112.049899;

				for (i = 0; i < DBGNO; i++) {
					printf("Input %f %f %f %f\n",in[i][0], in[i][1], in[i][2], in[i][3]);
					out_b2a_input((void *)&cx, out, in[i]);
					out_b2a_clut((void *)&cx, out, out);
					out_b2a_output((void *)&cx, out, out);
					printf("Output %f %f %f %f\n\n",out[0], out[1], out[2], out[3]);
				}
			}
#else /* !DEBUG_ONE */

#ifndef USE_LEASTSQUARES_APROX
			fprintf(stderr,"!!!!! profile/profout: USE_LEASTSQUARES_APROX undef !!!!!\n");
#endif
			if (icmSetMultiLutTables(
			        cx.ntables,
			        wo,
#ifdef USE_LEASTSQUARES_APROX
					ICM_CLUT_SET_APXLS,
#else
					0,
#endif
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

#ifdef WARN_CLUT_CLIPPING	/* Print warning if setting clut clips */
			/* Ignore clipping of the input table, because this happens */
			/* anyway due to Lab symetry adjustment. */
			if (wr_icco->warnc != 0 && wr_icco->warnc != 1) {
				warning("Values clipped in setting B2A LUT!");
			}
#endif /* WARN_CLUT_CLIPPING */
#endif /* !DEBUG_ONE */

			if (cx.abs_luo != NULL) {		/* Free up abstract transform */
				cx.abs_luo->del(cx.abs_luo);
				abs_xicc->del(abs_xicc);
				abs_icc->del(abs_icc);
				abs_fp->del(abs_fp);
			}

			if (cx.pmap != NULL)
				cx.pmap->del(cx.pmap), cx.pmap = NULL;
			if (cx.smap != NULL)
				cx.smap->del(cx.smap), cx.smap = NULL;
			if (cx.ixp != NULL)
				cx.ixp->del(cx.ixp), cx.ixp = NULL;
			if (cx.ox != NULL)
				cx.ox->del(cx.ox), cx.ox = NULL;

			if (src_xicc != NULL)
				src_xicc->del(src_xicc), src_xicc = NULL;
			if (src_icco != NULL)
				src_icco->del(src_icco), src_icco = NULL;

			if (verb)
				printf("Done B to A tables\n");
		}

		/* Set the ColorantTable PCS values */
		{
			unsigned int i;
			icmColorantTable *wo;
			double dv[MAX_CHAN];

			if ((wo = (icmColorantTable *)wr_icco->read_tag(
			           wr_icco, icSigColorantTableTag)) == NULL) 
				error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
	
			for (i = 0; i < wo->count; i++)
				dv[i] = 0.0;

			/* Lookup the colorant PCS values the recommended ICC way */
			for (i = 0; i < wo->count; i++) {
				dv[i] = 1.0;
				AtoB->lookup(AtoB, wo->data[i].pcsCoords, dv);
				dv[i] = 0.0;
			}
		}

#ifndef DISABLE_GAMUT_TAG
		/* Create Gamut clut for output type */
		/* This is not mandated for V2.4.0 Display profiles, but add it anyway */
		if (allintents) {	
			icmLut *wo;
			out_b2a_callback cx;
			double gres = 0.0;

			cx.verb = verb;
			cx.pcsspace = wantLab ? icSigLabData : icSigXYZData;
			cx.devspace = devspace;
			cx.x = (icxLuLut *)AtoB;		/* A2B icxLuLut */

			if (verb)
				printf("Creating gamut boundary table\n");

			/* Need to switch AtoB to be override Lab PCS */
			/* Do this the dirty way, by delving into xicclu and icclu. Alternatively */
			/* we could create an xlut set method, delete AtoB and recreate it, */
			/* or fix get_gamut to independently override convert to icmSigLabData */
			/* ~~~~~~~~999 should fix this !!! */
			cx.x->outs = icSigLabData;
			cx.x->pcs = icSigLabData;
			cx.x->plu->e_outSpace = icSigLabData;
			cx.x->plu->e_pcs = icSigLabData;
			cx.wantLab = wantLab;			/* Copy PCS flag over */

			if (oquality == 3) {	/* Ultra High */
	  	 		gres = 8.0;
			} else if (oquality == 2) {	/* High */
	  	 		gres = 9.0;
			} else if (oquality == 1) {	/* Medium */
	  	 		gres = 10.0;
			} else if (oquality == 0) {	/* Low quality */
	  	 		gres = 12.0;
			} else {					/* Extremely low */
	  	 		gres = 15.0;
			}

			/* Creat a gamut surface */
			if ((cx.gam = AtoB->get_gamut(AtoB, gres)) == NULL)
				error("Get_gamut failed: %d, %s",AtoB->pp->errc,AtoB->pp->err);
			
			if ((wo = (icmLut *)wr_icco->read_tag(
			           wr_icco, icSigGamutTag)) == NULL) 
				error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

			if (cx.verb) {
				unsigned int ui;
				cx.count = 0;
				cx.last = -1;
				for (cx.total = 1, ui = 0; ui < wo->inputChan; ui++, cx.total *= wo->clutPoints)
					; 
				printf(" 0%%"); fflush(stdout);
			}
#ifndef DEBUG_ONE	/* Skip this when debugging */
			if (wo->set_tables(wo,
					ICM_CLUT_SET_EXACT,
					&cx,				/* Context */
					cx.pcsspace,		/* Input color space */
					icSigGrayData,		/* Output color space */
					out_b2a_input,		/* Input transform PCS->PCS' */
					NULL, NULL,			/* Use default Lab' range */
					PCSp_bdist,			/* Lab' -> Boundary distance */
					NULL, NULL,			/* Use default Device' range */
					gamut_output) != 0)	/* Boundary distance to out of gamut value */
				error("Setting 16 bit PCS->Device Gamut Lut failed: %d, %s",wr_icco->errc,wr_icco->err);
#endif /* !DEBUG_ONE */
			if (cx.verb) {
				printf("\n");
			}
#ifdef WARN_CLUT_CLIPPING	/* Print warning if setting clut clips */
			if (wr_icco->warnc)
				warning("Values clipped in setting Gamut LUT");
#endif /* WARN_CLUT_CLIPPING */

			cx.gam->del(cx.gam);		/* Done with gamut object */
			cx.gam = NULL;

			if (verb)
				printf("Done gamut boundary table\n");
		}
#endif /* !DISABLE_GAMUT_TAG */

		/* Free up xicc stuff */
		AtoB->del(AtoB); /* Done with device to PCS lookup */
		wr_xicc->del(wr_xicc);

	}
	/* Gamma/Shaper + matrix profile */
	/* or XYZ cLUT with matrix as well. */
	if (!isLut || mtxtoo) {
		xicc *wr_xicc;			/* extention object */
		icxLuBase *xluo;		/* Forward ixcLu */
		int flags = 0;

		/* Wrap with an expanded icc */
		if ((wr_xicc = new_xicc(wr_icco)) == NULL)
			error("Creation of xicc failed");
		
		if (verb)
			flags |= ICX_VERBOSE;

		if (!mtxtoo)	/* Write matrix white/black/Luminance if no cLUT */
			flags |= ICX_WRITE_WBL;

		/* Setup Device -> XYZ conversion (Fwd) object from scattered data. */
		if ((xluo = wr_xicc->set_luobj(
		               wr_xicc, icmFwd, isdisp ? icmDefaultIntent : icRelativeColorimetric,
		               icmLuOrdRev,
			           flags | ICX_SET_WHITE | ICX_SET_BLACK, 		/* Compute white & black */
		               npat, tpat, dispLuminance, -1.0, smooth, avgdev,
		               NULL, oink, cal, iquality)) == NULL)
			error("%d, %s",wr_xicc->errc, wr_xicc->err);

		/* Free up xicc stuff */
		xluo->del(xluo);

		/* Set the ColorantTable PCS values */
		if (!mtxtoo) {
			unsigned int i;
			icmColorantTable *wo;
			double dv[MAX_CHAN];

			/* Get lookup object simply for fwd_relpcs_outpcs() */
			if ((xluo = wr_xicc->get_luobj(wr_xicc, ICX_CLIP_NEAREST, icmFwd,
			                     icRelativeColorimetric, icmSigDefaultData,
			                     icmLuOrdNorm, NULL, NULL)) == NULL)
				error ("%d, %s",wr_xicc->errc, wr_xicc->err);

			if ((wo = (icmColorantTable *)wr_icco->read_tag(
			           wr_icco, icSigColorantTableTag)) == NULL) 
				error("read_tag failed: %d, %s",wr_icco->errc,wr_icco->err);
	
			for (i = 0; i < wo->count; i++)
				dv[i] = 0.0;

			/* Lookup the colorant PCS values the recommended ICC way */
			for (i = 0; i < wo->count; i++) {
				dv[i] = 1.0;
				xluo->lookup(xluo, wo->data[i].pcsCoords, dv);
				/* Matrix profile can produce -ve values not representable by 16 bit XYZ */
				icmClipXYZ(wo->data[i].pcsCoords,wo->data[i].pcsCoords);
				dv[i] = 0.0;
			}
			xluo->del(xluo);
		}

		wr_xicc->del(wr_xicc);
	}

	/* We're done with any cal now */
	if (cal != NULL)
		cal->del(cal);

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
		double rerr = 0.0;
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

			rerr += mxd * mxd;
			aerr += mxd;
			nsamps++;
		}
		rerr = sqrt(rerr/nsamps);
		aerr /= nsamps;
		printf("profile check complete, peak err = %f, avg err = %f, RMS = %f\n",merr,aerr,rerr);

		/* Done with lookup object */
		luo->del(luo);

		/* Close the file */
		rd_icco->del(rd_icco);
		rd_fp->del(rd_fp);
	}

	free(tpat);
}

