
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
 * Based on the old iccXfm class.
 */

/*
 * This module expands the basic icclib functionality,
 * providing more functionality in exercising.
 * The implementation for the three different types
 * of profile representation, are in their own source files.
 */

/*
 * TTBD:
 *       Some of the error handling is crude. Shouldn't use
 *       error(), should return status.
 *
 */

#include <sys/types.h>
#include <string.h>
#include <ctype.h>
#ifdef __sun
#include <unistd.h>
#endif
#if defined(__IBMC__) && defined(_M_IX86)
#include <float.h>
#endif
#include "numlib.h"
#include "counters.h"
#include "plot.h"
#include "../h/sort.h"
#include "xicc.h"		/* definitions for this library */

#define USE_CAM			/* Use CIECAM02 for clipping and gamut mapping, else use Lab */

#undef DEBUG			/* Plot 1d Luts */

#ifdef DEBUG
#include "plot.h"
#endif

#define MAX_INVSOLN 4

static void xicc_free(xicc *p);
icxLuBase * xicc_get_luobj(xicc *p, int flags, icmLookupFunc func, icRenderingIntent intent,
                           icColorSpaceSignature pcsor, icmLookupOrder order,
                           icxViewCond *vc, icxInk *ink);
static icxLuBase *xicc_set_luobj(xicc *p, icmLookupFunc func, icRenderingIntent intent,
                            icmLookupOrder order, int flags, int no, cow *points,
                            double dispLuminance, double smooth, double avgdev,
                            icxViewCond *vc, icxInk *ink, int quality);
static void icxLutSpaces(icxLuBase *p, icColorSpaceSignature *ins, int *inn,
                         icColorSpaceSignature *outs, int *outn,
                         icColorSpaceSignature *pcs);
static void icxLuSpaces(icxLuBase *p, icColorSpaceSignature *ins, int *inn,
                        icColorSpaceSignature *outs, int *outn, 
                        icmLuAlgType *alg, icRenderingIntent *intt,
                        icmLookupFunc *fnc, icColorSpaceSignature *pcs);
static void icxLu_get_native_ranges (icxLuBase *p,
                              double *inmin, double *inmax, double *outmin, double *outmax);
static void icxLu_get_ranges (icxLuBase *p,
                              double *inmin, double *inmax, double *outmin, double *outmax);
static void icxLuRel_wh_bk_points(icxLuBase *p, double *wht, double *blk);
int xicc_get_viewcond(xicc *p, icxViewCond *vc);

/* The different profile types are in their own source filesm */
/* and are included to keep their functions private. (static) */
#include "xmono.c"
#include "xmatrix.c"
#include "xlut.c"		/* New xfit3 in & out optimising based profiles */
//#include "xlut1.c"	/* Old xfit1 device curve based profiles */

#ifdef NT		/* You'd think there might be some standards.... */
# ifndef __BORLANDC__
#  define stricmp _stricmp
# endif
#else
# define stricmp strcasecmp
#endif
/* Utilities */

/* Return a string description of the given enumeration value */
const char *icx2str(icmEnumType etype, int enumval) {

	if (etype == icmColorSpaceSignature) {
		if (((icColorSpaceSignature)enumval) == icxSigJabData)
			return "Jab";
		else if (((icColorSpaceSignature)enumval) == icxSigJChData)
			return "JCh";
		else if (((icColorSpaceSignature)enumval) == icxSigLChData)
			return "LCh";
	} else if (etype == icmRenderingIntent) {
		if (((icRenderingIntent)enumval) == icxAppearance)
			return "icxAppearance";
		else if (((icRenderingIntent)enumval) == icxAbsAppearance)
			return "icxAbsAppearance";
		else if (((icRenderingIntent)enumval) == icxPerceptualAppearance)
			return "icxPerceptualAppearance";
		else if (((icRenderingIntent)enumval) == icxAbsPerceptualAppearance)
			return "icxAbsPerceptualAppearance";
		else if (((icRenderingIntent)enumval) == icxSaturationAppearance)
			return "icxSaturationAppearance";
		else if (((icRenderingIntent)enumval) == icxAbsSaturationAppearance)
			return "icxAbsSaturationAppearance";
	}
	return icm2str(etype, enumval);
}

/* Common xicc stuff */

/* Return information about the native lut in/out colorspaces. */
/* Any pointer may be NULL if value is not to be returned */
static void
icxLutSpaces(
	icxLuBase *p,					/* This */
	icColorSpaceSignature *ins,		/* Return input color space */
	int *inn,						/* Return number of input components */
	icColorSpaceSignature *outs,	/* Return output color space */
	int *outn,						/* Return number of output components */
	icColorSpaceSignature *pcs		/* Return PCS color space */
) {
	p->plu->lutspaces(p->plu, ins, inn, outs, outn, pcs);
}

/* Return information about the overall lookup in/out colorspaces, */
/* including allowance for any PCS override. */
/* Any pointer may be NULL if value is not to be returned */
static void
icxLuSpaces(
	icxLuBase *p,					/* This */
	icColorSpaceSignature *ins,		/* Return input color space */
	int *inn,						/* Return number of input components */
	icColorSpaceSignature *outs,	/* Return output color space */
	int *outn,						/* Return number of output components */
	icmLuAlgType *alg,				/* Return type of lookup algorithm used */
    icRenderingIntent *intt,		/* Return the intent implemented */
    icmLookupFunc *fnc,				/* Return the profile function being implemented */
	icColorSpaceSignature *pcs		/* Return the effective PCS */
) {
    icmLookupFunc function;
	icColorSpaceSignature npcs;		/* Native PCS */

	p->plu->spaces(p->plu, NULL, inn, NULL, outn, alg, NULL, &function, &npcs);

	if (intt != NULL)
		*intt = p->intent;

	if (fnc != NULL)
		*fnc = function;

	if (ins != NULL)
		*ins = p->ins;

	if (outs != NULL)
		*outs = p->outs;

	if (pcs != NULL)
		*pcs = p->pcs;
}

/* Return the native (internaly visible) colorspace value ranges */
static void
icxLu_get_native_ranges (
icxLuBase *p,
double *inmin, double *inmax,		/* Return maximum range of inspace values */
double *outmin, double *outmax		/* Return maximum range of outspace values */
) {
	int i;
	if (inmin != NULL) {
		for (i = 0; i < p->inputChan; i++)
			inmin[i] = p->ninmin[i];
	}
	if (inmax != NULL) {
		for (i = 0; i < p->inputChan; i++)
			inmax[i] = p->ninmax[i];
	}
	if (outmin != NULL) {
		for (i = 0; i < p->outputChan; i++)
			outmin[i] = p->noutmin[i];
	}
	if (outmax != NULL) {
		for (i = 0; i < p->outputChan; i++)
			outmax[i] = p->noutmax[i];
	}
}

/* Return the effective (externaly visible) colorspace value ranges */
static void
icxLu_get_ranges (
icxLuBase *p,
double *inmin, double *inmax,		/* Return maximum range of inspace values */
double *outmin, double *outmax		/* Return maximum range of outspace values */
) {
	int i;
	if (inmin != NULL) {
		for (i = 0; i < p->inputChan; i++)
			inmin[i] = p->inmin[i];
	}
	if (inmax != NULL) {
		for (i = 0; i < p->inputChan; i++)
			inmax[i] = p->inmax[i];
	}
	if (outmin != NULL) {
		for (i = 0; i < p->outputChan; i++)
			outmin[i] = p->outmin[i];
	}
	if (outmax != NULL) {
		for (i = 0; i < p->outputChan; i++)
			outmax[i] = p->outmax[i];
	}
}

/* Return the relative media white and black points */
/* in the xlu effective PCS colorspace. Pointers may be NULL. */
static void icxLuRel_wh_bk_points(
struct _icxLuBase *p,
double *wht,
double *blk
) {
	icmXYZNumber whiteXYZ, blackXYZ;
	double white[3], black[3];

	/* Get the absolute XYZ profile black and white points */
	p->plu->wh_bk_points(p->plu, &whiteXYZ, &blackXYZ);

	white[0] = whiteXYZ.X;
	white[1] = whiteXYZ.Y;
	white[2] = whiteXYZ.Z;
	black[0] = blackXYZ.X;
	black[1] = blackXYZ.Y;
	black[2] = blackXYZ.Z;

	/* Convert to relative colorimetric if appropriate */
	switch (p->intent) {
		/* If it is relative */
		case icmDefaultIntent:				/* Shouldn't happen */
		case icPerceptual:
		case icRelativeColorimetric:
		case icSaturation:
			p->plu->XYZ_Abs2Rel(p->plu,white,white);
			p->plu->XYZ_Abs2Rel(p->plu,black,black);
			break;
		default:
			break;
	}

	/* Convert to the (possibly) override PCS */
	switch (p->pcs) {
		case icSigXYZData:
			break;								/* Don't have to do anyting */
		case icSigLabData:
			icmXYZ2Lab(&icmD50, white, white);	/* Convert from XYZ to Lab */
			icmXYZ2Lab(&icmD50, black, black);
			break;
		case icxSigJabData:
			p->cam->XYZ_to_cam(p->cam, white, white);	/* Convert from XYZ to Jab */
			p->cam->XYZ_to_cam(p->cam, black, black);
			break;
		default:
			break;
	}
	if (wht != NULL) {
		wht[0] = white[0];
		wht[1] = white[1];
		wht[2] = white[2];
	}

	if (blk != NULL) {
		blk[0] = black[0];
		blk[1] = black[1];
		blk[2] = black[2];
	}
}

/* Create an instance of an xicc object */
xicc *new_xicc(
icc *picc		/* icc we are expanding */
) {
	xicc *p;
	if ((p = (xicc *) calloc(1,sizeof(xicc))) == NULL)
		return NULL;
	p->pp = picc;
	p->del           = xicc_free;
	p->get_luobj     = xicc_get_luobj;
	p->set_luobj     = xicc_set_luobj;
	p->get_viewcond  = xicc_get_viewcond;
	return p;
}

/* Do away with the xicc (but not the icc!) */
static void xicc_free(
xicc *p
) {
	free (p);
}


/* Return an expanded lookup object, initialised */
/* from the icc. */
/* Return NULL on error, check errc+err for reason. */
/* Set the pcsor & intent to consistent and values if */
/* Jab and/or icxAppearance has been requested. */
/* Create the underlying icm lookup object that is used */
/* to create and implement the icx one. The icm will be used */
/* to translate from native to effective PCS, unless the */
/* effective PCS is Jab, in which case the icm will be set to */
/* have an effective PCS of XYZ. Since native<->effecive PCS conversion */
/* is done at the to/from_abs() stage, none of this affects the individual */
/* conversion steps, which will all talk the native PCS (unless merged). */
icxLuBase *xicc_get_luobj(
xicc *p,					/* this */
int flags,					/* clip, merge flags */
icmLookupFunc func,			/* Functionality */
icRenderingIntent intent,	/* Intent */
icColorSpaceSignature pcsor,/* PCS override (0 = def) */
icmLookupOrder order,		/* Search Order */
icxViewCond *vc,			/* Viewing Condition (may be NULL if pcsor is not CIECAM) */
icxInk *ink					/* inking details (NULL for default) */
) {
	icmLuBase *plu;
	icxLuBase *xplu;
	icmLuAlgType alg;
	icRenderingIntent n_intent = intent;			/* Native Intent to request */
	icColorSpaceSignature n_pcs = icmSigDefaultData;	/* Native PCS to request */

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	/* Ensure that appropriate PCS is slected for an appearance intent */
	if (intent == icxAppearance
	 || intent == icxAbsAppearance
	 || intent == icxPerceptualAppearance
	 || intent == icxAbsPerceptualAppearance
	 || intent == icxSaturationAppearance
	 || intent == icxAbsSaturationAppearance) {
		pcsor = icxSigJabData;

	/* Translate non-Jab intents to the equivalent appearance "intent" if pcsor == Jab. */
	/* This is how we get these when the UI's don't list all the apperances intents, */
	/* we select the analogous non-apperance intent with pcsor = Jab. */
	/* Note that Abs/non-abs selects between Apperance and AbsAppearance. */
	} else if (pcsor == icxSigJabData) {
		if (intent == icRelativeColorimetric)
			intent = icxAppearance;
		else if (intent == icAbsoluteColorimetric)
			intent = icxAbsAppearance;
		else if (intent == icPerceptual)
			intent = icxPerceptualAppearance;
		else if (intent == icmAbsolutePerceptual)
			intent = icxAbsPerceptualAppearance;
		else if (intent == icSaturation)
			intent = icxSaturationAppearance;
		else if (intent == icmAbsoluteSaturation)
			intent = icxAbsSaturationAppearance;
		else
			intent = icxAppearance;
	}

	/* Translate intent asked for into intent needed icclib */
	if      (intent == icxAppearance
	      || intent == icxAbsAppearance)
		n_intent = icAbsoluteColorimetric;
	else if (intent == icxPerceptualAppearance
	      || intent == icxAbsPerceptualAppearance)
		n_intent = icmAbsolutePerceptual;
	else if (intent == icxSaturationAppearance
	      || intent == icxAbsSaturationAppearance)
		n_intent = icmAbsoluteSaturation;

	if (pcsor != icmSigDefaultData)
		n_pcs = pcsor;			/* There is an icclib override */

	if (pcsor == icxSigJabData)	/* xicc override */
		n_pcs = icSigXYZData;	/* Translate to XYZ */

	/* Get icclib lookup object */
	if ((plu = p->pp->get_luobj(p->pp, func, n_intent, n_pcs, order)) == NULL) {
		p->errc = p->pp->errc;		/* Copy error */
		strcpy(p->err, p->pp->err);
		return NULL;
	}

	/* Figure out what the algorithm is */
	plu->spaces(plu, NULL, NULL, NULL, NULL, &alg, NULL, NULL, &n_pcs);

	/* make sure its "Abs CAM" */
	if (vc!= NULL
	 && (intent == icxAbsAppearance
	  || intent == icxAbsPerceptualAppearance
	  || intent == icxAbsSaturationAppearance)) {	/* make sure its "Abs CAM" */
		/* Set white point and flare color to D50 */
		vc->Wxyz[0] = icmD50.X/icmD50.Y;
		vc->Wxyz[1] = icmD50.Y/icmD50.Y;	// Normalise white reference to Y = 1 ?
		vc->Wxyz[2] = icmD50.Z/icmD50.Y;
	
		vc->Fxyz[0] = icmD50.X;
		vc->Fxyz[1] = icmD50.Y;
		vc->Fxyz[2] = icmD50.Z;
	}

	/* Call xiccLu wrapper creation */
	switch (alg) {
		case icmMonoFwdType:
			xplu = new_icxLuMono(p, flags, plu, func, intent, pcsor, vc, 0);
			break;
		case icmMonoBwdType:
			xplu = new_icxLuMono(p, flags, plu, func, intent, pcsor, vc, 1);
			break;
    	case icmMatrixFwdType:
			xplu = new_icxLuMatrix(p, flags, plu, func, intent, pcsor, vc, 0);
			break;
    	case icmMatrixBwdType:
			xplu = new_icxLuMatrix(p, flags, plu, func, intent, pcsor, vc, 1);
			break;
    	case icmLutType:
			xplu = new_icxLuLut(p, flags, plu, func, intent, pcsor, vc, ink);
			break;
		default:
			xplu = NULL;
			break;
	}

	return xplu;
}


/* Return an expanded lookup object, initialised */
/* from the icc, and then overwritten by a conversion */
/* created from the supplied scattered data points. */
/* The Lut is assumed to be a device -> native PCS profile. */
/* If the SET_WHITE and/or SET_BLACK flags are set, */
/* discover the white/black point, set it in the icc, */
/* and make the Lut relative to them. */
/* Return NULL on error, check errc+err for reason */
static icxLuBase *xicc_set_luobj(
xicc *p,					/* this */
icmLookupFunc func,			/* Functionality */
icRenderingIntent intent,	/* Intent */
icmLookupOrder order,		/* Search Order */
int flags,					/* white/black point, verbose flags etc. */
int no,						/* Number of points */
cow *points,				/* Array of input points */
double dispLuminance,		/* > 0.0 if display luminance value and is known */
double smooth,				/* RSPL smoothing factor, -ve if raw */
double avgdev,				/* reading Average Deviation as a proportion of the input range */
icxViewCond *vc,			/* Viewing Condition (NULL if not using CAM) */
icxInk *ink,				/* inking details (NULL for default) */
int quality					/* Quality metric, 0..3 */
) {
	icmLuBase *plu;
	icxLuBase *xplu = NULL;
	icmLuAlgType alg;

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	if (func != icmFwd) {
		p->errc = 1;
		sprintf(p->err,"Can only create Device->PCS profiles from scattered data.");
		xplu = NULL;
		return xplu;
	}

	/* Get icclib lookup object */
	if ((plu = p->pp->get_luobj(p->pp, func, intent, 0, order)) == NULL) {
		p->errc = p->pp->errc;		/* Copy error */
		strcpy(p->err, p->pp->err);
		return NULL;
	}

	/* Figure out what the algorithm is */
	plu->spaces(plu, NULL, NULL, NULL, NULL, &alg, NULL, NULL, NULL);

	/* Call xiccLu wrapper creation */
	switch (alg) {
		case icmMonoFwdType:
			p->errc = 1;
			sprintf(p->err,"Setting Monochrome Fwd profile from scattered data not supported.");
			plu->del(plu);
			xplu = NULL;		/* Not supported yet */
			break;

    	case icmMatrixFwdType:
			xplu = set_icxLuMatrix(p, plu, flags, no, points, dispLuminance, quality);
			break;

    	case icmLutType:
			/* ~~~ Should add check that it is a fwd profile ~~~ */
			xplu = set_icxLuLut(p, plu, intent, func, flags, no, points, dispLuminance, smooth, avgdev, vc, ink, quality);
			break;

		default:
			break;
	}

	return xplu;
}

/* ------------------------------------------------------ */
/* Viewing Condition Parameter stuff */

#ifdef NEVER	/* Not currently used */

/* Guess viewing parameters from the technology signature */
static void guess_from_techsig(
icTechnologySignature tsig,
double *Ybp
) {
double Yb = -1.0;

	switch (tsig) {
		/* These are all inputing either a representation of */
		/* a natural scene captured on another medium, or are assuming */
		/* that the medium is the original. A _good_ system would */
		/* let the user indicate which is the case. */
		case icSigReflectiveScanner:
		case icSigFilmScanner:
			Yb = 0.2;
			break;

		/* Direct scene to value devices. */
		case icSigDigitalCamera:
		case icSigVideoCamera:
			Yb = 0.2;
			break;

		/* Emmisive displays. */
		/* We could try tweaking the white point on the assumption */
		/* that the viewer will be adapted to a combination of both */
		/* the CRT white point, and the ambient light. */
		case icSigVideoMonitor:
		case icSigCRTDisplay:
		case icSigPMDisplay:
		case icSigAMDisplay:
			Yb = 0.2;
			break;

		/* Photo CD has its own viewing definitions */
		/* (It represents original scene colors) */
		case icSigPhotoCD:
			Yb = 0.2;
			break;

		/* Projection devices, either direct, or */
		/* via another intermediate medium. */
		case icSigProjectionTelevision:
			Yb = 0.1;		/* Assume darkened room, little background */
			break;
		case icSigFilmWriter:
			Yb = 0.0;		/* Assume a dark room - no background */
			break;

		/* Printed media devices. */
		case icSigInkJetPrinter:
		case icSigThermalWaxPrinter:
		case icSigElectrophotographicPrinter:
		case icSigElectrostaticPrinter:
		case icSigDyeSublimationPrinter:
		case icSigPhotographicPaperPrinter:
		case icSigPhotoImageSetter:
		case icSigGravure:
		case icSigOffsetLithography:
		case icSigSilkscreen:
		case icSigFlexography:
			Yb = 0.2;
			break;

		default:
			Yb = 0.2;
	}

	if (Ybp != NULL)
		*Ybp = Yb;
}

#endif /* NEVER */


/* See if we can read or guess the viewing conditions */
/* for an ICC profile. */
/* Return value 0 if it is well defined */
/* Return value 1 if it is a guess */
/* Return value 2 if it is not possible/appropriate */
int xicc_get_viewcond(
xicc *p,			/* Expanded profile we're working with */
icxViewCond *vc		/* Viewing parameters to return */
) {
	icc *pp = p->pp;	/* Base ICC */

	/* Numbers we're trying to find */
	ViewingCondition Ev = vc_none;
	double Wxyz[3] = {-1.0, -1.0, -1.0};	/* Adapting white color */
	double La = -1.0;		/* Adapting luminance */
	double Ixyz[3] = {-1.0, -1.0, -1.0};	/* Illuminant color */
	double Li = -1.0;						/* Illuminant luminance */
	double Lb = -1.0;		/* Backgrount luminance */
	double Yb = -1.0;		/* Background relative luminance to Lv */
	double Lve = -1.0;		/* Emissive device image luminance */
	double Lvr = -1.0;		/* Reflective device image luminance */
	double Lv = -1.0;		/* device image luminance */
	double Yf = -1.0;		/* Flare relative luminance to Lv */
	double Fxyz[3] = {-1.0, -1.0, -1.0};	/* Flare color */
    icTechnologySignature tsig = icMaxEnumTechnology; /* Technology Signature */
    icProfileClassSignature devc = icMaxEnumClass;
	int trans = -1;		/* Set to 0 if not transparency, 1 if it is */

	/* Collect all the information we can find */

	/* Emmisive devices image white luminance */
	{
		icmXYZArray *luminanceTag;

		if ((luminanceTag = (icmXYZArray *)pp->read_tag(pp, icSigLuminanceTag)) != NULL
         && luminanceTag->ttype == icSigXYZType && luminanceTag->size >= 1) {
			Lve = luminanceTag->data[0].Y;	/* Copy structure */
		} 
	}

	/* Flare: */
	{
		icmMeasurement *ro;

		if ((ro = (icmMeasurement *)pp->read_tag(pp, icSigMeasurementTag)) != NULL
		 && ro->ttype == icSigMeasurementType) {
	
    		Yf = ro->flare;
   			/* ro->illuminant ie D50, D65, D93, A etc. */
		}
	}

	/* Media White Point */
	{
		icmXYZArray *whitePointTag;

		if ((whitePointTag = (icmXYZArray *)pp->read_tag(pp, icSigMediaWhitePointTag)) != NULL
         && whitePointTag->ttype == icSigXYZType && whitePointTag->size >= 1) {
			Wxyz[0] = whitePointTag->data[0].X;
			Wxyz[1] = whitePointTag->data[0].Y;
			Wxyz[2] = whitePointTag->data[0].Z;
		}
	}

	/* ViewingConditions: */
	{
		icmViewingConditions *ro;

		if ((ro = (icmViewingConditions *)pp->read_tag(pp, icSigViewingConditionsTag)) != NULL
		 && ro->ttype == icSigViewingConditionsType) {
	
   			/* ro->illuminant.X */
   			/* ro->illuminant.Z */

    		Li = ro->illuminant.Y;

			/* Reflect illuminant off the media white */
    		Lvr = Li * Wxyz[1];

			/* Illuminant color */
			Ixyz[0] = ro->illuminant.X/ro->illuminant.Y;
			Ixyz[1] = 1.0;
			Ixyz[2] = ro->illuminant.Z/ro->illuminant.Y;
			
			/* Assume ICC surround is CICAM97 background */
    		/* ro->surround.X */
    		/* ro->surround.Z */
    		La = ro->surround.Y;

    		/* ro->stdIlluminant ie D50, D65, D93, A etc. */
		}
	}

	/* Stuff we might need */

	/* Technology: */
	{
		icmSignature *ro;

		/* Try and read the tag from the file */
		if ((ro = (icmSignature *)pp->read_tag(pp, icSigTechnologyTag)) != NULL
		 && ro->ttype != icSigSignatureType) {

		tsig = ro->sig;
		}
	}

    devc = pp->header->deviceClass;	/* Type of profile */
	if (devc == icSigLinkClass
	 || devc == icSigAbstractClass
	 || devc == icSigColorSpaceClass
	 || devc == icSigNamedColorClass)
		return 2;

	/*
    	icSigInputClass
    	icSigDisplayClass
    	icSigOutputClass
	*/

	if ((pp->header->flags & icTransparency) != 0)
		trans = 1;
	else
		trans = 0;


	/* figure Lv if we have the information */
	if (Lve >= 0.0)
		Lv = Lve;		/* Emmisive image white luminance */
	else
		Lv = Lvr;		/* Reflectance image white luminance */

	/* Fudge the technology signature */
	if (tsig == icMaxEnumTechnology) {
		if (devc == icSigDisplayClass)
			tsig = icSigCRTDisplay;
	}

#ifndef NEVER
	printf("Enumeration = %d\n", Ev);
	printf("Viewing Conditions:\n");
	printf("White adaptation color %f %f %f\n",Wxyz[0], Wxyz[1], Wxyz[2]);
	printf("Adapting Luminance La = %f\n",La);
	printf("Illuminant color %f %f %f\n",Ixyz[0], Ixyz[1], Ixyz[2]);
	printf("Illuminant Luminance Li = %f\n",Li);
	printf("Background Luminance Lb = %f\n",Lb);
	printf("Relative Background Yb = %f\n",Yb);
	printf("Emissive Image White Lve = %f\n",Lve);
	printf("Reflective Image White Lvr = %f\n",Lvr);
	printf("Device Image White Lv = %f\n",Lv);
	printf("Relative Flare Yf = %f\n",Yf);
	printf("Flare color %f %f %f\n",Fxyz[0], Fxyz[1], Fxyz[2]);
	printf("Technology = %s\n",tag2str(tsig));
	printf("deviceClass = %s\n",tag2str(devc));
	printf("Transparency = %d\n",trans);
#endif
	
	/* See if the viewing conditions are completely defined as ICC can do it */
	if (Wxyz[0] >= 0.0 && Wxyz[1] >= 0.0 && Wxyz[2] >= 0.0
	 && La >= 0.0
	 && Yb >= 0.0
	 && Lv >= 0.0
	 && Yf >= 0.0
	 && Fxyz[0] >= 0.0 && Fxyz[1] >= 0.0 && Fxyz[2] >= 0.0) {

		vc->Ev = vc_none;
		vc->Wxyz[0] = Wxyz[0];
		vc->Wxyz[1] = Wxyz[1];
		vc->Wxyz[2] = Wxyz[2];
		vc->La = La;
		vc->Yb = Yb;
		vc->Lv = Lv;
		vc->Yf = Yf;
		vc->Fxyz[0] = Fxyz[0];	
		vc->Fxyz[1] = Fxyz[1];	
		vc->Fxyz[2] = Fxyz[2];	
		return 0;
	}

	/* Hmm. We didn't get all the info an ICC can contain. */
	/* We will try to guess some reasonable defaults */

	/* Have we at least got an adaptation white point ? */
	if (Wxyz[0] < 0.0 || Wxyz[1] < 0.0 || Wxyz[2] < 0.0)
		return 2;						/* No */

	/* Have we got the technology ? */
	if (tsig == icMaxEnumTechnology)
		return 2;						/* Hopeless */

	/* Guess from the technology */
	switch (tsig) {

		/* This is inputing either a representation of */
		/* a natural scene captured on another a print medium, or */
		/* are is assuming that the medium is the original. */
		/* We will assume that the print is the original. */
		case icSigReflectiveScanner:
			{
				if (La < 0.0)	/* No adapting luminance */
					La = 34.0;	/* Use a practical print evaluation number */
				if (Yb < 0.0)	/* No background relative luminance */
					Yb = 0.2;	/* Assume grey world */
				if (Lv < 0.0)	/* No device image luminance */
					Ev = vc_average;	/* Assume average viewing conditions */
				if (Yf < 0.0)	/* No flare figure */
					Yf = 0.01;	/* Assume 1% flare */
				if (Fxyz[0] < 0.0 || Fxyz[1] < 0.0 || Fxyz[2] < 0.0)	/* No flare color */
					Fxyz[0] = Wxyz[0], Fxyz[1] = Wxyz[1], Fxyz[2] = Wxyz[2];
				break;
			}

		/* This is inputing either a representation of */
		/* a natural scene captured on another a photo medium, or */
		/* are is assuming that the medium is the original. */
		/* We will assume a compromise media original, natural scene */
		case icSigFilmScanner:
			{
				if (La < 0.0)	/* No adapting luminance */
					La = 50.0;	/* Use bright indoors, dull outdoors */
				if (Yb < 0.0)	/* No background relative luminance */
					Yb = 0.2;	/* Assume grey world */
				if (Lv < 0.0)	/* No device image luminance */
					Ev = vc_average;	/* Assume average viewing conditions */
				if (Yf < 0.0)	/* No flare figure */
					Yf = 0.005;	/* Assume 0.5% flare */
				if (Fxyz[0] < 0.0 || Fxyz[1] < 0.0 || Fxyz[2] < 0.0)	/* No flare color */
					Fxyz[0] = Wxyz[0], Fxyz[1] = Wxyz[1], Fxyz[2] = Wxyz[2];
				break;
			}

		/* Direct scene to value devices. */
		case icSigDigitalCamera:
		case icSigVideoCamera:
			{
				if (La < 0.0)	/* No adapting luminance */
					La = 110.0;	/* Use very bright indoors, usual outdoors */
				if (Yb < 0.0)	/* No background relative luminance */
					Yb = 0.2;	/* Assume grey world */
				if (Lv < 0.0)	/* No device image luminance */
					Ev = vc_average;	/* Assume average viewing conditions */
				if (Yf < 0.0)	/* No flare figure */
					Yf = 0.0;	/* Assume 0% flare */
				if (Fxyz[0] < 0.0 || Fxyz[1] < 0.0 || Fxyz[2] < 0.0)	/* No flare color */
					Fxyz[0] = Wxyz[0], Fxyz[1] = Wxyz[1], Fxyz[2] = Wxyz[2];
				break;
			}

		/* Emmisive displays. */
		/* Assume a video monitor is in a darker environment than a CRT */
		case icSigVideoMonitor:
			{
				if (La < 0.0)	/* No adapting luminance */
					La = 4.0;	/* Darkened work environment */
				if (Yb < 0.0)	/* No background relative luminance */
					Yb = 0.2;	/* Assume grey world */
				if (Lv < 0.0)	/* No device image luminance */
					Ev = vc_dim;	/* Assume dim viewing conditions */
				if (Yf < 0.0)	/* No flare figure */
					Yf = 0.01;	/* Assume 1% flare */
				if (Fxyz[0] < 0.0 || Fxyz[1] < 0.0 || Fxyz[2] < 0.0)	/* No flare color */
					Fxyz[0] = Wxyz[0], Fxyz[1] = Wxyz[1], Fxyz[2] = Wxyz[2];
				break;
			}


		/* Assume a typical work environment */
		case icSigCRTDisplay:
		case icSigPMDisplay:
		case icSigAMDisplay:
			{
				if (La < 0.0)	/* No adapting luminance */
					La = 33.0;	/* Typical work environment */
				if (Yb < 0.0)	/* No background relative luminance */
					Yb = 0.2;	/* Assume grey world */
				if (Lv < 0.0)	/* No device image luminance */
					Ev = vc_average;	/* Assume average viewing conditions */
				if (Yf < 0.0)	/* No flare figure */
					Yf = 0.02;	/* Assume 2% flare */
				if (Fxyz[0] < 0.0 || Fxyz[1] < 0.0 || Fxyz[2] < 0.0)	/* No flare color */
					Fxyz[0] = Wxyz[0], Fxyz[1] = Wxyz[1], Fxyz[2] = Wxyz[2];
				break;
			}

		/* Photo CD has its own viewing definitions */
		/* (It represents original scene colors) */
		case icSigPhotoCD:
			{
				if (La < 0.0)	/* No adapting luminance */
					La = 320.0;	/* Bright outdoors */
				if (Yb < 0.0)	/* No background relative luminance */
					Yb = 0.2;	/* Assume grey world */
				if (Lv < 0.0)	/* No device image luminance */
					Ev = vc_average;	/* Assume average viewing conditions */
				if (Yf < 0.0)	/* No flare figure */
					Yf = 0.00;	/* Assume 0% flare */
				if (Fxyz[0] < 0.0 || Fxyz[1] < 0.0 || Fxyz[2] < 0.0)	/* No flare color */
					Fxyz[0] = Wxyz[0], Fxyz[1] = Wxyz[1], Fxyz[2] = Wxyz[2];
				break;
			}

		/* Projection devices, either direct, or */
		/* via another intermediate medium. */
		/* Assume darkened room, little background */
		case icSigProjectionTelevision:
			{
				if (La < 0.0)	/* No adapting luminance */
					La = 7.0;	/* Dark environment */
				if (Yb < 0.0)	/* No background relative luminance */
					Yb = 0.1;	/* Assume little background */
				if (Lv < 0.0)	/* No device image luminance */
					Ev = vc_dim;	/* Dim environment */
				if (Yf < 0.0)	/* No flare figure */
					Yf = 0.01;	/* Assume 1% flare */
				if (Fxyz[0] < 0.0 || Fxyz[1] < 0.0 || Fxyz[2] < 0.0)	/* No flare color */
					Fxyz[0] = Wxyz[0], Fxyz[1] = Wxyz[1], Fxyz[2] = Wxyz[2];
				break;
			}
		/* Assume very darkened room, no background */
		case icSigFilmWriter:
			{
				if (La < 0.0)	/* No adapting luminance */
					La = 7.0;	/* Dark environment */
				if (Yb < 0.0)	/* No background relative luminance */
					Yb = 0.0;	/* Assume no background */
				if (Lv < 0.0)	/* No device image luminance */
					Ev = vc_dark;	/* Dark environment */
				if (Yf < 0.0)	/* No flare figure */
					Yf = 0.01;	/* Assume 1% flare */
				if (Fxyz[0] < 0.0 || Fxyz[1] < 0.0 || Fxyz[2] < 0.0)	/* No flare color */
					Fxyz[0] = Wxyz[0], Fxyz[1] = Wxyz[1], Fxyz[2] = Wxyz[2];
				break;
			}

		/* Printed media devices. */
		/* Assume a normal print viewing environment */
		case icSigInkJetPrinter:
		case icSigThermalWaxPrinter:
		case icSigElectrophotographicPrinter:
		case icSigElectrostaticPrinter:
		case icSigDyeSublimationPrinter:
		case icSigPhotographicPaperPrinter:
		case icSigPhotoImageSetter:
		case icSigGravure:
		case icSigOffsetLithography:
		case icSigSilkscreen:
		case icSigFlexography:
			{
				if (La < 0.0)	/* No adapting luminance */
					La = 40.0;	/* Use a practical print evaluation number */
				if (Yb < 0.0)	/* No background relative luminance */
					Yb = 0.2;	/* Assume grey world */
				if (Lv < 0.0)	/* No device image luminance */
					Ev = vc_average;	/* Assume average viewing conditions */
				if (Yf < 0.0)	/* No flare figure */
					Yf = 0.01;	/* Assume 1% flare */
				if (Fxyz[0] < 0.0 || Fxyz[1] < 0.0 || Fxyz[2] < 0.0)	/* No flare color */
					Fxyz[0] = Wxyz[0], Fxyz[1] = Wxyz[1], Fxyz[2] = Wxyz[2];
				break;
			}

		default:
			{
			return 2;
			}
	}

	return 1;
}

/* Write our viewing conditions to the underlying ICC profile, */
/* using a private tag. */
void xicc_set_viewcond(
xicc *p,			/* Expanded profile we're working with */
icxViewCond *vc		/* Viewing parameters to return */
) {
	//icc *pp = p->pp;	/* Base ICC */

	// ~~1	Not implemented yet
}



/* Return an enumerated viewing condition */
/* Return enumeration if OK, -999 if there is no such enumeration. */
int xicc_enum_viewcond(
xicc *p,			/* Expanded profile to get white point (May be NULL if desc NZ) */
icxViewCond *vc,	/* Viewing parameters to return, May be NULL if desc is nz */
int no,				/* Enumeration to return, -1 for default, -2 for none */
char *as,			/* String alias to number, NULL if none */
int desc			/* NZ - Just return a description of this enumeration in vc */
) {

	if (desc == 0) {	/* We're setting the viewing condition */
		icc *pp;		/* Base ICC */
		icmXYZArray *whitePointTag;

		if (p == NULL)
			return -999;
	
		pp = p->pp;
		if ((whitePointTag = (icmXYZArray *)pp->read_tag(pp, icSigMediaWhitePointTag)) != NULL
         && whitePointTag->ttype == icSigXYZType && whitePointTag->size >= 1) {
			vc->Wxyz[0] = whitePointTag->data[0].X;
			vc->Wxyz[1] = whitePointTag->data[0].Y;
			vc->Wxyz[2] = whitePointTag->data[0].Z;

			/* Set a default flare color */
			vc->Fxyz[0] = vc->Wxyz[0];
			vc->Fxyz[1] = vc->Wxyz[1];
			vc->Fxyz[2] = vc->Wxyz[2];
		} else {
			sprintf(p->err,"Enum VC: Failed to read Media White point");
			p->errc = 2;
			return -999;
		}
	}

	/*

	Typical adapting field luminances and white luminance in reflective setup:

	E = illuminance in Lux
	Lv = White luminance assuming 100% reflectance
	La = Adapting field luminance in cd/m^2, assuming 20% reflectance from surround
    
	    E     La     Lv   Condition 
	    11     0.7    4   Twilight
	    32     2     10   Subdued indoor lighting
	    64     4     20   Less than typical office light; sometimes recommended for
	                      display-only workplaces (sRGB)
	   350    22    111   Typical Office (sRGB annex D)
	   500    32    160   Practical print evaluationa (ISO-3664 P2)
	  1000    64    318   Good Print evaluation (CIE 116-1995)
	  1000    64    318   Television Studio lighting
	  1000    64    318   Overcast Outdoors
	  2000   127    637   Critical print evaluation (ISO-3664 P1)
	 10000   637   3183   Typical outdoors, full daylight 
	 50000  3185  15915   Bright summers day 

	*/

	if (no == -1
	 || (as != NULL && stricmp(as,"d") == 0)) {

		no = -1;
		if (vc != NULL) {
			vc->desc = "  d: Default Viewing Condition";
			vc->Ev = vc_average;	/* Average viewing conditions */
			vc->La = 50.0;			/* Practical to Good lighting */
			vc->Yb = 0.2;			/* Grey world */
			vc->Yf = 0.01;			/* 1% flare */
		}
	}
	else if (no == 0
	 || (as != NULL && stricmp(as,"pp") == 0)) {

		no = 0;
		if (vc != NULL) {
			vc->desc = " pp - Practical Reflection Print";
			vc->Ev = vc_average;	/* Average viewing conditions */
			vc->La = 32.0;			/* Use a practical print evaluation number */
			vc->Yb = 0.2;			/* Grey world */
			vc->Yf = 0.01;			/* 1% flare */
		}
	}
	else if (no == 1
	 || (as != NULL && stricmp(as,"pe") == 0)) {

		no = 1;
		if (vc != NULL) {
			vc->desc = " pe - Print evaluation environment";
			vc->Ev = vc_average;	/* Average viewing conditions */
			vc->La = 64.0;			/* Good */
			vc->Yb = 0.2;			/* Grey world */
			vc->Yf = 0.01;			/* 1% flare */
		}
	}
	else if (no == 2
	 || (as != NULL && stricmp(as,"mt") == 0)) {

		no = 2;
		if (vc != NULL) {
			vc->desc = " mt - Monitor in typical work environment";
			vc->Ev = vc_average;	/* Average viewing conditions */
			vc->La = 22.0;			/* Typical work environment */
			vc->Yb = 0.2;			/* Grey world */
			vc->Yf = 0.02;			/* 2% flare */
		}
	}
	else if (no == 3
	 || (as != NULL && stricmp(as,"mb") == 0)) {

		no = 3;
		if (vc != NULL) {
			vc->desc = " mb - Bright monitor in bright work environment";
			vc->Ev = vc_average;	/* Average viewing conditions */
			vc->La = 42.0;			/* Bright work environment */
			vc->Yb = 0.2;			/* Grey world */
			vc->Yf = 0.02;			/* 2% flare */
		}
	}
	else if (no == 4
	 || (as != NULL && stricmp(as,"md") == 0)) {

		no = 4;
		if (vc != NULL) {
			vc->desc = " md - Monitor in darkened work environment";
			vc->Ev = vc_dim;		/* Dim viewing conditions */
			vc->La = 4.0;			/* Darkened work environment */
			vc->Yb = 0.2;			/* Grey world */
			vc->Yf = 0.01;			/* 1% flare */
		}
	}
	else if (no == 5
	 || (as != NULL && stricmp(as,"jm") == 0)) {

		no = 5;
		if (vc != NULL) {
			vc->desc = " jm - Projector in dim environment";
			vc->Ev = vc_dim;		/* Dim viewing conditions */
			vc->La = 10.0;			/* Adaptation is from display */
			vc->Yb = 0.2;			/* Grey world */
			vc->Yf = 0.01;			/* 1% flare */
		}
	}
	else if (no == 6
	 || (as != NULL && stricmp(as,"jd") == 0)) {

		no = 6;
		if (vc != NULL) {
			vc->desc = " jd - Projector in dark environment";
			vc->Ev = vc_dark;		/* Dark viewing conditions */
			vc->La = 10.0;			/* Adaptation is from display */
			vc->Yb = 0.2;			/* Grey world */
			vc->Yf = 0.01;			/* 1% flare ? */
		}
	}
	else if (no == 7
	 || (as != NULL && stricmp(as,"pcd") == 0)) {

		no = 7;
		if (vc != NULL) {
			vc->desc = "pcd - Photo CD - original scene outdoors";
			vc->Ev = vc_average;	/* Average viewing conditions */
			vc->La = 320.0;			/* Typical outdoors, 1600 cd/m^2 */
			vc->Yb = 0.2;			/* Grey world */
			vc->Yf = 0.00;			/* 0% flare */
		}
	}
	else if (no == 8
	 || (as != NULL && stricmp(as,"ob") == 0)) {

		no = 8;
		if (vc != NULL) {
			vc->desc = " ob - Original scene - Bright Outdoors";
			vc->Ev = vc_average;	/* Average viewing conditions */
			vc->La = 2000.0;		/* Bright Outdoors */
			vc->Yb = 0.2;			/* Grey world */
			vc->Yf = 0.00;			/* 0% flare */
		}
	}
	else if (no == 9
	 || (as != NULL && stricmp(as,"cx") == 0)) {

		no = 9;
		if (vc != NULL) {
			vc->desc = " cx - Cut Sheet Transparencies on a viewing box";
			vc->Ev = vc_cut_sheet;	/* Cut sheet viewing conditions */
			vc->La = 53.0;			/* Dim, adapted to slide ? */
			vc->Yb = 0.2;			/* Grey world */
			vc->Yf = 0.01;			/* 1% flare ? */
		}
	}
	else {
		if (p != NULL) {
			sprintf(p->err,"Enum VC: Unrecognised enumeration %d",no);
			p->errc = 1;
		}
		return -999;
	}

	return no;
}

/* Debug: dump a Viewing Condition to standard out */
void xicc_dump_viewcond(
icxViewCond *vc
) {
	printf("Viewing Condition:\n");
	if (vc->Ev == vc_dark)
		printf("  Surround to Image: Dark\n");
	else if (vc->Ev == vc_dim)
		printf("  Surround to Image: Dim\n");
	else if (vc->Ev == vc_average)
		printf("  Surround to Image: Average\n");
	else if (vc->Ev == vc_cut_sheet)
		printf("  Transparency on Light box\n");
	printf("  Adapted white = %f %f %f\n",vc->Wxyz[0], vc->Wxyz[1], vc->Wxyz[2]);
	printf("  Adapted luminance = %f cd/m^2\n",vc->La);
	printf("  Background to image ratio = %f\n",vc->Yb);
	if (vc->Ev == vc_none)
		printf("  Image luminance = %f cd/m^2\n",vc->Lv);
	printf("  Flare to image ratio = %f\n",vc->Yf);
	printf("  Flare color = %f %f %f\n",vc->Fxyz[0], vc->Fxyz[1], vc->Fxyz[2]);
}


/* Debug: dump an Inking setup to standard out */
void xicc_dump_inking(icxInk *ik) {
	printf("Inking settings:\n");
	if (ik->tlimit < 0.0)
		printf("No total limit\n");
	else
		printf("Total limit = %f%%\n",ik->tlimit * 100.0);
	if (ik->klimit < 0.0)
		printf("No black limit\n");
	else
		printf("Black limit = %f%%\n",ik->klimit * 100.0);

	if (ik->k_rule == icxKvalue) {
		printf("Inking rule is a fixed K target\n");
	} if (ik->k_rule == icxKlocus) {
		printf("Inking rule is a fixed locus target\n");
	} if (ik->k_rule == icxKluma5 || ik->k_rule == icxKluma5k) {
		if (ik->k_rule == icxKluma5)
			printf("Inking rule is a 5 parameter locus function of L\n");
		else
			printf("Inking rule is a 5 parameter K function of L\n");
		printf("Ksmth = %f\n",ik->c.Ksmth);
		printf("Kstle = %f\n",ik->c.Kstle);
		printf("Kstpo = %f\n",ik->c.Kstpo);
		printf("Kenpo = %f\n",ik->c.Kenpo);
		printf("Kenle = %f\n",ik->c.Kenle);
		printf("Kshap = %f\n",ik->c.Kshap);
	} if (ik->k_rule == icxKl5l || ik->k_rule == icxKl5lk) {
		if (ik->k_rule == icxKl5l)
			printf("Inking rule is a 2x5 parameter locus function of L and K aux\n");
		else
			printf("Inking rule is a 2x5 parameter K function of L and K aux\n");
		printf("Min Ksmth = %f\n",ik->c.Ksmth);
		printf("Min Kstle = %f\n",ik->c.Kstle);
		printf("Min Kstpo = %f\n",ik->c.Kstpo);
		printf("Min Kenpo = %f\n",ik->c.Kenpo);
		printf("Min Kenle = %f\n",ik->c.Kenle);
		printf("Min Kshap = %f\n",ik->c.Kshap);
		printf("Max Ksmth = %f\n",ik->x.Ksmth);
		printf("Max Kstle = %f\n",ik->x.Kstle);
		printf("Max Kstpo = %f\n",ik->x.Kstpo);
		printf("Max Kenpo = %f\n",ik->x.Kenpo);
		printf("Max Kenle = %f\n",ik->x.Kenle);
		printf("Max Kshap = %f\n",ik->x.Kshap);
	} 
}

/* ------------------------------------------------------ */
/* Gamut Mapping Intent stuff */

/* Return an enumerated gamut mapping intent */
/* Return enumeration if OK, icxIllegalGMIntent if there is no such enumeration. */
int xicc_enum_gmapintent(
icxGMappingIntent *gmi,	/* Gamut Mapping parameters to return */
int no,					/* Enumeration selected, icxNoGMIntent for none */
char *as				/* Alias string selector, NULL for none */
) {
#ifdef USE_CAM
	int colccas = 0x2;	/* Use cas clipping for colorimetric style intents */
	int perccas = 0x1;	/* Use cas for perceptual style intents */
#else
	int colccas = 0x0;	/* Use Lab for colorimetric style intents */
	int perccas = 0x0;	/* Use Lab for perceptual style intents */
	fprintf(stderr,"!!!!!! Warning, USE_CAM is off in xicc.c !!!!!!\n");
#endif

	/* Assert default if no guidance given */
	if (no == icxNoGMIntent && as == NULL)
		no = icxDefaultGMIntent;

	if (no == 0
	 || no == icxAbsoluteGMIntent
	 || (as != NULL && stricmp(as,"a") == 0)) {
		/* Map Absolute Jab to Jab and clip out of gamut */
		no = 0;
		gmi->desc = " a - Absolute Colorimetric (in Jab) [ICC Absolute Colorimetric]";
		gmi->icci = icAbsoluteColorimetric;
		gmi->usecas  = colccas;		/* Use absolute appearance space */
		gmi->usemap  = 0;			/* Don't use gamut mapping */
		gmi->greymf  = 0.0;
		gmi->glumwcpf = 0.0;
		gmi->glumwexf = 0.0;
		gmi->glumbcpf = 0.0;
		gmi->glumbexf = 0.0;
		gmi->glumknf = 0.0;
		gmi->gamcpf  = 0.0;
		gmi->gamexf  = 0.0;
		gmi->gamknf  = 0.0;
		gmi->gampwf  = 0.0;
		gmi->gamswf  = 0.0;
		gmi->satenh  = 0.0;			/* No saturation enhancement */
	}
	else if (no == 1
	 || (as != NULL && stricmp(as,"aw") == 0)) {

		/* I'm not sure how often this intent is useful. It's less likely than */
		/* I though that a printer white point won't fit within the gamut */
		/* of a display profile, since the display white always has Y = 1.0, */
		/* and no paper has better than about 95% reflectance. */

		/* Map Absolute Jab to Jab and scale source to avoid clipping the white point */
		no = 1;
		gmi->desc = "aw - Absolute Colorimetric (in Jab) with scaling to fit white point";
		gmi->icci = icAbsoluteColorimetric;
		gmi->usecas = 0x100 | colccas;	/* Absolute Appearance space with scaling */
									/* to avoid clipping the source white point */
		gmi->usemap  = 0;			/* Don't use gamut mapping */
		gmi->greymf  = 0.0;
		gmi->glumwcpf = 0.0;
		gmi->glumwexf = 0.0;
		gmi->glumbcpf = 0.0;
		gmi->glumbexf = 0.0;
		gmi->glumknf = 0.0;
		gmi->gamcpf  = 0.0;
		gmi->gamexf  = 0.0;
		gmi->gamknf  = 0.0;
		gmi->gampwf  = 0.0;
		gmi->gamswf  = 0.0;
		gmi->satenh  = 0.0;			/* No saturation enhancement */
	}
	else if (no == 2
	 || (as != NULL && stricmp(as,"aa") == 0)) {

		/* Map Jab to Jab and clip out of gamut */
		no = 2;
		gmi->desc = "aa - Absolute Appearance";
		gmi->icci = icRelativeColorimetric;
		gmi->usecas  = perccas;		/* Appearance space */
		gmi->usemap  = 0;			/* Don't use gamut mapping */
		gmi->greymf  = 0.0;
		gmi->glumwcpf = 0.0;
		gmi->glumwexf = 0.0;
		gmi->glumbcpf = 0.0;
		gmi->glumbexf = 0.0;
		gmi->glumknf = 0.0;
		gmi->gamcpf  = 0.0;
		gmi->gamexf  = 0.0;
		gmi->gamknf  = 0.0;
		gmi->gampwf  = 0.0;
		gmi->gamswf  = 0.0;
		gmi->satenh  = 0.0;			/* No saturation enhancement */
	}
	else if (no == 3
	 || no == icxRelativeGMIntent
	 || (as != NULL && stricmp(as,"r") == 0)) {

		/* Map Jab to Jab and clip out of gamut */
		/* Linear transform of the neutral axis to align white */
		no = 3;
		gmi->desc = " r - White Point Matched Appearance [ICC Relative Colorimetric]";
		gmi->icci = icRelativeColorimetric;
		gmi->usecas  = perccas;		/* Appearance space */
		gmi->usemap  = 1;			/* Use gamut mapping */
		gmi->greymf  = 1.0;			/* And linearly map white point */
		gmi->glumwcpf = 1.0;
		gmi->glumwexf = 1.0;
		gmi->glumbcpf = 0.0;
		gmi->glumbexf = 0.0;
		gmi->glumknf = 0.0;
		gmi->gamcpf  = 0.0;
		gmi->gamexf  = 0.0;
		gmi->gamknf  = 0.0;
		gmi->gampwf  = 0.0;
		gmi->gamswf  = 0.0;
		gmi->satenh  = 0.0;			/* No saturation enhancement */
	}
	else if (no == 4
	 || (as != NULL && stricmp(as,"la") == 0)) {

		/* Map Jab to Jab, sigma map white/black points, and clip out of gamut */
		no = 4;
		gmi->desc = "la - Luminance axis matched Appearance";
		gmi->icci = icRelativeColorimetric;
		gmi->usecas  = perccas;		/* Appearance space */
		gmi->usemap  = 1;			/* Use gamut mapping */
		gmi->greymf  = 1.0;			/* Fully align grey axis */
		gmi->glumwcpf = 1.0;		/* Fully compress grey axis */
		gmi->glumwexf = 1.0;		/* Fully expand grey axis at black end */
		gmi->glumbcpf = 1.0;		/* Fully compress grey axis at black end */
		gmi->glumbexf = 1.0;		/* Fully expand grey axis */
		gmi->glumknf = 1.0;			/* Distort at white/black ends only */
		gmi->gamcpf  = 0.0;			/* No gamut compression */
		gmi->gamexf  = 0.0;			/* No gamut expansion */
		gmi->gamknf  = 0.0;			/* No knee in gamut compress/expand */
		gmi->gampwf  = 0.0;			/* No Perceptual surface weighting factor */
		gmi->gamswf  = 0.0;			/* No Saturation surface weighting factor */
		gmi->satenh  = 0.0;			/* No saturation enhancement */
	}
	else if (no == 5
	 || no == icxDefaultGMIntent
	 || no == icxPerceptualGMIntent
	 || (as != NULL && stricmp(as,"p") == 0)) {

		/* Map Jab to Jab, sigma map white/black, compress out of gamut */
		/* NOTE would like to be using sigma knee on gamut in these two, */
		/* but the current gamut mapping isn't linear enough to need extra knee :-) */
		no = 5;
		gmi->desc = " p - Perceptual (Preferred) [ICC Perceptual]";
		gmi->icci = icPerceptual;
		gmi->usecas  = perccas;		/* Appearance space */
		gmi->usemap  = 1;			/* Use gamut mapping */
		gmi->greymf  = 1.0;			/* Fully align grey axis */
		gmi->glumwcpf = 1.0;		/* Fully compress grey axis at white end */
		gmi->glumwexf = 1.0;		/* Fully expand grey axis at white end */
		gmi->glumbcpf = 1.0;		/* Fully compress grey axis at black end */
		gmi->glumbexf = 1.0;		/* Fully expand grey axis at black end */
		gmi->glumknf = 1.0;			/* Sigma knee in grey compress/expand */
		gmi->gamcpf  = 1.0;			/* Full gamut compression */
		gmi->gamexf  = 0.0;			/* No gamut expansion */
		gmi->gamknf  = 0.5;			/* Sigma knee in gamut compress/expand */
		gmi->gampwf  = 1.0;			/* Full Perceptual surface weighting factor */
		gmi->gamswf  = 0.0;			/* No Saturation surface weighting factor */
		gmi->satenh  = 0.0;			/* No saturation enhancement */
	}
	else if (no == 6
	 || (as != NULL && stricmp(as,"ms") == 0)) {

		/* Map Jab to Jab, sigma map whole gamut */
		no = 6;
		gmi->desc = "ms - Saturation";
		gmi->icci = icSaturation;
		gmi->usecas  = perccas;		/* Appearance space */
		gmi->usemap  = 1;			/* Use gamut mapping */
		gmi->greymf  = 1.0;			/* Fully align grey axis */
		gmi->glumwcpf = 1.0;		/* Fully compress grey axis at white end */
		gmi->glumwexf = 1.0;		/* Fully expand grey axis at white end */
		gmi->glumbcpf = 1.0;		/* Fully compress grey axis at black end */
		gmi->glumbexf = 1.0;		/* Fully expand grey axis at black end */
		gmi->glumknf = 1.0;			/* Sigma knee in grey compress/expand */
		gmi->gamcpf  = 0.95;		/* Almost full gamut compression */
		gmi->gamexf  = 1.2;			/* Excessive expansion to compensate for rspl effect */
		gmi->gamknf  = 1.0;			/* Sigma knee in gamut compress/expand */
		gmi->gampwf  = 0.0;			/* No Perceptual surface weighting factor */
		gmi->gamswf  = 1.0;			/* Full Saturation surface weighting factor */
		gmi->satenh  = 0.0;			/* No saturation enhancement */
	}
	else if (no == 7
	 || no == icxSaturationGMIntent
	 || (as != NULL && stricmp(as,"s") == 0)) {

		/* Map Jab to Jab, sigma map whole gamut */
		no = 7;
		gmi->desc = " s - Enhanced Saturation [ICC Saturation]";
		gmi->icci = icSaturation;
		gmi->usecas  = perccas;		/* Appearance space */
		gmi->usemap  = 1;			/* Use gamut mapping */
		gmi->greymf  = 1.0;			/* Fully align grey axis */
		gmi->glumwcpf = 1.0;		/* Fully compress grey axis at white end */
		gmi->glumwexf = 1.0;		/* Fully expand grey axis at white end */
		gmi->glumbcpf = 1.0;		/* Fully compress grey axis at black end */
		gmi->glumbexf = 1.0;		/* Fully expand grey axis at black end */
		gmi->glumknf = 1.0;			/* Sigma knee in grey compress/expand */
		gmi->gamcpf  = 0.90;		/* Almost full gamut compression */
		gmi->gamexf  = 1.3;			/* Excessive expansion to compensate for rspl effect */
		gmi->gamknf  = 1.0;			/* Sigma knee in gamut compress/expand */
		gmi->gampwf  = 0.0;			/* No Perceptual surface weighting factor */
		gmi->gamswf  = 1.0;			/* Full Saturation surface weighting factor */
		gmi->satenh  = 0.8;			/* Medium saturation enhancement */
	}
	else if (no == 8
	 || (as != NULL && stricmp(as,"al") == 0)) {

		/* Map Absolute Lab to Lab and clip out of gamut */
		no = 8;
		gmi->desc = "al - Absolute Colorimetric (Lab)";
		gmi->icci = icAbsoluteColorimetric;
		gmi->usecas  = 0x0;			/* Don't use appearance space, use L*a*b* */
		gmi->usemap  = 0;			/* Don't use gamut mapping */
		gmi->greymf  = 0.0;
		gmi->glumwcpf = 0.0;
		gmi->glumwexf = 0.0;
		gmi->glumbcpf = 0.0;
		gmi->glumbexf = 0.0;
		gmi->glumknf = 0.0;
		gmi->gamcpf  = 0.0;
		gmi->gamexf  = 0.0;
		gmi->gamknf  = 0.0;
		gmi->gampwf  = 0.0;
		gmi->gamswf  = 0.0;
		gmi->satenh  = 0.0;			/* No saturation enhancement */
	}
	else if (no == 9
	 || (as != NULL && stricmp(as,"rl") == 0)) {

		/* Map Lab to Lab and clip out of gamut */
		/* Linear transform of the neutral axis to align white */
		no = 3;
		gmi->desc = "rl - White Point Matched Appearance (Lab)";
		gmi->icci = icRelativeColorimetric;
		gmi->usecas  = 0x0;			/* Don't use appearance space, use L*a*b* */
		gmi->usemap  = 1;			/* Use gamut mapping */
		gmi->greymf  = 1.0;			/* And linearly map white point */
		gmi->glumwcpf = 1.0;
		gmi->glumwexf = 1.0;
		gmi->glumbcpf = 0.0;
		gmi->glumbexf = 0.0;
		gmi->glumknf = 0.0;
		gmi->gamcpf  = 0.0;
		gmi->gamexf  = 0.0;
		gmi->gamknf  = 0.0;
		gmi->gampwf  = 0.0;
		gmi->gamswf  = 0.0;
		gmi->satenh  = 0.0;			/* No saturation enhancement */
	}
	else {	/* icxIllegalGMIntent */
		return icxIllegalGMIntent;
	}

	return no;
}


/* Debug: dump a Gamut Mapping specification */
void xicc_dump_gmi(
icxGMappingIntent *gmi	/* Gamut Mapping parameters to return */
) {
	printf(" Gamut Mapping Specification:\n");
	if (gmi->desc != NULL)
		printf("  Description = '%s'\n",gmi->desc);
	printf("  Closest ICC intent = '%s'\n",icm2str(icmRenderingIntent,gmi->icci));

	if ((gmi->usecas & 0xff) == 0)
		printf("  Not using Color Apperance Space\n");
	else if ((gmi->usecas & 0xff) == 1)
		printf("  Using Color Apperance Space\n");
	else if ((gmi->usecas & 0xff) == 2)
		printf("  Using Absolute Color Apperance Space\n");

	if ((gmi->usecas & 0x100) != 0)
		printf("  Scaling source to avoid white point clipping\n");

	if (gmi->usemap == 0)
		printf("  Not using Mapping\n");
	else {
		printf("  Using Mapping with parameters:\n");
		printf("  Grey axis alignment   factor %f\n", gmi->greymf);
		printf("  Grey axis white compression factor %f\n", gmi->glumwcpf);
		printf("  Grey axis white expansion   factor %f\n", gmi->glumwexf);
		printf("  Grey axis black compression factor %f\n", gmi->glumbcpf);
		printf("  Grey axis black expansion   factor %f\n", gmi->glumbexf);
		printf("  Grey axis knee        factor %f\n", gmi->glumknf);
		printf("  Gamut compression factor %f\n", gmi->gamcpf);
		printf("  Gamut expansion   factor %f\n", gmi->gamexf);
		printf("  Gamut knee        factor %f\n", gmi->gamknf);
		printf("  Gamut Perceptual mapping weighting factor %f\n", gmi->gampwf);
		printf("  Gamut Saturation mapping weighting factor %f\n", gmi->gamswf);
		printf("  Saturation enhancement factor %f\n", gmi->satenh);
	}
}

/* ------------------------------------------------------ */

/* Utility function - given an open icc profile, */
/* estmate the total ink limit and black ink limit. */
/* Note that this is rather rough, because ICC profiles */
/* don't have a tag for this information, and ICC profiles */
/* don't have any straightforward way of identifying particular */
/* color channels for > 4 color. */
/* If there are no limits, or they are not discoverable or */
/* applicable, return values of -1.0 */

void icxGetLimits(icc *p, double *tlimit, double *klimit) {
	int nch;
	double max[MAX_CHAN];
	double total;

	total = p->get_tac(p, max);

	if (total < 0.0) {	/* Not valid */
		if (tlimit != NULL)
			*tlimit = -1.0;
		if (klimit != NULL)
			*klimit = -1.0;
		return;
	}

	nch = icmCSSig2nchan(p->header->colorSpace);

	/* No effective limit */
	if (tlimit != NULL) {
		if (total >= (double)nch) {
			*tlimit = -1.0;
		} else {
			*tlimit = total;
		}
	}

	if (klimit != NULL) {
		double black = 1.0;
		switch (p->header->colorSpace) {
	   	 	case icSigCmykData:
				black = max[3];
				break;
			/* Hmm. It would be nice to know the black channel of other spaces */
			/* This could be fixed by using other tags + heuristics on the profile */
			default:
				break;
		}

		if (black >= 1.0) {
			*klimit = -1.0;
		} else {
			*klimit = black;
		}
	}
}

/* Using the above function, set default total and black ink values */
void icxDefaultLimits(icc *p, double *tlout, double tlin, double *klout, double klin) {
	if (tlin < 0.0 || klin < 0.0) { 
		double tl, kl;

		icxGetLimits(p, &tl, &kl);

		if (tlin < 0.0)
			tlin = tl;

		if (klin < 0.0)
			klin = kl;
	}

	if (tlout != NULL)
		*tlout = tlin;

	if (klout != NULL)
		*klout = klin;
}

/* ------------------------------------------------------ */
/* Conversion and deltaE formular that include partial */
/* derivatives, for use within fit parameter optimisations. */

/* CIE XYZ to perceptual Lab with partial derivatives. */
void icxdXYZ2Lab(icmXYZNumber *w, double *out, double dout[3][3], double *in) {
	double wp[3], tin[3], dtin[3];
	int i;

	wp[0] = w->X, wp[1] = w->Y, wp[2] = w->Z;

	for (i = 0; i < 3; i++) {
		tin[i] = in[i]/wp[i];
		dtin[i] = 1.0/wp[i];

		if (tin[i] > 0.008856451586) {
			dtin[i] *= pow(tin[i], -2.0/3.0) / 3.0;
			tin[i] = pow(tin[i],1.0/3.0);
		} else {
			dtin[i] *= 7.787036979;
			tin[i] = 7.787036979 * tin[i] + 16.0/116.0;
		}
	}

	out[0] = 116.0 * tin[1] - 16.0;
	dout[0][0] = 0.0;
	dout[0][1] = 116.0 * dtin[1];
	dout[0][2] = 0.0;

	out[1] = 500.0 * (tin[0] - tin[1]);
	dout[1][0] = 500.0 * dtin[0];
	dout[1][1] = 500.0 * -dtin[1];
	dout[1][2] = 0.0;

	out[2] = 200.0 * (tin[1] - tin[2]);
	dout[2][0] = 0.0;
	dout[2][1] = 200.0 * dtin[1];
	dout[2][2] = 200.0 * -dtin[2];
}


/* Return the normal Delta E squared, given two Lab values, */
/* including partial derivatives. */
double icxdLabDEsq(double dout[2][3], double *Lab0, double *Lab1) {
	double rv = 0.0, tt;

	tt = Lab0[0] - Lab1[0];
	dout[0][0] =  2.0 * tt;
	dout[1][0] = -2.0 * tt;
	rv += tt * tt;
	tt = Lab0[1] - Lab1[1];
	dout[0][1] =  2.0 * tt;
	dout[1][1] = -2.0 * tt;
	rv += tt * tt;
	tt = Lab0[2] - Lab1[2];
	dout[0][2] =  2.0 * tt;
	dout[1][2] = -2.0 * tt;
	rv += tt * tt;
	return rv;
}

/* Return the CIE94 Delta E color difference measure, squared */
/* including partial derivatives. */
double icxdCIE94sq(double dout[2][3], double Lab0[3], double Lab1[3]) {
	double desq, _desq[2][3];
	double dlsq;
	double dcsq, _dcsq[2][2];	/* == [x][1,2] */
	double c12,   _c12[2][2];	/* == [x][1,2] */
	double dhsq, _dhsq[2][2];	/* == [x][1,2] */
	double rv;

	{
		double dl, da, db;

		dl = Lab0[0] - Lab1[0];
		dlsq = dl * dl;		/* dl squared */
		da = Lab0[1] - Lab1[1];
		db = Lab0[2] - Lab1[2];


		/* Compute normal Lab delta E squared */
		desq = dlsq + da * da + db * db;
		_desq[0][0] =  2.0 * dl;
		_desq[1][0] = -2.0 * dl;
		_desq[0][1] =  2.0 * da;
		_desq[1][1] = -2.0 * da;
		_desq[0][2] =  2.0 * db;
		_desq[1][2] = -2.0 * db;
	}

	{
		double c1, c2, dc, tt;

		/* Compute chromanance for the two colors */
		c1 = sqrt(Lab0[1] * Lab0[1] + Lab0[2] * Lab0[2]);
		c2 = sqrt(Lab1[1] * Lab1[1] + Lab1[2] * Lab1[2]);
		c12 = sqrt(c1 * c2);	/* Symetric chromanance */

		tt = 0.5 * (pow(c2, 0.5) + 1e-12)/(pow(c1, 1.5) + 1e-12);
		_c12[0][0] = Lab0[1] * tt;
		_c12[0][1] = Lab0[2] * tt;
		tt = 0.5 * (pow(c1, 0.5) + 1e-12)/(pow(c2, 1.5) + 1e-12);
		_c12[1][0] = Lab1[1] * tt;
		_c12[1][1] = Lab1[2] * tt;

		/* delta chromanance squared */
		dc = c2 - c1;
		dcsq = dc * dc;
		if (c1 < 1e-12 || c2 < 1e-12) {
			c1 += 1e-12;
			c2 += 1e-12;
		}
		_dcsq[0][0] = -2.0 * Lab0[1] * (c2 - c1)/c1;
		_dcsq[0][1] = -2.0 * Lab0[2] * (c2 - c1)/c1;
		_dcsq[1][0] =  2.0 * Lab1[1] * (c2 - c1)/c2;
		_dcsq[1][1] =  2.0 * Lab1[2] * (c2 - c1)/c2;
	}

	/* Compute delta hue squared */
	dhsq = desq - dlsq - dcsq;
	if (dhsq >= 0.0) {
		_dhsq[0][0] = _desq[0][1] - _dcsq[0][0];
		_dhsq[0][1] = _desq[0][2] - _dcsq[0][1];
		_dhsq[1][0] = _desq[1][1] - _dcsq[1][0];
		_dhsq[1][1] = _desq[1][2] - _dcsq[1][1];
	} else {
		dhsq = 0.0;
		_dhsq[0][0] = 0.0;
		_dhsq[0][1] = 0.0;
		_dhsq[1][0] = 0.0;
		_dhsq[1][1] = 0.0;
	}

	{
		double sc, scsq, scf;
		double sh, shsq, shf;

		/* Weighting factors for delta chromanance & delta hue */
		sc = 1.0 + 0.048 * c12;
		scsq = sc * sc;
		
		sh = 1.0 + 0.014 * c12;
		shsq = sh * sh;
	
		rv = dlsq + dcsq/scsq + dhsq/shsq;

		scf = 0.048 * -2.0 * dcsq/(scsq * sc);
		shf = 0.014 * -2.0 * dhsq/(shsq * sh);
		dout[0][0] = _desq[0][0];
		dout[0][1] = _dcsq[0][0]/scsq + _c12[0][0] * scf
		           + _dhsq[0][0]/shsq + _c12[0][0] * shf;
		dout[0][2] = _dcsq[0][1]/scsq + _c12[0][1] * scf
		           + _dhsq[0][1]/shsq + _c12[0][1] * shf;
		dout[1][0] = _desq[1][0];
		dout[1][1] = _dcsq[1][0]/scsq + _c12[1][0] * scf
		           + _dhsq[1][0]/shsq + _c12[1][0] * shf;
		dout[1][2] = _dcsq[1][1]/scsq + _c12[1][1] * scf
		           + _dhsq[1][1]/shsq + _c12[1][1] * shf;
		return rv;
	}
}

/* ------------------------------------------------------ */
/* Parameterized transfer/dot gain function. */
/* Used for device modelling. Including partial */
/* derivative for input and parameters. */

/* NOTE that clamping the input values seems to cause */
/* conjgrad() problems. */

/* Transfer function */
double icxTransFunc(
double *v,			/* Pointer to first parameter */
int    luord,		/* Number of parameters */
double vv			/* Source of value */
) {
	double g;
	int ord;

	/* Process all the shaper orders from low to high. */
	/* [These shapers were inspired by a Graphics Gem idea */
	/* (Gems IV, VI.3, "Fast Alternatives to Perlin's Bias and */
	/*  Gain Functions, pp 401). */
	/*  They have the nice properties that they are smooth, and */
	/*  are monotonic. The control parameter has been */
	/*  altered to have a range from -oo to +oo rather than 0.0 to 1.0 */
	/*  so that the search space is less non-linear. */
	for (ord = 0; ord < luord; ord++) {
		int nsec;			/* Number of sections */
		double sec;			/* Section */

		g = v[ord];			/* Parameter */

		nsec = ord + 1;		/* Increase sections for each order */

		vv *= (double)nsec;

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
		vv /= (double)nsec;
	}

	return vv;
}

/* Inverse transfer function */
double icxInvTransFunc(
double *v,			/* Pointer to first parameter */
int    luord,		/* Number of parameters */
double vv			/* Source of value */
) {
	double g;
	int ord;

	/* Process the shaper orders in reverse from high to low. */
	/* [These shapers were inspired by a Graphics Gem idea */
	/* (Gems IV, VI.3, "Fast Alternatives to Perlin's Bias and */
	/*  Gain Functions, pp 401). */
	/*  They have the nice properties that they are smooth, and */
	/*  are monotonic. The control parameter has been */
	/*  altered to have a range from -oo to +oo rather than 0.0 to 1.0 */
	/*  so that the search space is less non-linear. */
	for (ord = luord-1; ord >= 0; ord--) {
		int nsec;			/* Number of sections */
		double sec;			/* Section */

		g = -v[ord];		/* Inverse parameter */

		nsec = ord + 1;		/* Increase sections for each order */

		vv *= (double)nsec;

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
		vv /= (double)nsec;
	}

	return vv;
}

/* Transfer function with offset and scale */
double icxSTransFunc(
double *v,			/* Pointer to first parameter */
int    luord,		/* Number of parameters */
double vv,			/* Source of value */
double min,			/* Scale values */
double max
) {
	max -= min;

	vv = (vv - min)/max;
	vv = icxTransFunc(v, luord, vv);
	vv = (vv * max) + min;
	return vv;
}

/* Inverse Transfer function with offset and scale */
double icxInvSTransFunc(
double *v,			/* Pointer to first parameter */
int    luord,		/* Number of parameters */
double vv,			/* Source of value */
double min,			/* Scale values */
double max
) {
	max -= min;

	vv = (vv - min)/max;
	vv = icxInvTransFunc(v, luord, vv);
	vv = (vv * max) + min;
	return vv;
}

/* Transfer function with partial derivative */
/* with respect to the parameters. */
double icxdpTransFunc(
double *v,			/* Pointer to first parameter */
double *dv,			/* Return derivative wrt each parameter */
int    luord,		/* Number of parameters */
double vv			/* Source of value */
) {
	double g;
	int i, ord;

	/* Process all the shaper orders from high to low. */
	for (ord = 0; ord < luord; ord++) {
		double dsv;		/* del for del in g */
		double ddv;		/* del for del in vv */
		int nsec;		/* Number of sections */
		double sec;		/* Section */

		g = v[ord];			/* Parameter */

		nsec = ord + 1;		/* Increase sections for each order */

		vv *= (double)nsec;

		sec = floor(vv);
		if (((int)sec) & 1) {
			g = -g;				/* Alternate action in each section */
		}
		vv -= sec;
		if (g >= 0.0) {
			double tt = g - g * vv + 1.0;
			dsv = (vv * vv - vv)/(tt * tt);
			ddv = (g + 1.0)/(tt * tt);
			vv = vv/tt;
		} else {
			double tt = 1.0 - g * vv;
			dsv = (vv * vv - vv)/(tt * tt);
			ddv = (1.0 - g)/(tt * tt);
			vv = (vv - g * vv)/tt;
		}

		vv += sec;
		vv /= (double)nsec;
		dsv /= (double)nsec;
		if (((int)sec) & 1)
			dsv = -dsv;

		dv[ord] = dsv;
		for (i = ord - 1; i >= 0; i--)
			dv[i] *= ddv;
	}

	return vv;
}

/* Transfer function with offset and scale, and */
/* partial derivative with respect to the parameters. */
double icxdpSTransFunc(
double *v,			/* Pointer to first parameter */
double *dv,			/* Return derivative wrt each parameter */
int    luord,		/* Number of parameters */
double vv,			/* Source of value */
double min,			/* Scale values */
double max
) {
	int i;
	max -= min;

	vv = (vv - min)/max;
	vv = icxdpTransFunc(v, dv, luord, vv);
	vv = (vv * max) + min;
	for (i = 0; i < luord; i++)
		dv[i] *= max;
	return vv;
}


/* Transfer function with partial derivative */
/* with respect to the input value. */
double icxdiTransFunc(
double *v,			/* Pointer to first parameter */
double *pdin,		/* Return derivative wrt source value */
int    luord,		/* Number of parameters */
double vv			/* Source of value */
) {
	double g, din;
	int ord;

#ifdef NEVER
	if (vv < 0.0 || vv > 1.0) {
		if (vv < 0.0)
			vv = 0.0;
	 	else
			vv = 1.0;

		*pdin = 0.0;
		return vv;
	}
#endif
	din = 1.0;

	/* Process all the shaper orders from high to low. */
	for (ord = 0; ord < luord; ord++) {
		double ddv;		/* del for del in vv */
		int nsec;		/* Number of sections */
		double sec;		/* Section */

		g = v[ord];			/* Parameter */

		nsec = ord + 1;		/* Increase sections for each order */

		vv *= (double)nsec;

		sec = floor(vv);
		if (((int)sec) & 1) {
			g = -g;				/* Alternate action in each section */
		}
		vv -= sec;
		if (g >= 0.0) {
			double tt = g - g * vv + 1.0;
			ddv = (g + 1.0)/(tt * tt);
			vv = vv/tt;
		} else {
			double tt = 1.0 - g * vv;
			ddv = (1.0 - g)/(tt * tt);
			vv = (vv - g * vv)/tt;
		}

		vv += sec;
		vv /= (double)nsec;
		din *= ddv;
	}

	*pdin = din;
	return vv;
}

/* Transfer function with offset and scale, and */
/* partial derivative with respect to the input value. */
double icxdiSTransFunc(
double *v,			/* Pointer to first parameter */
double *pdv,		/* Return derivative wrt source value */
int    luord,		/* Number of parameters */
double vv,			/* Source of value */
double min,			/* Scale values */
double max
) {
	max -= min;

	vv = (vv - min)/max;
	vv = icxdiTransFunc(v, pdv, luord, vv);
	vv = (vv * max) + min;
	return vv;
}


/* Transfer function with partial derivative */
/* with respect to the parameters and the input value. */
double icxdpdiTransFunc(
double *v,			/* Pointer to first parameter */
double *dv,			/* Return derivative wrt each parameter */
double *pdin,		/* Return derivative wrt source value */
int    luord,		/* Number of parameters */
double vv			/* Source of value */
) {
	double g, din;
	int i, ord;

#ifdef NEVER
	if (vv < 0.0 || vv > 1.0) {
		if (vv < 0.0)
			vv = 0.0;
	 	else
			vv = 1.0;

		for (ord = 0; ord < luord; ord++)
			dv[ord] = 0.0;
		*pdin = 0.0;
		return vv;
	}
#endif
	din = 1.0;

	/* Process all the shaper orders from high to low. */
	for (ord = 0; ord < luord; ord++) {
		double dsv;		/* del for del in g */
		double ddv;		/* del for del in vv */
		int nsec;		/* Number of sections */
		double sec;		/* Section */

		g = v[ord];			/* Parameter */

		nsec = ord + 1;		/* Increase sections for each order */

		vv *= (double)nsec;

		sec = floor(vv);
		if (((int)sec) & 1) {
			g = -g;				/* Alternate action in each section */
		}
		vv -= sec;
		if (g >= 0.0) {
			double tt = g - g * vv + 1.0;
			dsv = (vv * vv - vv)/(tt * tt);
			ddv = (g + 1.0)/(tt * tt);
			vv = vv/tt;
		} else {
			double tt = 1.0 - g * vv;
			dsv = (vv * vv - vv)/(tt * tt);
			ddv = (1.0 - g)/(tt * tt);
			vv = (vv - g * vv)/tt;
		}

		vv += sec;
		vv /= (double)nsec;
		dsv /= (double)nsec;
		if (((int)sec) & 1)
			dsv = -dsv;

		dv[ord] = dsv;
		for (i = ord - 1; i >= 0; i--)
			dv[i] *= ddv;
		din *= ddv;
	}

	*pdin = din;
	return vv;
}

/* Transfer function with offset and scale, and */
/* partial derivative with respect to the */
/* parameters and the input value. */
double icxdpdiSTransFunc(
double *v,			/* Pointer to first parameter */
double *dv,			/* Return derivative wrt each parameter */
double *pdin,		/* Return derivative wrt source value */
int    luord,		/* Number of parameters */
double vv,			/* Source of value */
double min,			/* Scale values */
double max
) {
	int i;
	max -= min;

	vv = (vv - min)/max;
	vv = icxdpdiTransFunc(v, dv, pdin, luord, vv);
	vv = (vv * max) + min;
	for (i = 0; i < luord; i++)
		dv[i] *= max;
	return vv;
}

/* ------------------------------------------------------ */
/* Matrix cube interpolation, used for device modelling. */
/* Including partial derivative for input and parameters. */

/* Matrix cube interpolation - interpolate between 2^di output corner values. */
/* Parameters are assumed to be fdi groups of 2^di parameters. */
void icxCubeInterp(
double *v,			/* Pointer to first parameter */
int    fdi,			/* Number of output channels */
int    di,			/* Number of input channels */
double *out,		/* Resulting fdi values */
double *in			/* Input di values */
) {
	int e, f, g;
	double gw[1 << MXDI];		/* weight for each matrix grid cube corner */

	/* Compute corner weights needed for interpolation */
	gw[0] = 1.0;
	for (e = 0, g = 1; e < di; e++, g *= 2) {
		int i; 
		for (i = 0; i < g; i++) {
			gw[g+i] = gw[i] * in[e];
			gw[i] *= (1.0 - in[e]);
		}
	}

	/* Now compute the output values */
	for (f = 0; f < fdi; f++) {
		out[f] = 0.0;						/* For each output value */
		for (e = 0; e < (1 << di); e++) {	/* For all corners of cube */
			out[f] += gw[e] * *v;
			v++;
		}
	}
}


/* Matrix cube interpolation. with partial derivative */
/* with respect to the input and parameters. */
void icxdpdiCubeInterp(
double *v,			/* Pointer to first parameter value [fdi * 2^di] */
double *dv,			/* Return [2^di] deriv. wrt each parameter v */
double *din,		/* Return [fdi * di] deriv. wrt each input value */
int    fdi,			/* Number of output channels */
int    di,			/* Number of input channels */
double *out,		/* Resulting fdi values */
double *in			/* Input di values */
) {
	int e, ee, f, g;
	int dip2 = (1 << di);
	double gw[1 << MXDI];			/* weight for each matrix grid cube corner */

	/* Compute corner weights needed for interpolation */
	gw[0] = 1.0;
	for (e = 0, g = 1; e < di; e++, g *= 2) {
		int i; 
		for (i = 0; i < g; i++) {
			gw[g+i] = gw[i] * in[e];
			gw[i] *= (1.0 - in[e]);
		}
	}

	/* Now compute the output values */
	for (f = 0; f < fdi; f++) {
		out[f] = 0.0;						/* For each output value */
		for (ee = 0; ee < dip2; ee++) {	/* For all corners of cube */
			out[f] += gw[ee] * v[f * dip2 + ee];
		}
	}

	/* Copy del for parameter to return array */
	for (ee = 0; ee < dip2; ee++) {	/* For all other corners of cube */
		dv[ee] = gw[ee];					/* del from parameter */
	}

	/* Compute del from in[] value we want */
	for (e = 0; e < di; e++) {				/* For input we want del wrt */

		for (f = 0; f < fdi; f++)
			din[f * di + e] = 0.0;			/* Zero del ready of accumulation */

		for (ee = 0; ee < dip2; ee++) {		/* For all corners of cube weights, */
			int e2;							/* accumulate del from in[] we want. */
			double vv = 1.0;

			/* Compute in[] weighted cube corners for all except del of in[] we want */
			for (e2 = 0; e2 < di; e2++) {			/* const from non del inputs */
				if (e2 == e)
					continue;
				if (ee & (1 << e2))
					vv *= in[e2];
				else
					vv *= (1.0 - in[e2]);
			}

			/* Accumulate contribution of in[] we want for corner to out[] we want */
			if (ee & (1 << e)) {
				for (f = 0; f < fdi; f++)
					din[f * di + e] += v[f * dip2 + ee] * vv;
			} else {
				for (f = 0; f < fdi; f++)
					din[f * di + e] -= v[f * dip2 + ee] * vv;
			}
		}
	}
}

/* ------------------------------------------------------ */
/* Matrix cube simplex interpolation, used for device modelling. */

/* Matrix cube simplex interpolation - interpolate between 2^di output corner values. */
/* Parameters are assumed to be fdi groups of 2^di parameters. */
void icxCubeSxInterp(
double *v,			/* Pointer to first parameter */
int    fdi,			/* Number of output channels */
int    di,			/* Number of input channels */
double *out,		/* Resulting fdi values */
double *in			/* Input di values */
) {
	int    si[MAX_CHAN];		/* in[] Sort index, [0] = smalest */

// ~~999
//{
//	double tout[MXDO];
//
//	icxCubeInterp(v, fdi, di, tout, in);
//printf("\n~1 Cube interp result = %f\n",tout[0]);
//}

//printf("~1 icxCubeSxInterp: %f %f %f\n", in[0], in[1], in[2]);
	/* Do insertion sort on coordinates, smallest to largest. */
	{
		int ff, vf;
		unsigned int e;
		double v;
		for (e = 0; e < di; e++)
			si[e] = e;						/* Initial unsorted indexes */

		for (e = 1; e < di; e++) {
			ff = e;
			v = in[si[ff]];
			vf = ff;
			while (ff > 0 && in[si[ff-1]] > v) {
				si[ff] = si[ff-1];
				ff--;
			}
			si[ff] = vf;
		}
	}
//printf("~1 sort order %d %d %d\n", si[0], si[1], si[2]);
//printf("         from %f %f %f\n", in[si[0]], in[si[1]], in[si[2]]);

	/* Now compute the weightings, simplex vertices and output values */
	{
		unsigned int e, f;
		double w;		/* Current vertex weight */

		w = 1.0 - in[si[di-1]];		/* Vertex at base of cell */
		for (f = 0; f < fdi; f++) {
			out[f] = w * v[f * (1 << di)];
//printf("~1 out[%d] = %f = %f * %f\n",f,out[f],w,v[f * (1 << di)]);
		}

		for (e = di-1; e > 0; e--) {	/* Middle verticies */
			w = in[si[e]] - in[si[e-1]];
			v += (1 << si[e]);				/* Move to top of cell in next largest dimension */
			for (f = 0; f < fdi; f++) {
				out[f] += w * v[f * (1 << di)];
//printf("~1 out[%d] = %f += %f * %f\n",f,out[f],w,v[f * (1 << di)]);
			}
		}

		w = in[si[0]];
		v += (1 << si[0]);		/* Far corner from base of cell */
		for (f = 0; f < fdi; f++) {
			out[f] += w * v[f * (1 << di)];
//printf("~1 out[%d] = %f += %f * %f\n",f,out[f],w,v[f * (1 << di)]);
		}
	}
}

/* ------------------------------------------------------ */
/* Matrix multiplication, used for device modelling. */
/* Including partial derivative for input and parameters. */


/* 3x3 matrix multiplication, with the matrix in a 1D array */
/* with respect to the input and parameters. */
void icxMulBy3x3Parm(
	double out[3],			/* Return input multiplied by matrix */
	double mat[9],			/* Matrix organised in [slow][fast] order */
	double in[3]			/* Input values */
) {
	double *v = mat, ov[3];
	int e, f;

	/* Compute the output values */
	for (f = 0; f < 3; f++) {
		ov[f] = 0.0;						/* For each output value */
		for (e = 0; e < 3; e++) {
			ov[f] += *v++ * in[e];
		}
	}
	out[0] = ov[0];
	out[1] = ov[1];
	out[2] = ov[2];
}


/* 3x3 matrix multiplication, with partial derivatives */
/* with respect to the input and parameters. */
void icxdpdiMulBy3x3Parm(
	double out[3],			/* Return input multiplied by matrix */
	double dv[3][9],		/* Return deriv for each [output] with respect to [param] */
	double din[3][3],		/* Return deriv for each [output] with respect to [input] */
	double mat[9],			/* Matrix organised in [slow][fast] order */
	double in[3]			/* Input values */
) {
	double *v, ov[3];
	int e, f;

	/* Compute the output values */
	v = mat;
	for (f = 0; f < 3; f++) {
		ov[f] = 0.0;						/* For each output value */
		for (e = 0; e < 3; e++) {
			ov[f] += *v++ * in[e];
		}
	}

	/* Compute deriv. with respect to the matrix parameter % 3 */
	/* This is pretty simple for a matrix ... */
	for (f = 0; f < 3; f++) {
		for (e = 0; e < 9; e++) {
			if (e/3 ==  f)
				dv[f][e] = in[e % 3];
			else
				dv[f][e] = 0.0;
		}
	}
	
	/* Compute deriv. with respect to the input values */
	/* This is pretty simple for a matrix ... */
	v = mat;
	for (f = 0; f < 3; f++)
		for (e = 0; e < 3; e++)
			din[f][e] = *v++;

	out[0] = ov[0];
	out[1] = ov[1];
	out[2] = ov[2];
}

#undef stricmp























