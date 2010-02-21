
/* 
 * Color Correct a TIFF file, using an ICC Device link profile.
 * Version #2, that allows an arbitrary string of profiles, and
 * copes with TIFF L*a*b* input and output.
 *
 * Author:  Graeme W. Gill
 * Date:    29/5/2004
 * Version: 2.00
 *
 * Copyright 2000 - 2006 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* TTBD:

	Should make -d option default to the last device
	profile.

	Should add support for getting TIFF ink names from
	the ICC profile colorantTable if it exists,
	or guessing them if they don't, and
	then matching them to the incoming TIFF, and
	embedding them in the outgoing TIFF.
	Add flag to ignore inkname mismatches.


	Should add support for transfering any extra alpha
	planes from input to output, rather than simply ignoring them.


	Question: Should this be changed to also function as
	          a dedicated simple linker, capable of outputing
	          a device link formed from a sequence ?
	          If argyll functionality was properly modularized,
	          it would be possible to have a single arbitrary
	          smart link sequence for both purposes.


	There's the sugggestion that the CIELab and ICCLab encodings
	have different white points (D65 and D50 respecively -
	see <http://www.asmail.be/msg0055212264.html>). Should
    we convert the white point to D65 for CIELab, or make this
	an option ? Probably a bad idea if we regard the Lab in/out
	as relative colorimetric representation ??

	Ideally should automatically generate optimized per channel
	input and output curves, rather than depending on
	reasonable behaviour from the profiles.

	I don't think that the idea of an XYZ input/output space
	has been properly implemented, since the TIFF format
	doesn't have an encoding for it, and hence the l2y_curve
	and u2l_curve scaling is probably incorrect.

	Forcing XYZ or Lab input & output spaces to be [io]combined
	is also not actually necessary, if the necessary profile
	curves were applied and the colorspace ranges properly allowed for.

 */

/*
	This program is a framework that exercises the
	IMDI code, as well as a demonstration of profile linking.
	It can also do the raster data conversion using the
    floating point code in ICCLIB as a reference.

 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "tiffio.h"
#include "icc.h"
#include "xicc.h"
#include "imdi.h"

#undef DEBUG		/* Print detailed debug info */


void usage(char *diag, ...) {
	fprintf(stderr,"Color Correct a TIFF file using any sequence of ICC profiles or Calibrations, V%s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: cctiff [-options] { [-i intent] profile%s | calbrtn.cal ...} infile.tif outfile.tif\n",ICC_FILE_EXT);
	fprintf(stderr," -v              Verbose.\n");
	fprintf(stderr," -c              Combine linearisation curves into one transform.\n");
	fprintf(stderr," -p              Use slow precise correction.\n");
	fprintf(stderr," -k              Check fast result against precise, and report.\n");
	fprintf(stderr," -r n            Override the default CLUT resolution\n");
	fprintf(stderr," -t n            Choose TIFF output encoding from 1..n\n");
	fprintf(stderr," -a              Read and Write planes > 4 as alpha planes\n");
	fprintf(stderr," -I              Ignore any file or profile colorspace mismatches\n");
	fprintf(stderr," -D              Don't append or set the output TIFF description\n");
	fprintf(stderr," -e profile.[%s | tiff]  Optionally embed a profile in the destination TIFF file.\n",ICC_FILE_EXT_ND);
	fprintf(stderr,"\n");
	fprintf(stderr,"                 Then for each profile in sequence:\n");
	fprintf(stderr,"   -i intent       p = perceptual, r = relative colorimetric,\n");
	fprintf(stderr,"                   s = saturation, a = absolute colorimetric\n");
	fprintf(stderr,"   -o order        n = normal (priority: lut > matrix > monochrome)\n");
	fprintf(stderr,"                   r = reverse (priority: monochrome > matrix > lut)\n");
	fprintf(stderr,"   profile.[%s | tiff]  Device, Link or Abstract profile\n",ICC_FILE_EXT_ND);
	fprintf(stderr,"                   ( May be embedded profile in TIFF file)\n");
	fprintf(stderr,"                 or each calibration file in sequence:\n");
	fprintf(stderr,"   -d dir          f = forward cal. (default), b = backwards cal.\n");
	fprintf(stderr,"   calbrtn.cal     Device calibration file.\n");
	fprintf(stderr,"\n");
	fprintf(stderr," infile.tif      Input TIFF file in appropriate color space\n");
	fprintf(stderr," outfile.tif     Output TIFF file\n");
	exit(1);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Conversion functions from direct binary 0..n^2-1 == 0.0 .. 1.0 range */
/* to ICC luo input range, and the reverse. */
/* Note that all these functions are per-component, */
/* so the can be included in per-component input or output curves, */
/* if they PCS values were rescaled to be within the range 0.0 .. 1.0. */
/* Since we're not currently doing this, we always set i/ocpmbine */
/* if the input/output is PCS, so that real PCS values don't */
/* appear in the input/output curves. */

/* TIFF 8 bit CIELAB to standard L*a*b* */
/* Assume that a & b have been converted from signed to offset */
static void cvt_CIELAB8_to_Lab(double *out, double *in) {
	out[0] = in[0] * 100.0;
	out[1] = in[1] * 255.0 - 128.0;
	out[2] = in[2] * 255.0 - 128.0;
}

/* Standard L*a*b* to TIFF 8 bit CIELAB */
/* Assume that a & b will be converted from offset to signed */
static void cvt_Lab_to_CIELAB8(double *out, double *in) {
	out[0] = in[0] / 100.0;
	out[1] = (in[1] + 128.0) * 1.0/255.0;
	out[2] = (in[2] + 128.0) * 1.0/255.0;
}


/* TIFF 16 bit CIELAB to standard L*a*b* */
/* Assume that a & b have been converted from signed to offset */
static void cvt_CIELAB16_to_Lab(double *out, double *in) {
	out[0] = in[0] * 100.0;
	out[1] = (in[1] - 32768.0/65535.0) * 256.0;
	out[2] = (in[2] - 32768.0/65535.0) * 256.0;
}

/* Standard L*a*b* to TIFF 16 bit CIELAB */
/* Assume that a & b will be converted from offset to signed */
static void cvt_Lab_to_CIELAB16(double *out, double *in) {
	out[0] = in[0] / 100.0;
	out[1] = in[1]/256.0 + 32768.0/65535.0;
	out[2] = in[2]/256.0 + 32768.0/65535.0;
}


/* TIFF 8 bit ICCLAB to standard L*a*b* */
static void cvt_ICCLAB8_to_Lab(double *out, double *in) {
	out[0] = in[0] * 100.0;
	out[1] = (in[1] * 255.0) - 128.0;
	out[2] = (in[2] * 255.0) - 128.0;
}

/* Standard L*a*b* to TIFF 8 bit ICCLAB */
static void cvt_Lab_to_ICCLAB8(double *out, double *in) {
	out[0] = in[0] * 1.0/100.0;
	out[1] = (in[1] + 128.0) * 1.0/255.0;
	out[2] = (in[2] + 128.0) * 1.0/255.0;
}


/* TIFF 16 bit ICCLAB to standard L*a*b* */
static void cvt_ICCLAB16_to_Lab(double *out, double *in) {
	out[0] = in[0] * (100.0 * 65535.0)/65280.0;
	out[1] = (in[1] * (255.0 * 65535.0)/65280) - 128.0;
	out[2] = (in[2] * (255.0 * 65535.0)/65280) - 128.0;
}

/* Standard L*a*b* to TIFF 16 bit ICCLAB */
static void cvt_Lab_to_ICCLAB16(double *out, double *in) {
	out[0] = in[0] * 65280.0/(100.0 * 65535.0);
	out[1] = (in[1] + 128.0) * 65280.0/(255.0 * 65535.0);
	out[2] = (in[2] + 128.0) * 65280.0/(255.0 * 65535.0);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Convert a TIFF Photometric tag to an ICC colorspace. */
/* return 0 if not possible or applicable. */
icColorSpaceSignature 
TiffPhotometric2ColorSpaceSignature(
void (**ocvt)(double *out, double *in),	/* Return write conversion function, NULL if none */
void (**icvt)(double *out, double *in),	/* Return read conversion function, NULL if none */
int *smsk,			/* Return signed handling mask, 0x0 if none */
int pmtc,			/* Input TIFF photometric */
int bps,			/* Input Bits per sample */
int spp,			/* Input Samples per pixel */
int extra			/* Extra Samples per pixel, if any */
) {
	if (icvt != NULL)
		*icvt = NULL;		/* Default return values */
	if (ocvt != NULL)
		*ocvt = NULL;		/* Default return values */
	if (smsk != NULL)
		*smsk = 0x0;

//	if (extra > 0 && pmtc != PHOTOMETRIC_SEPARATED)
//		return 0x0;						/* We don't handle this */
		
	switch (pmtc) {
		case PHOTOMETRIC_MINISWHITE:	/* Subtractive Gray */
			return icSigGrayData;

		case PHOTOMETRIC_MINISBLACK:	/* Additive Gray */
			return icSigGrayData;

		case PHOTOMETRIC_RGB:
			return icSigRgbData;

		case PHOTOMETRIC_PALETTE:
			return 0x0;

		case PHOTOMETRIC_MASK:
			return 0x0;

		case PHOTOMETRIC_SEPARATED:
			/* Should look at the colorant names to figure out if this is CMY, CMYK */
			/* Should at least return both Cmy/3 or Cmyk/4 ! */
			switch(spp) {
				case 2:
					return icSig2colorData;
				case 3:
//					return icSig3colorData;
					return icSigCmyData;
				case 4:
//					return icSig4colorData;
					return icSigCmykData;
				case 5:
					return icSig5colorData;
				case 6:
					return icSig6colorData;
				case 7:
					return icSig7colorData;
				case 8:
					return icSig8colorData;
				case 9:
					return icSig9colorData;
				case 10:
					return icSig10colorData;
				case 11:
					return icSig11colorData;
				case 12:
					return icSig12colorData;
				case 13:
					return icSig13colorData;
				case 14:
					return icSig14colorData;
				case 15:
					return icSig15colorData;
			}

		case PHOTOMETRIC_YCBCR:
			return icSigYCbCrData;

		case PHOTOMETRIC_CIELAB:
			if (bps == 8) {
				if (icvt != NULL)
					*icvt = cvt_CIELAB8_to_Lab;
				if (ocvt != NULL)
					*ocvt = cvt_Lab_to_CIELAB8;
			} else {
				if (icvt != NULL)
					*icvt = cvt_CIELAB16_to_Lab;
				if (ocvt != NULL)
					*ocvt = cvt_Lab_to_CIELAB16;
			}
			*smsk = 0x6;				/* Treat a & b as signed */
			return icSigLabData;

		case PHOTOMETRIC_ICCLAB:
			if (bps == 8) {
				if (icvt != NULL)
					*icvt = cvt_ICCLAB8_to_Lab;
				if (ocvt != NULL)
					*ocvt = cvt_Lab_to_ICCLAB8;
			} else {
				if (icvt != NULL)
					*icvt = cvt_ICCLAB16_to_Lab;
				if (ocvt != NULL)
					*ocvt = cvt_Lab_to_ICCLAB16;
			}
			return icSigLabData;

		case PHOTOMETRIC_ITULAB:
			return 0x0;					/* Could add this with a conversion function */
										/* but have to allow for variable ITU gamut */
										/* (Tag 433, "Decode") */

		case PHOTOMETRIC_LOGL:
			return 0x0;					/* Could add this with a conversion function */

		case PHOTOMETRIC_LOGLUV:
			return 0x0;					/* Could add this with a conversion function */
	}
	return 0x0;
}


/* Convert an ICC colorspace to the corresponding possible TIFF Photometric tags. */
/* Return the number of matching tags, and 0 if there is no corresponding tag. */
int
ColorSpaceSignature2TiffPhotometric(
uint16 tags[10],				/* Pointer to return array, up to 10 */
icColorSpaceSignature cspace	/* Input ICC colorspace */
) {
	switch(cspace) {
		case icSigGrayData:
			tags[0] = PHOTOMETRIC_MINISBLACK;
			return 1;
		case icSigRgbData:
			tags[0] = PHOTOMETRIC_RGB;
			return 1;
		case icSigCmyData:
		case icSigCmykData:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;
		case icSigYCbCrData:
			tags[0] = PHOTOMETRIC_YCBCR;
			return 1;
		case icSigXYZData:
		case icSigLabData:
			tags[0] = PHOTOMETRIC_CIELAB;
			tags[1] = PHOTOMETRIC_ICCLAB;
#ifdef NEVER
			tags[2] = PHOTOMETRIC_ITULAB;
#endif
			return 2;

		case icSigLuvData:
		case icSigYxyData:		/* Could handle this with a conversion ?? */
		case icSigHsvData:
		case icSigHlsData:
			return 0;

		case icSig2colorData:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig3colorData:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig4colorData:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig5colorData:
		case icSigMch5Data:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig6colorData:
		case icSigMch6Data:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig7colorData:
		case icSigMch7Data:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig8colorData:
		case icSigMch8Data:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig9colorData:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig10colorData:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig11colorData:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig12colorData:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig13colorData:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig14colorData:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		case icSig15colorData:
			tags[0] = PHOTOMETRIC_SEPARATED;
			return 1;

		default:
			return 0;
	}
	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Compute the length of a double nul terminated string, including */
/* the nuls. */
static int zzstrlen(char *s) {
	int i;
	for (i = 0;; i++) {
		if (s[i] == '\000' && s[i+1] == '\000')
			return i+2;
	}
	return 0;
}

/* Convert an ICC colorspace to the corresponding TIFF Inkset tag */
/* return 0xffff if not possible or applicable. */
int
ColorSpaceSignature2TiffInkset(
icColorSpaceSignature cspace,
int *len,						/* Return length of ASCII inknames */
char **inknames					/* Return ASCII inknames if non NULL */
) {
	switch(cspace) {
		case icSigCmyData:
			return 0xffff;	// ~~9999
			if (inknames != NULL) {
				*inknames = "cyan\000magenta\000yellow\000\000";
				*len = zzstrlen(*inknames);
			}
			return 0;				/* Not CMYK */
		case icSigCmykData:
			if (inknames != NULL) {
				*inknames = NULL;	/* No inknames */
				*len = 0;
			}
			return INKSET_CMYK;

		case icSigGrayData:
		case icSigRgbData:
		case icSigYCbCrData:
		case icSigLabData:
		case icSigXYZData:
		case icSigLuvData:
		case icSigYxyData:
		case icSigHsvData:
		case icSigHlsData:
		case icSig2colorData:
		case icSig3colorData:
		case icSig4colorData:
		case icSig5colorData:
		case icSigMch5Data:
			return 0xffff;

		case icSig6colorData:
		case icSigMch6Data:
				/* This is a cheat and a hack. Should really make use of the */
				/* ColorantTable to determine the colorant names. */
				/* allowing cctiff to read it. */
			if (inknames != NULL) {
				*inknames = "cyan\000magenta\000yellow\000black\000orange\000green\000\000";
				*len = zzstrlen(*inknames);
			}
			return 0;				/* Not CMYK */

		case icSig7colorData:
		case icSigMch7Data:
			return 0xffff;

		case icSig8colorData:
		case icSigMch8Data:
				/* This is a cheat and a hack. Should really make use of the */
				/* ColorantTable to determine the colorant names. */
				/* allowing cctiff to read it. */
			if (inknames != NULL) {
				*inknames = "cyan\000magenta\000yellow\000black\000orange\000green\000lightcyan\000lightmagenta\000\000";
				*len = zzstrlen(*inknames);
			}
			return 0;				/* Not CMYK */
		case icSig9colorData:
		case icSig10colorData:
		case icSig11colorData:
		case icSig12colorData:
		case icSig13colorData:
		case icSig14colorData:
		case icSig15colorData:
		default:
			return 0xffff;
	}
	return 0xffff;
}

char *
Photometric2str(
int pmtc
) {
	static char buf[80];
	switch (pmtc) {
		case PHOTOMETRIC_MINISWHITE:
			return "Subtractive Gray";
		case PHOTOMETRIC_MINISBLACK:
			return "Additive Gray";
		case PHOTOMETRIC_RGB:
			return "RGB";
		case PHOTOMETRIC_PALETTE:
			return "Indexed";
		case PHOTOMETRIC_MASK:
			return "Transparency Mask";
		case PHOTOMETRIC_SEPARATED:
			return "Separated";
		case PHOTOMETRIC_YCBCR:
			return "YCbCr";
		case PHOTOMETRIC_CIELAB:
			return "CIELab";
		case PHOTOMETRIC_ICCLAB:
			return "ICCLab";
		case PHOTOMETRIC_ITULAB:
			return "ITULab";
		case PHOTOMETRIC_LOGL:
			return "CIELog2L";
		case PHOTOMETRIC_LOGLUV:
			return "CIELog2Luv";
	}
	sprintf(buf,"Unknown Photometric Tag %d",pmtc);
	return buf;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define YSCALE (2.0/1.3)

/* Extra non-linearity applied to BtoA XYZ PCS */
/* This distributes the LUT indexes more evenly in */
/* perceptual space, greatly improving the B2A accuracy of XYZ LUT */
/* Since typically XYZ doesn't use the full range of 0-2.0 allowed */
/* for in the encoding, we scale the cLUT index values to use the 0-1.3 range */

/* (For these functions the encoded XYZ 0.0 - 2.0 range is 0.0 - 1.0 ??) */

/* Y to L* */
static void y2l_curve(double *out, double *in, int isXYZ) {
	int i;
	double val;
	double isc = 1.0, osc = 1.0;

	/* Scale from 0.0 .. 1.999969 to 0.0 .. 1.0 and back */
	/* + range adjustment */
	if (isXYZ) {
		isc = 32768.0/65535.0 * YSCALE;
		osc = 65535.0/32768.0;
	}

	for (i = 0; i < 3; i++) {
		val = in[i] * isc;
		if (val > 0.008856451586)
			val = 1.16 * pow(val,1.0/3.0) - 0.16;
		else
			val = 9.032962896 * val;
		if (val > 1.0)
			val = 1.0;
		out[i] = val * osc;
	}
}

/* L* to Y */
static void l2y_curve(double *out, double *in, int isXYZ) {
	int i;
	double val;
	double isc = 1.0, osc = 1.0;

	/* Scale from 0.0 .. 1.999969 to 0.0 .. 1.0 and back */
	/* + range adjustment */
	if (isXYZ) {
		isc = 32768.0/65535.0;
		osc = 65535.0/32768.0 / YSCALE;
	}

	/* Use an L* like curve, scaled to the maximum XYZ value */
	for (i = 0; i < 3; i++) {
		val = in[i] * isc;
		if (val > 0.08)
			val = pow((val + 0.16)/1.16, 3.0);
		else
			val = val/9.032962896;
		out[i] = val * osc;
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Callbacks used to initialise imdi */


/* Information needed from a single profile */
struct _profinfo {
	char name[MAXNAMEL+1];
	icc *c;								/* If non-NULL, ICC profile. */
	xcal *cal;							/* if non-NULL, xcal rather than profile */

	/* Valid for both icc and xcal: */
	icColorSpaceSignature ins, outs;	/* Colorspace of conversion */
	int id, od;							/* Dimensions of conversion */

	/* Valid only for icc: */
	icmHeader *h;						
	icRenderingIntent intent;			/* Rendering intent chosen */
	icmLookupFunc func;					/* Type of function to use in lookup */
	icmLookupOrder order;				/* tag search order to use */
	icmLuAlgType alg;					/* Type of lookup algorithm used */
	int clutres;						/* If this profile uses a clut, what's it's res. ? */
	icColorSpaceSignature natpcs;		/* Underlying natural PCS */
	icmLuBase *luo;						/* Base Lookup type object */

}; typedef struct _profinfo profinfo;

/* Context for imdi setup callbacks */
typedef struct {
	/* Overall parameters */
	int verb;				/* Non-zero if verbose */
	icColorSpaceSignature ins, outs;	/* Input/Output spaces */
	int id, od, md;			/* Input/Output dimensions and max(id,od) */
	int width, height;		/* Width and heigh of raster */
	int isign_mask;			/* Input sign mask */
	int osign_mask;			/* Output sign mask */
	int icombine;			/* Non-zero if input curves are to be combined */
	int ocombine;			/* Non-zero if output curves are to be combined */
	int ilcurve;			/* 1 if input curves are to be concatenated with Y like ->L* curve */
							/* (2 if input curves are to be concatenated with Y->L* curve) */
	int olcurve;			/* 1 if output curves are to be concatenated with L*->Y like curve */
							/* (2 if output curves are to be concatenated with L*->Y curve) */
	void (*icvt)(double *out, double *in);	/* If non-NULL, Input format conversion */
	void (*ocvt)(double *out, double *in);	/* If non-NULL, Output format conversion */

	int nprofs;				/* Number of profiles in the sequence */
	int first, last;		/* index of first and last profiles/cals, 0 and nprofs-1 */ 
	int fclut, lclut;		/* first and last profiles/cals that are part of multi-d table */
	profinfo *profs;		/* Profile information */
} sucntx;


/* Input curve function */
static void input_curves(
	void *cntx,
	double *out_vals,
	double *in_vals
) {
	int i;
	sucntx *rx = (sucntx *)cntx;

//printf("~1 incurve in %f %f %f %f\n",in_vals[0],in_vals[1],in_vals[2],in_vals[3]);

	for (i = 0; i < rx->id; i++)
		out_vals[i] = in_vals[i];	/* Default to nothing */

	if (rx->icombine == 0) {		/* Not combined into multi-d table */

		/* TIFF input format conversion */
		if (rx->icvt != NULL) {					/* (Never used because PCS < 0.0 > 1.0) */
			rx->icvt(out_vals, out_vals);
		}

		/* Any concatinated input calibrations */
		for (i = rx->first; i < rx->fclut; i++) {
			if (rx->profs[i].func == icmFwd)
				rx->profs[i].cal->interp(rx->profs[i].cal, out_vals, out_vals);
			else
				rx->profs[i].cal->inv_interp(rx->profs[i].cal, out_vals, out_vals);
		}

		/* The input table of the first profile */
		/* (icombine is set if input is PCS) */
		if (rx->fclut <= rx->lclut) 
			rx->profs[rx->fclut].luo->lookup_in(rx->profs[rx->fclut].luo, out_vals, out_vals);
	}

	/* If input curve converts to Y like space, apply Y->L* curve */
	/* so as to index CLUT perceptually. */
	if (rx->ilcurve != 0) {
		y2l_curve(out_vals, out_vals, rx->ilcurve == 2);
	}
//printf("~1 incurve out %f %f %f %f\n",out_vals[0],out_vals[1],out_vals[2],out_vals[3]);
}

/* Multi-dim table function */
static void md_table(
void *cntx,
double *out_vals,
double *in_vals
) {
	sucntx *rx = (sucntx *)cntx;
	double vals[MAX_CHAN];
	icColorSpaceSignature prs;		/* Previous colorspace */	
	int i, j;
	
	for (i = 0; i < rx->id; i++)
		vals[i] = in_vals[i];		/* default is do nothing */

//printf("~1 md_table in %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);

	if (rx->ilcurve) {
		/* Apply L*->Y curve to compensate for curve applied after input curve */
		l2y_curve(vals, vals, rx->ilcurve == 2);
	}

	prs = rx->ins;

	/* If the input curves are being combined into clut: */
	if (rx->icombine != 0) {

		/* Any needed TIFF file format conversion */
		if (rx->icvt != NULL) {
			rx->icvt(vals, vals);
//printf("~1 md_table after icvt %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);
		}

		/* Any concatinated input calibrations */
		for (j = rx->first; j < rx->fclut; j++) {
			if (rx->profs[j].func == icmFwd)
				rx->profs[j].cal->interp(rx->profs[j].cal, vals, vals);
			else
				rx->profs[j].cal->inv_interp(rx->profs[j].cal, vals, vals);
		}
//printf("~1 md_table after in cals %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);
	}

	/* Do all the profile links in-between (if any) */
	for (j = rx->fclut; j <= rx->lclut; j++) {

		/* If it's a calibration */
		if (rx->profs[j].cal != NULL) {
			if (rx->profs[j].func == icmFwd)
				rx->profs[j].cal->interp(rx->profs[j].cal, vals, vals);
			else
				rx->profs[j].cal->inv_interp(rx->profs[j].cal, vals, vals);

		/* Else it's a profile */
		} else {
			/* Convert PCS for this profile */
			if (prs == icSigXYZData && rx->profs[j].ins == icSigLabData) {
				icmXYZ2Lab(&icmD50, vals, vals);
//printf("~1 md_table after XYZ2Lab %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);
			} else if (prs == icSigLabData && rx->profs[j].ins == icSigXYZData) {
				icmLab2XYZ(&icmD50, vals, vals);
//printf("~1 md_table after Lab2XYZ %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);
			}
		
			/* If first or last profile */
			if (j == rx->fclut || j == rx->lclut) {
				if (j != rx->fclut || rx->icombine) {
					rx->profs[j].luo->lookup_in(rx->profs[j].luo, vals, vals);
//printf("~1 md_table after input curve %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);
				}
				rx->profs[j].luo->lookup_core(rx->profs[j].luo, vals, vals);
//printf("~1 md_table after core %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);
				if (j != rx->lclut || rx->ocombine) {
					rx->profs[j].luo->lookup_out(rx->profs[j].luo, vals, vals);
				}
//printf("~1 md_table after output curve %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);
			/* Middle of chain */
			} else {
				rx->profs[j].luo->lookup(rx->profs[j].luo, vals, vals);
//printf("~1 md_table after middle %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);
			}
		}
		prs = rx->profs[j].outs;
	}

	/* convert last PCS to rx->outs PCS if needed */
	if (prs == icSigXYZData
	 && rx->outs == icSigLabData) {
		icmXYZ2Lab(&icmD50, vals, vals);
//printf("~1 md_table after out XYZ2Lab %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);
	} else if (prs == icSigLabData
	      && rx->outs == icSigXYZData) {
		icmLab2XYZ(&icmD50, vals, vals);
//printf("~1 md_table after out Lab2XYZ %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);
	}

	/* If the output curves are being combined into clut: */
	if (rx->ocombine != 0) {

		/* Any concatinated output calibrations */
		for (j = rx->lclut+1; j <= rx->last; j++) {
			if (rx->profs[j].func == icmFwd)
				rx->profs[j].cal->interp(rx->profs[j].cal, vals, vals);
			else
				rx->profs[j].cal->inv_interp(rx->profs[j].cal, vals, vals);
		}
//printf("~1 md_table after out cals %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);

		/* Any needed TIFF file format conversion */
		if (rx->ocvt != NULL) {
			rx->ocvt(vals, vals);
//printf("~1 md_table after out ocvt %f %f %f %f\n",vals[0],vals[1],vals[2],vals[3]);
		}

	}

	if (rx->olcurve) {
		/* Add Y->L* curve to cause interpolation in perceptual space */
		y2l_curve(vals, vals, rx->olcurve == 2);
	}

	for (i = 0; i < rx->od; i++)
		out_vals[i] = vals[i];
//printf("~1 md_table returns %f %f %f %f\n",out_vals[0],out_vals[1],out_vals[2],out_vals[3]);
}

/* Output curve function */
static void output_curves(
	void *cntx,
	double *out_vals,
	double *in_vals
) {
	sucntx *rx = (sucntx *)cntx;
	int i; 
	
//printf("~1 outurve in %f %f %f %f\n",in_vals[0],in_vals[1],in_vals[2],in_vals[3]);
	for (i = 0; i < rx->od; i++)
		out_vals[i] = in_vals[i];

	/* Apply L* -> Y curve to undo curve applied at CLUT output. */
	if (rx->olcurve != 0) {
		l2y_curve(out_vals, out_vals, rx->olcurve == 2);
	}

	if (rx->ocombine == 0) {	/* Not combined into multi-d table */

		/* The output table of the last profile */
		/* (ocombine is set if output is PCS) */
		if (rx->fclut <= rx->lclut)  {
			rx->profs[rx->lclut].luo->lookup_out(rx->profs[rx->lclut].luo, out_vals, out_vals);
//printf("~1 md_table after out curve %f %f %f %f\n",out_vals[0],out_vals[1],out_vals[2],out_vals[3]);
		}

		/* Any concatinated output calibrations */
		for (i = rx->lclut+1; i <= rx->last; i++) {
			if (rx->profs[i].func == icmFwd)
				rx->profs[i].cal->interp(rx->profs[i].cal, out_vals, out_vals);
			else
				rx->profs[i].cal->inv_interp(rx->profs[i].cal, out_vals, out_vals);
		}

		/* Any needed file format conversion */
		if (rx->ocvt != NULL) {					/* (Never used because PCS < 0.0 > 1.0) */
			rx->ocvt(out_vals, out_vals);
//printf("~1 md_table after out ocvt %f %f %f %f\n",out_vals[0],out_vals[1],out_vals[2],out_vals[3]);
		}
	}
//printf("~1 outurve out %f %f %f %f\n",out_vals[0],out_vals[1],out_vals[2],out_vals[3]);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


/* Check whether two colorspaces appear compatible */
/* return NZ if they match, Z if they don't. */
/* Compatible means any PCS == any PCS, or exact match */
int CSMatch(icColorSpaceSignature s1, icColorSpaceSignature s2) {
	if (s1 == s2)
		return 1;

	if ((s1 == icSigXYZData || s1 == icSigLabData)
	 && (s2 == icSigXYZData || s2 == icSigLabData))
		return 1;

	if ((s1 == icSig5colorData || s1 == icSigMch5Data)
	 && (s2 == icSig5colorData || s2 == icSigMch5Data))
		return 1;

	if ((s1 == icSig6colorData || s1 == icSigMch6Data)
	 && (s2 == icSig6colorData || s2 == icSigMch6Data))
		return 1;

	if ((s1 == icSig7colorData || s1 == icSigMch7Data)
	 && (s2 == icSig7colorData || s2 == icSigMch7Data))
		return 1;

	if ((s1 == icSig8colorData || s1 == icSigMch8Data)
	 && (s2 == icSig8colorData || s2 == icSigMch8Data))
		return 1;

	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

int
main(int argc, char *argv[]) {
	int fa, nfa;							/* argument we're looking at */
	char in_name[MAXNAMEL+1] = "";			/* Input raster file name */
	char out_name[MAXNAMEL+1] = "";			/* Output raster file name */
	char dst_pname[MAXNAMEL+1] = "";		/* Destination embedded profile file name */
	icc *deicc = NULL;						/* Destination embedded profile (if any) */
	icRenderingIntent next_intent;			/* Rendering intent for next profile */
	icmLookupOrder next_order;				/* tag search order for next profile */
	icmLookupFunc next_func;				/* Direction for next calibration */
	int last_dim;							/* Next dimentionality between conversions */
	icColorSpaceSignature last_colorspace;	/* Next colorspace between conversions */
	char *last_cs_file;		/* Name of the file the last colorspace came from */
	int doimdi = 1;			/* Use the fast overall integer conversion */
	int dofloat = 0;		/* Use the slow precice (float). */
	int check = 0;			/* Check fast (int) against slow (float) */
	int ochoice = 0;		/* Output photometric choice 1..n */
	int alpha = 0;			/* Use alpha for extra planes */
	int ignoremm = 0;		/* Ignore any colorspace mismatches */
	int nodesc = 0;			/* Don't append or set the description */
	int i, j, rv = 0;

	/* TIFF file info */
	TIFF *rh = NULL, *wh = NULL;

	int x, y, width, height;					/* Common size of image */
	uint16 bitspersample;						/* Bits per sample */
	uint16 resunits;
	float resx, resy;
	uint16 pconfig;								/* Planar configuration */

	uint16 rsamplesperpixel, wsamplesperpixel;	/* Channels per sample */
	uint16 rphotometric, wphotometric;			/* Photometrics */
	uint16 rextrasamples, wextrasamples;		/* Extra "alpha" samples */
	uint16 *rextrainfo, wextrainfo[MAX_CHAN];	/* Info about extra samples */
	char *rdesc = NULL;							/* Existing description */
	char *wdesc = "[ Color corrected by ArgyllCMS ]";	/* New description */

	tdata_t *inbuf = NULL, *outbuf = NULL, *precbuf = NULL;

	/* IMDI */
	imdi *s = NULL;
	sucntx su;				/* Setup context */
	unsigned char *inp[MAX_CHAN];
	unsigned char *outp[MAX_CHAN];
	int clutres = 0;		/* Default */

	/* Error check */
	int mxerr = 0;
	double avgerr = 0.0;
	double avgcount = 0.0;

	error_program = "cctiff";
	if (argc < 2)
		usage("Too few arguments");

	/* Set defaults */
	memset((void *)&su, 0, sizeof(sucntx));
	next_intent = icmDefaultIntent;
	next_func = icmFwd;
	next_order = icmLuOrdNorm;

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')	{	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage("Usage requested");

			/* Slow, Precise, not integer */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				doimdi = 0;
				dofloat = 1;
			}

			/* Combine per channel curves */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				su.icombine = 1;
				su.ocombine = 1;
			}

			/* Check curves */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				doimdi = 1;
				dofloat = 1;
				check = 1;
			}

			/* Use alpha planes for any over 4 */
			else if (argv[fa][1] == 'a' || argv[fa][1] == 'A') {
				alpha = 1;
			}

			/* CLUT resolution */
			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to -r flag");
				clutres = atoi(na);
				if (clutres < 2)
					usage("-r argument must be >= 2");
			}

			/* Output TIFF photometric choice */
			else if (argv[fa][1] == 't' || argv[fa][1] == 'T') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to -t flag");
				ochoice = atoi(na);
			}
			/* Destination TIFF embedded profile */
			else if (argv[fa][1] == 'e' || argv[fa][1] == 'E') {
				fa = nfa;
				if (na == NULL) usage("Expect profile name argument to -e flag");
				strncpy(dst_pname,na, MAXNAMEL); dst_pname[MAXNAMEL] = '\000';
			}

			/* Next profile Intent */
			else if (argv[fa][1] == 'i') {
				fa = nfa;
				if (na == NULL) usage("Missing argument to -i flag");
    			switch (na[0]) {
					case 'p':
					case 'P':
						next_intent = icPerceptual;
						break;
					case 'r':
					case 'R':
						next_intent = icRelativeColorimetric;
						break;
					case 's':
					case 'S':
						next_intent = icSaturation;
						break;
					case 'a':
					case 'A':
						next_intent = icAbsoluteColorimetric;
						break;
					default:
						usage("Unknown argument '%c' to -i flag",na[0]);
				}
			}

			/* Next profile search order */
			else if (argv[fa][1] == 'o') {
				fa = nfa;
				if (na == NULL) usage("Missing argument to -o flag");
    			switch (na[0]) {
					case 'n':
					case 'N':
						next_order = icmLuOrdNorm;
						break;
					case 'r':
					case 'R':
						next_order = icmLuOrdRev;
						break;
					default:
						usage("Unknown argument '%c' to -o flag",na[0]);
				}
			}

			/* Next calibraton direction */
			else if (argv[fa][1] == 'd') {
				fa = nfa;
				if (na == NULL) usage("Missing argument to -i flag");
    			switch (na[0]) {
					case 'f':
					case 'F':
						next_func = icmFwd;
						break;
					case 'b':
					case 'B':
						next_func = icmBwd;
						break;
					default:
						usage("Unknown argument '%c' to -d flag",na[0]);
				}
			}

			else if (argv[fa][1] == 'I')
				ignoremm = 1;

			else if (argv[fa][1] == 'D')
				nodesc = 1;

			/* Verbosity */
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				su.verb = 1;
			}

			else  {
				usage("Unknown flag '%c'",argv[fa][1]);
			}

		} else if (argv[fa][0] != '\000') {
			/* Get the next filename */
		
			if (su.nprofs == 0)
				su.profs = (profinfo *)malloc(sizeof(profinfo));
			else
				su.profs = (profinfo *)realloc(su.profs, (su.nprofs+1) * sizeof(profinfo));
			if (su.profs == NULL)
				error("Malloc failed in allocating space for profile info.");

			memset((void *)&su.profs[su.nprofs], 0, sizeof(profinfo));
			strncpy(su.profs[su.nprofs].name,argv[fa],MAXNAMEL);
			su.profs[su.nprofs].name[MAXNAMEL] = '\000';
			su.profs[su.nprofs].intent = next_intent;
			su.profs[su.nprofs].func = next_func;
			su.profs[su.nprofs].order = next_order;

			su.nprofs++;
			next_intent = icmDefaultIntent;
			next_func = icmFwd;
			next_order = icmLuOrdNorm;
		} else {
			break;
		}
	}

	/* The last two "profiles" are actually the input and output TIFF filenames */
	/* Unwind them */
	if (su.nprofs < 2)
		usage("Not enough arguments to specify input and output TIFF files");

	strncpy(out_name,su.profs[--su.nprofs].name, MAXNAMEL); out_name[MAXNAMEL] = '\000';
	strncpy(in_name,su.profs[--su.nprofs].name, MAXNAMEL); in_name[MAXNAMEL] = '\000';

	su.fclut = su.first = 0;
	su.lclut = su.last = su.nprofs-1;

	if (check && (!doimdi || !dofloat))
		error("Can't do check unless both integera and float processing are enabled");

/*

	Logic required:

	Discover input TIFF colorspace and set as (ICC) "next_space"
	Set any special input space encoding transform (ie. device, Lab flavour)

	For each profile:

		case abstract:
			set dir = fwd, intent = default
			check next_space == CIE
			next_space = CIE
	
		case dev link:
			set dir = fwd, intent = default
			check next_space == profile.in_devspace
			next_space = profile.out_devspace 

		case cal file:
			check next_space == cal.devspace
			next_space = cal.devspace 

		case colorspace/input/display/output:
			if colorspace
				set intent = default

			if next_space == CIE
				set dir = fwd
				next_space = profile.devspace
			else
				set dir = bwd
				check next_space == profile.devspace
				next_space = CIE
	
		create luo
	
	Make output TIFF colorspace match next_space

	Figure out how many calibrations can be concatinated into the input
	and output curves.

	Set any special output space encoding transform (ie. device, Lab flavour)

*/
	/* - - - - - - - - - - - - - - - */
	/* Open up input tiff file ready for reading */
	/* Discover input TIFF colorspace and set as (ICC) "last_colorspace" */
	/* Set any special input space encoding transform (ie. device, Lab flavour) */

	if ((rh = TIFFOpen(in_name, "r")) == NULL)
		error("error opening read file '%s'",in_name);

	TIFFGetField(rh, TIFFTAG_IMAGEWIDTH,  &width);
	TIFFGetField(rh, TIFFTAG_IMAGELENGTH, &height);

	TIFFGetField(rh, TIFFTAG_BITSPERSAMPLE, &bitspersample);
	if (bitspersample != 8 && bitspersample != 16) {
		error("TIFF Input file must be 8 or 16 bit/channel");
	}

	TIFFGetFieldDefaulted(rh, TIFFTAG_EXTRASAMPLES, &rextrasamples, &rextrainfo);
//	if (rextrasamples > 0 && alpha == 0)
//		error("TIFF Input file has extra samples per pixel - cctiff can't handle that");
		
	TIFFGetField(rh, TIFFTAG_PHOTOMETRIC, &rphotometric);
	TIFFGetField(rh, TIFFTAG_SAMPLESPERPIXEL, &rsamplesperpixel);

	/* Figure out how to handle the input TIFF colorspace */
	if ((su.ins = TiffPhotometric2ColorSpaceSignature(NULL, &su.icvt, &su.isign_mask, rphotometric,
	                                     bitspersample, rsamplesperpixel, rextrasamples)) == 0)
		error("Can't handle TIFF file photometric %s", Photometric2str(rphotometric));
	su.id = rsamplesperpixel;

	TIFFGetField(rh, TIFFTAG_PLANARCONFIG, &pconfig);
	if (pconfig != PLANARCONFIG_CONTIG)
		error ("TIFF Input file must be planar");

	TIFFGetField(rh, TIFFTAG_RESOLUTIONUNIT, &resunits);
	TIFFGetField(rh, TIFFTAG_XRESOLUTION, &resx);
	TIFFGetField(rh, TIFFTAG_YRESOLUTION, &resy);

	last_dim = su.id;
	last_colorspace = su.ins;
	last_cs_file = in_name;

	su.width = width;
	su.height = height;


	/* - - - - - - - - - - - - - - - */
	/* Check and setup the sequence of ICC profiles */

	/* For each profile in the sequence, configure it to transform the color */
	/* appropriately */
	for (i = su.first; i <= su.last; i++) {

		/* First see if it's a calibration file */
		if ((su.profs[i].cal = new_xcal()) == NULL)
			error("new_xcal failed");

		if ((su.profs[i].cal->read(su.profs[i].cal, su.profs[i].name)) == 0) {

			su.profs[i].ins = su.profs[i].outs = icx_colorant_comb_to_icc(su.profs[i].cal->devmask);
			if (su.profs[i].outs == 0)
				error ("Calibration file '%s' has unhandled device mask %s",su.profs[i].name,icx_inkmask2char(su.profs[i].cal->devmask,1));
			su.profs[i].id = su.profs[i].od = su.profs[i].cal->devchan;
			/* We use the user provided direction */

		/* else see if it's an ICC or embedded ICC */
		} else {
		
			su.profs[i].cal->del(su.profs[i].cal);	/* Clean up */
			su.profs[i].cal = NULL;

			if ((su.profs[i].c = read_embedded_icc(su.profs[i].name)) == NULL)
				error ("Can't read profile or calibration from file '%s'",su.profs[i].name);

			su.profs[i].h = su.profs[i].c->header;

			/* Deal with different profile classes, */
			/* and set the profile function and intent. */
			switch (su.profs[i].h->deviceClass) {
    		case icSigAbstractClass:
    		case icSigLinkClass:
					su.profs[i].func = icmFwd;
					su.profs[i].intent = icmDefaultIntent;
					break;

    		case icSigColorSpaceClass:
					su.profs[i].func = icmFwd;
					su.profs[i].intent = icmDefaultIntent;
					/* Fall through */

    		case icSigInputClass:
    		case icSigDisplayClass:
    		case icSigOutputClass:
					/* Note we don't handle an ambigious (both directions match) case. */
					/* We would need direction from the user to resolve this. */
					if (CSMatch(last_colorspace, su.profs[i].h->colorSpace)) {
						su.profs[i].func = icmFwd;
					} else {
						su.profs[i].func = icmBwd;		/* PCS -> Device */
					}
					break;
					/* Use the user provided intent */

				default:
					error("Can't handle deviceClass %s from file '%s'",
					     icm2str(icmProfileClassSignature,su.profs[i].h->deviceClass),
					     su.profs[i].c->err,su.profs[i].name);
			}

			/* Get a conversion object */
			if ((su.profs[i].luo = su.profs[i].c->get_luobj(su.profs[i].c, su.profs[i].func,
			              su.profs[i].intent, icmSigDefaultData, su.profs[i].order)) == NULL)
				error ("%d, %s from '%s'",su.profs[i].c->errc, su.profs[i].c->err, su.profs[i].name);
		
			/* Get details of conversion */
			su.profs[i].luo->spaces(su.profs[i].luo, &su.profs[i].ins, &su.profs[i].id,
			         &su.profs[i].outs, &su.profs[i].od, &su.profs[i].alg, NULL, NULL, NULL, NULL);

			/* Get native PCS space */
			su.profs[i].luo->lutspaces(su.profs[i].luo, NULL, NULL, NULL, NULL, &su.profs[i].natpcs);

			/* If this is a lut transform, find out its resolution */
			if (su.profs[i].alg == icmLutType) {
				icmLut *lut;
				icmLuLut *luluo = (icmLuLut *)su.profs[i].luo;		/* Safe to coerce */
				luluo->get_info(luluo, &lut, NULL, NULL, NULL);	/* Get some details */
				su.profs[i].clutres = lut->clutPoints;			/* Desired table resolution */
			} else 
				su.profs[i].clutres = 0;

		}

		/* Check that we can join to previous correctly */
		if (!ignoremm && !CSMatch(last_colorspace, su.profs[i].ins))
			error("Last colorspace %s from file '%s' doesn't match input space %s of profile %s",
		      icm2str(icmColorSpaceSignature,last_colorspace),
			  last_cs_file,
			  icm2str(icmColorSpaceSignature,su.profs[i].h->colorSpace),
			  su.profs[i].name);

		last_dim = icmCSSig2nchan(su.profs[i].outs);
		last_colorspace = su.profs[i].outs;
		last_cs_file = su.profs[i].name; 
	}
	
	su.od = last_dim;
	su.outs = last_colorspace;

	/* Go though the sequence again, and count the number of leading and */
	/* trailing calibrations that can be combined into the input and output */
	/* lookup curves */
	for (i = su.first; ; i++) {
		if (i > su.last || su.profs[i].c != NULL) {
			su.fclut = i;
			break;
		}
	}
	for (i = su.last; ; i--) {
		if (i < su.first || su.profs[i].c != NULL) {
			su.lclut = i;
			break;
		}
	}

	if (su.fclut > su.lclut) {	/* Hmm. All calibs, no profiles */
		su.fclut = su.first;	/* None at start */
		su.lclut = su.first-1;	/* All at the end */
	}
		
//printf("~1 first = %d, fclut = %d, lclut = %d, last = %d\n", su.first, su.fclut, su.lclut, su.last);

	su.md = su.id > su.od ? su.id : su.od;

	/* - - - - - - - - - - - - - - - */
	/* Open up the output TIFF file for writing */
	if ((wh = TIFFOpen(out_name, "w")) == NULL)
		error("Can\'t create TIFF file '%s'!",out_name);
	
	wsamplesperpixel = su.od;

	wextrasamples = 0;
	if (alpha && wsamplesperpixel > 4) {
		wextrasamples = wsamplesperpixel - 4;	/* Call samples > 4 "alpha" samples */
		for (j = 0; j < wextrasamples; j++)
			wextrainfo[j] = EXTRASAMPLE_UNASSALPHA;
	}

	/* Configure the output TIFF file appropriately */
	TIFFSetField(wh, TIFFTAG_IMAGEWIDTH,  width);
	TIFFSetField(wh, TIFFTAG_IMAGELENGTH, height);
	TIFFSetField(wh, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(wh, TIFFTAG_SAMPLESPERPIXEL, wsamplesperpixel);
	TIFFSetField(wh, TIFFTAG_BITSPERSAMPLE, bitspersample);
	TIFFSetField(wh, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

	TIFFSetField(wh, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	if (resunits) {
		TIFFSetField(wh, TIFFTAG_RESOLUTIONUNIT, resunits);
		TIFFSetField(wh, TIFFTAG_XRESOLUTION, resx);
		TIFFSetField(wh, TIFFTAG_YRESOLUTION, resy);
	}
		/* Perhaps the description could be more informative ? */
	if (TIFFGetField(rh, TIFFTAG_IMAGEDESCRIPTION, &rdesc) != 0) {
		char *ndesc = NULL;

		if ((ndesc = malloc(sizeof(char) * (strlen(rdesc) + strlen(wdesc) + 2))) == NULL)
			error("malloc failed on new desciption string");
		
		strcpy(ndesc, rdesc);
		if (nodesc == 0) {
			strcat(ndesc, " ");
			strcat(ndesc, wdesc);
		}
		TIFFSetField(wh, TIFFTAG_IMAGEDESCRIPTION, ndesc);
	
		free(ndesc);
	} else if (nodesc == 0) {
		TIFFSetField(wh, TIFFTAG_IMAGEDESCRIPTION, wdesc);
	}

	/* Lookup and decide what TIFF photometric suites the output colorspace */
	{
		int no_pmtc;					/* Number of possible photometrics */
		uint16 pmtc[10];				/* Photometrics of output file */
		if ((no_pmtc = ColorSpaceSignature2TiffPhotometric(pmtc,
		                                    last_colorspace)) == 0)
			error("TIFF file can't handle output colorspace '%s'!",
			      icm2str(icmColorSpaceSignature, last_colorspace));
	
		if (no_pmtc > 1) {		/* Need to choose a photometric */
			if (ochoice < 1 || ochoice > no_pmtc ) {
				printf("Possible photometrics (Output Encodings) for output colorspace %s are:\n",
				        icm2str(icmColorSpaceSignature,last_colorspace));
				for (i = 0; i < no_pmtc; i++)
					printf("%d: %s\n",i+1, Photometric2str(pmtc[i]));
				if (ochoice < 1 || ochoice > no_pmtc )
					error("An output photometric must be selected with the -t parameter");
				
			}
			wphotometric = pmtc[ochoice-1];
		} else {
			wphotometric = pmtc[0];
		}
	}

	/* Lookup what we need to handle this. */
	if ((su.outs = TiffPhotometric2ColorSpaceSignature(&su.ocvt, NULL, &su.osign_mask, wphotometric,
	                                     bitspersample, wsamplesperpixel, wextrasamples)) == 0)
		error("Can't handle TIFF file photometric %s", Photometric2str(wphotometric));
	TIFFSetField(wh, TIFFTAG_PHOTOMETRIC, wphotometric);

	if (alpha && wextrasamples > 0) {
		TIFFSetField(wh, TIFFTAG_EXTRASAMPLES, wextrasamples, wextrainfo);

	} else {

		if (wphotometric == PHOTOMETRIC_SEPARATED) {
			int iset;
			int inlen;
			char *inames;
			iset = ColorSpaceSignature2TiffInkset(su.outs, &inlen, &inames);
			if (iset != 0xffff && inlen > 0 && inames != NULL) {
				TIFFSetField(wh, TIFFTAG_INKSET, iset);
				if (inames != NULL) {
					TIFFSetField(wh, TIFFTAG_INKNAMES, inlen, inames);
				}
			}
		}
	}

	/* - - - - - - - - - - - - - - - */
	/* Setup any destination embedded profile */
	if (dst_pname[0] != '\000') {
		icmFile *fp;		/* Read fp for the profile */
		unsigned char *buf;
		int size;

		if ((deicc = read_embedded_icc(dst_pname)) == NULL)
			error("Unable to open profile for destination embedding '%s'",dst_pname);

		/* Check that it is compatible with the destination TIFF */
		if (deicc->header->deviceClass != icSigColorSpaceClass
		 && deicc->header->deviceClass != icSigInputClass
		 && deicc->header->deviceClass != icSigDisplayClass
		 && deicc->header->deviceClass != icSigOutputClass) {
			error("Destination embedded profile is wrong device class for embedding");
		}

		if (deicc->header->colorSpace != su.outs
		 || (deicc->header->pcs != icSigXYZData 
		  && deicc->header->pcs != icSigLabData)) {
			error("Destination embedded profile colorspaces don't match TIFF");
		}

		if ((fp = deicc->get_rfp(deicc)) == NULL)
			error("Failed to be able to read destination embedded profile");

		if ((size = fp->get_size(fp)) == 0)
			error("Failed to be able to get size of destination embedded profile");

		if ((buf = malloc(size)) == NULL)
			error("malloc failed on destination embedded profile size %d",size);

		if (fp->seek(fp,0))
			error("rewind on destination embedded profile failed");

		if (fp->read(fp, buf, 1, size) != size)
			error("reading destination embedded profile failed");

		/* (For iccv4 we would now fp->del(fp) because we got a reference) */

		if (TIFFSetField(wh, TIFFTAG_ICCPROFILE, size, buf) == 0)
			error("setting TIFF embedded ICC profiel field failed");

		free(buf);
		deicc->del(deicc);
	}

	/* - - - - - - - - - - - - - - - */
	if (su.fclut <= su.lclut
	 && ((su.profs[su.fclut].natpcs == icSigXYZData && su.profs[su.fclut].alg == icmMatrixFwdType))
	  || (su.profs[su.fclut].ins == icSigXYZData)) {
		su.ilcurve = 1;			/* Index CLUT with L* curve rather than Y */
	}

	/* Setup input/output curve use. */
	if (su.ins == icSigLabData || su.ins == icSigXYZData) {
		su.icombine = 1;		/* CIE can't be conveyed through 0..1 domain lookup */
	}

	if (su.fclut <= su.lclut
	 && ((su.profs[su.lclut].natpcs == icSigXYZData && su.profs[su.lclut].alg == icmMatrixBwdType))
	  || (su.profs[su.lclut].outs == icSigXYZData)) {
		su.olcurve = 1;			/* Interpolate in L* space rather than Y */
	}

	if (su.outs == icSigLabData || su.outs == icSigXYZData) {
		su.ocombine = 1;		/* CIE can't be conveyed through 0..1 domain lookup */
	}

	/* - - - - - - - - - - - - - - - */
	/* Report the connection sequence details */

	if (su.verb) {

		printf("Input TIFF file '%s'\n",in_name);
		printf("TIFF file colorspace is %s\n",icm2str(icmColorSpaceSignature,su.ins));
		printf("TIFF file photometric is %s\n",Photometric2str(rphotometric));
		printf("\n");

		printf("There are %d profiles/calibrations in the sequence:\n\n",su.nprofs);

		for (i = su.first; i <= su.last; i++) {
			if (su.profs[i].c != NULL) {
				icmFile *op;
				if ((op = new_icmFileStd_fp(stdout)) == NULL)
					error ("Can't open stdout");
				printf("Profile %d '%s':\n",i,su.profs[i].name);
				su.profs[i].h->dump(su.profs[i].h, op, 1);
				op->del(op);
				printf("Direction = %s\n",icm2str(icmTransformLookupFunc, su.profs[i].func));
				printf("Intent = %s\n",icm2str(icmRenderingIntent, su.profs[i].intent));
				printf("Algorithm = %s\n",icm2str(icmLuAlg, su.profs[i].alg));
			} else {
				printf("Calibration %d '%s':\n",i,su.profs[i].name);
				printf("Direction = %s\n",icm2str(icmTransformLookupFunc, su.profs[i].func));
				if (su.profs[i].cal->xpi.deviceMfgDesc != NULL)
					printf("Manufacturer: '%s'\n",su.profs[i].cal->xpi.deviceMfgDesc);
				if (su.profs[i].cal->xpi.modelDesc != NULL)
					printf("Model: '%s'\n",su.profs[i].cal->xpi.modelDesc);
				if (su.profs[i].cal->xpi.profDesc != NULL)
					printf("Description: '%s'\n",su.profs[i].cal->xpi.profDesc);
				if (su.profs[i].cal->xpi.copyright != NULL)
					printf("Copyright: '%s'\n",su.profs[i].cal->xpi.copyright);
			}

			if (i == 0 && su.icombine)
				printf("Input curves being combined\n");
			if (i == 0 && su.ilcurve)
				printf("Input curves being post-converted to L*\n");
			printf("Input space = %s\n",icm2str(icmColorSpaceSignature, su.profs[i].ins));
			printf("Output space = %s\n",icm2str(icmColorSpaceSignature, su.profs[i].outs));
			if (i == (su.last) && su.olcurve)
				printf("Output curves being pre-converted from L*\n");
			if (i == (su.last) && su.ocombine)
				printf("Output curves being combined\n");
			printf("\n");
		}

		printf("Output TIFF file '%s'\n",out_name);
		printf("TIFF file colorspace is %s\n",icm2str(icmColorSpaceSignature,su.outs));
		printf("TIFF file photometric is %s\n",Photometric2str(wphotometric));
		printf("\n");
	}

	/* - - - - - - - - - - - - - - - */
	/* Setup the imdi */

	if (check)
		doimdi = dofloat = 1;

	if (doimdi && su.nprofs > 0) {
		int aclutres = 0;	/* Automatically set res */
		imdi_options opts = opts_none;
	
		if (rextrasamples > 0) {		/* We need to skip the alpha */
			opts |= opts_istride;
		}
	
		/* Setup the imdi resolution */
		/* Choose the resolution from the highest lut resolution in the sequence, */
		/* or choose a default. */
		for (i = su.first; i <= su.last; i++) {
			if (su.profs[i].c != NULL
			 && su.profs[i].clutres > aclutres)
				aclutres = su.profs[i].clutres;
		}
		if (aclutres == 0) {
			aclutres = dim_to_clutres(su.id, 2);			/* High quality */

		} else if (aclutres < dim_to_clutres(su.id, 1)) {	/* Worse than medium */
			aclutres = dim_to_clutres(su.id, 1);
		}

		if (clutres == 0)
			clutres = aclutres;

		if (su.verb)
			printf("Using CLUT resolution %d\n",clutres);
	
		s = new_imdi(
			su.id,			/* Number of input dimensions */
			su.od,			/* Number of output dimensions */
							/* Input pixel representation */
			bitspersample == 8 ? pixint8 : pixint16,
							/* Output pixel representation */
			su.isign_mask,	/* Treat appropriate channels as signed */
			NULL,			/* No raster to callback channel mapping */
			prec_min,		/* Minimum of input and output precision */
			bitspersample == 8 ? pixint8 : pixint16,
			su.osign_mask,	/* Treat appropriate channels as signed */
			NULL,			/* No raster to callback channel mapping */
			clutres,		/* Desired table resolution */
			oopts_none,		/* Desired per channel output options */
			NULL,			/* Output channel check values */
			opts,			/* Desired processing direction and stride support */
			input_curves,	/* Callback functions */
			md_table,
			output_curves,
			(void *)&su		/* Context to callbacks */
		);
		
		if (s == NULL) {
	#ifdef NEVER
			printf("id = %d\n",su.id);
			printf("od = %d\n",su.od);
			printf("in bps = %d\n",bitspersample);
			printf("out bps = %d\n",bitspersample);
			printf("in signs = %d\n",su.isign_mask);
			printf("out signs = %d\n",su.osign_mask);
			printf("clutres = %d\n",clutres);
	#endif
			error("new_imdi failed");
		}
	}

	inbuf  = _TIFFmalloc(TIFFScanlineSize(rh));
	outbuf = _TIFFmalloc(TIFFScanlineSize(wh));
	inp[0] = (unsigned char *)inbuf;
	outp[0] = (unsigned char *)outbuf;
	if (dofloat || su.nprofs == 0)
		precbuf = _TIFFmalloc(TIFFScanlineSize(wh));


	/* - - - - - - - - - - - - - - - */
	/* Process colors to translate */
	/* (Should fix this to process a group of lines at a time ?) */

	for (y = 0; y < height; y++) {

		/* Read in the next line */
		if (TIFFReadScanline(rh, inbuf, y, 0) < 0)
			error ("Failed to read TIFF line %d",y);

		if (doimdi && su.nprofs > 0) {
			/* Do fast conversion */
			s->interp(s, (void **)outp, 0, (void **)inp, rsamplesperpixel, width);
		}
		
		if (dofloat || su.nprofs == 0) {
			/* Do floating point conversion into the precbuf[] */
			for (x = 0; x < width; x++) {
				int i;
				double in[MAX_CHAN], out[MAX_CHAN];
				
//printf("\n");
				if (bitspersample == 8) {
					for (i = 0; i < rsamplesperpixel; i++) {
						int v = ((unsigned char *)inbuf)[x * rsamplesperpixel + i];
//printf("~1 8 bit pixel value chan %d = %d\n",i,v);
						if (su.isign_mask & (1 << i))		/* Treat input as signed */
							v = (v & 0x80) ? v - 0x80 : v + 0x80;
//printf("~1 8 bit after treat as signed chan %d = %d\n",i,v);
						in[i] = v/255.0;
//printf("~1 8 bit fp chan %d value = %f\n",i,in[i]);
					}
				} else {
					for (i = 0; i < rsamplesperpixel; i++) {
						int v = ((unsigned short *)inbuf)[x * rsamplesperpixel + i];
//printf("~1 16 bit pixel value chan %d = %d\n",i,v);
						if (su.isign_mask & (1 << i))		/* Treat input as signed */
							v = (v & 0x8000) ? v - 0x8000 : v + 0x8000;
//printf("~1 16 bit after treat as signed chan %d = %d\n",i,v);
						in[i] = v/65535.0;
//printf("~1 16 bit fp chan %d value = %f\n",i,in[i]);
					}
				}

				if (su.nprofs > 0) {
					/* Apply the reference conversion */
					input_curves((void *)&su, out, in);
//for (i = 0; i < su.id; i++) printf("~1 after input curve chan %d = %f\n",i,out[i]);
					md_table((void *)&su, out, out);
//for (i = 0; i < su.od; i++) printf("~1 after md table chan %d = %f\n",i,out[i]);
					output_curves((void *)&su, out, out);
//for (i = 0; i < su.od; i++) printf("~1 after output curve chan %d = %f\n",i,out[i]);
				} else {
					for (i = 0; i < su.od; i++)
						 out[i] = in[i];
				}

				if (bitspersample == 8) {
					for (i = 0; i < su.od; i++) {
						int v = (int)(out[i] * 255.0 + 0.5);
//printf("~1 8 bit chan %d = %d\n",i,v);
						if (v < 0)
							v = 0;
						else if (v > 255)
							v = 255;
//printf("~1 8 bit after clip curve chan %d = %d\n",i,v);
						if (su.osign_mask & (1 << i))		/* Treat input as offset */
							v = (v & 0x80) ? v - 0x80 : v + 0x80;
//printf("~1 8 bit after treat as offset chan %d = %d\n",i,v);
						((unsigned char *)precbuf)[x * su.od + i] = v;
					}
				} else {
					for (i = 0; i < su.od; i++) {
						int v = (int)(out[i] * 65535.0 + 0.5);
//printf("~1 16 bit chan %d = %d\n",i,v);
						if (v < 0)
							v = 0;
						else if (v > 65535)
							v = 65535;
//printf("~1 16 bit after clip curve chan %d = %d\n",i,v);
						if (su.osign_mask & (1 << i))		/* Treat input as offset */
							v = (v & 0x8000) ? v - 0x8000 : v + 0x8000;
//printf("~1 16 bit after treat as offset chan %d = %d\n",i,v);
						((unsigned short *)precbuf)[x * su.od + i] = v;
					}
				}
			}

			if (check) {
				/* Compute the errors */
				for (x = 0; x < (width * su.od); x++) {
					int err;
					if (bitspersample == 8)
						err = ((unsigned char *)outbuf)[x] - ((unsigned char *)precbuf)[x];
					else
						err = ((unsigned short *)outbuf)[x] - ((unsigned short *)precbuf)[x];
					if (err < 0)
						err = -err;
					if (err > mxerr)
						mxerr = err;
					avgerr += (double)err;
					avgcount++;
				}
			}
		}
			
		if (dofloat || su.nprofs == 0) {	/* Use the results of the f.p. conversion */
			if (TIFFWriteScanline(wh, precbuf, y, 0) < 0)
				error ("Failed to write TIFF line %d",y);
		} else {
			if (TIFFWriteScanline(wh, outbuf, y, 0) < 0)
				error ("Failed to write TIFF line %d",y);
		}
	}

	if (check) {
		printf("Worst error = %d bits, average error = %f bits\n", mxerr, avgerr/avgcount);
		if (bitspersample == 8)
			printf("Worst error = %f%%, average error = %f%%\n",
			       mxerr/2.55, avgerr/(2.55 * avgcount));
		else
			printf("Worst error = %f%%, average error = %f%%\n",
			       mxerr/655.35, avgerr/(655.35 * avgcount));
	}


	if (inbuf != NULL) _TIFFfree(inbuf);
	if (outbuf != NULL) _TIFFfree(outbuf);
	if (precbuf != NULL) _TIFFfree(precbuf);

	/* Done with lookup object */
	if (s != NULL)
		s->del(s);

	TIFFClose(rh);		/* Close Input file */
	TIFFClose(wh);		/* Close Output file */

	/* Free up all the profiles etc. in the sequence. */
	for (i = 0; i < su.nprofs; i++) {
		if (su.profs[i].c != NULL) {				/* Has an ICC profile */
			su.profs[i].luo->del(su.profs[i].luo);	/* Lookup */
			su.profs[i].c->del(su.profs[i].c);	
		} else {
			su.profs[i].cal->del(su.profs[i].cal);	/* Calibration */
		}
	}

	return 0;
}

