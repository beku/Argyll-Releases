
/* 
 * Color Correct a TIFF file, using an ICC Device link profile.
 * Prototype of new version that allows an arbitrary string of profiles.
 *
 * Author:  Graeme W. Gill
 * Date:    29/5/2004
 * Version: 2.00
 *
 * Copyright 2000 - 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/*
 * Thanks to Neil Okamoto for the 16 bit TIFF mods.
 */

/* TTBD:

	Waiting on the implementation of the unified simple
	profile support ( Lut - Clut - Lut) to implement this.
	This has become entangled in ICCV4 support.

	Question: Should this be changed to also function as
	          a dedicated simple linker, capable of outputing
	          a device link ?
	          If argyll functionality was properly modularized,
	          it would be possible to have a single arbitrary
	          smart link chain for both purposes.
 */

/*
	This program is a framework that exercises the
	IMDI code, as well as a demonstration of simple
    profile linking.  It can also do the conversion using the
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
#include "tiffio.h"
#include "icc.h"
#include "imdi.h"

void error(char *fmt, ...), warning(char *fmt, ...);

void usage(char *mes) {
	fprintf(stderr,"Color Correct a TIFF file using an ICC device link profile, V%s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	if (mes != NULL)
		fprintf(stderr,"Error: '%s'\n",mes);
	fprintf(stderr,"usage: cctiff [-options] { [-i intent] profile.icm ...} infile.tif outfile.tif\n");
	fprintf(stderr," -v            Verbose.\n");
	fprintf(stderr," -c            Combine linearisation curves into one transform.\n");
	fprintf(stderr," -p            Use slow precise correction.\n");
	fprintf(stderr," -k            Check fast result against precise, and report.\n");
	fprintf(stderr," -o n          Choose TIFF output encoding from 1..n\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"               Then for each profile in chain:\n");
	fprintf(stderr,"   -i intent   p = perceptual, r = relative colorimetric,\n");
	fprintf(stderr,"               s = saturation, a = absolute colorimetric\n");
	fprintf(stderr,"   profile.icm Device, Link or Abstract profile.\n");
	fprintf(stderr,"\n");
	fprintf(stderr," infile.tif    Input TIFF file in appropriate color space\n");
	fprintf(stderr," outfile.tif   Output TIFF file\n");
	exit(1);
}

/* Conversion functions from direct binary 0..n^2-1 == 0.0 .. 1.0 range */
/* to ICC luo input range, and the reverse. */

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
	out[1] = (in[1] - 32768.0/65535.0)/256.0;
	out[2] = (in[2] - 32768.0/65535.0)/256.0;
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


/* Convert a TIFF Photometric tag to an ICC colorspace. */
/* return 0 if not possible or applicable. */
icColorSpaceSignature 
TiffPhotometric2ColorSpaceSignature(
void (**ocvt)(double *out, double *in),	/* Return write conversion function, NULL if none */
void (**icvt)(double *out, double *in),	/* Return read conversion function, NULL if none */
int *smsk,			/* Return signed handling mask, 0x0 if none */
int pmtc,			/* Input TIFF photometric */
int bps,			/* Input Bits per sample */
int spp				/* Input Samples per pixel */
) {
	if (icvt != NULL)
		*icvt = NULL;		/* Default return values */
	if (ocvt != NULL)
		*ocvt = NULL;		/* Default return values */
	if (smsk != NULL)
		*smsk = 0x0;

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
			switch(spp) {
				case 2:
					return icSig2colorData;
				case 3:
					return icSig3colorData;
				case 4:
					return icSig4colorData;
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
			*smsk = 0x6;				/* Tread a & b as signed */
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

/* Convert an ICC colorspace to the corresponding TIFF Inkset tag */
/* return 0xffff if not possible or applicable. */
int
ColorSpaceSignature2TiffInkset(
icColorSpaceSignature cspace,
char **inknames					/* Return ASCII inknames */
) {
	switch(cspace) {
		case icSigCmyData:
			if (inknames != NULL)
				*inknames = "cyan\000magenta\000yellow\000\000";
			return 0;				/* Not CMYK */
		case icSigCmykData:
			if (inknames != NULL)
				*inknames = NULL;	/* No inknames */
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
				/* This is a cheat and a hack. Should really make sure that icclink */
				/* transfers the right information from the destination */
				/* profile, and then copies it to the device profile, */
				/* allowing cctiff to read it. */
			if (inknames != NULL)
				*inknames = "cyan\000magenta\000yellow\000black\000orange\000green\000\000";
			return 0;				/* Not CMYK */

		case icSig7colorData:
		case icSigMch7Data:
			return 0xffff;

		case icSig8colorData:
		case icSigMch8Data:
				/* This is a cheat and a hack. Should really make sure that icclink */
				/* transfers the right information from the destination */
				/* profile, and then copies it to the device profile, */
				/* allowing cctiff to read it. */
			if (inknames != NULL)
				*inknames = "cyan\000magenta\000yellow\000black\000orange\000green\000lightcyan\000lightmagenta\000\000";
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
	sprintf(buf,"Unknown Tag %d",pmtc);
	return buf;
}

/* Callbacks used to initialise imdi */

/* Information needed from a profile */
struct _profinfo {
	char name[500];
	icmFile *fp;
	icc *c;
	icmHeader *h;
	icRenderingIntent intent;	/* Rendering intent chosen */
	icmLookupFunc func;			/* Type of function to use in lookup */
	icColorSpaceSignature ins, outs;	/* Colorspace of conversion */
	icColorSpaceSignature id, od;		/* Dimensions of conversion */
	icmLuAlgType alg;					/* Type of lookup algorithm used */
	icColorSpaceSignature natpcs;		/* Underlying natural PCS */

	icmLuBase *luo;				/* Base Lookup type object */
}; typedef struct _profinfo profinfo;

/* Context for imdi setup callbacks */
typedef struct {
	/* Overall parameters */
	int verb;			/* Non-zero if verbose */
	icColorSpaceSignature ins, outs;	/* Input/Output spaces */
	int id, od;			/* Input/Output dimensions */
	int isign_mask;		/* Input sign mask */
	int osign_mask;		/* Output sign mask */
	int icombine;		/* Non-zero if input curves are to be combined */
	int ocombine;		/* Non-zero if output curves are to be combined */
	void (*icvt)(double *out, double *in);	/* If non-NULL, Input format conversion */
	void (*ocvt)(double *out, double *in);	/* If non-NULL, Output format conversion */

	int nprofs;			/* Number of profiles in the chain */
	profinfo *profs;	/* Profile information */
} sucntx;

/* Input curve function */
double input_curve(
	void *cntx,
	int ch,
	double in_val
) {
	sucntx *rx = (sucntx *)cntx;
	int i;
	double vals[MAX_CHAN];

#ifdef NEVER
	if (rx->icombine)
		return in_val;

	if (rx->link) {

		for (i = 0; i < rx->id; i++)
			vals[i] = 0.0;
		vals[ch] = in_val;

		switch(rx->in.alg) {
		    case icmMonoFwdType: {
				icmLuMono *lu = (icmLuMono *)rx->in.luo;	/* Safe to coerce */
				lu->fwd_curve(lu, vals, vals);
				break;
			}
		    case icmMatrixFwdType: {
				icmLuMatrix *lu = (icmLuMatrix *)rx->in.luo;	/* Safe to coerce */
				lu->fwd_curve(lu, vals, vals);
				break;
			}
		    case icmLutType: {
				icmLuLut *lu = (icmLuLut *)rx->in.luo;	/* Safe to coerce */
				/* Since not PCS, in_abs and matrix cannot be valid, */
				/* so input curve on own is ok to use. */
				lu->input(lu, vals, vals);
				break;
			}
			default:
				error("Unexpected algorithm type in input curve");
		}
	} else {
		icmLuLut *lu = (icmLuLut *)rx->dev.luo;		/* Safe to coerce */
	
		for (i = 0; i < rx->id; i++)
			vals[i] = 0.0;
		vals[ch] = in_val;
	
		/* Since input not PCS, in_abs and matrix cannot be valid, */
		/* so input curve on own is ok to use. */
		lu->input(lu, vals, vals);
	
	}
#endif // NEVER
	return vals[ch];
}

/* Multi-dim table function */
void md_table(
void *cntx,
double *out_vals,
double *in_vals
) {
	sucntx *rx = (sucntx *)cntx;
	double pcsv[3];
	int i;
	
#ifdef NEVER
	if (rx->link) {
		double vals[MAX_CHAN];

		switch(rx->in.alg) {
		    case icmMonoFwdType: {
				icmLuMono *lu = (icmLuMono *)rx->in.luo;	/* Safe to coerce */
				if (rx->icombine) {
					lu->fwd_curve(lu, vals, in_vals);
					lu->fwd_map(lu, pcsv, vals);
				} else {
					lu->fwd_map(lu, pcsv, in_vals);
				}
				lu->fwd_abs(lu, pcsv, pcsv);
				break;
			}
		    case icmMatrixFwdType: {
				icmLuMatrix *lu = (icmLuMatrix *)rx->in.luo;	/* Safe to coerce */
				if (rx->icombine) {
					lu->fwd_curve(lu, vals, in_vals);
					lu->fwd_matrix(lu, pcsv, vals);
				} else {
					lu->fwd_matrix(lu, pcsv, in_vals);
				}
				lu->fwd_abs(lu, pcsv, pcsv);
				break;
			}
		    case icmLutType: {
				icmLuLut *lu = (icmLuLut *)rx->in.luo;	/* Safe to coerce */
				if (rx->icombine) {
					/* Since not PCS, in_abs and matrix cannot be valid, */
					/* so input curve on own is ok to use. */
					lu->input(lu, vals, in_vals);
					lu->clut(lu, pcsv, vals);
				} else {
					lu->clut(lu, pcsv, in_vals);
				}
				lu->output(lu, pcsv, pcsv);
				lu->out_abs(lu, pcsv, pcsv);
				break;
			}
			default:
				error("Unexpected algorithm type in clut lookup");
		}

		switch(rx->out.alg) {
		    case icmMonoBwdType: {
				icmLuMono *lu = (icmLuMono *)rx->out.luo;	/* Safe to coerce */
				lu->bwd_abs(lu, pcsv, pcsv);
				lu->bwd_map(lu, out_vals, pcsv);
				if (rx->ocombine) {
					lu->bwd_curve(lu, out_vals, out_vals);
				}
				break;
			}
		    case icmMatrixBwdType: {
				icmLuMatrix *lu = (icmLuMatrix *)rx->out.luo;	/* Safe to coerce */
				lu->bwd_abs(lu, pcsv, pcsv);
				lu->bwd_matrix(lu, out_vals, pcsv);
				if (rx->ocombine) {
					lu->bwd_curve(lu, out_vals, out_vals);
				}
				break;
			}
		    case icmLutType: {
				icmLuLut *lu = (icmLuLut *)rx->out.luo;	/* Safe to coerce */
				lu->in_abs(lu, pcsv, pcsv);
				lu->matrix(lu, pcsv, pcsv);
				lu->input(lu, pcsv, pcsv);
				lu->clut(lu, out_vals, pcsv);
				if (rx->ocombine) {
					lu->output(lu, out_vals, out_vals);
					/* Since not PCS, out_abs is never used */
				}
				break;
			}

			default:
				error("Unexpected algorithm type in clut lookup");
		}
	} else {
		icmLuLut *lu = (icmLuLut *)rx->dev.luo;		/* Safe to coerce */

		if (rx->icombine && rx->ocombine) {
			lu->lookup((icmLuBase *)lu, out_vals, in_vals);		/* Do everything here */
		} else {
			lu->clut(lu, out_vals, in_vals);
		}
	}
#endif // NEVER
}

/* Output curve function */
double output_curve(
void *cntx,
int ch,
double in_val
) {
	sucntx *rx = (sucntx *)cntx;
	int i;
	double vals[MAX_CHAN];
	
	if (rx->ocombine)
		return in_val;

#ifdef NEVER
	if (rx->link) {
		for (i = 0; i < rx->od; i++)
			vals[i] = 0.0;
		vals[ch] = in_val;
	
		switch(rx->out.alg) {
		    case icmMonoBwdType: {
				icmLuMono *lu = (icmLuMono *)rx->out.luo;	/* Safe to coerce */
				lu->bwd_curve(lu, vals, vals);
				break;
			}
		    case icmMatrixBwdType: {
				icmLuMatrix *lu = (icmLuMatrix *)rx->out.luo;	/* Safe to coerce */
				lu->bwd_curve(lu, vals, vals);
				break;
			}
		    case icmLutType: {
				icmLuLut *lu = (icmLuLut *)rx->out.luo;	/* Safe to coerce */
				lu->output(lu, vals, vals);
				/* Since not PCS, out_abs is never used */
				break;
			}
			default:
				error("Unexpected algorithm type in devop_devo()");
		}

	} else {
		icmLuLut *lu = (icmLuLut *)rx->dev.luo;		/* Safe to coerce */
	
		for (i = 0; i < rx->od; i++)
			vals[i] = 0.0;
		vals[ch] = in_val;
	
		/* Since output not PCS, out_abs cannot be valid, */
		lu->output(lu, vals, vals);
	
	}
#endif // NEVER
	return vals[ch];
}

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

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	char in_name[100];		/* Raster file name */
	char out_name[100];		/* Raster file name */
	icRenderingIntent next_intent;	/* Rendering intent for next profile */
	icColorSpaceSignature last_colorspace;		/* Next colorspace between conversions */
	int slow = 0;
	int check = 0;
	int ochoice = 0;		/* Output photometric choice */
	int i, rv = 0;

	/* TIFF file info */
	TIFF *rh = NULL, *wh = NULL;

	int x, y, width, height;					/* Common size of image */
	uint16 bitspersample;						/* Bits per sample */
	uint16 resunits;
	float resx, resy;
	uint16 pconfig;								/* Planar configuration */

	uint16 rsamplesperpixel, wsamplesperpixel;	/* Bits per sample */
	uint16 rphotometric, wphotometric;			/* Photometrics */

	tdata_t *inbuf, *outbuf, *checkbuf;

	/* IMDI */
	imdi *s = NULL;
	sucntx su;		/* Setup context */
	unsigned char *inp[MAX_CHAN];
	unsigned char *outp[MAX_CHAN];
	int clutres = 0;

	/* Error check */
	int mxerr = 0;
	double avgerr = 0.0;
	double avgcount = 0.0;

	if (argc < 2)
		usage("Too few arguments");

	su.verb = 0;
	su.icombine = 0;
	su.ocombine = 0;

	su.nprofs = 0;
	su.profs = NULL;
	next_intent = icmDefaultIntent;

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
		char mes[200];
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

			/* Slow, Precise */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				slow = 1;
			}

			/* Combine per channel curves */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				su.icombine = 1;
				su.ocombine = 1;
			}

			/* Check curves */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				check = 1;
			}

			/* Output photometric choice */
			else if (argv[fa][1] == 'o' || argv[fa][1] == 'O') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to -o flag");
				ochoice = atoi(na);
			}
			/* Next profile Intent */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
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
						sprintf(mes,"Unknown argument '%c' to -i flag",na[0]);
						usage(mes);
				}
			}

			/* Verbosity */
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				su.verb = 1;
			}

			else  {
				sprintf(mes,"Unknown flag '%c'",argv[fa][1]);
				usage(mes);
			}

		} else if (argv[fa][0] != '\000') {
			/* Get the next filename */
		
			if (su.nprofs == 0)
				su.profs = (profinfo *)malloc(sizeof(profinfo));
			else
				su.profs = (profinfo *)realloc(su.profs, (su.nprofs+1) * sizeof(profinfo));
			if (su.profs == NULL)
				error("Malloc failed in allocating space for profile info.");

			strcpy(su.profs[su.nprofs].name,argv[fa]);
			su.profs[su.nprofs].intent = next_intent;

			su.nprofs++;
			next_intent = icmDefaultIntent;
		} else {
			break;
		}
	}

	/* The last two "profiles" are actually the input and output TIFF filenames */
	/* Unwind them */

	if (su.nprofs < 3)
		usage("Not enough arguments to specify input and output TIFF files");

	strcpy(out_name,su.profs[--su.nprofs].name);
	strcpy(in_name,su.profs[--su.nprofs].name);

#ifdef NEVER
	/* Dump where we're at */
	printf("~1 There are %d profile in the chain\n",su.nprofs);
	for (i = 0; i < su.nprofs; i++) {
		printf("~1 Profile %d is '%s' intent 0x%x\n",i+1,su.profs[i].name,su.profs[i].intent);
	}
	printf("~1 Input TIFF is '%s'\n",in_name);
	printf("~1 Output TIFF is '%s'\n",out_name);
#endif

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

	TIFFGetField(rh, TIFFTAG_PHOTOMETRIC, &rphotometric);
	TIFFGetField(rh, TIFFTAG_SAMPLESPERPIXEL, &rsamplesperpixel);

	/* Figure out how to handle the input TIFF colorspace */
	if ((su.ins = TiffPhotometric2ColorSpaceSignature(NULL, &su.icvt, &su.isign_mask, rphotometric,
	                                     bitspersample, rsamplesperpixel)) == 0)
		error("Can't handle TIFF file photometric %s", Photometric2str(rphotometric));

	TIFFGetField(rh, TIFFTAG_PLANARCONFIG, &pconfig);
	if (pconfig != PLANARCONFIG_CONTIG)
		error ("TIFF Input file must be planar");

	TIFFGetField(rh, TIFFTAG_RESOLUTIONUNIT, &resunits);
	TIFFGetField(rh, TIFFTAG_XRESOLUTION, &resx);
	TIFFGetField(rh, TIFFTAG_YRESOLUTION, &resy);

	last_colorspace = su.ins;

	if (su.verb) {
		printf("TIFF file colorspace is %s\n",icm2str(icmColorSpaceSignature,su.ins));
	}

	/* - - - - - - - - - - - - - - - */
	/* Check and setup the chain of ICC profiles */

	/* For each profile in the chain, configure it to transform the color */
	/* appropriately */
	for (i = 0; i < su.nprofs; i++) {

		/* Open up the profile for reading */
		if ((su.profs[i].fp = new_icmFileStd_name(su.profs[i].name,"r")) == NULL)
			error ("Can't open profile '%s'",su.profs[i].name);
	
		if ((su.profs[i].c = new_icc()) == NULL)
			error ("Creation of ICC object for '%s' failed",su.profs[i].name);
	
		if ((rv = su.profs[i].c->read(su.profs[i].c, su.profs[i].fp, 0)) != 0)
			error ("%d, %s from '%s'",rv,su.profs[i].c->err,su.profs[i].name);
		su.profs[i].h = su.profs[i].c->header;

		/* Deal with different profile classes */
		switch (su.profs[i].h->deviceClass) {
    		case icSigAbstractClass:
    		case icSigLinkClass:
				su.profs[i].func = icmFwd;
				su.profs[i].intent = icmDefaultIntent;
				break;

    		case icSigColorSpaceClass:
				su.profs[i].intent = icmDefaultIntent;
				/* Fall through */

    		case icSigInputClass:
    		case icSigDisplayClass:
    		case icSigOutputClass:
				/* Note we don't handle an ambigious (both directions match) case. */
				/* We would need direction from the user to resolve this. */
				if (CSMatch(last_colorspace, su.profs[i].h->colorSpace))
					su.profs[i].func = icmFwd;
				else {
					su.profs[i].func = icmBwd;		/* PCS -> Device */
				}
				break;

			default:
				error("Can't handle deviceClass %s from file '%s'",
				     icm2str(icmProfileClassSignature,su.profs[i].h->deviceClass),
				     su.profs[i].c->err,su.profs[i].name);
		}

		/* Get a conversion object */
		if ((su.profs[i].luo = su.profs[i].c->get_luobj(su.profs[i].c, su.profs[i].func,
		                          su.profs[i].intent, icmSigDefaultData, icmLuOrdNorm)) == NULL)
			error ("%d, %s from '%s'",su.profs[i].c->errc, su.profs[i].c->err, su.profs[i].name);
	
		/* Get details of conversion */
		su.profs[i].luo->spaces(su.profs[i].luo, &su.profs[i].ins, &su.profs[i].id,
		         &su.profs[i].outs, &su.profs[i].od, &su.profs[i].alg, NULL, NULL, NULL);

		/* Get native PCS space */
		su.profs[i].luo->lutspaces(su.profs[i].luo, NULL, NULL, NULL, NULL, &su.profs[i].natpcs);

		if (su.verb) {
			icmFile *op;
			if ((op = new_icmFileStd_fp(stdout)) == NULL)
				error ("Can't open stdout");
			su.profs[i].h->dump(su.profs[i].h, op, 1);
			op->del(op);
		}

		/* Check that we can join to previous correctly */
		if (!CSMatch(last_colorspace, su.profs[i].ins))
			error("Last colorspace %s doesn't match input space %s of profile %s",
		      icm2str(icmColorSpaceSignature,last_colorspace),
				      icm2str(icmColorSpaceSignature,su.profs[i].h->colorSpace),
				      su.profs[i].name);
	}
	
	/* Setup input/output curve use.
	if (su.profs[0].natpcs == icSigXYZData)
		su.icombine = 1;			/* XYZ is to non-linear to be a benefit */
	if (su.profs[su.nprofs-1].natpcs == icSigXYZData)
		su.ocombine = 1;			/* XYZ is to non-linear to be a benefit */

	/* - - - - - - - - - - - - - - - */
	/* Open up the output file for writing */
	if ((wh = TIFFOpen(out_name, "w")) == NULL)
		error("Can\'t create TIFF file '%s'!",out_name);
	
	wsamplesperpixel = su.profs[su.nprofs-1].od;
	su.od = wsamplesperpixel;

	/* Lookup and decide what TIFF photometric suites the output colorspace */
	{
		int no_pmtc;					/* Number of possible photometrics */
		uint16 pmtc[10];				/* Photometrics of output file */
		if ((no_pmtc = ColorSpaceSignature2TiffPhotometric(pmtc,
		                                    last_colorspace)) == 0)
			error("TIFF file can't handle output colorspace '%s'!",
			      icm2str(icmColorSpaceSignature, last_colorspace));
	
		if (no_pmtc > 1) {		/* Choose photometric */
			if (su.verb) {
				printf("Possible photometrics for output colorspace %s are:\n",
				        icm2str(icmColorSpaceSignature,last_colorspace));
				for (i = 0; i < no_pmtc; i++)
					printf("%d: %s\n",i+1, Photometric2str(pmtc[i]));
			}
			if (ochoice < 1)
				ochoice = 1;
			else if (ochoice > no_pmtc)
				ochoice = no_pmtc;
			if (su.verb)
				printf("Using choice %d\n",ochoice);
		}
		wphotometric = pmtc[ochoice-1];
	}

	/* Lookup what we need to handle this. */
	if ((su.outs = TiffPhotometric2ColorSpaceSignature(&su.ocvt, NULL, &su.osign_mask, wphotometric,
	                                     bitspersample, wsamplesperpixel)) == 0)
		error("Can't handle TIFF file photometric %s", Photometric2str(wphotometric));

	/* Configure the output TIFF file appropriately */

	TIFFSetField(wh, TIFFTAG_PHOTOMETRIC, wphotometric);
	if (wphotometric == PHOTOMETRIC_SEPARATED) {
		int iset;
		char *inames;
		iset = ColorSpaceSignature2TiffInkset(su.outs, &inames);
		if (iset != 0xffff) {
			TIFFSetField(wh, TIFFTAG_INKSET, iset);
			if (inames != NULL) {
				TIFFSetField(wh, TIFFTAG_INKNAMES, inames);
			}
		}
	}

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
	TIFFSetField(wh, TIFFTAG_IMAGEDESCRIPTION, "Color corrected by Argyll");

	/* - - - - - - - - - - - - - - - */
	/* Setup the imdi */

	/* Setup the imdi resolution */
#ifdef NEVER
	if (only one profile in link)
		set res to res of first profile
	else
		select from table of link resolutions
#endif

	/* If one conversion, make the resolution the same */
	if (su.nprofs == 1 && su.profs[0].alg == icmLutType) {
		icmLut *lut;
		icmLuLut *luluo = (icmLuLut *)su.profs[0].luo;		/* Safe to coerce */
		luluo->get_info(luluo, &lut, NULL, NULL, NULL);	/* Get some details */
		clutres = lut->clutPoints;			/* Desired table resolution */
	}

	if (clutres <= 0) {

		/* Choose a "high quality" clut resolution. */
		/* ~~ should scan through all the conversions and */
		/* choose a level equal or greater than all conversions~~~~ */

		if ((clutres = dim_to_clutres(su.id, 2)) == 0)
			error ("Illegal number of input chanels");
	}
// ~~~9999

	if (!slow) {
		s = new_imdi(
			su.id,			/* Number of input dimensions */
			su.od,			/* Number of output dimensions */
							/* Input pixel representation */
			bitspersample == 8 ? pixint8 : pixint16,
							/* Output pixel representation */
			0x0,			/* Treat every channel as unsigned */
			bitspersample == 8 ? pixint8 : pixint16,
			0x0,			/* Treat every channel as unsigned */
			clutres,		/* Desired table resolution */
			input_curve,	/* Callback functions */
			md_table,
			output_curve,
			(void *)&su		/* Context to callbacks */
		);
	
		if (s == NULL)
			error("new_imdi failed");
	}

	/* - - - - - - - - - - - - - - - */
	/* Process colors to translate */
	/* (Should fix this to process a group of lines at a time ?) */

	inbuf  = _TIFFmalloc(TIFFScanlineSize(rh));
	outbuf = _TIFFmalloc(TIFFScanlineSize(wh));
	if (check)
		checkbuf = _TIFFmalloc(TIFFScanlineSize(wh));

	inp[0] = (unsigned char *)inbuf;
	outp[0] = (unsigned char *)outbuf;

	if (!slow) {		/* Fast */
		for (y = 0; y < height; y++) {

			/* Read in the next line */
			if (TIFFReadScanline(rh, inbuf, y, 0) < 0)
				error ("Failed to read TIFF line %d",y);

			/* Do fast conversion */
			s->interp(s, (void **)outp, (void **)inp, width);
			
			if (check) {
				/* Do floating point conversion */
				for (x = 0; x < width; x++) {
					int i;
					double in[MAX_CHAN], out[MAX_CHAN];
					
					if (bitspersample == 8)
						for (i = 0; i < su.id; i++)
							in[i] = ((unsigned char *)inbuf)[x * su.id + i]/255.0;
					else
						for (i = 0; i < su.id; i++)
							in[i] = ((unsigned short *)inbuf)[x * su.id + i]/65535.0;
					
					if (su.link) {
						if ((rv = su.in.luo->lookup(su.in.luo, out, in)) > 1)
							error ("%d, %s",su.dev.c->errc,su.dev.c->err);
						if ((rv = su.out.luo->lookup(su.out.luo, out, out)) > 1)
							error ("%d, %s",su.dev.c->errc,su.dev.c->err);
					} else {
						if ((rv = su.dev.luo->lookup(su.dev.luo, out, in)) > 1)
							error ("%d, %s",su.dev.c->errc,su.dev.c->err);
					}
	
					if (bitspersample == 8)
						for (i = 0; i < su.od; i++)
							((unsigned char *)checkbuf)[x * su.od + i] = (int)(out[i] * 255.0 + 0.5);
					else
						for (i = 0; i < su.od; i++)
							((unsigned short *)checkbuf)[x * su.od + i] = (int)(out[i] * 65535.0 + 0.5);
				}
				/* Compute the errors */
				for (x = 0; x < (width * su.od); x++) {
					int err;
					if (bitspersample == 8)
						err = ((unsigned char *)outbuf)[x] - ((unsigned char *)checkbuf)[x];
					else
						err = ((unsigned short *)outbuf)[x] - ((unsigned short *)checkbuf)[x];
					if (err < 0)
						err = -err;
					if (err > mxerr)
						mxerr = err;
					avgerr += (double)err;
					avgcount++;
				}
			}
				
			if (TIFFWriteScanline(wh, outbuf, y, 0) < 0)
				error ("Failed to write TIFF line %d",y);

		}

	} else {	/* Slow but precise */
		if (bitspersample == 8) {
			for (y = 0; y < height; y++) {

				/* Read in the next line */
				if (TIFFReadScanline(rh, inbuf, y, 0) < 0)
					error ("Failed to read TIFF line %d",y);

				/* Do floating point conversion */
				for (x = 0; x < width; x++) {
					int i;
					double in[MAX_CHAN], out[MAX_CHAN];
					
					for (i = 0; i < su.id; i++) {
						in[i] = ((unsigned char *)inbuf)[x * su.id + i]/255.0;
					}
					
					if (su.link) {
						if ((rv = su.in.luo->lookup(su.in.luo, out, in)) > 1)
							error ("%d, %s",su.dev.c->errc,su.dev.c->err);
						if ((rv = su.out.luo->lookup(su.out.luo, out, out)) > 1)
							error ("%d, %s",su.dev.c->errc,su.dev.c->err);
					} else {
						if ((rv = su.dev.luo->lookup(su.dev.luo, out, in)) > 1)
							error ("%d, %s",su.dev.c->errc,su.dev.c->err);
					}

					for (i = 0; i < su.od; i++) {
						((unsigned char *)outbuf)[x * su.od + i] = (int)(out[i] * 255.0 + 0.5);
					}
				}
				if (TIFFWriteScanline(wh, outbuf, y, 0) < 0)
					error ("Failed to write TIFF line %d",y);
			}
		} else if (bitspersample == 16) {
			for (y = 0; y < height; y++) {

				/* Read in the next line */
				if (TIFFReadScanline(rh, inbuf, y, 0) < 0)
					error ("Failed to read TIFF line %d",y);

				/* Do floating point conversion */
				for (x = 0; x < width; x++) {
					int i;
					double in[MAX_CHAN], out[MAX_CHAN];
					
					for (i = 0; i < su.id; i++) {
						in[i] = ((unsigned short *)inbuf)[x * su.id + i]/65535.0;
					}
					
					if (su.link) {
						if ((rv = su.in.luo->lookup(su.in.luo, out, in)) > 1)
							error ("%d, %s",su.dev.c->errc,su.dev.c->err);
						if ((rv = su.out.luo->lookup(su.out.luo, out, out)) > 1)
							error ("%d, %s",su.dev.c->errc,su.dev.c->err);
					} else {
						if ((rv = su.dev.luo->lookup(su.dev.luo, out, in)) > 1)
							error ("%d, %s",su.dev.c->errc,su.dev.c->err);
					}

					for (i = 0; i < su.od; i++) {
						((unsigned short *)outbuf)[x * su.od + i] = (int)(out[i] * 65535.0 + 0.5);
					}
				}
				if (TIFFWriteScanline(wh, outbuf, y, 0) < 0)
					error ("Failed to write TIFF line %d",y);
			}
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

	/* Done with lookup object */
	if (s != NULL)
		s->done(s);

	if (su.link) {
		su.in.luo->del(su.in.luo);
		su.in.c->del(su.in.c);
		su.in.fp->del(su.in.fp);
		su.out.luo->del(su.out.luo);
		su.out.c->del(su.out.c);
		su.out.fp->del(su.out.fp);
	} else {
		su.dev.luo->del(su.dev.luo);
		su.dev.c->del(su.dev.c);
		su.dev.fp->del(su.dev.fp);
	}

	_TIFFfree(inbuf);
	_TIFFfree(outbuf);
	if (check)
		_TIFFfree(checkbuf);

	TIFFClose(rh);		/* Close Input file */
	TIFFClose(wh);		/* Close Output file */

	return 0;
}


/* Basic printf type error() and warning() routines */

void
error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"cctiff: Error - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	exit (-1);
}

void
warning(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"cctiff: Warning - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}
