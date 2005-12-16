
/* 
 * Create a visualisation of the gamut of a TIFF file.
 *
 * Author:  Graeme W. Gill
 * Date:    00/9/20
 * Version: 1.00
 *
 * Copyright 2000, 2002 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the LICENCE.TXT file for licencing details.
 */

/*
 * TTBD:
 *
 *		Should record abs/relative intent, PCS space used together
 *		with viewing conditions, so that a mismatch within
 *		icclink can be detected or allowed for.
 *
 *       Need to cope with profile not having black point.
 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "tiffio.h"
#include "numlib.h"
#include "icc.h"
#include "gamut.h"
#include "xicc.h"

#ifdef NEVER
void error(char *fmt, ...), warning(char *fmt, ...);
#endif

void usage(void) {
	int i;
	fprintf(stderr,"Create VRML image of the gamut surface of a TIFF, V2.00\n");
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	fprintf(stderr,"usage: tiffgamut [-v level] profile.icm infile.tif\n");
	fprintf(stderr," -v            Verbose\n");
	fprintf(stderr," -i intent     p = perceptual, r = relative colorimetric,\n");
	fprintf(stderr,"               s = saturation, a = absolute, j = Appearance %s\n",icxcam_description(cam_default));
	fprintf(stderr," -o order      n = normal (priority: lut > matrix > monochrome)\n");
	fprintf(stderr,"               r = reverse (priority: monochrome > matrix > lut)\n");
	fprintf(stderr," -c viewcond   set viewing conditions for %s,\n",icxcam_description(cam_default));
	fprintf(stderr,"               either an enumerated choice, or a parameter:\n");
	for (i = 0; ; i++) {
		icxViewCond vc;
		if (xicc_enum_viewcond(NULL, &vc, i, 1))
			break;

		fprintf(stderr,"               %d: %s\n",i,vc.desc);
	}
	fprintf(stderr,"         s:surround    a = average, m = dim, d = dark,\n");
	fprintf(stderr,"                       c = transparency (default average)\n");
	fprintf(stderr,"         w:X:Y:Z       Adapted white point as XYZ (default media white)\n");
	fprintf(stderr,"         w:x:y         Adapted white point as x, y\n");
	fprintf(stderr,"         a:adaptation  Adaptation luminance in cd.m^2 (default 50.0)\n");
	fprintf(stderr,"         b:background  Background %% of image luminance (default 20)\n");
	fprintf(stderr,"         f:flare       Flare light %% of image luminance (default 1)\n");
	fprintf(stderr,"         f:X:Y:Z       Flare color as XYZ (default media white)\n");
	fprintf(stderr,"         f:x:y         Flare color as x, y\n");
	fprintf(stderr," -d sres       Surface resolution details 1.0 - 50.0\n");
	fprintf(stderr," -w            emit VRML .wrl file as well as CGATS .gam file\n");
	fprintf(stderr," -n            Don't add VRML axes or white/black point\n");
	exit(1);
}

#define GAMRES 15.0		/* Default surface resolution */

/* Convert an ICC colorspace to the corresponding TIFF Photometric tag */
/* return 0xffff if not possible. */

int
ColorSpaceSignature2TiffPhotometric(
icColorSpaceSignature cspace
) {
	switch(cspace) {
		case icSigGrayData:
			return PHOTOMETRIC_MINISBLACK;
		case icSigRgbData:
			return PHOTOMETRIC_RGB;
		case icSigCmykData:
			return PHOTOMETRIC_SEPARATED;
		case icSigYCbCrData:
			return PHOTOMETRIC_YCBCR;
		case icSigLabData:
			return PHOTOMETRIC_CIELAB;

		case icSigXYZData:
		case icSigLuvData:
		case icSigYxyData:
		case icSigHsvData:
		case icSigHlsData:
		case icSigCmyData:
		case icSig2colorData:
		case icSig3colorData:
		case icSig4colorData:
		case icSig5colorData:
		case icSigMch5Data:
		case icSig6colorData:
		case icSigMch6Data:
		case icSig7colorData:
		case icSigMch7Data:
		case icSig8colorData:
		case icSigMch8Data:
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
			return "CMYK";
		case PHOTOMETRIC_YCBCR:
			return "YCbCr";
		case PHOTOMETRIC_CIELAB:
			return "CIELab";
		case PHOTOMETRIC_LOGL:
			return "CIELog2L";
		case PHOTOMETRIC_LOGLUV:
			return "CIELog2Luv";
	}
	sprintf(buf,"Unknown Tag %d",pmtc);
	return buf;
}

int
main(int argc, char *argv[]) {
	int fa,nfa;					/* argument we're looking at */
	char prof_name[100];		/* Icc profile name */
	char in_name[100];			/* TIFF input file */
	char *xl, out_name[100];	/* VRML output file */

	icmFile *p_fp;
	icc *icco;
	xicc *xicco;
	icxViewCond vc;				/* Viewing Condition for CIECAM */
	int vc_e = -1;				/* Enumerated viewing condition */
	int vc_s = -1;				/* Surround override */
	double vc_wXYZ[3] = {-1.0, -1.0, -1.0};	/* Adapted white override in XYZ */
	double vc_wxy[2] = {-1.0, -1.0};		/* Adapted white override in x,y */
	double vc_a = -1.0;			/* Adapted luminance */
	double vc_b = -1.0;			/* Background % overid */
	double vc_f = -1.0;			/* Flare % overid */
	double vc_fXYZ[3] = {-1.0, -1.0, -1.0};	/* Flare color override in XYZ */
	double vc_fxy[2] = {-1.0, -1.0};		/* Flare color override in x,y */
	int verb = 0;
	int vrml = 0;
	int doaxes = 1;
	int rv = 0;

	TIFF *rh = NULL;
	int x, y, width, height;					/* Size of image */
	uint16 samplesperpixel, bitspersample;
	uint16 pconfig, photometric, pmtc;
	uint16 resunits;
	float resx, resy;
	tdata_t *inbuf;

	double gamres = GAMRES;				/* Surface resolution */
	gamut *gam;

	icxLuBase *luo;						/* Generic lookup object */
	icColorSpaceSignature ins, outs;	/* Type of input and output spaces */
	int inn, outn;						/* Number of components */
	icmLuAlgType alg;					/* Type of lookup algorithm */

	/* Lookup parameters */
	icmLookupFunc     func   = icmFwd;				/* Must be */
	icRenderingIntent intent = icAbsoluteColorimetric;	/* Default */
	icColorSpaceSignature pcsor = icSigLabData;		/* Default */
	icmLookupOrder    order  = icmLuOrdNorm;		/* Default */
	
	if (argc < 2)
		usage();

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
				usage();

			/* Verbosity */
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;
			}

			/* Intent */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL) usage();
    			switch (na[0]) {
					case 'p':
					case 'P':
						intent = icPerceptual;
						break;
					case 'r':
					case 'R':
						intent = icRelativeColorimetric;
						break;
					case 's':
					case 'S':
						intent = icSaturation;
						break;
					case 'a':
					case 'A':
						intent = icAbsoluteColorimetric;
						break;
					case 'j':
					case 'J':
						intent = icxAppearance;
						pcsor = icxSigJabData;
						break;
					default:
						usage();
				}
			}

			/* Search order */
			else if (argv[fa][1] == 'o' || argv[fa][1] == 'O') {
				fa = nfa;
				if (na == NULL) usage();
    			switch (na[0]) {
					case 'n':
					case 'N':
						order = icmLuOrdNorm;
						break;
					case 'r':
					case 'R':
						order = icmLuOrdRev;
						break;
					default:
						usage();
				}
			}

			/* Viewing conditions */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage();
				if (na[0] >= '0' && na[0] <= '9') {
					vc_e = atoi(na);
				} else if (na[0] == 's' || na[0] == 'S') {
					if (na[1] != ':')
						usage();
					if (na[2] == 'a' || na[2] == 'A') {
						vc_s = vc_average;
					} else if (na[2] == 'm' || na[2] == 'M') {
						vc_s = vc_dim;
					} else if (na[2] == 'd' || na[2] == 'D') {
						vc_s = vc_dark;
					} else if (na[2] == 'c' || na[2] == 'C') {
						vc_s = vc_cut_sheet;
					} else
						usage();
				} else if (na[0] == 'w' || na[0] == 'W') {
					double x, y, z;
					if (sscanf(na+1,":%lf:%lf:%lf",&x,&y,&z) == 3) {
						vc_wXYZ[0] = x; vc_wXYZ[1] = y; vc_wXYZ[2] = z;
					} else if (sscanf(na+1,":%lf:%lf",&x,&y) == 2) {
						vc_wxy[0] = x; vc_wxy[1] = y;
					} else
						usage();
				} else if (na[0] == 'a' || na[0] == 'A') {
					if (na[1] != ':')
						usage();
					vc_a = atof(na+2);
				} else if (na[0] == 'b' || na[0] == 'B') {
					if (na[1] != ':')
						usage();
					vc_b = atof(na+2);
				} else if (na[0] == 'f' || na[0] == 'F') {
					double x, y, z;
					if (sscanf(na+1,":%lf:%lf:%lf",&x,&y,&z) == 3) {
						vc_fXYZ[0] = x; vc_fXYZ[1] = y; vc_fXYZ[2] = z;
					} else if (sscanf(na+1,":%lf:%lf",&x,&y) == 2) {
						vc_fxy[0] = x; vc_fxy[1] = y;
					} else if (sscanf(na+1,":%lf",&x) == 1) {
						vc_f = x;
					} else
						usage();
				} else
					usage();
			}

			/* VRML output */
			else if (argv[fa][1] == 'w' || argv[fa][1] == 'W') {
				vrml = 1;
			}
			/* No axis output */
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N') {
				doaxes = 0;
			}
			/* Surface Detail */
			else if (argv[fa][1] == 'd' || argv[fa][1] == 'D') {
				fa = nfa;
				if (na == NULL) usage();
				gamres = atof(na);
			}

			else 
				usage();
		} else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(prof_name,argv[fa++]);

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(in_name,argv[fa++]);

	/* - - - - - - - - - - - - - - - - */
	/* Open up the profile for reading */
	if ((p_fp = new_icmFileStd_name(prof_name,"r")) == NULL)
		error ("Can't open file '%s'",prof_name);

	if ((icco = new_icc()) == NULL)
		error ("Creation of ICC object failed");

	if ((rv = icco->read(icco,p_fp,0)) != 0)
		error ("%d, %s",rv,icco->err);

	if (verb) {
		icmFile *op;
		if ((op = new_icmFileStd_fp(stdout)) == NULL)
			error ("Can't open stdout");
		icco->header->dump(icco->header, op, 1);
		op->del(op);
	}

	/* Check that the profile is appropriate */
	if (icco->header->deviceClass != icSigInputClass
	 && icco->header->deviceClass != icSigDisplayClass
	 && icco->header->deviceClass != icSigOutputClass)
		error("Profile isn't a device profile");

	/* Wrap with an expanded icc */
	if ((xicco = new_xicc(icco)) == NULL)
		error ("Creation of xicc failed");

	/* Setup the viewing conditions */
	if (xicc_enum_viewcond(xicco, &vc, -1, 0))
		error ("%d, %s",xicco->errc, xicco->err);

	if (vc_e >= 0)
		if (xicc_enum_viewcond(xicco, &vc, vc_e, 0))
			error ("%d, %s",xicco->errc, xicco->err);
	if (vc_s >= 0)
		vc.Ev = vc_s;
	if (vc_wXYZ[1] > 0.0) {
		/* Normalise it to current media white */
		vc.Wxyz[0] = vc_wXYZ[0]/vc_wXYZ[1] * vc.Wxyz[1];
		vc.Wxyz[2] = vc_wXYZ[2]/vc_wXYZ[1] * vc.Wxyz[1];
	} 
	if (vc_wxy[0] >= 0.0) {
		double x = vc_wxy[0];
		double y = vc_wxy[1];	/* If Y == 1.0, then X+Y+Z = 1/y */
		double z = 1.0 - x - y;
		vc.Wxyz[0] = x/y * vc.Wxyz[1];
		vc.Wxyz[2] = z/y * vc.Wxyz[1];
	}
	if (vc_a >= 0.0)
		vc.La = vc_a;
	if (vc_b >= 0.0)
		vc.Yb = vc_b/100.0;
	if (vc_f >= 0.0)
		vc.Yf = vc_f/100.0;
	if (vc_fXYZ[1] > 0.0) {
		/* Normalise it to current media white */
		vc.Fxyz[0] = vc_fXYZ[0]/vc_fXYZ[1] * vc.Fxyz[1];
		vc.Fxyz[2] = vc_fXYZ[2]/vc_fXYZ[1] * vc.Fxyz[1];
	}
	if (vc_fxy[0] >= 0.0) {
		double x = vc_fxy[0];
		double y = vc_fxy[1];	/* If Y == 1.0, then X+Y+Z = 1/y */
		double z = 1.0 - x - y;
		vc.Fxyz[0] = x/y * vc.Fxyz[1];
		vc.Fxyz[2] = z/y * vc.Fxyz[1];
	}

	/* Get a expanded color conversion object */
	if ((luo = xicco->get_luobj(xicco, 0
	   | ICX_CLIP_NEAREST
	           , func, intent, pcsor, order, &vc, NULL)) == NULL)
		error ("%d, %s",xicco->errc, xicco->err);

	luo->spaces(luo, &ins, &inn, &outs, &outn, &alg, NULL, NULL, NULL);

	/* - - - - - - - - - - - - - - - */
	/* Open up input tiff file ready for reading */
	/* Got arguments, so setup to process the file */
	if ((rh = TIFFOpen(in_name, "r")) == NULL)
		error("error opening read file '%s'",in_name);

	TIFFGetField(rh, TIFFTAG_IMAGEWIDTH,  &width);
	TIFFGetField(rh, TIFFTAG_IMAGELENGTH, &height);

	TIFFGetField(rh, TIFFTAG_BITSPERSAMPLE, &bitspersample);
	if (bitspersample != 8)
		error("TIFF Input file must be 8 bit/channel");

	TIFFGetField(rh, TIFFTAG_PHOTOMETRIC, &photometric);
	if  ((pmtc = ColorSpaceSignature2TiffPhotometric(ins)) == 0xffff)
		error("ICC  input colorspace '%s' can't be handled by a TIFF file!",
		      icm2str(icmColorSpaceSignature, ins));
	if (pmtc != photometric)
		error("ICC  input colorspace '%s' doesn't match TIFF photometric '%s'!",
		      icm2str(icmColorSpaceSignature,ins), Photometric2str(photometric));

	TIFFGetField(rh, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);
	if (inn != samplesperpixel)
		error ("TIFF Input file has %d input chanels mismatched to colorspace '%s'",
		       samplesperpixel, icm2str(icmColorSpaceSignature, ins));

	TIFFGetField(rh, TIFFTAG_PLANARCONFIG, &pconfig);
	if (pconfig != PLANARCONFIG_CONTIG)
		error ("TIFF Input file must be planar");

	TIFFGetField(rh, TIFFTAG_RESOLUTIONUNIT, &resunits);
	TIFFGetField(rh, TIFFTAG_XRESOLUTION, &resx);
	TIFFGetField(rh, TIFFTAG_YRESOLUTION, &resy);

	/* - - - - - - - - - - - - - - - */
	/* Creat a gamut surface */
	gam = new_gamut(gamres);

	/* - - - - - - - - - - - - - - - */
	/* Process colors to translate */
	/* (Should fix this to process a group of lines at a time ?) */

	inbuf  = _TIFFmalloc(TIFFScanlineSize(rh));

	for (y = 0; y < height; y++) {

		/* Read in the next line */
		if (TIFFReadScanline(rh, inbuf, y, 0) < 0)
			error ("Failed to read TIFF line %d",y);

		/* Do floating point conversion */
		for (x = 0; x < width; x++) {
			int i;
			double in[MAX_CHAN], out[MAX_CHAN];
			
			for (i = 0; i < inn; i++) {
				in[i] = ((unsigned char *)inbuf)[x * inn + i]/255.0;
			}
			
			if ((rv = luo->lookup(luo, out, in)) > 1)
				error ("%d, %s",icco->errc,icco->err);
			
			if (outs == icSigXYZData)	/* Convert to Lab */
				icmXYZ2Lab(&icco->header->illuminant, out, out);

			gam->expand(gam, out);
		}
	}

	/* Get White and Black points from the profile, and set them in the gamut */
	{
		double wp[3], bp[3];

		luo->rel_wh_bk_points(luo, wp, bp);
		gam->setwb(gam, wp, bp);
	}

	/* Done with lookup object */
	luo->del(luo);
	xicco->del(xicco);		/* Expansion wrapper */
	icco->del(icco);		/* Icc */
	p_fp->del(p_fp);

	TIFFClose(rh);		/* Close Input file */

	/* Create the VRML file */
	strcpy(out_name, in_name);
	if ((xl = strrchr(out_name, '.')) == NULL)	/* Figure where extention is */
		xl = out_name + strlen(out_name);

	strcpy(xl,".gam");
	if (gam->write_gam(gam,out_name))
		error ("write gamut failed on '%s'",out_name);

	if (vrml) {
		strcpy(xl,".wrl");
		if (gam->write_vrml(gam,out_name, doaxes))
			error ("write vrml failed on '%s'",out_name);
	}

	gam->del(gam);

	return 0;
}


#ifdef NEVER
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

#endif /* NEVER */
