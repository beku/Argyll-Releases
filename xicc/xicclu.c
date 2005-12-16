
/* 
 * xicc lookup/test utility
 *
 * This program is the analog of icclu, but allows reverse lookup
 * of transforms by making use of xicc interpolation code.
 * (Based on the old xfmlu.c)
 *
 * Author:  Graeme W. Gill
 * Date:    8/7/00
 * Version: 1.00
 *
 * Copyright 1999, 2000 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the LICENCE.TXT file for licencing details.
 */

/* TTBD:
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "numlib.h"
#include "plot.h"
#include "xicc.h"

#define USE_NEARCLIP		/* Our usual expectation */
#define XRES 128			/* Plotting resolution */

void usage(char *diag) {
	int i;
	fprintf(stderr,"Translate colors through an xicc, V1.00\n");
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	fprintf(stderr,"usage: xicclu [-options] profile\n");
	if (diag != NULL)
		fprintf(stderr,"Diagnostic: %s\n",diag);
	fprintf(stderr," -v            Verbose\n");
	fprintf(stderr," -g            Plot neutral axis instead of looking colors up.\n");
	fprintf(stderr," -f function   f = forward, b = backwards, g = gamut, p = preview\n");
	fprintf(stderr,"               if = inverted forward, ib = inverted backwards\n");
	fprintf(stderr," -i intent     p = perceptual, r = relative colorimetric,\n");
	fprintf(stderr,"               s = saturation, a = absolute, j = Appearance %s\n",icxcam_description(cam_default));
	fprintf(stderr,"               k = Absolute Appearance %s\n",icxcam_description(cam_default));
	fprintf(stderr," -o order      n = normal (priority: lut > matrix > monochrome)\n");
	fprintf(stderr,"               r = reverse (priority: monochrome > matrix > lut)\n");
	fprintf(stderr," -p oride      x = XYZ_PCS, l = Lab_PCS, y = Yxy, j = Jab, J = JCh\n");
	fprintf(stderr," -s scale      Scale device range 0.0 - scale rather than 0.0 - 1.0\n");
	fprintf(stderr," -k [zhxrlv]  Black generation: z = zero K,\n");
	fprintf(stderr,"               h = 0.5 K, x = max K, r = ramp K (def)\n");
	fprintf(stderr,"               l = extra PCS input is portion of K locus\n");
	fprintf(stderr,"               v = extra PCS input is K target value\n");
	fprintf(stderr," -k p stle stpo enpo enle shape\n");
	fprintf(stderr,"               stle: K level at White 0.0 - 1.0\n");
	fprintf(stderr,"               stpo: start point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"               enpo: End point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"               enle: K level at Black 0.0 - 1.0\n");
	fprintf(stderr,"               shape: 1.0 = straight, 0.0-1.0 concave, 1.0-2.0 convex\n");
	fprintf(stderr," -l tlimit     set total ink limit, 0 - 400%%\n");
	fprintf(stderr," -L klimit     set black ink limit, 0 - 100%%\n");
	fprintf(stderr," -a            show actual target values if clipped\n");
	fprintf(stderr," -m            merge output processing into clut\n");
	fprintf(stderr," -b            use CAM Jab for clipping\n");
	fprintf(stderr," -S            Use internal optimised separation for inverse 4d\n");

	fprintf(stderr," -c viewcond   set viewing conditions for CIECAM97s,\n");
	fprintf(stderr,"               either an enumerated choice, or a parameter:\n");
	for (i = 0; ; i++) {
		icxViewCond vc;
		if (xicc_enum_viewcond(NULL, &vc, i, 1))
			break;

		fprintf(stderr,"               %d: %s\n",i,vc.desc);
	}

	fprintf(stderr,"         s:surround    a = average, m = dim, d = dark,\n");
	fprintf(stderr,"                       c = transparency (default average)\n");
	fprintf(stderr,"         w:X:Y:Z       Adapted white point as XYZ (default media white, Abs: D50)\n");
	fprintf(stderr,"         w:x:y         Adapted white point as x, y\n");
	fprintf(stderr,"         a:adaptation  Adaptation luminance in cd.m^2 (default 50.0)\n");
	fprintf(stderr,"         b:background  Background %% of image luminance (default 20)\n");
	fprintf(stderr,"         f:flare       Flare light %% of image luminance (default 1)\n");
	fprintf(stderr,"         f:X:Y:Z       Flare color as XYZ (default media white, Abs: D50)\n");
	fprintf(stderr,"         f:x:y         Flare color as x, y\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"    The colors to be translated should be fed into standard in,\n");
	fprintf(stderr,"    one input color per line, white space separated.\n");
	fprintf(stderr,"    A line starting with a # will be ignored.\n");
	fprintf(stderr,"    A line not starting with a number will terminate the program.\n");
	exit(1);
}


int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	char prof_name[100];
	icmFile *fp;
	icc *icco;
	xicc *xicco;
	int doplot = 0;				/* Do grey axis plot */
	icxInk ink;					/* Ink parameters */
	int tlimit = -1;			/* Total ink limit as a % */
	int klimit = -1;			/* Black ink limit as a % */
	int intsep = 0;
	icxViewCond vc;				/* Viewing Condition for CIECAM97s */
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
	int actual = 0;
	int dual = 0;
	int merge = 0;
	int camclip = 0;
	int repYxy = 0;			/* Report Yxy */
	int repJCh = 0;			/* Report JCh */
	double scale = 0.0;		/* Device value scale factor */
	int rv = 0;
	char buf[200];
	double oin[MAX_CHAN], in[MAX_CHAN], out[MAX_CHAN];

	icxLuBase *luo, *aluo = NULL;
	icColorSpaceSignature ins, outs;	/* Type of input and output spaces */
	int inn, outn;						/* Number of components */
	icmLuAlgType alg;					/* Type of lookup algorithm */

	/* Lookup parameters */
	icmLookupFunc     func   = icmFwd;				/* Default */
	icRenderingIntent intent = icmDefaultIntent;	/* Default */
	icColorSpaceSignature pcsor = icmSigDefaultData;	/* Default */
	icmLookupOrder    order  = icmLuOrdNorm;		/* Default */
	int               inking = 3;					/* Default is ramp */
	double Kstle, Kstpo, Kenle, Kenpo, Kshap;		/* K curve params */
	int               invert = 0;
	
	error_program = argv[0];

	if (argc < 2)
		usage("Too few arguments");

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
				usage("Requested usage");

			/* Verbosity */
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;
			}
			/* Grey axis plot */
			else if (argv[fa][1] == 'g' || argv[fa][1] == 'G') {
				doplot = 1;
			}
			/* Actual target values */
			else if (argv[fa][1] == 'a' || argv[fa][1] == 'A') {
				actual = 1;
			}
			/* Merge output */
			else if (argv[fa][1] == 'm' || argv[fa][1] == 'M') {
				merge = 1;
			}
			/* Use CAM Jab for clipping on reverse lookup */
			else if (argv[fa][1] == 'b' || argv[fa][1] == 'B') {
				camclip = 1;
			}
			/* Use optimised internal separation */
			else if (argv[fa][1] == 'S') {
				intsep = 1;
			}
			/* Device scale */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -s");
				scale = atof(na);
				if (scale <= 0.0) usage("Illegal scale value");
			}
			/* function */
			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F') {
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -f");
    			switch (na[0]) {
					case 'f':
					case 'F':
						func = icmFwd;
						break;
					case 'b':
					case 'B':
						func = icmBwd;
						break;
					case 'g':
					case 'G':
						func = icmGamut;
						break;
					case 'p':
					case 'P':
						func = icmPreview;
						break;
					case 'i':
					case 'I':
						invert = 1;
						if (na[1] == 'f' || na[1] == 'F')
							func = icmFwd;
						else if (na[1] == 'b' || na[1] == 'B')
							func = icmBwd;
						else
							usage("Unknown parameter after flag -fi");
						break;
					default:
						usage("Unknown parameter after flag -f");
				}
			}

			/* Intent */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -i");
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
						intent = icxAppearance;
						repJCh = 0;
						break;
					case 'J':
						intent = icxAppearance;
						repJCh = 1;
						break;
					case 'k':
					case 'K':
						intent = icxAbsAppearance;
						break;
					default:
						usage("Unknown parameter after flag -i");
				}
			}

			/* PCS override */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -i");
    			switch (na[0]) {
					case 'x':
					case 'X':
						pcsor = icSigXYZData;
						repYxy = 0;
						break;
					case 'l':
					case 'L':
						pcsor = icSigLabData;
						repYxy = 0;
						break;
					case 'y':
					case 'Y':
						pcsor = icSigXYZData;
						repYxy = 1;
						break;
					case 'j':
						pcsor = icxSigJabData;
						repYxy = 0;
						repJCh = 0;
						break;
					case 'J':
						pcsor = icxSigJabData;
						repYxy = 0;
						repJCh = 1;
						break;
					default:
						usage("Unknown parameter after flag -i");
				}
			}

			/* Search order */
			else if (argv[fa][1] == 'o' || argv[fa][1] == 'O') {
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -o");
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
						usage("Unknown parameter after flag -o");
				}
			}

			/* Inking rule */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -k");
    			switch (na[0]) {
					case 'z':
					case 'Z':
						inking = 0;		/* Use minimum k */
						break;
					case 'h':
					case 'H':
						inking = 1;		/* Use 0.5 k */
						break;
					case 'x':
					case 'X':
						inking = 2;		/* Use maximum k */
						break;
					case 'r':
					case 'R':
						inking = 3;		/* Use ramping K */
						break;
					case 'l':
					case 'L':
						inking = 4;		/* Extra param is locus */
						break;
					case 'v':
					case 'V':
						inking = 5;		/* Extra param is K target */
						break;
					case 'p':
					case 'P':
						inking = 6;		/* Use curve parameter */
						++fa;
						if (fa >= argc) usage(NULL);
						Kstle = atof(argv[fa]);

						++fa;
						if (fa >= argc) usage(NULL);
						Kstpo = atof(argv[fa]);

						++fa;
						if (fa >= argc || argv[fa][0] == '-') usage(NULL);
						Kenpo = atof(argv[fa]);

						++fa;
						if (fa >= argc) usage(NULL);
						Kenle = atof(argv[fa]);

						++fa;
						if (fa >= argc || argv[fa][0] == '-') usage(NULL);
						Kshap = atof(argv[fa]);
						break;
					default:
						usage("Unknown parameter after flag -k");
				}
			}

			else if (argv[fa][1] == 'l') {
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -l");
				tlimit = atoi(na);
			}

			else if (argv[fa][1] == 'L') {
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -L");
				klimit = atoi(na);
			}

			/* Viewing conditions */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -c");
				if (na[0] >= '0' && na[0] <= '9') {
					vc_e = atoi(na);
				} else if (na[0] == 's' || na[0] == 'S') {
					if (na[1] != ':')
						usage("Unrecognised parameters after -cs");
					if (na[2] == 'a' || na[2] == 'A') {
						vc_s = vc_average;
					} else if (na[2] == 'm' || na[2] == 'M') {
						vc_s = vc_dim;
					} else if (na[2] == 'd' || na[2] == 'D') {
						vc_s = vc_dark;
					} else if (na[2] == 'c' || na[2] == 'C') {
						vc_s = vc_cut_sheet;
					} else
						usage("Unrecognised parameters after -cs:");
				} else if (na[0] == 'w' || na[0] == 'W') {
					double x, y, z;
					if (sscanf(na+1,":%lf:%lf:%lf",&x,&y,&z) == 3) {
						vc_wXYZ[0] = x; vc_wXYZ[1] = y; vc_wXYZ[2] = z;
					} else if (sscanf(na+1,":%lf:%lf",&x,&y) == 2) {
						vc_wxy[0] = x; vc_wxy[1] = y;
					} else
						usage("Unrecognised parameters after -cw");
				} else if (na[0] == 'a' || na[0] == 'A') {
					if (na[1] != ':')
						usage("Unrecognised parameters after -ca");
					vc_a = atof(na+2);
				} else if (na[0] == 'b' || na[0] == 'B') {
					if (na[1] != ':')
						usage("Unrecognised parameters after -cb");
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
						usage("Unrecognised parameters after -cf");
				} else
					usage("Unrecognised parameters after -c");
			}

			else 
				usage("Unknown flag");
		} else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage("Expecting profile file name");
	strcpy(prof_name,argv[fa]);

	if (doplot) {

		/* Force PCS to be Lab or Jab */
		repJCh = 0;
		if (pcsor != icxSigJabData)
			pcsor = icSigLabData;

		if ((invert == 0 && func != icmBwd)
		 || (invert != 0 && func != icmFwd))
			error("Must use -fb or -fif for grey axis plot");
	}

	/* Open up the profile for reading */
	if ((fp = new_icmFileStd_name(prof_name,"r")) == NULL)
		error ("Can't open file '%s'",prof_name);

	if ((icco = new_icc()) == NULL)
		error ("Creation of ICC object failed");

	if ((rv = icco->read(icco,fp,0)) != 0)
		error ("%d, %s",rv,icco->err);

	if (doplot) {
		if (icco->header->deviceClass != icSigInputClass
		 && icco->header->deviceClass != icSigDisplayClass
		 && icco->header->deviceClass != icSigOutputClass)
			error("Profile must be a device profile to plot neutral axis");
	}

	if (verb) {
		icmFile *op;
		if ((op = new_icmFileStd_fp(stdout)) == NULL)
			error ("Can't open stdout");
		icco->header->dump(icco->header, op, 1);
		op->del(op);
	}

	/* Wrap with an expanded icc */
	if ((xicco = new_xicc(icco)) == NULL)
		error ("Creation of xicc failed");

	/* Setup Inking and limit rules */
	if (tlimit >= 0)
		ink.tlimit = tlimit/100.0;	/* Set a total ink limit */
	else
		ink.tlimit = -1.0;			/* Don't use a limit */

	if (klimit >= 0)
		ink.klimit = klimit/100.0;	/* Set a black ink limit */
	else
		ink.klimit = -1.0;			/* Don't use a limit */

	if (inking == 0) {			/* Use minimum */
		ink.k_rule = icxKluma5;
		ink.c.Kstle = 0.0;
		ink.c.Kstpo = 0.0;
		ink.c.Kenpo = 1.0;
		ink.c.Kenle = 0.0;
		ink.c.Kshap = 1.0;
	} else if (inking == 1) {	/* Use 0.5  */
		ink.k_rule = icxKluma5;
		ink.c.Kstle = 0.5;
		ink.c.Kstpo = 0.0;
		ink.c.Kenpo = 1.0;
		ink.c.Kenle = 0.5;
		ink.c.Kshap = 1.0;
	} else if (inking == 2) {	/* Use maximum  */
		ink.k_rule = icxKluma5;
		ink.c.Kstle = 1.0;
		ink.c.Kstpo = 0.0;
		ink.c.Kenpo = 1.0;
		ink.c.Kenle = 1.0;
		ink.c.Kshap = 1.0;
	} else if (inking == 3) {	/* Use ramp  */
		ink.k_rule = icxKluma5;
		ink.c.Ksmth = ICXINKDEFSMTH;	/* Default curve smoothing */
		ink.c.Kstle = 0.0;
		ink.c.Kstpo = 0.0;
		ink.c.Kenpo = 1.0;
		ink.c.Kenle = 1.0;
		ink.c.Kshap = 1.0;
	} else if (inking == 4) {	/* Use locus  */
		ink.k_rule = icxKlocus;
	} else if (inking == 5) {	/* Use K target  */
		ink.k_rule = icxKvalue;
	} else {				/* Use specified curve */
		ink.k_rule = icxKluma5;
		ink.c.Ksmth = ICXINKDEFSMTH;	/* Default curve smoothing */
		ink.c.Kstle = Kstle;
		ink.c.Kstpo = Kstpo;
		ink.c.Kenpo = Kenpo;
		ink.c.Kenle = Kenle;
		ink.c.Kshap = Kshap;
	}

	/* Setup the viewing conditions */
	if (xicc_enum_viewcond(xicco, &vc, -1, 0))
		error ("%d, %s",xicco->errc, xicco->err);

//xicc_dump_viewcond(&vc);
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
//xicc_dump_viewcond(&vc);

	/* Get a expanded color conversion object */
	if ((luo = xicco->get_luobj(xicco, 0
#ifdef USE_NEARCLIP
	   | ICX_CLIP_NEAREST
#endif
	   | (intsep ? ICX_INT_SEPARATE : 0)
	   | (merge ? ICX_MERGE_CLUT : 0)
	   | (camclip ? ICX_CAM_CLIP : 0)
	                                  , func, intent, pcsor, order, &vc, &ink)) == NULL)
		error ("%d, %s",xicco->errc, xicco->err);

	/* Get details of conversion (Arguments may be NULL if info not needed) */
	if (invert)
		luo->spaces(luo, &outs, &outn, &ins, &inn, &alg, NULL, NULL, NULL);
	else
		luo->spaces(luo, &ins, &inn, &outs, &outn, &alg, NULL, NULL, NULL);

	/* If we can do check on clipped values */
	if (actual != 0) {
		if (invert == 0) {
			if (func == icmFwd || func == icmBwd) {
				if ((aluo = xicco->get_luobj(xicco, 0,
				     func == icmFwd ? icmBwd : icmFwd, intent, pcsor, order, &vc, &ink)) == NULL)
				error ("%d, %s",xicco->errc, xicco->err);
			}
		} else {
			aluo = luo;		/* We can use the same one */
		}
	}

	/* More information */
	if (verb) {
		int j;
		double inmin[MAX_CHAN], inmax[MAX_CHAN];
		double outmin[MAX_CHAN], outmax[MAX_CHAN];

		luo->get_native_ranges(luo, inmin, inmax, outmin,outmax);
		printf("Internal input value range: ");
		for (j = 0; j < inn; j++) {
			if (j > 0)
				fprintf(stdout," %f..%f",inmin[j], inmax[j]);
			else
				fprintf(stdout,"%f..%f",inmin[j], inmax[j]);
		}
		printf("\nInternal output value range: ");
		for (j = 0; j < outn; j++) {
			if (j > 0)
				fprintf(stdout," %f..%f",outmin[j], outmax[j]);
			else
				fprintf(stdout,"%f..%f",outmin[j], outmax[j]);
		}

		luo->get_ranges(luo, inmin, inmax, outmin,outmax);
		printf("\nInput value range: ");
		for (j = 0; j < inn; j++) {
			if (j > 0)
				fprintf(stdout," %f..%f",inmin[j], inmax[j]);
			else
				fprintf(stdout,"%f..%f",inmin[j], inmax[j]);
		}
		printf("\nOutput value range: ");
		for (j = 0; j < outn; j++) {
			if (j > 0)
				fprintf(stdout," %f..%f",outmin[j], outmax[j]);
			else
				fprintf(stdout,"%f..%f",outmin[j], outmax[j]);
		}
		printf("\n");
	}

	if (repYxy) {	/* report Yxy rather than XYZ */
		if (ins == icSigXYZData)
			ins = icSigYxyData; 
		if (outs == icSigXYZData)
			outs = icSigYxyData; 
	}

	if (repJCh) {	/* report JCh rather than Jab */
		if (ins == icxSigJabData)
			ins = icxSigJChData; 
		if (outs == icxSigJabData)
			outs = icxSigJChData; 
	}

	if (doplot) {
		int i, j;
		double xx[XRES];
		double yy[6][XRES];

		for (i = 0; i < XRES; i++) {
			double ival = (double)i/(XRES-1.0);

			/* Input is always Jab or Lab */
			in[0] = 100.0 * (1.0 - ival);
			in[1] = 0.0;
			in[2] = 0.0;
			in[3] = 0.0;

			/* Do the conversion */
			if (invert) {
				if ((rv = luo->inv_lookup(luo, out, in)) > 1)
					error ("%d, %s",xicco->errc,xicco->err);
			} else {
				if ((rv = luo->lookup(luo, out, in)) > 1)
					error ("%d, %s",xicco->errc,xicco->err);
			}

			xx[i] = in[0];
			for (j = 0; j < outn; j++)
				yy[j][i] = 100.0 * out[j];
		}

		/* plot order: Black Red Green Blue Yellow Purple */
		if (outs == icSigRgbData) {
			do_plot6(xx, NULL, yy[0], yy[1], yy[2], NULL, NULL, -XRES);

		} else if (outs == icSigCmykData) {
			do_plot6(xx, yy[3], yy[1], NULL, yy[0], yy[2], NULL, -XRES);

		} else {
		
			switch(outn) {
				case 1:
					do_plot6(xx, yy[0], NULL, NULL, NULL, NULL, NULL, -XRES);
					break;
				case 2:
					do_plot6(xx, yy[0], yy[1], NULL, NULL, NULL, NULL, -XRES);
					break;
				case 3:
					do_plot6(xx, yy[0], yy[1], yy[2], NULL, NULL, NULL, -XRES);
					break;
				case 4:
					do_plot6(xx, yy[0], yy[1], yy[2], yy[3], NULL, NULL, -XRES);
					break;
				case 5:
					do_plot6(xx, yy[0], yy[1], yy[2], yy[3], yy[4], NULL, -XRES);
					break;
				case 6:
					do_plot6(xx, yy[0], yy[1], yy[2], yy[3], yy[4], yy[5], -XRES);
					break;
			}
		}


	} else {
		/* Process colors to translate */
		for (;;) {
			int i,j;
			char *bp, *nbp;

			/* Read in the next line */
			if (fgets(buf, 200, stdin) == NULL)
				break;
			if (buf[0] == '#') {
				fprintf(stdout,"%s\n",buf);
				continue;
			}
			/* For each input number */
			for (bp = buf-1, nbp = buf, i = 0; i < MAX_CHAN; i++) {
				bp = nbp;
				out[i] = in[i] = oin[i] = strtod(bp, &nbp);
				if (nbp == bp)
					break;			/* Failed */
			}
			if (i == 0)
				break;

			/* If device data and scale */
			if(scale > 0.0
			 && ins != icxSigJabData
			 && ins != icxSigJChData
			 && ins != icSigXYZData
			 && ins != icSigLabData
			 && ins != icSigLuvData
			 && ins != icSigYCbCrData
			 && ins != icSigYxyData
			 && ins != icSigHsvData
			 && ins != icSigHlsData) {
				for (i = 0; i < MAX_CHAN; i++) {
					in[i] /= scale;
				}
			}

			if (repYxy && ins == icSigYxyData) {
				double Y = in[0];
				double x = in[1];
				double y = in[2];
				double z = 1.0 - x - y;
				double sum;
				if (y < 1e-6) {
					in[0] = in[1] = in[2] = 0.0;
				} else {
					sum = Y/y;
					in[0] = x * sum;
					in[1] = Y;
					in[2] = z * sum;
				}
			}

			/* JCh -> Jab */
			if (repJCh && ins == icxSigJChData) {
				double C = out[1];
				double h = out[2];
				out[1] = C * cos(3.14159265359/180.0 * h);
				out[2] = C * sin(3.14159265359/180.0 * h);
			}

			/* Do conversion */
			if (invert) {
				for (j = 0; j < MAX_CHAN; j++)
					out[j] = in[j];		/* Carry any auxiliary value to out for lookup */
				if ((rv = luo->inv_lookup(luo, out, in)) > 1)
					error ("%d, %s",xicco->errc,xicco->err);
			} else {
				if ((rv = luo->lookup(luo, out, in)) > 1)
					error ("%d, %s",xicco->errc,xicco->err);
			}

			if (repYxy && outs == icSigYxyData) {
				double X = out[0];
				double Y = out[1];
				double Z = out[2];
				double sum = X + Y + Z;
				if (sum < 1e-6) {
					out[0] = out[1] = out[2] = 0.0;
				} else {
					out[0] = Y;
					out[1] = X/sum;
					out[2] = Y/sum;
				}
			}

			/* Jab -> JCh */
			if (repJCh && outs == icxSigJChData) {
				double a = out[1];
				double b = out[2];
				out[1] = sqrt(a * a + b * b);
			    out[2] = (180.0/3.14159265359) * atan2(b, a);
				out[2] = (out[2] < 0.0) ? out[2] + 360.0 : out[2];
			}

			/* If device data and scale */
			if(scale > 0.0
			 && outs != icxSigJabData
			 && outs != icxSigJChData
			 && outs != icSigXYZData
			 && outs != icSigLabData
			 && outs != icSigLuvData
			 && outs != icSigYCbCrData
			 && outs != icSigYxyData
			 && outs != icSigHsvData
			 && outs != icSigHlsData) {
				for (i = 0; i < MAX_CHAN; i++) {
					out[i] *= scale;
				}
			}

			/* Output the results */
			for (j = 0; j < inn; j++) {
				if (j > 0)
					fprintf(stdout," %f",oin[j]);
				else
					fprintf(stdout,"%f",oin[j]);
			}
			printf(" [%s] -> %s -> ", icx2str(icmColorSpaceSignature, ins), icm2str(icmLuAlg, alg));

			for (j = 0; j < outn; j++) {
				if (j > 0)
					fprintf(stdout," %f",out[j]);
				else
					fprintf(stdout,"%f",out[j]);
			}
			printf(" [%s]", icx2str(icmColorSpaceSignature, outs));

			if (tlimit >= 0) {
				double tot;	
				for (tot = 0.0, j = 0; j < outn; j++) {
					tot += out[j];
				}
				printf(" Lim %f",tot);
			}
			if (rv == 0)
				fprintf(stdout,"\n");
			else {
				fprintf(stdout," (clip)\n");
				if (actual && aluo != NULL) {
					double cin[MAX_CHAN], de;
					if ((rv = aluo->lookup(aluo, cin, out)) > 1)
						error ("%d, %s",xicco->errc,xicco->err);

					for (de = 0.0, j = 0; j < inn; j++) {
						de += (cin[j] - in[j]) * (cin[j] - in[j]);
					}
					de = sqrt(de);
					printf("[Actual ");
					for (j = 0; j < inn; j++) {
						if (j > 0)
							fprintf(stdout," %f",cin[j]);
						else
							fprintf(stdout,"%f",cin[j]);
					}
					printf(", deltaE %f]\n",de);
				}
			}
		}
	}

	/* Done with lookup object */
	if (aluo != NULL && aluo != luo)
		luo->del(aluo);
	luo->del(luo);

	xicco->del(xicco);		/* Expansion wrapper */
	icco->del(icco);		/* Icc */
	fp->del(fp);

	return 0;
}

