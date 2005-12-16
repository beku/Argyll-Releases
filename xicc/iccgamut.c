
/* 
 * iccgamut
 *
 * Produce color surface gamut of an ICC profile.
 *
 * Author:  Graeme W. Gill
 * Date:    19/3/00
 * Version: 1.00
 *
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the LICENCE.TXT file for licencing details.
 */


/*
 * TTBD:
 *       To support CIACAM02 properly, need to cope with viewing parameters ?
 */

#define SURFACE_ONLY
#define GAMRES 10.0		/* Default surface resolution */

#define RGBRES 33	/* 33 */
#define CMYKRES 17	/* 17 */

#undef RANDOMPOINTS
#undef RANDCUBE

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "numlib.h"
#include "icc.h"
#include "xicc.h"
#include "gamut.h"

#ifdef NEVER
void error(char *fmt, ...), warning(char *fmt, ...);
#endif

void usage(char *diag) {
	int i;
	fprintf(stderr,"Create Lab/Jab gamut plot\n");
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	fprintf(stderr,"usage: iccgamut [options] profile\n");
	if (diag != NULL)
		fprintf(stderr,"Diagnostic: %s\n",diag);
	fprintf(stderr," -v            Verbose\n");
	fprintf(stderr," -d sres       Surface resolution details 1.0 - 50.0\n");
	fprintf(stderr," -w            emit VRML .wrl file as well as CGATS .gam file\n");
	fprintf(stderr," -n            Don't add VRML axes or white/black point\n");
	fprintf(stderr," -f function   f = forward*, b = backwards\n");
	fprintf(stderr," -i intent     p = perceptual, r = relative colorimetric,\n");
	fprintf(stderr,"               s = saturation, a = absolute*,\n");
	fprintf(stderr,"               j = Appearance %s, d = default\n",icxcam_description(cam_default));
	fprintf(stderr," -o order      n = normal (priority: lut > matrix > monochrome)\n");
	fprintf(stderr,"               r = reverse (priority: monochrome > matrix > lut)\n");
	fprintf(stderr," -l tlimit     set total ink limit, 0 - 400%% (default none)\n");
	fprintf(stderr," -L klimit     set black ink limit, 0 - 100%% (default none)\n");
	fprintf(stderr," -c viewcond   set viewing conditions for %s,\n",icxcam_description(cam_default));
	fprintf(stderr,"               either an enumerated choice, or a series of parameters:\n");
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
	fprintf(stderr,"\n");
	exit(1);
}

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	char prof_name[100];
	char *xl, out_name[100];
	icmFile *fp;
	icc *icco;
	xicc *xicco;
	gamut *gam;
	int verb = 0;
	int rv = 0;
	int vrml = 0;
	int doaxes = 1;
	double gamres = GAMRES;				/* Surface resolution */
	icxInk ink;							/* Ink parameters */
	int tlimit = -1;			/* Total ink limit as a % */
	int klimit = -1;			/* Black ink limit as a % */
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
	char buf[200];

	icxLuBase *luo;

	/* Lookup parameters */
	icmLookupFunc     func   = icmFwd;				/* Default */
	icRenderingIntent intent = icAbsoluteColorimetric;	/* Default */
	icColorSpaceSignature pcsor = icSigLabData;	/* Default */
	icmLookupOrder    order  = icmLuOrdNorm;		/* Default */
	
#ifdef NUMSUP_H
	error_program = "iccgamut";
#endif

	if (argc < 2)
		usage("Too few parameters");

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
				usage(NULL);

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
					default:
						usage("Unrecognised parameter after flag -f");
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
					case 'J':
						intent = icxAppearance;
						pcsor = icxSigJabData;
						break;
					case 'd':
					case 'D':
						intent = icmDefaultIntent;
						break;
					default:
						usage("Unrecognised parameter after flag -i");
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
						usage("Unrecognised parameter after flag -o");
				}
			}

			/* Verbosity */
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;
			}
			/* VRML output */
			else if (argv[fa][1] == 'w' || argv[fa][1] == 'W') {
				vrml = 1;
			}
			/* No axis output */
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N') {
				doaxes = 0;
			}
			/* Ink limit */
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


			/* Surface Detail */
			else if (argv[fa][1] == 'd' || argv[fa][1] == 'D') {
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -d");
				gamres = atof(na);
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

	if (fa >= argc || argv[fa][0] == '-') usage("Expected profile name");
	strcpy(prof_name,argv[fa]);

	/* Open up the profile for reading */
	if ((fp = new_icmFileStd_name(prof_name,"r")) == NULL)
		error ("Can't open file '%s'",prof_name);

	if ((icco = new_icc()) == NULL)
		error ("Creation of ICC object failed");

	if ((rv = icco->read(icco,fp,0)) != 0)
		error ("%d, %s",rv,icco->err);

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

	/* Setup ink limit */
	if (tlimit >= 0)
		ink.tlimit = tlimit/100.0;	/* Set a total ink limit */
	else
		ink.tlimit = -1.0;			/* Don't use a limit */

	if (klimit >= 0)
		ink.klimit = klimit/100.0;	/* Set a black ink limit */
	else
		ink.klimit = -1.0;			/* Don't use a limit */

	/* Setup a safe ink generation (not used) */
	ink.k_rule = icxKluma5;
	ink.c.Ksmth = ICXINKDEFSMTH;	/* Default smoothing */
	ink.c.Kstle = 0.0;		/* Min K at white end */
	ink.c.Kstpo = 0.0;		/* Start of transition is at white */
	ink.c.Kenle = 1.0;		/* Max K at black end */
	ink.c.Kenpo = 1.0;		/* End transition at black */
	ink.c.Kshap = 1.0;		/* Linear transition */

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
	if ((luo = xicco->get_luobj(xicco, 0, func, intent, pcsor, order, &vc, &ink)) == NULL)
		error ("%d, %s",xicco->errc, xicco->err);

	/* Creat a gamut surface */
	if ((gam = luo->get_gamut(luo, gamres)) == NULL)
		error ("%d, %s",xicco->errc, xicco->err);

	strcpy(out_name, prof_name);
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

	if (verb) {
		printf("Total volume of gamut is %f cubic colorspace units\n",gam->volume(gam));
	}
	gam->del(gam);

	luo->del(luo);			/* Done with lookup object */

	xicco->del(xicco);		/* Expansion wrapper */
	icco->del(icco);			/* Icc */
	fp->del(fp);


	return 0;
}


#ifdef NEVER
/* Basic printf type error() and warning() routines */

void
error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"icclu: Error - ");
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

	fprintf(stderr,"icclu: Warning - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}

#endif /* NEVER */