/* 
 * Argyll Color Correction System
 * Color Device profile generator.
 *
 * Author: Graeme W. Gill
 * Date:   15/2/97
 *
 * Copyright 1996-2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/*
 * This program takes in the scattered test chart
 * points, and interpolates them into a gridded
 * forward ICC device profile, as well as creating
 * backward conversions based on the forward grid.
 * 
 * Preview profiles are not currently generated.
 * 
 * The gamut clut should be implemented with xicc/rspl
 */

/*
 * TTBD:
 *      Add Argyll private tag to record ink limit etc. to automate link parameters.
 *      Estimate ink limit from B2A tables if no private tag ?
 *      Add used option for black relative
 *      Add used option for separate high res reverse tables
 *      Fix 400% assumptions for > 4 color devices ?
 *
 *      Should allow creating profiles from .MPP directly for <= 4 dev channels.
 *      Should allow creating profiles from existing ICC profiles (deprecate revfix ?)
 *
 *      Should allow creating profiles >4 channels by providing .MPP for input,
 *      dev link .icm for psudo-dev to device & .ti3 for Pseudo-dev to PCS.
 *		Note gamut should come from psudo-dev to PCS.
 */

#define DEBUG

#define verbo stdout

#include <stdio.h>
#include "copyright.h"
#include "config.h"
#include "cgats.h"
#include "xicc.h"
#include "prof.h"

#define DEFAVGDEV 0.25		/* Default average deviation percentage */

void usage(char *diag, ...) {
	int i;
	fprintf(stderr,"Create ICC profile, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: %s [-options] inoutfile\n",error_program);
	fprintf(stderr," -v         Verbose mode\n");
	fprintf(stderr," -A manufacturer Manufacturer description string\n");
	fprintf(stderr," -M model        Model description string\n");
	fprintf(stderr," -D description  Profile Description string (Default \"inoutfile\")\n");
	fprintf(stderr," -C copyright    Copyright string\n");
	fprintf(stderr," -q [lmhu]  Quality - Low, Medium (def), High, Ultra\n");
	fprintf(stderr," -b         B2A table will not be used (V. low quality, -B Extremely low)\n");
	fprintf(stderr," -y         Verify A2B profile\n");
	fprintf(stderr," -ni        Don't create input (Device) shaper curves\n");
	fprintf(stderr," -no        Don't create output (PCS) shaper curves\n");
	fprintf(stderr," -k [zhxr]  Black generation: z = zero K,\n");
	fprintf(stderr,"            h = 0.5 K (def), x = max K, r = ramp K\n");
	fprintf(stderr," -k p stle stpo enpo enle shape\n");
	fprintf(stderr,"            stle: K level at White 0.0 - 1.0\n");
	fprintf(stderr,"            stpo: start point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"            enpo: End point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"            enle: K level at Black 0.0 - 1.0\n");
	fprintf(stderr,"            shape: 1.0 = straight, 0.0-1.0 concave, 1.0-2.0 convex\n");
	fprintf(stderr," -l tlimit  set total ink limit, 0 - 400%% (default none)\n");
	fprintf(stderr," -L klimit  set black ink limit, 0 - 100%% (default none)\n");
	fprintf(stderr," -a [lxgs]  Algorithm type override\n");
	fprintf(stderr,"            l = Lab clut, x = XYZ lut\n");
	fprintf(stderr,"            g = gamma+matrix, s = shaper+matrix\n");
	fprintf(stderr,"            G = single gamma+matrix, S = single shaper+matrix\n");
	fprintf(stderr," -u         If Lut input profile, make it absolute (non-standard)\n");
	fprintf(stderr," -i illum   Choose illuminant for print/transparency spectral data:\n");
	fprintf(stderr,"            A, D50 (def.), D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ  Choose CIE Observer for spectral data:\n");
	fprintf(stderr,"            1931_2, 1964_10, S&B 1955_2, 1964_10c, shaw, J&V 1978_2 (def.)\n");
	fprintf(stderr," -f         Use Fluorescent Whitening Agent compensation\n");
	fprintf(stderr," -r avgdev  Average deviation of device+instrument readings as a percentage (default %3.1f%%)\n",DEFAVGDEV);
/*	fprintf(stderr," -rs smooth  RSPL suplimental optimised smoothing factor\n"); */
/*	fprintf(stderr," -rr smooth  RSPL raw underlying smoothing factor\n"); */
	fprintf(stderr," -s src.icc Apply gamut mapping to perceptual B2A table for given source space\n");
	fprintf(stderr," -S src.icc Apply gamut mapping to perceptual and saturation B2A table\n");
	fprintf(stderr," -g src.gam Use source image gamut as well for gamut mapping\n");
	fprintf(stderr," -p absprof Include abstract profile in output tables\n");
	fprintf(stderr," -t intent  Override gamut mapping intent for perceptual table:\n");
	fprintf(stderr," -T intent  Override gamut mapping intent for saturation table:\n");
	for (i = 0; ; i++) {
		icxGMappingIntent gmi;
		if (xicc_enum_gmapintent(&gmi, i))
			break;
		fprintf(stderr,"               %d: %s\n",i,gmi.desc);
	}
	fprintf(stderr," -c viewcond  set input viewing conditions for %s gamut mapping,\n",icxcam_description(cam_default));
	fprintf(stderr,"               either an enumerated choice, or a parameter\n");
	fprintf(stderr," -d viewcond  set output viewing conditions for %s gamut mapping\n",icxcam_description(cam_default));
	fprintf(stderr,"               either an enumerated choice, or a parameter\n");
	fprintf(stderr,"               Also enables out of gamut clipping CAM space.\n");
	fprintf(stderr,"               Enumerated Viewing Conditions:\n");
	for (i = 0; ; i++) {
		icxViewCond vc;
		if (xicc_enum_viewcond(NULL, &vc, i, 1))
			break;

		fprintf(stderr,"               %d: %s\n",i,vc.desc);
	}
	fprintf(stderr," inoutfile  Base name for input.ti3/output.icm file\n");
	exit(1);
}


main(int argc, char *argv[])
{
	int fa,nfa;					/* current argument we're looking at */
	int verb = 0;
	int iquality = 1;			/* A2B quality */
	int oquality = -1;			/* B2A quality same as A2B */
	int verify = 0;
	int noiluts = 0;			/* No input shaper luts */
	int nooluts = 0;			/* No output shaper luts */
	int nsabs = 0;				/* Make non-standard absolute input lut profile */
	int inking = 1;				/* Default 0.5 K */
	double Kstle, Kstpo, Kenle, Kenpo, Kshap;
	int tlimit = -1;			/* Total ink limit as a % */
	int klimit = -1;			/* Black ink limit as a % */
	int fwacomp = 0;			/* FWA compensation */
	double avgdev = DEFAVGDEV/100.0;	/* Average measurement deviation */
	double smooth = 1.0;		/* RSPL Smoothness factor (verification) */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType illum = icxIT_D50;	/* Spectral defaults */
	xspect cust_illum;			/* Custom illumination spectrum */
	icxObserverType observ = icxOT_Judd_Voss_2;
	char ipname[MAXNAMEL+1] = "\000";	/* Input icc profile - enables gamut map */
	char sgname[MAXNAMEL+1] = "\000";	/* Image source gamut name */
	char absname[MAXNAMEL+1] = "\000";	/* Abstract profile name */
	int sepsat = 0;				/* Create separate saturation B2A table */
	icxViewCond ivc_p;			/* Input Viewing Parameters for CAM */
	icxViewCond ovc_p;			/* Output Viewing Parameters for CAM (enables CAM clip) */
	int ivc_e = -1, ovc_e = -1;	/* Enumerated viewing condition */
	icxGMappingIntent pgmi;		/* default Perceptual gamut mapping intent */
	icxGMappingIntent sgmi;		/* default Saturation gamut mapping intent */
	char baname[MAXNAMEL+1] = "\000";	/* Input & Output base name */
	char inname[MAXNAMEL+1] = "\000";	/* Input cgats file base name */
	char outname[MAXNAMEL+1] = "\000";	/* Output cgats file base name */
	cgats *icg;			/* input cgats structure */
	int ti;				/* Temporary CGATs index */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	prof_atype ptype = prof_default;	/* Default for each type of device */
	profxinf xpi;		/* Extra profile information */

	error_program = argv[0];
	memset((void *)&xpi, 0, sizeof(profxinf));	/* Init extra profile info to defaults */

	/* Init VC overrides so that we know when the've been set */
	ivc_p.Ev = -1;
	ivc_p.Wxyz[0] = -1.0; ivc_p.Wxyz[1] = -1.0; ivc_p.Wxyz[2] = -1.0;
	ivc_p.Yb = -1.0;
	ivc_p.La = -1.0;
	ivc_p.Lv = -1.0;
	ivc_p.Yf = -1.0;
	ivc_p.Fxyz[0] = -1.0; ivc_p.Fxyz[1] = -1.0; ivc_p.Fxyz[2] = -1.0;

	ovc_p.Ev = -1;
	ovc_p.Wxyz[0] = -1.0; ovc_p.Wxyz[1] = -1.0; ovc_p.Wxyz[2] = -1.0;
	ovc_p.Yb = -1.0;
	ovc_p.La = -1.0;
	ovc_p.Lv = -1.0;
	ovc_p.Yf = -1.0;
	ovc_p.Fxyz[0] = -1.0; ovc_p.Fxyz[1] = -1.0; ovc_p.Fxyz[2] = -1.0;

	xicc_enum_gmapintent(&pgmi, icxPerceptualGMIntent);
	xicc_enum_gmapintent(&sgmi, icxSaturationGMIntent);

	if (argc < 2)
		usage("Too few arguments, got %d expect at least %d",argc-1,1);

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {

		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {	/* Look for any flags */
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

			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V')
				verb = 1;

			/* Manufacturer description string */
			else if (argv[fa][1] == 'A') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to manufacturer description flag -A");
				xpi.deviceMfgDesc = na;
			}

			/* Model description string */
			else if (argv[fa][1] == 'M') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to model description flag -M");
				xpi.modelDesc = na;
			}

			/* Profile Description */
			else if (argv[fa][1] == 'D'
			      || argv[fa][1] == 'E') {		/* Backwards compatibility with old versions */
				fa = nfa;
				if (na == NULL) usage("Expect argument to profile description flag -E");
				xpi.profDesc = na;
			}

			/* Copyright string */
			else if (argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to copyright flag -C");
				xpi.copyright = na;
			}

			/* Quality */
			else if (argv[fa][1] == 'q' || argv[fa][1] == 'Q') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to quality flag -q");
    			switch (na[0]) {
					case 'l':
					case 'L':
						iquality = 0;
						break;
					case 'm':
					case 'M':
						iquality = 1;
						break;
					case 'h':
					case 'H':
						iquality = 2;
						break;
					case 'u':
					case 'U':
						iquality = 3;
						break;
					default:
						usage("Unknown argument '%c' to quality flag -q",na[0]);
				}
			}
			else if (argv[fa][1] == 'b')
				oquality = 0;

			else if (argv[fa][1] == 'B')
				oquality = -2;

			else if (argv[fa][1] == 'y' || argv[fa][1] == 'Y')
				verify = 1;

			/* Disable input or output luts */
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N') {
				fa = nfa;
				if (na == NULL) {	/* Backwards compatible */
					nooluts = 1;
				} else {

					if (na[0] != 'i' && na[0] != 'I'
					 && na[0] != 'o' && na[0] != 'O')
						usage("Unknown argument '%c' to flag -n",na[0]);
	
					if (na[0] == 'i' || na[0] == 'I')
						noiluts = 1;
					if (na[0] == 'o' || na[0] == 'O')
						nooluts = 1;
				}
			}

			else if (argv[fa][1] == 'u' || argv[fa][1] == 'U')
				nsabs = 1;

			/* Inking rule */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to inking flag -k");
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
					case 'p':
					case 'P':
						inking = 4;		/* Use curve parameter */
						++fa;
						if (fa >= argc) usage("Too few arguments to inking flag -kp");
						Kstle = atof(argv[fa]);

						++fa;
						if (fa >= argc) usage("Too few arguments to inking flag -kp");
						Kstpo = atof(argv[fa]);

						++fa;
						if (fa >= argc || argv[fa][0] == '-') usage("Too few arguments to inking flag -kp");
						Kenpo = atof(argv[fa]);

						++fa;
						if (fa >= argc) usage("Too few arguments to inking flag -kp");
						Kenle = atof(argv[fa]);

						++fa;
						if (fa >= argc || argv[fa][0] == '-') usage("Too few arguments to inking flag -kp");
						Kshap = atof(argv[fa]);
						break;
					default:
						usage("Unknown inking rule (-k) argument '%c'",na[0]);
				}
			}

			/* Total Ink Limit */
			else if (argv[fa][1] == 'l') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to total ink limit flag -l");
				tlimit = atoi(na);
			}

			/* Black Ink Limit */
			else if (argv[fa][1] == 'L') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to black ink limit flag -L");
				klimit = atoi(na);
			}

			/* Algorithm type */
			else if (argv[fa][1] == 'a' || argv[fa][1] == 'A') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to algorithm flag -a");
    			switch (na[0]) {
					case 'l':
					case 'L':
						ptype = prof_clutLab;
						break;
					case 'x':
					case 'X':
						ptype = prof_clutXYZ;
						break;
					case 'g':
						ptype = prof_gammat;
						break;
					case 'G':
						ptype = prof_gam1mat;
						break;
					case 's':
						ptype = prof_shamat;
						break;
					case 'S':
						ptype = prof_sha1mat;
						break;
					default:
						usage("Unknown argument '%c' to algorithm flag -a",na[0] );
				}
			}
			/* Spectral Illuminant type */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to illuminant flag -i");
				if (strcmp(na, "A") == 0) {
					spec = 1;
					illum = icxIT_A;
				} else if (strcmp(na, "D50") == 0) {
					spec = 1;
					illum = icxIT_D50;
				} else if (strcmp(na, "D65") == 0) {
					spec = 1;
					illum = icxIT_D65;
				} else if (strcmp(na, "F5") == 0) {
					spec = 1;
					illum = icxIT_F5;
				} else if (strcmp(na, "F8") == 0) {
					spec = 1;
					illum = icxIT_F8;
				} else if (strcmp(na, "F10") == 0) {
					spec = 1;
					illum = icxIT_F10;
				} else {	/* Assume it's a filename */
					spec = 1;
					illum = icxIT_custom;
					if (read_xspect(&cust_illum, na) != 0)
						usage("Failed to read custom illuminant spectrum in file '%s'",na);
				}
			}

			/* Spectral Observer type */
			else if (argv[fa][1] == 'o' || argv[fa][1] == 'O') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to observer flag -o");
				if (strcmp(na, "1931_2") == 0) {			/* Classic 2 degree */
					spec = 1;
					observ = icxOT_CIE_1931_2;
				} else if (strcmp(na, "1964_10") == 0) {	/* Classic 10 degree */
					spec = 1;
					observ = icxOT_CIE_1964_10;
				} else if (strcmp(na, "1955_2") == 0) {		/* Stiles and Burch 1955 2 degree */
					spec = 1;
					observ = icxOT_Stiles_Burch_2;
				} else if (strcmp(na, "1978_2") == 0) {		/* Judd and Voss 1978 2 degree */
					spec = 1;
					observ = icxOT_Judd_Voss_2;
				} else if (strcmp(na, "shaw") == 0) {		/* Shaw and Fairchilds 1997 2 degree */
					spec = 1;
					observ = icxOT_Shaw_Fairchild_2;
				} else
					usage("Unrecognised argument '%s' to observer flag -o",na);
			}

			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F')
				fwacomp = 1;

			/* Average Deviation percentage */
			else if (argv[fa][1] == 'r') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to average deviation flag -r");
				if (na[0] == 's') {			/* (verification) */
					smooth = atof(na+1);
					if (smooth < 0.0)
						usage("Optimised smoothing factor argument to '-rs' must be over 0.0");
				} else if (na[0] == 'r') {	/* (verification) */
					smooth = atof(na+1);
					if (smooth < 0.0)
						usage("Raw smoothing factor argument to '-rr' must be over 0.0");
					smooth = -smooth;		/* Signal raw factor */
				} else {
					avgdev = 0.01 * atof(na);
					if (avgdev < 0.0 || avgdev > 1.0)
						usage("Average Deviation argument must be between 0.0 and 100.0");
				}
			}

			/* Percetual Source Gamut and Perceptual/Saturation Gamut Maping mode enable */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				if (argv[fa][1] == 'S')
					sepsat = 1;
				if (na == NULL)
					usage("Unrecognised argument to source gamut flag -%c",argv[fa][1]);

				fa = nfa;
				strncpy(ipname,na,MAXNAMEL); ipname[MAXNAMEL] = '\000';
			}

			/* Source image gamut */
			else if (argv[fa][1] == 'g' || argv[fa][1] == 'G') {
				if (na == NULL)
					usage("Unrecognised argument to source image gamut flag -g",argv[fa][1]);
				fa = nfa;
				strncpy(sgname,na,MAXNAMEL); sgname[MAXNAMEL] = '\000';
			}

			/* Abstract profile */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				if (na == NULL) usage("Expected abstract profile filename after -p");
				fa = nfa;
				strncpy(absname,na,MAXNAMEL); absname[MAXNAMEL] = '\000';
			}

			/* Perceptual Mapping intent override */
			else if (argv[fa][1] == 't') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to perceptul intent override flag -t");
				if (na[0] >= '0' && na[0] <= '9') {
					int i = atoi(na);
					if (xicc_enum_gmapintent(&pgmi, i) != 0)
						usage("Unrecognised intent '%s' to perceptual override flag -t",na);
				}
			}

			/* Saturation Mapping intent override */
			else if (argv[fa][1] == 'T') {
				int i = icxIllegalGMIntent;
				fa = nfa;
				if (na == NULL) usage("Expect argument to saturation intent override flag -T");
				if (na[0] >= '0' && na[0] <= '9') {
					i = atoi(na);
				} else {
	    			switch (na[0]) {
						case 'p':
						case 'P':
							i = icxPerceptualGMIntent;
							break;
						case 'r':
						case 'R':
							i = icxRelativeGMIntent;
							break;
						case 's':
						case 'S':
							i = icxSaturationGMIntent;
							break;
						case 'a':
						case 'A':
							i = icxAbsoluteGMIntent;
							break;
					}
				}
				if (xicc_enum_gmapintent(&sgmi, i) != 0)
					usage("Unrecognised intent '%s' to saturation override flag -T",na);
			}

			/* Viewing conditions */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'd') {
				icxViewCond *vc;

				if (argv[fa][1] == 'c') {
					vc = &ivc_p;
				} else {
					vc = &ovc_p;
				}

				fa = nfa;
				if (na == NULL) usage("Viewing conditions flag (-c) needs an argument");
				if (na[0] >= '0' && na[0] <= '9') {
					if (vc == &ivc_p)
						ivc_e = atoi(na);
					else
						ovc_e = atoi(na);
				} else if (na[0] == 's' || na[0] == 'S') {
					if (na[1] != ':')
						usage("Viewing conditions (-cs) missing ':'");
					if (na[2] == 'a' || na[2] == 'A') {
						vc->Ev = vc_average;
					} else if (na[2] == 'm' || na[2] == 'M') {
						vc->Ev = vc_dim;
					} else if (na[2] == 'd' || na[2] == 'D') {
						vc->Ev = vc_dark;
					} else if (na[2] == 'c' || na[2] == 'C') {
						vc->Ev = vc_cut_sheet;
					} else
						usage("Viewing condition (-c) unrecognised surround '%c'",na[2]);
				} else if (na[0] == 'w' || na[0] == 'W') {
					double x, y, z;
					if (sscanf(na+1,":%lf:%lf:%lf",&x,&y,&z) == 3) {
						vc->Wxyz[0] = x; vc->Wxyz[1] = y; vc->Wxyz[2] = z;
					} else if (sscanf(na+1,":%lf:%lf",&x,&y) == 2) {
						vc->Wxyz[0] = x; vc->Wxyz[1] = y;
					} else
						usage("Viewing condition (-cw) unrecognised white point '%s'",na+1);
				} else if (na[0] == 'a' || na[0] == 'A') {
					if (na[1] != ':')
						usage("Viewing conditions (-ca) missing ':'");
					vc->La = atof(na+2);
				} else if (na[0] == 'b' || na[0] == 'B') {
					if (na[1] != ':')
						usage("Viewing conditions (-cb) missing ':'");
					vc->Yb = atof(na+2)/100.0;
				} else if (na[0] == 'f' || na[0] == 'F') {
					double x, y, z;
					if (sscanf(na+1,":%lf:%lf:%lf",&x,&y,&z) == 3) {
						vc->Fxyz[0] = x; vc->Fxyz[1] = y; vc->Fxyz[2] = z;
					} else if (sscanf(na+1,":%lf:%lf",&x,&y) == 2) {
						vc->Fxyz[0] = x; vc->Fxyz[1] = y;
					} else if (sscanf(na+1,":%lf",&x) == 1) {
						vc->Yf = x/100.0;
					} else
						usage("Viewing condition (-cf) unrecognised flare '%s'",na+1);
				} else
					usage("Viewing condition (-c) unrecognised sub flag '%c'",na[0]);
			}

			else 
				usage("Unknown flag '%c'",argv[fa][1]);
		} else
			break;
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage("Missing input .ti3 and output ICC basename");
	strncpy(baname,argv[fa++],MAXNAMEL-4); baname[MAXNAMEL-4] = '\000';
	if (xpi.profDesc == NULL)
		xpi.profDesc = baname;	/* Default description */
	strcpy(inname,baname);
	strcpy(outname,baname);
	strcat(inname,".ti3");
#if defined (UNIX) || defined(__APPLE__)
	strcat(outname,".icc");
#else
	strcat(outname,".icm");
#endif

	if (fwacomp && spec == 0)
		error ("FWA compensation only works when viewer and/or illuminant selected");

	if (oquality == -1) {		/* B2A tables will be used */
		oquality = iquality;
	}

	/* Open and look at the .ti3 profile patches file */
	icg = new_cgats();			/* Create a CGATS structure */
	icg->add_other(icg, "CTI3"); 	/* our special input type is Calibration Target Information 3 */
	icg->add_other(icg, "CAL"); 	/* our special device Calibration state */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI3 format file");
	if (icg->ntables < 1)
		error ("Input file doesn't contain at least one table");

	/* See if CIE is actually available - some sources of .TI3 don't provide it */
	if (!spec
	 && icg->find_field(icg, 0, "LAB_L") < 0
	 && icg->find_field(icg, 0, "XYZ_X") < 0) {

		if (icg->find_kword(icg, 0, "SPECTRAL_BANDS") < 0)
			error ("Neither CIE nor spectral data found in file '%s'",inname);

		/* Switch to using spectral information */
		if (verb)
			printf("No CIE data found, switching to spectral with standard observer & D50\n");
		spec = 1;
		illum = icxIT_D50;
		observ = icxOT_CIE_1931_2;
	}
	
	/* If we requested spectral, check that it is available */
	if (spec) {
		if (icg->find_kword(icg, 0, "SPECTRAL_BANDS") < 0)
			error ("Requested spectral interpretation when data not available");
	}

	/* read the device class, and call function to create profile. */
	if ((ti = icg->find_kword(icg, 0, "DEVICE_CLASS")) < 0)
		error ("Input file doesn't contain keyword DEVICE_CLASS");

	if (strcmp(icg->t[0].kdata[ti],"OUTPUT") == 0) {
		icxInk ink;							/* Ink parameters */

		if ((ti = icg->find_kword(icg, 0, "TOTAL_INK_LIMIT")) >= 0) {
			double imax;
			imax = atof(icg->t[0].kdata[ti]);
			if (imax > 1e-4 && imax <= 400.0) {
				if (tlimit > 1e-4 && tlimit <= 400.0) {	/* User has specified limit as option */
					if (imax < tlimit) {
						warning("Ink limit greater than original chart! (%f > %f)",tlimit,imax);
					}
				} else {
					if (imax > 80.0)
						tlimit = imax - 10.0;	/* Rule of thumb - 10% below chart maximum */
					else
						tlimit = imax;
				}
			} else {
				tlimit = -1;
			}
		}

		if (tlimit >= 0) {
			if (verb)
				printf("Total ink limit being used is %d%%\n",tlimit);
			ink.tlimit = tlimit/100.0;	/* Set a total ink limit */
		} else {
			if (verb)
				printf("No total ink limit being used\n");
			ink.tlimit = -1.0;			/* Don't use a limit */
		}

		if (klimit >= 0) {
			if (verb)
				printf("Black ink limit being used is %d%%\n",klimit);
			ink.klimit = klimit/100.0;	/* Set a black ink limit */
		} else {
			if (verb)
				printf("No black ink limit being used\n");
			ink.klimit = -1.0;			/* Don't use a limit */
		}

		ink.c.Ksmth = ICXINKDEFSMTH;	/* default black curve smoothing */

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
			ink.c.Kstle = 0.0;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 1.0;
			ink.c.Kshap = 1.0;
		} else {				/* Use specified curve */
			ink.k_rule = icxKluma5;
			ink.c.Kstle = Kstle;
			ink.c.Kstpo = Kstpo;
			ink.c.Kenpo = Kenpo;
			ink.c.Kenle = Kenle;
			ink.c.Kshap = Kshap;
		}

		if (ptype == prof_default)
			ptype = prof_clutLab;
		else if (ptype != prof_clutLab && ptype != prof_clutXYZ) {
			error ("Output profile can only be a clut algorithm");
		}

		make_output_icc(ptype, verb, iquality, oquality, noiluts, nooluts, verify, &ink, outname,
		                icg, spec, illum, &cust_illum, observ, fwacomp, smooth, avgdev,
		                ipname[0] != '\000' ? ipname : NULL,
		                sgname[0] != '\000' ? sgname : NULL,
		                absname[0] != '\000' ? absname : NULL,
						sepsat, &ivc_p, &ovc_p, ivc_e, ovc_e,
						&pgmi, &sgmi, &xpi);

	} else if (strcmp(icg->t[0].kdata[ti],"INPUT") == 0) {

		if (ptype == prof_default)
			ptype = prof_clutLab;		/* For best possible quality */
		make_input_icc(ptype, verb, iquality, verify, nsabs, outname, icg, smooth, avgdev, &xpi);

	} else if (strcmp(icg->t[0].kdata[ti],"DISPLAY") == 0) {

		if (fwacomp)
			error ("FWA compensation isn't applicable to a display device");

		if (ptype == prof_default)
			ptype = prof_clutLab;		/* ?? or should it default to prof_shamat ?? */

		/* ICC V2.3 doesn't have display intents. */
		/* ICC V2.4 does. What should we do here ? */
		make_output_icc(ptype, verb, iquality, oquality, noiluts, nooluts, verify, NULL, outname,
		                icg, spec, illum, &cust_illum, observ, 0, smooth, avgdev,
		                NULL, NULL,
		                absname[0] != '\000' ? absname : NULL,
						0, &ivc_p, &ovc_p, ivc_e, ovc_e,
						&pgmi, &sgmi, &xpi);

	} else
		error ("Input file keyword DEVICE_CLASS has unknown value");

	icg->del(icg);		/* Clean up */

	return 0;
}










