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
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes in the scattered test chart
 * points, and interpolates them into a gridded
 * forward ICC device profile, as well as creating
 * backward conversions based on the forward grid.
 * 
 * Preview profiles are not currently generated.
 * 
 * The gamut cLUT should be implemented with xicc/rspl
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

#undef DEBUG
#undef DO_TIME			/* Time the operation */

#define verbo stdout

#include <stdio.h>
#include <string.h>
#include <time.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "cgats.h"
#include "xicc.h"
#include "prof.h"

#define DEFAVGDEV 0.5		/* Default average deviation percentage */
							/* This equates to a uniform added error of +/- 1% */

void usage(char *diag, ...) {
	int i;
	fprintf(stderr,"Create ICC profile, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: %s [-options] inoutfile\n",error_program);
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -A manufacturer Manufacturer description string\n");
	fprintf(stderr," -M model        Model description string\n");
	fprintf(stderr," -D description  Profile Description string (Default \"inoutfile\")\n");
	fprintf(stderr," -C copyright    Copyright string\n");
	fprintf(stderr," -q lmhu         Quality - Low, Medium (def), High, Ultra\n");
	fprintf(stderr," -b [lmhun]      Low quality B2A table - or specific B2A quality or none for input device\n");
	fprintf(stderr," -y              Verify A2B profile\n");
	fprintf(stderr," -ni             Don't create input (Device) shaper curves\n");
	fprintf(stderr," -np             Don't create input (Device) grid position curves\n");
	fprintf(stderr," -no             Don't create output (PCS) shaper curves\n");
	fprintf(stderr," -nc             Don't put the input .ti3 data in the profile\n");
	fprintf(stderr," -k zhxr         Black value target: z = zero K,\n");
	fprintf(stderr,"                 h = 0.5 K, x = max K, r = ramp K (def.)\n");
	fprintf(stderr," -k p stle stpo enpo enle shape\n");
	fprintf(stderr,"                 stle: K level at White 0.0 - 1.0\n");
	fprintf(stderr,"                 stpo: start point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"                 enpo: End point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"                 enle: K level at Black 0.0 - 1.0\n");
	fprintf(stderr,"                 shape: 1.0 = straight, 0.0-1.0 concave, 1.0-2.0 convex\n");
	fprintf(stderr," -K parameters   Same as -k, but target is K locus rather than K value itself\n");
	fprintf(stderr," -l tlimit       override total ink limit, 0 - 400%% (default from .ti3)\n");
	fprintf(stderr," -L klimit       override black ink limit, 0 - 100%% (default from .ti3)\n");
	fprintf(stderr," -a lxXgsmGS     Algorithm type override\n");
	fprintf(stderr,"                 l = Lab cLUT (def.), x = XYZ cLUT, X = display XYZ cLUT + matrix\n");
	fprintf(stderr,"                 g = gamma+matrix, s = shaper+matrix, m = matrix only,\n");
	fprintf(stderr,"                 G = single gamma+matrix, S = single shaper+matrix\n");
//  Development - not supported
//	fprintf(stderr," -I ver          Set ICC profile version > 2.2.0\n");
//	fprintf(stderr,"                 ver = 4, Enable ICC V4 creation\n");
	fprintf(stderr," -u              If Lut input profile, make it absolute (non-standard)\n");
	fprintf(stderr," -U scale        If input profile, scale media white point by scale\n");
	fprintf(stderr," -i illum        Choose illuminant for print/transparency spectral data:\n");
	fprintf(stderr,"                 A, C, D50 (def.), D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ       Choose CIE Observer for spectral data:\n");
	fprintf(stderr,"                 1931_2 (def), 1964_10, S&B 1955_2, shaw, J&V 1978_2\n");
	fprintf(stderr," -f              Use Fluorescent Whitening Agent compensation\n");
	fprintf(stderr," -r avgdev       Average deviation of device+instrument readings as a percentage (default %4.2f%%)\n",DEFAVGDEV);
/* Research options: */
/*	fprintf(stderr," -r sSMOOTH      RSPL suplimental optimised smoothing factor\n"); */
/*	fprintf(stderr," -r rSMOOTH      RSPL raw underlying smoothing factor\n"); */
	fprintf(stderr," -s src%s      Apply gamut mapping to output profile perceptual B2A table for given source space\n",ICC_FILE_EXT);
	fprintf(stderr," -S src%s      Apply gamut mapping to output profile perceptual and saturation B2A table\n",ICC_FILE_EXT);
	fprintf(stderr," -nP             Use colormetric source gamut to make output profile perceptual table\n");
	fprintf(stderr," -nS             Use colormetric source gamut to make output profile saturation table\n");
	fprintf(stderr," -g src.gam      Use source image gamut as well for output profile gamut mapping\n");
	fprintf(stderr," -p absprof      Incorporate abstract profile into output tables\n");
	fprintf(stderr," -t intent       Override gamut mapping intent for output profile perceptual table:\n");
	fprintf(stderr," -T intent       Override gamut mapping intent for output profile saturation table:\n");
	for (i = 0; ; i++) {
		icxGMappingIntent gmi;
		if (xicc_enum_gmapintent(&gmi, i, NULL) == icxIllegalGMIntent)
			break;
		fprintf(stderr,"              %s\n",gmi.desc);
	}
	fprintf(stderr," -c viewcond     set input viewing conditions for output profile %s gamut mapping,\n",icxcam_description(cam_default));
	fprintf(stderr,"                  either an enumerated choice, or a parameter\n");
	fprintf(stderr," -d viewcond     set output viewing conditions for output profile %s gamut mapping\n",icxcam_description(cam_default));
	fprintf(stderr,"                  either an enumerated choice, or a parameter\n");
	fprintf(stderr,"                  Also sets out of gamut clipping CAM space.\n");
	fprintf(stderr,"                  either an enumerated choice, or a series of parameters:value changes\n");
	for (i = 0; ; i++) {
		icxViewCond vc;
		if (xicc_enum_viewcond(NULL, &vc, i, NULL, 1, NULL) == -999)
			break;

		fprintf(stderr,"             %s\n",vc.desc);
	}
	fprintf(stderr," -P              Create gamut gammap_p.wrl and gammap_s.wrl diagostics\n");
	fprintf(stderr," -O outputfile   Override the default output filename.\n");
	fprintf(stderr," inoutfile       Base name for input.ti3/output%s file\n",ICC_FILE_EXT);
	exit(1);
}


int main(int argc, char *argv[]) {
	int fa,nfa,mfa;				/* current argument we're looking at */
#ifdef DO_TIME			/* Time the operation */
	clock_t stime, ttime;		/* Start and total times */
#endif
	int verb = 0;
	int iquality = 1;			/* A2B quality */
	int oquality = -1;			/* B2A quality same as A2B */
	int verify = 0;
	int noisluts = 0;			/* No input shaper luts */
	int noipluts = 0;			/* No input position luts */
	int nooluts = 0;			/* No output shaper luts */
	int nocied = 0;				/* No .ti3 CIE data in profile */
	int noptop = 0;				/* Use colormetric source gamut to make perceptual table */
	int nostos = 0;				/* Use colormetric source gamut to make saturation table */
	int gamdiag = 0;			/* Make gamut mapping diagnostic wrl plots */
	int nsabs = 0;				/* Make non-standard absolute input lut profile */
	double iwpscale = -1.0;		/* Input white point scale factor */
	int doinb2a = 1;			/* Create an input device B2A table */
	int inking = 3;				/* Default K target ramp K */
	int locus = 0;				/* Default K value target */
	double Kstle = 0.0, Kstpo = 0.0, Kenle = 0.0, Kenpo = 0.0, Kshap = 0.0;
	int tlimit = -1;			/* Total ink limit as a % */
	int klimit = -1;			/* Black ink limit as a % */
	int fwacomp = 0;			/* FWA compensation */
	double avgdev = DEFAVGDEV/100.0;	/* Average measurement deviation */
	double smooth = 1.0;		/* RSPL Smoothness factor (relative, for verification) */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType illum = icxIT_D50;	/* Spectral defaults */
	xspect cust_illum;			/* Custom illumination spectrum */
	icxObserverType observ = icxOT_CIE_1931_2;	/* The classic observer */
	char ipname[MAXNAMEL+1] = "";	/* Input icc profile - enables gamut map */
	char sgname[MAXNAMEL+1] = "";	/* Image source gamut name */
	char absname[MAXNAMEL+1] = "";	/* Abstract profile name */
	int sepsat = 0;				/* Create separate saturation B2A table */
	icxViewCond ivc_p;			/* Input Viewing Parameters for CAM */
	icxViewCond ovc_p;			/* Output Viewing Parameters for CAM (enables CAM clip) */
	int ivc_e = -1, ovc_e = -1;	/* Enumerated viewing condition */
	icxGMappingIntent pgmi;		/* default Perceptual gamut mapping intent */
	int pgmi_set = 0;			/* Set by user option */
	icxGMappingIntent sgmi;		/* default Saturation gamut mapping intent */
	int sgmi_set = 0;			/* Set by user option */
	char baname[MAXNAMEL+1] = "";	/* Input & Output base name */
	char inname[MAXNAMEL+1] = "";	/* Input cgats file base name */
	char outname[MAXNAMEL+1] = "";	/* Output cgats file base name */
	cgats *icg;					/* input cgats structure */
	int ti;						/* Temporary CGATs index */
	prof_atype ptype = prof_default;	/* Default for each type of device */
	int mtxtoo = 0;				/* NZ if matrix tags should be created for Display XYZ cLUT */
	icmICCVersion iccver = icmVersionDefault;	/* ICC profile version to create */
	profxinf xpi;		/* Extra profile information */

	
#ifdef DO_TIME			/* Time the operation */
	stime = clock();
#endif /* DO_TIME */
	error_program = argv[0];
	memset((void *)&xpi, 0, sizeof(profxinf));	/* Init extra profile info to defaults */

	/* Init VC overrides so that we know when the've been set */
	ivc_p.Ev = -1;
	ivc_p.Wxyz[0] = -1.0; ivc_p.Wxyz[1] = -1.0; ivc_p.Wxyz[2] = -1.0;
	ivc_p.La = -1.0;
	ivc_p.Yb = -1.0;
	ivc_p.Lv = -1.0;
	ivc_p.Yf = -1.0;
	ivc_p.Fxyz[0] = -1.0; ivc_p.Fxyz[1] = -1.0; ivc_p.Fxyz[2] = -1.0;

	ovc_p.Ev = -1;
	ovc_p.Wxyz[0] = -1.0; ovc_p.Wxyz[1] = -1.0; ovc_p.Wxyz[2] = -1.0;
	ovc_p.La = -1.0;
	ovc_p.Yb = -1.0;
	ovc_p.Lv = -1.0;
	ovc_p.Yf = -1.0;
	ovc_p.Fxyz[0] = -1.0; ovc_p.Fxyz[1] = -1.0; ovc_p.Fxyz[2] = -1.0;

	xicc_enum_gmapintent(&pgmi, icxPerceptualGMIntent, NULL);
	xicc_enum_gmapintent(&sgmi, icxSaturationGMIntent, NULL);

	if (argc <= 1)
		usage("Too few arguments, got %d expect at least %d",argc-1,1);

	/* Process the arguments */
	mfa = 1;		/* Minimum final arguments */
	for(fa = 1;fa < argc;fa++) {

		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1+mfa) < argc) {
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
			else if (argv[fa][1] == 'D') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to profile description flag -D");
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
			else if (argv[fa][1] == 'b') {
				if (na != NULL) {	/* Got a B2A quaiity */
					fa = nfa;
	    			switch (na[0]) {
						case 'l':
						case 'L':
							oquality = 0;
							break;
						case 'm':
						case 'M':
							oquality = 1;
							break;
						case 'h':
						case 'H':
							oquality = 2;
							break;
						case 'u':
						case 'U':
							oquality = 3;
							break;
						case 'n':				/* No B2A for input device */
						case 'N':
							oquality = -2;
							doinb2a = 0;
							break;
						default:
							usage("Unknown argument '%c' to quality flag -q",na[0]);
					}
				} else
					oquality = 0;
			}

			else if (argv[fa][1] == 'B') {
				oquality = -2;
				doinb2a = 0;
			}

			else if (argv[fa][1] == 'y' || argv[fa][1] == 'Y')
				verify = 1;

			/* Disable input or output luts */
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N') {
				fa = nfa;
				if (na == NULL) {	/* Backwards compatible */
					nooluts = 1;
				} else {
					if (na[0] == 'i')
						noisluts = 1;
					else if (na[0] == 'p')
						noipluts = 1;
					else if (na[0] == 'o')
						nooluts = 1;
					else if (na[0] == 'c')
						nocied = 1;
					else if (na[0] == 'P')
						noptop = 1;
					else if (na[0] == 'S')
						nostos = 1;
					else
						usage("Unknown argument '%c' to flag -n",na[0]);
				}
			}

			else if (argv[fa][1] == 'u') {
				nsabs = 1;
			}
			else if (argv[fa][1] == 'U') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to input white point scale flag -U");
				iwpscale = atof(na);
				if (iwpscale < 0.0 || iwpscale > 100.0)
					usage("Argument '%s' to flag -U out of range",na);
			}

			/* Inking rule */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to inking flag -k");
				if (argv[fa][1] == 'k')
					locus = 0;		/* Use K value target */
				else
					locus = 1;		/* Use K locus target */
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
						inking = 2;		/* Use maximum K */
						break;
					case 'r':
					case 'R':
						inking = 3;		/* Use ramping K */
						break;
					case 'p':
					case 'P':
						inking = 4;		/* Use parameter curve */
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
					case 'X':
						mtxtoo = 1;
						/* Fall though */
					case 'x':
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
					case 'm':
						ptype = prof_matonly;
						break;
					default:
						usage("Unknown argument '%c' to algorithm flag -a",na[0] );
				}
			}
			/* Profile version */
			else if (argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to version flag -I");
    			switch (na[0]) {
					case '4':
						iccver = icmVersion4_1;
						break;
					default:
						usage("Unknown argument '%c' to version flag -I",na[0] );
				}
			}
			/* Spectral Illuminant type */
			else if (argv[fa][1] == 'i') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to illuminant flag -i");
				if (strcmp(na, "A") == 0) {
					spec = 1;
					illum = icxIT_A;
				} else if (strcmp(na, "C") == 0) {
					spec = 1;
					illum = icxIT_C;
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
			else if (argv[fa][1] == 'o') {
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
				if (na[0] == 's') {			/* (relative, for verification) */
					smooth = atof(na+1);
					if (smooth < 0.0)
						usage("Optimised smoothing factor argument to '-rs' must be over 0.0");
				} else if (na[0] == 'r') {	/* (absolute, for testing) */
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
			else if (argv[fa][1] == 'p') {
				if (na == NULL) usage("Expected abstract profile filename after -p");
				fa = nfa;
				strncpy(absname,na,MAXNAMEL); absname[MAXNAMEL] = '\000';
			}

			/* Perceptual Mapping intent override */
			else if (argv[fa][1] == 't') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to perceptul intent override flag -t");
				if (xicc_enum_gmapintent(&pgmi, icxNoGMIntent, na) == icxIllegalGMIntent)
					usage("Unrecognised intent '%s' to perceptual override flag -t",na);
				pgmi_set = 1;
			}

			/* Saturation Mapping intent override */
			else if (argv[fa][1] == 'T') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to saturation intent override flag -T");
				if (xicc_enum_gmapintent(&sgmi, icxNoGMIntent, na) == icxIllegalGMIntent)
					usage("Unrecognised intent '%s' to saturation override flag -T",na);
				sgmi_set = 1;
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
				if (na[1] != ':') {
					if (vc == &ivc_p) {
						if ((ivc_e = xicc_enum_viewcond(NULL, NULL, -2, na, 1, NULL)) == -999)
							usage("Urecognised Enumerated Viewing conditions '%s'",na);
					} else {
						if ((ovc_e = xicc_enum_viewcond(NULL, NULL, -2, na, 1, NULL)) == -999)
							usage("Urecognised Enumerated Viewing conditions '%s'",na);
					}
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

			/* Gammut mapping diagnostic plots */
			else if (argv[fa][1] == 'P')
				gamdiag = 1;

			/* Output file name */
			else if (argv[fa][1] == 'O') {
				fa = nfa;
				if (na == NULL) usage("Output filename override (-O) needs an argument");
				strncpy(outname,na,MAXNAMEL); outname[MAXNAMEL] = '\000';
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
	strcat(inname,".ti3");
	if (outname[0] == '\000') {		/* If not overridden */
		strcpy(outname,baname);
		strcat(outname,ICC_FILE_EXT);
	}

	/* Issue some errors & warnings for strange combinations */
	if (fwacomp && spec == 0)
		error("FWA compensation only works when viewer and/or illuminant selected");

	if (pgmi_set && ipname[0] == '\000')
		warning("-t perceptual intent override only works if -s srcprof or -S srcprof is used");

	if (sgmi_set && ipname[0] == '\000')
		warning("-T saturation intent override only works if -S srcprof is used");

	if (sgmi_set && sepsat == 0) {	/* Won't do much otherwise */
		if (verb)
			printf("Saturation intent override was set, so adding saturation intent table\n");
		sepsat = 1;
	}

	if (sgname[0] != '\000' && ipname[0] == '\000')
		warning("-g srcgam will do nothing without -s srcprof or -S srcprof");

	if (oquality == -1) {		/* B2A tables will be used */
		oquality = iquality;
	}

	/* Open and look at the .ti3 profile patches file */
	icg = new_cgats();			/* Create a CGATS structure */
	icg->add_other(icg, "CTI3"); 	/* our special input type is Calibration Target Information 3 */
	icg->add_other(icg, "CAL"); 	/* our special device Calibration state */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->ntables == 0 || icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI3 format file");
	if (icg->ntables < 1)
		error ("Input file doesn't contain at least one table");

	/* See if CIE is actually available - some sources of .TI3 don't provide it */
	if (!spec
	 && icg->find_field(icg, 0, "LAB_L") < 0
	 && icg->find_field(icg, 0, "XYZ_X") < 0) {

		if (icg->find_kword(icg, 0, "SPECTRAL_BANDS") < 0)
			error("Neither CIE nor spectral data found in file '%s'",inname);

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
			int imax;
			imax = atoi(icg->t[0].kdata[ti]);
			if (imax > 0 && imax <= 400.0) {
				if (tlimit > 0 && tlimit <= 400.0) {	/* User has specified limit as option */
					if (imax < tlimit) {
						warning("Ink limit greater than original chart! (%d%% > %d%%)",tlimit,imax);
					}
				} else {
					if (imax > 80.0)
						tlimit = (int)imax - 10;	/* Rule of thumb - 10% below chart maximum */
					else
						tlimit = (int)imax;
				}
			}
		}

		/* (Note that this isn't set by any of the Argyll tools currently, */
		/*  but can be set manually.) */
		if ((ti = icg->find_kword(icg, 0, "BLACK_INK_LIMIT")) >= 0) {
			int kmax;
			kmax = atoi(icg->t[0].kdata[ti]);
			if (kmax > 0 && kmax <= 100.0) {
				if (klimit > 0 && klimit <= 100.0) {	/* User has specified limit as option */
					if (kmax < klimit) {
						warning("Black ink limit greater than original chart! (%d%% > %d%%)",klimit,kmax);
					}
				} else {
					klimit = (int)kmax;
				}
			}
		}

		if (tlimit >= 0 && tlimit < 400.0) {
			if (verb)
				printf("Total ink limit being used is %d%%\n",tlimit);
			ink.tlimit = tlimit/100.0;	/* Set a total ink limit */
		} else {
			if (verb)
				printf("No total ink limit being used\n");
			ink.tlimit = -1.0;			/* Don't use a limit */
		}

		if (klimit >= 0 && klimit < 100.0) {
			if (verb)
				printf("Black ink limit being used is %d%%\n",klimit);
			ink.klimit = klimit/100.0;	/* Set a black ink limit */
		} else {
			if (verb)
				printf("No black ink limit being used\n");
			ink.klimit = -1.0;			/* Don't use a limit */
		}

		ink.KonlyLmin = 0;				/* Use normal black Lmin for locus */
		ink.c.Ksmth = ICXINKDEFSMTH;	/* default black curve smoothing */
		ink.c.Kskew = ICXINKDEFSKEW;	/* default black curve skew */
		ink.x.Ksmth = ICXINKDEFSMTH;
		ink.x.Kskew = ICXINKDEFSKEW;

		if (inking == 0) {			/* Use minimum */
			ink.k_rule = locus ? icxKluma5 : icxKluma5k;
			ink.c.Kstle = 0.0;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 0.0;
			ink.c.Kshap = 1.0;
		} else if (inking == 1) {	/* Use 0.5  */
			ink.k_rule = locus ? icxKluma5 : icxKluma5k;
			ink.c.Kstle = 0.5;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 0.5;
			ink.c.Kshap = 1.0;
		} else if (inking == 2) {	/* Use maximum  */
			ink.k_rule = locus ? icxKluma5 : icxKluma5k;
			ink.c.Kstle = 1.0;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 1.0;
			ink.c.Kshap = 1.0;
		} else if (inking == 3) {	/* Use ramp  */
			ink.k_rule = locus ? icxKluma5 : icxKluma5k;
			ink.c.Kstle = 0.0;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 1.0;
			ink.c.Kshap = 1.0;
		} else {				/* Use specified curve */
			ink.k_rule = locus ? icxKluma5 : icxKluma5k;
			ink.c.Kstle = Kstle;
			ink.c.Kstpo = Kstpo;
			ink.c.Kenpo = Kenpo;
			ink.c.Kenle = Kenle;
			ink.c.Kshap = Kshap;
		}

		if (ptype == prof_default)
			ptype = prof_clutLab;
		else if (ptype != prof_clutLab && ptype != prof_clutXYZ) {
			error ("Output profile can only be a cLUT algorithm");
		}

		make_output_icc(ptype, 0, iccver, verb, iquality, oquality,
		                noisluts, noipluts, nooluts, nocied, noptop, nostos,
		                gamdiag, verify, &ink, inname, outname, icg, spec,
		                illum, &cust_illum, observ, fwacomp, smooth, avgdev,
		                ipname[0] != '\000' ? ipname : NULL,
		                sgname[0] != '\000' ? sgname : NULL,
		                absname[0] != '\000' ? absname : NULL,
						sepsat, &ivc_p, &ovc_p, ivc_e, ovc_e,
						&pgmi, &sgmi, &xpi);

	} else if (strcmp(icg->t[0].kdata[ti],"INPUT") == 0) {

		if (ptype == prof_default)
			ptype = prof_clutLab;		/* For best possible quality */
		make_input_icc(ptype, iccver, verb, iquality, oquality, noisluts, noipluts, nooluts, nocied,
		               verify, nsabs, iwpscale, doinb2a, inname, outname, icg,
		               spec, illum, &cust_illum, observ, smooth, avgdev, &xpi);

	} else if (strcmp(icg->t[0].kdata[ti],"DISPLAY") == 0) {

		if (fwacomp)
			error ("FWA compensation isn't applicable to a display device");

		if (ptype == prof_default)
			ptype = prof_clutLab;		/* ?? or should it default to prof_shamat ?? */

		/* If a source gamut is provided for a Display, then a V2.4.0 profile will be created */
		make_output_icc(ptype, mtxtoo, iccver, verb, iquality, oquality,
		                noisluts, noipluts, nooluts, nocied, noptop, nostos,
		                gamdiag, verify, NULL, inname, outname, icg, spec,
		                illum, &cust_illum, observ, 0, smooth, avgdev,
		                ipname[0] != '\000' ? ipname : NULL,
		                sgname[0] != '\000' ? sgname : NULL,
		                absname[0] != '\000' ? absname : NULL,
						sepsat, &ivc_p, &ovc_p, ivc_e, ovc_e,
						&pgmi, &sgmi, &xpi);

	} else
		error ("Input file keyword DEVICE_CLASS has unknown value");

	icg->del(icg);		/* Clean up */

#ifdef DO_TIME			/* Time the operation */
	ttime = clock() - stime;
	printf("Exectution time = %f seconds\n",(double)ttime/(double)CLOCKS_PER_SEC);
#endif /* DO_TIME */

	return 0;
}










