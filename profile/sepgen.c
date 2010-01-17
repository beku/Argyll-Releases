
/* 
 * Argyll Color Correction System
 * Optimised Seperation Generator
 *
 * Author: Graeme W. Gill
 * Date:   2002/10/11
 *
 * Copyright 2002 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes in a device ICC or MPP profile,
 * and creates an optimsed separation from either
 * RGB'/CMY' or CMYK' pseudo-device space to the real
 * device space.
 */

/*
 * TTBD:
 *      Again would like ink limit encoded in icc profile, or estimated from
 *      B2A table.
 */

#define DEBUG

#define verbo stdout

#include <stdio.h>
#include "copyright.h"
#include "config.h"
#include "cgats.h"
#include "numlib.h"
#include "xicc.h"
#include "prof.h"

void usage(char *diag, ...) {
	fprintf(stderr,"Create Optimsed separation, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: %s [-options] infile.[icm|mpp] outprofile.icm\n",error_program);
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -t              Create psuedo-CMY input device space,\n");
	fprintf(stderr," -T              Create psuedo-RGB input device space,\n");
	fprintf(stderr,"                 Default is psuedo-CMYK input device space,\n");
	fprintf(stderr," -e description  Description string\n");
	fprintf(stderr," -q [lmhu]       Quality - Low, Medium (def), High, Ultra\n");
	fprintf(stderr," -l tlimit       override total ink limit, 0 - 400+%%\n");
	fprintf(stderr," -L klimit       override black ink limit, 0 - 100%%\n");
	fprintf(stderr,"                 Black generation for pseudo-CMY or RGB:\n");
	fprintf(stderr," -k [zhxr]       z = zero K, h = 0.5 K,\n");
	fprintf(stderr,"                 x = max K, r = ramp K (def.)\n");
	fprintf(stderr," -k p stle stpo enpo enle shape\n");
	fprintf(stderr,"                 stle: K level at White 0.0 - 1.0\n");
	fprintf(stderr,"                 stpo: start point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"                 enpo: End point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"                 enle: K level at Black 0.0 - 1.0\n");
	fprintf(stderr,"                 shape: 1.0 = straight, 0.0-1.0 concave, 1.0-2.0 convex\n");
	fprintf(stderr," -K parameters   Same as -k, but target is K locus rather than K value itself\n");
	fprintf(stderr," -i illum        Choose illuminant for print/transparency spectral MPP profile:\n");
	fprintf(stderr,"                 A, C, D50 (def.), D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ       Choose CIE Observer for spectral MPP profile:\n");
	fprintf(stderr,"                 1931_2 (def), 1964_10, S&B 1955_2, shaw, J&V 1978_2\n");
	fprintf(stderr," -f              Use Fluorescent Whitening Agent compensation for spectral MPP profile\n");
	fprintf(stderr," infile.[icm|mpp] Name of device forward profile\n");
	fprintf(stderr," outfile.icm      Name for output.icm device to device separation profile file\n");
	exit(1);
	}


int main(int argc, char *argv[]) {
	int fa,nfa;					/* current argument we're looking at */
	int verb = 0;
	int quality = 1;			/* quality */
	int inking = 3;				/* Default ramp */
	int locus = 0;				/* Default K value target */
	double Kstle = 0.0, Kstpo = 0.0, Kenle = 0.0, Kenpo = 0.0, Kshap = 0.0;
	double tlimit = -1.0;		/* No total ink limit (0 .. inn)*/
	double klimit = -1.0;		/* No black ink limit (0 .. 1.0)*/
	int fwacomp = 0;			/* FWA compensation */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType illum = icxIT_D50;	/* Spectral defaults */
	xspect cust_illum;			/* Custom illumination spectrum */
	icxObserverType observ = icxOT_CIE_1931_2;
	static char inname[1000] = { 0 };		/* Input icc or mpp file base name */
	static char outname[1000] = { 0 };		/* Output icc file base name */
	profxinf xpi;				/* Extra profile information */

	icmFile *fp;				/* Input icc file */
	icc *icco = NULL;			/* Input icc */
	icmLuBase *luo = NULL;		/* Icc lookup object */
	mpp *mppo = NULL;			/* Input mpp */
	
	int inn = 4;				/* Input number of components */
	int iimask = ICX_CMYK;		/* Input ink mask */
	int outn;					/* Output number of components */
	inkmask oimask;				/* Output ink mask */
//	xsep *xsepo;				/* Separation object */
//	icc *sep_icco;				/* Output icc */

	icxInk ink;					/* K Inking parameters */

	int rv;

	error_program = argv[0];
	memset((void *)&xpi, 0, sizeof(profxinf));	/* Init extra profile info to defaults */

	if (argc < 3)
		usage("Too few arguments, got %d expect at least %d",argc-1,2);

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

			/* CMY rather than CMYK pseudo space */
			else if (argv[fa][1] == 't') {
				inn = 3;
				iimask = ICX_CMY;
			}

			/* RGB rather than CMYK pseudo space */
			else if (argv[fa][1] == 'T') {
				inn = 3;
				iimask = ICX_RGB;
			}

			/* Description */
			else if (argv[fa][1] == 'e' || argv[fa][1] == 'E') {
				fa = nfa;
				if (na == NULL)usage("Expect argument to description flag -e");
				xpi.profDesc = na;
			}

			/* Quality */
			else if (argv[fa][1] == 'q' || argv[fa][1] == 'Q') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to quality flag -q");
    			switch (na[0]) {
					case 'l':
					case 'L':
						quality = 0;
						break;
					case 'm':
					case 'M':
						quality = 1;
						break;
					case 'h':
					case 'H':
						quality = 2;
						break;
					case 'u':
					case 'U':
						quality = 3;
						break;
					default:
						usage("Unknown argument '%c' to quality flag -q",na[0]);
				}
			}

			/* Inking rule */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to inking flag -k");
				if (argv[fa][1] == 'k')
					locus = 0;			/* Use K value target */
				else
					locus = 1;			/* Use K locus target */
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
			else if (argv[fa][1] == 'l') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to total ink limit flag -l");
				tlimit = atof(na)/100.0;
			}
			else if (argv[fa][1] == 'L') {
				fa = nfa;
				if (na == NULL) usage("Expect argument to black ink limit flag -l");
				klimit = atof(na)/100.0;
			}

			/* Spectral Illuminant type */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
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

			else 
				usage("Unknown flag '%c'",argv[fa][1]);
		} else
			break;
	}

	/* Get the file name arguments */
	if (fa >= argc || argv[fa][0] == '-') usage("Missing input .icm or .mpp filename");
	strcpy(inname,argv[fa]);

	if (fa >= argc || argv[fa][0] == '-') usage("Missing output .icm filename");
	strcpy(outname,argv[fa]);

	if (fwacomp && spec == 0)
		error ("FWA compensation only works when viewer and/or illuminant selected");

	/* See if we have an icc input profile */
	if ((fp = new_icmFileStd_name(inname,"r")) == NULL)
		error ("Can't open file '%s'",inname);

	if ((icco = new_icc()) == NULL)
		error ("Creation of ICC object failed");

	if ((rv = icco->read(icco,fp,0)) == 0) {	/* It is an ICC */
		icColorSpaceSignature outs;

		/* Get a conversion object */
		if ((luo = icco->get_luobj(icco, icmFwd, icAbsoluteColorimetric,
		                           icSigXYZData, icmLuOrdNorm)) == NULL)
			error ("%d, %s",icco->errc, icco->err);
	
		/* Get details of conversion (Arguments may be NULL if info not needed) */
		luo->spaces(luo, &outs, &outn, NULL, NULL, NULL, NULL, NULL, NULL, NULL);

		switch(outs) {
			case icSigCmykData:
				oimask = ICX_CMYK;
				break;

			default:
				error ("Can't handle device space '%s' for separation yet",
				       icm2str(icmColorSpaceSignature, outs));
		}

		if (spec)
			error ("Requested spectral interpretation when data not available");

	} else {	/* Not an ICC, see if it's an MPP */
		int spec_n;
		double plimit;

		fp->del(fp);
		fp = NULL;
		icco->del(icco);
		icco = NULL;

		if ((mppo = new_mpp()) == NULL)
			error ("Creation of MPP object failed");

		if (mppo->read_mpp(mppo, inname))
			error("Read error : %s",mppo->err);

		/* Get various types of information about the mpp */
		mppo->get_info(mppo, &oimask, &outn, &plimit, &spec_n, NULL, NULL, NULL);

		if (spec && spec_n == 0)
			error ("Requested spectral interpretation when data not available");

		if (tlimit < 1e-4 || plimit < tlimit)	/* Set to device limit */
			tlimit = plimit;
	}


	if (verb) {
		char *ident = icx_inkmask2char(oimask, 1);
		printf("Creating separation from pseudo-%s to %s device spaces\n",
		        iimask == ICX_CMY ? "CMY" : iimask == ICX_RGB || iimask == ICX_IRGB ? "RGB" : "CMYK", ident);
		free(ident);
	}

	/* Configure ink limits */
	if (tlimit >= 0.0) {
		if (verb)
			printf("Total ink limit being used is %.0f%%\n",tlimit);
		ink.tlimit = tlimit;			/* Set a total ink limit */
	} else {
		if (verb)
			printf("No total ink limit being used\n");
		ink.tlimit = -1.0;			/* Don't use a limit */
	}

	if (klimit >= 0.0) {
		if (verb)
			printf("Black ink limit being used is %.0f%%\n",klimit);
		ink.klimit = klimit;		/* Set a black ink limit */
	} else {
		if (verb)
			printf("No black ink limit being used\n");
		ink.klimit = -1.0;			/* Don't use a black limit */
	}

	/* Configure black generation */
	/* This will be ignored if the pseudo input is CMYK */

	ink.KonlyLmin = 0;				/* Use normal black Lmin for locus */
	ink.c.Ksmth = ICXINKDEFSMTH;	/* Default curve smoothing */
	ink.c.Kskew = ICXINKDEFSKEW;	/* default curve skew */
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


	/* Generate the separation */
	// ~~999

	return 0;
}










