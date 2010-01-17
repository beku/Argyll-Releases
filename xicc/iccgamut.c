
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
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */


/*
 * TTBD:
 *       To support CIACAM02 properly, need to cope with viewing parameters ?
 */

#define SURFACE_ONLY
#define GAMRES 10.0		/* Default surface resolution */

#define USE_CAM_CLIP_OPT		/* Use CAM space to clip in */

#define RGBRES 33	/* 33 */
#define CMYKRES 17	/* 17 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "icc.h"
#include "xicc.h"
#include "gamut.h"
#include "counters.h"

static void diag_gamut(icxLuBase *p, double detail, int doaxes,
                       double tlimit, double klimit, char *outname);

void usage(char *diag) {
	int i;
	fprintf(stderr,"Create Lab/Jab gamut plot Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	fprintf(stderr,"usage: iccgamut [options] profile\n");
	if (diag != NULL)
		fprintf(stderr,"Diagnostic: %s\n",diag);
	fprintf(stderr," -v            Verbose\n");
	fprintf(stderr," -d sres       Surface resolution details 1.0 - 50.0\n");
	fprintf(stderr," -w            emit VRML .wrl file as well as CGATS .gam file\n");
	fprintf(stderr," -n            Don't add VRML axes or white/black point\n");
	fprintf(stderr," -k            Add VRML markers for prim. & sec. \"cusp\" points\n");
	fprintf(stderr," -f function   f = forward*, b = backwards\n");
	fprintf(stderr," -i intent     p = perceptual, r = relative colorimetric,\n");
	fprintf(stderr,"               s = saturation, a = absolute (default), d = profile default\n");
//  fprintf(stderr,"               P = absolute perceptual, S = absolute saturation\n");
	fprintf(stderr," -p oride      l = Lab_PCS (default), j = %s Appearance Jab\n",icxcam_description(cam_default),icxcam_description(cam_default));
	fprintf(stderr," -o order      n = normal (priority: lut > matrix > monochrome)\n");
	fprintf(stderr,"               r = reverse (priority: monochrome > matrix > lut)\n");
	fprintf(stderr," -l tlimit     set total ink limit, 0 - 400%% (estimate by default)\n");
	fprintf(stderr," -L klimit     set black ink limit, 0 - 100%% (estimate by default)\n");
	fprintf(stderr," -c viewcond   set viewing conditions for %s,\n",icxcam_description(cam_default));
	fprintf(stderr,"               either an enumerated choice, or a series of parameter:value changes\n");
	for (i = 0; ; i++) {
		icxViewCond vc;
		if (xicc_enum_viewcond(NULL, &vc, i, NULL, 1, NULL) == -999)
			break;

		fprintf(stderr,"           %s\n",vc.desc);
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
    fprintf(stderr," -s                    Create special cube surface topology plot\n");
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
	int docusps = 0;
	double gamres = GAMRES;		/* Surface resolution */
	int special = 0;			/* Special surface plot */
	int fl = 0;					/* luobj flags */
	icxInk ink;					/* Ink parameters */
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

	icxLuBase *luo;

	/* Lookup parameters */
	icmLookupFunc     func   = icmFwd;				/* Default */
	icRenderingIntent intent = -1;					/* Default */
	icColorSpaceSignature pcsor = icSigLabData;		/* Default */
	icmLookupOrder    order  = icmLuOrdNorm;		/* Default */
	
	error_program = argv[0];

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
					case 'd':
						intent = icmDefaultIntent;
						break;
					case 'a':
						intent = icAbsoluteColorimetric;
						break;
					case 'p':
						intent = icPerceptual;
						break;
					case 'r':
						intent = icRelativeColorimetric;
						break;
					case 's':
						intent = icSaturation;
						break;
					/* Argyll special intents to check spaces underlying */
					/* icxPerceptualAppearance & icxSaturationAppearance */
					case 'P':
						intent = icmAbsolutePerceptual;
						break;
					case 'S':
						intent = icmAbsoluteSaturation;
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

			/* PCS override */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -p");
    			switch (na[0]) {
					case 'l':
						pcsor = icSigLabData;
						break;
					case 'j':
						pcsor = icxSigJabData;
						break;
					default:
						usage("Unrecognised parameter after flag -p");
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
			/* No axis output in vrml */
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N') {
				doaxes = 0;
			}
			/* Do cusp markers in vrml */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				docusps = 1;
			}
			/* Special */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				special = 1;
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
#ifdef NEVER
				if (na[0] >= '0' && na[0] <= '9') {
					vc_e = atoi(na);
				} else
#endif
				if (na[1] != ':') {
					if ((vc_e = xicc_enum_viewcond(NULL, NULL, -2, na, 1, NULL)) == -999)
						usage("Urecognised Enumerated Viewing conditions");
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

	if (intent == -1) {
		if (pcsor == icxSigJabData)
			intent = icRelativeColorimetric;	/* Default to icxAppearance */
		else
			intent = icAbsoluteColorimetric;	/* Default to icAbsoluteColorimetric */
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

	/* Set the ink limits */
	icxDefaultLimits(xicco, &ink.tlimit, tlimit/100.0, &ink.klimit, klimit/100.0);

	if (verb) {
		if (ink.tlimit >= 0.0)
			printf("Total ink limit assumed is %3.0f%%\n",100.0 * ink.tlimit);
		if (ink.klimit >= 0.0)
			printf("Black ink limit assumed is %3.0f%%\n",100.0 * ink.klimit);
	}

	/* Setup a safe ink generation (not used) */
	ink.KonlyLmin = 0;		/* Use normal black Lmin for locus */
	ink.k_rule = icxKluma5k;
	ink.c.Ksmth = ICXINKDEFSMTH;	/* Default smoothing */
	ink.c.Kskew = ICXINKDEFSKEW;	/* default curve skew */
	ink.c.Kstle = 0.0;		/* Min K at white end */
	ink.c.Kstpo = 0.0;		/* Start of transition is at white */
	ink.c.Kenle = 1.0;		/* Max K at black end */
	ink.c.Kenpo = 1.0;		/* End transition at black */
	ink.c.Kshap = 1.0;		/* Linear transition */

	/* Setup the default viewing conditions */
	if (xicc_enum_viewcond(xicco, &vc, -1, NULL, 0, NULL) == -2)
		error ("%d, %s",xicco->errc, xicco->err);

	if (vc_e != -1)
		if (xicc_enum_viewcond(xicco, &vc, vc_e, NULL, 0, NULL) == -2)
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

	fl |= ICX_CLIP_NEAREST;		/* Don't setup rev uncessarily */

#ifdef USE_CAM_CLIP_OPT
	 fl |= ICX_CAM_CLIP;
#endif

#ifdef NEVER
	printf("~1 output space flags = 0x%x\n",fl);
	printf("~1 output space intent = %s\n",icx2str(icmRenderingIntent,intent));
	printf("~1 output space pcs = %s\n",icx2str(icmColorSpaceSignature,pcsor));
	printf("~1 output space viewing conditions =\n"); xicc_dump_viewcond(&vc);
	printf("~1 output space inking =\n"); xicc_dump_inking(&ink);
#endif

	strcpy(out_name, prof_name);
	if ((xl = strrchr(out_name, '.')) == NULL)	/* Figure where extention is */
		xl = out_name + strlen(out_name);

	strcpy(xl,".gam");

	/* Get a expanded color conversion object */
	if ((luo = xicco->get_luobj(xicco, fl, func, intent, pcsor, order, &vc, &ink)) == NULL)
		error ("%d, %s",xicco->errc, xicco->err);

	if (special) {
		if (func != icmFwd)
			error("Must be forward direction for special plot");
		strcpy(xl,".wrl");
		diag_gamut(luo, gamres, doaxes, tlimit/100.0, klimit/100.0, out_name); 
	} else {
		/* Creat a gamut surface */
		if ((gam = luo->get_gamut(luo, gamres)) == NULL)
			error ("%d, %s",xicco->errc, xicco->err);

		if (gam->write_gam(gam,out_name))
			error ("write gamut failed on '%s'",out_name);

		if (vrml) {
			strcpy(xl,".wrl");
			if (gam->write_vrml(gam,out_name, doaxes, docusps))
				error ("write vrml failed on '%s'",out_name);
		}

		if (verb) {
			printf("Total volume of gamut is %f cubic colorspace units\n",gam->volume(gam));
		}
		gam->del(gam);
	}

	luo->del(luo);			/* Done with lookup object */

	xicco->del(xicco);		/* Expansion wrapper */
	icco->del(icco);		/* Icc */
	fp->del(fp);


	return 0;
}

/* -------------------------------------------- */
/* Code for special gamut surface plot */

#define GAMUT_LCENT 50

/* Create a diagnostic gamut, illustrating */
/* device space "fold-over" */
static void diag_gamut(
icxLuBase *p,		/* Lookup object */
double detail,		/* Gamut resolution detail */
int doaxes,			/* Do Lab axes */
double tlimit,		/* Total ink limit */
double klimit,		/* K ink limit */
char *outname		/* Output VRML file */
) {
	int i, j;
	FILE *wrl;
	struct {
		double x, y, z;
		double wx, wy, wz;
		double r, g, b;
	} axes[5] = {
		{ 0, 0,  50-GAMUT_LCENT,  2, 2, 100,  .7, .7, .7 },	/* L axis */
		{ 50, 0,  0-GAMUT_LCENT,  100, 2, 2,   1,  0,  0 },	/* +a (red) axis */
		{ 0, -50, 0-GAMUT_LCENT,  2, 100, 2,   0,  0,  1 },	/* -b (blue) axis */
		{ -50, 0, 0-GAMUT_LCENT,  100, 2, 2,   0,  1,  0 },	/* -a (green) axis */
		{ 0,  50, 0-GAMUT_LCENT,  2, 100, 2,   1,  1,  0 },	/* +b (yellow) axis */
	};
	int vix;						/* Vertex index */
	DCOUNT(coa, MXDI, p->inputChan, 0, 0, 2);

	double col[1 << MXDI][3];		/* Color asigned to each major vertex */
	int res;

	if (tlimit < 0.0)
		tlimit = p->inputChan;
	if (klimit < 0.0)
		klimit = 1.0;

	/* Asign some colors to the combination nodes */
	for (i = 0; i < (1 << p->inputChan); i++) {
		int a, b, c, j;
		double h;

		j = (i ^ 0x5a5a5a5a) % (1 << p->inputChan);
		h = (double)j/((1 << p->inputChan)-1);

		/* Make fully saturated with chosen hue */
		if (h < 1.0/3.0) {
			a = 0;
			b = 1;
			c = 2;
		} else if (h < 2.0/3.0) {
			a = 1;
			b = 2;
			c = 0;
			h -= 1.0/3.0;
		} else {
			a = 2;
			b = 0;
			c = 1;
			h -= 2.0/3.0;
		}
		h *= 3.0;

		col[i][a] = (1.0 - h);
		col[i][b] = h;
		col[i][c] = d_rand(0.0, 1.0);
	}

	if (detail > 0.0)
		res = (int)(100.0/detail);	/* Establish an appropriate sampling density */
	else
		res = 4;

	if (res < 2)
		res = 2;

	if ((wrl = fopen(outname,"w")) == NULL)
		error("Error opening wrl output file '%s'",outname);

	/* Spit out a VRML 2 Object surface of gamut */
	fprintf(wrl,"#VRML V2.0 utf8\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"# Created by the Argyll CMS\n");
	fprintf(wrl,"Transform {\n");
  	fprintf(wrl,"children [\n");
    fprintf(wrl,"	NavigationInfo {\n");
	fprintf(wrl,"		type \"EXAMINE\"        # It's an object we examine\n");
	fprintf(wrl,"	} # We'll add our own light\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    DirectionalLight {\n");
	fprintf(wrl,"        direction 0 0 -1      # Light illuminating the scene\n");
	fprintf(wrl,"        direction 0 -1 0      # Light illuminating the scene\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    Viewpoint {\n");
	fprintf(wrl,"        position 0 0 340      # Position we view from\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");
	if (doaxes != 0) {
		fprintf(wrl,"# Lab axes as boxes:\n");
		for (i = 0; i < 5; i++) {
			fprintf(wrl,"Transform { translation %f %f %f\n", axes[i].x, axes[i].y, axes[i].z);
			fprintf(wrl,"\tchildren [\n");
			fprintf(wrl,"\t\tShape{\n");
			fprintf(wrl,"\t\t\tgeometry Box { size %f %f %f }\n",
			                  axes[i].wx, axes[i].wy, axes[i].wz);
			fprintf(wrl,"\t\t\tappearance Appearance { material Material ");
			fprintf(wrl,"{ diffuseColor %f %f %f} }\n", axes[i].r, axes[i].g, axes[i].b);
			fprintf(wrl,"\t\t}\n");
			fprintf(wrl,"\t]\n");
			fprintf(wrl,"}\n");
		}
		fprintf(wrl,"\n");
	}
	fprintf(wrl,"    Transform {\n");
	fprintf(wrl,"      translation 0 0 0\n");
	fprintf(wrl,"      children [\n");
	fprintf(wrl,"		Shape { \n");
	fprintf(wrl,"		    geometry IndexedFaceSet {\n");
	fprintf(wrl,"				solid FALSE\n");		/* Don't back face cull */
	fprintf(wrl,"				convex TRUE\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"		        coord Coordinate { \n");
	fprintf(wrl,"		            point [			# Verticy coordinates\n");


	/* Itterate over all the faces in the device space */
	/* generating the vertx positions. */
	DC_INIT(coa);
	vix = 0;
	while(!DC_DONE(coa)) {
		int e, m1, m2;
		double in[MXDI];
		double inl[MXDI];
		double out[3];
		double sum;

		/* Scan only device surface */
		for (m1 = 0; m1 < p->inputChan; m1++) {
			if (coa[m1] != 0)
				continue;

			for (m2 = m1 + 1; m2 < p->inputChan; m2++) {
				int x, y;

				if (coa[m2] != 0)
					continue;

				for (e = 0; e < p->inputChan; e++)
					in[e] = (double)coa[e];		/* Base value */

				/* Scan over 2D device space face */
				for (x = 0; x < res; x++) {				/* step over surface */
					in[m1] = x/(res - 1.0);
					for (y = 0; y < res; y++) {
						in[m2] = y/(res - 1.0);

						for (sum = 0.0, e = 0; e < p->inputChan; e++) {
							sum += inl[e] = in[e];
						}
						if (sum >= tlimit) {
							for (e = 0; e < p->inputChan; e++)
								inl[e] *= tlimit/sum;
						}
						if (p->inputChan >= 3 && inl[3] >= klimit)
							inl[3] = klimit;
						p->lookup(p, out, inl);
						fprintf(wrl,"%f %f %f,\n",out[1], out[2], out[0]-50.0);
						vix++;
					}
				}
			}
		}
		/* Increment index within block */
		DC_INC(coa);
	}

	fprintf(wrl,"					]\n");
	fprintf(wrl,"		        }\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"		        coordIndex [ 		# Indexes of poligon Verticies \n");

	/* Itterate over all the faces in the device space */
	/* generating the quadrilateral indexes. */
	DC_INIT(coa);
	vix = 0;
	while(!DC_DONE(coa)) {
		int e, m1, m2;
		double in[MXDI];

		/* Scan only device surface */
		for (m1 = 0; m1 < p->inputChan; m1++) {
			if (coa[m1] != 0)
				continue;

			for (m2 = m1 + 1; m2 < p->inputChan; m2++) {
				int x, y;

				if (coa[m2] != 0)
					continue;

				for (e = 0; e < p->inputChan; e++)
					in[e] = (double)coa[e];		/* Base value */

				/* Scan over 2D device space face */
				/* Only output quads under the total ink limit */
				/* Scan over 2D device space face */
				for (x = 0; x < res; x++) {				/* step over surface */
					for (y = 0; y < res; y++) {
						if (x < (res-1) && y < (res-1)) {
							fprintf(wrl,"%d, %d, %d, %d, -1\n", 
							vix, vix + 1, vix + 1 + res, vix + res);
						}
						vix++;
					}
				}
			}
		}
		/* Increment index within block */
		DC_INC(coa);
	}

	fprintf(wrl,"				]\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"				colorPerVertex TRUE\n");
	fprintf(wrl,"		        color Color {\n");
	fprintf(wrl,"		            color [			# RGB colors of each vertex\n");

	/* Itterate over all the faces in the device space */
	/* generating the vertx colors. */
	DC_INIT(coa);
	vix = 0;
	while(!DC_DONE(coa)) {
		int e, m1, m2;
		double in[MXDI];

		/* Scan only device surface */
		for (m1 = 0; m1 < p->inputChan; m1++) {
			if (coa[m1] != 0)
				continue;

			for (m2 = m1 + 1; m2 < p->inputChan; m2++) {
				int x, y;

				if (coa[m2] != 0)
					continue;

				for (e = 0; e < p->inputChan; e++)
					in[e] = (double)coa[e];		/* Base value */

				/* Scan over 2D device space face */
				for (x = 0; x < res; x++) {				/* step over surface */
					double xb = x/(res - 1.0);
					for (y = 0; y < res; y++) {
						int v0, v1, v2, v3;
						double yb = y/(res - 1.0);
						double rgb[3];

						for (v0 = 0, e = 0; e < p->inputChan; e++)
							v0 |= coa[e] ? (1 << e) : 0;		/* Binary index */

						v1 = v0 | (1 << m2);				/* Y offset */
						v2 = v0 | (1 << m2) | (1 << m1);	/* X+Y offset */
						v3 = v0 | (1 << m1);				/* Y offset */

						/* Linear interp between the main verticies */
						for (j = 0; j < 3; j++) {
							rgb[j] = (1.0 - yb) * (1.0 - xb) * col[v0][j]
							       +        yb  * (1.0 - xb) * col[v1][j]
							       + (1.0 - yb) *        xb  * col[v3][j]
							       +        yb  *        xb  * col[v2][j];
						}
						fprintf(wrl,"%f %f %f,\n",rgb[1], rgb[2], rgb[0]);
						vix++;
					}
				}
			}
		}
		/* Increment index within block */
		DC_INC(coa);
	}

	fprintf(wrl,"					] \n");
	fprintf(wrl,"		        }\n");
	fprintf(wrl,"		    }\n");
	fprintf(wrl,"		    appearance Appearance { \n");
	fprintf(wrl,"		        material Material {\n");
	fprintf(wrl,"					transparency 0.0\n");
	fprintf(wrl,"					ambientIntensity 0.3\n");
	fprintf(wrl,"					shininess 0.5\n");
	fprintf(wrl,"				}\n");
	fprintf(wrl,"		    }\n");
	fprintf(wrl,"		}	# end Shape\n");
	fprintf(wrl,"      ]\n");
	fprintf(wrl,"    }\n");

	fprintf(wrl,"\n");
	fprintf(wrl,"  ] # end of children for world\n");
	fprintf(wrl,"}\n");

	if (fclose(wrl) != 0)
		error("Error closing output file '%s'",outname);
}

