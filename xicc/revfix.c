
/* 
 * ICC reprocess to give true BtoA1 by inverting
 * the AtoB1 table, and also correct the neutral
 * axis of the BtoA0 and BtoA2 tables.
 *
 *
 * Author:  Graeme W. Gill
 * Date:    9/7/00
 * Version: 1.00
 *
 * Copyright 2000 Graeme W. Gill
 * Please refer to Licence.txt file for details.
 */

/* TTBD:
 *
 * Add support for proper gamut mapping, just like profile,
 * or deprecate revfix by modifying profile to work from
 * an existing profile ?
 * 
 * Remove the auxiliary fixup stuff when we have implemented
 * optimised separations.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "copyright.h"
#include "config.h"
#include "xicc.h"

#define USE_CAM_CLIP_OPT		/* Clip in CAM Jab space rather than Lab */
#undef DEBUG		/* Print each value changed */

void usage(void) {
	fprintf(stderr,"Invert AtoB1 to make BtoA1 for CMYK profiles, V1.20\n");
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	fprintf(stderr,"usage: revfix [-options] iccin iccout\n");
	fprintf(stderr," -v         Verbose\n");
	fprintf(stderr," -0         Process perceptual\n");
	fprintf(stderr," -1         Process absolute/relative colorimetric\n");
	fprintf(stderr," -2         Process saturation\n");
	fprintf(stderr," -r res     Override BtoA1 Clut res\n");
	fprintf(stderr," -k [ezhxr] e = same K as existing BtoA table (def)\n");
	fprintf(stderr,"            z = zero, h = 0.5 K, x = max K, r = ramp K\n");
	fprintf(stderr," -k p stle stpo enle enpo shape\n");
	fprintf(stderr,"            p = curve parameters\n");
	fprintf(stderr,"            stle: K level at White 0.0 - 1.0\n");
	fprintf(stderr,"            stpo: start point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"            enpo: End point of transition Wh 0.0 - Bk 1.0\n");
	fprintf(stderr,"            enle: K level at Black 0.0 - 1.0\n");
	fprintf(stderr,"            shape: 1.0 = straight, 0.0-1.0 concave, 1.0-2.0 convex\n");
	fprintf(stderr," -l tlimit  set total ink limit, 0 - 400%% (default none)\n");
	fprintf(stderr," -L klimit  set black ink limit, 0 - 100%% (default none)\n");
	fprintf(stderr," -p absprof Include abstract profile in output tables\n");
//	fprintf(stderr," -s         Use internal optimized separation for CMYK\n");
	exit(1);
}

/* ------------------------------------------- */
/* structure to support icc Lut initialisation calbacks */

/* ~~~Note that  we're not coping with a matrix or XYZ PCS properly here. */

struct _callback {
	int verb;			/* Verbosity */
	int total, count, last;	/* Progress count information */
	icColorSpaceSignature pcsspace;
	int inking;			/* k inking algorithm */
	icxLuLut *BtoA;		/* BtoA of table being processed */
	icxLuLut *AtoB;		/* AtoB of table being processed */
	icxLuLut *AtoB1;	/* AtoB of colorimetric table */

	icRenderingIntent abs_intent;	/* Desired abstract profile rendering intent */
	icxLuBase *abs_luo;				/* abstract profile tranform in PCS, NULL if none */

}; typedef struct _callback callback;


/* Utility to handle abstract profile application to PCS */
/* PCS in creating output table is always XYZ or Lab relative colorimetric, */
/* and abstract profile is absolute or relative, and will be */
/* XYZ if absolute, and PCS if relative. */
static void do_abstract(callback *p, double out[3], double in[3]) {
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];

//printf("~1 do_abstract got %f %f %f\n",in[0],in[1],in[2]);

	if (p->abs_intent == icAbsoluteColorimetric) {
		if (p->pcsspace == icSigLabData) {
			icmLab2XYZ(&icmD50, out, out);
//printf("~1 after Lab 2 XYZ got %f %f %f\n",out[0],out[1],out[2]);
		}
		p->AtoB1->plu->XYZ_Rel2Abs(p->AtoB1->plu, out, out);
//printf("~1 after rel to abs got %f %f %f\n",out[0],out[1],out[2]);
	}

	p->abs_luo->lookup(p->abs_luo, out, out);
//printf("~1 after abs_luo got %f %f %f\n",out[0],out[1],out[2]);

	if (p->abs_intent == icAbsoluteColorimetric) {
		p->AtoB1->plu->XYZ_Abs2Rel(p->AtoB1->plu, out, out);
//printf("~1 after abs2rel got %f %f %f\n",out[0],out[1],out[2]);
		if (p->pcsspace == icSigLabData) {
			icmXYZ2Lab(&icmD50, out, out);
//printf("~1 after XYZ to Lab got %f %f %f\n",out[0],out[1],out[2]);
		}
	}
//printf("~1 returning %f %f %f\n\n",out[0],out[1],out[2]);
}

/* - - - - - - - - - */
/* New input table */
void Lab_Labp(void *cntx, double out[3], double in[3]) {
	callback *p = (callback *)cntx;

#ifdef DEBUG
	printf("Got Lab %f %f %f\n",in[0],in[1],in[2]);
#endif
	if (p->AtoB != p->AtoB1) {
		/* Non-colorimetric, use existing input table */
		if (p->BtoA->input(p->BtoA, out, in) > 1)
			error ("%d, %s",p->BtoA->pp->errc,p->BtoA->pp->err);
	} else {
		/* Colorimetric, use inverse AtoB output */
		if (p->AtoB1->inv_output(p->AtoB1, out, in) > 1)
			error ("%d, %s",p->AtoB1->pp->errc,p->AtoB1->pp->err);
	}
#ifdef DEBUG
	printf("New Lab' %f %f %f\n",out[0],out[1],out[2]);
#endif
}

/* - - - - */
/*  clut  */

/* Normal CLUT routine */
void Labp_CMYKp(void *cntx, double out[4], double in[3]) {
	double temp[4], targetk;
	int rv;
	callback *p = (callback *)cntx;

#ifdef DEBUG
	printf("Got Lab' %f %f %f\n",in[0],in[1],in[2]);
#endif

	if (p->inking == 0) {	/* If we are copying existing K */
		
		/* Figure out what K value was previously here */
		if (p->AtoB != p->AtoB1) {
			/* Simple because BtoA input & output tables don't change */
			/* Figure out what DEV' K value the BtoA table currently has for this PCS' */
			if (p->BtoA->clut(p->BtoA, temp, in) > 1)
				error ("%d, %s",p->BtoA->pp->errc,p->BtoA->pp->err);
		} else {
			/* More complicated because old BtoA in/out tables are different */
			/* from the new ones. */
			/* We know that new BtoA in/out tables are inverse of AtoB in/out, */
			/* so we don't have to use BtoA1->inv_input, & BtoA1->inv_output */
			/* Convert PCS' to PCS */
			if (p->AtoB->output(p->AtoB, temp, in) > 1)
				error ("%d, %s",p->AtoB->pp->errc,p->AtoB->pp->err);
	
			/* Figure out what DEV K value the BtoA table currently has for this PCS */
			if (((icxLuBase *)p->BtoA)->lookup((icxLuBase *)p->BtoA, temp, temp) > 1)
				error ("%d, %s",p->BtoA->pp->errc,p->BtoA->pp->err);
	
			/* Convert DEV to DEV' */
			if (p->AtoB->input(p->AtoB, temp, temp) > 1)
				error ("%d, %s",p->AtoB->pp->errc,p->AtoB->pp->err);
		}

		targetk = temp[3];
#ifdef DEBUG
		printf("Got existing CMYK' %f %f %f %f\n",temp[0],temp[1],temp[2],temp[3]);
#endif
	}

	/* Copy the Lab in */
	temp[0] = in[0];
	temp[1] = in[1];
	temp[2] = in[2];

	if (p->AtoB != p->AtoB1) {
		double tt[4];

		/* Can't assume B2A in/out tables are inverses of AtoB */
		/* Convert PCS' -> PCS for this table */
		if (p->BtoA->inv_input(p->BtoA, temp, temp) > 1)
			error ("%d, %s",p->BtoA->pp->errc,p->BtoA->pp->err);
		/* Convert PCS -> PCS' for colorimetric */
		if (p->AtoB1->inv_output(p->AtoB1, temp, temp) > 1)
			error ("%d, %s",p->AtoB->pp->errc,p->AtoB->pp->err);

		/* Need to correct target k value */
		tt[0] = tt[1] = tt[2] = 0.5;
		tt[3] = targetk;
		/* DEV->DEV' for this table */
		p->BtoA->output(p->BtoA, tt, tt);
		/* DEV -> DEV' for colorimetric */
		p->AtoB1->input(p->AtoB1, tt, tt);
		targetk = tt[3];
	}

	/* Abstract profile applied before inversion */
	if (p->abs_luo != NULL) {
		do_abstract(p, temp, temp);
	}

	/* Invert AtoB1 clut, using set inking policy */
	out[3] = targetk;
	/* PCS' -> DEV' colorimetric */
	if ((rv = p->AtoB1->inv_clut(p->AtoB1, out, temp)) > 1)
		error ("%d, %s",p->AtoB1->pp->errc,p->AtoB1->pp->err);

	/* ~~~ Note that the ink limit will be wrong for non-colorimetric, */
	/* since AtoB1->inv_clut will be assuming A->toB1->inv_input as the output table, */
	/* while we will actually be using BtoA->output ~~~~ */
	/* Need to override ink limit computation function in icx for non-colorimetric. */

//#ifdef DEBUG
//		if (rv != 0) {
//			printf("Inversion clipped\n");
//		}
//#endif
	if (p->AtoB != p->AtoB1) {
		/* Can't assume B2A in/out tables are inverses of AtoB */
		/* Converts DEV' -> DEV colorimetric */
		if (p->AtoB1->inv_input(p->AtoB1, out, out) > 1)
			error ("%d, %s",p->AtoB->pp->errc,p->AtoB->pp->err);
		/* Convert DEV -> DEV' for this table */
		if (p->BtoA->inv_output(p->BtoA, out, out) > 1)
			error ("%d, %s",p->BtoA->pp->errc,p->BtoA->pp->err);
	}
#ifdef DEBUG
	printf("New CMYK' %f %f %f %f\n",out[0],out[1],out[2],out[3]);
	printf("\n");
#endif

	if (p->verb) {		/* Output percent intervals */
		int pc;
		p->count++;
		pc = p->count * 100.0/p->total + 0.5;
		if (pc != p->last) {
			printf("\r%2d%%",pc), fflush(stdout);
			p->last = pc;
		}
	}
}

/* - - - - - - - - - */
/* New output table */
void CMYKp_CMYK(void *cntx, double out[4], double in[4]) {
	callback *p = (callback *)cntx;

#ifdef DEBUG
	printf("Got CMYK' %f %f %f %f\n",in[0],in[1],in[2],in[3]);
#endif
	if (p->AtoB != p->AtoB1) {
		/* Non-colorimetric, use existing output table */
		if (p->BtoA->output(p->BtoA, out, in) > 1)
			error ("%d, %s",p->BtoA->pp->errc,p->BtoA->pp->err);
	} else {
		/* Colorimetric, use inverse AtoB input */
		if (p->AtoB1->inv_input(p->AtoB1, out, in) > 1)
			error ("%d, %s",p->AtoB1->pp->errc,p->AtoB1->pp->err);
	}
#ifdef DEBUG
	printf("New CMYK %f %f %f %f\n",out[0],out[1],out[2],out[3]);
#endif
}


/* ------------------------------------------- */

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	char in_name[MAXNAMEL+1];
	char out_name[MAXNAMEL+1];
	char abs_name[MAXNAMEL+1] = "\000";	/* Abstract profile name */
	icmFile *rd_fp, *wr_fp;
	icc *icco;
	int verb = 0;
	int clutres = 0;
	int do0 = 0;
	int do1 = 0;
	int do2 = 0;
	int inking = 0;			/* Default copy from existing */
	double Kstle, Kstpo, Kenle, Kenpo, Kshap;
	int tlimit = -1;		/* Total ink limit as a % */
	int klimit = -1;		/* Black ink limit as a % */
	int intsep = 0;			/* Not implimented in xicc yet ??? */
	int rv = 0;
	char buf[200];
	error_program = argv[0];

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
			else if (argv[fa][1] == '0') {
				do0 = 1;
			}
			else if (argv[fa][1] == '1') {
				do1 = 1;
			}
			else if (argv[fa][1] == '2') {
				do2 = 1;
			}
			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R') {
				fa = nfa;
				if (na == NULL) usage();
				clutres = atoi(na);
			}

			/* Inking rule */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				fa = nfa;
				if (na == NULL) usage();
    			switch (na[0]) {
					case 'e':
					case 'E':
						inking = 0;		/* Use existing table as guide */
						break;
					case 'z':
					case 'Z':
						inking = 1;		/* Use minimum k */
						break;
					case 'h':
					case 'H':
						inking = 2;		/* Use 0.5 k */
						break;
					case 'x':
					case 'X':
						inking = 3;		/* Use maximum k */
						break;
					case 'r':
					case 'R':
						inking = 4;		/* Use ramp */
						break;
					case 'p':
					case 'P':
						inking = 5;		/* Use curve parameter */
						++fa;
						if (fa >= argc) usage();
						Kstle = atof(argv[fa]);

						++fa;
						if (fa >= argc) usage();
						Kstpo = atof(argv[fa]);

						++fa;
						if (fa >= argc || argv[fa][0] == '-') usage();
						Kenpo = atof(argv[fa]);

						++fa;
						if (fa >= argc) usage();
						Kenle = atof(argv[fa]);

						++fa;
						if (fa >= argc || argv[fa][0] == '-') usage();
						Kshap = atof(argv[fa]);
						break;
					default:
						usage();
				}
			}
			else if (argv[fa][1] == 'l') {
				fa = nfa;
				if (na == NULL) usage();
				tlimit = atoi(na);
			}
			else if (argv[fa][1] == 'L') {
				fa = nfa;
				if (na == NULL) usage();
				klimit = atoi(na);
			}

			/* Use internal separation */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				intsep = 1;
			}

			/* Abstract profile */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				if (na == NULL) usage();
				fa = nfa;
				strncpy(abs_name,na,MAXNAMEL); abs_name[MAXNAMEL] = '\000';
			}

			else 
				usage();
		} else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(in_name,argv[fa++],MAXNAMEL); in_name[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(out_name,argv[fa++],MAXNAMEL); out_name[MAXNAMEL] = '\000';

	/* Open up the profile for reading */
	if ((rd_fp = new_icmFileStd_name(in_name,"r")) == NULL)
		error ("Can't open file '%s'",in_name);

	if ((icco = new_icc()) == NULL)
		error ("Creation of ICC object failed");

	/* Read header etc. */
	if ((rv = icco->read(icco,rd_fp,0)) != 0)
		error ("%d, %s",rv,icco->err);

	/* Read every tag */
	if (icco->read_all_tags(icco) != 0) {
		error("Unable to read all tags: %d, %s",icco->errc,icco->err);
	}

	rd_fp->del(rd_fp);

	/* ======================= */
	/* Check that it is a suitable icc */
	{
		icmHeader *rh = icco->header;

		if (rh->deviceClass != icSigOutputClass)
			error("Profile isn't an output device profile");

    	if (rh->colorSpace != icSigCmykData)
			error("Profile isn't for a CMYK device");

    	if (rh->pcs != icSigLabData)
			error("Profile is not using a PCS of Lab - can't cope with this yet");
	}

	if (verb && inking == 5) {
		double tL;
		printf("K parameters are are %f %f %f %f %f\n",Kstle, Kstpo, Kenpo, Kenle, Kshap);
		/* Plot expected K locus goal */
		for (tL = 100.0; tL >= 0.0; tL -= 10.0) {
			double rv, L;

			L = 0.01 * tL;

			/* Code from xlut.c: */
			/* Invert sense of L, so that 0.0 = white, 1.0 = black */
			L = 1.0 - L;

			if (L <= Kstpo && L <= Kenpo) {
				/* We are at white level */
				rv = Kstle;
			} else if (L >= Kstpo && L >= Kenpo) {
				/* We are at black level */
				rv = Kenle;
			} else {
				double g;
				/* We must be on the curve from start to end levels */

				if (Kstpo > Kenpo) {
					rv = (L - Kenpo)/(Kstpo - Kenpo);
				} else {
					rv = (L - Kstpo)/(Kenpo - Kstpo);
				}

				g = Kshap/2.0;

				/* A value of 0.5 will be tranlated to g */
				rv = rv/((1.0/g - 2.0) * (1.0 - rv) + 1.0);

				/* Transition between start end end levels */
				rv = rv * (Kenle - Kstle) + Kstle;
			}

			/* To be safe */
			if (rv < 0.0)
				rv = 0.0;
			else if (rv > 1.0)
				rv = 1.0;

			printf("L = %f, K locus = %f\n",tL, rv);
		}
	}

	/* ======================= */
	{
		int ii;						/* Current intent */
		icmLut *done[3];			/* record pointers to done Luts here */
		icRenderingIntent intent;
		icTagSignature sig;
		xicc *xicco;
		callback cb;						/* Callback support stucture */
		icxInk ink;							/* Ink parameters */
		icmLuAlgType alg;					/* Type of lookup algorithm */
		icmFile *abs_fp;	/* Abstract profile transform: */
		icc *abs_icc;
		xicc *abs_xicc;
		
		/* Wrap with an expanded icc */
		if ((xicco = new_xicc(icco)) == NULL)
			error ("Creation of xicc failed");

		/* Setup CMYK -> Lab conversion (Fwd) object */
		if (tlimit >= 0)
		ink.tlimit = tlimit/100.0;	/* Set a total ink limit */
		else
			ink.tlimit = -1.0;			/* Don't use a limit */

		if (klimit >= 0)
			ink.klimit = klimit/100.0;	/* Set a black ink limit */
		else
			ink.klimit = -1.0;			/* Don't use a limit */

		ink.c.Ksmth = ICXINKDEFSMTH;	/* Default curve smoothing */

		if (inking == 0) {
			ink.k_rule = icxKvalue;		/* K is auxiliary target */

		} else if (inking == 1) {		/* Use minimum */
			ink.c.Kstle = 0.0;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 0.0;
			ink.c.Kshap = 1.0;
		} else if (inking == 2) {		/* Use 0.5  */
			ink.k_rule = icxKluma5;
			ink.c.Kstle = 0.5;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 0.5;
			ink.c.Kshap = 1.0;
		} else if (inking == 3) {		/* Use maximum  */
			ink.k_rule = icxKluma5;
			ink.c.Kstle = 1.0;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 1.0;
			ink.c.Kshap = 1.0;
		} else if (inking == 4) {		/* Use ramp  */
			ink.k_rule = icxKluma5;
			ink.c.Kstle = 0.0;
			ink.c.Kstpo = 0.0;
			ink.c.Kenpo = 1.0;
			ink.c.Kenle = 1.0;
			ink.c.Kshap = 1.0;
		} else {						/* Use specified curve */
			ink.k_rule = icxKluma5;
			ink.c.Kstle = Kstle;
			ink.c.Kstpo = Kstpo;
			ink.c.Kenpo = Kenpo;
			ink.c.Kenle = Kenle;
			ink.c.Kshap = Kshap;
		}

		cb.verb = verb;
		cb.count = 0;
		cb.last = -1;
		cb.inking = inking;

		/* Setup our access to the device characteristic */
		if ((cb.AtoB1 = (icxLuLut *)xicco->get_luobj(xicco,
		                            ICX_CLIP_NEAREST
#ifdef USE_CAM_CLIP_OPT
									| ICX_CAM_CLIP
#endif
		                            | (intsep ? ICX_INT_SEPARATE : 0),
		                            icmFwd, icRelativeColorimetric,
		                            icmSigDefaultData, icmLuOrdNorm, NULL, &ink)) == NULL)
			error ("%d, %s",xicco->errc, xicco->err);

		cb.AtoB1->spaces((icxLuBase *)cb.AtoB1, NULL, NULL, NULL, NULL, &alg,
		                                              NULL, NULL, &cb.pcsspace);
		if (alg != icmLutType)
			error("Forward conversion is not a Lut");

		/* Open up the abstract profile if supplied, and setup luo */
		if (abs_name[0] != '\000') {
			if ((abs_fp = new_icmFileStd_name(abs_name,"r")) == NULL)
				error ("Can't open abstract profile file '%s'",abs_name);
			
			if ((abs_icc = new_icc()) == NULL)
				error ("Creation of Abstract profile ICC object failed");
	
			/* Read header etc. */
			if ((rv = abs_icc->read(abs_icc,abs_fp,0)) != 0)
				error ("%d, %s",rv,abs_icc->err);
			abs_icc->header;
	
			if (abs_icc->header->deviceClass != icSigAbstractClass)
				error("Abstract profile isn't an abstract profile");
	
			/* Take intended abstract intent from profile itself */
			if ((cb.abs_intent = abs_icc->header->renderingIntent) != icAbsoluteColorimetric)
				cb.abs_intent = icRelativeColorimetric;
	
			/* Wrap with an expanded icc */
			if ((abs_xicc = new_xicc(abs_icc)) == NULL)
				error ("Creation of abstract profile xicc failed");

			/* The abstract profile intent is assumed to determine how it gets applied. */
			/* Make abstract PCS XYZ if icAbsoluteColorimetric is needed. */
			if ((cb.abs_luo = abs_xicc->get_luobj(abs_xicc, 0, icmFwd, cb.abs_intent,
		        (cb.pcsspace == icSigLabData && cb.abs_intent == icRelativeColorimetric)
				             ? icSigLabData : icSigXYZData,
				icmLuOrdNorm, NULL, NULL)) == NULL)
				error ("%d, %s",abs_icc->errc, abs_icc->err);
		} else {
			cb.abs_luo = NULL;
		}
			
		/* for all intents, and not already processed */
		for (ii = 0; ii <= 2; ii++) {
			int i;
			icmLut *wo;

			switch(ii) {
				case 0:
					intent = icRelativeColorimetric;
					sig    = icSigBToA1Tag;
					if (do1 == 0)
						continue;
					break;
				case 1:
					intent = icPerceptual;
					sig    = icSigBToA0Tag;
					if (do0 == 0)
						continue;
					break;
				case 2:
					intent = icSaturation;
					sig    = icSigBToA2Tag;
					if (do2 == 0)
						continue;
					break;
			}

			if (intent == icRelativeColorimetric) {
				cb.AtoB = cb.AtoB1;
			} else {
				/* Setup CMYK -> Lab conversion (Fwd) object */
				if ((cb.AtoB = (icxLuLut *)xicco->get_luobj(xicco, 0, icmFwd, intent,
				                         icmSigDefaultData, icmLuOrdNorm, NULL, NULL)) == NULL)
					error ("%d, %s",xicco->errc, xicco->err);
	
				cb.AtoB->spaces((icxLuBase *)cb.AtoB, NULL, NULL, NULL, NULL, &alg, NULL, NULL, NULL);
				if (alg != icmLutType)
					error("Forwards conversion is not a Lut");
			}

			/* Setup Lab -> CMYK conversion (Bwd) object */
			if ((cb.BtoA = (icxLuLut *)xicco->get_luobj(xicco, 0, icmBwd, intent,
			                           icmSigDefaultData, icmLuOrdNorm, NULL, NULL)) == NULL)
				error ("%d, %s",xicco->errc, xicco->err);

			cb.BtoA->spaces((icxLuBase *)cb.BtoA, NULL, NULL, NULL, NULL, &alg, NULL, NULL, NULL);
			if (alg != icmLutType)
				error("Backwards conversion is not a Lut");

			/* Try and read the tag from the file */
			wo = (icmLut *)icco->read_tag(icco, sig);
			if (wo == NULL) 
				error("Can't find %s", icm2str(icmRenderingIntent, intent));

			/* Need to check that the cast is appropriate. */
			if (wo->ttype != icSigLut16Type && wo->ttype != icSigLut8Type)
				error("Lut table isn't Lut8 or Lut16 Type");

			/* Set reverse input table resolution to same as fwd output */
			wo->inputEnt = cb.AtoB->lut->outputEnt;

			/* Let user override for BtoA1 */
			if (sig == icSigBToA1Tag && clutres > 0) {
				if (verb)
					printf("Overriding existing clut resolution %d with %d\n",wo->clutPoints,clutres);
		    	wo->clutPoints = clutres;
			}
			/* If Lut8, make sure the input and output tables have 256 enries. */
			if (wo->ttype == icSigLut8Type) {
		    	wo->inputEnt = 256;
		    	wo->outputEnt = 256;
			}
			wo->allocate((icmBase *)wo);/* Allocate space */

			/* Make sure that we don't process a table twice */
			done[ii] = wo;
			for (i = ii-1; i >= 0; i--) {
				if (done[i] == wo) {
					break;
				}
			}
			if (i >= 0) {
				if (verb)
					printf("Skipping %s table - already done\n", icm2str(icmRenderingIntent, intent));

			} else {
		
				if (verb)
					printf("About to start processing %s\n", icm2str(icmRenderingIntent, intent));
	
				if (cb.verb) {
					for (cb.total = 1, i = 0; i < wo->inputChan; i++, cb.total *= wo->clutPoints)
						; 
					printf(" 0%%"), fflush(stdout);
				}
				/* Use helper function to do the hard work. */
				if (wo->set_tables(wo,
						&cb,						/* Context */
						icSigLabData, 				/* Input color space */
						icSigCmykData, 				/* Output color space */
						Lab_Labp,					/* Linear input transform Lab->Lab' */
						NULL, NULL,					/* Use default Lab' range */
						Labp_CMYKp,					/* Lab' -> CMYK' transfer function */
						NULL, NULL,					/* Use default CMYK' range */
						CMYKp_CMYK) != 0)			/* Output transfer function, CMYK'->CMYK (NULL = deflt) */
					error("Setting 16 bit Lab->CMYK Lut failed: %d, %s",icco->errc,icco->err);
	
				if (verb)
					printf("\nDone processing %s\n", icm2str(icmRenderingIntent, intent));
	
			}

			/* Done with this intents lookup object */
			cb.BtoA->del((icxLuBase *)cb.BtoA);
			if (cb.AtoB != cb.AtoB1)
				cb.AtoB->del((icxLuBase *)cb.AtoB);
		}
		/* Done with colorimetric intents AtoB1 lookup objects, and xicc */
		cb.AtoB1->del((icxLuBase *)cb.AtoB1);
		xicco->del(xicco);

		if (cb.abs_luo != NULL) {		/* Free up abstract transform */
			cb.abs_luo->del(cb.abs_luo);
			abs_xicc->del(abs_xicc);
			abs_icc->del(abs_icc);
			abs_fp->del(abs_fp);
		}

	}
	/* ======================================= */
	
	/* Open up the other profile for writing */
	if ((wr_fp = new_icmFileStd_name(out_name,"w")) == NULL)
		error ("Can't open file '%s'",out_name);

	if ((rv = icco->write(icco,wr_fp,0)) != 0)
		error ("Write file: %d, %s",rv,icco->err);

	wr_fp->del(wr_fp);
	icco->del(icco);

	return 0;
}

