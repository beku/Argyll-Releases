
/* 
 * icclink
 *
 * Link two device profiles to create a Device Link Profile
 *
 * Author:  Graeme W. Gill
 * Date:    25/11/00
 * Version: 2.10
 *
 * Copyright 2000 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* TTBD:
 *
 * Add input & output ink limit to allow for input CMYK profiles actual gamut.
 *     - should ink limit be stashed in profiles (yes). 
 *     - should ink limit be guessed if not in profiles (yes). 
 *
 * Abstract link support intent doesn't work properly for anyithing
 * other than absolute. This should really be fixed.
 */

/* NOTES:

	Normally the device side per channel curves are copied from
	the source profiles to the link profile on the assumption that
	the raw linearisation they do is good, and should be maintained
	for best overall transform accuracy. Since the intermediate
	link is done in an Lab like space, then this linearisation
	will even out quantisation error introduced by the Lut
	in perceptual space. In the case of a profile with a native
	XYZ PCS, these assumptions break, since the device curves
	would generally be linearising the device to Y, not a perceptual
	space. For this reason we do not transfer the device linearisation
	curves for XYZ based profiles.
	It would be possible to use device linearisation curves for XYZ based
    profiles, if the curve regeneration code worked as desired. This
    needs more investigation, development and testing.


	Colorspace representations are a bit of a mess. It's hard to know what space
	color is in at any point, and difficult to transform to match
	some other element. Putting the different colorspace support within
	the profile transforms is neat, but may not be flexible enough, since the
	needed information to do a transform (white point, viewing conditions)
	might be a bit too closely bound into the profile. An alternative might be
	to expand the colorspace definition from a tag to include all the
	other needed information, and create a general colorspace adapter
	to transform from one to the other, or to limit connections
	to a cannonical spaces such as XYZ.
	One big cause of confusion for user and implimentor is how to
	handle the conflicting intents. What does it mean to link
	a relative source to absolute abstract to CAM destination ???

	A non-absolute abstract profile is poorly defined (just like
	the device profiles) if it doesn't define the viewing conditions
	within which the transform is defined. 

	Having separated creating a gamut surface for a raster out
	(for better modularity), it then creates problems in applying
	an abstact transform. The abstract transform can't really be
	applied to a gamut surface, it really needs to be applied
	to the image data before gamut surfac extraction. This is yet
	another avenue for user & implementor confusion, with regard
	to intents, even without allowing the user to set the intent
	of the abstract profile.

 */

#undef XYZCURVE		/* Enable device linearisation curve regeneration (xdevlin) */
					/* stuff for XYZ profile */
					/* This seems to make accuracy worse with the current xdevlin.c */
					/* Need to re-implement xdevlin.c using xlut2.c optimisation code, */
					/* and test it again. */

#undef USE_MERGE_CLUT_OPT	/* When using inverse A2B table, merge the output luts */
							/* with the clut for faster operation, and clipping in */
							/* Jab space. Turned off because it affects the accuracy too much. */

#define USE_CAM_CLIP_OPT	/* Clip out of gamut in CAM space rather than XYZ or L*a*b* */

#define KHACKWIDTH 0		/* Number of grid points outside diagonal to force to K only */ 

#undef DEBUG		/* test a single value out */
#undef NEUTKDEBUG	/* print info about neutral L -> K mapping */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "copyright.h"
#include "config.h"
#include "icc.h"
#include "xicc.h"
#include "xdevlin.h"
#include "gamut.h"
#include "gammap.h"

void usage(char *diag, ...) {
	int i;
	fprintf(stderr,"Link ICC profiles, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licenced under the GPL\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: icclink [options] inprofile outprofile linkedprofile\n");
	fprintf(stderr," -v            Verbose\n");
	fprintf(stderr," -V            Verify existing profile, rather than link\n");
	fprintf(stderr," -q [lmhu]     Quality - Low, Medium (def), High, Ultra\n");
	fprintf(stderr," -r res        Override clut res. set by -q\n");
	fprintf(stderr," -n            Don't preserve device linearization curves in result\n");
	fprintf(stderr," -f            Special :- Force neutral colors to be K only output\n");
	fprintf(stderr," -F            Special :- Force all colors to be K only output\n");
	fprintf(stderr," -p absprof    Include abstract profile in link\n");
	fprintf(stderr," -s            Simple Mode (default)\n");
	fprintf(stderr," -g [src.gam]  Gamut Mapping Mode [optional source image gamut]\n");
	fprintf(stderr," -G [src.gam]  Gamut Mapping Mode using inverse outprofile A2B\n");
	fprintf(stderr," Simple Mode Options:\n");
	fprintf(stderr," -i in_intent  p = perceptual, r = relative colorimetric,\n");
	fprintf(stderr,"               s = saturation, a = absolute colorimetric\n");
	fprintf(stderr," -o out_intent p = perceptual, r = relative colorimetric,\n");
	fprintf(stderr,"               s = saturation, a = absolute colorimetric\n");
	fprintf(stderr," Mapping Mode Options:\n");
	fprintf(stderr," -i intent     set linking intent from the following choice:\n");
	for (i = 0; ; i++) {
		icxGMappingIntent gmi;
		if (xicc_enum_gmapintent(&gmi, i))
			break;

		fprintf(stderr,"               %d: %s\n",i,gmi.desc);
	}
	fprintf(stderr,"               (Can also use -ia, -ir, -ip or -is as aliases.)\n");
	fprintf(stderr," -w [J,a,b]    Use forced whitepoint hack\n");
//	fprintf(stderr," -WJ,a,b       Forced whitepoint adjustment by delta Jab\n");
	fprintf(stderr," -c viewcond   set input viewing conditions for %s,\n",icxcam_description(cam_default));
	fprintf(stderr,"               either an enumerated choice, or a parameter\n");
	fprintf(stderr," -d viewcond   set output viewing conditions for %s,\n",icxcam_description(cam_default));
	fprintf(stderr,"               either an enumerated choice, or a parameter\n");
	fprintf(stderr,"               Enumerated Viewing Conditions:\n");
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
	fprintf(stderr," Inverse outprofile A2B Options:\n");
	fprintf(stderr," -k [tzhxr]    CMYK Black generation\n");
	fprintf(stderr,"               t = transfer K from input to output\n");
	fprintf(stderr,"               z = zero K, h = 0.5 K (def), x = maximum K, r = ramp K\n");
	fprintf(stderr," -k p stle stpo enpo enle shape\n");
	fprintf(stderr,"               p = black locus generation curve parameters\n");
	fprintf(stderr," -k q stle0 stpo0 enpo0 enle0 shape0 stle2 stpo2 enpo2 enle2 shape2\n");
	fprintf(stderr,"               q = transfer input K to dual curve limits\n");
	fprintf(stderr," -l tlimit     set total ink limit, 0 - 400%% (default none)\n");
	fprintf(stderr," -L klimit     set black ink limit, 0 - 100%% (default none)\n");
	exit(1);
}

/* ------------------------------------------- */
/* structures to support icc calbacks */

/* Information needed from a profile */
struct _profinfo {
	/* Setup parameters */
	icRenderingIntent intent;	/* Selected ICC rendering intent */
	icxViewCond vc;				/* Viewing Condition for CAM */
	int inking;					/* k inking algorithm, 0 = input, 1 = min, */
								/* 2 = 0.5, 3 = max, 4 = curve, 5 = ramp */
	icxInk ink;					/* Ink parameters */

	/* Operational parameters */
	icmFile *fp;
	icc *c;
	icmHeader *h;
	xicc *x;
	icxLuBase *luo;				/* Base XLookup type object */
	icmLuAlgType alg;			/* Type of lookup algorithm */
	int chan;					/* Channels */
	int nocurve;				/* Flag to disable explicit device curve */
								/* This is use for native XYZ profiles */
	xdevlin *devlin;			/* Computed Device linearisation curve for XYZ profiles */
								/* Device -> linearised (Not properly implemented!) */
	double wp[3];				/* Lab white point for profile used by wphack */
}; typedef struct _profinfo profinfo;

/* Structure that holds all the color lookup information */
struct _link {
	/* Overall options */
	int verb;
	int total, count, last;	/* Progress count information */
	int mode;		/* 0 = simple mode, 1 = mapping mode, 2 = mapping mode with inverse A2B */
	int quality;	/* 0 = low, 1 = medium, 2 = high, 3 = ultra */
	int clutres;	/* 0 = quality default, !0 = override, then actual during link */

	int nhack;		/* 0 = off, 1 = hack to map input neutrals to output K only, 2 = all to K */
	rspl *pcs2k;	/* PCS to K lookup for nhack */
	int wphack;		/* 0 = off, 1 = hack to map input wp to output wp, 2 = to hwp[] */
	double hwp[3];	/* hack destination white point in PCS space */
	int wphacked;	/* Opperation flag, set nz if white point was translated */

	icmFile *abs_fp;	/* Abstract profile transform */
	icRenderingIntent abs_intent;
	icc *abs_icc;
	xicc *abs_xicc;
	icxLuBase *abs_luo;	/* NULL if none */

	icxGMappingIntent gmi;
	gammap *map;	/* Gamut mapping */
	
	/* Per profile setup information */
	profinfo in;
	profinfo out;

}; typedef struct _link link;


/* ------------------------------------------- */
/* Functions called back in setting up the transform table */

/* Input table, DevIn -> DevIn' */
void devi_devip(void *cntx, double *out, double *in) {
	int rv = 0;
	link *p = (link *)cntx;

#ifdef DEBUG
	printf("IP DEVIN %f %f %f\n",in[0], in[1], in[2]); fflush(stdout);
#endif

	if (p->in.nocurve) {	/* Use created linearisation curve */
		int i;
#ifdef XYZCURVE
		p->in.devlin->lin(p->in.devlin, out, in);		/* Linearisation curve */
#else
		for (i = 0; i < p->in.chan; i++)
			out[i] = in[i];
#endif
	} else {	/* We use an explicit curve */
		switch(p->in.alg) {
		    case icmMonoFwdType: {
				icxLuMono *lu = (icxLuMono *)p->in.luo;	/* Safe to coerce */
				rv |= lu->fwd_curve(lu, out, in);
				break;
			}
		    case icmMatrixFwdType: {
				icxLuMatrix *lu = (icxLuMatrix *)p->in.luo;	/* Safe to coerce */
				rv |= lu->fwd_curve(lu, out, in);
				break;
			}
		    case icmLutType: {
				icxLuLut *lu = (icxLuLut *)p->in.luo;	/* Safe to coerce */
				/* Since not PCS, in_abs and matrix cannot be valid, */
				/* so input curve on own is ok to use. */
				rv |= lu->input(lu, out, in);
				break;
			}
			default:
				error("Unexpected algorithm type %d in devi_devip()",p->in.alg);
		}
		if (rv >= 2)
			error("icc lookup failed: %d, %s",p->in.c->errc,p->in.c->err);
	}
#ifdef DEBUG
	printf("IP DEVOUT %f %f %f\n",out[0], out[1], out[2]); fflush(stdout);
#endif
}


/* - - - - - - - - - - - - */
/* clut, DevIn' -> DevOut' */
void devip_devop(void *cntx, double *out, double *in) {
	double pcsv[MAX_CHAN];			/* PCS intermediate value */
	double locus[MAX_CHAN];			/* Auxiliary locus values */
	double nin[3];		/* Diagnostic values */
	int wptrig = 0;		/* White point hack triggered */
	int ntrig = 0;		/* Neutral K only hack triggered */
	int rv = 0;
	link *p = (link *)cntx;

#ifdef DEBUG
	printf("DEVIN %f %f %f %f\n",in[0], in[1], in[2], in[3]); fflush(stdout);
#endif

	/* Handle neutral recognition/output K only hack */
	if (p->nhack == 1) {
		double thr = (KHACKWIDTH + 0.5)/(p->clutres-1.0); /* Match threshold */

		/* We want to see if the input colors are equal (Or a=b= 0.0 ??) */
		/* li.nhack should have set p->in.nocurve, so we should be getting raw */
		/* input space device values here. It also made sure that there are at */
		/* least 3 input channels. */

		if (fabs(in[0] - in[1]) < thr
		 && fabs(in[1] - in[2]) < thr
		 && fabs(in[2] - in[0]) < thr) {
			ntrig = 1;			/* Is neutral flag */
			nin[0] = in[0];		/* Save diagnostic values */
			nin[1] = in[1];
			nin[2] = in[2];
		}
	}

	/* Do DevIn' -> PCS */
	switch(p->in.alg) {
	    case icmMonoFwdType: {
			icxLuMono *lu = (icxLuMono *)p->in.luo;	/* Safe to coerce */

			if (p->in.nocurve) {	/* No explicit curve, so do implicit here */
#ifdef XYZCURVE
				p->in.devlin->invlin(p->in.devlin, pcsv, in);
				rv |= lu->fwd_curve(lu, pcsv, pcsv);
#else
				rv |= lu->fwd_curve(lu, pcsv, in);
#endif
				rv |= lu->fwd_map(lu, pcsv, pcsv);
			} else {
				rv |= lu->fwd_map(lu, pcsv, in);
			}
			rv |= lu->fwd_abs(lu, pcsv, pcsv);
			break;
		}
	    case icmMatrixFwdType: {
			icxLuMatrix *lu = (icxLuMatrix *)p->in.luo;	/* Safe to coerce */

			if (p->in.nocurve) {	/* No explicit curve, so do implicit here */
#ifdef XYZCURVE
				p->in.devlin->invlin(p->in.devlin, pcsv, in);
				rv |= lu->fwd_curve(lu, pcsv, pcsv);
#else
				rv |= lu->fwd_curve(lu, pcsv, in);
#endif
				rv |= lu->fwd_matrix(lu, pcsv, pcsv);
			} else {
				rv |= lu->fwd_matrix(lu, pcsv, in);
			}
			rv |= lu->fwd_abs(lu, pcsv, pcsv);
			break;
		}
	    case icmLutType: {
			icxLuLut *lu = (icxLuLut *)p->in.luo;	/* Safe to coerce */
			if (p->in.nocurve) {	/* No explicit curve, so do implicit here */
				/* Since not PCS, in_abs and matrix cannot be valid, */
				/* so input curve on won is ok to use. */
#ifdef XYZCURVE
				p->in.devlin->invlin(p->in.devlin, pcsv, in);
				rv |= lu->input(lu, pcsv, pcsv);
#else
				rv |= lu->input(lu, pcsv, in);
#endif
				rv |= lu->clut(lu, pcsv, pcsv);
			} else {
				rv |= lu->clut(lu, pcsv, in);
			}
			if (p->out.inking == 0) {
				lu->clut_locus(lu, locus, pcsv, in);	/* Compute possible locus values */
#ifdef DEBUG
printf("Got possible locus value %f %f %f %f\n",locus[0],locus[1],locus[2],locus[3]);
#endif
			}
			rv |= lu->output(lu, pcsv, pcsv);
			rv |= lu->out_abs(lu, pcsv, pcsv);
			break;
		}
		default:
			error("Unexpected algorithm type %d in devip of devip_devop()",p->in.alg);
	}

	/* At this point, the PCS is:
	 *
	 *	If not gamut mapped:
	 *		Lab in the intent selected for the source profile
	 *	If gamut mapped:
	 *		either
	 *			Absolute Lab
	 *		or
	 *			Jab derived from absolute XYZ via the in/out viewing conditions
	 */

#ifdef DEBUG
	printf("PCS before map %f %f %f\n",pcsv[0], pcsv[1], pcsv[2]); fflush(stdout);
#endif

	if (p->wphack) {
		int e;
		double dd = 0.0;
		for (e = 0; e < 3; e++) {
			double tt;
			tt = pcsv[e] - p->in.wp[e];
			dd += tt * tt;
		}
		dd = sqrt(dd);

		if (dd < 1.0) {		/* Triggered withing 1 delta E */
			p->wphacked++;
			wptrig = 1;
			if (p->wphack == 2) {
				for (e = 0; e < 3; e++)		/* Map input white to given white */
					pcsv[e] = p->hwp[e];
			} else {
				for (e = 0; e < 3; e++)		/* Map input white to output white */
					pcsv[e] = p->out.wp[e];
			}

#ifndef DEBUG
			if (p->verb)
#endif
			{
				printf("White point hack mapped %f %f %f to %f %f %f, hit withing %f\n",
    		    p->in.wp[0],p->in.wp[1],p->in.wp[2],pcsv[0], pcsv[1], pcsv[2],dd);
				fflush(stdout);
			}
		}
	}

	/* Do gamut mapping */
	if (wptrig == 0 && p->mode > 0 && p->gmi.usemap) {
		
		/* We've used pcsor to ensure PCS space is appropriate */
		p->map->domap(p->map, pcsv, pcsv);
	}

#ifdef DEBUG
	printf("PCS after map %f %f %f\n",pcsv[0], pcsv[1], pcsv[2]); fflush(stdout);
#endif

	/* Abstract profile transform, PCS -> PCS */
	/* pcsor -> abstract ->pcsor conversion */
	if (p->abs_luo != NULL) {
		/* Abstract profile is either absolute or relative.  */
		/* We need to convert the current PCS into something compatible. */
		/* This is more ugly than it really should be. */
		p->abs_luo->lookup(p->abs_luo, pcsv, pcsv);
#ifdef DEBUG
	printf("PCS after abstract %f %f %f\n",pcsv[0], pcsv[1], pcsv[2]); fflush(stdout);
#endif
	}

	/* Do PCS -> DevOut' */
	/* Neutral to K only hack has triggered */
	if (p->nhack == 2
	 || (p->nhack && ntrig)) {
		co pp;
		pp.p[0] = pcsv[0];		/* Input L value */
		p->pcs2k->interp(p->pcs2k, &pp);
		out[0] = out[1] = out[2] = 0.0;		/* We know output is CMYK */
		out[3] = pp.v[0];
		if (out[3] < 0.0)		/* rspl might have extrapolated */
			out[3] = 0.0;
		else if (out[3] > 1.0)
			out[3] = 1.0;

#ifndef DEBUG
		if (p->verb)
#endif
		if (ntrig) {
			printf("Neutral hack mapped %f %f %f to 0 0 0 %f\n",
   		    nin[0], nin[1], nin[2], out[3]); 
			fflush(stdout);
		}
		
	} else {
		switch(p->out.alg) {
		    case icmMonoBwdType: {
				icxLuMono *lu = (icxLuMono *)p->out.luo;	/* Safe to coerce */

				rv |= lu->bwd_abs(lu, pcsv, pcsv);
				rv |= lu->bwd_map(lu, out, pcsv);
				if (p->out.nocurve) {	/* No explicit curve, so do implicit here */
					rv |= lu->bwd_curve(lu, out, out);
#ifdef XYZCURVE
					p->out.devlin->invlin(p->out.devlin, out, out);
#endif
				}
				break;
			}
		    case icmMatrixBwdType: {
				icxLuMatrix *lu = (icxLuMatrix *)p->out.luo;	/* Safe to coerce */

				rv |= lu->bwd_abs(lu, pcsv, pcsv);
				rv |= lu->bwd_matrix(lu, out, pcsv);
				if (p->out.nocurve) {	/* No explicit curve, so do implicit here */
					rv |= lu->bwd_curve(lu, out, out);
#ifdef XYZCURVE
					p->out.devlin->invlin(p->out.devlin, out, out);
#endif
				}
				break;
			}
		    case icmLutType: {
				icxLuLut *lu = (icxLuLut *)p->out.luo;	/* Safe to coerce */

				if (p->mode < 2) {	/* Using B2A table */
					rv |= lu->in_abs(lu, pcsv, pcsv);
					rv |= lu->matrix(lu, pcsv, pcsv);
					rv |= lu->input(lu, pcsv, pcsv);
					rv |= lu->clut(lu, out, pcsv);
					if (p->out.nocurve) {	/* No explicit curve, so do implicit here */
						rv |= lu->output(lu, out, out);
#ifdef XYZCURVE
						p->out.devlin->invlin(p->out.devlin, out, out);
#endif
					}

				} else {	/* Use inverse A2B table */
					int i;
#ifdef USE_MERGE_CLUT_OPT
					/* Because we have used the ICX_MERGE_CLUT flag, we don't need */
					/* to call inv_out_abs() and inv_output() */
#else
					rv |= lu->inv_out_abs(lu, pcsv, pcsv);
					rv |= lu->inv_output(lu, pcsv, pcsv);
#endif

#ifdef DEBUG
printf("Calling inv_clut with locus %f %f %f %f and pcsv %f %f %f %f\n",
locus[0],locus[1],locus[2],locus[3],pcsv[0],pcsv[1],pcsv[2],pcsv[3]);
#endif

					rv |= lu->inv_clut(lu, locus, pcsv);
#ifdef DEBUG
printf("Got result %f %f %f %f\n", locus[0],locus[1],locus[2],locus[3]);
#endif

					/* Use locus as out[], in case we have auxiliary targets to give to inv_clut */
					for (i = 0; i < p->out.chan; i++)
						out[i] = locus[i];		/* Copy result */
					
					if (p->out.nocurve) {	/* No explicit curve, so do implicit here */
						rv |= lu->inv_input(lu, out, out);
#ifdef XYZCURVE
						p->out.devlin->invlin(p->out.devlin, out, out);
#endif
					}
				}
				break;
			}

			default:
				error("Unexpected algorithm type %d in devop of devip_devop()",p->out.alg);
		}
		if (rv >= 2)
			error("icc lookup failed: %d, %s",p->in.c->errc,p->in.c->err);
	}

	if (p->verb) {		/* Output percent intervals */
		int pc;
		p->count++;
		pc = p->count * 100.0/p->total + 0.5;
		if (pc != p->last) {
			printf("\r%2d%%",pc); fflush(stdout);
			p->last = pc;
		}
	}

#ifdef DEBUG
	printf("DEVOUT %f %f %f %f\n\n",out[0], out[1], out[2], out[3]); fflush(stdout);
#endif
}

/* - - - - - - - - - - - - - - - - */
/* Output table, DevOut' -> DevOut */
void devop_devo(void *cntx, double *out, double *in) {
	int rv = 0;
	link *p = (link *)cntx;

#ifdef DEBUG
	printf("OP DEVIN %f %f %f %f\n",in[0], in[1], in[2], in[4]); fflush(stdout);
#endif

	if (p->out.nocurve) {	/* Use linear curve */
		int i;
#ifdef XYZCURVE
		p->out.devlin->lin(p->out.devlin, out, in);		/* Generated linearisation curve */
#else
		for (i = 0; i < p->out.chan; i++)
			out[i] = in[i];
#endif

	} else {	/* We use an explicit curve */

		switch(p->out.alg) {
		    case icmMonoBwdType: {
				icxLuMono *lu = (icxLuMono *)p->out.luo;		/* Safe to coerce */
				rv |= lu->bwd_curve(lu, out, in);
				break;
			}
		    case icmMatrixBwdType: {
				icxLuMatrix *lu = (icxLuMatrix *)p->out.luo;	/* Safe to coerce */
				rv |= lu->bwd_curve(lu, out, in);
				break;
			}
		    case icmLutType: {
				if (p->mode < 2) {	/* Using B2A table */
					icxLuLut *lu = (icxLuLut *)p->out.luo;	/* Safe to coerce */
					rv |= lu->output(lu, out, in);
					/* Since not PCS, out_abs is never used */
					break;
				} else {	/* Use inverse A2B table */
					icxLuLut *lu = (icxLuLut *)p->out.luo;	/* Safe to coerce */
					rv |= lu->inv_input(lu, out, in);
					/* Since not PCS, inv_matrix and inv_in_abs is never used */
					break;
				}
			}
			default:
				error("Unexpected algorithm type in devop_devo()");
		}
		if (rv >= 2)
			error("icc lookup failed: %d, %s",p->in.c->errc,p->in.c->err);
	}
#ifdef DEBUG
	printf("OP DEVOUT %f %f %f %f\n",out[0], out[1], out[2], out[3]); fflush(stdout);
#endif
}

/* ------------------------------------------- */
/* Fixup L -> K only lookup table white and black values */

/* Context for fixup */
typedef struct {
	double kmax;
	double kmin;
} pcs2k_ctx;

/* Function to pass to rspl to re-set output values, */
/* to make them relative to the white and black points */
static void
fix_pcs2k_white(
	void *pp,			/* relativectx structure */
	double *out,		/* output value */
	double *in			/* input value */
) {
	double f;
	pcs2k_ctx *p = (pcs2k_ctx *)pp;

	/* Set white K value to zero, without affecting K max */
	f = 1.0/(p->kmax - p->kmin) * (out[0] - p->kmin);

	out[0] = f;
}

/* ------------------------------------------- */

int
main(int argc, char *argv[]) {
	int fa, nfa;				/* argument we're looking at */
	char in_name[MAXNAMEL+1];
	char sgam_name[MAXNAMEL+1] = "\000";	/* Source gamut name */
	char abs_name[MAXNAMEL+1] = "\000";		/* Abstract profile name */
	char out_name[MAXNAMEL+1];
	char link_name[MAXNAMEL+1];
	int verify = 0;				/* Do verify pass */
	int rv = 0;
	icxViewCond ivc, ovc;		/* Viewing Condition Overrides for in and out profiles */
	int ivc_e = -1, ovc_e = -1;	/* Enumerated viewing condition */
	link li;					/* Linking information structure */
	int i;

	error_program = argv[0];

	/* Set defaults */
	li.verb = 0;
	li.count = 0;
	li.last = -1;
	li.mode = 0;
	li.quality = 1;				/* Medium */
	li.clutres = 0;				/* No override */	
	li.nhack = 0;
	li.pcs2k = NULL;
	li.wphack = 0;
	li.wphacked = 0;
	li.abs_luo = NULL;				/* No abstract */
	li.hwp[0] = li.hwp[1] = li.hwp[2] = 0.0;
	li.map = NULL;
	li.in.intent  = icmDefaultIntent;	/* Default */
	li.in.inking  = 2;					/* Inking algorithm default */
	li.in.nocurve = 0;					/* Preserve device linearisation curve */
	li.in.devlin = NULL;				/* XYZ PCS device linearisation curve */
	li.out.intent = icmDefaultIntent;	/* Default */
	li.out.ink.tlimit = -1.0;			/* Default no total limit */
	li.out.ink.klimit = -1.0;			/* Default no black limit */
	li.out.inking  = 2;					/* Default 0.5 K */
	li.out.nocurve = 0;					/* Preserve device linearisation curve */
	li.out.devlin = NULL;				/* XYZ PCS device linearisation curve */

	xicc_enum_gmapintent(&li.gmi, -1);	/* Set default overall intent */

	/* Init VC overrides so that we know when the've been set */
	ivc.Ev = -1;
	ivc.Wxyz[0] = -1.0; ivc.Wxyz[1] = -1.0; ivc.Wxyz[2] = -1.0;
	ivc.Yb = -1.0;
	ivc.La = -1.0;
	ivc.Lv = -1.0;
	ivc.Yf = -1.0;
	ivc.Fxyz[0] = -1.0; ivc.Fxyz[1] = -1.0; ivc.Fxyz[2] = -1.0;

	ovc.Ev = -1;
	ovc.Wxyz[0] = -1.0; ovc.Wxyz[1] = -1.0; ovc.Wxyz[2] = -1.0;
	ovc.Yb = -1.0;
	ovc.La = -1.0;
	ovc.Lv = -1.0;
	ovc.Yf = -1.0;
	ovc.Fxyz[0] = -1.0; ovc.Fxyz[1] = -1.0; ovc.Fxyz[2] = -1.0;

	if (argc < 4)
		usage("Too few arguments, got %d expect at least 3",argc-1);

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
			else if (argv[fa][1] == 'v') {
				li.verb = 1;
			}

			/* Verify rather than link */
			else if (argv[fa][1] == 'V')
				verify = 1;

			/* Disable linearisation curve preservation */
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N') {
				li.in.nocurve = 1;
				li.out.nocurve = 1;
			}

			/* Hack to force input neutrals to K only output */
			else if (argv[fa][1] == 'f') {
				li.nhack = 1;
				li.in.nocurve = 1;		/* Disable input curve to preserve input equality */
			}

			/* Hack to force all input colors to K only output */
			else if (argv[fa][1] == 'F') {
				li.nhack = 2;
			}

			/* Quality */
			else if (argv[fa][1] == 'q' || argv[fa][1] == 'Q') {
				fa = nfa;
				if (na == NULL) usage("Quality flag (-q) needs an argument");
    			switch (na[0]) {
					case 'l':
					case 'L':
						li.quality = 0;
						break;
					case 'm':
					case 'M':
						li.quality = 1;
						break;
					case 'h':
					case 'H':
						li.quality = 2;
						break;
					case 'u':
					case 'U':
						li.quality = 3;
						break;
					default:
						usage("Unrecognised quality flag (-q) argument '%c'",na[0]);
				}
			}

			/* CLUT resolution override */
			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R') {
				int rr;
				fa = nfa;
				if (na == NULL) usage("Resolution flag (-r) needs an argument");
				rr = atoi(na);
				if (rr < 1 || rr > 100) usage("Resolution flag (-r) argument out of range (%d)",rr);
				li.clutres = rr;
			}

			/* Abstract profile */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				if (na == NULL) usage("Expected abstract profile filename after -a");
				fa = nfa;
				strncpy(abs_name,na,MAXNAMEL); abs_name[MAXNAMEL] = '\000';
			}

			/* Simple mode */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				li.mode = 0;
			}

			/* Maping mode */
			else if (argv[fa][1] == 'g' || argv[fa][1] == 'G') {
				li.mode = 1;
				if (argv[fa][1] == 'G') {
					li.mode = 2;
				}

				if (na != NULL) {	/* Found an optional source gamut */
					fa = nfa;
					strncpy(sgam_name,na,MAXNAMEL); sgam_name[MAXNAMEL] = '\000';
				}

			}
			/* White point hack */
			else if (argv[fa][1] == 'w' || argv[fa][1] == 'W') {
				li.wphack = 1;
				if (na != NULL) {		// To a particular white point
					fa = nfa;
					if (sscanf(na, " %lf , %lf , %lf ",&li.hwp[0], &li.hwp[1], &li.hwp[2]) == 3) {
						li.wphack = 2;
					} else
						usage("Couldn't parse hack white point (-w) value '%s'",na);
				}
				if (li.mode < 1)	/* Set minimum link mode */
					li.mode = 1;
			}
			/* Input profile Intent or Mapping mode intent */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL) usage("Input intent flag (-i) needs an argument");
    			switch (na[0]) {
					case 'p':
					case 'P':
						li.in.intent = icPerceptual;
						xicc_enum_gmapintent(&li.gmi, icxPerceptualGMIntent);
						break;
					case 'r':
					case 'R':
						li.in.intent = icRelativeColorimetric;
						xicc_enum_gmapintent(&li.gmi, icxRelativeGMIntent);
						break;
					case 's':
					case 'S':
						li.in.intent = icSaturation;
						xicc_enum_gmapintent(&li.gmi, icxSaturationGMIntent);
						break;
					case 'a':
					case 'A':
						li.in.intent = icAbsoluteColorimetric;
						xicc_enum_gmapintent(&li.gmi, icxAbsoluteGMIntent);
						break;
					default:
						if (na[0] >= '0' && na[0] <= '9') {
							int i = atoi(na);
							if (xicc_enum_gmapintent(&li.gmi, i) == 0) {
								if (li.mode < 1)	/* Set minimum link mode */
									li.mode = 1;
								break;
							}
						}
						usage("Input intent (-i) argument '%s' isn't recognised",na);
				}
			}

			/* Output profile Intent */
			else if (argv[fa][1] == 'o' || argv[fa][1] == 'O') {
				fa = nfa;
				if (na == NULL) usage("Output intent flag (-o) needs an argument");
    			switch (na[0]) {
					case 'p':
					case 'P':
						li.out.intent = icPerceptual;
						break;
					case 'r':
					case 'R':
						li.out.intent = icRelativeColorimetric;
						break;
					case 's':
					case 'S':
						li.out.intent = icSaturation;
						break;
					case 'a':
					case 'A':
						li.out.intent = icAbsoluteColorimetric;
						break;
					default:
						usage("Output intent (-o) argument '%s' not recognised",na);
				}
			}
			
			/* Viewing conditions */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C'
			      || argv[fa][1] == 'd' || argv[fa][1] == 'D') {
				icxViewCond *vc;

				if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
					vc = &ivc;
				} else {
					vc = &ovc;
				}

				fa = nfa;
				if (na == NULL) usage("Viewing conditions flag (-[cd]) needs an argument");
				if (na[0] >= '0' && na[0] <= '9') {
					if (vc == &ivc)
						ivc_e = atoi(na);
					else
						ovc_e = atoi(na);
				} else if (na[0] == 's' || na[0] == 'S') {
					if (na[1] != ':')
						usage("Viewing conditions (-[cd]s) missing ':'");
					if (na[2] == 'a' || na[2] == 'A') {
						vc->Ev = vc_average;
					} else if (na[2] == 'm' || na[2] == 'M') {
						vc->Ev = vc_dim;
					} else if (na[2] == 'd' || na[2] == 'D') {
						vc->Ev = vc_dark;
					} else if (na[2] == 'c' || na[2] == 'C') {
						vc->Ev = vc_cut_sheet;
					} else
						usage("Viewing condition (-[cd]) unrecognised surround '%c'",na[2]);
				} else if (na[0] == 'w' || na[0] == 'W') {
					double x, y, z;
					if (sscanf(na+1,":%lf:%lf:%lf",&x,&y,&z) == 3) {
						vc->Wxyz[0] = x; vc->Wxyz[1] = y; vc->Wxyz[2] = z;
					} else if (sscanf(na+1,":%lf:%lf",&x,&y) == 2) {
						vc->Wxyz[0] = x; vc->Wxyz[1] = y;
					} else
						usage("Viewing condition (-[cd]w) unrecognised white point '%s'",na+1);
				} else if (na[0] == 'a' || na[0] == 'A') {
					if (na[1] != ':')
						usage("Viewing conditions (-[cd]a) missing ':'");
					vc->La = atof(na+2);
				} else if (na[0] == 'b' || na[0] == 'B') {
					if (na[1] != ':')
						usage("Viewing conditions (-[cd]b) missing ':'");
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
						usage("Viewing condition (-[cd]f) unrecognised flare '%s'",na+1);
				} else
					usage("Viewing condition (-[cd]) unrecognised sub flag '%c'",na[0]);
				if (li.mode < 1)	/* Set minimum link mode */
					li.mode = 1;
			}

			/* Inking rule */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				fa = nfa;
				if (na == NULL) usage("Inking rule flag (-k) needs an argument");
    			switch (na[0]) {
					case 't':
					case 'T':
						li.out.inking = 0;		/* Use existing table as guide */
						break;
					case 'z':
					case 'Z':
						li.out.inking = 1;		/* Use minimum k */
						break;
					case 'h':
					case 'H':
						li.out.inking = 2;		/* Use half k */
						break;
					case 'x':
					case 'X':
						li.out.inking = 3;		/* Use maximum k */
						break;
					case 'r':
					case 'R':
						li.out.inking = 4;		/* Use ramp k */
						break;
					case 'p':
					case 'P':
					case 'q':
					case 'Q':
						li.out.inking = 5;		/* Use curve parameter */

						++fa;
						if (fa >= argc) usage("Inking rule (-kp) expects more parameters");
						li.out.ink.c.Kstle = atof(argv[fa]);

						++fa;
						if (fa >= argc) usage("Inking rule (-kp) expects more parameters");
						li.out.ink.c.Kstpo = atof(argv[fa]);

						++fa;
						if (fa >= argc || argv[fa][0] == '-') usage("Inking rule (-kp) expects more parameters");
						li.out.ink.c.Kenpo = atof(argv[fa]);

						++fa;
						if (fa >= argc || argv[fa][0] == '-') usage("Inking rule (-kp) expects more parameters");
						li.out.ink.c.Kenle = atof(argv[fa]);

						++fa;
						if (fa >= argc || argv[fa][0] == '-') usage("Inking rule (-kp) expects more parameters");
						li.out.ink.c.Kshap = atof(argv[fa]);

						if (na[0] == 'q' || na[0] == 'Q') {
							li.out.inking = 6;		/* Use transfer to dual curve parameter */

							++fa;
							if (fa >= argc) usage("Inking rule (-kq) expects more parameters");
							li.out.ink.x.Kstle = atof(argv[fa]);
	
							++fa;
							if (fa >= argc) usage("Inking rule (-kq) expects more parameters");
							li.out.ink.x.Kstpo = atof(argv[fa]);
	
							++fa;
							if (fa >= argc || argv[fa][0] == '-') usage("Inking rule (-kq) expects more parameters");
							li.out.ink.x.Kenpo = atof(argv[fa]);
	
							++fa;
							if (fa >= argc) usage("Inking rule (-kq) expects more parameters");
							li.out.ink.x.Kenle = atof(argv[fa]);
	
							++fa;
							if (fa >= argc || argv[fa][0] == '-') usage("Inking rule (-kq) expects more parameters");
							li.out.ink.x.Kshap = atof(argv[fa]);
	
						}
						break;
					default:
						usage("Inking rule (-k) unknown sub flag '%c'",na[0]);
				}
				if (li.mode < 2)	/* Set minimum link mode */
					li.mode = 2;
			}
			/* Ink limits */
			else if (argv[fa][1] == 'l') {
				int tlimit;
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -l");
				tlimit = atoi(na);
				if (tlimit >= 0)
					li.out.ink.tlimit = tlimit/100.0;
				else
					li.out.ink.tlimit = -1.0;
				if (li.mode < 2)	/* Set minimum link mode */
					li.mode = 2;
			}
			else if (argv[fa][1] == 'L') {
				int klimit;
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -L");
				klimit = atoi(na);
				if (klimit >= 0)
					li.out.ink.klimit = klimit/100.0;
				else
					li.out.ink.klimit = -1.0;
				if (li.mode < 2)	/* Set minimum link mode */
					li.mode = 2;
			}

			else 
				usage("Unknown flag '%c'",argv[fa][1]);
		} else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage("Missing input profile");
	strncpy(in_name,argv[fa++],MAXNAMEL); in_name[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage("Missing output profile");
	strncpy(out_name,argv[fa++],MAXNAMEL); out_name[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage("Missing result profile");
	strncpy(link_name,argv[fa++],MAXNAMEL); link_name[MAXNAMEL] = '\000';

	if (li.verb)
		printf("Got options\n");

	/* - - - - - - - - - - - - - - - - - - - */
	/* Open up the input device profile for reading */
	if ((li.in.fp = new_icmFileStd_name(in_name,"r")) == NULL)
		error ("Can't open file '%s'",in_name);
	
	if ((li.in.c = new_icc()) == NULL)
		error ("Creation of Input profile ICC object failed");

	/* Read header etc. */
	if ((rv = li.in.c->read(li.in.c,li.in.fp,0)) != 0)
		error ("%d, %s",rv,li.in.c->err);
	li.in.h = li.in.c->header;

	/* Check that it is a suitable device input icc */
	if (li.in.h->deviceClass != icSigInputClass
	 && li.in.h->deviceClass != icSigDisplayClass
	 && li.in.h->deviceClass != icSigOutputClass
	 && li.in.h->deviceClass != icSigColorSpaceClass)	/* For sRGB etc. */
		error("Input profile isn't a device profile");

	/* Wrap with an expanded icc */
	if ((li.in.x = new_xicc(li.in.c)) == NULL)
		error ("Creation of input profile xicc failed");

	/* - - - - - - - - - - - - - - - - - - - */
	/* Open up the abstract profile if requested */
	if (abs_name[0] != '\000') {
		if ((li.abs_fp = new_icmFileStd_name(abs_name,"r")) == NULL)
			error ("Can't open abstract profile file '%s'",abs_name);
		
		if ((li.abs_icc = new_icc()) == NULL)
			error ("Creation of Abstract profile ICC object failed");

		/* Read header etc. */
		if ((rv = li.abs_icc->read(li.abs_icc,li.abs_fp,0)) != 0)
			error ("%d, %s",rv,li.abs_icc->err);
		li.abs_icc->header;

		if (li.abs_icc->header->deviceClass != icSigAbstractClass)
			error("Abstract profile isn't an abstract profile");

		/* Take intended abstract intent from profile itself */
		if ((li.abs_intent = li.abs_icc->header->renderingIntent) != icAbsoluteColorimetric)
			li.abs_intent = icRelativeColorimetric;

		/* Wrap with an expanded icc */
		if ((li.abs_xicc = new_xicc(li.abs_icc)) == NULL)
			error ("Creation of abstract profile xicc failed");
	}
	/* - - - - - - - - - - - - - - - - - - - */
	/* Open up the output device output profile for reading */
	if ((li.out.fp = new_icmFileStd_name(out_name,"r")) == NULL)
		error ("Can't open file '%s'",out_name);
	
	if ((li.out.c = new_icc()) == NULL)
		error ("Creation of Output profile ICC object failed");

	/* Read header etc. */
	if ((rv = li.out.c->read(li.out.c,li.out.fp,0)) != 0)
		error ("%d, %s",rv,li.out.c->err);
	li.out.h = li.out.c->header;

	if (li.out.h->deviceClass != icSigInputClass
	 && li.out.h->deviceClass != icSigDisplayClass
	 && li.out.h->deviceClass != icSigOutputClass
	 && li.out.h->deviceClass != icSigColorSpaceClass)	/* For sRGB etc. */
		error("Output profile isn't a device profile");

	/* Wrap with an expanded icc */
	if ((li.out.x = new_xicc(li.out.c)) == NULL)
		error ("Creation of output profile xicc failed");

	/* - - - - - - - - - - - - - - - - - - - */
	/* Deal with setting up all the options */

	/* deal with output black generation. */
	/* Ink limits will have been set in option parsing */

	li.out.ink.c.Ksmth = ICXINKDEFSMTH;		/* default curve smoothing */
	li.out.ink.x.Ksmth = ICXINKDEFSMTH;

	switch (li.out.inking) {
		case 0:			/* Use input profile K locus */
			/* Sanity check */
			if (li.in.h->colorSpace != li.in.h->colorSpace)
				error("Can't transfer black ink in & out unless the same colorspaces");
			li.out.ink.k_rule = icxKlocus;		/* Given as aux parameter in PCS -> Device */
			break;
		case 1:			/* Minimum K */
			li.out.ink.k_rule = icxKluma5;
			li.out.ink.c.Kstle = 0.0;
			li.out.ink.c.Kstpo = 0.0;
			li.out.ink.c.Kenpo = 1.0;
			li.out.ink.c.Kenle = 0.0;
			li.out.ink.c.Kshap = 1.0;
			break;
		case 2:			/* 0.5 K */
			li.out.ink.k_rule = icxKluma5;
			li.out.ink.c.Kstle = 0.5;
			li.out.ink.c.Kstpo = 0.0;
			li.out.ink.c.Kenpo = 1.0;
			li.out.ink.c.Kenle = 0.5;
			li.out.ink.c.Kshap = 1.0;
			break;
		case 3:			/* Maximum K */
			li.out.ink.k_rule = icxKluma5;
			li.out.ink.c.Kstle = 1.0;
			li.out.ink.c.Kstpo = 0.0;
			li.out.ink.c.Kenpo = 1.0;
			li.out.ink.c.Kenle = 1.0;
			li.out.ink.c.Kshap = 1.0;
			break;
		case 4:			/* Ramp K */
			li.out.ink.k_rule = icxKluma5;
			li.out.ink.c.Kstle = 0.0;
			li.out.ink.c.Kstpo = 0.0;
			li.out.ink.c.Kenpo = 1.0;
			li.out.ink.c.Kenle = 1.0;
			li.out.ink.c.Kshap = 1.0;
			break;
		case 5:			/* Curve */
			li.out.ink.k_rule = icxKluma5;
			break;		/* Other params already set by options */
		case 6:			/* Use input profile K locus + dual curve limits */
			/* Sanity check */
			if (li.in.h->colorSpace != li.in.h->colorSpace)
				error("Can't transfer black ink in & out unless the same colorspaces");
			li.out.ink.k_rule = icxKl5l;		/* Given as aux parameter in PCS -> Device */
			break;		/* Other params already set by options */
	}

	for (i = 0; i < 2; i++) {
		xicc *x;
		icxViewCond *v, *vc;
		int es;

		if (i == 0) {
			v = &ivc;			/* Override parameters */
			vc = &li.in.vc;		/* Target parameters */
			es = ivc_e;
			x = li.in.x;		/* xicc */
		} else {
			v = &ovc;			/* Override parameters */
			vc = &li.out.vc;	/* Target parameters */
			es = ovc_e;
			x = li.out.x;		/* xicc */
		}
		
		/* Set the default */
		xicc_enum_viewcond(x, vc, -1, 0);

		/* Override the viewing conditions */
		/* (?? Could move this code into xicc_enum_viewcond() as an option ??) */
		if (es >= 0)
			if (xicc_enum_viewcond(x, vc, es, 0))
				error ("%d, %s",x->errc, x->err);
		if (v->Ev >= 0)
			vc->Ev = v->Ev;
		if (v->Wxyz[0] >= 0.0 && v->Wxyz[1] > 0.0 && v->Wxyz[2] >= 0.0) {
			/* Normalise XYZ to current media white */
			vc->Wxyz[0] = v->Wxyz[0]/v->Wxyz[1] * vc->Wxyz[1];
			vc->Wxyz[2] = v->Wxyz[2]/v->Wxyz[1] * vc->Wxyz[1];
		} 
		if (v->Wxyz[0] >= 0.0 && v->Wxyz[1] >= 0.0 && v->Wxyz[2] < 0.0) {
			/* Convert Yxy to XYZ */
			double x = v->Wxyz[0];
			double y = v->Wxyz[1];	/* If Y == 1.0, then X+Y+Z = 1/y */
			double z = 1.0 - x - y;
			vc->Wxyz[0] = x/y * vc->Wxyz[1];
			vc->Wxyz[2] = z/y * vc->Wxyz[1];
		}
		if (v->La >= 0.0)
			vc->La = v->La;
		if (v->Yb >= 0.0)
			vc->Yb = v->Yb;
		if (v->Yf >= 0.0)
			vc->Yf = vc->Yf;
		if (v->Fxyz[0] >= 0.0 && v->Fxyz[1] > 0.0 && v->Fxyz[2] >= 0.0) {
			/* Normalise XYZ to current media white */
			vc->Fxyz[0] = v->Fxyz[0]/v->Fxyz[1] * vc->Fxyz[1];
			vc->Fxyz[2] = v->Fxyz[2]/v->Fxyz[1] * vc->Fxyz[1];
		}
		if (v->Fxyz[0] >= 0.0 && v->Fxyz[1] >= 0.0 && v->Fxyz[2] < 0.0) {
			/* Convert Yxy to XYZ */
			double x = v->Fxyz[0];
			double y = v->Fxyz[1];	/* If Y == 1.0, then X+Y+Z = 1/y */
			double z = 1.0 - x - y;
			vc->Fxyz[0] = x/y * vc->Fxyz[1];
			vc->Fxyz[2] = z/y * vc->Fxyz[1];
		}
	}

	/* If we are using the gamut map mode, then the */
	/* ICC profile intent we want is AbsoluteColorimetric */
	if (li.mode > 0) {
		li.in.intent = icAbsoluteColorimetric;
		li.out.intent = icAbsoluteColorimetric;
		li.abs_intent = icAbsoluteColorimetric;
	}

	if (li.verb)
		printf("Configured options\n");

	/* - - - - - - - - - - - - - - - - - - - */
	/* Setup the profile color lookup information */
	{
		icColorSpaceSignature pcsor;	/* PCS to use between in & out profiles */
		icmLuAlgType oalg;				/* Native output algorithm */
		icColorSpaceSignature natpcs;	/* Underlying native output PCS */
		
		pcsor = icSigLabData;			/* Default use Lab as PCS */

		if (li.mode > 0 && li.gmi.usecas) {
			pcsor = icxSigJabData;		/* Use CAM as PCS */

			if (li.gmi.usecas == 2) {	/* Absolute Appearance space */
				double mxw;

				/* Make absolute common white point average between the two */
				li.in.vc.Wxyz[0] = 0.5 * (li.in.vc.Wxyz[0] + li.out.vc.Wxyz[0]);
				li.in.vc.Wxyz[1] = 0.5 * (li.in.vc.Wxyz[1] + li.out.vc.Wxyz[1]);
				li.in.vc.Wxyz[2] = 0.5 * (li.in.vc.Wxyz[2] + li.out.vc.Wxyz[2]);

				/* And scale it Y to be equal to 1.0 */
				mxw = 1.0/li.in.vc.Wxyz[1];
				li.in.vc.Wxyz[0] *= mxw;
				li.in.vc.Wxyz[1] *= mxw;
				li.in.vc.Wxyz[2] *= mxw;

				/* Set the output vc to be the same as the input */
				li.out.vc = li.in.vc;	/* Structure copy */
			}
		}

		if (li.verb)
			printf("Loading input A2B table\n");

		/* Get an input profile xicc conversion object */
//xicc_dump_viewcond(&li.in.vc);

		if ((li.in.luo = li.in.x->get_luobj(li.in.x, 0, icmFwd, li.in.intent,
		                                    pcsor, icmLuOrdNorm, &li.in.vc, NULL)) == NULL) {
			error("get xlookup object failed: %d, %s",li.in.x->errc,li.in.x->err);
		}
	
		/* Get details of overall conversion */
		li.in.luo->spaces(li.in.luo, NULL, &li.in.chan, NULL, NULL, &li.in.alg,
		                  NULL, NULL, NULL);

		/* Grab the white point in case the wphack needs it */
		li.in.luo->rel_wh_bk_points(li.in.luo, li.in.wp, NULL);

		/* Get native PCS space */
		li.in.luo->lutspaces(li.in.luo, NULL, NULL, &natpcs, NULL, NULL);
	
		if (natpcs == icSigXYZData) {
			double min[MAX_CHAN], max[MAX_CHAN];
			
			li.in.luo->get_ranges(li.in.luo, min, max, NULL, NULL);	/* Get input space range */
			li.in.nocurve = 1;		/* Don't use input profiles input curve */
#ifdef XYZCURVE
			li.in.devlin = new_xdevlin(li.in.chan, min, max, li.in.luo,
			                    (void (*) (void *, double *, double *))li.in.luo->lookup); 
#endif
			
			if (li.verb)
				printf("No explicit curve on input profile\n");
		}


		/* Setup any abstract profile to match the chosen PCS */
		/* We aren't checking whether the input/abstract/output profile */
		/* intents really make any sense. It's assumed at the moment */
		/* that the user knows what they're doing! */
		if (abs_name[0] != '\000') {

			if ((li.abs_luo = li.abs_xicc->get_luobj(li.abs_xicc, 0, icmFwd, li.abs_intent,
			        pcsor, icmLuOrdNorm, &li.out.vc, NULL)) == NULL)
					error ("%d, %s",li.abs_icc->errc, li.abs_icc->err);
		}

		// Figure out whether the output profile is a Lut profile or not */
		{
			icmLuBase *plu;

			/* Get temporary icm lookup object */
			/* (Use Fwd just in case profile is missing B2A !!!!) */
			if ((plu = li.out.c->get_luobj(li.out.c, icmFwd, icmDefaultIntent, icmSigDefaultData,
			                            icmLuOrdNorm)) == NULL) {
				error("get icm lookup object failed: on '%s' %d, %s",out_name,li.out.c->errc,li.out.c->err);
			}

			/* Check what the algorithm is */
			plu->spaces(plu, NULL, NULL, NULL, NULL, &oalg, NULL, NULL, NULL);

			/* release the icm lookup */
			plu->del(plu);

		}

		if (oalg != icmLutType || li.mode < 2) {	/* Using B2A table or inv. mono/matrix */
			if (li.verb)
				printf("Loading output B2A table\n");

			if ((li.out.luo = li.out.x->get_luobj(li.out.x, 0, icmBwd, li.out.intent,
			                                      pcsor, icmLuOrdNorm, &li.out.vc, NULL)) == NULL) {
				error("get xlookup object failed: %d, %s",li.out.x->errc,li.out.x->err);
			}
		
			/* Get details of overall conversion */
			li.out.luo->spaces(li.out.luo, NULL, NULL, NULL, &li.out.chan, &li.out.alg,
			                   NULL, NULL, NULL);

			/* Grab the white point in case the wphack needs it */
			li.out.luo->rel_wh_bk_points(li.out.luo, li.out.wp, NULL);

			/* Get native PCS space */
			li.out.luo->lutspaces(li.out.luo, &natpcs, NULL, NULL, NULL, NULL);

		} else {	/* Using inverse A2B Lut for output conversion */
			if (li.verb)
				printf("Loading output inverse A2B table\n");

//xicc_dump_viewcond(&li.out.vc);

			if ((li.out.luo = li.out.x->get_luobj(li.out.x, ICX_CLIP_NEAREST
#ifdef USE_MERGE_CLUT_OPT
			                  | ICX_MERGE_CLUT
#endif
#ifdef USE_CAM_CLIP_OPT
			                  | ICX_CAM_CLIP
#endif
			                  , icmFwd, li.out.intent, pcsor, icmLuOrdNorm, &li.out.vc,
			                  &li.out.ink)) == NULL) {
				error("get xlookup object failed: %d, %s",li.out.x->errc,li.out.x->err);
			}
		
			/* Get details of overall conversion */
			li.out.luo->spaces(li.out.luo, NULL, &li.out.chan, NULL, NULL, &li.out.alg,
			                   NULL, NULL, NULL);

			/* Grab the white point in case the wphack needs it */
			li.out.luo->rel_wh_bk_points(li.out.luo, li.out.wp, NULL);

			/* Get native PCS space */
			li.out.luo->lutspaces(li.out.luo, NULL, NULL, &natpcs, NULL, NULL);
		}

		/* If we need an PCS'->K' mapping for the neutral axis to K hack */
		if (li.nhack) {
			icxLuBase *luo;		/* Base XLookup type object */
			icmLuAlgType alg;	/* Type of lookup algorithm */
			co ips[256];		/* Initialisation points */
			datai glow;			/* Grid low scale */
			datai ghigh;		/* Grid high scale */
			datao vlow;			/* Data value low normalize */
			datao vhigh;		/* Data value high normalize */
			double Lmax, Lmin;	/* Max and Min L values that result */
			int gres;

			if (li.out.h->colorSpace != icSigCmykData)
				error("Netral Axis K only requested with non CMYK output profile");

			if (li.in.chan < 3)
				error("Netral Axis K only requested with input profile with less than 3 channels");

			if ((li.pcs2k = new_rspl(1, 1)) == NULL) {
				error("Failed to create an rspl object");
			}

			/* Get a device to PCS lookup object to use to lookup K->PCS */
			if ((luo = li.out.x->get_luobj(li.out.x, 0,
			                  icmFwd, li.out.intent, pcsor, icmLuOrdNorm, &li.out.vc,
			                  NULL)) == NULL) {
				error("get xlookup object failed: %d, %s",li.out.x->errc,li.out.x->err);
			}
			/* Get details of overall conversion */
			li.out.luo->spaces(li.out.luo, NULL, NULL, NULL, NULL, &alg, NULL, NULL, NULL);
			if (alg != icmLutType)
				error ("Unexpected algorithm type for CMYK output profile");

			/* Setup the initialisation points */
			Lmax = -100.0;
			Lmin = 1000.0;
			for (i = 0; i < 256; i++) {
				icxLuLut *lu = (icxLuLut *)luo;	/* Safe to coerce */
				double in[4], pcsv[3];
				in[0] = in[1] = in[2] = 0.0;
				in[3] = i/(255.0);

				/* Want to do dev' -> PCS conversion to match */
				/* the normal inverse PCS-> dev' */
				if (li.out.nocurve) {	/* No explicit curve, so do implicit here */
					/* Since not PCS, in_abs and matrix cannot be valid, */
					/* so input curve on won is ok to use. */
#ifdef XYZCURVE
					li.out.devlin->invlin(li.out.devlin, pcsv, in);
					lu->input(lu, pcsv, pcsv);
#else
					lu->input(lu, pcsv, in);
#endif
					lu->clut(lu, pcsv, pcsv);
				} else {
					lu->clut(lu, pcsv, in);
				}
				lu->output(lu, pcsv, pcsv);
				lu->out_abs(lu, pcsv, pcsv);
		
				/* We force the rspl to be a forward conversion by swapping K and PCS */
				ips[i].p[0] = pcsv[0];		/* PCS as input */
				ips[i].v[0] = in[3];		/* K as output */
#ifdef NEUTKDEBUG
				printf("L %f -> K %f\n",pcsv[0], in[3]);
#endif /* NEUTKDEBUG */

				if (pcsv[0] > Lmax)			/* Track min and max values */
					Lmax = pcsv[0];
				if (pcsv[0] < Lmin)
					Lmin = pcsv[0];
			}

			glow[0] = 0.0;
			ghigh[0] = 100.0;
			vlow[0] = 0.0;
			vhigh[0] = 1.0;
			gres = 256;

			li.pcs2k->fit_rspl(li.pcs2k, 0, ips, 256, glow, ghigh, &gres, vlow, vhigh, 1.0, 0.005);

			/* Fixup the white point for neutral axis to K hack */
			{
				pcs2k_ctx cx;		/* White point fixup context */
				co pp;				/* Lookup the min and max K values */

				pp.p[0] = Lmax;
				li.pcs2k->interp(li.pcs2k, &pp);
				cx.kmin = pp.v[0];
				pp.p[0] = Lmin;
				li.pcs2k->interp(li.pcs2k, &pp);
				cx.kmax = pp.v[0];

#ifdef NEUTKDEBUG
				printf("Lmax %f, Lmin %f, Kmin %f, Kmax %f\n",Lmax, Lmin, cx.kmin, cx.kmax);
#endif /* NEUTKDEBUG */

				li.pcs2k->re_set_rspl(li.pcs2k, 0, (void *)&cx, fix_pcs2k_white);
			}

#ifdef NEUTKDEBUG
			{
				co pp;
			
				printf("\n");
				pp.p[0] = 94.925292;		/* Input L value */
				li.pcs2k->interp(li.pcs2k, &pp);
				printf("L %f -> K %f\n",pp.p[0], pp.v[0]);
			
				pp.p[0] = 94.925180;		/* Input L value */
				li.pcs2k->interp(li.pcs2k, &pp);
				printf("L %f -> K %f\n",pp.p[0], pp.v[0]);
			
				pp.p[0] = 17.699968;		/* Input L value */
				li.pcs2k->interp(li.pcs2k, &pp);
				printf("L %f -> K %f\n",pp.p[0], pp.v[0]);
			}
#endif /* NEUTKDEBUG */

		}	/* end if neutral axis to K hack */

		if (natpcs == icSigXYZData) {
			double min[MAX_CHAN], max[MAX_CHAN];
			
			li.out.luo->get_ranges(li.out.luo, NULL, NULL, min, max);	/* Get output space range */
			li.out.nocurve = 1;		/* Don't use output profiles output curve */
#ifdef XYZCURVE
			li.out.devlin = new_xdevlin(li.out.chan, min, max, li.out.luo,
			    (void (*) (void *, double *, double *))		/* cast */ 
				(li.mode < 2 ? li.out.luo->inv_lookup : li.out.luo->lookup));
#endif
			
			if (li.verb)
				printf("No explicit curve on output profile\n");
		}
	}

	/* - - - - - - - - - - - - - - - - - - - */
	/* Setup the gamut mapping */
// ~~~~ need to account for possible abstract profile after source !!!!
// ~~~~ also need to fix tiffgamut to allow for abstract profile !!!!

	if (li.verb)
		printf("Gamut mapping mode is '%s'\n",li.mode == 0 ? "Simple" : li.mode == 1 ? "Mapping" : "Mapping inverse A2B");
	if (li.verb && li.mode > 0)
		printf("Gamut mapping intent is '%s'\n",li.gmi.desc);

	/* In gamut mapping mode, the PCS used will always be absolute */
	/* intent from the input and output profiles, and either */
	/* lab or Jab space, with the given in/out viewing conditions */
	/* for the latter. The xluo->get_gamut work in the set pcsor */
	/* for each xluo. */
	if (li.mode > 0  && li.gmi.usemap) {
		gamut *csgam, *igam, *ogam;
		double gres;			/* Gamut surface feature resolution */
		int    mapres;			/* Mapping rspl resolution */

		if (li.verb)
			printf("Creating Gamut Mapping\n");

		if (li.quality == 0) {			/* Low quality */
  	 		gres = 12.0;
  	 		mapres = 17;
		} else if (li.quality == 1) {	/* Medium */
  	 		gres = 10.0;
  	 		mapres = 25;
		} else if (li.quality == 2) {	/* High */
  	 		gres = 9.0;
  	 		mapres = 33;
		} else { 						/* Ultra High */
  	 		gres = 8.0;
  	 		mapres = 43;
		}

		/* Creat the source colorspace gamut surface */
		if (li.verb)
			printf(" Finding Source Colorspace Gamut\n");

		/* Creat the source image gamut surface in the selected pcsor space */
		if ((csgam = li.in.luo->get_gamut(li.in.luo, gres)) == NULL)
			error ("%d, %s",li.in.x->errc, li.in.x->err);

		/* Grab a given source image gamut. Note that this interface is */
		/* rather loose and dangerous, because the PCS the imported gamut */
		/* is in, is not marked. If gamut was fixed to include the required */
		/* information, this problem could be detected and/or allowed for. */
		if (sgam_name[0] != '\000') {		/* Optional source gamut - ie. from an images */

			if (li.verb)
				printf(" Loading Image Source Gamut '%s'\n",sgam_name);

			igam = new_gamut(gres);

			if (igam->read_gam(igam, sgam_name))
				error("Reading source gamut '%s' failed",sgam_name);
		} else {
			igam = NULL;	/* NULL signals no source image gamut */
		}

		/* Creat the destination gamut surface */
		if (li.verb)
			printf(" Finding Destination Gamut\n");

		if ((ogam = li.out.luo->get_gamut(li.out.luo, gres)) == NULL)
			error ("%d, %s",li.out.x->errc, li.out.x->err);

		if (li.verb)
			printf(" Creating Gamut match\n");

		li.map = new_gammap(li.verb, csgam, igam, ogam, li.gmi.greymf,
		                    li.gmi.glumwcpf, li.gmi.glumwexf, li.gmi.glumbcpf, li.gmi.glumbexf,
		                    li.gmi.glumknf,
		                    li.gmi.gamcpf, li.gmi.gamexf, li.gmi.gamknf,
		                    li.gmi.gampwf, li.gmi.gamswf, li.gmi.satenh,
		                    mapres, NULL, NULL);
		if (li.map == NULL)
			error ("Failed to make gamut map transform");

		ogam->del(ogam);
		if (igam != NULL)
			igam->del(igam);
		csgam->del(csgam);
	}

	/* - - - - - - - - - - - - - - - - - - - */
	/* Verify the given link, assuming all the options are the same */
	if (verify) {
		icmFile *rd_fp;
		icc *rd_icc;
		icmLuBase *luo;

		icColorSpaceSignature ins, outs;	/* Type of input and output spaces */
		int inn, outn;						/* Number of components */
		icmLuAlgType alg;					/* Type of lookup algorithm */

		/* Lookup parameters */
		icmLookupFunc     func   = icmFwd;				/* Default */
		icRenderingIntent intent = icmDefaultIntent;	/* Default */
		icmLookupOrder    order  = icmLuOrdNorm;		/* Default */

		int gc[MAX_CHAN];								/* Grid counter */
		int vres = 8;		//~~9
		double imin[MAX_CHAN], imax[MAX_CHAN];			/* Range of input values */
		double omin[MAX_CHAN], omax[MAX_CHAN];			/* Range of output values */
		double in[MAX_CHAN];							/* Input value */
		double ref[MAX_CHAN];							/* Reference output value */
		double out[MAX_CHAN];							/* Output value */
		double aerr, perr;								/* Average, Peak error */
		double nerr;
		int count, total;
		int pc, lastpc;
		double pin[MAX_CHAN], pref[MAX_CHAN], pout[MAX_CHAN];	/* Peak error values */
		
		if (li.verb)
			printf("Setting up to verify the link\n");

		/* Open up the link file for reading */
		if ((rd_fp = new_icmFileStd_name(link_name,"r")) == NULL)
			error ("Verify: Can't open file '%s'",link_name);

		if ((rd_icc = new_icc()) == NULL)
			error ("Verify: Creation of ICC object failed");

		if ((rv = rd_icc->read(rd_icc,rd_fp,0)) != 0)
			error ("%d, %s",rv,rd_icc->err);

		/* Check that the profile is appropriate */
		if (rd_icc->header->deviceClass != icSigLinkClass)
			error("Profile isn't a device link profile");

		/* Get a conversion object */
		if ((luo = rd_icc->get_luobj(rd_icc, func, intent, icmSigDefaultData, order)) == NULL)
			error ("%d, %s",rd_icc->errc, rd_icc->err);

		/* Get details of conversion (Arguments may be NULL if info not needed) */
		luo->spaces(luo, &ins, &inn, &outs, &outn, &alg, NULL, NULL, NULL);

		if (alg != icmLutType)
			error ("DeviceLink profile doesn't have Lut !");

		/* Get the icm value ranges */
		luo->get_ranges(luo, imin, imax, omin, omax);

		/* Init the grid counter */
		for (i = 0; i < inn; i++)
			gc[i] = 0;

		for (total = 1, i = 0; i < inn; i++, total *= vres)
				; 

		count = 0;
		lastpc = 0;
		nerr = aerr = perr = 0.0;
		for(i = 0; i < inn; ) {
			int j;
			double err;

			/* Create and scale input */
			for (j = 0; j < inn; j++) {
				in[j] = gc[j]/(vres-1.0);
				in[j] = in[j] * (imax[j] - imin[j]) + imin[j];
			}

//			printf("Input %f %f %f %f\n",in[0], in[1], in[2], in[3]);

			/* Create the reference output value */
			devi_devip((void *)&li, ref, in);
			devip_devop((void *)&li, ref, ref);
			devop_devo((void *)&li, ref, ref);

			/* Lookup the icm output value */
			if ((rv = luo->lookup(luo, out, in)) > 1)
				error ("%d, %s",rd_icc->errc,rd_icc->err);
		
//			printf("Output %f %f %f %f\n",out[0], out[1], out[2], out[3]);
//			printf("Ref    %f %f %f %f\n",ref[0], ref[1], ref[2], ref[3]);
//			printf("\n");

			/* Unscale output and compare the results */
			for (err = 0.0, j = 0; j < outn; j++) {
				double o,r;
				o = (out[j] - omin[j])/(omax[j] - omin[j]);
				r = (ref[j] - omin[j])/(omax[j] - omin[j]);
				err += (o - r) * (o - r);
			}
			err = sqrt(err);

			aerr += err;
			nerr++;
			if (err > perr) {
				perr = err;
				for (j = 0; j < outn; j++) {
					pin[j]  = in[j];
					pref[j] = ref[j];
					pout[j] = out[j];
				}
			}

			count++;
			pc = count * 100.0/total + 0.5;
			if (pc != lastpc) {
				printf("\r%2d%%",pc); fflush(stdout);
				lastpc = pc;
			}

			/* Increment the grid counter */
			for (i = 0; i < inn; i++) {
				if (++gc[i] < vres)
					break;			/* No carry */
				gc[i] = 0;			/* Reset digit */
			}
		}

		if (li.verb)
			printf("Finished verfication\n");

		printf("Average error = %f%%, peak error = %f%%\n",aerr * 100.0/nerr, perr * 100.0);
		printf("Input %f %f %f %f\n",pin[0], pin[1], pin[2], pin[3]);
		printf("Output %f %f %f %f\n",pout[0], pout[1], pout[2], pout[3]);
		printf("Ref    %f %f %f %f\n",pref[0], pref[1], pref[2], pref[3]);

		luo->del(luo);
		rd_icc->del(rd_icc);
		rd_fp->del(rd_fp);

	/* - - - - - - - - - - - - - - - - - - - */
	/* Create the link profile */
	} else {
		icmFile *wr_fp;
		icc *wr_icc;

		if (li.verb)
			printf("Creating link profile\n");

		/* Open up the link file for writing */
		if ((wr_fp = new_icmFileStd_name(link_name,"w")) == NULL)
			error ("Write: Can't open file '%s'",link_name);

		if ((wr_icc = new_icc()) == NULL)
			error ("Write: Creation of ICC object failed");

		/* Add all the tags required */

		/* The header: */
		{
			icmHeader *wh = wr_icc->header;

			/* Values that must be set before writing */
			wh->deviceClass     = icSigLinkClass;			/* We are creating a link ! */
			wh->colorSpace      = li.in.h->colorSpace;		/* Input profile device space */
			wh->pcs             = li.out.h->colorSpace;		/* Output profile device space */
			wh->renderingIntent = li.out.h->renderingIntent;	/* Assume output dominates */

			/* Values that should be set before writing */
			wh->manufacturer = 0;
			wh->model        = 0;
			wh->attributes.l = 0;
			wh->flags        = 0;

		}
		/* Profile Description Tag: */
		{
			icmTextDescription *wo;
			char *dst = "Device Link profile - See ProfileSequenceDescTag for more infornmation";
			if ((wo = (icmTextDescription *)wr_icc->add_tag(
			           wr_icc, icSigProfileDescriptionTag,	icSigTextDescriptionType)) == NULL) 
				error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

			wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
			wo->allocate((icmBase *)wo);/* Allocate space */
			strcpy(wo->desc, dst);		/* Copy the string in */
		}
		/* Copyright Tag: */
		{
			icmText *wo;
			char *crt = "Argyll is Open Source Software";
			if ((wo = (icmText *)wr_icc->add_tag(
			           wr_icc, icSigCopyrightTag,	icSigTextType)) == NULL) 
				error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

			wo->size = strlen(crt)+1; 	/* Allocated and used size of text, inc null */
			wo->allocate((icmBase *)wo);/* Allocate space */
			strcpy(wo->data, crt);		/* Copy the text in */
		}
		/* ProfileSequenceDescTag: */
		{
			int i;
			icmProfileSequenceDesc *wo;
			if ((wo = (icmProfileSequenceDesc *)wr_icc->add_tag(
			           wr_icc, icSigProfileSequenceDescTag, icSigProfileSequenceDescType)) == NULL) 
				return 1;

			wo->count = 2; 		/* Number of descriptions in sequence */
			wo->allocate((icmBase *)wo);	/* Allocate space for all the DescStructures */

			/* Fill in each description structure in sequence */

			/* Input profile */
			for (i = 0; i < wo->count; i++) {
				icc *iccs;
				icmHeader *sh;
				icmSignature *tsig;
				icmTextDescription *ddesc;
				icmTextDescription *mdesc;

				if (i == 0) {
					iccs = li.in.c;		/* Input profile */
					sh = li.in.h;		/* Input profile header */
				} else if (i == (wo->count-1)) {
					iccs = li.out.c;	/* Output profile */
					sh = li.out.h;		/* Output profile header */
				} else {
					error("Abstract profiles in link not implemented yet!");
				}

				/* Try and read the technology tag */
				if ((tsig = (icmSignature *)iccs->read_tag(iccs, icSigTechnologyTag)) != NULL) {
					if (tsig->ttype != icSigSignatureType)	/* oops */
						tsig = NULL;
				}

				/* Try and read the Device Manufacturers Description Tag */
				if ((ddesc = (icmTextDescription *)iccs->read_tag(
				                                      iccs, icSigDeviceMfgDescTag)) != NULL) {
					if (ddesc->ttype != icSigTextDescriptionType)	/* oops */
						ddesc = NULL;
				}
			
				/* Try and read the Model Manufacturers Description Tag */
				if ((mdesc = (icmTextDescription *)iccs->read_tag(
				                                      iccs, icSigDeviceModelDescTag)) != NULL) {
					if (mdesc->ttype != icSigTextDescriptionType)	/* oops */
						mdesc = NULL;
				}
			
				/* Header information */
				wo->data[i].deviceMfg   = sh->manufacturer;
				wo->data[i].deviceModel = sh->model;
				wo->data[i].attributes  = sh->attributes;

				/* Technology signature */
				if (tsig != NULL)
					wo->data[i].technology = tsig->sig; 

				if (ddesc != NULL) {
					wo->data[i].device.size = ddesc->size;
					wo->data[i].allocate(&wo->data[i]); 			/* Allocate space */
					strcpy(wo->data[i].device.desc, ddesc->desc);

					wo->data[i].device.ucLangCode = ddesc->ucLangCode;
					wo->data[i].device.ucSize = ddesc->ucSize;
					wo->data[i].allocate(&wo->data[i]);					/* Allocate space */
					memcpy(wo->data[i].device.ucDesc, ddesc->ucDesc, 2 * ddesc->ucSize);

					wo->data[i].device.scCode = ddesc->scCode;	
					wo->data[i].device.scSize = ddesc->scSize;		
					strcpy(wo->data[i].device.scDesc, ddesc->scDesc);
				}

				/* model Text description */
				if (mdesc != NULL) {
					wo->data[i].model.size = mdesc->size;
					wo->data[i].allocate(&wo->data[i]); 			/* Allocate space */
					strcpy(wo->data[i].model.desc, mdesc->desc);

					wo->data[i].model.ucLangCode = mdesc->ucLangCode;
					wo->data[i].model.ucSize = mdesc->ucSize;
					wo->data[i].allocate(&wo->data[i]);					/* Allocate space */
					memcpy(wo->data[i].model.ucDesc, mdesc->ucDesc, 2 * mdesc->ucSize);

					wo->data[i].model.scCode = mdesc->scCode;	
					wo->data[i].model.scSize = mdesc->scSize;		
					strcpy(wo->data[i].model.scDesc, mdesc->scDesc);
				}
			}
		}
		/* 16 bit input device -> output device lut: */
		{
			icmLut *wo;

			/* Link Lut = AToB0 */
			if ((wo = (icmLut *)wr_icc->add_tag(
			           wr_icc, icSigAToB0Tag,	icSigLut16Type)) == NULL) 
				error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

			wo->inputChan = li.in.chan;
			wo->outputChan = li.out.chan;

			if (li.quality >= 3)
	    		wo->inputEnt = 4096;
			else if (li.quality == 2)
	    		wo->inputEnt = 2048;
			else
	    		wo->inputEnt = 256;
			
			/* See discussion in imdi/imdi_gen.c for ideal numbers */
			switch (li.in.chan) {
				case 0:
					error ("Illegal number of input chanels");
				case 1:
					if (li.quality >= 3)
			  		  	wo->clutPoints = 4370;
					else if (li.quality == 2)
			  		  	wo->clutPoints = 772;
					else
			  		  	wo->clutPoints = 256;
					break;

				case 2:
					if (li.quality >= 2)
		  		  		wo->clutPoints = 256;
					else
		  		  		wo->clutPoints = 86;
					break;
				case 3:
					if (li.quality >= 3)
		  		  		wo->clutPoints = 52;
					else if (li.quality == 2)
		  		  		wo->clutPoints = 33;
					else if (li.quality == 1)
		  		  		wo->clutPoints = 17;
					else
		  		  		wo->clutPoints = 9;
					break;
				case 4:
					if (li.quality >= 3)
		  		  		wo->clutPoints = 33;
					else if (li.quality == 2)
		  		  		wo->clutPoints = 18;
					else if (li.quality == 1)
		  		  		wo->clutPoints = 9;
					else
		  		  		wo->clutPoints = 6;
					break;
				case 5:
					if (li.quality >= 3)
		  		  		wo->clutPoints = 18;
					else if (li.quality == 2)
		  		  		wo->clutPoints = 16;
					else 
						wo->clutPoints = 9;
					break;
				case 6:
					if (li.quality >= 3)
		  		  		wo->clutPoints = 12;
					else if (li.quality == 2)
		  		  		wo->clutPoints = 9;
					else 
						wo->clutPoints = 6;
					break;
				case 7:
					if (li.quality >= 3)
		  		  		wo->clutPoints = 8;
					else if (li.quality == 2)
		  		  		wo->clutPoints = 7;
					else 
						wo->clutPoints = 6;
					break;
				case 8:
					if (li.quality >= 3)
		  		  		wo->clutPoints = 7;
					else if (li.quality == 2)
		  		  		wo->clutPoints = 6;
					else 
						wo->clutPoints = 5;
					break;
				deault: /* > 8 chan */
					wo->clutPoints = 3;
					break;
			}	
			if (li.quality >= 3)
	    		wo->outputEnt = 4096;
			else if (li.quality == 2)
	    		wo->outputEnt = 2048;
			else
	    		wo->outputEnt = 256;
			if (li.clutres > 0)		/* clut resolution override */
  		  		wo->clutPoints = li.clutres;
			li.clutres = wo->clutPoints;	/* Actual resolution */
			wo->allocate((icmBase *)wo);/* Allocate space */

			/* Special case if input profile is Lut with matrix */
			if (li.in.alg == icmLutType) {
				icxLuLut *lu = (icxLuLut *)li.in.luo;
				lu->get_matrix(lu, wo->e);		/* Copy it across */
			}

			if (li.verb)
				printf("Filling in Lut table\n");
#ifdef DEBUG
#define DBGNO 2		/* Up to 10 */

#ifdef NEVER
			/* Test a single given cmyk -> cmyk value */
			{
				double in[10][MAX_CHAN];
				double out[MAX_CHAN];
				in[0][0] = 0.0;
				in[0][1] = 0.0;
				in[0][2] = 0.0;
				in[0][3] = 0.0;

				in[0][0] = 1.0;
				in[0][1] = 1.0;
				in[0][2] = 1.0;
				in[0][3] = 0.0;

				for (i = 0; i < DBGNO; i++) {
					printf("Input %f %f %f %f\n",in[i][0], in[i][1], in[i][2], in[i][3]);
					devi_devip((void *)&li, out, in[i]);
					devip_devop((void *)&li, out, out);
					devop_devo((void *)&li, out, out);
					printf("Output %f %f %f %f\n\n",out[0], out[1], out[2], out[3]);
				}
			}
#else
			/* Test a single given rgb -> rgb value */
			{
				double in[10][MAX_CHAN];
				double out[MAX_CHAN];
				in[0][0] = 0.0;
				in[0][1] = 0.0;
				in[0][2] = 0.0;

				in[0][0] = 1.0;
				in[0][1] = 1.0;
				in[0][2] = 1.0;

				for (i = 0; i < DBGNO; i++) {
					printf("Input %f %f %f\n",in[i][0], in[i][1], in[i][2]);
					devi_devip((void *)&li, out, in[i]);
					devip_devop((void *)&li, out, out);
					devop_devo((void *)&li, out, out);
					printf("Output %f %f %f\n\n",out[0], out[1], out[2]);
				}
			}

#endif
#else
			/* Use helper function to do the hard work. */
			if (li.verb) {
				for (li.total = 1, i = 0; i < wo->inputChan; i++, li.total *= wo->clutPoints)
					; 
				printf(" 0%%"); fflush(stdout);
			}
			if (wo->set_tables(wo,
				&li,						/* Context */
				li.in.h->colorSpace,		/* Input color space */
				li.out.h->colorSpace,		/* Output color space */
				devi_devip,					/* Input transfer tables devi->devi' */
				NULL, NULL,					/* Use default input colorspace range */
				devip_devop,				/* devi' -> devo' transfer function */
				NULL, NULL,					/* Default output colorspace range */
				devop_devo					/* Output transfer tables, devo'->devo */
			) != 0)
				error("Setting 16 bit Lut failed: %d, %s",wr_icc->errc,wr_icc->err);
			if (li.verb) {
				printf("\n");
			}
#endif
		}

		if (li.verb && li.wphack && li.wphacked == 0)
			printf("Warning :- white point hack didn't trigger!\n");
		if (li.verb && li.wphack && li.wphacked > 1)
			printf("Warning :- white point hack trigger more than once! (%d)\n",li.wphacked);

		if (li.verb)
			printf("Writing out file\n");

		/* Write the file out */
		if ((rv = wr_icc->write(wr_icc,wr_fp,0)) != 0)
			error ("Write file: %d, %s",rv,wr_icc->err);
		
		wr_icc->del(wr_icc);
		wr_fp->del(wr_fp);
	}

	/* - - - - - - - - - - - - - - - - - - - */
	/* Cleanup source profiles and exit */

	if (li.pcs2k != NULL)			/* Free up PCS->K lookup for neutral hack */
		li.pcs2k->del(li.pcs2k);

	if (li.map != NULL)
		li.map->del(li.map);

	if (li.abs_luo != NULL) {		/* Free up abstract transform */
		li.abs_luo->del(li.abs_luo);
		li.abs_xicc->del(li.abs_xicc);
		li.abs_icc->del(li.abs_icc);
		li.abs_fp->del(li.abs_fp);
	}

	li.in.luo->del(li.in.luo);
	li.in.x->del(li.in.x);
	li.in.c->del(li.in.c);
	li.in.fp->del(li.in.fp);
	if (li.in.devlin != NULL)
		li.in.devlin->del(li.in.devlin);

	li.out.luo->del(li.out.luo);
	li.out.x->del(li.out.x);
	li.out.c->del(li.out.c);
	li.out.fp->del(li.out.fp);
	if (li.out.devlin != NULL)
		li.out.devlin->del(li.out.devlin);

	return 0;
}







