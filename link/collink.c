
/* 
 * collink
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
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* TTBD:
 *
 * Device curve resolution should be taken from the device profiles,
 * rather than depending on the quality settings (scRGB compression curves)
 *
 * Abstract link support intent doesn't work properly for anything
 * other than absolute. This should really be fixed.
 */

/* NOTES:

	Normally the device side per channel curves are copied from
	the source profiles to the link profile on the assumption that
	the raw linearisation they do is good, and should be maintained
	for best overall transform accuracy. Since the intermediate
	link is done in an Lab like space, then this linearisation
	will even out quantisation error introduced by the Lut
	in perceptual space. In the case of a Matrix profile with a native
	XYZ PCS, these assumptions break, since the device curves
	would generally be linearising the device to Y, not a perceptual
	space. For this reason we add a Y to L* type curve to the input
	linearisation curves. For XYZ CLUT based profiles, we just have
	to hope that the profile has been created well, and that the
	input curves distribute the indexes reasonably perceptually
	evenly. For an output (XYZ) Matrix profile, we also
	use an L* style mixing space.  

	In general the per channel curves should be:

	  No curve, or

	  An optimized curve, or

      The source profile per channel curve plus
	  a Y to L* curve it's a Matrix profile.


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
	to the image data before gamut surface extraction. This is yet
	another avenue for user & implementor confusion, with regard
	to intents, even without allowing the user to set the intent
	of the abstract profile.

 */

#undef USE_MERGE_CLUT_OPT	/* [Undef] When using inverse A2B table, merge the output luts */
							/* with the clut for faster operation, and clipping in */
							/* Jab space. Turned off because it affects the accuracy too much, */
							/* and xicc handles Jab clip without this now. */

#define USE_CAM_CLIP_OPT	/* [Define] Clip out of gamut in CAM space rather than XYZ or L*a*b* */
#define ENKHACK				/* [Define] Enable K hack code */
#undef WARN_CLUT_CLIPPING	/* [Undef] Print warning if setting clut clips */

#undef DEBUG		/* Report values of each sample transformed */
#undef DEBUGC 		/* ie "if (tt)" */		/* Debug condition */
#undef DEBUG_ONE	/* test a single value out. Look for DBGNO to set value. */
#undef NEUTKDEBUG	/* print info about neutral L -> K mapping */


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
#include "gammap.h"
// ~~~99
#include "vrml.h"

void usage(char *diag, ...) {
	int i;
	fprintf(stderr,"Link ICC profiles, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"  Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: collink [options] srcprofile dstprofile linkedprofile\n");
	fprintf(stderr," -v              Verbose\n");
	fprintf(stderr," -A manufacturer Manufacturer description string\n");
	fprintf(stderr," -M model        Model description string\n");
	fprintf(stderr," -D description  Profile Description string (Default \"inoutfile\")\n");
	fprintf(stderr," -C copyright    Copyright string\n");
	fprintf(stderr," -V              Verify existing profile, rather than link\n");
	fprintf(stderr," -q lmhu         Quality - Low, Medium (def), High, Ultra\n");
	fprintf(stderr," -r res          Override clut res. set by -q\n");
	fprintf(stderr," -n              Don't preserve device linearization curves in result\n");
	fprintf(stderr," -f              Special :- Force neutral colors to be K only output\n");
	fprintf(stderr," -fk             Special :- Force K only neutral colors to be K only output\n");
	fprintf(stderr," -F              Special :- Force all colors to be K only output\n");
	fprintf(stderr," -fcmy           Special :- Force 100%% C,M or Y only to stay pure \n");
	fprintf(stderr," -p absprof      Include abstract profile in link\n");
	fprintf(stderr," -s              Simple Mode (default)\n");
	fprintf(stderr," -g [src.gam]    Gamut Mapping Mode [optional source image gamut]\n");
	fprintf(stderr," -G [src.gam]    Gamut Mapping Mode using inverse outprofile A2B\n");
	fprintf(stderr," Simple Mode Options:\n");
	fprintf(stderr," -i in_intent    p = perceptual, r = relative colorimetric,\n");
	fprintf(stderr,"                 s = saturation, a = absolute colorimetric\n");
	fprintf(stderr," -o out_intent   p = perceptual, r = relative colorimetric,\n");
	fprintf(stderr,"                 s = saturation, a = absolute colorimetric\n");
	fprintf(stderr," Gamut Mapping   Mode Options:\n");
	fprintf(stderr," -i intent       set linking intent from the following choice:\n");
	for (i = 0; ; i++) {
		icxGMappingIntent gmi;
		if (xicc_enum_gmapintent(&gmi, i, NULL) == icxIllegalGMIntent)
			break;

		fprintf(stderr,"              %s\n",gmi.desc);
	}
	fprintf(stderr," -w [J,a,b]      Use forced whitepoint hack [optional target point]\n");
//	fprintf(stderr," -W J,a,b        Forced whitepoint adjustment by delta Jab\n");
	fprintf(stderr," -c viewcond     set source viewing conditions for %s,\n",icxcam_description(cam_default));
	fprintf(stderr,"                 either an enumerated choice, or a parameter\n");
	fprintf(stderr," -d viewcond     set destination viewing conditions for %s,\n",icxcam_description(cam_default));
	fprintf(stderr,"                 either an enumerated choice, or parameter:value changes\n");
	for (i = 0; ; i++) {
		icxViewCond vc;
		if (xicc_enum_viewcond(NULL, &vc, i, NULL, 1, NULL) == -999)
			break;

		fprintf(stderr,"            %s\n",vc.desc);
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
	fprintf(stderr," -t tlimit       set source total ink limit, 0 - 400%% (estimate by default)\n");
	fprintf(stderr," -T klimit       set source black ink limit, 0 - 100%% (estimate by default)\n");
	fprintf(stderr," Inverse outprofile A2B Options:\n");
	fprintf(stderr," -k tezhxr       CMYK Black generation\n");
	fprintf(stderr,"                 t = transfer K from source to destination, e = retain K of destination B2A table\n");
	fprintf(stderr,"                 z = zero K, h = 0.5 K, x = maximum K, r = ramp K (default)\n");
	fprintf(stderr," -k p stle stpo enpo enle shape\n");
	fprintf(stderr,"                 p = black target generation curve parameters\n");
	fprintf(stderr," -k q stle0 stpo0 enpo0 enle0 shape0 stle2 stpo2 enpo2 enle2 shape2\n");
	fprintf(stderr,"                 q = transfer source K to dual curve limits\n");
	fprintf(stderr," -K parameters   Same as -k, but target is K locus rather than K value itself\n");
	fprintf(stderr," -l tlimit       set destination total ink limit, 0 - 400%% (estimate by default)\n");
	fprintf(stderr," -L klimit       set destination black ink limit, 0 - 100%% (estimate by default)\n");
	fprintf(stderr," -P              Create gamut gammap_p.wrl and gammap_s.wrl diagostics\n");
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
								/* 2 = 0.5, 3 = max, 4 = ramp, 5 = curve, 6 = dual curve */
								/* 7 = outpupt profile K value */
	int locus;					/* 0 = K target value, 1 = K locus value */
	icxInk ink;					/* Ink parameters */

	/* Operational parameters */
	icc *c;
	icmHeader *h;
	xicc *x;
	icxLuBase *luo;				/* Base XLookup type object */
	icmLuAlgType alg;			/* Type of lookup algorithm */
	icColorSpaceSignature csp;	/* Colorspace */
	int chan;					/* Channels */
	int nocurve;				/* NZ to not use device curve in per channel curve */
	int lcurve;					/* 1 to apply a Y like to L* curve for XYZ Matrix profiles */
								/* 2 to apply a Y to L* curve for XYZ space */
	double wp[3];				/* Lab/Jab white point for profile used by wphack & xyzscale */
	icxLuBase *b2aluo;			/* B2A lookup for inking == 7 */
}; typedef struct _profinfo profinfo;

/* Structure that holds all the color lookup information */
struct _link {
	/* Overall options */
	int verb;
	int gamdiag;	/* nz, create gammap.wrl diagnostic */
	int total, count, last;	/* Progress count information */
	int mode;		/* 0 = simple mode, 1 = mapping mode, 2 = mapping mode with inverse A2B */
	int quality;	/* 0 = low, 1 = medium, 2 = high, 3 = ultra */
	int clutres;	/* 0 = quality default, !0 = override, then actual during link */
	int src_kbp;	/* nz = Use K only black point as src gamut black point */
	int dst_kbp;	/* nz = Use K only black point as dst gamut black point */
	int dst_cmymap;	/* masks C = 1, M = 2, Y = 4 to force 100% cusp map */

	icColorSpaceSignature pcsor;	/* PCS to use between in & out profiles */

	int nhack;		/* 0 = off, 1 = hack to map input neutrals to output K only, */
					/* 2 = map 000K to output K only, 3 = all to K */
	int cmyhack;	/* CMY 100% colorant map though hack, 1 = C, 2 = M, 4 = Y */
	rspl *pcs2k;	/* PCS L to K' lookup for nhack */
	int wphack;		/* 0 = off, 1 = hack to map input wp to output wp, 2 = to hwp[] */
	double hwp[3];	/* hack destination white point in PCS space */
	int wphacked;	/* Operation flag, set nz if white point was translated */
	int rel_oride;	/* Relative override flag */

	icmFile *abs_fp;	/* Abstract profile transform */
	icRenderingIntent abs_intent;
	icc *abs_icc;
	xicc *abs_xicc;
	icxLuBase *abs_luo;	/* NULL if none */

					/* (We current assume that xyzscale can't be used with gmi) */
	double xyzscale;	/* < 1.0 if Y is to be scaled in destination XYZ space */
	double swxyz[3];	/* Source white point in XYZ */

	icxGMappingIntent gmi;
	gammap *map;	/* Gamut mapping */
	gammap *Kmap;	/* Gamut mapping K in to K out nhack == 2 and K in to K out */
	

	/* Per profile setup information */
	profinfo in;
	profinfo out;

}; typedef struct _link link;


/* ------------------------------------------- */
//#define YSCALE 1.0
#define YSCALE (2.0/1.3)

/* Extra non-linearity applied to BtoA XYZ PCS */
/* This distributes the LUT indexes more evenly in */
/* perceptual space, greatly improving the B2A accuracy of XYZ LUT */
/* Since typically XYZ doesn't use the full range of 0-2.0 allowed */
/* for in the encoding, we scale the cLUT index values to use the 0-1.3 range */

/* (For these functions the encoded XYZ 0.0 - 2.0 range is 0.0 - 1.0 ??) */

/* Y to L* */
static void y2l_curve(double *out, double *in, int isXYZ) {
	int i;
	double val;
	double isc = 1.0, osc = 1.0;

	/* Scale from 0.0 .. 1.999969 to 0.0 .. 1.0 and back */
	/* + range adjustment */
	if (isXYZ) {
		isc = 32768.0/65535.0 * YSCALE;
		osc = 65535.0/32768.0;
	}

	for (i = 0; i < 3; i++) {
		val = in[i] * isc;
		if (val > 0.008856451586)
			val = 1.16 * pow(val,1.0/3.0) - 0.16;
		else
			val = 9.032962896 * val;
		if (val > 1.0)
			val = 1.0;
		out[i] = val * osc;
	}
}

/* L* to Y */
static void l2y_curve(double *out, double *in, int isXYZ) {
	int i;
	double val;
	double isc = 1.0, osc = 1.0;

	/* Scale from 0.0 .. 1.999969 to 0.0 .. 1.0 and back */
	/* + range adjustment */
	if (isXYZ) {
		isc = 32768.0/65535.0;
		osc = 65535.0/32768.0 / YSCALE;
	}

	/* Use an L* like curve, scaled to the maximum XYZ value */
	for (i = 0; i < 3; i++) {
		val = in[i] * isc;
		if (val > 0.08)
			val = pow((val + 0.16)/1.16, 3.0);
		else
			val = val/9.032962896;
		out[i] = val * osc;
	}
}

/* ------------------------------------------- */
/* Functions called back in setting up the transform table */

#ifdef DEBUGC

static int tt = 0;

#endif	/* DEBUGC */

/* Input table, DevIn -> DevIn' */
void devi_devip(void *cntx, double *out, double *in) {
	int rv = 0;
	link *p = (link *)cntx;

#ifdef DEBUGC
	if (in[0] == 1.0 && in[1] == 1.0 && in[2] == 1.0 && in[3])
		tt = 1;
#endif

#ifdef DEBUG
#ifdef DEBUGC
	DEBUGC
#endif
	printf("DevIn->DevIn' got %f %f %f %f\n",in[0], in[1], in[2], in[3]); fflush(stdout);
#endif

	if (p->in.nocurve) {	/* Don't use profile per channel curves */
		int i;
		for (i = 0; i < p->in.chan; i++)
			out[i] = in[i];
	} else {				/* Use input profile per channel curves */
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

	if (p->in.lcurve) { 		/* Apply Y to L* */
//printf("~1 y2l_curve got %f %f %f, isXYZ %d\n",in[0],in[1],in[2],p->in.lcurve == 2);
		y2l_curve(out, out, p->in.lcurve == 2);
	}

#ifdef DEBUG
#ifdef DEBUGC
	DEBUGC
#endif
	printf("DevIn->DevIn' ret %f %f %f %f\n",out[0], out[1], out[2], in[3]); fflush(stdout);
#endif
}


/* - - - - - - - - - - - - */
/* clut, DevIn' -> DevOut' */
void devip_devop(void *cntx, double *out, double *in) {
	double win[MAX_CHAN];			/* working input values */
	double pcsv[MAX_CHAN];			/* PCS intermediate value, pre-gamut map */
	double pcsvm[MAX_CHAN];			/* PCS intermediate value, post-gamut map */
	double locus[MAX_CHAN];			/* Auxiliary locus values */
	int wptrig = 0;				/* White point hack triggered */
	double konlyness = 0.0;		/* Degree of K onlyness */
	int ntrig = 0;				/* K only output hack triggered */
	int cmytrig = 0;			/* CMY output hack triggered */
	int i, rv = 0;
	link *p = (link *)cntx;

#ifdef DEBUG
#ifdef DEBUGC
	DEBUGC
#endif
	printf("DevIn'->DevOut' got %f %f %f %f\n",in[0], in[1], in[2], in[3]); fflush(stdout);
#endif

#ifdef ENKHACK
	/* Handle neutral recognition/output K only hack */
	/* (see discussion at top of file for generalization of this idea) */
	if (p->nhack == 1 || p->nhack == 2) {
		double thr = (0.5)/(p->clutres-1.0); /* Match threshold */

		if (p->nhack == 1) {
			/* We want to see if the input colors are equal (Or a=b= 0.0 ??) */
			/* li.nhack should have set p->in.nocurve, so we should be getting raw */
			/* input space device values here. It also made sure that there are at */
			/* least 3 input channels. */

		    if (fabs(in[0] - in[1]) < thr
		     && fabs(in[1] - in[2]) < thr
		     && fabs(in[2] - in[0]) < thr)
				ntrig = 1;			/* K only output triggered flag */

		} else if (p->nhack == 2) {
			double maxcmy;		/* Comute a degree of source "K onlyness" */
			double maxcmyk;

			maxcmy = in[0];			/* Compute minimum of CMY */
			if (in[1] > maxcmy)
				maxcmy = in[1];
			if (in[2] > maxcmy)
				maxcmy = in[2];

			maxcmyk = maxcmy;		/* Compute minimum of all inks */
			if (in[3] > maxcmyk)
				maxcmyk = in[3];

//printf("~1 maxcmy = %f, maxcmyk = %f, in[3] = %f\n",maxcmy,maxcmyk,in[3]);
			if (in[3] <= 0.0 || maxcmy > in[3]) {
				konlyness = 0.0;
			} else {
				konlyness = (in[3] - maxcmy)/in[3];	
			}

			/* As we approach no colorant, blend towards no Konlyness */
			if (maxcmyk < 0.2)
				konlyness *= maxcmyk/0.2;

			/* We want to see if the input colors are exactly K only. */
		    if (in[0] < thr
		     && in[1] < thr
		     && in[2] < thr)
				ntrig = 1;			/* K only output triggered flag */
#ifdef DEBUG
#ifdef DEBUGC
			DEBUGC
#endif
			printf("konlyness set to %f\n",konlyness);
#endif

		}
	}
	/* Handle 100% CMY hack */
	if (p->cmyhack != 0) {
		double thr = (0.5)/(p->clutres-1.0); /* Match threshold */

		if (p->cmyhack & 1) {
		    if (in[0] > (1.0 - thr)
		     && in[1] < thr
		     && in[2] < thr
		     && (p->in.chan < 4 || in[3] < thr))
				cmytrig |= 1;
		}
		if (p->cmyhack & 2) {
		    if (in[0] < thr
		     && in[1] > (1.0 - thr)
		     && in[2] < thr
		     && (p->in.chan < 4 || in[3] < thr))
				cmytrig |= 2;
		}
		if (p->cmyhack & 4) {
		    if (in[0] < thr
		     && in[1] < thr
		     && in[2] > (1.0 - thr)
		     && (p->in.chan < 4 || in[3] < thr))
				cmytrig |= 4;
		}
	}
#endif /* ENKHACK */

	if (p->in.lcurve) {	/* Apply L* to Y */
		l2y_curve(win, in, p->in.lcurve == 2);
#ifdef DEBUG
#ifdef DEBUGC
	DEBUGC
#endif
	printf("win[] set to L* value %f %f %f %f\n",win[0], win[1], win[2], win[3]); fflush(stdout);
#endif

	} else {
		for (i = 0; i < p->in.chan; i++)
			win[i] = in[i];
#ifdef DEBUG
#ifdef DEBUGC
	DEBUGC
#endif
	printf("win[] set to in[] value %f %f %f %f\n",win[0], win[1], win[2], win[3]); fflush(stdout);
#endif
	}

	/* Do DevIn' -> PCS */
	switch(p->in.alg) {
	    case icmMonoFwdType: {
			icxLuMono *lu = (icxLuMono *)p->in.luo;	/* Safe to coerce */

			if (p->in.nocurve) {	/* No explicit curve, so do implicit here */
				rv |= lu->fwd_curve(lu, pcsv, win);
				rv |= lu->fwd_map(lu, pcsv, pcsv);
			} else {
				rv |= lu->fwd_map(lu, pcsv, win);
			}
			rv |= lu->fwd_abs(lu, pcsv, pcsv);
			break;
		}
	    case icmMatrixFwdType: {
			icxLuMatrix *lu = (icxLuMatrix *)p->in.luo;	/* Safe to coerce */

			if (p->in.nocurve) {	/* No explicit curve, so do implicit here */
				rv |= lu->fwd_curve(lu, pcsv, win);
				rv |= lu->fwd_matrix(lu, pcsv, pcsv);
			} else {
				rv |= lu->fwd_matrix(lu, pcsv, win);
			}
			rv |= lu->fwd_abs(lu, pcsv, pcsv);
			break;
		}
	    case icmLutType: {
			icxLuLut *lu = (icxLuLut *)p->in.luo;	/* Safe to coerce */
			if (p->in.nocurve) {	/* No explicit curve, so we've got Dev */
				/* Since not PCS, in_abs and matrix cannot be valid, */
				/* so input curve on own is ok to use. */
				rv |= lu->input(lu, pcsv, win);		/* Dev  -> Dev' */
				rv |= lu->clut(lu, pcsv, pcsv);		/* Dev' -> PCS' */
			} else {	/* We've got Dev' */
				rv |= lu->clut(lu, pcsv, win);		/* Dev' -> PCS' */
			}
			/* We've got the input profile PCS' at this point. */

			/* If we're transfering the K value from the input profile to the */
			/* output, copy it into locus[], which will be given to the inverse */
			/* lookup function, else the inverse lookup will generate a K using */
			/* the curve parameters. */  
//printf("~1 out.inking = %d\n",p->out.inking);
			if (p->out.inking == 0 || p->out.inking == 6) {
				if (p->out.locus) {
					/* Converts PCS' to K locus proportion */
					lu->clut_locus(lu, locus, pcsv, win);	/* Compute possible locus values */
//printf("~1 looked up locus value\n");
				} else {
					for (i = 0; i < p->in.chan; i++)		/* Target is K input value */
						locus[i] = win[i];
					/* Convert K' to K value ready for aux target */
					if (!p->in.nocurve) {	/* we have an input curve, so convert Dev' -> Dev */
						lu->inv_input(lu, locus, locus);
					}
//printf("~1 copied win to locus\n");
				}
#ifdef DEBUG
#ifdef DEBUGC
				DEBUGC
#endif
				printf("Got possible K %s of %f %f %f %f\n",p->out.locus ? "locus" : "value", locus[0],locus[1],locus[2],locus[3]);
#endif
			}
			rv |= lu->output(lu, pcsv, pcsv);	/* PCS' ->     */
			rv |= lu->out_abs(lu, pcsv, pcsv);	/*         PCS */
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
     *
     *	and locus[] contains any auxiliar target values if the
	 *  auxiliary is not being created by a rule applied to the PCS.
	 */

	/*
	 * The order to do this intermediate processing is hard to figure out,
	 * as is the interaction between such elements. How should the
	 * abstract profile be properly handled ?
	 * what should we do if the wphack is on and Y scaling is on ?
	 */

#ifdef DEBUG
#ifdef DEBUGC
	DEBUGC
#endif
	printf("PCS before map %f %f %f\n",pcsv[0], pcsv[1], pcsv[2]); fflush(stdout);
#endif

	if (p->wphack) {
		int e;
		double dd = 0.0;
		for (e = 0; e < 3; e++) {			/* Does this match the input white point ? */
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

	/* Do luminence scaling if requested */
	if (wptrig == 0 && p->xyzscale < 1.0) {
		double xyz[3];

//printf("~1 got xyzscale = %f\n",p->xyzscale);
//printf("PCS %f %f %f\n",pcsv[0], pcsv[1], pcsv[2]);

		/* Convert our PCS to XYZ */
		if (p->pcsor == icxSigJabData) {
			/* We're being bad in delving inside the xluo, but we'll fix it latter */
			p->out.luo->cam->cam_to_XYZ(p->out.luo->cam, xyz, pcsv);
		} else
			error("Internal :- not setup to handle Y scaling and non-Jab PCS");

//printf("XYZ %f %f %f\n",xyz[0], xyz[1], xyz[2]);
		/* Scale it */
		xyz[0] *= p->xyzscale;
		xyz[1] *= p->xyzscale;
		xyz[2] *= p->xyzscale;

//printf("scaled XYZ %f %f %f\n",xyz[0], xyz[1], xyz[2]);
		/* Convert back to PCS */
		if (p->pcsor == icxSigJabData) {
			/* We're being bad in delving inside the xluo, but we'll fix it latter */
			p->out.luo->cam->XYZ_to_cam(p->out.luo->cam, pcsv, xyz);
		} else
			error("Internal :- not setup to handle Y scaling and non-Jab PCS");

//printf("scaled PCS %f %f %f\n",pcsv[0], pcsv[1], pcsv[2]);
#ifdef DEBUG
#ifdef DEBUGC
	DEBUGC
#endif
		printf("PCS after Y scale %f %f %f\n",pcsv[0], pcsv[1], pcsv[2]); fflush(stdout);
#endif
	}


	/* Do gamut mapping */
	if (wptrig == 0 && p->mode > 0 && p->gmi.usemap) {
		/* We've used pcsor to ensure PCS space is appropriate */
		
		/* Doing XXXK -> XXXK */
		if (p->nhack == 2) {
			/* Ideally we would create a 4D PCSK -> PCSK gamut mapping */
			/* to smoothly and accurately cope with the changing source */
			/* and destination gamuts acording to their degree of "K onlyness". */
			/* In practice we're going to simply interpolated between */
			/* two extremes: unrestricted gamut and K only black gamut. */
			double map0[3], map1[3];

			/* Compute blend of normal gamut map and Konly to Konly gamut map */	
			{
				p->map->domap(p->map, map0, pcsv);
				p->Kmap->domap(p->Kmap, map1, pcsv);
				icmBlend3(pcsvm, map0, map1, konlyness);
			}

#ifdef DEBUG
#ifdef DEBUGC
	DEBUGC
#endif
			printf("PCS after map0 %f %f %f map1 %f %f %f\n", map0[0], map0[1], map0[2], map1[0], map1[1], map1[2]);
#endif

		/* Normal gamut mapping */
		} else {
			{
				p->map->domap(p->map, pcsvm, pcsv);
			}
		}


#ifdef DEBUG
#ifdef DEBUGC
		DEBUGC
#endif
		printf("PCS after map %f %f %f\n",pcsvm[0], pcsvm[1], pcsvm[2]); fflush(stdout);
#endif
	} else {
		pcsvm[0] = pcsv[0];
		pcsvm[1] = pcsv[1];
		pcsvm[2] = pcsv[2];
	}

	/* Gamut mapped PCS value is now in pcsvm[] */

	/* Abstract profile transform, PCS -> PCS */
	/* pcsor -> abstract -> pcsor conversion */
	/* We're applying any abstract profile after gamut mapping, */
	/* on the assumption is is primarily being used to "correct" the */
	/* output device. Ideally the gamut mapping should take the change */
	/* the abstract profile has on the output device into account, but */
	/* currently we're not doing this.. */
	if (wptrig == 0 && p->abs_luo != NULL) {
		/* Abstract profile is either absolute or relative.  */
		/* We need to convert the current PCS into something compatible. */
		/* This is more ugly than it really should be, so we're ignoring it. */
		/* We should really run the source through the abstract profile before */
		/* creating the gamut mapping, to be able to use abstract with gamut */
		/* mapping properly. */
		p->abs_luo->lookup(p->abs_luo, pcsvm, pcsvm);
#ifdef DEBUG
#ifdef DEBUGC
	DEBUGC
#endif
		printf("PCS after abstract %f %f %f\n",pcsvm[0], pcsvm[1], pcsvm[2]); fflush(stdout);
#endif
	}

	/* If we're using the existing B2A inking to determine K, */
	/* lookup the output profiles K value for this PCS */
	if (p->mode >= 2 && p->out.inking == 7) {
		double tdevv[MAX_CHAN];	

//printf("~1 dealing with out.inking = %d\n",p->out.inking);
		if (p->out.alg != icmLutType || p->out.c->header->colorSpace != icSigCmykData)
			error ("Attempting to use non-CMYK output profile to determine K inking");

		/* Lookup PCS in B2A of output profile to get target K value */  
//printf("~1 looking up pcs %f %f %f in B2A\n", pcsvm[0], pcsvm[1], pcsvm[2]);
		p->out.b2aluo->lookup(p->out.b2aluo, tdevv, pcsvm);
//printf("~1 resulting dev %f %f %f %f\n", tdevv[0], tdevv[1], tdevv[2], tdevv[3]);

		if (p->out.locus) {
			double tpcsv[MAX_CHAN];	
			icxLuLut *lu = (icxLuLut *)p->out.luo;		/* Safe to coerce */

			/* Convert PCS to PCS' ready for locus lookup */
			lu->in_abs(lu, tpcsv, pcsvm);
			lu->matrix(lu, tpcsv, tpcsv);
			lu->input(lu, tpcsv, tpcsv);
			lu->clut_locus(lu, locus, tpcsv, tdevv);	/* Compute locus values */
		} else {
			for (i = 0; i < p->out.chan; i++)		/* Target is K value */
				locus[i] = tdevv[i];
		}
#ifdef DEBUG
#ifdef DEBUGC
		DEBUGC
#endif
		printf("Got possible K %s of %f %f %f %f\n",p->out.locus ? "locus" : "value", locus[0],locus[1],locus[2],locus[3]);
#endif
	}

	/* Do PCS -> DevOut' */
	if (p->nhack == 3		/* All to K only */
	 || ntrig	 			/* Neutral or K only to K only hack has triggered */
	 || cmytrig) {			/* 100% CMY rough hack has triggered */

		if (p->nhack == 3 || ntrig) { /* Neutral to K only hack has triggered */
			co pp;
			pp.p[0] = pcsvm[0];					/* Input L value */
			p->pcs2k->interp(p->pcs2k, &pp);	/* L -> K' */
			if (pp.v[0] < 0.0)		/* rspl might have extrapolated */
				pp.v[0] = 0.0;
			else if (pp.v[0] > 1.0)
				pp.v[0] = 1.0;
			out[0] = out[1] = out[2] = 0.0;		/* We know output is CMYK' */
			out[3] = pp.v[0];

#ifndef DEBUG
			if (p->verb)
#endif
			if (ntrig) {
				printf("Neutral hack mapped %s to 0 0 0 %f\n", icmPdv(p->in.chan,win), out[3]); 
				fflush(stdout);
			}
		} else if (cmytrig) { /* 100% CMY rough hack has triggered */
			if (cmytrig & 1) {
				out[0] = 1.0;
				out[1] = out[2] = out[3] = 0.0;
			}
			if (cmytrig & 2) {
				out[1] = 1.0;
				out[0] = out[2] = out[3] = 0.0;
			}
			if (cmytrig & 4) {
				out[2] = 1.0;
				out[0] = out[1] = out[3] = 0.0;
			}

#ifndef DEBUG
			if (p->verb)
#endif
			if (cmytrig != 0) {
				if (p->in.chan == 4) 
					printf("CMY hack mapped %s to %s\n",icmPdv(p->in.chan, win), icmPdv(p->out.chan, out));
				fflush(stdout);
			}
		}
	} else {	/* Neutral to K hack has NOT triggered */
		switch(p->out.alg) {
		    case icmMonoBwdType: {
				icxLuMono *lu = (icxLuMono *)p->out.luo;	/* Safe to coerce */

				rv |= lu->bwd_abs(lu, pcsvm, pcsvm);
				rv |= lu->bwd_map(lu, out, pcsvm);
				if (p->out.nocurve) {	/* No explicit curve, so do implicit here */
					rv |= lu->bwd_curve(lu, out, out);
				}
				break;
			}
		    case icmMatrixBwdType: {
				icxLuMatrix *lu = (icxLuMatrix *)p->out.luo;	/* Safe to coerce */

				rv |= lu->bwd_abs(lu, pcsvm, pcsvm);
				rv |= lu->bwd_matrix(lu, out, pcsvm);
				if (p->out.nocurve) {	/* No explicit curve, so do implicit here */
					rv |= lu->bwd_curve(lu, out, out);
				}
				break;
			}
		    case icmLutType: {
				icxLuLut *lu = (icxLuLut *)p->out.luo;	/* Safe to coerce */

				if (p->mode < 2) {	/* Using B2A table */
					rv |= lu->in_abs(lu, pcsvm, pcsvm);
					rv |= lu->matrix(lu, pcsvm, pcsvm);
					rv |= lu->input(lu, pcsvm, pcsvm);
					rv |= lu->clut(lu, out, pcsvm);
					if (p->out.nocurve) {	/* No explicit curve, so do implicit here */
						rv |= lu->output(lu, out, out);
					}

				} else {	/* Use inverse A2B table */
					int i;
#ifdef USE_MERGE_CLUT_OPT
					/* Because we have used the ICX_MERGE_CLUT flag, we don't need */
					/* to call inv_out_abs() and inv_output() */
#else
					rv |= lu->inv_out_abs(lu, pcsvm, pcsvm);
					rv |= lu->inv_output(lu, pcsvm, pcsvm);
#endif

#ifdef DEBUG
#ifdef DEBUGC
					DEBUGC
#endif
					printf("Calling inv_clut with K aux targets %f %f %f %f and pcsvm %f %f %f %f\n",
					locus[0],locus[1],locus[2],locus[3],pcsvm[0],pcsvm[1],pcsvm[2],pcsvm[3]);
#endif

					/* locus[] contains possible K target or locus value, */
					/* so copy it to out[] so that inv_clut will use it. */
					for (i = 0; i < p->out.chan; i++)
						out[i] = locus[i];

					rv |= lu->inv_clut(lu, out, pcsvm);
#ifdef DEBUG
#ifdef DEBUGC
					DEBUGC
#endif
					printf("Got result %f %f %f %f\n", out[0],out[1],out[2],out[3]);
#endif

					
					if (p->out.nocurve) {	/* No explicit curve, so do implicit here */
						rv |= lu->inv_input(lu, out, out);
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

	if (p->out.lcurve) 		/* Apply Y to L* */
		y2l_curve(out, out, p->out.lcurve == 2);

#ifdef DEBUG
#ifdef DEBUGC
	DEBUGC
#endif
	printf("DevIn'->DevOut' ret %f %f %f %f\n\n",out[0], out[1], out[2], out[3]); fflush(stdout);
#endif


	if (p->verb) {		/* Output percent intervals */
		int pc;
		p->count++;
		pc = (int)(p->count * 100.0/p->total + 0.5);
		if (pc != p->last) {
			printf("\r%2d%%",pc); fflush(stdout);
			p->last = pc;
		}
	}
}

/* - - - - - - - - - - - - - - - - */
/* Output table, DevOut' -> DevOut */
void devop_devo(void *cntx, double *out, double *in) {
	int rv = 0;
	link *p = (link *)cntx;
	int i;

#ifdef DEBUG
#ifdef DEBUGC
	DEBUGC
#endif
	printf("DevOut'->DevOut got %f %f %f %f\n",in[0], in[1], in[2], in[4]); fflush(stdout);
#endif

	for (i = 0; i < p->out.chan; i++)
		out[i] = in[i];

	if (p->out.lcurve) 		/* Apply L* to Y */
		l2y_curve(out, out, p->out.lcurve == 2);

	if (p->out.nocurve == 0) {	/* Using per channel curves */

		switch(p->out.alg) {
		    case icmMonoBwdType: {
				icxLuMono *lu = (icxLuMono *)p->out.luo;		/* Safe to coerce */
				rv |= lu->bwd_curve(lu, out, out);
				break;
			}
		    case icmMatrixBwdType: {
				icxLuMatrix *lu = (icxLuMatrix *)p->out.luo;	/* Safe to coerce */
				rv |= lu->bwd_curve(lu, out, out);
				break;
			}
		    case icmLutType: {
				if (p->mode < 2) {	/* Using B2A table */
					icxLuLut *lu = (icxLuLut *)p->out.luo;	/* Safe to coerce */
					rv |= lu->output(lu, out, out);
					/* Since not PCS, out_abs is never used */
					break;
				} else {	/* Use inverse A2B table */
					icxLuLut *lu = (icxLuLut *)p->out.luo;	/* Safe to coerce */
					rv |= lu->inv_input(lu, out, out);
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
#ifdef DEBUGC
	DEBUGC
#endif
	printf("DevOut'->DevOut ret %f %f %f %f\n",out[0], out[1], out[2], out[3]); fflush(stdout);
#endif
#ifdef DEBUGC
	tt = 0;
#endif
}

/* ------------------------------------------- */
/* Fixup L -> K only lookup table white and black values, */
/* to compensate for inexact rspl fitting */

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

	/* Scale so that kmin->hmax becomes 0 to 1 */
	f = (out[0] - p->kmin)/(p->kmax - p->kmin);

	out[0] = f;
}

/* ------------------------------------------- */
/* powell() callback to set XYZ scaling factor */

static double xyzoptfunc(void *cntx, double *v) {
	link *p = (link *)cntx;
	double swxyz[3], jab[3], dev[MAX_CHAN];
	double rv;
	int rc = 0;

	rv = 2.0 - v[0];	/* Make Y as large as possible */

	/* If we wanted to use this function to maximise the brightness */
	/* we would not limit the scale to 1.0 */
	if (v[0] > 1.0) {
		rv += 1000.0;
		return rv;
	}
	if (v[0] < 0.0) {
		rv += 100.0;
		return rv;
	}
	swxyz[0] = v[0] * p->swxyz[0];
	swxyz[1] = v[0] * p->swxyz[1];
	swxyz[2] = v[0] * p->swxyz[2];

//printf("~1 scaled white XYZ = %f %f %f\n", swxyz[0], swxyz[1], swxyz[2]);

	if (p->pcsor == icxSigJabData) {
		/* We're being bad in delving inside the xluo, but we'll fix it latter */
		p->out.luo->cam->XYZ_to_cam(p->out.luo->cam, jab, swxyz);
	} else
		error("Internal :- not setup to handle Y scaling and non-Jab PCS");

//printf("~1 scaled white Jab = %f %f %f\n", jab[0], jab[1], jab[2]);

	/* Run the target PCS backwards through the output space to see if it clips */
	switch(p->out.alg) {
	    case icmMonoBwdType: {
			icxLuMono *lu = (icxLuMono *)p->out.luo;	/* Safe to coerce */

			rc = lu->bwd_lookup(p->out.luo, dev, jab);
			break;
		}
	    case icmMatrixBwdType: {
			icxLuMatrix *lu = (icxLuMatrix *)p->out.luo;	/* Safe to coerce */

			rc = lu->bwd_lookup(p->out.luo, dev, jab);
			break;
		}
	    case icmLutType: {
			icxLuLut *lu = (icxLuLut *)p->out.luo;	/* Safe to coerce */

			if (p->mode < 2)	/* Using B2A table */
				rc = lu->lookup(p->out.luo, dev, jab);
			else				/* Use inverse A2B table */
				rc = lu->inv_lookup(p->out.luo, dev, jab);
			break;
		}
		default:
			error("Unexpected algorithm type %d in devop of devip_devop()",p->out.alg);
	}
//printf("~1 device = %f %f %f, rc = %d\n", dev[0], dev[1], dev[2], rc);
	if (rc != 0)
		rv += 500.0;

//printf("~1 xyzoptfunc rv %f from xyzscale %f\n\n",rv,v[0]);
	return rv;
}

/* ------------------------------------------- */

int
main(int argc, char *argv[]) {
	int fa, nfa, mfa;				/* argument we're looking at */
	char in_name[MAXNAMEL+1];
	char sgam_name[MAXNAMEL+1] = "\000";	/* Source gamut name */
	char abs_name[MAXNAMEL+1] = "\000";		/* Abstract profile name */
	char out_name[MAXNAMEL+1];
	char link_name[MAXNAMEL+1];
	int verify = 0;				/* Do verify pass */
	int outinkset = 0;			/* The user specfied an output inking */
	int intentset = 0;			/* The user specified an intent */
	int vcset = 0;				/* Viewing conditions were set by user */
	int modeset = 0;			/* The gamut mapping mode was set by the user */
	int rv = 0;
	icxViewCond ivc, ovc;		/* Viewing Condition Overrides for in and out profiles */
	int ivc_e = -1, ovc_e = -1;	/* Enumerated viewing condition */
	link li;					/* Linking information structure */
	int isJab = 0;				/* (Derived from li.mode & li.gmi) NZ if Jab link space */
	int in_curve_res = 0;		/* Input profile A2B input curve resolution (if known) */
	int out_curve_res = 0;		/* Output profile B2A output curve resolution (if known) */
	profxinf xpi;				/* Extra profile information */
	int i;

	error_program = argv[0];
	memset((void *)&xpi, 0, sizeof(profxinf));	/* Init extra profile info to defaults */
	memset((void *)&li, 0, sizeof(link));

	/* Set defaults */
	li.verb = 0;
	li.count = 0;
	li.last = -1;
	li.mode = 0;						/* Default simple link mode */
	li.quality = 1;						/* Medium quality */
	li.clutres = 0;						/* No resolution override */	
	li.nhack = 0;
	li.cmyhack = 0;						/* Mask for 100% purity through mapping of CMY */
	li.pcs2k = NULL;
	li.wphack = 0;
	li.wphacked = 0;
	li.abs_luo = NULL;					/* No abstract */
	li.xyzscale = 1.0;					/* No XYZ scaling */
	li.hwp[0] = li.hwp[1] = li.hwp[2] = 0.0;
	li.map = NULL;
	li.Kmap = NULL;
	li.in.intent  = icmDefaultIntent;	/* Default */
	li.in.ink.tlimit = -1.0;			/* Default no total limit */
	li.in.ink.klimit = -1.0;			/* Default no black limit */
	li.in.inking  = 4;					/* Inking algorithm default = ramp */
	li.in.locus   = 0;					/* Default K value target */
	li.in.nocurve = 0;					/* Preserve device linearisation curve */
	li.in.lcurve = 0;					/* Don't apply a Y to L* curve after device curve */
	li.out.intent = icmDefaultIntent;	/* Default */
	li.out.ink.tlimit = -1.0;			/* Default no total limit */
	li.out.ink.klimit = -1.0;			/* Default no black limit */
	li.out.ink.KonlyLmin = 0;			/* Use normal black Lmin for locus */
	li.out.ink.c.Ksmth = ICXINKDEFSMTH;	/* default curve smoothing */
	li.out.ink.c.Kskew = ICXINKDEFSKEW;	/* default curve skew */
	li.out.ink.x.Ksmth = ICXINKDEFSMTH;
	li.out.ink.x.Kskew = ICXINKDEFSKEW;
	li.out.inking  = 4;					/* Default ramp K */
	li.out.locus   = 0;					/* Default K value target */
	li.out.nocurve = 0;					/* Preserve device linearisation curve */
	li.out.lcurve = 0;					/* Don't apply an L* to Y curve before device curve */
	li.out.b2aluo = NULL;				/* B2A lookup for inking == 7 */

	xicc_enum_gmapintent(&li.gmi, icxDefaultGMIntent, NULL);	/* Set default overall intent */

	/* Init VC overrides so that we know when the've been set */
	ivc.Ev = -1;
	ivc.Wxyz[0] = -1.0; ivc.Wxyz[1] = -1.0; ivc.Wxyz[2] = -1.0;
	ivc.La = -1.0;
	ivc.Yb = -1.0;
	ivc.Lv = -1.0;
	ivc.Yf = -1.0;
	ivc.Fxyz[0] = -1.0; ivc.Fxyz[1] = -1.0; ivc.Fxyz[2] = -1.0;

	ovc.Ev = -1;
	ovc.Wxyz[0] = -1.0; ovc.Wxyz[1] = -1.0; ovc.Wxyz[2] = -1.0;
	ovc.La = -1.0;
	ovc.Yb = -1.0;
	ovc.Lv = -1.0;
	ovc.Yf = -1.0;
	ovc.Fxyz[0] = -1.0; ovc.Fxyz[1] = -1.0; ovc.Fxyz[2] = -1.0;

	if (argc < 4)
		usage("Too few arguments, got %d expect at least 3",argc-1);

	/* Process the arguments */
	mfa = 3;        /* Minimum final arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */

		if (argv[fa][0] == '-')	{	/* Look for any flags */
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
				usage("Requested usage");

			/* Verbosity */
			else if (argv[fa][1] == 'v') {
				li.verb = 1;
			}

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

			/* Verify rather than link */
			else if (argv[fa][1] == 'V')
				verify = 1;

			/* Disable profile per channel curve use in device link output */
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N') {
				li.in.nocurve = 1;
				li.out.nocurve = 1;
			}

			/* Hack to force input neutrals to K only output */
			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F') {

				if (argv[fa][1] == 'f') {
					if (na != NULL) {		/* XXXK -> XXXK hack */
						int j;
						fa = nfa;
						for (j = 0; ; j++) {
							if (na[j] == '\000')
								break;
							if (na[j] == 'k' || na[j] == 'K')
								li.nhack = 2;
							else if (na[j] == 'c' || na[j] == 'C')
								li.cmyhack |= 0x1;
							else if (na[j] == 'm' || na[j] == 'M')
								li.cmyhack |= 0x2;
							else if (na[j] == 'y' || na[j] == 'Y')
								li.cmyhack |= 0x4;
							else
								usage("Unexpected argument '%c' to -f flag",na[j]);
						}

					} else {				/* Neutral -> 000K hack */
						li.nhack = 1;
						li.in.nocurve = 1;	/* Disable input curve to preserve input equality */
					}
				} else {
					li.nhack = 3;			/* All -> 000K Hack */
				}
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
			else if (argv[fa][1] == 'p') {
				if (na == NULL) usage("Expected abstract profile filename after -a");
				fa = nfa;
				strncpy(abs_name,na,MAXNAMEL); abs_name[MAXNAMEL] = '\000';
			}

			/* Simple mode */
			else if (argv[fa][1] == 's') {
				li.mode = 0;
				modeset = 1;
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
				modeset = 1;
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
			}
			/* Input profile Intent or Mapping mode intent */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL) usage("Input intent flag (-i) needs an argument");
				/* Record it for simple mode */
    			switch (na[0]) {
					case 'p':
					case 'P':
						li.in.intent = icPerceptual;
						break;
					case 'r':
					case 'R':
						li.in.intent = icRelativeColorimetric;
						break;
					case 's':
					case 'S':
						li.in.intent = icSaturation;
						break;
					case 'a':
					case 'A':
						li.in.intent = icAbsoluteColorimetric;
						break;
					default:
						li.in.intent = icMaxEnumIntent;		/* Detect error later */
				}
				/* Record it for gamut mapping mode */
				if (xicc_enum_gmapintent(&li.gmi, icxNoGMIntent, na) == -999)
					usage("Input intent (-i) argument '%s' isn't recognised",na);
				intentset = 1;
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
#ifdef NEVER
				if (na[0] >= '0' && na[0] <= '9') {
					if (vc == &ivc)
						ivc_e = atoi(na);
					else
						ovc_e = atoi(na);
				} else
#endif
				if (na[1] != ':') {
					if (vc == &ivc) {
						if ((ivc_e = xicc_enum_viewcond(NULL, NULL, -2, na, 1, NULL)) == -999)
							usage("Unrecognised viewing condition enumeration '%s'",na);
					} else {
						if ((ovc_e = xicc_enum_viewcond(NULL, NULL, -2, na, 1, NULL)) == -999)
							usage("Unrecognised viewing condition enumeration '%s'",na);
					}
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
				vcset = 1;			/* Viewing conditions were set by user */
			}

			/* Inking rule */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				fa = nfa;
				if (na == NULL) usage("Inking rule flag (-k) needs an argument");
				if (argv[fa][1] == 'k')
					li.out.locus = 0;			/* Use K value target */
				else
					li.out.locus = 1;			/* Use K locus target */
    			switch (na[0]) {
					case 't':
					case 'T':
						li.out.inking = 0;		/* Use input K value for output */
						break;
					case 'e':
					case 'E':
						li.out.inking = 7;		/* Use output K value as guide */
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
				outinkset = 1;		/* The user set an inking */
			}
			/* Input ink limits */
			else if (argv[fa][1] == 't') {
				int tlimit;
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -t");
				tlimit = atoi(na);
				if (tlimit >= 0)
					li.in.ink.tlimit = tlimit/100.0;
				else
					li.in.ink.tlimit = -1.0;
			}
			else if (argv[fa][1] == 'T') {
				int klimit;
				fa = nfa;
				if (na == NULL) usage("No parameter after flag -T");
				klimit = atoi(na);
				if (klimit >= 0)
					li.in.ink.klimit = klimit/100.0;
				else
					li.in.ink.klimit = -1.0;
			}
			/* Output ink limits */
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

			/* Gammut mapping diagnostic plots */
			else if (argv[fa][1] == 'P')
				li.gamdiag = 1;

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

	if (xpi.profDesc == NULL)
		xpi.profDesc = link_name;	/* Default description */

	if (li.verb)
		printf("Got options\n");

	/* - - - - - - - - - - - - - - - - - - - */
#ifndef ENKHACK				/* Enable K hack code */
	warning("!!!!!! linkl/collink.c ENKHACK not enabled !!!!!!");
#endif
	/* - - - - - - - - - - - - - - - - - - - */
	/* Sanity checking/defaulting of options */

	/* Deal with options that need link mode -g */
	if (li.mode < 1
	 && (li.in.intent == icMaxEnumIntent		/* User set a smart linking intent */
	  || vcset								/* Viewing conditions were set by user */
	  || li.wphack)) {
		if (modeset) {
			if (li.in.intent == icMaxEnumIntent)
				warning("Complex intent can't work with -s linking mode");
			else if (vcset)
				warning("Viewing conditions are ignored with -s linking mode");
			else if (li.wphack)
				warning("White point hack is ignored with -s linking mode");
		} else  {
			if (li.verb) {
				if (li.in.intent == icMaxEnumIntent)
					printf("Setting -g to enable Gamut Mapping mode intent\n");
				else if (vcset)
					printf("Setting -g to enable viewing conditions\n");
				else if (li.wphack)
					printf("Setting -g to enable white point hack\n");
			}
			li.mode = 1;
		}
	}

	/* Deal with options that need link mode -G */
	if (li.mode < 2
	  && (outinkset						/* The user set a K inking rule */
	   || li.out.ink.tlimit >= 0.0		/* The user set an output total limit */
	   || li.out.ink.klimit >= 0.0)) {	/* The user set an output black limit */
		if (modeset) {
			if (outinkset)
				warning("Black inking can't work with -s or -g linking mode");
			else if (li.out.ink.tlimit >= 0.0 || li.out.ink.klimit >= 0.0)
				warning("Ink limiting can't work with -s linking mode");
		} else {
			if (li.verb) {
				if (outinkset)
					printf("Setting -G to enable black inking\n");
				else if (li.out.ink.tlimit >= 0.0 || li.out.ink.klimit >= 0.0)
					printf("Setting -G to enable ink limiting\n");
			}
			li.mode = 2;
		}
	}

	/* Deal with options that complement -f -F */
	if (li.nhack || li.cmyhack) {

		/* Ideally we need to set K inking and map to K only black point, which require -G mode */
		if (li.mode < 2) {
			if (li.nhack == 1) {		/* All neutrals to K only */
				if (modeset) {
					warning("-f will give best result with -G mode");
				} else {
					if (li.verb)
						printf("Setting -G mode to complement -f option\n");
					li.mode = 2;
				}
			} else if (li.nhack == 2) {	/* K only in to K only out */
				if (modeset) {
					warning("For better results use -G mode with -fk option");
				} else {
					if (li.verb)
						printf("Setting -G mode to complement -fk option\n");
					li.mode = 2;
				}
			} else if (li.nhack == 3) {	/* All to K only out */
				if (modeset) {
					warning("For better results use -G mode with -F option");
				} else {
					if (li.verb)
						printf("Setting -G mode to complement -F option\n");
					li.mode = 2;
				}
			}
			if (li.cmyhack != 0) {	/* Map pure 100% CMY to pure CMY */
				if (modeset) {
					warning("For better results use -G mode with -fcmy options");
				} else {
					if (li.verb)
						printf("Setting -G mode to complement -fcmy options\n");
					li.mode = 2;
				}
			}
		}

		/* Ideally we should use an appropriate K inking */ 
		if (li.mode >= 2) {											/* Gammut mapping mode */
			if (li.nhack == 1 && li.out.inking != 3) {		/* All neutrals to K only */
				if (outinkset) {
					warning("For better results use -kx with -f option");
				} else {
					if (li.verb)
						printf("Setting -kx to complement -f option\n");
					li.out.inking = 3;			/* Use maximum K */
				}
			} else if (li.nhack == 2 && li.out.inking != 0) {	/* K only in to K only out */
				if (outinkset) {
					warning("For better results use -kt with -fk option");
				} else {
					if (li.verb)
						printf("Setting -kt to complement -fk option\n");
					li.out.inking = 0;		/* Use input K value for output */
				}
			} else if (li.nhack == 3 && li.out.inking != 3) {		/* All colors to K only */
				if (modeset) {
					warning("For better results use -kx with -F option");
				} else {
					if (li.verb)
						printf("Setting -kx to complement -f option\n");
					li.out.inking = 3;			/* Use maximum K */
				}
			}
		}

		/* Ideally we should use an appropriate gamut mapping */
		if (li.mode >= 1) {											/* Gammut mapping mode */

			if (li.gmi.usemap == 0 || li.gmi.greymf < 1.0			/* Not mapping black point */
		     || li.gmi.glumbcpf < 1.0 || li.gmi.glumbexf < 1.0) {
				if (li.nhack == 1) {		/* All neutrals to K only */
					if (intentset) {
						warning("For better results use an intent that maps black point with -f option");
					} else {
						if (li.verb)
							printf("Setting -ip intent to complement -f option\n");
						if (xicc_enum_gmapintent(&li.gmi, icxNoGMIntent, "p") == -999)
							usage("Internal, intent 'p' isn't recognised");
						li.dst_kbp = 1;				/* Map to K only black point */
						li.out.ink.KonlyLmin = 1;	/* Use K only black Lmin for locus */
					}
				} else if (li.nhack == 3) {	/* All to K only out */
					if (intentset) {
						warning("For better results use an intent that maps black point with -F option");
					} else {
						if (li.verb)
							printf("Setting -ip intent to complement -F option\n");
						if (xicc_enum_gmapintent(&li.gmi, icxNoGMIntent, "p") == -999)
							usage("Internal, intent 'p' isn't recognised");
						li.dst_kbp = 1;				/* Map to K only black point */
						li.out.ink.KonlyLmin = 1;	/* Use K only black Lmin for locus */
					}
				}

			/* Got an appropriate intent, so set mapping to K only black point */
			} else if (li.nhack == 1 || li.nhack == 3) {
				li.dst_kbp = 1;				/* Map to K only black point */
				li.out.ink.KonlyLmin = 1;	/* Use K only black Lmin for locus */
			}
			if (li.cmyhack != 0) {	/* Map pure 100% CMY to pure CMY */
				if (intentset) {
					if (strcmp(li.gmi.as, "s") != 0)
						warning("For better results use -is with -fcmy options");
				} else {
					if (li.verb)
						printf("Setting -is intent to complement -fcmy options\n");
					if (xicc_enum_gmapintent(&li.gmi, icxNoGMIntent, "s") == -999)
						usage("Internal, intent 's' isn't recognised");
					li.dst_cmymap = li.cmyhack;
				}
			}
		}
	}

	if (li.mode == 0) {
		if (li.in.intent == icMaxEnumIntent)
			usage("Input intent (-i) argument isn't recognised for simple mapping mode");
	}

	if (li.wphack && (li.gmi.usecas & 0x100) != 0)
		usage("Can't use 'white point hack' and Luminence scaling intent together");

	/* - - - - - - - - - - - - - - - - - - - */
	/* Open up the input device profile for reading, and read header etc. */
	if ((li.in.c = read_embedded_icc(in_name)) == NULL)
		error ("Can't open file '%s'",in_name);
	li.in.h = li.in.c->header;

	/* Check that it is a suitable device input icc */
	if (li.in.h->deviceClass != icSigInputClass
	 && li.in.h->deviceClass != icSigDisplayClass
	 && li.in.h->deviceClass != icSigOutputClass
	 && li.in.h->deviceClass != icSigColorSpaceClass)	/* For sRGB etc. */
		error("Input profile '%s' isn't a device profile",in_name);

	/* Wrap with an expanded icc */
	if ((li.in.x = new_xicc(li.in.c)) == NULL)
		error ("Creation of input profile xicc failed");

	/* Set the default ink limits if not set on command line */
	icxDefaultLimits(li.in.x, &li.in.ink.tlimit, li.in.ink.tlimit, &li.in.ink.klimit, li.in.ink.klimit);

	if (li.verb) {
		if (li.in.ink.tlimit >= 0.0)
			printf("Input total ink limit assumed is %3.0f%%\n",100.0 * li.in.ink.tlimit);
		if (li.in.ink.klimit >= 0.0)
			printf("Input black ink limit assumed is %3.0f%%\n",100.0 * li.in.ink.klimit);
	}

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
	/* Open up the output device output profile for reading, and read header etc. */
	if ((li.out.c = read_embedded_icc(out_name)) == NULL)
		error ("Can't open file '%s'",out_name);
	li.out.h = li.out.c->header;

	if (li.out.h->deviceClass != icSigInputClass
	 && li.out.h->deviceClass != icSigDisplayClass
	 && li.out.h->deviceClass != icSigOutputClass
	 && li.out.h->deviceClass != icSigColorSpaceClass)	/* For sRGB etc. */
		error("Output profile isn't a device profile");

	/* Wrap with an expanded icc */
	if ((li.out.x = new_xicc(li.out.c)) == NULL)
		error ("Creation of output profile xicc failed");

	/* Set the default ink limits if not set on command line */
	icxDefaultLimits(li.out.x, &li.out.ink.tlimit, li.out.ink.tlimit, &li.out.ink.klimit, li.out.ink.klimit);

	if (li.verb) {
		if (li.out.ink.tlimit >= 0.0)
			printf("Output total ink limit assumed is %3.0f%%\n",100.0 * li.out.ink.tlimit);
		if (li.out.ink.klimit >= 0.0)
			printf("Output black ink limit assumed is %3.0f%%\n",100.0 * li.out.ink.klimit);
	}

	/* deal with output black generation. */
	/* Ink limits will have been set in option parsing */

	switch (li.out.inking) {
		case 0:			/* Use input profile K level or locus */
			/* Sanity check */
			if (li.in.h->colorSpace != li.out.h->colorSpace)
				error("Can't transfer black ink in & out unless the same colorspaces");
			li.out.ink.k_rule = li.out.locus ? icxKlocus : icxKvalue;	/* Given as aux parameter in PCS -> Device */
			break;
		case 7:			/* Use output profile K level or locus */
			li.out.ink.k_rule = li.out.locus ? icxKlocus : icxKvalue;	/* Given as aux parameter in PCS -> Device */
			break;
		case 1:			/* Minimum K */
			li.out.ink.k_rule = li.out.locus ? icxKluma5 : icxKluma5k;
			li.out.ink.c.Kstle = 0.0;
			li.out.ink.c.Kstpo = 0.0;
			li.out.ink.c.Kenpo = 1.0;
			li.out.ink.c.Kenle = 0.0;
			li.out.ink.c.Kshap = 1.0;
			break;
		case 2:			/* 0.5 K */
			li.out.ink.k_rule = li.out.locus ? icxKluma5 : icxKluma5k;
			li.out.ink.c.Kstle = 0.5;
			li.out.ink.c.Kstpo = 0.0;
			li.out.ink.c.Kenpo = 1.0;
			li.out.ink.c.Kenle = 0.5;
			li.out.ink.c.Kshap = 1.0;
			break;
		case 3:			/* Maximum K */
			li.out.ink.k_rule = li.out.locus ? icxKluma5 : icxKluma5k;
			li.out.ink.c.Kstle = 1.0;
			li.out.ink.c.Kstpo = 0.0;
			li.out.ink.c.Kenpo = 1.0;
			li.out.ink.c.Kenle = 1.0;
			li.out.ink.c.Kshap = 1.0;
			break;
		case 4:			/* Ramp K */
			li.out.ink.k_rule = li.out.locus ? icxKluma5 : icxKluma5k;
			li.out.ink.c.Kstle = 0.0;
			li.out.ink.c.Kstpo = 0.0;
			li.out.ink.c.Kenpo = 1.0;
			li.out.ink.c.Kenle = 1.0;
			li.out.ink.c.Kshap = 1.0;
			break;
		case 5:			/* Curve */
			li.out.ink.k_rule = li.out.locus ? icxKluma5 : icxKluma5k;
			break;		/* Other params already set by options */
		case 6:			/* Use input profile K locus + dual curve limits */
			/* Sanity check */
			if (li.in.h->colorSpace != li.out.h->colorSpace)
				error("Can't transfer black ink in & out unless the same colorspaces");
			li.out.ink.k_rule = li.out.locus ? icxKl5l : icxKl5lk;	/* Aux param in PCS -> Device */
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
		xicc_enum_viewcond(x, vc, -1, NULL, 0, NULL);

		/* Override the viewing conditions */
		/* (?? Could move this code into xicc_enum_viewcond() as an option ??) */
		if (es != -1) {
			if (xicc_enum_viewcond(x, vc, es, NULL, 0, NULL) == -999)
				error ("%d, %s",x->errc, x->err);
		}
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
			vc->Yf = v->Yf;
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

	if (li.verb)
		printf("Configured options\n");

	/* - - - - - - - - - - - - - - - - - - - */
	/* Setup the profile color lookup information */
	{
		icmLuAlgType oalg;				/* Native output algorithm */
		icColorSpaceSignature natpcs;	/* Underlying native output PCS */
		int flb = 0, fl = 0;			/* luobj flags */
		
		li.pcsor = icSigLabData;			/* Default use Lab as PCS */

		/* If we are using the gamut map mode, then setup */
		/* the intents and pcsor appropriately. */
		if (li.mode > 0) {

			if ((li.gmi.usecas & 0xff) != 0) {
				li.pcsor = icxSigJabData;		/* Use CAM as PCS */
				isJab = 1;

				if ((li.gmi.usecas & 0xff) == 0x2) {	/* Absolute Appearance space */
					double mxw;
	
					li.in.intent = li.out.intent = li.abs_intent = icxAbsAppearance;

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
				} else {
					/* Not Abs Appearance space */
					li.in.intent = li.out.intent = li.abs_intent = icxAppearance;
				}
			} else {
				/* Not Appearance space */
				li.in.intent = li.out.intent = li.abs_intent = icAbsoluteColorimetric;
			}
		}

		if (li.verb)
			printf("Loading input A2B table\n");

		/* default flags for all xicc luobj's */
		flb = ICX_CLIP_NEAREST;
		if (li.verb)
			flb |= ICX_VERBOSE;

		/* Get an input profile xicc conversion object */
		fl = flb;
#ifdef USE_MERGE_CLUT_OPT
		fl |= ICX_MERGE_CLUT;
#endif

#ifdef NEVER
		printf("~1 input space flags = 0x%x\n",fl);
		printf("~1 input space intent = %s\n",icx2str(icmRenderingIntent,li.in.intent));
		printf("~1 input space pcs = %s\n",icx2str(icmColorSpaceSignature,li.pcsor));
		printf("~1 input space viewing conditions =\n"); xicc_dump_viewcond(&li.in.vc);
		printf("~1 input space inking =\n"); xicc_dump_inking(&li.in.ink);
#endif
		if ((li.in.luo = li.in.x->get_luobj(li.in.x, fl, icmFwd, li.in.intent,
		                                    li.pcsor, icmLuOrdNorm, &li.in.vc, &li.in.ink)) == NULL) {
			error("get xlookup object failed: %d, %s",li.in.x->errc,li.in.x->err);
		}
	
		/* Get details of overall conversion */
		li.in.luo->spaces(li.in.luo, &li.in.csp, &li.in.chan, NULL, NULL, &li.in.alg,
		                  NULL, NULL, NULL);

		/* Get the input profile A2B input curve resolution */
		/* (This is pretty rough - this should work for non LUT types as well!!) */
		{
			if (li.in.alg== icmLutType) {
				icmLut *lut;
				icxLuLut *luluo = (icxLuLut *)li.in.luo;		/* Safe to coerce */
				luluo->get_info(luluo, &lut, NULL, NULL, NULL);	/* Get some details */
				in_curve_res = lut->inputEnt;
			}
		}

		/* Grab the white point in case the wphack or xyzscale needs it */
		li.in.luo->efv_wh_bk_points(li.in.luo, li.in.wp, NULL, NULL);

		/* Get native PCS space */
		li.in.luo->lutspaces(li.in.luo, NULL, NULL, &natpcs, NULL, NULL);
	
		if (li.in.nocurve == 0 && natpcs == icSigXYZData
		  && (li.in.alg == icmMatrixFwdType || li.in.alg == icmMatrixBwdType
		      || li.in.csp == icSigXYZData)) {
			li.in.lcurve = 1;		/* Use Y to L* and L* to Y for input */
		    if (li.in.csp == icSigXYZData) {
				li.in.lcurve = 2;		/* Use real Y to L* and L* to Y for input */
				li.in.nocurve = 1;		/* Don't trust the curve that comes with it */
			}
			if (li.verb)
				printf("Using Y to L* and L* to Y curves for input\n");
		}

		/* Setup any abstract profile to match the chosen PCS */
		/* We aren't checking whether the input/abstract/output profile */
		/* intents really make any sense. It's assumed at the moment */
		/* that the user knows what they're doing! */
		if (abs_name[0] != '\000') {

			if ((li.abs_luo = li.abs_xicc->get_luobj(li.abs_xicc, flb, icmFwd, li.abs_intent,
			        li.pcsor, icmLuOrdNorm, &li.out.vc, NULL)) == NULL)
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
			plu->spaces(plu, NULL, NULL, NULL, NULL, &oalg, NULL, NULL, NULL, NULL);

			/* release the icm lookup */
			plu->del(plu);

		}

		if (oalg != icmLutType || li.mode < 2) {	/* Using B2A table or inv. mono/matrix */
			if (li.verb)
				printf("Loading output B2A table\n");

			if ((li.out.luo = li.out.x->get_luobj(li.out.x, flb, icmBwd, li.out.intent,
			                              li.pcsor, icmLuOrdNorm, &li.out.vc, &li.out.ink)) == NULL) {
				error("get xlookup object failed: %d, %s",li.out.x->errc,li.out.x->err);
			}
			/* Get details of overall conversion */
			li.out.luo->spaces(li.out.luo, NULL, NULL, NULL, &li.out.chan, &li.out.alg,
			                   NULL, NULL, NULL);

			/* Get the output profile B2A output curve resolution */
			/* (This is pretty rough - this should work for non LUT types as well!!) */
			{
				if (li.out.alg== icmLutType) {
					icmLut *lut;
					icxLuLut *luluo = (icxLuLut *)li.out.luo;		/* Safe to coerce */
					luluo->get_info(luluo, &lut, NULL, NULL, NULL);	/* Get some details */
					out_curve_res = lut->outputEnt;
				}
			}
	
			/* Grab the white point in case the wphack or xyzscale needs it */
			li.out.luo->efv_wh_bk_points(li.out.luo, li.out.wp, NULL, NULL);

			/* Get native PCS space */
			li.out.luo->lutspaces(li.out.luo, &natpcs, NULL, NULL, NULL, NULL);

		} else {	/* Using inverse A2B Lut for output conversion */

			fl = flb;
#ifdef USE_MERGE_CLUT_OPT
			fl |= ICX_MERGE_CLUT;
#endif
#ifdef USE_CAM_CLIP_OPT
			fl |= ICX_CAM_CLIP;
#endif
			if (li.verb)
				printf("Loading output inverse A2B table\n");

#ifdef NEVER
			printf("~1 output space flags = 0x%x\n",fl);
			printf("~1 output space intent = %s\n",icx2str(icmRenderingIntent,li.out.intent));
			printf("~1 output space pcs = %s\n",icx2str(icmColorSpaceSignature,li.pcsor));
			printf("~1 output space viewing conditions =\n"); xicc_dump_viewcond(&li.out.vc);
			printf("~1 output space inking =\n"); xicc_dump_inking(&li.out.ink);
#endif

			if ((li.out.luo = li.out.x->get_luobj(li.out.x, fl, icmFwd,
			                  li.out.intent, li.pcsor, icmLuOrdNorm, &li.out.vc,
			                  &li.out.ink)) == NULL) {
				error("get xlookup object failed: %d, %s",li.out.x->errc,li.out.x->err);
			}
		
			/* Get details of overall conversion */
			li.out.luo->spaces(li.out.luo, &li.out.csp, &li.out.chan, NULL, NULL, &li.out.alg,
			                   NULL, NULL, NULL);

			/* Get the output profile A2B input curve resolution */
			/* (This is pretty rough - this should work for non LUT types as well!!) */
			{
				if (li.out.alg== icmLutType) {
					icmLut *lut;
					icxLuLut *luluo = (icxLuLut *)li.out.luo;		/* Safe to coerce */
					luluo->get_info(luluo, &lut, NULL, NULL, NULL);	/* Get some details */
					out_curve_res = lut->inputEnt;
				}
			}

			/* Grab the white point in case the wphack or xyzscale needs it */
			li.out.luo->efv_wh_bk_points(li.out.luo, li.out.wp, NULL, NULL);

			/* Get native PCS space */
			li.out.luo->lutspaces(li.out.luo, NULL, NULL, &natpcs, NULL, NULL);

			/* If we need a B2A lookup to get the existing K */
			if (li.out.inking == 7) {
				if ((li.out.b2aluo = li.out.x->get_luobj(li.out.x, flb, icmBwd,
					li.out.intent, li.pcsor, icmLuOrdNorm, &li.out.vc, NULL)) == NULL) {
					error("get B2A xlookup object failed: %d, %s",li.out.x->errc,li.out.x->err);
				}
			}
		}
		
		/* If we need an PCS->K' mapping for the neutral axis to K hack. */
		/* What we do is lookup the L for K values from 0 to 1, */
		/* and then invert this to create an L to K lookup. */
		/* If the gamut mapping is set to map to the K only black point, */
		/* it should all work well... */
		if (li.nhack) {
			icxLuBase *luo;		/* Base XLookup type object */
			icmLuAlgType alg;	/* Type of lookup algorithm */
			co ips[256];		/* Initialisation points */
			datai glow;			/* Grid low scale */
			datai ghigh;		/* Grid high scale */
			datao vlow;			/* Data value low normalize */
			datao vhigh;		/* Data value high normalize */
			double Lmax, Lmin;	/* Max and Min L values that result */
			int grres;
			double avgdev[MXDO];

			if (li.out.h->colorSpace != icSigCmykData)
				error("Neutral Axis K only requested with non CMYK output profile");

			if (li.in.chan < 3)
				error("Neutral Axis K only requested with input profile with less than 3 channels");

			if (li.nhack == 2 && li.in.h->colorSpace != icSigCmykData)
				error("Neutral Axis 000K only requested with input profile that is not CMYK");

			if ((li.pcs2k = new_rspl(RSPL_NOFLAGS, 1, 1)) == NULL) {
				error("Failed to create an rspl object");
			}

			/* Get a device to PCS lookup object to use to lookup K->PCS */
			if ((luo = li.out.x->get_luobj(li.out.x, flb,
			                  icmFwd, li.out.intent, li.pcsor, icmLuOrdNorm, &li.out.vc,
			                  NULL)) == NULL) {
				error("get xlookup object failed: %d, %s",li.out.x->errc,li.out.x->err);
			}
			/* Get details of overall conversion */
			luo->spaces(luo, NULL, NULL, NULL, NULL, &alg, NULL, NULL, NULL);
			if (alg != icmLutType)
				error ("Unexpected algorithm type for CMYK output profile");

			/* Setup the initialisation points */
			Lmax = -100.0;
			Lmin = 1000.0;
			for (i = 0; i < 256; i++) {
				icxLuLut *lu = (icxLuLut *)luo;	/* Safe to coerce */
				double in[4], pcsv[4];
				in[0] = in[1] = in[2] = 0.0;
				in[3] = i/(255.0);

				/* Want to do dev' -> PCS conversion to match the */
				/* normal inverse PCS-> dev' used in devip_devop() */
				if (li.out.nocurve) {	/* No explicit curve, so do implicit here */
					/* Since not PCS, in_abs and matrix cannot be valid, */
					/* so input curve on own is ok to use. */
					lu->input(lu, pcsv, in);
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
				printf("L %f -> K' %f\n",pcsv[0], in[3]);
#endif /* NEUTKDEBUG */

				if (pcsv[0] > Lmax)			/* Track min and max L values */
					Lmax = pcsv[0];
				if (pcsv[0] < Lmin)
					Lmin = pcsv[0];
			}

			glow[0] = 0.0;
			ghigh[0] = 100.0;
			vlow[0] = 0.0;
			vhigh[0] = 1.0;
			grres = 256;
			avgdev[0] = 0.005;

			li.pcs2k->fit_rspl(li.pcs2k, 0, ips, 256, glow, ghigh, &grres, vlow, vhigh, 1.0, avgdev, NULL);

			/* Fixup the white and black points for neutral axis to K hack. */
			/* This is to make sure that they exactly match the fwd mapping */
			/* after the rspl is fitted. */
			{
				pcs2k_ctx cx;		/* White point fixup context */
				co pp;				/* Lookup the min and max K values */

				pp.p[0] = Lmax;
				li.pcs2k->interp(li.pcs2k, &pp);
				cx.kmin = pp.v[0];	/* Ideally would be 0 */
				pp.p[0] = Lmin;
				li.pcs2k->interp(li.pcs2k, &pp);
				cx.kmax = pp.v[0];	/* Ideally would be 1 */

#ifdef NEUTKDEBUG
				printf("Before fix: Lmax %f, Lmin %f, Kmin %f, Kmax %f\n",Lmax, Lmin, cx.kmin, cx.kmax);
#endif /* NEUTKDEBUG */

				li.pcs2k->re_set_rspl(li.pcs2k, 0, (void *)&cx, fix_pcs2k_white);
#ifdef NEUTKDEBUG
				pp.p[0] = Lmax;
				li.pcs2k->interp(li.pcs2k, &pp);
				cx.kmin = pp.v[0];	/* Ideally would be 0 */
				pp.p[0] = Lmin;
				li.pcs2k->interp(li.pcs2k, &pp);
				cx.kmax = pp.v[0];	/* Ideally would be 1 */
				printf("After fix: Lmax %f, Lmin %f, Kmin %f, Kmax %f\n",Lmax, Lmin, cx.kmin, cx.kmax);
#endif /* NEUTKDEBUG */
			}

		}	/* end if neutral axis to K hack */

		if (li.cmyhack != 0) {
			if (li.in.h->colorSpace != icSigCmyData
			 && li.in.h->colorSpace != icSigCmykData)
				error("100% CMY mapping requested with non CMY or CMYK input profile");

			if (li.out.h->colorSpace != icSigCmyData
			 && li.out.h->colorSpace != icSigCmykData)
				error("100% CMY mapping requested with non CMY or CMYK output profile");
		}

		if (li.out.nocurve == 0 && natpcs == icSigXYZData
		  && (li.out.alg == icmMatrixFwdType || li.out.alg == icmMatrixBwdType
		      || li.out.csp == icSigXYZData)) {
			li.out.lcurve = 1;		/* Use Y to L* and L* to Y for output */
		    if (li.out.csp == icSigXYZData) {
				li.out.lcurve = 2;		/* Use real Y to L* and L* to Y for output */
				li.out.nocurve = 1;		/* Don't trust the curve that comes with it */
			}
			if (li.verb)
				printf("Using Y to L* and L* to Y curves for output\n");
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
	/* for the latter. The xluo->get_gamut work in the set li.pcsor */
	/* for each xluo. */
	if (li.mode > 0 && li.gmi.usemap) {
		gamut *csgam, *igam, *ogam;
		double sgres;			/* Source gamut surface feature resolution */
		double dgres;			/* Destination gamut surface feature resolution */
		int    mapres;			/* Mapping rspl resolution */

		if (li.verb)
			printf("Creating Gamut Mapping\n");

		/* Gamut mapping will extend given grid res to encompas */
		/* source gamut by a margin. */
		if (li.quality == 3) {			/* Ultra High */
  	 		sgres = 7.0;
  	 		dgres = 7.0;
  	 		mapres = 41;
		} else if (li.quality == 2) {	/* High */
  	 		sgres = 8.0;
  	 		dgres = 8.0;
  	 		mapres = 33;
		} else if (li.quality == 1) {	/* Medium */
  	 		sgres = 10.0;
  	 		dgres = 10.0;
  	 		mapres = 25;
		} else {						/* Low quality */
  	 		sgres = 12.0;
  	 		dgres = 12.0;
  	 		mapres = 17;
		}

		/* Creat the source colorspace gamut surface */
		if (li.verb)
			printf(" Finding Source Colorspace Gamut with res %f\n",sgres);

		/* Creat the source image gamut surface in the selected li.pcsor space */
		if ((csgam = li.in.luo->get_gamut(li.in.luo, sgres)) == NULL)
			error ("%d, %s",li.in.x->errc, li.in.x->err);

		/* Grab a given source image gamut. */
		if (sgam_name[0] != '\000') {		/* Optional source gamut - ie. from an images */

			if (li.verb)
				printf(" Loading Image Source Gamut '%s'\n",sgam_name);

			igam = new_gamut(sgres, isJab, 0);	/* isJab will be overriden by gamut file */

			if (igam->read_gam(igam, sgam_name))
				error("Reading source gamut '%s' failed",sgam_name);

			if (igam->getisjab(igam) != isJab) {
				/* Should really convert to/from Jab here! */
				warning("Image gamut is wrong colorspace for link (Lab != Jab)");

				/* This will actually error in the gamut mapping code */
				/* Note that we're not checking relative/absolute colorspace here. */
				/* At the moment it's up to the user to get this right. */
			}

		} else {
			igam = NULL;	/* NULL signals no source image gamut */
		}

		/* Creat the destination gamut surface */
		if (li.verb)
			printf(" Finding Destination Gamut with res %f\n",dgres);

		if ((ogam = li.out.luo->get_gamut(li.out.luo, dgres)) == NULL)
			error ("%d, %s",li.out.x->errc, li.out.x->err);

		if (li.verb)
			printf(" Creating Gamut match\n");

		li.map = new_gammap(li.verb, csgam, igam, ogam, &li.gmi,
		                    li.src_kbp, li.dst_kbp, li.cmyhack, li.rel_oride,
		                    mapres, NULL, NULL, li.gamdiag ? "gammap.wrl" : NULL
		);
		if (li.map == NULL)
			error ("Failed to make gamut map transform");

		if (li.nhack == 2) {
			if (li.verb)
				printf(" Creating K only black to K only black Gamut match\n");

			li.Kmap = new_gammap(li.verb, csgam, igam, ogam, &li.gmi,
			                    1, 1, li.cmyhack, li.rel_oride,
			                    mapres, NULL, NULL, li.gamdiag ? "gammap.wrl" : NULL
			);
			if (li.Kmap == NULL)
				error ("Failed to make K only gamut map transform");
		}

		ogam->del(ogam);
		if (igam != NULL)
			igam->del(igam);
		csgam->del(csgam);
	}

	/* If we've got a request for Absolute Appearance mode with scaling */
	/* to avoid clipping the source white point, compute the needed XYZ scaling factor. */
	/* We assume that the white point hack can't be used at the same time. */
	if (li.mode > 0 && li.wphack == 0 && (li.gmi.usecas & 0x100) != 0) {
		double xyzscale[1], sa[1];

		/* We already have the source space white point in li.in.wp[] */

		/* Convert it to destination XYZ */
		if (li.pcsor == icxSigJabData) {
			/* We're being bad in delving inside the xluo, but we'll fix it latter */
			li.out.luo->cam->cam_to_XYZ(li.out.luo->cam, li.swxyz, li.in.wp);
		} else
			error("Internal :- not setup to handle Y scaling and non-Jab PCS");

//printf("~1 Source white Jab = %f %f %f\n", li.in.wp[0], li.in.wp[1], li.in.wp[2]);
//printf("~1 Source white XYZ = %f %f %f\n", li.swxyz[0], li.swxyz[1], li.swxyz[2]);

		/* Compute the bigest scale factor less than or equal to 1.0, */
		/* that doesn't clip the li.swxyz[] on the destination gamut */
		sa[0] = 0.1;
		xyzscale[0] = 0.5;
		if (powell(NULL, 1, xyzscale, sa, 1e-6, 2000, xyzoptfunc, (void *)&li, NULL, NULL) != 0) {
			warning("set_icxLuLut: XYZ scale powell failed to converge - set scale to 1.0");
		} else {
			li.xyzscale = xyzscale[0];
			if (li.verb)
				printf("Set XYZ scale factor to %f\n",li.xyzscale);
		}
	}

	/* Create the link profile */
	if (verify == 0) {
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
			if (li.mode > 0) {
				wh->renderingIntent = li.gmi.icci;			/* Closest ICC intent */
			} else {
				wh->renderingIntent = li.out.intent;		/* Output intent chosen */
			}

			/* Values that should be set before writing */
			if (xpi.manufacturer != 0L)
				wh->manufacturer = xpi.manufacturer;
			else
				wh->manufacturer = str2tag("????");
	
			if (xpi.model != 0L)
				wh->model = xpi.model;
			else
		    	wh->model        = str2tag("????");
	
			/* Values that may be set before writing */
			if (xpi.creator != 0L)
				wh->creator = xpi.creator;

			wh->attributes.l = 0;
			wh->flags        = 0;
#ifdef NT
			wh->platform = icSigMicrosoft;
#endif
#ifdef __APPLE__
			wh->platform = icSigMacintosh;
#endif
#if defined(UNIX) && !defined(__APPLE__)
			wh->platform = icmSig_nix;
#endif
		}
		/* Profile Description Tag: */
		{
			icmTextDescription *wo;
			char *dst, dstm[200];			/* description */
	
			if (xpi.profDesc != NULL)
				dst = xpi.profDesc;
			else {
				dst = "Device Link profile - See ProfileSequenceDescTag for more information";
				dst = dstm;
			}
	
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
			char *crt;
	
			if (xpi.copyright != NULL)
				crt = xpi.copyright;
			else
				crt = "Copyright, the creator of this profile";
	
			if ((wo = (icmText *)wr_icc->add_tag(
			           wr_icc, icSigCopyrightTag,	icSigTextType)) == NULL) 
				error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);
	
			wo->size = strlen(crt)+1; 	/* Allocated and used size of text, inc null */
			wo->allocate((icmBase *)wo);/* Allocate space */
			strcpy(wo->data, crt);		/* Copy the text in */
		}
		/* Device Manufacturers Description Tag: */
		if (xpi.deviceMfgDesc != NULL) {
			icmTextDescription *wo;
			char *dst = xpi.deviceMfgDesc;
	
			if ((wo = (icmTextDescription *)wr_icc->add_tag(
			           wr_icc, icSigDeviceMfgDescTag,	icSigTextDescriptionType)) == NULL) 
				error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);
	
			wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
			wo->allocate((icmBase *)wo);/* Allocate space */
			strcpy(wo->desc, dst);		/* Copy the string in */
		}
		/* Model Description Tag: */
		if (xpi.modelDesc != NULL) {
			icmTextDescription *wo;
			char *dst = xpi.modelDesc;
	
			if ((wo = (icmTextDescription *)wr_icc->add_tag(
			           wr_icc, icSigDeviceModelDescTag,	icSigTextDescriptionType)) == NULL) 
				error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);
	
			wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
			wo->allocate((icmBase *)wo);/* Allocate space */
			strcpy(wo->desc, dst);		/* Copy the string in */
		}
		/* ProfileSequenceDescTag: */
		{
			unsigned int i;
			icmProfileSequenceDesc *wo;
			if ((wo = (icmProfileSequenceDesc *)wr_icc->add_tag(
			           wr_icc, icSigProfileSequenceDescTag, icSigProfileSequenceDescType)) == NULL) 
				return 1;

			wo->count = 2; 		/* Number of descriptions in sequence */
			if (wo->allocate((icmBase *)wo) != 0)	/* Allocate space for all the DescStructures */
				error("allocate failed: %d, %s",wr_icc->errc,wr_icc->err);

			/* Fill in each description structure in sequence */

			/* For each profile in the chain */
			for (i = 0; i < wo->count; i++) {
				icc *iccs = NULL;
				icmHeader *sh = NULL;
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
					if (wo->data[i].allocate(&wo->data[i]) != 0)	/* Allocate space */
						error("allocate failed: %d, %s",wr_icc->errc,wr_icc->err);
					strcpy(wo->data[i].device.desc, ddesc->desc);

					wo->data[i].device.ucLangCode = ddesc->ucLangCode;
					wo->data[i].device.ucSize = ddesc->ucSize;
					if (wo->data[i].allocate(&wo->data[i]) != 0)	/* Allocate space */
						error("allocate failed: %d, %s",wr_icc->errc,wr_icc->err);
					memcpy(wo->data[i].device.ucDesc, ddesc->ucDesc, 2 * ddesc->ucSize);

					wo->data[i].device.scCode = ddesc->scCode;	
					wo->data[i].device.scSize = ddesc->scSize;		
					strcpy((char *)wo->data[i].device.scDesc, (char *)ddesc->scDesc);
				}

				/* model Text description */
				if (mdesc != NULL) {
					wo->data[i].model.size = mdesc->size;
					if (wo->data[i].allocate(&wo->data[i])!= 0)		/* Allocate space */
						error("allocate failed: %d, %s",wr_icc->errc,wr_icc->err);
					strcpy(wo->data[i].model.desc, mdesc->desc);

					wo->data[i].model.ucLangCode = mdesc->ucLangCode;
					wo->data[i].model.ucSize = mdesc->ucSize;
					if (wo->data[i].allocate(&wo->data[i]) != 0) 	/* Allocate space */
						error("allocate failed: %d, %s",wr_icc->errc,wr_icc->err);
					memcpy(wo->data[i].model.ucDesc, mdesc->ucDesc, 2 * mdesc->ucSize);

					wo->data[i].model.scCode = mdesc->scCode;	
					wo->data[i].model.scSize = mdesc->scSize;		
					strcpy((char *)wo->data[i].model.scDesc, (char *)mdesc->scDesc);
				}
			}
		}
		/* ColorantTable: */
		{
			int i;
			unsigned int j;
			int repclip = 0;

			/* For the first and last profile in the chain */
			/* (Note that we're assuming that the link output PCS is always Lab) */ 
			for (i = 0; i < 2; i++) {
				icc *iccs;
				icmHeader *sh;
				icmColorantTable *ro;
				icmColorantTable *wo;
				icTagSignature cts;

				if (i == 0) {
					iccs = li.in.c;		/* Input profile */
					sh = li.in.h;		/* Input profile header */
					cts = icSigColorantTableTag; 
				} else {
					iccs = li.out.c;	/* Output profile */
					sh = li.out.h;		/* Output profile header */
					cts = icSigColorantTableOutTag;
				}

				/* Try and read the input ColorantTable */
				if ((ro = (icmColorantTable *)iccs->read_tag(
			   	          iccs, icSigColorantTableTag)) != NULL) { 

					/* Create a ColorantTable in the output device link */
					if ((wo = (icmColorantTable *)wr_icc->add_tag(
			           wr_icc, cts, icSigColorantTableType)) == NULL) 
						error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

					/* Copy everything across */
			   		wo->count = ro->count;
					if (wo->allocate((icmBase *)wo) != 0)
						error("allocate failed: %d, %s",wr_icc->errc,wr_icc->err);

					for (j = 0; j < wo->count; j++) {
						strcpy(wo->data[j].name, ro->data[j].name);
						if (sh->pcs != icSigLabData) {
							icmXYZ2Lab(&icmD50, wo->data[j].pcsCoords, ro->data[j].pcsCoords);
							/* For device links the colorant table must be Lab PCS, */
							/* but embarassingly, XYZ profiles can have colorant values */
							/* not representable in the Lab PCS range. */
							if (icmClipLab(wo->data[j].pcsCoords, wo->data[j].pcsCoords)) {
								if (repclip)
									warning("Colorant Tag Lab value was clipped");
								repclip = 1;
							}
						} else {
							icmAry2Ary(wo->data[j].pcsCoords, ro->data[j].pcsCoords);
						}
					}

				} else {	/* Do this the hard way */
					icmLuBase *luo;
					unsigned int count;
					double dv[MAX_CHAN];
					double cvals[MAX_CHAN][3];
					inkmask imask;

					/* Get a lookup to read colorant values */
					if ((luo = iccs->get_luobj(iccs, icmFwd, icRelativeColorimetric,
					                           icSigLabData, icmLuOrdNorm)) == NULL)
						goto skip_coloranttable;

			   		count = icmCSSig2nchan(sh->colorSpace);
					for (j = 0; j < count; j++)
						dv[j] = 0.0;

					/* Lookup the colorant Lab values the recommended ICC way */
					for (j = 0; j < count; j++) {
						dv[j] = 1.0;
						luo->lookup(luo, cvals[j], dv);
						/* For device links the colorant table must be Lab PCS, */
						/* but embarassingly, XYZ profiles can have colorant values */
						/* not representable in the Lab PCS range. */
						if (icmClipLab(cvals[j], cvals[j])) {
							if (repclip)
								warning("Colorant Tag Lab value was clipped");
							repclip = 1;
						}
						dv[j] = 0.0;
					}
					luo->del(luo);

					/* Lookup colorant names */
					if ((imask = icx_icc_cv_to_colorant_comb(sh->colorSpace, iccs->header->deviceClass, cvals)) == 0)
						goto skip_coloranttable;

					/* Create a ColorantTable in the output device link */
					if ((wo = (icmColorantTable *)wr_icc->add_tag(
			           wr_icc, cts, icSigColorantTableType)) == NULL) 
						error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

			   		wo->count = count;
					if (wo->allocate((icmBase *)wo) != 0)
						error("allocate failed: %d, %s",wr_icc->errc,wr_icc->err);
					
					for (j = 0; j < count; j++) {
						inkmask iimask;			/* Individual ink mask */
						char *name;
						
						iimask = icx_index2ink(imask, j);
						name = icx_ink2string(iimask); 
						if (strlen(name) > 31)
							error("Internal: colorant name exceeds 31 characters");
						strcpy(wo->data[j].name, name);
						wo->data[j].pcsCoords[0] = cvals[j][0];
						wo->data[j].pcsCoords[1] = cvals[j][1];
						wo->data[j].pcsCoords[2] = cvals[j][2];
					}
				}
				/* Jump to here if we can't figure out what to put in ColorantTag */
				skip_coloranttable:;
			}
		}
		/* 16 bit input device -> output device lut: */
		{
			int inputEnt, outputEnt, clutPoints;
			icmLut *wo;

			/* Setup the cLUT resolutions */
			if (li.quality >= 3)
	    		inputEnt = 4096;
			else if (li.quality == 2)
	    		inputEnt = 2048;
			else
	    		inputEnt = 256;

			/* Make sure that we have at least the number of input entries as the */
			/* input profile. */
			if (in_curve_res > inputEnt)
				inputEnt = in_curve_res;

			/* See discussion in imdi/imdi_gen.c for ideal numbers */
			switch (li.in.chan) {
				case 0:
					error ("Illegal number of input chanels");
				case 1:
					if (li.quality >= 3)
			  		  	clutPoints = 255;
					else if (li.quality == 2)
			  		  	clutPoints = 255;
					else
			  		  	clutPoints = 255;
					break;

				case 2:
					if (li.quality >= 2)
		  		  		clutPoints = 255;
					else
		  		  		clutPoints = 86;
					break;
				case 3:
					if (li.quality >= 3)
		  		  		clutPoints = 52;
					else if (li.quality == 2)
		  		  		clutPoints = 33;
					else if (li.quality == 1)
		  		  		clutPoints = 17;
					else
		  		  		clutPoints = 9;
					break;
				case 4:
					if (li.quality >= 3)
		  		  		clutPoints = 33;
					else if (li.quality == 2)
		  		  		clutPoints = 18;
					else if (li.quality == 1)
		  		  		clutPoints = 9;
					else
		  		  		clutPoints = 6;
					break;
				case 5:
					if (li.quality >= 3)
		  		  		clutPoints = 18;
					else if (li.quality == 2)
		  		  		clutPoints = 16;
					else 
						clutPoints = 9;
					break;
				case 6:
					if (li.quality >= 3)
		  		  		clutPoints = 12;
					else if (li.quality == 2)
		  		  		clutPoints = 9;
					else 
						clutPoints = 6;
					break;
				case 7:
					if (li.quality >= 3)
		  		  		clutPoints = 8;
					else if (li.quality == 2)
		  		  		clutPoints = 7;
					else 
						clutPoints = 6;
					break;
				case 8:
					if (li.quality >= 3)
		  		  		clutPoints = 7;
					else if (li.quality == 2)
		  		  		clutPoints = 6;
					else 
						clutPoints = 5;
					break;
				default: /* > 8 chan */
					clutPoints = 3;
					break;
			}	

			if (li.clutres > 0)		/* clut resolution override */
  		  		clutPoints = li.clutres;
			li.clutres = clutPoints;	/* Actual resolution */

			if (li.quality >= 3)
	    		outputEnt = 4096;
			else if (li.quality == 2)
	    		outputEnt = 2048;
			else
	    		outputEnt = 256;

			/* Make sure that we have at least the number of input entries as the */
			/* output profile. */
			if (out_curve_res > outputEnt)
				outputEnt = out_curve_res;


			/* Link Lut = AToB0 */
			if ((wo = (icmLut *)wr_icc->add_tag(
			           wr_icc, icSigAToB0Tag,	icSigLut16Type)) == NULL) 
				error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

			wo->inputChan = li.in.chan;
			wo->outputChan = li.out.chan;

			/* Setup the tables resolutions */
    		wo->inputEnt = inputEnt;
  		  	wo->clutPoints = clutPoints;
	  		wo->outputEnt = outputEnt;

			if (wo->allocate((icmBase *)wo) != 0)	/* Allocate space */
				error("allocate failed: %d, %s",wr_icc->errc,wr_icc->err);

			/* Special case if input profile is Lut with matrix */
			/* (Does this do anything since input is not XYZ ?) */
			if (li.in.alg == icmLutType && li.in.nocurve == 0) {
				icxLuLut *lu = (icxLuLut *)li.in.luo;
				lu->get_matrix(lu, wo->e);		/* Copy it across */
			}



			if (li.verb)
				printf("Filling in Lut table\n");
#ifdef DEBUG_ONE
#define DBGNO 1		/* Up to 10 */

#ifndef NEVER
			/* Test a single given rgb/cmyk -> cmyk value */
			{
				double in[10][MAX_CHAN];
				double out[MAX_CHAN];
				in[0][0] = 1.0;
				in[0][1] = 1.0;
				in[0][2] = 1.0;
				in[0][3] = 0.0;

				in[1][0] = 1.0;
				in[1][1] = 1.0;
				in[1][2] = 1.0;
				in[1][3] = 0.0;

				for (i = 0; i < DBGNO; i++) {
					printf("Input %f %f %f %f\n",in[i][0], in[i][1], in[i][2], in[i][3]);
					devi_devip((void *)&li, out, in[i]);
					devip_devop((void *)&li, out, out);
					devop_devo((void *)&li, out, out);
					printf("Output %f %f %f %f\n\n",out[0], out[1], out[2], out[3]);
				}
			}
#endif /* NEVER */

#else	/* !DEBUG_ONE */
			/* Use helper function to do the hard work. */
			if (li.verb) {
				unsigned int ui;
				int itotal;
				for (itotal = 1, ui = 0; ui < li.in.chan; ui++, itotal *= clutPoints)
					; 
				li.total = itotal;
				/* Allow for extra lookups due to ICM_CLUT_SET_APXLS */
				for (itotal = 1, ui = 0; ui < li.in.chan; ui++, itotal *= (clutPoints-1))
					; 
				li.total += itotal;
				li.count = 0;
				printf(" 0%%"); fflush(stdout);
			}
			if (icmSetMultiLutTables(
				1,
				&wo,
				ICM_CLUT_SET_APXLS,			/* Use aproximate least squares */
				&li,						/* Context */
				li.in.h->colorSpace,		/* Input color space */
				li.out.h->colorSpace,		/* Output color space */
				devi_devip,					/* Input transfer tables devi->devi' */
				NULL, NULL,					/* Use default input colorspace range */
				devip_devop,				/* devi' -> devo' transfer function */
				NULL, NULL,					/* Default output colorspace range */
				devop_devo					/* Output transfer tables, devo'->devo */
			) != 0) {
				error("Setting 16 bit Lut failed: %d, %s",wr_icc->errc,wr_icc->err);
			}
			if (li.verb) {
				printf("\n");
			}
#ifdef WARN_CLUT_CLIPPING
			if (wr_icc->warnc)
				warning("Values clipped in setting device link LUT");
#endif /* WARN_CLUT_CLIPPING */

#endif	/* !DEBUG_ONE */

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

	/* - - - - - - - - - - - - - - - - - - - */
	/* Verify the given link, assuming all the options are the same */
	} else {
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
		luo->spaces(luo, &ins, &inn, &outs, &outn, &alg, NULL, NULL, NULL, NULL);

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
			pc = (int)(count * 100.0/total + 0.5);
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
	}

	/* - - - - - - - - - - - - - - - - - - - */
	/* Cleanup source profiles and exit */

	if (li.pcs2k != NULL)			/* Free up PCS->K lookup for neutral hack */
		li.pcs2k->del(li.pcs2k);

	if (li.map != NULL)
		li.map->del(li.map);
	if (li.Kmap != NULL)
		li.Kmap->del(li.Kmap);

	if (li.abs_luo != NULL) {		/* Free up abstract transform */
		li.abs_luo->del(li.abs_luo);
		li.abs_xicc->del(li.abs_xicc);
		li.abs_icc->del(li.abs_icc);
		li.abs_fp->del(li.abs_fp);
	}

	li.in.luo->del(li.in.luo);
	li.in.x->del(li.in.x);
	li.in.c->del(li.in.c);

	if (li.out.b2aluo != NULL)
		li.out.b2aluo->del(li.out.b2aluo);
	li.out.luo->del(li.out.luo);
	li.out.x->del(li.out.x);
	li.out.c->del(li.out.c);

	return 0;
}







