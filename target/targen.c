/* 
 * Argyll Color Correction System
 * Test target chart Generator.
 *
 * Author: Graeme W. Gill
 * Date:   28/9/96
 *
 * Copyright 1996 - 2004, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* This program generates a CGATS.5 compatibe file, that */
/* containing device color test patch values. */

/* TTBD:

	Should add an option to quantize the device values -
	ie. to 8 bits, so that quantizing error isn't introduced
	through a reproduction path. Should this be a default ?

	Using adaptive patch creation for grey colorspace is broken.
	This should really be disabled for grey space, or fixed.

	Adaptive only does a first go placement, itterative improvement
	is turned off, because the algorithm doesn't work.
	The useful level for the -A parameter (degree of adapatation)
	hasn't been determined.

 */

/* Description:

   >> THIS NEEDS REVISION <<

   Nearly all current Color correction systems generate test charts (or
   device characterisation target charts) by laying out a regular rectangular
   grid of test points in device space (Targen will do this if you feed it a non-zero
   m option). On some consideration, this approach is far from optimal. Not only
   is a regular grid inefficent in packing the multidimentional device space,
   but if the points are spaced evenly in device space, they will be poorly
   spaced in human perceptual space, and errors in perceptual space are
   the ultimate arbiter of the end profiles accuracy. Some commercial
   color systems tackle the latter problem by "pre-linearising" the device,
   which amounts to distorting the regular device space grid points with
   a perceptual inverse per device chanel lookup curve.

   The approach I have taken with Argyll, is a little different. By
   using an iterative sphere packing algorithm, I constrain the given
   number of test points to the devices physical gamut (including an
   ink limit for a printer device), and then try and pack the points
   evenly in human perceptual space, or even space them to minimise
   curvature approximation errors. Because the packing is a stocastic
   process, the resulting points are distributed without evident
   patterns.

#ifdef NEVER
   For higher dimensional spaces, where the aim is to create a
   more aproximate device profile, I've used a "perfect simplex
   latice" generator to lay out perfectly packed sample points
   in perceptual space. The latice spacing is sized by an
   iterative search to (hopefully) create the right number of
   test points.
#else
   For higher dimensional spaces, where the aim is to create a
   more aproximate device profile, I've used an "incremental far
   point" point generator, that for each added point, locates
   the device values that result in a percetual value farthest
   from any existing points in the test set.
#endif

   Another issue with laying test points out in regular grids, is
   that this means that the device response is poorly sampled
   (since the grids are usually coarse), and this can make it
   impossible to create detailed device linearisation "shaper"
   curves from the resulting data !
   Ideally, in any colorspace (input or output), when viewed from
   any possible angle, none of the test data points should appear
   to line up. The Argyll target generator seems to acheive this goal.

 */

#undef DEBUG

#define VRML_DIAG		/* Enable option to dump a VRML of the resulting full spread points */
#define ADDRECCLIPPOINTS	/* Add ink limited clipping points to regular grid */
#define EMPH_NEUTRAL	/* Emphasise neutral axis, like CIE94 does */
#define NEF 1.0			/* Amount to emphasis neutral axis, 0 = none, 1 = CIE94 */
#define DEFANGLE 0.3333	/* For simdlat and simplat */
#define SIMDLAT_TYPE SIMDLAT_BCC	/* Simdlat geometry type */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "vrml.h"
#include "rspl.h"
#include "cgats.h"
#include "icc.h"
#include "xicc.h"
#include "targen.h"
//#include "ppoint.h"
#include "ofps.h"
#include "ifarp.h"
#include "simplat.h"
#include "simdlat.h"
#include "prand.h"

#include <stdarg.h>

#define min2(a,b) ((a) < (b) ? (a) : (b))
#define min3(a,b,c)  (min2((a), min2((b),(c))))
#define max2(a,b) ((a) > (b) ? (a) : (b))
#define max3(a,b,c)  (max2((a), max2((b),(c))))

/* 32 bit pseudo random sequencer */
#define PSRAND32(S) (((S) & 0x80000000) ? (((S) << 1) ^ 0xa398655d) : ((S) << 1))

/* ---------------------------- */
/* The perception function data */
/* (Used for test point distribution) */
struct _pcpt {
/* public: */
	void (*del)(struct _pcpt *s);	/* We're done with it */

	int (*is_specific)(struct _pcpt *s);	/* Is a specific model, not defaulte */

	/* Conversions */
	void (*dev_to_perc)(struct _pcpt *s, double *out, double *in);	/* N-chan Perceptual */
	void (*dev_to_XYZ)(struct _pcpt *s, double *out, double *in);	/* Absolute XYZ */
	void (*dev_to_rLab)(struct _pcpt *s, double *out, double *in);	/* Relative Lab */
	void (*den_to_dev)(struct _pcpt *s, double *out, double *in);	/* Density to device */
	void (*rLab_to_dev)(struct _pcpt *s, double *out, double *in);	/* Lab to device */

/* private: */
	inkmask mask;		/* xcolorants inkmask */
	int di;				/* Number of Device dimensions */
	

	/* ICC profile based */
	icmFile *fp;
	icc   *icco;
	icmLuBase *luo;		/* Device -> rLab conversion */
	icmLuBase *luo2;	/* Device -> XYZ conversion */

	/* MPP profile based */
	mpp *mlu;			/* Device -> XYZ */

	/* Xcolorants model based */
	icxColorantLu *clu;	/* Device -> CIE */

	rspl *nlin[MXTD - 3];	/* Perceptual linearisation for other chanels */
	int e;					/* Chanel being set */

	/* Reverse lookup support */
	double ilimit;			/* Ink limit (scale 1.0) */
	double den[3];			/* Target density or Lab */
	double uniform;			/* NZ if target is uniform */	
	int kchan;				/* Set to the K chanel (-1 if none) */

}; typedef struct _pcpt pcpt;


/* Absolute XYZ conversion function */
/* Device values 0.0 - 1.0 are converted into XYZ values */
/* (Used for downstream checking) */
static void
pcpt_to_XYZ(pcpt *s, double *out, double *in) {

	if (s->luo2 != NULL)
		s->luo2->lookup(s->luo2, out, in);
	else if (s->mlu != NULL)
		s->mlu->lookup(s->mlu, out, in);
	else  if (s->clu != NULL)
		s->clu->dev_to_XYZ(s->clu, out, in);
	else {	/* Linear conversion */
		out[0] = 100.0 * in[0];
		out[1] = 100.0 * in[1] - 50.0;
		out[2] = 100.0 * in[2] - 50.0;
		icmLab2XYZ(&icmD50, out, out);
	}
}


/* Relative Lab conversion function */
/* Device values 0.0 - 1.0 are converted into Lab values */
/* (Used for VRML visualisation checking) */
static void
pcpt_to_rLab(pcpt *s, double *out, double *in) {

	if (s->luo != NULL)
		s->luo->lookup(s->luo, out, in);
	else if (s->mlu != NULL) {
		s->mlu->lookup(s->mlu, out, in);
		icmXYZ2Lab(&icmD50, out, out);
	} else if (s->clu != NULL) 
		s->clu->dev_to_rLab(s->clu, out, in);
	else {	/* Linear conversion */
		out[0] = 100.0 * in[0];
		out[1] = 100.0 * in[1] - 50.0;
		out[2] = 100.0 * in[2] - 50.0;
	}
		
}

/* Perceptual conversion function */
/* Device values 0.0 - 1.0 are converted into perceptually uniform 0.0 - 100.0 */
static void
pcpt_to_nLab(pcpt *s, double *out, double *in) {
	int e;

	/* If we have some sort of perceptual conversion */
	if (s->luo != NULL || s->mlu != NULL || s->clu != NULL) {
		double lab[3];

		if (s->luo != NULL)
			s->luo->lookup(s->luo, lab, in);
		else if (s->mlu != NULL) {
			s->mlu->lookup(s->mlu, lab, in);
			icmXYZ2Lab(&icmD50, lab, lab);
		} else 
			s->clu->dev_to_rLab(s->clu, lab, in);

#ifdef EMPH_NEUTRAL		/* Emphasise neutral axis, like CIE94 does */
		{
			double bl = 0.35 * NEF;	/* Strength of neutral axis emphasis - about 2:1 for NEF = 1 */
			double c;		/* Chromanance */

			c = sqrt(lab[1] * lab[1] + lab[2] * lab[2]);	/* Compute chromanance */

			c = 2.6624 / (1.0 + 0.013 * c);		/* Full strength scale factor */
			c = 1.0 + bl * (c - 1.0);			/* Reduced strength scale factor */

			lab[1] *= c;			/* scale a & b */
			lab[2] *= c;
		}
#endif
		/* Copy Lab values to output */
		for (e = 0; e < (s->di < 3 ? s->di : 3); e++)
			out[e] = lab[e];

		/* Lookup perceptual linearised auxiliary values */
		for (e = 0; e < (s->di-3); e++) {
			co cc;
			cc.p[0] = in[3 + e];
			s->nlin[e]->interp(s->nlin[e], &cc);
			out[3 + e] = cc.v[0];
		}

	} else {
		/* Default linear in Device space */

		for (e = 0; e < s->di; e++)
			out[e] = 100.0 * in[e];
			if (e == 1 || e == 2)
				out[e] -= 50.0;		/* Make it Lab like */
	}
}


/* Return the largest distance of the point outside the device gamut. */
/* This will be 0 if inside the gamut, and > 0 if outside.  */
static double
pcpt_in_dev_gamut(pcpt *s, double *d) {
	int e;
	int di = s->di;
	double tt, dd = 0.0;
	double ss = 0.0;

	for (e = 0; e < di; e++) {
		ss += d[e];

		tt = 0.0 - d[e];
		if (tt > 0.0) {
			if (tt > dd)
				dd = tt;
		}
		tt = d[e] - 1.0; 
		if (tt > 0.0) {
			if (tt > dd)
				dd = tt;
		}
	}
	tt = ss - s->ilimit;
	if (tt > 0.0) {
		if (tt > dd)
			dd = tt;
	}
	return dd;
}

/* optimisation function to find device values */
/* for a target density value. */
static double efunc(void *edata, double p[]) {
	pcpt *s = (pcpt *)edata;
	int e, di = s->di;
	double rv, xyz[3], den[4];

//printf("~1 efunc got dev %f %f %f %f\n",p[0],p[1],p[2],p[3]);

	pcpt_to_XYZ(s, xyz, p);		/* Convert device to XYZ */
//printf("~1 efunc got XYZ %f %f %f\n",xyz[0],xyz[1],xyz[2]);
	icx_XYZ2Tdens(den, xyz);	/* Convert XYZ to approx statusT density */
//printf("~1 efunc got density %f %f %f %f\n",den[0],den[1],den[2],den[3]);

//printf("~1 efunc got in_dev_gamut %f\n",pcpt_in_dev_gamut(s, p));

	/* Penalise for being out of gamut */
	rv = 5000.0 * pcpt_in_dev_gamut(s, p);

	/* Error to density target */
	{
		double ss = 0.0;
		for (e = 0; e < 3; e++) {
			double tt;
			tt = s->den[e] - den[e];
			ss += tt * tt;
		}
		rv += ss;
//printf("~1 efunc target den %f %f %f, err = %f, toterr %f\n",s->den[0],s->den[1],s->den[2],ss, rv);
	}

	{
		/* Minimise all channels beyond the */
		/* (assumed) first primary 3, but don't count black. */
		/* Minimise all channels except black if nchan >= 4 and uniform target */
		double ss = 0.0;

		for (e = 0; e < di; e++) {
			double tt = 0.0;

			if (di >= 4 && s->uniform && e < 3 && e != s->kchan)
				tt = p[e];			/* Minimise primary non-black if uniform */
//			else if (!s->uniform && (e < 3 || e == s->kchan))
//				tt = p[e];			/* Minimise sum of primaries & black if uniform */
			else if (e >= 3 && e != s->kchan)
				tt = 3.0 * p[e]; 	/* Supress non-primary, non-black */

			ss += tt;
		}
		rv += 1.5 * ss * ss;
//printf("~1 efunc sum err = %f, toterr %f\n",ss, rv);
	}

//printf("~1 returning %f\n",rv);
	return rv;
}

/* Given target CMY densities, return a suitable device value */ 
static void
pcpt_den_to_dev(pcpt *s, double *out, double *in) {
	int e, di = s->di;
	double tt, sr[MXTD];	/* Search radius */

//printf("\n");
//printf("~1 targen density = %f %f %f\n",in[0],in[1],in[2]);
//printf("~1 di = %d, ilimit = %f\n",s->di,s->ilimit);

	for (e = 0; e < 3; e++)
		s->den[e] = in[e];

	for (e = 0; e < di; e++) {
		sr[e] = 0.5;			/* Device space search radius */
		out[e] = 0.5;
	}

	if (fabs(in[0] - in[1]) < 0.1
	 && fabs(in[0] - in[2]) < 0.1
	 && fabs(in[1] - in[2]) < 0.1) {
		s->uniform = 1;
//printf("~1 uniform set\n");
	} else
		s->uniform = 0;

	if (powell(&tt, di, out, sr,  0.00001, 1000, efunc, (void *)s) != 0 || tt >= 50000.0) {
		error("targen: powell failed, tt = %f\n",tt);
	}

	/* Filter out silly values */
	for (e = 0; e < di; e++) {
		if (out[e] < 0.001)
			out[e] = 0.0;
		else if (out[e] > 0.999)
			out[e] = 1.0;
	}
//printf("~1 returning device values %f %f %f\n",out[0],out[1],out[2]);
}

/* Optimisation function to find device values */
/* for a target Lab value. */
static double efunc2(void *edata, double p[]) {
	pcpt *s = (pcpt *)edata;
	int e, di = s->di;
	double rv, lab[3];

//printf("~1 efunc2 got dev %f %f %f %f\n",p[0],p[1],p[2],p[3]);
//printf("~1 efunc2 got dev %f %f %f %f %f %f\n",p[0],p[1],p[2],p[3],p[4],p[5]);

	pcpt_to_rLab(s, lab, p);		/* Convert device to rLab */
//printf("~1 efunc2 got Lab %f %f %f\n",lab[0],lab[1],lab[2]);

//printf("~1 efunc2 got in_dev_gamut %f\n",pcpt_in_dev_gamut(s, p));

	/* Penalise for being out of gamut */
	rv = 5000.0 * pcpt_in_dev_gamut(s, p);

	/* Error to Lab target */
	{
		double ss = 0.0;
		for (e = 0; e < 3; e++) {
			double tt;
			tt = s->den[e] - lab[e];
			ss += tt * tt;
		}
		rv += ss;
//printf("~1 efunc2 target Lab %f %f %f, err = %f, toterr %f\n",s->den[0],s->den[1],s->den[2],ss, rv);
	}

	{
		int f;

		/* Minimise all channels except K, and especially any */
		/* beyond the first primary 3 or 4. */
		double ss = 0.0;

		if ((s->mask & ICX_CMYK) == ICX_CMYK)
			f = 4;
		else 
			f = 3;
		for (e = 0; e < di; e++) {
			if (e >= f)
				ss += 10.0 * p[e]; 	/* Supress non-primary */
			else if (e < 3)
				ss += 0.05 * p[e]; 	/* Supress first 3 primary slightly */
		}
		rv += ss * ss;
//printf("~1 efunc2 sum err = %f, toterr %f\n",ss, rv);
	}

//printf("~1 efunc2 returning %f\n\n",rv);
	return rv;
}

/* Given target Lab densities, return a suitable device value */ 
static void
pcpt_rLab_to_dev(pcpt *s, double *out, double *in) {
	int e, di = s->di;
	double tt, sr[MXTD];	/* Search radius */

//printf("\n");
//printf("#######################3\n");
//printf("~1 targen Lab = %f %f %f\n",in[0],in[1],in[2]);
//printf("~1 di = %d, ilimit = %f\n",s->di,s->ilimit);

	for (e = 0; e < 3; e++)
		s->den[e] = in[e];

	for (e = 0; e < di; e++) {
		sr[e] = 0.5;			/* Device space search radius */
		out[e] = 0.5;
	}

	if (powell(&tt, di, out, sr,  0.00001, 1000, efunc2, (void *)s) != 0 || tt >= 50000.0) {
		error("targen: powell failed, tt = %f\n",tt);
	}

	/* Filter out silly values */
	for (e = 0; e < di; e++) {
		if (out[e] <= 0.02)
			out[e] = 0.0;
		else if (out[e] >= 0.98)
			out[e] = 1.0;
	}
//printf("~1 returning device values %f %f %f\n",out[0],out[1],out[2]);
}

/* Callback to setup s->nlin[e] mapping */
static void set_nlin(void *cbntx, double *out, double *in) {
	pcpt *s = (pcpt *)cbntx;	/* Object we're setting up from */
	int e, di = s->di;
	double dev[MXTD];
	double lab[3];

	/* Just input extra channel into perceptual type lookup */
	for (e = 0; e < di; e++) 
		dev[e] = 0.0;
	dev[3 + s->e] = in[0];

	if (s->luo != NULL) {
		s->luo->lookup(s->luo, lab, dev);
	} else if (s->mlu != NULL) {
		s->mlu->lookup(s->mlu, lab, dev);
		icmXYZ2Lab(&icmD50, lab, lab);
	} else  if (s->clu != NULL) {
		s->clu->dev_to_rLab(s->clu, lab, dev);
	} else {
		lab[0] = 100.0 * in[0];
	}

	/* ~~~ should we make this delta lab along locus, rather than L value ??? */
	out[0] = lab[0];
}

/* Is a specific model, not default */
int pcpt_is_specific(pcpt *s) {
	if (s->luo2 != NULL || s->mlu != NULL)
		return 1;
	return 0;
}

/* Free the pcpt */
static void pcpt_del(pcpt *s) {

	if (s != NULL) {
		int e;

		if (s->luo != NULL) {
			s->luo->del(s->luo);
			s->luo2->del(s->luo2);
			s->icco->del(s->icco);
			s->fp->del(s->fp);
		}
		if (s->mlu != NULL) {
			s->mlu->del(s->mlu);
		}
		if (s->clu != NULL) {
			s->clu->del(s->clu);
		}
		for (e = 0; e < (s->di-3); e++) {
			if (s->nlin[e] != NULL)
				s->nlin[e]->del(s->nlin[e]);
		}
		
		free(s);
	}
}

/* Create a pcpt conversion class */
pcpt *new_pcpt(
char *profName,			/* ICC or MPP profile path, NULL for default, "none" for linear */
inkmask mask,			/* xcolorants mask */
double *ilimit			/* ink sum limit % input and return, -1 if default */
) {
	int e;
	pcpt *s;

	if ((s = (pcpt *)calloc(1, sizeof(pcpt))) == NULL) {
		fprintf(stderr,"targen: malloc failed allocating pcpt object\n");
		exit(-1);
	}
	
	/* Initialise methods */
	s->del          = pcpt_del;
	s->is_specific  = pcpt_is_specific;
	s->dev_to_perc  = pcpt_to_nLab;
	s->dev_to_XYZ   = pcpt_to_XYZ;
	s->dev_to_rLab  = pcpt_to_rLab;
	s->den_to_dev   = pcpt_den_to_dev;
	s->rLab_to_dev  = pcpt_rLab_to_dev;

	s->mask = mask;
	s->di = icx_noofinks(mask);
	/* See if we have a profile */
	if (profName != NULL
	 && profName[0] != '\000'
	 && strcmp(profName, "none") != 0
	 && strcmp(profName, "NONE") != 0) {
		int rv = 0;
		
		/* Try and open the file as an ICC profile */
		if ((s->fp = new_icmFileStd_name(profName,"r")) == NULL)
			error ("Can't open device profile '%s'",profName);
	
	
		if ((s->icco = new_icc()) == NULL)
			error ("Creation of ICC object failed");
	
		if ((rv = s->icco->read(s->icco,s->fp,0)) == 0) {
			icColorSpaceSignature ins, outs;	/* Type of input and output spaces */

			/* Get a conversion object for relative Lab */
			if ((s->luo = s->icco->get_luobj(s->icco, icmFwd, icRelativeColorimetric,
			                                 icSigLabData, icmLuOrdNorm)) == NULL) {
				if ((s->luo = s->icco->get_luobj(s->icco, icmFwd, icmDefaultIntent,
				                                 icSigLabData, icmLuOrdNorm)) == NULL) {
					error ("%d, %s",s->icco->errc, s->icco->err);
				}
			}

			/* Get a conversion object for absolute XYZ */
			if ((s->luo2 = s->icco->get_luobj(s->icco, icmFwd, icAbsoluteColorimetric,
			                                 icSigXYZData, icmLuOrdNorm)) == NULL) {
				if ((s->luo2 = s->icco->get_luobj(s->icco, icmFwd, icmDefaultIntent,
				                                 icSigXYZData, icmLuOrdNorm)) == NULL) {
					error ("%d, %s",s->icco->errc, s->icco->err);
				}
			}
		
			/* Get details of conversion (Arguments may be NULL if info not needed) */
			s->luo->spaces(s->luo, &ins, NULL, &outs, NULL, NULL, NULL, NULL, NULL);

			if (icx_colorant_comb_match_icc(mask, ins) == 0) {
				s->luo->del(s->luo);
				error("ICC profile doesn't match device!");
			}

		/* Set the default ink limits if not set by user */
		if (*ilimit < 0.0) {
			icxDefaultLimits(s->icco, ilimit, *ilimit/100.0, NULL, -1.0);
			*ilimit *= 100.0;
			if (*ilimit >= 0.0)
				*ilimit += 10.0;
		}

		} else {	/* Not a valid ICC */
			/* Close out the ICC profile */
			s->icco->del(s->icco);
			s->icco = NULL;
			s->fp->del(s->fp);
			s->fp = NULL;
		}

		/* If we don't have an ICC lookup object, look for an MPP */
		if (s->luo == NULL) {
			int imask;
			double dlimit;

			if ((s->mlu = new_mpp()) == NULL)
				error ("Creation of MPP object failed");

			if ((rv = s->mlu->read_mpp(s->mlu, profName)) != 0)
				error ("%d, %s",rv,s->mlu->err);

			s->mlu->get_info(s->mlu, &imask, NULL, &dlimit, NULL, NULL, NULL, NULL);

			if (mask != imask) {
				s->mlu->del(s->mlu);
				error("MPP profile doesn't match device!");
			}
			if (*ilimit < 0.0)	 {/* If not user specified, use MPP inklimit */
				*ilimit = 100.0 * dlimit + 10.0;
			}
		}
	}

	/* Fall back on an xcolorants model */
	if (s->luo == NULL
	 && s->mlu == NULL
	 && strcmp(profName, "none") != 0
	 && strcmp(profName, "NONE") != 0) {
		if ((s->clu = new_icxColorantLu(mask)) == NULL)
			error ("Creation of xcolorant lu object failed");
	}
	/* else leave pointers NULL */

	if (*ilimit < 0.0)
		s->ilimit = (double)s->di;	/* Default to no limit */
	else
		s->ilimit = *ilimit/100.0;

	if (s->di > 1)
		s->kchan = icx_ink2index(mask, ICX_BLACK);
	else
		s->kchan = -1;

	/* Create extra chanel linearisation lookups */
	for (e = 0; e < (s->di-3); e++) {
		double inmin = 0.0, inmax = 1.0;
		double outmax = 100.0;
		int gres = 256;

		if ((s->nlin[e] = new_rspl(1, 1)) == NULL)
			error("RSPL creation failed");

		s->e = e;	/* Chanel to set */
		s->nlin[e]->set_rspl(s->nlin[e], 0, s, set_nlin,
		                     &inmin, &inmax, &gres, &inmax, &outmax);
	}

	return s;
}

/* ------------------------------------ */

void
usage(int level) {
	int i;
	fprintf(stderr,"Generate Target deviceb test chart color values, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	fprintf(stderr,"usage: targen [options] outfile\n");
	fprintf(stderr," -v [level]      Verbose mode [optional level 1..N]\n");
	fprintf(stderr," -d col_comb     choose colorant combination from the following:\n");
	for (i = 0; ; i++) {
		char *desc; 
		if (icx_enum_colorant_comb(i, &desc) == 0)
			break;
		fprintf(stderr,"                 %d: %s\n",i,desc);
	}
	fprintf(stderr," -D colorant     Add or delete colorant from combination:\n");
	if (level == 0)
		fprintf(stderr,"                 (Use -?? to list known colorants)\n");
	else {
		fprintf(stderr,"                 %d: %s\n",0,"Additive");
		for (i = 0; ; i++) {
			char *desc; 
			if (icx_enum_colorant(i, &desc) == 0)
				break;
			fprintf(stderr,"                 %d: %s\n",i+1,desc);
		}
	}
	fprintf(stderr," -e patches      White test patches (default 4)\n");
	fprintf(stderr," -s steps        Single channel steps (default grey 50, color 0)\n");
	fprintf(stderr," -g steps        Grey axis RGB or CMY steps (default 0)\n");
	fprintf(stderr," -m steps        Multidimensional cube steps (default 2)\n");
	fprintf(stderr," -f patches      Add iterative & adaptive full spread patches to total (default grey 0, color 836)\n");
	fprintf(stderr,"                 Default is Optimised Farthest Point Sampling (OFPS)\n");
	fprintf(stderr,"  -t             Use incremental far point for full spread\n");
	fprintf(stderr,"  -r             Use device space random for full spread\n");
	fprintf(stderr,"  -R             Use perceptual space random for full spread\n");
	fprintf(stderr,"  -q             Use device space-filling quasi-random for full spread\n");
	fprintf(stderr,"  -i             Use device space body centered cubic grid for full spread\n");
	fprintf(stderr,"  -I             Use perceptual space body centered cubic grid for full spread\n");
	fprintf(stderr,"  -a angle       Simplex grid angle 0.0 - 0.5 for B.C.C. grid, default %f\n",DEFANGLE);
	fprintf(stderr,"  -A adaptation  Degree of adaptation for preconditioning in OFPS 0.0 - 1.0 (default 0.0)\n");
	fprintf(stderr," -l ilimit       Total ink limit in %%(default = none) \n");
	fprintf(stderr," -c profile      Optional device ICC or MPP pre-conditioning profile filename\n");
	fprintf(stderr,"                 (Use \"none\" to turn off any conditioning)\n");
	fprintf(stderr," -F L,a,b,rad    Filter out samples outside Lab sphere.\n");
#ifdef VRML_DIAG
	fprintf(stderr," -w              Dump diagnostic outfile.wrl file (Lab locations)\n");
	fprintf(stderr," -W              Dump diagnostic outfile.wrl file (Device locations)\n");
#endif /* VRML_DIAG */
	fprintf(stderr," outfile         Base name for output(.ti1)\n");
	exit(1);
}

/* Test if outside filter sphere. Return nz if it is */
int dofilt(
	pcpt *pdata,		/* Perceptual conversion routine */
	double *filt,		/* Filter sphere definition */
	double *dev			/* Device values to check */
) {
	int i;
	double Lab[3], rr;
	pdata->dev_to_rLab(pdata, Lab, dev);
	for (rr = 0.0, i = 0; i < 3; i++) {
		double tt = Lab[i] - filt[i];
		rr += tt * tt;
	}
	if (rr > (filt[3] * filt[3])) {
//printf("~1 rejecting rad %f of %f %f %f <=> %f %f %f\n",sqrt(rr),Lab[0],Lab[1],Lab[2],filt[0],filt[1],filt[2]);
		return 1;
	}
	return 0;
}

int main(int argc, char *argv[]) {
	int i, j, k;
	int fa, nfa, mfa;		/* current argument we're looking at */
	int verb = 0;			/* Verbose flag */
#ifdef VRML_DIAG
	int dumpvrml = 0;		/* Dump diagnostic .wrl file */
#endif /* VRML_DIAG */
	inkmask nmask = 0;		/* Ink mask combination */
	int di = 0;				/* Output dimensions */
	char *ident;			/* Ink combination identifier */
	int esteps = 4;			/* White color patches */
	int ssteps = -1;		/* Single channel steps */
	int gsteps = 0;			/* Composite grey wedge steps */
	int msteps = 2;			/* Regular grid multidimensional channel steps */
	int fsteps = -1;		/* Fitted Multidimensional patches */
	int uselat = 0;			/* Use incremental far point alg. for full spread points */
	int userand = 0;		/* Use random for full spread points */
	int useqrand = 0;		/* Use sobol for full spread points */
	int usedsim = 0;		/* Use device space simplex grid */
	int usepsim = 0;		/* Use perceptual space simplex grid */
	double simangle = DEFANGLE;	/* BCC grid angle */
	double dadapt = 0.0;	/* Degree of iterative adaptation */
	double ilimit = -1.0;	/* Ink limit (default none) */
	int filter = 0;			/* Filter values */
	double filt[4] = { 50,0,0,0 };	
	static char fname[200] = { 0 };		/* Output file base name */
	static char pname[200] = { 0 };		/* Device profile base */
	static char wname[200] = { 0 };		/* Diagnostic .wrl name */
	char buf[200];			/* Genaral use text buffer */
	int id = 1;				/* Sample ID */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	cgats *pp;				/* cgats structure */
	long stime,ttime;
	pcpt *pdata;			/* Space linearisation callback struct */
	fxpos *fxlist = NULL;	/* Fixed point list for full spread */
	int fxlist_a = 0;		/* Fixed point list allocation */
	int fxno = 0;			/* The number of fixed points */

	if (argc <= 1)
		usage(0);

#ifdef NUMSUP_H
	error_program = "targen";
#endif

	/* Process the arguments */
	mfa = 1;        /* Minimum final arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {		/* Look for any flags */
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

			if (argv[fa][1] == '?' || argv[fa][1] == '-') {
				if (argv[fa][2] == '?' || argv[fa][2] == '-')
					usage(1);
				usage(0);
			}

			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;
				if (na != NULL && na[0] >= '0' && na[0] <= '9') {
					verb = atoi(na);
					fa = nfa;
				}
			}

			/* Select the ink enumeration */
			else if (argv[fa][1] == 'd') {
				fa = nfa;
				if (na == NULL) usage(0);
				i = atoi(na);
				if (i == 0 && na[0] != '0')
					usage(0);
				if ((nmask = icx_enum_colorant_comb(i, NULL)) == 0)
					usage(0);
			}
			/* Toggle the colorant in ink combination */
			else if (argv[fa][1] == 'D') {
				int tmask;
				fa = nfa;
				if (na == NULL) usage(0);
				i = atoi(na);
				if (i == 0 && na[0] != '0')
					usage(0);
				if (i == 0)
					tmask = ICX_ADDITIVE;
				else
					if ((tmask = icx_enum_colorant(i-1, NULL)) == 0)
						usage(0);
				nmask ^= tmask;
			}
			/* White color patches */
			else if (argv[fa][1] == 'e' || argv[fa][1] == 'E') {
				int tt;
				fa = nfa;
				if (na == NULL) usage(0);
				if ((tt = atoi(na)) >= 0)
					esteps = tt;
			}
			/* Individual chanel steps */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				int tt;
				fa = nfa;
				if (na == NULL) usage(0);
				if ((tt = atoi(na)) >= 0)
					ssteps = tt;
			}
			/* RGB or CMY grey wedge steps */
			else if (argv[fa][1] == 'g' || argv[fa][1] == 'G') {
				int tt;
				fa = nfa;
				if (na == NULL) usage(0);
				if ((tt = atoi(na)) >= 0)
					gsteps = tt;
			}
			/* Multidimentional cube steps */
			else if (argv[fa][1] == 'm' || argv[fa][1] == 'M') {
				int tt;
				fa = nfa;
				if (na == NULL) usage(0);
				if ((tt = atoi(na)) >= 0) {
					msteps = tt;
					if (msteps == 1)
						msteps = 2;
				}
			}
			/* Full even spread Multidimentional patches */
			else if (argv[fa][1] == 'f') {
				int tt;
				fa = nfa;
				if (na == NULL) usage(0);
				if ((tt = atoi(na)) >= 0)
					fsteps = tt;
			}

			/* Use incremental far point algorithm for full spread */
			else if (argv[fa][1] == 't' || argv[fa][1] == 'T') {
				uselat = 1;
				userand = 0;
				useqrand = 0;
				usedsim = 0;
				usepsim = 0;
			}

			/* Random requested */
			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R') {
				uselat = 0;
				if (argv[fa][1] == 'R')
					userand = 2;
				else
					userand = 1;
				useqrand = 0;
				usedsim = 0;
				usepsim = 0;
			}

			/* Space filling quasi-random requested */
			else if (argv[fa][1] == 'q' || argv[fa][1] == 'Q') {
				uselat = 0;
				userand = 0;
				useqrand = 1;
				usedsim = 0;
				usepsim = 0;
			}


			/* Device simplex grid requested */
			else if (argv[fa][1] == 'i') {
				uselat = 0;
				userand = 0;
				useqrand = 0;
				usedsim = 1;
				usepsim = 0;
			}

			/* Perceptual simplex grid requested */
			else if (argv[fa][1] == 'I') {
				uselat = 0;
				userand = 0;
				useqrand = 0;
				usedsim = 0;
				usepsim = 1;
			}

			/* Simplex grid angle */
			else if (argv[fa][1] == 'a') {
				fa = nfa;
				if (na == NULL) usage(0);
				simangle = atof(na);
			}

			/* Degree of iterative adaptation */
			else if (argv[fa][1] == 'A') {
				fa = nfa;
				if (na == NULL) usage(0);
				dadapt = atof(na);
			}

			/* Ink limit percentage */
			else if (argv[fa][1] == 'l' || argv[fa][1] == 'L') {
				double tt;
				fa = nfa;
				if (na == NULL) usage(0);
				if ((tt = atof(na)) > 0.0)
					ilimit = tt;
				}

			/* ICC profile for perceptual linearisation */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage(0);
				strcpy (pname, na);
				}

			/* Filter out samples outside given sphere */
			else if (argv[fa][1] == 'F') {
				fa = nfa;
				if (na == NULL) usage(0);
				if (sscanf(na, " %lf,%lf,%lf,%lf ",&filt[0], &filt[1], &filt[2], &filt[3]) != 4)
					usage(0);
				filter = 1;
			}

#ifdef VRML_DIAG
			else if (argv[fa][1] == 'w')
				dumpvrml = 1;
			else if (argv[fa][1] == 'W')
				dumpvrml = 2;
#endif /* VRML_DIAG */
			else 
				usage(0);
		}
		else
			break;
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage(0);
	strcpy(fname,argv[fa]);
	strcat(fname,".ti1");

	strcpy(wname,argv[fa]);
	strcat(wname,".wrl");

	/* Set default colorant combination as CMYK */
	if (nmask == 0)
		nmask = ICX_CMYK;
	
	ident = icx_inkmask2char(nmask); 
	di = icx_noofinks(nmask);	/* Lookup number of dimensions */
	stime = clock();

	/* Implement some defaults */
	if (di == 1) {
		if (ssteps < 0)
			ssteps = 50;
		if (fsteps < 0)
			fsteps = 0;
	} else {
		if (ssteps < 0)		/* Defaults */
			ssteps = 0;
		if (fsteps < 0)
			fsteps = 836;
	}

	/* Do some sanity checking */
	if (di == 1) {
		if (ssteps == 0 && fsteps == 0 && msteps == 0)
			error ("Must have some Gray steps");
		if (gsteps > 0) {
			warning ("Composite grey steps ignored for monochrome output");
			gsteps = 0;
		}
	} else if (di == 3) {
		if (ssteps == 0 && fsteps == 0 && msteps == 0 && gsteps == 0)
			error ("Must have some single or multi dimensional RGB or CMY steps");
	} else {
		if (ssteps == 0 && fsteps == 0 && msteps == 0 && gsteps == 0)
			error ("Must have some single or multi dimensional steps");
	}

	/* Deal with ICC, MPP or fallback profile */
	if ((pdata = new_pcpt(pname, nmask, &ilimit)) == NULL) {
		error("Perceptual lookup object creation failed");
	}

	if (verb) {
		char *ss;
	
		if ((ss = icx_inkmask2char(nmask)) == NULL)
			error("Failed to find short string description for inkmask 0x%x\n",nmask);
		printf("%s test chart\n",ss);

		if (ssteps > 0)
			printf("Single channel steps = %d\n",ssteps);
		if (gsteps > 0)
			printf("Compostie Grey steps = %d\n",gsteps);
		if (fsteps > 0)
			printf("Full spread patches = %d\n",fsteps);
		if (msteps > 0)
			printf("Multi-dimention steps = %d\n",msteps);
		if (ilimit > 0.0)
			printf("Ink limit = %f%%\n",ilimit);
		if (filter) {
			printf("Filtering out samples outside sphere at %f %f %f radius %f\n",
			        filt[0], filt[1], filt[2], filt[3]);
		}
	}
	pp = new_cgats();	/* Create a CGATS structure */
	pp->add_other(pp, "CTI1"); 	/* our special type is Calibration Target Information 1 */

	pp->add_table(pp, tt_other, 0);	/* Add the first table for target points */
	pp->add_table(pp, tt_other, 0);	/* Add the second table for density pre-defined device values */
	pp->add_table(pp, tt_other, 0);	/* Add the second table for device pre-defined device values */
	pp->add_kword(pp, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 1",NULL);
	pp->add_kword(pp, 1, "DESCRIPTOR", "Argyll Calibration Target chart information 1",NULL);
	pp->add_kword(pp, 2, "DESCRIPTOR", "Argyll Calibration Target chart information 1",NULL);
	pp->add_kword(pp, 0, "ORIGINATOR", "Argyll targen", NULL);
	pp->add_kword(pp, 1, "ORIGINATOR", "Argyll targen", NULL);
	pp->add_kword(pp, 2, "ORIGINATOR", "Argyll targen", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	pp->add_kword(pp, 0, "CREATED",atm, NULL);

	/* Make available the aproximate white point to allow relative */
	/* interpretation of the aproximate XYZ values */
	{
		int e;
		double val[MXTD], XYZ[3];

		/* Setup device white */
		if (nmask & ICX_ADDITIVE)
			for (e = 0; e < di; e++)
				val[e] = 1.0;
		else
			for (e = 0; e < di; e++)
				val[e] = 0.0;
		pdata->dev_to_XYZ(pdata, XYZ, val);		/* Lookup white XYZ */

		sprintf(buf,"%f %f %f", 100.0  * XYZ[0], 100.0 * XYZ[1], 100.0 * XYZ[2]);
		pp->add_kword(pp, 0, "APPROX_WHITE_POINT", buf, NULL);
	}

	pp->add_field(pp, 0, "SAMPLE_ID", cs_t);
	pp->add_field(pp, 1, "INDEX", i_t);			/* Index no. 0..7 in second table */
	pp->add_field(pp, 2, "INDEX", i_t);			/* Index no. 0..7 in third table */

	/* Setup CGATS fields */
	{
		int j;
		char c_ilimit[20];

		for (j = 0; j < di; j++) {
			int imask;
			char fname[100];

			imask = icx_index2ink(nmask, j);
			sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : ident,
			                      icx_ink2char(imask));

			pp->add_field(pp, 0, fname, r_t);
			pp->add_field(pp, 1, fname, r_t);
			pp->add_field(pp, 2, fname, r_t);
		}

		pp->add_kword(pp, 0, "COLOR_REP", ident, NULL);

		if (ilimit > 0.0) {
			sprintf(c_ilimit,"%5.1f",ilimit);
			pp->add_kword(pp, 0, "TOTAL_INK_LIMIT", c_ilimit, NULL);
		}

	}

	/* ilimit is assumed to be in a valid range from here on */
	if (ilimit < 0.0)
		ilimit = 100.0 * di;	/* default is no limit */
	ilimit *= 0.01;		/* Convert ilimit from % to device sum */

	/* Add expected XYZ values to aid previews, scan recognition & strip recognition */
	pp->add_field(pp, 0, "XYZ_X", r_t);
	pp->add_field(pp, 0, "XYZ_Y", r_t);
	pp->add_field(pp, 0, "XYZ_Z", r_t);
	pp->add_field(pp, 1, "XYZ_X", r_t);
	pp->add_field(pp, 1, "XYZ_Y", r_t);
	pp->add_field(pp, 1, "XYZ_Z", r_t);
	pp->add_field(pp, 2, "XYZ_X", r_t);
	pp->add_field(pp, 2, "XYZ_Y", r_t);
	pp->add_field(pp, 2, "XYZ_Z", r_t);

	/* Note if the expected values are expected to be accurate */
	if (pdata->is_specific(pdata))
		pp->add_kword(pp, 0, "ACCURATE_EXPECTED_VALUES", "true", NULL);

	/* Only use optimsed full spread if <= 4 dimensions, else use ifarp */
	if (di > 4 
	 && userand == 0		/* Not other high D useful method */
	 && useqrand == 0
	 && usedsim == 0
	 && usepsim == 0)
		uselat = 1;

	/* Allocate space to record fixed steps */
	{
		fxlist_a = 4;
		if ((fxlist = (fxpos *)malloc(sizeof(fxpos) * fxlist_a)) == NULL)
			error ("Failed to malloc fxlist");
	}

	/* White color patches */
	if (esteps > 0)	{
		int j, e;

		sprintf(buf,"%d",esteps);
		pp->add_kword(pp, 0, "WHITE_COLOR_PATCHES",buf, NULL);
	
		for (j = 0; j < esteps; j++) {
			double val[MXTD], XYZ[3];
			cgats_set_elem ary[1 + MXTD + 3];
	
			if (nmask & ICX_ADDITIVE) {
				for (e = 0; e < di; e++) {
					val[e] = 1.0;			/* White is full colorant */
				}
			} else {
				for (e = 0; e < di; e++) {
					val[e] = 0.0;			/* White is no colorant */
				}
			}
	
			/* Apply general filter */
			if (filter && dofilt(pdata, filt, val))
				continue;

			sprintf(buf,"%d",id++);
			ary[0].c = buf;
			for (e = 0; e < di; e++)
				ary[1 + e].d = 100.0 * val[e];
			pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
			ary[1 + di + 0].d = 100.0 * XYZ[0];
			ary[1 + di + 1].d = 100.0 * XYZ[1];
			ary[1 + di + 2].d = 100.0 * XYZ[2];
	
			pp->add_setarr(pp, 0, ary);
	
			if (fxlist != NULL) {		/* Note in fixed list */
				if (fxno >= fxlist_a) {
					fxlist_a *= 2;
					if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
						error ("Failed to malloc fxlist");
				}
				for (e = 0; e < di; e++)
					fxlist[fxno].p[e] = val[e];
				fxno++;
			}
		}
	}

	/* Primary wedge steps */
	if (ssteps > 0)	{
		sprintf(buf,"%d",ssteps);
		pp->add_kword(pp, 0, "SINGLE_DIM_STEPS",buf, NULL);
		for (j = 0; j < di; j++) {
			for (i = 0; i < ssteps; i++) {
				int addp, e;
				double val[MXTD];
				cgats_set_elem ary[1 + MXTD + 3];

				addp = 1;			/* Default add the point */

				for (e = 0; e < di; e++) {
					if (e == j)
						val[e] = (double)i/(ssteps-1);
					else
						val[e] = 0.0;
				}

				/* See if it is already in the fixed list */
				if (fxlist != NULL) {
					int k;
					for (k = 0; k < fxno; k++) {
						for (e = 0; e < di; e++) {
							double tt;
							tt = fabs(fxlist[k].p[e] - val[e]);
							if (tt > 0.001)
								break;			/* Not identical */
						}
						if (e >= di)
							break;				/* Was identical */
					}
					if (k < fxno)				/* Found an identical patch */
						addp = 0;				/* Don't add the point */
				}

				/* Apply general filter */
				if (filter && dofilt(pdata, filt, val))
					addp = 0;

				if (addp) {
					double XYZ[3];
		
					sprintf(buf,"%d",id++);
					ary[0].c = buf;
					for (e = 0; e < di; e++)
						ary[1 + e].d = 100.0 * val[e];
					pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
					ary[1 + di + 0].d = 100.0 * XYZ[0];
					ary[1 + di + 1].d = 100.0 * XYZ[1];
					ary[1 + di + 2].d = 100.0 * XYZ[2];
	
					pp->add_setarr(pp, 0, ary);
		
					if (fxlist != NULL) {		/* Note in fixed list */
						if (fxno >= fxlist_a) {
							fxlist_a *= 2;
							if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
								error ("Failed to malloc fxlist");
						}
						for (e = 0; e < di; e++)
							fxlist[fxno].p[e] = val[e];
						fxno++;
					}
				}
			}
		}
	}

	/* Composite gray wedge steps */
	if (gsteps > 0) {
		int cix[3];		/* Composite indexes */ 
		sprintf(buf,"%d",gsteps);
		pp->add_kword(pp, 0, "COMP_GREY_STEPS",buf, NULL);


		if (nmask & ICX_ADDITIVE) {	/* Look for the RGB */
			cix[0] = icx_ink2index(nmask, ICX_RED);
			cix[1] = icx_ink2index(nmask, ICX_GREEN);
			cix[2] = icx_ink2index(nmask, ICX_BLUE);

		} else {					/* Look for the CMY */
			cix[0] = icx_ink2index(nmask, ICX_CYAN);
			cix[1] = icx_ink2index(nmask, ICX_MAGENTA);
			cix[2] = icx_ink2index(nmask, ICX_YELLOW);
		}
		if (cix[0] < 0 || cix[1] < 0 || cix[2] < 0)
			error("Composite grey wedges aren't appropriate for %s device\n",ident);

		for (i = 0; i < gsteps; i++) {
			int addp, e;
			double sum, val[MXTD];
			cgats_set_elem ary[1 + MXTD + 3];

			addp = 1;			/* Default add the point */
			sum = 0.0;

			for (e = 0; e < di; e++) {
				if (e == cix[0] || e == cix[1] || e == cix[2])
					sum += val[e] = (double)i/(gsteps-1);
				else
					val[e] = 0.0;
			}

			/* See if it is already in the fixed list */
			if (fxlist != NULL) {
				int k;
				for (k = 0; k < fxno; k++) {
					for (e = 0; e < di; e++) {
						double tt;
						tt = fabs(fxlist[k].p[e] - val[e]);
						if (tt > 0.001)
							break;			/* Not identical */
					}
					if (e >= di)
						break;				/* Was identical */
				}
				if (k < fxno)				/* Found an identical patch */
					addp = 0;				/* Don't add the point */
			}

			if (sum > ilimit)
				addp = 0;

			/* Apply general filter */
			if (filter && dofilt(pdata, filt, val))
				addp = 0;

			if (addp) {
				double XYZ[3];
	
				sprintf(buf,"%d",id++);
				ary[0].c = buf;
				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * val[e];
				pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
				ary[1 + di + 0].d = 100.0 * XYZ[0];
				ary[1 + di + 1].d = 100.0 * XYZ[1];
				ary[1 + di + 2].d = 100.0 * XYZ[2];

				pp->add_setarr(pp, 0, ary);
	
				if (fxlist != NULL) {		/* Note in fixed list */
					if (fxno >= fxlist_a) {
						fxlist_a *= 2;
						if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
							error ("Failed to malloc fxlist");
					}
					for (e = 0; e < di; e++)
						fxlist[fxno].p[e] = val[e];
					fxno++;
				}
			}
		}
	}

	/* Regular Gridded Multi dimension steps */
	if (msteps > 0) {
		int gc[MXTD];			/* Grid coordinate */

		sprintf(buf,"%d",msteps);
		pp->add_kword(pp, 0, "MULTI_DIM_STEPS",buf, NULL);

		for (j = 0; j < di; j++)
			gc[j] = 0;			/* init coords */
			
		for (;;) {	/* For all grid points */
			double val[MXTD];
			double sum;
			int addp, e;

			addp = 1;			/* Default add the point */

			for (sum = 0.0, e = 0; e < di; e++)
				sum += val[e] = (double)gc[e]/(msteps-1);
			
			if (sum > ilimit)
				addp = 0;		/* Don't add patches over ink limit */

			/* See if it is already in the fixed list */
			if (addp && fxlist != NULL) { 
				int k;
				for (k = 0; k < fxno; k++) {
					for (e = 0; e < di; e++) {
						double tt;
						tt = fabs(fxlist[k].p[e] - val[e]);
						if (tt > 0.001)
							break;			/* Not identical */
					}
					if (e >= di)
						break;				/* Was identical */
				}
				if (k < fxno)				/* Found an identical patch */
					addp = 0;				/* Don't add the point */
			}

			/* Apply general filter */
			if (filter && dofilt(pdata, filt, val))
				addp = 0;

			/* Add patch to list if OK */
			if (addp) {
				double XYZ[3];
				cgats_set_elem ary[1 + MXTD + 3];

				sprintf(buf,"%d",id++);
				ary[0].c = buf;

				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * val[e];
				pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
				ary[1 + di + 0].d = 100.0 * XYZ[0];
				ary[1 + di + 1].d = 100.0 * XYZ[1];
				ary[1 + di + 2].d = 100.0 * XYZ[2];

				pp->add_setarr(pp, 0, ary);

				if (fxlist != NULL) {		/* Note in fixed list */
					if (fxno >= fxlist_a) {
						fxlist_a *= 2;
						if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
							error ("Failed to malloc fxlist");
					}
					for (e = 0; e < di; e++)
						fxlist[fxno].p[e] = val[e];
					fxno++;
				}
			}

			/* Increment grid index and position */
			for (j = 0; j < di; j++) {
				gc[j]++;
				if (gc[j] < msteps)
					break;	/* No carry */
				gc[j] = 0;
			}
			if (j >= di)
				break;		/* Done grid */
		}

#ifdef ADDRECCLIPPOINTS
		/* Add extra points that intersect */
		/* grid, and lie on ink limit plane */
		if (ilimit < (di * 100.0)) {
			double val[MXTD], tv;
			for (k = 0; k < di; k++) {	/* dimension not on grid */
				for (j = 0; j < di; j++)
					gc[j] = 0;			/* init coords */
					
				for (;;) {	/* Until done */
					for (tv = 0.0, j = 0; j < di; j++)
						if (j != k)
							tv += val[j] = (double)gc[j]/(msteps-1);
					if (tv <= ilimit && (tv + 1.0) >= ilimit) {	/* Will intersect */
						double fr;
						val[k] = ilimit - tv; 	/* Point of intersection */
						fr = fmod((val[k] * msteps), 1.0);
						if (fr > 0.05 && fr < 0.95) {		/* Not within 5% of a grid point */ 
							int addp, e;
							cgats_set_elem ary[1 + MXTD + 3];
			
							addp = 1;			/* Default add the point */

							/* See if it is already in the fixed list */
							if (fxlist != NULL) {
								int k;
								for (k = 0; k < fxno; k++) {
									for (e = 0; e < di; e++) {
										double tt;
										tt = fabs(fxlist[k].p[e] - ary[1 + e].d);
										if (tt > 0.001)
											break;			/* Not identical */
									}
									if (e >= di)
										break;				/* Was identical */
								}
								if (k < fxno)				/* Found an identical patch */
									addp = 0;				/* Don't add the point */
							}

							/* Apply general filter */
							if (filter && dofilt(pdata, filt, val))
								addp = 0;

							if (addp) {
								double XYZ[3];
								sprintf(buf,"%d",id++);
								ary[0].c = buf;
				
								for (e = 0; e < di; e++)
									ary[1 + e].d = 100.0 * val[e];
								pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
								ary[1 + di + 0].d = 100.0 * XYZ[0];
								ary[1 + di + 1].d = 100.0 * XYZ[1];
								ary[1 + di + 2].d = 100.0 * XYZ[2];

								pp->add_setarr(pp, 0, ary);
	
								if (fxlist != NULL) {		/* Note in fixed list */
									if (fxno >= fxlist_a) {
										fxlist_a *= 2;
										if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
											error ("Failed to malloc fxlist");
									}
									for (e = 0; e < di; e++)
										fxlist[fxno].p[e] = val[e];
									fxno++;
								}
							}
						}
					}

					/* Increment grid index and position */
					for (j = 0; j < di; j++) {
						gc[j]++;
						if (j != k && gc[j] < msteps)
							break;	/* No carry */
						gc[j] = 0;
					}
					if (j >= di)
						break;		/* ALL done */
				}
			}
		}
#endif /* ADDRECCLIPPOINTS */
	}

	if (fsteps > fxno) { /* Top up with full spread (perceptually even) patches */

		if (userand == 1 || useqrand) {	/* Generate random numbers. Don't check for duplicates */
			int i, j;
			sobol *sl = NULL;

			if (useqrand) {
				if ((sl = new_sobol(di)) == NULL)
					error("Creating sobol sequence generator failed");
			}

			/* Create more points up to fsteps */
			if (verb)
				printf("\n");
			for (j = 0, i = fxno; i < fsteps;) {
				int e;
				double sum;
				double val[MXTD], XYZ[3];
				cgats_set_elem ary[1 + MXTD + 3];

				if (sl != NULL) {
					if (sl->next(sl, val))
						error("Run out of sobol random numbers!");

					for (sum = 0.0, e = 0; e < di; e++)
						sum += val[e];
				} else {	/* else uniform random distribution */
					for (sum = 0.0, e = 0; e < di; e++)
						sum += val[e] = d_rand(0.0, 1.0);
				}

				if (sum > ilimit)
					continue;

				/* Apply general filter */
				if (filter && dofilt(pdata, filt, val))
					continue;

				sprintf(buf,"%d",id++);
				ary[0].c = buf;
				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * val[e];
				pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
				ary[1 + di + 0].d = 100.0 * XYZ[0];
				ary[1 + di + 1].d = 100.0 * XYZ[1];
				ary[1 + di + 2].d = 100.0 * XYZ[2];
	
				pp->add_setarr(pp, 0, ary);

				if (fxlist != NULL) {		/* Note in fixed list to allow stats later */
					if (fxno >= fxlist_a) {
						fxlist_a *= 2;
						if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
							error ("Failed to malloc fxlist");
					}
					for (e = 0; e < di; e++)
						fxlist[fxno].p[e] = val[e];
					fxno++;
				}

				if (verb) {
					printf("\rAdded %d/%d",i+1,fxno); fflush(stdout);
				}
				i++, j++;
			}
			if (verb)
				printf("\n");

			sprintf(buf,"%d",j);
			if (sl != NULL)
				pp->add_kword(pp, 0, "SPACEFILING_RANDOM_PATCHES", buf, NULL);
			else
				pp->add_kword(pp, 0, "RANDOM_PATCHES", buf, NULL);

			if (sl != NULL)
				sl->del(sl);

		} else {
//			ppoint *s = NULL;
			ofps *s = NULL;
			ifarp *t = NULL;
			simdlat *dx = NULL;
			simplat *px = NULL;
			prand *rx = NULL;

			if (uselat)	 {
				t = new_ifarp(di, ilimit, fsteps, fxlist, fxno,
				                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);
				sprintf(buf,"%d",fsteps - fxno);
				pp->add_kword(pp, 0, "INC_FAR_PATCHES", buf, NULL);
			} else if (usedsim) {
				dx = new_simdlat(di, ilimit, fsteps, fxlist, fxno, SIMDLAT_TYPE, simangle,
			                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);
				sprintf(buf,"%d",fsteps - fxno);
				pp->add_kword(pp, 0, "SIMPLEX_DEVICE_PATCHES", buf, NULL);
			} else if (usepsim) {
				px = new_simplat(di, ilimit, fsteps, fxlist, fxno, simangle,
			                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);
				sprintf(buf,"%d",fsteps - fxno);
				pp->add_kword(pp, 0, "SIMPLEX_PERCEPTUAL_PATCHES", buf, NULL);
			} else if (userand == 2) {
				rx = new_prand(di, ilimit, fsteps, fxlist, fxno,
			                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);
				sprintf(buf,"%d",fsteps - fxno);
				pp->add_kword(pp, 0, "ERROR_OPTIMISED_PATCHES", buf, NULL);

			} else {		/* Default full spread algorithm */
// 				/* Old one */
//				s = new_ppoint(di, ilimit, fsteps, fxlist, fxno,
//			                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);

				/* New one, Optimised Farthest Point Sampling */
				s = new_ofps(di, ilimit, fsteps, dadapt, fxlist, fxno,
			                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);
				sprintf(buf,"%d",fsteps - fxno);
				pp->add_kword(pp, 0, "ERROR_OPTIMISED_PATCHES", buf, NULL);
			}

	
			for (;;) {
				int e;
				double XYZ[3], val[MXTD];
				cgats_set_elem ary[1 + MXTD + 3];
				if (( s ? s->read(s, val, NULL) :
				      t ? t->read(t, val, NULL) :
				     dx ? dx->read(dx, val, NULL) :
				     rx ? rx->read(rx, val, NULL) :
				          px->read(px, val, NULL)))
					break;
	
				/* Filter out silly values from ppoint */
				for (e = 0; e < di; e++) {
					if (val[e] < 0.001)
						val[e] = 0.0;
					else if (val[e] > 0.999)
						val[e] = 1.0;
				}
			
				/* Apply general filter */
				if (filter && dofilt(pdata, filt, val))
					continue;

				sprintf(buf,"%d",id++);
				ary[0].c = buf;
				
				for (e = 0; e < di; e++)
					ary[1 + e].d = 100.0 * val[e];
				pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */
				ary[1 + di + 0].d = 100.0 * XYZ[0];
				ary[1 + di + 1].d = 100.0 * XYZ[1];
				ary[1 + di + 2].d = 100.0 * XYZ[2];
	
				pp->add_setarr(pp, 0, ary);

				if (fxlist != NULL) {		/* Note in fixed list to allow stats later */
					if (fxno >= fxlist_a) {
						fxlist_a *= 2;
						if ((fxlist = (fxpos *)realloc(fxlist, sizeof(fxpos) * fxlist_a)) == NULL)
							error ("Failed to malloc fxlist");
					}
					for (e = 0; e < di; e++)
						fxlist[fxno].p[e] = val[e];
					fxno++;
				}
			}
			(s ? s->del(s) : t ? t->del(t) : dx ? dx->del(dx) : rx ? rx->del(rx) : px->del(px));
		}
	}

	/* Use ofps to measure the stats of the points */
	if (verb > 1
     && di <= 4
	 && (userand || useqrand || usedsim || usepsim || uselat)) {
		ofps *s;
		printf("Computing device space point stats:\n");
		s = new_ofps(di, ilimit, fxno, 0.0, fxlist, fxno,
	                (void(*)(void *, double *, double *))pdata->dev_to_perc, (void *)pdata);
		s->stats(s);
		printf("Max distance stats: Min = %f, Average = %f, Max = %f\n",s->mn,s->av,s->mx);

		s->del(s);
	}

	/* Add the eight entries in the second table. */
	/* These are legal device values that we think may */
	/* give all combinations of min/max CMY density values. */
	/* These are typically used for DTP51 and DTP41 patch separators. */
	{
		int i;

		pp->add_kword(pp, 1, "DENSITY_EXTREME_VALUES", "8", NULL);

		for (i = 0; i < 8; i++) {
			int e;
			double den[4], val[MXTD], XYZ[3];
			cgats_set_elem ary[1 + MXTD + 3];

			/* Setup target density combination */
			for (e = 0; e < 3; e++) {
				if (i & (1 << e)) 
					den[e] = 2.5;
				else
					den[e] = -0.5;
			}

			/* Lookup device values for target density */
			pdata->den_to_dev(pdata, val, den);	
			pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */

			ary[0].i = i;
			for (e = 0; e < di; e++)
				ary[1 + e].d = 100.0 * val[e];
			ary[1 + di + 0].d = 100.0 * XYZ[0];
			ary[1 + di + 1].d = 100.0 * XYZ[1];
			ary[1 + di + 2].d = 100.0 * XYZ[2];

			pp->add_setarr(pp, 1, ary);

		}
	}

	/* Add the nine entries in the third table. */
	/* These are legal device values that we calculate */
	/* give all combinations of typical CMY device values + 50% CMY */
	/* These are typically use for DTP20 bar coding. */
	{
		int i;
	
		icxColorantLu *ftarg = NULL;

		/* If not possible to use native space, use fake CMY */
		if ((nmask & ICX_CMYK) != ICX_CMYK
		 && (nmask & ICX_CMY) != ICX_CMY
		 && (nmask & ICX_RGB) != ICX_RGB) {
			if ((ftarg = new_icxColorantLu(ICX_CMY)) == NULL)
				error ("Creation of xcolorant lu object failed");
		}
		
		pp->add_kword(pp, 2, "DEVICE_COMBINATION_VALUES", "9", NULL);

		for (i = 0; i < 9; i++) {
			int e;
			double val[MXTD], lab[3], XYZ[3];
			cgats_set_elem ary[1 + MXTD + 3];

			for (e = 0; e < di; e++)
				val[e] = 0.0;

			/* Setup target device combination */
			/* Order must be White, Cyan, Magenta, Blue Yellow Green Red Black */
			if (ftarg != NULL || (nmask & ICX_CMY) == ICX_CMY) {
				for (e = 0; e < 3; e++) {
					if (i & (1 << e)) 
						val[e] = 1.0;
					else
						val[e] = 0.0;
				}
				if (i == 7)
					val[3] = 1.0;
				if (i == 8) {
					val[0] = val[1] = val[2] = 0.4;
					val[3] = 0.0;
				}

			} else {	/* RGB like */
				for (e = 0; e < 3; e++) {
					if (i & (1 << e)) 
						val[e] = 0.0;
					else
						val[e] = 1.0;
				}
				if (i == 8)
					val[0] = val[1] = val[2] = 0.6;
			}
			
			/* If target space isn't something we recognise, convert it */
			if (ftarg != NULL) {
				ftarg->dev_to_rLab(ftarg, lab, val);
				pdata->rLab_to_dev(pdata, val, lab);	

			} else if (ilimit < (double)di) {	/* Do a simple ink limit */
				double tot = 0.0;
				for (e = 0; e < di; e++)
					tot += val[e];
				if (tot > ilimit) {
					for (e = 0; e < di; e++)
						val[e] *= ilimit/tot;
				}
			}
			pdata->dev_to_XYZ(pdata, XYZ, val);		/* Add expected XYZ */

			ary[0].i = i;
			for (e = 0; e < di; e++)
				ary[1 + e].d = 100.0 * val[e];
			ary[1 + di + 0].d = 100.0 * XYZ[0];
			ary[1 + di + 1].d = 100.0 * XYZ[1];
			ary[1 + di + 2].d = 100.0 * XYZ[2];

			pp->add_setarr(pp, 2, ary);
		}

		if (ftarg != NULL)
			ftarg->del(ftarg);
	}

	ttime = clock() - stime;
	if (verb) {
		printf("Total number of patches = %d\n",id-1);
		if (id < (1 + (1 << di)))
			printf("WARNING : not enough patches for %d channels, need at least %d\n",di,(1 + (1 << di)));
		printf("Execution time = %f seconds\n",ttime/(double)CLOCKS_PER_SEC);
		// ~~99
		fprintf(stderr,"Execution time = %f seconds\n",ttime/(double)CLOCKS_PER_SEC);
	}

	if (pp->write_name(pp, fname))
		error("Write error : %s",pp->err);

#ifdef VRML_DIAG		/* Dump a VRML of the resulting points */
	if (dumpvrml != 0) {
		vrml *wrl;
		int nsets = pp->t[0].nsets;
		double rad;
		double dev[MXTD], Lab[3], col[3];

		wrl = new_vrml(wname, dumpvrml == 2 ? 0 : 1);

		/* Fudge sphere diameter */
		rad = 15.0/pow(nsets, 1.0/(double)(di <= 3 ? di : 3));

		for (i = 0; i < nsets; i++) {
			for (j = 0; j < di; j++)
				dev[j] = 0.01 * *((double *)pp->t[0].fdata[i][j + 1]);
			pdata->dev_to_rLab(pdata, Lab, dev);
			wrl->Lab2RGB(wrl, col, Lab);

			/* Fudge device values into Lab space */
			if (dumpvrml == 2) {
				Lab[0] = 100.0 * dev[0];
				Lab[1] = 100.0 * dev[1] - 50.0;
				Lab[2] = 100.0 * dev[2] - 50.0;
			}
			wrl->add_marker(wrl, Lab, col, rad);
		}
		wrl->del(wrl);		/* Write file and delete */
	}
#endif /* VRML_DIAG */

	pdata->del(pdata);	/* Cleanup perceptual conversion */

	free(ident);
	pp->del(pp);		/* Clean up */
	if (fxlist != NULL)
		free(fxlist);

	return 0;
}






