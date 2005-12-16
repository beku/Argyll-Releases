
/* 
 * cam02
 *
 * Color Appearance Model, based on
 * CIECAM02, "The CIECAM02 Color Appearance Model"
 * by Nathan Moroney, Mark D. Fairchild, Robert W.G. Hunt, Changjun Li,
 * M. Ronnier Luo and Todd Newman, IS&T/SID Tenth Color Imaging
 * Conference, with the addition of the Viewing Flare
 * model described on page 487 of "Digital Color Management",
 * by Edward Giorgianni and Thomas Madden, and the
 * Helmholtz-Kohlraush effect, using the equation
 * the Bradford-Hunt 96C model as detailed in Mark Fairchilds
 * book "Color Appearance Models". 
 *
 * Author:  Graeme W. Gill
 * Date:    17/1/2004
 * Version: 1.00
 *
 * This file is based on cam97s3.c by Graeme Gill.
 *
 * Copyright 2004 Graeme W. Gill
 * Please refer to COPYRIGHT file for details.
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the LICENCE.TXT file for licencing details.
 */


/* Note that XYZ values are normalised to 1.0 consistent */
/* with the ICC convention (not 100.0 as assumed by the CIECAM spec.) */
/* Note that all whites are assumed to be normalised (ie. Y = 1.0) */

/* Various minor changes have been made to allow the CAM conversions to */
/* function over a much greater range of XYZ and Jab values than */
/* the functions described in the above references. This is */
/* because such values arise in the process of gamut mapping, and */
/* in scanning through the grid of PCS values needed to fill in */
/* the A2B table of an ICC profile. Such values have no correlation */
/* to a real color value, but none the less need to be handled without */
/* causing an exception, in a geometrically consistent and reversible */
/* fashion. */ 

#include <copyright.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "xcam.h"
#include "cam02.h"

#undef DIAG				/* Print internal value diagnostics for each conversion */
#undef DISABLE_MATRIX	/* Debug - wire XYZ to rgba */

#define CAM02_PI 3.14159265359

/* Utility function */
/* Return a viewing condition enumeration from the given Ambient and */
/* Adapting/Surround Luminance. */
static ViewingCondition cam02_Ambient2VC(
double La,		/* Ambient Luminance (cd/m^2) */
double Lv		/* Luminance of white in the Viewing/Scene/Image field (cd/m^2) */
) {
	double r;

	if (fabs(La) < 1e-10) 		/* Hmm. */
		r = 1.0;
	else
		r = La/Lv;

	if (r < 0.01)
		return vc_dark;
	if (r < 0.2)
		return vc_dim;
	return vc_average;
}

static void cam_free(cam02 *s);
static int set_view(struct _cam02 *s, ViewingCondition Ev, double Wxyz[3],
                    double Yb, double La, double Lv, double Yf, double Fxyz[3],
					int hk);
static int XYZ_to_cam(struct _cam02 *s, double *Jab, double *xyz);
static int cam_to_XYZ(struct _cam02 *s, double *xyz, double *Jab);

/* Create a cam02 conversion object, with default viewing conditions */
cam02 *new_cam02(void) {
	cam02 *s;
	double D50[3] = { 0.9642, 1.0000, 0.8249 };

	if ((s = (cam02 *)calloc(1, sizeof(cam02))) == NULL) {
		fprintf(stderr,"cam02: malloc failed allocating object\n");
		exit(-1);
	}
	
	/* Initialise methods */
	s->del      = cam_free;
	s->set_view = set_view;
	s->XYZ_to_cam = XYZ_to_cam;
	s->cam_to_XYZ = cam_to_XYZ;

	/* Set a default viewing condition ?? */
	/* set_view(s, vc_average, D50, 0.2, 33.0, 0.0, 0.0, D50, 0); */

	return s;
}

static void cam_free(cam02 *s) {
	if (s != NULL)
		free(s);
}

/* A version of the pow() function that preserves the */
/* sign of its first argument. */
static double spow(double x, double y) {
	return x < 0.0 ? -pow(-x,y) : pow(x,y);
}

static int set_view(
cam02 *s,
ViewingCondition Ev,	/* Enumerated Viewing Condition */
double Wxyz[3],	/* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
double Yb,		/* Relative Luminance of Background to reference white */
double La,		/* Adapting/Surround Luminance cd/m^2 */
double Lv,		/* Luminance of white in the Viewing/Scene/Image field (cd/m^2) */
				/* Ignored if Ev is set to other than vc_none */
double Yf,		/* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
double Fxyz[3],	/* The Flare white coordinates (typically the Ambient color) */
int hk			/* Flag, NZ to use Helmholtz-Kohlraush effect */
) {
	double tt;

	if (Ev == vc_none)		/* Compute enumerated viewing condition */
		Ev = cam02_Ambient2VC(La, Lv);
	/* Transfer parameters to the object */
	s->Ev = Ev;
	s->Wxyz[0] = Wxyz[0];
	s->Wxyz[1] = Wxyz[1];
	s->Wxyz[2] = Wxyz[2];
	s->Yb = Yb > 0.005 ? Yb : 0.005;	/* Set minimum to avoid divide by 0.0 */
	s->La = La;
	s->Yf = Yf;
	s->Fxyz[0] = Fxyz[0];
	s->Fxyz[1] = Fxyz[1];
	s->Fxyz[2] = Fxyz[2];
	s->hk = hk;

	/* Compute the internal parameters by category */
	switch(s->Ev) {
		case vc_dark:
			s->C = 0.525;
			s->Nc = 0.8;
			s->F = 0.8;
			break;
		case vc_dim:
			s->C = 0.59;
			s->Nc = 0.95;
			s->F = 0.9;
			break;
		case vc_cut_sheet:
			s->C = 0.41;
			s->Nc = 0.8;
			s->F = 0.8;
			break;
		default:	/* average */
			s->C = 0.69;
			s->Nc = 1.0;
			s->F = 1.0;
			break;
	}

	/* Compute values that only change with viewing parameters */

	/* Figure out the Flare contribution to the flareless XYZ input */
	tt = s->Yf * s->Wxyz[1]/s->Fxyz[1];
	s->Fsxyz[0] = tt * s->Fxyz[0];
	s->Fsxyz[1] = tt * s->Fxyz[1];
	s->Fsxyz[2] = tt * s->Fxyz[2];

	/* Rescale so that the sum of the flare and the input doesn't exceed white */
	s->Fsc = s->Wxyz[1]/(s->Fsxyz[1] + s->Wxyz[1]);
	s->Fsxyz[0] *= s->Fsc;
	s->Fsxyz[1] *= s->Fsc;
	s->Fsxyz[2] *= s->Fsc;
	s->Fisc = 1.0/s->Fsc;

	/* Sharpened cone response white values */
	s->rgbW[0] =  0.7328 * s->Wxyz[0] + 0.4296 * s->Wxyz[1] - 0.1624 * s->Wxyz[2];
	s->rgbW[1] = -0.7036 * s->Wxyz[0] + 1.6975 * s->Wxyz[1] + 0.0061 * s->Wxyz[2];
	s->rgbW[2] =  0.0030 * s->Wxyz[0] + 0.0136 * s->Wxyz[1] + 0.9834 * s->Wxyz[2];

	/* Degree of chromatic adaptation */
	s->D = s->F * (1.0 - exp((-s->La - 42.0)/92.0)/3.6);

	/* Precompute Chromatic transform values */
	s->Drgb[0] = s->D * (s->Wxyz[1]/s->rgbW[0]) + 1.0 - s->D;
	s->Drgb[1] = s->D * (s->Wxyz[1]/s->rgbW[1]) + 1.0 - s->D;
	s->Drgb[2] = s->D * (s->Wxyz[1]/s->rgbW[2]) + 1.0 - s->D;

	/* Chromaticaly transformed white value */
	s->rgbcW[0] = s->Drgb[0] * s->rgbW[0];
	s->rgbcW[1] = s->Drgb[1] * s->rgbW[1];
	s->rgbcW[2] = s->Drgb[2] * s->rgbW[2];
	
	/* Transform from spectrally sharpened, to Hunt-Pointer_Estevez cone space */
	s->rgbpW[0] =  0.7409790970135310 * s->rgbcW[0]
	            +  0.2180251556757356 * s->rgbcW[1]
	            +  0.0410057473107336 * s->rgbcW[2];
	s->rgbpW[1] =  0.2853532916858800 * s->rgbcW[0]
	            +  0.6242015741188160 * s->rgbcW[1]
	            +  0.0904451341953042 * s->rgbcW[2];
	s->rgbpW[2] = -0.0096276087384294 * s->rgbcW[0]
	            -  0.0056980312161134 * s->rgbcW[1]
	            +  1.0153256399545427 * s->rgbcW[2];

	/* Background induction factor */
    s->n = s->Yb/ s->Wxyz[1];
	s->nn = pow(1.64 - pow(0.29, s->n), 0.73);	/* Pre computed value */

	/* Lightness contrast factor ?? */
	{
		double k;

		k = 1.0 / (5.0 * s->La + 1.0);
		s->Fl = 0.2 * pow(k , 4.0) * 5.0 * s->La
		      + 0.1 * pow(1.0 - pow(k , 4.0) , 2.0) * pow(5.0 * s->La , 1.0/3.0);
	}

	/* Background and Chromatic brightness induction factors */
	s->Nbb   = 0.725 * pow(1.0/s->n, 0.2);
	s->Ncb   = s->Nbb;

	/* Base exponential nonlinearity */
	s->z = 1.48 + pow(s->n , 0.5);

	/* Post-adapted cone response of white */
	tt = pow(s->Fl * s->rgbpW[0], 0.42);
	s->rgbaW[0] = (400.1 * tt + 2.713) / (tt + 27.13);
	tt = pow(s->Fl * s->rgbpW[1], 0.42);
	s->rgbaW[1] = (400.1 * tt + 2.713) / (tt + 27.13);
	tt = pow(s->Fl * s->rgbpW[2], 0.42);
	s->rgbaW[2] = (400.1 * tt + 2.713) / (tt + 27.13);

	/* Achromatic response of white */
    s->Aw = (2.0 * s->rgbaW[0] + s->rgbaW[1] + (1.0/20.0) * s->rgbaW[2] - 0.305) * s->Nbb;

#ifdef DIAG
	printf("Scene parameters:\n");
	printf("Viewing condition Ev = %d\n",s->Ev);
	printf("Ref white Wxyz = %f %f %f\n", s->Wxyz[0], s->Wxyz[1], s->Wxyz[2]);
	printf("Relative liminance of background Yb = %f\n", s->Yb);
	printf("Adapting liminance La = %f\n", s->La);
	printf("Flare Yf = %f\n", s->Yf);
	printf("Flare color Fxyz = %f %f %f\n", s->Fxyz[0], s->Fxyz[1], s->Fxyz[2]);

	printf("Internal parameters:\n");
	printf("Surround Impact C = %f\n", s->C);
	printf("Chromatic Induction Nc = %f\n", s->Nc);
	printf("Adaptation Degree F = %f\n", s->F);

	printf("Pre-computed values\n");
	printf("Sharpened cone white rgbW = %f %f %f\n", s->rgbW[0], s->rgbW[1], s->rgbW[2]);
	printf("Degree of chromatic adaptation D = %f\n", s->D);
	printf("Chromatic transform values Drgb = %f %f %f\n", s->Drgb[0], s->Drgb[1], s->Drgb[2]);
	printf("Chromatically transformed white rgbcW = %f %f %f\n", s->rgbcW[0], s->rgbcW[1], s->rgbcW[2]);
	printf("Hunter-P-E cone response white rgbpW = %f %f %f\n", s->rgbpW[0], s->rgbpW[1], s->rgbpW[2]);
	printf("Background induction factor n = %f\n", s->n);
	printf("Lightness contrast factor Fl = %f\n", s->Fl);
	printf("Background brightness induction factor Nbb = %f\n", s->Nbb);
	printf("Chromatic brightness induction factor Ncb = %f\n", s->Ncb);
	printf("Base exponential nonlinearity z = %f\n", s->z);
	printf("Post adapted cone response white rgbaW = %f %f %f\n", s->rgbaW[0], s->rgbaW[1], s->rgbaW[2]);
	printf("Achromatic response of white Aw = %f\n", s->Aw);
	printf("\n");
#endif
	return 0;
}

/* Conversions */
static int XYZ_to_cam(
struct _cam02 *s,
double Jab[3],
double XYZ[3]
) {
	int i;
	double xyz[3], rgb[3], rgbp[3], rgba[3], rgbaW[3], rgbc[3], rgbcW[3];
	double a, b, nab, J, C, h, e, A, ss;
	double ttd, ttA, ttS, tt;

//printf("\n~1 Forward conversion:\n");
//printf("~1 XYZ %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]);

#ifdef DISABLE_MATRIX

	rgba[0]= XYZ[0];
	rgba[1]= XYZ[1];
	rgba[2]= XYZ[2];

#else /* !DISABLE_MATRIX */

	/* Add in flare */
	xyz[0] = s->Fsc * XYZ[0] + s->Fsxyz[0];
	xyz[1] = s->Fsc * XYZ[1] + s->Fsxyz[1];
	xyz[2] = s->Fsc * XYZ[2] + s->Fsxyz[2];

	/* Spectrally sharpened cone responses */
	rgb[0] =  0.7328 * xyz[0] + 0.4296 * xyz[1] - 0.1624 * xyz[2];
	rgb[1] = -0.7036 * xyz[0] + 1.6975 * xyz[1] + 0.0061 * xyz[2];
	rgb[2] =  0.0030 * xyz[0] + 0.0136 * xyz[1] + 0.9834 * xyz[2];
	
	/* Chromaticaly transformed sample value */
	rgbc[0] = s->Drgb[0] * rgb[0];
	rgbc[1] = s->Drgb[1] * rgb[1];
	rgbc[2] = s->Drgb[2] * rgb[2];
	
	/* Transform from spectrally sharpened, to Hunt-Pointer_Estevez cone space */
	rgbp[0] =  0.7409790970135310 * rgbc[0]
	        +  0.2180251556757356 * rgbc[1]
	        +  0.0410057473107336 * rgbc[2];
	rgbp[1] =  0.2853532916858800 * rgbc[0]
	        +  0.6242015741188160 * rgbc[1]
	        +  0.0904451341953042 * rgbc[2];
	rgbp[2] = -0.0096276087384294 * rgbc[0]
	        -  0.0056980312161134 * rgbc[1]
	        +  1.0153256399545427 * rgbc[2];

	/* Post-adapted cone response of sample. */
	/* rgba[] has a minimum value of 0.1 for XYZ[] = 0 and no flare. */
	/* We add a symetric negative compression region, plus linear segments at */
	/* the ends of this conversion to allow numerical handling of a */
	/* very wide range of values. */
	for (i = 0; i < 3; i++) {
		if (rgbp[i] < 0.0) {
			tt = pow(s->Fl * -rgbp[i], 0.42);
			if (tt < 10756.69231)
				rgba[i] = (2.713 - 397.387 * tt) / (tt + 27.13);
			else 
				rgba[i] = (2.713 - tt) / 27.13;

		} else {
			tt = pow(s->Fl * rgbp[i], 0.42);
			if (tt < 10824.87)
				rgba[i] = (400.1 * tt + 2.713) / (tt + 27.13);
			else
				rgba[i] = (tt + 2.713) / 27.13;
		}
	}

#endif /* !DISABLE_MATRIX */

	/* Deal with the core rgb to A, S & h conversion: */

//printf("~1 rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);

	/* Note that the minimum values of rgba[] for XYZ = 0 is 0.1, */
	/* hence magic 0.305 below comes from the following weighting of rgba[], */
	/* to base A at 0.0 */
	ttA = 2.0 * rgba[0] + rgba[1] + (1.0/20.0) * rgba[2];

	/* Avoid the inverse blowing up due to ttA == 0 */
	/* Note that the three threshold values (Fwd ttA, Fwd ttd and Bwd ttA) */
	/* have bee carefuly juggled to minimise errors in typical situations. */
	if (fabs(ttA) < 0.00002) {
		double shift = ttA < 0.0 ? -0.00002 : 0.00002;
//printf("~1 ttA %f being limited to %f\n",ttA,shift);
		rgba[0] += (shift-ttA) * (20.0 * 40.0)/(40.0 * 40.0 + 20.0 * 20.0 + 1.0);
		rgba[1] += (shift-ttA) * (20.0 * 20.0)/(40.0 * 40.0 + 20.0 * 20.0 + 1.0);
		rgba[2] += (shift-ttA) * (20.0 * 1.0)/(40.0 * 40.0 + 20.0 * 20.0 + 1.0);
		ttA = 2.0 * rgba[0] + rgba[1] + (1.0/20.0) * rgba[2];
//printf("~1 After adjustment, ttA = %e\n",ttA);
	}

	/* Preliminary Saturation denominator */
	ttd = rgba[0] + rgba[1] + (21.0/20.0) * rgba[2];
	
	/* Avoid Preliminary Saturation denominator blowing things up */
	if (fabs(ttd) < 0.00000001) {
		double shift = ttd < 0.0 ? -0.00000001 : 0.00000001;
//printf("~1 ttd %f being limited to %f\n",ttd,shift);
		rgba[0] += shift * 20.0/61;
		rgba[1] += shift * 20.0/61;
		rgba[2] += shift * 21.0/40;
		ttd = rgba[0] + rgba[1] + (21.0/20.0) * rgba[2];
		ttA = 2.0 * rgba[0] + rgba[1] + (1.0/20.0) * rgba[2];
//printf("~1 After adjustment, ttA = %e, ttd = %e\n",ttA,ttd);
	}

	/* Preliminary red-green & yellow-blue opponent dimensions */
	a     = rgba[0] - 12.0 * rgba[1]/11.0 + rgba[2]/11.0;
    b     = (1.0/9.0) * (rgba[0] + rgba[1] - 2.0 * rgba[2]);
	nab   = sqrt(a * a + b * b);		/* Normalised a, b */

	/* Deal with sign flip in ttd */
	if (ttd < 0.0) {
//printf("~1 Experimental ttd sign flip triggered\n");
		ttd = -ttd;
		a = -a;
		b = -b;
//printf("~1 After sign flip: a = %f, b = %f, ttd = %f\n",a,b,ttd);
	}

	/* Preliminary Saturation - always +ve */
	ttS = nab/ttd;

	/* Final hue angle */
    h = (180.0/CAM02_PI) * atan2(b,a);
	h = (h < 0.0) ? h + 360.0 : h;

//printf("~1 ttd = %f, ttA = %f, ttS = %f, h = %f\n",ttd,ttA,ttS,h);

	/* Finish of conversion from core values to Jab: */

	/* Eccentricity factor */
	e = (cos(h * CAM02_PI/180.0 + 2.0) + 3.8);

	/* Achromatic response */
	A = (ttA - 0.305) * s->Nbb;

	/* Lightness */
	/* Derived directly from Acromatic response. */
	J = spow(A/s->Aw, s->C * s->z);		/* J/100  - keep Sign */

	/* Preliminary saturation always +ve */
	ss = (12500.0/13.0 * s->Nc * s->Ncb * e) * ttS;

	/* Croma - always +ve */
	tt = fabs(J);
	C = pow(ss, 0.9) * spow(tt, 0.5) * s->nn;
	
//printf("~1 ss = %f, A = %f, J = %f, C = %f\n",ss,A,J,C);

 	/* Helmholtz-Kohlraush effect */
	if (s->hk) {
		double kk = C/300.0 * sin(CAM02_PI * fabs(0.5 * (h - 90.0))/180.0);
		if (kk > 0.9)		/* Limit kk to a reasonable range */
			kk = 0.9;
		J = J + (1.0 - J) * kk;
	}

	J *= 100.0;		/* Scale J */

	/* Compute Jab value */
	Jab[0] = J;
	if (nab > 1e-10) {
		Jab[1] = C * a/nab;
		Jab[2] = C * b/nab;
	} else {
		Jab[1] = 0.0;
		Jab[2] = 0.0;
	}

//printf("~1 Jab %f %f %f\n",Jab[0], Jab[1], Jab[2]);

#ifdef DIAG
	printf("Processing:\n");
	printf("XYZ = %f %f %f\n", XYZ[0], XYZ[1], XYZ[2]);
	printf("Including flare XYZ = %f %f %f\n", xyz[0], xyz[1], xyz[2]);
	printf("Sharpened cone sample rgb = %f %f %f\n", rgb[0], rgb[1], rgb[2]);
	printf("Chromatically transformed sample value rgbc = %f %f %f\n", rgbc[0], rgbc[1], rgbc[2]);
	printf("Hunt-P-E cone space rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]);
	printf("Post adapted cone response rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);
	printf("Prelim red green a = %f, b = %f\n", a, b);
	printf("Hue angle h = %f\n", h);
	printf("Eccentricity factor e = %f\n", e);
	printf("Achromatic response A = %f\n", A);
	printf("Lightness J = %f\n", J);
	printf("Prelim Saturation ss = %f\n", ss);
	printf("Chroma C = %f\n", C);
	printf("Jab = %f %f %f\n", Jab[0], Jab[1], Jab[2]);
	printf("\n");
#endif
	return 0;
}

static int cam_to_XYZ(
struct _cam02 *s,
double XYZ[3],
double Jab[3]
) {
	int i;
	double xyz[3], rgb[3], rgbp[3], rgba[3], rgbaW[3], rgbc[3], rgbcW[3];
	double ja, jb, aa, ab, a, b, J, C, h, e, A, ss;
	double tt, ttA, ttS;

//printf("\n~1 Reverse conversion\n");
//printf("~1 Jab %f %f %f\n",Jab[0], Jab[1], Jab[2]);

	J = Jab[0] * 0.01;	/* J/100 */
	ja = Jab[1];
	jb = Jab[2];

	/* Convert Jab to core A, S & h values: */

	/* Compute hue angle */
    h = (180.0/CAM02_PI) * atan2(jb, ja);
	h = (h < 0.0) ? h + 360.0 : h;
	
	/* Compute chroma value */
	C = sqrt(ja * ja + jb * jb);		/* Must be Always +ve */

 	/* Helmholtz-Kohlraush effect */
	if (s->hk) {
		double kk = C/300.0 * sin(CAM02_PI * fabs(0.5 * (h - 90.0))/180.0);
		if (kk > 0.9)		/* Limit kk to a reasonable range */
			kk = 0.9;
		J = (J - kk)/(1.0 - kk);
	}

	/* Eccentricity factor */
	e = (cos(h * CAM02_PI/180.0 + 2.0) + 3.8);

	/* Achromatic response */
	A = spow(J, 1.0/(s->C * s->z)) * s->Aw;

	/* Preliminary Saturation - always +ve */
	tt = fabs(J);
	ss = pow(C/(pow(tt, 0.5) * s->nn), 1.0/0.9);	/* keep +ve */

//printf("~1 ss = %f, A = %f, J = %f, C = %f\n",ss,A,J,C);

	ttS = ss/(12500.0/13.0 * e * s->Nc * s->Ncb);	/* Simplified variable, +ve */
    ttA = (A/s->Nbb)+0.305;							/* Simplified variable, +ve */

//printf("~1 ttA = %f, ttS = %f, h = %f\n",ttA,ttS,h);

	/* Invert A, S & h into rgb values */

	/* Attempt to avoid loosing saturation/hue information */
	if (fabs(ttA) < 0.00002) {
//printf("~1 ttA being limited to +-0.00002\n");
		if (ttA < 0.0)
			ttA = -0.00002;
		else
			ttA = 0.00002;
	}
	
	/* Solve for preliminary a & b, taking care of numerical problems */
	aa = fabs(ja);
	ab = fabs(jb);

	if (aa < 1e-10 && ab < 1e-10) {
		double sign = ttA < 0.0 ? -1.0 : 1.0;
		a = sign * 1e-8 * cos(CAM02_PI * h/180.0);
		b = sign * 1e-8 * sin(CAM02_PI * h/180.0);
//printf("~1 ttS nearly zero, preserve hue angle\n");

	} else if (aa > ab) {
		double ttanh = jb/ja;
		double sign = (h > 90.0 && h <= 270.0) ? -1.0 : 1.0;
		double tsden;		/* Solution denominator */

		tsden = (108.0/23.0 * ttanh + 11.0/23.0) * ttS + sign * sqrt(ttanh * ttanh + 1.0);
//printf("~1 aa > ab, ttanh = %f, signA = %f, tsden = %f\n",ttanh, sign, tsden);

#ifdef NEVER
		if (fabs(tsden) < 1e-10)  {		/* Happens if ttA == 0 */
			printf("~1 #### Warning tsden is %e, composed of %f + %f\n",tsden,(108.0/23.0 * ttanh + 11.0/23.0) * ttS,sign * sqrt(ttanh * ttanh + 1.0));
		}
#endif

		a = (ttA * ttS) / tsden;
		b = a * ttanh;

	} else {	/* ab > aa */
		double itanh = ja/jb;
		double sign = (h > 180.0 && h <= 360.0) ? -1.0 : 1.0;
		double isden;		/* Solution denominator */

		isden = (108.0/23.0 + 11.0/23.0 * itanh) * ttS + sign * sqrt(itanh * itanh + 1.0);
//printf("~1 ab > aa, itanh = %f, signB = %f, isden = %f\n",itanh,sign, isden);

#ifdef NEVER
		if (fabs(isden) < 1e-10)		/* Happens if ttA == 0 */
			printf("~1 #### Warning isden is %e, composed of %f + %f\n",isden,(108.0/23.0 + 11.0/23.0 * itanh) * ttS,sign * sqrt(itanh * itanh + 1.0));
#endif

		b = (ttA * ttS) / isden;
		a = b * itanh;
	}

	/* Post-adapted cone response of sample */
    rgba[0] = (20.0/61.0) * ttA
	        + ((41.0 * 11.0)/(61.0 * 23.0)) * a
	        + ((288.0 * 1.0)/(61.0 * 23.0)) * b;
    rgba[1] = (20.0/61.0) * ttA
	        - ((81.0 * 11.0)/(61.0 * 23.0)) * a
	        - ((261.0 * 1.0)/(61.0 * 23.0)) * b;
    rgba[2] = (20.0/61.0) * ttA
	        - ((20.0 * 11.0)/(61.0 * 23.0)) * a
	        - ((20.0 * 315.0)/(61.0 * 23.0)) * b;

//printf("~1 rgba %f %f %f\n",rgba[0], rgba[1], rgba[2]);

#ifdef DISABLE_MATRIX

	XYZ[0] = rgba[0];
	XYZ[1] = rgba[1];
	XYZ[2] = rgba[2];

#else /* !DISABLE_MATRIX */

	/* Hunt-Pointer_Estevez cone space */
	/* (with linear segments at the ends) */
	tt = 1.0/s->Fl;
	for (i = 0; i < 3; i++) {
		if (rgba[i] < 0.1) {
			double ta = rgba[i] > -396.387 ? rgba[i] : -396.387;
			rgbp[i] = -tt * pow((2.713 - 27.13 * rgba[i] )/(397.387 + ta), 1.0/0.42);
		} else {
			double ta = rgba[i] < 399.1 ? rgba[i] : 399.1;
			rgbp[i] =  tt * pow((27.13 * rgba[i] -2.713)/(400.1 - ta), 1.0/0.42);
		}
	}

	/* Chromaticaly transformed sample value */
	rgbc[0] =  1.5591523979049677 * rgbp[0]
	        -  0.5447226796590880 * rgbp[1]
	        -  0.0144453097698588 * rgbp[2];
	rgbc[1] = -0.7143267176368630 * rgbp[0]
	        +  1.8503099728895096 * rgbp[1]
	        -  0.1359761119854705 * rgbp[2];
	rgbc[2] =  0.0107755117023383 * rgbp[0]
	        +  0.0052187662221759 * rgbp[1]
	        +  0.9840056143203700 * rgbp[2];

	/* Spectrally sharpened cone responses */
	rgb[0]  =  rgbc[0]/s->Drgb[0];
	rgb[1]  =  rgbc[1]/s->Drgb[1];
	rgb[2]  =  rgbc[2]/s->Drgb[2];

	/* XYZ values */
	xyz[0] =  1.0961238208355140 * rgb[0]
	       -  0.2788690002182872 * rgb[1]
	       +  0.1827451793827730 * rgb[2];
	xyz[1] =  0.4543690419753590 * rgb[0]
	       +  0.4735331543074120 * rgb[1]
	       +  0.0720978037172291 * rgb[2];
	xyz[2] = -0.0096276087384294 * rgb[0]
	       -  0.0056980312161134 * rgb[1]
	       +  1.0153256399545427 * rgb[2];

	/* Subtract flare */
	XYZ[0] = s->Fisc * (xyz[0] - s->Fsxyz[0]);
	XYZ[1] = s->Fisc * (xyz[1] - s->Fsxyz[1]);
	XYZ[2] = s->Fisc * (xyz[2] - s->Fsxyz[2]);

#endif /* !DISABLE_MATRIX */

//printf("~1 XYZ = %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]);

#ifdef DIAG
	printf("Processing:\n");
	printf("Jab = %f %f %f\n", Jab[0], Jab[1], Jab[2]);
	printf("Chroma C = %f\n", C);
	printf("Preliminary Saturation ss = %f\n", ss);
	printf("Lightness J = %f\n", J * 100.0);
	printf("Achromatic response A = %f\n", A);
	printf("Eccentricity factor e = %f\n", e);
	printf("Hue angle h = %f\n", h);
	printf("Prelim red green a = %f, b = %f\n", a, b);
	printf("Post adapted cone response rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);
	printf("Hunt-P-E cone space rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]);
	printf("Chromatically transformed sample value rgbc = %f %f %f\n", rgbc[0], rgbc[1], rgbc[2]);
	printf("Sharpened cone sample rgb = %f %f %f\n", rgb[0], rgb[1], rgb[2]);
	printf("Including flare XYZ = %f %f %f\n", xyz[0], xyz[1], xyz[2]);
	printf("XYZ = %f %f %f\n", XYZ[0], XYZ[1], XYZ[2]);
	printf("\n");
#endif
	return 0;
}




