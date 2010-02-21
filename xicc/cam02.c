
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
 * Helmholtz-Kohlraush effect, using the equation from
 * the Bradford-Hunt 96C model as detailed in Mark Fairchild's
 * book "Color Appearance Models". 
 * The Slight modification to the Hunt-Pointer-Estevez matrix
 * recommended by Kuo, Zeis and Lai in IS&T/SID 14th Color Imaging
 * Conference "Robust CIECAM02 Implementation and Numerical
 * Experiment within an ICC Workflow", page 215, together with
 * their matrix formulation of inversion has been adopted.
 *
 * In addition the algorithm has unique extensions to allow
 * it to operate reasonably over an unbounded domain.
 *
 * Author:  Graeme W. Gill
 * Date:    17/1/2004
 * Version: 1.00
 *
 * This file is based on cam97s3.c by Graeme Gill.
 *
 * Copyright 2004 - 2009 Graeme W. Gill
 * Please refer to COPYRIGHT file for details.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */


/* Note that XYZ values are normalised to 1.0 consistent */
/* with the ICC convention (not 100.0 as assumed by the CIECAM spec.) */
/* Note that all whites are assumed to be normalised (ie. Y = 1.0) */

/*
	Various additions and changes have been made to allow the CAM
	conversions to and from an unbounded range of XYZ and Jab values,
	in a (somewhat) geometrically consistent maner. This is because
	such values arise in the process of gamut mapping, and in scanning through
	the grid of PCS values needed to fill in the A2B table of an
	ICC profile. Such values have no correlation to a real color
	value, but none the less need to be handled without causing an
	exception, in a geometrically consistent and reversible
	fashion.

	The changes are:

	The Hunt-Pointer-Estevez matrix has been increased in precission by
	1 digit [Kuo, Zeise & Lai, IS&T/SID 14th Color Imaging Conference pp. 215-219]

	The sharpened cone response matrix third row has been changed from
	[0.0030, 0.0136, 0.9834] to [0.0000, 0.0000, 1.0000]
	[Susstrunk & Brill, IS&T/SID 14th Color Imaging Conference Late Breaking News sublement]

	The post-adaptation non-linearity response has had a straight
	line negative extension added, that preserves negative values
	but expands them more than the positive characteristic on conversion
	to the CAM, and therefore compresses such values more on
	coversion from the CAM back to XYZ. This helps minimise 
	shifts in resulting color for CAM values near or just below zero.
	The positive region has a tangential linear extension added, so
	that the range of values is not limited.
		
	Re-arrange ss equation into separated effects,
	k1, k2, k3, and then limit these to avoid divide by zero
	and sign change errors in fwd and bwd converson.

	To avoid chroma and hue angle information being lost when the
	J value used to scale the chroma is 0, and to ensure
	that J = 0, a,b != 0 values have a valid mapping into
	XYZ space, the J value used to multiply Chroma, is limited
	to be equivalent to not less than A == 0.1.

	The Helmholtz-Kohlraush effect is crafted to have resonable
	effects over the range of J from 0 to 100, and have more
	moderate effects outside this range. 
*/

#include <copyright.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "xcam.h"
#include "cam02.h"
#include "numlib.h"

#define ENABLE_DDL			/* Enable k1,k2,k3 overall ss limit values (seems to be the best scheme) */
#undef ENABLE_SS			/* Disable overall ss limit values (not the scheme used) */

#undef ENTRACE				/* Enable internal value runtime tracing if s->trace != 0 */
#undef DIAG1				/* Print internal value diagnostics for conditions setup */
#undef DIAG2				/* Print internal value diagnostics for each conversion */
#undef TRACKMINMAX			/* Track min/max DD & SS limit values (run with locus cam02test) */
#undef DISABLE_MATRIX		/* Debug - wire XYZ to rgba */
#undef DISABLE_NONLIN		/* Debug - wire rgbp to rgba */
#undef DISABLE_TTD			/* Debug - disable ttd vector 'tilt' */
#undef DISABLE_HHKR			/* Debug - disable Helmholtz-Kohlraush */
#undef DISABLE_BLUECLIP		/* Debug - disable blue clipping  */

#define BLUECLIMIT 1.3		/* Blue clipping max ratio of z to x+y */

#ifdef NEVER
#define NLDLIMIT 0.005		/* Non-linearity lower crossover to straight line */
#define NLDLSLOPE 400.0		/* Lower slope if NLDLIMIT == 0.0 */
#define NLDLPOW 2.0			/* Lower power. < 1 improves black->blue problem */

#else	/* Another approach to improve -ve black goes to blue problem ?? */

#define NLDLIMIT 0.0		/* Non-linearity lower crossover to straight line */
#define NLDLSLOPE 0.01		/* Lower slope if NLDLIMIT == 0.0 */
							/* (Smaller improves black->blue problem) */
#define NLDLPOW 2.0			/* Lower power. < 1 improves black->blue problem ? */

#endif	/* NEVER */

#define NLULIMIT 1e5		/* Non-linearity upper crossover to straight line */

#ifdef ENABLE_SS
#define SSLLIMIT 0.22		/* Overall ab scale lower limit */
#define SSULIMIT 580.0      /* Overall ab scale upper limit */
#endif /* ENABLE_SS */

#define DDLLIMIT 0.55		/* ab component -k3:k2 ratio limit (must be < 1.0) */
// NEED 0.9993 to match reference CIECAM02 behaviour within the spectrum locus. */
// The value 0.60 limits the blue turning black problem though. */
//#define DDULIMIT 0.9993		/* ab component k3:k1 ratio limit (must be < 1.0) */
#define DDULIMIT 0.90		/* ab component k3:k1 ratio limit (must be < 1.0) */
#define SSMINcJ 0.005		/* ab scale cJ minimum value */
#define JLIMIT 0.005		/* J encoding cutover point straight line (0 - 1.0 range) */
#define HKLIMIT 0.5			/* Maximum Helmholtz-Kohlraush lift */

#ifdef TRACKMINMAX
double minss = 1e60;
double maxss = -1e60;
double minlrat = 0.0;
double maxurat = 0.0;
#define noslots 103
double slotsd[noslots];
double slotsu[noslots];
double minj = 1e38, maxj = -1e38;
#endif /* TRACKMINMAX */

#ifdef ENTRACE
#define TRACE(xxxx) if (s->trace) printf xxxx ;
#else
#define TRACE(xxxx) 
#endif

static void cam_free(cam02 *s);
static int set_view(struct _cam02 *s, ViewingCondition Ev, double Wxyz[3],
	                double La, double Yb, double Lv, double Yf, double Fxyz[3],
					int hk, int noclip);
static int XYZ_to_cam(struct _cam02 *s, double *Jab, double *xyz);
static int cam_to_XYZ(struct _cam02 *s, double *xyz, double *Jab);

static double spow(double val, double pp) {
	if (val < 0.0)
		return -pow(-val, pp);
	else
		return pow(val, pp);
}

/* Create a cam02 conversion object, with default viewing conditions */
cam02 *new_cam02(void) {
	cam02 *s;
//	double D50[3] = { 0.9642, 1.0000, 0.8249 };

	if ((s = (cam02 *)calloc(1, sizeof(cam02))) == NULL) {
		fprintf(stderr,"cam02: malloc failed allocating object\n");
		exit(-1);
	}
	
	/* Initialise methods */
	s->del      = cam_free;
	s->set_view = set_view;
	s->XYZ_to_cam = XYZ_to_cam;
	s->cam_to_XYZ = cam_to_XYZ;

	/* Set default range handling limits */
	s->nldlimit = NLDLIMIT;
	s->nldxslope = NLDLSLOPE;
	s->nldpow = NLDLPOW;
	s->nlulimit = NLULIMIT;
	s->ddllimit = DDLLIMIT;
	s->ddulimit = DDULIMIT;
	s->ssmincj = SSMINcJ;
	s->jlimit = JLIMIT;
	s->hklimit = 1.0 / HKLIMIT;

	/* Set a default viewing condition ?? */
	/* set_view(s, vc_average, D50, 33, 0.2, 0.0, 0.0, D50, 0); */

#ifdef TRACKMINMAX
	{
		int i;
		for (i = 0; i < noslots; i++) {
			slotsd[i] = 0.0;
			slotsu[i] = 0.0;
		}
	}
#endif /* TRACKMINMAX */

	return s;
}

static void cam_free(cam02 *s) {

#ifdef TRACKMINMAX
	{
		int i;
		for (i = 0; i < noslots; i++) {
			printf("Slot %d = %f, %f\n",i,slotsd[i], slotsu[i]);
		}
		printf("minj = %f, maxj = %f\n",minj,maxj);
	
		printf("minss = %f\n",minss);
		printf("maxss = %f\n",maxss);
		printf("minlrat = %f\n",minlrat);
		printf("maxurat = %f\n",maxurat);
	}
#endif /* TRACKMINMAX */

	if (s != NULL)
		free(s);
}

/* Return value is always 0 */
static int set_view(
cam02 *s,
ViewingCondition Ev,	/* Enumerated Viewing Condition */
double Wxyz[3],	/* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
double La,		/* Adapting/Surround Luminance cd/m^2 */
double Yb,		/* Relative Luminance of Background to reference white (range 0.0 .. 1.0) */
double Lv,		/* Luminance of white in the Viewing/Scene/Image field (cd/m^2) */
				/* Ignored if Ev is set to other than vc_none */
double Yf,		/* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
double Fxyz[3],	/* The Flare white coordinates (typically the Ambient color) */
int hk,			/* Flag, NZ to use Helmholtz-Kohlraush effect */
int noclip		/* Flag, NZ to not clip to useful gamut before XYZ_to_cam() */
) {
	double tt, t1, t2;

	if (Ev == vc_none) {
		/* Compute the internal parameters by interpolation */
		int i;
		double r, bf;
		/* Dark, dim, average, above average */
		double t_C[4]  = { 0.525, 0.59, 0.69, 1.0 }; 
		double t_Nc[4] = { 0.800, 0.95, 1.00, 1.0 }; 
		double t_F[4]  = { 0.800, 0.90, 1.00, 1.0 }; 

		if (La < 1e-10) 		/* Hmm. */
			La = 1e-10;
		r = La/Lv;
		if (r < 0.0)
			r = 0.0;
		else if (r > 1.0)
			r = 1.0;
	
		if (r < 0.1) {			/* Dark to Dim */
			i = 0;
			bf = r/0.1;
		} else if (r < 0.2) {	/* Dim to Average */
			i = 1;
			bf = (r - 0.1)/0.1;
		} else {				/* Average to above average */
			i = 2;
			bf = (r - 0.2)/0.8;
		}
		s->C  = t_C[i] * (1.0 - bf)  + t_C[i+1] * bf;
		s->Nc = t_Nc[i] * (1.0 - bf) + t_Nc[i+1] * bf;
		s->F  = t_F[i] * (1.0 - bf)  + t_F[i+1] * bf;
	} else {
		/* Compute the internal parameters by category */
		switch(Ev) {
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
	}

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
	s->noclip = noclip;

	/* The rgba vectors */
	s->Va[0] = 1.0;
	s->Va[1] = -12.0/11.0;
	s->Va[2] = 1.0/11.0;

	s->Vb[0] = 1.0/9.0;
	s->Vb[1] = 1.0/9.0;
	s->Vb[2] = -2.0/9.0;

	s->VttA[0] = 2.0;
	s->VttA[1] = 1.0;
	s->VttA[2] = 1.0/20.0;

	/* (Not used) */
	s->Vttd[0] = 1.0;
	s->Vttd[1] = 1.0;
	s->Vttd[2] = 21.0/20.0;

	/* Vttd in terms of the VttA, Va and Vb vectors */
	s->dcomp[0] = 1.0;
	s->dcomp[1] = -11.0/23.0;
	s->dcomp[2] = -108.0/23.0;

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
	s->rgbW[2] =  0.0000 * s->Wxyz[0] + 0.0000 * s->Wxyz[1] + 1.0000 * s->Wxyz[2];

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
	s->rgbpW[0] =  0.7409744840453773 * s->rgbcW[0]
	            +  0.2180245944753982 * s->rgbcW[1]
	            +  0.0410009214792244 * s->rgbcW[2];
	s->rgbpW[1] =  0.2853532916858801 * s->rgbcW[0]
	            +  0.6242015741188157 * s->rgbcW[1]
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
	s->rgbaW[0] = 400.0 * tt / (tt + 27.13) + 0.1;
	tt = pow(s->Fl * s->rgbpW[1], 0.42);
	s->rgbaW[1] = 400.0 * tt / (tt + 27.13) + 0.1;
	tt = pow(s->Fl * s->rgbpW[2], 0.42);
	s->rgbaW[2] = 400.0 * tt / (tt + 27.13) + 0.1;

	/* Achromatic response of white */
	s->Aw = (s->VttA[0] * s->rgbaW[0]
	      +  s->VttA[1] * s->rgbaW[1]
	      +  s->VttA[2] * s->rgbaW[2] - 0.305) * s->Nbb;

	/* Non-linearity lower crossover value */
	tt = pow(s->Fl * s->nldlimit, 0.42);
	s->nldxval = 400.0 * tt / (tt + 27.13) + 0.1;

	/* Non-linearity lower crossover slope, passes through 0 */
	if (s->nldlimit > 0.0) {
		s->nldxslope = (s->nldxval-0.1)/s->nldlimit;
	}
//printf("~1 nldxval = %f, nlxdslope = %f\n",s->nldxval,s->nldxslope);

	/* Non-linearity upper crossover value */
	tt = pow(s->Fl * s->nlulimit, 0.42);
	s->nluxval = 400.0 * tt / (tt + 27.13) + 0.1;

	/* Non-linearity upper crossover slope, set to asymtope */
	t1 = s->Fl * s->nlulimit;
	t2 = pow(t1, 0.42) + 27.13;
	s->nluxslope = 0.42 * s->Fl * 400.0 * 27.13 / (pow(t1,0.58) * t2 * t2);


	/* Limited A value at J = JLIMIT */
	s->lA = pow(s->jlimit, 1.0/(s->C * s->z)) * s->Aw;

#ifdef DIAG1
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
	printf("                            nn = %f\n", s->nn);
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

/* Conversions. Return values are always 0 */
static int XYZ_to_cam(
struct _cam02 *s,
double Jab[3],
double XYZ[3]
) {
	int i;
	double xyz[3], rgb[3], rgbp[3], rgba[3], rgbc[3];
	double a, b, ja, jb, J, JJ, C, h, e, A, ss;
	double ttA, rS, cJ, tt;
	double k1, k2, k3;
	int clip = 0;			// ~~999

#ifdef DIAG2		/* Incase in == out */
	double XYZi[3];

	XYZi[0] = XYZ[0];
	XYZi[1] = XYZ[1];
	XYZi[2] = XYZ[2];
#endif

	TRACE(("\nForward conversion:\n"))
	TRACE(("XYZ %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]))

#ifdef DISABLE_MATRIX

	rgba[0] = XYZ[0];
	rgba[1] = XYZ[1];
	rgba[2] = XYZ[2];

#else /* !DISABLE_MATRIX */

	/* Add in flare */
	xyz[0] = s->Fsc * XYZ[0] + s->Fsxyz[0];
	xyz[1] = s->Fsc * XYZ[1] + s->Fsxyz[1];
	xyz[2] = s->Fsc * XYZ[2] + s->Fsxyz[2];

	/* Try and prevent crazy blue behaviour */
#ifndef DISABLE_BLUECLIP
	if (s->noclip == 0) {
		double lz, lxy;

		lz = xyz[2];
		lxy = 0.7328/0.1624 * xyz[0] + 0.4296/0.1624 * xyz[1];

		if (lxy > 1e-6 && lz > 1e-6) {
			lxy = BLUECLIMIT * lxy;		/* Set limit */
			lxy = 1.0/(lxy * lxy);
			lz = 1.0/sqrt(lxy + 1.0/(lz * lz));
		}
		
		/* Spectrally sharpened cone responses */
		rgb[0] =  0.7328 * xyz[0] + 0.4296 * xyz[1] - 0.1624 * lz;
		rgb[1] = -0.7036 * xyz[0] + 1.6975 * xyz[1] + 0.0061 * lz;
		rgb[2] =  0.0000 * xyz[0] + 0.0000 * xyz[1] + 1.0000 * lz;

	} else
#endif /* !DISABLE_BLUECLIP */
	{
		/* Spectrally sharpened cone responses */
		rgb[0] =  0.7328 * xyz[0] + 0.4296 * xyz[1] - 0.1624 * xyz[2];
		rgb[1] = -0.7036 * xyz[0] + 1.6975 * xyz[1] + 0.0061 * xyz[2];
		rgb[2] =  0.0000 * xyz[0] + 0.0000 * xyz[1] + 1.0000 * xyz[2];
	}

	
	/* Chromaticaly transformed sample value */
	rgbc[0] = s->Drgb[0] * rgb[0];
	rgbc[1] = s->Drgb[1] * rgb[1];
	rgbc[2] = s->Drgb[2] * rgb[2];
	
	/* Transform from spectrally sharpened, to Hunt-Pointer_Estevez cone space */
	rgbp[0] =  0.7409744840453773 * rgbc[0]
	        +  0.2180245944753982 * rgbc[1]
	        +  0.0410009214792244 * rgbc[2];
	rgbp[1] =  0.2853532916858801 * rgbc[0]
	        +  0.6242015741188157 * rgbc[1]
	        +  0.0904451341953042 * rgbc[2];
	rgbp[2] = -0.0096276087384294 * rgbc[0]
	        -  0.0056980312161134 * rgbc[1]
	        +  1.0153256399545427 * rgbc[2];

	TRACE(("rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]))
#ifdef DISABLE_NONLIN
	rgba[0] = 400.0/27.13 * rgbp[0];
	rgba[1] = 400.0/27.13 * rgbp[1];
	rgba[2] = 400.0/27.13 * rgbp[2];
#else	/* !DISABLE_NONLIN */
	/* Post-adapted cone response of sample. */
	/* rgba[] has a minimum value of 0.1 for XYZ[] = 0 and no flare. */
	/* We add a negative compression region, plus a linear segment at */
	/* the end of the +ve conversion to allow numerical handling of a */
	/* very wide range of values. */

	for (i = 0; i < 3; i++) {
		if (rgbp[i] < s->nldlimit) {
			rgba[i] = 0.1 + s->nldxslope * spow(rgbp[i], s->nldpow);
		} else {
			if (rgbp[i] <= s->nlulimit) {
				tt = pow(s->Fl * rgbp[i], 0.42);
				rgba[i] = 400.0 * tt / (tt + 27.13) + 0.1;
			} else {
				rgba[i] = s->nluxval + s->nluxslope * (rgbp[i] - s->nlulimit);
			}
		}
	}
#endif	/* !DISABLE_NONLIN */

#endif /* !DISABLE_MATRIX */

	/* Note that the minimum values of rgba[] for XYZ = 0 is 0.1, */
	/* hence magic 0.305 below comes from the following weighting of rgba[], */
	/* to base A at 0.0 */
	/* Deal with the core rgb to A, S & h conversion: */

	TRACE(("rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]))

	/* Preliminary Acromatic response */
	ttA = s->VttA[0] * rgba[0] + s->VttA[1] * rgba[1] + s->VttA[2] * rgba[2];

	/* Achromatic response */
	A = (ttA - 0.305) * s->Nbb;

	/* Preliminary red-green & yellow-blue opponent dimensions */
	/* a, b & ttd form an (almost) orthogonal coordinate set. */
	/* ttA is in a different direction */
	a = s->Va[0] * rgba[0] + s->Va[1] * rgba[1] + s->Va[2] * rgba[2];
	b = s->Vb[0] * rgba[0] + s->Vb[1] * rgba[1] + s->Vb[2] * rgba[2];

	/* restricted Saturation to stop divide by zero */
	/* (The exact value isn't important because the numerator dominates as a,b aproach 0 */
	rS = sqrt(a * a + b * b);
	if (rS < DBL_EPSILON)
		rS = DBL_EPSILON;
	TRACE(("ttA = %f, a = %f, b = %f, rS = %f, A = %f\n", ttA,a,b,rS,A))

	/* Lightness J, Derived directly from Acromatic response. */
	/* Cuttover to a straight line segment when J < 0.005, */
#ifdef NEVER		/* Cut to a straight line */
	if (A >= s->lA) {
		J = pow(A/s->Aw, s->C * s->z);		/* J/100  - keep Sign */
	} else {
		J = s->jlimit/s->lA * A;			/* Straight line */
		TRACE(("limited Acromatic to straight line\n"))
	}
#else			/* Symetric */
	if (A >= 0.0) {
		J = pow(A/s->Aw, s->C * s->z);		/* J/100  - keep Sign */
	} else {
		J = -pow(-A/s->Aw, s->C * s->z);		/* J/100  - keep Sign */
		TRACE(("symetric Acromatic\n"))
	}
#endif

	/* Constraied (+ve, non-zero) J */
	if (A > 0.0) {
		cJ = pow(A/s->Aw, s->C * s->z);
		if (cJ < s->ssmincj)
			cJ = s->ssmincj;
	} else {
		cJ = s->ssmincj;
	}

	TRACE(("J = %f, cJ = %f\n",J,cJ))

	/* Final hue angle */
	h = (180.0/DBL_PI) * atan2(b,a);
	h = (h < 0.0) ? h + 360.0 : h;

	/* Eccentricity factor */
	e = (12500.0/13.0 * s->Nc * s->Ncb * (cos(h * DBL_PI/180.0 + 2.0) + 3.8));

	/* ab scale components */
	k1 = pow(s->nn, 1.0/0.9) * e * pow(cJ, 1.0/1.8)/pow(rS, 1.0/9.0);
	k2 = pow(cJ, 1.0/(s->C * s->z)) * s->Aw/s->Nbb + 0.305;
	k3 = s->dcomp[1] * a + s->dcomp[2] * b;

	TRACE(("Raw k1 = %f, k2 = %f, k3 = %f, raw ss = %f\n",k1, k2, k3, pow(k1/(k2 + k3), 0.9)))

#ifdef TRACKMINMAX
	{
		int sno;
		double lrat, urat;

		ss = pow(k1/(k2 + k3), 0.9);

		if (ss < minss)
			minss = ss;
		if (ss > maxss)
			maxss = ss;

		lrat = -k3/k2;
		urat =  k3 * pow(ss, 10.0/9.0) / k1; 
		if (lrat > minlrat)
			minlrat = lrat;
		if (urat > maxurat)
			maxurat = urat;

		/* Record distribution of ss min/max vs. J for */
		/* regions outside a,b == 0 */

		sno = (int)(J * 100.0 + 0.5);
		if (sno < 0) {
			sno = 101;
			if (J < minj)
				minj = J;
		} else if (sno > 100) {
			sno = 102;
			if (J > maxj)
				maxj = J;
		}
		if (slotsd[sno] < lrat)
			slotsd[sno] = lrat;
		if (slotsu[sno] < urat)
			slotsu[sno] = urat;
	}
#endif /* TRACKMINMAX */

#ifdef DISABLE_TTD 
	ss = pow((k1/k2), 0.9);

#else	/* !TRACKMINMAX */

#ifdef ENABLE_DDL	

#ifndef NEVER
	/* Limit ratio of k3 to k2 to stop zero or -ve ss */
	if (k3 < -k2 * s->ddllimit) {
		k3 = -k2 * s->ddllimit;
		TRACE(("k3 limited to %f due to k3:k2 ratio, ss = %f\n",k3,pow(k1/(k2 + k3), 0.9)))
		clip = 1;
	}
#else
	/* Alternative to above ?? */
	/* Limit ratio of k2 to k3 to stop zero or -ve ss */
	if (k2 < -k3/s->ddllimit) {
		k2 = -k3/s->ddllimit;
		TRACE(("k2 limited to %f due to k3:k2 ratio, ss = %f\n",k2,pow(k1/(k2 + k3), 0.9)))
		clip = 1;
	}
#endif 	/* NEVER */

#ifndef NEVER
	/* See if there is going to be a problem in bwd, and limit k3 if there is */
	if (k3 > (k2 * s->ddulimit/(1.0 - s->ddulimit))) {
		k3 = (k2 * s->ddulimit/(1.0 - s->ddulimit));
		TRACE(("k3 limited to %f to allow for bk3:bk1 bwd limit\n",k3, ss))
		clip = 1;
	}
#else
	/* Alternative to above ?? */
	NO = should limit k1 as alternative !
	/* See if there is going to be a problem in bwd, and limit k2 if ther is */
	if (k2 < k3 * (1.0 - s->ddulimit)/s->ddulimit) {
		k2 = k3 * (1.0 - s->ddulimit)/s->ddulimit;
		TRACE(("k2 limited to %f to allow for bk3:bk1 bwd limit\n",k2, ss))
		clip = 1;
	}
#endif

#endif /* ENABLE_DDL */	

	/* Compute the ab scale factor */
	ss = pow(k1/(k2 + k3), 0.9);

#endif	/* !TRACKMINMAX */

#ifdef ENABLE_SS
	if (ss < SSLLIMIT)
		ss = SSLLIMIT;
	else if (ss > SSULIMIT)
		ss = SSULIMIT;
#endif /* ENABLE_SS */

#ifdef NEVER
// -------------------------------------------
	/* Show ss components ~~999 */
	if (s->retss) {
		Jab[0] = 1.0;
		if (clip)
			Jab[0] = 0.0;
		Jab[1] = (log10(ss) - +1.5)/(2.0 - +1.5);
		Jab[2] = (log10(ss) - +1.0)/(2.5 - +1.0);

		for (i = 0; i < 3; i++) {
			if (Jab[i] < 0.0)
				Jab[i] = 0.0;
			else if (Jab[i] > 1.0)
				Jab[i] = 1.0;
		}
		return 0;
	}
// -------------------------------------------
#endif /* NEVER */

	ja = a * ss;
	jb = b * ss;

	/* Chroma - always +ve, used just for HHKR */
	C = sqrt(ja * ja + jb * jb);
	
	TRACE(("ss = %f, A = %f, J = %f, C = %f, h = %f\n",ss,A,J,C,h))

	JJ = J;
#ifndef DISABLE_HHKR
 	/* Helmholtz-Kohlraush effect */
	if (s->hk && J < 1.0) {
		double kk = C/300.0 * sin(DBL_PI * fabs(0.5 * (h - 90.0))/180.0);
		if (kk > 1e-6) 	/* Limit kk to a reasonable range */
			kk = 1.0/(s->hklimit + 1.0/kk);
		JJ = J + (1.0 - (J > 0.0 ? J : 0.0)) * kk;
		TRACE(("JJ = %f from J = %f, kk = %f\n",JJ,J,kk))
	}
#endif /* DISABLE_HHKR */

	JJ *= 100.0;	/* Scale J */

	/* Compute Jab value */
	Jab[0] = JJ;
	Jab[1] = ja;
	Jab[2] = jb;

#ifdef NEVER	/* Brightness/Colorfulness */
	{
		double M, Q;

		Q = (4.0/s->C) * sqrt(JJ/100.0) * (s->Aw + 4.0) * pow(s->Fl, 0.25);
		M = C * pow(s->Fl, 0.25);

		printf("Lightness = %f, Chroma = %f\n",J,C);
		printf("Brightness = %f, Colorfulness = %f\n",M,Q);
	}
#endif /* NEVER */

	TRACE(("Jab %f %f %f\n",Jab[0], Jab[1], Jab[2]))
	TRACE(("\n"))

#ifdef DIAG2
	printf("Processing:\n");
	printf("XYZ = %f %f %f\n", XYZi[0], XYZi[1], XYZi[2]);
	printf("Including flare XYZ = %f %f %f\n", xyz[0], xyz[1], xyz[2]);
	printf("Sharpened cone sample rgb = %f %f %f\n", rgb[0], rgb[1], rgb[2]);
	printf("Chromatically transformed sample value rgbc = %f %f %f\n", rgbc[0], rgbc[1], rgbc[2]);
	printf("Hunt-P-E cone space rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]);
	printf("Post adapted cone response rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]);
	printf("Prelim red green a = %f, b = %f\n", a, b);
	printf("Hue angle h = %f\n", h);
	printf("Eccentricity factor e = %f\n", e);
	printf("Achromatic response A = %f\n", A);
	printf("Lightness J = %f, H.K. Lightness = %f\n", J * 100, JJ);
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
	double xyz[3], rgb[3], rgbp[3], rgba[3], rgbc[3];
	double a, b, ja, jb, J, JJ, C, rC, h, e, A, ss;
	double tt, cJ, ttA;
	double k1, k2, k3;

#ifdef DIAG2		/* Incase in == out */
	double Jabi[3];

	Jabi[0] = Jab[0];
	Jabi[1] = Jab[1];
	Jabi[2] = Jab[2];

#endif

	TRACE(("\nReverse conversion:\n"))
	TRACE(("Jab %f %f %f\n",Jab[0], Jab[1], Jab[2]))

	JJ = Jab[0] * 0.01;	/* J/100 */
	ja = Jab[1];
	jb = Jab[2];

	/* Convert Jab to core A, S & h values: */

	/* Compute hue angle */
	h = (180.0/DBL_PI) * atan2(jb, ja);
	h = (h < 0.0) ? h + 360.0 : h;
	
	/* Compute chroma value */
	C = sqrt(ja * ja + jb * jb);	/* Must be Always +ve, Can be NZ even if J == 0 */

	/* Preliminary Restricted chroma, always +ve and NZ */
	/* (The exact value isn't important because the numerator dominates as a,b aproach 0) */
	rC = C;
	if (rC < DBL_EPSILON)
		rC = DBL_EPSILON;

	J = JJ;
#ifndef DISABLE_HHKR
 	/* Undo Helmholtz-Kohlraush effect */
	if (s->hk && J < 1.0) {
		double kk = C/300.0 * sin(DBL_PI * fabs(0.5 * (h - 90.0))/180.0);
		if (kk > 1e-6) 	/* Limit kk to a reasonable range */
			kk = 1.0/(s->hklimit + 1.0/kk);
		J = (JJ - kk)/(1.0 - kk);
		if (J < 0.0)
			J = JJ - kk;
		TRACE(("J = %f from JJ = %f, kk = %f\n",J,JJ,kk))
	}
#endif /* DISABLE_HHKR */

	/* Achromatic response */
#ifdef NEVER	/* cut to straight line */
	if (J >= s->jlimit) {
		A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
	} else {	/* In the straight line segment */
		A = s->lA/s->jlimit * J;
		TRACE(("Undo Acromatic straight line\n"))
	}
#else			/* Symetric */
	if (J >= 0.0) {
		A = pow(J, 1.0/(s->C * s->z)) * s->Aw;
	} else {	/* In the straight line segment */
		A = -pow(-J, 1.0/(s->C * s->z)) * s->Aw;
		TRACE(("Undo symetric Acromatic\n"))
	}
#endif

	/* Preliminary Acromatic response +ve */ 
	ttA = (A/s->Nbb)+0.305;

	if (A > 0.0) {
		cJ = pow(A/s->Aw, s->C * s->z);
		if (cJ < s->ssmincj)
			cJ = s->ssmincj;
	} else {
		cJ = s->ssmincj;
	}
	TRACE(("C = %f, A = %f from J = %f, cJ = %f\n",C, A,J,cJ))

	/* Eccentricity factor */
	e = (12500.0/13.0 * s->Nc * s->Ncb * (cos(h * DBL_PI/180.0 + 2.0) + 3.8));

	/* ab scale components */
	k1 = pow(s->nn, 1.0/0.9) * e * pow(cJ, 1.0/1.8)/pow(rC, 1.0/9.0);
	k2 = pow(cJ, 1.0/(s->C * s->z)) * s->Aw/s->Nbb + 0.305;
	k3 = s->dcomp[1] * ja + s->dcomp[2] * jb;

	TRACE(("Raw k1 = %f, k2 = %f, k3 = %f, raw ss = %f\n",k1, k2, k3, (k1 - k3)/k2))

#ifdef DISABLE_TTD 
	ss = k1/k2;
#else	/* !DISABLE_TTD */
#ifdef ENABLE_DDL	

#ifndef NEVER
	/* Limit ratio of k3 to k1 to stop zero or -ve ss */
	if (k3 > (k1 * s->ddulimit)) {
		k3 = k1 * s->ddulimit;
		TRACE(("k3 limited to %f due to k3:k1 ratio, ss = %f\n",k3,(k1 - k3)/k2))
	}
#else
	/* Alternate to above ?? */
	/* Limit ratio of k2 to k3 to stop zero or -ve ss */
	if (k1 < (k3 / s->ddulimit)) {
		k1 = (k3 / s->ddulimit);
		TRACE(("k1 limited to %f due to k3:k1 ratio, ss = %f\n",k1,(k1 - k3)/k2))
	}
#endif

#ifndef NEVER
	/* See if there is going to be a problem in fwd */
	if (k3 < -k1 * s->ddllimit/(1.0 - s->ddllimit)) {
		/* Adjust ss to allow for fwd limitd computation */
		k3 = -k1 * s->ddllimit/(1.0 - s->ddllimit);
		TRACE(("k3 set to %f to allow for fk3:fk2 fwd limit\n",k3))
	}
#else
	/* Alternate to above ?? */
	/* See if there is going to be a problem in fwd */
	NO, should limit k2 as alternate to above
	if (k1 < -k3 * (1.0 - s->ddllimit)/s->ddllimit) {
		/* Adjust ss to allow for fwd limitd computation */
		k1 = -k3 * (1.0 - s->ddllimit)/s->ddllimit;
		TRACE(("k1 set to %f to allow for fk3:fk2 fwd limit\n",k1))
	}
#endif
#endif /* ENABLE_DDL */	

	/* Compute the ab scale factor */
	ss = (k1 - k3)/k2;

#endif	/* !DISABLE_TTD */

#ifdef ENABLE_SS
	if (ss < SSLLIMIT)
		ss = SSLLIMIT;
	else if (ss > SSULIMIT)
		ss = SSULIMIT;
#endif /* ENABLE_SS */

	/* Unscale a and b */
	a = ja / ss;
	b = jb / ss;

	TRACE(("ss = %f, ttA = %f, a = %f, b = %f\n",ss,ttA,a,b))

	/* Solve for post-adapted cone response of sample */
	/* using inverse matrix on ttA, a, b */
	rgba[0] = (20.0/61.0) * ttA
	        + ((41.0 * 11.0)/(61.0 * 23.0)) * a
	        + ((288.0 * 1.0)/(61.0 * 23.0)) * b;
	rgba[1] = (20.0/61.0) * ttA
	        - ((81.0 * 11.0)/(61.0 * 23.0)) * a
	        - ((261.0 * 1.0)/(61.0 * 23.0)) * b;
	rgba[2] = (20.0/61.0) * ttA
	        - ((20.0 * 11.0)/(61.0 * 23.0)) * a
	        - ((20.0 * 315.0)/(61.0 * 23.0)) * b;

	TRACE(("rgba %f %f %f\n",rgba[0], rgba[1], rgba[2]))
	
#ifdef DISABLE_MATRIX

	XYZ[0] = rgba[0];
	XYZ[1] = rgba[1];
	XYZ[2] = rgba[2];

#else /* !DISABLE_MATRIX */

#ifdef DISABLE_NONLIN
	rgbp[0] = 27.13/400.0 * rgba[0];
	rgbp[1] = 27.13/400.0 * rgba[1];
	rgbp[2] = 27.13/400.0 * rgba[2];
#else	/* !DISABLE_NONLIN */

	/* Hunt-Pointer_Estevez cone space */
	/* (with linear segment at the +ve end) */
	for (i = 0; i < 3; i++) {
		if (rgba[i] < s->nldxval) {
			rgbp[i] = spow((rgba[i] - 0.1)/s->nldxslope, 1.0/s->nldpow);
		} else if (rgba[i] <= s->nluxval) {
			tt = rgba[i] - 0.1;
			rgbp[i] = pow((27.13 * tt)/(400.0 - tt), 1.0/0.42)/s->Fl;
		} else {
			rgbp[i] = s->nlulimit + (rgba[i] - s->nluxval)/s->nluxslope;
		}
	}
#endif /* !DISABLE_NONLIN */

	TRACE(("rgbp %f %f %f\n",rgbp[0], rgbp[1], rgbp[2]))

	/* Chromaticaly transformed sample value */
	rgbc[0] =  1.5591630679450694 * rgbp[0]
	        -  0.5447249392120921 * rgbp[1]
	        -  0.0144381287329769 * rgbp[2];
	rgbc[1] = -0.7143316061228973 * rgbp[0]
	        +  1.8503110081052354 * rgbp[1]
	        -  0.1359794019823380 * rgbp[2];
	rgbc[2] =  0.0107755854444190 * rgbp[0]
	        +  0.0052187506061015 * rgbp[1]
	        +  0.9840056639494795 * rgbp[2];

	/* Spectrally sharpened cone responses */
	rgb[0]  =  rgbc[0]/s->Drgb[0];
	rgb[1]  =  rgbc[1]/s->Drgb[1];
	rgb[2]  =  rgbc[2]/s->Drgb[2];

	/* XYZ values */
	xyz[0] =  1.0978566630062390 * rgb[0]
	       -  0.2778434300014611 * rgb[1]
	       +  0.1799867669952221 * rgb[2];
	xyz[1] =  0.4550526940154284 * rgb[0]
	       +  0.4739377688665520 * rgb[1]
	       +  0.0710095371180196 * rgb[2];
	xyz[2] =  0.0000000000000000 * rgb[0]
	       +  0.0000000000000000 * rgb[1]
	       +  1.0000000000000000 * rgb[2];

	/* Subtract flare */
	XYZ[0] = s->Fisc * (xyz[0] - s->Fsxyz[0]);
	XYZ[1] = s->Fisc * (xyz[1] - s->Fsxyz[1]);
	XYZ[2] = s->Fisc * (xyz[2] - s->Fsxyz[2]);

#endif /* !DISABLE_MATRIX */

	TRACE(("XYZ = %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]))
	TRACE(("\n"))

#ifdef DIAG2
	printf("Processing:\n");
	printf("Jab = %f %f %f\n", Jabi[0], Jabi[1], Jabi[2]);
	printf("Chroma C = %f\n", C);
	printf("Preliminary Saturation ss = %f\n", ss);
	printf("Lightness J = %f, H.K. Lightness = %f\n", J * 100, JJ * 100);
	printf("Achromatic response A = %f\n", A);
	printf("Eccentricity factor e = %f\n", e);
	printf("Hue angle h = %f\n", h);
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

