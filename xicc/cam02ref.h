
/* 
 * cam02, unoptimised, reference version for testing.
 *
 * Color Appearance Model.
 *
 * Author:  Graeme W. Gill
 * Date:    17/1/2004
 * Version: 1.00
 *
 * This file is based on cam97s3.c by Graeme Gill.
 *
 * Copyright 2004, 2007 Graeme W. Gill
 * Please refer to COPYRIGHT file for details.
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */


/* Note that XYZ values are normalised to 1.0 consistent */
/* with the ICC convention (not 100.0 as assumed by the CIECAM spec.) */
/* Note that all whites are assumed to be normalised (ie. Y = 1.0) */

#undef DIAG			/* Print internal value diagnostics for each conversion */

/* ---------------------------------- */

struct _cam02ref {
/* Public: */
	void (*del)(struct _cam02ref *s);	/* We're done with it */

	int (*set_view)(
		struct _cam02ref *s,
		ViewingCondition Ev,	/* Enumerated Viewing Condition */
		double Wxyz[3],	/* Reference/Adapted White XYZ (Y scale 1.0) */
		double La,		/* Adapting/Surround Luminance cd/m^2 */
		double Yb,		/* Luminance of Background relative to reference white (range 0.0 .. 1.0) */
		double Lv,		/* Luminance of white in the Viewing/Scene/Image field (cd/m^2) */
						/* Ignored if Ev is set */
		double Yf,		/* Flare as a fraction of the reference white (range 0.0 .. 1.0) */
		double Fxyz[3],	/* The Flare white coordinates (typically the Ambient color) */
		int hk			/* Flag, NZ to use Helmholtz-Kohlraush effect */
	);

	/* Conversions */
	int (*XYZ_to_cam)(struct _cam02ref *s, double *out, double *in);

/* Private: */
	/* Scene parameters */
	ViewingCondition Ev;	/* Enumerated Viewing Condition */
	double Wxyz[3];	/* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
	double Yb;		/* Relative Luminance of Background to reference white (Y range 0.0 .. 1.0) */
	double La;		/* Adapting/Surround Luminance cd/m^2 */
	double Yf;		/* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
	double Fxyz[3];	/* The Flare white coordinates (typically the Ambient color) */

	/* Internal parameters */
	double  C;		/* Surround Impact */
	double Nc;		/* Chromatic Induction */
	double  F;		/* Adaptation Degree */

	/* Pre-computed values */
	double Fsc;			/* Flare scale */
	double Fisc;		/* Inverse flare scale */
	double Fsxyz[3];	/* Scaled Flare color coordinates */
	double rgbW[3];		/* Sharpened cone response white values */
	double D;			/* Degree of chromatic adaption */
	double Drgb[3];		/* Chromatic transformation value */
	double rgbcW[3];	/* Chromatically transformed white value */
	double rgbpW[3];	/* Hunt-Pointer-Estevez cone response space white */
	double n;			/* Background induction factor */
	double nn;			/* Precomuted function of n */
	double Fl;			/* Lightness contrast factor ?? */
	double Nbb;			/* Background brightness induction factors */
	double Ncb;			/* Chromatic brightness induction factors */
	double z;			/* Base exponential nonlinearity */
	double rgbaW[3];	/* Post-adapted cone response of white */
	double Aw;			/* Achromatic response of white */

	/* Option flags */
	int hk;				/* Use Helmholtz-Kohlraush effect */
	int trace;			/* Trace internal values */
	double nnlimit;		/* Return error if nlinear is less than this */
	double jlimit;		/* Return error if J is less than this */

}; typedef struct _cam02ref cam02ref;

/* ---------------------------------- */

/* Utility function */
/* Return a viewing condition enumeration from the given Ambient and */
/* Adapting/Surround Luminance. */
static ViewingCondition cam02ref_Ambient2VC(
double La,		/* Ambient Luminence (cd/m^2) */
double Lv		/* Luminence of white in the Viewing/Scene/Image field (cd/m^2) */
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

static void cam02ref_free(cam02ref *s);
static int cam02ref_set_view(struct _cam02ref *s, ViewingCondition Ev, double Wxyz[3],
                    double Yb, double La, double Lv, double Yf, double Fxyz[3], int hk);
static int cam02ref_XYZ_to_cam(struct _cam02ref *s, double *Jab, double *xyz);

/* Create a cam02 conversion object, with default viewing conditions */
cam02ref *new_cam02ref(void) {
	cam02ref *s;
	double D50[3] = { 0.9642, 1.0000, 0.8249 };

	if ((s = (cam02ref *)malloc(sizeof(cam02ref))) == NULL) {
		fprintf(stderr,"cam02: malloc failed allocating object\n");
		exit(-1);
	}
	
	/* Initialise methods */
	s->del      = cam02ref_free;
	s->set_view = cam02ref_set_view;
	s->XYZ_to_cam = cam02ref_XYZ_to_cam;

	return s;
}

static void cam02ref_free(cam02ref *s) {
	if (s != NULL)
		free(s);
}

static int cam02ref_set_view(
cam02ref *s,
ViewingCondition Ev,	/* Enumerated Viewing Condition */
double Wxyz[3],	/* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
double Yb,		/* Relative Luminence of Background to reference white */
double La,		/* Adapting/Surround Luminance cd/m^2 */
double Lv,		/* Luminence of white in the Viewing/Scene/Image field (cd/m^2) */
				/* Ignored if Ev is set to other than vc_none */
double Yf,		/* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
double Fxyz[3],	/* The Flare white coordinates (typically the Ambient color) */
int hk			/* Flag, NZ to use Helmholtz-Kohlraush effect */
) {
	double tt;

	if (Ev == vc_none)		/* Compute enumerated viewing condition */
		Ev = cam02ref_Ambient2VC(La, Lv);
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
	s->rgbW[2] =  0.0000 * s->Wxyz[0] + 0.0000 * s->Wxyz[1] + 1.0000 * s->Wxyz[2];

	/* Degree of chromatic adaption */
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
	printf("Adaption Degree F = %f\n", s->F);

	printf("Pre-computed values\n");
	printf("Sharpened cone white rgbW = %f %f %f\n", s->rgbW[0], s->rgbW[1], s->rgbW[2]);
	printf("Degree of chromatic adaption D = %f\n", s->D);
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


#define REFTRACE(xxxx) if (s->trace) printf xxxx ;

/* Conversion. Returns NZ and -1, -1, -1 if input is out of range */
static int cam02ref_XYZ_to_cam(
struct _cam02ref *s,
double Jab[3],
double XYZ[3]
) {
	int i;
	double xyz[3], rgb[3], rgbp[3], rgba[3], rgbaW[3], rgbc[3], rgbcW[3];
	double a, b, rS, J, C, h, e, A, ss;
	double ttd, tt;

	REFTRACE(("\nReference forward conversion:\n"))
	REFTRACE(("XYZ %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]))

	/* Add in flare */
	xyz[0] = s->Fsc * XYZ[0] + s->Fsxyz[0];
	xyz[1] = s->Fsc * XYZ[1] + s->Fsxyz[1];
	xyz[2] = s->Fsc * XYZ[2] + s->Fsxyz[2];

	/* Spectrally sharpened cone responses */
	rgb[0] =  0.7328 * xyz[0] + 0.4296 * xyz[1] - 0.1624 * xyz[2];
	rgb[1] = -0.7036 * xyz[0] + 1.6975 * xyz[1] + 0.0061 * xyz[2];
	rgb[2] =  0.0000 * xyz[0] + 0.0000 * xyz[1] + 1.0000 * xyz[2];
	
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

	REFTRACE(("rgbp = %f %f %f\n", rgbp[0], rgbp[1], rgbp[2]))

	/* Post-adapted cone response of sample. */
	/* rgba[] has a minimum value of 0.1 for XYZ[] = 0 and no flare. */
	/* We add a symetric negative compression region, plus linear segments at */
	/* the ends of this conversion to allow numerical handling of a */
	/* very wide range of values. */
	for (i = 0; i < 3; i++) {
		if (rgbp[i] < 0.0 || rgbp[i] < s->nnlimit) {
			Jab[0] = Jab[1] = Jab[2] = -1.0;
			return 1; 
		} else {
			tt = pow(s->Fl * rgbp[i], 0.42);
			if (tt < 10824.87)
				rgba[i] = (400.1 * tt + 2.713) / (tt + 27.13);
			else
				rgba[i] = (tt + 2.713) / 27.13;
		}
	}

	REFTRACE(("rgba = %f %f %f\n", rgba[0], rgba[1], rgba[2]))

	/* Preliminary red-green & yellow-blue opponent dimensions */
	a     = rgba[0] - 12.0 * rgba[1]/11.0 + rgba[2]/11.0;
    b     = (1.0/9.0) * (rgba[0] + rgba[1] - 2.0 * rgba[2]);
	rS   = sqrt(a * a + b * b);		/* Normalised a, b */

	/* Preliminary Saturation */
	/* Note that the minimum values for rgba[] for XYZ = 0 is 0.1 */
	/* Hence magic 0.305 below comes from the following weighting of rgba[] */
	ttd = rgba[0] + rgba[1] + (21.0/20.0) * rgba[2];

	/* Achromatic response */
	/* Note that the minimum values of rgba[] for XYZ = 0 is 0.1, */
	/* hence magic 0.305 below comes from the following weighting of rgba[], */
	/* to base A at 0.0 */
	A = (2.0 * rgba[0] + rgba[1] + (1.0/20.0) * rgba[2] - 0.305) * s->Nbb;

	REFTRACE(("a = %f, b = %f, ttd = %f, rS = %f, A = %f\n", a, b, ttd, rS, A))

	/* Lightness */
	J = pow(A/s->Aw, s->C * s->z);		/* J/100  - keep Sign */

	/* Hue angle */
    h = (180.0/DBL_PI) * atan2(b,a);
	h = (h < 0.0) ? h + 360.0 : h;

	/* Eccentricity factor */
	e = (cos(h * DBL_PI/180.0 + 2.0) + 3.8);

	if (J < DBL_EPSILON || J < s->jlimit || ttd < DBL_EPSILON) {
		REFTRACE(("J = %f, ttd = %f, exit with error\n", J, ttd))
		Jab[0] = Jab[1] = Jab[2] = -1.0;
		return 1; 
	} 

	ss = (12500.0/13.0 * s->Nc * s->Ncb * rS * e) / ttd;

	/* Chroma - Keep C +ve and make sure J doesn't force it to 0  */
	C = pow(ss, 0.9) * sqrt(J) * s->nn;
	
	REFTRACE(("ss = %f, C = %f\n", ss, C))

 	/* Helmholtz-Kohlraush effect */
	if (s->hk && J < 1.0) {
		double JJ, kk = C/300.0 * sin(DBL_PI * fabs(0.5 * (h - 90.0))/180.0);
		if (kk > 0.9)		/* Limit kk to a reasonable range */
			kk = 0.9;
		JJ = J + (1.0 - J) * kk;
		REFTRACE(("JJ = %f from J = %f, kk = %f\n",JJ,J,kk))
		J = JJ;
	}

	J *= 100.0;		/* Scale J */

	/* Compute Jab value */
	Jab[0] = J;
	if (rS >= DBL_EPSILON) {
		Jab[1] = C * a/rS;
		Jab[2] = C * b/rS;
	} else {
		Jab[1] = 0.0;
		Jab[2] = 0.0;
	}

	REFTRACE(("Returning Jab %f %f %f\n", Jab[0],Jab[1],Jab[2]))

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



