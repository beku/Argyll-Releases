/* 
 * International Color Consortium color transform expanded support
 *
 * Author:  Graeme W. Gill
 * Date:    2/7/00
 * Version: 1.00
 *
 * Copyright 2000, 2003 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * Based on the old iccXfm class.
 */

/*
 * This module provides the expands icclib functionality
 * for matrix profiles.
 * This file is #included in xicc.c, to keep its functions private.
 */

/*
 * TTBD:
 *
 *      Some of the error handling is crude. Shouldn't use
 *      error(), should return status.
 *
 *		Should allow for offset in curves - this will greatly improve
 *		profile quality on non-calibrated displays. See spectro/dispcal.c
 *		spectro/moncurve.c. Use conjgrad() instead of powell() to speed things up.
 *      Note that if curves have scale, the scale will have to be
 *      normalized back to zero by scaling the matrix before storing
 *      the result in the ICC profile.
 */

#define USE_CIE94_DE	/* Use CIE94 delta E measure when creating fit */

/* Weights in shaper parameters, to minimise unconstrained "wiggles" */
#define MXNORDERS 30			/* Maximum shaper harmonic orders to use */
#define XSHAPE_MAG  5000.0		/* Overall shaper parameter magnitide */
#define XSHAPE_BASE  0.00001	/* 0 & 1 harmonic weight */
#define XSHAPE_HBASE 0.0001		/* 2nd and higher additional weight */

#undef DEBUG			/* Extra printfs */
#undef DEBUG_PLOT		/* Plot curves */

/* ========================================================= */
/* Forward and Backward Matrix type conversion */
/* Return 0 on success, 1 if clipping occured, 2 on other error */

/* Individual components of Fwd conversion: */
static int
icxLuMatrixFwd_curve (
icxLuMatrix *p,		/* This */
double *out,		/* Vector of output values */
double *in			/* Vector of input values */
) {
	return ((icmLuMatrix *)p->plu)->fwd_curve((icmLuMatrix *)p->plu, out, in);
}

static int
icxLuMatrixFwd_matrix (
icxLuMatrix *p,		/* This */
double *out,		/* Vector of output values */
double *in			/* Vector of input values */
) {
	return ((icmLuMatrix *)p->plu)->fwd_matrix((icmLuMatrix *)p->plu, out, in);
}

static int
icxLuMatrixFwd_abs (
icxLuMatrix *p,		/* This */
double *out,		/* Vector of output values */
double *in			/* Vector of input values */
) {
	int rv = 0;
	rv |= ((icmLuMatrix *)p->plu)->fwd_abs((icmLuMatrix *)p->plu, out, in);

	if (p->pcs == icxSigJabData) {
		p->cam->XYZ_to_cam(p->cam, out, out);
	}
	return rv;
}


/* Overall Fwd conversion */
static int
icxLuMatrixFwd_lookup (
icxLuBase *pp,		/* This */
double *out,		/* Vector of output values */
double *in			/* Vector of input values */
) {
	int rv = 0;
	icxLuMatrix *p = (icxLuMatrix *)pp;
	rv |= icxLuMatrixFwd_curve(p, out, in);
	rv |= icxLuMatrixFwd_matrix(p, out, out);
	rv |= icxLuMatrixFwd_abs(p, out, out);
	return rv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Given a relative XYZ or Lab PCS value, convert in the fwd direction into */ 
/* the nominated output PCS (ie. Absolute, Jab etc.) */
/* (This is used in generating gamut compression in B2A tables) */
void icxLuMatrix_fwd_relpcs_outpcs(
icxLuBase *pp,
icColorSpaceSignature is,		/* Input space, XYZ or Lab */
double *out, double *in) {
	icxLuMatrix *p = (icxLuMatrix *)pp;

	if (is == icSigLabData && p->natpcs == icSigXYZData) {
		icmLab2XYZ(&icmD50, out, in);
		icxLuMatrixFwd_abs(p, out, out);
	} else if (is == icSigXYZData && p->natpcs == icSigLabData) {
		icmXYZ2Lab(&icmD50, out, in);
		icxLuMatrixFwd_abs(p, out, out);
	} else {
		icxLuMatrixFwd_abs(p, out, in);
	}
}

/* - - - - - - - - - - - - - - - - - - - - - */
/* Individual components of Bwd conversion: */

static int
icxLuMatrixBwd_abs (
icxLuMatrix *p,		/* This */
double *out,		/* Vector of output values */
double *in			/* Vector of input values */
) {
	int rv = 0;

	if (p->pcs == icxSigJabData) {
		p->cam->cam_to_XYZ(p->cam, out, in);
		rv |= ((icmLuMatrix *)p->plu)->bwd_abs((icmLuMatrix *)p->plu, out, out);
	} else {
		rv |= ((icmLuMatrix *)p->plu)->bwd_abs((icmLuMatrix *)p->plu, out, in);
	}
	return rv;
}

static int
icxLuMatrixBwd_matrix (
icxLuMatrix *p,		/* This */
double *out,		/* Vector of output values */
double *in			/* Vector of input values */
) {
	return ((icmLuMatrix *)p->plu)->bwd_matrix((icmLuMatrix *)p->plu, out, in);
}

static int
icxLuMatrixBwd_curve (
icxLuMatrix *p,		/* This */
double *out,		/* Vector of output values */
double *in			/* Vector of input values */
) {
	return ((icmLuMatrix *)p->plu)->bwd_curve((icmLuMatrix *)p->plu, out, in);
}

/* Overall Bwd conversion */
static int
icxLuMatrixBwd_lookup (
icxLuBase *pp,		/* This */
double *out,		/* Vector of output values */
double *in			/* Vector of input values */
) {
	int rv = 0;
	icxLuMatrix *p = (icxLuMatrix *)pp;
	rv |= icxLuMatrixBwd_abs(p, out, in);
	rv |= icxLuMatrixBwd_matrix(p, out, out);
	rv |= icxLuMatrixBwd_curve(p, out, out);
	return rv;
}

static void
icxLuMatrix_free(
icxLuBase *p
) {
	p->plu->del(p->plu);
	if (p->cam != NULL)
		p->cam->del(p->cam);
	free(p);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Given a nominated output PCS (ie. Absolute, Jab etc.), convert it in the bwd */
/* direction into a relative XYZ or Lab PCS value */
/* (This is used in generating gamut compression in B2A tables) */
void icxLuMatrix_bwd_outpcs_relpcs(
icxLuBase *pp,
icColorSpaceSignature os,		/* Output space, XYZ or Lab */
double *out, double *in) {
	icxLuMatrix *p = (icxLuMatrix *)pp;

	icxLuMatrixFwd_abs(p, out, in);
	if (os == icSigXYZData && p->natpcs == icSigLabData) {
		icmLab2XYZ(&icmD50, out, out);
	} else if (os == icSigXYZData && p->natpcs == icSigLabData) {
		icmXYZ2Lab(&icmD50, out, out);
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static gamut *icxLuMatrixGamut(icxLuBase *plu, double detail); 

/* Do the basic icxLuMatrix creation and initialisation */
static icxLuMatrix *
alloc_icxLuMatrix(
	xicc                  *xicp,
	icmLuBase             *plu,			/* Pointer to Lu we are expanding (ours) */
	int                   dir,			/* 0 = fwd, 1 = bwd */
	int                   flags			/* clip, merge flags */
) {
	icxLuMatrix *p;

	if ((p = (icxLuMatrix *) calloc(1,sizeof(icxLuMatrix))) == NULL)
		return NULL;

	p->pp                = xicp;
	p->plu               = plu;
	p->del               = icxLuMatrix_free;
	p->lutspaces         = icxLutSpaces;
	p->spaces            = icxLuSpaces;
	p->get_native_ranges = icxLu_get_native_ranges;
	p->get_ranges        = icxLu_get_ranges;
	p->efv_wh_bk_points  = icxLuEfv_wh_bk_points;
	p->get_gamut         = icxLuMatrixGamut;
	p->fwd_relpcs_outpcs = icxLuMatrix_fwd_relpcs_outpcs;
	p->bwd_outpcs_relpcs = icxLuMatrix_bwd_outpcs_relpcs;

	p->nearclip = 0;				/* Set flag defaults */
	p->mergeclut = 0;
	p->noisluts = 0;
	p->noipluts = 0;
	p->nooluts = 0;

	p->intsep = 0;

	p->fwd_lookup = icxLuMatrixFwd_lookup;
	p->fwd_curve  = icxLuMatrixFwd_curve;
	p->fwd_matrix = icxLuMatrixFwd_matrix;
	p->fwd_abs    = icxLuMatrixFwd_abs;
	p->bwd_lookup = icxLuMatrixBwd_lookup;
	p->bwd_abs    = icxLuMatrixBwd_abs;
	p->bwd_matrix = icxLuMatrixBwd_matrix;
	p->bwd_curve  = icxLuMatrixBwd_curve;

	if (dir) {		/* Bwd */
		p->lookup     = icxLuMatrixBwd_lookup;
		p->inv_lookup = icxLuMatrixFwd_lookup;
	} else {
		p->lookup     = icxLuMatrixFwd_lookup;
		p->inv_lookup = icxLuMatrixBwd_lookup;
	}

	/* There are no matrix specific flags */
	p->flags = flags;

	/* Get details of internal, native color space */
	p->plu->lutspaces(p->plu, &p->natis, NULL, &p->natos, NULL, &p->natpcs);

	/* Get other details of conversion */
	p->plu->spaces(p->plu, NULL, &p->inputChan, NULL, &p->outputChan, NULL, NULL, NULL, NULL, NULL);

	return p;
}

/* We setup valid fwd and bwd component conversions, */
/* but setup only the asked for overal conversion. */
static icxLuBase *
new_icxLuMatrix(
xicc                  *xicp,
int                   flags,		/* clip, merge flags */
icmLuBase             *plu,			/* Pointer to Lu we are expanding */
icmLookupFunc         func,			/* Functionality requested */
icRenderingIntent     intent,		/* Rendering intent */
icColorSpaceSignature pcsor,		/* PCS override (0 = def) */
icxViewCond           *vc,			/* Viewing Condition (NULL if pcsor is not CIECAM) */
int                   dir			/* 0 = fwd, 1 = bwd */
) {
	icxLuMatrix *p;

	/* Do basic creation and initialisation */
	if ((p = alloc_icxLuMatrix(xicp, plu, dir, flags)) == NULL)
		return NULL;

	p->func = func;

	/* Init the CAM model */
	if (pcsor == icxSigJabData) {
		if (vc != NULL)		/* One has been provided */
			p->vc  = *vc;		/* Copy the structure */
		else
			xicc_enum_viewcond(xicp, &p->vc, -1, NULL, 0, NULL);	/* Use a default */
		p->cam = new_icxcam(cam_default);
		p->cam->set_view(p->cam, p->vc.Ev, p->vc.Wxyz, p->vc.La, p->vc.Yb, p->vc.Lv, p->vc.Yf,
		                 p->vc.Fxyz, XICC_USE_HK, XICC_NOCAMCL);
	} else {
		p->cam = NULL;
	}

	/* Remember the effective intent */
	p->intent = intent;

	/* Get the effective spaces */
	plu->spaces(plu, &p->ins, NULL, &p->outs, NULL, NULL, NULL, NULL, &p->pcs, NULL);

	/* Override with pcsor */
	if (pcsor == icxSigJabData) {
		p->pcs = pcsor;		
		if (func == icmBwd || func == icmGamut || func == icmPreview)
			p->ins = pcsor;
		if (func == icmFwd || func == icmPreview)
			p->outs = pcsor;
	}

	/* In general the native and effective ranges of the icx will be the same as the */
	/* underlying icm lookup object. */
	p->plu->get_lutranges(p->plu, p->ninmin, p->ninmax, p->noutmin, p->noutmax);
	p->plu->get_ranges(p->plu, p->inmin,  p->inmax,  p->outmin,  p->outmax);

	/* If we have a Jab PCS override, reflect this in the effective icx range. */
	/* Note that the ab ranges are nominal. They will exceed this range */
	/* for colors representable in L*a*b* PCS */
	if (p->ins == icxSigJabData) {
		p->inmin[0] = 0.0;		p->inmax[0] = 100.0;
		p->inmin[1] = -128.0;	p->inmax[1] = 128.0;
		p->inmin[2] = -128.0;	p->inmax[2] = 128.0;
	} else if (p->outs == icxSigJabData) {
		p->outmin[0] = 0.0;		p->outmax[0] = 100.0;
		p->outmin[1] = -128.0;	p->outmax[1] = 128.0;
		p->outmin[2] = -128.0;	p->outmax[2] = 128.0;
	} 

	return (icxLuBase *)p;
}

/* ========================================================== */
/* xicc creation code                                         */
/* ========================================================== */

/* Context for figuring input curves */
typedef struct {
	rspl *r;			/* Device -> PCS rspl */
	int linear;			/* Flag */
	double nmin, nmax;	/* PCS End points to linearise */
	double min, max;	/* device End points to linearise */
} mxinctx;

#define NPARMS (15 + 3 * MXNORDERS)

/* Context for optimising matrix */
typedef struct {
	int verb;				/* Verbose */
	int optdim;				/* Optimisation dimensions */
	int isGamma;			/* NZ if gamma + matrix, else shaper */
	int isShTRC;			/* NZ if shared TRC */
	int norders;			/* Number of shaper orders */
	double v[NPARMS];		/* Holder for parameters */
	double sa[NPARMS];		/* Initial search area */
							/* Rest are matrix : */
							/* 0  1  2   R   X   */
							/* 3  4  5 * G = Y   */
							/* 6  7  8   B   Z   */
							/* For Gamma: */
							/* 9, 10, 11 are gamma */
							/* Else for shaper only: */
							/* 9,  10, 11 are Offset */
							/* 12, 13, 14 are 0th harmonics */
							/* 15, 16, 17 are 1st harmonics */
							/* 18, 19, 20 are 2nd harmonics */
							/* 24, 25, 26 etc. */
	cow *points;			/* List of test points as dev->Lab */
	int nodp;				/* Number of data points */
} mxopt;

/* Per chanel function being optimised */
static void mxmfunc1(mxopt *p, int j, double *v, double *out, double *in) {
	double vv, g;
	int ps = 3;		/* Parameter spacing */

	vv = *in;

	if (p->isShTRC) {
		j = 0;
		ps = 1;				/* Only one channel */
	}


	if (p->isGamma) {		/* Pure Gamma */
		/* Apply gamma */
		g = v[9 + j];
		if (g <= 0.0)
			vv = 1.0;
		else {
			if (vv >= 0.0) {
				vv = pow(vv, g);
			} else {
				vv = -pow(-vv, g);
			}
		}
	} else {		/* Add extra shaper parameters */
		int ord;

		/* Apply gamma as order 0 */
		g = v[9 + ps + j];
		if (g <= 0.0)
			vv = 1.0;
		else {
			if (vv >= 0.0) {
				vv = pow(vv, g);
			} else {
				vv = -pow(-vv, g);
			}
		}

		/* Process all the shaper orders from high to low. */
		/* [These shapers were inspired by a Graphics Gem idea */
		/* (Gems IV, VI.3, "Fast Alternatives to Perlin's Bias and */
		/*  Gain Functions, pp 401). */
		/*  They have the nice properties that they are smooth, and */
		/*  can't be non-monotonic. The control parameter has been */
		/*  altered to have a range from -oo to +oo rather than 0.0 to 1.0 */
		/*  so that the search space is less non-linear. ] */
		for (ord = 1; ord < p->norders; ord++) {
			int nsec;			/* Number of sections */
			double sec;			/* Section */
			g = v[9 + ps + ord * ps + j];	/* Parameter */
	
			nsec = ord + 1;		/* Increase sections for each order */
	
			vv *= (double)nsec;
	
			sec = floor(vv);
			if (((int)sec) & 1)
				g = -g;				/* Alternate action in each section */
			vv -= sec;
			if (g >= 0.0) {
				vv = vv/(g - g * vv + 1.0);
			} else {
				vv = (vv - g * vv)/(1.0 - g * vv);
			}
			vv += sec;
			vv /= (double)nsec;
		}

		/* Apply offset */
		g = v[9 + j];		/* Offset value */

		if (g >= 1.0) {
			vv = 1.0;
		} else if (g > 0.0) {	/* g > 0 && g < 1.0 */
			vv = g + ((1.0 - g) * vv);
		}

		if (vv < 0.0)
			vv = 0.0;
		else if (vv > 1.0)
			vv = 1.0;
	}

	*out = vv;
}

/* Function being optimised */
static void mxmfunc(mxopt *p, double *v, double *xyz, double *in) {
	int j;
	double rgb[3];

	/* Apply per channel processing */
	for (j = 0; j < 3; j++)
		mxmfunc1(p, j, v, &rgb[j], &in[j]);

	/* Apply matrix */
	xyz[0] = v[0] * rgb[0] + v[1] * rgb[1] + v[2] * rgb[2];
	xyz[1] = v[3] * rgb[0] + v[4] * rgb[1] + v[5] * rgb[2];
	xyz[2] = v[6] * rgb[0] + v[7] * rgb[1] + v[8] * rgb[2];
}

/* return the sum of the squares of the current shaper parameters */
static double xshapmag(
mxopt  *p,			/* Base of optimisation structure */
double *v			/* Pointer to parameters */
) {
	double tt, w, tparam = 0.0;
	int f, g;

	if (p->isGamma) {		/* Pure Gamma only */
		return 0.0;
	}

	if (p->isShTRC) {
		/* Offset value */
		tt = v[9];
		tt *= tt;
		w = XSHAPE_BASE;
		tparam += w * tt;

		/* Shaper values */
		for (f = 0; f < p->norders; f++) {
			tt = v[12 + f];
			tt *= tt;
			if (f <= 1)	/* First or second curves */
				w = XSHAPE_BASE;
			else
				w = XSHAPE_BASE + (f-1) * XSHAPE_HBASE;
			tparam += w * tt;
		}
		return XSHAPE_MAG * tparam;
	}

	/* Offset value */
	for (g = 0; g < 3; g++) {
		tt = v[9 + g];
		tt *= tt;
		w = XSHAPE_BASE;
		tparam += w * tt;
	}
	/* Shaper values */
	for (f = 0; f < p->norders; f++) {
		for (g = 0; g < 3; g++) {
			tt = v[12 + 3 * f + g];
			tt *= tt;
			if (f <= 1)	/* First or second curves */
				w = XSHAPE_BASE;
			else
				w = XSHAPE_BASE + (f-1) * XSHAPE_HBASE;
			tparam += w * tt;
		}
	}
	return XSHAPE_MAG * tparam/3.0;
}

/* Matrix optimisation function handed to powell() */
static double mxoptfunc(void *edata, double *v) {
	mxopt *p = (mxopt *)edata;
	double rv = 0.0, smv;
	double xyz[3], lab[3];
	int i;

	for (i = 0; i < p->nodp; i++) {

		/* Apply our function */
//printf("%f %f %f -> %f %f %f\n", p->points[i].p[0], p->points[i].p[1], p->points[i].p[2],
//xyz[0], xyz[1], xyz[2]);
		mxmfunc(p, v, xyz, p->points[i].p);
	
		/* Convert to Lab */
		icmXYZ2Lab(&icmD50, lab, xyz);
//printf("%f %f %f -> %f %f %f, target %f %f %f\n", p->points[i].p[0], p->points[i].p[1], p->points[i].p[2],
//lab[0], lab[1], lab[2], p->points[i].v[0], p->points[i].v[1], p->points[i].v[2]);
	
		/* Accumulate total delta E squared */
#ifdef USE_CIE94_DE
		rv += p->points[i].w * icmCIE94sq(lab, p->points[i].v);
#else
		rv += p->points[i].w * icmLabDEsq(lab, p->points[i].v);
#endif
	}

	/* Normalise error to be an average delta E squared */
	rv /= (double)p->nodp;

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	smv = xshapmag(p, v);
	rv += smv;

#ifdef DEBUG
printf("~9(%f)mxoptfunc returning %f\n",smv,rv);
#endif

	return rv;
}

/* Matrix progress function handed to powell() */
static void mxprogfunc(void *pdata, int perc) {
	mxopt *p = (mxopt *)pdata;

	if (p->verb) {
		printf("\r% 3d%%",perc); 
		if (perc == 100)
			printf("\n");
		fflush(stdout);
	}
}


/* Create icxLuMatrix and undelying tone reproduction curves and */
/* colorant tags, initialised from the icc, and then overwritten */
/* by a conversion created from the supplied scattered data points. */

/* The scattered data is assumed to map Device -> native PCS (ie. dir = Fwd) */
/* NOTE:- in theory once this icxLuMatrix is setup, it can be */
/* called to translate color values. In practice I suspect */
/* that the icxLuMatrix hasn't been setup completely enough to allows this. */
/* Might be easier to close it and re-open it ? */
static icxLuBase *
set_icxLuMatrix(
xicc               *xicp,
icmLuBase          *plu,			/* Pointer to Lu we are expanding (ours) */	
int                flags,			/* white/black point flags */
int                nodp,			/* Number of points */
cow                *ipoints,		/* Array of input points in XYZ space */
double             dispLuminance,	/* > 0.0 if display luminance value and is known */
double             wpscale,			/* > 0.0 if input white point is to be scaled */
int                quality			/* Quality metric, 0..3 */
) {
	icxLuMatrix *p;						/* Object being created */
	icc *icco = xicp->pp;				/* Underlying icc object */
	icmLuMatrix *pmlu = (icmLuMatrix *)plu;	/* icc matrix lookup object */
	int luflags = 0;					/* icxLuMatrix alloc clip, merge flags */
	int isGamma = 0;					/* NZ if gamma rather than shaper */
	int isShTRC = 0;					/* NZ if shared TRCs */
	int inputChan = 3;					/* Must be RGB like */
	int outputChan = 3;					/* Must be the PCS */
	icmHeader *h = icco->header;		/* Pointer to icc header */
	int rsplflags = 0;					/* Flags for scattered data rspl */
	int e, f, i, j;
	int maxits = 200;					/* Optimisation stop params */
	double stopon = 0.01;				/* Absolute delta E change to stop on */
	cow *points;		/* Copy of ipoints */
	mxopt os;			/* Optimisation information */
	double rerr;

#ifdef DEBUG_PLOT
	#define	XRES 100
	double xx[XRES];
	double y1[XRES];
#endif /* DEBUG_PLOT */

	if (flags & ICX_VERBOSE)
		rsplflags |= RSPL_VERBOSE;

	luflags = flags;		/* Transfer straight though ? */

	/* Check out some things about the profile */
	{
		icmCurve *wor, *wog, *wob;
		wor = pmlu->redCurve;
		wog = pmlu->greenCurve;
		wob = pmlu->blueCurve;

		if (wor == wog) {
			if (wog != wob) {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_matrix: TRC sharing is inconsistent");
				return NULL;
			}
			isShTRC = 1;
		}
		if (wor->flag != wog->flag || wog->flag != wob->flag) {
			xicp->errc = 1;
			sprintf(xicp->err,"icx_set_matrix: TRC type is inconsistent");
			return NULL;
		}
		if (wor->flag == icmCurveGamma) {
			isGamma = 1;
		}
	}

	/* Do basic icxLu creation and initialisation */
	if ((p = alloc_icxLuMatrix(xicp, plu, 0, luflags)) == NULL)
		return NULL;

	p->func = icmFwd;		/* Assumed by caller */

	/* Get the effective spaces of underlying icm, and set icx the same */
	plu->spaces(plu, &p->ins, NULL, &p->outs, NULL, NULL, &p->intent, NULL, &p->pcs, NULL);

	/* For set_icx the effective pcs has to be the same as the native pcs */

	/* Sanity check for matrix */
	if (p->pcs != icSigXYZData) {
		p->pp->errc = 1;
		sprintf(p->pp->err,"Can't create matrix profile with PCS of Lab !");
		p->del((icxLuBase *)p);
		return NULL;
	}

	/* In general the native and effective ranges of the icx will be the same as the */
	/* underlying icm lookup object. */
	p->plu->get_lutranges(p->plu, p->ninmin, p->ninmax, p->noutmin, p->noutmax);
	p->plu->get_ranges(p->plu, p->inmin,  p->inmax,  p->outmin,  p->outmax);

	/* If we have a Jab PCS override, reflect this in the effective icx range. */
	/* Note that the ab ranges are nominal. They will exceed this range */
	/* for colors representable in L*a*b* PCS */
	if (p->ins == icxSigJabData) {
		p->inmin[0] = 0.0;		p->inmax[0] = 100.0;
		p->inmin[1] = -128.0;	p->inmax[1] = 128.0;
		p->inmin[2] = -128.0;	p->inmax[2] = 128.0;
	} else if (p->outs == icxSigJabData) {
		p->outmin[0] = 0.0;		p->outmax[0] = 100.0;
		p->outmin[1] = -128.0;	p->outmax[1] = 128.0;
		p->outmin[2] = -128.0;	p->outmax[2] = 128.0;
	} 

	/* ------------------------------- */

	/* Allocate the array passed to fit_rspl() */
	if ((points = (cow *)malloc(sizeof(cow) * nodp)) == NULL) {
		p->pp->errc = 2;
		sprintf(p->pp->err,"Allocation of scattered coordinate array failed");
		p->del((icxLuBase *)p);
		return NULL;
	}

	/* Setup points ready for optimisation */
	for (i = 0; i < nodp; i++) {
		for (e = 0; e < inputChan; e++)
			points[i].p[e] = ipoints[i].p[e];
		
		for (f = 0; f < outputChan; f++)
			points[i].v[f] = ipoints[i].v[f];

		points[i].w = ipoints[i].w;

		/* Make sure its Lab for delta E calculation */
		icmXYZ2Lab(&icmD50, points[i].v, points[i].v);
	}

	/* Setup for optimising run */
	if (flags & ICX_VERBOSE)
		os.verb = 1;
	else
		os.verb = 0;
	os.points = points;
	os.nodp   = nodp;
	os.isShTRC = 0;

	/* Set quality/effort  factors */
	if (quality >= 3) {			/* Ultra high */
		os.norders = 20;
		maxits = 5000;
		stopon = 0.000001;
	} else if (quality == 2) {	/* High */
		os.norders = 15;
		maxits = 2000;
		stopon = 0.00001;
	} else if (quality == 1) {	/* Medium */
		os.norders = 10;
		maxits = 1000;
		stopon = 0.0001;
	} else {					/* Low */
		os.norders = 5;
		maxits = 500;
		stopon = 0.001;
	}
	if (os.norders > MXNORDERS)
		os.norders = MXNORDERS;
	
	/* Set initial optimisation values */
	os.v[0] = 0.4;  os.v[1] = 0.4;  os.v[2] = 0.2;		/* Matrix */
	os.v[3] = 0.2;  os.v[4] = 0.8;  os.v[5] = 0.1;
	os.v[6] = 0.02; os.v[7] = 0.15; os.v[8] = 1.3;

	if (isGamma) {		/* Just gamma curve */
		os.isGamma = 1;
		os.optdim = 12;
		os.v[9] = os.v[10] = os.v[11] = 2.4;					/* Gamma */ 
	} else {		/* Creating input curves */
		os.isGamma = 0;
		os.optdim = 12 + 3 * os.norders;
		os.v[9] = os.v[10] = os.v[11] = 0.0;	/* Offset */
		os.v[12] = os.v[13] = os.v[14] = 2.0; 	/* 0th Harmonic */
		for (i = 15; i < os.optdim; i++)
			os.v[i] = 0.0; 						/* Higher orders */
	}

	/* Set search area to starting values */
	for (j = 0; j < os.optdim; j++)
		os.sa[j] = 0.2;					/* Matrix, Gamma, Offsets, harmonics */

	if (isShTRC) {						/* Adjust things for shared */
		os.isShTRC = 1;

		/* Pack red paramenters down */
		for (i = 9; i < os.optdim; i++) {
			os.v[i] = os.v[(i - 9) * 3 + 9];
			os.sa[i] = os.sa[(i - 9) * 3 + 9];
		}
		os.optdim = ((os.optdim - 9)/3) + 9;
	}

	if (os.verb)
		printf("Creating matrix and curves...\n"); 

	if (powell(&rerr, os.optdim, os.v, os.sa, stopon, maxits,
	           mxoptfunc, (void *)&os, mxprogfunc, (void *)&os) != 0)
		warning("Powell failed to converge, residual error = %f",rerr);

#ifdef NEVER
printf("Matrix = %f %f %f\n",os.v[0], os.v[1], os.v[2]);
printf("         %f %f %f\n",os.v[3], os.v[4], os.v[5]);
printf("         %f %f %f\n",os.v[6], os.v[7], os.v[8]);
if (isGamma) {		/* Creating input curves */
if (isShTRC) 
	printf("Gamma = %f\n",os.v[9]);
else
	printf("Gamma = %f %f %f\n",os.v[9], os.v[10], os.v[11]);
} else {		/* Creating input curves */
	if (isShTRC) 
		printf("Offset = %f\n",os.v[9]);
	else
		printf("Offset = %f %f %f\n",os.v[9], os.v[10], os.v[11]);
	for (j = 0; j < os.norders; j++) {
		if (isShTRC) 
			printf("%d harmonics = %f\n",j, os.v[10 + j]);
		else
			printf("%d harmonics = %f %f %f\n",j, os.v[12 + j * 3], os.v[13 + j * 3],
			                                           os.v[14 + j * 3]);
	}
}
#endif /* NEVER */

	/* Deal with white/black points */
	if (flags & (ICX_SET_WHITE | ICX_SET_BLACK)) {
		double wp[3];	/* Absolute White point in XYZ */
		double bp[3];	/* Absolute Black point in XYZ */

		if (flags & ICX_VERBOSE)
			printf("Find white & black points\n");

		icmXYZ2Ary(wp, icmD50); 		/* Set a default value - D50 */
		icmXYZ2Ary(bp, icmBlack); 		/* Set a default value - absolute black */

		/* Figure out the device values for white */
		if (h->deviceClass == icSigInputClass) {
			double dwhite[MXDI], dblack[MXDI];	/* Device white and black values */
			double Lmax = -1e60;
			double Lmin = 1e60;

			/* We assume that the input target is well behaved, */
			/* and that it includes a white and black point patch, */
			/* and that they have the extreme L values */
	
			/* Discover the white and black points */
			for (i = 0; i < nodp; i++) {
				if (points[i].v[0] > Lmax) {
					Lmax = 
					wp[0] = points[i].v[0];
					wp[1] = points[i].v[1];
					wp[2] = points[i].v[2];
					dwhite[0] = points[i].p[0];
					dwhite[1] = points[i].p[1];
					dwhite[2] = points[i].p[2];
				}
				if (points[i].v[0] < Lmin) {
					Lmin = 
					bp[0] = points[i].v[0];
					bp[1] = points[i].v[1];
					bp[2] = points[i].v[2];
					dblack[0] = points[i].p[0];
					dblack[1] = points[i].p[1];
					dblack[2] = points[i].p[2];
				}
			}

			/* Lookup device white/black values in model */
			mxmfunc(&os, os.v, wp, dwhite);
			mxmfunc(&os, os.v, bp, dblack);

			/* If we were given an input white point scale factor, apply it */
			if (wpscale >= 0.0) {
				wp[0] *= wpscale;
				wp[1] *= wpscale;
				wp[2] *= wpscale;
			}

		} else {	/* Assume Monitor class */

			switch(h->colorSpace) {

				case icSigCmyData: {
					double cmy[3];

					/* Lookup white value */
					for (e = 0; e < inputChan; e++)
						cmy[e] = 0.0;

					mxmfunc(&os, os.v, wp, cmy);

					if (flags & ICX_VERBOSE)
						printf("Initial white point = %f %f %f\n",wp[0],wp[1],wp[2]);

					/* Lookup black value */
					for (e = 0; e < inputChan; e++)
						cmy[e] = 1.0;

					mxmfunc(&os, os.v, bp, cmy);

					if (flags & ICX_VERBOSE)
						printf("Initial black point = %f %f %f\n",bp[0],bp[1],bp[2]);
					break;
				}

				case icSigRgbData: {
					double rgb[3];

					/* Lookup white value */
					for (e = 0; e < inputChan; e++)
						rgb[e] = 1.0;

					mxmfunc(&os, os.v, wp, rgb);

					if (flags & ICX_VERBOSE)
						printf("Initial white point = %f %f %f\n",wp[0],wp[1],wp[2]);

					/* Lookup black value */
					for (e = 0; e < inputChan; e++)
						rgb[e] = 0.0;

					mxmfunc(&os, os.v, bp, rgb);
					
					if (flags & ICX_VERBOSE)
						printf("Initial black point = %f %f %f\n",bp[0],bp[1],bp[2]);
					break;
				}

				default: {
					xicp->errc = 1;
					sprintf(xicp->err,"set_icxLuMatrix: can't handle color space %s",
					                           icm2str(icmColorSpaceSignature, h->colorSpace));
					p->del((icxLuBase *)p);
					return NULL;
					break;
				}
			}
		}

		/* If this is a display, adjust the white point to be */
		/* exactly Y = 1.0, and compensate the matrix, dispLuminance */
		/* and black point accordingly. */
		if (h->deviceClass == icSigDisplayClass) {
			double scale = 1.0/wp[1];
			int i;

			for (i = 0; i < 9; i++) {
				os.v[i] *= scale;
			}

			dispLuminance *= wp[1];

			for (i = 0; i < 3; i++) {
				wp[i] *= scale;
				bp[i] *= scale;
			}
		}

		/* Absolute luminance tag */
		if (flags & ICX_WRITE_WBL
		 && h->deviceClass == icSigDisplayClass
		 && dispLuminance > 0.0) {
			icmXYZArray *wo;
			if ((wo = (icmXYZArray *)icco->read_tag(
			           icco, icSigLuminanceTag)) == NULL)  {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_luminance: couldn't find luminance tag");
				p->del((icxLuBase *)p);
				return NULL;
			}
			if (wo->ttype != icSigXYZArrayType) {
				xicp->errc = 1;
				sprintf(xicp->err,"luminance: tag has wrong type");
				p->del((icxLuBase *)p);
				return NULL;
			}

			wo->size = 1;
			wo->allocate((icmBase *)wo);	/* Allocate space */
			wo->data[0].X = 0.0;
			wo->data[0].Y = dispLuminance;
			wo->data[0].Z = 0.0;

			if (flags & ICX_VERBOSE)
				printf("Display Luminance = %f\n", wo->data[0].Y);
		}

		/* Write white and black tags */
		if ((flags & ICX_WRITE_WBL)
		 && (flags & ICX_SET_WHITE)) { /* White Point Tag: */
			icmXYZArray *wo;
			if ((wo = (icmXYZArray *)icco->read_tag(
			           icco, icSigMediaWhitePointTag)) == NULL)  {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_white_black: couldn't find white tag");
				p->del((icxLuBase *)p);
				return NULL;
			}
			if (wo->ttype != icSigXYZArrayType) {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_white_black: white tag has wrong type");
				p->del((icxLuBase *)p);
				return NULL;
			}

			wo->size = 1;
			wo->allocate((icmBase *)wo);	/* Allocate space */
			wo->data[0].X = wp[0];
			wo->data[0].Y = wp[1];
			wo->data[0].Z = wp[2];

			if (flags & ICX_VERBOSE)
				printf("White point XYZ = %f %f %f\n",wp[0],wp[1],wp[2]);
		}
		if ((flags & ICX_WRITE_WBL)
		 && (flags & ICX_SET_BLACK)) { /* Black Point Tag: */
			icmXYZArray *wo;
			if ((wo = (icmXYZArray *)icco->read_tag(
			           icco, icSigMediaBlackPointTag)) == NULL)  {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_white_black: couldn't find black tag");
				p->del((icxLuBase *)p);
				return NULL;
				}
			if (wo->ttype != icSigXYZArrayType) {
				xicp->errc = 1;
				sprintf(xicp->err,"icx_set_white_black: black tag has wrong type");
				p->del((icxLuBase *)p);
				return NULL;
			}

			wo->size = 1;
			wo->allocate((icmBase *)wo);	/* Allocate space */
			wo->data[0].X = bp[0];
			wo->data[0].Y = bp[1];
			wo->data[0].Z = bp[2];

			if (flags & ICX_VERBOSE)
				printf("Black point XYZ = %f %f %f\n",bp[0],bp[1],bp[2]);
		}

		if (flags & ICX_VERBOSE)
			printf("Fixup matrix for white point\n");

		/* Fix matrix to be relative to D50 white point, rather than absolute */
		{
			icmXYZNumber swp;
			double mat[3][3];
			
			icmAry2XYZ(swp, wp);

			/* Transfer from parameter to matrix */
			mat[0][0] = os.v[0]; mat[0][1] = os.v[1]; mat[0][2] = os.v[2];
			mat[1][0] = os.v[3]; mat[1][1] = os.v[4]; mat[1][2] = os.v[5];
			mat[2][0] = os.v[6]; mat[2][1] = os.v[7]; mat[2][2] = os.v[8];

			/* Adapt matrix */
			icmChromAdaptMatrix(ICM_CAM_MULMATRIX | ICM_CAM_BRADFORD, icmD50, swp, mat);

			/* Transfer back to parameters */ 
			os.v[0] = mat[0][0]; os.v[1] = mat[0][1]; os.v[2] = mat[0][2];
			os.v[3] = mat[1][0]; os.v[4] = mat[1][1]; os.v[5] = mat[1][2];
			os.v[6] = mat[2][0]; os.v[7] = mat[2][1]; os.v[8] = mat[2][2];
			if (flags & ICX_VERBOSE) {
				printf("After white point adjust:\n");
				printf("Matrix = %f %f %f\n",os.v[0], os.v[1], os.v[2]);
				printf("         %f %f %f\n",os.v[3], os.v[4], os.v[5]);
				printf("         %f %f %f\n",os.v[6], os.v[7], os.v[8]);
			}
		}
	}

	if (flags & ICX_VERBOSE)
		printf("Done gamma/shaper and matrix creation\n");

	/* Write the gamma/shaper and matrix to the icc memory structures */
	if (!isGamma) {		/* Creating input curves */
		unsigned int ui;
		icmCurve *wor, *wog, *wob;
		wor = pmlu->redCurve;
		wog = pmlu->greenCurve;
		wob = pmlu->blueCurve;

		for (ui = 0; ui < wor->size; ui++) {
			double in, rgb[3];

			for (j = 0; j < 3; j++) {

				in = (double)ui / (wor->size - 1.0);
	
				mxmfunc1(&os, j, os.v, &rgb[j], &in);

			}
			wor->data[ui] = rgb[0];	/* Curve values 0.0 - 1.0 */
			if (!isShTRC) {
				wog->data[ui] = rgb[1];
				wob->data[ui] = rgb[2];
			}
		}
#ifdef DEBUG_PLOT
		/* Display the result fit */
		for (j = 0; j < 3; j++) {
			for (i = 0; i < XRES; i++) {
				double x, y;
				xx[i] = x = i/(double)(XRES-1);
				mxmfunc1(&os, j, os.v, &y, &x);
				y1[i] = y;
			}
			do_plot(xx,y1,NULL,NULL,XRES);
		}
#endif /* DEBUG_PLOT */


	} else {			/* Gamma */
		icmCurve *wor, *wog, *wob;
		wor = pmlu->redCurve;
		wog = pmlu->greenCurve;
		wob = pmlu->blueCurve;
		wor->data[0] = os.v[9];	/* Gamma values */
		if (!isShTRC) {
			wog->data[0] = os.v[10];
			wob->data[0] = os.v[11];
		}
	}

	/* Matrix values */
	{
		icmXYZArray *wor, *wog, *wob;
		wor = pmlu->redColrnt;
		wog = pmlu->greenColrnt;
		wob = pmlu->blueColrnt;
		wor->data[0].X = os.v[0]; wor->data[0].Y = os.v[3]; wor->data[0].Z = os.v[6];
		wog->data[0].X = os.v[1]; wog->data[0].Y = os.v[4]; wog->data[0].Z = os.v[7];
		wob->data[0].X = os.v[2]; wob->data[0].Y = os.v[5]; wob->data[0].Z = os.v[8];

		/* Load into pmlu matrix and inverse ??? */
	}

	/* Free the coordinate lists */
	free(points);

	if (flags & ICX_VERBOSE)
		printf("Profile done\n");

	return (icxLuBase *)p;
}

/* ========================================================= */

/* Given an xicc lookup object, returm a gamut object. */
/* Note that the PCS must be Lab or Jab */
/* Return NULL on error, check errc+err for reason */
static gamut *icxLuMatrixGamut(
icxLuBase   *plu,		/* this */
double       detail		/* gamut detail level, 0.0 = def */
) {
	xicc     *p = plu->pp;				/* parent xicc */
	icxLuMatrix *lumat = (icxLuMatrix *)plu;	/* Lookup xMatrix type object */
	icColorSpaceSignature pcs;
	icmLookupFunc func;
	double white[3], black[3], kblack[3];
	gamut *gam;
	int res;		/* Sample point resolution */
	int i, e;

	if (detail == 0.0)
		detail = 10.0;

	/* get some details */
	plu->spaces(plu, NULL, NULL, NULL, NULL, NULL, NULL, &func, &pcs);

	if (func != icmFwd && func != icmBwd) {
		p->errc = 1;
		sprintf(p->err,"Creating Gamut surface for anything other than Device <-> PCS is not supported.");
		return NULL;
	}

	if (pcs != icSigLabData && pcs != icxSigJabData) {
		p->errc = 1;
		sprintf(p->err,"Creating Gamut surface PCS of other than Lab or Jab is not supported.");
		return NULL;
	}

	gam = new_gamut(detail, pcs == icxSigJabData);

	/* Explore the gamut by itterating through */
	/* it with sample points in device space. */

	res = (int)(600.0/detail);	/* Establish an appropriate sampling density */

	if (res < 40)
		res = 40;

	/* Since matrix profiles can't be non-monotonic, */
	/* just itterate through the surface colors. */
	for (i = 0; i < 3; i++) {
		int co[3];
		int ep[3];
		int co_e = 0;

		for (e = 0; e < 3; e++) {
			co[e] = 0;
			ep[e] = res;
		}
		ep[i] = 2;

		while (co_e < 3) {
			double in[3];
			double out[3];
	
			for (e = 0; e < 3; e++) 		/* Convert count to input value */
				in[e] = co[e]/(ep[e]-1.0);
	
			/* Always use the device->PCS conversion */
			if (lumat->fwd_lookup((icxLuBase *)lumat, out, in) > 1)
				error ("%d, %s",p->errc,p->err);
	
			gam->expand(gam, out);
		
			/* Increment the counter */
			for (co_e = 0; co_e < 3; co_e++) {
				co[co_e]++;
				if (co[co_e] < ep[co_e])
					break;	/* No carry */
				co[co_e] = 0;
			}
		}
	}

#ifdef NEVER
	/* Try it twice */
	for (i = 0; i < 3; i++) {
		int co[3];
		int ep[3];
		int co_e = 0;

		for (e = 0; e < 3; e++) {
			co[e] = 0;
			ep[e] = res;
		}
		ep[i] = 2;

		while (co_e < 3) {
			double in[3];
			double out[3];
	
			for (e = 0; e < 3; e++) 		/* Convert count to input value */
				in[e] = co[e]/(ep[e]-1.0);
	
			/* Always use the device->PCS conversion */
			if (lumat->fwd_lookup((icxLuBase *)lumat, out, in) > 1)
				error ("%d, %s",p->errc,p->err);
	
			gam->expand(gam, out);
		
			/* Increment the counter */
			for (co_e = 0; co_e < 3; co_e++) {
				co[co_e]++;
				if (co[co_e] < ep[co_e])
					break;	/* No carry */
				co[co_e] = 0;
			}
		}
	}
#endif

#ifdef NEVER	// (doesn't seem to make much difference) 
	/* run along the primary ridges in more detail too */
	/* just itterate through the surface colors. */
	for (i = 0; i < 3; i++) {
		int j;
		double in[3];
		double out[3];
		
		res *= 4;

		for (j = 0; j < res; j++) {
			double vv = i/(res-1.0);

			in[0] = in[1] = in[2] = vv;
			in[i] = 0.0;

			if (lumat->fwd_lookup((icxLuBase *)lumat, out, in) > 1)
				error ("%d, %s",p->errc,p->err);
			gam->expand(gam, out);

			in[0] = in[1] = in[2] = 0.0;
			in[i] = vv;

			if (lumat->fwd_lookup((icxLuBase *)lumat, out, in) > 1)
				error ("%d, %s",p->errc,p->err);
			gam->expand(gam, out);
		}
	}
#endif

	/* Put the white and black points in the gamut */
	plu->efv_wh_bk_points(plu, white, black, kblack);
	gam->setwb(gam, white, black, kblack);

	/* set the cusp points by itterating through the 0 & 100% colorant combinations */
	{
		DCOUNT(co, 3, 3, 0, 0, 2);

		gam->setcusps(gam, 0, NULL);
		DC_INIT(co);
		while(!DC_DONE(co)) {
			int e;
			double in[3];
			double out[3];
	
			if (!(co[0] == 0 && co[1] == 0 && co[2] == 0)
			 && !(co[0] == 1 && co[1] == 1 && co[2] == 1)) {	/* Skip white and black */
				for (e = 0; e < 3; e++)
					in[e] = (double)co[e];
	
				/* Always use the device->PCS conversion */
				if (lumat->fwd_lookup((icxLuBase *)lumat, out, in) > 1)
					error ("%d, %s",p->errc,p->err);
				gam->setcusps(gam, 3, out);
			}

			DC_INC(co);
		}
		gam->setcusps(gam, 2, NULL);
	}

#ifdef NEVER		/* Not sure if this is a good idea ?? */
	gam->getwb(gam, NULL, NULL, white, black);	/* Get the actual gamut white and black points */
	gam->setwb(gam, white, black);				/* Put it back as colorspace one */
#endif

	return gam;
}

#ifdef DEBUG
#undef DEBUG
#endif
