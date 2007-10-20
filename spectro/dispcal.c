
/* 
 * Argyll Color Correction System
 * Display callibrator.
 *
 * Author: Graeme W. Gill
 * Date:   14/10/2005
 *
 * Copyright 1996 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* This program displays test patches, and takes readings from a display device */
/* in order to create a RAMDAC calibration curve (usually stored in the ICC vcgt tag) */

/* This is the third version of the program. */

/* TTBD

	We could improve the robustnes against LCD warm up
	effects by reading the white every N readings,
	and then linearly interpolating the white readings,
	and scaling them back to a reference white.

	Need to add flare measure/subtract, to improve
	projector calibration ? - need to add to dispread too.

	It may be an improvement to automatically set the black
	ajustment level using a heuristic based on the black level
	and the change in brightness that would result from 100%
	black point grey adjustment.

 */

#ifdef __MINGW32__
# define WINVER 0x0500
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include <stdarg.h>
#if defined (NT)
#include <conio.h>
#endif
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "xicc.h"
#include "xspect.h"
#include "xcolorants.h"
#include "cgats.h"
#include "insttypes.h"
#include "icoms.h"
#include "inst.h"
#include "dispwin.h"
#include "dispsup.h"
#include "rspl.h"
#include "moncurve.h"
#include "targen.h"
#include "ofps.h"
#include "icc.h"
#include "sort.h"

#include "spyd2setup.h"			/* Enable Spyder 2 access */

#undef RES_TEST		/* See if we can detect better than 8 bit operation */

#undef DEBUG
#undef DEBUG_OFFSET		/* Keep test window out of the way */
#undef DEBUG_PLOT		/* Plot curve each time around */
#undef CHECK_MODEL		/* Do readings to check the accuracy of our model */

/* Invoke with -dfake for testing with a fake device. */
/* Will use a fake.icm profile if present, or a built in fake */
/* device behaviour if not. */

#define COMPORT 1			/* Default com port 1..4 */
#define REFINE_GAIN 0.80	/* Refinement correction damping/gain */
#define MAX_RPTS 15			/* Maximum tries at making a sample meet the current threshold */
#define VER_RES 100			/* Verification resolution */
#define NEUTRAL_BLEND_RATE 8.0		/* Suddenness of transition for -k factor < 1.0 */
#define ADJ_JACOBIAN		/* Adjust the Jacobian predictor matrix each time */
#define JAC_COR_FACT 0.5	/* Amount to correct Jacobian by (to filter noise) */
#define MOD_DIST_POW 1.4	/* Power used to distribute test samples for model building */
#define REFN_DIST_POW 1.4	/* Power used to distribute test samples for grey axis refinement */
#define CHECK_DIST_POW 1.4	/* Power used to distribute test samples for grey axis refinement */
#define CLIP				/* Clip RGB during refinement */

#define VERBOUT stdout

#ifdef DEBUG_PLOT
#include "plot.h"
#endif

#if defined(__APPLE__)
#define DEF_GAMMA 1.8
#else
#define DEF_GAMMA 2.2
#endif

#define VLENGTH(aa) (sqrt(aa[0] * aa[0] + aa[1] * aa[1] + aa[2] * aa[2]))

/* Compute approximate technical gamma from black/50% grey/white readings, */
/* (assumes an input offset gamma curve shape) */
static double tech_gamma(double bY, double gY, double wY) {
	int i;
	double grat, brat, gioff, gvv, gamma = 2.2;

	grat = gY/wY;
	brat = bY/wY;
	
	/* Iterative solution */
	for (i = 0; i < 5; i++) {
		gioff = pow(brat, 1.0/gamma);
		gvv = 0.5 + (1.0 - 0.5) * gioff;
		gamma = log(grat) / log(gvv);
	}
	return gamma;
}

/* Compute approximate popular gamma from black/50% grey/white readings, */
/* (assumes a zero based gamma curve shape) */
static double pop_gamma(double bY, double gY, double wY) {
	int i;
	double grat, brat, gioff, gvv, gamma = 2.2;

	grat = gY/wY;
	brat = bY/wY;
	
	gamma = log(grat) / log(0.5);
	return gamma;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Sample points used in initial device model optimisation */

typedef struct {
	double dev[3];		/* Device values */
	double lab[3];		/* Read value */
	double w;			/* Weighting */
} optref;

/* - - - - - - - - - - - - - - - - - - - */
/* device RGB inverse solution code */

/* Context for calibration solution */
typedef struct {
	double wh[3];		/* White absolute XYZ value */
	double bk[3];		/* Black absolute XYZ value */

	/* Target model */
	int gammat;			/* Gamma type, 0 = number, 1 = Lab, 2 = sRGB */
	double tgamma;		/* Technical Gamma target */
	double gioff;		/* Gamma curve input zero offset */
	double tooff;		/* Target output offset */
	int nat;			/* Flag - nz if native white target */

	double twh[3];		/* Target white absolute XYZ value */
	icmXYZNumber twN;	/* Same as above as XYZNumber */

	double tbk[3];		/* Target black point color */
	double tbL[3];		/* Same as above as Lab */
	icmXYZNumber tbN;	/* Same as above as XYZNumber */

	double fm[3][3];	/* Forward, aprox. linear RGB -> XYZ */
	double bm[3][3];	/* Backwards, aprox. XYZ -> linear RGB */
	mcv *dcvs[3];		/* Device RGB channel to linearised RGB curves */

	mcv *rdac[3];		/* Current RGB to RGB ramdac curves */

	double xyz[3];		/* Target xyz value */

	/* optimisation information */
	int np;				/* Total number of optimisation parameters */
	int co[3];			/* Offset in the parameters to each curve offset */
	int nc[3];			/* Number of for each curve */
	int nrp;			/* Total number of reference points */
	optref *rp;			/* reference points */
	double *dtin_iv;	/* Temporary array :- dp for input curves */
} cctx;

/* Return the xyz that is predicted by our aproximate device model */
/* by the given device RGB. */
static void fwddev(cctx *x, double xyz[3], double rgb[3]) {
	double lrgb[3];
	int j;

//printf("~1 fwddev called with rgb %f %f %f\n",rgb[0],rgb[1],rgb[2]);

	/* Convert device RGB into linear light RGB via curves */
	for (j = 0; j < 3; j++)
		lrgb[j] = x->dcvs[j]->interp(x->dcvs[j], rgb[j]);

//printf("~1 fwddev got linear RGB %f %f %f\n",lrgb[0],lrgb[1],lrgb[2]);

	/* Convert linear light RGB into XYZ via the matrix */
	icmMulBy3x3(xyz, x->fm, lrgb);

//printf("~1 fwddev got final xyz %f %f %f\n",xyz[0],xyz[1],xyz[2]);
}

/* Return the closest device RGB predicted by our aprox. device model */
/* to generate the given xyz. */
/* Return 0 on exact match, 1 on nearest match */
static int invdev(cctx *x, double rgb[3], double xyz[3]) {
	double lrgb[3];
	int j;

//printf("~1 invdev called with xyz %f %f %f\n",xyz[0],xyz[1],xyz[2]);

	/* Convert XYZ to linear light RGB via the inverse matrix */
	icmMulBy3x3(lrgb, x->bm, xyz);
//printf("~1 invdev; lin light rgb = %f %f %f\n",lrgb[0],lrgb[1],lrgb[2]);

	/* Convert linear light RGB to device RGB via inverse curves */
	for (j = 0; j < 3; j++) {
		rgb[j] = x->dcvs[j]->inv_interp(x->dcvs[j], lrgb[j]);
#ifdef CLIP
		if (rgb[j] < 0.0)
			rgb[j] = 0.0;
		else if (rgb[j] > 1.0)
			rgb[j] = 1.0;
#endif
	}
//printf("~1 invdev; inverse curves rgb = %f %f %f\n",rgb[0],rgb[1],rgb[2]);

	return 0;
}


/* Overall optimisation support */

/* Set the optimsation parameter number and offset values in cctx, */
/* and return an array filled in with the current parameters. */
/* Allocate temporary arrays */
static double *dev_get_params(cctx *x) {
	double *p, *tp;
	int i, j;

	x->np = 9;
	for (i = 0; i < 3; i++)
		x->np += x->dcvs[i]->luord;

	if ((p = (double *)malloc(x->np * sizeof(double))) == NULL)
		error("dev_params malloc failed");
	
	tp = p;
	
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			*tp++ = x->fm[i][j];

	for (i = 0; i < 3; i++) {
		x->co[i] = tp - p;				/* Offset to start */
		for (j = 0; j < x->dcvs[i]->luord; j++)
			*tp++ = x->dcvs[i]->pms[j];
		x->nc[i] = (tp - p) - x->co[i];	/* Number */
	}
		
	if ((x->dtin_iv = (double *)malloc(x->np * sizeof(double))) == NULL)
		error("dev_params malloc failed");

	return p;
}

/* Given a set of parameters, put them back into the model */
static void dev_put_params(cctx *x, double *p) {
	int i, j;

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			x->fm[i][j] = *p++;
	for (i = 0; i < 3; i++)
		for (j = 0; j < x->dcvs[i]->luord; j++)
			x->dcvs[i]->pms[j] = *p++;
}

/* Device model optimisation function handed to powell() */
static double dev_opt_func(void *edata, double *v) {
	cctx *x = (cctx *)edata;
	int i, j;
	double tw = 0.0;
	double rv, smv;

#ifdef NEVER
	printf("params =");
	for (i = 0; i < x->np; i++)
		printf(" %f",v[i]);
	printf("\n");
#endif

	/* For all our data points */
	rv = 0.0;
	for (i = 0; i < x->nrp; i++) {
		double lrgb[3];		/* Linear light RGB */
		double xyz[3], lab[3];

		/* Convert through device curves */
		for (j = 0; j < 3; j++)
			lrgb[j] = x->dcvs[j]->interp_p(x->dcvs[j], v + x->co[j], x->rp[i].dev[j]);

		/* Convert linear light RGB into XYZ via the matrix */
		icxMulBy3x3Parm(xyz, v, lrgb);

		/* Convert to Lab */
		icmXYZ2Lab(&x->twN, lab, xyz);
	
		/* Compute delta E squared */
//printf("~1 point %d: Lab is %f %f %f, should be %f %f %f\n",
//i, lab[0], lab[1], lab[2], x->rp[i].lab[0], x->rp[i].lab[1], x->rp[i].lab[2]);
		rv += icmCIE94sq(lab, x->rp[i].lab) * x->rp[i].w;
		tw += x->rp[i].w;
	}

	/* Normalise error to be a weighted average delta E squared and scale smoothing */
	rv /= (tw * 1200.0);

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	smv = 0.0;
	for (j = 0; j < 3; j++)
		smv += x->dcvs[j]->shweight_p(x->dcvs[j], v + x->co[j], 1.0);
	rv += smv;

#ifdef NEVER
	printf("rv = %f (%f)\n",rv, smv);
#endif
	return rv;
}


/* Device model optimisation function handed to conjgrad() */
static double dev_dopt_func(void *edata, double *dv, double *v) {
	cctx *x = (cctx *)edata;
	int i, j, k;
	int f, ee, ff, jj;
	double tw = 0.0;
	double rv, smv;

	double dmato_mv[3][9];		/* Del in mat out due to del in matrix param vals */
	double dmato_tin[3][3];		/* Del in mat out due to del in matrix input values */
	double dout_lab[3][3];		/* Del in out due to XYZ to Lab conversion */
	double de_dout[2][3];		/* Del in delta E due to input Lab values */

#ifdef NEVER
	printf("params =");
	for (i = 0; i < x->np; i++)
		printf(" %f",v[i]);
	printf("\n");
#endif

	/* Zero the accumulated partial derivatives */
	for (i = 0; i < x->np; i++)
		dv[i] = 0.0;

	/* For all our data points */
	rv = 0.0;
	for (i = 0; i < x->nrp; i++) {
		double lrgb[3];		/* Linear light RGB */
		double xyz[3], lab[3];

		/* Apply the input channel curves */
		for (j = 0; j < 3; j++)
			lrgb[j] = x->dcvs[j]->dinterp_p(x->dcvs[j], v + x->co[j],
			                         x->dtin_iv + x->co[j], x->rp[i].dev[j]);

		/* Convert linear light RGB into XYZ via the matrix */
		icxdpdiMulBy3x3Parm(xyz, dmato_mv, dmato_tin, v, lrgb);
		
		/* Convert to Lab */
		icxdXYZ2Lab(&x->twN, lab, dout_lab, xyz);
	
		/* Compute delta E squared */
//printf("~1 point %d: Lab is %f %f %f, should be %f %f %f\n",
//i, lab[0], lab[1], lab[2], x->rp[i].lab[0], x->rp[i].lab[1], x->rp[i].lab[2]);
		rv += icxdCIE94sq(de_dout, lab, x->rp[i].lab) * x->rp[i].w;
		de_dout[0][0] *= x->rp[i].w;
		de_dout[0][1] *= x->rp[i].w;
		de_dout[0][2] *= x->rp[i].w;
		tw += x->rp[i].w;

		/* Compute and accumulate partial difference values for each parameter value */

		/* Input channel curves */
		for (ee = 0; ee < 3; ee++) {				/* Parameter input chanel */
			for (k = 0; k < x->nc[ee]; k++) {		/* Param within channel */
				double vv = 0.0;
				jj = x->co[ee] + k;					/* Overall input curve param */

				for (ff = 0; ff < 3; ff++) {		/* Lab channels */
					for (f = 0; f < 3; f++) {		/* XYZ channels */
						vv += de_dout[0][ff] * dout_lab[ff][f]
						    * dmato_tin[f][ee] * x->dtin_iv[jj];
					}
				}
				dv[jj] += vv;
			}
		}

		/* Matrix parameters */
		for (k = 0; k < 9; k++) {				/* Matrix parameter */
			double vv = 0.0;

			for (ff = 0; ff < 3; ff++) {		/* Lab channels */
				for (f = 0; f < 3; f++) {		/* XYZ channels */
					vv += de_dout[0][ff] * dout_lab[ff][f]
					    * dmato_mv[f][k];
				}
			}
			dv[k] += vv;
		}
	}

	/* Normalise error to be a weighted average delta E squared and scale smoothing */
	rv /= (tw * 1200.0);
	for (i = 0; i < x->np; i++)
		dv[i] /= (tw * 900.0);

	/* Sum with shaper parameters squared, to */
	/* minimise unsconstrained "wiggles" */
	smv = 0.0;
	for (j = 0; j < 3; j++)
		smv += x->dcvs[j]->dshweight_p(x->dcvs[j], v + x->co[j], x->dtin_iv + x->co[j], 1.0);
	rv += smv;

#ifdef NEVER
	printf("drv = %f (%f)\n",rv, smv);
#endif
	return rv;
}

#ifdef NEVER
/* Check partial derivative function within dev_opt_func() using powell() */

static double dev_opt_func(void *edata, double *v) {
	cctx *x = (cctx *)edata;
	int i;
	double dv[2000];
	double rv, drv;
	double trv;
	
	rv = dev_opt_func_(edata, v);
	drv = dev_dopt_func(edata, dv, v);

	if (fabs(rv - drv) > 1e-6) {
		printf("######## RV MISMATCH is %f should be %f ########\n",drv, rv);
		exit(0);
	}

	/* Check each parameter delta */
	for (i = 0; i < x->np; i++) {
		double del;

		v[i] += 1e-7;
		trv = dev_opt_func_(edata, v);
		v[i] -= 1e-7;
		
		/* Check that del is correct */
		del = (trv - rv)/1e-7;
		if (fabs(dv[i] - del) > 1.0) {
//printf("~1 del = %f from (trv %f - rv %f)/0.1\n",del,trv,rv);
			printf("######## EXCESSIVE at v[%d] is %f should be %f ########\n",i,dv[i],del);
			exit(0);
		}
	}
	return rv;
}
#endif

/* =================================================================== */
/* Structure to save aproximate model readings in */
typedef struct {
	double v;			/* Input value */
	double xyz[3];		/* Reading */
} sxyz;

/* ------------------------------------------------------------------- */
#ifdef __APPLE__

/* Workaround for a ppc gcc 3.3 optimiser bug... */
static int gcc_bug_fix(int i) {
	static int nn;
	nn += i;
	return nn;
}
#endif	/* APPLE */


/* =================================================================== */
/* Calibration sample point support. This allows the successive */
/* refinement of our neutral sample points */

/* A sample point */
typedef struct {
	double v;			/* Desired input value */
	double rgb[3];		/* Input value through calibration curves */
	double tXYZ[3];		/* Target XYZ */
	double XYZ[3];		/* Read XYZ */
	double deXYZ[3];	/* Delta XYZ wanted to target */
	double de;			/* Delta Lab to target */
	double dc;			/* Delta XYZ to target */

	double pXYZ[3];		/* Previous measured XYZ */
	double pdXYZ[3];	/* Delta XYZ intended from previous measure */
	double pdrgb[3];	/* Delta rgb made to previous */

	double dXYZ[3];		/* Actual delta XYZ resulting from previous delta rgb */

	double j[3][3];		/* Aproximate Jacobian (del RGB -> XYZ) */
	double ij[3][3];	/* Aproximate inverse Jacobian (del XYZ-> del RGB) */
} csp;

/* All the sample points */
typedef struct {
	int no;				/* Number of samples */
	int _no;			/* Allocation */
	csp *s;			/* List of samples */
} csamp;

static void free_alloc_csamp(csamp *p) {
	if (p->s != NULL)
		free(p->s);
	p->s = NULL;
}

/* Initialise v values */
static void init_csamp_v(csamp *p, cctx *x, int psrand) {
	int i, j;
	sobol *so = NULL;

	if (psrand != 0) {	/* Use pseudo random distribution for verification */
		if ((so = new_sobol(1)) == NULL)
			error("New sobol failed");
	}
	
	/* Generate the sample points */
	for (i = 0; i < p->no; i++) {
		double vv;

#ifdef __APPLE__
		gcc_bug_fix(i);
#endif
		if (so != NULL) {
			if (i == 0)
				vv = 1.0;
			else if (i == 1)
				vv = 0.0;
			else
				so->next(so, &vv);
		} else
			vv = i/(p->no - 1.0);
		vv = pow(vv, REFN_DIST_POW);	/* Skew sample points to be slightly perceptual */
		p->s[i].v = vv;
	}

	if (so != NULL) {
		/* Sort it so white is last */
#define HEAP_COMPARE(A,B) (A.v < B.v) 
	HEAPSORT(csp,p->s,p->no)
#undef HEAP_COMPARE
		so->del(so);
	}
}

/* Initialise txyz values from v values */
static void init_csamp_txyz(csamp *p, cctx *x) {
	int i, j;

	/* Set the sample points targets */
	for (i = 0; i < p->no; i++) {
		double y, vv;
		double Lab[3];
		double bl;

		vv = p->s[i].v;

		/* Compute target Y value for this device inpute */
		if (x->gammat == 0) {			/* Power curve */
			/* Allow for zero input offset */
			y = x->gioff + (1.0 - x->gioff) * vv;
			y = pow(y, x->tgamma);
		} else if (x->gammat == 2) {	/* sRGB curve */
			if (vv <= 0.03928)
				y = vv/12.92;
			else
				y = pow((0.055 + vv)/1.055, 2.4);
			/* Allow for zero output offset */
			y = x->tooff + (1.0 - x->tooff) * y;
		} else {	/* L* curve */
			y = 100.0 * vv;				/* L* target value */

			/* Convert L* to Y */
			y = (y + 16.0)/116.0;
			if (y > 24.0/116.0)
				y = pow(y,3.0);
			else
				y = (y - 16.0/116.0)/7.787036979;
			/* Allow for zero output offset */
			y = x->tooff + (1.0 - x->tooff) * y;
		}
		/* Convert Y to L* */
		if (y > 0.008856451586)
			y = pow(y,1.0/3.0);
		else
			y = 7.787036979 * y + 16.0/116.0;
		Lab[0] = 116.0 * y - 16.0;
		Lab[1] = Lab[2] = 0.0;	/* Target is neutral */

		/* Compute blended neutral target a* b* */
		bl = pow((1.0 - vv), NEUTRAL_BLEND_RATE);		/* Crossover near the black */
		Lab[1] = bl * x->tbL[1];
		Lab[2] = bl * x->tbL[2];

		icmLab2XYZ(&x->twN, p->s[i].tXYZ, Lab);		/* XYZ Value to aim for */

//printf("~1 %d: target XYZ %.2f %.2f %.2f, Lab %.2f %.2f %.2f\n",i, p->s[i].tXYZ[0],p->s[i].tXYZ[1],p->s[i].tXYZ[2], Lab[0],Lab[1],Lab[2]);
	}
}


/* Allocate the sample points and initialise them with the */
/* target device and XYZ values, and first cut device values. */
static void init_csamp(csamp *p, cctx *x, int doupdate, int verify, int psrand, int no) {
	int i, j;
	
	p->_no = p->no = no;

	if ((p->s = (csp *)malloc(p->_no * sizeof(csp))) == NULL)
		error("csamp malloc failed");

	/* Compute v and txyz */
	init_csamp_v(p, x, psrand);
	init_csamp_txyz(p, x);

	/* Generate the sample points */
	for (i = 0; i < no; i++) {
		double dd, vv;

#ifdef __APPLE__
		gcc_bug_fix(i);
#endif
		vv = p->s[i].v;

		if (verify == 2) {		/* Verifying through installed curve */
			/* Make RGB values the input value */
			p->s[i].rgb[0] = p->s[i].rgb[1] = p->s[i].rgb[2] = vv;

		} else if (doupdate) {	/* Start or verify through current cal curves */
			for (j = 0; j < 3; j++) {
				p->s[i].rgb[j] = x->rdac[j]->interp(x->rdac[j], vv);
#ifdef CLIP
				if (p->s[i].rgb[j] < 0.0)
					p->s[i].rgb[j] = 0.0;
				else if (p->s[i].rgb[j] > 1.0)
					p->s[i].rgb[j] = 1.0;
#endif
			}
		} else {		/* we have model */
			/* Lookup an initial device RGB for that target by inverting */
			/* the approximate forward device model */
			p->s[i].rgb[0] = p->s[i].rgb[1] = p->s[i].rgb[2] = vv;
			invdev(x, p->s[i].rgb, p->s[i].tXYZ);
		}
		/* Force white to be native if native flag set */
		if (x->nat && i == (no-1)) {
//printf("~1 Forcing white rgb to be 1,1,1\n");
			p->s[i].rgb[0] = p->s[i].rgb[1] = p->s[i].rgb[2] = 1.0;
		}
		
//printf("~1 Inital point %d rgb %f %f %f\n",i,p->s[i].rgb[0],p->s[i].rgb[1],p->s[i].rgb[2]);

		/* Compute the approximate inverse Jacobian at this point */
		/* by taking the partial derivatives wrt to each device */
		/* channel of our aproximate forward model */
		if (verify != 2) {
			double refXYZ[3], delXYZ[3];
			fwddev(x, refXYZ, p->s[i].rgb);
			if (vv < 0.5)
				dd = 0.02;
			else
				dd = -0.02;
			/* Matrix organization is J[XYZ][RGB] for del RGB->del XYZ*/
			for (j = 0; j < 3; j++) {
				p->s[i].rgb[j] += dd;
				fwddev(x, delXYZ, p->s[i].rgb);
				p->s[i].j[0][j] = (delXYZ[0] - refXYZ[0]) / dd;
				p->s[i].j[1][j] = (delXYZ[1] - refXYZ[1]) / dd;
				p->s[i].j[2][j] = (delXYZ[2] - refXYZ[2]) / dd;
				p->s[i].rgb[j] -= dd;
			}
			if (icmInverse3x3(p->s[i].ij, p->s[i].j)) {
				warning("dispcal: inverting Jacobian failed (1)");
				*((char *)0) = 55;
			}
		}
	}
}

/* Return a linear XYZ interpolation */
static void csamp_interp(csamp *p, double xyz[3], double v) {
	int i, j;
	double b;

	if (p->no < 2)
		error("Calling csamp_interp with less than two existing samples");

	/* Locate the pair surrounding our input value */
	for (i = 0; i < (p->no-1); i++) {
		if (v >= p->s[i].v && v <= p->s[i+1].v)
			break;
	}
	if (i >= (p->no-1))
		error("csamp_interp out of range");
	
	b = (v - p->s[i].v)/(p->s[i+1].v - p->s[i].v);

	for (j = 0; j < 3; j++) {
		xyz[j] = b * p->s[i+1].XYZ[j] + (1.0 - b) * p->s[i].XYZ[j];
	}
}

/* Re-initialise a CSP with a new number of points. */
/* Interpolate the device values and jacobian. */
/* Set the current rgb from the current RAMDAC curves if not verifying */
static void reinit_csamp(csamp *p, cctx *x, int verify, int psrand, int no) {
	csp *os;			/* Old list of samples */
	int ono;			/* Old number of samples */
	int i, j, k, m;
	
	if (no == p->no)
		return;			/* Nothing has changed */

	os = p->s;			/* Save the existing per point information */
	ono = p->no;

	init_csamp(p, x, 0, 2, psrand, no);

	p->_no = p->no = no;

	/* Interpolate the current device values */
	for (i = 0; i < no; i++) {
		double vv, b;

		vv = p->s[i].v;

		/* Locate the pair surrounding our target value */
		for (j = 0; j < ono-1; j++) {
			if (vv >= os[j].v && vv <= os[j+1].v)
				break;
		}
		if (j >= (ono-1))
			error("csamp interp. out of range");
		
		b = (vv - os[j].v)/(os[j+1].v - os[j].v);
	
		for (k = 0; k < 3; k++) {
			if (verify == 2) {

				p->s[i].rgb[k] = b * os[j+1].rgb[k] + (1.0 - b) * os[j].rgb[k];

			} else {	/* Lookup rgb from current calibration curves */
				for (m = 0; m < 3; m++) {
					p->s[i].rgb[m] = x->rdac[m]->interp(x->rdac[m], vv);
#ifdef CLIP
					if (p->s[i].rgb[m] < 0.0)
						p->s[i].rgb[m] = 0.0;
					else if (p->s[i].rgb[m] > 1.0)
						p->s[i].rgb[m] = 1.0;
#endif
				}
			}
			p->s[i].XYZ[k] = b * os[j+1].XYZ[k] + (1.0 - b) * os[j].XYZ[k];
			p->s[i].deXYZ[k] = b * os[j+1].deXYZ[k] + (1.0 - b) * os[j].deXYZ[k];
			p->s[i].pXYZ[k] = b * os[j+1].pXYZ[k] + (1.0 - b) * os[j].pXYZ[k];
			p->s[i].pdrgb[k] = b * os[j+1].pdrgb[k] + (1.0 - b) * os[j].pdrgb[k];
			p->s[i].dXYZ[k] = b * os[j+1].dXYZ[k] + (1.0 - b) * os[j].dXYZ[k];
			for (m = 0; m < 3; m++) 
				p->s[i].j[k][m] = b * os[j+1].j[k][m] + (1.0 - b) * os[j].j[k][m];

		}
		if (icmInverse3x3(p->s[i].ij, p->s[i].j)) {
			warning("dispcal: inverting Jacobian failed (2)");
			*((char *)0) = 55;
		}
		/* Compute expected delta XYZ using new Jacobian */
		icmMulBy3x3(p->s[i].pdXYZ, p->s[i].j, p->s[i].pdrgb);

		p->s[i].de = b * os[j+1].de + (1.0 - b) * os[j].de;
		p->s[i].dc = b * os[j+1].dc + (1.0 - b) * os[j].dc;
	}

	free(os);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef NEVER
/* Do a linear interp of the ramdac */
static void interp_ramdac(double cal[256][3], double drgb[3], double rgb[3]) {
	int i, j;
	int gres = 256;
	double w;

	/* For r,g & b */
	for (j = 0; j < 3; j++) {
		int mi, gres_1 = gres-1;
		double t, vv = rgb[j];
		t = gres * vv;
		mi = (int)floor(t);			/* Grid coordinate */
		if (mi < 0)					/* Limit to valid cube base index range */
			mi = 0;
		else if (mi >= gres_1)
			mi = gres_1-1;
		w = t - (double)mi;	 		/* 1.0 - weight */

		drgb[j] = (1.0 - w) * cal[mi][j] + w * cal[mi+1][j];
	}
}
#endif	/* NEVER */

/* Given an XYZ, compute the color temperature and the delta E to the locus */
static double comp_ct(
	double *de,		/* If non-NULL, return CIEDE2000 to locus */
	double lxyz[3],	/* If non-NULL, return normalised XYZ on locus */
	int plank,		/* NZ if Plankian locus, 0 if Daylight locus */
	int dovct,		/* NZ if visual match, 0 if traditional correlation */
	double xyz[3]	/* Color to match */
) {
	double ct_xyz[3];		/* XYZ on locus */
	double nxyz[3];			/* Normalised input color */
	double ct, ctde;		/* Color temperature & delta E to Black Body locus */
	icmXYZNumber wN;

	if ((ct = icx_XYZ2ill_ct(ct_xyz, plank != 0 ? icxIT_Ptemp : icxIT_Dtemp,
	                         icxOT_CIE_1931_2, NULL, xyz, NULL, dovct)) < 0)
		error ("Got bad color temperature conversion\n");

	if (de != NULL) {
		icmAry2XYZ(wN, ct_xyz);
		icmAry2Ary(nxyz, xyz);
		nxyz[0] /= xyz[1];
		nxyz[2] /= xyz[1];
		nxyz[1] /= xyz[1];
		ctde = icmXYZCIE2K(&wN, nxyz, ct_xyz);
		*de = ctde;
	}
	if (lxyz != NULL) {
		icmAry2Ary(lxyz, ct_xyz);
	}
	return ct;
}
	
/* =================================================================== */

void usage(char *diag, ...) {
	disppath **dp;
	icoms *icom;

	fprintf(stderr,"Calibrate a Display, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	if (setup_spyd2() == 2)
		fprintf(stderr,"WARNING: This file contains a proprietary firmware image, and may not be freely distributed !\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: dispcal [options] outfile\n");
	fprintf(stderr," -v [n]               Verbose mode\n");
#if defined(UNIX) && !defined(__APPLE__)
	fprintf(stderr," -display displayname Choose X11 display name\n");
	fprintf(stderr," -d n[,m]             Choose the display n from the following list (default 1)\n");
	fprintf(stderr,"                      Optionally choose different display m for VideoLUT access\n"); 
#else
	fprintf(stderr," -d n                 Choose the display from the following list (default 1)\n");
#endif
//	fprintf(stderr," -d fake              Use a fake display device for testing, fake.icm if present\n");
	dp = get_displays();
	if (dp == NULL || dp[0] == NULL)
		fprintf(stderr,"    ** No displays found **\n");
	else {
		int i;
		for (i = 0; ; i++) {
			if (dp[i] == NULL)
				break;
			fprintf(stderr,"    %d = '%s'\n",i+1,dp[i]->description);
		}
	}
	free_disppaths(dp);
	fprintf(stderr," -c listno            Set communication port from the following list (default %d)\n",COMPORT);
	if ((icom = new_icoms()) != NULL) {
		icompath **paths;
		if ((paths = icom->get_paths(icom)) != NULL) {
			int i;
			for (i = 0; ; i++) {
				if (paths[i] == NULL)
					break;
				if (strlen(paths[i]->path) >= 8
				  && strcmp(paths[i]->path+strlen(paths[i]->path)-8, "Spyder2)") == 0
				  && setup_spyd2() == 0)
					fprintf(stderr,"    %d = '%s' !! Disabled - no firmware !!\n",i+1,paths[i]->path);
				else
					fprintf(stderr,"    %d = '%s'\n",i+1,paths[i]->path);
			}
		} else
			fprintf(stderr,"    ** No ports found **\n");
		icom->del(icom);
	}
	fprintf(stderr," -r                   Report on the calibrated display then exit\n");
	fprintf(stderr," -R                   Report on the uncalibrated display then exit\n");
	fprintf(stderr," -m                   Skip adjustment of the monitor controls\n");
	fprintf(stderr," -u [profile.icm]     Update previous calibration [update profile with new calibration]\n");
	fprintf(stderr," -q [lmh]             Quality - Low, Medium (def), High\n");
	fprintf(stderr," -y c|l               Display type, c = CRT, l = LCD\n");
	fprintf(stderr," -t [temp]            White Daylight locus target, optional target temperaturee in deg. K (deflt.)\n");
	fprintf(stderr," -T [temp]            White Black Body locus target, optional target temperaturee in deg. K\n");
	fprintf(stderr," -w x,y        	      Set the target white point as chromaticity coordinates\n");
#ifdef NEVER	/* Not worth confusing people about this */
	fprintf(stderr," -o                   Show CCT/CDT rather than VCT/VDT during native white point adjustment\n");
#endif
	fprintf(stderr," -b bright            Set the target brightness in cd/m^2\n");
	fprintf(stderr," -g gamma             Set the target response curve gamma (Def. %3.1f)\n",DEF_GAMMA);
	fprintf(stderr,"                      Use \"-gl\" for L*a*b* curve\n");
	fprintf(stderr,"                      Use \"-gs\" for sRGB curve\n");
	fprintf(stderr," -k factor            Amount to try and correct black point. Default 1.0, LCD default 0.0\n");
	fprintf(stderr," -e [n]               Run n verify passes on final curves\n");
	fprintf(stderr," -E                   Run only verify pass on installed calibration curves\n");
	fprintf(stderr," -p ho,vo,ss          Position test window and scale it\n");
	fprintf(stderr,"                      ho,vi: 0.0 = left/top, 0.5 = center, 1.0 = right/bottom etc.\n");
	fprintf(stderr,"                      ss: 0.5 = half, 1.0 = normal, 2.0 = double etc.\n");
#if defined(UNIX) && !defined(__APPLE__)
	fprintf(stderr," -n                   Don't set override redirect on test window\n");
#endif
	fprintf(stderr," -K                   Run instrument calibration first (used rarely)\n");
	fprintf(stderr," -N                   Disable auto calibration of instrument\n");
	fprintf(stderr," -H                   Use high resolution spectrum mode (if available)\n");
	fprintf(stderr," -D [level]           Print debug diagnostics to stderr\n");
	fprintf(stderr," outfile              Base name for output .cal file\n");
	exit(1);
	}

int main(int argc, char *argv[])
{
	int i, j, k;
	int fa, nfa, mfa;					/* current argument we're looking at */
	disppath *disp = NULL;				/* Display being used */
	double patsize = 100.0;				/* size of displayed color test patch */
	double patscale = 1.0;				/* scale factor for test patch size */
	double ho = 0.0, vo = 0.0;			/* Test window offsets, -1.0 to 1.0 */
	int verb = 0;
	int debug = 0;
	int fake = 0;						/* Use the fake device for testing */
	int override = 1;					/* Override redirect on X11 */
	int docalib = 0;					/* Do a manual instrument calibration */
	int doreport = 0;					/* 1 = Report the current uncalibrated display response */
										/* 2 = Report the current calibrated display response */
	int docontrols = 1;					/* Do adjustment of the display controls */
	int doupdate = 0;				    /* Do an update rather than a fresh calbration */
	char iccname[MAXNAMEL+1] = { 0 };	/* Optional icc profile to update */
	int comport = COMPORT;				/* COM port used */
	instType itype = instUnknown;		/* Default target instrument - none */
	int dtype = 0;						/* Display kind, 0 = default, 1 = CRT, 2 = LCD */
	int nocal = 0;						/* Disable auto calibration */
	int highres = 0;					/* Use high res mode if available */
	double temp = 0.0;					/* Color temperature (0 = native) */
	int planckian = 0;					/* 0 = Daylight, 1 = Planckian color locus */
	int dovct = 1;						/* Show VXT rather than CXT for adjusting white point */
	double wpx = 0.0, wpy = 0.0;		/* White point xy (native) */
	double tbright = 0.0;				/* Target brightness ( 0.0 == max)  */
	double gamma = DEF_GAMMA;			/* Popular Gamma target */
	double bkcorrect = -1.0;			/* Level of black point correction */ 
	int quality = -99;					/* Quality level, -2 = v, -1 = l, 0 = m, 1 = h, 2 = u */
	int isteps = 22;					/* Initial measurement steps/3 (medium) */
	int rsteps = 64;					/* Refinement measurement steps (medium) */
	double errthr = 1.5;				/* Error threshold for refinement steps (medium) */
	int mxits = 3;						/* maximum iterations (medium) */
	int verify = 0;						/* Do a verify after last refinement, 2 = do only verify. */
	int nver = 0;						/* Number of verify passes after refinement */
	char outname[MAXNAMEL+1] = { 0 };	/* Output cgats file base name */
	int spectral = 0;					/* Want spectral data from instrument */
	disprd *dr = NULL;					/* Display patch read object */
	cctx x;								/* Context for calibration solution */
	csamp asgrey;						/* Main calibration loop test points */
	int it;								/* verify & refine iteration */
	int rv;
	int fitord = 30;					/* More seems to make curves smoother */
	int donat;							/* Set native device colors ? (else via hw lut) */
	int errc;							/* Return value from new_disprd() */

	set_exe_path(argv[0]);				/* Set global exe_path and error_program */
	setup_spyd2();						/* Load firware if available */

	x.gammat = 0 ;						/* Default gamma type */
	x.tgamma = 0.0;						/* Default technical gamma */

#ifdef DEBUG_OFFSET
	ho = 0.8;
	vo = -0.8;
#endif

#if defined(DEBUG) || defined(DEBUG_OFFSET) || defined(DEBUG_PLOT)
	printf("!!!!!! Debug turned on !!!!!!\n");
#endif

	if (argc <= 1)
		usage("Too few arguments");

	/* Process the arguments */
	mfa = 1;        /* Minimum final arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')		/* Look for any flags */
			{
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1+mfa) < argc) {
					if (argv[fa+1][0] != '-')
						{
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?' || argv[fa][1] == '-') {
				usage("Usage requested");

			} else if (argv[fa][1] == 'v') {
				verb = 1;
				if (na != NULL && na[0] >= '0' && na[0] <= '9') {
					verb = atoi(na);
					fa = nfa;
				}

			/* Display number */
			} else if (argv[fa][1] == 'd') {
#if defined(UNIX) && !defined(__APPLE__)
				int ix, iv;

				if (strcmp(&argv[fa][2], "isplay") == 0 || strcmp(&argv[fa][2], "ISPLAY") == 0) {
					if (++fa >= argc || argv[fa][0] == '-') usage("Parameter expected following -display");
					setenv("DISPLAY", argv[fa], 1);
				} else {
					if (na == NULL) usage("Parameter expected following -d");
					fa = nfa;
					if (strcmp(na,"fake") == 0) {
						fake = 1;
					} else {
						if (sscanf(na, "%d,%d",&ix,&iv) != 2) {
							ix = atoi(na);
							iv = 0;
						}
						if ((disp = get_a_display(ix-1)) == NULL)
							usage("-d parameter %d out of range",ix);
						if (iv > 0)
							disp->rscreen = iv-1;
					}
				}
#else
				int ix;
				if (na == NULL) usage("Parameter expected following -d");
				fa = nfa;
				if (strcmp(na,"fake") == 0) {
					fake = 1;
				} else {
					ix = atoi(na);
					if ((disp = get_a_display(ix-1)) == NULL)
						usage("-d parameter %d out of range",ix);
				}
#endif

			} else if (argv[fa][1] == 'K') {
				docalib = 1;

			} else if (argv[fa][1] == 'N') {
				nocal = 1;

			/* High res mode */
			} else if (argv[fa][1] == 'H') {
				highres = 1;

			} else if (argv[fa][1] == 'D') {
				debug = 1;
				if (na != NULL && na[0] >= '0' && na[0] <= '9') {
					debug = atoi(na);
					fa = nfa;
				}

			} else if (argv[fa][1] == 'k') {
				fa = nfa;
				if (na == NULL) usage("Paramater expected following -k");
				bkcorrect = atof(na);
				if (bkcorrect < 0.0 || bkcorrect > 1.0) usage ("-k parameter must be between 0.0 and 1.0");
			} else if (argv[fa][1] == 'e') {
				verify = 1;
				nver = 1;
				if (na != NULL && na[0] >= '0' && na[0] <= '9') {
					nver = atoi(na);
					fa = nfa;
				}

			} else if (argv[fa][1] == 'E') {
				verify = 2;
				mfa = 0;

#if defined(UNIX) && !defined(__APPLE__)
			} else if (argv[fa][1] == 'n') {
				override = 0;
#endif /* UNIX */
			/* COM port  */
			} else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage("Paramater expected following -c");
				comport = atoi(na);
				if (comport < 1 || comport > 50) usage("-c parameter %d out of range",comport);

			} else if (argv[fa][1] == 'r' || argv[fa][1] == 'R') {
				if (argv[fa][1] == 'R')
					doreport = 1;		/* raw */
				else
					doreport = 2;		/* Calibrated */
				mfa = 0;

			} else if (argv[fa][1] == 'm' || argv[fa][1] == 'M') {
				docontrols = 0;

			} else if (argv[fa][1] == 'u' || argv[fa][1] == 'U') {
				doupdate = 1;
				docontrols = 0;

				if (na != NULL) {	/* Found an optional icc profile */
					fa = nfa;
					strncpy(iccname,na,MAXNAMEL); iccname[MAXNAMEL] = '\000';
				}

			/* Quality */
			} else if (argv[fa][1] == 'q' || argv[fa][1] == 'Q') {
				fa = nfa;
				if (na == NULL) usage("Parameter expected following -q");
    			switch (na[0]) {
					case 'L':			/* Test value */
						quality = -3;
						break;
					case 'v':
						quality = -2;
						break;
					case 'l':
						quality = -1;
						break;
					case 'm':
					case 'M':
						quality = 0;
						break;
					case 'h':
					case 'H':
						quality = 1;
						break;
					case 'u':
					case 'U':
						quality = 2;
						break;
					default:
						usage("-q parameter '%c' not recognised",na[0]);
				}

			/* Display type */
			} else if (argv[fa][1] == 'y' || argv[fa][1] == 'Y') {
				fa = nfa;
				if (na == NULL) usage("Parameter expected after -y");
				if (na[0] == 'c' || na[0] == 'C')
					dtype = 1;
				else if (na[0] == 'l' || na[0] == 'L')
					dtype = 2;
				else
					usage("-y parameter '%c' not recognised",na[0]);

			/* Daylight color temperature */
			} else if (argv[fa][1] == 't' || argv[fa][1] == 'T') {
				if (argv[fa][1] == 'T')
					planckian = 1;
				else
					planckian = 0;
				if (na != NULL) {
					fa = nfa;
					temp = atof(na);
					if (temp < 1000.0 || temp > 15000.0) usage("-%c parameter %f out of range",argv[fa][1], temp);
				}

			/* White point as x, y */
			} else if (argv[fa][1] == 'w' || argv[fa][1] == 'W') {
				fa = nfa;
				if (na == NULL) usage("Parameter expected after -w");
				if (sscanf(na, " %lf,%lf ", &wpx, &wpy) != 2)
					usage("-w parameter '%s' not recognised",na);

			/* Show CXT rather than VXT when adjusting native white point */
			} else if (argv[fa][1] == 'o' || argv[fa][1] == 'O') {
				dovct = 0;

			/* Brightness */
			} else if (argv[fa][1] == 'b' || argv[fa][1] == 'B') {
				fa = nfa;
				if (na == NULL) usage("Parameter expected after -b");
				tbright = atof(na);
				if (tbright <= 0.0 || tbright > 100000.0) usage("-b parameter %f out of range",tbright);

			/* Gamma */
			} else if (argv[fa][1] == 'g' || argv[fa][1] == 'G') {
				fa = nfa;
				if (na == NULL) usage("Parameter expected after -g");
				if ((na[0] == 'l' || na[0] == 'L') && na[1] == '\000')
					x.gammat = 1;
				else if ((na[0] == 's' || na[0] == 'S') && na[1] == '\000')
					x.gammat = 2;
				else {
					gamma = atof(na);
					if (gamma <= 0.0 || gamma > 10.0) usage("-g parameter %f out of range",gamma);
					x.gammat = 0;
				}
			/* Test patch offset and size */
			} else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				fa = nfa;
				if (na == NULL) usage("Parameter expected after -p");
				if (sscanf(na, " %lf,%lf,%lf ", &ho, &vo, &patscale) != 3)
					usage("-p parameter '%s' not recognised",na);
				if (ho < 0.0 || ho > 1.0
				 || vo < 0.0 || vo > 1.0
				 || patscale <= 0.0 || patscale > 50.0)
					usage("-p parameters %f %f %f out of range",ho,vo,patscale);
				ho = 2.0 * ho - 1.0;
				vo = 2.0 * vo - 1.0;

			} else 
				usage("Flag '-%c' not recognised",argv[fa][1]);
		} else
			break;
	}

	if (bkcorrect < 0.0) {
		if (dtype == 2)		/* LCD */
			bkcorrect = 0.0;	/* Default to no black point correction for LCD */
		else
			bkcorrect = 1.0;	/* Default to full black point correction for other displays */
	}

	patsize *= patscale;

	if (!fake && disp == NULL && (disp = get_a_display(0)) == NULL) {
		error("Unable to open the display");
	}

	if (docalib) {
		if ((rv = disprd_calibration(itype, comport, dtype, nocal, disp, override, patsize, ho, vo, verb, debug)) != 0) {
			error("docalibration failed with return value %d\n",rv);
		}
	}

	if (verify != 2 && doreport == 0) {
		/* Get the file name argument */
		if (fa >= argc || argv[fa][0] == '-') usage("Output filname parameter not found");
		strncpy(outname,argv[fa++],MAXNAMEL-4); outname[MAXNAMEL-4] = '\000';
		strcat(outname,".cal");
	}

	if (verify == 2) {
		if (doupdate)
			warning("Update flag ignored because we're doing a verify only");
		doupdate = 0;
		docontrols = 0;
	}

	if (doreport != 0) {
		if (verify == 2)
			warning("Verify flag ignored because we're doing a report only");
		verify = 0;
	}

	/* Get ready to do some readings */
	donat = 1;		/* Normally calibrate against native response */
	if (verify == 2 || doreport == 2)
		donat = 0;	/* But measure calibrated response of verify or report calibrated */ 
	if ((dr = new_disprd(&errc, itype, fake ? -99 : comport, dtype, nocal, highres, donat,
	                     NULL, disp, override, patsize, ho, vo,
	                     spectral, verb, VERBOUT, debug, "fake.icm")) == NULL)
		error("dispread failed with '%s'\n",disprd_err(errc));

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if (doreport) {
		col tcols[3] = {	/* Base set of test colors */
			{ 0.0, 0.0, 0.0 },
			{ 0.5, 0.5, 0.5 },
			{ 1.0, 1.0, 1.0 }
		};
		double cct, cct_de;
		double cdt, cdt_de;
		double vct, vct_de;
		double vdt, vdt_de;
		double gamma, w[3], wp[2];

		if ((rv = dr->read(dr, tcols, 3, 1, 3, 1, 0)) != 0) {
			dr->del(dr);
			error("display read failed with '%s'\n",disprd_err(rv));
		} 

//printf("~1 Got black = %f, half = %f, white = %f\n",tcols[0].aXYZ[1],tcols[1].aXYZ[1],tcols[2].aXYZ[1]);
		/* Normalised XYZ white point */
		w[0] = tcols[2].aXYZ[0]/tcols[2].aXYZ[1];
		w[1] = tcols[2].aXYZ[1]/tcols[2].aXYZ[1];
		w[2] = tcols[2].aXYZ[2]/tcols[2].aXYZ[1];

		/* White point chromaticity coordinates */
		wp[0] = w[0]/(w[0] + w[1] + w[2]);
		wp[1] = w[1]/(w[0] + w[1] + w[2]);

		cct = comp_ct(&cct_de, NULL, 1, 0, w);	/* Compute CCT */
		cdt = comp_ct(&cdt_de, NULL, 0, 0, w);	/* Compute CDT */
		vct = comp_ct(&vct_de, NULL, 1, 1, w);	/* Compute VCT */
		vdt = comp_ct(&vdt_de, NULL, 0, 1, w);	/* Compute VDT */

		/* Approximate Gamma - use the gross curve shape for robustness */
		gamma = pop_gamma(tcols[0].aXYZ[1], tcols[1].aXYZ[1], tcols[2].aXYZ[1]);

		if (doreport == 2)
			printf("Current calibration response:\n");
		else
			printf("Uncalibrated response:\n");
		printf("Black level = %.2f cd/m^2\n",tcols[0].aXYZ[1]);
		printf("White level = %.2f cd/m^2\n",tcols[2].aXYZ[1]);
		printf("Aprox. gamma = %.2f\n",gamma);
		printf("Contrast ratio = %.0f:1\n",tcols[2].aXYZ[1]/tcols[0].aXYZ[1]);
		printf("White chromaticity coordinates %.4f, %.4f\n",wp[0],wp[1]);
		printf("White    Correlated Color Temperature = %.0fK, DE to locus = %4.1f\n",cct,cct_de);
		printf("White Correlated Daylight Temperature = %.0fK, DE to locus = %4.1f\n",cdt,cdt_de);
		printf("White        Visual Color Temperature = %.0fK, DE to locus = %4.1f\n",vct,vct_de);
		printf("White     Visual Daylight Temperature = %.0fK, DE to locus = %4.1f\n",vdt,vdt_de);
		dr->del(dr);
		exit(0);
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* If we're updating, retrieve the previously used settings, */
	/* and device model */
	if (doupdate) {
		cgats *icg;						/* output cgats structure */
		int nsamp;
		mcvco *rdv[3];					/* Scattered data for ramdac curves */
		int fi;							/* Field index */
		int si[4];						/* Set fields */

		if (verb)
			printf("Updating previous calibration\n");

		icg = new_cgats();				/* Create a CGATS structure */
		icg->add_other(icg, "CAL"); 	/* our special type is Calibration file */
		
		if (icg->read_name(icg, outname)) {
			dr->del(dr);
			error("Can't update '%s' - read error : %s",outname, icg->err);
		}

		if (icg->ntables == 0
		 || icg->t[0].tt != tt_other || icg->t[0].oi != 0
		 || icg->t[1].tt != tt_other || icg->t[1].oi != 0) {
			dr->del(dr);
			error("Can't update '%s' - wrong type of file",outname);
		}
		if (icg->ntables < 2) {
			dr->del(dr);
			error("Can't update '%s' - there aren't two tables",outname);
		}

//printf("~1 reading previous cal, got 2 tables\n");

		/* Read in the setup, user and model values */

		dtype = 0;		/* Override any user setting */
		if ((fi = icg->find_kword(icg, 0, "DEVICE_TYPE")) >= 0) {
			if (strcmp(icg->t[0].kdata[fi], "CRT") == 0)
				dtype = 1;
			else if (strcmp(icg->t[0].kdata[fi], "LCD") == 0)
				dtype = 2;
			else {
				dr->del(dr);
				error ("Can't update '%s' - field 'DISPLAY_TYPE' has unrecognised value '%s'",
				        outname,icg->t[0].kdata[fi]);
			}
		}
//printf("~1 dealt with device type\n");
	
		if ((fi = icg->find_kword(icg, 0, "TARGET_WHITE_XYZ")) < 0) {
			dr->del(dr);
			error ("Can't update '%s' - can't find field 'TARGET_WHITE_XYZ'",outname);
		}
		if (sscanf(icg->t[0].kdata[fi], "%lf %lf %lf", &x.twh[0], &x.twh[1], &x.twh[2]) != 3) {
			dr->del(dr);
			error ("Can't update '%s' - reading field 'TARGET_WHITE_XYZ' failed",outname);
		}
		if ((fi = icg->find_kword(icg, 0, "NATIVE_TARGET_WHITE")) >= 0) {
			wpx = wpy = 0.0;
			temp = 0.0;
			tbright = 0.0;
		} else {
			wpx = wpy = 0.0001;
			temp = 1.0;
			tbright = 1.0;
		}
//printf("~1 dealt with target white\n");

		if ((fi = icg->find_kword(icg, 0, "TARGET_GAMMA")) < 0) {
			dr->del(dr);
			error ("Can't update '%s' - can't find field 'TARGET_GAMMA'",outname);
		}
		if (strcmp(icg->t[0].kdata[fi], "L_STAR") == 0)
			x.gammat = 1;
		else if (strcmp(icg->t[0].kdata[fi], "sRGB") == 0)
			x.gammat = 2;
		else {
			x.gammat = 0;
			gamma = atof(icg->t[0].kdata[fi]);
			if (gamma < 0.1 || gamma > 5.0) {
				dr->del(dr);
				error ("Can't update '%s' - field 'TARGET_GAMMA' has bad value %f",outname,gamma);
			}
		}
//printf("~1 dealt with target gamma\n");

		if ((fi = icg->find_kword(icg, 0, "BLACK_POINT_CORRECTION")) < 0) {
			dr->del(dr);
			error ("Can't update '%s' - can't find field 'BLACK_POINT_CORRECTION'",outname);
		}
		bkcorrect = atof(icg->t[0].kdata[fi]);
		if (bkcorrect < 0.0 || bkcorrect > 1.0) {
			dr->del(dr);
			error ("Can't update '%s' - field 'BLACK_POINT_CORRECTION' has bad value %f",outname,bkcorrect);
		}
//printf("~1 dealt with black point correction\n");

		if ((fi = icg->find_kword(icg, 0, "QUALITY")) < 0) {
			dr->del(dr);
			error ("Can't update '%s' - can't find field 'QUALITY'",outname);
		}
		if (quality < -50) {	/* User hasn't overridden quality */
			if (strcmp(icg->t[0].kdata[fi], "ultra low") == 0)
				quality = -3;
			else if (strcmp(icg->t[0].kdata[fi], "very low") == 0)
				quality = -2;
			else if (strcmp(icg->t[0].kdata[fi], "low") == 0)
				quality = -1;
			else if (strcmp(icg->t[0].kdata[fi], "medium") == 0)
				quality = 0;
			else if (strcmp(icg->t[0].kdata[fi], "high") == 0)
				quality = 1;
			else if (strcmp(icg->t[0].kdata[fi], "ultra high") == 0)
				quality = 2;
			else {
				dr->del(dr);
				error ("Can't update '%s' - field 'QUALITY' has unrecognised value '%s'",
				        outname,icg->t[0].kdata[fi]);
			}
		}
//printf("~1 dealt with quality\n");

		/* Read in the last set of calibration curves used */
		if ((nsamp = icg->t[0].nsets) < 2) {
			dr->del(dr);
			error("Can't update '%s' - %d not enough data points in calibration curves",
			        outname,nsamp);
		}
//printf("~1 got %d points in calibration curves\n",nsamp);

		for (k = 0; k < 3; k++) {
			if ((x.rdac[k] = new_mcv()) == NULL) {
				dr->del(dr);
				error("new_mcv x.rdac[%d] failed",k);
			}
			if ((rdv[k] = malloc(sizeof(mcvco) * nsamp)) == NULL) {
				dr->del(dr);
				error ("Malloc of scattered data points failed");
			}
		}
//printf("~1 allocated calibration curve objects\n");

		/* Read the current calibration curve points (usually 256 of them) */
		for (k = 0; k < 4; k++) {
			char *fnames[4] = { "RGB_I", "RGB_R", "RGB_G", "RGB_B" };

			if ((si[k] = icg->find_field(icg, 0, fnames[k])) < 0) {
				dr->del(dr);
				error ("Can't updata '%s' - can't find field '%s'",outname,fnames[k]);
			}
			if (icg->t[0].ftype[si[k]] != r_t) {
				dr->del(dr);
				error ("Can't updata '%s' - field '%s' is wrong type",outname,fnames[k]);
			}
		}
//printf("~1 Found calibration curve fields\n");

		for (i = 0; i < nsamp; i++) {
			rdv[0][i].p =
			rdv[1][i].p =
			rdv[2][i].p =
			                *((double *)icg->t[0].fdata[i][si[0]]);
			for (k = 0; k < 3; k++) {		/* RGB */
				rdv[k][i].v = *((double *)icg->t[0].fdata[i][si[k + 1]]);
			}
			rdv[0][i].w = rdv[1][i].w = rdv[2][i].w = 1.0;
		}
//printf("~1 Read calibration curve data points\n");
		for (k = 0; k < 3; k++) {
			x.rdac[k]->fit(x.rdac[k], 0, fitord, rdv[k], nsamp, 1.0);
			free (rdv[k]);
		}
//printf("~1 Fitted calibration curves\n");

		/* Read in the per channel forward model curves */
		for (k = 0; k < 3; k++) {
			char *fnames[3] = { "R_P", "G_P", "B_P" };
			double *pp;

//printf("~1 Reading device curve channel %d\n",k);
			if ((si[k] = icg->find_field(icg, 1, fnames[k])) < 0) {
				dr->del(dr);
				error ("Can't updata '%s' - can't find field '%s'",outname,fnames[k]);
			}
			if (icg->t[1].ftype[si[k]] != r_t) {
				dr->del(dr);
				error ("Can't updata '%s' - field '%s' is wrong type",outname,fnames[k]);
			}
			/* Create the model curves */
			if ((pp = (double *)malloc(icg->t[1].nsets * sizeof(double))) == NULL) {
				dr->del(dr);
				error ("Malloc of device curve parameters");
			}
			for (i = 0; i < icg->t[1].nsets; i++)
				pp[i] = *((double *)icg->t[1].fdata[i][si[k]]);

			if ((x.dcvs[k] = new_mcv_p(pp, icg->t[1].nsets)) == NULL) {
				dr->del(dr);
				error("new_mcv x.dcvs[%d] failed",k);
			}
			free(pp);
		}
		
		icg->del(icg);
//printf("~1 read in previous settings and device model\n");
	}

	/* Be nice - check we can read the iccprofile before calibrating the display */
	if (doupdate && iccname[0] != '\000') {
		icmFile *ic_fp;
		icc *icco;

		if ((icco = new_icc()) == NULL) {
			dr->del(dr);
			error ("Creation of ICC object to read profile '%s' failed",iccname);
		}

		/* Open up the profile for reading */
		if ((ic_fp = new_icmFileStd_name(iccname,"r")) == NULL) {
			dr->del(dr);
			error ("Can't open file '%s'",iccname);
		}

		/* Read header etc. */
		if ((rv = icco->read(icco,ic_fp,0)) != 0) {
			dr->del(dr);
			error ("Reading profile '%s' failed with %d, %s",iccname, rv,icco->err);
		}

		ic_fp->del(ic_fp);

		if (icco->find_tag(icco, icSigVideoCardGammaTag) != 0) {
			dr->del(dr);
			error("Can't find VideoCardGamma tag in file '%s': %d, %s",
			      iccname, icco->errc,icco->err);
		}
		icco->del(icco);
	}
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	/* Convert quality level to iterations etc. */
	/* Note that final tollerance is often double the */
	/* final errth, because one more corrections is always */
	/* performed after the last reading. */
    switch (quality) {
		case -3:				/* Test value */
			isteps = 3;
			rsteps = 9;
			mxits = 1;
			errthr = 2.5;
			break;
		case -2:				/* Very low */
			isteps = 10;
			rsteps = 16;
			errthr = 2.0;
			if (doupdate)
				mxits = 1;
			else
				mxits = 1;
			break;
		case -1:				/* Low */
			isteps = 12;
			rsteps = 32;
			errthr = 1.2;
			if (doupdate)
				mxits = 1;
			else
				mxits = 2;
			break;
		default:
		case 0:					/* Medum */
			quality = 0;		/* In case it wasn't set */
			isteps = 16;
			rsteps = 64;
			errthr = 0.8;
			if (doupdate)
				mxits = 1;
			else
				mxits = 3;
			break;
		case 1:					/* High */
			isteps = 20;
			rsteps = 96;
			errthr = 0.5;
			if (doupdate)
				mxits = 1;
			else
				mxits = 4;
			break;
		case 2:					/* Ultra */
			isteps = 24;
			rsteps = 128;
			errthr = 0.3;
			if (doupdate)
				mxits = 1;
			else
				mxits = 5;
			break;
	}

	/* Set native white target flag in cctx so that other things can play the game.. */
	if (wpx == 0.0 && wpy == 0.0 && temp == 0.0 && tbright == 0.0)
		x.nat = 1;
	else
		x.nat = 0;

	if (verb) {
		if (dtype == 1)
			printf("Display type is CRT\n");
		else if (dtype == 2)
			printf("Display type is LCD\n");

		if (doupdate) {
			if (x.nat)
				printf("Target white = native white point & brightness\n");
			else
				printf("Target white = XYZ %f %f %f\n",
				       x.twh[0], x.twh[1], x.twh[2]);
		} else {
			if (wpx > 0.0 || wpy > 0.0)
				printf("Target white = xy %f %f\n",wpx,wpy);
			else if (temp > 0.0) {
				if (planckian)
					printf("Target white = %f degrees kelvin Planckian (black body) spectrum\n",temp);
				else
					printf("Target white = %f degrees kelvin Daylight spectrum\n",temp);
			} else
				printf("Target white = native white point\n");
	
			if (tbright > 0.0)
				printf("Target brightness = %f cd/m^2\n",tbright);
			else
				printf("Target brightness = native brightness\n");
		}

		if (x.gammat == 0)
			printf("Target gamma = %f\n",gamma);
		else if (x.gammat == 1)
			printf("Target gamma = L* curve\n");
		else if (x.gammat == 2)
			printf("Target gamma = sRGB curve\n");
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Go through the procedure of adjusting monitor controls */
	if (docontrols) {
		int rgbch = 0;			/* Got RBG Yxy ? */
		double rgbXYZ[3][3];	/* The RGB XYZ */

		/* Until the user is done */
		printf("\nDisplay adjustment menu:");
		for (;;) {
			int c;

			/* Print the menue of adjustments */
			printf("\nPress 1 .. 7\n");
			printf("1) Black level (CRT: Brightness)\n");
			printf("2) White point (Color temperature, R,G,B, Gain)\n");
			printf("3) White level (CRT: Contrast, LCD: Brightness)\n");
			printf("4) Black point (R,G,B, Offset)\n");
			printf("5) Check all\n");
			printf("6) Continue on to calibration\n");
			printf("7) Exit\n");

			empty_con_chars();
			c = next_con_char();

			/* Black level adjustment */
			/* Due to the possibility of the channel offsets not being even, */
			/* we use the largest of the XYZ values after they have been */
			/* scaled to be even acording to the white XYZ balance. */
			/* It's safer to set the black level a bit low, and then the */
			/* calibration curves can bump the low ones up. */
			if (c == '1') {
				col tcols[3] = {	/* Base set of test colors */
					{ 0.0, 0.0, 0.0 },
					{ 0.5, 0.5, 0.5 },		/* And 1% values */
					{ 1.0, 1.0, 1.0 }
				};
				int ff;
				double mgamma, tar1, dev1;
		
				printf("Doing some initial measurements\n");
				/* Do an initial set of readings to set 1% output mark */
				if ((rv = dr->read(dr, tcols, 3, 0, 0, 1, 0)) != 0) {
					dr->del(dr);
					error("display read failed with '%s'\n",disprd_err(rv));
				} 

				/* Popular Gamma - Gross curve shape */
				mgamma = pop_gamma(tcols[0].aXYZ[1], tcols[1].aXYZ[1], tcols[2].aXYZ[1]);

				dev1 = pow(0.01, 1.0/mgamma);
//printf("~1 device level for 1%% output = %f\n",dev1);
				tcols[1].r = tcols[1].g = tcols[1].b = dev1;
				tar1 = 0.01 * tcols[2].aXYZ[1];

				printf("\nAdjust CRT brightness to get target level. Press space when done.\n");
				printf("   Target %.2f\n",tar1);
				for (ff = 0;; ff ^= 1) {
					double dir;			/* Direction to adjust brightness */
					double sv[3], val1;				/* Scaled values */
					if ((rv = dr->read(dr, tcols+1, 1, 0, 0, 1, ' ')) != 0) {
						if (rv == 4)
							break;			/* User is done with this adjustment */
						dr->del(dr);
						error("display read failed with '%s'\n",disprd_err(rv));
					} 
					sv[0] = tcols[1].aXYZ[0] * tcols[2].aXYZ[1]/tcols[2].aXYZ[0];
					sv[1] = tcols[1].aXYZ[1];
					sv[2] = tcols[1].aXYZ[2] * tcols[2].aXYZ[2]/tcols[2].aXYZ[2];
//printf("~1 scaled readings = %f %f %f\n",sv[0],sv[1],sv[2]);
					val1 = sv[1];
					if (sv[0] > val1)
						val1 = sv[0];
					if (sv[2] > val1)
						val1 = sv[2];
					dir = tar1 - val1;
					if (fabs(dir) < 0.01)
						dir = 0.0;
					printf("\r%c Current %.2f  %c",
					       ff == 0 ? '/' : '\\',
					       val1,
					       dir < 0.0 ? '-' : dir > 0.0 ? '+' : '=');
					fflush(stdout);
					c = poll_con_char();
					if (c == ' ')
						break;
				}
				printf("\n");

			/* White point adjustment */
			} else if (c == '2') {
				int nat = 0;		/* NZ if using native white as target */
				col tcols[1] = {	/* Base set of test colors */
					{ 1.0, 1.0, 1.0 }
				};
				int ff;
				double tYxy[3];				/* Target white chromaticities */
				icmXYZNumber tXYZ;			/* Target white as XYZ */
				double tLab[3];				/* Target white as Lab or UCS */
				double tarw;				/* Target brightness */
				double Lab[3];				/* Last measured point Lab or UCS */
				double ct = 0.0, ct_de;		/* Color temperature & delta E to white locus */
		
				printf("Doing some initial measurements\n");

				if (rgbch == 0) {	/* Figure the RGB chromaticities */
					col ccols[3] = {
						{ 1.0, 0.0, 0.0 },
						{ 0.0, 1.0, 0.0 },
						{ 0.0, 0.0, 1.0 }
					};
					if ((rv = dr->read(dr, ccols, 3, 0, 0, 1, 0)) != 0) {
						dr->del(dr);
						error("display read failed with '%s'\n",disprd_err(rv));
					}
					for (i = 0; i < 3; i++)
						icmAry2Ary(rgbXYZ[i], ccols[i].aXYZ);
					rgbch = 1;
				}
				/* Do an initial set of readings to set full output mark */
				if ((rv = dr->read(dr, tcols, 1, 0, 0, 1, 0)) != 0) {
					dr->del(dr);
					error("display read failed with '%s'\n",disprd_err(rv));
				}
				
				/* Figure out the target white chromaticity */
				if (wpx > 0.0 || wpy > 0.0) {	/* xy coordinates */
					tYxy[0] = 1.0;
					tYxy[1] = wpx;
					tYxy[2] = wpy;
				} else if (temp > 0.0) {		/* Daylight color temperature */
					double XYZ[3];
					if (planckian)
						rv = icx_ill_sp2XYZ(XYZ, icxOT_default, NULL, icxIT_Ptemp, temp, NULL);
					else
						rv = icx_ill_sp2XYZ(XYZ, icxOT_default, NULL, icxIT_Dtemp, temp, NULL);
					if (rv != 0)
						error("Failed to compute XYZ of target color temperature %f\n",temp);
					icmXYZ2Yxy(tYxy, XYZ);
				} else {						/* Native white */
					icmXYZ2Yxy(tYxy, tcols[0].aXYZ);
					nat = 1;
				}

				/* Figure out the target white brightness */
				/* Note we're not taking the device gamut into account here */
				if (tbright > 0.0)			/* Given brightness */
					tarw = tbright;
				else						/* Native/maximum brightness */
					tarw = tcols[0].aXYZ[1];

				if (!nat) {	/* Target is a specified white */
					printf("\nAdjust R,G & B gain to get target x,y. Press space when done.\n");
					printf("   Target B %.2f, x %.4f, y %.4f\n",
					        tarw, tYxy[1],tYxy[2]);

				} else {	/* Target is native white */
					printf("\nAdjust R,G & B gain to desired white point. Press space when done.\n");
					/* Compute the CT and delta E to white locus of target */
					ct = comp_ct(&ct_de, NULL, planckian, dovct, tcols[0].aXYZ);
					printf("  Initial B %.2f, x %.4f, y %.4f, %c%cT %4.0fK DE %4.1f\n",
					        tarw, tYxy[1],tYxy[2],
				            dovct ? 'V' : 'C', planckian ? 'C' : 'D', ct,ct_de);
				}
				for (ff = 0;; ff ^= 1) {
					double dir;			/* Direction to adjust brightness */
					double Yxy[3];		/* Yxy of current reading */
					double rgbdir[3];	/* Direction to adjust RGB */
					double rgbxdir[3];	/* Biggest to move */
					double bdir, terr;
					int bx = 0;
					
					if ((rv = dr->read(dr, tcols, 1, 0, 0, 1, ' ')) != 0) {
						if (rv == 4)
							break;			/* User is done with this adjustment */
						dr->del(dr);
						error("display read failed with '%s'\n",disprd_err(rv));
					} 
					dir = tarw - tcols[0].aXYZ[1];
					if (fabs(dir) < 0.01)
						dir = 0.0;

					icmXYZ2Yxy(Yxy, tcols[0].aXYZ);

					if (!nat) {	/* Target is a specified white */
						/* Compute values we need for delta E and RGB direction */
						icmYxy2XYZ(tLab, tYxy);
						tLab[0] /= tLab[1];
						tLab[2] /= tLab[1];
						tLab[1] /= tLab[1];
						icmAry2XYZ(tXYZ, tLab);					/* Lab white reference */
						icmXYZ2Lab(&tXYZ, tLab, tLab);			/* Target Lab */
	
						icmAry2Ary(Lab, tcols[0].aXYZ);
						Lab[0] /= Lab[1];
						Lab[2] /= Lab[1];
						Lab[1] /= Lab[1];
						icmXYZ2Lab(&tXYZ, Lab, Lab);			/* Current Lab */

					} else {	/* Target is native white */
						double lxyz[3];	/* Locus XYZ */
						ct = comp_ct(&ct_de, lxyz, planckian, dovct, tcols[0].aXYZ);

						icmXYZ2Yxy(tYxy, lxyz);
						/* lxyz is already normalised */
						icmAry2XYZ(tXYZ, lxyz);					/* Lab white reference */
						if (dovct)
							icmXYZ2Lab(&tXYZ, tLab, lxyz);			/* Target Lab */
						else
							icmXYZ2UCS(&tXYZ, tLab, lxyz);			/* Target UCS */
	
						icmAry2Ary(Lab, tcols[0].aXYZ);
						Lab[0] /= Lab[1];
						Lab[2] /= Lab[1];
						Lab[1] /= Lab[1];
						if (dovct)
							icmXYZ2Lab(&tXYZ, Lab, Lab);			/* Current Lab */
						else
							icmXYZ2UCS(&tXYZ, Lab, Lab);			/* Current UCS */
					}

					/* Compute dot products */
					bdir = 0.0;
					for (i = 0; i < 3; i++) {
						double rgbLab[3];
						
						if (dovct)
							icmXYZ2Lab(&tXYZ, rgbLab, rgbXYZ[i]);
						else
							icmXYZ2UCS(&tXYZ, rgbLab, rgbXYZ[i]);
						rgbdir[i] = (tLab[1] - Lab[1]) * (rgbLab[1] - Lab[1])
						          + (tLab[2] - Lab[2]) * (rgbLab[2] - Lab[2]);
						rgbxdir[i] = 0.0;
						if (fabs(rgbdir[i]) > fabs(bdir)) {
							bdir = rgbdir[i];
							bx = i;
						}
					}

					/* See how close to the target we are */
					terr = sqrt((tLab[1] - Lab[1]) * (tLab[1] - Lab[1])
					          + (tLab[2] - Lab[2]) * (tLab[2] - Lab[2]));
					if (terr < 0.1)
						rgbdir[0] = rgbdir[1] = rgbdir[2] = 0.0;
					rgbxdir[bx] = rgbdir[bx];

	
					if (!nat) {
						printf("\r%c Current B %.2f, x %.4f, y %.4f  DE %4.1f  R%c%c G%c%c B%c%c ",
					       ff == 0 ? '/' : '\\',
					       tcols[0].aXYZ[1], Yxy[1], Yxy[2],
					       icmCIE2K(tLab, Lab),
					       rgbdir[0] < 0.0 ? '-' : rgbdir[0] > 0.0 ? '+' : '=',
					       rgbxdir[0] < 0.0 ? '-' : rgbxdir[0] > 0.0 ? '+' : ' ',
					       rgbdir[1] < 0.0 ? '-' : rgbdir[1] > 0.0 ? '+' : '=',
					       rgbxdir[1] < 0.0 ? '-' : rgbxdir[1] > 0.0 ? '+' : ' ',
					       rgbdir[2] < 0.0 ? '-' : rgbdir[2] > 0.0 ? '+' : '=',
					       rgbxdir[2] < 0.0 ? '-' : rgbxdir[2] > 0.0 ? '+' : ' ');
					} else {
						printf("\r%c Current B %.2f, x %.4f, y %.4f  %c%cT %4.0fK DE %4.1f  R%c%c G%c%c B%c%c ",
					       ff == 0 ? '/' : '\\',
					       tcols[0].aXYZ[1], Yxy[1], Yxy[2],
				           dovct ? 'V' : 'C', planckian ? 'C' : 'D', ct,ct_de,
					       rgbdir[0] < 0.0 ? '-' : rgbdir[0] > 0.0 ? '+' : '=',
					       rgbxdir[0] < 0.0 ? '-' : rgbxdir[0] > 0.0 ? '+' : ' ',
					       rgbdir[1] < 0.0 ? '-' : rgbdir[1] > 0.0 ? '+' : '=',
					       rgbxdir[1] < 0.0 ? '-' : rgbxdir[1] > 0.0 ? '+' : ' ',
					       rgbdir[2] < 0.0 ? '-' : rgbdir[2] > 0.0 ? '+' : '=',
					       rgbxdir[2] < 0.0 ? '-' : rgbxdir[2] > 0.0 ? '+' : ' ');
					}
					fflush(stdout);
					c = poll_con_char();
					if (c == ' ')
						break;
				}
				printf("\n");

			/* White level adjustment */
			} else if (c == '3') {
				col tcols[1] = {	/* Base set of test colors */
					{ 1.0, 1.0, 1.0 }
				};
				int ff;
				double tarw;
		
				printf("Doing some initial measurements\n");
				/* Do an initial set of readings to set full output mark */
				if ((rv = dr->read(dr, tcols, 1, 0, 0, 1, 0)) != 0) {
					dr->del(dr);
					error("display read failed with '%s'\n",disprd_err(rv));
				}
				
				/* Figure out the target white brightness */
				/* Note we're not taking the device gamut into account here */
				if (tbright > 0.0)			/* Given brightness */
					tarw = tbright;
				else						/* Native/maximum brightness */
					tarw = tcols[0].aXYZ[1];

				if (tbright > 0.0) {
					printf("\nAdjust CRT Contrast or LCD Brightness to get target level. Press space when done.\n");
					printf("   Target %.2f\n", tarw);
				} else {
					printf("\nAdjust CRT Contrast or LCD Brightness to desired level. Press space when done.\n");
					printf("  Initial %.2f\n", tarw);
				}
				for (ff = 0;; ff ^= 1) {
					double dir;			/* Direction to adjust brightness */
					
					if ((rv = dr->read(dr, tcols, 1, 0, 0, 1, ' ')) != 0) {
						if (rv == 4)
							break;			/* User is done with this adjustment */
						dr->del(dr);
						error("display read failed with '%s'\n",disprd_err(rv));
					} 
					dir = tarw - tcols[0].aXYZ[1];
					if (fabs(dir) < 0.01)
						dir = 0.0;

					if (tbright > 0.0)
						printf("\r%c Current %.2f  %c",
					       ff == 0 ? '/' : '\\',
					       tcols[0].aXYZ[1],
					       dir < 0.0 ? '-' : dir > 0.0 ? '+' : '=');
					else
						printf("\r%c Current %.2f   ",
					       ff == 0 ? '/' : '\\',
					       tcols[0].aXYZ[1]);
					fflush(stdout);
					c = poll_con_char();
					if (c == ' ')
						break;
				}
				printf("\n");

			/* Black point adjustment */
			} else if (c == '4') {
				col tcols[3] = {	/* Base set of test colors */
					{ 0.0, 0.0, 0.0 },
					{ 0.5, 0.5, 0.5 },		/* And 1% values */
					{ 1.0, 1.0, 1.0 }
				};
				int ff;
				double tYxy[3];				/* Target white chromaticities */
				icmXYZNumber tXYZ;			/* Target white as XYZ */
				double tLab[3];				/* Target white as Lab or UCS */
				double mgamma, tar1, dev1;
				double Lab[3];				/* Last measured point Lab */
		
				printf("Doing some initial measurements\n");

				if (rgbch == 0) {	/* Figure the RGB chromaticities */
					col ccols[3] = {
						{ 1.0, 0.0, 0.0 },
						{ 0.0, 1.0, 0.0 },
						{ 0.0, 0.0, 1.0 }
					};
					if ((rv = dr->read(dr, ccols, 3, 0, 0, 1, 0)) != 0) {
						dr->del(dr);
						error("display read failed with '%s'\n",disprd_err(rv));
					}
					for (i = 0; i < 3; i++)
						icmAry2Ary(rgbXYZ[i], ccols[i].aXYZ);
					rgbch = 1;
				}
				/* Do an initial set of readings to set 1% output mark */
				if ((rv = dr->read(dr, tcols, 3, 0, 0, 1, 0)) != 0) {
					dr->del(dr);
					error("display read failed with '%s'\n",disprd_err(rv));
				}
				
				/* Popular Gamma - Gross curve shape */
				mgamma = pop_gamma(tcols[0].aXYZ[1], tcols[1].aXYZ[1], tcols[2].aXYZ[1]);

				dev1 = pow(0.01, 1.0/mgamma);
				tcols[1].r = tcols[1].g = tcols[1].b = dev1;
				tar1 = 0.01 * tcols[2].aXYZ[1];

				/* Figure out the target white chromaticity */
				if (wpx > 0.0 || wpy > 0.0) {	/* xy coordinates */
					tYxy[0] = 1.0;
					tYxy[1] = wpx;
					tYxy[2] = wpy;
		
				} else if (temp > 0.0) {		/* Daylight color temperature */
					double XYZ[3];
					if (planckian)
						rv = icx_ill_sp2XYZ(XYZ, icxOT_default, NULL, icxIT_Ptemp, temp, NULL);
					else
						rv = icx_ill_sp2XYZ(XYZ, icxOT_default, NULL, icxIT_Dtemp, temp, NULL);
					if (rv != 0)
						error("Failed to compute XYZ of target color temperature %f\n",temp);
					icmXYZ2Yxy(tYxy, XYZ);
				} else {						/* Native white */
					icmXYZ2Yxy(tYxy, tcols[2].aXYZ);
				}

				printf("\nAdjust R,G & B offsets to get target x,y. Press space when done.\n");
				printf("   Target B %.2f, x %.4f, y %.4f\n", tar1, tYxy[1],tYxy[2]);
				for (ff = 0;; ff ^= 1) {
					double dir;			/* Direction to adjust brightness */
					double sv[3], val1;				/* Scaled values */
					double Yxy[3];		/* Yxy of current reading */
					double rgbdir[3];	/* Direction to adjust RGB */
					double rgbxdir[3];	/* Biggest to move */
					double bdir, terr;
					int bx = 0;
					
					if ((rv = dr->read(dr, tcols+1, 1, 0, 0, 1, ' ')) != 0) {
						if (rv == 4)
							break;			/* User is done with this adjustment */
						dr->del(dr);
						error("display read failed with '%s'\n",disprd_err(rv));
					} 
				
					/* Compute 1% direction */
					sv[0] = tcols[1].aXYZ[0] * tcols[2].aXYZ[1]/tcols[2].aXYZ[0];
					sv[1] = tcols[1].aXYZ[1];
					sv[2] = tcols[1].aXYZ[2] * tcols[2].aXYZ[2]/tcols[2].aXYZ[2];
					val1 = sv[1];
					if (sv[0] > val1)
						val1 = sv[0];
					if (sv[2] > val1)
						val1 = sv[2];

					dir = tar1 - val1;
					if (fabs(dir) < 0.01)
						dir = 0.0;

					/* Compute numbers for black point error and direction */
					icmYxy2XYZ(tLab, tYxy);
					tLab[0] /= tLab[1];
					tLab[2] /= tLab[1];
					tLab[1] /= tLab[1];
					icmAry2XYZ(tXYZ, tLab);				/* Lab white reference */
					icmXYZ2Lab(&tXYZ, tLab, tLab);

					icmXYZ2Yxy(Yxy, tcols[1].aXYZ);
					icmAry2Ary(Lab, tcols[1].aXYZ);
					Lab[0] /= Lab[1];
					Lab[2] /= Lab[1];
					Lab[1] /= Lab[1];
					icmXYZ2Lab(&tXYZ, Lab, Lab);
	
					/* Compute dot products */
					bdir = 0.0;
					for (i = 0; i < 3; i++) {
						double rgbLab[3];
						
						icmXYZ2Lab(&tXYZ, rgbLab, rgbXYZ[i]);
						rgbdir[i] = (tLab[1] - Lab[1]) * (rgbLab[1] - Lab[1])
						          + (tLab[2] - Lab[2]) * (rgbLab[2] - Lab[2]);
						rgbxdir[i] = 0.0;
						if (fabs(rgbdir[i]) > fabs(bdir)) {
							bdir = rgbdir[i];
							bx = i;
						}
					}

					/* See how close to the target we are */
					terr = sqrt((tLab[1] - Lab[1]) * (tLab[1] - Lab[1])
					          + (tLab[2] - Lab[2]) * (tLab[2] - Lab[2]));
					if (terr < 0.1)
						rgbdir[0] = rgbdir[1] = rgbdir[2] = 0.0;
					rgbxdir[bx] = rgbdir[bx];

			 		printf("\r%c Current B %.2f, x %.4f, y %.4f  DE %4.1f  R%c%c G%c%c B%c%c ",
					       ff == 0 ? '/' : '\\',
					       val1, Yxy[1], Yxy[2],
					       icmCIE2K(tLab, Lab),
					       rgbdir[0] < 0.0 ? '-' : rgbdir[0] > 0.0 ? '+' : '=',
					       rgbxdir[0] < 0.0 ? '-' : rgbxdir[0] > 0.0 ? '+' : ' ',
					       rgbdir[1] < 0.0 ? '-' : rgbdir[1] > 0.0 ? '+' : '=',
					       rgbxdir[1] < 0.0 ? '-' : rgbxdir[1] > 0.0 ? '+' : ' ',
					       rgbdir[2] < 0.0 ? '-' : rgbdir[2] > 0.0 ? '+' : '=',
					       rgbxdir[2] < 0.0 ? '-' : rgbxdir[2] > 0.0 ? '+' : ' ');
					fflush(stdout);
					c = poll_con_char();
					if (c == ' ')
						break;
				}
				printf("\n");

			/* Report on how well we current meet the targets */
			} else if (c == '5') {
				int nat = 0;		/* NZ if using native white as target */
				col tcols[4] = {	/* Set of test colors */
					{ 0.0, 0.0, 0.0 },
					{ 0.5, 0.5, 0.5 },
					{ 1.0, 1.0, 1.0 },
					{ 0.0, 0.0, 0.0 }	/* 1% test value */
				};
				double tYxy[3];				/* Target white chromaticities */
				double sv[3], val1;			/* Scaled values */
				double mgamma, tarw, tar1, dev1, tarh;
				icmXYZNumber tXYZ;
				double tLab[3], wYxy[3], wLab[3], bYxy[3], bLab[3];
				double ct, ct_de;			/* Color temperature & delta E to white locus */
		
				printf("Doing check measurements\n");

				/* Do an initial set of readings to set 1% output mark */
				if ((rv = dr->read(dr, tcols, 3, 0, 0, 1, 0)) != 0) {
					dr->del(dr);
					error("display read failed with '%s'\n",disprd_err(rv));
				}
				
				/* Approximate Gamma - use the gross curve shape for robustness */
				mgamma = pop_gamma(tcols[0].aXYZ[1], tcols[1].aXYZ[1], tcols[2].aXYZ[1]);

				dev1 = pow(0.01, 1.0/mgamma);
				tcols[3].r = tcols[3].g = tcols[3].b = dev1;
				tar1 = 0.01 * tcols[2].aXYZ[1];

				/* Read the 1% value */
				if ((rv = dr->read(dr, tcols+3, 1, 0, 0, 1, 0)) != 0) {
					dr->del(dr);
					error("display read failed with '%s'\n",disprd_err(rv));
				}
				sv[0] = tcols[3].aXYZ[0] * tcols[2].aXYZ[1]/tcols[2].aXYZ[0];
				sv[1] = tcols[3].aXYZ[1];
				sv[2] = tcols[3].aXYZ[2] * tcols[2].aXYZ[2]/tcols[2].aXYZ[2];
				val1 = sv[1];
				if (sv[0] > val1)
					val1 = sv[0];
				if (sv[2] > val1)
					val1 = sv[2];

				/* Figure out the target white brightness */
				/* Note we're not taking the device gamut into account here */
				if (tbright > 0.0)			/* Given brightness */
					tarw = tbright;
				else						/* Native/maximum brightness */
					tarw = tcols[2].aXYZ[1];

				/* Figure out the target white chromaticity */
				if (wpx > 0.0 || wpy > 0.0) {	/* xy coordinates */
					tYxy[0] = 1.0;
					tYxy[1] = wpx;
					tYxy[2] = wpy;
		
				} else if (temp > 0.0) {		/* Daylight color temperature */
					double XYZ[3];
					if (planckian)
						rv = icx_ill_sp2XYZ(XYZ, icxOT_default, NULL, icxIT_Ptemp, temp, NULL);
					else
						rv = icx_ill_sp2XYZ(XYZ, icxOT_default, NULL, icxIT_Dtemp, temp, NULL);
					if (rv != 0)
						error("Failed to compute XYZ of target color temperature %f\n",temp);
					icmXYZ2Yxy(tYxy, XYZ);
				} else {						/* Native white */
					icmXYZ2Yxy(tYxy, tcols[2].aXYZ);
					nat = 1;
				}

				/* Figure out the target 50% device output value */
				x.tooff = tcols[0].aXYZ[1]/tcols[2].aXYZ[1];	/* Aprox. black output offset */
				if (x.gammat == 0 || x.gammat == 2) {
					if (x.gammat == 0) {
						tarh = pow(0.5, gamma);
					} else {
						if (0.5 <= 0.03928)
							tarh = 0.5/12.92;
						else
							tarh = pow((0.055 + 0.5)/1.055, 2.4);
						/* Allow for zero output offset */
						tarh = x.tooff + (1.0 - x.tooff) * tarh;
					}

				} else {	/* L* curve */
					tarh = 50.0;
	
					/* Convert L* to Y */
					tarh = (tarh + 16.0)/116.0;
					if (tarh > 24.0/116.0)
						tarh = pow(tarh,3.0);
					else
						tarh = (tarh - 16.0/116.0)/7.787036979;
					/* Allow for zero output offset */
					tarh = x.tooff + (1.0 - x.tooff) * tarh;
				}
				/* Convert from Y fraction to absolute Y */
				tarh = tarh * tcols[2].aXYZ[1];

				/* Compute various white point values */
				icmYxy2XYZ(tLab, tYxy);
				tLab[0] /= tLab[1];
				tLab[2] /= tLab[1];
				tLab[1] /= tLab[1];
				icmAry2XYZ(tXYZ, tLab);
				icmXYZ2Lab(&tXYZ, tLab, tLab);

				icmXYZ2Yxy(wYxy, tcols[2].aXYZ);
				icmAry2Ary(wLab, tcols[2].aXYZ);
				wLab[0] /= wLab[1];
				wLab[2] /= wLab[1];
				wLab[1] /= wLab[1];
				icmXYZ2Lab(&tXYZ, wLab, wLab);

				icmXYZ2Yxy(bYxy, tcols[3].aXYZ);
				icmAry2Ary(bLab, tcols[3].aXYZ);
				bLab[0] /= bLab[1];
				bLab[2] /= bLab[1];
				bLab[1] /= bLab[1];
				icmXYZ2Lab(&tXYZ, bLab, bLab);

				/* And color temperature */
				ct = comp_ct(&ct_de, NULL, planckian, dovct, tcols[2].aXYZ);

				printf("\n");

				if (tbright > 0.0)			/* Given brightness */
					printf("  Target Brightness = %.2f, Current = %5.2f, error = % .1f%%\n",
				       tarw, tcols[2].aXYZ[1], 
				       100.0 * (tcols[2].aXYZ[1] - tarw)/tarw);
				else
					printf("  Current Brightness = %.2f\n", tcols[2].aXYZ[1]);
				
				printf("  Target 50%% Level  = %.2f, Current = %5.2f, error = % .1f%%\n",
				       tarh, tcols[1].aXYZ[1], 
				       100.0 * (tcols[1].aXYZ[1] - tarh)/tarw);
				
				printf("  Target Near Black = %5.2f, Current = %5.2f, error = % .1f%%\n",
				       tar1, val1, 
				       100.0 * (val1 - tar1)/tarw);

				if (!nat)
					printf("  Target white = x %.4f, y %.4f, Current = x %.4f, y %.4f, error = %5.2f DE\n",
				       tYxy[1], tYxy[2], wYxy[1], wYxy[2], icmCIE2K(tLab, wLab));
				else
					printf("  Current white = x %.4f, y %.4f, %c%cT %4.0fK DE %4.1f\n",
				       wYxy[1], wYxy[2], dovct ? 'V' : 'C', planckian ? 'C' : 'D', ct,ct_de);

				printf("  Target black = x %.4f, y %.4f, Current = x %.4f, y %.4f, error = %5.2f DE\n",
				       tYxy[1], tYxy[2], bYxy[1], bYxy[2], icmCIE2K(tLab, bLab));

			} else if (c == '6') {
				if (!verb) {		/* Tell user command has been accepted */
					if (verify == 2)
						printf("Commencing device verification\n");
					else
						printf("Commencing device calibration\n");
				}
				break;
			} else if (c == '7' || c == 0x03 || c == 0x1b) {
				printf("Exiting\n");
				dr->del(dr);
				exit(0);
			}
		}
	}
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

	/* Take a small number of readings, and compute basic */
	/* informations such as black & white, white target, */
	/* aproximate matrix based display forward and reverse model. */
#ifdef RES_TEST
	printf("Running 8 bit res test.\n");
	/* See if we can detect whether the LUT's have better than 8 bit precision */
	{
		col ttt[3] = {	/* Base set of test colors */
			{ 220.0/255.0, 220.0/255.0, 220.0/255.0 },
			{ 220.5/255.0, 220.5/255.0, 220.5/255.0 },
			{ 221.0/255.0, 221.0/255.0, 221.0/255.0 }
		};
		double L1 = 0.0, L_5 = 0.0, L_1 = 0.0;
		int n;
		for (n = 1; n < 30; n++) {
			if ((rv = dr->read(dr, ttt, 3, 0, 1, 0)) != 0) {
				dr->del(dr);
				error("display read failed with '%s'\n",disprd_err(rv));
			} 
			L_1 += ttt[0].aXYZ[1];
			L_5 += ttt[1].aXYZ[1];
			L1  += ttt[2].aXYZ[1];
			printf("Measured %f %f %f, average over %d = %f, %f, %f\n",
			        ttt[0].aXYZ[1], ttt[1].aXYZ[1], ttt[2].aXYZ[1], n, L_1/n, L_5/n, L1/n);
		}
	}
#endif

	/* Read the base test set */
	{
		icmXYZNumber mrd;		/* Number for matrix */
		icmXYZNumber mgn;
		icmXYZNumber mbl;
		icmXYZNumber mwh;

		col base[6] = {	/* Base set of test colors */
			{ 0.0, 0.0, 0.0 },		/* 0 - Black */
			{ 1.0, 0.0, 0.0 },		/* 1 - Red */
			{ 0.0, 1.0, 0.0 },		/* 2 - Green */
			{ 0.0, 0.0, 1.0 },		/* 3 - Blue */
			{ 1.0, 1.0, 1.0 },		/* 4 - White */
			{ 0.0, 0.0, 0.0 }		/* 5 - Black */
		};

		if (verb) {
			if (verify == 2)
				printf("Commencing device verification\n");
			else
				printf("Commencing device calibration\n");
		}
		if ((rv = dr->read(dr, base, 6, 1, 6, 1, 0)) != 0) {
			dr->del(dr);
			error("display read failed with '%s'\n",disprd_err(rv));
		} 

		if (base[0].aXYZ_v == 0) {
			dr->del(dr);
			error("Failed to get an absolute XYZ value from the instrument!\n");
		}

		/* Average black relative from 2 readings */
		x.bk[0] = 0.5 * (base[0].aXYZ[0] + base[5].aXYZ[0]);
		x.bk[1] = 0.5 * (base[0].aXYZ[1] + base[5].aXYZ[1]);
		x.bk[2] = 0.5 * (base[0].aXYZ[2] + base[5].aXYZ[2]);

		if (verb) {
			printf("Black = XYZ %6.2f %6.2f %6.2f\n",x.bk[0],x.bk[1],x.bk[2]);
			printf("Red   = XYZ %6.2f %6.2f %6.2f\n",base[1].aXYZ[0], base[1].aXYZ[1], base[1].aXYZ[2]);
			printf("Green = XYZ %6.2f %6.2f %6.2f\n",base[2].aXYZ[0], base[2].aXYZ[1], base[2].aXYZ[2]);
			printf("Blue  = XYZ %6.2f %6.2f %6.2f\n",base[3].aXYZ[0], base[3].aXYZ[1], base[3].aXYZ[2]);
			printf("White = XYZ %6.2f %6.2f %6.2f\n",base[4].aXYZ[0], base[4].aXYZ[1], base[4].aXYZ[2]);
		}

		/* Copy other readings into place */
		icmAry2Ary(x.bk, base[0].aXYZ);
		icmAry2Ary(x.wh, base[4].aXYZ);
		icmAry2XYZ(x.twN, x.wh);	/* Use this as Lab reference white until we establish target */

		icmAry2XYZ(mrd, base[1].aXYZ);
		icmAry2XYZ(mgn, base[2].aXYZ);
		icmAry2XYZ(mbl, base[3].aXYZ);
		icmAry2XYZ(mwh, base[4].aXYZ);

		/* Setup forward matrix */
		if (icmRGBprim2matrix(mwh, mrd, mgn, mbl, x.fm)) {
			dr->del(dr);
			error("Aprox. fwd matrix unexpectedly singular\n");
		}

#ifdef DEBUG
		if (verb) {
			printf("Forward matrix is:\n");
			printf("%f %f %f\n", x.fm[0][0], x.fm[0][1], x.fm[0][2]);
			printf("%f %f %f\n", x.fm[1][0], x.fm[1][1], x.fm[1][2]);
			printf("%f %f %f\n", x.fm[2][0], x.fm[2][1], x.fm[2][2]);
		}
#endif

		/* Compute bwd matrix */
		if (icmInverse3x3(x.bm, x.fm)) {
			dr->del(dr);
			error("Inverting aprox. fwd matrix failed");
		}
	}

	/* Now do some more readings, to compute the basic per channel */
	/* transfer characteristics, and then a device model. */
	if (verify != 2 && !doupdate) {
		col set[4];				/* Variable to read up to 4 values from the display */
		sxyz *asrgb[4];			/* samples for r, g, b & w */

		for (j = 0; j < 4; j++) {
			if ((asrgb[j] = (sxyz *)malloc(isteps * sizeof(sxyz))) == NULL) {
				dr->del(dr);
				error ("Malloc of array of readings failed");
			}
		}

		for (i = 0; i < isteps; i++) {
			double vv, v[4];

#ifdef __APPLE__
			gcc_bug_fix(i);
#endif
			vv = i/(isteps - 1.0);
			vv = pow(vv, MOD_DIST_POW);
			for (j = 0; j < 4; j++) {
				v[j] = vv;
				set[j].r = set[j].g = set[j].b = 0.0;
				if (j == 0)
					set[j].r = v[j];
				else if (j == 1)
					set[j].g = v[j];
				else if (j == 2)
					set[j].b = v[j];
				else
					set[j].r = set[j].g = set[j].b = v[j];
			}

			if ((rv = dr->read(dr, set, 4, i * 4+1, isteps * 4, 1, 0)) != 0) {
				dr->del(dr);
				error("display read failed with '%s'\n",disprd_err(rv));
			} 
			for (j = 0; j < 4; j++) {
//printf("~1 R = %f, G = %f, B = %f, XYZ = %f %f %f\n",
//set[j].r, set[j].g, set[j].b, set[j].aXYZ[0], set[j].aXYZ[1], set[j].aXYZ[2]);
				asrgb[j][i].v = v[j];
				asrgb[j][i].xyz[0] = set[j].aXYZ[0];
				asrgb[j][i].xyz[1] = set[j].aXYZ[1];
				asrgb[j][i].xyz[2] = set[j].aXYZ[2];
			}
		}

		/* Convert RGB channel samples to curves */
		{
			mcvco *sdv;		/* Points used to create cvs[], RGB */
			double blrgb[3];
			double *op;		/* Parameters to optimise */
			double *sa;		/* Search area */
			double re;		/* Residual error */

			/* Transform measured  black back to linearised RGB values */
			icmMulBy3x3(blrgb, x.bm, x.bk);
//printf("~1 model black should be %f %f %f\n", x.bk[0], x.bk[1], x.bk[2]);
//printf("~1 linearised RGB should be %f %f %f\n", blrgb[0], blrgb[1], blrgb[2]);

			if ((sdv = malloc(sizeof(mcvco) * isteps)) == NULL) {
				dr->del(dr);
				error ("Malloc of scattered data points failed");
			}
			for (k = 0; k < 3; k++) {	/* Create the model curves */
				for (i = 0; i < isteps; i++) {
					sdv[i].p = asrgb[k][i].v;
					sdv[i].v = VLENGTH(asrgb[k][i].xyz);
					sdv[i].w = 1.0;
//printf("~1 chan %d, entry %d, p = %f, v = %f from XYZ %f %f %f\n",
//k,i,x.sdv[k][i].p,x.sdv[k][i].v, asrgb[k][i].xyz[0], asrgb[k][i].xyz[1], asrgb[k][i].xyz[2]);
				}
				if ((x.dcvs[k] = new_mcv()) == NULL) {
					dr->del(dr);
					error("new_mcv x.dcvs[%d] failed",k);
				}
				x.dcvs[k]->fit(x.dcvs[k], 0, fitord, sdv, isteps, 5.0);

				/* Scale the whole curve so the output is scaled to 1.0 */
				x.dcvs[k]->force_scale(x.dcvs[k], 1.0);

				/* Force curves to produce this lrgb for 0.0 */
				x.dcvs[k]->force_0(x.dcvs[k], blrgb[k]);
			}
			free(sdv);

			/* Setup list of reference points ready for optimisation */
			x.nrp = 4 * isteps;
			if ((x.rp = (optref *)malloc(sizeof(optref) * x.nrp)) == NULL) {
				dr->del(dr);
				error ("Malloc of measurement reference points failed");
			}
			for (k = 0; k < 4; k++) {
				for (i = 0; i < isteps; i++) {
					int ii = k * isteps + i;
					double v[3];
				
					v[0] = v[1] = v[2] = 0.0;
					if (k == 0)
						v[k] = asrgb[k][i].v;
					else if (k == 1)
						v[k] = asrgb[k][i].v;
					else if (k == 2)
						v[k] = asrgb[k][i].v;
					else
						v[0] = v[1] = v[2] = asrgb[k][i].v;
					icmAry2Ary(x.rp[ii].dev, v);
					icmXYZ2Lab(&x.twN, x.rp[ii].lab, asrgb[k][i].xyz);
					if (k == 3)		/* White */
						x.rp[ii].w = 0.5;
					else
						x.rp[ii].w = 0.16667;
				}
			}

			/* Get parameters and setup for optimisation */
			op = dev_get_params(&x);
			if ((sa = malloc(x.np * sizeof(double))) == NULL)
				error ("Malloc of scattered data points failed");
			
			for (i = 0; i < x.np; i++)
				sa[i] = 0.1;

			/* Do optimisation */
#ifdef NEVER
			if (powell(&re, x.np, op, sa, 1e-5, 2000, dev_opt_func, (void *)&x) != 0)
				error ("Model powell failed, re = %f",re);
#else
			if (conjgrad(&re, x.np, op, sa, 1e-5, 2000, dev_opt_func, dev_dopt_func, (void *)&x) != 0) {
				if (re > 1e-2)
					error("Model conjgrad failed, residual error = %f",re);
				else
					warning("Model conjgrad failed, residual error = %f",re);
			}
#endif

			/* Put optimised parameters in place */
			dev_put_params(&x, op);

			free(x.rp);
			free(x.dtin_iv);		/* Free temporary arrays */
			x.rp = NULL;
			x.nrp = 0;
			free(sa);
			free(op);
		}

#ifdef DEBUG_PLOT
		/* Plot the current calc curves */
		{
			#define	XRES 256
			double xx[XRES];
			double yy[3][XRES];
			double xyz[3];
			for (i = 0; i < XRES; i++) {
				xx[i] = i/(XRES-1.0);
				for (j = 0; j < 3; j++)
					yy[j][i] = x.dcvs[j]->interp(x.dcvs[j], xx[i]);
			}
			printf("Channel curves\n");
			do_plot(xx,yy[0],yy[1],yy[2],XRES);
			#undef XRES
		}
#endif

		/* We're done with asrgb[] */
		for (j = 0; j < 4; j++)
			free(asrgb[j]);
	}

#ifdef CHECK_MODEL
	/* Check how well our fwd model agrees with the device */
	if (verify != 2) {
		col set[3];					/* Variable to read up to 3 values from the display */
		int nn = 27;
		double alab[3], mxyz[3], mlab[3];	/* Actual and model Lab */
		double mnerr;		/* Maximum neutral error */
		double mnv;			/* Value where maximum error is */
		double anerr;		/* Average neutral error */

		mnerr = anerr = 0.0;
		for (i = 0; i < nn; i++) {
			double vv, v[3];
			double de;

			vv = i/(nn - 1.0);
			vv = pow(vv, CHECK_DIST_POW);
			v[0] = v[1] = v[2] = vv;
			set[0].r = set[0].g = set[0].b = vv;

			if ((rv = dr->read(dr, set, 1, i+1, nn, 1, 0)) != 0) {
				dr->del(dr);
				error("display read failed with '%s'\n",disprd_err(rv));
			} 
			icmXYZ2Lab(&x.twN, alab, set[0].aXYZ);
		
			fwddev(&x, mxyz, v);
			icmXYZ2Lab(&x.twN, mlab, mxyz);

			de = icmCIE2K(mlab, alab);
			if (de > mnerr) {
				mnerr = de;
				mnv = vv;
			}
			anerr += de;
			
			printf("RGB %.3f -> Lab %.2f %.2f %.2f, model %.2f %.2f %.2f, DE %f\n",
			vv, alab[0], alab[1], alab[2], mlab[0], mlab[1], mlab[2],de);
		}
		anerr /= (double)nn;
		printf("Maximum error (@ %f) = %f deltaE\n",mnv, mnerr);
		printf("Average error = %f deltaE\n",anerr);
	}
#endif /* CHECK_MODEL */

	/* Figure out our calibration curve parameter targets */
	if (!doupdate) {
		double _mat[3][3];
		double *mat[3];		/* Forward matrix in numlib format */
		double xx[3];		/* RHS and solution */
		double kk;			/* Smallest k */
		mat[0] = _mat[0], mat[1] = _mat[1], mat[2] = _mat[2];

		/* Figure out the target white point */
		if (wpx > 0.0 || wpy > 0.0) {	/* xy coordinates */
			double Yxy[3];
			Yxy[0] = 1.0;
			Yxy[1] = wpx;
			Yxy[2] = wpy;
			icmYxy2XYZ(x.twh, Yxy);

		} else if (temp > 0.0) {		/* Daylight color temperature */
			if (planckian)
				rv = icx_ill_sp2XYZ(x.twh, icxOT_default, NULL, icxIT_Ptemp, temp, NULL);
			else
				rv = icx_ill_sp2XYZ(x.twh, icxOT_default, NULL, icxIT_Dtemp, temp, NULL);
			if (rv != 0)
				error("Failed to compute XYZ of target color temperature %f\n",temp);
//printf("~1 Raw target from temp %f XYZ = %f %f %f\n",temp,x.twh[0],x.twh[1],x.twh[2]);
		} else {						/* Native white */
			x.twh[0] = x.wh[0]/x.wh[1];
			x.twh[1] = x.wh[1]/x.wh[1];
			x.twh[2] = x.wh[2]/x.wh[1];
		}

		/* Convert it to absolute white target */
		if (tbright > 0.0) {			/* Given brightness */
			x.twh[0] *= tbright;
			x.twh[1] *= tbright;
			x.twh[2] *= tbright;
		} else {						/* Native/maximum brightness */
			x.twh[0] *= x.wh[1];
			x.twh[1] *= x.wh[1];
			x.twh[2] *= x.wh[1];
			if (verb)
				printf("Initial native brightness target = %f cd/m^2\n", x.twh[1]);
		}

		/* Now see if target white will fit in gamut. */

		/* Setup for intersection solution */
		kk = 1e6;
		for (i = 0; i < 3; i++) {	/* Intersect with upper 3 faces of cube */

			/* Setup equation to solve */
			for (j = 0; j < 3; j++) {
				for (k = 0; k < 3; k++) {
					if (k == i)
						mat[j][k] = -x.twh[j];		/* Target white in column */
					else
						mat[j][k] = x.fm[j][k];		/* Matrix */
				}
			}
			xx[0] = -x.fm[0][i];	/* RHS is matrix column */
			xx[1] = -x.fm[1][i];
			xx[2] = -x.fm[2][i];

#ifdef DEBUG
			if (verb) {
				printf("A.X = B to solve is:\n");
				printf("%f %f %f . X = %f\n", mat[0][0], mat[0][1], mat[0][2], xx[0]);
				printf("%f %f %f . X = %f\n", mat[1][0], mat[1][1], mat[1][2], xx[1]);
				printf("%f %f %f . X = %f\n", mat[2][0], mat[2][1], mat[2][2], xx[2]);
			}
#endif
			/* Solve for scale factor and other 2 device values */
			if (solve_se(mat, xx, 3)) {
				warning ("White point vector intesection failure");
			} else {
#ifdef DEBUG
				if (verb)
					printf("Got solution for index %d of %f %f %f\n",i,xx[0],xx[1],xx[2]);
#endif
				if (xx[0] > 0.0 && xx[1] > 0.0 && xx[2] > 0.0 && xx[i] < kk) {
					kk = xx[i];
				}
			}

#ifdef DEBUG
			{	/* Verify the solution */
				double rgb[3], xyz[3];
				
				rgb[0] = xx[0];
				rgb[1] = xx[1];
				rgb[2] = xx[2];
				rgb[i] = 1.0;
				
				xyz[0] = x.fm[0][0] * rgb[0] + x.fm[0][1] * rgb[1] + x.fm[0][2] * rgb[2];
				xyz[1] = x.fm[1][0] * rgb[0] + x.fm[1][1] * rgb[1] + x.fm[1][2] * rgb[2];
				xyz[2] = x.fm[2][0] * rgb[0] + x.fm[2][1] * rgb[1] + x.fm[2][2] * rgb[2];

				xyz[0] /= xx[i];
				xyz[1] /= xx[i];
				xyz[2] /= xx[i];

				printf("Verify xyz is %f %f %f, should be %f %f %f\n",
				xyz[0], xyz[1], xyz[2], x.twh[0], x.twh[1], x.twh[2]);
			}
#endif
		}

		if (kk < 0.9999) {
			x.twh[0] *= kk;		/* Scale brightness to fit */
			x.twh[1] *= kk;
			x.twh[2] *= kk;
			if (verb)
				printf("Had to scale brightness from %f to %f to fit within gamut\n",x.twh[1]/kk, x.twh[1]);
//printf("~1 Scaled target from temp %f XYZ = %f %f %f\n",temp,x.twh[0],x.twh[1],x.twh[2]);
		}

		if (verb)
			printf("Target white value is XYZ %f %f %f\n",x.twh[0],x.twh[1],x.twh[2]);

	}
	icmAry2XYZ(x.twN, x.twh);		/* Need this for Lab conversions */

	{
		double tbkLab[3]; 

		icmXYZ2Lab(&x.twN, tbkLab, x.bk);	/* Convert measured black to Lab */
	
		/* Now blend the a* b* with that of the target white point */
		x.tbL[0] = tbkLab[0];
		x.tbL[1] = bkcorrect * 0.0   + (1.0 - bkcorrect) * tbkLab[1];
		x.tbL[2] = bkcorrect * 0.0   + (1.0 - bkcorrect) * tbkLab[2];
	
		/* And make this the black hue to aim for */
		icmLab2XYZ(&x.twN, x.tbk, x.tbL);
		icmAry2XYZ(x.tbN, x.tbk);
		if (verb)
			printf("Adjusted target black XYZ %.2f %.2f %.2f, Lab %.2f %.2f %.2f\n",
	        x.tbk[0], x.tbk[1], x.tbk[2], x.tbL[0], x.tbL[1], x.tbL[2]);
	}

	/* Figure out the gamma curve black offset value */
	/* that will give us the black level we actually have. */
	{
		double yy, tby;			/* Target black y */
		double thyr;			/* Target 50 % y ratio */
		
		/* Make target black Y as high as necessary */
		/* to get the black point hue */
		// ~~~~9999 should do this by increasing L* until XYZ > x.bk ??? */
		tby = x.bk[1];
//printf("Target Y from Y = %f\n",tby);
		yy = x.bk[0] * x.tbk[1]/x.tbk[0];
//printf("Target Y from X = %f\n",yy);
		if (yy > tby)
			tby = yy;
		yy = x.bk[2] * x.tbk[1]/x.tbk[2];
//printf("Target Y from Z = %f\n",yy);
		if (yy > tby)
			tby = yy;

		if (verb) {
			double tbp[3], tbplab[3];
	        tbp[0] = x.tbk[0] * tby/x.tbk[1];
			tbp[1] = tby;
			tbp[2] = x.tbk[2] * tby/x.tbk[1];
			icmXYZ2Lab(&x.twN, tbplab, tbp);
			printf("Target black after min adjust: XYZ %.2f %.2f %.2f, Lab %.2f %.2f %.2f\n",
			        tbp[0], tbp[1], tbp[2], tbplab[0], tbplab[1], tbplab[2]);
		}

		/* Figure out the x.gioff and tgamma needed to get this x.tooff and gamma */
		/* (x.gioff and tgamma are only used for x.gammat == 0, ie. a power) */ 
		x.tooff = tby / x.twh[1];			/* Convert to relative */
		thyr = pow(0.5, gamma);
		x.tgamma = gamma;

//printf("~1 x.tooff = %f, thyr = %f\n",x.tooff,thyr);
		/* Iterative solution */
		for (i = 0; i < 6; i++) {
			double gvv;
			x.gioff = pow(x.tooff, 1.0/x.tgamma);
			gvv = 0.5 + (1.0 - 0.5) * x.gioff;
			x.tgamma = log(thyr) / log(gvv);
		}
//printf("~1 x.gioff = %f, x.tgamma = %f\n",x.gioff,x.tgamma);

	}

	if (verb)
		printf("Current gamma curve offset = %f, Gamma curve power = %f\n",x.gioff,x.tgamma);

	/* - - - - - - - - - - - - - - - - - - - - - */
	/* Start with a scaled down number of test points and refine threshold, */
	/* and double/halve these on each iteration. */
	if (verb)
		printf("Total Iteration %d, Final Samples = %d Final Repeat threshold = %f\n",
		        mxits, rsteps, errthr);
	if (verify == 2) {
		rsteps = VER_RES;
		errthr = 0.0;
	} else {
		rsteps /= (1 << (mxits-1));
		errthr *= (double)(1 << (mxits-1));
	}

	/* Setup the initial calibration test point values */
	init_csamp(&asgrey, &x, doupdate, verify, verify == 2 ? 1 : 0, rsteps);
	
	/* Calculate the initial calibration curve values */
	if (verify != 2 && !doupdate) {
		int nsamp = 128;
		mcvco *sdv[3];				/* Scattered data for creating mcv */

		for (j = 0; j < 3; j++) {
			if ((x.rdac[j] = new_mcv()) == NULL) {
				dr->del(dr);
				error("new_mcv x.rdac[%d] failed",j);
			}
		}

		for (j = 0; j < 3; j++) {
			if ((sdv[j] = malloc(sizeof(mcvco) * rsteps)) == NULL) {
				dr->del(dr);
				error ("Malloc of scattered data points failed");
			}
		}

		if (verb)
			printf("Creating initial calibration curves...\n");

		/* Copy the sample points */
		for (i = 0; i < rsteps; i++) {
			for (j = 0; j < 3; j++) {
				sdv[j][i].p = asgrey.s[i].v;
				sdv[j][i].v = asgrey.s[i].rgb[j];
				sdv[j][i].w = 1.0;
			}
		}
		if (x.nat)		/* Make curve go thought white if possible */
			sdv[0][rsteps-1].w = sdv[1][rsteps-1].w = sdv[2][rsteps-1].w = 10.0;

		/* Create an initial set of RAMDAC curves */
		for (j = 0; j < 3; j++)
			x.rdac[j]->fit(x.rdac[j], 0, fitord, sdv[j], rsteps, 1.0);

		/* Make sure that if we are using native brightness and white point, */
		/* that the curves go to a perfect 1.0 ... */
		if (x.nat) {
			for (j = 0; j < 3; j++)
				x.rdac[j]->force_1(x.rdac[j], 1.0);
		}

		for (j = 0; j < 3; j++)
			free (sdv[j]);
	}

#ifdef DEBUG_PLOT
	/* Plot the current curves */
	if (verify != 2) {
		#define	XRES 255
		double xx[XRES];
		double y1[XRES];
		double y2[XRES];
		double y3[XRES];
		double rgb[3];
		for (i = 0; i < XRES; i++) {
			double drgb[3], rgb[3];
			xx[i] = i/(XRES-1.0);
			rgb[0] = rgb[1] = rgb[2] = xx[i];
			for (j = 0; j < 3; j++)
				drgb[j] = x.rdac[j]->interp(x.rdac[j], rgb[j]);
			y1[i] = drgb[0];
			y2[i] = drgb[1];
			y3[i] = drgb[2];
		}
		printf("Current ramdac curves\n");
		do_plot(xx,y1,y2,y3,XRES);
		#undef XRES
	}
#endif

	/* Now we go into the main verify & refine loop */
	for (it = verify == 2 ? mxits : 0; it < mxits || verify != 0; rsteps *= 2, errthr /= 2.0, it++)  {

		col set[2];				/* Variable to read one or two value from the display */

		/* Verify pass */
		if (it >= mxits)
			rsteps = VER_RES;		/* Fixed verification resolution */

		/* re-init asgrey if the number of test points has changed */
		reinit_csamp(&asgrey, &x, verify, (verify == 2 || it >= mxits) ? 1 : 0, rsteps);

		if (verb) {
			if (it >= mxits)
				printf("Doing verify pass with %d sample points\n",rsteps);
			else
				printf("Doing iteration %d with %d sample points and repeat threshold of %f DE\n",
				                                                             it+1,rsteps, errthr);
		}
		/* Read and adjust each step */
		/* Do this white to black to track drift in native white point */
		for (i = rsteps-1; i >= 0; i--) {
			double rpt;
			double bestrgb[3];		/* In case we fail */
			double bestxyz[3];
			double bestde = 1e7;
			double bestdc = 1e7;
			double rgain = REFINE_GAIN;		/* Scale down if lots of repeats */

			/* Until we meet the necessary accuracy or give up */
			for (rpt = 0;rpt < MAX_RPTS; rpt++) {
				int gworse = 0;		/* information flag */
				double wde;			/* informational */

				set[0].r = asgrey.s[i].rgb[0];
				set[0].g = asgrey.s[i].rgb[1];
				set[0].b = asgrey.s[i].rgb[2];
				set[0].id = NULL;

				/* Read patches (no auto cr in case we repeat last patch) */
				if ((rv = dr->read(dr, set, 1, rsteps-i, rsteps, 0, 0)) != 0) {
					dr->del(dr);
					error("display read failed with '%s'\n",disprd_err(rv));
				} 
	
				icmAry2Ary(asgrey.s[i].pXYZ, asgrey.s[i].XYZ);	/* Remember previous XYZ */
				icmAry2Ary(asgrey.s[i].XYZ, set[0].aXYZ);		/* Transfer current reading */

				/* If native white and we've just measured it, */
				/* and we're not doing a verification, */
				/* adjust all the other point targets txyz to track the white. */
				if (x.nat && i == (rsteps-1) && it < mxits && asgrey.s[i].v == 1.0) {

					icmAry2Ary(x.twh, asgrey.s[i].XYZ);	/* Set current white */
					icmAry2XYZ(x.twN, x.twh);			/* Need this for Lab conversions */
														/* Recompute txyz */
					init_csamp_txyz(&asgrey, &x);
//printf("~1 Just reset target white to native white\n");
				}

				/* Compute the current change wanted to hit target */
				icmSub3(asgrey.s[i].deXYZ, asgrey.s[i].tXYZ, asgrey.s[i].XYZ);
				asgrey.s[i].de = icmXYZLabDE(&x.twN, asgrey.s[i].tXYZ, asgrey.s[i].XYZ);
				/* Eudclidian difference of XYZ, because this doesn't always track Lab */
				asgrey.s[i].dc = icmLabDE(asgrey.s[i].tXYZ, asgrey.s[i].XYZ);

				/* Compute change from last XYZ */
				icmSub3(asgrey.s[i].dXYZ, asgrey.s[i].XYZ, asgrey.s[i].pXYZ);


#ifdef DEBUG
				printf("\n\nTest point %d, v = %f\n",rsteps - i,asgrey.s[i].v);
				printf("Current rgb %f %f %f -> XYZ %f %f %f, de %f, dc %f\n", 
				asgrey.s[i].rgb[0], asgrey.s[i].rgb[1], asgrey.s[i].rgb[2],
				asgrey.s[i].XYZ[0], asgrey.s[i].XYZ[1], asgrey.s[i].XYZ[2],
				asgrey.s[i].de, asgrey.s[i].dc);
				printf("Target XYZ %f %f %f, delta needed %f %f %f\n", 
				asgrey.s[i].tXYZ[0], asgrey.s[i].tXYZ[1], asgrey.s[i].tXYZ[2],
				asgrey.s[i].deXYZ[0], asgrey.s[i].deXYZ[1], asgrey.s[i].deXYZ[2]);
				if (rpt > 0) {
					printf("Intended XYZ change %f %f %f, actual change %f %f %f\n", 
					asgrey.s[i].pdXYZ[0], asgrey.s[i].pdXYZ[1], asgrey.s[i].pdXYZ[2],
					asgrey.s[i].dXYZ[0], asgrey.s[i].dXYZ[1], asgrey.s[i].dXYZ[2]);
				}
#endif

				if (it < mxits) {		/* Not verify, apply correction */
					int impj = 0;		/* We improved the Jacobian */
					int dclip = 0;		/* We clipped the new RGB */

#ifdef ADJ_JACOBIAN
					int isclipped = 0;

#ifndef CLIP		/* Check for cliping */
					/* Don't try and update the Jacobian if the */	
					/* device values are going out of gamut, */
					/* and being clipped without Jac correction being aware. */
					for (j = 0; j < 3; j++) {
						if (asgrey.s[i].rgb[j] <= 0.0 || asgrey.s[i].rgb[j] >= 1.0) {
							isclipped = 1;
							break;
						}
					}
#endif /* !CLIP */

					/* Compute a correction to the Jacobian if we can. */
					/* (Don't do this unless we have a solid previous */
					/* reading for this patch) */ 
					if (rpt > 0 && isclipped == 0) {
						double nsdrgb;			/* Norm squared of pdrgb */
						double spdrgb[3];		/* Scaled previous delta rgb */
						double dXYZerr[3];		/* Error in previous prediction */
						double jadj[3][3];		/* Adjustment to Jacobian */
						double tj[3][3];		/* Temp Jacobian */
						double itj[3][3];		/* Temp inverse Jacobian */

						/* Use Broyden's Formula */
						icmSub3(dXYZerr, asgrey.s[i].dXYZ, asgrey.s[i].pdXYZ);
//printf("~1 Jacobian error = %f %f %f\n", dXYZerr[0], dXYZerr[1], dXYZerr[2]);
						nsdrgb = icmNorm3sq(asgrey.s[i].pdrgb);
						if (nsdrgb > 1e-4) {
							icmScale3(spdrgb, asgrey.s[i].pdrgb, 1.0/nsdrgb);
							icmTensMul3(jadj, dXYZerr, spdrgb);

#ifdef DEBUG
							/* Check that new Jacobian predicts previous delta XYZ */
							{
								double eXYZ[3];
	
								/* Make a full adjustment to temporary Jac */
								icmAdd3x3(tj, asgrey.s[i].j, jadj);
								icmMulBy3x3(eXYZ, tj, asgrey.s[i].pdrgb);
								icmSub3(eXYZ, eXYZ, asgrey.s[i].dXYZ);
								printf("Jac check resid %f %f %f\n", eXYZ[0], eXYZ[1], eXYZ[2]);
							}
#endif	/* DEBUG */

							/* Add part of our correction to actual Jacobian */
							icmScale3x3(jadj, jadj, JAC_COR_FACT);
							icmAdd3x3(tj, asgrey.s[i].j, jadj);
							if (icmInverse3x3(itj, tj) == 0) {		/* Invert OK */
								icmCpy3x3(asgrey.s[i].j, tj);		/* Use adjusted */
								icmCpy3x3(asgrey.s[i].ij, itj);
								impj = 1;
							}
//else printf("~1 ij failed\n");
						}
//else printf("~1 nsdrgb was 0\n");
					}
//else if (isclipped) printf("~1 no j update: rgb is clipped\n");
#endif	/* ADJ_JACOBIAN */

					/* Track the best solution we've found */
					if (asgrey.s[i].de <= bestde) {
						bestde = asgrey.s[i].de;
						bestdc = asgrey.s[i].dc;
						bestrgb[0] = asgrey.s[i].rgb[0];
						bestrgb[1] = asgrey.s[i].rgb[1];
						bestrgb[2] = asgrey.s[i].rgb[2];
						bestxyz[0] = asgrey.s[i].XYZ[0];
						bestxyz[1] = asgrey.s[i].XYZ[1];
						bestxyz[2] = asgrey.s[i].XYZ[2];

					} else if (asgrey.s[i].dc > bestdc) {
						/* we got worse in Lab and XYZ ! */

						wde = asgrey.s[i].de;

						/* If we've wandered too far, return to best we found */
						if (asgrey.s[i].de > (3.0 * bestde)) {
							asgrey.s[i].de = bestde;
							asgrey.s[i].dc = bestdc;
							asgrey.s[i].rgb[0] = bestrgb[0];
							asgrey.s[i].rgb[1] = bestrgb[1];
							asgrey.s[i].rgb[2] = bestrgb[2];
							asgrey.s[i].XYZ[0] = bestxyz[0];
							asgrey.s[i].XYZ[1] = bestxyz[1];
							asgrey.s[i].XYZ[2] = bestxyz[2];
							icmSub3(asgrey.s[i].deXYZ, asgrey.s[i].tXYZ, asgrey.s[i].XYZ);
						}

						/* If the Jacobian hasn't changed, moderate the gain */
						if (impj == 0)
							rgain *= 0.8;		/* We might be overshooting */
						gworse = 1;
					}

					/* See if we need to repeat */
					if (asgrey.s[i].de <= errthr) {	
						if (verb > 1)
							printf("Point %d Delta E %f, OK\n",rsteps - i,asgrey.s[i].de);
						break;	/* No more retries */
					}
					if ((rpt+1) >= MAX_RPTS) {
						asgrey.s[i].de = bestde;			/* Restore to best we found */
						asgrey.s[i].dc = bestdc;
						asgrey.s[i].rgb[0] = bestrgb[0];
						asgrey.s[i].rgb[1] = bestrgb[1];
						asgrey.s[i].rgb[2] = bestrgb[2];
						if (verb > 1)
							printf("Point %d Delta E %f, Fail\n",rsteps - i,asgrey.s[i].de);
						break;	/* No more retries */
					}
					if (verb > 1) {
						if (gworse)
							printf("Point %d Delta E %f, Repeat (got worse)\n", rsteps - i, wde);
						else
							printf("Point %d Delta E %f, Repeat\n", rsteps - i,asgrey.s[i].de);
					}
				
					/* Compute refinement of rgb */
					icmMulBy3x3(asgrey.s[i].pdrgb, asgrey.s[i].ij, asgrey.s[i].deXYZ);
//printf("~1 delta needed %f %f %f -> delta RGB %f %f %f\n",
//asgrey.s[i].deXYZ[0], asgrey.s[i].deXYZ[1], asgrey.s[i].deXYZ[2],
//asgrey.s[i].pdrgb[0], asgrey.s[i].pdrgb[1], asgrey.s[i].pdrgb[2]);

					/* Gain scale */
					icmScale3(asgrey.s[i].pdrgb, asgrey.s[i].pdrgb, rgain);
//printf("~1 delta RGB after gain scale %f %f %f\n",
//asgrey.s[i].pdrgb[0], asgrey.s[i].pdrgb[1], asgrey.s[i].pdrgb[2]);

#ifdef CLIP
					/* Component wise clip */
					for (j = 0; j < 3; j++) {		/* Check for clip */
						if ((-asgrey.s[i].pdrgb[j]) > asgrey.s[i].rgb[j]) {
							asgrey.s[i].pdrgb[j] = -asgrey.s[i].rgb[j];
							dclip = 1;
						}
						if (asgrey.s[i].pdrgb[j] > (1.0 - asgrey.s[i].rgb[j])) {
							asgrey.s[i].pdrgb[j] = (1.0 - asgrey.s[i].rgb[j]);
							dclip = 1;
						}
					}
#ifdef DEBUG
					if (dclip) printf("delta RGB after clip %f %f %f\n",
					       asgrey.s[i].pdrgb[0], asgrey.s[i].pdrgb[1], asgrey.s[i].pdrgb[2]);
#endif	/* DEBUG */
#endif	/* CLIP */
					/* Compute next on the basis of this one RGB */
					icmAdd3(asgrey.s[i].rgb, asgrey.s[i].rgb, asgrey.s[i].pdrgb);

					/* Save expected change in XYZ */
					icmMulBy3x3(asgrey.s[i].pdXYZ, asgrey.s[i].j, asgrey.s[i].pdrgb);

#ifdef DEBUG
					printf("New rgb %f %f %f from expected del XYZ %f %f %f\n",
					       asgrey.s[i].rgb[0], asgrey.s[i].rgb[1], asgrey.s[i].rgb[2],
					       asgrey.s[i].pdXYZ[0], asgrey.s[i].pdXYZ[1], asgrey.s[i].pdXYZ[2]);
#endif
				} else {	/* Verification, so no repeat */
					break;
				}

			}	/* Next repeat */
		}		/* Next resolution step */
		if (verb)
			printf("\n");			/* Final return for patch count */

#ifdef DEBUG_PLOT
		/* Plot the measured response XYZ */
		{
			#define	XRES 256
			double xx[XRES];
			double yy[3][XRES];
			double xyz[3];
			for (i = 0; i < XRES; i++) {
				xx[i] = i/(XRES-1.0);
				csamp_interp(&asgrey, xyz, xx[i]);
				for (j = 0; j < 3; j++)
					yy[j][i] = xyz[j];
			}
			printf("Measured neutral axis XYZ\n",k);
			do_plot(xx,yy[0],yy[1],yy[2],XRES);
			#undef XRES
		}
#endif

		/* Check out the accuracy of the results: */
		{
			double ctwh[3];		/* Current target white */
			icmXYZNumber ctwN;	/* Same as above as XYZNumber */
			double brerr;		/* Brightness error */
			double cterr;		/* Color temperature delta E */
			double mnerr;		/* Maximum neutral error */
			double mnv = 0.0;	/* Value where maximum error is */
			double anerr;		/* Average neutral error */
			double lab1[3], lab2[3];
			
			/* Brightness */
			brerr = asgrey.s[asgrey.no-1].XYZ[1] - x.twh[1];
		
			/* Compensate for brightness error */
			for (j = 0; j < 3; j++)
				ctwh[j] = x.twh[j] * asgrey.s[asgrey.no-1].XYZ[1]/x.twh[1];
			icmAry2XYZ(ctwN, ctwh);		/* Need this for Lab conversions */
			
			/* Color temperature error */
			icmXYZ2Lab(&ctwN, lab1, ctwh);		/* Should be 100,0,0 */
			icmXYZ2Lab(&ctwN, lab2, asgrey.s[asgrey.no-1].XYZ);
			cterr = icmLabDE(lab1, lab2);

			/* check delta E of all the sample points */
			/* We're checking against our given brightness and */
			/* white point target. */
			mnerr = anerr = 0.0;
			for (i = 0; i < asgrey.no; i++) {
				double err = asgrey.s[i].de;

//printf("RGB %.3f -> Lab %.2f %.2f %.2f, target %.2f %.2f %.2f, DE %f\n",
//asgrey.s[i].v, lab2[0], lab2[1], lab2[2], lab1[0], lab1[1], lab1[2], err);
				if (err > mnerr) {
					mnerr = err;
					mnv = asgrey.s[i].v;
				}
				anerr += err;
			}
			anerr /= (double)asgrey.no;

			if (verb || it >= mxits) {
				if (it >= mxits)
					printf("Verification results:\n");
				printf("Brightness error = %f cd/m^2\n",brerr);
				printf("White point error = %f deltaE\n",cterr);
				printf("Maximum neutral error (@ %f) = %f deltaE\n",mnv, mnerr);
				printf("Average neutral error = %f deltaE\n",anerr);
			}
		}

		/* Verify loop exit */
		if (it >= (mxits + nver -1)) {
			break;
		}

		/* Convert our test points into calibration curves. */
		/* The call to reinit_csamp() will then convert the */
		/* curves back to current test point values. */
		/* This applies some level of cohesion between the test points, */
		/* as well as forcing monotomicity */
		if (it < mxits) {
			mcvco *sdv[3];				/* Scattered data for mcv */

			for (j = 0; j < 3; j++) {
				if ((sdv[j] = malloc(sizeof(mcvco) * asgrey.no)) == NULL) {
					dr->del(dr);
					error ("Malloc of scattered data points failed");
				}
			}

			if (verb)
				printf("Computing update to calibration curves...\n");

			for (i = 0; i < asgrey.no; i++) {
				/* Use fixed rgb's */
				for (j = 0; j < 3; j++) {
					sdv[j][i].p = asgrey.s[i].v;
					sdv[j][i].v = asgrey.s[i].rgb[j];
					sdv[j][i].w = 1.0;
				}
			}
			if (x.nat)		/* Make curve go thought white if possible */
				sdv[0][rsteps-1].w = sdv[1][rsteps-1].w = sdv[2][rsteps-1].w = 10.0;

			for (j = 0; j < 3; j++)
				x.rdac[j]->fit(x.rdac[j], 0, fitord, sdv[j], asgrey.no, 1.0);

			/* Make sure that if we are using native brightness and white point, */
			/* that the curves go to a perfect 1.0 ... */
			if (x.nat) {
				for (j = 0; j < 3; j++)
					x.rdac[j]->force_1(x.rdac[j], 1.0);
			}

			for (j = 0; j < 3; j++)
				free(sdv[j]);
#ifdef DEBUG_PLOT
			/* Plot the current curves */
			{
				#define	XRES 255
				double xx[XRES];
				double y1[XRES];
				double y2[XRES];
				double y3[XRES];
				double rgb[3];
				for (i = 0; i < XRES; i++) {
					double drgb[3], rgb[3];
					xx[i] = i/(XRES-1.0);
					rgb[0] = rgb[1] = rgb[2] = xx[i];
					for (j = 0; j < 3; j++)
						drgb[j] = x.rdac[j]->interp(x.rdac[j], rgb[j]);
					y1[i] = drgb[0];
					y2[i] = drgb[1];
					y3[i] = drgb[2];
				}
				printf("Current ramdac curves\n");
				do_plot(xx,y1,y2,y3,XRES);
				#undef XRES
			}
#endif
		}
	}	/* Next refine/verify loop */

	free_alloc_csamp(&asgrey);		/* We're done with test points */
	dr->del(dr);	/* Now we're done with test window */

	/* Write out the resulting calibration file */
	if (verify != 2) {
		int calres = 256;				/* 256 steps in calibration */
		cgats *ocg;						/* output cgats structure */
		time_t clk = time(0);
		struct tm *tsp = localtime(&clk);
		char *atm = asctime(tsp);		/* Ascii time */
		cgats_set_elem *setel;			/* Array of set value elements */
		int ncps;						/* Number of curve parameters */
		double *cps[3];					/* Arrays of curve parameters */
		char *bp = NULL, buf[100];		/* Buffer to sprintf into */

		ocg = new_cgats();				/* Create a CGATS structure */
		ocg->add_other(ocg, "CAL"); 	/* our special type is Calibration file */

		ocg->add_table(ocg, tt_other, 0);	/* Add a table for RAMDAC values */
		ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Device Calibration State",NULL);
		ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll dispcal", NULL);
		atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
		ocg->add_kword(ocg, 0, "CREATED",atm, NULL);

		ocg->add_kword(ocg, 0, "DEVICE_CLASS","DISPLAY", NULL);

		/* Put the target parameters in the CGATS file too */
		if (dtype != 0)
			ocg->add_kword(ocg, 0, "DEVICE_TYPE", dtype == 1 ? "CRT" : "LCD", NULL);
		
		if (wpx == 0.0 && wpy == 0.0 && temp == 0.0 && tbright == 0.0)
			ocg->add_kword(ocg, 0, "NATIVE_TARGET_WHITE","", NULL);

		sprintf(buf,"%f %f %f", x.twh[0], x.twh[1], x.twh[2]);
		ocg->add_kword(ocg, 0, "TARGET_WHITE_XYZ",buf, NULL);

		sprintf(buf,"%f", gamma);
		ocg->add_kword(ocg, 0, "TARGET_GAMMA",x.gammat == 0 ? buf : x.gammat == 1 ? "L_STAR" : "sRGB", NULL);
		sprintf(buf,"%f", bkcorrect);
		ocg->add_kword(ocg, 0, "BLACK_POINT_CORRECTION", buf, NULL);

		/* Write rest of setup */
	    switch (quality) {
			case -3:				/* Test value */
				bp = "ultra low";
				break;
			case -2:				/* Very low */
				bp = "very low";
				break;
			case -1:				/* Low */
				bp = "low";
				break;
			case 0:					/* Medum */
				bp = "medium";
				break;
			case 1:					/* High */
				bp = "high";
				break;
			case 2:					/* Ultra */
				bp = "ultra high";
				break;
			default:
				error("unknown quality level %d",quality);
		}
		ocg->add_kword(ocg, 0, "QUALITY",bp, NULL);

		ocg->add_field(ocg, 0, "RGB_I", r_t);
		ocg->add_field(ocg, 0, "RGB_R", r_t);
		ocg->add_field(ocg, 0, "RGB_G", r_t);
		ocg->add_field(ocg, 0, "RGB_B", r_t);

		if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * 4)) == NULL)
			error("Malloc failed!");

		/* Write the video lut curve values */
		for (i = 0; i < calres; i++) {
			double vv, rgb[3];

#ifdef __APPLE__
			gcc_bug_fix(i);
#endif
			vv = i/(calres-1.0);
			for (j = 0; j < 3; j++) {
				double cc;
				cc = x.rdac[j]->interp(x.rdac[j], vv);
				if (cc < 0.0)
					cc = 0.0;
				else if (cc > 1.0)
					cc = 1.0;
				rgb[j] = cc;
			}

			setel[0].d = vv;
			setel[1].d = rgb[0];
			setel[2].d = rgb[1];
			setel[3].d = rgb[2];

			ocg->add_setarr(ocg, 0, setel);
		}

		free(setel);

		/* Write some of the device model information to a second */
		/* table, so that we can update the calibration latter on without */
		/* having to read R,G & B curves. */

		ocg->add_table(ocg, tt_other, 0);	/* Add a second table for setup and model */
		ocg->add_kword(ocg, 1, "DESCRIPTOR", "Argyll Calibration options and model",NULL);
		ocg->add_kword(ocg, 1, "ORIGINATOR", "Argyll dispcal", NULL);
		atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
		ocg->add_kword(ocg, 1, "CREATED",atm, NULL);


		/* Write device model curves */
		ocg->add_field(ocg, 1, "R_P", r_t);
		ocg->add_field(ocg, 1, "G_P", r_t);
		ocg->add_field(ocg, 1, "B_P", r_t);

		if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * 3)) == NULL)
			error("Malloc failed!");

		ncps = -1;
		for (i = 0; i < 3; i++) {
			int nn;
			nn = x.dcvs[i]->get_params(x.dcvs[i], &cps[i]);
			if (ncps != -1 && ncps != nn)
				error("Expect device model linearisation curves to have the same order");
			ncps = nn;
		}

		for (i = 0; i < ncps; i++) {
			setel[0].d = cps[0][i];
			setel[1].d = cps[1][i];
			setel[2].d = cps[2][i];
			ocg->add_setarr(ocg, 1, setel);
		}

		for (i = 0; i < 3; i++)
			free(cps[i]);
		free(setel);

		if (ocg->write_name(ocg, outname))
			error("Write error : %s",ocg->err);

		ocg->del(ocg);		/* Clean up */
	}

	/* Update the ICC file with the new LUT curves */
	if (verify != 2 && doupdate && iccname[0] != '\000') {
		icmFile *ic_fp;
		icc *icco;
		int j, i;
		icmVideoCardGamma *wo;

		if ((icco = new_icc()) == NULL)
			error ("Creation of ICC object to read profile '%s' failed",iccname);

		/* Open up the profile for reading */
		if ((ic_fp = new_icmFileStd_name(iccname,"r")) == NULL)
			error ("Can't open file '%s'",iccname);

		/* Read header etc. */
		if ((rv = icco->read(icco,ic_fp,0)) != 0)
			error ("Reading profile '%s' failed with %d, %s",iccname, rv,icco->err);

		/* Read every tag */
		if (icco->read_all_tags(icco) != 0) {
			error("Unable to read all tags from '%s': %d, %s",iccname, icco->errc,icco->err);
		}

		ic_fp->del(ic_fp);

		wo = (icmVideoCardGamma *)icco->read_tag(icco, icSigVideoCardGammaTag);
		if (wo == NULL)
			error("Can't find VideoCardGamma tag in file '%s': %d, %s",
			      iccname, icco->errc,icco->err);

		wo->tagType = icmVideoCardGammaTableType;
		wo->u.table.channels = 3;             /* rgb */
		wo->u.table.entryCount = 256;         /* full lut */
		wo->u.table.entrySize = 2;            /* 16 bits */
		wo->allocate((icmBase*)wo);
		for (j = 0; j < 3; j++) {
			for (i = 0; i < 256; i++) {
				double cc, vv = i/(256-1.0);
				cc = x.rdac[j]->interp(x.rdac[j], vv);
				((unsigned short*)wo->u.table.data)[256 * j + i] = (int)(cc * 65535.0 + 0.5);
			}
		}

		/* Open up the profile again writing */
		if ((ic_fp = new_icmFileStd_name(iccname,"w")) == NULL)
			error ("Can't open file '%s' for writing",iccname);

		if ((rv = icco->write(icco,ic_fp,0)) != 0)
			error ("Write to file '%s' failed: %d, %s",iccname, rv,icco->err);

		ic_fp->del(ic_fp);
		icco->del(icco);
	}

	if (verify != 2) {
		for (j = 0; j < 3; j++)
			x.rdac[j]->del(x.rdac[j]);
	
		for (k = 0; k < 3; k++)
			x.dcvs[k]->del(x.dcvs[k]);
	}

	free_a_disppath(disp);

	return 0;
}


