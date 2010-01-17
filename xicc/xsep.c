
/* 
 * International Color Consortium color transform expanded support
 *
 * Author:  Graeme W. Gill
 * Date:    1/5/01
 * Version: 1.00
 *
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 */

/*
 * This class handles the creation of an optimised color separation.
 * This is a middle layer that has a repertoire of ink separation
 * rules that can be applied to suite the particular colorspace and
 * wishes of the user. 
 *
 * !!! Under Development - not finished yet !!!
 */

#include "numlib.h"
#include "rspl_imp.h"
#include "xicc.h"
#include "xsep.h"

/* Return the separation given the pseudo-device values */
/* (Assumes that we have loaded or computed the separation into a */
/*  separation object) */
static void xsep_lookup(
xsep *p,
double out[MXDO],	/* Actual device output (ie. CMYKOG) */
double in[MXDI]		/* Pseudo device CMYx input */
) {
	co tc;
	int i;

	for (i = 0; i < p->pdi; i++)
		tc.p[i] = in[i];

	p->sep->interp(p->sep, &tc); 

	for (i = 0; i < p->ddi; i++)
		out[i] = tc.v[i];
}

static void xsep_del(
xsep *p
) {

	if (p->sep != NULL)
		p->sep->del(p->sep);
		
	free(p);
}


#ifdef NEVER	/* Not used yet */

/* The optimisation function passed to opt_rspl */
/* Returns value is the "error" for this point. */
static double optfunc(
void *fdata,
double *inout,	/* Pointers to fdi+pdi+adi values for the grid point being optimised. */
				/* corresponding to  */
double *surav,	/* Pointers to fdi+pdi values which are the average of the */
				/* neighbors of this grid point. Pointer will NULL if this */
				/* is a surface grid point. */
int first,		/* Flag, NZ if this is the first optimisation of this point. */
double cw		/* The (grid resolution) curvature weighting factor */
) {
	xsep *p = (xsep *)fdata;

	/* If first time at this grid resolution, compute needed locus extreme values */
	/* and place them in the additional values allocation */
	if (first) {
	}

	if (surav == NULL) {	/* Must be a surface point */
		/* Optimise to maximise distance along vector direction */
		/* together with device value smoothness, inking rule */
		/* targets and ink limit. */

		/* Compute surface smoothness value */

	} else {				/* Must be an interior point */
		/* Compute Lab target as mean of surrounding targets, */
		/* and optimise towards it, together with device value */
		/* smoothness, inking rule targets and ink limit. */

		/* Compute interior smoothness value */
	}

	return 0.0;
}

#endif /* NEVER */


/* default target vectors (Labx) for CMYx input */
static double def_tvect[16][4] = {
	{ 100.0,   0.0,    0.0,  100.0 },		/* 0 0 0 0 */
	{  50.0, -60.0, -100.0,  100.0 },		/* C 0 0 0 */
	{  50.0, 100.0,    0.0,  100.0 },		/* 0 M 0 0 */
	{  25.0,  60.0, -100.0,  100.0 },		/* C M 0 0 */
	{  90.0,   0.0,  100.0,  100.0 },		/* 0 0 Y 0 */
	{  50.0,-100.0,   60.0,  100.0 },		/* C 0 Y 0 */
	{  50.0, 100.0,  100.0,  100.0 },		/* 0 M Y 0 */
	{   0.0,   0.0,    0.0,  100.0 },		/* C M Y 0 */
	{ 100.0,   0.0,    0.0,    0.0 },		/* 0 0 0 x */
	{  50.0, -60.0, -100.0,    0.0 },		/* C 0 0 x */
	{  50.0, 100.0,    0.0,    0.0 },		/* 0 M 0 x */
	{  25.0,  60.0, -100.0,    0.0 },		/* C M 0 x */
	{  90.0,   0.0,  100.0,    0.0 },		/* 0 0 Y x */
	{  50.0,-100.0,   60.0,    0.0 },		/* C 0 Y x */
	{  50.0, 100.0,  100.0,    0.0 },		/* 0 M Y x */
	{   0.0,   0.0,    0.0,    0.0 }		/* C M Y x */
};

/* Create the optimised separation - return NULL on error */
/* We assume that we want to create a separation from pseudo-device */
/* channels CMY' to the real device channels. There may also be auxiliary */
/* channels that come after CMY' and Lab. The schident[] is used to align */
/* the pseudo chanels with their real counterparts. If there are no */
/* counterparts, then xsep will use a default. The schident[3] is assumed */
/* to be the chanel to use for the black (C+M+Y) chanel as an alternative */
/* to the C+M+Y direction. */
xsep *new_xsep(
int pdi,			/* pseudo device CMY'/CMYK', and target Lab/LabL (PCS) dimensions */
int ddi, 			/* Device Dimensions/Channels */
int (*dev2lab) (void *d2lcntx, double *out, double *in),	/* Device to Lab callback */
void *d2lcntx,		/* dev2lab callback context */
double (*dev2ink) (void *d2icntx, double *out, double *in),	/* Device to ink callback */
					/* Return 2 or 3 totals: Total ink, CMY sup. group, K sup. */
void *d2icntx,		/* dev2ink callback context */
double totlimit,	/* Value not to be exceeded by dev2ink() Total ink, 0.0 .. N.0 */
double smoothness,	/* Device value smoothness, nominally 1.0 */
int schident[3],	/* Standard channel ident for psudo CMY' chanels. */
					/* -1 = no mapping else index of corresponding device channel. */
icxScat cat[MXDO]	/* Device channel control and function category */
){
	xsep *p;		/* this */
	rspl *sep;		/* Mapping we will create */
	int i, j;
//	double **vdata;	/* pdi^2 array of function, target and additional values to init */
	double *vdatap[16];

	/* Sanity check */
	if (pdi != 3 && pdi != 4) {
		return NULL;
	}

	if (ddi < 1 || ddi > MXDO) {
		return NULL;
	}

	if ((p = (xsep *) calloc(1,sizeof(xsep))) == NULL)
		return NULL;

	p->lookup      = xsep_lookup;
	p->del         = xsep_del;

	/* We need to setup the initial CMYx' colorant target vectors */
	/* This will be from the device counterparts (ie. like CMYK), */
	/* or using a default set of directions. */
	
	for (i = 0; i < pdi; i++) {
		if (schident[i] < 0)
			break;				/* Not all defined */
		for (j = 0; j < i; j++)
			if (schident[j] == schident[i])
				break;			/* Two are the same */
		if (j < i)
			break;
	}
	if (i < pdi) {				/* Fall back on default */
		for (i = 0; i < (1 << pdi); i++) {
			vdatap[i] = def_tvect[i];	// ~~~99
		}
	} else {
	}

	/* !!! Stuff goes here !!! */
	/* This will rely on the opt_rspl method of the rspl class */

	/* Create the rspl that maps CMYx -> Device values */
	if ((sep = new_rspl(RSPL_NOFLAGS, pdi, ddi)) == NULL) {
		free(p);
		return NULL;
	}

#ifdef NEVER
	/* Set values by multi-grid optimisation using the provided function. */
	/* The input values to the rspl are typically the CMY' or CMYx' color space. */
	/* The output values are the device values */
	/* The target values are typically Lab or Labx (vector dir on surface) */
	/* The additional value are typically the device auxiliary chanel locus limits, */
	/* computed on the first call at a given grid resolution. */
	/* The func needs to compute the smoothness value, and weight it */
	/* by the curvature value. */
	int sep->opt_rspl(
		sep,			/* this */
		int flags,		/* Combination of flags */
		int pdi,		/* Dimensionality of target data (typically Lab or Labx) */
		int adi,		/* Additional grid point data allowance */
		double **vdata,	/* pdi^2 array of function, target and additional values to init */
						/* array corners with. */
		double (*func)(void *fdata, double *inout, double *surav, int first, double cw),
						/* Optimisation function */
		void *fdata,	/* Opaque data needed by function */
		datai glow,		/* Grid low scale - NULL = default 0.0 */
		datai ghigh,	/* Grid high scale - NULL = default 1.0 */
		int gres,		/* Spline grid resolution */
		datao vlow,		/* Data value low normalize, NULL = default 0.0 */
		datao vhigh		/* Data value high normalize - NULL = default 1.0 */
	);
#endif /* NEVER */

	return p;
}























