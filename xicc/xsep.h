#ifndef XSEP_H
#define XSEP_H
/* 
 * International Color Consortium color transform expanded support
 *
 * Author:  Graeme W. Gill
 * Date:    1/5/01
 * Version: 1.10
 *
 * Copyright 2000, 2002 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 */

/*
 * This class handles the creation of an optimised color separation.
 * This may be used when the ICX_OPT_K flag is set when creating an xicc
 * lookup object, or when an explicit separation is needed for
 * a more than 4 color device.
 *
 * This class is independent of other icc or icx classes.
 *
 * The pseudo-device channels are almost always CMY, or CMYK,
 * but this is not essential.
 */

/* The color separation class */
struct _xsep {

  /* Private: */
	rspl *sep;			/* separation rspl */
	int pdi;			/* pseudo device CMYx', and target Labx (PCS) dimensions */
	int ddi; 			/* Device Dimensions */

	/* More private stuff goes here */


  /* Public: */
	void (*del)(struct _xsep *p);

	/* Return the separation given the pseudo-device values */
	void (*lookup)(struct _xsep *p, double *dev, double *pdev);

}; typedef struct _xsep xsep;

/*
 * Black generation curve
 *                         E.P.
 *                          _______   End level
 *                         /
 *                        /
 *                       /
 * Start level    ------/
 *                     S.P.
 */
//typedef struct {
//	double Kstle;		/* K start level at white end (0.0 - 1.0)*/
//	double Kstpo;		/* K start point as prop. of L locus (0.0 - 1.0) */
//	double Kenpo;		/* K end point as prop. of L locus (0.0 - 1.0) */
//	double Kenle;		/* K end level at Black end (0.0 - 1.0) */
//	double Kshap;		/* K transition shape, 0.0-1.0 concave, 1.0-2.0 convex */
//} icxKrule;

/* Per device channel separation control and function category */
typedef enum {
    icxS_CMY_Prim   = 0x0,	/* CMY control group, Primary colorant */
    icxS_CMY_Sup    = 0x1,	/* CMY control group, Suplimental colorant */
    icxS_K_Prim     = 0x2,	/* K control group, Primary colorant */
    icxS_K_Sup      = 0x3	/* K control group, Suplimental colorant */
} icxScat;


/* 
  Setup information :

	Overall want:
		CMY or CMYK
		Ink limit

		Ink or direction for CMYK vectors

	For each ink want:

		K group or other
		Device value aim point

  Operation:
	For CMYK, treat K completely separately,
	and treat CMY->Device-K separately.

	For CMY, need K generation curve.
	Problem is making this compatible with other Argyll black curves -
	really need min and max K map for the whole gamut.
	Either compute this each time we do optimisation, or compute
	it once for each grid point. Complicated by composite light+dark black.
	How does one measure "amount" of black ???a - suplimental lookup
	of kK -> equivalent K ??
*/

/* Create the optimised separation - return NULL on error */
/* We assume that we want to create a separation from pseudo-device */
/* channels CMY' or CMYK' to the real device channels. */
/* The schident[] is used to align the pseudo chanels with their real */
/* counterparts in Lab vector direction. If there are no */
/* counterparts, then xsep will use a default. Since pseudo K' is */
/* only used when there are corresponding device K channel(s), there is */
/* no need to specify its corresponding channel. */
xsep *new_xsep(
	int pdi,			/* pseudo device CMY'/CMYK', and target Lab/LabL (PCS) dimensions */
	int ddi, 			/* Device Dimensions/Channels */
	int (*dev2lab) (void *d2lcntx, double *out, double *in),	/* Device to Lab callback */
	void *d2lcntx,		/* dev2lab callback context */
	double (*dev2ink) (void *d2icntx, double *out, double *in),	/* Device to ink callback */
						/* Return Total ink */
	void *d2icntx,		/* dev2ink callback context */
	double totlimit,	/* Value not to be exceeded by dev2ink() Total ink, 0.0 .. N.0 */
	double smoothness,	/* Device value smoothness, nominally 1.0 */
	int schident[3],	/* Standard channel ident for psudo CMY' chanels. */
						/* -1 = no mapping else index of corresponding device channel. */
	icxScat cat[MXDO]	/* Device channel control and function category */
);

#endif /* XSEP_H */



































