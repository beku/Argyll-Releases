
#ifndef MCV_H
#define MCV_H

/* 
 * Argyll Color Correction System
 * Monotonic curve class for display calibration.
 *
 * Author: Graeme W. Gill
 * Date:   30/10/2005
 *
 * Copyright 2005 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the LICENCE.TXT file for licencing details.
 *
 * This is based on the monotonic curve equations used elsewhere,
 * but currently intended to support the display calibration process.
 * moncurve is not currently general, missing:
 *
 *  input scaling
 *  output scaling
 *
 * The nominal input and output ranges are 0.0 to 1.0,
 */

/* A test patch value */
typedef struct {
	double p;		/* Position */
	double v;		/* Value */
} mcvco;

struct _mcv {

  /* Public: */
	void (*del)(struct _mcv *p);

	/* Fit the curve to the given points */
	void (*fit) (struct _mcv *p,
	            int verb,		/* Vebosity level, 0 = none */
	            int order,		/* Number of curve orders, 1..MCV_MAXORDER */
		        mcvco *d,		/* Array holding scattered initialisation data */
		        int ndp,		/* Number of data points */
	            double smooth	/* Degree of smoothing, 1.0 = normal */			
	);

	/* Scale the the output so that the value for input 1.0, */
	/* is the given value. */
	void (*force_scale) (struct _mcv *p,
	            double target	/* Target output value */
	);

	/* Translate a value through the curve */
	double (*interp) (struct _mcv *p,
	                  double in);	/* Input value */

	/* Translate a value through backwards the curve */
	double (*inv_interp) (struct _mcv *p,
	                  double in);	/* Input value */

  /* Private: */
	int verb;				/* Verbose */
	int luord;				/* Lookup order including scale */
	double *pms;			/* Allocated curve parameters */
	double *upms;			/* Pointer to parameters to use */
	double *dv;				/* Work space for dv's */

	mcvco *d;				/* Array holding scattered initialisation data */
	int ndp;				/* Number of data points */
	double smooth;			/* Smoothing factor */

}; typedef struct _mcv mcv;

/* Create a new, uninitialised mcv */
mcv *new_mcv(void);

#endif /* MCV */

