
#ifndef XFIT1_H
#define XFIT1_H

/* 
 * International Color Consortium color transform expanded support
 *
 * Author:  Graeme W. Gill
 * Date:    28/6/00
 * Version: 1.00
 *
 * Copyright 2000 - 2007 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 */

#define MXLUORD 10		/* Shaper harmonic orders to use */
#define PCSD 3			/* PCS dimensions */

#define MXPARMS (((1 << MXDI) * MXDO) + (MXDI + MXDO) * (MXLUORD))

/* Optimisation mask legal combinations */
/* (Reserve flag 2 for oc2_s) */
typedef enum {
	oc_i   = 1,			/* Input */
	oc_m   = 8,			/* Matrix */
	oc_im  = 1 + 8,		/* Input and matrix */
	oc_o   = 4,			/* Output */
	oc_io  = 1 + 4,		/* Input and Output */
	oc_mo  = 8 + 4,		/* Matrix and output */
	oc_imo = 1 + 8 + 4	/* Input, matrix and output */
} optcomb;

/* Flag values */
#define XFIT_VERB     0x0001		/* Verbose output during fitting */

#define XFIT_FM_INPUT   0x0002		/* Use input space for fit metric (default is output) */

#define XFIT_FM_XYZ     0x0010		/* Metric space is XYZ */
#define XFIT_FM_LAB     0x0020		/* Metric space is Lab */
#define XFIT_FM_DEV     0x0030		/* Metric space is Device */
#define XFIT_FM_LU      0x0040		/* Metric space should use dev->Lab lookup */
#define XFIT_FM_MASK    0x00F0		/* Mask for metric handling */

#define XFIT_OUT_DEV  0x0002		/* The output is DEV so scale matrix search by 100 */ 
#define XFIT_IN_ZERO  0x0100		/* Adjust input curves 1 & 2 for zero ~~not implemented */
#define XFIT_OUT_ZERO 0x0200		/* Adjust output curves 1 & 2 for zero */

/* Context for optimising input and output luts */
struct _xfit1 {
	int verb;				/* Verbose */
	int flags;				/* Behaviour flags */
	int di, fdi;			/* Dimensionaluty of input and output */
	optcomb tcomb;			/* Target elements to fit */

	void *cntx2;			/* Context of callback */
	double (*to_de2)(void *cntx2, double *in1, double *in2);	
							/* callback to convert in or out value to fit metric squared */


	int iluord[MXDI];		/* Input Shaper order actualy used (must be <= MXLUORD) */
	int oluord[MXDO];		/* Output Shaper order actualy used (must be <= MXLUORD) */
	double in_min[MXDI];	/* Input value scaling minimum */
	double in_max[MXDI];	/* Input value scaling maximum */
	double out_min[MXDO];	/* Output value scaling minimum */
	double out_max[MXDO];	/* Output value scaling maximum */

	int in_off;				/* Input  parameters offset */
	int in_offs[MXDI];		/* Input  parameter offsets for each channel from v[0] */
	int in_cnt;				/* Input  parameters count */
	int mat_off;			/* Matrix parameters offset from v[0] */
	int mat_cnt;			/* Matrix parameters count */
	int out_off;			/* Output parameters offset from v[0] */
	int out_offs[MXDO];		/* Output  parameter offsets for each channel from v[0] */
	int out_cnt;			/* Output parameters count */
	int tot_cnt;			/* Total parameter count */

	double *v;				/* Holder for parameters */
							/* Optimisation parameters are layed out:         */
							/*                                                */
							/* Input curves:, di groups of iluord[e] parameters */
							/*                                                */
							/* Matrix: fdi groups of 2 ^ di parameters        */
							/*                                                */
							/* Output curves:, fdi groups of oluord[f] parameters */
							/*                                                */

	int nodp;				/* Number of data points                          */
	cow *ipoints;			/* List of test points as in->out                */
	cow *points;			/* Test points as in->Lab */

	/* Optimisation state */
	optcomb opt_msk;		/* Optimisation mask: 3 = i+m, 2 = m, 6 = m+o, 7 = i+m+o */
	int opt_off;			/* Optimisation parameters offset */
	int opt_cnt;			/* Optimisation parameters count */
	int symch;				/* Output channel being adjusted for symetry */
	double *wv;				/* Paraneters being optimised */
	double *sa;				/* Search area */

	/* Methods */
	void (*del)(struct _xfit1 *p);

	/* Do the fitting. Return nz on error */ 
	int (*fit)(
		struct _xfit1 *p, 
		int flags,				/* Flag values */
		int di,					/* Input dimensions */
		int fdi,				/* Output dimensions */
		int nodp,				/* Number of data points */
		cow *points,			/* Array of data points to fit */
		int gres[MXDI],			/* clut resolutions being optimised for */
		double in_min[MXDI],	/* Input value scaling minimum */
		double in_max[MXDI],	/* Input value scaling maximum */
		double out_min[MXDO],	/* Output value scaling minimum */
		double out_max[MXDO],	/* Output value scaling maximum */
		int iord[],				/* Order of input shaper curve for each dimension */
		int dumy[],				/* Dummy array - not used */
		int oord[],				/* Order of output shaper curve for each dimension */
		optcomb tcomb,			/* Target elements to fit. */
		void *cntx2,			/* Context of callback */
		double (*to_de)(void *cntx2, double *in1, double *in2)	
								/* callback to convert in or out value to fit metric squared */
	);

	/* Lookup a value though an input curve */
	double (*incurve)(struct _xfit1 *p, double in, int chan);

	/* Inverse Lookup a value though an input curve */
	double (*invincurve)(struct _xfit1 *p, double in, int chan);

	/* Lookup a value though an output curve */
	double (*outcurve)(struct _xfit1 *p, double in, int chan);

	/* Inverse Lookup a value though an output curve */
	double (*invoutcurve)(struct _xfit1 *p, double in, int chan);

}; typedef struct _xfit1 xfit1;

xfit1 *new_xfit1();

#endif /* XFIT1_H */



































