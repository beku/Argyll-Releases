
#ifndef XFIT2_H
#define XFIT2_H

/* 
 * Clut per channel curve fitting
 *
 * Author:  Graeme W. Gill
 * Date:    27/5/2007
 * Version: 1.00
 *
 * Copyright 2000 - 2007 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * Based on the rspl.c and xlut.c code.
 */

#define MXLUORD 10		/* Shaper harmonic orders to use */
#define PCSD 3			/* PCS dimensions */

#define MXPARMS (((1 << MXDI) * MXDO) + (MXDI + MXDO) * (MXLUORD))

/* Optimisation mask legal combinations */
typedef enum {
	oc2_i   = 1,			/* Input */
	oc2_s   = 2,			/* Input sub-grid */
	oc2_o   = 4,			/* Output */
	oc2_iso = 1 + 2 + 4,	/* Input, Sub-grid and Output */
} optcomb2;

/* Flag values */
#define XFIT_VERB       0x0001		/* Verbose output during fitting */

#define XFIT_FM_INPUT   0x0002		/* Use input space for fit metric (default is output) */

#define XFIT_FM_XYZ     0x0010		/* Metric space is XYZ */
#define XFIT_FM_LAB     0x0020		/* Metric space is Lab */
#define XFIT_FM_DEV     0x0030		/* Metric space is Device */
#define XFIT_FM_LU      0x0040		/* Metric space should use dev->Lab lookup */
#define XFIT_FM_MASK    0x00F0		/* Mask for metric handling */

#define XFIT_IN_ZERO    0x0100		/* Adjust input curves 1 & 2 for zero (~~not implemented) */
#define XFIT_OUT_ZERO   0x0200		/* Adjust output curves 1 & 2 for zero (~~ will delete) */


#define XFIT_OPTGRID_RANGE 0x0080	/* Optimize inner grid around usaed range */

/* A set of curve test points for a single axis and single grid span */
typedef struct {
	int ns, ntp;	/* Number of sets, number of points per set */
	double x0, x1;	/* Span */
	co *lup;		/* in->out reference point values lup[ns * np] */
	co *lupd;		/* in'->out' point values lupd[ns * np] */
} xfit2_ctp;

/* Context for optimising input and output luts */
struct _xfit2 {
	int verb;				/* Verbose */
	int flags;				/* Behaviour flags */
	int di, fdi;			/* Dimensionaluty of input and output */
	optcomb2 tcomb;			/* Target elements to fit */
	int gres[MXDI];			/* clut resolutions being optimised for */

	rspl *iaxs[MXDI];		/* rspl's for each input coord */
	double imin[MXDI], imax[MXDI];	/* Range iaxs use on input */
	int orthres;			/* iaxs rspl resolution orthogonal to input coord */

	void *cntx2;			/* Context of callback */
	double (*to_de2)(void *cntx, double *in1, double *in2);	
							/* callback to convert in or out value to fit metric squared */

	int iluord[MXDI];		/* Input Shaper order actualy used (must be <= MXLUORD) */
	int sluord[MXDI];		/* Sub-grid shaper order */
	int oluord[MXDO];		/* Output Shaper order actualy used (must be <= MXLUORD) */
	double in_min[MXDI];	/* Input value scaling minimum */
	double in_max[MXDI];	/* Input value scaling maximum */
	double out_min[MXDO];	/* Output value scaling minimum */
	double out_max[MXDO];	/* Output value scaling maximum */

	int in_off;				/* Input  parameters offset */
	int in_offs[MXDI];		/* Input  parameter offsets for each channel from v[0] */
	int sub_cnt;			/* Sub-grid  parameters count */
	int sub_off;			/* Sub-grid  parameters offset */
	int sub_offs[MXDI];		/* Sub-grid  parameter offsets for each channel from v[0] */
							/* There are sluord[] * (gres[]-1) parameters per channel */
	int in_cnt;				/* Input  parameters count */
	int out_off;			/* Output parameters offset from v[0] */
	int out_offs[MXDO];		/* Output  parameter offsets for each channel from v[0] */
	int out_cnt;			/* Output parameters count */
	int tot_cnt;			/* Total parameter count */

	double *v;				/* Holder for parameters */
							/* Optimisation parameters are layed out:             */
							/*                                                    */
							/* Input curves:, di groups of iluord[e] parameters   */
							/*                                                    */
							/* Output curves:, fdi groups of oluord[f] parameters */
							/*                                                    */

	int nodp;				/* Number of data points                              */
	cow *ipoints;			/* Reference to test points as in->out                */

	/* Optimisation state */
	int optt;				/* 0 = input, 1 = sub-grid, 2 = output curves */  
	int och;				/* In or output channel being optimized */
	int osp;				/* Input sub-grid span being optimised */
	int opt_off;			/* Optimisation parameters offset from v[0] */
	int opt_cnt;			/* Optimisation parameters count */
	double *wv;				/* Parameters being optimised */
	double *sa;				/* Search area */
	double *uerrv;			/* Uniform error values for each span */

	xfit2_ctp **ctp; 		/* Pre-computed test/ref points ctp[di][gres-1] */

	/* Methods */
	void (*del)(struct _xfit2 *p);

	/* Do the fitting. Return nz on error */ 
	int (*fit)(
		struct _xfit2 *p, 
		int flags,				/* Flag values */
		int di,					/* Input dimensions */
		int fdi,				/* Output dimensions */
		int nodp,				/* Number of data points */
		cow *ipoints,			/* Array of data points to fit - referece taken */
		int gres[MXDI],			/* clut resolutions being optimised for */
		double in_min[MXDI],	/* Input value scaling/domain minimum */
		double in_max[MXDI],	/* Input value scaling/domain maximum */
		double out_min[MXDO],	/* Output value scaling/range minimum */
		double out_max[MXDO],	/* Output value scaling/range maximum */
		int iord[],				/* Order of input shaper curve for each dimension */
		int sord[],				/* Order of input sub-grid shaper curve for each dimension */
		int oord[],				/* Order of output shaper curve for each dimension */
		optcomb2 tcomb,			/* Target elements to fit. */
		void *cntx2,			/* Context of callback */
		double (*to_de2)(void *cntx, double *in1, double *in2)	
								/* callback to convert in or out value to fit metric squared */
	);

	/* Lookup a value though an input curve */
	double (*incurve)(struct _xfit2 *p, double in, int chan);

	/* Inverse Lookup a value though an input curve */
	double (*invincurve)(struct _xfit2 *p, double in, int chan);

	/* Lookup a value though an output curve */
	double (*outcurve)(struct _xfit2 *p, double in, int chan);

	/* Inverse Lookup a value though an output curve */
	double (*invoutcurve)(struct _xfit2 *p, double in, int chan);

}; typedef struct _xfit2 xfit2;

xfit2 *new_xfit2();

#endif /* XFIT2_H */



































