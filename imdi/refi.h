
/* Reference floating point interpolator constructed out of rspl's */
/* This provides imdi functionality in floating point */
/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include "../rspl/rspl.h"

/* ------------------------------------------------ */

typedef struct {
	int id, od;		/* Input and output dimensions */
	int inres;		/* Desired input table resolution */
	int clutres;	/* Desired clut table resolution */
	int outres;		/* Desired output table resolution */
	rspl *in[MXDI];
	rspl *clut;
	rspl *out[MXDO];

	double (*input_curve) (void *cntx, int ch, double in_val);
	void   (*md_table)    (void *cntx, double *out_vals, double *in_vals);
	double (*output_curve)(void *cntx, int ch, double in_val);
	void *cntx;		/* Context to callbacks */
	int chan;		/* Current callback channel */
} refi;

refi *new_refi(
	int id,			/* Number of input dimensions */
	int od,			/* Number of output dimensions */
	int inres,		/* Desired input table resolution */
	int clutres,	/* Desired clut table resolution */
	int outres,		/* Desired output table resolution */

	/* Callbacks to lookup the table values */
	double (*input_curve) (void *cntx, int ch, double in_val),
	void   (*md_table)    (void *cntx, double *out_vals, double *in_vals),
	double (*output_curve)(void *cntx, int ch, double in_val),
	void *cntx		/* Context to callbacks */
);

void refi_free(refi *r);


/* Component interpolations */
double refi_input(void *cntx, int ch, double in_val);
void refi_clut(void *cntx, double *out_vals, double *in_vals);
double refi_output(void *cntx, int ch, double in_val);

