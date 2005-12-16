/* Test support code */
/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include "stdio.h"
#include "stdlib.h"
#include "refi.h"

/* Callbackes used to setup rspl's */
static void inputlu(
void *cbctx,
double *out,
double *in
) {
	refi *r = (refi *)cbctx;

	*out = r->input_curve(r->cntx, r->chan, *in);
}

static void clutlu(
void *cbctx,
double *out,
double *in
) {
	refi *r = (refi *)cbctx;

	r->md_table(r->cntx, out, in);
}

static void outputlu(
void *cbctx,
double *out,
double *in
) {
	refi *r = (refi *)cbctx;

	*out = r->output_curve(r->cntx, r->chan, *in);
}


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
) {
	refi *r;
	int e;
	int gres[MXDI];

	if ((r = (refi *)malloc(sizeof(refi))) == NULL) {
		fprintf(stderr,"Malloc of refi failed\n");
		exit (-1);
	}

	r->id = id;
	r->od = od;
	r->inres = inres;
	r->clutres = clutres;
	r->outres = outres;
	r->input_curve  = input_curve;
	r->md_table     = md_table;
	r->output_curve = output_curve;
	r->cntx		    = cntx;

	/* Create some input interpolations */
	for (e = 0; e < id; e++) {
		if ((r->in[e] = new_rspl(1, 1)) == NULL) {
			fprintf(stderr,"new_rspl failed\n");
			exit (-1);
		}
		r->chan = e;
		r->in[e]->set_rspl(r->in[e], 0, (void *)r, inputlu, NULL, NULL, &inres, NULL, NULL);
	}

	/* Clut */
	if ((r->clut = new_rspl(id, od)) == NULL) {
		fprintf(stderr,"new_rspl failed\n");
		exit (-1);
	}
	for (e = 0; e < id; e++)
		gres[e] = clutres;
	r->clut->set_rspl(r->clut, 0, (void *)r, clutlu, NULL, NULL, gres, NULL, NULL);
	
	/* Create some output interpolations */
	for (e = 0; e < od; e++) {
		if ((r->out[e] = new_rspl(1, 1)) == NULL) {
			fprintf(stderr,"new_rspl failed\n");
			exit (-1);
		}
		r->chan = e;
		r->out[e]->set_rspl(r->out[e], 0, (void *)r, outputlu, NULL, NULL, &outres, NULL, NULL);
	}

	return r;
}

/* Run an interpolation through an input table */
double refi_input(
void *cntx,
int ch,
double in_val
) {
	refi *r = (refi *)cntx;
	co vals;				/* Input and output values */

	vals.p[0] = in_val;
	r->in[ch]->interp(r->in[ch], &vals);
	return vals.v[0];
}

/* Run an interpolation through an clut table */
void refi_clut(
void *cntx,
double *out_vals,
double *in_vals
) {
	refi *r = (refi *)cntx;
	int e;
	co vals;				/* Input and output values */

	for (e = 0; e < r->id; e++) 
		vals.p[e] = in_vals[e];
	r->clut->interp(r->clut, &vals);
	for (e = 0; e < r->od; e++)
		out_vals[e] = vals.v[e];
}

/* Run an interpolation through an output table */
double refi_output(
void *cntx,
int ch,
double in_val
) {
	refi *r = (refi *)cntx;
	co vals;				/* Input and output values */

	vals.p[0] = in_val;
	r->out[ch]->interp(r->out[ch], &vals);
	return vals.v[0];
}

void
refi_free(
refi *r
) {
	int e;

	for (e = 0; e < r->id; e++) {
		r->in[e]->del(r->in[e]);
	}

	r->clut->del(r->clut);

	for (e = 0; e < r->od; e++) {
		r->out[e]->del(r->out[e]);
	}
}
























