#ifndef NUMSUP_H
#define NUMSUP_H

/* Numerical routine general support declarations */
/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

/* Some default math limits */
#ifndef DBL_EPSILON
#define DBL_EPSILON     2.2204460492503131e-016		/* 1.0+DBL_EPSILON != 1.0 */
#endif
#ifndef DBL_MIN
#define DBL_MIN         2.2250738585072014e-308		/* IEEE 64 bit min value */
#endif
#ifndef DBL_MAX
#define DBL_MAX         1.7976931348623158e+308		/* IEEE 64 bit max value */
#endif

extern char *error_program;
extern FILE *error_out;
extern FILE *warn_out;
extern FILE *verbose_out;
extern int verbose_level;

#include <stdarg.h>
extern void error(char *fmt, ...), warning(char *fmt, ...);

/* Numerical recipes vector/matrix support functions */
float *fvector(int nl,int nh);
float *fvectorz(int nl,int nh);
void free_fvector(float *v,int nl,int nh);

float **fmatrix(int nrl, int nrh, int ncl, int nch);
float **fmatrixz(int nrl, int nrh, int ncl, int nch);
void free_fmatrix(float **m, int nrl, int nrh, int ncl, int nch);

double *dvector(int nl,int nh);
double *dvectorz(int nl,int nh);
void free_dvector(double *v,int nl,int nh);

double **dmatrix(int nrl, int nrh, int ncl, int nch);
double **dmatrixz(int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);

double **dhmatrix(int nrl, int nrh, int ncl, int nch);
double **dhmatrixz(int nrl, int nrh, int ncl, int nch);
void free_dhmatrix(double **m, int nrl, int nrh, int ncl, int nch);

void copy_dmatrix(double **dst, double **src, int nrl, int nrh, int ncl, int nch);

double **convert_dmatrix(double *a,int nrl,int nrh,int ncl,int nch);
void free_convert_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);

int *ivector(int nl,int nh);
int *ivectorz(int nl,int nh);
void free_ivector(int *v,int nl,int nh);

int **imatrix(int nrl,int nrh,int ncl,int nch);
int **imatrixz(int nrl,int nrh,int ncl,int nch);
void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);

#endif /* NUMSUP_H */
