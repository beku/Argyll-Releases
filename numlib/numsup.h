#ifndef NUMSUP_H
#define NUMSUP_H

/* Numerical routine general support declarations */
/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <math.h>

#ifdef __cplusplus
	extern "C" {
#endif

/* Some default math limits and constants */
#ifndef DBL_EPSILON
#define DBL_EPSILON     2.2204460492503131e-016		/* 1.0+DBL_EPSILON != 1.0 */
#endif
#ifndef DBL_MIN
#define DBL_MIN         2.2250738585072014e-308		/* IEEE 64 bit min value */
#endif
#ifndef DBL_MAX
#define DBL_MAX         1.7976931348623158e+308		/* IEEE 64 bit max value */
#endif
#ifndef DBL_PI
#define DBL_PI         3.1415926535897932384626433832795
#endif

/* Globals used to hold certain information */
extern char *exe_path;			/* Malloce'd - won't be freed on exit (intended leak) */
extern char *error_program;
extern FILE *error_out;
extern FILE *warn_out;
extern FILE *verbose_out;
extern int verbose_level;

extern void set_exe_path(char *arg0);
//extern void error(char *fmt, ...), warning(char *fmt, ...);
extern void (*error)(char *fmt, ...);
extern void (*warning)(char *fmt, ...);
extern void (*verbose)(int level, char *fmt, ...);

/* Numerical recipes vector/matrix support functions */

/* Double */
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


/* Float */
float *fvector(int nl,int nh);
float *fvectorz(int nl,int nh);
void free_fvector(float *v,int nl,int nh);

float **fmatrix(int nrl, int nrh, int ncl, int nch);
float **fmatrixz(int nrl, int nrh, int ncl, int nch);
void free_fmatrix(float **m, int nrl, int nrh, int ncl, int nch);


/* Int */
int *ivector(int nl,int nh);
int *ivectorz(int nl,int nh);
void free_ivector(int *v,int nl,int nh);

int **imatrix(int nrl,int nrh,int ncl,int nch);
int **imatrixz(int nrl,int nrh,int ncl,int nch);
void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);


/* Short */
short *svector(int nl, int nh);
short *svectorz(int nl, int nh);
void free_svector(short *v, int nl, int nh);

short **smatrix(int nrl,int nrh,int ncl,int nch);
short **smatrixz(int nrl,int nrh,int ncl,int nch);
void free_smatrix(short **m,int nrl,int nrh,int ncl,int nch);


/* Cast a double to an IEEE754 single precision value, */
/* in a platform independent fashion. */
unsigned int doubletoIEEE754(double d);

/* Cast a an IEEE754 single precision value to a double, */
/* in a platform independent fashion. */
double IEEE754todouble(unsigned int ip);


#ifdef __cplusplus
	}
#endif

#endif /* NUMSUP_H */
