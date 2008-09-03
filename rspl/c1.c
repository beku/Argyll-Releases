
/************************************************/
/* Investigate various curve approximations     */
/************************************************/

/* Discrete regularized spline versions */
/* Standard test + Random testing */

/* Author: Graeme Gill
 * Date:   4/10/95
 * Date:   5/4/96
 *
 * Copyright 1995, 1996 Graeme W. Gill
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#undef DIAG
#undef DIAG2
#undef GLOB_CHECK
#define RES2			/* Do multiple test at various resolutions */
#define AVGDEV 0.0		/* Average deviation of function data */

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "plot.h"
#include "rspl.h"

double lin();
void usage(void);

#define TRIALS 20		/* Number of random trials */
#define SKIP 0		/* Number of random trials to skip */

#define MIN_PNTS 5
#define MAX_PNTS 40

#define MIN_RES 20
#define MAX_RES 2000

double xa[MAX_PNTS];
double ya[MAX_PNTS];

#define XRES 100

#define PNTS 10
#define GRES 400
//#define GRES 800
double t1xa[PNTS] = { 0.2, 0.25, 0.30, 0.35,  0.40,  0.44, 0.48, 0.51,  0.64,  0.75  };
double t1ya[PNTS] = { 0.3, 0.35, 0.4,  0.41,  0.42,  0.46, 0.5,  0.575, 0.48,  0.75  };
co test_points[MAX_PNTS];

double lin(double x, double xa[], double ya[], int n);

void usage(void) {
	fprintf(stderr,"Test 1D rspl interpolation\n");
	fprintf(stderr,"Author: Graeme W. Gill\n");
	fprintf(stderr,"usage: c1 [options]\n");
	fprintf(stderr," -s smooth     Use given smoothness (default 1.0)\n");
	fprintf(stderr," -x            Use extra fitting\n");
	exit(1);
}

int main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	int i,j, n;
	double x;
	double xx[XRES];
	double yy[6][XRES];
	rspl *rss;		/* incremental solution version */
	datai low,high;
	int gres[MXDI];
	double smooth = 1.0;
	int extra = 0;
	double avgdev[MXDO];

	low[0] = 0.0;
	high[0] = 1.0;
	avgdev[0] = AVGDEV;

	error_program = "c1";

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')	{	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage();

			/* smoothness */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				fa = nfa;
				if (na == NULL) usage();
				smooth = atof(na);
			}

			else if (argv[fa][1] == 'x' || argv[fa][1] == 'X') {
				extra = 1;

			} else 
				usage();
		} else
			break;
	}

	for (n = 0; n < TRIALS; n++) {
		double lrand = 0.0;		/* Amount of level randomness */
		int pnts;
		int fres;

		if (n == 0) {	/* Standard versions */
#ifdef NEVER		/* Doubled up points */
			pnts = 2 * PNTS;
			fres = GRES; 
			for (i = 0; i < pnts; i++) {
				xa[i * 2 + 0] = t1xa[i] - 0.01;
				ya[i * 2 + 0] = t1ya[i];
				xa[i * 2 + 1] = t1xa[i] + 0.01;
				ya[i * 2 + 1] = t1ya[i];
			}
#else
			pnts = PNTS;
			fres = GRES; 
			for (i = 0; i < pnts; i++) {
				xa[i] = t1xa[i];
				ya[i] = t1ya[i];
			}
#endif
			printf("Trial %d, points = %d, res = %d, level randomness = %f\n",n,pnts,fres,lrand);
		} else {	/* Random versions */
			lrand = d_rand(0.0,0.1);		/* Amount of level randomness */
			pnts = i_rand(MIN_PNTS,MAX_PNTS);
			fres = i_rand(MIN_RES,MAX_RES);

			printf("Trial %d, points = %d, res = %d, level randomness = %f\n",n,pnts,fres,lrand);

			/* Create X values */
			xa[0] = d_rand(0.5,1.0);
			for (i = 1; i < pnts; i++)
				xa[i] = xa[i-1] + d_rand(0.5,1.0);
			for (i = 0; i < pnts; i++)	/* Divide out */
				xa[i] = (xa[i]/xa[pnts-1]);

			/* Create y values */
			ya[0] = xa[0];
			for (i = 0; i < pnts; i++)
				ya[i] = ya[i-1] + d_rand(0.2,1.0) + d_rand(-0.2,0.3) + d_rand(-0.2,0.3);
			for (i = 0; i < pnts; i++)	/* Divide out */
				ya[i] = (ya[i]/ya[pnts-1]);
		}

		if (n < SKIP)
			continue;

		/* Create the object */
		rss =  new_rspl(RSPL_NOFLAGS,
		                1,				/* di */
		                1);				/* fdi */

		for (i = 0; i < pnts; i++) {
			test_points[i].p[0] = xa[i];
			test_points[i].v[0] = ya[i];
		}
		gres[0] = fres;

#ifdef RES2
		if (n != 0) {
#endif
		/* Fit to scattered data */
		rss->fit_rspl(rss,
		           0 | extra ? RSPL_EXTRAFIT : 0 ,
		           test_points,			/* Test points */
		           pnts,	/* Number of test points */
		           low, high, gres,		/* Low, high, resolution of grid */
		           NULL, NULL,			/* Default data scale */
		           smooth,				/* Smoothing */
		           avgdev,				/* Average deviation */
		           NULL);				/* iwidth */


		/* Display the result */
		for (i = 0; i < XRES; i++) {
			co tp;	/* Test point */
			x = i/(double)(XRES-1);
			xx[i] = x;
			yy[0][i] = lin(x,xa,ya,pnts);
			tp.p[0] = x;
			rss->interp(rss, &tp);
			yy[1][i] = tp.v[0];
			if (yy[1][i] < -0.2)
				yy[1][i] = -0.2;
			else if (yy[1][i] > 1.2)
				yy[1][i] = 1.2;
		}
		
		do_plot(xx,yy[0],yy[1],NULL,XRES);

#ifdef RES2
		} else {	/* Multiple resolution version */
			int gresses[5];
			for (j = 0; j < 5; j++) {
#ifndef NEVER
				if (j == 0)
					gres[0] = fres/8;
				else if (j == 1)
					gres[0] = fres/4;
				else if (j == 2)
					gres[0] = fres/2;
				else if (j == 3)
					gres[0] = fres;
				else 
					gres[0] = fres * 2;
#else 	/* Check sensitivity to griding of data points */
				if (j == 0)
					gres[0] = 192;
				else if (j == 1)
					gres[0] = 193;
				else if (j == 2)
					gres[0] = 194;
				else if (j == 3)
					gres[0] = 195;
				else 
					gres[0] = 196;
#endif
				gresses[j] = gres[0];
	
				rss->fit_rspl(rss,
#ifdef EXTRAFIT
		           RSPL_EXTRAFIT |		/* Extra fit flag */
#endif
			           0,
			           test_points,			/* Test points */
			           pnts,	/* Number of test points */
			           low, high, gres,		/* Low, high, resolution of grid */
			           NULL, NULL,			/* Default data scale */
			           smooth,				/* Smoothing */
			           avgdev,				/* Average deviation */
			           NULL);				/* iwidth */
	
				/* Get the result */
				for (i = 0; i < XRES; i++) {
					co tp;	/* Test point */
					x = i/(double)(XRES-1);
					xx[i] = x;
					yy[0][i] = lin(x,xa,ya,pnts);
					tp.p[0] = x;
					rss->interp(rss, &tp);
					yy[1+j][i] = tp.v[0];
					if (yy[1+j][i] < -0.2)
						yy[1+j][i] = -0.2;
					else if (yy[1+j][i] > 1.2)
						yy[1+j][i] = 1.2;
				}
			}
	
		printf("Black = lin, Red = %d, Green = %d, Blue = %d, Yellow = %d, Purple = %d\n",
		       gresses[0], gresses[1], gresses[2], gresses[3], gresses[4]);
		do_plot6(xx,yy[0],yy[1],yy[2],yy[3],yy[4],yy[5],XRES);
	}
#endif /* RES2 */
	}	/* next trial */
	return 0;
}


double
lin(
double x,
double xa[],
double ya[],
int n) {
	int i;
	double y;

	if (x < xa[0])
		return ya[0];
	else if (x > xa[n-1])
		return ya[n-1];

	for (i = 0; i < (n-1); i++)
		if (x >=xa[i] && x <= xa[i+1])
			break;

	x = (x - xa[i])/(xa[i+1] - xa[i]);

	y = ya[i] + (ya[i+1] - ya[i]) * x;
	
	return y;
}


/******************************************************************/
/* Error/debug output routines */
/******************************************************************/

/* Next u function done with optimization */

/* Structure to hold data for optimization function */
struct _edatas {
	rspl *rss;
	int j;
	}; typedef struct _edatas edatas;

#ifdef GLOB_CHECK
/* Overall Global optimization method */
/* Definition of the optimization function handed to powell() */
double efunc2(void *edata, double p[])
	{
	int j;
	double rv;
	rspl *rss = (rspl *)edata;
	for (j = 0; j < rss->nig; j++)	/* Ugg */
		rss->u[j].v = p[j];
	rv = rss->efactor(rss);
#ifdef DIAG2
	/* printf("\r%e",rv); */
	printf("%e\n",rv);
#endif
	return rv;
	}

solveu(rss)
rspl *rss;
	{
	int j;
	double *cp;
	double *s;

	cp = dvector(0,rss->nig);
	s = dvector(0,rss->nig);
	for (j = 0; j < rss->nig; j++)	/* Ugg */
		{
		cp[j] = rss->u[j].v;
		s[j] = 0.1;
		}
	powell(rss->nig,cp,s,0.0000001,1000,efunc2,(void *)rss);
	}
#endif /* GLOB_CHECK */
