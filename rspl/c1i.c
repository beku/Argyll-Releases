
/************************************************/
/* Investigate various curve approximations     */
/************************************************/

/* Discrete regularized spline versions */
/* Standard test + Random testing */
/* Test incremental version */

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
#undef RES2			/* Do second resolution test */

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
#define MAX_RES 500

#define XRES 100
#define GRES 100

#define PNTS 10
#define PNTS2 8
double t1xa[PNTS] = { 0.2, 0.25, 0.30, 0.35,  0.40,  0.44, 0.48, 0.75,  0.51,  0.64  };
double t1ya[PNTS] = { 0.3, 0.35, 0.4,  0.41,  0.42,  0.46, 0.5,  0.75,  0.575, 0.48  };
co tpts[MAX_PNTS];

double lin(double x, double xa[], double ya[], int n);

int main() {
	int i, j;
	double x;
	double xx[XRES];		/* Up to 3 graphs */
	double y1[XRES];		/* Full */
	double y2[XRES];		/* 8 points */
	double y3[XRES];		/* 8 + 2 points */
	rspl *rss;		/* incremental solution version */
	datai low,high;
	int gres[MXDI];
	double avgdev[MXDO];

	low[0] = 0.0;
	high[0] = 1.0;
	avgdev[0] = 0.0;

	error_program = "Curve1";

	{
		int pnts;

		/* Standard versions with all the points */
		pnts = PNTS;
		for (i = 0; i < pnts; i++) {
			tpts[i].p[0] = t1xa[i];
			tpts[i].v[0] = t1ya[i];
		}
		gres[0] = GRES;

		/* Create the object */
		rss =  new_rspl(RSPL_NOFLAGS, 1, 1);		/* di */

		/* Fit to scattered data */
		rss->fit_rspl(rss,
		           /* RSPL_EXTRAFIT | */	/* Extra fit flag */
		           0,
		           tpts,			/* Test points */
		           pnts,	/* Number of test points */
		           low, high, gres,	/* Low, high, resolution of grid */
		           NULL, NULL,			/* Default data scale */
		           1.0,					/* Smoothing */
		           avgdev,				/* Average deviation */
		           NULL);				/* iwidth */

		/* Save the result */
		for (i = 0; i < XRES; i++) {
			co tp;	/* Test point */
			x = i/(double)(XRES-1);
			xx[i] = x;
			tp.p[0] = x;
			rss->interp(rss, &tp);
			y1[i] = tp.v[0];
			if (y1[i] < -0.2)
				y1[i] = -0.2;
			else if (y1[i] > 1.2)
				y1[i] = 1.2;
		}

		/* Done with that one */
		rss->del(rss);

		/* Create the object */
		rss = new_rspl(RSPL_NOFLAGS, 1, 1);				/* di */

		/* Do second run with 8 points */
		pnts = PNTS2;
		for (i = 0; i < pnts; i++) {
			tpts[i].p[0] = t1xa[i];
			tpts[i].v[0] = t1ya[i];
		}
		gres[0] = GRES;

		/* Fit to scattered data */
		rss->fit_rspl(rss,
		           RSPL_INCREMENTAL |
		           /* RSPL_EXTRAFIT | */	/* Extra fit flag */
		           0,
		           tpts,			/* Test points */
		           pnts,	/* Number of test points */
		           low, high, gres,/* Low, high, resolution of grid */
		           NULL, NULL,			/* Default data scale */
		           1.0,					/* Smoothing */
		           avgdev,				/* Average deviation */
		           NULL);				/* iwidth */

		/* Save the result */
		for (i = 0; i < XRES; i++) {
			co tp;	/* Test point */
			x = i/(double)(XRES-1);
			xx[i] = x;
			tp.p[0] = x;
			rss->interp(rss, &tp);
			y2[i] = tp.v[0];
			if (y2[i] < -0.2)
				y2[i] = -0.2;
			else if (y2[i] > 1.2)
				y2[i] = 1.2;
		}

#ifndef NEVER
		/* Do second run adding 2 points */
		pnts = PNTS2;
		for (j = 0, i = PNTS2; i < PNTS; i++, j++) {
			tpts[j].p[0] = t1xa[i];
			tpts[j].v[0] = t1ya[i];
		}

#ifdef NEVER
/* Extra points same as first two */
tpts[0].p[0] = 0.2;
tpts[0].v[0] = 0.3;
tpts[1].p[0] = 0.30;
tpts[1].v[0] = 0.35;
#endif

#ifdef NEVER
/* Very slightly different last 2 points to distinguish the two graphs */
tpts[0].p[0] = 0.51;
tpts[0].v[0] = 0.575;
tpts[1].p[0] = 0.64;
tpts[1].v[0] = 0.40;
#endif

		/* Add the last 2 points */
		rss->add_rspl(rss,
		           0,
		           tpts,			/* Test points */
		           2);					/* Number of test points */

		/* Save the result */
		for (i = 0; i < XRES; i++) {
			co tp;	/* Test point */
			x = i/(double)(XRES-1);
			xx[i] = x;
			tp.p[0] = x;
			rss->interp(rss, &tp);
			y3[i] = tp.v[0];
			if (y3[i] < -0.2)
				y3[i] = -0.2;
			else if (y3[i] > 1.2)
				y3[i] = 1.2;
		}

		rss->del(rss);

		printf("Black = Std all, Red = minus2, Green = incremental\n");
		do_plot(xx,y1,y2,y3,XRES);
#else
		printf("Black = Std all, Red = minus2\n");
		do_plot(xx,y1,y2,NULL,XRES);
#endif
	}
		
	return 0;
}


double
lin(
double x,
double xa[],
double ya[],
int n)
	{
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

void
usage(void) {
	printf("1D rspl curve incremental  test by Graeme Gill, Version %s\n",ARGYLL_VERSION_STR);
	printf("usage: c1i\n");
	exit(1);
	}

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
