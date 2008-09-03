
/* Multi-dimentional minizer using Powell or Conjugate Gradient methods */
/* This is good for smoother, well behaved functions. */

/* Code is an original expression of the algorithms decsribed in */
/* "Numerical Recipes in C", by W.H.Press, B.P.Flannery, */
/* S.A.Teukolsky & W.T.Vetterling. */

/*
 * Copyright 2000, 2006 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* TTBD:
   Fix error handling to return status (malloc, excessive itters)
   Create to "safe" library ?
   Make standalone - ie remove numsup ?
 */

/* Note that all arrays are indexed from 0 */

#include "numsup.h"
#include "powell.h"

#undef SLOPE_SANITY_CHECK		/* exermental */
#undef DEBUG					/* Some debugging printfs (not comprehensive) */

#ifdef DEBUG
#undef DBG
#define DBG(xxx) printf xxx ;
#else
#undef DBG
#define DBG(xxx) 
#endif

static double linmin(double p[], double xi[], int n, double ftol,
	double (*func)(void *fdata, double tp[]), void *fdata);

/* Standard interface for powell function */
/* return 0 on sucess, 1 on failure due to excessive itterions */
/* Result will be in cp */
int powell(
double *rv,				/* If not NULL, return the residual error */
int di,					/* Dimentionality */
double cp[],			/* Initial starting point */
double s[],				/* Size of initial search area */
double ftol,			/* Tollerance of error change to stop on */
int maxit,				/* Maximum iterations allowed */
double (*func)(void *fdata, double tp[]),		/* Error function to evaluate */
void *fdata				/* Opaque data needed by function */
) {
	int i;
	double **dmtx;		/* Direction vector */
	double *spt;		/* Sarting point before exploring all the directions */
	double *xpt;		/* Extrapolated point */
	double *svec;		/* Search vector */
	int    iter;
	double retv; 		/* Returned function value at p */

	dmtx = dmatrixz(0, di-1, 0, di-1);	/* Zero filled */
	spt  = dvector(0,di-1);
	xpt  = dvector(0,di-1);
	svec = dvector(0,di-1);

	/* Create initial direction matrix by */
	/* placing search start on diagonal */
	for (i = 0; i < di; i++)
		dmtx[i][i] = s[i];
	
	/* Save the starting point */
	for (i = 0; i < di; i++)
		spt[i] = cp[i];

	/* Initial function evaluation */
	retv = (*func)(fdata, cp);

	/* Itterate untill we converge on a solution, or give up. */
	for (iter = 1; iter < maxit; iter++) {
		int j;
		double lretv;			/* Last function return value */
		int    ibig = 0;		/* Index of biggest delta */
		double del = 0.0;		/* Biggest function value decrease */
		double pretv;			/* Previous function return value */

		pretv = retv;			/* Save return value at top of itteration */

		/* Loop over all directions in the set */
		for (i = 0; i < di; i++) {

			for (j = 0; j < di; j++)	/* Extract this direction to make search vector */
				svec[j] = dmtx[j][i];

			/* Minimize in that direction */
			lretv = retv;
			retv = linmin(cp, svec, di, ftol, func, fdata);

			/* Record bigest function decrease, and dimension it occured on */
			if (fabs(lretv - retv) > del) {
				del = fabs(lretv - retv);
				ibig = i;
			}
		}

		/* If we have reached a suitable tollerance, then finish */
		if (2.0 * fabs(pretv - retv) <= ftol * (fabs(pretv) + fabs(retv))) {
			break;
		}

		for (i = 0; i < di; i++) {
			svec[i] = cp[i] - spt[i];	/* Average direction moved after minimization round */
			xpt[i]  = cp[i] + svec[i];	/* Extrapolated point after round of minimization */
			spt[i]  = cp[i];			/* New start point for next round */
		}

		/* Function value at extrapolated point */
		lretv = (*func)(fdata, xpt);

		if (lretv < pretv) {			/* If extrapolation is an improvement */
			double t, t1, t2;

			t1 = pretv - retv - del;
			t2 = pretv - lretv;
			t = 2.0 * (pretv -2.0 * retv + lretv) * t1 * t1 - del * t2 * t2;
			if (t < 0.0) {
				/* Move to the minimum of the new direction */
				retv = linmin(cp, svec, di, ftol, func, fdata);

				for (i = 0; i < di; i++) 		/* Save the new direction */
					dmtx[i][ibig] = svec[i];
			}
		}
	}

//printf("~1 iters = %d\n",iter);
	/* Free up all the temporary vectors and matrix */
	free_dvector(svec,0,di-1);
	free_dvector(xpt,0,di-1);
	free_dvector(spt,0,di-1);
	free_dmatrix(dmtx, 0, di-1, 0, di-1);

	if (rv != NULL)
		*rv = retv;

	if (iter < maxit)
		return 0;

	return 1;		/* Failed due to execessive itterations */
}

/* -------------------------------------- */
/* Conjugate Gradient optimiser */
/* return 0 on sucess, 1 on failure due to excessive itterions */
/* Result will be in cp */
/* Note that we could use gradient in line minimiser, */
/* but haven't bothered yet. */
int conjgrad(
double *rv,				/* If not NULL, return the residual error */
int di,					/* Dimentionality */
double cp[],			/* Initial starting point */
double s[],				/* Size of initial search area */
double ftol,			/* Tollerance of error change to stop on */
int maxit,				/* Maximum iterations allowed */
double (*func)(void *fdata, double tp[]),		/* Error function to evaluate */
double (*dfunc)(void *fdata, double dp[], double tp[]),		/* Gradient function to evaluate */
void *fdata				/* Opaque data needed by function */
) {
	int i, iter;
	double *svec;		/* Search vector */
	double *ssvec;		/* Scaled search vector */
	double *gvec;		/* G direction vector */
	double *hvec;		/* H direction vector */
	double retv; 		/* Returned function value at p */

	svec = dvector(0,di-1);
	ssvec = dvector(0,di-1);
	gvec  = dvector(0,di-1);
	hvec  = dvector(0,di-1);

	/* Initial function evaluation */
	retv = (*dfunc)(fdata, svec, cp);

	/* Initial vector setup */
	for (i = 0; i < di; i++) {
		gvec[i] = hvec[i] = svec[i] = -svec[i];	/* Inverse gradient */
		ssvec[i] = s[i] * svec[i];
	}

	/* Itterate untill we converge on a solution, or give up. */
	for (iter = 1; iter < maxit; iter++) {
		double gamden, gamnum, gam;
		double pretv;			/* Previous function return value */

		DBG(("conjrad: about to do linmin\n"))
		pretv = retv;
		retv = linmin(cp, ssvec, di, 5.0 * ftol, func, fdata);

		/* If we have reached a suitable tollerance, then finish */
		/* (I don't understand why this tends to stop sooner than powell) */
		if (20.0 * fabs(pretv - retv) <= ftol * (fabs(pretv) + fabs(retv))) {
			break;
		}

		DBG(("conjrad: recomputing direction\n"))
		(*dfunc)(fdata, svec, cp);

		/* Compute gamma */
		for (gamnum = gamden = 0.0, i = 0; i < di; i++) {
			gamnum += svec[i] * (gvec[i] + svec[i]);
			gamden += gvec[i] * gvec[i];
		}

		if (gamden == 0.0) {		/* Gradient is exactly zero */
			DBG(("conjrad: gradient is exactly zero\n"))
			break;
		}

		gam = gamnum/gamden;
		DBG(("conjrad: gamma = %f = %f/%f\n",gam,gamnum,gamden))

		/* Adjust seach direction */
		for (i = 0; i < di; i++) {
			gvec[i] = -svec[i];
			svec[i] = hvec[i] = gvec[i] + gam * hvec[i];
			ssvec[i] = s[i] * svec[i];
		}

	}
	/* Free up all the temporary vectors and matrix */
	free_dvector(hvec,0,di-1);
	free_dvector(gvec,0,di-1);
	free_dvector(svec,0,di-1);
	free_dvector(ssvec,0,di-1);

	if (rv != NULL)
		*rv = retv;

	if (iter < maxit)
		return 0;

	return 1;		/* Failed due to execessive itterations */
}

/*------------------------------*/
#define POWELL_GOLD 1.618034
#define POWELL_CGOLD 0.3819660
#define POWELL_MAXIT 100

/* Line bracketing and minimisation routine. */
/* Return value at minimum. */
static double linmin(
double cp[],		/* Start point, and returned value */
double xi[],		/* Search vector */
int di,				/* Dimensionality */
double ftol,		/* Tolerance to stop on */
double (*func)(void *fdata, double tp[]),		/* Error function to evaluate */
void *fdata)		/* Opaque data for func() */
{
	int i;
	double ax, xx, bx;	/* Search vector multipliers */
	double af, xf, bf;	/* Function values at those points */
	double *xt, XT[10];	/* Trial point */

	if (di <= 10)
		xt = XT;
	else
		xt = dvector(0, di-1);			/* Vector for trial point */

	/* -------------------------- */
	/* First bracket the solution */

	DBG(("linmin: Bracketing solution\n"))

	/* The line is measured as startpoint + offset * search vector */
	/* Start with ax being vector offset 0.0 */
	ax = 0.0;
	for (i = 0; i < di; i++)
		xt[i] = cp[i] + ax * xi[i];
	af = (*func)(fdata, xt);

	/* xx being vector offset 0.618 */
	xx =  1.0/POWELL_GOLD;
	for (i = 0; i < di; i++)
		xt[i] = cp[i] + xx * xi[i];
	xf = (*func)(fdata, xt);

	DBG(("linmin: Initial points a:%f:%f -> b:%f:%f\n",ax,af,xx,xf))

	/* Fix it so that we are decreasing from point a -> x */
	if (xf > af) {
		double tt;
		tt = ax; ax = xx; xx = tt;
		tt = af; af = xf; xf = tt;
	}
	DBG(("linmin: Ordered Initial points a:%f:%f -> b:%f:%f\n",ax,af,xx,xf))

	bx = xx + POWELL_GOLD * (xx-ax);	/* Guess b beyond a -> x */
	for (i = 0; i < di; i++)
		xt[i] = cp[i] + bx * xi[i];
	bf = (*func)(fdata, xt);

	DBG(("linmin: Initial bracket a:%f:%f x:%f:%f b:%f:%f\n",ax,af,xx,xf,bx,bf))

#ifdef SLOPE_SANITY_CHECK
	/* If we're not seeing a slope indicitive of progress */
	/* of order ftol, give up straight away */
	if (2000.0 * fabs(xf - bf) <= ftol * (fabs(xf) + fabs(bf))
	 && 2000.0 * fabs(af - xf) <= ftol * (fabs(af) + fabs(xf))) {
		DBG(("linmin: giving up because slope is too shallow\n"))
		if (xt != XT)
			free_dvector(xt,0,di-1);

		if (bf < xf) {
			xf = bf;
			xx = bx;
		}

		/* Compute solution vector */
		for (i = 0; i < di; i++) 
			cp[i] += xx * xi[i];
		return xf;
	}
#endif /* SLOPE_SANITY_CHECK */


	/* While not bracketed */
	while (xf > bf) {
		double ulim, ux, uf;
		double tt, r, q;

		DBG(("linmin: Not bracketer a:%f:%f x:%f%f b:%f:%f\n",ax,af,xx,xf,bx,bf))
//		DBG(("linmin: Not bracketed because xf %f > bf %f\n",xf, bf))
//		DBG(("        ax = %f, xx = %f, bx = %f\n",ax,xx,bx))

		/* Compute ux by parabolic interpolation from a, x & b */
		q = (xx - bx) * (xf - af);
		r = (xx - ax) * (xf - bf);
		tt = q - r;
		if (tt >= 0.0 && tt < 1e-20)				/* If +ve too small */
			tt = 1e-20;
		else if (tt <= 0.0 && tt > -1e-20)		/* If -ve too small */
			tt = -1e-20;
		ux = xx - ((xx - bx) * q - (xx - ax) * r) / (2.0 * tt);
		ulim = xx + 100.0 * (bx - xx);			/* Extrapolation limit */

//printf("~1 ux = %f, ulim = %f\n",ux,ulim);
		if ((xx - ux) * (ux - bx) > 0.0) {		/* u is between x and b */

			for (i = 0; i < di; i++)			/* Evaluate u */
				xt[i] = cp[i] + ux * xi[i];
			uf = (*func)(fdata, xt);

//printf("~1 u is between x and b, uf = %f\n",uf);

			if (uf < bf) {						/* Minimum is between x and b */
//printf("~1 min is between x and b\n");
				ax = xx; af = xf;
				xx = ux; xf = uf;
				break;
			} else if (uf > xf) {				/* Minimum is between a and u */
//printf("~1 min is between a and u\n");
				bx = ux; bf = uf;
				break;
			}

			/* Parabolic fit didn't work, look further out in direction of b */
			ux = bx + POWELL_GOLD * (bx-xx);
//printf("~1 parabolic fit didn't work,look further in direction of b (%f)\n",ux);

		} else if ((bx - ux) * (ux - ulim) > 0.0) {	/* u is between b and limit */
			for (i = 0; i < di; i++)			/* Evaluate u */
				xt[i] = cp[i] + ux * xi[i];
			uf = (*func)(fdata, xt);

//printf("~1 u is between b and limit uf = %f\n",uf);
			if (uf > bf) {						/* Minimum is between x and u */
//printf("~1 min is between x and uf\n");
				ax = xx; af = xf;
				xx = bx; xf = bf;
				bx = ux; bf = uf;
				break;
			}
			xx = bx; xf = bf;					/* Continue looking */
			bx = ux; bf = uf;
			ux = bx + POWELL_GOLD * (bx - xx);	/* Test beyond b */
//printf("~1 continue looking beyond b (%f)\n",ux);

		} else if ((ux - ulim) * (ulim - bx) >= 0.0) {	/* u is beyond limit */
			ux = ulim;
//printf("~1 use limit\n");
		} else {							/* u is to left side of x ? */
			ux = bx + POWELL_GOLD * (bx-xx);
//printf("~1 look gold beyond b (%f)\n",ux);
		}
		/* Evaluate u, and move into place at b */
		for (i = 0; i < di; i++)
			xt[i] = cp[i] + ux * xi[i];
		uf = (*func)(fdata, xt);
//printf("~1 lookup ux %f value uf = %f\n",ux,uf);
		ax = xx; af = xf;
		xx = bx; xf = bf;
		bx = ux; bf = uf;
//printf("~1 move along to the right (a<-x, x<-b, b-<u)\n");
	}
	DBG(("linmin: Got bracket a:%f:%f x:%f:%f b:%f:%f\n",ax,af,xx,xf,bx,bf))
	/* Got bracketed minimum between a -> x -> b */
//printf("~1 got bracketed minimum at %f (%f), %f (%f), %f (%f)\n",ax,af,xx,xf,bx,bf);

	/* --------------------------------------- */
	/* Now use brent minimiser bewteen a and b */
	{
		/* a and b bracket solution */
		/* x is best function value so far */
		/* w is second best function value so far */
		/* v is previous second best, or third best */
		/* u is most recently tested point */
		double wx, vx, ux;			/* Search vector multipliers */
		double wf, vf = 0.0, uf;	/* Function values at those points */
		int iter;
		double de = 0.0;	/* Distance moved on previous step */
		double e = 0.0;		/* Distance moved on 2nd previous step */

		/* Make sure a and b are in ascending order */
		if (ax > bx) {
			double tt;
			tt = ax; ax = bx; bx = tt;
			tt = af; af = bf; bf = tt;
		}

		wx = vx = xx;	/* Initial values of other center points */
		wf = xf = xf;

		for (iter = 1; iter <= POWELL_MAXIT; iter++) {
			double mx = 0.5 * (ax + bx);		/* m is center of bracket values */
			double tol1 = ftol * fabs(xx) + 1e-10;
			double tol2 = 2.0 * tol1;

			DBG(("linmin: Got bracket a:%f:%f x:%f:%f b:%f:%f\n",ax,af,xx,xf,bx,bf))

			/* See if we're done */
//printf("~1 linmin check %f <= %f\n",fabs(xx - mx), tol2 - 0.5 * (bx - ax));
			if (fabs(xx - mx) <= (tol2 - 0.5 * (bx - ax))) {
				DBG(("linmin: We're done because %f <= %f\n",fabs(xx - mx), tol2 - 0.5 * (bx - ax)))
				break;
			}

			if (fabs(e) > tol1) {			/* Do a trial parabolic fit */
				double te, p, q, r;
				r = (xx-wx) * (xf-vf);
				q = (xx-vx) * (xf-wf);
				p = (xx-vx) * q - (xx-wx) * r;
				q = 2.0 * (q - r);
				if (q > 0.0)
					p = -p;
				else
					q = -q;
				te = e;				/* Save previous e value */
				e = de;				/* Previous steps distance moved */

				DBG(("linmin: Trial parabolic fit\n" ))

				if (fabs(p) >= fabs(0.5 * q * te) || p <= q * (ax-xx) || p >= q * (bx-xx)) {
					/* Give up on the parabolic fit, and use the golden section search */
					e = ((xx >= mx) ? ax-xx : bx-xx);	/* Override previous distance moved */
					de = POWELL_CGOLD * e;
					DBG(("linmin: Moving to golden section search\n" ))
				} else {	/* Use parabolic fit */
					de = p/q;			/* Change in xb */
					ux = xx + de;		/* Trial point according to parabolic fit */
					if ((ux - ax) < tol2 || (bx - ux) < tol2) {
						if ((mx - xx) > 0.0)	/* Don't use parabolic, use tol1 */
							de = tol1;			/* tol1 is +ve */
						else
							de = -tol1;
					}
					DBG(("linmin: Using parabolic fit\n" ))
				}
			} else {	/* Keep using the golden section search */
				e = ((xx >= mx) ? ax-xx : bx-xx);	/* Override previous distance moved */
				de = POWELL_CGOLD * e;
				DBG(("linmin: Continuing golden section search\n" ))
			}

			if (fabs(de) >= tol1) {		/* If de moves as much as tol1 would */
				ux = xx + de;			/* use it */
				DBG(("linmin: ux = %f = xx %f + de %f\n",ux,xx,de))
			} else {					/* else move by tol1 in direction de */
				if (de > 0.0) {
					ux = xx + tol1;
					DBG(("linmin: ux = %f = xx %f + tol1 %f\n",ux,xx,tol1))
				} else {
					ux = xx - tol1;
					DBG(("linmin: ux = %f = xx %f - tol1 %f\n",ux,xx,tol1))
				}
			}

			/* Evaluate function */
			for (i = 0; i < di; i++)
				xt[i] = cp[i] + ux * xi[i];
			uf = (*func)(fdata, xt);

			if (uf <= xf) {					/* Found new best solution */
				if (ux >= xx) {	
					ax = xx; af = xf;		/* New lower bracket */
				} else {
					bx = xx; bf = xf;		/* New upper bracket */
				}
				vx = wx; vf = wf;			/* New previous 2nd best solution */
				wx = xx; wf = xf;			/* New 2nd best solution from previous best */
				xx = ux; xf = uf;			/* New best solution from latest */
				DBG(("linmin: found new best solution\n"))
			} else {						/* Found a worse solution */
				if (ux < xx) {
					ax = ux; af = uf;		/* New lower bracket */
				} else {
					bx = ux; bf = uf;		/* New upper bracket */
				}
				if (uf <= wf || wx == xx) {	/* New 2nd best solution, or equal best */
					vx = wx; vf = wf;		/* New previous 2nd best solution */
					wx = ux; wf = uf;		/* New 2nd best from latest */
				} else if (uf <= vf || vx == xx || vx == wx) {	/* New 3rd best, or equal 1st & 2nd */
					vx = ux; vf = uf;		/* New previous 2nd best from latest */
				}
				DBG(("linmin: found new worse solution\n"))
			}
		}
		/* !!! should do something if iter > POWELL_MAXIT !!!! */
		/* Solution is at xx, xf */

		/* Compute solution vector */
		for (i = 0; i < di; i++) 
			cp[i] += xx * xi[i];
	}

	if (xt != XT)
		free_dvector(xt,0,di-1);
	// printf("~~~ line minimizer returning %e\n",xf);
	return xf;
}

#undef POWELL_GOLD
#undef POWELL_CGOLD
#undef POWELL_MAXIT

/**************************************************/
