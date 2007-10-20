/***************************************************/
/* Linear Simultaeous equation solver */
/***************************************************/

/* General simultaneous equation solver. */
/* Code was inspired by the algorithm decsribed in section 2.3 */
/* of "Numerical Recipes in C", by W.H.Press, B.P.Flannery, */
/* S.A.Teukolsky & W.T.Vetterling. */

/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include "numsup.h"
#include "ludecomp.h"

#undef DO_POLISH
#undef DO_CHECK

/* Solve the simultaneous linear equations A.X = B */
/* Return 1 if the matrix is singular, 0 if OK */
int
solve_se(
double **a,	/* A[][] input matrix, returns LU decomposition of A */
double  *b,	/* B[]   input array, returns solution X[] */
int      n	/* Dimensionality */
) {
	double rip;		/* Row interchange parity */
	int *pivx, PIVX[10];
#if defined(DO_POLISH) || defined(DO_CHECK)
	double **sa;	/* save input matrix values */
	double *sb;		/* save input vector values */
#endif

	if (n <= 10)
		pivx = PIVX;
	else
		pivx = ivector(0, n-1);

#if defined(DO_POLISH) || defined(DO_CHECK)
	sa = dmatrix(0, n-1, 0, n-1);
	sb = dvector(0, n-1);

	/* Copy input matrix and vector values */
	for (i = 0; i < n; i++) {
		sb[i] = b[i];
		for (j = 0; j < n; j++)
			sa[i][j] = a[i][j];
	}
#endif

	if (lu_decomp(a, n, pivx, &rip)) {
#if defined(DO_POLISH) || defined(DO_CHECK)
		free_dvector(sb, 0, n-1);
		free_dmatrix(sa, 0, n-1, 0, n-1);
		if (pivx != PIVX)
			free_ivector(pivx, 0, n-1);
#endif
		return 1;
	}

	lu_backsub(a, n, pivx, b);

#ifdef DO_POLISH
	lu_polish(n, sa, a, pivx, sb, b);	/* Improve the solution */
#endif

#ifdef DO_CHECK
	/* Check that the solution is correct */
	for (i = 0; i < n; i++) {
		double sum, temp;
		sum = 0.0;
		for (j = 0; j < n; j++)
			sum += sa[i][j] * b[j];
		temp = fabs(sum - sb[i]);
		if (temp > 1e-6) {
			free_dvector(sb, 0, n-1);
			free_dmatrix(sa, 0, n-1, 0, n-1);
			if (pivx != PIVX)
				free_ivector(pivx, 0, n-1);
			return 2;
		}
	}
#endif
#if defined(DO_POLISH) || defined(DO_CHECK)
	free_dvector(sb, 0, n-1);
	free_dmatrix(sa, 0, n-1, 0, n-1);
#endif
	if (pivx != PIVX)
		free_ivector(pivx, 0, n-1);
	return 0;
}

/* Decompose the square matrix A[][] into lower and upper triangles */
/* Return 1 if the matrix is singular. */
int
lu_decomp(
double **a,		/* A input array, output upper and lower triangles. */
int      n,		/* Dimensionality */
int     *pivx,	/* Return pivoting row permutations record */
double  *rip	/* Row interchange parity, +/- 1.0, used for determinant */
) {
	int    i, j;
	double *rscale, RSCALE[10];		/* Implicit scaling of each row */

	if (n <= 10)
		rscale = RSCALE;
	else
		rscale = dvector(0, n-1);

	/* For each row */
	for (i = 0; i < n; i++) {
		double big;
		/* For each column in row */
		for (big = 0.0, j=0; j < n; j++) {
			double temp;
			temp = fabs(a[i][j]);
			if (temp > big)
				big = temp;
		}
		if (fabs(big) <= DBL_MIN) {
			if (rscale != RSCALE)
				free_dvector(rscale, 0, n-1);
			return 1;		/* singular matrix */
		}
		rscale[i] = 1.0/big;	/* Save the scaling */
	}

	/* For each column (Crout's method) */
	for (*rip = 1.0, j = 0; j < n; j++) {
		double big;
		int k, bigi = 0;

		/* For each row */
		for (i = 0; i < j; i++) {
			double sum;
			sum = a[i][j];
			for (k = 0; k < i; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
		}

		/* Find largest pivot element */
		for (big = 0.0, i = j; i < n; i++) {
			double sum, temp;

			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;

			temp = rscale[i] * fabs(sum);	/* Figure of merit */
			if (temp >= big) {
				big = temp;		/* Best so far */
				bigi = i;		/* Remember index */
			}
		}
		
		/* If we need to interchange rows */
		if (j != bigi) {
			{	/* Take advantage of matrix storage to swap pointers to rows */
				double *temp;
				temp = a[bigi];
				a[bigi] = a[j];
				a[j] = temp;
			}
			*rip = -(*rip);				/* Another row interchange */
			rscale[bigi] = rscale[j];	/* Interchange scale factor */
		}
		
		pivx[j] = bigi;					/* Record pivot */
		if (fabs(a[j][j]) <= DBL_MIN) {
			if (rscale != RSCALE)
				free_dvector(rscale, 0, n-1);
			return 1; 					/* Pivot element is zero, so matrix is singular */
		}

		/* Divide by the pivot element */
		if (j != (n-1)) {
			double temp;
			temp = 1.0/a[j][j];
			for (i = j+1; i < n; i++)
				a[i][j] *= temp;
		}
	}
	if (rscale != RSCALE)
		free_dvector(rscale, 0, n-1);
	return 0;
}

/* Solve a set of simultaneous equations from the */
/* LU decomposition, by back substitution. */
void
lu_backsub(
double **a,		/* A[][] LU decomposed matrix */
int      n,		/* Dimensionality */
int     *pivx,	/* Pivoting row permutations record */
double  *b		/* Input B[] vecor, return X[] */
) {
	int i, j;
	int nvi;		/* When >= 0, indicates non-vanishing B[] index */

	/* Forward substitution, undo pivoting on the fly */
	for (nvi = -1, i = 0; i < n; i++) {
		int px;
		double sum;

		px = pivx[i];
		sum = b[px];
		b[px] = b[i];
		if (nvi >= 0) {
			for (j = nvi; j < i; j++)
				sum -= a[i][j] * b[j];
		} else {
			if (sum != 0.0)
				nvi = i;			/* Found non-vanishing element */
		}
		b[i] = sum;
	}

	/* Back substitution */
	for (i = (n-1); i >= 0; i--) {
		double sum;
		sum = b[i];
		for (j = i+1; j < n; j++)
			sum -= a[i][j] * b[j];
		b[i] = sum/a[i][i];
	}
}


/* Improve a solution of equations */
void
lu_polish(
double **a,			/* Original A[][] matrix */
double **lua,		/* LU decomposition of A[][] */
int      n,			/* Dimensionality */
double  *b,			/* B[] vector of equation */
double  *x,			/* X[] solution to be polished */
int     *pivx		/* Pivoting row permutations record */
) {
	int i, j;
	double *r, R[10];		/* Residuals */

	if (n <= 10)
		r = R;
	else
		r = dvector(0, n-1);

	/* Accumulate the residual error */
	for (i = 0; i < n; i++) {
		double sum;
		sum = -b[i];
		for (j = 0; j < n; j++)
			sum += a[i][j] * x[j];
		r[i] = sum;
	}

	/* Solve for the error */
	lu_backsub(lua, n, pivx, r);

	/* Subtract error from the old solution */
	for (i = 0; i < n; i++)
		x[i] -= r[i];
	
	if (r != R)
		free_dvector(r, 0, n-1);
}













