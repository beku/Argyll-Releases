
#ifndef LUDECOMP_H
#define LUDECOMP_H

/*
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#ifdef __cplusplus
	extern "C" {
#endif

/* NOTE:- lu decomp rearanges the rows of the matrix */
/* by swapping pointers rather than exchanging data, */
/* so the matrix must be addressed by the **pointer */
/* if it is re-used after an ludecomp!!! */

/* Solve the simultaneous linear equations A.X = B */
/* Return 1 if the matrix is singular, 0 if OK */
int
solve_se(
double **a,	/* A[][] input matrix, returns LU decimposition of A */
double  *b,	/* B[]   input array, returns solution X[] */
int      n	/* Dimensionality */
);

/* Decompose the square matrix A[][] into lower and upper triangles */
/* Return 1 if the matrix is singular. */
int
lu_decomp(
double **a,		/* A input array, output upper and lower triangles. */
int      n,		/* Dimensionality */
int     *pivx,	/* Return pivoting row permutations record */
double  *rip	/* Row interchange parity, +/- 1.0, used for determinant */
);

/* Solve a set of simultaneous equations from the */
/* LU decomposition, by back substitution. */
void
lu_backsub(
double **a,		/* A[][] LU decomposed matrix */
int      n,		/* Dimensionality */
int     *pivx,	/* Pivoting row permutations record */
double  *b		/* Input B[] vecor, return X[] */
);

/* Polish a solution for equations */
void
lu_polish(
double **a,			/* Original A[][] matrix */
double **lua,		/* LU decomposition of A[][] */
int      n,			/* Dimensionality */
double  *b,			/* B[] vector of equation */
double  *x,			/* X[] solution to be polished */
int     *pivx		/* Pivoting row permutations record */
);

#ifdef __cplusplus
	}
#endif

#endif /* LUDECOMP_H */
