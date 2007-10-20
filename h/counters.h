
/* 
 * Argyll Color Correction System
 * Multi-dimensional counter macros.
 *
 * Author: Graeme W. Gill
 * Date:   28/9/96
 *
 * Copyright 1996 - 2006, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#ifndef COUNTERS_H

/* ------------------------------------------------------- */
/* Macros for a multi-dimensional counter. */
/* Declare the counter name nn, maximum di mxdi, dimensions di, & count */

#define DCOUNT(nn, mxdi, di, start, reset, count) 				\
	int nn[mxdi];	/* counter value */						\
	int nn##_di = (di);		/* Number of dimensions */		\
	int nn##_stt = (start);	/* start count value */			\
	int nn##_rst = (reset);	/* reset on carry value */		\
	int nn##_res = (count); /* last count +1 */				\
	int nn##_e				/* dimension index */

/* Set the counter value to 0 */
#define DC_INIT(nn) 								\
{													\
	for (nn##_e = 0; nn##_e < nn##_di; nn##_e++)	\
		nn[nn##_e] = nn##_stt;						\
	nn##_e = 0;										\
}

/* Increment the counter value */
#define DC_INC(nn)									\
{													\
	for (nn##_e = 0; nn##_e < nn##_di; nn##_e++) {	\
		nn[nn##_e]++;								\
		if (nn[nn##_e] < nn##_res)					\
			break;	/* No carry */					\
		nn[nn##_e] = nn##_rst;						\
	}												\
}

/* After increment, expression is TRUE if counter is done */
#define DC_DONE(nn)									\
	(nn##_e >= nn##_di)
	
/* (Do we need a version of the above that tracks the actual input coords ?) */
/* ------------------------------------------------------- */
/* Same as above, but allows for variable resolution on each axis */

#define ECOUNT(nn, mxdi, di, count)				 				\
	int nn[mxdi];	/* counter value */						\
	int nn##_di = (di);		/* Number of dimensions */		\
	int *nn##_res = (count);/* last count +1 */				\
	int nn##_e				/* dimension index */

/* Set the counter value to 0 */
#define EC_INIT(nn) 								\
{													\
	for (nn##_e = 0; nn##_e < nn##_di; nn##_e++)	\
		nn[nn##_e] = 0;								\
	nn##_e = 0;										\
}

/* Increment the counter value */
#define EC_INC(nn)									\
{													\
	for (nn##_e = 0; nn##_e < nn##_di; nn##_e++) {	\
		nn[nn##_e]++;								\
		if (nn[nn##_e] < nn##_res[nn##_e])			\
			break;	/* No carry */					\
		nn[nn##_e] = 0;								\
	}												\
}

/* After increment, expression is TRUE if counter is done */
#define EC_DONE(nn)									\
	(nn##_e >= nn##_di)
	
/* (Do we need a version of the above that tracks the actual input coords ?) */

/* ------------------------------------------------------- */
/* Macros combination counter */
/* Declare the counter name nn, combinations out of total */
/* Maximum combinations is DI+2 */

#define COMBO(nn, mxdi, comb, total) 				\
	int nn[mxdi+2];			/* counter value */				\
	int nn##_cmb = (comb);	/* number of combinations*/		\
	int nn##_tot = (total);	/* out of total possible */		\
	int nn##_e				/* dimension index */

/* Set total to new setting */
#define CB_SETT(nn, total)					 		\
	nn##_tot = (total)	/* total possible */

/* Set combinations to new setting */
#define CB_SETC(nn, comb)					 		\
	nn##_cmb = (comb)	/* number of combinations*/

/* Set the counter to its initial value */
#define CB_INIT(nn) 								\
{													\
	for (nn##_e = 0; nn##_e < nn##_cmb; nn##_e++)	\
		nn[nn##_e] = nn##_cmb-nn##_e-1;				\
	nn##_e = 0;										\
}

/* Increment the counter value */
#define CB_INC(nn)									\
{													\
	for (nn##_e = 0; nn##_e < nn##_cmb; nn##_e++) {	\
		nn[nn##_e]++;								\
		if (nn[nn##_e] < (nn##_tot-nn##_e)) {		\
			int nn##_ee;		/* No carry */		\
			for (nn##_ee = nn##_e-1; nn##_ee >= 0; nn##_ee--)	\
				nn[nn##_ee] = nn[nn##_ee+1] + 1;	\
			break;									\
		}											\
	}												\
}

/* After increment, expression is TRUE if counter is done */
#define CB_DONE(nn)									\
	(nn##_e >= nn##_cmb)
	
/* ------------------------------------------------------- */
/* Macros simplex combination counter. */
/* Based on COMBO, but skips invalid simplex combinations */

#define XCOMBO(nn, mxdi, comb, total) 						\
		 COMBO(nn, mxdi, comb, total)

/* Set total to new setting */
#define XCB_SETT(nn, total)					 			\
         CB_SETT(nn, total)

/* Set combinations to new setting */
#define XCB_SETC(nn, comb)					 			\
         CB_SETC(nn, comb)


/* Set the counter to its initial value */
#define XCB_INIT(nn) 									\
{														\
	int nn##_ii;										\
														\
	for (nn##_e = 0; nn##_e < nn##_cmb; nn##_e++)		\
		nn[nn##_e] = nn##_cmb-nn##_e-1;					\
	for (nn##_ii = 1; nn##_ii < nn##_cmb; nn##_ii++) {	\
		if ((nn[nn##_ii-1] ^ nn[nn##_ii]) & nn[nn##_ii])\
			break;	/* Went from 0 to 1 */				\
	}													\
	if (nn##_ii < nn##_cmb)	{ /* Fix invalid combination */	\
		XCB_INC(nn);									\
	}													\
	nn##_e = 0;											\
}

/* Increment the counter value */
#define XCB_INC(nn)										\
{														\
	int nn##_ii = 0;									\
														\
	while (nn##_ii < nn##_cmb) {						\
		for (nn##_e = 0; nn##_e < nn##_cmb; nn##_e++) {	\
			nn[nn##_e]++;								\
			if (nn[nn##_e] < (nn##_tot-nn##_e)) {		\
				int nn##_ee;		/* No carry */		\
				for (nn##_ee = nn##_e-1; nn##_ee >= 0; nn##_ee--)	\
					nn[nn##_ee] = nn[nn##_ee+1] + 1;	\
				break;									\
			}											\
		}												\
		if (nn##_e >= nn##_cmb)							\
			break;		/* Done */						\
														\
		/* Reject invalid combinations */				\
		for (nn##_ii = 1; nn##_ii < nn##_cmb; nn##_ii++) {		\
			if ((nn[nn##_ii-1] ^ nn[nn##_ii]) & nn[nn##_ii]) 	\
				break;	/* Went from 0 to 1 */			\
		}												\
	}													\
}

/* After increment, expression is TRUE if counter is done */
#define XCB_DONE(nn)									\
         CB_DONE(nn)
	
/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

#define COUNTERS_H
#endif /* COUNTERS_H */
