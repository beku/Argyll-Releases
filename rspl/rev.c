/* 
 * Argyll Color Correction System
 * Multi-dimensional regularized spline data structure
 *
 * Reverse interpolation support code.
 *
 * Author: Graeme W. Gill
 * Date:   30/1/00
 *
 * Copyright 1999 - 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 *
 * Latest simplex/linear equation version.
 */

/* TTBD:
	Get rid of error() calls - return status instead

	Need to add a hefty overview and explanation of
	how all this works, before I forget it !

	ie:

	  Basic function requirements:  exact, auxil, locus, clip

	  Fwd cell - reverse cell list lookup

	  Basic layout di -> fdi + auxils + ink limit

	  Basic search strategy

	  Sub Simplex decomposition & properties

	  How each type of function finds solutions
		Sub-simplex dimensionality & dof + target dim & dof
		Linear algebra choices.
		
	  How final solutions are chosen

 */

/* TTBD:

	Allow function callback to set auxiliary values for 
	flag RSPL_AUXLOCUS. 
	How to pass enough info back to aux_compute() ?

	Should auxil return multiple solutions if it finds them ???

 */


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <memory.h>
#include <time.h>

#undef SSS				/* ~~~ TEST exaustinve nn search */

#define INKSCALE 5000.0	/* For ink limit weighting to fudge SVD least squares solution */

#include "rspl_imp.h"
#include "numlib.h"
#include "../h/sort.h"		/* Heap sort */

//#define DMALLOC_GLOBALS
//#include "dmalloc.h"
//#undef DMALLOC_GLOBALS

#undef DEBUG

#define DOSORT				/* Cell sort */
#define ALLOC_ALL_CACHE		/* Allocate the maximum cell cache at the start */

/* Print a vectors value */
#define DBGVI(text, dim, out, vec, end)			\
{	int pveci;									\
	printf("%s",text);							\
	for (pveci = 0 ; pveci < (dim); pveci++)		\
		printf(out,(vec)[pveci]);				\
	printf(end);								\
}

/* Print a matrix value */
#define DBGMI(text, rows, cols, out, mat, end)		\
{	int pveci, pvecr;								\
	printf("%s",text);								\
	for (pvecr = 0 ; pvecr < (rows); pvecr++) {		\
		for (pveci = 0 ; pveci < (cols); pveci++)		\
			printf(out,(mat)[pvecr][pveci]);		\
	if ((pvecr+1) < (rows))							\
		printf("\n");								\
	}												\
	printf(end);									\
}

/* Do an arbitrary printf */
#define DBGI(text) printf text ;

#undef DBG
#undef DBGV
#undef DBGM

#undef NEVER
#define ALWAYS

#undef DEBUG	// ~~1

#ifdef DEBUG
#undef DBG
#undef DBGV
#undef DBGM
#define DBG(xxx) DBGI(xxx)
#define DBGV(xxx) DBGVI xxx
#define DBGM(xxx) DBGMI xxx
#else
#undef DBG
#undef DBGV
#undef DBGM
#define DBG(xxx) 
#define DBGV(xxx) 
#define DBGM(xxx) 
#endif

/* Convention is to use:
   i to index grid points u.a
   n to index data points d.a
   e to index position dimension di
   f to index output function dimension fdi
   j misc and cube corners
   k misc
 */

#define	EPS (1e-10)			/* Allowance for numeric error */

extern void error(char *fmt, ...), warning(char *fmt, ...);

static void make_rev(rspl *s);

static cell *get_rcell(schbase *b, int ix);
static void uncache_rcell(revcache *r, cell *cp);
#define unget_rcell(r, cp) uncache_rcell(r, cp)		/* These are the same */
static void invalidate_revcache(revcache *rc);

/* ====================================================== */

static schbase *init_search(rspl *s, int flags, double *av, int *auxm,
                        double *v, double *cdir, co *cpp, int mxsoln, enum ops op);
static void adjust_search(rspl *s, int flags, double *av, enum ops op);
static schbase *set_search_limit(rspl *s, double (*limit)(void *vcntx, double *in),
                                 void *lcntx, double limitv);
static void set_lsearch(rspl *s, int e);
static void free_search(schbase *b);

static int *calc_fwd_cell_list(rspl *s, double *v);

static void init_line_eq(schbase *b, double st[MXRO], double de[MXRO]);
static int *init_line(rspl *s, line *l, double st[MXRO], double de[MXRO]);
static int *next_line_cell(line *l);

#ifdef SSS
/* Structure to hold clip sphere state information */
typedef struct {
	struct _rspl *s;	/* Pointer to parent rspl */
	int ix[MXRO];		/* Offset counter */
	int c;				/* layer */
	int base[MXRO];		/* base grid location */
	double off[MXRO];	/* offset within base cell */
	int loc[MXRO];		/* Current location */
	double cdist;		/* Best possible cell distance */
	double ldist;		/* Best possible layer distance */
} sphere;

static int *init_sphere(rspl *s, sphere *l, double st[MXRO]);
static int *next_sphere_cell(sphere *l);
#endif /* SSS */

static int *get_nearest_list(rspl *s, double *v);

static void search_list(schbase *b, int *rip, unsigned int tcount);

static void clear_limitv(rspl *s);

#ifdef STATS
static char *opnames[6] = { "exact", "clipv", "clipn", "auxil", "locus" };
#endif /* STATS */

#define INF_DIST 1e38		/* Stands for infinite "current best" distance */

/* ====================================================== */
/* Set the ink limit information for any reverse interpolation. */
/* Calling this will clear the reverse interpolaton cache. */
static void
rev_set_limit_rspl(
	rspl *s,		/* this */
	double (*limit)(void *lcntx, double *in),	/* Optional input space limit function. Function */
					/* should evaluate in[0..di-1], and return number that is not to exceed */
					/* limitv. NULL if not used */
	void *lcntx,	/* Context passed to limit() */
	double limitv	/* Value that limit() is not to exceed */
) {
	schbase *b;

	/* This is a restricted size function */
	if (s->di > MXRI)
		error("rspl: rev_set_limit can't handle di = %d",s->di);
	if (s->fdi > MXRO)
		error("rspl: rev_set_limit can't handle fdi = %d",s->fdi);

	b = set_search_limit(s, limit, lcntx, limitv);	/* Init and set limit info */

	if (s->rev.inited) {		/* If cache has been allocated */
		invalidate_revcache(s->rev.cache);		/* Invalidate the reverse cache */
	}

	/* Invalidate any ink limit values cached with the grid data */
	clear_limitv(s);
}

/* Get the ink limit information for any reverse interpolation. */
static void
rev_get_limit_rspl(
	rspl *s,		/* this */
	double (**limit)(void *lcntx, double *in),	/* Return pointer to function of NULL if not set */
	void **lcntx,	/* return context pointer */
	double *limitv	/* Return limit value */
) {
	schbase *b = s->rev.sb;

	/* This is a restricted size function */
	if (s->di > MXRI)
		error("rspl: rev_get_limit can't handle di = %d",s->di);
	if (s->fdi > MXRO)
		error("rspl: rev_get_limit can't handle fdi = %d",s->fdi);

	if (b == NULL) {
		*limit = NULL;
		*lcntx = NULL;
		*limitv = 0.0;
	} else {
		*limit = b->limit;
		*lcntx = b->cntx;
		*limitv = b->limitv/ INKSCALE;
	}
}

#define RSPL_CERTAIN 0x80000000 						/* WILLCLIP hint is certain */
#define RSPL_WILLCLIP2 (RSPL_CERTAIN | RSPL_WILLCLIP)	/* Clipping will certainly be needed */

/* Do reverse interpolation given target output values and (optional) auxiliary target */
/* input values. Return number of results and clipping flag. If return value == mxsoln, */
/* then there might be more results. The target values returned will correspond to the */
/* actual (posssibly clipped) point. The return value is the number of solutions + */
/* a clipped flag. Properly set hint flags improve performance, but a correct result should */
/* be returned if the RSPL_NEARCLIP is set, even if they are not set correctly. */
static int
rev_interp_rspl(
	rspl *s,		/* this */
	int flags,		/* Hint flag */
	int mxsoln,		/* Maximum number of solutions allowed for */
	int *auxm,		/* Array of di mask flags, !=0 for valid auxliaries (NULL if no auxiliaries) */
	double cdir[MXRO],	/* Clip vector direction wrt to cpp[0].v and length - NULL if not used */
	co *cpp			/* Given target output space value in cpp[0].v[] +  */
					/* target input space auxiliaries in cpp[0].p[], return */
					/* input space solutions in cpp[0..retval-1].p[], and */
) {
	int e, di = s->di;
	int f, fdi = s->fdi;
	int i, *rip = NULL;
	cell *cells;			/* Candidate cells */
	schbase *b = NULL;		/* Base search information */
	double auxv[MXRI];		/* Locus proportional auxiliary values */
	int didclip = 0;		/* flag - set if we clipped the target */
	
	DBGV(("\nrev interp called with out targets", fdi, " %f", cpp[0].v, "\n"));

	/* This is a restricted size function */
	if (di > MXRI)
		error("rspl: rev_interp can't handle di = %d",di);
	if (fdi > MXRO)
		error("rspl: rev_interp can't handle fdi = %d",fdi);

	if (auxm != NULL) {
		double ax[MXRI];
		for (i = 0; i < di; i++) {
			if (auxm[i] != 0)
				ax[i] = cpp[0].p[i];
			else
				ax[i] = 0.0;
		}
		DBGV(("                  auxiliaries mask", di, " %d", auxm, "\n"));
		DBGV(("                auxiliaries values", di, " %f", ax, "\n"));
	}
	DBG(("di = %d, fdi = %d\n",di, fdi));
	DBG(("flags = 0x%x\n",flags));

	mxsoln &= RSPL_NOSOLNS;		/* Prevent silliness */

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Auxiliary is proportion of locus, so we need to find locus extent */	
	if (flags & RSPL_AUXLOCUS) {
		DBG(("rev interp has aux targets as proportion of locus\n"));

		flags &= ~RSPL_WILLCLIP;		/* Reset hint flag, as we will figure it out */

		/* For each valid auxiliary */
		for (e = 0; e < di; e++) {
			if (auxm[e] == 0)
				continue;			/* Skip unsused auxiliaries */
	
			/* Do search for min and max */
			DBG(("rev locus searching for aux %d min/max\n", e));
			if (b == NULL) {
				b = init_search(s, flags, cpp[0].p, auxm, cpp[0].v, cdir, cpp, mxsoln, locus);
#ifdef STATS
				s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
			} else
				set_lsearch(s, e);		/* Reset locus search for next auxiliary */

			if (rip == NULL) {		/* Not done this yet */
				rip = calc_fwd_cell_list(s, cpp[0].v); /* Reverse grid index for out target */
				if (rip == NULL) {
					DBG(("Got NULL list (point outside range) for auxiliary locus search\n"));
					flags |= RSPL_WILLCLIP2;
					break;
				}
			}
	
			search_list(b, rip, s->get_next_touch(s)); /* Setup, sort and search the list */
	
			if (b->min > b->max) {			/* Failed to find locus */
				DBG(("rev interp failed to find locus for aux %d, so expect clip\n",e));
				flags |= RSPL_WILLCLIP2;
				break;
			}
			auxv[e] = (cpp[0].p[e] * (b->max - b->min)) + b->min;
		}

		DBG(("rev interp got all locuses, so expect exact result\n",e));
		if (!(flags & RSPL_WILLCLIP)) {
			flags |= RSPL_EXACTAUX;				/* Got locuses, so expect exact result */
		}
	}

	/* Init the search information */
	if (b == NULL)
		b = init_search(s, flags, cpp[0].p, auxm, cpp[0].v, cdir, cpp, mxsoln, exact);
	else
		adjust_search(s, flags, auxv, exact);		/* Using proportion of locus aux */
	
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* If hinted that we will not need to clip, look for exact solution. */
	if (!(flags & RSPL_WILLCLIP)) {
		DBG(("Trying exact search\n"));

		/* First do an exact search (init will select auxil if requested) */
		adjust_search(s, flags, NULL, exact);
	
		/* Figure out the reverse grid index appropriate for this request */
		if (rip == NULL)	/* Not done this yet */
			rip = calc_fwd_cell_list(s, cpp[0].v);
	
#ifdef STATS
		s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
		if (rip != NULL) {
			/* Setup, sort and search the list */
			search_list(b, rip, s->get_next_touch(s));
		} else {
			DBG(("Got NULL list (point outside range) for first exact reverse cell\n"));
		}
	
		/* If we selected exact aux, but failed to find a solution, relax expectation */
		if (b->nsoln == 0 & b->naux > 0 && (flags & RSPL_EXACTAUX)) {
			DBG(("Searching for exact match to auxiliary target failed, so try again\n"));
			adjust_search(s, flags & ~RSPL_EXACTAUX, NULL, exact);

#ifdef STATS
			s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
			/* Candidate cell list should be the same */
			if (rip != NULL) {
				/* Setup, sort and search the list */
				search_list(b, rip, s->get_next_touch(s));
			} else {
				DBG(("Got NULL list (point outside range) for nearest search reverse cell\n"));
			}
		}
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* If the exact search failed, and we should look for a nearest solution */
	if (b->nsoln == 0 && (flags & RSPL_NEARCLIP)) {
		unsigned int tcount;	/* grid touch count for this opperation */

		adjust_search(s, flags, NULL, clipn);
		tcount = s->get_next_touch(s);		/* Get next grid touched generation count */

#ifdef STATS
		s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
		/* We get returned a list of cube base indexes of all cubes that have */
		/* the closest valid vertex value to the target value. */
		/* (This may not result in the true closest point if the geometry of */
		/* the vertex values is localy non-smooth or self interesecting, */
		/* but seems to return a good result in most realistic situations ?) */

#ifdef NEVER
		rip = calc_fwd_cell_list(s, cpp[0].v);		/* Try proximity of target */
		if (rip != NULL)
			search_list(b, rip, tcount);
#endif /* NEVER */
	
#ifdef SSS
		{
		sphere sh;				/* Structure to hold sphere context */
		unsigned int tcount;	/* grid touch count for this opperation */

		DBG(("Starting a clipping nearest search now!!\n"));

		adjust_search(s, flags, NULL, clipn);

		tcount = s->get_next_touch(s);		/* Get next grid touched generation count */

#ifdef STATS
		s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
		/* We need to search reverse grid cells in an expanding sphere, */
		/* until we run out of cells, or we have looked at all cells that */
		/* are as close or closer to the current best solution. */

#ifdef NEVER	// ~2
		rip = init_sphere(s, &sh, cpp[0].v);	/* Init the sphere cell searcher */
		for (; sh.ldist <= b->cdist ;
		                    rip = next_sphere_cell(&sh)) {
			if (rip == NULL || sh.cdist > b->cdist) {
				DBG(("Got NULL list or worse possible for this reverse cell\n"));
				continue;
			}

			/* Setup, sort and search the list */
			search_list(b, rip, tcount);
		}
#else
		{	/* Exaustive search */
		int **rpp;
		for (rpp = s->rev.rev; rpp < (s->rev.rev + s->rev.no); rpp++) {
			rip = *rpp;
			if (rip != NULL) {
				search_list(b, rip, tcount);
			}
		}
		}
#endif

		}
#else /* !SSS */
		/* Get list of cells enclosing nearest vertex */
		if ((rip = get_nearest_list(s, cpp[0].v)) != NULL)
			search_list(b, rip, tcount); /* Setup, sort and search the list */
#endif

		if (b->nsoln > 0)
			didclip = RSPL_DIDCLIP;
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* If we still don't have a solution, do a vector direction clip */
	if (b->nsoln == 0 && b->canvecclip) {
		/* Find clipping solution in vector direction */
		line ln;				/* Structure to hold line context */
		unsigned int tcount;	/* grid touch count for this opperation */

		DBG(("Starting a clipping vector search now!!\n"));

		adjust_search(s, flags, NULL, clipv);

		tcount = s->get_next_touch(s);		/* Get next grid touched generation count */

#ifdef STATS
		s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
		init_line_eq(b, b->v, cdir);				/* Init the implicit line equation */
		rip = init_line(s, &ln, cpp[0].v, cdir);	/* Init the line cell dda */
//~~1 HACK!!! should be <= 1.0 !!!
		for (; ln.t <= 2.0;
		                    rip = next_line_cell(&ln)) {
			if (rip == NULL) {
				DBG(("Got NULL list for this reverse cell\n"));
				continue;
			}

			/* Setup, sort and search the list */
			search_list(b, rip, tcount);

			/* If we have found a solution, then abort the search - */
			/* this line will be taking us away from the best solution. */
			if (b->nsoln > 0)
				break;
		}
		if (b->nsoln > 0)
			didclip = RSPL_DIDCLIP;
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* If the clipped solution seems to have been jumping to conclusions, */
	/* search for an exact solution. */
	if (didclip && (flags & RSPL_WILLCLIP && !(flags & RSPL_CERTAIN))
	 && (b->cdist/s->get_out_scale(s)) < 0.002) {
		co c_cpp       = b->cpp[0];	/* Save clip solution in case we want it */
		double c_idist = b->idist;	
		int c_nsoln    = b->nsoln;
		int c_pauxcell = b->pauxcell;
		int c_cdist    = b->cdist;
		int c_iclip    = b->iclip;

		DBG(("Trying exact search again\n"));

		/* Do an exact search (init will select auxil if requested) */
		adjust_search(s, flags & ~RSPL_WILLCLIP, NULL, exact);
	
		/* Figure out the reverse grid index appropriate for this request */
		rip = calc_fwd_cell_list(s, cpp[0].v);
	
#ifdef STATS
		s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
		if (rip != NULL) {
			/* Setup, sort and search the list */
			search_list(b, rip, s->get_next_touch(s));
		} else {
			DBG(("Got NULL list (point outside range) for first exact reverse cell\n"));
		}
	
		/* If we selected exact aux, but failed to find a solution, relax expectation */
		if (b->nsoln == 0 & b->naux > 0 && (flags & RSPL_EXACTAUX)) {
			DBG(("Searching for exact match to auxiliary target failed, so try again\n"));
			adjust_search(s, flags & ~RSPL_EXACTAUX, NULL, exact);

#ifdef STATS
			s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
			/* Candidate cell list should be the same */
			if (rip != NULL) {
				/* Setup, sort and search the list */
				search_list(b, rip, s->get_next_touch(s));
			} else {
				DBG(("Got NULL list (point outside range) for nearest search reverse cell\n"));
			}
		}

		/* If we did get an exact solution */
		if (b->nsoln > 0) {
			DBG(("Deciding to return exact solution after finding clipped\n"));
			didclip = 0;		/* Reset did-clip and return exact solution */

		} else {
			DBG(("keeping clipped solution\n"));
			/* Restore the clipped solution */
			b->cpp[0] = c_cpp;
			b->idist = c_idist;	
			b->nsoln = c_nsoln;
			b->pauxcell = c_pauxcell;
			b->cdist = c_cdist;
			b->iclip = c_iclip;
		}
	}

	if (b->nsoln > 0) {
		DBGV(("rev interp returning 1st soln: ",di," %f", cpp[0].p, "\n"));
	}
	DBG(("rev interp returning %d solutions%s\n",b->nsoln, didclip ? " [clip]" : ""));

	return b->nsoln | didclip;
}

/* ------------------------------------------------------------------------------------ */
/* Do reverse search for the auxiliary min/max ranges of the solution locus for the */
/* given target output values. */
/* Return number of locus segments found, up to mxsoln. 0 will be returned if no solutions */
/* are found. */

static int
rev_locus_segs_rspl (
	rspl *s,		/* this */
	int *auxm,		/* Array of di mask flags, !=0 for valid auxliaries (NULL if no auxiliaries) */
	co *cpp,		/* Input value in cpp[0].v[] */
	int mxsoln,		/* Maximum number of solutions allowed for */
	double min[][MXRI],	/* Array of min[MXRI] to hold return segment minimum values. */
	double max[][MXRI]	/* Array of max[MXRI] to hold return segment maximum values. */
) {
	int e, di = s->di;
	int f, fdi = s->fdi;
	int six;		/* solution index */
	int *rip = NULL;
	int rv = 1;				/* Return value */
	schbase *b = NULL;		/* Base search information */
	
	DBGV(("rev locus called with out targets", fdi, " %f", cpp[0].v, "\n"));
	
	/* This is a restricted size function */
	if (di > MXRI)
		error("rspl: rev_locus_segs can't handle di = %d",di);
	if (fdi > MXRO)
		error("rspl: rev_locus_segs can't handle fdi = %d",fdi);

	if (mxsoln < 1)
		return 0;			/* Guard against silliness */

	if (auxm != NULL) {
		int i;
		double ax[MXRI];
		for (i = 0; i < di; i++) {
			if (auxm[i] != 0)
				ax[i] = cpp[0].p[i];
			else
				ax[i] = 0.0;
		}
		DBGV(("                  auxiliaries mask", di, " %d", auxm, "\n"));
		DBGV(("                auxiliaries values", di, " %f", ax, "\n"));
	}

	/* Init default return values */
	for (six = 0; six < mxsoln; six++) {
		for (e = 0; e < di; e++) {
			if (auxm[e] == 0) {
				min[six][e] = max[six][e] = 0;	/* Return 0 for unused auxiliaries */
			} else {
				min[six][e] = 1.0;			/* max < min indicates invalid range */
				max[six][e] = 0.0;
			}
		}
	}

	/* For each valid auxiliary */
	for (e = 0; e < di; e++) {
		if (auxm[e] == 0)
			continue;			/* Skip unsused auxiliaries */

		/* Do search for min and max */
		DBG(("rev locus searching for aux %d min/max\n", e));
		if (b == NULL)
			b = init_search(s, 0, cpp[0].p, auxm, cpp[0].v, NULL, cpp, mxsoln, locus);
		else
			set_lsearch(s, e);		/* Reset locus search for next auxiliary */

		if (rip == NULL) {		/* Not done this yet */
			rip = calc_fwd_cell_list(s, cpp[0].v); /* Reverse grid index for this request */
			if (rip == NULL) {
				DBG(("Got NULL list (point outside range) for auxiliary locus search\n"));
				rv = 0;
				break;
			}
		}

		search_list(b, rip, s->get_next_touch(s)); /* Setup, sort and search the list */

		if (b->min > b->max) {
			rv = 0;				/* Failed to find a result */
			break;
		}

		if (b->asegs == 0) {		/* Overall min max only */

			min[0][e] = b->min;		/* Save single result */
			max[0][e] = b->max;

		} else {				/* Tracking auxiliary segments */
			int si;					/* Start i */
			int i, j, ff;

			/* Sort the segment list */
#define 	HEAP_COMPARE(A,B) (A.xval < B.xval)
			HEAPSORT(axisec, b->axisl, b->axisln)
#undef 		HEAP_COMPARE

#ifdef NEVER
for (i = 0; i < b->axisln; i++) {
printf("~2 xval = %f, verts = ",b->axisl[i].xval);
for (f = 0; f < b->axisl[i].nv; f++)
printf(" %d", b->axisl[i].vix[f]);
printf("\n");
}
#endif
			/* Find the segments by finding common verticies */
			six = si = i = 0;

			min[six][e] = b->axisl[i].xval;

			for (i++; i < (b->axisln-1); i++) {
				/* Check if any i and i-1 to j are connected */
				for (j = i-1; j >= si; j--) {
					for (f = 0; f < b->axisl[j].nv; f++) {
						for (ff = 0; ff < b->axisl[i].nv; ff++) {
							if (b->axisl[j].vix[f] == b->axisl[i].vix[ff])
								break;		/* Found a link */
						}
						if (ff < b->axisl[i].nv)
							break;
					}
					if (f < b->axisl[j].nv)
						break;
				}
				if (j < si) {	/* Wasn't linked */
					int ii, jj;
					/* Think we found a break. Check that all the rest of */
					/* the entries don't have any links to the previous group */

					/* This could be rather a slow way of checking ! (On^2) */
					for (ii = i+1; ii < (b->axisln); ii++) {
						for (jj = i-1; jj >= si; jj--) {
							for (f = 0; f < b->axisl[jj].nv; f++) {
								for (ff = 0; ff < b->axisl[ii].nv; ff++) {
									if (b->axisl[jj].vix[f] == b->axisl[ii].vix[ff])
										break;		/* Found a link */
								}
								if (ff < b->axisl[ii].nv)
									break;
							}
							if (f < b->axisl[jj].nv)
								break;
						}
						if (jj >= si)
							break;
					}
					if (ii >= b->axisln) {	/* Wasn't forward linked */
						/* Nothing ahead links to last group */
						max[six][e] = b->axisl[i-1].xval;

						/* If we run out of solution space */
						/* merge the last segments */
						if ((six+1) < mxsoln) {
							six++;
							min[six][e] = b->axisl[i].xval;
						}
					}
				}
			}
			max[six++][e] = b->axisl[i].xval;

			if (six > rv)
				rv = six;
		}
	}

#ifdef STATS
	s->rev.st[b->op].searchcalls++;
#endif	/* STATS */
	if (rv) {
		for (six = 0; six < rv; six++) {
			DBG(("rev locus returning:\n"));
			DBGV(("     min", di, " %f", min[six], "\n"));
			DBGV(("     max", di, " %f", max[six], "\n"));
		}
	}

	DBG(("rev locus returning status %d\n",rv));
	return rv;
}

/* ------------------------------------------------------------------------------------ */
typedef double mxdi_ary[MXRI];

/* Do reverse search for the locus of the auxiliary input values given a target output. */
/* Return 1 on finding a valid solution, and 0 if no solutions are found. */
static int
rev_locus_rspl(
	rspl *s,		/* this */
	int *auxm,		/* Array of di mask flags, !=0 for valid auxliaries (NULL if no auxiliaries) */
	co *cpp,		/* Input value in cpp[0].v[] */
	double min[MXRI],/* Return minimum auxiliary values */
	double max[MXRI] /* Return maximum auxiliary values */
) {

	/* Use segment routine to compute oveall locus */
	return rev_locus_segs_rspl (s, auxm, cpp, 1, (mxdi_ary *)min, (mxdi_ary *)max);
}

/* ------------------------------------------------------------------------------------ */

#undef DEBUG	// ~~1

#ifdef DEBUG
#undef DBG
#undef DBGV
#undef DBGM
#define DBG(xxx) DBGI(xxx)
#define DBGV(xxx) DBGVI xxx
#define DBGM(xxx) DBGMI xxx
#else
#undef DBG
#undef DBGV
#undef DBGM
#define DBG(xxx) 
#define DBGV(xxx) 
#define DBGM(xxx) 
#endif

/* ------------------------------------------------ */
/* subroutines of top level reverse lookup routine */

static int exact_setsort(schbase *b, cell *c);
static int exact_compute(schbase *b, simplex *x);

static int auxil_setsort(schbase *b, cell *c);
static int auxil_check(schbase *b, cell *c);
static int auxil_compute(schbase *b, simplex *x);

static int locus_setsort(schbase *b, cell *c);
static int locus_check(schbase *b, cell *c);
static int locus_compute(schbase *b, simplex *x);

static int clipv_setsort(schbase *b, cell *c);
static int clipv_check(schbase *b, cell *c);
static int clipv_compute(schbase *b, simplex *x);

static int clipn_setsort(schbase *b, cell *c);
static int clipn_check(schbase *b, cell *c);
static int clipn_compute(schbase *b, simplex *x);

/* Allocate the search base structure */
static schbase *
alloc_sb(rspl *s) {
	schbase *b;
	if ((b = s->rev.sb = (schbase *)calloc(1, sizeof(schbase))) == NULL)
		error("rspl malloc failed - rev.sb structure");

	b->s     = s;				/* rsp */
	b->pauxcell =				/* Previous solution cell indexes */
	b->plmaxcell = 
	b->plmincell = -1;
	return b;
}

/* Do the basic search type independent initialization */
static schbase *	/* Return pointer to base search information */
init_search(
	rspl *s,		/* rsp; */
	int flags,		/* Hint flag */

	double *av,		/* Auxiliary input values - may be NULL */
	int *auxm,		/* Array of di mask flags, !=0 for valid auxliaries (NULL if no auxiliaries) */
					/* Locus search will search for max/min of first valid auxlilary */
	double *v,		/* Output value target, NULL if none */
	double *cdir,	/* Clip vector direction, NULL if none */
	co *cpp,		/* Array that hold solutions, NULL if none. */
	int mxsoln,		/* Maximum number of solutions allowed for */
	enum ops op		/* Type of reverse search operation requested */
) {
	schbase *b = NULL;		/* Pointer to search base information structure */
	int e, di = s->di;
	int f, fdi = s->fdi;

	DBG(("Initializing search\n"));

	if (s->rev.inited == 0) 		/* Compute reverse info if it doesn't exist */
		make_rev(s);

	/* If first time initialisation (Third section init) */
	if ((b = s->rev.sb) == NULL) {
		b = alloc_sb(s);
	}

	/* Init some basic search info */
	b->op    = op;				/* operation */
	b->flags = flags;			/* hint flags */
	b->canvecclip = 0;			/* Assume invalid clip direction */

	b->ixc = (1<<di)-1;			/* Cube index of corner that holds maximum input values */

	/* Figure out if auxiliaries have been requested */
	b->naux = 0;
	b->auxbm = 0;
	if (auxm != NULL) {
		unsigned bm;

		if (mxsoln > 1)
			b->asegs = 1;		/* Find all segments */
		else
			b->asegs = 0;		/* Find only overall aux locus range */

		for (e = di-1, bm = 1 << e; e >= 0; e--, bm >>= 1) {	/* Record auxiliary mask bits */
			if (av != NULL)
				b->av[e] = av[e];	/* Auxiliary target values */
			b->auxm[e] = auxm[e];	/* Auxiliary mask */
			if (auxm[e] != 0) {
				b->auxbm |= bm;			/* Auxiliary bit mask */
				b->auxi[b->naux++] = e;	/* Index of next auxiliary input to be used */
				/* Auxiliary locus extent */
				b->lxi = e;			/* Assume first one */
				b->max = -INF_DIST;	/* In case searching for max */
				b->min =  INF_DIST;	/* In case searching for minimum */
				b->axisln = 0;		/* No intersects in list */
			}
		}
	}

	/* Figure out if the clip direction is meaningfull */
	/* Check that the clip vector makes sense */
	if (cdir != NULL) {	/* Clip vector is specified */
		double ss;
		for (ss = 0.0, f = 0; f < fdi; f++) {
			double tt = cdir[f];
			b->cdir[f] = tt;
			ss += tt * tt;
		}

		if (ss > 1e-6) {
			b->canvecclip = 1;	/* It has a non-zero length */
			ss = sqrt(ss);
			/* Compute normalised clip vector direction */
			for (f = 0; f < fdi; f++) {
				b->ncdir[f] = b->cdir[f]/ss;
			}
		}
	}

	if (di <= fdi)		/* Only allow auxiliaries if di > fdi */
		b->naux = 0;

	/* Switch to appropriate operation */
	if (b->op == exact && (b->naux > 0 || di != fdi)) {
		b->op = auxil;
	} else if (b->op == auxil && b->naux == 0 && di == fdi) {
		b->op = exact;
	}

	/* Set appropriate functions for type of operation */
	switch (b->op) {
		case exact:
			b->setsort = exact_setsort;
			b->check   = NULL;
			b->compute = exact_compute;
			b->snsdi = b->ensdi = di;	/* Search full dimension simplex, expect point soln. */
			break;
		case auxil:
			b->setsort = auxil_setsort;
			b->check   = auxil_check;
			b->compute = auxil_compute;
			b->snsdi = di;				/* Start here DOF = di-fdi locus solutions */
			b->ensdi = fdi;				/* End with DOF = 0 for point solutions */
			break;
		case locus:
			b->setsort = locus_setsort;
			b->check   = locus_check;
			b->compute = locus_compute;
			b->snsdi = b->ensdi = fdi;	/* Search for point solutions */
			break;
		case clipv:
			b->setsort = clipv_setsort;
			b->check   = clipv_check;
			b->compute = clipv_compute;
											/* Clip vector 1 dimension in output space, */
			b->snsdi = b->ensdi = fdi-1;	/* search planes for combined point solution */
			break;
		case clipn:
			b->setsort = clipn_setsort;
			b->check   = clipn_check;
			b->compute = clipn_compute;
			b->snsdi = 0;				/* Start with DOF = 0 for point solutions */
			b->ensdi = fdi-1;			/* End on DOF = di-fdi-1 on surfaces of simplexes */
			break;
		default:
			error("init_search: Unknown operation %d\n",b->op);
	}

	if (v != NULL) {
		for (f = 0; f < fdi; f++)	/* Record target output values */
			b->v[f] = v[f];
		b->v[fdi] = b->limitv;		/* Limitvalue is output target for limit clip subsimplexes */
	}

	b->mxsoln = mxsoln;				/* Allow solutions to be returned */
	b->cpp    = cpp;				/* Put solutions here */
	b->nsoln = 0;					/* No solutions at present */
	b->iclip = 0;					/* Default solution isn't above ink limit */

	if (flags & RSPL_EXACTAUX)		/* Expect to be able to match auxiliary target exactly */
		b->idist = 0.00001;			/* Best input distance to beat - helps sort/triage */
	else
		b->idist = INF_DIST;		/* Best input distance to beat. */

	b->cdist = INF_DIST;			/* Best clip distance to beat. */

	DBG(("Search initialized\n"));

	return b;
}

/* Adjust the search */
static void
adjust_search(
	rspl *s,		/* rsp; */
	int flags,		/* Hint flag */
	double *av,		/* Auxiliary input values - may be NULL */
	enum ops op		/* Type of reverse search operation requested */
) {
	schbase *b = s->rev.sb;		/* Pointer to search base information structure */
	int e, di = s->di;
	int f, fdi = s->fdi;

	DBG(("Adjusting search\n"));

	b->op    = op;				/* operation */
	b->flags = flags;			/* hint flags */

	/* Switch to appropriate operation */
	if (b->op == exact && (b->naux > 0 || di != fdi)) {
		b->op = auxil;
	} else if (b->op == auxil && b->naux == 0 && di == fdi) {
		b->op = exact;
	}

	/* Update auxiliary target values */
	if (av != NULL) {
		for (e = 0; e < b->naux; e++) {
			int ee = b->auxi[e];
			b->av[ee] = av[ee];
		}
	}

	/* Set appropriate functions for type of operation */
	switch (b->op) {
		case exact:
			b->setsort = exact_setsort;
			b->check   = NULL;
			b->compute = exact_compute;
			b->snsdi = b->ensdi = di;		/* Expect point solution */
			break;
		case auxil:
			b->setsort = auxil_setsort;
			b->check   = auxil_check;
			b->compute = auxil_compute;
			b->snsdi = di;				/* Start here DOF = di-fdi locus solutions */
			b->ensdi = fdi;				/* End with DOF = 0 for point solutions, */
			break;						/* will early exit DOF. */
		case locus:
			b->setsort = locus_setsort;
			b->check   = locus_check;
			b->compute = locus_compute;
			b->snsdi = b->ensdi = fdi;	/* Search for point solutions */
			break;
		case clipv:
			b->setsort = clipv_setsort;
			b->check   = clipv_check;
			b->compute = clipv_compute;
											/* Clip vector 1 dimension in output space, */
			b->snsdi = b->ensdi = fdi-1;	/* so the intersection with the simplex is a point. */
			break;
		case clipn:
			b->setsort = clipn_setsort;
			b->check   = clipn_check;
			b->compute = clipn_compute;
			b->snsdi = 0;				/* Start with DOF = 0 for point solutions */
			b->ensdi = fdi-1;			/* End on DOF = di-fdi-1 on surfaces of simplexes */
			break;						/* Will go through all DOF */
		default:
			error("init_search: Unknown operation %d\n",b->op);
	}

	b->nsoln = 0;					/* No solutions at present */

	if (flags & RSPL_EXACTAUX)		/* Expect to be able to match auxiliary target exactly */
		b->idist = 0.00001;			/* Best input distance to beat - helps sort/triage */
	else
		b->idist = INF_DIST;		/* Best input distance to beat. */

	b->cdist = INF_DIST;			/* Best clip distance to beat. */

	DBG(("Search adjusted\n"));
}

/* Adjust existing locus search for a different auxiliary */
static void
set_lsearch(
rspl *s,
int e			/* Next auxiliary */
) {
	schbase *b = s->rev.sb;		/* Pointer to search base information structure */

	b->lxi = e;			/* Assume first one */
	b->max = -INF_DIST;	/* In case searching for max */
	b->min =  INF_DIST;	/* In case searching for minimum */
	b->axisln = 0;		/* No intersects in list */
}

/* Set the limit search information */
/* Note this doesn't create or init the main rev information. */
static schbase *	/* Return pointer to base search information */
set_search_limit(
	rspl *s,		/* rsp; */
	double (*limit)(void *vcntx, double *in),	/* Optional input space limit function. Function */
					/* should evaluate in[0..di-1], and return number that is not to exceed */
					/* limitv. NULL if not used */
	void *lcntx,	/* Context passed to limit() */
	double limitv	/* Value that limit() is not to exceed */
) {
	schbase *b = NULL;		/* Pointer to search base information structure */

	/* If sb info needs initialising (Third section init) */
	if ((b = s->rev.sb) == NULL) {
		b = alloc_sb(s);
	}

	b->limit = limit;			/* Input limit function */
	b->cntx  = lcntx; 			/* Context passed to limit() */
	b->limitv= INKSCALE * limitv; 	/* Context passed to values not to be exceedded by limit() */
	if (limit != NULL)
		b->limiten = 1;				/* enable limiting by default */
	else
		b->limiten = 0;				/* No limit function, so limiting not enabled. */

	return b;
}

/* Free any search specific data */
static void
free_search(
schbase *b	/* Base search information */
) {
	DBG(("Freeing search\n"));

	/* Clip line implicit equation (incuding space for ink target) */
	if (b->cla != NULL) {
		int fdi = b->s->fdi;
		free_dmatrix(b->cla, 0, fdi-1, 0, fdi);
		b->cla = NULL;
	}

	/* Auxiliary segment list */
	if (b->axisl > 0) {
		free (b->axisl);
		b->axisl = NULL;
		b->axislz = 0;
		b->axisln = 0;
	}

	/* Sorted cell list */
	if (b->lclistz > 0) {
		free (b->lclist);
		b->lclist = NULL;
		b->lclistz = 0;
	}

	/* Simplex filter list */
	if (b->lsxfilt > 0) {
		free (b->sxfilt);
		b->sxfilt = NULL;
		b->lsxfilt = 0;
	}

	/* nn cell list */
	if (b->nnlistz > 0) {
		free (b->nnlist);
		b->nnlist = NULL;
		b->nnlistz = 0;
	}
	free(b);
}

/* Return the pointer to the list of fwd cells given */
/* the target output values. Return NULL if none in list. */
static int *
calc_fwd_cell_list(
	rspl *s,		/* this */
	double *v		/* Output values */
) {
	int f, fdi = s->fdi;
	int **rpp;
	int rgres_1 = s->rev.res - 1;

	for (rpp = s->rev.rev, f = 0; f < fdi; f++) {
		int mi;
		double gw = s->rev.gw[f];
		double t = (v[f] - s->rev.gl[f])/gw;
		mi = (int)floor(t);				/* Grid coordinate */
		if (mi < 0 || mi > rgres_1) { 	/* If outside valid reverse range */
			return NULL;
		}
		rpp += mi * s->rev.coi[f];	/* Accumulate reverse grid pointer */
	}
	return *rpp;
}

static int rev_sch(schbase *b, cell *c, int level);
void alloc_simplexes(cell *c, int nsdi);

/* Given a pointer to a list of fwd cells, cull cells that
/* cannot contain or improve the solution, sort the list, */
/* and then compute the final best solution. */
static void
search_list(
schbase *b,				/* Base search information */
int     *rip,			/* Pointer to base of cell list, entry 0 = allocated space */
unsigned int tcount		/* grid touch count for this opperation */
) {
	rspl *s = b->s;
	int e, di = s->di;
	int nsdi;
	int i;
	int nilist;			/* Number in cell list */
	
	DBG(("search_list called\n"));

	/* (rip[0] contains number of fwd cells in the list) */
	if (b->lclistz < rip[0]) {	/* Allocate more space if needed */

		if (b->lclistz > 0) {	/* Free old space before allocating new */
			free (b->lclist);
		}
		b->lclistz = 0;
		/* Allocate enough space for all the candidate cells */
		if ((b->lclist = (cell **)malloc(rip[0] * sizeof(cell *))) == NULL)
			error("rev: malloc failed - candidate cell list, count %d",rip[0]);
		b->lclistz = rip[0];	/* Current allocated space */
	}
		
	/* For each chunk of the list that we can fit in the rcache: */
	for(rip++; *rip != -1;)  {

		/* Go through all the candidate fwd cells, and build up the list of search cells */
		for(nilist = 0; *rip != -1; rip++)  {
			int ix = *rip;				/* Fwd cell index */
			float *fcb = s->g.a + ix * s->g.pss;	/* Pointer to base float of fwd cell */
			cell *c;

			if (TOUCHF(fcb) >= tcount) {	/* If we have visited this cell before */
				DBG((" Already touched cell index %d\n",ix));
				continue;
			}
			/* Get pointers to cells from cache, and lock it in the cache */
			if ((c = get_rcell(b, ix)) == NULL) {
				static warned = 0;
				if (!warned) {
					warning("Reverse Cell Cache exausted, processing in chunks");
					warned = 1;
				}
				DBG(("revcache is exausted, do search in chunks\n"));
				break;		/* cache has run out of room, so abandon, and do it next time */
			}

			DBG(("checking out cell %d\n",ix));
			TOUCHF(fcb) = tcount;			/* Touch it */

			/* Check mandatory conditions, and compute search key */
			if (!b->setsort(b, c)) {
				DBG(("cell %d rejected from list\n",ix));
				unget_rcell(s->rev.cache, c);
				continue;
			}
			DBG(("cell %d accepted into list\n",ix));

			b->lclist[nilist++] = c; /* Cell is accepted as recursion candidate */
		}

		if (nilist == 0) {
			DBG(("List was empty\n"));
		}

#ifdef DOSORT
		/* If appropriate, sort child cells into best order */
		/* == sort key smallest to largest */
		switch (b->op) {
			case locus:
				{	/* Special case, adjust sort values */
					double min = INF_DIST, max = -INF_DIST;
					for (i = 0; i < nilist; i++) {
						cell *c = b->lclist[i];
						if (c->sort < min)
							min = c->sort;
						if (c->sort > max)
							max = c->sort;
					}
					max = min + max;	/* Total of min/max */
					min = 0.5 * max;	/* Average sort value */
					for (i = 0; i < nilist; i++) {
						cell *c = b->lclist[i];
						if (c->ix == b->plmincell || c->ix == b->plmaxcell) {
							c->sort = -1.0;		/* Put previous solution cells at head of list */
						} else if (c->sort > min) {
							c->sort = max - c->sort;	/* Reflect about average */
						}
					}
				}
				/* Fall through to sort */
			case auxil:
			case clipv:
			case clipn:
#define 	HEAP_COMPARE(A,B) (A->sort < B->sort)
				HEAPSORT(cell *,b->lclist, nilist)
#undef 		HEAP_COMPARE
				break;
		}
#endif /* DOSORT */

		DBG(("List sorted, about to search\n"));
#ifdef NEVER
		printf("\n~1 Op = %s, Cell sort\n",opnames[b->op]);
		for (i = 0; i < nilist; i++) {
			printf("~1 List %d, cell %d, sort = %f\n",i,b->lclist[i]->ix,b->lclist[i]->sort);
		}
#endif /* NEVER */

		/* 
			Tried reversing the "for each cell" and "for each level" loops,
			but it made a negligible difference to the performance.
		 */

		/* For each dimensionality of sub-simplexes, in given order */
		DBG(("Searching from level %d to level %d\n",b->snsdi, b->ensdi));
		for (nsdi = b->snsdi;;) {

			DBG(("\n******************\n"));
			DBG(("Searching level %d\n",nsdi));

			/* For each cell in the list */
			for (i = 0; i < nilist; i++) {
				cell *c = b->lclist[i];
				int j, nospx;		/* Number of simplexes in cell */

				/* For those searches that have an optimisation goal, */
				/* re-check the cell to see if the goal can still improve on. */
				if (b->check != NULL && !b->check(b, c))
					continue;

#ifdef STATS
				s->rev.st[b->op].csearched++;
#endif /* STATS */

				if (c->sx[nsdi] == NULL) {
					alloc_simplexes(c, nsdi);	/* Do level 1 initialisation for nsdi */
				}
				/* For each simplex in a cell */
				nospx = c->sxno[nsdi];	/* Number of nsdi simplexes */
				for (j = 0; j < nospx; j++) {
					simplex *x = &c->sx[nsdi][j];

					if (b->limiten == 0) {
						if (x->flags & SPLX_CLIPSX)		/* If limiting is disabled, we're */
							continue;					/* not interested in clip plane simplexes */
					} else  {
						if (x->flags & SPLX_OVLIMIT)	/* If limiting is enabled, we're */
							continue;					/* not interested in simplexes over the limit */
					}
		
#ifdef STATS
					s->rev.st[b->op].ssearched++;
#endif /* STATS */
					if (b->compute(b, x)) {
						DBG(("search aborted by compute\n"));
						break;					/* Found enough solutions */
					}
				}	/* Next Simplex */
			}	/* Next cell */

			if (nsdi == b->ensdi)
				break;

			/* Next Simplex dimensionality */
			if (b->ensdi < b->snsdi) {
				if (nsdi == b->snsdi && b->nsoln > 0)	/* Don't continue though decreasing */
					break;			/* sub-simplex dimensions if we found a solution at */
									/* the highest dimension level. */
				nsdi--;
			} else if (b->ensdi > b->snsdi) {
				nsdi++;				/* Continue through increasing sub-simplex dimenionality */
			}						/* until we get to the top. */
		}

		/* Unlock the cache cells now that we're done with them */
		for (i = 0; i < nilist; i++) {
			unget_rcell(s->rev.cache, b->lclist[i]);
		}
	}	/* Next chunk */

	DBG(("search_list complete\n"));
	return;
}

/* ------------------------------------- */
/* Vector search in output space support */

/* Setup the line, and fetch the first cell */
/* Return the pointer to the list of fwd cells, NULL if none in list. */
static int *
init_line(
	rspl *s,			/* this */
	line *l,			/* line structure */
	double st[MXRO],	/* start of line */
	double de[MXRO]		/* line direction and length */
) {
	int f, fdi = s->fdi;
	int **rpp;
	int rgres_1 = s->rev.res - 1;
	int nvalid = 0;		/* Flag set if outside reverse grid range */

	DBGV(("Line from ", fdi, " %f", st, "\n"));
	DBGV(("In dir    ", fdi, " %f", de, "\n"));
	DBGV(("gl        ", fdi, " %f", s->rev.gl, "\n"));
	DBGV(("gh        ", fdi, " %f", s->rev.gh, "\n"));
	DBGV(("gw        ", fdi, " %f", s->rev.gw, "\n"));
	
	/* Init */
	l->s = s;
	for (f = 0; f < fdi; f++) {
		l->st[f] = st[f] - s->rev.gl[f];
		l->de[f] = de[f];
		if (de[f] > 0.0)
			l->di[f] = 1;	/* Axis increments */
		else if (de[f] < 0.0)
			l->di[f] = -1;
		else
			l->di[f] = 0;
	}
	l->t = 0.0;
	DBGV(("increments =", fdi, " %d", l->di, "\n"));

	/* Figure out the starting cell */
	for (rpp = s->rev.rev, f = 0; f < fdi; f++) {
		double gw = s->rev.gw[f];
		double t = l->st[f]/gw;
		l->ci[f] = (int)floor(t);					/* Grid coordinate */
		if (l->ci[f] < 0 || l->ci[f] > rgres_1) 	/* If outside valid reverse range */
			nvalid = 1;
		rpp += l->ci[f] * s->rev.coi[f];	/* Accumulate reverse grid pointer */
	}
	DBGV(("current line cell = ", fdi, " %d", l->ci, "")); DBG((",  t = %f, nvalid = %d\n",l->t,nvalid));
#ifdef DEBUG
{
int ii;
double tt;
printf("Current cell = ");
for (ii = 0; ii < fdi; ii++) {
	tt = l->ci[ii] * s->rev.gw[ii] + s->rev.gl[ii];
	printf(" %f - %f",tt,tt+s->rev.gw[ii]);
}
printf("\n");
}
#endif	/* DEBUG */
	if (nvalid)
		return NULL;
	return *rpp;
}

/* Get the next cell on the line. */
/* Return the pointer to the list of fwd cells, NULL if none in list. */
static int *
next_line_cell(
	line *l		/* line structure */
) {
	rspl *s = l->s;
	int bf = 0, f, fdi = s->fdi;
	int **rpp;
	int rgres_1 = s->rev.res - 1;
	double bt = 100.0;	/* Best (smalest +ve) parameter value to move */

	/* See which axis cell crossing we will hit next */
	for (f = 0; f < fdi; f++) {
		double t;
		double gw = s->rev.gw[f];
		if (l->de[f] != 0) {
			t = ((l->ci[f] + l->di[f]) * gw - l->st[f])/l->de[f];
			DBG(("t for dim %d = %f\n",f,t));
			if (t < bt) {
				bt = t;
				bf = f;		/* Best direction to move */
			}
		}
	}

	/* Move to the next reverse grid coordinate */
	l->ci[bf] += l->di[bf];
	l->t = bt;

	DBGV(("current line cell =", fdi, " %d", l->ci, "")); DBG((",  t = %f\n",l->t));

#ifdef DEBUG
{
int ii;
double tt;
printf("Current cell = ");
for (ii = 0; ii < fdi; ii++) {
	tt = l->ci[ii] * s->rev.gw[ii] + s->rev.gl[ii];
	printf(" %f - %f",tt,tt+s->rev.gw[ii]);
}
printf("\n");
}
#endif	/* DEBUG */

	/* Compute reverse cell index */
	for (rpp = s->rev.rev, f = 0; f < fdi; f++) {
		if (l->ci[f] < 0 || l->ci[f] > rgres_1) { 	/* If outside valid reverse range */
			DBG(("Outside list on dim %d, 0 <= %d <= %d\n", f, l->ci[f],rgres_1));
			return NULL;
		}
		rpp += l->ci[f] * s->rev.coi[f];	/* Accumulate reverse grid pointer */
	}
	return *rpp;
}

/* ------------------------------------- */
/* Clip nearest support. */

int rev_nnfind(rspl *s, double rad, double *v);

/* Using the reverse nearest neighbor vertex search function, */
/* construct a reverse search list */
static int *get_nearest_list(
rspl *s,			/* RSPL that we are part of */
double *v			/* Target output values */
) {
	double vv[MXRO];			/* Closest vertex output value */
	int f, fdi = s->fdi;
	float *fcb;					/* Pointer to base float of fwd cell */
	int ix;

	if ((ix = rev_nnfind(s, 0.0, v)) < 0) 
		return NULL;

	fcb = s->g.a + ix * s->g.pss;	/* Pointer to grid output value */

	for (f = 0; f < fdi; f++)
		vv[f] = fcb[f];				/* Convert from float to double */

	return calc_fwd_cell_list(s, vv);	/* return list of cells enclosing this */
}

#ifdef SSS
/* "Spherical" search in output space support */

/* Initialise the "sphere" counter to the first location */
static int *init_sphere(
rspl *s,			/* RSPL that we are part of */
sphere *sh,			/* Pointer to sphere counter object */
double *base		/* Real base grid cell index */
) {
	int fdi  = s->fdi;
	int gres = s->rev.res;
	int **rpp;
	int f;

	sh->s     = s;
	sh->c     = 0;		/* Start at layer 0 */

	sh->cdist = 0.0;		/* We are in the target cell */
	sh->ldist = 0.0;

	/* Figure out the starting cell */
	for (rpp = s->rev.rev, f = 0; f < fdi; f++) {
		double gw = s->rev.gw[f];
		double t = (base[f] - s->rev.gl[f])/gw;
		sh->loc[f] = sh->base[f] = (int)floor(t);	/* Grid coordinate */
		sh->off[f] = t - (double)sh->base[f];		/* Offset within grid */
		sh->ix[f] = 0;
		if (sh->loc[f] < 0 || sh->loc[f] >= gres) {	/* If outside valid reverse range */
			sh->cdist = 2 * INF_DIST;		/* Give up - can't cope */
			sh->ldist = 2 * INF_DIST;
			DBG(("init_sphere search failed because target is outside grid\n"));
			return NULL;
		}
		rpp += sh->loc[f] * s->rev.coi[f];	/* Accumulate reverse grid pointer */
	}
	DBGV(("first sphere cell = ", fdi, " %d", sh->loc, "\n")); 
	return *rpp;
}

/* Return the next reverse list pointer. */
/* set ridiculously large sh.ldist if we reached end */
static int *
next_sphere_cell(sphere *sh) {
	rspl *s  = sh->s;
	int fdi  = s->fdi;
	int gres = s->rev.res;
	int **rpp;
	int f;

	do {	/* while locations are not within grid, or we're not at the end */

		/* Increment counter */
		for (f = 0; f < fdi; f++) {
			int j;

			if (f == 0) {
				/* See if any non-ls counters are +/- c */
				for (j = f+1; j < fdi; j++) {
					if (sh->ix[j] <= -sh->c || sh->ix[j] >= sh->c)
						break;			/* Found one that is */
				}
				if (j < fdi)					/* Found one that is */
					sh->ix[f]++;				/* Itterate through center */
				else
					sh->ix[f] += 2 * sh->c;		/* Skip center section */
			} else {
				sh->ix[f]++;					/* All non ls increment steadily */
			}
			sh->loc[f] = sh->base[f] + sh->ix[f];

			if (sh->ix[f] <= sh->c)
				break;						/* No carry yet */ 
			sh->ix[f] = -sh->c;
			sh->loc[f] = sh->base[f] + sh->ix[f];
		}
		
		if (f >= fdi) {	/* Carry to layer */

			/* Increment layer */
			sh->c++;

			if (sh->c >= gres) {
				sh->cdist = 2 * INF_DIST;		/* Give up - can't cope */
				sh->ldist = 2 * INF_DIST;
				return NULL;			/* We're at the end */
			}

			for (f = 0; f < fdi; f++) {
				sh->ix[f] = -sh->c;
				sh->loc[f] = sh->base[f] + sh->ix[f];
			}

			/* Compute the best possible distance for this layer */
			for (f = 0; f < fdi; f++) {
				double tt;
				tt = sh->off[f] < 0.5 ? sh->off[f] : 1.0 - sh->off[f];	/* Closest axis */
				tt += (double)sh->c - 1.0;								/* plus other layers */
				tt = tt * s->rev.gw[f] + s->rev.gl[f];		/* Scale back to output space coords */
				if (f == 0 || tt < sh->ldist)
					sh->ldist = tt;
			}
		
		}

		/* Check if this offset is within the grid */
		for (f = 0; f < fdi; f++) {
			if (sh->loc[f] < 0 || sh->loc[f] >= gres)
				break;			/* No its not */
		}
	} while(f < fdi);		/* Until we get one in the grid */

	/* Compute the best possible distance for this cell */
	sh->cdist = 0.0;
	for (f = 0; f < fdi; f++) {

		double tt;
		if (sh->ix[f] > 0)
			tt = (double)sh->ix[f] - sh->off[f];
		else if (sh->ix[f] < 0)
			tt = sh->off[f] - (double)sh->ix[f] - 1.0;
		else	/* == 0 */
			tt = 0.0;
		tt = tt * s->rev.gw[f] + s->rev.gl[f];		/* Scale back to output space coords */
		sh->cdist += tt * tt;
	}
	sh->cdist = sqrt(sh->cdist);

	DBGV(("current sphere cell = ", fdi, " %d", sh->loc, "\n")); 

	/* Compute reverse cell index */
	for (rpp = s->rev.rev, f = 0; f < fdi; f++) {
		rpp += sh->loc[f] * s->rev.coi[f];	/* Accumulate reverse grid pointer */
	}

	return *rpp;
}

#endif /* SSS */

/* =================================================== */
/* The cell and simplex solver top level routines */

static int add_lu_svd(simplex *x);
static int add_locus(schbase *b, simplex *x);
static int add_auxil_lu_svd(schbase *b, simplex *x);
static int within_simplex(simplex *x, double *p);
static void simplex_to_abs(simplex *x, double *in, double *out);

static int auxil_solve(schbase *b, simplex *x, double *xp);

/* ---------------------- */
/* Exact search functions */
/* Return non-zero if cell is acceptable */
static int exact_setsort(schbase *b, cell *c) {
	int ee, di  = b->s->di;
	int f, fdi  = b->s->fdi;
	double ss;

	DBG(("Reverse exact search, evaluate and set sort key on cell\n"));

	/* Check that the target lies within the cell bounding sphere */
	for (ss = 0.0, f = 0; f < fdi; f++) {
		double tt = c->bcent[f] - b->v[f];
		ss += tt * tt;
	}
	if (ss > c->bradsq) {
		DBG(("Cell is rejected - bounding sphere\n"));
		return 0;
	}

	if (b->limiten != 0 && c->limmin > b->limitv) {
		DBG(("Cell is rejected - ink limit, min = %f, limit = %f\n",c->limmin,b->limitv));
		return 0;
	}

	/* Sort can't be used, because we return all solutions */
	c->sort = 0.0;

	DBG(("Cell is accepted\n"));

	return 1;
}

/* Compute a solution for a given sub-simplex (if there is one) */
/* Return 1 if search should be aborted */
static int exact_compute(schbase *b, simplex *x) {
	rspl *s     = b->s;
	cell *c     = x->c;
	int e, di = s->di, sdi  = x->sdi;
	int f, fdi  = s->fdi;
	int i;
	datai xp;	/* solution in simplex relative coord order */
	datai p;	/* absolute solution */
	int wsrv;	/* Within simplex return value */

	DBG(("\nExact: computing possible solution\n"));

#ifdef DEBUG
	/* Sanity check */
	if (sdi != fdi || sdi != di || x->efdi != fdi) {
printf("di = %d, fdi = %d\n",di,fdi);
printf("sdi = %d, efdi = %d\n",sdi,x->efdi);
		error("rspl exact reverse interp called with sdi != fdi, sdi != di, efdi != fdi");
		/* !!! could switch to SVD solution if di != fdi ?? !!! */
	}
#endif

	/* This may not be worth it here since it may not filter out */
	/* many more simplexes than the cube check did. */
	/* This is due to full dimension simplexes all sharing the main */
	/* diagonal axis. */

	/* Check that the target lies within the simplex bounding cube */
	for (f = 0; f < fdi; f++) {
		if (b->v[f] < *x->min[f] || b->v[f] > *x->max[f]) {
			DBG(("Simplex is rejected - bounding cube\n"));
			return 0;
		}
	}

	/* Create the LU decomp needed to exactly solve */
	if (add_lu_svd(x)) {
		DBG(("LU decomp was singular, skip simplex\n"));
		return 0;
	}

	/* Init the RHS B[] vector (note di == fdi) */
	for (f = 0; f < fdi; f++) {
		xp[f] = b->v[f] - x->v[di][f];
	}

	/* Compute the solution (in simplex space) */
	lu_backsub(x->d_u, sdi, (int *)x->d_w, xp);

	/* Check that the solution is within the simplex */
	if ((wsrv = within_simplex(x, xp)) == 0) {
		DBG(("Solution rejected because not in simplex\n"));
		return 0;
	}

	/* Convert solution from simplex relative to absolute space */
	simplex_to_abs(x, p, xp);

	/* Check if a very similiar input solution has been found before */
	for (i = 0; i < b->nsoln; i++) {
		double tt;
		for (e = 0; e < di; e++) {
			tt = b->cpp[i].p[e] - p[e];
			if (fabs(tt) > (2 * EPS))
				break;	/* Mismatch */
		}
		if (e >= di)	/* Found good match */
			break;
	}

	/* Probably alias caused by solution lying close to a simplex boundary */
	if (i < b->nsoln) {
		DBG(("Another solution has been found before - index %d\n",i));
		return 0;		/* Skip this, since betters been found before */
	}

	/* Check we haven't overflowed space */
	if (i >= b->mxsoln) {
		DBG(("Run out of space for new solution\n"));
		return 1;		/* Abort */
	}

	DBG(("######## Accepting new solution\n"));

	/* Put solution in place */
	for (e = 0; e < di; e++)
		b->cpp[i].p[e] = p[e];
	for (f = 0; f < fdi; f++)
		b->cpp[i].v[f] = b->v[f];	/* Assumed to be an exact solution */
	if (i == b->nsoln)
		b->nsoln++;
	if (wsrv == 2)					/* Is above (disabled) ink limit */
		b->iclip = 1;
	return 0;
}

/* -------------------------- */
/* Auxiliary search functions */
static int auxil_setsort(schbase *b, cell *c) {
	int f, fdi  = b->s->fdi;
	int ee, ixc = b->ixc;
	double ss, sort;

	DBG(("Reverse auxiliary search, evaluate and set sort key on cell\n"));

	if (b->s->di <= fdi) {	/* Assert ~1 */
		error("rspl auxiliary reverse interp called with di <= fdi (%d %d)", b->s->di, fdi);
	}

	/* Check that the target lies within the cell bounding sphere */
	for (ss = 0.0, f = 0; f < fdi; f++) {
		double tt = c->bcent[f] - b->v[f];
		ss += tt * tt;
	}
	if (ss > c->bradsq) {
		DBG(("Cell is rejected - bounding sphere\n"));
		return 0;
	}

	if (b->limiten != 0 && c->limmin > b->limitv) {
		DBG(("Cell is rejected - ink limit, min = %f, limit = %f\n",c->limmin,b->limitv));
		return 0;
	}

	/* Check if this cell could possible improve b->idist */
	/* and compute sort key as the distance to auxilliary target */
	/* (We may have a non INF_DIST idist before commencing the */
	/* search if we already know that the auxiliary target is */
	/* within gamut - the usual usage case!) */
	for (sort = 0.0, ee = 0; ee < b->naux; ee++) {
		double tt;
		int ei = b->auxi[ee];
		if (c->p[0][ei]   >= (b->av[ei] + b->idist)
		 || c->p[ixc][ei] <= (b->av[ei] - b->idist)) {
			DBG(("Doesn't contain solution that will be closer to auxiliary goal\n"));
			return 0;
		}
		tt = (c->p[0][ei] + c->p[ixc][ei]) - b->av[ei];
		sort += tt * tt;
	}
	c->sort = sort + 0.01 * ss;

	if (c->ix == b->pauxcell)
		c->sort = -1.0;			/* Put previous calls solution cell at top of sort list */

	DBG(("Cell is accepted\n"));
	return 1;
}

/* Re-check whether it's worth searching cell */
static int auxil_check(schbase *b, cell *c) {
	int ee, ixc = b->ixc;

	DBG(("Reverse auxiliary search, re-check cell\n"));

	/* Check if this cell could possible improve b->idist */
	/* and compute sort key as the distance to auxilliary target */
	for (ee = 0; ee < b->naux; ee++) {
		int ei = b->auxi[ee];
		if (c->p[0][ei]   >= (b->av[ei] + b->idist)
		 || c->p[ixc][ei] <= (b->av[ei] - b->idist)) {
			DBG(("Doesn't contain solution that will be closer to auxiliary goal\n"));
			return 0;
		}
	}
	DBG(("Cell is still ok\n"));
	return 1;
}

/* Compute a solution for a given simplex (if there is one) */
/* Return 1 if search should be aborted */
static int auxil_compute(schbase *b, simplex *x) {
	rspl *s     = b->s;
	cell *c     = x->c;
	int e, di   = s->di;
	int f, fdi  = s->fdi;
	int ixc = b->ixc;
	datai xp;		/* solution in simplex relative coord order */
	datai p;		/* absolute solution */
	double idist;	/* Auxiliary input distance */
	int wsrv;		/* Within simplex return value */

	DBG(("\nAuxil: computing possible solution\n"));

	/* Check that the target lies within the simplex bounding cube */
	for (f = 0; f < fdi; f++) {
		if (b->v[f] < *x->min[f] || b->v[f] > *x->max[f]) {
			DBG(("Simplex is rejected - bounding cube\n"));
			return 0;
		}
	}

	/* Check if this cell could possible improve b->idist */
	for (e = 0; e < b->naux; e++) {
		int ei = b->auxi[e];					/* pmin/max[] is indexed in input space */
		if (*x->pmin[ei] >= (b->av[ei] + b->idist)
		 || *x->pmax[ei] <= (b->av[ei] - b->idist)) {
			DBG(("Simplex doesn't contain solution that will be closer to auxiliary goal\n"));
			return 0;
		}
	}

//printf("~~ About to create svd decomp\n");
	/* Create the SVD or LU decomp needed to compute solution or locus */
	if (add_lu_svd(x)) {
		DBG(("SVD decomp failed, skip simplex\n"));
		return 0;
	}

//printf("~~ About to solve locus for aux target\n");
	/* Now solve for locus parameter that minimises */
	/* distance to auxliary target. */
	if ((wsrv = auxil_solve(b, x, xp)) == 0) {
		DBG(("Target auxiliary along locus is outside simplex,\n"));
		DBG(("or computation failed, skip simplex\n"));
		return 0;
	}

//printf("~~ About to convert solution to absolute space\n");
	/* Convert solution from simplex relative to absolute space */
	simplex_to_abs(x, p, xp);

//printf("~~ soln = %f %f %f %f\n",p[0],p[1],p[2],p[3]);
//printf("~~ About to compute auxil distance\n");
	/* Compute distance to auxiliary target */
	for (idist = 0.0, e = 0; e < b->naux; e++) {
		int ei = b->auxi[e];
		double tt = b->av[ei] - p[ei];
		idist += tt * tt;
	}
	idist = sqrt(idist);

	/* We want the smallest error from auxiliary target */
	if (b->nsoln != 0 && idist >= b->idist) {	/* Equal or worse auxiliary solution */
		DBG(("idist = %f, better solution has been found before\n",idist));
		return 0;
	}

	/* Solution is accepted */
	DBG(("######## Accepting new solution with idist %f\n",idist));
	for (e = 0; e < di; e++)
		b->cpp[0].p[e] = p[e];
	for (f = 0; f < fdi; f++)
		b->cpp[0].v[f] = b->v[f];	/* Assumed to be an exact solution */
	b->idist = idist;
	b->nsoln = 1;
	b->pauxcell = c->ix;
	if (wsrv == 2)					/* Is above (disabled) ink limit */
		b->iclip = 1;

	return 0;
}

/* ------------------------------------ */
/* Locus range search functions */

static int locus_setsort(schbase *b, cell *c) {
	int f, fdi  = b->s->fdi;
	int lxi = b->lxi;	/* Auxiliary we are finding min/max of */
	int ixc = b->ixc;
	int ee, p2di = (1<<b->s->di);
	double sort, ss;

	DBG(("Reverse locus evaluate and set sort key on cell\n"));

#ifdef DEBUG
	if (b->s->di <= fdi) {	/* Assert ~1 */
		error("rspl auxiliary locus interp called with di <= fdi");
	}
#endif /* DEBUG */

	/* Check that the target lies within the cell bounding sphere */
	for (ss = 0.0, f = 0; f < fdi; f++) {
		double tt = c->bcent[f] - b->v[f];
		ss += tt * tt;
	}
	if (ss > c->bradsq) {
		DBG(("Cell is rejected - bounding sphere\n"));
		return 0;
	}

	if (b->limiten != 0 && c->limmin > b->limitv) {
		DBG(("Cell is rejected - ink limit, min = %f, limit = %f\n",c->limmin,b->limitv));
		return 0;
	}

	/* Check if this cell could possible improve the locus min/max */
	if (b->asegs == 0) {	/* If we aren't find all segments of the locus */
		if (c->p[0][lxi] >= b->min && c->p[ixc][lxi] <= b->max ) {
			DBG(("Doesn't contain solution that will expand the locus\n"));
			return 0;
		}
	}

	/* Compute sort index from average of auxiliary values */
	sort = (c->p[0][b->lxi] + c->p[ixc][b->lxi]);
	
	c->sort = sort + 0.01 * ss;

	DBG(("Cell is accepted\n"));
	return 1;
}

/* Re-check whether it's worth searching simplexes */
static int locus_check(schbase *b, cell *c) {
	int lxi  = b->lxi;	/* Auxiliary we are finding min/max of */
	int ixc = b->ixc;

	DBG(("Reverse locus re-check\n"));

	/* Check if this cell could possible improve the locus min/max */
	if (b->asegs == 0) {	/* If we aren't find all segments of the locus */
		if (c->p[0][lxi] >= b->min && c->p[ixc][lxi] <= b->max ) {
			DBG(("Doesn't contain solution that will expand the locus\n"));
			return 0;
		}
	}

	DBG(("Cell is still ok\n"));
	return 1;
}

static int auxil_locus(schbase *b, simplex *x);

/* We expect to be given a sub-simplex with no DOF, to give an exact solution */
static int locus_compute(schbase *b, simplex *x) {
	rspl *s  = b->s;
	cell *c  = x->c;
	int di   = s->di;
	int f, fdi  = s->fdi;
	int lxi  = b->lxi;	/* Auxiliary we are finding min/max of */

	DBG(("\nLocus: computing possible solution\n"));

	/* Check that the target lies within the simplex bounding cube */
	for (f = 0; f < fdi; f++) {
		if (b->v[f] < *x->min[f] || b->v[f] > *x->max[f]) {
			DBG(("Simplex is rejected - bounding cube\n"));
			return 0;
		}
	}

	/* Check if simplex could possible improve the locus min/max */
	if (b->asegs == 0) {	/* If we aren't find all segments of the locus */
		if (*x->pmin[lxi] >= b->min && *x->pmax[lxi] <= b->max ) {
			DBG(("Simplex doesn't contain solution that will expand the locus\n"));
			return 0;
		}
	}

//printf("~~ About to create svd decomp\n");
	/* Create the SVD decomp needed to compute solution extreme points */
	if (add_lu_svd(x)) {
		DBG(("SVD decomp failed, skip simplex\n"));
		return 0;
	}

//printf("~~ About to solve locus for aux extremes\n");
	/* Now solve for locus parameter that are at the extremes */
	/* of the axiliary we are interested in. */
	if (!auxil_locus(b, x)) {
		DBG(("Target auxiliary is outside simplex,\n"));
		DBG(("or computation failed, skip simplex\n"));
		return 0;
	}

	return 0;
}

/* ------------------- */
/* Vector clipping search functions */
static int clipv_setsort(schbase *b, cell *c) {
	int f, fdi  = b->s->fdi;
	double ss, dp;

	DBG(("Reverse clipping search evaluate cell\n"));

//printf("~~sphere center = %f %f %f, radius %f\n",c->bcent[0],c->bcent[1],c->bcent[2],sqrt(c->bradsq));
	/* Check if the clipping line intersects the bounding sphere */
	/* First compute dot product cdir . (bcent - v) */
	/* == distance to center of sphere in direction of clip vector */
	for (dp = 0.0, f = 0; f < fdi; f++) {
		dp += b->ncdir[f] * (c->bcent[f] - b->v[f]);
	}

	if (b->limiten != 0 && c->limmin > b->limitv) {
		DBG(("Cell is rejected - ink limit, min = %f, limit = %f\n",c->limmin,b->limitv));
		return 0;
	}

//printf("~~ dot product = %f\n",dp);
	/* Now compute closest distance to sphere center */
	for (ss = 0.0, f = 0; f < fdi; f++) {
		double tt = b->v[f] + dp * b->ncdir[f] - c->bcent[f];
		ss += tt * tt;
	}

//printf("~~ distance to sphere center = %f\n",sqrt(ss));
	if (ss > c->bradsq) {
		DBG(("Cell is rejected - wrong direction or bounding sphere\n"));
		return 0;
	}
	c->sort = dp;		/* May be -ve if beyond clip target point ? */

	DBG(("Cell is accepted\n"));
	return 1;
}

/* Clipping check functions */
/* Note that we don't bother with this check in setsort(), */
/* because we assume that nothing will set a small cdist */
/* before the search commences (unlike auxil). */
/* Note that line search loop exits on finding any solution. */
static int clipv_check(schbase *b, cell *c) {

	DBG(("Reverse clipping re-check\n"));

	if (b->cdist < INF_DIST) {	/* If some clip solution has been found */
		int f, fdi = b->s->fdi;
		double dist, brad;
		/* Compute a conservative "best possible solution clip distance" */
		for (dist = 0.0, f = 0; f < fdi ; f++) {
			double tt = (c->bcent[f] - b->v[f]);
			dist += tt * tt;
		}
		dist = sqrt(dist); /* Target distance to bounding */

		if (dist >= (c->brad + b->cdist)) {	/* Equal or worse clip solution */
			DBG(("Cell best possible solution worse than current\n"));
			return 0;
		}
	}

	DBG(("Cell is still ok\n"));
	return 1;
}

static int vnearest_clip_solve(schbase *b, simplex *x, double *xp, double *xv, double *err);

/* Compute a clip solution */
static int clipv_compute(schbase *b, simplex *x) {
	rspl   *s  = b->s;
	int f, fdi = s->fdi;
	datai p;				/* Input space solution */
	datao v;				/* Output space solution */
	double err;				/* output error of solution */
	int wsrv;	/* Within simplex return value */

	DBG(("Clips: computing possible solution\n"));

	/* Compute a solution value */
	if ((wsrv = vnearest_clip_solve(b, x, p, v, &err)) == 0) {
		DBG(("Doesn't contain a solution\n"));
		return 0;
	}

	/* We want the smallest clip error */
	/* (Should we reject points in -ve vector direction ??) */
	if (err >= b->cdist) {	/* Equal or worse clip solution */
		DBG(("better solution has been found before\n"));
		return 0;
	}

	simplex_to_abs(x, b->cpp[0].p, p);	/* Convert to abs. space & copy */

	DBG(("######## Accepting new clipv solution with error %f\n",err));
#ifdef DEBUG
	if (b->limiten != 0) {
		DBG(("######## Ink value = %f, limit %f\n",b->limit(b->cntx, b->cpp[0].p), b->limitv));
	}
#endif

	/* Put solution in place */
	for (f = 0; f < fdi; f++)
		b->cpp[0].v[f] = v[f];
	b->cdist = err;
	b->nsoln = 1;
	if (wsrv == 2)					/* Is above (disabled) ink limit */
		b->iclip = 1;

	return 0;
}

/* ------------------- */
/* Nearest clipping search functions */
static int clipn_setsort(schbase *b, cell *c) {
	int f, fdi  = b->s->fdi;
	double ss;

	DBG(("Reverse nearest clipping search evaluate cell\n"));

	/* Compute a conservative "best possible solution clip distance" */
	for (ss = 0.0, f = 0; f < fdi ; f++) {
		double tt = (c->bcent[f] - b->v[f]);
		ss += tt * tt;
	}
	ss = sqrt(ss); /* Target distance to bounding sphere */
	ss -= c->brad;
	if (ss < 0.0)
		ss = 0.0;

	/* Check that the cell could possibly improve the solution */
	if (b->cdist < INF_DIST) {	/* If some clip solution has been found */
		if (ss >= b->cdist) {	/* Equal or worse clip solution */
			DBG(("Cell best possible solution worse than current\n"));
			return 0;
		}
	}

	if (b->limiten != 0 && c->limmin > b->limitv) {
		DBG(("Cell is rejected - ink limit, min = %f, limit = %f\n",c->limmin,b->limitv));
		return 0;
	}

	c->sort = ss;		/* May be -ve if beyond clip target point ? */

	DBG(("Cell is accepted\n"));
	return 1;
}

/* Clipping check functions */
static int clipn_check(schbase *b, cell *c) {

	DBG(("Reverse nearest clipping re-check\n"));

	if (b->cdist < INF_DIST) {	/* If some clip solution has been found */
		/* re-use sort value, best possible distance to solution */
		if (c->sort >= b->cdist) {	/* Equal or worse clip solution */
			DBG(("Cell best possible solution worse than current\n"));
			return 0;
		}
	}

	DBG(("Cell is still ok\n"));
	return 1;
}

static int nnearest_clip_solve(schbase *b, simplex *x, double *xp, double *xv, double *err);

/* Compute a clip solution */
static int clipn_compute(schbase *b, simplex *x) {
	rspl   *s  = b->s;
	int f, fdi = s->fdi;
	datai p;				/* Input space solution */
	datao v;				/* Output space solution */
	double err;				/* output error of solution */
	int wsrv;	/* Within simplex return value */

	DBG(("Clipn: computing possible solution\n"));

	/* Compute a solution value */
	if ((wsrv = nnearest_clip_solve(b, x, p, v, &err)) == 0) {
		DBG(("Doesn't contain a solution\n"));
		return 0;
	}

	/* We want the smallest clip error */
	if (err >= b->cdist) {	/* Equal or worse clip solution */
		DBG(("better solution has been found before\n"));
		return 0;
	}

	DBG(("######## Accepting new clipn solution with error %f\n",err));

	simplex_to_abs(x, b->cpp[0].p, p);	/* Convert to abs. space & copy */

	/* Put solution in place */
	for (f = 0; f < fdi; f++)
		b->cpp[0].v[f] = v[f];
	b->cdist = err;
	b->nsoln = 1;
	if (wsrv == 2)					/* Is above (disabled) ink limit */
		b->iclip = 1;

	return 0;
}

/* -------------------------------------------------------- */
/* Cell/simplex solver middle level code */

/* Find the point on this sub-simplexes solution locus that is */
/* closest to the target auxiliary values, and return it in xp[] */
/* Return zero if this point canot be calculated, */
/* or it lies outside the simplex. */
/* Return 1 normally, and 2 if the solution would be over the ink limit */
static int
auxil_solve(
schbase *b,
simplex *x,
double *xp		/* Return solution xp[sdi] */
) {
	rspl *s = b->s;
	cell *c  = x->c;
	int ee, e, di = s->di, sdi = x->sdi; 
	int f,    fdi = s->fdi, efdi = x->efdi; 
	int dof = sdi-efdi;			 /* Degree of freedom of simplex locus */
	int *icomb = x->psxi->icomb; /* abs -> simplex coordinate translation */
	double auxt[MXRI];			/* Simplex relative auxiliary targets */
	double bb[MXRI];
	int wsrv;	/* Within simplex return value */

	DBG(("axuil_solve called\n"));

	if (dof < 0)
		error("Error - auxil_solve got sdi < efdi (%d < %d) - don't know how to handle this",sdi, efdi);

	/* If there is no locus, compute an exact solution */
	if (dof == 0) {
		DBG(("axuil_solve dof = zero\n"));

		/* Init the RHS B[] vector (note sdi == efdi) */
		for (f = 0; f < efdi; f++) {
			xp[f] = b->v[f] - x->v[sdi][f];
		}

		/* Compute the solution (in simplex space) */
		lu_backsub(x->d_u, sdi, (int *)x->d_w, xp);

		if ((wsrv = within_simplex(x, xp)) != 0) {
			DBG(("Got solution\n"));
			return wsrv;				/* OK, got solution */
		}

		DBG(("No solution\n"));
		return 0;
	}

	/* There is a locus, so find solution nearest auxiliaries */

	/* Compute locus for target function values (if sdi > efdi) */
	if (add_locus(b, x)) {
		DBG(("Locus computation failed, skip simplex\n"));
		return 0;
	}

	/* Convert aux targets from absolute space to simplex relative */
	for (e = 0; e < di; e++) {	/* For abs coords */
		int ei = icomb[e];		/* Simplex coord */

		if (ei >= 0 &&  b->auxm[e] != 0) {
			auxt[ei] = (b->av[e] - c->p[0][e])/s->g.w[e];	/* Only sets those needed */
		}
	}

	if (dof == 1 && b->naux == 1) {		/* Special case, because it's common and easy! */
		int ei = icomb[b->auxi[0]];		/* Simplex relative auxiliary index */
		double tt;

		DBG(("axuil_solve dof = naux = 1\n"));
		if (ei < 0)
			return 0;					/* Not going to find solution */
		if ((tt = x->lo_l[ei][0]) == 0.0)
			return 0;
		tt = (auxt[ei] - x->lo_bd[ei])/tt;	/* Parameter solution for target auxiliary */

		/* Back substitute parameter */
		for (e = 0; e < sdi; e++) {
			xp[e] = x->lo_bd[e] + tt * x->lo_l[e][0];
		}
		if ((wsrv = within_simplex(x, xp)) != 0) {
			DBG(("Got solution\n"));
			return wsrv;				/* OK, got solution */
		}
		DBG(("No solution\n"));
		return 0;
	}

	/* Compute the locus decompositions needed (info #5) */
	if (add_auxil_lu_svd(b, x)) {	/* Will set x->naux */
		DBG(("LU/SVD decomp failed\n"));
		return 0;
	}

	/* Setup B[], equation RHS  */
	for (e = ee = 0; ee < b->naux; ee++) {
		int ei = icomb[b->auxi[ee]];		/* Simplex relative auxiliary index */
		if (ei >= 0)						/* Usable auxiliary on this sub simplex */ 
			bb[e++] = auxt[ei] - x->lo_bd[ei];
	}
	if (e != x->naux)	/* Assert */
		error("Internal error - auxil_solve got mismatching number of auxiliaries");

	if (x->naux == dof) {			/* Use LU decomp to solve */
		DBG(("axuil_solve using LU\n"));
		lu_backsub(x->ax_u, dof, (int *)x->ax_w, bb);

	} else if (x->naux > 0) {	/* Use SVD to solve least squares */
		DBG(("axuil_solve using SVD\n"));
		svdbacksub(x->ax_u, x->ax_w, x->ax_v, bb, bb, x->naux, dof);

	} else {	/* x->naux == 0 */
		DBG(("axuil_solve  naux = 0\n"));
		for (f = 0; f < dof; f++)
			bb[f] = 0.0;		/* Use base solution ?? */
	}

	/* Now back substitute the locus parameters */
	/* to calculate the solution point (in simplex space) */
	for (e = 0; e < sdi; e++) {
		double tt;
		for (tt = 0.0, f = 0; f < dof; f++) {
			tt += bb[f] * x->lo_l[e][f];
		}
		xp[e] = x->lo_bd[e] + tt;
	}

	if ((wsrv = within_simplex(x, xp)) != 0) {
		DBG(("Got solution\n"));
		return wsrv;				/* OK, got solution */
	}
	DBG(("No solution\n"));
	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Compute the min/max values for the current auxiliary of interest. */
/* Return zero if this point canot be calculated, */
/* or it lies outside the simplex. */
/* Return 1 normally, 2 if it would be outside the simplex if limting was enabled */
/* We expect to get a sub-simplex that will give an exact solution. */
static int
auxil_locus(
schbase *b,
simplex *x
) {
	rspl *s = b->s;
	cell *c  = x->c;
	int ee, e, di = s->di, sdi = x->sdi; 
	int f,    fdi = s->fdi, efdi = x->efdi; 
	double pp[MXRI];
	int wsrv;	/* Within simplex return value */

	DBG(("axuil_locus called\n"));

	if (sdi != efdi)
		error("Internal error - auxil_locus got sdi != efdi (%d < %d)",sdi, efdi);

	/* Init the RHS B[] vector (note sdi == efdi) */
	for (f = 0; f < efdi; f++) {
		pp[f] = b->v[f] - x->v[sdi][f];
	}

	/* Compute the solution (in simplex space) */
	lu_backsub(x->d_u, sdi, (int *)x->d_w, pp);

	/* Check that the solution is within the simplex */
	if ((wsrv = within_simplex(x, pp)) != 0) {
		double xval;
		int lxi = b->lxi;	/* Auxiliary we are finding min/max of (Abs space) */
		int xlxi = x->psxi->icomb[lxi];	/* Auxiliary we are finding min/max of (simplex space) */

		DBG(("Got locus solution within simplex\n"));

		/* Compute auxiliary value for this solution (absolute space) */
		xval = c->p[0][lxi];
		if (xlxi >= 0)				/* Simplex param value */
			xval += s->g.w[lxi] * pp[xlxi];
		else if (xlxi == -2)		/* 1 value */
			xval += s->g.w[lxi];
									/* Else 0 value */

		if (b->asegs != 0) {		/* Tracking auxiliary segments */
			if (b->axisln >= b->axislz) {	/* Need some more space in list */
				if (b->axislz == 0) {
					b->axislz = 10;
					if ((b->axisl = (axisec *)malloc(b->axislz * sizeof(axisec))) == NULL)
						error("rev: malloc failed - Auxiliary intersect list size %d",b->axislz);
				} else {
					b->axislz *= 2;
					if ((b->axisl = (axisec *)realloc(b->axisl, b->axislz * sizeof(axisec)))
					    == NULL)
						error("rev: realloc failed - Auxiliary intersect list size %d",b->axislz);
				}
			}
			b->axisl[b->axisln].xval = xval;
			b->axisl[b->axisln].nv = x->sdi + 1;
			for (f = 0; f <= x->sdi; f++) {
				b->axisl[b->axisln].vix[f] = c->ix + s->g.hi[x->psxi->offs[f]];
			}
			b->axisln++;
		}

		/* If this solution is expands the min or max, save it */
		if (xval < b->min) {
			DBG(("######## Improving minimum to %f\n",xval));
			b->min = xval;
			b->plmincell = c->ix;
		}
		if (xval > b->max) {
			DBG(("######## Improving maximum to %f\n",xval));
			b->max = xval;
			b->plmaxcell = c->ix;
		}
	} else {
		DBG(("Solution wasn't within the simplex\n"));
		return 0;
	}

	return wsrv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Find the point on the clip line locus and simplexes */
/* valid surface, that is closest to the target output value. */
/* We expect to be given a sub simplex with sdi = fdi-1, and efdi = fdi */
/* or a limit sub-simplex with sdi = fdi, and efdi = fdi+1 */
/* Return zero if solution canot be calculated, */
/* return 1 normally, 2 if solution would be above the (disabled) ink limit */
static int
vnearest_clip_solve(
schbase *b,
simplex *x,
double *xp,		/* Return solution (simplex parameter space) */
double *xv,		/* Return solution (output space) */
double *err		/* Output error distance at solution point */
) {
	cell *c = x->c;
	rspl *s = b->s;
	int ee, e, di = s->di, sdi = x->sdi; 
	int f, fdi = s->fdi, efdi = x->efdi; 
	int i, g;
	int wsrv;	/* Within simplex return value */

	double *ta[MXRO], TA[MXRO][MXRO];
	double tb[MXRO];

	DBG(("Vector nearest clip solution called, cell %d, splx %d\n", c->ix, x->si));

	/* Setup temporary matricies */
	for (f = 0; f < sdi; f++) {
		ta[f] = TA[f];
	}

	/* Substitute simplex equation for output values V */
	/* in terms of sub-simplex parameters P, */
	/* into  clip line implicit equation in V, to give */
	/* clip line simplex implicit equation in terms of P (simplex input space) */
	/* If this is a limit sub-simlex, the ink limit part of the clip vector */
	/* equations will be used. */

	/* LHS: ta[sdi][sdi] = cla[sdi][efdi] * vv[efdi][sdi] */
	/* RHS: tb[sdi] = clb[sdi] - cla[sdi][efdi] * vv_di[efdi] */
	for (f = 0; f < sdi; f++) {
		double tt;
		for (e = 0; e < sdi; e++) {
			for (tt = 0.0, g = 0; g < efdi; g++)
				tt += b->cla[f][g] * (x->v[e][g] - x->v[e+1][g]);
			ta[f][e] = tt;
		}
		for (tt = 0.0, g = 0; g < efdi; g++)
			tt += b->cla[f][g] * x->v[sdi][g];
		tb[f] = b->clb[f] - tt;
	}

	/* Compute the solution */
	if (gen_solve_se(ta, tb, sdi, sdi)) {
		DBG(("Equation solution failed!\n"));
		return 0;		/* No solution */
	}

	/* Check that the solution is within the simplex */
	if ((wsrv = within_simplex(x, tb)) != 0) {
		double vv[MXRO];			/* space for trial solution */
		double dist;				/* distance to clip target */

		DBG(("Got solution within simplex\n"));

		/* Compute the output space solution point */
		for (f = 0; f < fdi; f++) {
			double tt = 0.0;
			for (e = 0; e < sdi; e++) {
				tt += (x->v[e][f] - x->v[e+1][f]) * tb[e];
			}
			xv[f] = tt + x->v[sdi][f];
		}

		/* Copy to return array */
		for (e = 0; e < sdi; e++)
			xp[e] = tb[e];

		/* Compute distance to clip target */
		for (dist = 0.0, f = 0; f < fdi ; f++) {
			double tt = (b->v[f] - xv[f]);
			dist += tt * tt;
		}
		DBGV(("Vector clip output soln: ",fdi," %f", xv, "\n"));

		/* Return the solution in xp[]m xv[] and *err */
		*err = sqrt(dist);

		DBG(("Vector clip returning a solution with error %f\n",*err));
		return wsrv;
	}

	DBG(("Vector clip solution not in simplex\n"));
	return 0;		/* No solution */
}

/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Find the point on the simplexes valid surface, that is closest */
/* to the target output value. */
/* We expect to be given a sub simplex with sdi = fdi-1, and efdi = fdi */
/* or a limit sub-simplex with sdi = fdi, and efdi = fdi+1 */
/* Return zero if solution canot be calculated, */
/* return 1 normally, 2 if solution would be above the (disabled) ink limit */
static int
nnearest_clip_solve(
schbase *b,
simplex *x,
double *xp,		/* Return solution (simplex parameter space) */
double *xv,		/* Return solution (output space) */
double *err		/* Output error distance at solution point */
) {
	cell *c = x->c;
	rspl *s = b->s;
	int ee, e, di = s->di, sdi = x->sdi; 
	int f, fdi = s->fdi, efdi = x->efdi; 
	int i, g;
	double tb[MXRO];		/* RHS & Parameter solution */
	double dist;			/* distance to clip target */
	int wsrv = 0;			/* Within simplex return value */

	DBG(("Nearest clip solution called, cell %d, splx %d\n", c->ix, x->si));

	if (sdi == 0) {		/* Solution is vertex */
		wsrv = 1;
		for (f = 0; f < efdi; f++)
			xv[f] = x->v[sdi][f]; 		/* Copy vertex value */
		if (x->v[sdi][fdi] > b->limitv) {
			if (b->limiten)				/* (is this needed ? simplex will have SPLX_OVLIMIT ?) */
				return 0;				/* Over ink limit - no good */
			wsrv = 2;					/* Would be over */
		}
		DBG(("Got assumed vertex solution\n"));
	} else {
#ifdef NEVER	/* Don't specialise ink limit version - use INKSCALE fudge instead */
		if (!(x->flags & SPLX_CLIPSX)) {	/* Not an ink limited plane simplex */
		
#endif
			/* Create the SVD decomp needed for least squares solution */
			if (add_lu_svd(x)) {
				DBG(("SVD decomp failed, skip simplex\n"));
				return 0;
			}
		
			/* Setup RHS to solve */
			for (f = 0; f < efdi; f++)
				tb[f] = b->v[f] - x->v[sdi][f]; 

			/* Find least squares solution */
			svdbacksub(x->d_u, x->d_w, x->d_v, tb, tb, efdi, sdi);
	
			/* Check that the solution is within the simplex */
			if ((wsrv = within_simplex(x, tb)) == 0) {
				DBG(("Nearest clip solution not in simplex\n"));
				return 0;		/* No solution */
			}
	
			DBG(("Got solution within simplex\n"));
	
			/* Compute the output space solution point */
			for (f = 0; f < fdi; f++) {
				double tt = 0.0;
				for (e = 0; e < sdi; e++) {
					tt += (x->v[e][f] - x->v[e+1][f]) * tb[e];
				}
				xv[f] = tt + x->v[sdi][f];
			}
#ifdef NEVER /* ~~1 Haven't figured out equations to make this a special case. */
			 /* Content to use INKSCALE fudge and rely on SVD least squares. */
		} else {
			/* We can't use the given equations, because we want the solution */
			/* to lie exactly on the ink limit plane, and be least squares to the */
			/* other target parameters. */
			/* Extract the ink limit parameters, and transform them into */
			/* a parameterised surface for this simplex. */
			/* Substitute the ink plane equation into the remaining target */
			/* parameter equations, and solve for least squares. */

		}
#endif
	}

	/* Copy to return array */
	for (e = 0; e < sdi; e++)
		xp[e] = tb[e];

	/* Compute distance to clip target */
	for (dist = 0.0, f = 0; f < fdi ; f++) {
		double tt = (b->v[f] - xv[f]);
		dist += tt * tt;
	}
	DBGV(("Nearest clip output soln: ",fdi," %f", xv, "\n"));

	/* Return the solution in xp[]m xv[] and *err */
	*err = sqrt(dist);

	DBG(("Nearest clip returning a solution with error %f\n",*err));
	return wsrv;
}


#ifdef NEVER
/* Utility to convert an implicit ink limit plane equation */
/* (held at the end of the simplex output value equations), */
/* into a parameterized surface equation. */
static void
compute_param_limit_surface(
schbase *b,
simplex *x
) {
	rspl *s = b->s;
	int ff, f, fdi = s->fdi;
	int i, p;
	double lgst;

double st[MXRO],	/* Start point */
double de[MXRO]		/* Delta */
	DBG(("Computing clipping line implicit equation, dim = %d\n", fdi));
	
	/* Pick a pivot element - the smallest */
	for (lgst = -1.0, p = -1, f = 0; f < fdi; f++) {
		double tt = de[f];
		b->cdir[f] = tt;		/* Stash this away */
		tt = fabs(tt);
		if (tt > lgst) {
			lgst = tt;
			p = f;
		}
	}
	if (p < 0)	/* Shouldn't happen */
		error("rspl rev, internal, trying to cope with zero length clip line\n");
	
	if (b->cla == NULL)
		b->cla = dmatrix(0, fdi-1, 0, fdi);	/* Allow for ink limit supliment */

	for (i = ff = 0;  ff < fdi; ff++) {	/* For the input rows */
		if (ff == p) {
			continue;					/* Skip pivot row */
		}
		for (f = 0; f < fdi; f++) {		/* For input & output columns */
			if (f == p) {
				b->cla[i][f] = -de[ff];	/* Last column is -ve delta value */
			} else if (f == ff) {
				b->cla[i][f] = de[p];	/* Diagonal is pivot value */
			} else {
				b->cla[i][f] = 0.0;		/* Else zero */
			}
		}
		b->clb[i] = de[p] * st[ff] - de[ff] * st[p];
		i++;
	}

	/* Add ink limit target equation - */
	/* interpolated ink value == target */
	if (b->limit != NULL) {
		for (i = 0;  i < (fdi-1); i++)
			b->cla[i][fdi] = 0.0;

		for (f = 0; f < fdi; f++) 
			b->cla[fdi-1][f] = 0.0;
		
		b->cla[fdi-1][fdi] = 1.0;
		b->clb[fdi-1] = b->limitv;
	}

#ifdef NEVER
/* Verify that the implicit equation is correct */
{
	double pnt[MXRO], v[MXRO];
	double pa;	/* Parameter */
	for (pa = 0.0; pa <= 1.0; pa += 0.125) {
		for (f = 0; f < fdi; f++) {
			pnt[f] = st[f] + pa * de[f];
		}

		/* Verify the implicit equation */
		for (ff = 0; ff < (fdi-1); ff++) {
			v[ff] = 0.0;
			for (f = 0; f < fdi; f++) {
				v[ff] += b->cla[ff][f] * pnt[f];
			}
			v[ff] -= b->clb[ff];
			if (v[ff] < 0.0)
				v[ff] = -v[ff];
			if (v[ff] > 0.000001) {
				printf("Point on clip line = %f %f %f\n",pnt[0],pnt[1],pnt[2]);
				printf("Implicit %d error of = %f\n",ff, v[ff]);
			}
		}
	}
}
#endif /* NEVER */

}

#endif




/* -------------------------------------------------------- */
/* Cell/simplex object lower level code */

/* Utility to get or calculate a vertexes ink limit value */
static double get_limitv(
schbase *b,			/* Base search information */
int ix,				/* fwd index of cell */
float *fcb,			/* Pointer to base of vertex value array (can be NULL) */
double *p			/* Array of input values (can be NULL to compute) */
) {
	rspl *s = b->s;
	float *base = fcb;
	double lv;
	if (base == NULL)
		base = b->s->g.a + ix * b->s->g.pss;
	lv = base[-1];					/* Fetch existing ink limit function value */
	if ((float)lv == L_UNINIT) {			/* Not been computed yet */
		if (p != NULL) {
			base[-1] = lv = INKSCALE * b->limit(b->cntx, p);	/* Do it */
		} else {
			int e, di = s->di;
			double pp[MXRI];			/* Copy from float to double */
			int tix;					/* Temp fwd cell index */

			for (tix = ix, e = 0; e < di; e++) {
				int dix;
				dix = tix % s->g.res[e];
				tix /= s->g.res[e];
				pp[e] = s->g.l[e] + (double)dix * s->g.w[e];	/* Base point */
			}
			base[-1] = lv = INKSCALE * b->limit(b->cntx, pp);	/* Do it */
		}
		s->g.limitv_cached = 1;			/* At least one limit value is cached */
	}
	return lv;
}

/* Utility to invalidate all the ink limit values */
/* cached in the main rspl array */
static void clear_limitv(
rspl *s
) {
	int i;
	float *gp;		/* Grid point pointer */

	if (s->g.limitv_cached != 0) {	/* If any have been set */
		/* Unset them all */
		for (i = 0, gp = s->g.a; i < s->g.no; i++, gp += s->g.pss) {
			gp[-1] = L_UNINIT;
		}
		s->g.limitv_cached = 0;
	}
}

/* Cell code */

static void free_cell_info(cell *c);
static cell *cache_rcell(revcache *r, int ix);
static void uncache_rcell(revcache *r, cell *cp);

/* Return a pointer to an appropriate reverse cell */
/* cache structure. None of the sub simplex lists will */
/* be initialised. */
/* NOTE: must unget_cell() (== uncache_rcell()) when cell */
/* is no longer needed */
/* Return NULL if we ran out of room in the cache */
static cell *get_rcell(
schbase *b,			/* Base search information */
int ix				/* fwd index of cell */
) {
	rspl *s = b->s;
	int ee, e, di = s->di;
	int p2di = (1<<di);
	int ff, f, fdi = s->fdi;
	int dof;
	cell *c;

	c = cache_rcell(s->rev.cache, ix);		/* Fetch it from the cache and lock it */
	if (c == NULL)
		return NULL;

	if (!(c->flags & CELL_FLAG_1)) {			/* Have to (re)initialize cell & simplexes */
		int tix;								/* Temp fwd cell index */
		float *fcb = s->g.a + ix * s->g.pss;	/* Pointer to base float of fwd cell */

		/* Free any existing simplexes, since the number or use migh change */
		free_cell_info(c);

		/* Setup lower edge flags */
		c->lwredge = 0;
		for (e = 0; e < di; e++) {		/* Transfer cell verticy values from grid */
			if (G_FL(fcb, e) == 4)		/* Cell is on lower edge of grid */
				c->lwredge |= (1 << e);	/* Set flag */
		}

		/* Compute basic Cell info and vertex output values */
		for (ee = 0; ee < p2di; ee++) {
			float *vp = fcb + s->g.fhi[ee];
			for (f = 0; f < fdi; f++)		/* Transfer cell verticy values from grid */
				c->v[ee][f] = vp[f];

			/* ~~ reset any other cell info that will be stale */
		}

		/* Convert from cell index, to absolute fwd coord base values */
		c->limmin = INF_DIST;				/* and min/max values */
		c->limmax = -INF_DIST;
		for (tix = ix, e = 0; e < di; e++) {
			int dix;
			dix = tix % s->g.res[e];
			tix /= s->g.res[e];
			c->p[0][e] = s->g.l[e] + (double)dix * s->g.w[e];	/* Base point */
		}
		if (b->limit != NULL) {			/* Compute ink limit values at base verticy */
			double lv = get_limitv(b, ix, fcb, c->p[0]); /* Fetch or generate limit value */
			c->v[0][fdi] = lv;
			if (lv < c->limmin)	/* And min/max for this cell */
				c->limmin = lv;
			if (lv > c->limmax)
				c->limmax = lv;
		}
			
		/* Setup cube verticy input position values, and ink limit values */
		for (ee = 1; ee < p2di; ee++) {
			for (e = 0; e < di; e++) {
				c->p[ee][e] = c->p[0][e];
				if (ee & (1 << e))
					c->p[ee][e] += s->g.w[e];		/* In input space offset */
			}
			if (b->limit != NULL) {			/* Compute ink limit values at cell verticies */
				double lv = get_limitv(b, ix, fcb + s->g.fhi[ee], c->p[ee]);
				c->v[ee][fdi] = lv;
				if (lv < c->limmin)	/* And min/max for this cell */
					c->limmin = lv;
				if (lv > c->limmax)
					c->limmax = lv;
			}
		}
		
		/* Compute the output bounding sphere for fast rejection testing */
		{
			double *min[MXRO], *max[MXRO];	/* Pointers to points with min/max values */
			double radsq = -1.0;			/* Span/radius squared */
			double rad;
			int spf;
			
			/* Find verticies of cell that have min and max values in output space */
			for (f = 0; f < fdi; f++)
				min[f] = max[f] = NULL;

			for (ee = 0; ee < p2di; ee++) {
				double *vp = c->v[ee];
				for (f = 0; f < fdi; f++) {
					if (min[f] == NULL || min[f][f] > vp[f])
						min[f] = vp;
					if (max[f] == NULL || max[f][f] < vp[f])
						max[f] = vp;
				}
			}

			/* Find the pair of points with the largest span (diameter) in output space */
			for (ff = 0; ff < fdi; ff++) {
				double ss;
				for (ss = 0.0, f = 0; f < fdi; f++) {
					double tt;
					tt = max[ff][f] - min[ff][f];
					ss += tt * tt;
				}
				if (ss > radsq) {
					radsq = ss;
					spf = ff;		/* Output dimension max was in */
				}
			}

			/* Set initial bounding sphere */
			for (f = 0; f < fdi; f++) {
				c->bcent[f] = (max[spf][f] + min[spf][f])/2.0;
			}
			radsq /= 4.0;			/* diam^2 -> rad^2 */
			c->bradsq = radsq;
			rad = c->brad = sqrt(radsq);
			
			/* Go though all the points again, expanding sphere if necessary */
			for (ee = 0; ee < p2di; ee++) {
				double ss;
				double *vp = c->v[ee];

				/* Compute distance squared of point to bounding shere */
				for (ss = 0.0, f = 0; f < fdi; f++) {
					double tt = vp[f] - c->bcent[f];
					ss += tt * tt;
				}
				if (ss > radsq) {
					double tt;
					/* DBG(("Expanding bounding sphere by %f\n",sqrt(ss) - rad)); */

					ss = sqrt(ss) + EPS;			/* Radius to point */
					rad = (rad + ss)/2.0;
					c->bradsq = radsq = rad * rad;
					tt = ss - rad;
					for (f = 0; f < fdi; f++) {
						c->bcent[f] = (rad * c->bcent[f] + tt * vp[f])/ss;
					}

				} else {
					/* DBG(("Bounding sphere encloses by %f\n",rad - sqrt(ss))); */
				}
			}
		}
		c->flags = CELL_FLAG_1;
	}

	return c;
}

void free_simplex_info(cell *c, int dof);

/* Free up any allocated information in a cell */
static void
free_cell_info(
cell *c
) {
	int nsdi;
	
	/* Free up all the simplexes */
	for (nsdi = 0; nsdi <= c->s->di; nsdi++) {
		if (c->sx[nsdi] != NULL) {
			free_simplex_info(c, nsdi);
			c->sx[nsdi] = NULL;
		}
	}
	/* ~~ free any other cell information */
}

/* - - - - - -  */
/* Simplex code */

/* Allocate and do the basic initialisation for a DOF list of simplexes */
void alloc_simplexes(
cell *c,
int nsdi			/* Non limited sub simplex dimensionality */
) {
	rspl *s = c->s;
	schbase *b = s->rev.sb;
	int ee, e, di = s->di;
	int f, fdi = s->fdi;
	int lsdi;			/* Ink limited Sub-simplex sdi */
	int tsxno;			/* Total number of DOF simplexes */
	int nsxno;			/* Number of non-ink limited DOF simplexes */
	int si, so;			/* simplex index in and out */

	DBG(("Allocating level %d sub simplexes in cell %d\n",nsdi,c->ix));
	if (c->sx[nsdi] != NULL)
		error("rspl rev, internal, trying allocate already allocated simplexes\n");

	/* Figure out how many simplexes will be at this nsdi */
	lsdi = nsdi + 1;	/* Limit simplexes sdi */

	tsxno = nsxno = s->rev.sspxi[nsdi].nospx;
	if (b->limit != NULL && lsdi <= di)
		tsxno += s->rev.sspxi[lsdi].nospx;

	/* Make sure there is enough space in temp simplex filter list */
	if (b->lsxfilt < tsxno) {	/* Allocate more space if needed */

		if (b->lsxfilt > 0) {	/* Free old space before allocating new */
			free (b->sxfilt);
		}
		b->lsxfilt = 0;
		/* Allocate enough space for all the candidate cells */
		if ((b->sxfilt = (char *)malloc(tsxno * sizeof(char))) == NULL)
			error("rev: malloc failed - temp simplex filter list, count %d",tsxno);
		b->lsxfilt = tsxno;	/* Current allocated space */
	}
		
	/* Figure out the number of simplexes that will actually be needed */
	for (si = so = 0; si < tsxno; si++) {
		psxinfo *psxi;
		int *icomb, *offs;
		int sdi = nsdi;
		int efdi = fdi;
		int ssi = si;
		int isclip = 0;
		if (si >= nsxno) {				/* If limit boundary simplex */
			sdi++;						/* One more dimension */
			efdi++;						/* One more constraint */
			ssi -= nsxno;				/* In second half of list */
			isclip++;					/* Limit clipped simplex */
		}
		psxi = &s->rev.sspxi[sdi].spxi[ssi];
		icomb = psxi->icomb;
		offs  = psxi->offs;

		b->sxfilt[si] = 0;				/* Assume simplex won't be used */

		/* See if simplex can be ignored because it is accounted for in another cell */
		for (e = 0; e < di; e++) {
			if (icomb[e] == -1 && !(c->lwredge & (1 << e))) {
				break;	/* Simplex is on lower edge of cell,  */
						/* while cell is not on lower edge of grid */
			}
		}
		if (e < di) {
			continue;	/* Skip this simplex */
		}

		/* Check if ink limit is relevant */
		if (b->limit != NULL) {
			double max = -INF_DIST;
			double min =  INF_DIST;

			for (e = 0; e <= sdi; e++) {		/* For all the simplex verticies */
				int i = offs[e];
				double vv = c->v[i][fdi];		/* Ink limit value */
				if (vv < min)
					min = vv;
				if (vv > max)
					max = vv;
			}

			
			if (isclip) {
				if (max < b->limitv || min > b->limitv) {
					continue;	/* Discard this simplex - it can't straddle the ink limit */
				}
			} else {
				if (min > b->limitv) {
					continue;				/* Discard this simplex - it is above the ink limit */
				}
			}
		}

		b->sxfilt[si] |= 1;		/* This cell will be OK */
		so++;
	}

	DBG(("There are %d level %d sub simplexes\n",so, nsdi));
	/* Allocate space for all the DOF simplexes that will be used */
	if (so > 0 && (c->sx[nsdi] = (simplex *) malloc(so * sizeof(simplex))) == NULL)
		error("rspl malloc failed - reverse cell simplexes");

	/* Setup SPLX_FLAG_1 level information in the simplex */
	for (si = so = 0; si < tsxno; si++) {
		simplex *x;
		psxinfo *psxi;
		int *icomb;
		int sdi, efdi;
		int ssi;

		if (b->sxfilt[si] == 0)		/* Don't use this one */
			continue;

#ifdef STATS
		s->rev.st[b->op].sinited++;
#endif /* STATS */

		x = &c->sx[nsdi][so];
		x->flags = 0;

		sdi = nsdi;
		efdi = fdi;
		ssi = si;
		if (si >= nsxno) {				/* If limit boundary simplex */
			sdi++;						/* One more dimension */
			efdi++;						/* One more constraint */
			ssi -= nsxno;				/* In second half of list */
			x->flags |= SPLX_CLIPSX;	/* Limit clipped simplex */
		}

		if (b->sxfilt[si] == 3)		/* Don't use this one */
			x->flags |= SPLX_OVLIMIT;	/* This one is over the ink limit */

		psxi = &s->rev.sspxi[sdi].spxi[ssi];
		icomb = psxi->icomb;

		/* Fill in the other simplex details */
		x->psxi = psxi;					/* Pointer to constant per simplex info */
		x->c    = c;					/* Parent cell */
		x->si   = so;					/* Diagnostic, simplex offset in list */
		x->sdi  = sdi;					/* Copy of simplex dimensionaity */
		x->efdi = efdi;					/* Copy of effective output dimensionality */

		/* Compute pointers to simplex vertex output and limit values */

		for (e = 0; e <= sdi; e++) {		/* For all the simplex verticies */
			int i = x->psxi->offs[e];
			x->v[e] = c->v[i];

			/* Setup output bounding box value pointers (the hard way) */
			if (e == 0) {					/* Init to first vertex of simplex */
				for (f = 0; f <= fdi; f++)				/* Output space */
					x->min[f] = x->max[f] = &c->v[i][f];
			} else {
				for (f = 0; f <= fdi; f++) {			/* Output space + ink sum */
					double vv, *vp;
					if (f == fdi && b->limit == NULL)
						continue;			/* Skip ink */
					vp = &c->v[i][f];
					vv = *vp;
					if (vv < *x->min[f])
						x->min[f] = vp;
					else if (vv > *x->max[f])
						x->max[f] = vp;
				}
			}
		}

		/* Setup input bounding box value pointers (the easy way) */
		for (ee = 0; ee < di; ee++) {
			x->pmin[ee] = &c->p[x->psxi->pmino[ee]][ee];
			x->pmax[ee] = &c->p[x->psxi->pmaxo[ee]][ee];
		}

		x->flags |= SPLX_FLAG_1;		/* vv & iv done, nothing else */

		x->aloc2 = x->aloc5 = NULL;		/* Matrix allocations not done yet */
		so++;
	}
	c->sxno[nsdi] = so;				/* Record actual number in list */
	c->flags |= CELL_FLAG_2;		/* Note that cell now has simplexes */
}

/* Free up any allocated for a list of sub-simplexes */
void
free_simplex_info(
cell *c,
int nsdi			/* non limit sub simplex dimensionaity */
) {
	int si, sxno = c->sxno[nsdi];	/* Number of simplexes */

	for (si = 0; si < sxno; si++) { /* For all the simplexes */
		simplex *x = &c->sx[nsdi][si];
		int sdi  = x->sdi;
		int efdi = x->efdi;
		int dof  = sdi - efdi;	/* Extra DOF */
		int naux = x->naux;

		if (x->aloc2 != NULL)
			free(x->aloc2);

		if (x->aloc5 != NULL)
			free(x->aloc5);

		/* ~~ free any other simplex information */
	}
	free(c->sx[nsdi]);
	c->sx[nsdi] = NULL;

	/* ~~ free any other cell information */
}

/* - - - - - - - - - - - - */
/* Check that an input space vector is within a given simplex, */
/* and that it meets any ink limit. */
/* Return zero if outside the simplex, */
/* 1 normally if within the simplex, */
/* and 2 if it would be over the ink limit if limit was enabled. */
static int
within_simplex(
simplex *x,				/* Simplex */
double *p				/* Input coords in simplex space */
) {
	cell *c = x->c;
	rspl *s = c->s;
	schbase *b = s->rev.sb;
	int    fdi = s->fdi;
	int e, sdi = x->sdi;		/* simplex dimensionality */
	double cp, lp;
	int rv = 1;
	/* EPS is allowance for numeric error */
	/* (Don't want solutions falling down */
	/* the numerical cracks between the simplexes) */

	/* Check we are within baricentric limits */
	for (lp = 0.0, e = 0; e < sdi; e++) {
		cp = p[e];
		if ((cp+EPS) < lp)		/* Outside baricentric or not in correct */
			return 0;			/* order for this simplex  */
		lp = cp;
	}
	if ((1.0+EPS) < lp)			/* outside baricentric range */
		return 0;

	/* Compute limit using interp. - assume simplex would have been trivially rejected */
	if (b->limit != NULL) {
		double sum = 0.0;			/* Might be over the limit */
		for (e = 0; e < sdi; e++)
			sum += p[e] * (x->v[e][fdi] - x->v[e+1][fdi]);
		sum += x->v[sdi][fdi];
		if (sum > b->limitv) {
			if (b->limiten != 0) {
	 			return 0;			/* Exceeds ink limit */
			} else {
				rv = 2;				/* would have exceeded limit */
			}
		}
	}

//~~1 is this needed ?????
	/* Constrain to legal values */
	for (e = 0; e < sdi; e++) {
		cp = p[e];
		if (cp < 0.0)
			p[e] = 0.0;
		else if (cp > 1.0)
			p[e] = 1.0;
	}
	return rv;
}

/* Convert vector from simplex space to absolute space */
static void simplex_to_abs(
simplex *x,
double *out,	/* output in absolute space */
double *in		/* Input in simplex space */
) {
	cell *c     = x->c;
	rspl *s     = c->s;
	int e, di   = s->di;
	int *icomb  = x->psxi->icomb;	/* Coord combination order */

	for (e = 0; e < di; e++) {		/* For each absolute coord */
		double ov = c->p[0][e];		/* Base value */
		int ee = icomb[e];			/* Simplex param index */
		if (ee >= 0)				/* Simplex param value */
			ov += s->g.w[e] * in[ee];
		else if (ee == -2)			/* 1 value */
			ov += s->g.w[e];
									/* Else 0 value */
		out[e] = ov;
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Given the parametric clip line equation, compute the */
/* implicit equation in terms of the absolute output space. */
/* Pad equation with target ink limit in case it is use */
/* with CLIPSX sub-simplexes. */
/* Note that no line equation values are returned if fdi = 1, */
/* since there is no such thing as an implicit line equation. */
static void
init_line_eq(
schbase *b,
double st[MXRO],	/* Start point */
double de[MXRO]		/* Delta */
) {
	rspl *s = b->s;
	int ff, f, fdi = s->fdi;
	int i, p;
	double lgst;

	DBG(("Computing clipping line implicit equation, dim = %d\n", fdi));
	
	/* Pick a pivot element */
	for (lgst = -1.0, p = -1, f = 0; f < fdi; f++) {
		double tt = de[f];
		b->cdir[f] = tt;		/* Stash this away */
		tt = fabs(tt);
		if (tt > lgst) {
			lgst = tt;
			p = f;
		}
	}
	if (p < 0)	/* Shouldn't happen */
		error("rspl rev, internal, trying to cope with zero length clip line\n");
	
	if (b->cla == NULL)
		b->cla = dmatrix(0, fdi-1, 0, fdi);	/* Allow for ink limit supliment */

	for (i = ff = 0;  ff < fdi; ff++) {	/* For the input rows */
		if (ff == p) {
			continue;					/* Skip pivot row */
		}
		for (f = 0; f < fdi; f++) {		/* For input & output columns */
			if (f == p) {
				b->cla[i][f] = -de[ff];	/* Last column is -ve delta value */
			} else if (f == ff) {
				b->cla[i][f] = de[p];	/* Diagonal is pivot value */
			} else {
				b->cla[i][f] = 0.0;		/* Else zero */
			}
		}
		b->clb[i] = de[p] * st[ff] - de[ff] * st[p];
		i++;
	}

	/* Add ink limit target equation - */
	/* interpolated ink value == target */
	if (b->limit != NULL) {
		for (i = 0;  i < (fdi-1); i++)
			b->cla[i][fdi] = 0.0;

		for (f = 0; f < fdi; f++) 
			b->cla[fdi-1][f] = 0.0;
		
		b->cla[fdi-1][fdi] = 1.0;
		b->clb[fdi-1] = b->limitv;
	}

#ifdef NEVER
/* Verify that the implicit equation is correct */
{
	double pnt[MXRO], v[MXRO];
	double pa;	/* Parameter */
	for (pa = 0.0; pa <= 1.0; pa += 0.125) {
		for (f = 0; f < fdi; f++) {
			pnt[f] = st[f] + pa * de[f];
		}

		/* Verify the implicit equation */
		for (ff = 0; ff < (fdi-1); ff++) {
			v[ff] = 0.0;
			for (f = 0; f < fdi; f++) {
				v[ff] += b->cla[ff][f] * pnt[f];
			}
			v[ff] -= b->clb[ff];
			if (v[ff] < 0.0)
				v[ff] = -v[ff];
			if (v[ff] > 0.000001) {
				printf("Point on clip line = %f %f %f\n",pnt[0],pnt[1],pnt[2]);
				printf("Implicit %d error of = %f\n",ff, v[ff]);
			}
		}
	}
}
#endif /* NEVER */

}

/* - - - - - -  */
/* Simpex solution info #2 */

/* Create the LU or SVD decomp needed to compute solution or locus. */
/* Return non-zero if it cannot be created */
static int
add_lu_svd(simplex *x) {

	if (x->flags & SPLX_FLAG_2F) {		/* Previously failed */
		return 1;
	}
	if (!(x->flags & SPLX_FLAG_2)) {
		int ee, e, sdi = x->sdi; 
		int f, efdi = x->efdi; 
		int dof = sdi-efdi;		/* Degree of freedom of locus, or -ve over specification */
		int adof = dof >= 0 ? dof : 0;		/* Allocation dof */
		int i;

		if (x->aloc2 == NULL) {	/* Allocate space for matricies and arrays */
			/* Do this in one hit to minimise malloc overhead */
			if (dof == 0) {
				int i;
				char *mem;
				int asize = sizeof(double) * (efdi * sdi)
				          + sizeof(double *) * efdi 
				          + sizeof(int) * sdi;

				if ((x->aloc2 = mem = (char *) malloc(asize)) == NULL)
					error("rspl malloc failed - reverse cell sub-simplex matricies");

				/* Allocate biggest to smallest (double, pointers, ints) */
				/* to make sure that items lie on the natural boundaries. /

				/* Reserve matrix doubles */
				mem += efdi * sdi * sizeof(double);

				/* Allocate pointers */
				x->d_u = (double **)mem, mem += efdi * sizeof(double *);

				/* Allocate ints */
				x->d_w = (double *)mem, mem += sdi * sizeof(int);

#ifdef DEBUG
				if (mem != (x->aloc2 + asize))
					error("~1 aloc2a assert failed! Is %d, should be %d\n",mem - x->aloc2,asize);
#endif /* DEBUG */

				/* Reset and allocate matrix doubles */
				mem = x->aloc2; 
				for (i = 0; i < efdi; i++)
					x->d_u[i] = (double *)mem,	mem += sdi * sizeof(double);

			} else {
				int i;
				char *mem;
				int asize = sizeof(double) * (sdi * (efdi + sdi + adof + 2) + efdi)
				          + sizeof(double *) * (efdi + 2 * sdi);

				if ((x->aloc2 = mem = (char *) malloc(asize)) == NULL)
					error("rspl malloc failed - reverse cell sub-simplex matricies");

				/* Allocate biggest to smallest (double, pointers, ints) */
				/* to make sure that items lie on the natural boundaries. /

				/* Reserve matrix doubles */
				mem += sdi * (efdi + sdi + adof) * sizeof(double);

				/* Allocate doubles */
				x->lo_xb = (double *)mem, mem += efdi * sizeof(double);
				x->lo_bd = (double *)mem; mem += sdi * sizeof(double);
				x->d_w = (double *)mem, mem += sdi * sizeof(double);

				/* Allocate pointers */
				x->d_u = (double **)mem, mem += efdi * sizeof(double *);
				x->d_v = (double **)mem, mem += sdi * sizeof(double *);
				x->lo_l = (double **)mem, mem += sdi * sizeof(double *);

#ifdef DEBUG
				if (mem != (x->aloc2 + asize))
					error("~1 aloc2b assert failed! Is %d, should be %d\n",mem - x->aloc2,asize);
#endif /* DEBUG */

				/* Reset and allocate matrix doubles */
				mem = x->aloc2;
				for (i = 0; i < efdi; i++)
					x->d_u[i] = (double *)mem,	mem += sdi * sizeof(double);
				for (i = 0; i < sdi; i++)
					x->d_v[i] = (double *)mem,	mem += sdi * sizeof(double);
				for (i = 0; i < sdi; i++)
					x->lo_l[i] = (double *)mem,	mem += adof * sizeof(double);

				/* Init any values that will be read before being written to. */
				for (f = 0; f < efdi; f++)
					x->lo_xb[f] = 1e100;		/* Silly value */
			}
		}

		/* Setup matrix from vertex values */
		for (f = 0; f < efdi; f++)
			for (e = 0; e < sdi; e++)
				x->d_u[f][e] = x->v[e][f] - x->v[e+1][f];

		if (dof == 0) {	/* compute LU */
			double rip;
#ifdef STATS
			x->c->s->rev.st[x->c->s->rev.sb->op].sinited2a++;
#endif /* STATS */
			if (lu_decomp(x->d_u, sdi, (int *)x->d_w, &rip)) {
				x->flags |= SPLX_FLAG_2F;	/* Failed */
				return 1;
			}
		} else {
//printf("~~ Creating SVD decomp, sdi = %d, efdi = %d\n", sdi, efdi);

#ifdef STATS
			x->c->s->rev.st[x->c->s->rev.sb->op].sinited2b++;
#endif /* STATS */
			if (svdecomp(x->d_u, x->d_w, x->d_v, efdi, sdi)) {
				x->flags |= SPLX_FLAG_2F;	/* Failed */
				return 1;
			}
	
			/* Threshold the singular values W[] */ 
			svdthresh(x->d_w, sdi);
	
			if (dof >= 0) {		/* If we expect a locus */
//printf("~~ got dif %d locus from SVD\n",dof);
				/* copy the locus direction coefficients out */
				for (i = e = 0; e < sdi; e++) {
					if (x->d_w[e] == 0.0) {		/* Found a zero W[] */
						if (i < dof) {
							for (ee = 0; ee < sdi; ee++) {	/* Copy column of V[][] */
								x->lo_l[ee][i] = x->d_v[ee][e];
							}
						}
						i++;
					}
				}
				if (i != dof) {
//printf("~~ got unexpected dof in svd\n");
					x->flags |= SPLX_FLAG_2F;	/* Failed */
					return 1;					/* Didn't get expected d.o.f. */
				}
			}
		}
		x->flags |= SPLX_FLAG_2;	/* Set flag so that it isn't attempted again */
	}
	return 0;
}

/* - - - - - -  */
/* Simplex solution info #4 */

/* Calculate the solution locus equation for this simplex and target */
/* (The direction was calculated by add_svd(), but now calculate */
/* the base solution point for this particular reverse lookup) */
/* Return non-zero if this point canot be calculated */
/* We are assuming that sdi > efdi */
static int
add_locus(
schbase *b,
simplex *x
) {
	int e,  sdi = x->sdi; 
	int f, efdi = x->efdi; 
	int doback = 0;

#ifdef STATS
	x->c->s->rev.st[x->c->s->rev.sb->op].sinited4++;
#endif /* STATS */
	/* Use output of svdcmp() to solve overspecified and/or */
	/* singular equation A.x = b */

	/* Init the RHS B[] vector, and check if it doesn't match */
	/* that used to compute base value last time. */
	for (f = 0; f < efdi; f++) {
		double xb = b->v[f] - x->v[sdi][f];
		if (x->lo_xb[f] != xb) {
			x->lo_xb[f] = xb;
			doback = 1;			/* RHS differs, so re-compute */
		}
	}
	
#ifdef STATS
	if (doback && (x->flags & SPLX_FLAG_4))
		x->c->s->rev.st[x->c->s->rev.sb->op].sinited4i++;
#endif /* STATS */

	/* Compute locus */
	if (doback || !(x->flags & SPLX_FLAG_4))
		svdbacksub(x->d_u, x->d_w, x->d_v, x->lo_xb, x->lo_bd, efdi, sdi);
	
	x->flags |= SPLX_FLAG_4;

	return 0;
}

/* - - - - - -  */
/* Simplex solution info #5 */

/* Compute LU or SVD decomp of lo_l */
/* Return non-zero if this canot be calculated. */
static int
add_auxil_lu_svd(
schbase *b,
simplex *x
) {
	int ee, e, sdi = x->sdi; 
	int f, efdi = x->efdi; 
	int dof = sdi-efdi;		/* Degree of freedom of locus */
	int naux = b->naux;		/* Number of auxiliaries actually available */

#ifdef STATS
	if (x->aaux != b->naux || x->auxbm != b->auxbm)
		x->c->s->rev.st[x->c->s->rev.sb->op].sinited5i++;
#endif /* STATS */

	if (x->aaux != b->naux) {	/* Number of auxiliaries has changed */
		if (x->aloc5 != NULL) {
			free(x->aloc5);
			x->aloc5 = NULL;
		}
		x->flags &= ~(SPLX_FLAG_5 | SPLX_FLAG_5F);	/* Force recompute */
	}
	
	if (x->auxbm != b->auxbm) {	/* Different selection of auxiliaries */
		x->flags &= ~(SPLX_FLAG_5 | SPLX_FLAG_5F);	/* Force recompute */
	}

	if (x->flags & SPLX_FLAG_5F) {		/* Previously failed */
		return 1;
	}
	if (!(x->flags & SPLX_FLAG_5)) {
		int *icomb = x->psxi->icomb; /* abs -> simplex coordinate translation */

		if (x->aloc5 == NULL) {	/* Allocate space for matricies and arrays */
			/* Do this in one hit to minimise malloc overhead */
			if (naux == dof) {
				int i;
				char *mem;
				int asize = sizeof(double *) * naux
				          + sizeof(double) * (naux * dof)
				          + sizeof(int) * dof;

				if ((x->aloc5 = mem = (char *) malloc(asize)) == NULL)
					error("rspl malloc failed - reverse cell sub-simplex matricies");

				/* Allocate biggest to smallest (double, pointers, ints) */
				/* to make sure that items lie on the natural boundaries. /

				/* Reserve matrix doubles */
				mem += naux * dof * sizeof(double);

				/* Allocate pointers and ints */
				x->d_u = (double **)mem, mem += naux * sizeof(double *);
				x->d_w = (double *)mem, mem += dof * sizeof(int);

#ifdef DEBUG
				if (mem != (x->aloc5 + asize))
					error("aloc5a assert failed! Is %d, should be %d\n",mem - x->aloc5,asize);
#endif /* DEBUG */

				/* Reset and allocate matrix doubles */
				mem = x->aloc5;
				for (i = 0; i < naux; i++)
					x->d_u[i] = (double *)mem,	mem += dof * sizeof(double);
			} else {
				int i;
				char *mem;
				int asize = sizeof(double *) * (naux + dof) 
				          + sizeof(double) * (dof * (naux + dof + 1));

				if ((x->aloc5 = mem = (char *) malloc(asize)) == NULL)
					error("rspl malloc failed - reverse cell sub-simplex matricies");

				/* Allocate biggest to smallest (double, pointers, ints) */
				/* to make sure that items lie on the natural boundaries. /

				/* Reserve matrix doubles */
				mem += dof * (naux + dof) * sizeof(double);

				/* Allocate doubles */
				x->ax_w = (double *)mem, mem += dof * sizeof(double);

				/* Allocate pointers, ints */
				x->ax_u = (double **)mem, mem += naux * sizeof(double *);
				x->ax_v = (double **)mem, mem += dof * sizeof(double *);

#ifdef DEBUG
				if (mem != (x->aloc5 + asize))
					error("aloc5b assert failed! Is %d, should be %d\n",mem - x->aloc5,asize);
#endif /* DEBUG */

				/* Reset and allocate matrix doubles */
				mem = x->aloc5;
				for (i = 0; i < naux; i++)
					x->ax_u[i] = (double *)mem,	mem += dof * sizeof(double);
				for (i = 0; i < dof; i++)
					x->ax_v[i] = (double *)mem, mem += dof * sizeof(double);
			}
			x->aaux = naux;				/* Number of auxiliaries allocated for */
		}
	
		/* Setup A[][] matrix to decompose, and figure number of auxiliaries actually needed */
		for (ee = naux = 0; ee < b->naux; ee++) {
			int ei = icomb[b->auxi[ee]];		/* Simplex relative auxiliary index */
			if (ei < 0)
				continue;		/* aux corresponds with fixed input value for this simplex */
			for (f = 0; f < dof; f++)
				x->ax_u[naux][f] = x->lo_l[ei][f];
			naux++;
		}
		x->naux = naux;					/* Number of auxiliaries actually available */
		x->auxbm = b->auxbm;			/* Mask of auxiliaries used */

		if (naux == dof) {				/* Use LU decomp to solve exactly */
			double rip;

#ifdef STATS
			x->c->s->rev.st[x->c->s->rev.sb->op].sinited5a++;
#endif /* STATS */
			if (lu_decomp(x->ax_u, dof, (int *)x->ax_w, &rip)) {
				x->flags |= SPLX_FLAG_5F;
				return 1;
			}

		} else if (naux > 0) {			/* Use SVD to solve least squares */

#ifdef STATS
			x->c->s->rev.st[x->c->s->rev.sb->op].sinited5b++;
#endif /* STATS */
			if (svdecomp(x->ax_u, x->ax_w, x->ax_v, naux, dof)) {
				x->flags |= SPLX_FLAG_5F;
				return 1;
			}
	
			/* Threshold the singular values W[] */ 
			svdthresh(x->ax_w, dof);
		} /* else naux == 0, don't setup anything */

		x->flags |= SPLX_FLAG_5;
	}
	return 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */

/* Initialise the given static sub-simplex verticy information table in s->r.sspxi[] */
static void
init_ssimplex_info(
rspl *s,
int sdi					/* Sub-simplex dimensionality (range 0 - di) */
) {
	int e, di = s->di;		/* Dimensionality */
	ssxinfo *xip;			/* Pointer to sub-simplex info structure */
	int vi, nospx;			/* Number of sub-simplexes */
	XCOMBO(vcmb, sdi+1, 1 << di);/* Simplex dimension sdi out of cube dimension di counter */

	xip = &s->rev.sspxi[sdi];
	if (xip->spxi != NULL)	/* Assert */
		error("rspl rev, internal, init_ssimplex_info called on already init'd\n");

//printf("~~ init_ssimplex_info called with sdi = %d\n",sdi);
	/* First count the number of sub-simplexes */
	nospx = 0;
	XCB_INIT(vcmb);
	while (!XCB_DONE(vcmb)) {
		nospx++;
		XCB_INC(vcmb);
	}

	xip->sdi = sdi;
	xip->nospx = nospx;
	if ((xip->spxi = (psxinfo *) calloc(nospx, sizeof(psxinfo))) == NULL)
		error("rspl malloc failed - reverse cell sub-simplex info array");
	
//printf("~~ no subsimplex = %d\n",nospx);
	/* For all sub-simplexes */
	XCB_INIT(vcmb);
	for (vi = 0; vi < nospx; vi++) {
		psxinfo *x = &xip->spxi[vi];
		int i;

		/* XCOMB generates verticies in order from max to min offest */

		/* Compute Absolute -> Parameter mapping */
		for (e = 0; e < di; e++) {				/* For each absolute axis */

			if ((vcmb[sdi] & (1<<e)) != 0) {
				x->icomb[e] = -2;	/* This abs is always '1' */

			} else if ((vcmb[0] & (1<<e)) == 0) {
				x->icomb[e] = -1;	/* This abs is always '0' */

			} else {
				for (i = 0; i < sdi; i++) {	/* For each verticy in large to small order (!first) */
					if ((vcmb[i]   & (1<<e)) != 0 && 
					    (vcmb[i+1] & (1<<e)) == 0) {/* Transition from offset 1 to 0 */
						x->icomb[e] = i;	/* This is parameter */
						break;
					}
				}
			}
		}
		
		/* Compute fwd grid offsets for each simplex vertex in baricentric order */
		for (i = 0; i <= sdi; i++) {	/* For each verticy */
			int pmin[MXRI], pmax[MXRI];
			x->offs[i]  = vcmb[i];
			x->foffs[i] = s->g.fhi[vcmb[i]];

			/* Setup input coordinate bounding box value offsets */
			if (i == 0) {								/* Init to first vertex of simplex */
				for (e = 0; e < di; e++) {				/* Input space */
					x->pmino[e] = x->pmaxo[e] = vcmb[i];
					pmin[e] = pmax[e] = vcmb[i] & (1<<e);
				}
			} else {
				for (e = 0; e < di; e++) {			/* Input space */
					double vv = vcmb[i] & (1<<e);
					if (vv < pmin[e]) {				/* Adjust min/max offsets */
						x->pmino[e] = vcmb[i];
						pmin[e] = vv;
					} else if (vv > pmax[e]) {
						x->pmaxo[e] = vcmb[i];
						pmax[e] = vv;
					}
				}
			}
		}

#ifdef NEVER
printf("~~Verticies   = ");
for (i = 0; i <= sdi; i++)
	printf("%d ",vcmb[i]);
printf("\n");

printf("~~Abs -> Parm = ");
for (e = 0; e < di; e++)
	printf("%d ",x->icomb[e]);
printf("\n");

printf("~~Offset      = ");
for (e = 0; e <= sdi; e++)
	printf("%d ",x->foffs[e]);
printf("\n");
printf("\n");
#endif /* NEVER */

		/* Increment the counter value */
		XCB_INC(vcmb);
	}
}

/* Free the given sub-simplex verticy information */
static void
free_ssimplex_info(
rspl *s,
int sdi					/* Sub-simplex dimensionality (range 1 - di) */
) {
	ssxinfo *xip;	/* Pointer to sub-simplex info structure */

	xip = &s->rev.sspxi[sdi];
	if (xip->spxi == NULL)	/* Assert */
		return;

	free(xip->spxi);
	xip->spxi = 0;
}

/* ====================================================== */
/* Reverse cell cache code                                */

static int increase_revcache(revcache *rc);

/* Allocate and initialise the reverse cell cache */
static revcache *
alloc_revcache(
rspl *s
) {
	revcache *rc;
	char *ev;
	int i;

	if ((rc = (revcache *) calloc(1, sizeof(revcache))) == NULL)
		error("rspl malloc failed - reverse cell cache");
	
	rc->s = s;		/* For stats */

	/* Check for environment variable tweak */
	rc->rae = REV_ALLOC_ENTRIES;
	if ((ev = getenv("ARGYLL_REV_CACHE_MULT")) != NULL) {
		double mm;
		mm = atof(ev);
		if (mm > 0.1 && mm < 20.0)
			rc->rae = (int)(rc->rae * mm + 0.5);
	}
//printf("~1 reverse cache entries = %d\n",rc->rae);

	/* Allocate the initial number of blocks of entries */
	if ((rc->cells = (cell **)malloc(rc->rae * sizeof(cell *))) == NULL)
		error("rspl malloc failed - reverse cell cache cells array");
	rc->ncells = rc->rae;

	rc->lrutop = rc->lrubot = NULL;		/* Nothing in lru list */

#ifdef ALLOC_ALL_CACHE
	while (rc->nextalloc < rc->rae)
		increase_revcache(rc);
#endif

	return rc;
}

/* Free the reverse cell cache */
static void
free_revcache(revcache *rc) {
	int i;
	cell *cp;

	/* Free any stuff allocated in the cell contents */
	for (cp = rc->lrubot; cp != NULL; cp = cp->lruup) {
		free_cell_info(cp);
	}

	/* Free cell allocations */
	for (i = 0; i < rc->nextalloc; i++) {
		free(rc->cells[i]);
	}

	free(rc->cells);
	free(rc);
}

/* Allocate another block of cells, and add them to the cache. */
/* Return non-zero if we ran out of room to lock another cell */
static int
increase_revcache(
revcache *rc
) {
	cell *nxcells;
	int i;

	if (rc->nextalloc >= rc->ncells) {
		return 1;

#ifdef NEVER
/*		error("rspl cell cache exausted!"); */
		/* We need to re-allocate some more. */
		rc->ncells += rc->rae;
		if ((rc->cells = (cell **)realloc(rc->cells, rc->ncells * sizeof(cell *))) == NULL)
			error("rspl realloc failed - reverse cell cache cells array");
#endif
	}

	if ((nxcells = (cell *) calloc(REV_CACHE_INC, sizeof(cell))) == NULL)
		error("rspl malloc failed - reverse cache cells");

	/* Add cells to the cache lru linked list */
	if (rc->lrutop == NULL)				/* List was empty */
		rc->lrutop = &nxcells[0];
	else {
		rc->lrubot->lrudown = &nxcells[0];	/* Splice into bottom */
		nxcells[0].lruup = rc->lrubot;
	}

	nxcells[0].s = rc->s;
	for (i = 1; i < REV_CACHE_INC; i++) {
		nxcells[i-1].lrudown = &nxcells[i];
		nxcells[i].lruup     = &nxcells[i-1];
		nxcells[i].s = rc->s;
	}
	rc->lrubot = &nxcells[REV_CACHE_INC-1];

	rc->cells[rc->nextalloc++] = nxcells;
	
//printf("~1 cache is now %d cells\n",rc->nextalloc * REV_CACHE_INC);
	return 0;
}

#define HASH(xx) ((xx) % REV_HASH_SIZE)

/* Return a pointer to an appropriate reverse cell */
/* cache structure. cell->flags will be 0 if the cell */
/* has been reallocated. cell contents will be 0 if */
/* never used before. */
/* The cell reference count is incremented, so that it */
/* can't be thrown out of the cache. The cell must be */
/* released with uncache_rcell() when it's no longer needed. */
/* return NULL if we ran out of room in the cache */
static cell *cache_rcell(
revcache *r,		/* Reverse cache structure */
int ix				/* fwd index of cell */
) {
	int hit = 0;
	int hash;
	cell *cp;
	
	hash = HASH(ix);		/* Compute hash of fwd cell index */

	/* See if we get a cache hit */
	for (cp = r->hashtop[hash]; cp != NULL; cp = cp->hlink) {
		if (ix == cp->ix) {	/* Hit */
			hit = 1;
#ifdef STATS
			r->s->rev.st[r->s->rev.sb->op].chits++;
#endif /* STATS */
			break;
		}
	}
	if (!hit) {			/* Re-use the least recently used cell */
		int ohash;
		cell *c;

		for (cp = NULL; cp == NULL; ) {
			/* Find the least recently used, non-referenced cell */
			for (cp = r->lrubot; cp != NULL && cp->refcount > 0; cp = cp->lruup)
				;
		
			if (cp == NULL) {
				if (increase_revcache(r))	/* Allocate more cache entries */
					return NULL;			/* Can't lock any more into cache */
			}
		}

#ifdef STATS
		r->s->rev.st[r->s->rev.sb->op].cmiss++;
#endif /* STATS */

		/* Remove from current hash index */
		ohash = HASH(cp->ix);			/* Old hash */
		if (r->hashtop[ohash] == cp) {
			r->hashtop[ohash] = cp->hlink;
		} else {
			for (c = r->hashtop[ohash]; c != NULL && c->hlink != cp; c = c->hlink)
				;
			if (c != NULL)
				c->hlink = cp->hlink;
		}

		/* Add this to hash index */
		cp->hlink = r->hashtop[hash];
		r->hashtop[hash] = cp;	/* Add to hash table and list */

		cp->ix = ix;
		cp->flags = 0;			/* Contents needs re-initializing */
	}
	
	/* Move slected cell to the top of the lru list */
	if (cp->lruup != NULL) {		/* This one wasn't already at top */
		cp->lruup->lrudown = cp->lrudown;
		if (cp->lrudown == NULL)	/* This was bottom */
			r->lrubot = cp->lruup;	/* New bottom */
		else
			cp->lrudown->lruup = cp->lruup;
		/* Put this one at the top */
		r->lrutop->lruup = cp;
		cp->lrudown = r->lrutop;
		r->lrutop = cp;
		cp->lruup = NULL;
	}
	cp->refcount++;

	return cp;
}

/* Invalidate the whole cache */
static void
invalidate_revcache(
revcache *rc)
{
	int i;
	cell *cp;

	/* Free any stuff allocated in the cell contents */
	for (cp = rc->lrubot; cp != NULL; cp = cp->lruup) {
		free_cell_info(cp);
		cp->refcount = 0;		/* Make sure they can now be reused */
		cp->ix = 0;
	}

	/* Clear the hash table so they can't be hit */
	for (i = 0; i < REV_HASH_SIZE; i++) {
		rc->hashtop[i] = NULL;
	}
}


/* Tell the cache that we aren't using this cell anymore */
static void uncache_rcell(
revcache *r,		/* Reverse cache structure */
cell *cp
) {
	if (cp->refcount > 0)
		cp->refcount--;
	else
		warning("rspl cell cache assert: refcount overdecremented!");
}

/* ====================================================== */
/* Reverse rspl setup functions                           */

static void free_rev_nn(rspl *s);

/* Initialise all three sections of the rev data in rspl. */
/* Note that reverse cell lookup tables are not */
/* allocated & created until the first call */
/* to a reverse interpolation function */
/* Init reverse interp elements in rspl */
void
init_rev(rspl *s) {

	/* First section */
	s->rev.inited = 0;
	s->rev.res = 0;
	s->rev.no = 0;
	s->rev.rev = NULL;
	s->rev.cache = NULL;

	/* Second section */
	s->rev.touch = NULL;

	/* Third section */
	s->rev.sb = NULL;

	s->rev_set_limit   = rev_set_limit_rspl;
	s->rev_get_limit   = rev_get_limit_rspl;
	s->rev_interp      = rev_interp_rspl;
	s->rev_locus       = rev_locus_rspl;
	s->rev_locus_segs  = rev_locus_segs_rspl;
}

/* Free up all the three sections of the reverse interpolation info */
void free_rev(
rspl *s		/* Pointer to rspl grid */
) {
	int e, di = s->di;
	int i,**rpp;
		
	/* If first section has been initialised */
	if (s->rev.inited != 0)	 {

#ifdef STATS
		{
		int totcalls = 0;
		for (i = 0; i < 5; i++) {
			totcalls += s->rev.st[i].searchcalls;
		}

		printf("\n===============================\n");
		printf("di = %d, do = %d\n",s->di, s->fdi);
		for (i = 0; i < 5; i++) {
			int calls = s->rev.st[i].searchcalls;
			if (calls == 0) 
				continue;
			printf("\n- - - - - - - - - - - - - - - -\n");
			printf("Operation %s\n",opnames[i]);
			printf("Search calls = %d = %f%%\n",s->rev.st[i].searchcalls,
			100.0 * s->rev.st[i].searchcalls/totcalls);
			printf("Cells searched/call = %f\n",s->rev.st[i].csearched/(double)calls);
			printf("Simplexes searched/call = %f\n",s->rev.st[i].ssearched/(double)calls);
			printf("Simplexes inited level 1/call = %f\n",s->rev.st[i].sinited/(double)calls);
			printf("Simplexes inited level 2 (LU)/call = %f\n",s->rev.st[i].sinited2a/(double)calls);
			printf("Simplexes inited level 2 (SVD)/call = %f\n",s->rev.st[i].sinited2b/(double)calls);
			printf("Simplexes invalidated level 4/call = %f\n",s->rev.st[i].sinited4i/(double)calls);
			printf("Simplexes inited level 4/call = %f\n",s->rev.st[i].sinited4/(double)calls);
			printf("Simplexes invalidated level 5/call = %f\n",s->rev.st[i].sinited5i/(double)calls);
			printf("Simplexes inited level 5 (LU)/call = %f\n",s->rev.st[i].sinited5a/(double)calls);
			printf("Simplexes inited level 5 (SVD)/call = %f\n",s->rev.st[i].sinited5b/(double)calls);
			if ((s->rev.st[i].chits + s->rev.st[i].cmiss) == 0)
				printf("No cache calls\n");
			else
				printf("Cell hit rate = %f%%\n",
					100.0 * s->rev.st[i].chits/(double)(s->rev.st[i].chits + s->rev.st[i].cmiss));
		}
		printf("\n===============================\n");
		}
#endif /* STATS */

		/* Sub-simplex information */
		for (e = 0; e <= di; e++) {
			free_ssimplex_info(s, e);
		}

		if (s->rev.rev != NULL) {
			/* Free arrays at grid points */
			for (rpp = s->rev.rev; rpp < (s->rev.rev + s->rev.no); rpp++) {
				if (*rpp != NULL)
					free(*rpp);
			}
			free(s->rev.rev);
			s->rev.rev = NULL;
		}
		s->rev.res = 0;
		s->rev.no = 0;

		free_revcache(s->rev.cache);	/* Reverse cell cache */
		s->rev.cache = NULL;

		s->rev.inited = 0;
	}

	/* Free up second section */
	free_rev_nn(s);					/* Free reverse nearest neighbor stuff */

	/* Free up third section */
	if (s->rev.sb != NULL) {
		free_search(s->rev.sb);
		s->rev.sb = NULL;
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Make reverse index list (rev First section init) */
/* This is called by a reverse interpolation call */
/* that discovers that the reverse index list haven't */
/* been initialised. */
static void make_rev(
rspl *s
) {
	int e,f;
	int di = s->di;
	int fdi = s->fdi;
	int rgno, gno = s->g.no;
	int rgres;
	int rgres_1;	/* rgres -1 == maximum base coord value */
	datao rgmin, rgmax;
	int i;			/* Index of fwd grid point */
	float *gp;		/* Pointer to fwd grid points */
	
	DBG(("make_rev called\n"));

	/* Reverse cell cache allocation */
	s->rev.cache = alloc_revcache(s);

	/* Sub-simplex information */
	for (e = 0; e <= di; e++) {
		init_ssimplex_info(s, e);
	}

	/* The reverse lookup relies on a search of the fwd interpolation tables. */
	/* To eliminate out of gamut points quickly, to provide */
	/* a starting point for the search, and to guarantee */
	/* that all possible reverse solutions are discovered, a spatial indexing */
	/* structure is used to provide a list of starting candidate forward indexes */
	/* for a given output value. */
	/* The reverse structure contains an fdi dimensional cell grid, each element of the */
	/* cell grid holding the indexes of the forward interpolation grid, which intersect */
	/* that ranges of output values. */
	/* Note that unlike the forward grid which us composed of verticies, */
	/* this grid is composed of cells. */

	s->get_out_range(s, rgmin, rgmax);	/* overall output min/max */

	/* Expand out range to encompass declared range */
	for (f = 0; f < fdi; f++) {
		if ((s->d.vl[f] + s->d.vw[f]) > rgmax[f])
				rgmax[f] = s->d.vl[f] + s->d.vw[f];
		if (s->d.vl[f] < rgmin[f])
				rgmin[f] = s->d.vl[f];
	}

	/* Expand out range slightly to allow for out of gamut points */
	for (f = 0; f < fdi; f++) {
		double del = (rgmax[f] - rgmin[f]) * 0.10;	/* Expand by +/- 10% */
		rgmax[f] += del;
		rgmin[f] -= del;
	}
//printf("~~got output range\n");

	/* Heuristic - main grid resolution ? */
	if ((rgres = s->g.mres + 1) < 4)
		rgres = 4;
	s->rev.res = rgres;			/* == number of cells per side */
	rgres_1 = rgres-1;

	for (rgno = 1, f = 0; f < fdi; f++, rgno *= rgres);	/* Number of elements in the rev.grid */
	s->rev.no = rgno;

	/* Compute coordinate increments */
	s->rev.coi[0] = 1;
	for (f = 1; f < fdi; f++)
		s->rev.coi[f] = s->rev.coi[f-1] * rgres;

	for (f = 0; f < fdi; f++) {
		s->rev.gl[f] = rgmin[f];
		s->rev.gh[f] = rgmax[f];
		s->rev.gw[f] = (rgmax[f] - rgmin[f])/(double)rgres;
	}

	if ((s->rev.rev = (int **) calloc(rgno, sizeof(int *))) == NULL)
		error("rspl malloc failed - rev.grid points");

	/*
	 * The grid contains pointers to lists of grid cube base indexes.
	 * If the pointer is NULL, then there are no base indexes in that list.
	 * A non NULL list uses the first element to indicate the alocation
	 * of the list. The last used entry for the list is always -1.
	 */

//printf("~~filling in reverse grid\n");
	/* For all grid points, form the cube with that point at its base, */
	/* and determine the bounding box of the output values that could */
	/* in that cube. */
	for (gp = s->g.a, i = 0; i < gno; gp += s->g.pss, i++) {
		int ee;
		datao min,max;
		int imin[MXRO], imax[MXRO], gc[MXRO];

//printf("~~~ i = %d/%d\n",i,gno);
		/* Skip cubes that are on the outside edge of the grid */
		for (e = 0; e < di; e++) {
			if(G_FL(gp, e) == 0)		/* At the top edge */
				break;
		}
		if (e < di) {	/* Top edge - skip this cube */
			continue;
		}

		/* Find the output value bounding box values for this grid cell */
		for (f = 0; f < fdi; f++)	/* Init output min/max */
			min[f] = max[f] = gp[f];
	
		/* For all other grid points in the cube */
		for (ee = 1; ee < (1 << di); ee++) {
			float *gt = gp + s->g.fhi[ee];	/* Pointer to cube vertex */
			
			/* Update bounding box for this grid point */
			for (f = 0; f < fdi; f++) {
				if (min[f] > gt[f])	
					 min[f] = gt[f];
				if (max[f] < gt[f])
					 max[f] = gt[f];
			}
		}

		/* Figure out intersection range in reverse grid */
		for (f = 0; f < fdi; f++) {
			double t;
			int mi;
			double gw = s->rev.gw[f];
			t = (min[f] - s->rev.gl[f])/gw;
			mi = (int)floor(t);			/* Grid coordinate */
			if (mi < 0)					/* Limit to valid cube base index range */
				mi = 0;
			else if (mi > rgres_1)
				mi = rgres_1;
			imin[f] = mi;	
			t = (max[f] - s->rev.gl[f])/gw;
			mi = (int)floor(t);			/* Grid coordinate */
			if (mi < 0)					/* Limit to valid cube base index range */
				mi = 0;
			else if (mi > rgres_1)
				mi = rgres_1;
			imax[f] = mi;	
		}

/* ~1 */
/*
printf("Scanning over grid:\n");
for (f = 0; f < fdi; f++) {
printf("Min[%d] = %d -> Max[%d] = %d\n",f,imin[f],f,imax[f]);
}
*/
		/* Now register forward index and vector with all the reverse grid cells */
		for (f = 0; f < fdi; f++)
			gc[f] = imin[f];	/* init coords */

		for (f=0; f < fdi;) {	/* For all of intersect cube */
			int **rpp, *rp;
			
			/* Compute pointer to grid cell */
			for (rpp = s->rev.rev, f = 0; f < fdi; f++)
				rpp += gc[f] * s->rev.coi[f];
			rp = *rpp;

/* ~1 */
/*
printf("Currently at grid:\n");
for (f = 0; f < fdi; f++) {
printf("gc[%d] = %d\n",f,gc[f]);
}
*/
			if (rp == NULL) {
				if ((rp = (int *) malloc(sizeof(int) * 5)) == NULL)
					error("rspl malloc failed - rev.grid entry");
				*rpp = rp;
				rp[0] = 5;		/* Actual allocation */
				rp[1] = i;
				rp[2] = -1;
			} else {
				int z, ll = rp[0];
				for (z = 1; z < ll; z++)	/* Find place to put next entry */
					if (rp[z] == -1)
						break;
				if (z >= (ll-1)) {			/* Not enough space */
					ll *= 2;
					if ((rp = (int *) realloc(rp, sizeof(int) * ll)) == NULL)
						error("rspl realloc failed - rev.grid entry");
					*rpp = rp;
					rp[0] = ll;
				}
				rp[z] = i;
				rp[z+1] = -1;
			}
			/* Increment index */
			for (f = 0; f < fdi; f++) {
				gc[f]++;
				if (gc[f] <= imax[f])
					break;	/* No carry */
				gc[f] = imin[f];
			}
		}	/* Next reverse grid point in intersecting cube */
	}	/* Next base grid point */
//printf("~~filled in reverse grid\n");

	s->rev.inited = 1;

	DBG(("make_rev finished\n"));
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* Make reverse nearest neighbor acceleration structures */
static void make_rev_nn(
rspl *s
) {
	/* The reverse clip nearest lookup relies on being able to quicly find the */
	/* nearest valid vertex value. We use an acceleration structure that minimised */
	/* the number of verticies examined in any detail. */
	/* We create fdi sorted axis projection lists, that allow us to */
	/* sweep an expanding cubic window centered on the target output value. */
	/* The touch bookkeeping tracks which verticy value are in the cube. */
	/* Hash randomization to improve hit rate: */
	static unsigned char key[] = { 0x43, 0x6F, 0x70, 0x79, 0x72, 0x69, 0x67, 0x68,
	                               0x74, 0x20, 0x32, 0x30, 0x30, 0x34, 0x20, 0x47,
	                               0x72, 0x61, 0x65, 0x6D, 0x65, 0x20, 0x57, 0x2E,
        	                       0x20, 0x47, 0x69, 0x6C, 0x6C, 0x00 };

	rev_struct *rs = &s->rev;
	int f, fdi = s->fdi;
	float *gbase = s->g.a;		/* Grid base address */
	int    gpss  = s->g.pss;	/* Grid increment between values */
	int gno = s->g.no;
	int i;			/* Index of fwd grid point */
	
	DBG(("make_rev_nn called\n"));

	/* Allocate the arrays spaces */
	for (f = 0; f < fdi; f++) {
		if ((rs->sax[f] = (float **)malloc(sizeof(float *) * gno)) == NULL)
			error("rspl rev malloc failed - make_rev_nn allocate sorted index array");
	}
	if ((rs->touch = (unsigned int *)calloc(sizeof(unsigned int), gno)) == NULL)
		error("rspl rev malloc failed - make_rev_nn allocate touch array");

	/* Fill out the axis arrays */
	for (i = 0; i < gno; i++) {
		for (f = 0; f < fdi; f++) {
			rs->sax[f][i] = gbase + i * gpss;	/* Address of point vector */
		}
	}

	/* Sort them */
	for (f = 0; f < fdi; f++) {
#define 	HEAP_COMPARE(A,B) (A[f] < B[f])
			HEAPSORT(float *, &rs->sax[f][0], gno)
#undef 		HEAP_COMPARE
	}

	DBG(("make_rev_nn finished\n"));
}

#define NN_INF 1e307

/* Return index of closest point to target */
/* Return -1 if not found */
int rev_nnfind(
rspl *s,
double rad,		/* If 0.0, return closest. If > 0.0 return next within radius squared */
double *q) {	/* Target point */
	rev_struct *rs = &s->rev;
	schbase *b = rs->sb;
	int f, fdi = s->fdi;
	float *gbase = s->g.a;		/* Grid base address */
	int    gpss  = s->g.pss;	/* Grid increment between values */
	int    gno = s->g.no;
	int e, i;

	if (rs->touch == NULL) {
		make_rev_nn(s);			/* Init search structure */
	}
	if (rad == 0.0 || rad != b->crad) {	/* Init search */

		rs->tbase += fdi;			/* Increment after last search */

		if ((rs->tbase + fdi) < rs->tbase) {	/* Overflow of touch count */
			memset((void *)rs->touch, 0, sizeof(unsigned int) * gno);
			rs->tbase = 0;
		}

		/* Find starting indexes within each axis array, */
		/* using binary search. */
		for (f = 0; f < fdi; f++) {
			int    i0, i1, i2;
			double v0, v1, v2;
			double ww;
	
			i0 = 0;
			i2 = gno - 1;
			v0 = rs->sax[f][i0][f];
			v2 = rs->sax[f][i2][f];
			if (q[f] <= v0) {
				i2 = i0;
				v2 = v0;
			}
			else if (q[f] >= v2) {
				i0 = i2;
				v0 = v2;
			}
			for (; (i2 - i0) > 1; ) {
				i1 = (i2 + i0)/2;		/* Trial point */
				v1 = rs->sax[f][i1][f];	/* Value at trial */
				if (v1 < q[f]) {
					i0 = i1;			/* Take top half */
					v0 = v1;
				} else {
					i2 = i1;			/* Take bottom half */
					v2 = v1;
				}
			}
	
			b->wex[f * 2 + 0] = i0;
			b->wex[f * 2 + 1] = i2;
	
			ww = v0 - q[f];
			b->wed[f * 2 + 0] = ww * ww;
			ww = v2 - q[f];
			b->wed[f * 2 + 1] = ww * ww;
		}
		b->bw = 0.0;		/* Current window distance */
		b->crad = rad;
	}

	/* Expand a fdi dimenstional cube centered on the target point, */
	/* jumping to the next nearest point on any axis, discovering */
	/* any points within the expanding window using the touch bookkeeping. */

	/* The first point found establishes the initial best distance. */
	/* When the window expands beyond the best distance, stop */

	{
		double bdist = 1e308;
		int bpoint = -1;

		/* Until we're done */
		for (;;) {
			int ee;			/* Axis & direction */
			int ff;			/* Axis */
			int ii;			/* Index of chosen point */
			float *fcb;		/* Pointer to chosen point */

			/* find next window increment */
			ee = 0;
			ii = b->wex[ee];
			b->bw = b->wed[ee];
			for (e = 1; e < (2 * fdi); e++) {
				if (b->wed[e] < b->bw) {
					ee = e;
					ii = b->wex[e];
					b->bw = b->wed[e];
				}
			}

			if (b->bw == NN_INF)
				break;				/* Can't expand any more */

			/* Chosen point on ee axis/direction, index ii */
			ff = ee / 2;			/* Axis only */

			/* Touch this new window boundary point */
			fcb = rs->sax[ff][ii];
			i = (fcb - gbase)/gpss;
			rs->touch[i] = (rs->touch[i] < rs->tbase ? rs->tbase : rs->touch[i]) + 1;

			/* Check the point out */
			if (rs->touch[i] == (rs->tbase + fdi)) {	/* Is within window on all axes */
				double tdist;

				for (tdist = 0.0, f = 0; f < fdi; f++) {
					double ww = rs->sax[ff][ii][f] - q[f];
					tdist += ww * ww;
				}
				if (tdist < bdist		/* Possible new best point */
				 && (b->limiten == 0 || get_limitv(b, i, fcb, NULL) <= b->limitv)) {
					bpoint = i;
					bdist = tdist;
					if (rad > 0.0 && tdist <= rad)		/* Return this next point */
						break;
				}
			}

			if (b->bw > bdist) {
				break;					/* Any further points will be worse, so done */
			}

			/* Increment next window edge candidate, and figure new edge distance */
			if (ee & 1) {					/* Top */
				if (++b->wex[ee] >= gno) {
					b->wed[ee] = NN_INF;
				} else {
					double ww = rs->sax[ff][b->wex[ee]][ff] - q[ff];
					b->wed[ee] = ww * ww;
				}
			} else {
				if (--b->wex[ee] < 0) {
					b->wed[ee] = NN_INF;
				} else {
					double ww = rs->sax[ff][b->wex[ee]][ff] - q[ff];
					b->wed[ee] = ww * ww;
				}
			}
		}

		return bpoint;
	}
}

static void free_rev_nn(
rspl *s
) {
	rev_struct *rs = &s->rev;
	int f, fdi = s->fdi;

	if (rs->touch != NULL) {		/* Has been initialised */

		for (f = 0; f < fdi; f++) {
			if (rs->sax[f] != NULL) {
				free(rs->sax[f]);
				rs->sax[f] = NULL;
			}
		}

		free(rs->touch);
		rs->touch = NULL;
	}
}

/* ====================================================== */

#undef DBGV
#undef DBG
#define DBGV(xxx)
#define DBG(xxx) 






