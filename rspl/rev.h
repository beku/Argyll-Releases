
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
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 * Latest simplex/linear equation version.
 */

/* Data structures used by reverse lookup code. */
/* Note that the reverse lookup code only supports a more limited */
/* dimension range than the general rspl code. */

/*
 * Note on simplex parameter space. 
 *
 * Simplex interpolation is normaly done in Baricentric space, where there
 * is one more baricentric coordinate than dimensions, and the sum of all
 * the baricentric coordinates must be 1. 
 *
 * To simplify things, we work in a "Simplex parameter" space, in which
 * there are only dimension parameters, and each directly corresponds
 * with a cartesian coordinate, but the parameter order corresponds with
 * the baricentric order.
 *
 * For example, given cartesian coordinates D0, D1, D2
 * these should be sorted from smallest to largest, thereby
 * choosing a particular simplex within a cube, and allowing
 * a correspondence to the parameter coordinates, ie:
 *
 * D1 D0 D2		Smallest -> Largest cartesian sort
 * P2 P1 P0 	Corresponding Parameter coordinates (note reverse order!)
 * 
 * B0 = P0		Conversion to Baricentric coordinates
 * B1 = P1 - P0
 * B2 = P2 - P1
 * B3 = 1  - P2
 *
 * The vertex values directly correspond to Baricentric coordinates,
 * giving the usual interpolation equation of:
 *
 *   VV0 * B0
 * + VV1 * B1
 * + VV2 * B2
 * + VV3 * B3
 *
 * Reversing the Parameter -> Baricentric equations gives the
 * following interpolation equation using Parameter coordinates:
 *
 *   VV0 - VV1 * P0
 * + VV1 - VV2 * P1
 * + VV2 - VV3 * P2
 * + VV3
 *
 * It is this which is used in rev.c for solving the reverse
 * interpolation problem.
 */

#undef STATS			/* Turn on rev interp stats */

/* A structure to hold per simplex coordinate and vertex mapping */
typedef struct {
	int icomb[MXRI];	/* Absolute[di]->Parameter[sdi], -1 == value 0, -2 == value 1 */
						/* Note that many Abs can map to one Param. */
	int offs[MXRI+1];	/* Offsets to simplex verticies within cube */
	int foffs[MXRI+1];	/* Fwd grid floating offsets to simplex verticies from cube base */
	int pmino[MXRI], pmaxo[MXRI]; /* Cube verticy offsets to setup simplex */
						/* pmin[] and pmax[] pointers. */
} psxinfo;

/* Sub simplexes of a cube information structure */
typedef struct {
	int sdi;			/* Sub-simplex dimensionality */
	int nospx;			/* Number of sub-simplexs per cube */
	psxinfo *spxi;		/* Per sub-simplex info array, NULL if not initialised */
} ssxinfo;

/* - - - - - - - - - - - - - - - - - - - - - */
/* Simplex definition. Each top level fwd interpolation cell, */
/* is decomposed into sub-simplexes. Sub-simplexes are of equal or */
/* lower dimensionality (ie. faces, edges, verticies) to the cube. */
typedef struct {
	struct _cell *c;			/* Pointer to parent cell */
	int si;						/* Diagnostic - simplex number within level */
	int sdi;					/* Sub-simplex dimensionality */
	int efdi;					/* Effective fdi. This will be = fdi for a non clip */
								/* plane simplex, and fdi+1 for a clip plane simplex */
								/* The DOF (Degress of freedom) withing this simplex = sdi - efdi */

	psxinfo *psxi;				/* Generic per simplex info */
	short flags;				/* Various flags */

#define SPLX_CLIPSX  0x01		/* This is a clip plane simplex */
#define SPLX_OVLIMIT 0x02		/* This non-clip plane simplex is over the ink limit */

#define SPLX_FLAG_1  0x04		/* v, linmin/max  initialised */
#define SPLX_FLAG_2  0x08		/* lu/svd initialised */
#define SPLX_FLAG_2F 0x10		/* lu/svd init. failed */
#define SPLX_FLAG_4  0x20		/* locus found */
#define SPLX_FLAG_5  0x40		/* auxiliary lu/svd initialised */
#define SPLX_FLAG_5F 0x80		/* auxiliary lu/svd init. failed */


#define SPLX_FLAGS  (SPLX_FLAG_1 | SPLX_FLAG_2 | SPLX_FLAG_2F \
                   | SPLX_FLAG_4 | SPLX_FLAG_5 | SPLX_FLAG_5F)

	double *v[MXRI+1];			/* Pointers to Vertex values for this simplex */
								/* v[0..sdi][0..fdi-1] are the output interpolation values */
								/* v[0..sdi][fdi] are the ink limit interpolation values */

								/* Baricentric vv[x][y] = (v[y][x] - v[y+1][x]) */
								/* and         vv[x][sdi] = v[sdi][x]           */

								/* Note that #num indicates appropriate flag number */
								/* and *num indicates a validator */

	double *pmin[MXRI], *pmax[MXRI]; /* Pointers to simplex vertex  */
								/* input space min and max values [di] */

	double *min[MXRO+1], *max[MXRO+1]; /* Pointers to simplex vertex [fdi+1] */
								/* and ink limit bounding values at *minmax[fdi] */

	/* If sdi == efdi, this holds the LU decomposition */
	/* else this holds the SVD and solution locus info */

	char *aloc2;		/* Memory allocation for #2 & #4 */

	/* double **d_u;	   LU decomp of vv, U[0..efdi-1][0..sdi-1]		#2 */
	/* int *d_w;		   LU decomp of vv, W[0..sdi-1]					#2 */

	double **d_u;		/* SVD decomp of vv, U[0..efdi-1][0..sdi-1]		#2 */
	double *d_w;		/* SVD decomp of vv, W[0..sdi-1]				#2 */
	double **d_v;		/* SVD decomp of vv, V[0..sdi-1][0..sdi-1]		#2 */

						/* Degrees of freedom = dof = sdi - efdi */
	double **lo_l;		/* Locus coefficients,  [0..sdi-1][0..dof-1]	#2 */

	double *lo_xb;		/* RHS used to compute lo_bd [0..efdi-1]		*4 */
	double *lo_bd;		/* Locus base solution, [0..sdi-1]				#4 */

	unsigned auxbm;		/* aux bitmap mask for ax_lu and ax_svd 		*5 */
	int      aaux;		/* naux count for allocation					*5 */
	int      naux;		/* naux for calculation (may be < aaux ?)		*5 */

	/* if (sdi-efdi = dof) == naux this holds LU of lo_l */
	/* else this holds the SVD of lo_l */

	char *aloc5;		/* Memory allocation for #5 */

	/* double **ax_u;	   LU decomp of lo_l							#5 */
	/* int *ax_w;		   Pivot record for ax_lu decomp				#5 */

	double **ax_u;		/* SVD decomp of lo_l, U[0..naux-1][0..dof-1]	#5 */
	double *ax_w;		/* SVD decomp of lo_l, W[0..dof-1]				#5 */
	double **ax_v;		/* SVD decomp of lo_l, V[0..dof-1][0..dof-1]	#5 */

} simplex;

/* A candidate search cell (cell cache entry structure) */
struct _cell {
	struct _rspl *s;		/* Pointer to parent rspl */

	/* Cache information */
	int ix;					/* Fwd cell index */
	struct _cell *hlink;	/* Link to other cells with this hash */
	struct _cell *lrudown;	/* Links to next least recently used cell */
	struct _cell *lruup;
	int refcount;			/* Reference count */
	int flags;				/* Non-zero if the cell has been initialised */
#define CELL_FLAG_1  0x01	/* Basic initialisation */
#define CELL_FLAG_2  0x02	/* Simplex information initialised */

	/* Use information */
	double sort;			/* Sort key */

	int    lwredge;			/* Lower edge flags, bits correspond with dimension */
	double limmin, limmax;	/* limit() min/max for this cell */
    double bcent[MXRO];		/* Output value bounding shere center */
    double brad;			/* Output value bounding shere radius */
    double bradsq;			/* Output value bounding shere radius squared */

	double p[POW2MXRI][MXRI]; /* Vertex input positions for this cube. */
							/* Used to by x->pmin/pmax[] & ink limit */

	double v[POW2MXRI][MXRO+1]; /* Vertex data for this cube. Used to by x->v[] */
							/* v[][fdi] is the ink limit values, if relevant */

	simplex *sx[MXRI+1];	/* Lists of simplexes that make up this cell. */
							/* Each list indexed by the non-limited simplex */
							/* dimensionality (similar to sspxi[]) */
							/* Each list will be NULL if it hasn't been created yet */
	int sxno[MXRI+1];		/* Corresponding count of each list */

}; typedef struct _cell cell;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */


/* The structure where cells are allocated and cached. */
/* Enough space is needed to cache all the cells/simplexes */
/* for a full aux. locus, or the query will be processed in */
/* several "chunks". */
/* This sets the basic memory usage of the rev code. Ideally */
/* it should set itself depending on the available RAM in the system. */
/* An environment variable "ARGYLL_REV_CACHE_MULT" can be set to a value > 1 to */
/* increase the ALLOC_ENTRIES. Note a value > 7 exceeds hash index size! */
#define REV_HASH_SIZE		39679	/* A prime is good (preferably >= CACHE_INC * ALLOC_ENTRIES) */
#define REV_CACHE_INC		512		/* Batch of cells to add to cache at a time */
#define REV_ALLOC_ENTRIES	10		/* Limit of allocation batches - sets base memory usage */
//#define REV_CACHE_INC		32		/* Test value:- Batch of cells to add to cache at a time */
//#define REV_ALLOC_ENTRIES	2		/* Test value:- Limit of allocation batches */

/* Holds the cell cache specific information */
typedef struct {
	struct _rspl *s;				/* Pointer to parent rspl */
	int rae;						/* possibly tweaked REV_ALLOC_ENTRIES */
	cell **cells;					/* The reverse cell allocations */
	int ncells;						/* Number of cells allocated */
	int  nextalloc;					/* Next cells to allocate */
	cell *hashtop[REV_HASH_SIZE];	/* Top of hash list */
	cell *lrutop, *lrubot;			/* Top and bottom pointers of LRU list */
} revcache;

/* common search information */
/* Type of (internal) reverse search */
enum ops {
	exact  = 0,		/* Search for all input values that exactly map to given output value */
	clipv  = 1,		/* Search for input values that map to outermost solution along a vector */
	clipn  = 2,		/* Search for input values that map to closest solution */
	auxil  = 3,		/* Search for input values that map to given output, and closest to auxiliary target */
	locus  = 4 };	/* Return range of auxiliary values that contains solution */

/* + possible exact with clip */

/* Structure to hold clip line state information */
typedef struct {
	struct _rspl *s;	/* Pointer to parent rspl */
	double st[MXRO];	/* start of line - reverse grid base value */
	double de[MXRO];	/* direction of line */
	int    di[MXRO];	/* incerement in line direction */
	int    ci[MXRO];	/* current rev grid coordinate */
	double t;			/* Parameter 0.0 - 1.0, line finished if t > 1.0 */
} line;


/* Structure to hold aux value of an intersection of a */
/* solution locus with a sub-simplex. Used when asegs flag is set */
typedef struct {
	double xval;		/* Auxiliary value */
	int nv;				/* Number of verticies valid */
	int vix[MXRI+1];	/* Verticy indexes of sub-simplex involved */
} axisec;

/* -------------------------------------------- */
/* Information needed/cached for reverse lookup */
struct _schbase {
	struct _rspl *s;	/* Pointer to parent rspl */

	int flags;			/* Hint flags */
	enum ops op;		/* Type of reverse search operation */
	int ixc;			/* Cube index of corner that holds maximum input values */

	int snsdi, ensdi;	/* Start and end extra sub-simplex dimensionality */
	int (*setsort)(struct _schbase *b, cell *c);	/* Function to triage & set cube sort index */
	int (*check)(struct _schbase *b, cell *c);		/* Function to recheck cube worth keeping */
	int (*compute)(struct _schbase *b, simplex *x);	/* Function to compute a simplex solution */

	double v[MXRO+1];	/* Target output value, + ink limit */
	double av[MXRI];	/* Target auxiliary values */
	int auxm[MXRI];		/* aux mask flags */
	unsigned auxbm;		/* aux bitmap mask */
	int naux;			/* Number of auxiliary target input values */
	int auxi[MXRI];		/* aux list of auxiliary target input values */
	double idist;		/* best input distance auxiliary target found (smaller is better) */

	double (*limit)(void *cntx, double *in);	/* Optional input space qualifier function. */
	void *cntx;			/* Context passed to limit() */
	double limitv;		/* Value not to be exceeded by limit() */
	int limiten;		/* Flag - limiting is enabled */

	int canvecclip;		/* Non-zero if vector clip direction usable */
	double cdir[MXRO];	/* Clip vector direction and length wrt. v[] */
	double ncdir[MXRO];	/* Normalised clip vector */
	double **cla;		/* Clip vector LHS implicit equation matrix [fdi][fdi+1] (inc. ink tgt.) */
	double clb[MXRO+1];	/* Clip vector RHS implicit equation vector [fdi+1] (inc. ink tgt.) */
	double cdist;		/* Best clip locus distance found (aim is min +ve) */
	int    iclip;		/* NZ if result is above (disabled) ink limit */

	int mxsoln;			/* Maximum number of solutions that we want */
	int nsoln;			/* Current number of solutions found */
	co *cpp;			/* Store solutions here */

	int   lxi;			/* Locus search axiliary index */
	double min, max;	/* current extreme locus values for locus search */
	int    asegs;		/* flag - find all search segments */
	int    axisln;		/* Number of elements used in axisl[] */
	int    axislz;		/* Space allocated to axisl[] */
	axisec *axisl;		/* Auxiliary intersections */

	int lclistz;		/* Allocated space to cell sort list */
	cell **lclist;		/* Sorted list of pointers to candidate cells */

	int pauxcell;		/* Indexe of previous call solution cell, -1 if not relevant */
	int plmaxcell;		/* Indexe of previous call solution cell, -1 if not relevant */
	int plmincell;		/* Indexe of previous call solution cell, -1 if not relevant */

	int lsxfilt;		/* Allocated space of simplex filter list */
	char *sxfilt;		/* Flag for simplexes that should be in a cell */

	double crad;			/* nn current radius distance */
	double bw;				/* nn current window distance */
	int wex[MXRO * 2];		/* nn current window edge indexes */
	double wed[MXRO * 2];	/* nn current window edge distances */

	int nnlistz;		/* space allocated to nnlist */
	int *nnlist;		/* list used to return nn cell list */

}; typedef struct _schbase schbase;

/* ----------------------------------------- */

#ifdef STATS
struct _stats {
	int 	searchcalls;/* Number of top level searches */
	int 	csearched;	/* Cells searched */
	int 	ssearched;	/* Simplexes searched */
	int		sinited;	/* Simplexes initialised to base level */
	int		sinited2a;	/* Simplexes initialised to 2nd level with LU */
	int		sinited2b;	/* Simplexes initialised to 2nd level with SVD */
	int		sinited4i;	/* Simplexes invalidated at 4th level */
	int		sinited4;	/* Simplexes initialised to 4th level */
	int		sinited5i;	/* Simplexes invalidated at 5th level */
	int		sinited5a;	/* Simplexes initialised to 5th level with LU */
	int		sinited5b;	/* Simplexes initialised to 5th level with SVD */
	int		chits;		/* Cells hit in cache */
	int		cmiss;		/* Cells misses in cache */
}; typedef struct _stats stats;
#endif /* STATS */

/* ----------------------------------------- */
/* Reverse info stored in main rspl function */
struct _rev_struct {

	/* First section */
	/* Has been initialised if inited != 0 */

	int inited;			/* Non-zero if this section has been initialised */

	/* Reverse grid lookup information */
	int res;			/* Reverse grid resolution == ncells on a side */
	int no;				/* Total number of points in reverse grid = rev.res ^ fdi */
	int coi[MXRO];		/* Coordinate increments */
	datao gl,gh,gw;		/* Reverse grid low, high, grid cell width */
	int  **rev;			/* rev.no pointers to lists of fwd grid indexes. */
						/* First int has number + allocation length */
						/* Last int is -1 */

	/* Sub-dimension simplex information */
	ssxinfo sspxi[MXRI+1];/* One per sub dimenstionality at offset sdi */

	revcache *cache;	/* Where cells are allocated and cached */

#ifdef STATS
	stats st[5];	/* Set of stats info indexed by enum ops */
#endif	/* STATS */


	/* Second section */
	/* Has been initialised if touch != NULL */

	/* Reverse nearest neighbor lookup information */
	float **sax[MXRO];	/* Sorted axis arrays, NULL if not inited. */
	unsigned int tbase;	/* Touch base value for this search */
	unsigned int *touch;/* Per value touch count */

	/* Third section */
	/* Has been initialise if sb != NULL */
	schbase *sb;		/* Structure holding calculated per-search call information */

}; typedef struct _rev_struct rev_struct;
















