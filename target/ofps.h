
#ifndef OFPS_H

/* 
 * Argyll Color Correction System
 *
 * Optimised Farthest Point Sampling
 *
 * Author: Graeme W. Gill
 * Date:   6/9/2004
 *
 * Copyright 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#ifndef MXPD
#define MXPD 4			/* Maximum ofps dimentionality */
#define POW2MXPD 16		/* 2 ^ MXPD */
#define POW3MXPD 81		/* 3 ^ MXPD */
#endif

/* A plane equation. The planes make up the surfaces of the */
/* Voronoi polyhedra. (Plane equations are duplicated, not shared) */
struct _peq {
	double pe[MXPD+1];	/* Vertex plane equation. First di elements are normalized */
						/* outward pointing normal (from first sample point), last element */
						/* is constant. If point . pe > 0, then point is outside surface */
	double cp[MXPD];	/* Closest point possible on plane to poi */
	int poi;			/* Sample point of interest */
	int ix;				/* Index of other sample point forming plane, -ve for gamut surface */
	int del;			/* Marked for deletion ? */
	double rads;		/* Radius squared of cp[] to sample point */
	struct _peq *link;	/* Linked list of free peq's */
}; typedef struct _peq peq;

/* A vertex. This is a point of intersection of the various planes, */
/* and are the vericies of the Voronoi polyhedra */
/* (Non gamut boundary verticies are shared) */
struct _vtx {
	double p[MXPD];		/* Vertex location */
	int nix[MXPD+1];	/* di+1 Sample point node indexes involved in vertex */
	double rads;		/* Radius squared to sample points on circimsphere */
	int vvalid;			/* Flag indicating v[], iv[] & eserr are valid */
	double v[MXPD];		/* Subjective value at vertex (Labj) ? */
	double iv[MXPD];	/* Interpolated subjective value */
	double rserr;		/* Raw (non padapt) estimated sampling error */
	double eserr;		/* Estimated sampling error */
	int del;			/* Marked for deletion ? */
	int refc;			/* Reference count */
	struct _vtx *link;	/* Linked list of free/used vtx's */
	struct _vtx **plp;	/* Pointer to link pointer in used list */

	struct _vtx *n;		/* Next in acceleration list */
	struct _vtx **pn;	/* Pointer to link pointer in acceleration list */
}; typedef struct _vtx vtx;

/* Pointer to vertex with pointers to planes involved */
typedef struct {
	vtx *p;
	peq *pp[MXPD];		/* Pointers to planes involved in vertex */
} vtxp;

/* A sample point node. */
/* Sample points are the points around which the Voronoi polyhedra */
/* are constructed. If a network is constructed between nodes that */
/* form a plane, the netork will be the Delaunay tesselation. */
struct _node {
	int    ix;			/* Index of node in s->n[] */
	int    fx;			/* nz if point is fixed */
	double p[MXPD];		/* Device coordinate position */
	double v[MXPD];		/* Subjective value (Labk) ? */

	int nvp;			/* Number of Voronoi surface planes */ 
	int _nvp;			/* Number allocated */
	peq **vp;			/* List of Voronoi surface planes */

	int nvv;			/* Number of Voronoi surface verticies */ 
	int _nvv;			/* Number allocated */
	vtxp *vv;			/* List of Voronoi surface verticies */

	double dmxs;		/* double maxiumum distance to sample point squared (for fast reject) */
	double op[MXPD];	/* Previous device coordinates during opt */

	struct _node *n;	/* Next in acceleration list */
}; typedef struct _node node;

/* An acceleration structure cube */
struct _acell {
	unsigned flag;		/* Access flag */
	node *head;			/* List of nodes inside acceleration cell */
	vtx  *vhead;		/* List of verticies with all real nodes inside acceleration cell */
	double p[MXPD];		/* Center of cube location */
}; typedef struct _acell acell;

/* A spiral list element */
struct _spir {
	int cid;			/* cell index displacement */
	double mpd;			/* Minimum possible distance squared */
}; typedef struct _spir spir;

/* Main sample point object */
struct _ofps {
/* private: */
	int di;			/* Point dimensionality */
	double ilimit;	/* Ink limit - limit on sum of p[] */
	double imin[MXPD];	/* Ink limit - limit on min of p[], must be >= 0.0 */
	double imax[MXPD];	/* Ink limit - limit on min of p[], must be <= 1.0  */
	int fnp;		/* Number of existing fixed points in list */
	int tinp;		/* target number of total points in list, including fixed points */
	int maxits;		/* Maximum itterative improvement passes */
	double padapt;	/* Perceptual adaptation, 0.0 = none */

	int np;			/* Number of points currently in list */
	node *n;		/* tinp list of points */

	/* Perceptual function handed in */
	void (*percept)(void *od, double *out, double *in);
	void *od;		/* Opaque data for perceptual point */
	
	/* Other info */
	int rix;			/* Next read index */
	double mn,mx,av;	/* Sampling distance stats */
	double pmn,pmx,pav;	/* Perceptually weighted sampling distance stats */
	double mxmvsq;		/* Maximum movement during optimisation */
	
	/* Gamut surface definition */
	node gn;			/* Voronoi surface for gamut */

	/* Currently used vertex */
	struct _vtx *uvtx;	/* Linked list of used vtx's */

	/* Free Plane and Vertex */
	struct _peq *fpeq;	/* Linked list of free peq's */
	struct _vtx *fvtx;	/* Linked list of free vtx's */

	/* Acceleration structure */
	int agres;		/* Acceleration grid resolution */
	double gw; 		/* Grid cell width */
	int gim[MXPD];	/* Grid index multiplier */
	double ghd;		/* Grid half diagonal */
	int nig;		/* Number of cells in grid */
	acell *grid;	/* Pointer to array of grid structures */
	unsigned flag;	/* Access flag */
	int nis;		/* Number of cells in spiral list */
	spir *spiral;	/* Spiral list */

/* public: */
	/* return non-zero if the perceptual point is within the device gammut */
	int (*pig)(struct _ofps *s, double *p);

	/* Initialise, ready to read out all the points */
	void (*reset)(struct _ofps *s);

	/* Read the next set of non-fixed points values */
	/* return non-zero when no more points */
	/* f may be NULL */
	int (*read)(struct _ofps *s, double *p, double *f);

	/* Calculate and print stats */
	void (*stats)(struct _ofps *s);

	/* Destroy ourselves */
	void (*del)(struct _ofps *s);

	}; typedef struct _ofps ofps;


/* Constructor */
extern ofps *new_ofps(int di, double ilimit, int npoints, double dadaptation,
	fxpos *fxlist, int fxno, 
	void (*percept)(void *od, double *out, double *in), void *od);

/* Extended constructor */
ofps *new_ofps_ex(
int di,					/* Dimensionality of device space */
double ilimit,			/* Ink limit (sum of device coords max) */
double *imin,			/* Ink limit - limit on min of p[], normally >= 0.0 */
double *imax,			/* Ink limit - limit on min of p[], normally <= 1.0  */
int tinp,				/* Total number of points to generate, including fixed */
double dadaptation,		/* Degree of adaptation to device characteristic 0.0 - 1.0 */
fxpos *fxlist,			/* List of existing fixed points (may be NULL) */
int fxno,				/* Number of existing fixes points */
void (*percept)(void *od, double *out, double *in),		/* Perceptual lookup func. */
void *od				/* context for Perceptual function */
);

/* ------------------------------------------------------- */
/* Macros for a di dimensional counter */
/* Declare the counter name nn, dimensions di, & count */

#define DCOUNT(nn, di, start, reset, count) 				\
	int nn[MXPD];	/* counter value */						\
	int nn##_di = (di);		/* Number of dimensions */		\
	int nn##_stt = (start);	/* start count value */			\
	int nn##_rst = (reset);	/* reset on carry value */		\
	int nn##_res = (count);	/* last count +1 */				\
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
	
/* ------------------------------------------------------- */
/* Macros combination counter */
/* Declare the counter name nn, combinations out of total */
/* Maximum combinations is DI+2 */

#define COMBO(nn, comb, total) 				\
	int nn[MXPD+2];			/* counter value */				\
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
	

#define OFPS_H
#endif /* OFPS_H */
