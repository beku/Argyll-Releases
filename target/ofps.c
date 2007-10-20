
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

/* TTBD:

   Itterative improvement for adapative doesn't work, so it's turned off,
   and you only get the initial placement if dadaptation > 0.0.

   Need to figure out how to get itterative improvement of adaptive working
   properly. The current itterative code depends on the space being Euclidean,
   so distorting the space metric for adaptive breaks it.

	See ~9
 */

/*
	Description:

		We build a Voronoi volume description around each sample point,
	so that we can locate the furthest point (vertex) from all existing sample
	points. Initially, we add sampling points at the largest distance verticies.
	This gives us an optimal distribution within a tollerance of 2:1
	We then iteratively improve the distribution by moving each point to be at
	the center of the minimal bounding sphere touching it's surrounding vertex points.

	Ideas for improving speed:

		Share planes to reduce perceptual computation.
		Share planes to allow vertexes to be computed once (reduce matrix soln time).
		Improve search for next location to seed point.
		Improve search for next point in neighborhood when computing Voronoi for added point.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include "numlib.h"
#include "sort.h"
#include "plot.h"
#include "icc.h"
#include "xcolorants.h"
#include "targen.h"
#include "ofps.h"

//#include <iperf.h>

#undef DEBUG

#ifndef NEVER	// Real settings

# define MAXITS 20		/* Number of optimisation itterations (0 to disable optimisation) */
# define AMAXITS 0		/* Adaptive number of itterations */
# define STOP_TOL 0.0005	/* Stopping tollerance */
# define SA_ADAPT 0.5	/* Standalone test, adaptation level */

#else			// development settings

# define MAXITS 1000		/* Number of optimisation itterations (0 to disable optimisation) */
# define AMAXITS 1000		/* Adaptive number of itterations */
# define STOP_TOL 0.00001	/* Stopping tollerance */

# define DUMP_PLOT		/* Show on screen plot */
# define DUMP_VTX 1		/* Display the vertex locations too */
# define PERC_PLOT 0	/* Emit perceptive space plots */
# define DO_WAIT 1		/* Wait for user key after each plot */
# define SA_ADAPT 0.9	/* Standalone test, adaptation level */

#endif /* NEVER */

#define NUMTOL 1e-10	/* Numerical tollerance */
#define DOACCEL			/* Use near cell grid acceleration */
#define TNPAGRID 1.0	/* Target nodes per acceleration grid cell */
#define FASTREJECT		/* Fast reject test for acceleration */
#define RANDOM_PERTERB	/* Perpterb initial placement to break up patterns */
#define PERTERB_AMOUNT 0.01
#define OPT_INITIAL_OVERSHOOT 1.0	/* Optimisation movement initial overshoot (1.5) */

#define ALWAYS
#undef NEVER

#ifdef NEVER
#ifdef	__STDC__
#include <stdarg.h>
void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);
#else
#include <varargs.h>
void error(), warning(), verbose();
#endif
#endif	/* NEVER */

#ifdef DUMP_PLOT
static void dump_image(ofps *s, int pcp, int dwt, int vtx);
#endif

static void ofps_stats(ofps *s);
static int ofps_point2cell(ofps *s, double *p);
static void ofps_add_vacc(ofps *s, vtx *vx);
static void ofps_rem_vacc(ofps *s, vtx *vx);

/* --------------------------------------------------- */
/* Default convert the nodes device coordinates into approximate perceptual coordinates */
/* (usually overriden by caller supplied function) */
static void
default_ofps_to_percept(void *od, double *p, double *d) {
	ofps *s = (ofps *)od;
	int e;

#ifndef NEVER
	/* Default Do nothing - copy device to perceptual. */
	for (e = 0; e < s->di; e++) {
		double tt = d[e];
		if (e == 0)
			tt = pow(tt, 2.0);
		else
			tt = pow(tt, 0.5);
		p[e] = tt * 100.0;
	}
#else
	for (e = 0; e < s->di; e++) {
		double tt = d[e];
		/* Two slopes with a sharp turnover in X */
		if (e == 0) {
			if (tt < 0.5)
				tt = tt * 0.3/0.5;
			else
				tt = 0.3 + ((tt-0.5) * 0.7/0.5);
		}
		p[e] = tt * 100.0;
	}
#endif
}

#ifdef NEVER /* Not currently used */
/* Return the distance of the device value from the device gamut */
/* This will be -ve if the point is outside */
/* If bvp is non-null, the index of the closest dim times 2 */
/* will be returned for the 0.0 boundary, dim * 2 + 1 for the 1.0 */
/* boundary, and di * 2 for the ink limit boundary. */
static double
ofps_in_dev_gamut(ofps *s, double *d, int *bvp) {
	int e;
	int di = s->di;
	double tt;
	double dd = 100.0;				/* Worst distance outside */
	double ss = 0.0;				/* Sum of values */
	int bv = di;
	for (e = 0; e < di; e++) {
		tt = d[e] - s->imin[e];
		if (tt < dd) {
			dd = tt;
			bv = e * 2;
		}
		tt = s->imax[e] - d[e];
		if (tt < dd) {
			dd = tt;
			bv = e * 2 + 1;
		}
		ss += d[e];					/* Track sum */
	}
	ss = (s->ilimit-ss)/di;	/* Axis aligned distance to ink limit */
	tt = sqrt((double)di) * ss;	/* Diagonal distance to ink limit */
	if (tt < dd) {
		dd = tt;
		bv = di * 2;
	}
	if (bvp != NULL)
		*bvp = bv;
	return dd;
}
#endif /* NEVER */

/* Given the new intended device coordinates, */
/* clip the new position to the device gamut edge */
/* return non-zero if the point was clipped */
static int
ofps_clip_point(ofps *s, double *d) {
	int e;
	double ss = 0.0;
	int rv = 0;
	for (e = 0; e < s->di; e++) {
		if (d[e] < s->imin[e]) {
			d[e] = s->imin[e];
			rv |= 1;
		} else if (d[e] > s->imax[e]) {
			d[e] = s->imax[e];
			rv |= 1;
		}
		ss += d[e];
	}
	if (ss > s->ilimit) {
		ss = (ss - s->ilimit)/s->di;
		for (e = 0; e < s->di; e++)
			d[e] -= ss;
		rv |= 1;
	}
	return rv;
}

/* --------------------------------------------------- */
/* Plane and vertex alloc/free support */

static peq *new_peq(ofps *s) {
	peq *p;

	if (s->fpeq != NULL) {	/* re-use one we've got */
		p = s->fpeq;
		s->fpeq = p->link;
		memset((void *)p, 0, sizeof(peq));

	} else {
		if ((p = (peq *)calloc(sizeof(peq), 1)) == NULL)
			error ("ofps: malloc failed on new vertex");
	}

	return p;
}

static void del_peq(ofps *s, peq *p) {
	p->link = s->fpeq;
	s->fpeq = p;
}

/* Get a new vertex */
static vtx *new_vtx(ofps *s) {
	vtx *p;

	if (s->fvtx != NULL) {	/* re-use one we've got */
		p = s->fvtx;
		s->fvtx = p->link;
		memset((void *)p, 0, sizeof(vtx));

	} else {
		if ((p = (vtx *)calloc(sizeof(vtx), 1)) == NULL)
			error ("ofps: malloc failed on new vertex");
	}

	/* Link vertex to currently used list */
	p->link = s->uvtx;
	if (s->uvtx != NULL)
		s->uvtx->plp = &p->link;
	s->uvtx = p;
	p->plp = &s->uvtx;

	return p;
}

/* Remove a vertx from the used list */
/* (Used for gamut boundary vertexes) */
static void remu_vtx(ofps *s, vtx *p) {
	if (p->plp != NULL) {		/* If is on used list, remove it */
		*p->plp = p->link;
		if (p->link != NULL)
			p->link->plp = p->plp;
	}
	p->plp = NULL;
	p->link = NULL;
}

static void del_vtx(ofps *s, vtx *p) {
	if (--p->refc <= 0) {

		if (p->plp != NULL) {		/* If is on used list, remove it */
			*p->plp = p->link;
			if (p->link != NULL)
				p->link->plp = p->plp;
		}

		p->link = s->fvtx;		/* Add to free list */
		p->plp = NULL;
		s->fvtx = p;
	}
}

/* Copy the contents of a vertex, but leave the used list pointers intact */
static void copy_vtx(vtx *dst, vtx *src) {
	vtx *link = dst->link;
	vtx **plp = dst->plp;

	*dst = *src;
	
	dst->link = link;
	dst->plp = plp;
}

/* --------------------------------------------------- */
/* Node basic support functions */

/* Free any allocated content of a node, but not the node itself. */
static void node_free(ofps *s, node *p) {
	int i;

	/* Free up list of Voronoi surface planes */
	if (p->vp != NULL) {
		for (i = 0; i < p->nvp; i++)
			del_peq(s, p->vp[i]);
		free(p->vp);
	}

	/* Free up list of Voronoi verticies */
	if (p->vv != NULL) {
		for (i = 0; i < p->nvv; i++)
			del_vtx(s, p->vv[i].p);
		free(p->vv);
	}
}

/* Add a plane to the node */
/* (Computes radius squared to node) */
static void node_add_plane(ofps *s, node *p, peq *pl) {
	int e, di = s->di;
	
	if (p->_nvp == 0) {
		p->_nvp = 4;		/* Initial allocation */
		if ((p->vp = (peq **)malloc(sizeof(peq *) * p->_nvp)) == NULL)
			error ("ofps: malloc failed on node plane pointers");
	} else if (p->nvp >= p->_nvp) {
		p->_nvp *= 2;		/* Double allocation */
		if ((p->vp = (peq **)realloc(p->vp, sizeof(peq *) * p->_nvp)) == NULL)
			error ("ofps: realloc failed on node plane pointers");
	}
#ifdef DEBUG
	{
		int e, di = s->di;
		printf("~1 Node 0x%x add pln 0x%x @ %d: ",p,pl,p->nvp);
		for (e = 0; e <= di; e++)
			printf("%f ",pl->pe[e]);
		printf("\n");
	}
#endif
	p->vp[p->nvp++] = pl;

	/* Compute the cp[]'s radius squared to the node position */
	for (pl->rads = 0.0, e = 0; e < di; e++) {
		double tt = p->p[e] - pl->cp[e];	/* sample point to Voronoi plane center */
		pl->rads += tt * tt;
	}
}

/* Delete any marked planes from the node (and free them) */
static void node_del_planes(ofps *s, node *p) {
	int i, j;
	
	for (j = i = 0; i < p->nvp; i++) {
		if (p->vp[i]->del == 0)
			p->vp[j++] = p->vp[i];
		else {
#ifdef DEBUG
		printf("~1 Node 0x%x deleted pln 0x%x @ %d\n",p,p->vp[i],i);
#endif
			del_peq(s, p->vp[i]);
		}
	}
	p->nvp = j;
}

/* Add a vertex pointer & vertex to the node */
/* (Increments reference count on vertex) */
static void node_add_vertex(ofps *s, node *pp, vtxp *vxp) {
	vtx *vx = vxp->p;

	if (pp->_nvv == 0) {
		pp->_nvv = 4;		/* Initial allocation */
		if ((pp->vv = (vtxp *)malloc(sizeof(vtxp) * pp->_nvv)) == NULL)
			error ("ofps: malloc failed on node vertex pointers");
	} else if (pp->nvv >= pp->_nvv) {
		pp->_nvv *= 2;		/* Double allocation */
		if ((pp->vv = (vtxp *)realloc(pp->vv, sizeof(vtxp) * pp->_nvv)) == NULL)
			error ("ofps: realloc failed on node vertex pointers");
	}
#ifdef DEBUG
	{
		int e, di = s->di;
		printf("~1 Node 0x%x add vtx @ %d: ",p,pp->nvv);
		for (e = 0; e < di; e++)
			printf("%f ",vx->p[e]);
		printf("from pln: ");
		for (e = 0; e < di; e++)
			printf("0x%x ",vxp->pp[e]);
		printf("\n");
	}
#endif

	vx->refc++;
	pp->vv[pp->nvv++] = *vxp;		/* Copy into place */

	/* Add vertex to acceleration structure if it is not there already. */
	if (vx->refc == 1)
		ofps_add_vacc(s, vx);
}

/* Delete any marked verticies from the node (and free them) */
static void node_del_verticies(ofps *s, node *p) {
	int i, j;
	
	for (j = i = 0; i < p->nvv; i++) {
		if (p->vv[i].p->del == 0) {
			p->vv[j++] = p->vv[i];

		} else {	/* Decrement reference count and delete if zero */
#ifdef DEBUG
			printf("~1 Node 0x%x deleted vtx 0x%x @ %d\n",p,p->vv[i],i);
#endif
			/* Remove vertex from acceleration structure if it will be deleted. */
			if (p->vv[i].p->refc == 1)
				ofps_rem_vacc(s, p->vv[i].p);
			del_vtx(s, p->vv[i].p);
		}
	}
	p->nvv = j;
}

/* Sort a vertex node index array. */
/* This is to speed up searching for a match */
/* Sort largest to smallest (so fake gamut nodes are last) */
static void sort_nix(ofps *s, int *nix) {
	int i, j, t;
	int di = s->di;		/* There are di+1 nodes */

	/* Do a really simple exchange sort */
	for (i = 0; i < di; i++) {
		for (j = i+1; j <= di; j++) {
			if (nix[i] < nix[j]) {
				t = nix[j]; 
				nix[j] = nix[i];
				nix[i] = t;
			}
		}
	}
}

/* Locate an existing vertex within a node, that has the same */
/* nodes that that determine it */
static vtx *locate_vertex(ofps *s, node *pp, int *nix) {
	int i;
	int e, di = s->di;

	for (i = 0; i < pp->nvv; i++) {
		vtx *vx = pp->vv[i].p;

		for (e = 0; e <= di; e++) {
			if (vx->nix[e] != nix[e])
				break;
		}
		if (e > di)
			break;		/* Found match */
	}
	if (i >= pp->nvv)
		return NULL;

	return pp->vv[i].p;
}

/* Given a target vertex and a base vertex, */
/* determine the extrapolated perceptual values */
/* of the first vertex. */
/* return nz if this can't be computed */
static int extrap_vtx(ofps *s, vtx *vx, vtx *vb) {
	int i, e, di = s->di;
	double **ta, *TTA[8], TA[8][8];
	double *tb, TB[8];
	node *np, *nn;
	int rv = 0;

	nn = &s->n[vb->nix[di]];		/* Last node */

	if (di <= 8) {
		for (e = 0; e < di; e++)
			TTA[e] = TA[e];
		ta = TTA;
		tb = TB;
	} else {
		ta = dmatrix(0, di-1, 0, di-1);
		tb = dvector(0, di-1);
	}

	/* Setup to compute the baricentric coordinates that */
	/* places vx->p in simplex defined by nodes around vb */
	for (i = 0; i < di; i++) {
		np = &s->n[vb->nix[i]];		/* node being delt with */
		for (e = 0; e < di; e++)
			ta[e][i] = np->p[e] - nn->p[e];
		tb[i] = vx->p[i] - nn->p[i];
	}
	/* Solve the simultaneous linear equations A.x = B */
	/* Return 1 if the matrix is singular, 0 if OK */
	if (solve_se(ta, tb, di)) {
		rv = 1;
	} else {
		double bn;

		/* Compute last baricentric weighting */
		for (bn = 1.0, i = 0; i < di; i++)
			bn -= tb[i];

		/* Use baricentric values to compute interpolated perceptual values */
		for (e = 0; e < di; e++)
			vx->iv[e] = bn * nn->v[e];

		for (i = 0; i < di; i++) {
			np = &s->n[vb->nix[i]];		/* node being delt with */
			for (e = 0; e < di; e++)
				vx->iv[e] += tb[i] * np->v[e];
		}
	}

	if (ta != TTA) {
		free_dmatrix(ta, 0, di-1, 0, di-1);
		free_dvector(tb, 0, di-1);
	}

	return rv;
}

/* Compute a vertexes estimates sampling error, and return it */
static double vtx_eserr(ofps *s, vtx *vx) {
	int i, e, di = s->di;
	double de;
	node *np;

	if (vx->vvalid)
		return vx->eserr;

	/* If not adaptive, then error is simple the radius of the circumsphere, so */
	if (s->padapt == 0.0) {
		vx->eserr = vx->rads;
		vx->vvalid = 1;
		return vx->eserr;
	}

	/* Lookup perceptual value at vertex location */
	s->percept(s->od, vx->v, vx->p);
	
	/* If this vertex is the center of a real simplex (and not on the gamut boundary) */
	/* Compute average of real sample points perceptual value associated with vertex */
	if (vx->nix[s->di] >= 0) {

		/* We use an equal weighting, because a vertex is at the intersection */
		/* of bisectors, and hence is at the circumcenter of the nodes. */
		for (e = 0; e < di; e++)
			vx->iv[e] = 0.0;
		for (i = 0; i <= di; i++) {
			np = &s->n[vx->nix[i]];
			for (e = 0; e < di; e++)
				vx->iv[e] += np->v[e];
		}
		for (e = 0; e < di; e++)
			vx->iv[e] /= (double)(di+1);

	} else {
		/* Find the closest non-gamut boundary simplex, and use it to extrapolate */

		int pci, tci, sli;		/* Point/Test cell/spiral list index */
		acell *cp;				/* Acceleration cell */
		double bdist = 1e80;	/* Best distance squared */
		vtx *bvx = NULL;		/* Best vertex */

		s->flag++;							/* Marker flag */
		pci = ofps_point2cell(s, vx->p);	/* Index of cell of interest */
	
		/* Search by searching grid cells in a spiral distance order for verticies */
		for (sli = 0; sli < s->nis; sli++) {
			vtx *v1;
			if (s->spiral[sli].mpd > bdist) {
				break;	/* This cell couln't possibly give closer vertex, so end search */
			}
	
			tci = pci + s->spiral[sli].cid;
			if (tci < 0 || tci >= s->nig)
				continue;					/* Outside boundary */
			cp = &s->grid[tci];
			if (cp->flag == s->flag)
				continue;					/* Already visited, must be clip wrap around */
			cp->flag = s->flag;				/* Mark this one as searched */
			for (v1 = cp->vhead; v1 != NULL; v1= v1->n) {
				double dd;
				for (dd = 0.0, e = 0; e < di; e++) {
					double tt = vx->p[e] - v1->p[e];
					dd += tt * tt;
				}
				if (dd < bdist) {
					bdist = dd;
					bvx = v1;
				}
			}
		}

		if (bvx != NULL) {
			/* Compute extrapolated value for vx from nodes around bvx */
			if (extrap_vtx(s, vx, bvx)) {
				/* Fail - fake zero error result */
				for (e = 0; e < di; e++)
					vx->iv[e] = vx->v[e];
			}

		} else {
			/* Failed to find any full verticies */
			/* Fake zero error result */
			for (e = 0; e < di; e++)
				vx->iv[e] = vx->v[e];
		}

		/* Fake result */
		for (e = 0; e < di; e++)
			vx->iv[e] = vx->v[e];
	}
	
	/* Compute estimated interpolation error */
	for (de = 0.0, e = 0; e < di; e++) {
		double tt = vx->v[e] - vx->iv[e];
		de += tt * tt;
	}
// ~~99
//	de *= 1.0/(100.0 * 100.0);		/* Scale perceptual to device levels */
	de *= 1.0/(100.0);		/* Scale perceptual to device levels */

	/* Blend distances squared acording to degree of adaptation */
	vx->rserr = de;
	vx->eserr = (1.0 - s->padapt) * vx->rads + s->padapt * de;
	vx->vvalid = 1;

	return vx->eserr;
}

/* --------------------------------------------------- */
/* Plane/vertex computation routines */

/* Given two sample point indexes, compute the plane between them. */
/* The first point is the point of interest. */
/* (This will fail with a divide by zero error if two points are coincident) */
static void comp_peq(ofps *s, peq *vp, int poi, int ix) {
	node *p0 = &s->n[poi], *p1 = &s->n[ix];
	int e, di = s->di;
	double sum = 0.0;

	/* Compute plane normal from poi to ix */
	for (e = 0; e < di; e++) {
		double tt = p1->p[e] - p0->p[e];
		vp->pe[e] = tt;
		sum += tt * tt;
	}
	sum = sqrt(sum);

	/* Normalise it */
	for (e = 0; e < di; e++)
		vp->pe[e] /= sum;

	/* Compute mid point */
	for (e = 0; e < di; e++) 
		vp->cp[e] = 0.5 * (p1->p[e] + p0->p[e]);

	/* Compute the plane equation constant */
	for (vp->pe[di] = 0.0, e = 0; e < di; e++) 
		vp->pe[di] -= vp->pe[e] * vp->cp[e];

	vp->poi = poi;		/* Sample point of interest */
	vp->ix = ix;		/* Index of sample point forming this plane with poi */
}

/* Compute a plane equation with the opposite orientation to */
/* the given plane. */
static void inv_peq(ofps *s, peq *ivp, peq *vp) {
	int e, di = s->di;

	for (e = 0; e <= di; e++)
		ivp->pe[e] = -vp->pe[e];

	for (e = 0; e < di; e++)
		ivp->cp[e] = vp->cp[e];

	ivp->poi = vp->ix;
	ivp->ix = vp->poi;
}

/* Given a vertex with it's di * pointers to planes set, */
/* compute the intersection point. */
/* return nz if there is no intersection */
static int comp_vtx(ofps *s, vtxp *vxp) {
	vtx *vx = vxp->p;
	int i, e, di = s->di;
	double **ta, *TTA[8], TA[8][8];
	int rv = 0;

	if (di <= 8) {
		for (e = 0; e < di; e++)
			TTA[e] = TA[e];
		ta = TTA;
	} else
		ta = dmatrix(0, di-1, 0, di-1);

	for (i = 0; i < di; i++) {
		for (e = 0; e < di; e++)
			ta[i][e] = vxp->pp[i]->pe[e];	/* Plane normal becomes row of matrix */
		vx->p[i] = -vxp->pp[i]->pe[di];		/* Plane constant becomes target */
	}
	/* Solve the simultaneous linear equations A.x = B */
	/* Return 1 if the matrix is singular, 0 if OK */
	if (solve_se(ta, vx->p, di))
		rv = 1;

	if (ta != TTA)
		free_dmatrix(ta, 0, di-1, 0, di-1);

	return rv;
}

/* Add a plane to a Voronoi surface around a test point */
/* return nz if it was added to the surface */
/* Note that numerical inacuracy can lead to strange results here, ie. */
/* planes being accepted but no corresponding verticies etc. */
static int add_to_vsurf(ofps *s,
int share,		/* NZ if resulting vertexes should be shared */
node *pp,		/* Node in question */
peq *vp			/* Plane to add */
) {
	int e, di = s->di;
	int i, j;
	vtxp vv;					/* vertex being added */
	int oneabove;
	COMBO(co, di-1, di);		/* di-1 out of di combination counter */

#ifdef DEBUG
	printf("~1 Adding plane 0x%x to node 0x%x\n",vp, pp);
	printf("   plane eqn = ");
	for (e = 0; e <= di; e++)
		printf("%f ",vp->pe[e]);
	printf("\n");
#endif

//PERF("Entering add_to_vsurf");

	vp->del = 0;

	/* Setup as if this plane will be accepted */
	/* Mark all existing planes for deletion */
	for (j = 0; j < pp->nvp; j++)
		pp->vp[j]->del = 1;

//PERF("Marked planes");
	/* Check all the existing verticies against the new plane, */
	/* and allow any planes that don't have a vertex below the */
	/* new plane to be marked for deletion. */
	for (oneabove = 0, i = 0; i < pp->nvv; i++) {
		vtxp *vxp = &pp->vv[i];
		vtx *vx = vxp->p;
		double v;
	
		for (v = vp->pe[di], e = 0; e < di; e++) 	/* Compute relation to plane equation */
			v += vp->pe[e] * vx->p[e];
		if (v > NUMTOL) {		/* Point is above plane */
			vx->del = 1;		/* Mark vertex for deletion */
			oneabove = 1;		/* Plane will be added */
		} else {
			vx->del = 0;		/* Vertex will remain */
			for (e = 0; e < di; e++)
				vxp->pp[e]->del = 0;		/* This plane is still being used */
		}
	}
//PERF("Checked existing verticies");
	if (oneabove == 0) {
#ifdef DEBUG
		printf("~1 Plane 0x%x doesn't intersect current surface",vp);
#endif
		return 0;				/* Plane wasn't added. caller should free it. */
	}

	/* delete vertex references to planes that will be deleted */
	for (i = 0; i < pp->nvv; i++) {
		vtxp *vxp = &pp->vv[i];
	
		for (e = 0; e < di; e++) {
			if (vxp->pp[e]->del != 0)
				vxp->pp[e] = NULL;
		}
	}

//PERF("Deleted plane references");

	/* Delete any planes that we don't want anymore */
	node_del_planes(s, pp);

//PERF("Deleted planes");

	/* Add out new plane to node */
	node_add_plane(s, pp, vp);
//PERF("Added plane");

	/* Compute new verticies due to the new plane */
	for (i = 0; i < pp->nvv; i++) {
		vtxp *vxp = &pp->vv[i];
		vtx *vx = vxp->p;

		if (vx->del == 0)		/* Vertex won't be deleted, so not cutoff by new plane */
			continue;

#ifdef DEBUG
		printf("~1 Check for new verts arnd del vtx 0x%x: ",vx);
		for (e = 0; e < di; e++)
			printf("%f ",vx->p[e]);
		printf("\n");
#endif
		/* Combine di-1 planes that formed to be deleted vertex */
		/* with out new plane, to compute potential new vertex */
		CB_INIT(co);

		while (!CB_DONE(co)) {
			int nix[MXPD+1];

			/* Setup to find existing vertex */
			for (e = 0; e < (di-1); e++) {
				if (vxp->pp[co[e]] == NULL)
					goto next_combo;			/* Involves a deleted plane */
				vv.pp[e] = vxp->pp[co[e]];		/* Existing planes */
				nix[e] = vxp->pp[co[e]]->ix;	/* Planes other sample point */
			}
			vv.pp[di-1] = vp;				/* And our new plane */
			nix[di-1] = vp->ix;				/* New planes other sample point */
			nix[di] = vp->poi;				/* This sample point */
			
			sort_nix(s, nix);

#ifdef DEBUG
			printf("~1 Testing for new vertex 0x%x from planes:\n",vx);
			for (j = 0; j < di; j++) {
				printf("   0x%x: ",vv.pp[j]);
				for (e = 0; e <= di; e++)
					printf("%f ",vv.pp[j]->pe[e]);
				printf("\n");
			}
#endif
			if (share) {
				/* Locate existing vertex that is circumcircle of the same nodes */

				if ((vv.p = locate_vertex(s, &s->n[vp->ix], nix)) != NULL) {
					node_add_vertex(s, pp, &vv);
				} else
					goto create_new;

			} else {
				/* We create a new vertex */
create_new:

				/* Allocate space for potential new vertex */
				vv.p = new_vtx(s);

				/* Setup new vertex */
				for (e = 0; e <= di; e++)
					vv.p->nix[e] = nix[e];
				
				/* Compute intersection point of di planes */
				if(comp_vtx(s, &vv)) {
#ifdef DEBUG
					printf("~1 No resulting vertex\n");
#endif
					del_vtx(s, vv.p);
					goto next_combo;		/* No intersection */
				}
		
#ifdef DEBUG
				printf("~1 Resulting location is: ");
				for (e = 0; e < di; e++)
					printf("%f ",vv.p->p[e]);
				printf("\n");
#endif
				/* Check if new vertex is within Voronoi surface */
				for (j = 0; j < pp->nvp; j++) {
					double v;
					peq *pl = pp->vp[j];
		
					/* Compute relation of new vertex to Voronoi planes */
					for (v = pl->pe[di], e = 0; e < di; e++)
						v += pl->pe[e] * vv.p->p[e];

#ifdef DEBUG
					printf("~1 Testing against plane 0x%x: ",pl);
					for (e = 0; e <= di; e++)
						printf("%f ",pl->pe[e]);
					printf(" result %f\n",v);
#endif
					if (v > NUMTOL) {			/* Point is above plane */
#ifdef DEBUG
						printf("~1 Vertex is outside Voronoi volume\n");
#endif
						del_vtx(s, vv.p);
						goto next_combo;	/* Not in Voronoi volume */
					}
				}

				/* Compute the vertexes radius squared to the node position */
				for (vv.p->rads = 0.0, e = 0; e < di; e++) {
					double tt = pp->p[e] - vv.p->p[e];
					vv.p->rads += tt * tt;
				}

				/* Add this new vertex */
				node_add_vertex(s, pp, &vv);
			}

			vxp = &pp->vv[i]; 	/* pp->vv[] may have been re-alloced after adding */
			vx = vxp->p;		/* vertex, so re-get pointers. */
next_combo:
			CB_INC(co);
		}
	}
//PERF("Added verticies");

	/* Now delete all the verticies that are no longer within Voronoi volume */
	node_del_verticies(s, pp);

//PERF("Deleted verticies");

//PERF("Exiting add_to_vsurf");
	return 1;					/* Plane was added. caller should lose reference to it */
}

/* Compute the Voronoi surface around the given sample point node from scratch. */
static void node_voronoi(
ofps *s,
int poi		/* Index of sample point to update/create Voronoi surface */
) {
	node *pp = &s->n[poi];		/* Node in question */
	peq *vp = NULL;				/* Next plane to check */
	int e, di = s->di;
	int i, j;
#ifdef DOACCEL
	int pci, tci, sli;		/* Point/Test cell/spiral list index */
	acell *cp;				/* Acceleration cell */
#endif /* DOACCEL */

//int nt, nr, nf, na;		// Stats

	/* Deconstruct the current Voronoi information associated with this sample node */
	while (pp->nvp > 0) {	/* Planes */
		del_peq(s, pp->vp[--pp->nvp]);
	};
	while (pp->nvv > 0) {	/* Verticies */
		del_vtx(s, pp->vv[--pp->nvv].p);
	};

	/* Misc other init */
	pp->ix = poi;

	/* Setup initial surface around this node by cloning the gamut surface */
	for (i = 0; i < s->gn.nvp; i++) {		/* Clone planes */
		peq *vp;							/* plane being added */

		vp = new_peq(s);
		*vp = *s->gn.vp[i];					/* Copy plane from gamut */
		vp->poi = poi;						/* Fix point of interest for cloned plane */
		node_add_plane(s, pp, vp);			/* Add plane to node */
	}
	pp->dmxs = -1e80;
	for (i = 0; i < s->gn.nvv; i++) {		/* Clone verticies */
		vtxp vv;							/* vertex being added */

		vv.p = new_vtx(s);
		copy_vtx(vv.p, s->gn.vv[i].p);		/* Copy vertex from gamut */
		for (e = 0; e < di; e++) {			/* Fix pointers to planes */
			vv.pp[e] = pp->vp[-s->gn.vv[i].pp[e]->ix-1];	/* Address of cloned plane */
			vv.p->nix[e] = vv.pp[e]->ix;				/* Planes other sample point */
		}
		vv.p->nix[di] = poi;				/* Fix point of interest in vertex */
		sort_nix(s, vv.p->nix);
		vv.p->refc = 0;

		/* Compute the vertexes radius squared to the node position */
		for (vv.p->rads = 0.0, e = 0; e < di; e++) {
			double tt = pp->p[e] - vv.p->p[e];
			vv.p->rads += tt * tt;
		}
		if (vv.p->rads > pp->dmxs)		/* Track furthest vertex distance for fast reject */
			pp->dmxs = vv.p->rads;

		node_add_vertex(s, pp, &vv);	/* Add vertex to node */
	}
	pp->dmxs *= 4.0;					/* Double squared */

//nt = nr = nf = na = 0.0;
	/* Go through the all other sample points accumulating their bisecting planes. */
#ifdef DOACCEL
	s->flag++;							/* Marker flag */
	pci = ofps_point2cell(s, pp->p);	/* Index of cell of interest */

	/* Accelerate search by searching grid cells in a spiral distance order */
	for (sli = 0; sli < s->nis; sli++) {
		node *p1;
		if (s->spiral[sli].mpd > pp->dmxs) {
			break;	/* This cell couln't possibly give closer point, so end search */
		}

		tci = pci + s->spiral[sli].cid;
		if (tci < 0 || tci >= s->nig)
			continue;					/* Outside boundary */
		cp = &s->grid[tci];
		if (cp->flag == s->flag)
			continue;					/* Already visited, must be clip wrap around */
		cp->flag = s->flag;				/* Mark this one as searched */
		for (p1 = cp->head; p1 != NULL; p1= p1->n) {
			int ix = p1->ix;			/* Index of point we're checking */

#else	/* !DOACCEL */
	/* Go through every cell */
	for (i = 0; i < s->np; i++) {
		int ix = i;
		node *p1 = &s->n[ix];		/* Node being considered */

#endif	/* !DOACCEL */
		double sum;
//nt++;
		if (ix == poi)
			continue;		/* Skip point of interest */

#ifdef FASTREJECT
		/* Fast reject. Reject if proposed node is further */
		/* away that double the furthest current vertex. */
		for (sum = 0.0, e = 0; e < di; e++)  {
			double tt = p1->p[e] - pp->p[e];
			sum += tt * tt;
		}
		if (sum > pp->dmxs) {
//nr++;
			continue;
		}
#endif /* FASTREJECT */

		if (vp == NULL)
			vp = new_peq(s);

//nf++;
		/* Compute the dividing plane between them */
		comp_peq(s, vp, poi, ix);

		/* Add plane to current surface */
		if (add_to_vsurf(s, 0, pp, vp)) {
//na++;
			/* Keep track of furthest vertex to assist fast reject */
			pp->dmxs = -1e80;
			for (j = 0; j < pp->nvv; j++) {
				vtx *vx = pp->vv[j].p;
				if (vx->rads > pp->dmxs)
					pp->dmxs = vx->rads;
			}
			pp->dmxs *= 4.0;	/* Double squared */

			vp = NULL;			/* Plane got used */
		}
	}
#ifdef DOACCEL
	}
#endif /* DOACCEL */
	if (vp != NULL)		/* Didn't use last plane */
		del_peq(s, vp);

//printf("~1 of %d, fast rej %d, full test %d, accepted %d\n",nt,nr,nf,na);
}

/* Update the current Voroni surface with one extra plane/node */
static void node_add_neigbor(
ofps *s,
peq *ovp		/* Inverse of plane to add */
) {
	node *pp = &s->n[ovp->ix];		/* Node of interest to add plane to */
	peq *vp;						/* Plane to add */

	vp = new_peq(s);

	/* Compute the properly oriented dividing plane between them */
	inv_peq(s, vp, ovp);

	/* Add plane to current surface */
	if (add_to_vsurf(s, 1, pp, vp) == 0) {
		del_peq(s, vp);		/* Didn't get used */

	} else {
		int i;

		/* Update the dmxs allowing for this new plane */
		pp->dmxs = -1e80;
		for (i = 0; i < pp->nvv; i++) {		/* For all the verticies */
			vtx *vx = pp->vv[i].p;
			if (vx->rads > pp->dmxs)
				pp->dmxs = vx->rads;
		}
		pp->dmxs *= 4.0;					/* Double squared */
	}
}

/* --------------------------------------------------- */

#ifdef NEVER
/* Structure to hold data for optimization function */
struct _edatas {
	ofps *s;			/* ppoint structure */
	node *pp;			/* Node in question */
	double *c;			/* Circumcenter */
}; typedef struct _edatas edatas;

/* Definition of the optimization functions handed to powell(.) */
/* We want to minimise the differences between the distance */
/* weighted vertex eserr values, in an attempt to make them equal. */
static double efunc1(void *edata, double p[]) {
	edatas *ep = (edatas *)edata;
	ofps *s = ep->s;
	int e, di = s->di;
	node *pp = ep->pp;
	double mnerr = 1e80;
	double mxerr = -1e80;
	double avg = 0.0;
	int i;

	/* For each vertex around this sampling node */
	for (i = 0; i < pp->nvv; i++) {
		vtx *vv = pp->vv[i].p;
		double sum;

		/* Compute distance to vertex */
		for (sum = 0.0, e = 0; e < di; e++) {
			double tt = p[e] - vv->p[e];
			sum += tt * tt;
		}
		sum = sqrt(sum);
		avg += sum;
	}
	avg /= (double)pp->nvv;
	
	/* For each vertex around this sampling node */
	for (i = 0; i < pp->nvv; i++) {
		vtx *vv = pp->vv[i].p;
		double sum, werr;

		/* Compute distance to vertex */
		for (sum = 0.0, e = 0; e < di; e++) {
			double tt = p[e] - vv->p[e];
			sum += tt * tt;
		}
		sum = sqrt(sum);

		/* Computed distance weighted estimated sampling error */
		werr = sum/avg * vv->rserr;
//		werr = sum;				/* Just distance = center */
		if (werr > mxerr)
			mxerr = werr;
		if (werr < mnerr)
			mnerr = werr;
	}
//printf("~9 efunc returning %f\n",mxerr - mnerr);
	return mxerr - mnerr;
}
#endif

/* Move the node to optimise it location amongst the surrounding */
/* verticies. */
/* When there is no adaptation, we compute a minimal bounding sphere */
/* that encloses the vertex points. */
/* */
/* When we are using adaptation, find a node location that minimises the */
/* maximum inverse distance weighted vertex eserr value. (This doesn't work!) */
// ~~99
/* Use -1 to compute center for gamut "node" */
static void comp_opt(ofps *s, int poi, double oshoot) {
	node *pp;						/* Node in question */
	int e, di = s->di;
	double radsq = -1.0;			/* Span/radius squared */
	double rad;
	double sum;
	int i, j;
	int bi = 0, bj = 0;

	if (poi < 0)
		pp = &s->gn;			/* gamut "node" */
	else
		pp = &s->n[poi];		/* Node in question */

	for (e = 0; e < di; e++)
		pp->op[e] = pp->p[e];		/* record previous position */

	if (pp->nvv > 0) {

		/* Locate a center point that minimises the maximum distance of an enclosing sphere */
		/* Find the two verticies that are farthest apart. Brute force search */
		for (i = 0; i < (pp->nvv-1); i++) {
			for (j = i+1; j < pp->nvv; j++) {
				for (sum = 0.0, e = 0; e < di; e++) {
					double tt = pp->vv[i].p->p[e] - pp->vv[j].p->p[e];
					sum += tt * tt;
				}
				if (sum > radsq) {
					radsq = sum;
					bi = i;
					bj = j;
				}
			}
		}
		
		/* Set initial bounding sphere */
		for (e = 0; e < di; e++)
			pp->p[e] = 0.5 * (pp->vv[bi].p->p[e] + pp->vv[bj].p->p[e]);
		radsq /= 4.0;			/* diam^2 -> rad^2 */

		rad = sqrt(radsq);
		
		/* Go though all the points again, expanding sphere if necessary */
		for (i = 0; i < pp->nvv; i++) {

			if (i == bi || i == bj)
				continue;

			/* Compute distance squared of vertex to bounding sphere center */
			for (sum = 0.0, e = 0; e < di; e++) {
				double tt = pp->vv[i].p->p[e] - pp->p[e];
				sum += tt * tt;
			}
			if (sum > radsq) {
				double tt;

				sum = sqrt(sum) + 1e-10;			/* Radius to point */
				rad = 0.5 * (rad + sum);
				radsq = rad * rad;
				tt = sum - rad;
				for (e = 0; e < di; e++)
					pp->p[e] = (rad * pp->p[e] + tt * pp->vv[i].p->p[e])/sum;
			}
		}

		/* If adaptive, move towards vertex with highest eserr */
		if (s->padapt > 0.0 && poi >= 0) {
			double sum;

			/* For each surrounding vertex, compute the distance weighted eserr, */
			/* and bias the enclosing sphere center point */

			/* Make sure eserr is valid for all the verticies */
			for (sum = 0.0, i = 0; i < pp->nvv; i++) {
				vtx *vv = pp->vv[i].p;
				vtx_eserr(s, vv);		/* Estimated sampling error */
				sum += vv->rserr;		/* Track sum of weights */
			}

			if (sum > 0.0) {
				double wp[MXPD];

				for (e = 0; e < di; e++)
					wp[e] = 0.0;

				/* Apply weighted vector towards verticies */
				for (i = 0; i < pp->nvv; i++) {
					vtx *vv = pp->vv[i].p;
			
					for (e = 0; e < di; e++) {
						wp[e] += vv->rserr/sum * vv->p[e];
					}
				}
				for (e = 0; e < di; e++)
					pp->p[e] = 0.5 * pp->p[e] + 0.5 * wp[e];
			}
		}
	}

	/* Apply overshoot */
	if (oshoot > 1.0) {
		for (e = 0; e < di; e++)
			pp->p[e] = (pp->p[e] - pp->op[e]) * oshoot + pp->op[e];
	}

	/* Clip the new location */
	ofps_clip_point(s, pp->p);

	/* Compute how far the point has moved */
	for (sum = 0.0, e = 0; e < di; e++) {
		double tt = pp->op[e] - pp->p[e];
		sum += tt * tt;
	}
	if (sum > s->mxmvsq)	/* Track maximum movement */
		s->mxmvsq = sum;
}

/* --------------------------------------------------- */
/* Setup gamut Voronoi surface */
static void ofps_binit(ofps *s) {
	int e, di = s->di;
	int doink = 0;
	node *gn = &s->gn;		/* Gamut surface "node" */
	peq *vp;				/* plane being added */
	vtxp vv;				/* vertex being added */
	int i;					/* plane combination */
	DCOUNT(co, di, 0, 0, 2);	/* Count through corner verticies */

	if (s->ilimit < (double)di) 	/* Ink limit is active */
		doink = 1;

	/* Init the plane equations */
	for (i = 0; i < (2 * di); i++) {	/* unit cell at 0 */
		int ii = i >> 1;				/* Dimension */

		vp = new_peq(s);
		for (e = 0; e < di; e++)
			vp->pe[e] = 0.0;
		vp->pe[ii] = i & 1 ? 1.0 : -1.0;					/* Normal */
		vp->pe[di] = i & 1 ? -s->imax[ii] : s->imin[ii];	/* Constant */

		/* Number the planes to aid cloning plus sample point index */
		vp->poi = -999;				/* fake sample point node of interest */
		vp->ix = -i-1;				/* -1 to -2di fake other nodes */

		node_add_plane(s, gn, vp);			/* Add plane to the "node" */
	}
	/* Add initial verticies. Create vertex for each di combination of planes */
	DC_INIT(co);

	while (!DC_DONE(co)) {
		vv.p = new_vtx(s);

		/* Compute vertex location */
		for (e = 0; e < di; e++) {
			if (co[e] != 0)
				vv.p->p[e] = s->imax[e];
			else
				vv.p->p[e] = s->imin[e];
		}

		/* Compute planes and nodes involved */
		for (e = 0; e < di; e++) {
			if (co[e] != 0) {
				vv.pp[e]     = gn->vp[(e << 1) + 1];
				vv.p->nix[e] = gn->vp[(e << 1) + 1]->ix;
			} else {
				vv.pp[e]     = gn->vp[(e << 1) + 0];
				vv.p->nix[e] = gn->vp[(e << 1) + 0]->ix;
			}
		}
		vv.p->nix[di] = -999;

		node_add_vertex(s, gn, &vv);		/* Add vertex to the "node" */

		DC_INC(co);
	}

	if (doink) {	/* Ink limit plane is orthogonal to diagonal */
		double len;
		vp = new_peq(s);
		len = 1.0/sqrt((double)di);	/* Normalised length */
		for (e = 0; e < di; e++)
			vp->pe[e] = len;
		vp->pe[di] = -s->ilimit * len;

		/* Number the planes to aid cloning plus sample point index */
		vp->poi = -999;
		vp->ix = -2 * di -1;

		/* Add plane to current surface, recomputing verticies etc. */
		if (add_to_vsurf(s, 0, gn, vp) == 0)
			del_peq(s, vp);		/* Didn't get used */
	}

	/* Remove the gamut surface vertexes from used list, */
	/* so they don't get scanned for next node location. */
	for (i = 0; i < gn->nvv; i++) {
		remu_vtx(s, gn->vv[i].p);
	}
}

/* --------------------------------------------------- */
/* Setup the acceleration structure */
static void
ofps_init_acc(ofps *s) {
	int i, e;
	int di = s->di;

	/* Create acceleration grid array */

	/* Choose a grid resolution that aims for aproximately TNPAGRID nodes per grid */
	s->agres = (int)(pow(s->tinp/TNPAGRID, 1.0/di) + 0.5);
	if (s->agres < 1)
		s->agres = 1;

#ifdef DEBUG
	printf("~1 accel grid res = %d\n",s->agres);
#endif

	/* Cell width */
	s->gw = 1.0/s->agres;

	/* Compute grid index multipliers */
	for (s->gim[0] = 1, e = 1; e < di; s->gim[e] = s->gim[e-1] * s->agres, e++)
		;
	
	/* Compute half cell diagonal */
	s->ghd = sqrt((double)di * s->gw * s->gw/4.0);

	/* Compute number of cells in grid */
	for (s->nig = 1, e = 0; e < di; e++)
		s->nig *= s->agres;

	/* Allocate grid */
	if ((s->grid = (acell *)malloc(sizeof(acell) * s->nig)) == NULL)
		error ("ofps: malloc failed for acceleration grid");

	s->flag = 0;

	/* Initialise grid */
	{
		DCOUNT(co, di, 0, 0, s->agres);

		i = 0;
		DC_INIT(co);
		while (!DC_DONE(co)) {
			acell *cp = &s->grid[i];

			cp->flag = 0;
			cp->head = NULL;
			cp->vhead = NULL;

			for (e = 0; e < di; e++)	/* Center of cell */
				cp->p[e] = (co[e] + 0.5) * s->gw;

			DC_INC(co);
			i++;
		}
	}

	/* Create the spiral list */

	/* Sriral list is hypercube of side 2 * s->agres -1 */
	for (s->nis = 1, e = 0; e < di; s->nis *= (2 * s->agres -1), e++)
		;
	s->nis++;		/* One more for guard point */

	if ((s->spiral = (spir *)malloc(sizeof(spir) * s->nis)) == NULL)
		error ("ofps: malloc failed on spiral list");

	/* Initialise list from cube */
	{
		DCOUNT(co, di, -s->agres+1, -s->agres+1, s->agres);
		
		i = 0;
		DC_INIT(co);
		while (!DC_DONE(co)) {
			int io;
			double ss, tt;

			/* Calc index offset */
			for (io = 0, ss = 0.0, e = 0; e < di; e++) {
				io += s->gim[e] * co[e];
				tt = (double)co[e];
				if (tt > 1.0)
					tt -= 1.0;
				else if (tt < -1.0)
					tt += 1.0;
				else
					tt = 0.0;
				ss += tt * tt; /* Axis distance from corners of center to current cell */
			}
			s->spiral[i].cid = io;
			s->spiral[i].mpd = ss * s->gw * s->gw;

			DC_INC(co);
			i++;
		}

		/* Add guard point */
		s->spiral[i].cid = 0;
		s->spiral[i].mpd = 1e100;
	}

	/* Sort the list in order of mpd */
#define HEAP_COMPARE(A,B) ((A).mpd < (B).mpd)
	HEAPSORT(spir, s->spiral, s->nis);
#undef HEAP_COMPARE

#ifdef NEVER
	for (i = 0; i < s->nis;i++) {
		printf("Spr[%d] = %d, %f\n",i,s->spiral[i].cid,s->spiral[i].mpd);
	}
#endif /* NEVER */
}

/* Convert a location to an acceleration cell index */
static int
ofps_point2cell(ofps *s, double *p) {
	int i, e, di = s->di;
	int agres = s->agres;

	for (i = e = 0; e < di; e++) {
		int t;
		t = (int)floor(agres * p[e]);
		if (t < 0)
			t = 0;
		else if (t >= agres)
			t = (agres-1);
		i += s->gim[e] * t;
	}
	return i;
}

/* Add a node to the acceleration structure */
static void
ofps_add_nacc(ofps *s, node *n) {
	acell *cp = &s->grid[ofps_point2cell(s, n->p)];
	n->n = cp->head;
	cp->head = n;
}

/* Add a vertex to the acceleration structure (if not boundary vertex) */
static void
ofps_add_vacc(ofps *s, vtx *vx) {
	if (vx->nix[s->di] < 0) {
		vx->pn = NULL;
		vx->n = NULL;
		return;		
	} else {
		acell *cp = &s->grid[ofps_point2cell(s, vx->p)];
		vx->n = cp->vhead;
		if (cp->vhead != NULL)
			cp->vhead->pn = &vx->n;
		cp->vhead = vx;
		vx->pn = &cp->vhead;
	}
}

/* Remove a vertex from the acceleration structure */
static void
ofps_rem_vacc(ofps *s, vtx *vx) {
	if (vx->pn != NULL) {		/* If is on acceleration list, remove it */
		*vx->pn = vx->n;
		if (vx->n != NULL)
			vx->n->pn = vx->pn;
	}
	vx->pn = NULL;
	vx->n = NULL;
}

/* Clear the acceleration structure */
static void
ofps_reset_acc(ofps *s) {
	int i;

	s->flag = 0;
	for (i = 0; i < s->nig; i++) {
		acell *cp = &s->grid[i];
		cp->flag = 0;
		cp->head = NULL;
		cp->vhead = NULL;
	}
}


/* --------------------------------------------------- */

/* Seed the object with any fixed points */
static void
ofps_add_fixed(
ofps *s,
fxpos *fxlist,			/* List of existing fixed points */
int fxno				/* Number in fixed list */
) {
	int e, di = s->di;
	int i, j, ii;

	/* Add fixed points if there are any */
	if (fxno > 0) {

		for (i = 0; (i < fxno) && (i < s->tinp); i++) {
			node *p = &s->n[s->np];	/* Destination for point */

			/* Reject any duplicate points, or Voronoi will get confused.. */
			for (ii = 0; ii < s->np; ii++) {
				for (e = 0; e < di; e++) {
					if (fabs(s->n[ii].p[e] - fxlist[i].p[e]) > 1e-5)
						break;		/* Not a match */
				}
				if (e >= di)
					break;			/* Is a match */
			}
			if (ii < s->np)
				continue;			/* Skip adding this point */

			s->fnp = ++s->np;

			for (e = 0; e < di; e++)
				p->op[e] = p->p[e] = fxlist[i].p[e];

			p->fx = 1;			/* is a fixed point */
			s->percept(s->od, p->v, p->p);

			/* Compute the Voronoi for it */
			node_voronoi(s, s->np-1);

			/* Add bisecting plane due to new node to all of its neighbors */
			for (j = 0; j < p->nvp; j++) {
				peq *vp = p->vp[j];				/* plane being added */
				if (vp->ix >= 0)				/* If not gamut plane */
					node_add_neigbor(s, vp);	/* add new point as neighbor */
			}

			/* Add node to acceleration structure */
			ofps_add_nacc(s, p);
			printf("\rAdded fixed %d/%d",i,fxno); fflush(stdout);
		}
		/* Adjust target points for number of duplicate fixed points */
		s->tinp -= (fxno - s->fnp);
	}
}

/* Seed the object with any movable incremental farthest points. */
static void
ofps_seed(ofps *s) {
	int e, di = s->di;
	int i, j;
	double rerr;

	printf("\n");

	/* Seed the non-fixed points */
	for (j = 0, i = s->fnp; i < s->tinp; i++, j++) {
		node *p = &s->n[i];		/* New node */

		if (i == 0) {			/* No initial fixed points */
			/* If there are no fixed points, place first point */
			/* in the center of the gamuts Voronoi verticies. */
			comp_opt(s, -1, 1.0);

			/* Add the new point */
			for (e = 0; e < di; e++)
				p->op[e] = p->p[e] = s->gn.p[e];

			s->percept(s->od, p->v, p->p);
			s->np = i+1;

			/* Compute the Voronoi for it */
			node_voronoi(s, i);

		} else {
			int k;

			double mx = -1e80;
			vtx *vx, *bvx = NULL;

			/* Locate the Voronoi vertex with the greatest distance to a sampling points */
			/* (Might be speedup if keep vertex worst distance sorted ?) */
			for (vx = s->uvtx; vx != NULL; vx = vx->link) {
				vtx_eserr(s, vx);		/* Estimated sampling error */
				if (vx->eserr > mx) {
					mx = vx->eserr;
					bvx = vx;
				}
			}

			/* Add the new point */
#ifdef RANDOM_PERTERB
			rerr = PERTERB_AMOUNT * sqrt(bvx->rads);
			for (e = 0; e < di; e++)
				p->p[e] = bvx->p[e] + d_rand(-rerr, rerr);
			ofps_clip_point(s, p->p);
			for (e = 0; e < di; e++)
				p->op[e] = p->p[e];
#else /* !RANDOM_PERTERB */
			for (e = 0; e < di; e++)
				p->op[e] = p->p[e] = bvx->p[e];
#endif /* !RANDOM_PERTERB */

			s->percept(s->od, p->v, p->p);
			s->np = i+1;

			/* Compute the Voronoi for it */
			node_voronoi(s, i);

			/* Add bisecting plane due to new node to all of its neighbors */
			for (k = 0; k < p->nvp; k++) {
				peq *vp = p->vp[k];				/* plane being added */
				if (vp->ix >= 0)				/* If not gamut plane */
					node_add_neigbor(s, vp);	/* add new point as neighbor */
			}
		}

		/* Add node to acceleration structure */
		ofps_add_nacc(s, p);

		if (j == 11 || i == (s->tinp-1)) {
			printf("\rAdded %d/%d",s->np,s->tinp); fflush(stdout);
			j = 0;
		}
	}
	printf("\n");
}

/* Recreate the Voronoi diagram with the current point positions */
/* (Could avoid adding fixed points again, by making a copy and restoring it ??) */
static void
ofps_redo_voronoi(
ofps *s
) {
	int j;

	/* Clear out the acceleration structure */
	ofps_reset_acc(s);

	/* Add the points in again */
	for (s->np = 0; s->np < s->tinp; s->np++) {
		node *p = &s->n[s->np];	/* Destination for point */

		/* Compute the Voronoi for it */
		node_voronoi(s, s->np);

		/* Add bisecting plane due to new node to all of its neighbors */
		for (j = 0; j < p->nvp; j++) {
			peq *vp = p->vp[j];				/* plane being added */
			if (vp->ix >= 0)				/* If not gamut plane */
				node_add_neigbor(s, vp);	/* add new point as neighbor */
		}
		/* Add node to acceleration structure */
		ofps_add_nacc(s, p);
	}
}

/* --------------------------------------------------- */
/* Optimise the seeded sample points */

static void
ofps_optimise(
ofps *s
) {
	double oshoot, ioshoot = OPT_INITIAL_OVERSHOOT;
	int i, j;

	oshoot = ioshoot;
	for (j = 0; j < s->maxits; j++) {	/* Up to maximum number of itterations */
#ifdef DUMP_PLOT
		dump_image(s, 0, 0, DUMP_VTX);		/* Device, No wait, no verticies */
#endif /* DUMP_PLOT */

		s->mxmvsq = 0.0;
		for (i = s->fnp; i < s->np; i++) {
			comp_opt(s, i, oshoot);
		}
		if (j < 10) {
			oshoot = (ioshoot - 1.0) * (10.0 - j)/10.0 + 1.0;
		} else {
			oshoot = 1.0;
		}

		/* Rebuild Voronoi surface from scratch */
		ofps_redo_voronoi(s);

		ofps_stats(s);
		printf("Maxmv = %f, stats Min = %f, Average = %f, Max = %f\n",sqrt(s->mxmvsq),s->mn,s->av,s->mx);
		printf("                Percep: Min = %f, Average = %f, Max = %f\n",s->pmn,s->pav,s->pmx);
		if (sqrt(s->mxmvsq) < STOP_TOL)
			break;
	}
}


/* --------------------------------------------------- */
/* Statistics */

static void ofps_stats(ofps *s) {
	int i, j;

	s->pmn = s->mn = 1e80;
	s->pmx = s->mx = -1e80;
	s->pav = s->av = 0.0;

	/* Figure out the stats for each sample point, and overall stats */
	for (i = 0; i < s->np; i++) {
		double mx, mn;
		double pmx, pmn;

		mx = -1e80;
		mn = 1e80;
		pmx = -1e80;
		pmn = 1e80;

		for (j = 0; j < s->n[i].nvv; j++) { 	/* For all the verticies */
			double es;
			vtx *vx = s->n[i].vv[j].p;
			if (vx->rads < mn)
				mn = vx->rads;
			if (vx->rads > mx)
				mx = vx->rads;
			es = vtx_eserr(s, vx);
			if (es < pmn)
				pmn = es;
			if (es > pmx)
				pmx = es;
		}
		for (j = 0; j < s->n[i].nvp; j++) {		/* For all the planes */
			peq *pl = s->n[i].vp[j];
			if (pl->rads < mn)
				mn = pl->rads;
		}
		mn = sqrt(mn);
		mx = sqrt(mx);
		pmn = sqrt(pmn);
		pmx = sqrt(pmx);

		if (mn < s->mn)
			s->mn = mn;
		if (mx > s->mx)
			s->mx = mx;
		s->av += mn;
		s->av += mx;
		if (pmn < s->pmn)
			s->pmn = pmn;
		if (pmx > s->pmx)
			s->pmx = pmx;
		s->pav += pmn;
		s->pav += pmx;
	}

	if (s->np > 0) {
		s->av /= (2.0 * s->np);
		s->pav /= (2.0 * s->np);
	}
}

/* --------------------------------------------------- */
/* Support accessing the list of generated sample points */

/* Reset the read index */
static void
ofps_reset(ofps *s) {
	s->rix = 0;
}

/* Read the next non-fixed point value */
/* Return nz if no more */
static int
ofps_read(ofps *s, double *p, double *f) {
	int e;

	/* Advance to next non-fixed point */
	while(s->rix < s->np && s->n[s->rix].fx)
		s->rix++;
	
	if (s->rix >= s->np)
		return 1;

	/* Return point info to caller */
	for (e = 0; e < s->di; e++) {
		if (p != NULL)
			p[e] = s->n[s->rix].p[e];
		if (f != NULL)
			f[e] = s->n[s->rix].v[e];
	}
	s->rix++;

	return 0;
}

/* --------------------------------------------------- */
/* Main object creation/destruction */

/* Destroy ourselves */
static void
ofps_del(ofps *s) {
	int i;

	/* Free our nodes */
	for (i = 0; i < s->np; i++) {
		node_free(s, &s->n[i]);
	} 
	free(s->n);

	/* gamut "node" */
	node_free(s, &s->gn);

	/* Any free planes or vertexes */
	while (s->fpeq != NULL) {
		peq *p = s->fpeq;
		s->fpeq = p->link;
		free(p);
	}
	while (s->fvtx != NULL) {
		vtx *p = s->fvtx;
		s->fvtx = p->link;
		free(p);
	}
	free (s);
}

/* Constructor */
ofps *new_ofps(
int di,					/* Dimensionality of device space */
double ilimit,			/* Ink limit (sum of device coords max) */
int tinp,				/* Total number of points to generate, including fixed */
double dadaptation,		/* Degree of adaptation to device characteristic 0.0 - 1.0 */
fxpos *fxlist,			/* List of existing fixed points (may be NULL) */
int fxno,				/* Number of existing fixes points */
void (*percept)(void *od, double *out, double *in),		/* Perceptual lookup func. */
void *od				/* context for Perceptual function */
) {
	int i;
	double imin[MXPD];	/* Ink limit - limit on min of p[] */
	double imax[MXPD];	/* Ink limit - limit on min of p[] */

	for (i = 0; i < di; i++) {
		imin[i] = 0.0;
		imax[i] = 1.0;
	}

	return new_ofps_ex(di, ilimit, imin, imax, tinp, dadaptation, fxlist, fxno, percept, od);
}

/* Extended constructor */
ofps *new_ofps_ex(
int di,					/* Dimensionality of device space */
double ilimit,			/* Ink limit (sum of device coords max) */
double *imin,			/* Ink limit - limit on min of p[], usually >= 0.0 */
double *imax,			/* Ink limit - limit on min of p[], usually <= 1.0  */
int tinp,				/* Total number of points to generate, including fixed */
double dadaptation,		/* Degree of adaptation to device characteristic 0.0 - 1.0 */
fxpos *fxlist,			/* List of existing fixed points (may be NULL) */
int fxno,				/* Number of existing fixes points */
void (*percept)(void *od, double *out, double *in),		/* Perceptual lookup func. */
void *od				/* context for Perceptual function */
) {
	int i;
	ofps *s;

	if ((s = (ofps *)calloc(sizeof(ofps), 1)) == NULL)
		error ("ofps: malloc failed");

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	if (di > MXPD)
		error ("ofps: Can't handle di %d",di);

	if (dadaptation < 0.0)
		dadaptation = 0.0;
	else if (dadaptation > 1.0)
		dadaptation = 1.0;

	s->di = di;

	if (tinp < fxno)	/* Make sure we return at least the fixed points */
		tinp = fxno;

	s->tinp = tinp;		/* Target total number of points */
	s->ilimit = ilimit;
	for (i = 0; i < di; i++) {
		s->imin[i] = imin[i];
		s->imax[i] = imax[i];
	}

	/* If no perceptual function given, use default */
	if (percept == NULL) {
		s->percept = default_ofps_to_percept;
		s->od = s;
		s->padapt = 0.0;
	} else {
		s->percept = percept;
		s->od = od;
		s->padapt = dadaptation;	/* Degree of perceptual adaptation */
	}

	if (s->padapt > 0.0)		/* Because this doesn't work yet */
		s->maxits = AMAXITS;
	else
		s->maxits = MAXITS;

	/* Init method pointers */
	s->reset = ofps_reset;
	s->read  = ofps_read;
	s->stats = ofps_stats;
	s->del   = ofps_del;

	/* Allocate the space for the target number of points */
	if ((s->n = (node *)calloc(sizeof(node), s->tinp)) == NULL)
		error ("ofps: malloc failed on sample nodes");
	s->np = s->fnp = 0;

	/* Setup gamut boundary surfaces */
	ofps_binit(s);

	/* Setup acceleration structures */
	ofps_init_acc(s);

	/* Setup the fixed points */
	ofps_add_fixed(s, fxlist, fxno);

	if (fxno > 0) {
		ofps_stats(s);
		printf("Stats after fixed points: Min = %f, Average = %f, Max = %f\n",s->mn,s->av,s->mx);
		printf("                  Percep: Min = %f, Average = %f, Max = %f\n",s->pmn,s->pav,s->pmx);
	}

	if (tinp > fxno) {
		/* Create the moveable points */
		ofps_seed(s);

		ofps_stats(s);
		printf("Stats after seeding points: Min = %f, Average = %f, Max = %f\n",s->mn,s->av,s->mx);
		printf("                    Percep: Min = %f, Average = %f, Max = %f\n",s->pmn,s->pav,s->pmx);
	
		/* Optimise the points */
		ofps_optimise(s);
	
		/* Print some stats */
		ofps_stats(s);
		printf("Stats after opt: Min = %f, Average = %f, Max = %f\n",s->mn,s->av,s->mx);
		printf("         Percep: Min = %f, Average = %f, Max = %f\n",s->pmn,s->pav,s->pmx);
	}

	ofps_reset(s);		/* Reset read index */

	return s;
}

/* =================================================== */

#ifdef STANDALONE_TEST

/* Graphics Gems curve */
static double gcurve(double vv, double g) {
	if (g >= 0.0) {
		vv = vv/(g - g * vv + 1.0);
	} else {
		vv = (vv - g * vv)/(1.0 - g * vv);
	}
	return vv;
}

#ifdef NEVER
static void sa_percept(void *od, double *out, double *in) {
	double lab[3];
	
	clu->dev_to_rLab(clu, lab, in);

	out[0] = lab[0];
//	out[1] = (lab[1]+100.0)/2.0;
	out[1] = (lab[2]+100.0)/2.0;
}
#else

static void sa_percept(void *od, double *p, double *d) {

#ifdef NEVER
	/* Default two curves with some interaction */
	p[0] = 100.0 * gcurve(d[0], -4.5);
	p[1] = 100.0 * gcurve(d[1], 2.8);
	p[1] = 0.8 * p[1] + 0.2 * p[0];

#else
	int e;
	for (e = 0; e < 2; e++) {
		double tt = d[e];
		/* Two slopes with a sharp turnover in X */
		if (e == 0) {
			if (tt < 0.5)
				tt = tt * 0.3/0.5;
			else
				tt = 0.3 + ((tt-0.5) * 0.7/0.5);
		}
		p[e] = tt * 100.0;
	}
#endif
}
#endif


int
main(argc,argv)
int argc;
char *argv[];
{
	int npoints = 50;
	ofps *s;
	long stime,ttime;
	fxpos fx[4];		/* Any fixed points */

	error_program = argv[0];

	printf("Standalone test of ofps, argument is number of points, default %d\n",npoints);

	if (argc > 1)
		npoints = atoi(argv[1]);

	fx[0].p[0] = 0.1;
	fx[0].p[1] = 0.2;

	fx[1].p[0] = 0.5;
	fx[1].p[1] = 0.3;

	fx[2].p[0] = 0.3;
	fx[2].p[1] = 0.6;

	fx[3].p[0] = 0.9;
	fx[3].p[1] = 0.2;

	/* Create the required points */
	stime = clock();
	s = new_ofps(2, 1.5, npoints, SA_ADAPT, fx, 0, sa_percept, (void *)NULL);

	ttime = clock() - stime;
	printf("Execution time = %f seconds\n",ttime/(double)CLOCKS_PER_SEC);

#ifdef DUMP_PLOT
//	printf("Perceptual plot:\n");
//	dump_image(s, 1, DO_WAIT, DUMP_VTX);

	printf("Device plot:\n");
	dump_image(s, 0, DO_WAIT, DUMP_VTX);
#endif /* DUMP_PLOT */

	s->del(s);

	return 0;
}

#ifdef NEVER
/* Basic printf type error() and warning() routines */
#ifdef	__STDC__
void
error(char *fmt, ...)
#else
void
error(va_alist) 
va_dcl
#endif
{
	va_list args;
#ifndef	__STDC__
	char *fmt;
#endif

	fprintf(stderr,"ofps: Error - ");
#ifdef	__STDC__
	va_start(args, fmt);
#else
	va_start(args);
	fmt = va_arg(args, char *);
#endif
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	fflush(stdout);
	exit (-1);
}
#endif /* NEVER */
#endif /* STANDALONE_TEST */


#if defined(DEBUG) || defined(DUMP_PLOT)

/* Dump the current point positions to a plot window file */
static void
dump_image(ofps *s, int pcp, int dwt, int dvx) {
	int i, j;
	double minx, miny, maxx, maxy;
	static double *x1a = NULL;		/* Current sample locations */
	static double *y1a = NULL;
	static double *x2a = NULL;		/* Previous sample locations */
	static double *y2a = NULL;
	static int _n3 = 0;				/* Current Voronoi verticies */
	static double *x3a = NULL;
	static double *y3a = NULL;
	int n3 = 0;

	if (pcp != 0) {	/* Perceptual range */
		minx = 0.0;	/* Assume */
		maxx = 100.0;
		miny = 0.0;
		maxy = 100.0;
	} else {
		minx = 0.0;	/* Assume */
		miny = 0.0;
		maxx = 1.0;
		maxy = 1.0;
	}
	
	if (x1a == NULL) {
		if ((x1a = (double *)malloc(s->tinp * sizeof(double))) == NULL)
			error ("ofps: malloc failed");
		if ((y1a = (double *)malloc(s->tinp * sizeof(double))) == NULL)
			error ("ofps: malloc failed");
		if ((x2a = (double *)malloc(s->tinp * sizeof(double))) == NULL)
			error ("ofps: malloc failed");
		if ((y2a = (double *)malloc(s->tinp * sizeof(double))) == NULL)
			error ("ofps: malloc failed");
	}

	/* Add sample node location */
	for (i = 0; i < s->np; i++) {
		node *p = &s->n[i];
		
		if (pcp != 0) {
			x2a[i] = x1a[i] = p->v[0];
			y2a[i] = y1a[i] = p->v[1];
		} else {
			x1a[i] = p->p[0];
			y1a[i] = p->p[1];
			x2a[i] = p->op[0];
			y2a[i] = p->op[1];
		}
	}

	if (dvx) {
		vtx *vx;

		if (x3a == NULL) {		/* Initial allocation */
			_n3 = s->np * 4;
			if ((x3a = (double *)malloc(_n3 * sizeof(double))) == NULL)
				error ("ofps: malloc failed");
			if ((y3a = (double *)malloc(_n3 * sizeof(double))) == NULL)
				error ("ofps: malloc failed");
		}

		/* Add Voronoi verticies */
		for (vx = s->uvtx; vx != NULL; vx = vx->link) {

			if (n3 >= _n3) {		/* need more space */
				_n3 *= 2;
				if ((x3a = (double *)realloc(x3a, _n3 * sizeof(double))) == NULL)
					error ("ofps: malloc failed");
				if ((y3a = (double *)realloc(y3a, _n3 * sizeof(double))) == NULL)
					error ("ofps: malloc failed");
			}
			x3a[n3]   = vx->p[0];
			y3a[n3++] = vx->p[1];
		}
	}

	/* Plot the vectors */
	do_plot_vec(minx, maxx, miny, maxy, 
				x1a, y1a, x2a, y2a, s->np, dwt, x3a, y3a, dvx ? n3 : 0);
}

#endif /* DEBUG || DUMP_PLOT */















