
/* 
 * gamut
 *
 * Gamut support routines.
 *
 * Author:  Graeme W. Gill
 * Date:    9/3/2000
 * Version: 1.00
 *
 * Copyright 2000 - 2006 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* TTBD:

*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "icc.h"
#include "gamut.h"
#include "cgats.h"
#include "numlib.h"
#include "sort.h"			/* ../h sort macro */
#include "counters.h"		/* ../h counter macros */
#include "xlist.h"			/* ../h expandable list macros */

#define COLORED_VRML

#define DO_TWOPASS			/* Second pass with adjustment based on first pass */

#define USE_FAKE_SEED		/* Use fake seed points */
#define FAKE_SEED_SIZE 0.1	/* Usually 0.1 */
#define TRIANG_TOL 1e-10	/* Triangulation tollerance, usually 1e-10 */

#undef TEST_CONVEX_HULL		/* Use pure convex hull, not log hull */

#undef DEBUG_TRIANG			/* Enable detailed triangulation debugging */
#undef DEBUG_TRIANG_VRML	/* Create debug.wrl for each step of triangulation */
							/* (Only on second pass if #define DO_TWOPASS) */
#undef DEBUG_TRIANG_VRML_STEP	/* Wait for return after each step */

#undef DEBUG_SPLIT_VRML	/* Create debug.wrl for each step of triangle plane split */

#undef TEST_LOOKUP
#undef TEST_NEAREST
#undef TEST_NOSORT			/* Turn off sorted insersion of verticies */

#undef SHOW_BUCKETS			/* Show vertex buckets as surface */
#undef SHOW_SPHERE			/* Show surface on sphere */
#undef SHOW_HULL_PNTS		/* Show log() length convex hull points */

#undef ASSERTS				/* Do internal checking */

/* These routines support:

   representing the 3D gamut boundary of a device or image as
   radial surface height, described by a triangular poligon hull.

   Interogate the surface to find the point lying on the hull
   in the same radial direction as the query point.

   Interogate the surface to find the point lying on the hull
   that is the closest to the query point.

   Save the gamut as a vrml format viewable file.
   Save the gamut as a CGATS format .gam file.

 */

/* TTBD:
 *
 *  Would be nice to take the exact colorspace specification
 * (ie. Lab vs. Jab + viewing conditions), and store them
 *  in the .gam file, so that a warning can be issued if
 *  the gamut colorspace is a mismatch in icclink, or to be able
 *  to translate the verticies into the correct colorespace.
 *
 *  Add interface to transform all the nodes, while
 *  keeping structure (For use within gamut map, or
 *  transforming existing gamut from Lab to Jab etc.)
 *
 *  Add inteface to fetch the triangle information ?
 *
 *  Need to cleanup error handling. We just exit() at the moment.
 *
 *  Replace BSP tree optmisation with ball tree, to speedup
 *  radial, nearest, and vector search ?
 *
 *  The log surface stuff is a compromise, that ends up with
 *  some dings and nicks, and a not fully detailed/smooth surface.
 *  The fundamental limitation is the use of the Delaunay triangulation
 *  criteria, and the triangulation algorithm dependence on it for consistency.
 *  Want to switch to triangulation algorithm that doesn't
 *  depend on this, and can triangulate concave objects,
 *  so that something like alpha-shapes criteria can be
 *  used to filter out non surface points. Inserting verticies
 *  from largest radius to smallest seems to do the right thing
 *  with cusp ridges, and this property needs to be retained.
 *
 */


#ifndef M_PI
#define M_PI (3.1415926535897932384626433832795)
#endif

static void triangulate(gamut *s);
static void del_gamut(gamut *s);
static void expand_gamut(gamut *s, double in[3]);
static double getsres(gamut *s);
static int getisjab(gamut *s);
static int compatible(gamut *s, gamut *t);
static int nrawverts(gamut *s);
static int getrawvert(gamut *s, double pos[3], int ix);
static int nraw0verts(gamut *s);
static int getraw0vert(gamut *s, double pos[3], int ix);
static int nverts(gamut *s);
static int getvert(gamut *s, double *rad, double pos[3], int ix);
static void startnexttri(gamut *s);
static int getnexttri(gamut *s, int v[3]);
static double volume(gamut *s);
static int write_trans_vrml(gamut *s, char *filename, int doaxes, int docusps,
	void (*transform)(void *cntx, double out[3], double in[3]), void *cntx);
static int write_vrml(gamut *s, char *filename, int doaxes, int docusps);
static int write_gam(gamut *s, char *filename);
static int read_gam(gamut *s, char *filename);
static double radial(gamut *s, double out[3], double in[3]);
static double nradial(gamut *s, double out[3], double in[3]);
static void nearest(gamut *s, double out[3], double in[3]);
static void setwb(gamut *s, double *wp, double *bp);
static int getwb(gamut *s, double *cswp, double *csbp, double *gawp, double *gabp);
static void setcusps(gamut *s, int flag, double in[3]);
static int getcusps(gamut *s, double cusps[6][3]);
static int compute_vector_isect(gamut *s, double *p1, double *p2, double *min, double *max, double *mint, double *maxt);
static double log_scale(double ss);
static int intersect(gamut *s, gamut *s1, gamut *s2);
static int vect_intersect(gamut *s, double *rvp, double *ip, double *p1, double *p2, gtri *t);

/* in isecvol.c: */
extern double isect_volume(gamut *s1, gamut *s2);

/* ------------------------------------ */

/* Hue directions in degrees for Lab and Jab */
/* Must be in increasing order */
double gam_hues[2][7] = {
	{
		/* Lab */
		36.0,			/* Red */
		101.0,			/* Yellow */
		149.0,			/* Green */
		225.0,			/* Cyan */
		300.0,			/* Blue */
		337.0,			/* Magenta */
		36.0 + 360.0	/* Red */
	},
	{
		/* Jab */
		28.0,			/* Red */
		101.0,			/* Yellow */
		148.0,			/* Green */
		211.0,			/* Cyan */
		269.0,			/* Blue */
		346.0,			/* Magenta */
		28.0 + 360.0	/* Red */
	}
};


/* ------------------------------------ */
static
gvert *new_gvert(
gamut *s,
gquad *p,			/* Parent quad (may be NULL) */
int i,				/* Intended node in quad */
int f,				/* Flag value to be OR'ed */
double pp[3],		/* Point in xyz rectangular coordinates, absolute */
double rr[3],		/* Radial coordinates */
double lrr0,		/* log scaled rr[0] */
double sp[3],		/* Point mapped to surface of unit sphere, relative to center */
double ch[3]		/* Point mapped for convex hull testing, relative to center */
) {
	gvert *v;

	if (s->doingfake == 0 && s->ul != NULL) {	/* There is an unused one available */
		v = s->ul;
		s->ul = v->ul;

	} else {				/* Allocate a new one */

		if (s->nv >= s->na) {	/* We need a new slot in the list */
			if (s->na == 0) {
				s->na = 5;
				if ((s->verts = (gvert **)malloc(s->na * sizeof(gvert *))) == NULL) {
					fprintf(stderr,"gamut: malloc failed on %d gvert pointer\n",s->na);
					exit (-1);
				}
			} else {
				s->na *= 2;
				if ((s->verts = (gvert **)realloc(s->verts, s->na * sizeof(gvert *))) == NULL) {
					fprintf(stderr,"gamut: realloc failed on %d gvert pointer\n",s->na);
					exit (-1);
				}
			}
		}
	
		if ((v = (gvert *)calloc(1, sizeof(gvert))) == NULL) {
			fprintf(stderr,"gamut: malloc failed on gvert object\n");
			exit (-1);
		}
		s->verts[s->nv] = v;
		v->n = s->nv++;
	}
	v->tag = 1;

	if  (p != NULL) {
		v->w = 0.5 * p->w;
		v->h = 0.5 * p->h;

		v->hc = p->hc;
		if (i & 1)
			v->hc += 0.5 * v->w;
		else
			v->hc -= 0.5 * v->w;
	
		v->vc = p->vc;
		if (i & 2)
			v->vc += 0.5 * v->h;
		else
			v->vc -= 0.5 * v->h;
	} else {
		v->w = 0.0;
		v->h = 0.0;
		v->hc = 0.0;
		v->vc = 0.0;
	}
		
	v->f = GVERT_NONE | f;
	v->ul = NULL;
	v->rc = 1;

	v->p[0] = pp[0];
	v->p[1] = pp[1];
	v->p[2] = pp[2];
	v->r[0] = rr[0];
	v->r[1] = rr[1];
	v->r[2] = rr[2];
	v->lr0  = lrr0;
	v->sp[0] = sp[0];
	v->sp[1] = sp[1];
	v->sp[2] = sp[2];
	v->ch[0] = ch[0];
	v->ch[1] = ch[1];
	v->ch[2] = ch[2];

	return v;
}

/* Set the size of gvert angular segment */
static
void set_gvert(
gvert *v,		/* gvert to set */
gquad *p,		/* Parent quad */
int i			/* Intended node in quad */
) {
	v->w = 0.5 * p->w;
	v->h = 0.5 * p->h;

	v->hc = p->hc;
	if (i & 1)
		v->hc += 0.5 * v->w;
	else
		v->hc -= 0.5 * v->w;

	v->vc = p->vc;
	if (i & 2)
		v->vc += 0.5 * v->h;
	else
		v->vc -= 0.5 * v->h;
}

/* Increment the reference count on a gvert */
static gvert *inc_gvert(gamut *s, gvert *v) {
	if (v == NULL)
		return v;
	v->rc++;
//printf("~1 incremented count on 0x%x to %d\n",v,v->rc);
	return v;
}

/* Decrement the reference count on a gvert */
/* Copes with NULL gvert. */
/* If the reference count goes to 0, place the */
/* vert on the unused list. */
static void dec_gvert(gamut *s, gvert *v) {
	if (v == NULL)
		return;

	v->rc--;
//printf("~1 decremended count on 0x%x to %d\n",v,v->rc);

#ifdef ASSERTS
	if (v->tag != 1)
		error("Assert: doing decremented on gquad node");

	if (v->rc < 0)
		error("Assert: decremented gvert ref count too far");
#endif
	
	if (v->rc <= 0) {	/* Add it to the unused list */
		memset((void *)v, 0, sizeof(gvert));
		v->ul = s->ul;
		s->ul = v;
	}
}

/* Delete all the gverts */
static void del_gverts(gamut *s) {
	int i;
	for (i = 0; i < s->nv; i++) {
		free(s->verts[i]);
	}
	if (s->verts != NULL) {
		free(s->verts);
		s->na = 0;
		s->nv = 0;
	}
}

/* ------------------------------------ */

static
gquad *new_gquad(
gquad *p,			/* Parent quad */
int i				/* Intended node in quad */
) {
	gquad *q;
	if ((q = (gquad *)calloc(1, sizeof(gquad))) == NULL) {
		fprintf(stderr,"gamut: calloc failed on gquad object\n");
		exit (-1);
	}
	q->tag = 2;

	q->w = 0.5 * p->w;
	q->h = 0.5 * p->h;

	q->hc = p->hc;
	if (i & 1)
		q->hc += 0.5 * q->w;
	else
		q->hc -= 0.5 * q->w;

	q->vc = p->vc;
	if (i & 2)
		q->vc += 0.5 * q->h;
	else
		q->vc -= 0.5 * q->h;
		
	return q;
}

/* Same as above, but create with explicit size */
static
gquad *new_gquad2(
double l,		/* Left border */
double r,		/* Right border */
double b,		/* Top border */
double t		/* Bottom border */
) {
	gquad *q;
	if ((q = (gquad *)calloc(1, sizeof(gquad))) == NULL) {
		fprintf(stderr,"gamut: calloc failed on gquad object\n");
		exit (-1);
	}
	q->tag = 2;
	q->w = r - l;
	q->h = t - b;
	q->hc = (l + r) * 0.5;
	q->vc = (t + b) * 0.5;
	return q;
}

static void
del_gquad(gquad *q) {
	int i;

	if (q == NULL)
		return;
	for (i = 0; i < 4; i++) {
		gnode *n = (gnode *)q->qt[i][0];
		if (n != NULL && n->tag == 2)
			del_gquad((gquad *)n);		/* Recurse */
	}
	free(q);
}

/* Helper functions */

/* Given a gquad and a location, decide which quandrant we're in */
static int gquad_quadrant(gquad *q, double p[3]) {
	int i;

#ifdef ASSERTS
	if (p[1] < (q->hc - q->w * 0.5 - 1e-10)
			 || p[1] > (q->hc + q->w * 0.5 + 1e-10)
			 || p[2] < (q->vc - q->h * 0.5 - 1e-10)
	 || p[2] > (q->vc + q->h * 0.5 + 1e-10)) {
		fprintf(stderr,"error! point doesn't fall into bucket chosen for it!!!\n");
		fprintf(stderr,"point x: %f < %f\n", p[1], (q->hc - q->w * 0.5));
		fprintf(stderr,"point x: %f > %f\n", p[1], (q->hc + q->w * 0.5));
		fprintf(stderr,"point y: %f < %f\n", p[2], (q->vc - q->h * 0.5));
		fprintf(stderr,"point y: %f > %f\n", p[2], (q->vc + q->h * 0.5));
		exit(-1);
	}
#endif

	i = 0;
	if (p[1] >= q->hc)
		i |= 1;
	if (p[2] >= q->vc)
		i |= 2;

	return i;
}

/* ------------------------------------ */
/* Allocate a BSP decision node structure */
gbspn *new_gbspn(void) {
	gbspn *t;
	static int n = 0; 		/* Serial number */
	if ((t = (gbspn *) calloc(1, sizeof(gbspn))) == NULL) {
		fprintf(stderr,"gamut: malloc failed - bspn node\n");
		exit(-1);
	}
	t->tag = 1;		/* bspn node */
	t->n = n++;

	return t;
}

/* Delete a BSP decision node struture */
void del_gbspn(gbspn *t) {
	free(t);
}

/* ------------------------------------ */
/* Allocate a BSP tree triangle list node structure */
gbspl *new_gbspl(
int nt,			/* Number of triangles in the list */
gtri **t		/* List of triangles to copy into structure */
) {
	gbspl *l;
	int i;
	static int n = 0; 		/* Serial number */
	if ((l = (gbspl *) calloc(1, sizeof(gbspl) + nt * sizeof(gtri *))) == NULL) {
		fprintf(stderr,"gamut: malloc failed - bspl triangle tree node\n");
		exit(-1);
	}
	l->tag = 3;		/* bspl node */
	l->n = n++;
	l->nt = nt;
	for (i = 0; i < nt; i++)
		l->t[i] = t[i];

	return l;
}

/* Delete a BSP tree triangle list structure */
void del_gbspl(gbspl *l) {
	free(l);
}

/* ------------------------------------ */
/* Allocate a triangle structure */
gtri *new_gtri(void) {
	gtri *t;
	static int n = 0; 		/* Serial number */
	if ((t = (gtri *) calloc(1, sizeof(gtri))) == NULL) {
		fprintf(stderr,"gamut: malloc failed - gamut surface triangle\n");
		exit(-1);
	}
	t->tag = 2;		/* Triangle */
	t->n = n++;

	return t;
}

/* Delete a triangle struture */
void del_gtri(gtri *t) {
	free(t);
}

/* ------------------------------------ */
/* Allocate an edge structure */
gedge *new_gedge(void) {
	gedge *t;
	static int n = 0; 		/* Serial number */
	if ((t = (gedge *) calloc(1, sizeof(gedge))) == NULL) {
		fprintf(stderr,"gamut: malloc failed - triangle edge\n");
		exit(-1);
	}
	t->n = n++;
	return t;
}

/* Delete an edge struture */
void del_gedge(gedge *t) {
	free(t);
}

/* ------------------------------------ */
	
/* Create a standard gamut map */
gamut *new_gamut(
double sres,			/* Resolution (in rect coord units) of surface triangles */
						/* 0.0 = default */
int isJab				/* Flag indicating Jab space */
) {
	gamut *s;

#ifdef ASSERTS
	fprintf(stderr,">>>>>>> ASSERTS ARE COMPILED INTO GAMUT.C <<<<<<<\n");
#endif /* ASSERTS */

	if ((s = (gamut *)calloc(1, sizeof(gamut))) == NULL) {
		fprintf(stderr,"gamut: calloc failed on gamut object\n");
		exit (-1);
	}

	if (sres <= 0.0)
		sres = 10.0;	/* default */
	if (sres > 15.0)	/* Anything less is very poor */
		sres = 15.0;
	s->sres = sres;

	if (isJab != 0)
		s->isJab = 1;

	/* Center point for radial values, surface creation etc. */
	/* To compare two gamuts using radial values, their cent must */
	/* be the same. */
	s->cent[0] = 50.0;
	s->cent[1] = 0.0;
	s->cent[2] = 0.0;

	/* Create top level quadtree nodes */
	s->tl = new_gquad2(-M_PI, 0.0, -M_PI/2.0, M_PI/2.0);	/* Left one */
	s->tr = new_gquad2(0.0, M_PI, -M_PI/2.0, M_PI/2.0);	/* Right one */

	INIT_LIST(s->tris);			/* Init triangle list */
	INIT_LIST(s->edges);		/* Init edge list (?) */
	s->read_inited = 0;
	s->lu_inited = 0;
	s->ne_inited = 0;
	s->cswbset = 0;
	s->gawbset = 0;

	/* Setup methods */
	s->del         = del_gamut;
	s->expand      = expand_gamut;
	s->getsres     = getsres;
	s->getisjab    = getisjab;
	s->compatible  = compatible;
	s->nrawverts   = nrawverts;
	s->getrawvert  = getrawvert;
	s->nraw0verts  = nraw0verts;
	s->getraw0vert = getraw0vert;
	s->nverts      = nverts;
	s->getvert     = getvert;
	s->startnexttri = startnexttri;
	s->getnexttri = getnexttri;
	s->getvert     = getvert;
	s->volume      = volume;
	s->intersect   = intersect;
	s->radial      = radial;
	s->nradial     = nradial;
	s->nearest     = nearest;
	s->vector_isect = compute_vector_isect;
	s->setwb       = setwb;
	s->getwb       = getwb;
	s->setcusps    = setcusps;
	s->getcusps    = getcusps;
	s->write_vrml  = write_vrml;
	s->write_trans_vrml = write_trans_vrml;
	s->write_gam   = write_gam;
	s->read_gam    = read_gam;

	return s;
}

static void del_gnn(gnn *p);
static void del_gbsp(gbsp *n);

/* Free and clear the triangulation structures, */
/* and clear the triangulation vertex flags. */
static void del_triang(gamut *s) {
	int i;
	gtri *tp;		/* Triangle pointer */
	gedge *ep;

	if (s->tris != NULL) {
		tp = s->tris; 					/* Delete all the triangles */
		FOR_ALL_ITEMS(gtri, tp) {
			DEL_LINK(s->tris, tp);
		} END_FOR_ALL_ITEMS(tp);

		INIT_LIST(s->tris);			/* Init triangle list */
	}
	
	if (s->edges != NULL) {
		ep = s->edges; 					/* Delete all the edges */
		FOR_ALL_ITEMS(gedge, ep) {
			DEL_LINK(s->edges, ep);
		} END_FOR_ALL_ITEMS(ep);
		INIT_LIST(s->edges);		/* Init edge list */
	}
	
	/* Free radial lookup acceleration structures */
	if (s->lutree != NULL) {
		del_gbsp(s->lutree);
		s->lutree = NULL;
	}
	s->lu_inited = 0;

	if (s->nns != NULL) {
		del_gnn(s->nns);
		s->nns = NULL;
	}
	s->ne_inited = 0;

	/* Reset the vertex flags triangulation changes */
	for (i = 0; i < s->nv; i++) {
		s->verts[i]->f &= ~GVERT_TRI;
		s->verts[i]->f &= ~GVERT_INSIDE;
	}
}

static void
del_gamut(gamut *s) {
	del_gquad(s->tl);
	del_gquad(s->tr);

	del_triang(s);
	del_gverts(s);

	free(s);
}

/* ===================================================== */
/* Segmented maxima filter code */
/* ===================================================== */

/*
	This implementation has two twists on the usual
    segmented maxima filtering:

	* Rather than using a uniform radial grid, it
      uses an adaptive, quadtree structure, so that
      rather than having a constant angular resolution
      for the segments, the sements are chosen so as to
      aproximately correspond to a constant surface detail
      level. [ A subtle problem with this approach is that
      some points will be discarded early on, that wouldn't
      be discarded later, when the quadtree is finer. A hack
      is to run the points throught twice.]

    * Rather than keep the single sample with the
      largest radius in each radial segment,
      four samples are kept, each being the largest
      in a different direction. This helps avoid
      "knicks" in edges, where a non edge sample
      displaces an edge sample within a segment.

 */

/* Helper function that returns nz if v1 should replace v2 */
static int smreplace(
gamut *s,
int d,			/* Slot number */
gvert *v1,		/* Candidate vertex */
gvert *v2		/* Existing vertex */
) {
	double xx, w[3], c[3];
	int j;
	if (v2 == NULL)
		return 1;

	/* Filter out any points that are almost identical */
	/* This can cause numerical problems in the triangle BSP tree creation. */ 
	for (xx = 0.0, j = 0; j < 3; j++) {
		double tt = v1->p[j] - v2->p[j];
		xx += tt * tt;
	}
	if (xx < (1e-4 * 1e-4))
		return 0;

	c[0] = s->cent[0];
	c[1] = s->cent[1];
	c[2] = s->cent[2];

	/* Set L, a & b weighting depending on the slot */
	switch(d) {
		case 1:
			w[0] = 0.5;
			w[1] = 1.0;
			w[2] = 1.0;
			break;
		case 2:
			w[0] = 2.0;
			w[1] = 1.0;
			w[2] = 1.0;
			break;
		case 3:
			w[0] = 4.0;
			w[1] = 1.0;
			w[2] = 1.0;
			break;
		case 4:
			w[0] = 0.0;
			w[1] = 1.0;
			w[2] = 0.1;
			break;
		case 5:
			w[0] = 0.0;
			w[1] = 0.1;
			w[2] = 1.0;
			break;
		default:
			w[0] = 1.0;
			w[1] = 1.0;
			w[2] = 1.0;
			break;
	}
	w[0] *= w[0];		/* Because we're weighting the squares */
	w[1] *= w[1];
	w[2] *= w[2];
	return ( w[0] * (v1->p[0] - c[0]) * (v1->p[0] - c[0])
	       + w[1] * (v1->p[1] - c[1]) * (v1->p[1] - c[1])
	       + w[2] * (v1->p[2] - c[2]) * (v1->p[2] - c[2]))
	     > ( w[0] * (v2->p[0] - c[0]) * (v2->p[0] - c[0])
	       + w[1] * (v2->p[1] - c[1]) * (v2->p[1] - c[1])
	       + w[2] * (v2->p[2] - c[2]) * (v2->p[2] - c[2]));
}


/* Expand the gamut by adding a point */
static void expand_gamut(
gamut *s,
double pp[3]		/* rectangular coordinate of point */
) {
	gnode *n;		/* Current node */
	gvert *nv, *ov;	/* new vertex, old vertex */
	gquad *q;		/* Parent quad */
	int i;			/* Sub element within quad */
	int k;			/* Index of direction slot */
	double rr[3];	/* Radial coordinate version of pp[] */
	double sp[3];	/* Unit shere mapped version of pp[] relative to center */
	double ch[3];	/* Convex hull testing mapped version of pp[] relative to center */
	double lrr0;	/* log scaled rr[0] */
	double hang, vang;	/* Critical angles for this points depth */
	double aa;
	int j;

	if (s->tris != NULL || s->read_inited || s->lu_inited || s->ne_inited) {
		fprintf(stderr,"Can't add points to gamut now!\n");
		exit(-1);
	}

	if (s->doingfake == 0)
		s->cu_inited = 0;		/* Invalidate cust info */

	/* Convert to radial coords */
	gamut_rect2radial(s, rr, pp);

	if (rr[0] < 1e-6) 		/* Ignore a point right at the center */
		return;

	/* Figure log scaled radius */
	lrr0 = log_scale(rr[0]);

	/* Compute unit shere mapped location */
	aa = 1.0/rr[0];				/* Adjustment to put in on unit sphere */
	for (j = 0; j < 3; j++)
		sp[j] = (pp[j] - s->cent[j]) * aa;

	/* Compute hull testing mapped version */
	for (j = 0; j < 3; j++)
		ch[j] = sp[j] * lrr0;

	/* compute twice angle resolution required (will compare to parent size) */
	hang = pow(rr[0], 1.01) * fabs(cos(rr[2]));
	if (hang < 1e-9)
		hang = 1e-9;
	hang = 4.0 * s->sres/hang;
	vang = 4.0 * s->sres/pow(rr[0], 1.01);

//printf("~1 Point at %f %f %f, radial %f %f %f\n", pp[0], pp[1], pp[2], rr[0], rr[1], rr[2]);
//printf("~1     shere at %f %f %f, log %f %f %f, vhang %f %f\n", sp[0], sp[1], sp[2], ch[0], ch[1], ch[2], vang, hang);

	/* If nofilter flag is set, add all points as verticies */
	if (s->nofilter) {

		/* Filter out any points that are almost identical. */
		/* This can cause numerical problems in the triangle BSP tree creation. */ 
		for (i = 0; i < s->nv; i++) {
			double xx;
	
			for (xx = 0.0, j = 0; j < 3; j++) {
				double tt = pp[j] - s->verts[i]->p[j];
				xx += tt * tt;
			}
			if (xx < (1e-4 * 1e-4)) {
				return;			/* Ignore duplicate */
			}
		}
		
		/* Create a vertex for the point we're possibly adding */
		new_gvert(s, NULL, 0, GVERT_SET | (s->doingfake ? GVERT_FAKE : 0), pp, rr, lrr0, sp, ch);

		return;
	}
	/* else filter using adaptive segmented maxima */

	/* Start by looking at the top level quads */
	if (rr[1] >= 0.0) {
		q = s->tr;
	} else {
		q = s->tl;
	}
	n = (gnode *)q;
//printf("~1 Starting with quad 0x%x, width %f, height %f\n",q, q->w, q->h);

	/* Now recurse until we have a virtex at the right location and depth */
	for (;;) {
		/* Recurse into quad node n */
		q = (gquad *)n;					/* Parent node */
		i = gquad_quadrant(q, rr);		/* Quadrand of parent */
		n = q->qt[i][0];				/* Child node in quadrant */

//printf("~1 Current quad 0x%x, width %f, height %f\n",q, q->w, q->h);
//printf("~1 Current child in quadrant %d, node 0x%x, type %d\n", i, n, n != NULL ? n->tag : 0);

		/* If we're at the right depth to create a vertex, break out of decent loop. */

		if (n == NULL) {			/* Create new node */
			if (q->w <= hang && q->h <= vang) {
//printf("~1 We're at the right depth to add vertex\n");
				break;
			}
			/* Else create a new quad */
			n = (gnode *)new_gquad(q, i);
			q->qt[i][0] = n;
//printf("~1 Empty child node not deep enough, creating new quad node 0x%x\n",n);

		/* If we've found verticies at this node */
		} else if (n->tag == 1) {
			int j;
			gquad *qq;				/* New child quad */
			gvert *vv[NSLOTS];		/* Existing verticies at this level */

			if (q->w <= hang && q->h <= vang) {
//printf("~1 We're at the right depth to replace vertex\n");
				break;
			}
//printf("~1 deepening verticies\n");

			for (k = 0; k < NSLOTS; k++)
				vv[k] = (gvert *)q->qt[i][k];	/* Save pointers to current verticies */
			
//printf("~1 existing verticies are 0x%x, 0x%x, 0x%x, 0x%x\n", vv[0], vv[1], vv[2], vv[3]);

			/* make a quad to replace the current verticies */
			qq = new_gquad(q, i);
			n = (gnode *)qq;
			q->qt[i][0] = n;
			for (k = 1; k < NSLOTS; k++)
				q->qt[i][k] = NULL;

//printf("~1 added quad 0x%x to quadrant %d\n",i,q);

			/* Distribute verticies that were here, into new quad */
			for (j = 0; j < NSLOTS; j++) {			/* For all existing verticies */

				if (vv[j] == NULL)
					continue;

//printf("~1 re-distributing verticy 0x%x\n",vv[j]);
				i = gquad_quadrant(qq, vv[j]->r);	/* Quadrant for existing vertex */

				set_gvert(vv[j], qq, i);			/* Update vertex node location */

				nv = vv[j];
				for (k = 0; nv != NULL && k < NSLOTS; k++) {	/* For direction slot */
					ov = (gvert *)qq->qt[i][k];
					if (smreplace(s, k, nv, ov)) {
						if (k == 0) {			/* Track points that are in k == 0 direction */
							if (ov != NULL && ov->k0 > 0)
								ov->k0--;
							nv->k0++;
						}
#ifndef NEVER
						qq->qt[i][k] = (gnode *)inc_gvert(s, nv);
						dec_gvert(s, ov);
#else
						/* Use slots for best, 2nd best, etc */
						qq->qt[i][k] = (gnode *)inc_gvert(s, nv);
						dec_gvert(s, nv);
						nv = ov;
#endif
//printf("Node 0x%x rc %d at %f %f %f is replacing\n",nv, nv->rc, nv->p[0], nv->p[1], nv->p[2]);
//if (ov != NULL) printf(" replacing node 0x%x rc %d at %f %f %f\n",ov, ov->rc, ov->p[0], ov->p[1], ov->p[2]);
//else printf(" NULL\n");
					}
				}
				dec_gvert(s, nv);
			}
		}
		/* Else it's a quad, and we will decend into it */

	}	/* keep decending until we find right depth */

//printf("~1 Got parent quad 0x%x, quadrant %d, vertex 0x%x\n", q, i, n);

	/* Create a vertex for the point we're possibly adding */
	nv = new_gvert(s, q, i, GVERT_SET, pp, rr, lrr0, sp, ch);

	/* Replace any existing gverts with this one */
	for (k = 0; k < NSLOTS; k++) {		/* For direction slot */
		ov = (gvert *)q->qt[i][k];
		if (smreplace(s, k, nv, ov)) {
			if (k == 0) {			/* Track points that are in k == 0 direction */
				if (ov != NULL && ov->k0 > 0)
					ov->k0--;
				nv->k0++;
			}
#ifndef NEVER
			q->qt[i][k] = (gnode *)inc_gvert(s, nv);
			dec_gvert(s, ov);
#else
			/* Use slots for best, 2nd best, etc */
			q->qt[i][k] = (gnode *)inc_gvert(s, nv);
			dec_gvert(s, nv);
			nv = ov;
#endif
		}
	}
	dec_gvert(s, nv);	 /* Make sure it's reclaimed if wasn't used */

//printf("~1 Point is done\n\n");
}

/* ------------------------------------ */
/* Initialise this gamut with the intersection of the */
/* the two given gamuts. Return NZ on error. */
/* Return 1 if gamuts are not compatible */
/* (We assume that the gamut is currently empty) */
static int intersect(gamut *s, gamut *sa, gamut *sb) {
	int i, j, k;
	gamut *s1, *s2;

	if (sa->compatible(sa, sb) == 0)
		return 1;

	if IS_LIST_EMPTY(sa->tris)
		triangulate(sa);
	if IS_LIST_EMPTY(sb->tris)
		triangulate(sb);

	s->isJab = sa->isJab;
	for (j = 0; j < 3; j++)
		s->cent[j] = sa->cent[j];

	/* Clear some flags */
	s->cswbset = 0;
	s->cswbset = 0;
	s->dcuspixs = 0;

	/* Don't filter the points (gives a more accurate result ?) */
	s->nofilter = 1;

	/* Add each source gamuts verticies that lie within */
	/* the other gamut */
	for (k = 0; k < 2; k++) {
		gtri *tp1, *tp2;		/* Triangle pointers */

		if (k == 0) {
			s1 = sa;
			s2 = sb;
		} else {
			s1 = sb;
			s2 = sa;
		}
		for (i = 0; i < s1->nv; i++) {
			double pl;
	
			if (!(s1->verts[i]->f & GVERT_TRI))
				continue;
	
			pl = s2->nradial(s2, NULL, s1->verts[i]->p);
			if (pl <= (1.0 + 1e-9)) {
				expand_gamut(s, s1->verts[i]->p);
				s1->verts[i]->f &= ~GVERT_ISOS;
			} else {
				s1->verts[i]->f |= GVERT_ISOS;
			}
		}

		/* Now find the edges that intersect the other gamut */
		tp1 = s1->tris;
		FOR_ALL_ITEMS(gtri, tp1) {

			for (j = 0; j < 3; j++) {
				if ((tp1->e[j]->v[0]->f ^ tp1->e[j]->v[1]->f) & GVERT_ISOS) {

					/* Exhaustive search of other triangles */
					tp2 = s2->tris;
					FOR_ALL_ITEMS(gtri, tp2) {
						double pv;
						double tt[3];

						/* Do a min/max intersection elimination test */
						for (i = 0; i < 3; i++) {
							if (tp2->mix[1][i] < tp1->mix[0][i]
							 || tp2->mix[0][i] > tp1->mix[1][i])
								break;			/* min/max don't overlap */
						}
						if (i < 3)
							continue;			/* Skip this triangle, it can't intersect */

						if (vect_intersect(s1, &pv, tt, tp1->e[j]->v[0]->p, tp1->e[j]->v[1]->p, tp2)
						 && pv >= (0.0 - 1e-10) && pv <= (1.0 + 1e-10)) {
							expand_gamut(s, tt);
						}
					} END_FOR_ALL_ITEMS(tp2);
				}
			}

		} END_FOR_ALL_ITEMS(tp1);
	}

	s->nofilter = 0;

	return 0;
}

/* ------------------------------------ */
/* Locate the vertices most likely to correspond to the */
/* primary and secondary colors (cusps) */
/*
 * Notes:
 *
 *	To better support gamuts of devices with more than 4 colorants,
 *  it may be necessary to add another flag type that expands
 *  the cusps to lie on the actual gamut surface, as for
 *  some devices this lies outside the pure colorant combinations.
 */ 

static void setcusps(gamut *s, int flag, double in[3]) {
	int i, j;

	if (flag == 0) {	/* Reset */
		for (j = 0; j < 6; j++) {
			s->cusps[j][0] = 0.0;		/* Marker values */
			s->cusps[j][1] = 0.0;
			s->cusps[j][2] = 0.0;
		}
		s->dcuspixs = 0;
		s->cu_inited = 0;
		return;

	} else if (flag == 2) {	/* Finalize */

		if (s->dcuspixs > 0) {
			double JCh[3];
			double hues[6];
			int r, br = 0.0;
			double berr;

			/* Figure out where to put the ones we got */
			for (j = 0; j < 6; j++) {	/* Compute the hues */
				icmLab2LCh(JCh, s->dcusps[j]);
				hues[j] = JCh[2];
			}

			/* Sort them into hue order */
			for (j = 0; j < 5; j++) {
				for (i = j+1; i < 6; i++) {
					if (hues[j] > hues[i]) {
						double tt;
						tt = hues[j]; hues[j] = hues[i]; hues[i] = tt;
						tt = s->dcusps[j][0]; s->dcusps[j][0] = s->dcusps[i][0]; s->dcusps[i][0] = tt;
						tt = s->dcusps[j][1]; s->dcusps[j][1] = s->dcusps[i][1]; s->dcusps[i][1] = tt;
						tt = s->dcusps[j][2]; s->dcusps[j][2] = s->dcusps[i][2]; s->dcusps[i][2] = tt;
					}
				}
			}

			/* Figure out which is the best match by rotation */
			berr = 1e6;
			for (r = 0; r < 6; r++) {
				double terr = 0.0;
				for (j = 0; j < 6; j++) {
					double tt;

					tt = fabs(gam_hues[s->isJab][j] - hues[(j + r) % 6]);
					if (tt > 180.0)
						tt = 360.0 - tt;
					terr += tt;
				}
				if (terr < berr) {
					br = r;
					berr = terr;
				}
			}
			/* Place them at that rotation */
			for (j = 0; j < 6; j++) {	/* Compute the hues */
//printf("~1 placing hue %f ix %d into hue %f ix %d\n", hues[(j + br) % 6],(j + br) % 6, gam_hues[s->isJab][j] ,j);

				s->cusps[j][0] = s->dcusps[(j + br) % 6][0];
				s->cusps[j][1] = s->dcusps[(j + br) % 6][1];
				s->cusps[j][2] = s->dcusps[(j + br) % 6][2];
			}
		}

		/* Check we've got a cusp */
		for (j = 0; j < 6; j++) {
			if (s->cusps[j][0] == 0.0
			 && s->cusps[j][1] == 0.0
			 && s->cusps[j][2] == 0.0) {
//printf("~1 not enough cusps were found = slot %d hue %f is empty\n",j,gam_hues[s->isJab][j]);
				s->cu_inited = 0;
				return;			/* Not all have been set */
			}
		}
		s->cu_inited = 1;
		return;

	} else if (flag == 3) {	/* Definite 1/6 cusp */

		if (s->dcuspixs >= 6) {
//printf("~1 too many cusp values added\n");
			return;						/* Error - extra cusp ignored */
		}
		s->dcusps[s->dcuspixs][0] = in[0];
		s->dcusps[s->dcuspixs][1] = in[1];
		s->dcusps[s->dcuspixs++][2] = in[2];

	} else {	/* Consider another point as a cusp point */
		double JCh[3];
		int bj = 0, sbj = 0;
		double bh = 1e6, sbh = 1e6;
		double ns, es;

		icmLab2LCh(JCh, in);

		/* See which hue it is closest and 2nd closet to cusp hue. */
		for (j = 0; j < 6; j++) {
			double tt;

			tt = fabs(gam_hues[s->isJab][j] - JCh[2]);
			if (tt > 180.0)
				tt = 360.0 - tt;

			if (tt < bh) {
				if (bh < sbh) {
					sbh = bh;
					sbj = bj;
				}
				bh = tt;
				bj = j;
			} else if (tt < sbh) {
				sbh = tt;
				sbj = j;
			}
		}

		/* Compute distance of existing and new */
		es = s->cusps[bj][1] * s->cusps[bj][1] + s->cusps[bj][2] * s->cusps[bj][2];
		ns = in[1] * in[1] + in[2] * in[2];
//printf("~1 chroma dist of chexisting %f, new %f\n",es,ns);
		if (ns > es) {
//printf("~1 New closest\n");
			s->cusps[bj][0] = in[0];
			s->cusps[bj][1] = in[1];
			s->cusps[bj][2] = in[2];

		} else {	/* If 2nd best has no entry, use this to fill it */
			if (s->cusps[sbj][0] == 0.0
			 && s->cusps[sbj][1] == 0.0
			 && s->cusps[sbj][2] == 0.0) {
//printf("~1 Fill with 2nd closest\n");
				s->cusps[sbj][0] = in[0];
				s->cusps[sbj][1] = in[1];
				s->cusps[sbj][2] = in[2];
			}
		}
	}
}


/* Get the cusp values for red, yellow, green, cyan, blue & magenta */
/* Return nz if there are no cusps available */
static int getcusps(
gamut *s,
double cusps[6][3]
) {
	int i, j;

	if (s->cu_inited == 0)
		return 1;

	for (i = 0; i < 6; i++)
		for (j = 0; j < 3; j++)
			cusps[i][j] = s->cusps[i][j];

	return 0;
}

/* ===================================================== */
/* Triangulation code */
/* ===================================================== */

/* Given three points, compute the normalised plane equation */
/* of a surface through them. */
/* Return non-zero on error */
static int plane_equation(
double *eq,		/* Return equation parameters */
double *p0,		/* The three points */
double *p1,
double *p2
) {
	double ll, v1[3], v2[3];

	/* Compute vectors along edges */
	v1[0] = p1[0] - p0[0];
	v1[1] = p1[1] - p0[1];
	v1[2] = p1[2] - p0[2];

	v2[0] = p2[0] - p0[0];
	v2[1] = p2[1] - p0[1];
	v2[2] = p2[2] - p0[2];

	/* Compute cross products v1 x v2, which will be the normal */
	eq[0] = v1[1] * v2[2] - v1[2] * v2[1];
	eq[1] = v1[2] * v2[0] - v1[0] * v2[2];
	eq[2] = v1[0] * v2[1] - v1[1] * v2[0];

	/* Normalise the equation */
	ll = sqrt(eq[0] * eq[0] + eq[1] * eq[1] + eq[2] * eq[2]);
	if (ll < 1e-10) {
		return 1;
	}
	eq[0] /= ll;
	eq[1] /= ll;
	eq[2] /= ll;

	/* Compute the plane equation constant */
	eq[3] = - (eq[0] * p0[0])
	        - (eq[1] * p0[1])
	        - (eq[2] * p0[2]);

#ifdef NEVER
	/* Veritify the plane equation */
	{
		double c;
		c = eq[0] * p0[0]
	      + eq[1] * p0[1]
	      + eq[2] * p0[2]
		  + eq[3];
		if (fabs(c) > 1e-10) {
			printf("Plane equation check 0 failed by %f\n",c);
		}
		c = eq[0] * p1[0]
	      + eq[1] * p1[1]
	      + eq[2] * p1[2]
		  + eq[3];
		if (fabs(c) > 1e-10) {
			printf("Plane equation check 1 failed by %f\n",c);
		}
		c = eq[0] * p2[0]
	      + eq[1] * p2[1]
	      + eq[2] * p2[2]
		  + eq[3];
		if (fabs(c) > 1e-10) {
			printf("Plane equation check 2 failed by %f\n",c);
		}
	}
#endif /* NEVER */
	return 0;
}

/* Compute the log surface plane equation for the triangle */
/* and other triangle attributes. (Doesn't depend on edge info.) */
void
comptriattr(
gamut *s,
gtri *t
) {
	static double v0[3] = {0.0, 0.0, 0.0};

	/* Compute the plane equation for the absolute triangle. */
	/* This is used for testing if a point is inside the gamut hull. */
	plane_equation(t->pe, t->v[0]->p, t->v[1]->p, t->v[2]->p);

	/* Compute the plane equation for the triangle */
	/* based on the log compressed convex hull verticies. */
	/* This is used for convex hull construction. */
	plane_equation(t->che, t->v[0]->ch, t->v[1]->ch, t->v[2]->ch);

	/* Compute the plane equation for the triangle */
	/* mapped to the surface of the sphere */
	/* This can be used for point in triangle testing ?? */
	plane_equation(t->spe, t->v[0]->sp, t->v[1]->sp, t->v[2]->sp);

	/* Compute the plane equations of the spherical mapped vertex */
	/* values with regard to the center of the sphere, so that */
	/* a point in triangle test can be performed, and baricentric, */
	/* coordinates can be computed. */
	plane_equation(t->ee[0], v0, t->v[1]->sp, t->v[2]->sp);
	plane_equation(t->ee[1], v0, t->v[2]->sp, t->v[0]->sp);
	plane_equation(t->ee[2], v0, t->v[0]->sp, t->v[1]->sp);

#ifdef NEVER // ???
#ifdef ASSERTS
	{
	int j;
	double tt[3];	/* Triangle test point */
	double ds;
	for (j = 0; j < 3; j++) {
		tt[j] = (t->v[0]->p[j] + t->v[1]->p[j] + t->v[2]->p[j])/3.0;
		tt[j] -= s->cent[j];			/* Make it center relative */
	}
	for (j = 0; j < 3; j++) {
		ds = t->ee[j][0] * tt[0]	/* Point we know is inside */
		    + t->ee[j][1] * tt[1]
	        + t->ee[j][2] * tt[2]
		    + t->ee[j][3];
		if (ds > 1e-8)
			break;		/* Not within triangle */
	}
	if (j < 3) {
		fprintf(stderr,"Assert: point expected to be within triangle %d (vx %d %d %d) is not\n",
		        t->n, t->v[0]->n, t->v[1]->n, t->v[2]->n);
		fprintf(stderr,"Known point is %f, expect -ve\n",ds);
		exit(-1);
		}
	}
#endif /* ASSERTS */
#endif /* NEVER */

}

/* By using the pow() or log() of the radial distance, */
/* blended with a sphere surface, we try and strike a compromise */
/* between a pure convex hull surface, and a pure Delaunay triangulation */
/* the latter which would show dings and nicks from points */
/* that fall withing the "real" gamut. */
static double log_scale(double ss) {
	double aa;

#ifdef TEST_CONVEX_HULL
	return ss;
#else
#ifdef NEVER		/* (Not using this version) */
	aa = (2.0 + ss)/3.0;	/* Blend with sphere */
	aa = log(aa);			/* Allow for concave slope */
	if (aa < 0.0)			/* but constrain to be +ve */
		aa = 0.0;
#else				/* (Using simpler version) */
	aa = 20.0 * pow(ss, 0.25);
#endif
#endif	/* TEST_CONVEX_HULL */

	return aa;
}

/* Comput r[], lr0, sp[] and ch[] from p[] for all verticies */
/* (Note that lr0 will be the first cut, log_scale() value) */ 
static void
compute_vertex_coords(
gamut *s
) {
	int i, j;

	for (i = 0; i < s->nv; i++) {
		gamut_rect2radial(s, s->verts[i]->r, s->verts[i]->p);

		if (s->verts[i]->r[0] < 1e-6) { 		/* Ignore a point right at the center */
			s->verts[i]->lr0 = 0.0;
			for (j = 0; j < 3; j++) {
				s->verts[i]->sp[j] = 0.0;
				s->verts[i]->ch[j] = 0.0;
			}
		} else {
			double aa;

			/* Figure log scaled radius */
			s->verts[i]->lr0 = log_scale(s->verts[i]->r[0]);

			/* Compute unit shere mapped location */
			aa = 1.0/s->verts[i]->r[0];		/* Adjustment to put in on unit sphere */
			for (j = 0; j < 3; j++)
				s->verts[i]->sp[j] = (s->verts[i]->p[j] - s->cent[j]) * aa;
		
			/* Compute hull testing mapped version */
			for (j = 0; j < 3; j++)
				s->verts[i]->ch[j] = s->verts[i]->p[j] * s->verts[i]->lr0;
		}
	}
}

/* Sort the verticies from maximum radius */
static void sort_verticies(
gamut *s
) {
	int i;

#ifndef TEST_NOSORT

	/* Sort them */
//#define 	HEAP_COMPARE(A,B) (A->r[0] > B->r[0])
#define 	HEAP_COMPARE(A,B) (A->lr0 > B->lr0)
			HEAPSORT(gvert *, s->verts, s->nv)
#undef HEAP_COMPARE

#endif /* !TEST_NOSORT */

	/* Renumber them */
	for (i = 0; i < s->nv; i++) {
		s->verts[i]->n = i;
	}
}

/* Number just the verticies that have been set, */
/* and those that have been used in the convex hull */
static void renumber_verticies(
gamut *s
) {
	int i, j;

	for (j = i = 0; i < s->nv; i++) {
		if (!(s->verts[i]->f & GVERT_SET))
			continue;

		s->verts[i]->sn = j;
		j++;
	}
	s->nsv = j;

	for (j = i = 0; i < s->nv; i++) {
		if (!(s->verts[i]->f & GVERT_TRI))
			continue;

		s->verts[i]->tn = j;
		j++;
	}
	s->ntv = j;
}

/* Given a list of target convex hull locations, */
/* return pointers to the closest verticies in those directions */
/* (Used to setup initial convex hull) */
void closest_verticies(
gamut *s,			/* Gamut with convex hull points to search */
gvert **tvs,		/* Return pointers to best verticies */
double **tgts,		/* Target vectors */		
int nn				/* Number of points */
) {
	int i, j, m;
	double *bestd;

	if ((bestd = (double *) malloc(nn * sizeof(double))) == NULL) {
		fprintf(stderr,"gamut: malloc failed - closest verticies best distance\n");
		exit(-1);
	}

	/* We are using a brute force search */
	for (m = 0; m < nn; m++) {
		bestd[m] = -1e38;
		tvs[m] = NULL;
	}

	for (m = 0; m < nn; m++) {			/* For all target vectors */
		double ss;

		/* Compute dot product of target vector and verticy vectors */
		for (i = 0; i < s->nv; i++) {			/* For all the search points */

			if (!(s->verts[i]->f & GVERT_SET))
				continue;

			if (!(s->verts[i]->f & GVERT_FAKE))
				continue;

			/* If this vector is already been used for a target, skip it */
			for (j = 0; j < m; j++) {
				if (tvs[j] == s->verts[i])
					break;
			}
			if (j < m)
				continue;

			for (ss = 0.0, j = 0; j < 3; j++) {
				double t1, t2;
				t1 = tgts[m][j];
				t2 = s->verts[i]->ch[j];
				ss += t1 * t2;
			}

			if (ss > bestd[m]) {	/* Found better candidate */
				bestd[m] = ss;
				tvs[m] = s->verts[i];
			}
		}
	}

	/* Sanity check */
	for (m = 0; m < nn; m++) {
		if (tvs[m] == NULL) {
			fprintf(stderr,"gamut: internal error - enough seed verticies not found\n");
			exit(-1);
		}
	}

#ifdef DEBUG_TRIANG
	printf("Returning closest verticies %d %d %d %d\n",
	        tvs[0]->n, tvs[1]->n, tvs[2]->n, tvs[3]->n);
	for (m = 0; m < nn; m++) {
		printf("Vert %d ch = %f %f %f\n", m, tvs[m]->ch[0], tvs[m]->ch[1], tvs[m]->ch[2]);
	}
#endif
	free(bestd);
}

#ifdef ASSERTS

/* Diagnpostic aid */
/* Check that the triangulation adjacenty info is OK */
static void check_triangulation(gamut *s, int final) {
	int i, j;
	gtri *tp;		/* Triangle pointer */
	gedge *ep;		/* Edge pointer */
	int failed = 0;

	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {

		/* Check verticies for duplication */
		for (i = 0; i < 3; i++) {
			for (j = i+1; j < 3; j++) {
				if (tp->v[i] == tp->v[j]) {
					failed = 1;
	printf("Validation failed - duplicate verticies:\n");
	printf("Triangle %d, has verticies %d %d %d\n", tp->n, tp->v[0]->n, tp->v[1]->n, tp->v[2]->n);
	fflush(stdout);
				}
			}
		}

		/* Check edges for duplication */
		for (i = 0; i < 3; i++) {
			for (j = i+1; j < 3; j++) {
				if (tp->e[i] == tp->e[j]) {
					failed = 1;
	printf("Validation failed - duplicate connectivity:\n");
	printf("Triangle %d, has verticies %d %d %d\n", tp->n, tp->v[0]->n, tp->v[1]->n, tp->v[2]->n);
	printf("Triangle %d, has edges %d %d %d\n", tp->n, tp->e[0]->n, tp->e[1]->n, tp->e[2]->n);
	fflush(stdout);
				}
			}
		}

		/* Check connectivity */
		for (i = 0; i < 3; i++) {
			gtri *t1, *t2;
			gedge *e;
			int ei1, ei2;
			int tei;						/* Edges index for this triangle [0..1] */

			e = tp->e[i];					/* The edge in question */
			tei = tp->ei[i];
			ei1 = e->ti[tei];				/* Edges record of edge index withing this triangle */

			/* Check that the edges reconing of what index edge it is */
			/* for this triangle is correct */
			if (ei1 != i) {
				failed = 1;
	printf("Validation failed - triangle edge index doesn't match record withing edge:\n");
	printf("Triangle %d, edge index %d edge %d has record %d\n", tp->n, i, e->n, ei1);
	fflush(stdout);
			}

			/* Check that the edges pointer to the triangle is this triangle */
			if (tp != e->t[tei]) {
				failed = 1;
	printf("Validation failed - edge doesn't point back to triangle:\n");
	printf("Triangle %d, edge index %d is edge %d\n",tp->n, i, e->n);
	printf("Edge     %d, triangle index %d is triangle %d\n", e->n, tei, e->t[tei]->n);
	printf("Edge     %d, triangle index %d is triangle %d\n", e->n, tei^1, e->t[tei^1]->n);
	fflush(stdout);
			}

			/* Check the verticies for this edge match edge record */
			if ((e->v[0] != tp->v[i] || e->v[1] != tp->v[(i+1) % 3])
			 && (e->v[1] != tp->v[i] || e->v[0] != tp->v[(i+1) % 3])) {
				failed = 1;
	printf("Validation failed - edge doesn't have same verticies as triangle expects:\n");
	printf("Triangle %d, has verticies %d %d\n", tp->n, tp->v[i]->n, tp->v[(i+1) % 3]->n);
	printf("Edge     %d, has verticies %d %d\n", e->n, e->v[0]->n, e->v[1]->n);
	fflush(stdout);
			}

			t2 = e->t[tei ^ 1];		/* The other triangle */
			ei2 = e->ti[tei ^ 1];		/* Edges index number withing triangle t2 */

			if (t2 == tp) {
				failed = 1;
	printf("Validation failed - connects to itself:\n");
	printf("Triangle %d, has edges %d %d %d\n", tp->n, tp->e[0]->n, tp->e[1]->n, tp->e[2]->n);
	fflush(stdout);
			}

			/* Check that the connection is reflective */
			if (e != t2->e[ei2]) {
				failed = 1;
	printf("Validation failed - connectivity not reflected:\n");
	printf("Triangle %d, edge index %d points to edge %d\n",tp->n, i, e->n);
	printf("Triangle %d, edge index %d points to edge %d\n",t2->n, ei2, t2->e[ei2]->n);
	fflush(stdout);
			}
		}
	} END_FOR_ALL_ITEMS(tp);
	if (failed) {
		exit(-1);
	}

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Check that every point is part of a triangle and edge */

	for (i = 0; i < s->nv; i++) {		/* Reset the assert flag */
		gvert *v = s->verts[i];

		v->as = 0;

		/* Check out the flags */
		if (!(v->f & GVERT_SET)) {
			if ((v->f & GVERT_TRI)
			 || (v->f & GVERT_INSIDE)) {
	printf("Validation failed - vertex %d has strange flags 0x%x\n",i, v->f);
	fflush(stdout);
				failed = 1;
			}
		} else {
			if ((v->f & GVERT_TRI) && (v->f & GVERT_INSIDE)) {
	printf("Validation failed - vertex %d has strange flags 0x%x\n",i, v->f);
	fflush(stdout);
				failed = 1;
			}
		}
	}

	FOR_ALL_ITEMS(gtri, tp) {
		for (i = 0; i < 3; i++)
			tp->v[i]->as |= 1;		/* Vertex is in a triangle */
	} END_FOR_ALL_ITEMS(tp);

	ep = s->edges; 
	FOR_ALL_ITEMS(gedge, ep) {
		ep->v[0]->as |= 2;			/* Vertex is in an edge */
		ep->v[1]->as |= 2;
		ep->as = 0;					/* Reset the assert flag */
	} END_FOR_ALL_ITEMS(ep);

	for (i = 0; i < s->nv; i++) {
		if (s->verts[i]->f & GVERT_TRI) {
			if ((s->verts[i]->as & 1) == 0) {
	printf("Validation failed - vertex %d is not in any triangles\n",i);
	fflush(stdout);
				failed = 1;
			}
			if ((s->verts[i]->as & 2) == 0) {
	printf("Validation failed - vertex %d is not in any edge\n",i);
	fflush(stdout);
				failed = 1;
			}
		}
	}


	/* - - - - - - - - - - - - - - - - - - - - - - */
	/* Check that every edge is part of a triangle */

	/* as flag in triangle was reset above */
	FOR_ALL_ITEMS(gtri, tp) {
		for (i = 0; i < 3; i++)
			tp->e[i]->as |= 1;		/* Mark edge used in triangle */
	} END_FOR_ALL_ITEMS(tp);

	ep = s->edges; 
	FOR_ALL_ITEMS(gedge, ep) {
		if (ep->as != 1) {
	printf("Validation failed - edge %d is not in any triangle\n",ep->n);
	fflush(stdout);
			failed = 1;
		}
	} END_FOR_ALL_ITEMS(ep);

	if (failed) {
		exit(-1);
	}
}

#endif /* ASSERTS */

/* -------------------------------------- */
/* Add a face to the hit list, if it is not a duplicate. */
static void add_to_hit_list(
gamut *s, 
gtri **hlp,		/* Hit list */
gtri *cf		/* Face to be added (triangle verts 0, 1) */
) {
	gtri *tp;		/* Triangle pointer */
	gvert *c0 = cf->v[0];
	gvert *c1 = cf->v[1];

//printf("Adding face to hit list %d: %d %d\n",
//cf->n, cf->v[0]->n, cf->v[1]->n);

	tp = *hlp; 
	/* Search current faces in hit list */
	FOR_ALL_ITEMS(gtri, tp) {
		gvert *v0 = tp->v[0];
		gvert *v1 = tp->v[1];
		if ((c0 == v0 && c1 == v1)			/* Same face from other side */
		 || (c0 == v1 && c1 == v0)) {
			/* Duplicate found */
//printf("Duplicate found %d: %d %d\n",
//tp->n, tp->v[0]->n, tp->v[1]->n);
			DEL_LINK(*hlp, tp);		/* Delete from the hit list */

			/* Check face is common */
			if (cf->e[0] != tp->e[0]) {
				fprintf(stderr,"gamut: internal error - face match inconsistency\n");
				exit(-1);
			}
			/* Delete edge */
			DEL_LINK(s->edges, cf->e[0]);
			del_gedge(cf->e[0]);

			/* Delete the two faces */
			del_gtri(tp);
			del_gtri(cf);
			return;
		}
	} END_FOR_ALL_ITEMS(tp);

	/* Safe to add it to face hit list */
	/* This removes triangle from triangles list ? */
	ADD_ITEM_TO_BOT(*hlp, cf);
//printf("Face added\n");
}

/* Add a triangles faces to the hit list. */
static void add_tri_to_hit_list(
gamut *s, 
gtri **hlp,		/* Hit list */
gtri *tp		/* Triangle faces to be added */
) {
	int j;
	gtri *t1, *t2;

	/* In case some verticies disapear below the log surface, */
	/* and don't remain part of the triangulation, we mark them off. */
	for (j = 0; j < 3 ; j++) {
		tp->v[j]->f &= ~GVERT_TRI;
		tp->v[j]->f |= GVERT_INSIDE;
	}

	/* Decompose the triangle into three faces, each face being stored */
	/* into a triangle created on the hit list, using verticices 0, 1. */
	/* The edges adjacency info remains valid for the three faces, */
	/* as does the edge plane equation. */
	DEL_LINK(s->tris, tp);		/* Delete it from the triangulation list */
	t1 = new_gtri();
	t1->v[0] = tp->v[1];		/* Duplicate with rotated faces */
	t1->v[1] = tp->v[2];
	t1->e[0] = tp->e[1];		/* Edge adjacency for this edge */
	t1->ei[0] = tp->ei[1];		/* Edge index of this triangle */
	t1->e[0]->t[t1->ei[0]] = t1; /* Fixup reverse adjacency for valid edge */
	t1->e[0]->ti[t1->ei[0]] = 0; /* Rotated index of new triangles edge */
	t1->e[1] = t1->e[2] = NULL;	/* be safe */
	for (j = 0; j < 4; j++)		/* Copy edge plane equation */
		t1->ee[2][j] = tp->ee[0][j];

	t2 = new_gtri();
	t2->v[0] = tp->v[2];		/* Duplicate with rotated faces */
	t2->v[1] = tp->v[0];
	t2->e[0] = tp->e[2];		/* Edge adjacency for this edge */
	t2->ei[0] = tp->ei[2];		/* Edge index of this triangle */
	t2->e[0]->t[t2->ei[0]] = t2; /* Fixup reverse adjacency for valid edge */
	t2->e[0]->ti[t2->ei[0]] = 0; /* Rotated index of new triangles edge */
	t2->e[1] = t2->e[2] = NULL;	/* be safe */
	for (j = 0; j < 4; j++)		/* Copy edge plane equation */
		t2->ee[2][j] = tp->ee[1][j];

	tp->e[1] = tp->e[2] = NULL;	/* be safe */
	add_to_hit_list(s, hlp, tp);	/* Add edge 0 to hit list as is */
	add_to_hit_list(s, hlp, t1);	/* Add edge 1 to hit list */
	add_to_hit_list(s, hlp, t2);	/* Add edge 2 to hit list */
}

#ifdef DEBUG_TRIANG
	typedef struct {
		int tix[3];		/* Triangle indexes */
		int type;		/* 0 = hit, 1 = added */
	} tidxs;
#endif

/* Insert a vertex into the triangulation */
static void insert_vertex(
gamut *s, 
gvert *v		/* Vertex to insert */
) {
	gtri *tp, *tp2;		/* Triangle pointers */
	gtri *hl;			/* Triangle face hit list (polygon faces) */
	double tol = TRIANG_TOL;
	int hit = 0;		/* Vertex expands hull flag */
#ifdef DEBUG_TRIANG
	int intri = 0;		/* Vertex landed in a triangle */
	XLIST(tidxs, hittris)
	tidxs xxs;

	XLIST_INIT(tidxs, &hittris);
#endif

#ifdef DEBUG_TRIANG
	printf("Adding vertex to triang %d: %f %f %f\n", v->n, v->p[0], v->p[1], v->p[2]);
#endif

	/* First we search the current triangles, and convert */
	/* any trianges that are visible from the new point, */
	/* into a list of faces stored on the face */
	/* hit list. */
	/* We are using a brute force search, which will make the */
	/* algorithm speed proportional to n^2. For better performance */
	/* with a large number of verticies, an acceleration structure */
	/* should be used to speed circumradius hit detection. */
	v->f &= ~GVERT_INSIDE;	/* Reset flags */
	v->f &= ~GVERT_TRI;
	INIT_LIST(hl);
	hit = 0;
	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		double c;

		/* Check the depth out compared to this triangle log plane equation */
		c = tp->che[0] * v->ch[0]
	      + tp->che[1] * v->ch[1]
	      + tp->che[2] * v->ch[2]
		  + tp->che[3];

		/* If vertex is above the log hull surface, add triangle to the hit list. */
		if (c < 0.0) {
#ifdef DEBUG_TRIANG
			int j;
#endif
			hit = 1;

#ifdef DEBUG_TRIANG
			printf("Got a hit on triangle %d: %d %d %d\n",
			           tp->n, tp->v[0]->n, tp->v[1]->n, tp->v[2]->n);
#endif

#ifdef DEBUG_TRIANG
			for (j = 0; j < 3; j++) {
				double ds;
				ds = tp->ee[j][0] * v->ch[0]
				   + tp->ee[j][1] * v->ch[1]
		   	    + tp->ee[j][2] * v->ch[2]
				   + tp->ee[j][3];
				if (ds > tol) {
					printf("Vertex is not in triangle by %e\n",ds);
					break;
				}
			}
			if (j >= 3)
				intri = 1;		/* Landed in this triangle */

			xxs.tix[0] = tp->v[0]->n, xxs.tix[1] = tp->v[1]->n, xxs.tix[2] = tp->v[2]->n;
			xxs.type = 0;
			XLIST_ADD(&hittris, xxs)
#endif
			add_tri_to_hit_list(s, &hl, tp);
		}
	} END_FOR_ALL_ITEMS(tp);
	
	if (hit == 0) {

//printf("No hits - must be inside the log hull\n");
		v->f |= GVERT_INSIDE;	/* This point is inside the log hull */
		v->f &= ~GVERT_TRI;
	} else {
		int changed = 1;

#ifdef DEBUG_TRIANG
if (!intri) printf("~1 ###### vertex didn't land in any triangle! ########\n");
#endif

//printf("Checking out hit polygon edges:\n");
		/* Now we must make a pass though the hit list, checking that each */
		/* hit list face will make a correctly oriented, non-sliver triangle */
		/* when joined to the vertex. */
		/* Do this check until there are no changes */
		for (;changed != 0 ;) {
			tp = hl; 
			changed = 0;
			FOR_ALL_ITEMS(gtri, tp) {
				/* Check which side of the edge our vertex is */
				double ds;
				ds = tp->ee[2][0] * v->ch[0]
				   + tp->ee[2][1] * v->ch[1]
			       + tp->ee[2][2] * v->ch[2]
				   + tp->ee[2][3];
//printf("Vertex margin to edge = %e\n",ds);
				/* If vertex is not to the right of this edge by tol */
				/* add associated triangle to the hit list. */
				if (ds > -tol) {
					gtri *xtp;
//printf("~1 ###### vertex on wrong side by %e - expand hit list  ######\n",ds);
					if (tp->e[0]->t[0] != tp)
						xtp = tp->e[0]->t[0];
					else
						xtp = tp->e[0]->t[1];
//printf("Got a hit on triangle %d: %d %d %d\n", xtp->n, xtp->v[0]->n, xtp->v[1]->n, xtp->v[2]->n);

#ifdef DEBUG_TRIANG
					xxs.tix[0] = xtp->v[0]->n, xxs.tix[1] = xtp->v[1]->n, xxs.tix[2] = xtp->v[2]->n;
					xxs.type = 1;
					XLIST_ADD(&hittris, xxs)
#endif

					add_tri_to_hit_list(s, &hl, xtp);
					changed = 1;
					break;
				}
			} END_FOR_ALL_ITEMS(tp);
		}

#ifdef DEBUG_TRIANG_VRML
		write_diag_vrml(s, v->ch, hittris.no, hittris.list, hl);
#endif

//printf("About to turn polygon faces into triangles\n");
		/* Turn all the faces that made it to the */
		/* hit list, into triangles using the new vertex. */
		tp = hl; 
		FOR_ALL_ITEMS(gtri, tp) {
			tp->v[2] = v;				/* Add third vertex to face to make triangle */
			comptriattr(s, tp);			/* Compute triangle attributes */

			/* Find the new adjacent triangles from the triangles being formed, */
			/* to maintain edge adjacency information. */
			/* Do only one edge at a time, since each adjacency */
			/* will be visited twice. */
			tp2 = hl; 
			FOR_ALL_ITEMS(gtri, tp2) {
				if (tp2->v[0] == tp->v[1]) {	/* Found 1/2 tp/tp2 edge adjacency */
					gedge *e;
					e = new_gedge();
					ADD_ITEM_TO_BOT(s->edges, e);	/* Append to edge list */
					tp->e[1] = e;			/* Point to edge */
					tp->ei[1] = 0;			/* edges 0th triangle */
					e->t[0] = tp;			/* triangles 1st edge */
					e->ti[0] = 1;			/* triangles 1st edge */
					tp2->e[2] = e;			/* Point to edge */
					tp2->ei[2] = 1;			/* edges 1st triangle */
					e->t[1] = tp2;			/* Triangles 2nd edge */
					e->ti[1] = 2;			/* Triangles 2nd edge */
					e->v[0] = v;			/* Add the two verticies */
					e->v[1] = tp->v[1];
				}
			} END_FOR_ALL_ITEMS(tp2);

//printf("~1 Creating new triangle %d: %d %d %d\n", tp->n, tp->v[0]->n, tp->v[1]->n, tp->v[2]->n);
		} END_FOR_ALL_ITEMS(tp);
		
#ifdef DEBUG_TRIANG_VRML
		tp = hl; 
		hittris.no = 0;
		FOR_ALL_ITEMS(gtri, tp) {
			xxs.tix[0] = tp->v[0]->n, xxs.tix[1] = tp->v[1]->n, xxs.tix[2] = tp->v[2]->n;
			xxs.type = 2;
			XLIST_ADD(&hittris, xxs)
		} END_FOR_ALL_ITEMS(tp);
		write_diag_vrml(s, v->ch, hittris.no, hittris.list, NULL);
#ifdef DEBUG_TRIANG_VRML_STEP
#ifdef DO_TWOPASS
		if (s->pass > 0)
#endif	/* DO_TWOPASS */
		getchar();
#endif
#endif

		/* Move them to the triangulation. */
		tp = hl; 
		FOR_ALL_ITEMS(gtri, tp) {
			int j;
			DEL_LINK(hl, tp);				/* Gone from the hit list */
			ADD_ITEM_TO_BOT(s->tris, tp);	/* Append to triangulation list */
			for (j = 0; j < 3 ; j++) {		/* Verticies weren't dropped from triangulation */
				tp->v[j]->f |= GVERT_TRI;
				tp->v[j]->f &= ~GVERT_INSIDE;
			}
		} END_FOR_ALL_ITEMS(tp);

		v->f |= GVERT_TRI;		/* This vertex has been added to triangulation */
		v->f &= ~GVERT_INSIDE;	/* and it's not inside */
	}

#ifdef DEBUG_TRIANG
	XLIST_FREE(&hittris);
#endif
}

/* - - - - - - - - - - - - - - - - */

/* Create the convex hull surface triangulation */
static void triangulate_ch(
gamut *s
) {
	/* Establish the base triangulation */
	{
		int i, j;
		gvert *tvs[4];	/* Initial verticies */
		gtri  *tr[4];	/* Initial triangles */
		gedge *ed[6];	/* Initial edges */

#ifdef USE_FAKE_SEED
		double fsz = FAKE_SEED_SIZE;		/* Initial tetra size */
		double ff[3];
		int onf;
		static double foffs[4][3] = {		/* tetrahedral offsets */
			{  0.0,   0.0,       1.0 },
			{  0.0,   0.80254,  -0.5 },
			{  0.75, -0.330127, -0.5 },
			{ -0.75, -0.330127, -0.5 }
		};

		/* Delete any current fake verticies */
		for (j = i = 0; i < s->nv; i++) {
			if (!(s->verts[i]->f & GVERT_FAKE))
				s->verts[j++] = s->verts[i];
			else
				free(s->verts[i]);
		}
		s->nv = j;

		/* Re-create fake points on each pass */
		onf = s->nofilter;
		s->nofilter = 1;			/* Turn off filtering */
		s->doingfake = 1;			/* Adding fake points */

		for (i = 0; i < 4; i++) {
			ff[0] = fsz * foffs[i][2] + s->cent[0];
			ff[1] = fsz * foffs[i][0] + s->cent[1];
			ff[2] = fsz * foffs[i][1] + s->cent[2];
			expand_gamut(s, ff);
		}

		s->nofilter = onf;
		s->doingfake = 0;

		/* Locate them for initial triangulation */
		for (j = i = 0; i < s->nv && j < 4; i++) {
			if (s->verts[i]->f & GVERT_FAKE) {
				tvs[j++] = s->verts[i];
			}
		}
#else	/* !USE_FAKE_SEED */
		/* Tetrahedral target points, outside log surface */
		/* (Assuming input values are < 1000) */
		double tgvals[4][3] = {
			{  0.0e4,   0.0e4,        1.0e4 },
			{  0.0e4,   0.8660254e4, -0.5e4 },
			{  0.75e4, -0.4330127e4, -0.5e4 },
			{ -0.75e4, -0.4330127e4, -0.5e4 }
		};
		double *tgts[4] = {
			tgvals[0], tgvals[1], tgvals[2], tgvals[3]
		};

		/* Find the closest verticies to the ideal tetrahedral */
		closest_verticies(s, tvs, tgts, 4);
#endif	/* !USE_FAKE_SEED */
		
#ifdef NEVER
		printf("Initial verticies:\n");
		for (i = 0; i < 4; i++) {
			printf(" %d: %f %f %f\n",tvs[i]->n, tvs[i]->p[0], tvs[i]->p[1], tvs[i]->p[2]);
		}
#endif
		/* Setup the initial triangulation */
		for (i = 0; i < 4; i++) {
			tr[i] = new_gtri();
		}

		for (i = 0; i < 6; i++) {
			ed[i] = new_gedge();
			ADD_ITEM_TO_BOT(s->edges, ed[i]);
		}

		/* Enter the edge verticies */
		ed[0]->v[0] = tvs[0];
		ed[0]->v[1] = tvs[1];
		ed[1]->v[0] = tvs[1];
		ed[1]->v[1] = tvs[2];
		ed[2]->v[0] = tvs[0];
		ed[2]->v[1] = tvs[2];
		ed[3]->v[0] = tvs[0];
		ed[3]->v[1] = tvs[3];
		ed[4]->v[0] = tvs[1];
		ed[4]->v[1] = tvs[3];
		ed[5]->v[0] = tvs[2];
		ed[5]->v[1] = tvs[3];

		/* Triangle facing in the +x, +y +z direction */
		tr[0]->v[0] = tvs[0];
		tr[0]->v[1] = tvs[1];
		tr[0]->v[2] = tvs[2];

		tr[0]->e[0] = ed[0];		/* Should make edge joining a function ? */
		tr[0]->ei[0] = 0;
		ed[0]->t[0] = tr[0];
		ed[0]->ti[0] = 0;
	
		tr[0]->e[1] = ed[1];
		tr[0]->ei[1] = 0;
		ed[1]->t[0] = tr[0];
		ed[1]->ti[0] = 1;
	
		tr[0]->e[2] = ed[2];
		tr[0]->ei[2] = 0;
		ed[2]->t[0] = tr[0];
		ed[2]->ti[0] = 2;
	
		comptriattr(s, tr[0]);				/* Compute triangle attributes */
		ADD_ITEM_TO_BOT(s->tris, tr[0]);	/* Append to list */
		
		/* Triangle facing in the -x, +y +z direction */
		tr[1]->v[0] = tvs[0];
		tr[1]->v[1] = tvs[3];
		tr[1]->v[2] = tvs[1];

		tr[1]->e[0] = ed[3];
		tr[1]->ei[0] = 0;
		ed[3]->t[0] = tr[1];
		ed[3]->ti[0] = 0;
	
		tr[1]->e[1] = ed[4];
		tr[1]->ei[1] = 0;
		ed[4]->t[0] = tr[1];
		ed[4]->ti[0] = 1;
	
		tr[1]->e[2] = ed[0];
		tr[1]->ei[2] = 1;
		ed[0]->t[1] = tr[1];
		ed[0]->ti[1] = 2;
	
		comptriattr(s, tr[1]);				/* Compute triangle attributes */
		ADD_ITEM_TO_BOT(s->tris, tr[1]);	/* Append to list */
		
		/* Triangle facing in the -y +z direction */
		tr[2]->v[0] = tvs[0];
		tr[2]->v[1] = tvs[2];
		tr[2]->v[2] = tvs[3];

		tr[2]->e[0] = ed[2];
		tr[2]->ei[0] = 1;
		ed[2]->t[1] = tr[2];
		ed[2]->ti[1] = 0;
	
		tr[2]->e[1] = ed[5];
		tr[2]->ei[1] = 0;
		ed[5]->t[0] = tr[2];
		ed[5]->ti[0] = 1;
	
		tr[2]->e[2] = ed[3];
		tr[2]->ei[2] = 1;
		ed[3]->t[1] = tr[2];
		ed[3]->ti[1] = 2;
	
		comptriattr(s, tr[2]);					/* Compute triangle attributes */
		ADD_ITEM_TO_BOT(s->tris, tr[2]);	/* Append to list */
		
		/* Triangle facing in the -z direction */
		tr[3]->v[0] = tvs[1];
		tr[3]->v[1] = tvs[3];
		tr[3]->v[2] = tvs[2];

		tr[3]->e[0] = ed[4];
		tr[3]->ei[0] = 1;
		ed[4]->t[1] = tr[3];
		ed[4]->ti[1] = 0;
	
		tr[3]->e[1] = ed[5];
		tr[3]->ei[1] = 1;
		ed[5]->t[1] = tr[3];
		ed[5]->ti[1] = 1;
	
		tr[3]->e[2] = ed[1];
		tr[3]->ei[2] = 1;
		ed[1]->t[1] = tr[3];
		ed[1]->ti[1] = 2;

		comptriattr(s, tr[3]);				/* Compute triangle attributes */
		ADD_ITEM_TO_BOT(s->tris, tr[3]);	/* Append to list */

		/* The four used verticies are now part of the triangulation */
		for (i = 0; i < 4; i++) {
			tvs[i]->f |= GVERT_TRI;
//printf("Base triangle %d: %d %d %d (Verticies 0x%x, 0x%x, 0x%x, 0x%x)\n",
//tr[i]->n, tr[i]->v[0]->n, tr[i]->v[1]->n, tr[i]->v[2]->n, tr[i]->v[0], tr[i]->v[1], tr[i]->v[2]);
		}
#ifdef ASSERTS
		check_triangulation(s, 0);
#endif
	}
	
	/* Sort the verticies from maximum radius, */
	/* to make our log convex hull logic work */
	sort_verticies(s);

	{
		int i;
		/* Complete the triangulation by adding all the remaining verticies */
		/* in order of decreasing radius, so that those below the log */
		/* convex hull get discarded. */
		for (i = 0; i < s->nv; i++) {
			if (!(s->verts[i]->f & GVERT_SET)
			 || (s->verts[i]->f & GVERT_TRI)
			 || (s->verts[i]->f & GVERT_INSIDE)) {
				continue;
			}

			insert_vertex(s, s->verts[i]);
#ifdef ASSERTS
			check_triangulation(s, 0);
#endif
		}
	}

	/* Number the used verticies */
	renumber_verticies(s);

#ifdef ASSERTS
	check_triangulation(s, 1);
#endif
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Compute a new convex hull mapping radius for the sample points, */
/* on the basis of the initial mapping. */
static void compute_smchrad(
gamut *s
) {
	int i, j;
	double zz[3] = { 0.0, 0.0, 1.0 };
	double rot[3][3];
	double ssr = 0.5 * s->sres;		/* Sample size radius in delta E */
//	double ssr = 10.0;
	int res = 4;			/* resolution of sample grid */

//printf("~1 computing smoothed chrads\n");

	/* Compute the new log surface value */
	for (i = 0; i < s->nv; i++) {
		double pr;			/* Current surface radius at this vertex */ 
		double rad, rw;		/* Smoothed radius, weight */
		double out[3];		/* Point on surface for this vertex */
		int x, y;

		if (!(s->verts[i]->f & GVERT_SET)) {
			continue;
		}

//printf("~1 vertex %d, %f %f %f\n",i, s->verts[i]->p[0], s->verts[i]->p[1], s->verts[i]->p[2]);

		/* Find the average surface level near this point */
		pr = s->radial(s, out, s->verts[i]->p);

//printf("~1 surface radius %f, location %f %f %f\n",pr, out[0],out[1],out[2]);

		/* Compute a rotation that lines Z up with the radial direction */
		out[0] -= s->cent[0];		/* Radial vector through point to surface */
		out[1] -= s->cent[1];
		out[2] -= s->cent[2];
		zz[2] = pr;					/* Z vector of same length */
		icmRotMat(rot, zz, out);	/* Compute vector from Z to radial */
		out[0] += s->cent[0];
		out[1] += s->cent[1];
		out[2] += s->cent[2];
		rad = rw = 0.0;
	
		/* Sample a rectangular array orthogonal to radial vector, */
		/* and weight samples appropriately */
		for (x = 0; x < res; x++) {
			for (y = 0; y < res; y++) {
				double tt, off[3], rv, in[3];

				off[0] = 2.0 * (x/(res-1.0) - 0.5);		/* -1.0 to 1.0 */
				off[1] = 2.0 * (y/(res-1.0) - 0.5);		/* -1.0 to 1.0 */
				off[2] = 0.0;

				rv = off[0] * off[0] + off[1] * off[1];
				if (rv > 1.0)
					continue;			/* Outside circle */
#ifdef NEVER
				rv = 1.0 - sqrt(rv);	/* Radius from center */
				rv = rv * rv * (3.0 - 2.0 * rv);	/* Spline weight it */
#else
				rv = 1.0;				/* Moving average weighting */
#endif

				off[0] *= ssr;			/* Scale offset to sampling radius */
				off[1] *= ssr;

				/* Rotate offset to be orthogonal to radial vector */
				icmMulBy3x3(off, rot, off);
			
				/* Add offset to surface point over current vertex */
				in[0] = out[0] + off[0];
				in[1] = out[1] + off[1];
				in[2] = out[2] + off[2];

//printf("~1 grid %d %d, weight %f, offset %f %f %f\n",x,y,rv,off[0],off[1],off[2]);

				/* Sum weighted radius at sample point */
				tt = s->radial(s, NULL, in);
				tt = log_scale(tt);
				rad += rv * tt;
				rw += rv;

			}
		}
		/* Compute sample filtered radius at the sample point */
		rad /= rw;
//printf("~1 sampled radius = %f\n\n",rad);

		/* Now compute new hull mapping radius for this point, */
		/* based on dividing out the sampled radius */

//		s->verts[i]->lr0 = 40.0 + s->verts[i]->lr0 - rad;
		s->verts[i]->lr0 = 40.0 + log_scale(s->verts[i]->r[0]) - rad;
		/* Prevent silliness */
		if (s->verts[i]->lr0 < (2.0 * FAKE_SEED_SIZE))
			s->verts[i]->lr0 = (2.0 * FAKE_SEED_SIZE);

//printf("~1 new lr0 = %f\n\n",s->verts[i]->lr0);

		/* recompute ch[] for new lr0 */
		for (j = 0; j < 3; j++)
			s->verts[i]->ch[j] = s->verts[i]->sp[j] * s->verts[i]->lr0;
	}
}

/* ===================================================== */
/* Overall triangulation */
static void triangulate(
gamut *s
) {

	/* Create the convex hull */
	triangulate_ch(s);

#if defined(DO_TWOPASS) && !defined(TEST_CONVEX_HULL)
#ifdef DEBUG_TRIANG
	printf("############ Starting second pass ###################\n");
#endif
	compute_smchrad(s);
	del_triang(s);
	s->pass++;
	triangulate_ch(s);

	/* Three passes is slightly better, but slower... */
//	compute_smchrad(s);
//	del_triang(s);
//	s->pass++;
//	triangulate_ch(s);
#endif /* DO_TWOPASS && !TEST_CONVEX_HULL */
}

/* ===================================================== */
/* Code that makes use of the triangulation */
/* ===================================================== */

/* return the current surface resolution */
static double getsres(
gamut *s
) {
	return s->sres;
}

/* return the isJab flag value */
static int getisjab(
gamut *s
) {
	return s->isJab;
}

/* return nz if the two gamut are compatible */
static int compatible(
gamut *s, struct _gamut *t) {
	int j;

	/* The same colorspace ? */
	if ((s->isJab && !t->isJab)
	 || (!s->isJab && t->isJab)) {
		return 0;
	}

	/* The same gamut center ? */
	for (j = 0; j < 3; j++) {
		if (fabs(s->cent[j] - t->cent[j]) > 1e-9) {
			return 0;
		}
	}
	return 1;
}


/* Return the number of raw verticies used to construct surface */
static int nrawverts(
gamut *s
) {
	int i, nrv = 0;

	/* Sort them so that triangulate doesn't mess indexing up */
	sort_verticies(s);

	/* Count them */
	for (i = 0; i < s->nv; i++) {
		if (s->verts[i]->f & GVERT_SET)
			nrv++;
	}

	return nrv;
}

/* Return the raw (triangle and non-triangle surface) verticies */
/* location given its index. */
/* return the next (sparse) index, or -1 if beyond last */
static int getrawvert(
gamut *s,
double pos[3],		/* Return absolute position */
int ix				/* Input index */
) {
	if (ix < 0)
		return -1;

	/* Find then next used in the triangulation */
	for (; ix < s->nv; ix++) {
		if (!(s->verts[ix]->f & GVERT_SET))
			continue;
		break;
	}

	if (ix >= s->nv)
		return -1;

	pos[0] = s->verts[ix]->p[0];
	pos[1] = s->verts[ix]->p[1];
	pos[2] = s->verts[ix]->p[2];

	return ix+1;
}

/* Return the number of raw direction 0 verticies used */
/* to construct surface. (Direction 0 is radial direction maxima) */
static int nraw0verts(
gamut *s
) {
	int i, nrv = 0;

	/* Sort them so that triangulate doesn't mess indexing up */
	sort_verticies(s);

	/* Count them */
	for (i = 0; i < s->nv; i++) {
		if ((s->verts[i]->f & GVERT_SET)
		 && (s->verts[i]->k0 > 0))
			nrv++;
	}

	return nrv;
}

/* Return the raw (triangle and non-triangle surface)  direction 0 */
/* verticies location given its index. (Direction 0 is radial direction maxima) */
/* return the next (sparse) index, or -1 if beyond last */
static int getraw0vert(
gamut *s,
double pos[3],		/* Return absolute position */
int ix				/* Input index */
) {
	if (ix < 0)
		return -1;

	/* Find then next used in the triangulation and direction 0 */
	for (; ix < s->nv; ix++) {
		if (!(s->verts[ix]->f & GVERT_SET)
		 || !(s->verts[ix]->k0 > 0))
			continue;
		break;
	}

	if (ix >= s->nv)
		return -1;

	pos[0] = s->verts[ix]->p[0];
	pos[1] = s->verts[ix]->p[1];
	pos[2] = s->verts[ix]->p[2];

	return ix+1;
}

/* Return the number of verticies in the triangulated surface */
static int nverts(
gamut *s
) {
	if IS_LIST_EMPTY(s->tris)
		triangulate(s);

	return s->ntv;
}

/* Return the verticies location and radius given its index. */
/* return the next (sparse) index, or -1 if beyond last */
static int getvert(
gamut *s,
double *rad,		/* Return radial radius */
double pos[3],		/* Return absolute position */
int ix				/* Input index */
) {
	if (ix >= s->nv)
		return -1;

	/* Find then next used in the triangulation */
	for (; ix < s->nv; ix++) {
		if (!(s->verts[ix]->f & GVERT_TRI))
			continue;
		break;
	}

	if (rad != NULL)
		*rad   = s->verts[ix]->r[0];
	if (pos != NULL) {
		pos[0] = s->verts[ix]->p[0];
		pos[1] = s->verts[ix]->p[1];
		pos[2] = s->verts[ix]->p[2];
	}

	return ix+1;
}


/* Reset indexing through triangles for getnexttri() */
static void startnexttri(gamut *s) {
	if IS_LIST_EMPTY(s->tris)
		triangulate(s);

	s->nexttri = NULL;
}

/* Return the next surface triange, nz on no more */
static int getnexttri(
gamut *s,
int v[3]		/* Return indexes for same order as getvert() */
) {
	if IS_LIST_EMPTY(s->tris)
		triangulate(s);

	if (s->nexttri == NULL) {
		s->nexttri = s->tris;
		if (s->nexttri == NULL)
			return 1;
	} else {	
		s->nexttri = NEXT_FWD(s->nexttri);
		if (s->nexttri == s->tris)
			return 1;
	}

	v[0] = s->nexttri->v[0]->tn;
	v[1] = s->nexttri->v[1]->tn;
	v[2] = s->nexttri->v[2]->tn;
	return 0;
}

/* ===================================================== */

/* Return the total volume of the gamut */
/* [Could we subtract the volumes of two gamuts by substracting */
/*  the dotproductarea of the two sufaces ??] */
static double volume(
gamut *s
) {
	int i, j;
	gtri *tp;		/* Triangle pointer */
	double vol;		/* Gamut volume */

	if IS_LIST_EMPTY(s->tris)
		triangulate(s);

	/* Compute the area of each triangle in the list, */
	/* and accumulate the gamut volume. */
	tp = s->tris; 
	vol = 0.0;
	FOR_ALL_ITEMS(gtri, tp) {
		double sp, ss[3];		/* Triangle side lengths */
		double area;			/* Area of this triangle */
		double dp;				/* Dot product of point in triangle and normal */
		
		for (i = 0; i < 3; i++) {	/* For each edge */
			for (ss[i] = 0.0, j = 0; j < 3; j++) {
				double dd = tp->e[i]->v[1]->p[j] - tp->e[i]->v[0]->p[j];
				ss[i] += dd * dd;
			}
			ss[i] = sqrt(ss[i]);
		}

		/* semi-perimeter */
		sp = 0.5 * (ss[0] + ss[1] + ss[2]);

		/* Area of triangle */
		area = sqrt(sp * (sp - ss[0]) * (sp - ss[1]) * (sp - ss[2]));
		
		/* Dot product between first vertex in triangle and the unit normal vector */
		dp = tp->v[0]->p[0] * tp->pe[0]
		   + tp->v[0]->p[1] * tp->pe[1]
		   + tp->v[0]->p[2] * tp->pe[2];

		/* Accumulate gamut volume */
		vol += dp * area;

	} END_FOR_ALL_ITEMS(tp);

	vol = fabs(vol)/3.0;

	return vol;
}

/* ===================================================== */
/* ===================================================== */
/* Given a point, */
/* return the distance to the gamut surface. */

static void init_lu(gamut *s);
static gtri *radial_point_triang(gamut *s, gbsp *np, double in[3]);
static double radial_point(gamut *s, gbsp *np, double in[3]);

/* Implementation for following two functions: */
/* Given a point, return the point in that direction */
/* that lies on the gamut surface. */
/* Return the radial length of the input and radial length of result */
static void
_radial(
gamut *s,
double *ir,		/* return input radius (may be NULL) */
double *or,		/* return output radius (may be NULL) */
double *out,	/* result point (absolute) (may be NULL) */
double *in		/* input point (absolute)*/
) {
	int j;
	double ss, rv;
	double nin[3];	/* Normalised input vector */

	if IS_LIST_EMPTY(s->tris)
		triangulate(s);
		
	/* We have to find out which triangle the point is in */
	if (s->lu_inited == 0) {
		init_lu(s);				/* Init BSP search tree */
	}
//if (trace) printf("~1 radial called with %f %f %f\n", in[0], in[1], in[2]);

	for (j = 0; j < 3; j++)
		nin[j] = in[j] - s->cent[j]; /* relative to gamut center */

	for (ss = 0.0, j = 0; j < 3; j++)
		ss += nin[j] * nin[j];
	ss = sqrt(ss);
	if (ss > 1e-9) {				/* Normalise to 1.0 */
		for (j = 0; j < 3; j++)
			nin[j] /= ss;
	} else {
		nin[0] = 1.0;
		nin[1] = nin[2] = 0.0;
	}

//if (trace) printf("~1 Normalised in = %f %f %f\n", nin[0], nin[1], nin[2]);
	rv = radial_point(s, s->lutree, nin);

	if (rv < 0.0) {
		fprintf(stderr,"gamut: radial internal error - failed to find triangle\n");
		exit (-1);
	}

	if (out != NULL) {
		for (j = 0; j < 3; j++)
			out[j] = nin[j] * rv + s->cent[j];		/* Scale out to surface length, absolute */
//if (trace) printf("~1 result = %f %f %f\n",out[0], out[1], out[2]);
	}

	if (ir != NULL) {
//if (trace) printf("~1 input radius res = %f\n",ss);
		*ir = ss;
	}

	if (or != NULL) {
//if (trace) printf("~1 output radius res = %f\n",rv);
		*or = rv;
	}
}

/* Given a point, return the point in that direction */
/* that lies on the gamut surface */
/* Return the normalised radial radius to the surface point */
static double
nradial(
gamut *s,
double *out,	/* result point (absolute) (May be NULL) */
double *in		/* input point (absolute)*/
) {
	double ss, rv;
	
	_radial(s, &ss, &rv, out, in);
	return ss/rv;
}

/* Given a point, return the point in that direction */
/* that lies on the gamut surface */
/* Return the radial radius to the surface point */
static double
radial(
gamut *s,
double *out,	/* result point (absolute) (May be NULL) */
double *in		/* input point (absolute)*/
) {
	double ss, rv;
	
	_radial(s, &ss, &rv, out, in);
	return rv;
}

void lu_split(gamut *s, gbsp **np, int rdepth, gtri **list, int llen);

/* Setup the radial lookup function acceleration structure */
static void
init_lu(
gamut *s
) {
	static double v0[3] = {0.0, 0.0, 0.0};
	static 
	gedge *ep;		/* Edge pointer */
	gtri *tp;		/* Triangle pointer */
	gtri **tlist;
	int ntris;

//printf("~1 init_lu called\n");

	/* Create mean angle dividing plane equations */
	ep = s->edges; 
	FOR_ALL_ITEMS(gedge, ep) {
		plane_equation(ep->re, v0, ep->v[0]->sp, ep->v[1]->sp);
	} END_FOR_ALL_ITEMS(ep);

	/* Create the initial triangle list */
	/* First count them */
	ntris = 0;
	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		ntris++;
	} END_FOR_ALL_ITEMS(tp);

	/* Allocate a list */
	if ((tlist = (gtri **) malloc(ntris * sizeof(gtri *))) == NULL) {
		fprintf(stderr,"gamut: malloc failed - top level triangle list (%d entries)\n",ntris);
		exit(-1);
	}

	/* Then add them to the list */
	ntris = 0;
	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		tlist[ntris] = tp;
		ntris++;
	} END_FOR_ALL_ITEMS(tp);

	/* Recursively split them */
	lu_split(s, &s->lutree, 0, tlist, ntris);

	free(tlist);

//printf("~1 init_lu done\n");
	s->lu_inited = 1;
}

/* Recursive routine to choose a partition plane, */
/* and then split the triangle list between the */
/* +ve and -ve sides. */
void
lu_split(
gamut *s,
gbsp **np,		/* Address of node pointer to be set */
int rdepth,		/* Current recursion depth */
gtri **list,	/* Current triangle list */
int llen		/* Number of triangles in the list */
) {
	int ii, jj;			/* Progress through edges */
	int pcount;			/* Current best scored try */
	int ncount;
	int bcount;
	int mcount;
	double peqs[4] = { 0.0, 0.0, 0.0, 0.0 };
	gtri **plist, **nlist;	/* New sub-lists */
	int pix, nix;			/* pos/ned sublist indexes */
	gbspn *bspn;			/* BSP decision node */

//printf("~1\nlu_split called at depth %d with %d triangles\n",rdepth, llen);
#ifdef DEBUG
	if (llen <= 3) {
		int i;
		for (i = 0; i < llen; i++) {
			printf("Triang index %d = %d\n",i, list[i]->n);
			printf("Triang verts %d %d %d\n", 
			list[i]->v[0]->tn, list[i]->v[1]->tn, list[i]->v[2]->tn);
			printf("Vert 0 at %.18f %.18f %.18f\n",list[i]->v[0]->sp[0], list[i]->v[0]->sp[1], list[i]->v[0]->sp[2]);
			printf("Vert 1 at %.18f %.18f %.18f\n",list[i]->v[1]->sp[0], list[i]->v[1]->sp[1], list[i]->v[1]->sp[2]);
			printf("Vert 2 at %.18f %.18f %.18f\n",list[i]->v[2]->sp[0], list[i]->v[2]->sp[1], list[i]->v[2]->sp[2]);
		}
	}
#endif /* DEBUG */

	if ((rdepth+1) >= BSPDEPTH) {	/* Oops */
		printf("gamut internal error: ran out of recursion depth in BSP\n");
		exit (-1);
	}

	pcount = ncount = bcount = -1;
	mcount = 0;
	/* test every edge in turn */
	for (ii = jj = 0;ii < llen;) {
		double eqs[4];
		int i;
		gedge *ep;			/* Edge pointer */
		int pc, nc, bc;		/* Score a try, postive count, negative count, both count */
		int mc;				/* Minumum count */

		ep = list[ii]->e[jj];
		eqs[0] = ep->re[0];	/* Use this edge */
		eqs[1] = ep->re[1];
		eqs[2] = ep->re[2];
		eqs[3] = ep->re[3];
		if (++jj > 2) {
			jj = 0;
			ii++;
		}

		/* Do the trial split */
		pc = nc = bc = 0;
		for (i = 0; i < llen; i++) {
			int j;
			int po, ne;

			/* Compute distance from plane of all verticies in triangle */
			po = ne = 0;
			for (j = 0; j < 3; j++) {	/* For triangle verticies */
				double ds;
				/* Compute distance to dividing plane of this vertex */
				ds = eqs[0] * list[i]->v[j]->sp[0]
				   + eqs[1] * list[i]->v[j]->sp[1]
	    		   + eqs[2] * list[i]->v[j]->sp[2]
				   + eqs[3];
				/* Figure if the verticies are clearly to one side of the plane */
				if (ds > 1e-10) {
					po++;
				} else if (ds < -1e-10) {
					ne++;
				}
			}
			/* Score this split */
			if (po) {
				pc++;
				if (ne) {
					nc++;
					bc++;
					list[i]->sort = 3;	/* Both */
				} else {
					list[i]->sort = 1;	/* +ve */
				}
			} else if (ne) {
				nc++;
				list[i]->sort = 2;	/* -ve */
			} else {				/* Hmm. Neither */
				bc++;
				list[i]->sort = 3;	/* Assume both */
			}
		}
		mc = pc < nc ? pc : nc;		/* Size of smallest group */
		mc -= bc;
//printf("~1 lu_split trial %d, mc %d, pc %d, nc %d, bc %d\n",ii * 3 + jj, mc, pc, nc, bc);
		if (mc > mcount) {			/* New largest small group */
			mcount = mc;
			pcount = pc;
			ncount = nc;
			bcount = bc;
			peqs[0] = eqs[0];
			peqs[1] = eqs[1];
			peqs[2] = eqs[2];
			peqs[3] = eqs[3];
//printf("~1 new best - plane mc = %d, %f %f %f %f\n",mc, peqs[0], peqs[1], peqs[2], peqs[3]);
			for (i = 0; i < llen; i++) {
				list[i]->bsort = list[i]->sort;
			}
		}
	}

#ifdef DEBUG_SPLIT_VRML
	write_split_diag_vrml(s, list, llen);
	getchar();
#endif /* DEBUG_SPLIT_VRML */

	if (ii >= llen && bcount < 0) {	/* We failed to find a split plane. */
		/* This is usually a result of the list being 2 or more triangles */
		/* that do not share any edges (disconected from each other), and */
		/* lying so that any split plane formed from an edge of one, */
		/* intersects one of the others. */

		/* Leave our list of triangles as the leaf node, */
		/* and let the search algorithms deal with this. */

		*np = (gbsp *)new_gbspl(llen, list);
//printf("~1 lu_split returning with a non split list of %d triangles\n",llen);
		return;
	}

	/* Divide the triangles into two lists */

	bspn = new_gbspn();				/* Next node */
	*np = (gbsp *)bspn;				/* Put it in place */
	bspn->pe[0] = peqs[0];			/* Plane equation */
	bspn->pe[1] = peqs[1];
	bspn->pe[2] = peqs[2];
	bspn->pe[3] = peqs[3];

	/* Allocate the sub lists */
	if ((plist = (gtri **) malloc(pcount * sizeof(gtri *))) == NULL) {
		fprintf(stderr,"gamut: malloc failed - pos sub-list\n");
		exit(-1);
	}
	if ((nlist = (gtri **) malloc(ncount * sizeof(gtri *))) == NULL) {
		fprintf(stderr,"gamut: malloc failed - neg sub-list\n");
		exit(-1);
	}

	/* Fill them in */
	for (pix = nix = ii = 0; ii < llen; ii++) {
		if (list[ii]->bsort & 1) {	/* Positive */
			plist[pix] = list[ii];
			pix++;
		}
		if (list[ii]->bsort & 2) {	/* Negative */
			nlist[nix] = list[ii];
			nix++;
		}
	}

	/* Recurse if there are more triangles to split */
	if (pix == 1) {
		bspn->po = (gbsp *)plist[0];			/* leaf node */ 
//printf("~1 pos leaf with triangle %d\n",plist[0]->n);
	} else if (pix > 1) {
//printf("~1 About to recurse on positive with list of %d\n",pix);
		lu_split(s, &bspn->po, rdepth+1, plist, pix);
	}

	if (nix == 1) {
//printf("~1 neg leaf with triangle %d\n",nlist[0]->n);
		bspn->ne = (gbsp *)nlist[0];			/* leaf node */ 
	} else if (nix > 1) {
//printf("~1 About to recurse on negative with list of %d\n",nix);
		lu_split(s, &bspn->ne, rdepth+1, nlist, nix);
	}

	free(plist);
	free(nlist);
//printf("~1 lu_split returning\n");
}

/* Given a point and a node in the BSP tree, recurse down */
/* the correct side of the tree, or return the triangle or */
/* triangles that could be on the gamut surface. Return NULL */
/* if it wasn't in triangle. */
static gtri *radial_point_triang(
gamut *s,
gbsp *np,		/* BSP node pointer we're at */
double *nin		/* Normalised center relative point */
) {
//if (trace) printf("~1 rad_pnt_tri: BSP 0x%x tag = %d, point %f %f %f\n", np,np->tag,nin[0],nin[1],nin[2]);
	if (np->tag == 1) {		/* It's a BSP node */
		gbspn *n = (gbspn *)np;
		double ds;

		ds = n->pe[0] * nin[0]
		   + n->pe[1] * nin[1]
	       + n->pe[2] * nin[2]
		   + n->pe[3];

//if (trace) printf("~1 checking against BSP plane, ds = %e\n",ds);
		if (ds >= 0)
			return radial_point_triang(s, n->po, nin);
		else
			return radial_point_triang(s, n->ne, nin);

	} else if (np->tag == 2) {	/* It's a triangle */
//if (trace) printf("~1 returning triangle %d\n",((gtri *)np)->n);
		return (gtri *)np;

	} else {			/* It's a triangle list */	
		gbspl *n = (gbspl *)np;
		int i, j;
//if (trace) printf("~1 got triangle list - radial failed to split triangles!\n");

		/* Go through the list and stop at the first triangle */
		/* that the node lies in. */
		for (i = 0; i < n->nt; i++) {
			/* Check if the point is within this triangle */
			for (j = 0; j < 3; j++) {
				double ds;
				ds = n->t[i]->ee[j][0] * nin[0]
				   + n->t[i]->ee[j][1] * nin[1]
			       + n->t[i]->ee[j][2] * nin[2]
				   + n->t[i]->ee[j][3];
				if (ds > 1e-10)
					break;			/* Not within triangle */
			}
			if (j >= 3) {
//if (trace) printf("~1 located triangle from list that we're in %d\n",n->t[i]->n);
				return n->t[i];
			}
		}
	}

//if (trace) printf("~1 failed to find a triangle\n");
	return NULL;
}

/* Return the location on the surface of the triangle */
/* that is intersected by the vector in the direction */
/* of the given relative point.  Return the distance to */
/* the gamut surface. Return -1 if it wasn't in triangle. */
static double radial_point(
gamut *s,
gbsp *np,		/* BSP node pointer we're at */
double *nin		/* Normalised center relative point */
) {
	gtri *t;
	double rv;

//if (trace) printf("~1 radial_point: BSP 0x%x tag = %d, point %f %f %f\n", np,np->tag,nin[0],nin[1],nin[2]);

	t = radial_point_triang(s, np, nin);

	/* Due to numerical issues, the BSP accelerated result may (very occassionally) */
	/* be in error. Double check that we've landed in the correct triangle. */
	if (t != NULL) {
		double ds;
		int j;

		/* Check if the closest point is within this triangle */
		for (j = 0; j < 3; j++) {
			double ds;
			ds = t->ee[j][0] * nin[0]
			   + t->ee[j][1] * nin[1]
		       + t->ee[j][2] * nin[2]
			   + t->ee[j][3];
//			if (ds > 1e-10) {
			if (1) {
//if (trace) printf("radial: lookup point wasn't within its BSP triangle (%f) !!\n",ds);
				t = NULL;
				break;
			}
		}

	}

	/* If we failed to find a triangle, or the result was incorrect, do a */
	/* brute force search to be sure of the result. */
	if (t == NULL) {
		gtri *tp, *btp = NULL;		/* Triangle pointer */
		double bds = 1e38;			/* Track smallest "out of side" */
		int j;

		tp = s->tris;
		FOR_ALL_ITEMS(gtri, tp) {
			double tds = -1e38;		/* Track bigest "out of side" */

			/* Track the triangle that the point most comfortably lands in */
			for (j = 0; j < 3; j++) {
				double ds;
				ds = tp->ee[j][0] * nin[0]
				   + tp->ee[j][1] * nin[1]
			       + tp->ee[j][2] * nin[2]
				   + tp->ee[j][3];
				if (ds > tds) {
					tds = ds;
				}
			}
			if (tds < bds) {
				bds = tds;
				btp = tp;
//if (trace) printf("radial: found new best triangle %d bds %f\n",tp->n,bds);
			}
		} END_FOR_ALL_ITEMS(tp);

		if (bds < 1e-10)	/* Found an acceptable result */
			t = btp;
	}

	if (t == NULL)
		return -1.0;

	/* Compute the intersection of the input vector with the triangle plane */
	/* (Since nin[] is already relative, we don't need to subtract cent[] from it) */
	rv = -(t->pe[0] * s->cent[0] + t->pe[1] * s->cent[1] + t->pe[2] * s->cent[2] + t->pe[3])/
			      (t->pe[0] * nin[0] + t->pe[1] * nin[1] + t->pe[2] * nin[2]);

#ifdef ASSERTS
	/* check the result */
	{
		double tt[3];
		double ds;
		int j;
		for (j = 0; j < 3; j++)			/* Compute result, absolute */
			tt[j] = nin[j] * rv + s->cent[j];
		
		ds = t->pe[0] * tt[0]
		   + t->pe[1] * tt[1]
		   + t->pe[2] * tt[2]
		   + t->pe[3];
		
		if (fabs(ds) > 1e-6) {
			fprintf(stderr,"radial: distance to plane not zero! %e\n",ds);
			exit(-1);
		}
		
		/* Check if the closest point is within this triangle */
		for (j = 0; j < 3; j++) {
			double ds;
			ds = t->ee[j][0] * (tt[0] - s->cent[0])
			   + t->ee[j][1] * (tt[1] - s->cent[1])
		       + t->ee[j][2] * (tt[2] - s->cent[2])
			   + t->ee[j][3];
			if (ds > 1e-8) {
				fprintf(stderr,"radial: lookup point wasn't within its triangle (%f) !!\n",ds);
				exit(-1);
			}
		}
	}
#endif /* ASSERTS */

//if (trace) printf("~1 radial_point: rv = %f\n",rv);
	return rv;
}

/* Recursively free a gbsp node and all its children */
static void del_gbsp(gbsp *n) {
	int tag = n->tag;

	if (tag == 1) {				/* Another decision node */
		gbspn *dn = (gbspn *)n;
		del_gbsp(dn->po);		/* Delete children */
		del_gbsp(dn->ne);
		del_gbspn(dn);			/* And itself */

	} else if (tag == 3) {		/* If a triangle list */
		gbspl *dl = (gbspl *)n;
		del_gbspl(dl);			/* Delete itself */
	}

	/* Don't delete triangles since they have their own linked list. */
}

/* =================================== */
/* Given a point, */
/* return the nearest point on the gamut surface. */

#define GNN_INF 1e307
static void init_ne(gamut *s);
static double ne_point_on_tri(gamut *s, gtri *t, double *out, double *in);

/* Given an absolute point, return the point on the gamut */
/* surface that is closest to it. */
static void
nearest(
gamut *s,
double *rout,	/* result point (absolute) */
double *q		/* Target point (absolute) */
) {
	gnn *p;				/* Pointer to nearest neighbor structure */
	int e, i;
	double r[3] = {0.0, 0.0, 0.0 };	/* Possible solution point */
	double out[3];		/* Current best output value */
	int wex[3 * 2];		/* Current window edge indexes */
	double wed[3 * 2];	/* Current window edge distances */
						/* Indexes are axis * 2 +0 for lower edge, */
						/* +1 for upper edge of search box. */
						/* We are comparing lower edge of search box */
						/* with upper edge of bounding box etc. */ 

//printf("~1 nearest called\n");
	if IS_LIST_EMPTY(s->tris)
		triangulate(s);

	/* We have to find out which triangle the point will be nearest */
	if (s->ne_inited == 0) {
		init_ne(s);				/* Init nn structure */
	}
	p = s->nns;

	if ((p->tbase + 3) < p->tbase) {	/* Overflow of touch count */
		for (i = 0; i < p->n; i++)
			p->sax[0][i]->touch = 0;		/* reset it in all the objects */
		p->tbase = 0;
	}
	p->ttarget = p->tbase + 3;		/* Target touch value */

//printf("\n");
//printf("Query point is %f %f %f\n",q[0], q[1], q[2]);

	/* Find starting indexes within axis arrays */
	for (e = 0; e < (2 * 3); e++) {	/* For all axes min & max */
		int f = e/2;			/* Axis */
		int mm = (e ^ 1) & 1;	/* Min/Max index used for edges */
		int i0, i1, i2;
		double v0, v1, v2;
		double qf, ww;

		/* Binary search this edge */
		qf = q[f]; 		/* strength reduced q[f] */

//printf("\n");
//printf("isearching axis %d %s for %f\n",f, e & 1 ? "max" : "min", qf);
		i0 = 0;
		i2 = p->n - 1;
		v0 = p->sax[e][i0]->mix[mm][f];
		v2 = p->sax[e][i2]->mix[mm][f];
//printf("start points %d - %d, bound %f - %f\n",i0, i2, v0, v2);

		if (qf <= v0) {
			i2 = i0;
			v2 = v0;
		} else if (qf >= v2) {
			i0 = i2;
			v0 = v2;
		} else {
			do {
				i1 = (i2 + i0)/2;		/* Trial point */
				v1 = p->sax[e][i1]->mix[mm][f];	/* Value at trial */
				if (v1 < qf) {
					i0 = i1;			/* Take top half */
					v0 = v1;
				} else {
					i2 = i1;			/* Take bottom half */
					v2 = v1;
				}
//printf("current point %d - %d, bound %f - %f\n",i0, i2, v0, v2);
			} while ((i2 - i0) > 1);
		}

		if (e & 1) {			/* Max side of window */
			int tc;				/* total object count */

			ww = v2 - qf;
			wed[e] = fabs(ww) * ww;
			wex[e] = i2;

			/* Check that min and max together will cover at least p->n objects */
			tc = p->n - i2 + wex[e ^ 1] + 1;
//printf("got %d, expected %d\n",tc, p->n);

			/* (I don't really understand why this works!) */
			if (tc < p->n) {		/* We haven't accounted for all the objects */
				int el = e ^ 1;		/* Low side sax */
				int ti0, ti2;
				double tv0, tv2;

				ti0 = wex[el];
				ti2 = i2;
//printf("We have straddling objects, initial indexes are %d - %d\n",ti0, ti2);

				/* While straddling objects remain undiscovered: */
				while (tc < p->n) {
					tv0 =  GNN_INF;		/* Guard values */
					tv2 = -GNN_INF;

					/* Increment low side until we find a straddler */
					while (ti0 < (p->n-1)) {
						ww = p->sax[el][++ti0]->mix[0][f];	/* Position of the other end */
						if (ww < qf) {
//printf("found low object %d at index %d that straddles\n",p->sax[el][ti0]-p->base,ti0);
							tv0 = qf - p->sax[el][ti0]->mix[1][f];
							break;
						}
					}

					/* Decrement high side until we find a straddler */
					while (ti2 > 0) {
						ww = p->sax[e][--ti2]->mix[1][f];	/* Position of the other end */
						if (ww > qf) {
//printf("found high object %d at index %d that straddles\n",p->sax[e][ti2]-p->base,ti2);
							tv2 = p->sax[e][ti2]->mix[0][f] - qf;
							break;
						}
					}
					/* Choose the closest */
					if (tv0 > tv2) {
						wed[el] = fabs(tv0) * tv0;
						wex[el] = ti0;
						tc++;
					} else {
						wed[e] = fabs(tv2) * tv2;
						wex[e] = ti2;
						tc++;
					}
				}
//printf("After correction we have %d - %d\n",wex[e^1], wex[e]);
			}
		} else {				/* Min side of window */
			ww = q[f] - v0;
			wed[e] = fabs(ww) * ww;
			wex[e] = i0;
		}
	}

	/* Expand a 3 dimenstional cube centered on the target point, */
	/* jumping to the next nearest point on any axis, discovering */
	/* any bounding boxes that are within the expanding window */
	/* by checking their touch count. */

	/* The first point found establishes the initial best distance. */
	/* When the window expands beyond the point where it can improve */
	/* the best distance, stop */

	{
		double bw = 0.0;		/* Current window distance */
		double bdist = 1e308;	/* Best possible distance to an object outside the window */
		gtri *bobj = NULL;
		int ptested = 0;		/* Stats */
		int pcalced = 0;		/* Stats */

		/* Until we're done */
		for (;;ptested++) {
			int ee;				/* Axis & expanding box edge */
			int ff;				/* Axis */
			int ii;				/* Index of chosen point */
			gtri *ob;			/* Current object */
			unsigned int ctv;	/* Current touch value */
//printf("\n");
//printf("wwidth = %f, bdist = %f, window = %d-%d, %d-%d, %d-%d\n",
//bw, bobj == NULL ? 0.0 : bdist, wex[0], wex[1], wex[2], wex[3], wex[4], wex[5]);
//printf("window edge distances are = %f-%f, %f-%f, %f-%f\n",
//wed[0], wed[1], wed[2], wed[3], wed[4], wed[5]);

			/* find next (smallest) window increment axis and direction */
			ee = 0;
			ii = wex[ee];
			bw = wed[ee];
			for (e = 1; e < (2 * 3); e++) {
				if (wed[e] < bw) {
					ee = e;
					ii = wex[e];
					bw = wed[e];
				}
			}
//printf("Next best is axisdir %d, object %d, axis index %d, best possible dist %f\n",
//ee, p->sax[ee][ii] - p->base, ii, bw);

			if (bw == GNN_INF || bw > bdist) {
				break;		/* Can't got any further, or further points will be worse */
			}

#ifdef ASSERTS
if (ii < 0 || ii >= p->n) {
printf("Assert: went out of bounds of sorted axis array\n");
exit(0);
}
#endif
			/* Chosen point on ee axis/direction, index ii */
			ff = ee / 2;			/* Axis only */

			ob = p->sax[ee][ii];

			/* Touch value of current object */
			ctv = ob->touch;

			if (ctv < p->ttarget) {		/* Not been dealt with before */

				/* Touch this new window boundary point */
				ob->touch = ctv = ((ctv < p->tbase) ? p->tbase : ctv) + 1;

//printf("New touch count on %d is %d, target %d\n",
//ob - p->base, p->sax[ee][ii]->touch, p->ttarget);

				/* Check the point out */
				if (ctv == (p->tbase + 3)) {	/* Is within window on all axes */
					double tdist;

					pcalced++;		/* Stats */

					/* Compute distance from query point to this object */
					tdist = ne_point_on_tri(s, ob, r, q);

//printf("Got new best point %d, dist %f\n",i,tdist);
					if (tdist < bdist) {	/* New best point */
						bobj = ob;
						bdist = tdist;
						out[0] = r[0];
						out[1] = r[1];
						out[2] = r[2];
					}
				}
			}

			/* Increment next window edge candidate, and figure new edge distance */
			if (ee & 1) {					/* Top */
				if (++wex[ee] >= p->n) {
					wed[ee] = GNN_INF;
					wex[ee]--;
				} else {
					double ww = p->sax[ee][wex[ee]]->mix[0][ff] - q[ff];
					wed[ee] = fabs(ww) * ww;
				}
			} else {
				if (--wex[ee] < 0) {
					wed[ee] = GNN_INF;
					wex[ee]++;
				} else {
					double ww = q[ff] - p->sax[ee][wex[ee]]->mix[1][ff];
					wed[ee] = fabs(ww) * ww;
				}
			}
		}

//printf("Searched %d points out of %d = %f%%\n",ptested, p->n, 100.0 * ptested/p->n);

		p->tbase += 3;		/* Next touch */

		rout[0] = out[0];	/* Copy results to output */
		rout[1] = out[1];
		rout[2] = out[2];
		return;
	}
}


/* Setup the nearest function acceleration structure */
static void
init_ne(
gamut *s
) {
	gnn *p;
	int i, k;
	gtri *tp;		/* Triangle pointer */
	int ntris;

//printf("~1 init_ne called\n");

	/* Allocate the nearest neighbor acceleration structure */
	if ((s->nns = p = (gnn *) calloc(1, sizeof(gnn))) == NULL) {
		fprintf(stderr,"gamut: calloc failed - gnn structure\n");
		exit(-1);
	}

	/* Count triangles */
	ntris = 0;
	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		ntris++;
	} END_FOR_ALL_ITEMS(tp);

	p->n = ntris;
	p->tbase = 0;		/* Initialse touch flag */

	/* Allocate the arrays spaces */
	for (k = 0; k < (3 * 2); k++) {
		if ((p->sax[k] = (gtri **)malloc(sizeof(gtri *) * ntris)) == NULL)
			error("Failed to allocate sorted index array");
	}

	/* For each triangle, create the triangle bounding box values, */
	/* and add them to the axis lists. */
	tp = s->tris; 
	i = 0;
	FOR_ALL_ITEMS(gtri, tp) {
		int j;
		for (j = 0; j < 3; j++) {	/* Init */
			tp->mix[0][j] = 1e38;
			tp->mix[1][j] = -1e38;
		}
		for (k = 0; k < 3; k++) {
			for (j = 0; j < 3; j++) {
				if (tp->v[k]->p[j] < tp->mix[0][j])		/* New min */
					tp->mix[0][j] = tp->v[k]->p[j];
				if (tp->v[k]->p[j] > tp->mix[1][j])		/* New max */
					tp->mix[1][j] = tp->v[k]->p[j];
			}
			p->sax[k * 2 + 0][i] = tp;
			p->sax[k * 2 + 1][i] = tp;
		}
		i++;
	} END_FOR_ALL_ITEMS(tp);


	/* Sort the axis arrays */
	for (k = 0; k < 3; k++) {

			/* Sort upper edges of bounding box */
#define 	HEAP_COMPARE(A,B) (A->mix[1][k] < B->mix[1][k])
			HEAPSORT(gtri *, &p->sax[k * 2 + 0][0], ntris)
#undef HEAP_COMPARE

			/* Sort lower edges of bounding box */
#define 	HEAP_COMPARE(A,B) (A->mix[0][k] < B->mix[0][k])
			HEAPSORT(gtri *, &p->sax[k * 2 + 1][0], ntris)
#undef HEAP_COMPARE
	}
	s->ne_inited = 1;

//printf("~1 init_ne done\n");
}

/* Free everything */
static void del_gnn(gnn *p) {
	int k;

	for (k = 0; k < (3 * 2); k++) {
		free (p->sax[k]);
	}

	free(p);
}

/* ===================================================== */
/* Define the colorspaces white and black point. May be NULL if unknown. */
/* Note that as in all of the gamut library, we assume that we are in */
/* an L*a*b* or Jab type color space. */
static void setwb(
gamut *s,
double *wp,
double *bp
) {
	if (wp != NULL) {
		s->cs_wp[0] = wp[0];
		s->cs_wp[1] = wp[1];
		s->cs_wp[2] = wp[2];
	} else {
		s->cs_wp[0] = 100.0;
		s->cs_wp[1] = 0.0;
		s->cs_wp[2] = 0.0;
	}

	if (bp != NULL) {
		s->cs_bp[0] = bp[0];
		s->cs_bp[1] = bp[1];
		s->cs_bp[2] = bp[2];
	} else {
		s->cs_bp[0] = 0.0;
		s->cs_bp[1] = 0.0;
		s->cs_bp[2] = 0.0;
	}

	s->cswbset = 1;
}


/* Compute the gamut white/black points, assuming */
/* that the colorspace white/black points have been set */
/* The gamut white/black are the points on the colorspace */
/* white/black axis that have the same L values as the */
/* extremes within the gamut. */
static void compgawb(gamut *s) {
	int i;
	double ff, Lmax, Lmin;

	if (s->cswbset == 0 || s->gawbset != 0)
		return;		/* Nothing to do */

	Lmax = -1000.0;
	Lmin =  1000.0;

	/* Discover min and max L values */
	for (i = 0; i < s->nv; i++) {
		if ((s->verts[i]->f & GVERT_SET) == 0 )
			continue;
	
		if (s->verts[i]->p[0] > Lmax)
			Lmax = s->verts[i]->p[0];
		if (s->verts[i]->p[0] < Lmin)
			Lmin = s->verts[i]->p[0];
	}

	if (Lmax > s->cs_wp[0])		/* Strange */
		Lmax = s->cs_wp[0];
	if (Lmin < s->cs_bp[0])		/* Also strange */
		Lmin = s->cs_bp[0];

	/* Locate points along colorspace grey axis */
	/* that correspond to the L extremes */
	ff = (Lmax - s->cs_bp[0])/(s->cs_wp[0] - s->cs_bp[0]);
	s->ga_wp[0] = Lmax;
	s->ga_wp[1] = ff * (s->cs_wp[1] - s->cs_bp[1]) + s->cs_bp[1];
	s->ga_wp[2] = ff * (s->cs_wp[2] - s->cs_bp[2]) + s->cs_bp[2];

	ff = (Lmin - s->cs_bp[0])/(s->cs_wp[0] - s->cs_bp[0]);
	s->ga_bp[0] = Lmin;
	s->ga_bp[1] = ff * (s->cs_wp[1] - s->cs_bp[1]) + s->cs_bp[1];
	s->ga_bp[2] = ff * (s->cs_wp[2] - s->cs_bp[2]) + s->cs_bp[2];
	
	s->gawbset = 1;
}

/* Get the colorspace and gamut white & black points. */
/* Return pointers may be NULL */
/* Return non-zero if not possible. */
static int getwb(
gamut *s,
double *cswp,
double *csbp,
double *gawp,
double *gabp
) {
	if (s->cswbset == 0) {
		return 1;
	}

	if (cswp != NULL) {
		cswp[0] = s->cs_wp[0];
		cswp[1] = s->cs_wp[1];
		cswp[2] = s->cs_wp[2];
	}
	
	if (csbp != NULL) {
		csbp[0] = s->cs_bp[0];
		csbp[1] = s->cs_bp[1];
		csbp[2] = s->cs_bp[2];
	}
	
	if (gawp != NULL || gabp != NULL)
		compgawb(s);		/* make sure we have gamut white/black available */

	if (gawp != NULL) {
		gawp[0] = s->ga_wp[0];
		gawp[1] = s->ga_wp[1];
		gawp[2] = s->ga_wp[2];
	}

	if (gabp != NULL) {
		gabp[0] = s->ga_bp[0];
		gabp[1] = s->ga_bp[1];
		gabp[2] = s->ga_bp[2];
	}

	return 0;
}


/* ---------------------------------------------------- */
/* Per-triangle primitives used to compute radial & vector */
/* intersection, and nearest point. */
/* See if the given triangle intersect the given vector. */
/* Return 1 if it doesm 0 if it doesn't */
static int vect_intersect(
gamut *s,
double *rvp,	/* parameter, 0.0 = p1, 1.0 = p2 */
double *ip,		/* return intersection point */
double *p1,		/* First point of vector (ie black) */
double *p2,		/* Second point of vector (ie white) */
gtri *t			/* Triangle in question */
) {
	double rv;			/* Axis parameter value */
	double gv[3];		/* Grey axis vector */
	double ival[3];		/* Intersection value */
	double den;
	int j;

	gv[0] = p2[0] - p1[0];
	gv[1] = p2[1] - p1[1];
	gv[2] = p2[2] - p1[2];

	den = t->pe[0] * gv[0] + t->pe[1] * gv[1] + t->pe[2] * gv[2];
	if (fabs(den) < 1e-10) {
		return 0;
	}

	/* Compute the intersection of the grey axis vector with the triangle plane */
	rv = -(t->pe[0] * p1[0] + t->pe[1] * p1[1] + t->pe[2] * p1[2] + t->pe[3])/den;

	/* Compute the actual intersection point */
	ival[0] = p1[0] + rv * gv[0];
	ival[1] = p1[1] + rv * gv[1];
	ival[2] = p1[2] + rv * gv[2];

	/* Check if the intersection point is within the triangle */
	for (j = 0; j < 3; j++) {
		double ds;
		ds = t->ee[j][0] * (ival[0] - s->cent[0])	/* Convert to relative for edge check */
		   + t->ee[j][1] * (ival[1] - s->cent[1])
	       + t->ee[j][2] * (ival[2] - s->cent[2])
		   + t->ee[j][3];
		if (ds > 1e-8) {
			return 0;		/* Not within triangle */
		}
	}

	/* Got an intersection point */
	ip[0] = ival[0];
	ip[1] = ival[1];
	ip[2] = ival[2];

	*rvp = rv;

	return 1;
}

/* Given a point and a triangle, return the closest point on */
/* the triangle closest to point. Also return the distance squared */
static double ne_point_on_tri(
gamut *s,
gtri *t,		/* Triangle to use */
double *out,	/* Absolute output point */
double *in		/* Absolute input point */
) {
	int j;
	double rv;
	double bdist;

	/* Compute the point on the triangles plane, that is orthogonal */
	/* (closest) to the target point. */
	rv = (t->pe[0] * in[0] + t->pe[1] * in[1] + t->pe[2] * in[2] + t->pe[3])/
	     (t->pe[0] * t->pe[0] + t->pe[1] * t->pe[1] + t->pe[2] * t->pe[2]);

	out[0] = in[0] - rv * t->pe[0];
	out[1] = in[1] - rv * t->pe[1];
	out[2] = in[2] - rv * t->pe[2];

	/* Check if the closest point is within this triangle */
	for (j = 0; j < 3; j++) {
		double ds;
		ds = t->ee[j][0] * (out[0] - s->cent[0])	/* Convert to relative for edge check */
		   + t->ee[j][1] * (out[1] - s->cent[1])
	       + t->ee[j][2] * (out[2] - s->cent[2])
		   + t->ee[j][3];
		if (ds > 1e-8) {
			break;		/* Not within triangle */
		}
	}
	if (j >= 3) {	/* It's OK */
		return rv * rv;	/* rv is distance since pe length is 1.0 */
	}

	/* Not in triangle, so find closest point along any edge, */
	/* or at the verticies. */
	bdist = 1e38;
	for (j = 0; j < 3; j++) {	/* For each edge */
		gedge *e = t->e[j];
		int k;
		double nu, de, ds;
		for (de = 0.0, k = 0; k < 3; k++) {
			double tt = e->v[1]->p[k] - e->v[0]->p[k];
			de += tt * tt;
		}
		for (nu = 0.0, k = 0; k < 3; k++)
			nu += (e->v[1]->p[k] - e->v[0]->p[k]) * (in[k] - e->v[0]->p[k]);

		ds = nu/de;

		if (ds >= 0.0 && ds <= 1.0) {	/* Valid edge */
			double tout[3], ss;
			for (ss = 0.0, k = 0; k < 3; k++) {
				tout[k] = e->v[0]->p[k] + ds * (e->v[1]->p[k] - e->v[0]->p[k]);
				ss += (in[k] - tout[k]) * (in[k] - tout[k]); 
			}
			if (ss < bdist) {
				bdist = ss;
				out[0] = tout[0];
				out[1] = tout[1];
				out[2] = tout[2];
			}
		} 
	}

	for (j = 0; j < 3; j++) {	/* For each vertex */
		int k;
		double ss;
		for (ss = 0.0, k = 0; k < 3; k++) {
			double tt;
			tt = in[k] - t->v[j]->p[k];
			ss += tt * tt;
		}

		if (ss < bdist) {
			bdist = ss;
			out[0] = t->v[j]->p[0];
			out[1] = t->v[j]->p[1];
			out[2] = t->v[j]->p[2];
		} 
	}

	return bdist;
}

/* ----------------------------------- */

#ifdef NEVER	/* Alternate brute force radial() */
/* Given a point, return the point in that direction */
/* that lies on the gamut surface */
/* Return the radial radius to the surface point */
/* We use a simple exaustive search, so this is not very fast. */
static double
radial(
gamut *s,
double *out,	/* result point (absolute)*/
double *in		/* input point (absolute)*/
) {
	gtri *tp;
	int j;
	double ss, rv;
	double nin[3];	/* Normalised input vector */

//printf("~1 radial called with %f %f %f\n", in[0], in[1], in[2]);
	if IS_LIST_EMPTY(s->tris)
		triangulate(s);

	/* Compute vector length to center point */
	for (ss = 0.0, j = 0; j < 3; j++)
		ss += (in[j] - s->cent[j]) * (in[j] - s->cent[j]);
	ss = 1.0/sqrt(ss);				/* Normalising factor */
	for (ss = 0.0, j = 0; j < 3; j++)
		nin[j] = s->cent[j] + (in[j] - s->cent[j]) * ss;

	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		if (vect_intersect(s, &rv, out, s->cent, nin, tp)) {
			if (rv > 0.0)	/* Expect only one intersection */
				break;
		}
	} END_FOR_ALL_ITEMS(tp);

//printf("~1 result = %f %f %f\n",out[0], out[1], out[2]);

	return rv;
}
#endif /* NEVER */


#ifdef NEVER	/* Alternate brute force nearest()  */
/* Given an absolute point, return the point on the gamut */
/* surface that is closest to it. */
static void
nearest(
gamut *s,
double *out,	/* result point (absolute) */
double *q		/* Target point (absolute) */
) {
	gtri *tp;
	double bdist = 1e308;	/* Best possible distance to an object outside the window */


	if IS_LIST_EMPTY(s->tris)
		triangulate(s);

	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		double r[3];		/* Possible solution point */
		double tdist;

		/* Compute distance from query point to this object */
		tdist = ne_point_on_tri(s, tp, r, q);

		if (tdist < bdist) {	/* New best point */
			bdist = tdist;
			out[0] = r[0];
			out[1] = r[1];
			out[2] = r[2];
		}
	} END_FOR_ALL_ITEMS(tp);
}
#endif /* NEVER */

/* Given a vector, find the two extreme intersection with */
/* the gamut surface. */
/* We use a simple exuastive search, so this is not very fast. */
/* Return 0 if there is no intersection */
static int compute_vector_isect(
gamut *s,
double *p1,		/* First point (ie black) */
double *p2,		/* Second point (ie white) */
double *omin,	/* Return gamut surface points, min = closest to p1 */
double *omax,	/* max = farthest from p1 */
double *omnt,	/* Return parameter values for p1 and p2, 0 being at p1, */
double *omxt	/* and 1 being at p2 */
) {
	gtri *tp;
	double ip[3], min[3], max[3];
	double mint, maxt;
	int j;

	if IS_LIST_EMPTY(s->tris)
		triangulate(s);

	maxt = -1e68;	/* Setup to find min/max */
	mint =  1e68;

	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		double rv;
		if (vect_intersect(s, &rv, ip, p1, p2, tp)) {
			if (rv > maxt) {
				max[0] = ip[0];
				max[1] = ip[1];
				max[2] = ip[2];
				maxt = rv;
			}
			if (rv < mint) {
				min[0] = ip[0];
				min[1] = ip[1];
				min[2] = ip[2];
				mint = rv;
			}
		}
	} END_FOR_ALL_ITEMS(tp);

	if (((omax != NULL || omxt != NULL) && maxt == -1e68)
	 || ((omin != NULL || omnt != NULL) && mint == 1e68)) {
		return 0;
	}

	if (omax != NULL)
		for (j = 0; j < 3; j++)
			omax[j] = max[j];

	if (omin != NULL)
		for (j = 0; j < 3; j++)
			omin[j] = min[j];

	if (omxt != NULL)
		*omxt = maxt;

	if (omnt != NULL)
		*omnt = mint;

	return 1;
}


/* ===================================================== */
/* ===================================================== */
/* ----------------------------------- */
/* Write to a VRML .wrl file */
/* Return non-zero on error */
static int write_vrml(
gamut *s,
char *filename,
int doaxes,			/* Non-zero if axes are to be written */
int docusps			/* Non-zero if cusp points are to be marked */
) {
	return write_trans_vrml(s, filename, doaxes, docusps, NULL, NULL);
}

/* Write to a VRML .wrl file */
/* Return non-zero on error */
static int write_trans_vrml(
gamut *s,
char *filename,
int doaxes,			/* Non-zero if axes are to be written */
int docusps,		/* Non-zero if cusp points are to be marked */
void (*transform)(void *cntx, double out[3], double in[3]),	/* Optional transformation callback */
void *cntx
) {
	int i;
	gtri *tp;		/* Triangle pointer */
	FILE *wrl;
	struct {
		double x, y, z;
		double wx, wy, wz;
		double r, g, b;
	} axes[5] = {
		{ 0 - s->cent[1], 0 - s->cent[2],  50 - s->cent[0], 2, 2, 100, .7, .7, .7 },
																			/* L axis */
		{ 50 - s->cent[1], 0 - s->cent[2],  0 - s->cent[0], 100, 2, 2,  1,  0,  0 },
																			/* +a (red) axis */
		{ 0 - s->cent[1], -50 - s->cent[2], 0 - s->cent[0], 2, 100, 2,  0,  0,  1 },
																			/* -b (blue) axis */
		{ -50 - s->cent[1], 0 - s->cent[2], 0 - s->cent[0], 100, 2, 2,  0,  1,  0 },
																			/* -a (green) axis */
		{ 0 - s->cent[1],  50 - s->cent[2], 0 - s->cent[0], 2, 100, 2,  1,  1,  0 },
																			/* +b (yellow) axis */
	};

	/* Define the labels */
	struct {
		double x, y, z;
		double size;
		char *string;
		double r, g, b;
	} labels[6] = {
		{ -2 - s->cent[1], 2 - s->cent[2],  - s->cent[0] + 100 + 10, 10, "+L*",  .7, .7, .7 },
																			/* Top of L axis */
		{ -2 - s->cent[1], 2 - s->cent[2],  - s->cent[0] - 10,      10, "0",    .7, .7, .7 },
																			/* Bottom of L axis */
		{ 100 + 5 - s->cent[1], -3 - s->cent[2],  0 - s->cent[0],  10, "+a*",  1,  0,  0 },
																			/* +a (red) axis */
		{ -5 - s->cent[1], -100 - 10 - s->cent[2], 0 - s->cent[0],  10, "-b*",  0,  0,  1 },
																			/* -b (blue) axis */
		{ -100 - 15 - s->cent[1], -3 - s->cent[2], 0 - s->cent[0],  10, "-a*",  0,  0,  1 },
																			/* -a (green) axis */
		{ -5 - s->cent[1],  100 + 5 - s->cent[2], 0 - s->cent[0],  10, "+b*",  1,  1,  0 },
																			/* +b (yellow) axis */
	};

	if IS_LIST_EMPTY(s->tris)
		triangulate(s);
		
	if ((wrl = fopen(filename,"w")) == NULL) {
		fprintf(stderr,"Error opening output file '%s'\n",filename);
		return 2;
	}

	/* Spit out a VRML 2 Object surface of gamut */

	fprintf(wrl,"#VRML V2.0 utf8\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"# Created by the Argyll CMS\n");
	fprintf(wrl,"Transform {\n");
  	fprintf(wrl,"children [\n");
    fprintf(wrl,"	NavigationInfo {\n");
	fprintf(wrl,"		type \"EXAMINE\"        # It's an object we examine\n");
	fprintf(wrl,"	} # We'll add our own light\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    DirectionalLight {\n");
	fprintf(wrl,"        intensity 0.2\n");
	fprintf(wrl,"        ambientIntensity 0.1\n");
	fprintf(wrl,"        direction -1 -1 -1\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"    DirectionalLight {\n");
	fprintf(wrl,"        intensity 0.6\n");
	fprintf(wrl,"        ambientIntensity 0.2\n");
	fprintf(wrl,"        direction 1 1 1\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    Viewpoint {\n");
	fprintf(wrl,"        position 0 0 250      # Position we view from\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");
	if (doaxes != 0) {
		fprintf(wrl,"# Lab axes as boxes:\n");
		for (i = 0; i < 5; i++) {
			fprintf(wrl,"Transform { translation %f %f %f\n", axes[i].x, axes[i].y, axes[i].z);
			fprintf(wrl,"\tchildren [\n");
			fprintf(wrl,"\t\tShape {\n");
			fprintf(wrl,"\t\t\tgeometry Box { size %f %f %f }\n",
			                  axes[i].wx, axes[i].wy, axes[i].wz);
			fprintf(wrl,"\t\t\tappearance Appearance { material Material ");
			fprintf(wrl,"{ diffuseColor %f %f %f} }\n", axes[i].r, axes[i].g, axes[i].b);
			fprintf(wrl,"\t\t}\n");
			fprintf(wrl,"\t]\n");
			fprintf(wrl,"}\n");
		}
		fprintf(wrl,"# Axes identification:\n");
		for (i = 0; i < 6; i++) {
			fprintf(wrl,"Transform { translation %f %f %f\n", labels[i].x, labels[i].y, labels[i].z);
			fprintf(wrl,"\tchildren [\n");
			fprintf(wrl,"\t\tShape {\n");
			fprintf(wrl,"\t\t\tgeometry Text { string [\"%s\"]\n",labels[i].string);
			fprintf(wrl,"\t\t\t\tfontStyle FontStyle { family \"SANS\" style \"BOLD\" size %f }\n",
			                                            labels[i].size);
			fprintf(wrl,"\t\t\t\t}\n");
			fprintf(wrl,"\t\t\tappearance Appearance { material Material ");
			fprintf(wrl,"{ diffuseColor %f %f %f} }\n", labels[i].r, labels[i].g, labels[i].b);
			fprintf(wrl,"\t\t}\n");
			fprintf(wrl,"\t]\n");
			fprintf(wrl,"}\n");
		}
		fprintf(wrl,"\n");
	}
	fprintf(wrl,"    Transform {\n");
	fprintf(wrl,"      translation 0 0 0\n");
	fprintf(wrl,"      children [\n");
	fprintf(wrl,"		Shape { \n");
	fprintf(wrl,"		    geometry IndexedFaceSet {\n");
	fprintf(wrl,"				ccw FALSE\n");
	fprintf(wrl,"				convex TRUE\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"		        coord Coordinate { \n");
	fprintf(wrl,"		            point [			# Verticy coordinates\n");

	/* Spit out the point values, in order. */
	/* Note that a->x, b->y, L->z */
	for (i = 0; i < s->nv; i++) {
		double out[3];

#ifdef SHOW_BUCKETS		/* Show vertex buckets as surface */
		if (!(s->verts[i]->f & GVERT_SET))
#else
		if (!(s->verts[i]->f & GVERT_TRI))
#endif
			continue;

#ifdef SHOW_BUCKETS		/* Show vertex buckets as surface */
		{
			double cc[3], rr[3];
# ifdef SHOW_SPHERE			/* Show surface on sphere */
			rr[0] = 50.0;	/* Sphere radius */
# else
			rr[0] = s->verts[i]->r[0],
# endif /* SHOW_SPHERE */

			rr[1] = s->verts[i]->hc - 0.5 * s->verts[i]->w;
			rr[2] = s->verts[i]->vc - 0.5 * s->verts[i]->h;
			gamut_radial2rect(s, cc, rr);
			fprintf(wrl,"%f %f %f,\n",cc[1], cc[2], cc[0]);
	
			rr[1] = s->verts[i]->hc - 0.5 * s->verts[i]->w;
			rr[2] = s->verts[i]->vc + 0.5 * s->verts[i]->h;
			gamut_radial2rect(s, cc, rr);
			fprintf(wrl,"%f %f %f,\n",cc[1], cc[2], cc[0]);
	
			rr[1] = s->verts[i]->hc + 0.5 * s->verts[i]->w;
			rr[2] = s->verts[i]->vc + 0.5 * s->verts[i]->h;
			gamut_radial2rect(s, cc, rr);
			fprintf(wrl,"%f %f %f,\n",cc[1], cc[2], cc[0]);
	
			rr[1] = s->verts[i]->hc + 0.5 * s->verts[i]->w;
			rr[2] = s->verts[i]->vc - 0.5 * s->verts[i]->h;
			gamut_radial2rect(s, cc, rr);
			fprintf(wrl,"%f %f %f,\n",cc[1], cc[2], cc[0]);
		}

#else	/* Show point data */

# ifdef SHOW_SPHERE			/* Show surface on sphere */
		fprintf(wrl,"%f %f %f,\n",s->verts[i]->sp[1], s->verts[i]->sp[2],
		                           s->verts[i]->sp[0]);
# else
# ifdef SHOW_HULL_PNTS
		fprintf(wrl,"%f %f %f,\n",s->verts[i]->ch[1], s->verts[i]->ch[2],
		                           s->verts[i]->ch[0]);
# else
		/* Show normal gamut surface */
		out[0] = s->verts[i]->p[0];
		out[1] = s->verts[i]->p[1];
		out[2] = s->verts[i]->p[2];

		if (transform)
			transform(cntx, out, out);		/* Do transform */

		fprintf(wrl,"%f %f %f,\n",out[1]-s->cent[1], out[2]-s->cent[2], out[0]-s->cent[0]);

# endif /* SHOW_HULL_PNTS */
# endif /* SHOW_SPHERE */

#endif /* SHOW_BUCKETS */

	}
	fprintf(wrl,"					]\n");
	fprintf(wrl,"		        }\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"		        coordIndex [ 		# Indexes of poligon Verticies \n");

#ifdef SHOW_BUCKETS		/* Show vertex buckets as surface */
	for (i = 0; i < s->nv; i++) {
		int j = s->verts[i]->sn;
		if (!(s->verts[i]->f & GVERT_SET))
			continue;
		fprintf(wrl,"%d, %d, %d, %d, -1\n", j * 4, j * 4 + 1, j * 4 + 2, j * 4 + 3);
	}
#else	/* Show gamut triangular surface */
	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		fprintf(wrl,"%d, %d, %d, -1\n", tp->v[0]->tn, tp->v[1]->tn, tp->v[2]->tn);
	} END_FOR_ALL_ITEMS(tp);
#endif /* SHOW_BUCKETS */

	fprintf(wrl,"				]\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"				colorPerVertex TRUE\n");
	fprintf(wrl,"		        color Color {\n");
	fprintf(wrl,"		            color [			# RGB colors of each vertex\n");

	/* Spit out the colors for each vertex */
	for (i = 0; i < s->nv; i++) {
		double rgb[3];
#ifdef SHOW_BUCKETS		/* Show vertex buckets as surface */
		if (!(s->verts[i]->f & GVERT_SET))
#else
		if (!(s->verts[i]->f & GVERT_TRI))
#endif
			continue;

#ifdef COLORED_VRML
		gamut_Lab2RGB(rgb, s->verts[i]->p);
#else
		rgb[0] = rgb[1] = rgb[2] = 1.0;
#endif
		fprintf(wrl,"%f %f %f,\n", rgb[0], rgb[1], rgb[2]);
#ifdef SHOW_BUCKETS		/* Show vertex buckets as surface */
		fprintf(wrl,"%f %f %f,\n", rgb[0], rgb[1], rgb[2]);
		fprintf(wrl,"%f %f %f,\n", rgb[0], rgb[1], rgb[2]);
		fprintf(wrl,"%f %f %f,\n", rgb[0], rgb[1], rgb[2]);
#endif /* SHOW_BUCKETS */
	}
	fprintf(wrl,"					] \n");
	fprintf(wrl,"		        }\n");
	fprintf(wrl,"		    }\n");
	fprintf(wrl,"		    appearance Appearance { \n");
	fprintf(wrl,"		        material Material {\n");
	fprintf(wrl,"					transparency 0.0\n");
	fprintf(wrl,"					ambientIntensity 0.3\n");
	fprintf(wrl,"					shininess 0.5\n");
	fprintf(wrl,"				}\n");
	fprintf(wrl,"		    }\n");
	fprintf(wrl,"		}	# end Shape\n");
	fprintf(wrl,"      ]\n");
	fprintf(wrl,"    }\n");


	if (s->gawbset && doaxes) {

		/* Show the gamut white and black points */
		fprintf(wrl,"\n");
		fprintf(wrl,"    Transform {\n");
		fprintf(wrl,"      translation %f %f %f\n",s->ga_wp[1]-s->cent[1], s->ga_wp[2]-s->cent[2], s->ga_wp[0]-s->cent[0]);
		fprintf(wrl,"      children [\n");
		fprintf(wrl,"		Shape { \n");
		fprintf(wrl,"		 geometry Sphere { radius 2.0 }\n");
		fprintf(wrl,"        appearance Appearance { material Material { diffuseColor 0.9 0.9 0.9 } }\n");
		fprintf(wrl,"		} \n");
		fprintf(wrl,"      ]\n");
		fprintf(wrl,"    }\n");
		fprintf(wrl,"\n");
		fprintf(wrl,"    Transform {\n");
		fprintf(wrl,"      translation %f %f %f\n",s->ga_bp[1]-s->cent[1], s->ga_bp[2]-s->cent[2], s->ga_bp[0]-s->cent[0]);
		fprintf(wrl,"      children [\n");
		fprintf(wrl,"		Shape { \n");
		fprintf(wrl,"		 geometry Sphere { radius 2.0 }\n");
		fprintf(wrl,"        appearance Appearance { material Material { diffuseColor 0.9 0.9 0.9 } }\n");
		fprintf(wrl,"		} \n");
		fprintf(wrl,"      ]\n");
		fprintf(wrl,"    }\n");
	}

	if (docusps && s->cu_inited != 0) {
		double ccolors[6][3] = {
			{ 1.0, 0.1, 0.1 },		/* Red */
			{ 1.0, 1.0, 0.1 },		/* Yellow */
			{ 0.1, 1.0, 0.1 },		/* Green */
			{ 0.1, 1.0, 1.0 },		/* Cyan */
			{ 0.1, 0.1, 1.0 },		/* Blue */
			{ 1.0, 0.1, 1.0 }		/* Magenta */
		};

		for (i = 0; i < 6; i++) {
			fprintf(wrl,"\n");
			fprintf(wrl,"    Transform {\n");
			fprintf(wrl,"      translation %f %f %f\n",s->cusps[i][1]-s->cent[1], s->cusps[i][2]-s->cent[2], s->cusps[i][0]-s->cent[0]);
			fprintf(wrl,"      children [\n");
			fprintf(wrl,"		Shape { \n");
			fprintf(wrl,"		 geometry Sphere { radius 2.0 }\n");
			fprintf(wrl,"        appearance Appearance { material Material { diffuseColor %f %f %f } }\n", ccolors[i][0],ccolors[i][1],ccolors[i][2]);
			fprintf(wrl,"		} \n");
			fprintf(wrl,"      ]\n");
			fprintf(wrl,"    }\n");
		}
	}

#ifdef TEST_LOOKUP
	{
		int i, j;
		double in[3], out[3];

		fprintf(wrl,"\n");
		fprintf(wrl,"Shape {\n");
		fprintf(wrl,"  geometry PointSet { \n");
		fprintf(wrl,"    coord Coordinate { \n");
		fprintf(wrl,"	   point [\n");

		for (i = 0; i < 10; i++) {
			double ss;
			/* Create random vector relative to center, absolute */
			in[0] = (rand() / (double)RAND_MAX) - 0.5 + s->cent[0];
			in[1] = (rand() / (double)RAND_MAX) - 0.5 + s->cent[1];
			in[2] = (rand() / (double)RAND_MAX) - 0.5 + s->cent[2];
				
			s->radial(s, out, in);	/* Lookup point on gamut surface */

			out[0] = (out[0] - s->cent[0]) * 1.01 + s->cent[0];
			out[1] = (out[1] - s->cent[1]) * 1.01 + s->cent[1];
			out[2] = (out[2] - s->cent[2]) * 1.01 + s->cent[2];
			fprintf(wrl,"%f %f %f,\n",out[1]-s->cent[1], out[2]-s->cent[2], out[0]-s->cent[0]);
		}
		fprintf(wrl,"      ]\n");
		fprintf(wrl,"    }\n");
		fprintf(wrl,"  }\n");
		fprintf(wrl,"} # end shape\n");

	}
#endif	/* TEST_LOOKUP */

#ifdef TEST_NEAREST
	{
#define NTPTS 500
		int i, j;
		double in[3], out[3];

		fprintf(wrl,"\n");
		fprintf(wrl,"Shape {\n");
		fprintf(wrl,"  geometry IndexedLineSet { \n");
		fprintf(wrl,"    coord Coordinate { \n");
		fprintf(wrl,"	   point [\n");

		for (i = 0; i < NTPTS; i++) {
			double ss;
			/* Create random vector relative to center */
			in[0] = (rand() / (double)RAND_MAX) - 0.5;
			in[1] = (rand() / (double)RAND_MAX) - 0.5;
			in[2] = (rand() / (double)RAND_MAX) - 0.5;

#ifndef NEVER	/* Make points just above surface */
			in[0] += s->cent[0];			/* Make absolute */
			in[1] += s->cent[1];
			in[2] += s->cent[2];
			s->radial(s, in, in);	/* Lookup point on gamut surface */
			in[0] = (in[0] - s->cent[0]) * 1.20 + s->cent[0];	/* Extend by 10% */
			in[1] = (in[1] - s->cent[1]) * 1.20 + s->cent[1];
			in[2] = (in[2] - s->cent[2]) * 1.20 + s->cent[2];
#else
			/* Make distance 150 */
			ss = sqrt(in[0] * in[0] + in[1] * in[1] + in[2] * in[2]);
			in[0] = 60.0/ss * in[0] + s->cent[0];
			in[1] = 60.0/ss * in[1] + s->cent[1];
			in[2] = 60.0/ss * in[2] + s->cent[2];
#endif

//			s->radial(s, out, in);	/* Lookup point on gamut surface */

			s->nearest(s, out, in);	/* Nearest point on gamut surface */

			fprintf(wrl,"%f %f %f,\n",in[1]-s->cent[1], in[2]-s->cent[2], in[0]-s->cent[0]);
			fprintf(wrl,"%f %f %f,\n",out[1]-s->cent[1], out[2]-s->cent[2], out[0]-s->cent[0]);
		}

		fprintf(wrl,"      ]\n");
		fprintf(wrl,"    }\n");
		fprintf(wrl,"  coordIndex [\n");

		for (i = 0; i < NTPTS; i++) {
			fprintf(wrl,"%d, %d, -1,\n", i * 2, i * 2 + 1);
		}
		fprintf(wrl,"    ]\n");
		fprintf(wrl,"  }\n");
		fprintf(wrl,"} # end shape\n");

	}
#endif	/* TEST_NEAREST */

	fprintf(wrl,"\n");
	fprintf(wrl,"  ] # end of children for world\n");
	fprintf(wrl,"}\n");

	if (fclose(wrl) != 0) {
		fprintf(stderr,"Error closing output file '%s'\n",filename);
		return 2;
	}

	return 0;
}


/* ----------------------------------- */
/* Write to a CGATS .gam file */
/* Return non-zero on error */
static int write_gam(
gamut *s,
char *filename
) {
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	int i;
	gtri *tp;		/* Triangle pointer */
	cgats *gam;
	char buf[100];

	if IS_LIST_EMPTY(s->tris)
		triangulate(s);
		
	gam = new_cgats();	/* Create a CGATS structure */
	gam->add_other(gam, "GAMUT");

	gam->add_table(gam, tt_other, 0);	/* Start the first table as type "GAMUT" */

	gam->add_kword(gam, 0, "DESCRIPTOR", "Argyll Gamut surface poligon data", NULL);
	gam->add_kword(gam, 0, "ORIGINATOR", "Argyll CMS gamut library", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	gam->add_kword(gam, 0, "CREATED",atm, NULL);

#ifdef NEVER
	/* would be nice to add extra info like description, source (ie icc filename) etc. */
	gam->add_kword(gam, 0, "DEVICE_CLASS","INPUT", NULL);	/* What sort of device this is */
#endif
	if (s->isJab)
		gam->add_kword(gam, 0, "COLOR_REP","JAB", NULL);
	else
		gam->add_kword(gam, 0, "COLOR_REP","LAB", NULL);

	sprintf(buf,"%f %f %f", s->cent[0], s->cent[1], s->cent[2]);
	gam->add_kword(gam, 0, "GAMUT_CENTER",buf, NULL);

	/* If the white and black points are known, put them in the file */
	if (s->cswbset) {

		compgawb(s);		/* make sure we have gamut white/black available */

		sprintf(buf,"%f %f %f", s->cs_wp[0], s->cs_wp[1], s->cs_wp[2]);
		gam->add_kword(gam, 0, "CSPACE_WHITE",buf, NULL);

		sprintf(buf,"%f %f %f", s->ga_wp[0], s->ga_wp[1], s->ga_wp[2]);
		gam->add_kword(gam, 0, "GAMUT_WHITE",buf, NULL);

		sprintf(buf,"%f %f %f", s->cs_bp[0], s->cs_bp[1], s->cs_bp[2]);
		gam->add_kword(gam, 0, "CSPACE_BLACK",buf, NULL);

		sprintf(buf,"%f %f %f", s->ga_bp[0], s->ga_bp[1], s->ga_bp[2]);
		gam->add_kword(gam, 0, "GAMUT_BLACK",buf, NULL);
	}

	/* If cusp values are known, put them in the file */
	if (s->cu_inited != 0) {
		char buf1[50], buf2[100];
		char *cnames[6] = { "RED", "YELLOW", "GREEN", "CYAN", "BLUE", "MAGENTA" };

		for (i = 0; i < 6; i++) {
			sprintf(buf1,"CUSP_%s", cnames[i]);
			sprintf(buf2,"%f %f %f", s->cusps[i][0], s->cusps[i][1], s->cusps[i][2]);
			gam->add_kword(gam, 0, buf1, buf2, NULL);
		}
	}

	gam->add_kword(gam, 0, NULL, NULL, "First come the triangle verticy location");

	gam->add_field(gam, 0, "VERTEX_NO", i_t);
	gam->add_field(gam, 0, "LAB_L", r_t);
	gam->add_field(gam, 0, "LAB_A", r_t);
	gam->add_field(gam, 0, "LAB_B", r_t);

	/* Spit out the vertex values, in order. */
	for (i = 0; i < s->nv; i++) {
		if (!(s->verts[i]->f & GVERT_TRI))
			continue;
		gam->add_set(gam, 0, s->verts[i]->tn,
		             s->verts[i]->p[0], s->verts[i]->p[1], s->verts[i]->p[2]);
	}

	gam->add_table(gam, tt_other, 0);	/* Start the second table */
	gam->set_table_flags(gam, 1, 1, 1, 0);	/* Suppress id & kwords */
	gam->add_kword(gam, 1, NULL, NULL, "And then come the triangles");

	gam->add_field(gam, 1, "VERTEX_0", i_t);
	gam->add_field(gam, 1, "VERTEX_1", i_t);
	gam->add_field(gam, 1, "VERTEX_2", i_t);

	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		gam->add_set(gam, 1, tp->v[0]->tn, tp->v[1]->tn, tp->v[2]->tn);
	} END_FOR_ALL_ITEMS(tp);


	if (gam->write_name(gam, filename)) {
		fprintf(stderr,"Error writing to file '%s' : '%s'\n",filename, gam->err);
		return 2;
	}

	gam->del(gam);		/* Clean up */
	return 0;
}

/* ----------------------------------- */
/* Read from a CGATS .gam file */
/* Return non-zero on error */
static int read_gam(
gamut *s,
char *filename
) {
	int i;
	cgats *gam;
	gtri *tp;
	int nverts;
	int ntris;
	int Lf, af, bf;			/* Fields holding L, a & b data */
	int v0f, v1f, v2f;		/* Fields holding verticies 0, 1 & 2 */
	int cw, cb;				/* Colorspace white, black keyword indexes */
	int gw, gb;				/* Gamut white, black keyword indexes */

	if (s->tris != NULL || s->read_inited || s->lu_inited || s->ne_inited) {
		fprintf(stderr,"Can't add read into gamut after it is initialised!\n");
		return 1;
	}

	gam = new_cgats();	/* Create a CGATS structure */

	gam->add_other(gam, "GAMUT");		/* Setup to cope with a gamut file */

	if (gam->read_name(gam, filename)) {
		fprintf(stderr,"Input file '%s' error : %s",filename, gam->err);
		return 1;
	}

	if (gam->t[0].tt != tt_other || gam->t[0].oi != 0) {
		fprintf(stderr,"Input file isn't a GAMUT format file");
		return 1;
	}
	if (gam->ntables != 2) {
		fprintf(stderr,"Input file doesn't contain exactly two tables");
		return 1;
	}

	/* Figure the basic colorspace information */
	s->isJab = 0;
	if ((cw = gam->find_kword(gam, 0, "COLOR_REP")) >= 0) {
		if (strcmp(gam->t[0].kdata[cw], "JAB") == 0)
			s->isJab = 1;
	}

	/* If we can find the the colorspace white and black points, add them to the gamut */
	cw = gam->find_kword(gam, 0, "CSPACE_WHITE");
	cb = gam->find_kword(gam, 0, "CSPACE_BLACK");
	if (cw >= 0 && cb >= 0) {
		int ok = 1;
		if (sscanf(gam->t[0].kdata[cw], "%lf %lf %lf",
		           &s->cs_wp[0], &s->cs_wp[1], &s->cs_wp[2]) != 3) {
			ok = 0;
		}

		if (sscanf(gam->t[0].kdata[cb], "%lf %lf %lf",
		           &s->cs_bp[0], &s->cs_bp[1], &s->cs_bp[2]) != 3) {
			ok = 0;
		}

		if (ok) {
			s->cswbset = 1;
		}
	}

	/* If we can find the the gamut white and black points, add them to the gamut */
	gw = gam->find_kword(gam, 0, "GAMUT_WHITE");
	gb = gam->find_kword(gam, 0, "GAMUT_BLACK");
	if (gw >= 0 && gb >= 0) {
		int ok = 1;
		if (sscanf(gam->t[0].kdata[gw], "%lf %lf %lf",
		           &s->ga_wp[0], &s->ga_wp[1], &s->ga_wp[2]) != 3) {
			ok = 0;
		}

		if (sscanf(gam->t[0].kdata[gb], "%lf %lf %lf",
		           &s->ga_bp[0], &s->ga_bp[1], &s->ga_bp[2]) != 3) {
			ok = 0;
		}

		if (ok) {
			s->gawbset = 1;
		}
	}

	/* See if there are cusp values */
	{
		int kk;
		char buf1[50];
		char *cnames[6] = { "RED", "YELLOW", "GREEN", "CYAN", "BLUE", "MAGENTA" };

		for (i = 0; i < 6; i++) {
			sprintf(buf1,"CUSP_%s", cnames[i]);
			if ((kk = gam->find_kword(gam, 0, buf1)) < 0)
				break;

			if (sscanf(gam->t[0].kdata[kk], "%lf %lf %lf",
		           &s->cusps[i][0], &s->cusps[i][1], &s->cusps[i][2]) != 3) {
				break;
			}
		}
		if (i >= 6)
			s->cu_inited = 1;
	}


	if ((nverts = gam->t[0].nsets) <= 0) {
		fprintf(stderr,"No verticies");
		return 1;
	}
	if ((ntris = gam->t[1].nsets) <= 0) {
		fprintf(stderr,"No triangles");
		return 1;
	}

	/* Get ready to read the verticy data */
	if ((Lf = gam->find_field(gam, 0, "LAB_L")) < 0) {
		fprintf(stderr,"Input file doesn't contain field LAB_L");
		return 1;
	}
	if (gam->t[0].ftype[Lf] != r_t) {
		fprintf(stderr,"Field LAB_L is wrong type");
		return 1;
	}
	if ((af = gam->find_field(gam, 0, "LAB_A")) < 0) {
		fprintf(stderr,"Input file doesn't contain field LAB_A");
		return 1;
	}
	if (gam->t[0].ftype[af] != r_t) {
		fprintf(stderr,"Field LAB_A is wrong type");
		return 1;
	}
	if ((bf = gam->find_field(gam, 0, "LAB_B")) < 0) {
		fprintf(stderr,"Input file doesn't contain field LAB_B");
		return 1;
	}
	if (gam->t[0].ftype[bf] != r_t) {
		fprintf(stderr,"Field LAB_B is wrong type");
		return 1;
	}

	/* Allocate an array to point at the verts */
	if ((s->verts = (gvert **)malloc(nverts * sizeof(gvert *))) == NULL) {
		fprintf(stderr,"gamut: malloc failed on gvert pointer\n");
		return 2;
	}
	s->nv = s->na = nverts;
	
	for (i = 0; i < nverts; i++) {
		gvert *v;

		/* Allocate and fill in each verticies basic information */
		if ((v = (gvert *)calloc(1, sizeof(gvert))) == NULL) {
			fprintf(stderr,"gamut: malloc failed on gvert object\n");
			return 2;
		}
		s->verts[i] = v;
		v->tag = 1;
		v->tn = v->n = i;
		v->f = GVERT_SET | GVERT_TRI;		/* Will be part of the triangulation */

		v->p[0] = *((double *)gam->t[0].fdata[i][Lf]);
		v->p[1] = *((double *)gam->t[0].fdata[i][af]);
		v->p[2] = *((double *)gam->t[0].fdata[i][bf]);

		gamut_rect2radial(s, v->r, v->p);
	}
	s->ntv = i;

	/* Compute the other vertex values */
	compute_vertex_coords(s);

	/* Get ready to read the triangle data */
	if ((v0f = gam->find_field(gam, 1, "VERTEX_0")) < 0) {
		fprintf(stderr,"Input file doesn't contain field VERTEX_0");
		return 1;
	}
	if (gam->t[1].ftype[v0f] != i_t) {
		fprintf(stderr,"Field VERTEX_0 is wrong type");
		return 1;
	}
	if ((v1f = gam->find_field(gam, 1, "VERTEX_1")) < 0) {
		fprintf(stderr,"Input file doesn't contain field VERTEX_1");
		return 1;
	}
	if (gam->t[1].ftype[v1f] != i_t) {
		fprintf(stderr,"Field VERTEX_1 is wrong type");
		return 1;
	}
	if ((v2f = gam->find_field(gam, 1, "VERTEX_2")) < 0) {
		fprintf(stderr,"Input file doesn't contain field VERTEX_2");
		return 1;
	}
	if (gam->t[1].ftype[v2f] != i_t) {
		fprintf(stderr,"Field VERTEX_2 is wrong type");
		return 1;
	}

	/* Create all the triangles */
	for (i = 0; i < ntris; i++) {
		gtri *t;
		int v0, v1, v2;

		t = new_gtri();
		ADD_ITEM_TO_BOT(s->tris, t);	/* Append to triangulation list */

		v0 = *((int *)gam->t[1].fdata[i][v0f]);
		v1 = *((int *)gam->t[1].fdata[i][v1f]);
		v2 = *((int *)gam->t[1].fdata[i][v2f]);

		t->v[0] = s->verts[v0];
		t->v[1] = s->verts[v1];
		t->v[2] = s->verts[v2];

		comptriattr(s, t);		/* Compute triangle attributes */
	}

	/* Connect edge information */
	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		int en;

		for (en = 0; en < 3; en++) {	/* For each edge */
			gedge *e;
			gvert *v0, *v1;				/* The two verticies of the edge */
			gtri *tp2;					/* The other triangle */
			int em;						/* The other edge */
			gvert *w0, *w1;				/* The other verticies */
			
			v0 = tp->v[en];
			v1 = tp->v[en < 2 ? en+1 : 0];
		
			if (v0->n > v1->n)
				continue;				/* Skip every other edge */

			/* Find the corresponding edge of the other triangle */
			w0 = w1 = NULL;
			tp2 = s->tris; 
			FOR_ALL_ITEMS(gtri, tp2) {
				for (em = 0; em < 3; em++) {	/* For each edge */
					w0 = tp2->v[em];
					w1 = tp2->v[em < 2 ? em+1 : 0];
					if (v0 == w1 && v1 == w0)	/* Found it */
						break;
				}
				if (em < 3)
					break;
			} END_FOR_ALL_ITEMS(tp2);
			if (w0 == NULL) {
				/* Should clean up ? */
				fprintf(stderr,".gam file triangle data is not consistent\n");
				return 1;
			}

			if (tp->e[en] != NULL
			 || tp2->e[em] != NULL) {
				fprintf(stderr,".gam file triangle data is not consistent\n");
				fprintf(stderr,"tp1->e[%d] = 0x%lx, tp2->e[%d]= 0x%lx\n",en,
				        (unsigned long)tp->e[en],em,(unsigned long)tp2->e[em]);
				return 1;
			}

			/* Creat the edge structure */
			e = new_gedge();
			ADD_ITEM_TO_BOT(s->edges, e);	/* Append to edge list */
			tp->e[en] = e;			/* This edge */
			tp->ei[en] = 0;			/* 0th triangle in edge */
			e->t[0] = tp;			/* 0th triangle is tp */
			e->ti[0] = en;			/* 0th triangles en edge */
			tp2->e[em] = e;			/* This edge */
			tp2->ei[em] = 1;		/* 1st triangle in edge */
			e->t[1] = tp2;			/* 1st triangle is tp2 */
			e->ti[1] = em;			/* 1st triangles em edge */
			e->v[0] = v0;			/* The two verticies */
			e->v[1] = v1;
		}
	} END_FOR_ALL_ITEMS(tp);

	gam->del(gam);			/* Clean up */

	s->read_inited = 1;			/* It's now valid */

#ifdef ASSERTS
	check_triangulation(s, 1);	/* Check out our work */
#endif

	return 0;
}

/* ===================================================== */
/* ===================================================== */

/* Convert from rectangular to radial coordinates */
void
gamut_rect2radial(
gamut *s,
double out[3],			/* Radius, longitude, lattitude out */
double in[3]			/* Lab in */
) {
	double L, a, b;	/* Lab values */
	double R, g, t;	/* Radial value */
	double c;	/* Chromatic length */
	

	L = in[0] - s->cent[0];	/* Offset value */
	a = in[1] - s->cent[1];
	b = in[2] - s->cent[2];
	c = a * a + b * b;
	R = c + L * L;
	c = sqrt(c);		/* Saturation */
	R = sqrt(R);		/* Vector length */

	if (R < 1e-6) {	/* Hmm, a point at the center */
		g = t = 0.0;
	} else {
	
		/* Figure out the longitude, -pi to +pi */
		if (c < 1e-6) {
			g = 0.0;
		} else {
			g = asin(b/c);
			if (a < 0.0) {
				if (b >= 0.0)
					g = M_PI - g;
				else
					g = -g - M_PI;
			}
		}
	
		/* Figure out the lattitude, -pi/2 to +pi/2 */
		t = asin(L/R);
	}
	out[0] = R;
	out[1] = g;
	out[2] = t;
}

/* Convert from radial to rectangular coordinates */
void
gamut_radial2rect(
gamut *s,
double out[3],			/* Lab out */
double in[3]			/* Radius, longitude, lattitude in */
) {
	double R, g, t;	/* Radial value */
	double L, a, b;	/* Lab values */
	double c;		/* Chromatic length */

	R = in[0];
	g = in[1];
	t = in[2];

	L = R * sin(t);
	c = R * cos(t);

	a = c * cos(g);
	b = c * sin(g);

	out[0] = L + s->cent[0];
	out[1] = a + s->cent[1];
	out[2] = b + s->cent[2];
}


/* -------------------------------------------------- */

/* Convert a gamut Lab value to an RGB value for display purposes */
void
gamut_Lab2RGB(double *out, double *in) {
	double L = in[0], a = in[1], b = in[2];
	double x,y,z,fx,fy,fz;
	double R, G, B;

	/* Scale so that black is visible */
	L = L * (100 - 40.0)/100.0 + 40.0;

	/* First convert to XYZ using D50 white point */
	if (L > 8.0) {
		fy = (L + 16.0)/116.0;
		y = pow(fy,3.0);
	} else {
		y = L/903.2963058;
		fy = 7.787036979 * y + 16.0/116.0;
	}

	fx = a/500.0 + fy;
	if (fx > 24.0/116.0)
		x = pow(fx,3.0);
	else
		x = (fx - 16.0/116.0)/7.787036979;

	fz = fy - b/200.0;
	if (fz > 24.0/116.0)
		z = pow(fz,3.0);
	else
		z = (fz - 16.0/116.0)/7.787036979;

	x *= 0.9642;	/* Multiply by white point, D50 */
	y *= 1.0;
	z *= 0.8249;

	/* Now convert to sRGB values */
	R = x * 3.2410  + y * -1.5374 + z * -0.4986;
	G = x * -0.9692 + y * 1.8760  + z * 0.0416;
	B = x * 0.0556  + y * -0.2040 + z * 1.0570;

	if (R < 0.0)
		R = 0.0;
	else if (R > 1.0)
		R = 1.0;

	if (G < 0.0)
		G = 0.0;
	else if (G > 1.0)
		G = 1.0;

	if (B < 0.0)
		B = 0.0;
	else if (B > 1.0)
		B = 1.0;

	R = pow(R, 1.0/2.2);
	G = pow(G, 1.0/2.2);
	B = pow(B, 1.0/2.2);

	out[0] = R;
	out[1] = G;
	out[2] = B;
}


/* -------------------------------------------------- */

#ifdef DEBUG_TRIANG

/* Write a surface contrsuction diagnostic VRML .wrl file */
static int write_diag_vrml(
gamut *s,
double vv[3],		/* Vertex being added */
int nh,				/* Number of hit triangles */
tidxs *hixs,		/* verticy indexes of hit triangles */
gtri *hl			/* Edge hit list (may be NULL) */
) {
	char *filename;
	int i, j;
	gtri *tp;		/* Triangle pointer */
	FILE *wrl;

	if (hl)
		filename = "diag1.wrl";		/* Triangles hit */
	else
		filename = "diag2.wrl";		/* Triangles formed */

	if ((wrl = fopen(filename,"w")) == NULL) {
		fprintf(stderr,"Error opening output file '%s'\n",filename);
		return 2;
	}

	/* Spit out a VRML 2 Object surface of gamut */

	fprintf(wrl,"#VRML V2.0 utf8\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"# Created by the Argyll CMS\n");
	fprintf(wrl,"Transform {\n");
  	fprintf(wrl,"children [\n");
    fprintf(wrl,"	NavigationInfo {\n");
	fprintf(wrl,"		type \"EXAMINE\"        # It's an object we examine\n");
	fprintf(wrl,"	} # We'll add our own light\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    DirectionalLight {\n");
	fprintf(wrl,"        direction 0 0 -1      # Light illuminating the scene\n");
	fprintf(wrl,"        direction 0 -1 0      # Light illuminating the scene\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    Viewpoint {\n");
	fprintf(wrl,"        position 0 0 200      # Position we view from\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");

	fprintf(wrl,"    Transform {\n");
	fprintf(wrl,"      translation 0 0 0\n");
	fprintf(wrl,"      children [\n");
	fprintf(wrl,"		Shape { \n");
	fprintf(wrl,"		    geometry IndexedFaceSet {\n");
	fprintf(wrl,"				ccw FALSE\n");
	fprintf(wrl,"				convex TRUE\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"		        coord Coordinate { \n");
	fprintf(wrl,"		            point [			# Verticy coordinates\n");

	/* Spit out the vertex values, in order. */
	/* Note that a->x, b->y, L->z */
	for (i = 0; i < s->nv; i++) {
		double out[3];

		/* Show normal gamut surface */
		out[0] = s->verts[i]->ch[0];
		out[1] = s->verts[i]->ch[1];
		out[2] = s->verts[i]->ch[2];

		fprintf(wrl,"%f %f %f,\n",out[1], out[2], out[0]);
	}
	fprintf(wrl,"					]\n");
	fprintf(wrl,"		        }\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"		        coordIndex [ 		# Indexes of poligon Verticies \n");

	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		fprintf(wrl,"%d, %d, %d, -1\n", tp->v[0]->n, tp->v[1]->n, tp->v[2]->n);
	} END_FOR_ALL_ITEMS(tp);

	for (i = 0; i < nh; i++) { 
		fprintf(wrl,"%d, %d, %d, -1\n", hixs[i].tix[0], hixs[i].tix[1], hixs[i].tix[2]);
	}

	fprintf(wrl,"				]\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"				colorPerVertex FALSE\n");
	fprintf(wrl,"		        color Color {\n");
	fprintf(wrl,"		            color [			# RGB colors of each vertex\n");

	/* Spit out the colors for each face */
	tp = s->tris; 
	FOR_ALL_ITEMS(gtri, tp) {
		fprintf(wrl,"%f %f %f,\n", 0.7, 0.7, 0.7);
	} END_FOR_ALL_ITEMS(tp);
	for (i = 0; i < nh; i++) { 
		if (hixs[i].type == 0)
			fprintf(wrl,"%f %f %f,\n", 0.4, 1.0, 0.4);	/* Green for hit */
		else if (hixs[i].type == 1)
			fprintf(wrl,"%f %f %f,\n", 0.4, 0.4, 1.0);	/* Blue for extra */
		else
			fprintf(wrl,"%f %f %f,\n", 0.8, 0.8, 0.2);	/* Yellow for new */
	}
	fprintf(wrl,"					] \n");
	fprintf(wrl,"		        }\n");
	fprintf(wrl,"		    }	# end IndexedFaceSet\n");

	fprintf(wrl,"		    appearance Appearance { \n");
	fprintf(wrl,"		        material Material {\n");
	fprintf(wrl,"					transparency 0.0\n");
	fprintf(wrl,"					ambientIntensity 0.3\n");
	fprintf(wrl,"					shininess 0.5\n");
	fprintf(wrl,"				}\n");
	fprintf(wrl,"		    }\n");
	fprintf(wrl,"		}	# end Shape\n");
	fprintf(wrl,"      ]\n");
	fprintf(wrl,"    } # end of transform\n");

	/* center of gamut */
	fprintf(wrl,"\n");
	fprintf(wrl,"    Transform {\n");
	fprintf(wrl,"      translation %f %f %f\n",0.0, 0.0, 0.0);
	fprintf(wrl,"      children [\n");
	fprintf(wrl,"		Shape { \n");
	fprintf(wrl,"		 geometry Sphere { radius 1.5 }\n");
	fprintf(wrl,"        appearance Appearance { material Material { diffuseColor %f %f %f } }\n", 1.0, 1.0, 0.0);
	fprintf(wrl,"		} \n");
	fprintf(wrl,"      ]\n");
	fprintf(wrl,"    }\n");

	/* vertex being added */
	fprintf(wrl,"\n");
	fprintf(wrl,"    Transform {\n");
	fprintf(wrl,"      translation %f %f %f\n",vv[1], vv[2], vv[0]);
	fprintf(wrl,"      children [\n");
	fprintf(wrl,"		Shape { \n");
	fprintf(wrl,"		 geometry Sphere { radius 1.5 }\n");
	fprintf(wrl,"        appearance Appearance { material Material { diffuseColor %f %f %f } }\n", 1.0, 0.0, 0.0);
	fprintf(wrl,"		} \n");
	fprintf(wrl,"      ]\n");
	fprintf(wrl,"    }\n");


	/* Verticies for Polygon edges */
	if (hl != NULL) {
		double base[3] = { 0.0, 0.0, 1.0 };		/* Default orientation of cone is b axis */
		tp = hl; 
		FOR_ALL_ITEMS(gtri, tp) {
			double len;
			double loc[3];
			double vec[3];
			double axis[3];		/* Axis to rotate around */
			double rot;			/* In radians */

//printf("~1 edge vert %d to %d\n",tp->v[0]->n, tp->v[1]->n);
//printf("~1 edge %f %f %f to %f %f %f\n",
//tp->v[0]->ch[0], tp->v[0]->ch[1], tp->v[0]->ch[2],
//tp->v[1]->ch[0], tp->v[1]->ch[1], tp->v[1]->ch[2]);

			icmAdd3(loc, tp->v[1]->ch, tp->v[0]->ch);
			icmScale3(loc, loc, 0.5);
			icmSub3(vec, tp->v[1]->ch, tp->v[0]->ch);
			len = icmNorm3(vec);

			if (len < 1.0)
				len = 1.0;

			icmNormalize3(base, base, 1.0);
			icmNormalize3(vec, vec, 1.0);
			icmCross3(axis, base, vec);
			rot = icmDot3(base, vec);
//printf("~1 Axis = %f %f %f\n",axis[0],axis[1],axis[2]);
			if (icmNorm3sq(axis) < 1e-10) {		/* 0 or 180 degrees */
				double base2[3];
				int mxi = 0;
				base2[0] = vec[1];		/* Comute vector in a different direction */
				base2[1] = vec[2];
				base2[2] = vec[0];
				for (j = 1; j < 3; j++) {
					if (fabs(base2[j]) > fabs(base2[mxi])) 
						mxi = j;
				}
				base2[mxi] = -base2[mxi];
					
				icmCross3(axis, base2, vec);
				if (icmNorm3sq(axis) < 1e-10) {		/* 0 or 180 degrees */
					error("VRML rotate axis still too small");
				}
				if (rot < 0.0)
					rot = 3.1415926;
				else
					rot = 0.0;			
			} else {
				rot = acos(rot);
//printf("~1 rotation %f\n",rot);
			}

			fprintf(wrl,"\n");
			fprintf(wrl,"    Transform {\n");
			fprintf(wrl,"      rotation %f %f %f %f\n",axis[1], axis[2], axis[0], rot);
			fprintf(wrl,"      translation %f %f %f\n",loc[1], loc[2], loc[0]);
			fprintf(wrl,"      children [\n");
			fprintf(wrl,"		Shape { \n");
			fprintf(wrl,"		 geometry Cone { bottomRadius 0.5 height %f }\n",len);
			fprintf(wrl,"        appearance Appearance { material Material { diffuseColor 0.7 0.0 1.0 } }\n");
			fprintf(wrl,"		} \n");
			fprintf(wrl,"      ]\n");
			fprintf(wrl,"    }\n");
		} END_FOR_ALL_ITEMS(tp);
	}

	fprintf(wrl,"\n");
	fprintf(wrl,"  ] # end of children for world\n");
	fprintf(wrl,"}\n");

	if (fclose(wrl) != 0) {
		fprintf(stderr,"Error closing output file '%s'\n",filename);
		return 2;
	}

	return 0;
}

#endif /* DEBUG_TRIANG */


#ifdef DEBUG_SPLIT_VRML

/* Write a triangle split diagnostic VRML .wrl file */
static int write_split_diag_vrml(
gamut *s,
gtri **list,	/* Triangle list */
int llen		/* Number of triangles in the list */
) {
	char *filename;
	int i, j;
	FILE *wrl;

	filename = "diag3.wrl";		/* Triangles split */

	if ((wrl = fopen(filename,"w")) == NULL) {
		fprintf(stderr,"Error opening output file '%s'\n",filename);
		return 2;
	}

	/* Spit out a VRML 2 Object surface of gamut */

	fprintf(wrl,"#VRML V2.0 utf8\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"# Created by the Argyll CMS\n");
	fprintf(wrl,"Transform {\n");
  	fprintf(wrl,"children [\n");
    fprintf(wrl,"	NavigationInfo {\n");
	fprintf(wrl,"		type \"EXAMINE\"        # It's an object we examine\n");
	fprintf(wrl,"	} # We'll add our own light\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    DirectionalLight {\n");
	fprintf(wrl,"        ambientIntensity 0.3           # Ambient light illuminating the scene\n");
	fprintf(wrl,"        direction 0 0 -1      # Light illuminating the scene\n");
	fprintf(wrl,"        direction 0 -1 0      # Light illuminating the scene\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    Viewpoint {\n");
	fprintf(wrl,"        position 0 0 5      # Position we view from\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");

	fprintf(wrl,"    Transform {\n");
	fprintf(wrl,"      translation 0 0 0\n");
	fprintf(wrl,"      children [\n");
	fprintf(wrl,"		Shape { \n");
	fprintf(wrl,"		    geometry IndexedFaceSet {\n");
	fprintf(wrl,"				ccw FALSE\n");
	fprintf(wrl,"				convex TRUE\n");
	fprintf(wrl,"				solid FALSE\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"		        coord Coordinate { \n");
	fprintf(wrl,"		            point [			# Verticy coordinates\n");

	/* Spit out the vertex values, in order. */
	for (i = 0; i < llen; i++) {
		fprintf(wrl,"%f %f %f,\n",list[i]->v[0]->sp[0], list[i]->v[0]->sp[1], list[i]->v[0]->sp[2]);
		fprintf(wrl,"%f %f %f,\n",list[i]->v[1]->sp[0], list[i]->v[1]->sp[1], list[i]->v[1]->sp[2]);
		fprintf(wrl,"%f %f %f,\n",list[i]->v[2]->sp[0], list[i]->v[2]->sp[1], list[i]->v[2]->sp[2]);
	}
	fprintf(wrl,"					]\n");
	fprintf(wrl,"		        }\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"		        coordIndex [ 		# Indexes of poligon Verticies \n");

	for (i = 0; i < llen; i++) {
		fprintf(wrl,"%d, %d, %d, -1\n", i * 3 + 0, i * 3 + 1, i * 3 + 2);
	}
	fprintf(wrl,"				]\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"				colorPerVertex FALSE\n");
	fprintf(wrl,"		        color Color {\n");
	fprintf(wrl,"		            color [			# RGB colors of each vertex\n");

	/* Spit out the colors for each face */
	for (i = 0; i < llen; i++) {
		if (list[i]->bsort == 1) {	/* Positive */
			fprintf(wrl,"%f %f %f,\n", 1.0, 0.3, 0.3);		/* Red */
		} else if (list[i]->bsort == 2) {	/* Negative */
			fprintf(wrl,"%f %f %f,\n", 0.3, 1.0, 0.3);		/* Green */
		} else if (list[i]->bsort == 3) {	/* Both */
			fprintf(wrl,"%f %f %f,\n", 1.0, 1.0, 0.3);		/* Yellow */
		} else {	/* Neither */
			fprintf(wrl,"%f %f %f,\n", 0.3, 0.3, 1.0);		/* Blue */
		} 
	}
	fprintf(wrl,"					] \n");
	fprintf(wrl,"		        }\n");
	fprintf(wrl,"		    }	# end IndexedFaceSet\n");

	fprintf(wrl,"		    appearance Appearance { \n");
	fprintf(wrl,"		        material Material {\n");
	fprintf(wrl,"					transparency 0.0\n");
	fprintf(wrl,"					ambientIntensity 0.3\n");
	fprintf(wrl,"					shininess 0.5\n");
	fprintf(wrl,"				}\n");
	fprintf(wrl,"		    }\n");
	fprintf(wrl,"		}	# end Shape\n");
	fprintf(wrl,"      ]\n");
	fprintf(wrl,"    } # end of transform\n");

	/* center of gamut */
	fprintf(wrl,"\n");
	fprintf(wrl,"    Transform {\n");
	fprintf(wrl,"      translation %f %f %f\n",0.0, 0.0, 0.0);
	fprintf(wrl,"      children [\n");
	fprintf(wrl,"		Shape { \n");
	fprintf(wrl,"		 geometry Sphere { radius 0.05 }\n");
	fprintf(wrl,"        appearance Appearance { material Material { diffuseColor %f %f %f } }\n", 1.0, 1.0, 0.0);
	fprintf(wrl,"		} \n");
	fprintf(wrl,"      ]\n");
	fprintf(wrl,"    }\n");

	fprintf(wrl,"\n");
	fprintf(wrl,"  ] # end of children for world\n");
	fprintf(wrl,"}\n");

	if (fclose(wrl) != 0) {
		fprintf(stderr,"Error closing output file '%s'\n",filename);
		return 2;
	}

	return 0;
}

#endif /* DEBUG_SPLIT_VRML */




















































