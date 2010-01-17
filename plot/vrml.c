
/*
 * Simple diagnostic VRML function library for debugging
 *
 * Copyright 2005 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "numlib.h"
#include "icc.h"
#include "gamut.h"
#include "vrml.h"

#define MAKE_SOLID

/* Add a shere at the given location. */
/* if col[] is NULL, use natural color. */
/* Need to do this before or after start_line_set()/dd_vertex()/make_lines() ! */
static void add_marker(vrml *s, double pos[3], double col[3], double rad) {
	double rgb[3];

	if (rad <= 0.0)
		rad = 1.0;

	if (col == NULL)
		s->Lab2RGB(s, rgb, pos);
	else {
		rgb[0] = col[0];
		rgb[1] = col[1];
		rgb[2] = col[2];
	}

	fprintf(s->fp,"    # Shere\n");
	fprintf(s->fp,"    Transform { translation %f %f %f\n", pos[1], pos[2], pos[0]-s->lcent);
	fprintf(s->fp,"      children [\n");
	fprintf(s->fp,"        Shape{\n");
	fprintf(s->fp,"          geometry Sphere { radius %f}\n", rad);
	fprintf(s->fp,"          appearance Appearance { material Material ");
	fprintf(s->fp,"{ diffuseColor %f %f %f} }\n", rgb[0], rgb[1], rgb[2]);
	fprintf(s->fp,"        }\n");
	fprintf(s->fp,"      ]\n");
	fprintf(s->fp,"    }\n");
}

/* Add a cone marker to the plot. col == NULL for natural color  */
/* Need to do this before or after start_line_set()/dd_vertex()/make_lines() ! */
static void add_cone(vrml *s, double pp0[3], double pp1[3], double col[3], double rad) {
	double rgb[3];
	double p0[3], p1[3];

	icmAry2Ary(p0, pp0);
	icmAry2Ary(p1, pp1);

//printf("~1 cone %f %f %f -> %f %f %f rad %f\n", p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], rad);

	if (rad <= 0.0)
		rad = 1.0;

	if (col == NULL) {
		icmAdd3(rgb, p1, p0);
		icmScale3(rgb, rgb, 0.5);		/* Compute half way value */
		s->Lab2RGB(s, rgb, rgb);
	} else {
		rgb[0] = col[0];
		rgb[1] = col[1];
		rgb[2] = col[2];
	}

	p0[0] -= s->lcent;
	p1[0] -= s->lcent;

	{
		double base[3] = { 0.0, 0.0, 1.0 };		/* Default orientation of cone is b axis */
		double len;
		double loc[3];
		double vec[3];
		double axis[3];		/* Axis to rotate around */
		double rot;			/* In radians */
		int j;

//printf("~1 edge vert %d to %d\n",tp->v[0]->n, tp->v[1]->n);
//printf("~1 edge %f %f %f to %f %f %f\n",
//tp->v[0]->ch[0], tp->v[0]->ch[1], tp->v[0]->ch[2],
//tp->v[1]->ch[0], tp->v[1]->ch[1], tp->v[1]->ch[2]);

		icmAdd3(loc, p1, p0);
		icmScale3(loc, loc, 0.5);		/* Compute half way value */
		icmSub3(vec, p1, p0);
		len = icmNorm3(vec);
//printf("~1 loc = %f %f %f\n", loc[0], loc[1], loc[2]);
//printf("~1 vec = %f %f %f\n", vec[0], vec[1], vec[2]);
//printf("~1 len = %f\n", len);

		if (len < 0.1)
			len = 0.1;

		icmNormalize3(base, base, 1.0);
		icmNormalize3(vec, vec, 1.0);
		icmCross3(axis, base, vec);
		rot = icmDot3(base, vec);
//printf("~1 base = %f %f %f\n", base[0], base[1], base[2]);
//printf("~1 vec = %f %f %f\n", vec[0], vec[1], vec[2]);
//printf("~1 axis = %f %f %f, rot = %f\n",axis[0],axis[1],axis[2],rot);
		if (icmNorm3sq(axis) < 1e-10) {		/* 0 or 180 degrees */
			double base2[3];
			int mxi = 0;
//printf("~1 computing a different axis\n");
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

		fprintf(s->fp,"\n");
		fprintf(s->fp,"    # Cone\n");
		fprintf(s->fp,"    Transform {\n");
		fprintf(s->fp,"      rotation %f %f %f %f\n",axis[1], axis[2], axis[0], rot);
		fprintf(s->fp,"      translation %f %f %f\n",loc[1], loc[2], loc[0]);
		fprintf(s->fp,"      children [\n");
		fprintf(s->fp,"		Shape { \n");
		fprintf(s->fp,"		 geometry Cone { bottomRadius %f height %f }\n",rad,len);
		fprintf(s->fp,"        appearance Appearance { material Material { diffuseColor %f %f %f } }\n",rgb[0],rgb[1],rgb[2]);
		fprintf(s->fp,"		} \n");
		fprintf(s->fp,"      ]\n");
		fprintf(s->fp,"    }\n");
	}
}

/* Add a text marker to the plot. col == NULL for natural color  */
/* (Need to do this before or after start_line_set()/dd_vertex()/make_lines() !) */
static void add_text(vrml *s, char *text, double p[3], double col[3], double size) {
	double rgb[3];

	if (size <= 0.0)
		size = 1.0;

	if (col == NULL) {
		s->Lab2RGB(s, rgb, p);
	} else {
		rgb[0] = col[0];
		rgb[1] = col[1];
		rgb[2] = col[2];
	}
	fprintf(s->fp,"    # Text\n");
	fprintf(s->fp,"    Transform { translation %f %f %f\n", p[1], p[2], p[0]-s->lcent);
	fprintf(s->fp,"      children [\n");
	fprintf(s->fp,"        Shape{\n");
	fprintf(s->fp,"          geometry Text { string [\"%s\"]\n",text);
	fprintf(s->fp,"            fontStyle FontStyle { family \"SANS\" style \"BOLD\" size %f }\n",
	                                                                                      size);
	fprintf(s->fp,"                        }\n");
	fprintf(s->fp,"          appearance Appearance { material Material ");
	fprintf(s->fp,"{ diffuseColor %f %f %f} }\n", rgb[0], rgb[1], rgb[2]);
	fprintf(s->fp,"        }\n");
	fprintf(s->fp,"      ]\n");
	fprintf(s->fp,"    }\n");
}

/* Start building up verticies that will be converted to lines */
/* Set can be from 0 - 9 */
static void start_line_set(vrml *s, int set) {
	if (set < 0 || set > 9)
		error("vrml start_line_set set %d out of range",set);
	s->set[set].npoints = 0;
}

/* Add a verticy with color */
static void add_col_vertex_l(vrml *s, int set, double pos[3], double col[3], int last) {

	if (set < 0 || set > 9)
		error("vrml start_line_set set %d out of range",set);

	if (s->set[set].npoints >= s->set[set].paloc) {
		s->set[set].paloc = (s->set[set].paloc + 10) * 2;
		if (s->set[set].pary == NULL)
			s->set[set].pary = malloc(s->set[set].paloc * 6 * (sizeof(double) + sizeof(int)));
		else
			s->set[set].pary = realloc(s->set[set].pary, s->set[set].paloc * 6 * (sizeof(double) + sizeof(int)));

		if (s->set[set].pary == NULL)
			error("VRML malloc failed at count %d\n",s->set[set].paloc);
	}
	s->set[set].pary[s->set[set].npoints].pp[0] = pos[0];
	s->set[set].pary[s->set[set].npoints].pp[1] = pos[1];
	s->set[set].pary[s->set[set].npoints].pp[2] = pos[2];
	s->set[set].pary[s->set[set].npoints].cc[0] = col[0];
	s->set[set].pary[s->set[set].npoints].cc[1] = col[1];
	s->set[set].pary[s->set[set].npoints].cc[2] = col[2];
	s->set[set].pary[s->set[set].npoints].last = last;
	s->set[set].npoints++;
}

/* Add a verticy with color */
static void add_col_vertex(vrml *s, int set, double pos[3], double col[3]) {

	add_col_vertex_l(s, set, pos, col, 0);
}

/* Add a color verticy */
static void add_vertex(vrml *s, int set, double pos[3]) {
	double col[3] = { -1.0, -1.0, -1.0 };

	add_col_vertex_l(s, set, pos, col, 0);
}

/* Turn the last added vertex into the last vertex of the line */
static void make_last_vertex(vrml *s, int set) {

	if (set < 0 || set > 9)
		error("vrml start_line_set set %d out of range",set);

	if (s->set[set].npoints <= 0)
		warning("vrml plot: tried to set last point with no points added!\n");
	else
		s->set[set].pary[s->set[set].npoints-1].last = 1;
}

/* Convert the verticies to lines, ppset verticies per line (or .last flag) */
static void make_lines(vrml *s, int set, int ppset) {
	int i, j;

	if (set < 0 || set > 9)
		error("vrml start_line_set set %d out of range",set);

	fprintf(s->fp,"\n");
	fprintf(s->fp,"    # Lines\n");
	fprintf(s->fp,"    Shape {\n");
	fprintf(s->fp,"      geometry IndexedLineSet { \n");
	fprintf(s->fp,"        coord Coordinate { \n");
	fprintf(s->fp,"          point [\n");

	for (i = 0; i < s->set[set].npoints; i++) {
		fprintf(s->fp,"            %f %f %f,\n",s->set[set].pary[i].pp[1], s->set[set].pary[i].pp[2],
		                            s->set[set].pary[i].pp[0] - s->lcent);
	}
	
	fprintf(s->fp,"          ]\n");
	fprintf(s->fp,"        }\n");
	fprintf(s->fp,"        coordIndex [\n");

	for (i = 0; i < s->set[set].npoints;) {
		fprintf(s->fp,"          ");
		for (j = 0; i < s->set[set].npoints  && j < ppset; j++) {
			fprintf(s->fp,"%d, ", i++);
			if (s->set[set].pary[i-1].last != 0)
				break;
		}
		fprintf(s->fp,"-1,\n");
	}
	fprintf(s->fp,"        ]\n");

	/* Color */
	fprintf(s->fp,"        colorPerVertex TRUE\n");
	fprintf(s->fp,"        color Color {\n");
	fprintf(s->fp,"          color [			# RGB colors of each vertex\n");

	for (i = 0; i < s->set[set].npoints; i++) {
		double rgb[3], Lab[3];

		if (s->set[set].pary[i].cc[0] < 0.0) {
			Lab[0] = s->set[set].pary[i].pp[0];
			Lab[1] = s->set[set].pary[i].pp[1];
			Lab[2] = s->set[set].pary[i].pp[2];
			s->Lab2RGB(s, rgb, Lab);
			fprintf(s->fp,"            %f %f %f,\n", rgb[0], rgb[1], rgb[2]);
		} else {
			fprintf(s->fp,"            %f %f %f,\n", s->set[set].pary[i].cc[0], s->set[set].pary[i].cc[1], s->set[set].pary[i].cc[2]);
		}
	}
	fprintf(s->fp,"          ] \n");
	fprintf(s->fp,"        }\n");
	/* End color */

	fprintf(s->fp,"      }\n");
	fprintf(s->fp,"    } # end shape\n");
}

/* Convert the verticies to triangles */
static void make_triangles_imp(
vrml *s,
int set,
double trans,	/* Transparency level */
int ixcol,		/* NZ for using index color */
double cc[3]	/* Surface color, cc == NULL or cc[0] < 0.0 for natural color */
) {
	int i, nverts, ix;
	int v[3];

	if (set < 0 || set > 9)
		error("vrml start_line_set set %d out of range",set);

	fprintf(s->fp,"    # Triangles\n");
	fprintf(s->fp,"    Transform {\n");
	fprintf(s->fp,"      translation 0 0 0\n");
	fprintf(s->fp,"      children [\n");
	fprintf(s->fp,"        Shape { \n");
	fprintf(s->fp,"          geometry IndexedFaceSet {\n");
//	fprintf(s->fp,"            ccw FALSE\n");
	fprintf(s->fp,"            convex TRUE\n");
#ifdef MAKE_SOLID
	fprintf(s->fp,"            solid FALSE\n");	/* If we want them visible from both sides */
#endif
	fprintf(s->fp,"\n");
	fprintf(s->fp,"            coord Coordinate { \n");
	fprintf(s->fp,"              point [			# Verticy coordinates\n");

	/* Spit out the point values, in order. */
	/* Note that a->x, b->y, L->z */
	for (i = 0; i < s->set[set].npoints; i++) {
		fprintf(s->fp,"                %f %f %f,\n",s->set[set].pary[i].pp[1], s->set[set].pary[i].pp[2],
		                            s->set[set].pary[i].pp[0] - s->lcent);
	}
	fprintf(s->fp,"              ]\n");
	fprintf(s->fp,"            }\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"            coordIndex [ 		# Indexes of poligon Verticies \n");

	for (i = 0; i < s->ntris; i++) {
		if (s->tary[i].set == set)
			fprintf(s->fp,"              %d, %d, %d, -1\n", s->tary[i].ix[0], s->tary[i].ix[1], s->tary[i].ix[2]);
	}

	fprintf(s->fp,"            ]\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"            colorPerVertex TRUE\n");
	fprintf(s->fp,"            color Color {\n");
	fprintf(s->fp,"            color [			# RGB colors of each vertex\n");

	/* Spit out the colors for each vertex */
	for (i = 0; i < s->set[set].npoints; i++) {
		double out[3];
		double rgb[3];

		if (ixcol) {
			fprintf(s->fp,"              %f %f %f,\n",s->set[set].pary[i].cc[0], s->set[set].pary[i].cc[1], s->set[set].pary[i].cc[2]);
		} else {
			if (cc == NULL || cc[0] < 0.0) {
				s->Lab2RGB(s, rgb, s->set[set].pary[i].pp);
				fprintf(s->fp,"              %f %f %f,\n", rgb[0], rgb[1], rgb[2]);
			} else {
				fprintf(s->fp,"              %f %f %f,\n", cc[0], cc[1], cc[2]);
			}
		}
	}
	fprintf(s->fp,"              ] \n");
	fprintf(s->fp,"            }\n");
	fprintf(s->fp,"          }\n");
	fprintf(s->fp,"          appearance Appearance { \n");
	fprintf(s->fp,"            material Material {\n");
	fprintf(s->fp,"              transparency %f\n",trans);
	fprintf(s->fp,"              ambientIntensity 0.3\n");
	fprintf(s->fp,"              shininess 0.5\n");
	fprintf(s->fp,"            }\n");
	fprintf(s->fp,"          }\n");
	fprintf(s->fp,"        }	# end Shape\n");
	fprintf(s->fp,"      ]\n");
	fprintf(s->fp,"    }\n");
}

/* Convert the verticies to triangles with vertex color */
static void make_triangles_vc(
vrml *s,
int set,
double trans	/* Transparency level */
) {
	make_triangles_imp(s, set, trans, 1, NULL);
}

/* Convert the verticies to triangles with color */
static void make_triangles(
vrml *s,
int set,
double trans,	/* Transparency level */
double cc[3]	/* Surface color, cc == NULL or cc[0] < 0.0 for natural color */
) {
	make_triangles_imp(s, set, trans, 0, cc);
}

/* Add a triangle */
static void add_triangle(vrml *s, int set, int ix[3]) {

	if (set < 0 || set > 9)
		error("vrml start_line_set set %d out of range",set);

	if (s->ntris >= s->taloc) {
		s->taloc = (s->taloc + 10) * 2;
		if (s->tary == NULL)
			s->tary = malloc(s->taloc * 4 * sizeof(int));
		else
			s->tary = realloc(s->tary, s->taloc * 4 * sizeof(int));

		if (s->tary == NULL)
			error("VRML malloc failed at count %d\n",s->taloc);
	}
	s->tary[s->ntris].set = set;
	s->tary[s->ntris].ix[0] = ix[0];
	s->tary[s->ntris].ix[1] = ix[1];
	s->tary[s->ntris].ix[2] = ix[2];
	s->ntris++;
}

/* Create a gamut surface solid or wireframe from the given gamut. */
/* Use the given transparency level. */
/* Display in natural colors if c[0] < 0.0, */
/* or the given color otherwise */
static void make_gamut_surface_2(
vrml *s,
gamut *g,
double trans,	/* Transparency level */
int wire,		/* Z for solid, NZ for wireframe */
double cc[3]	/* Surface color, cc[0] < 0.0 for natural color */
) {
	int i, nverts, ix;
	int v[3];

	nverts = g->nverts(g);

	if (nverts == 0)
		return;

	fprintf(s->fp,"    # Gamut surface\n");
	fprintf(s->fp,"    Transform {\n");
	fprintf(s->fp,"      translation 0 0 0\n");
	fprintf(s->fp,"      children [\n");
	fprintf(s->fp,"        Shape { \n");
	if (wire) {
		fprintf(s->fp,"          geometry IndexedLineSet {\n");
	} else {
		fprintf(s->fp,"          geometry IndexedFaceSet {\n");
//		fprintf(s->fp,"            ccw FALSE\n");
		fprintf(s->fp,"            convex TRUE\n");
#ifdef MAKE_SOLID
		fprintf(s->fp,"            solid FALSE\n");	/* If we want them visible from both sides */
#endif
	}
	fprintf(s->fp,"\n");
	fprintf(s->fp,"            coord Coordinate { \n");
	fprintf(s->fp,"              point [			# Verticy coordinates\n");

	/* Spit out the point values, in order. */
	/* Note that a->x, b->y, L->z */
	for (ix = i = 0; ix >= 0 && i < nverts; i++) {
		double out[3];

		ix = g->getvert(g, NULL, out, ix);
		fprintf(s->fp,"              %f %f %f,\n",out[1], out[2], out[0]-s->lcent);
	}
	fprintf(s->fp,"              ]\n");
	fprintf(s->fp,"            }\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"          coordIndex [ 		# Indexes of poligon Verticies \n");

	g->startnexttri(g);
	while (g->getnexttri(g, v) == 0) {
		if (wire) {
			if (v[0] < v[1])				/* Only output 1 wire of two on an edge */
				fprintf(s->fp,"              %d, %d, -1\n", v[0], v[1]);
			if (v[1] < v[2])
				fprintf(s->fp,"              %d, %d, -1\n", v[1], v[2]);
			if (v[2] < v[0])
				fprintf(s->fp,"              %d, %d, -1\n", v[2], v[0]);
		} else {
			fprintf(s->fp,"            %d, %d, %d, -1\n", v[0], v[1], v[2]);
		}
	}
	fprintf(s->fp,"              ]\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"            colorPerVertex TRUE\n");
	fprintf(s->fp,"            color Color {\n");
	fprintf(s->fp,"              color [			# RGB colors of each vertex\n");

	/* Spit out the colors for each vertex */
	for (ix = i = 0; ix >= 0 && i < nverts; i++) {
		double out[3];
		double rgb[3];

		ix = g->getvert(g, NULL, out, ix);

		if (cc == NULL || cc[0] < 0.0) {
			s->Lab2RGB(s, rgb, out);
			fprintf(s->fp,"                %f %f %f,\n", rgb[0], rgb[1], rgb[2]);
		} else {
			fprintf(s->fp,"                %f %f %f,\n", cc[0], cc[1], cc[2]);
		}
	}
	fprintf(s->fp,"              ] \n");
	fprintf(s->fp,"            }\n");
	fprintf(s->fp,"          }\n");
	fprintf(s->fp,"          appearance Appearance { \n");
	fprintf(s->fp,"            material Material {\n");
	fprintf(s->fp,"              transparency %f\n",trans);
	fprintf(s->fp,"              ambientIntensity 0.3\n");
	fprintf(s->fp,"              shininess 0.5\n");
	fprintf(s->fp,"            }\n");
	fprintf(s->fp,"          }\n");
	fprintf(s->fp,"        }	# end Shape\n");
	fprintf(s->fp,"      ]\n");
	fprintf(s->fp,"    }\n");
}

/* Create a gamut surface from the given gamut. */
/* Use the given transparency level. */
/* Display in natural colors if c[0] < 0.0, */
/* or the given color otherwise */
static void make_gamut_surface(
vrml *s,
gamut *g,
double trans,	/* Transparency level */
double cc[3]	/* Surface color, cc[0] < 0.0 for natural color */
) {
	s->make_gamut_surface_2(s, g, trans, 0, cc);
}

/* Add cusp markers from a gamut surface */
/* Use the given transparency level. */
/* Display in natural colors if c[0] < 0.0, */
/* or the given color otherwise */
static void add_cusps(
vrml *s,
gamut *g,
double trans,	/* Transparency level */
double cc[3]	/* Surface color, cc[0] < 0.0 for natural color, NULL for default */
) {
	double cusps[6][3];
	double ccolors[6][3] = {
		{ 1.0, 0.1, 0.1 },		/* Red */
		{ 1.0, 1.0, 0.1 },		/* Yellow */
		{ 0.1, 1.0, 0.1 },		/* Green */
		{ 0.1, 1.0, 1.0 },		/* Cyan */
		{ 0.1, 0.1, 1.0 },		/* Blue */
		{ 1.0, 0.1, 1.0 }		/* Magenta */
	};
	double rgb[3];
	double *cv = NULL;
	int i;
	int v[3];

	if (g->getcusps(g, cusps) != 0)
		return;

	fprintf(s->fp,"    # Cusps\n");
	for (i = 0; i < 6; i++) {
		if (cc == NULL) {
			cv = ccolors[i];
		} else if (cc[0] < 0.0) {
			s->Lab2RGB(s, rgb, cusps[i]);
			cv = rgb;
		} else {
			cv = cc;
		}
		fprintf(s->fp,"\n");
		fprintf(s->fp,"    Transform {\n");
		fprintf(s->fp,"      translation %f %f %f\n",cusps[i][1], cusps[i][2], cusps[i][0]-s->lcent);
		fprintf(s->fp,"      children [\n");
		fprintf(s->fp,"		   Shape { \n");
		fprintf(s->fp,"		    geometry Sphere { radius 2.0 }\n");
		fprintf(s->fp,"         appearance Appearance { \n");
		fprintf(s->fp,"           material Material {\n");
		fprintf(s->fp,"             transparency %f\n",trans);
		fprintf(s->fp,"             ambientIntensity 0.3\n");
		fprintf(s->fp,"             shininess 0.5\n");
		fprintf(s->fp,"             diffuseColor %f %f %f\n", cv[0],cv[1],cv[2]);
		fprintf(s->fp,"           }\n");
		fprintf(s->fp,"         }\n");
		fprintf(s->fp,"		  } \n");
		fprintf(s->fp,"      ]\n");
		fprintf(s->fp,"    }\n");
	}
}

/* Clear verticies and triangles */
static void clear(vrml *s) {
	int i;

	for (i = 0; i < 10; i++) {
		if (s->set[i].pary != NULL)
			free(s->set[i].pary);
		s->set[i].pary = NULL;
		s->set[i].npoints = s->set[i].paloc = 0;
	}
	if (s->tary != NULL)
		free(s->tary);
	s->tary = NULL;
	s->ntris = s->taloc = 0;
}
	
/* Helper :- convert a Lab value to RGB for display purposes */
static void Lab2RGB(vrml *s, double *out, double *in) {
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

	/* For a black background: */
//	R = R * 0.85 + 0.15;
//	G = G * 0.85 + 0.15;
//	B = B * 0.85 + 0.15;

	/* For a white background: */
	R = R * 0.70 + 0.05;
	G = G * 0.70 + 0.05;
	B = B * 0.70 + 0.05;

	out[0] = R;
	out[1] = G;
	out[2] = B;
}

static void del_vrml(vrml *s);

/* Constructor */
vrml *new_vrml(char *name, int doaxes) {
	vrml *s;
	double gamut_lcent = 50;

	/* Axes definition */
	struct {
		double x, y, z;
		double wx, wy, wz;
		double r, g, b;
	} axes[5] = {
		{ 0, 0,  50-gamut_lcent,  2, 2, 100,  .7, .7, .7 },	/* L axis */
		{ 50, 0,  0-gamut_lcent,  100, 2, 2,   1,  0,  0 },	/* +a (red) axis */
		{ 0, -50, 0-gamut_lcent,  2, 100, 2,   0,  0,  1 },	/* -b (blue) axis */
		{ -50, 0, 0-gamut_lcent,  100, 2, 2,   0,  1,  0 },	/* -a (green) axis */
		{ 0,  50, 0-gamut_lcent,  2, 100, 2,   1,  1,  0 },	/* +b (yellow) axis */
	};
	int i;

	if ((s = (vrml *)calloc(1, sizeof(vrml))) == NULL) {
		return NULL;
	}

	s->del                 = del_vrml;
	s->add_marker          = add_marker;
	s->add_cone            = add_cone;
	s->add_text            = add_text;
	s->start_line_set      = start_line_set;
	s->add_vertex          = add_vertex;
	s->add_col_vertex      = add_col_vertex;
	s->make_last_vertex    = make_last_vertex;
	s->add_triangle        = add_triangle;
	s->make_lines          = make_lines;
	s->make_triangles      = make_triangles;
	s->make_triangles_vc   = make_triangles_vc;
	s->make_gamut_surface  = make_gamut_surface;
	s->make_gamut_surface_2  = make_gamut_surface_2;
	s->add_cusps           = add_cusps;
	s->clear               = clear;
	s->Lab2RGB             = Lab2RGB;

	s->lcent = gamut_lcent;

	if ((s->fp = fopen(name,"w")) == NULL) {
		free(s);
		error("Malloc of vrml plot object failed");
		return NULL;
	}

	fprintf(s->fp,"#VRML V2.0 utf8\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"# Created by the Argyll CMS\n");
	fprintf(s->fp,"Transform {\n");
	fprintf(s->fp,"  children [\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"    NavigationInfo {\n");
	fprintf(s->fp,"      type \"EXAMINE\"        # It's an object we examine\n");
	fprintf(s->fp,"    } # We'll add our own light\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"    DirectionalLight {\n");
	fprintf(s->fp,"      direction 0 0 -1      # Light illuminating the scene\n");
	fprintf(s->fp,"      direction 0 -1 0      # Light illuminating the scene\n");
	fprintf(s->fp,"    }\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"    Viewpoint {\n");
	fprintf(s->fp,"      position 0 0 340      # Position we view from\n");
	fprintf(s->fp,"    }\n");
	fprintf(s->fp,"\n");
	if (doaxes != 0) {
		fprintf(s->fp,"    # Lab axes as boxes:\n");
		for (i = 0; i < 5; i++) {
			fprintf(s->fp,"    Transform { translation %f %f %f\n", axes[i].x, axes[i].y, axes[i].z);
			fprintf(s->fp,"      children [\n");
			fprintf(s->fp,"        Shape {\n");
			fprintf(s->fp,"          geometry Box { size %f %f %f }\n",
			                             axes[i].wx, axes[i].wy, axes[i].wz);
			fprintf(s->fp,"          appearance Appearance {");
			fprintf(s->fp,"            material Material { diffuseColor %f %f %f }\n", axes[i].r, axes[i].g, axes[i].b);
			fprintf(s->fp,"          }\n");
			fprintf(s->fp,"        }\n");
			fprintf(s->fp,"      ]\n");
			fprintf(s->fp,"    }\n");
		}
		fprintf(s->fp,"\n");
	}

	return s;
}

/* Finish writing the file and free ourselves */
static void del_vrml(vrml *s) {
	int i;

	fprintf(s->fp,"\n");
	fprintf(s->fp,"  ] # end of children for world\n");
	fprintf(s->fp,"}\n");

	fflush(s->fp);
	if (fclose(s->fp) != 0)
		error("VRML: Error closing VRML file\n");

	for (i = 0; i < 10; i++) {
		if (s->set[i].pary)
			free(s->set[i].pary);
	}
	if (s->tary)
		free(s->tary);
    free(s);
}

