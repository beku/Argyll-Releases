
/*
 * Simple diagnostic VRML function library for debugging
 *
 * Copyright 2005 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include "gamut.h"
#include "vrml.h"

/* Finish writing the file and free ourselves */
static void del_vrml(vrml *s) {

	fprintf(s->fp,"\n");
	fprintf(s->fp,"  ] # end of children for world\n");
	fprintf(s->fp,"}\n");

	if (fclose(s->fp) != 0)
		error("VRML: Error closing VRML file\n");
}

/* Add a shere at the given location. */
/* if col[] is NULL, use natural color. */
/* Need to do this before or after start_line_set()/dd_vertex()/make_lines() ! */
static void add_marker(vrml *s, double pos[3], double col[3]) {
	int j;
	double rad = 1.0;
	double rgb[3];

	if (col == NULL)
		s->Lab2RGB(s, rgb, pos);
	else {
		rgb[0] = col[0];
		rgb[1] = col[1];
		rgb[2] = col[2];
	}

	fprintf(s->fp,"Transform { translation %f %f %f\n", pos[1], pos[2], pos[0]-s->lcent);
	fprintf(s->fp,"\tchildren [\n");
	fprintf(s->fp,"\t\tShape{\n");
	fprintf(s->fp,"\t\t\tgeometry Sphere { radius %f}\n", rad);
	fprintf(s->fp,"\t\t\tappearance Appearance { material Material ");
	fprintf(s->fp,"{ diffuseColor %f %f %f} }\n", rgb[0], rgb[1], rgb[2]);
	fprintf(s->fp,"\t\t}\n");
	fprintf(s->fp,"\t]\n");
	fprintf(s->fp,"}\n");
}

/* Start building up verticies that will be converted to lines */
static void start_line_set(vrml *s) {
	s->npoints = 0;

	fprintf(s->fp,"\n");
	fprintf(s->fp,"Shape {\n");
	fprintf(s->fp,"  geometry IndexedLineSet { \n");
	fprintf(s->fp,"    coord Coordinate { \n");
	fprintf(s->fp,"	   point [\n");
}

/* Add a verticy with color */
static void add_col_vertex(vrml *s, double pos[3], double col[3]) {

	fprintf(s->fp,"%f %f %f,\n",pos[1], pos[2], pos[0] - s->lcent);
	
	if (s->npoints >= s->paloc) {
		s->paloc = (s->paloc + 10) * 2;
		if (s->pary == NULL)
			s->pary = malloc(s->paloc * 6 * sizeof(double));
		else
			s->pary = realloc(s->pary, s->paloc * 6 * sizeof(double));

		if (s->pary == NULL)
			error("VRML malloc failed");
	}
	s->pary[s->npoints].pp[0] = pos[0];
	s->pary[s->npoints].pp[1] = pos[1];
	s->pary[s->npoints].pp[2] = pos[2];
	s->pary[s->npoints].cc[0] = col[0];
	s->pary[s->npoints].cc[1] = col[1];
	s->pary[s->npoints].cc[2] = col[2];
	s->npoints++;
}

/* Add a color verticy */
static void add_vertex(vrml *s, double pos[3]) {
	double col[3] = { -1.0, -1.0, -1.0 };

	add_col_vertex(s, pos, col);
}

/* Convert the verticies to lines, ppset verticies per line */
static void make_lines(vrml *s, int ppset) {
	int i, j;

	fprintf(s->fp,"      ]\n");
	fprintf(s->fp,"    }\n");
	fprintf(s->fp,"  coordIndex [\n");

	for (i = 0; i < s->npoints;) {
		for (j = 0; j < ppset; j++, i++) {
			fprintf(s->fp,"%d, ", i);
		}
		fprintf(s->fp,"-1,\n");
	}
	fprintf(s->fp,"    ]\n");

	/* Color */
	fprintf(s->fp,"            colorPerVertex TRUE\n");
	fprintf(s->fp,"            color Color {\n");
	fprintf(s->fp,"              color [			# RGB colors of each vertex\n");

	for (i = 0; i < s->npoints; i++) {
		double rgb[3], Lab[3];

		if (s->pary[i].cc[0] < 0.0) {
			Lab[0] = s->pary[i].pp[0];
			Lab[1] = s->pary[i].pp[1];
			Lab[2] = s->pary[i].pp[2];
			s->Lab2RGB(s, rgb, Lab);
			fprintf(s->fp,"                %f %f %f,\n", rgb[0], rgb[1], rgb[2]);
		} else {
			fprintf(s->fp,"                %f %f %f,\n", s->pary[i].cc[0], s->pary[i].cc[1], s->pary[i].cc[2]);
		}
	}
	fprintf(s->fp,"              ] \n");
	fprintf(s->fp,"            }\n");
	/* End color */

	fprintf(s->fp,"  }\n");
	fprintf(s->fp,"} # end shape\n");
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
	int i, nverts, ix;
	int v[3];

	nverts = g->nverts(g);

	if (nverts == 0)
		return;

	fprintf(s->fp,"    Transform {\n");
	fprintf(s->fp,"      translation 0 0 0\n");
	fprintf(s->fp,"      children [\n");
	fprintf(s->fp,"		Shape { \n");
	fprintf(s->fp,"		    geometry IndexedFaceSet {\n");
	fprintf(s->fp,"				ccw FALSE\n");
	fprintf(s->fp,"				convex TRUE\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"		        coord Coordinate { \n");
	fprintf(s->fp,"		            point [			# Verticy coordinates\n");

	/* Spit out the point values, in order. */
	/* Note that a->x, b->y, L->z */
	for (ix = i = 0; ix >= 0 && i < nverts; i++) {
		double out[3];

		ix = g->getvert(g, NULL, out, ix);
		fprintf(s->fp,"%f %f %f,\n",out[1], out[2], out[0]-s->lcent);
	}
	fprintf(s->fp,"					]\n");
	fprintf(s->fp,"		        }\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"		        coordIndex [ 		# Indexes of poligon Verticies \n");

	g->startnexttri(g);
	while (g->getnexttri(g, v) == 0) {
		fprintf(s->fp,"%d, %d, %d, -1\n", v[0], v[1], v[2]);
	}
	fprintf(s->fp,"				]\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"				colorPerVertex TRUE\n");
	fprintf(s->fp,"		        color Color {\n");
	fprintf(s->fp,"		            color [			# RGB colors of each vertex\n");

	/* Spit out the colors for each vertex */
	for (ix = i = 0; ix >= 0 && i < nverts; i++) {
		double out[3];
		double rgb[3];

		ix = g->getvert(g, NULL, out, ix);

		if (cc[0] < 0.0) {
			s->Lab2RGB(s, rgb, out);
			fprintf(s->fp,"%f %f %f,\n", rgb[0], rgb[1], rgb[2]);
		} else {
			fprintf(s->fp,"%f %f %f,\n", cc[0], cc[1], cc[2]);
		}
	}
	fprintf(s->fp,"					] \n");
	fprintf(s->fp,"		        }\n");
	fprintf(s->fp,"		    }\n");
	fprintf(s->fp,"		    appearance Appearance { \n");
	fprintf(s->fp,"		        material Material {\n");
	fprintf(s->fp,"					transparency %f\n",trans);
	fprintf(s->fp,"					ambientIntensity 0.3\n");
	fprintf(s->fp,"					shininess 0.5\n");
	fprintf(s->fp,"				}\n");
	fprintf(s->fp,"		    }\n");
	fprintf(s->fp,"		}	# end Shape\n");
	fprintf(s->fp,"      ]\n");
	fprintf(s->fp,"    }\n");
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

	out[0] = R;
	out[1] = G;
	out[2] = B;
}

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

	s->del                = del_vrml;
	s->add_marker         = add_marker;
	s->start_line_set     = start_line_set;
	s->add_vertex         = add_vertex;
	s->add_col_vertex     = add_col_vertex;
	s->make_lines         = make_lines;
	s->make_gamut_surface = make_gamut_surface;
	s->Lab2RGB            = Lab2RGB;

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
	fprintf(s->fp,"children [\n");
	fprintf(s->fp,"	NavigationInfo {\n");
	fprintf(s->fp,"		type \"EXAMINE\"        # It's an object we examine\n");
	fprintf(s->fp,"	} # We'll add our own light\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"    DirectionalLight {\n");
	fprintf(s->fp,"        direction 0 0 -1      # Light illuminating the scene\n");
	fprintf(s->fp,"        direction 0 -1 0      # Light illuminating the scene\n");
	fprintf(s->fp,"    }\n");
	fprintf(s->fp,"\n");
	fprintf(s->fp,"    Viewpoint {\n");
	fprintf(s->fp,"        position 0 0 340      # Position we view from\n");
	fprintf(s->fp,"    }\n");
	fprintf(s->fp,"\n");
	if (doaxes != 0) {
		fprintf(s->fp,"# Lab axes as boxes:\n");
		for (i = 0; i < 5; i++) {
			fprintf(s->fp,"Transform { translation %f %f %f\n", axes[i].x, axes[i].y, axes[i].z);
			fprintf(s->fp,"\tchildren [\n");
			fprintf(s->fp,"\t\tShape{\n");
			fprintf(s->fp,"\t\t\tgeometry Box { size %f %f %f }\n",
			                  axes[i].wx, axes[i].wy, axes[i].wz);
			fprintf(s->fp,"\t\t\tappearance Appearance { material Material ");
			fprintf(s->fp,"{ diffuseColor %f %f %f} }\n", axes[i].r, axes[i].g, axes[i].b);
			fprintf(s->fp,"\t\t}\n");
			fprintf(s->fp,"\t]\n");
			fprintf(s->fp,"}\n");
		}
		fprintf(s->fp,"\n");
	}

	return s;
}

