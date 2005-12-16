
/* 
 * viewgam
 *
 * Gamut support routines.
 *
 * Author:  Graeme W. Gill
 * Date:    4/10/00
 * Version: 1.00
 *
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "gamut.h"
#include "cgats.h"

/* 
	This program reads one or more CGATS format triangular gamut
	surface descriptions, and combines them into a VRML file,
	so that the gamuts can be visually compared.

 */

/* TTBD:
 *
 */

#undef HALF_HACK /* 27.0 */

#ifdef NEVER
void error(char *fmt, ...), warning(char *fmt, ...);
#endif

void usage(char *diag, ...) {
	fprintf(stderr,"View gamuts Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: viewgam { [-c color] [-t trans] [-w|s] infile.gam } ... outfile.wrl\n");
	fprintf(stderr," -c color      Color to make gamut, r = red, g = green, b = blue\n"); 
	fprintf(stderr,"               c = cyan, m = magenta, y = yellow, e = grey, w = white\n"); 
	fprintf(stderr,"               n = natural color\n"); 
	fprintf(stderr," -t trans      Set transparency from 0.0 (opaque) to 1.0 (invisible)\n"); 
	fprintf(stderr," -w            Show as a wireframe\n");
	fprintf(stderr," -s            Show as a solid surace\n");
	fprintf(stderr," infile.gam    Name of .gam file\n");
	fprintf(stderr,"               Repeat above for each input file\n\n");
	fprintf(stderr," -n            Don't add Lab axes\n");
	fprintf(stderr," outfile.wrl   Name of output .wrl file\n");
	fprintf(stderr,"\n");
	exit(1);
}

#define MXGAMTS 5		/* Maximum number of gamuts that can be displayed */

#define GCENT 50.0

typedef enum {
	gam_red       = 0,
	gam_green     = 1,
	gam_blue      = 2,
	gam_cyan      = 3,
	gam_magenta   = 4,	
	gam_yellow    = 5,
	gam_grey      = 6,
	gam_white     = 7,
	gam_natural   = 8
} gam_colors;

struct {
	double r, g, b;
} color_rgb[8] = {
	{ 1, 0, 0 },	/* gam_red */
	{ 0, 1, 0 },	/* gam_green */
	{ 0, 0, 1 },	/* gam_blue */
	{ 0, 1, 1 },	/* gam_cyan */
	{ 1, 0, 1 },	/* gam_magenta */
	{ 1, 1, 0 },	/* gam_yellow */
	{ .1, .1, .1 },	/* gam_grey */
	{ .7, .7, .7 }	/* gam_white */
};

typedef enum {
	gam_solid     = 0,
	gam_wire      = 1,
	gam_points    = 2
} gam_reps;

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	int n, ng = 0;			/* Current number of input gamuts */
	char in_name[MXGAMTS+1][100];
	gam_colors in_colors[MXGAMTS+1];	/* Color enum for each input */
	double in_trans[MXGAMTS+1];		/* Transparency for each input */
	gam_reps in_rep[MXGAMTS+1];		/* Representation enum for each input */
	int doaxes = 1;
	FILE *wrl;
	char *xl, out_name[100];
	FILE *fp;
	int verb = 0;
	int rv = 0;
	struct {
		double x, y, z;
		double wx, wy, wz;
		double r, g, b;
	} axes[5] = {
		{ 0, 0,  50-GCENT,  2, 2, 100,  .7, .7, .7 },	/* L axis */
		{ 50, 0,  0-GCENT,  100, 2, 2,   1,  0,  0 },	/* +a (red) axis */
		{ 0, -50, 0-GCENT,  2, 100, 2,   0,  0,  1 },	/* -b (blue) axis */
		{ -50, 0, 0-GCENT,  100, 2, 2,   0,  1,  0 },	/* -a (green) axis */
		{ 0,  50, 0-GCENT,  2, 100, 2,   1,  1,  0 },	/* +b (yellow) axis */
	};

	if (argc < 3)
		usage("Too few arguments, got %d expect at least 2",argc-1);

	/* Set some defaults */
	for (n = 0; n < MXGAMTS; n++) {
		switch(n) {
			case 0:
				in_colors[0] = gam_natural;
				in_rep[0]   = gam_solid;
				in_trans[0]  = 0.0;
				break;
			case 1:
				in_colors[1] = gam_white;
				in_rep[1]   = gam_wire;
				in_trans[1]  = 0.0;
				break;
			case 2:
				in_colors[2] = gam_red;
				in_rep[2]   = gam_wire;
				in_trans[3]  = 0.0;
				break;
			case 3:
				in_colors[3] = gam_cyan;
				in_rep[3]   = gam_wire;
				in_trans[3]  = 0.0;
				break;
			case 4:
				in_colors[4] = gam_yellow;
				in_rep[4]   = gam_wire;
				in_trans[4]  = 0.2;
				break;
			case 5:
				in_colors[5] = gam_green;
				in_rep[5]   = gam_wire;
				in_trans[5]  = 0.4;
				break;
			default:
				in_colors[n] = gam_blue;
				in_rep[n]   = gam_wire;
				in_trans[n]  = 0.6;
				break;
		}
	}

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')	{	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage("Usage requested");

			/* Color */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage("Expect argument after flag -c");
    			switch (na[0]) {
					case 'r':
					case 'R':
						in_colors[ng] = gam_red;
						break;
					case 'g':
					case 'G':
						in_colors[ng] = gam_green;
						break;
					case 'b':
					case 'B':
						in_colors[ng] = gam_blue;
						break;
					case 'c':
					case 'C':
						in_colors[ng] = gam_cyan;
						break;
					case 'm':
					case 'M':
						in_colors[ng] = gam_magenta;
						break;
					case 'y':
					case 'Y':
						in_colors[ng] = gam_yellow;
						break;
					case 'e':
					case 'E':
						in_colors[ng] = gam_grey;
						break;
					case 'w':
					case 'W':
						in_colors[ng] = gam_white;
						break;
					case 'n':
					case 'N':
						in_colors[ng] = gam_natural;
						break;
					default:
						usage("Unknown argument after flag -c '%c'",na[0]);
				}
			}

			/* Transparency */
			else if (argv[fa][1] == 't' || argv[fa][1] == 'T') {
				double v;
				fa = nfa;
				if (na == NULL) usage("Expect argument after flag -t");
				v = atof(na);
				if (v < 0.0)
					v = 0.0;
				else if (v > 1.0)
					v = 1.0;
				in_trans[ng] = v;
			}

			/* Solid output */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				in_rep[ng] = gam_solid;
			}

			/* Wireframe output */
			else if (argv[fa][1] == 'w' || argv[fa][1] == 'W') {
				in_rep[ng] = gam_wire;
			}

			/* No axis output */
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N') {
				doaxes = 0;
			}

			else 
				usage("Unknown flag '%c'",argv[fa][1]);

		} else if (argv[fa][0] != '\000') { /* Got a non-flag */
			strcpy(in_name[ng], argv[fa]);
			ng++;
			if (ng >= MXGAMTS)
				break;

		} else {
			break;
		}
	}

	/* The last "gamut" is actually the output VRML filename, */
	/* so unwind it. */

	if (ng < 2)
		usage("Not enough arguments to specify output VRML files");

	strcpy(out_name, in_name[--ng]);


#ifdef NEVER
	for (n = 0; n < ng; n++) {
		printf("Input file %d is '%s'\n",n,in_name[n]);
		printf("Input file %d has color %d\n",n,in_colors[n]);
		printf("Input file %d has rep %d\n",n,in_rep[n]);
		printf("Input file %d has trans %f\n",n,in_trans[n]);
		
	}
	printf("Output file is '%s'\n",out_name);
#endif	/* NEVER */

	/* Open up the output file */
	if ((wrl = fopen(out_name,"w")) == NULL)
		error("Error opening output file '%s'\n",out_name);
	
	/* Write the header info */

	fprintf(wrl,"#VRML V2.0 utf8\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"# Created by the Argyll CMS\n");
	fprintf(wrl,"Transform {\n");
  	fprintf(wrl,"  children [\n");
    fprintf(wrl,"    NavigationInfo {\n");
	fprintf(wrl,"      type \"EXAMINE\"        # It's an object we examine\n");
	fprintf(wrl,"    } # We'll add our own light\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    DirectionalLight {\n");
	fprintf(wrl,"      direction 0 0 -1      # Light illuminating the scene\n");
	fprintf(wrl,"      direction 0 -1 0      # Light illuminating the scene\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");
	fprintf(wrl,"    Viewpoint {\n");
	fprintf(wrl,"      position 0 0 340      # Position we view from\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");
	if (doaxes) {
		fprintf(wrl,"    # Lab axes as boxes:\n");
		for (n = 0; n < 5; n++) {
			fprintf(wrl,"    Transform { translation %f %f %f\n", axes[n].x, axes[n].y, axes[n].z);
			fprintf(wrl,"      children [\n");
			fprintf(wrl,"        Shape{\n");
			fprintf(wrl,"          geometry Box { size %f %f %f }\n",
			                       axes[n].wx, axes[n].wy, axes[n].wz);
			fprintf(wrl,"          appearance Appearance { material Material ");
			fprintf(wrl,"{ diffuseColor %f %f %f} }\n", axes[n].r, axes[n].g, axes[n].b);
			fprintf(wrl,"        }\n");
			fprintf(wrl,"      ]\n");
			fprintf(wrl,"    }\n");
		}
	}

#ifdef NEVER

Add in labels:

# Axe descriptions L* a* b* and O

    Transform { translation -5.0 115.0 -50.0
            children[
                  Shape{
                         geometry Text {        string ["+B*"]
                             fontStyle FontStyle { family "SANS"
                                                   style "BOLD"
                                                   size 10.0 } }
                         appearance Appearance { material Material { diffuseColor 0.7 0.7 0.7 } }
                        }
                    ]
               }

    Transform { translation -5.0 -120.0 -50.0
            children[
                  Shape{
                         geometry Text {        string ["-B*"]
                             fontStyle FontStyle { family "SANS"
                                                   style "BOLD"
                                                   size 10.0 } }
                         appearance Appearance { material Material { diffuseColor 0.7 0.7 0.7 } }
                        }
                    ]
               }

    Transform { translation 115.0 -3.0 -50.0
            children[
                  Shape{
                         geometry Text {        string ["+A*"]
                             fontStyle FontStyle { family "SANS"
                                                   style "BOLD"
                                                   size 10.0 } }
                         appearance Appearance { material Material { diffuseColor 0.7 0.7 0.7 } }
                        }
                    ]
               }

    Transform { translation -125.0 -3.0 -50.0
            children[
                  Shape{
                         geometry Text {        string ["-A*"]
                             fontStyle FontStyle { family "SANS"
                                                   style "BOLD"
                                                   size 10.0 } }
                         appearance Appearance { material Material { diffuseColor 0.7 0.7 0.7 } }
                        }
                    ]
               }

    Transform { translation -2.0 2.0 56.0
                rotation 0.85 -0.35 -0.35 1.76269
            children[
                  Shape{
                         geometry Text {        string ["+L*"]
                             fontStyle FontStyle { family "SANS"
                                                   style "BOLD"
                                                   size 10.0 } }
                         appearance Appearance { material Material { diffuseColor 0.7 0.7 0.7 } }
                        }
                    ]
               }

    Transform { translation -2.0 2.0 -60.0
                rotation 0.85 -0.35 -0.35 1.76269
            children[
                  Shape{
                         geometry Text {        string ["0"]
                             fontStyle FontStyle { family "SANS"
                                                   style "BOLD"
                                                   size 10.0 } }
                         appearance Appearance { material Material { diffuseColor 0.7 0.7 0.7 } }
                        }
                    ]
               }


#endif /* NEVER */

	/* Read each input in turn */
	for (n = 0; n < ng; n++) {
		int i;
		cgats *pp;
		int nverts;
		int ntris;
		int Lf, af, bf;			/* Fields holding L, a & b data */
		int v0f, v1f, v2f;		/* Fields holding verticies 0, 1 & 2 */

		pp = new_cgats();	/* Create a CGATS structure */
	
		/* Setup to cope with a gamut file */
		pp->add_other(pp, "GAMUT");
	
		if (pp->read_name(pp, in_name[n]))
			error("Input file '%s' error : %s",in_name[n], pp->err);
	
		if (pp->t[0].tt != tt_other || pp->t[0].oi != 0)
			error("Input file isn't a GAMUT format file");
		if (pp->ntables != 2)
			error("Input file doesn't contain exactly two tables");

		if ((nverts = pp->t[0].nsets) <= 0)
			error("No verticies");
		if ((ntris = pp->t[1].nsets) <= 0)
			error("No triangles");

		if ((Lf = pp->find_field(pp, 0, "LAB_L")) < 0)
			error("Input file doesn't contain field LAB_L");
		if (pp->t[0].ftype[Lf] != r_t)
			error("Field LAB_L is wrong type");
		if ((af = pp->find_field(pp, 0, "LAB_A")) < 0)
			error("Input file doesn't contain field LAB_A");
		if (pp->t[0].ftype[af] != r_t)
			error("Field LAB_A is wrong type");
		if ((bf = pp->find_field(pp, 0, "LAB_B")) < 0)
			error("Input file doesn't contain field LAB_B");
		if (pp->t[0].ftype[bf] != r_t)
			error("Field LAB_B is wrong type");

		/* Write the vertexes out */
		fprintf(wrl,"\n");
		fprintf(wrl,"    Transform {\n");
		fprintf(wrl,"      translation 0 0 0\n");
		fprintf(wrl,"      children [\n");
		fprintf(wrl,"        Shape { \n");
		if (in_rep[n] == gam_wire) {
			fprintf(wrl,"          geometry IndexedLineSet {\n");
		} else {
			fprintf(wrl,"          geometry IndexedFaceSet {\n");
			fprintf(wrl,"            ccw FALSE\n");
			fprintf(wrl,"            convex TRUE\n");
		}
		fprintf(wrl,"\n");
		fprintf(wrl,"            coord Coordinate { \n");
		fprintf(wrl,"              point [			# Verticy coordinates\n");

		/* Spit out the point values, in order. */
		/* Note that a->x, b->y, L->z */
		for (i = 0; i < nverts; i++) {
			double L, a, b;
			L = *((double *)pp->t[0].fdata[i][Lf]);
			a = *((double *)pp->t[0].fdata[i][af]);
			b = *((double *)pp->t[0].fdata[i][bf]);
			fprintf(wrl,"                %f %f %f,\n",a, b, L - GCENT);
		}
		fprintf(wrl,"              ]\n");
		fprintf(wrl,"            }\n");
		fprintf(wrl,"\n");

		/* Write the triangles/wires out */
		if ((v0f = pp->find_field(pp, 1, "VERTEX_0")) < 0)
			error("Input file doesn't contain field VERTEX_0");
		if (pp->t[1].ftype[v0f] != i_t)
			error("Field VERTEX_0 is wrong type");
		if ((v1f = pp->find_field(pp, 1, "VERTEX_1")) < 0)
			error("Input file doesn't contain field VERTEX_1");
		if (pp->t[1].ftype[v1f] != i_t)
			error("Field VERTEX_1 is wrong type");
		if ((v2f = pp->find_field(pp, 1, "VERTEX_2")) < 0)
			error("Input file doesn't contain field VERTEX_2");
		if (pp->t[1].ftype[v2f] != i_t)
			error("Field VERTEX_2 is wrong type");

		fprintf(wrl,"            coordIndex [ 		# Indexes of poligon Verticies \n");

		for (i = 0; i < ntris; i++) {
			int v0, v1, v2;
			v0 = *((int *)pp->t[1].fdata[i][v0f]);
			v1 = *((int *)pp->t[1].fdata[i][v1f]);
			v2 = *((int *)pp->t[1].fdata[i][v2f]);

#ifdef HALF_HACK 
			if (*((double *)pp->t[0].fdata[v0][Lf]) < HALF_HACK
			 || *((double *)pp->t[0].fdata[v1][Lf]) < HALF_HACK
			 || *((double *)pp->t[0].fdata[v2][Lf]) < HALF_HACK)
				continue;
#endif /* HALF_HACK */

			if (in_rep[n] == gam_wire) {
				if (v0 < v1)				/* Only output 1 wire of two on an edge */
					fprintf(wrl,"              %d, %d, -1\n", v0, v1);
				if (v1 < v2)
					fprintf(wrl,"              %d, %d, -1\n", v1, v2);
				if (v2 < v0)
					fprintf(wrl,"              %d, %d, -1\n", v2, v0);
			} else {
				fprintf(wrl,"              %d, %d, %d, -1\n", v0, v1, v2);
			}
		}
		fprintf(wrl,"            ]\n");
		fprintf(wrl,"\n");

		/* Write the colors out */
		if (in_colors[n] == gam_natural) {
			fprintf(wrl,"            colorPerVertex TRUE\n");
			fprintf(wrl,"            color Color {\n");
			fprintf(wrl,"              color [			# RGB colors of each vertex\n");

			for (i = 0; i < nverts; i++) {
				double rgb[3], Lab[3];
				Lab[0] = *((double *)pp->t[0].fdata[i][Lf]);
				Lab[1] = *((double *)pp->t[0].fdata[i][af]);
				Lab[2] = *((double *)pp->t[0].fdata[i][bf]);
				gamut_Lab2RGB(rgb, Lab);
				fprintf(wrl,"                %f %f %f,\n", rgb[0], rgb[1], rgb[2]);
			}
			fprintf(wrl,"              ] \n");
			fprintf(wrl,"            }\n");
		}
		fprintf(wrl,"          }\n");
		fprintf(wrl,"          appearance Appearance { \n");
		fprintf(wrl,"            material Material {\n");
		if (in_trans[n] > 0.0) {
			fprintf(wrl,"              transparency %f\n", in_trans[n]);
		}
		fprintf(wrl,"              ambientIntensity 0.3\n");
		fprintf(wrl,"              shininess 0.5\n");
		if (in_colors[n] != gam_natural) {
			fprintf(wrl,"              emissiveColor %f %f %f\n",
			   color_rgb[in_colors[n]].r, color_rgb[in_colors[n]].g, color_rgb[in_colors[n]].b);
		}
		fprintf(wrl,"            }\n");
		fprintf(wrl,"          }\n");
		fprintf(wrl,"        }	# end Shape\n");
		fprintf(wrl,"      ] # end children\n");
		fprintf(wrl,"    } # end Transform\n");
		fprintf(wrl,"\n");

		pp->del(pp);		/* Clean up */
	}

	/* Write the trailer */
	fprintf(wrl,"  ] # end of children for world\n");
	fprintf(wrl,"}\n");

	/* Close the file */
	fclose(wrl);

	return 0;
}

#ifdef NEVER
/* Basic printf type error() and warning() routines */

void
error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"icclu: Error - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	exit (-1);
}

void
warning(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"icclu: Warning - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}

#endif /* NEVER */
