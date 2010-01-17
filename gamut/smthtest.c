
/* 
 * nearsmth test code. Test the smoothed nearpoint routine.
 *
 * Author:  Graeme W. Gill
 * Date:    17/1/2002
 * Version: 1.00
 *
 * Copyright 2002, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* TTBD:
 *
 */

#undef DEBUG		/* test a single value out */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#if defined(__IBMC__) && defined(_M_IX86)
#include <float.h>
#endif

#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "gamut.h"
#include "nearsmth.h"

double m21po[3] = { 2.0, 1.0, 2.0 };    /* Many to 1 filter mixing power LCh (theoretically 2) */

/* Mapping weights */
gammapweights weights[] = {
	{
		gmm_default,

		{		/* Weighting of absolute error of destination from source */
			1.0,	/* Absolute error overall weight */
			{
				1.0,	/* Absolute luminance error weight */
				1.0,	/* Absolute chroma error weight */
				1.0		/* Absolute hue error weight */
			}
		},
		{		/* Weighting of relative error of destination points to each */
				/* other, compared to source points to each other. */
			1.0,	/* Relative error overall weight */
			{
				1.0,	/* Relative luminance error weight */
				1.0,	/* Relative chroma error weight */
				1.0		/* Relative hue error weight */
			}
		},
		{		/* Weighting of error between destination point and source */
				/* point radially mapped to destination. */
			0.0,	/* Radial error overall weight */
			{
				1.0,	/* Radial luminance error weight */
				1.0,	/* Radial chroma error weight */
				1.0		/* Radial hue error weight */
			}
		},
	
		0.0
	}
};

#define OVERSHOOT 1.0

void usage(void) {
	fprintf(stderr,"Create smoothed near mapping between two gamuts, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	fprintf(stderr,"usage: smthtest [options] ingamut outgamut diag_vrml\n");
	fprintf(stderr," -v            Verbose\n");
//	fprintf(stderr," -s nearf      Absolute delta E weighting\n");
	exit(1);
}

FILE *start_vrml(char *name, int doaxes);
void start_line_set(FILE *wrl);
void add_vertex(FILE *wrl, double pp[3]);
void make_lines(FILE *wrl, int ppset);
void end_vrml(FILE *wrl);

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	char *xl;
	char in_name[100];
	char out_name[100];
	char diag_name[100];
	int verb = 0;
	double nearf = 1.0;		/* Absolute delta E weightign */

	gamut *gin, *gout;		/* Input and Output gamuts */
	nearsmth *nsm;			/* Returned list of near smooth points */
	int nnsm;				/* Number of near smoothed points */
	FILE *wrl;				/* VRML output file */

	gammapweights xweights[14];

	int i;

#if defined(__IBMC__) && defined(_M_IX86)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
#endif

	error_program = argv[0];

	if (argc < 3)
		usage();

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
				usage();

			/* Verbosity */
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;
			}

			/* Smoothing factor */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				fa = nfa;
				if (na == NULL) usage();
				nearf = atof(na);
			}
			else 
				usage();
		} else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(in_name,argv[fa++]);

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(out_name,argv[fa++]);

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(diag_name,argv[fa++]);

	/* - - - - - - - - - - - - - - - - - - - */
	/* read the input device gamut */

	gin = new_gamut(0.0, 0);

	if ((xl = strrchr(in_name, '.')) == NULL) {	/* Add .gam extention if there isn't one */
		xl = in_name + strlen(in_name);
		strcpy(xl,".gam");
	}

	if (gin->read_gam(gin, in_name))
		error("Reading input gamut failed");

	/* - - - - - - - - - - - - - - - - - - - */
	/* read the output device gamut */

	gout = new_gamut(0.0, 0);

	if ((xl = strrchr(out_name, '.')) == NULL) { /* Add .gam extention if there isn't one */
		xl = out_name + strlen(out_name);
		strcpy(xl,".gam");
	}

	if (gout->read_gam(gout, out_name))
		error("Reading output gamut failed");

	/* - - - - - - - - - - - - - - - - - - - */

	
	/* Convert from compact to explicit hextant weightings */
	expand_weights(xweights, weights);

	/* Create the near point mapping */
	nsm = near_smooth(verb, &nnsm, gin, gin, gout, 0, 0, NULL, xweights, 0.1, 0.1, 1, 1, 2.0, 17, 0.0);
	if (nsm == NULL)
		error("Creating smoothed near points failed");

	/* Output the src to smoothed near point vectors */
	if ((xl = strrchr(diag_name, '.')) == NULL) { /* Add .wrl extention if there isn't one */
		xl = diag_name + strlen(diag_name);
		strcpy(xl,".wrl");
	}

	wrl = start_vrml(diag_name, 1);
	start_line_set(wrl);

	for (i = 0; i < nnsm; i++) {
		add_vertex(wrl, nsm[i].sv);			/* Source gamut point */
		add_vertex(wrl, nsm[i].dv);			/* Smoother destination value */

//		add_vertex(wrl, nsm[i].drv);		/* Radial points */
	} 
	make_lines(wrl, 2);
	end_vrml(wrl);

	/* Clean up */
	free_nearsmth(nsm, nnsm);

	gout->del(gout);
	gin->del(gin);

	return 0;
}

/* ------------------------------------------------ */
/* Some simple functions to do basic VRML work */

#ifndef GAMUT_LCENT
#define GAMUT_LCENT 50.0
#endif
static int npoints = 0;
static int paloc = 0;
static struct { double pp[3]; } *pary;

static void Lab2RGB(double *out, double *in);

FILE *start_vrml(char *name, int doaxes) {
	FILE *wrl;
	struct {
		double x, y, z;
		double wx, wy, wz;
		double r, g, b;
	} axes[5] = {
		{ 0, 0,  50-GAMUT_LCENT,  2, 2, 100,  .7, .7, .7 },	/* L axis */
		{ 50, 0,  0-GAMUT_LCENT,  100, 2, 2,   1,  0,  0 },	/* +a (red) axis */
		{ 0, -50, 0-GAMUT_LCENT,  2, 100, 2,   0,  0,  1 },	/* -b (blue) axis */
		{ -50, 0, 0-GAMUT_LCENT,  100, 2, 2,   0,  1,  0 },	/* -a (green) axis */
		{ 0,  50, 0-GAMUT_LCENT,  2, 100, 2,   1,  1,  0 },	/* +b (yellow) axis */
	};
	int i;
	
	if ((wrl = fopen(name,"w")) == NULL)
		error("Error opening VRML file '%s'\n",name);

	npoints = 0;

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
	fprintf(wrl,"        position 0 0 340      # Position we view from\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"\n");
	if (doaxes != 0) {
		fprintf(wrl,"# Lab axes as boxes:\n");
		for (i = 0; i < 5; i++) {
			fprintf(wrl,"Transform { translation %f %f %f\n", axes[i].x, axes[i].y, axes[i].z);
			fprintf(wrl,"\tchildren [\n");
			fprintf(wrl,"\t\tShape{\n");
			fprintf(wrl,"\t\t\tgeometry Box { size %f %f %f }\n",
			                  axes[i].wx, axes[i].wy, axes[i].wz);
			fprintf(wrl,"\t\t\tappearance Appearance { material Material ");
			fprintf(wrl,"{ diffuseColor %f %f %f} }\n", axes[i].r, axes[i].g, axes[i].b);
			fprintf(wrl,"\t\t}\n");
			fprintf(wrl,"\t]\n");
			fprintf(wrl,"}\n");
		}
		fprintf(wrl,"\n");
	}

	return wrl;
}

void
start_line_set(FILE *wrl) {

	fprintf(wrl,"\n");
	fprintf(wrl,"Shape {\n");
	fprintf(wrl,"  geometry IndexedLineSet { \n");
	fprintf(wrl,"    coord Coordinate { \n");
	fprintf(wrl,"	   point [\n");
}

void add_vertex(FILE *wrl, double pp[3]) {

	fprintf(wrl,"%f %f %f,\n",pp[1], pp[2], pp[0]-GAMUT_LCENT);
	
	if (paloc < (npoints+1)) {
		paloc = (paloc + 10) * 2;
		if (pary == NULL)
			pary = malloc(paloc * 3 * sizeof(double));
		else
			pary = realloc(pary, paloc * 3 * sizeof(double));

		if (pary == NULL)
			error ("Malloc failed");
	}
	pary[npoints].pp[0] = pp[0];
	pary[npoints].pp[1] = pp[1];
	pary[npoints].pp[2] = pp[2];
	npoints++;
}


void make_lines(FILE *wrl, int ppset) {
	int i, j;

	fprintf(wrl,"      ]\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"  coordIndex [\n");

	for (i = 0; i < npoints;) {
		for (j = 0; j < ppset; j++, i++) {
			fprintf(wrl,"%d, ", i);
		}
		fprintf(wrl,"-1,\n");
	}
	fprintf(wrl,"    ]\n");

	/* Color */
	fprintf(wrl,"            colorPerVertex TRUE\n");
	fprintf(wrl,"            color Color {\n");
	fprintf(wrl,"              color [			# RGB colors of each vertex\n");

	for (i = 0; i < npoints; i++) {
		double rgb[3], Lab[3];
		Lab[0] = pary[i].pp[0];
		Lab[1] = pary[i].pp[1];
		Lab[2] = pary[i].pp[2];
		Lab2RGB(rgb, Lab);
		fprintf(wrl,"                %f %f %f,\n", rgb[0], rgb[1], rgb[2]);
	}
	fprintf(wrl,"              ] \n");
	fprintf(wrl,"            }\n");
	/* End color */

	fprintf(wrl,"  }\n");
	fprintf(wrl,"} # end shape\n");

}

void end_vrml(FILE *wrl) {

	fprintf(wrl,"\n");
	fprintf(wrl,"  ] # end of children for world\n");
	fprintf(wrl,"}\n");

	if (fclose(wrl) != 0)
		error("Error closing VRML file\n");
}


/* Convert a gamut Lab value to an RGB value for display purposes */
static void
Lab2RGB(double *out, double *in) {
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


