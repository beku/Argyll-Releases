
/* 
 * Argyll.
 * 
 * Check for B2A table PCS->Device interpolation faults
 *
 * Author:  Graeme W. Gill
 * Date:    2000/12/11
 * Version: 1.00
 *
 * Copyright 2000 Graeme W. Gill
 * Please refer to License.txt file for details.
 */

/* TTBD:
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "numlib.h"
#include "copyright.h"
#include "config.h"
#include "icc.h"

void usage(void) {
	fprintf(stderr,"Check PCS->Device Interpolation faults of ICC file, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	fprintf(stderr,"usage: icheck [-v] [-w] infile\n");
	fprintf(stderr," -v        verbose\n");
	fprintf(stderr," -w        create VRML visualisation\n");
	fprintf(stderr," -x        Use VRML axies\n");
	exit(1);
}

FILE *start_vrml(char *name, int doaxes);
void start_line_set(FILE *wrl);
void add_vertex(FILE *wrl, double pp[3]);
void make_lines(FILE *wrl, int ppset);
void end_vrml(FILE *wrl);

int
main(
	int argc,
	char *argv[]
) {
	int fa,nfa;				/* argument we're looking at */
	int verb = 0;
	int dovrml = 0;
	int doaxes = 0;
	char in_name[100];
	char out_name[100], *xl;
	icmFile *rd_fp;
	icc *rd_icco;
	int rv = 0;

	/* Check variables */
	icmLuBase *luof, *luob;		/* A2B and B2A table lookups */
	icmLuLut *lluof, *lluob;	/* Lookup Lut type object */
	int gres;					/* Grid resolution of B2A */
	icColorSpaceSignature ins, outs;	/* Type of input and output spaces */
	int inn;							/* Number of input chanels */
	icmLuAlgType alg;
	FILE *wrl = NULL;
	
	error_program = argv[0];

	if (argc < 2)
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

			/* Verbosity */
			if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;
			}
			/* VRML */
			else if (argv[fa][1] == 'w' || argv[fa][1] == 'W') {
				dovrml = 1;
			}
			/* Axes */
			else if (argv[fa][1] == 'x' || argv[fa][1] == 'X') {
				doaxes = 1;
			}
			else if (argv[fa][1] == '?')
				usage();
			else 
				usage();
		}
		else
			break;
	}

	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(in_name,argv[fa]);

	strcpy(out_name, in_name);
	if ((xl = strrchr(out_name, '.')) == NULL)	/* Figure where extention is */
		xl = out_name + strlen(out_name);
	strcpy(xl,".wrl");

	/* Open up the file for reading */
	if ((rd_fp = new_icmFileStd_name(in_name,"r")) == NULL)
		error ("Read: Can't open file '%s'",in_name);

	if ((rd_icco = new_icc()) == NULL)
		error ("Read: Creation of ICC object failed");

	/* Read the header and tag list */
	if ((rv = rd_icco->read(rd_icco,rd_fp,0)) != 0)
		error ("Read: %d, %s",rv,rd_icco->err);

	/* Get a Device to PCS conversion object */
	if ((luof = rd_icco->get_luobj(rd_icco, icmFwd, icAbsoluteColorimetric, icSigLabData, icmLuOrdNorm)) == NULL) {
		if ((luof = rd_icco->get_luobj(rd_icco, icmFwd, icmDefaultIntent, icSigLabData, icmLuOrdNorm)) == NULL)
			error ("%d, %s",rd_icco->errc, rd_icco->err);
	}

	/* Get a PCS to Device conversion object */
	if ((luob = rd_icco->get_luobj(rd_icco, icmBwd, icAbsoluteColorimetric, icSigLabData, icmLuOrdNorm)) == NULL) {
		if ((luob = rd_icco->get_luobj(rd_icco, icmBwd, icmDefaultIntent, icSigLabData, icmLuOrdNorm)) == NULL)
			error ("%d, %s",rd_icco->errc, rd_icco->err);
	}

	/* Get details of conversion (for B2A direction) */
	luob->spaces(luob, &outs, NULL, &ins, &inn, &alg, NULL, NULL, NULL, NULL);

	if (alg != icmLutType) {
		error("Expecting Lut based profile");
	}

	if (outs != icSigLabData) {
		error("Expecting Lab PCS");
	}

	lluof = (icmLuLut *)luof;	/* Lookup Lut type object */
	lluob = (icmLuLut *)luob;	/* Lookup Lut type object */

	gres = lluob->lut->clutPoints;

	if (dovrml) {
		wrl = start_vrml(out_name, doaxes);
		start_line_set(wrl);
	}

	{
		double aerr = 0.0;
		double ccount = 0.0;
		double merr = 0.0;
		double tcount = 0.0;
		int co[3];	/* PCS grid counter */

		/* Itterate throught the PCS clut grid cells */
		for (co[2] = 0; co[2] < (gres-1); co[2]++) {
			for (co[1] = 0; co[1] < (gres-1); co[1]++) {
				for (co[0] = 0; co[0] < (gres-1); co[0]++) {
					int j, k, m;
					int cc[3];				/* Cube corner offsets */
					double pcs[8][3], wpcsd;
					double apcs[3];
					double adev[MAX_CHAN];
					double check[3];		/* Check PCS */
					double ier;				/* Interpolation error */

					apcs[0] = apcs[1] = apcs[2] = 0.0;
					for (k = 0; k < inn; k++)
						adev[k] = 0.0;

					/* For each corner of the PCS grid based at the current point, */
					/* average the PCS and Device values */
					m = 0;
					for (cc[2] = 0; cc[2] < 2; cc[2]++, m++) {
						for (cc[1] = 0; cc[1] < 2; cc[1]++) {
							for (cc[0] = 0; cc[0] < 2; cc[0]++) {
								double dev[MAX_CHAN];

								pcs[m][0] = (co[0] + cc[0])/(gres - 1.0);
								pcs[m][1] = (co[1] + cc[1])/(gres - 1.0);
								pcs[m][2] = (co[2] + cc[2])/(gres - 1.0);

								/* Match icclib settable() range */
								pcs[m][0] = pcs[m][0] * 100.0;
								pcs[m][1] = (pcs[m][1] * 254.0) - 127.0;
								pcs[m][2] = (pcs[m][2] * 254.0) - 127.0;

//printf("Input PCS %f %f %f\n", pcs[m][0], pcs[m][1], pcs[m][2]);

								/* PCS to (cliped) Device */
								if ((rv = lluob->clut(lluob, dev, pcs[m])) > 1)
									error ("%d, %s",rd_icco->errc,rd_icco->err);

								/* (clipped) Device to (clipped) PCS */
								if ((rv = lluof->clut(lluof, pcs[m], dev)) > 1)
									error ("%d, %s",rd_icco->errc,rd_icco->err);

								apcs[0] += pcs[m][0];
								apcs[1] += pcs[m][1];
								apcs[2] += pcs[m][2];

//printf("Corner PCS %f %f %f -> %f %f %f %f\n",
//pcs[m][0], pcs[m][1], pcs[m][2], dev[0], dev[1], dev[2], dev[3]);

								for (k = 0; k < inn; k++)
									adev[k] += dev[k];
							}
						}
					}

					for (j = 0; j < 3; j++)
						apcs[j] /= 8.0;

					for (k = 0; k < inn; k++)
						adev[k] /= 8.0;
	
					/* Compute worst case distance of PCS corners to average PCS */
					wpcsd = 0.0;
					for (m = 0; m < 8; m++) {
						double ss;
						for (ss = 0.0, j = 0; j < 3; j++) {
							double tt = pcs[m][j] - apcs[j];
							ss += tt * tt;
						}
						ss = sqrt(ss);
						if (ss > wpcsd)
							wpcsd = ss;
					}
					wpcsd *= 0.75;		/* Set threshold at 75% of most distant corner */
					/* Set a worst case */
					if (wpcsd < 1.0)
						wpcsd = 1.0;

//					else if (wpcsd > 3.0)
//						wpcsd = 3.0;


//printf("Average PCS %f %f %f, Average Device %f %f %f %f\n",
//apcs[0], apcs[1], apcs[2], adev[0], adev[1], adev[2], adev[3]);

					/* Average Device to PCS */
					if ((rv = lluof->clut(lluof, check, adev)) > 1)
						error ("%d, %s",rd_icco->errc,rd_icco->err);

//printf("Check PCS %f %f %f\n",
//check[0], check[1], check[2]);

					/* Compute error in PCS vs. Device interpolation */
					for (ier = 0.0, j = 0; j < 3; j++) {
						double tt = apcs[j] - check[j];
						ier += tt * tt;
					}
					ier = sqrt(ier);

//printf("Average PCS %f %f %f, Check PCS %f %f %f, error %f\n",
//apcs[0], apcs[1], apcs[2], check[0], check[1], check[2], ier);

					aerr += ier;
					ccount++;
					if (ier > merr)
						merr = ier;

					if (ier > wpcsd) {
						tcount++;

						printf("ier = %f, Dev = %f %f %f %f\n",
						       ier, adev[0], adev[1], adev[2], adev[3]);
						if (dovrml) {
							add_vertex(wrl, apcs);
							add_vertex(wrl, check);
						}
					}

//printf("~1 ier = %f\n",ier);
//printf("\n");


					if (verb)
						printf("."), fflush(stdout);
				}
			}
		}

		if (dovrml) {
			make_lines(wrl, 2);
			end_vrml(wrl);
		}
	
		aerr /= ccount;
	
		printf("Average interpolation error %f, maximum %f\n",aerr, merr);
		printf("Number outside corner radius = %f%%\n",tcount * 100.0/ccount);
	}

	/* Done with lookup objects */
	luof->del(luof);
	luob->del(luob);

	rd_icco->del(rd_icco);
	rd_fp->del(rd_fp);

	return 0;
}

/* ------------------------------------------------ */
/* Some simple functions to do basix VRML work */

#define GAMUT_LCENT 50.0
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

