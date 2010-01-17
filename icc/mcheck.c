
/* 
 * International Color Consortium Format Library (icclib)
 * Check the device chanel to PCS monotonicity.
 *
 * Author:  Graeme W. Gill
 * Date:    2000/12/11
 * Version: 2.12
 *
 * Copyright 2000 - 2005 Graeme W. Gill
 *
 * This material is licensed with an "MIT" free use license:-
 * see the License.txt file in this directory for licensing details.
 */

/* TTBD:
 *
 * Make general device input, not just CMYK
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "icc.h"

void error(char *fmt, ...), warning(char *fmt, ...);

void usage(void) {
	fprintf(stderr,"Check device to PCS monotonicity of a CMYK ICC file, V%s\n",ICCLIB_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill\n");
	fprintf(stderr,"usage: kcheck [-v] [-w] infile\n");
	fprintf(stderr," -v        verbose\n");
	fprintf(stderr," -c        Check just Cyan monotonicity\n");
	fprintf(stderr," -m        Check just Magenta monotonicity\n");
	fprintf(stderr," -y        Check just Yellow monotonicity\n");
	fprintf(stderr," -k        Check just Black monotonicity\n");
	fprintf(stderr," -w        create VRML visualisation\n");
	exit(1);
}

#define MGR 50		/* Maximum grid resolution handled */


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
	int cchan = -1;			/* default all */
	int dovrml = 0;
	int doaxes = 0;
	char in_name[500];
	char out_name[500], *xl;
	icmFile *rd_fp;
	icc *wr_icco, *rd_icco;		/* Keep object separate */
	int rv = 0;

	/* Check variables */
	icmLuBase *luo;
	icmLuLut *luluto;	/* Lookup xLut type object */
	int gres;			/* Grid resolution */
	icColorSpaceSignature ins, outs;	/* Type of input and output spaces */
	int inn;							/* Number of input chanels */
	icmLuAlgType alg;
	FILE *wrl;
	int dx[4];			/* Device index mapping */
	int chan, cs, ce;
	
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
			/* Cyan */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				cchan = 0;
			}
			/* Magenta */
			else if (argv[fa][1] == 'm' || argv[fa][1] == 'M') {
				cchan = 1;
			}
			/* Yellow */
			else if (argv[fa][1] == 'y' || argv[fa][1] == 'Y') {
				cchan = 2;
			}
			/* Black */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				cchan = 3;
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
	if ((luo = rd_icco->get_luobj(rd_icco, icmFwd, icRelativeColorimetric, icSigLabData, icmLuOrdNorm)) == NULL) {
		if ((luo = rd_icco->get_luobj(rd_icco, icmFwd, icmDefaultIntent, icSigLabData, icmLuOrdNorm)) == NULL)
			error ("%d, %s",rd_icco->errc, rd_icco->err);
	}
	/* Get details of conversion */
	luo->spaces(luo, &ins, &inn, &outs, NULL, &alg, NULL, NULL, NULL);

	if (alg != icmLutType) {
		error("Expecting Lut based profile");
	}

	if (ins != icSigCmykData) {
		error("Expecting CMYK device");
	}

	if (outs != icSigLabData) {
		error("Expecting Lab PCS");
	}

	luluto = (icmLuLut *)luo;	/* Lookup xLut type object */

	gres = luluto->lut->clutPoints;
	if (gres > MGR) {
		error("Can't handle grid resolution greater than %d\n",MGR);
	}

	if (dovrml) {
		wrl = start_vrml(out_name, doaxes);
		start_line_set(wrl);
	}

	/* For all the device chanels chosen */
	if (cchan < 0) {
		cs = 0;
		ce = inn;
	} else {
		cs = cchan;
		ce = cs + 1;
	}
	for (chan = cs; chan < ce; chan++) {

		/* Check the monotonicity of the output for a given device input */
		int co[4];
		if (chan == 0) {
			dx[0] = 1;
			dx[1] = 2;
			dx[2] = 3;
			dx[3] = 0;		/* Cyan is variable */
		} else if (chan == 1) {
			dx[0] = 0;
			dx[1] = 2;
			dx[2] = 3;
			dx[3] = 1;		/* Magenta is variable */
		} else if (chan == 2) {
			dx[0] = 0;
			dx[1] = 1;
			dx[2] = 3;
			dx[3] = 2;		/* Yellow is variable */
		} else if (chan == 3) {
			dx[0] = 0;
			dx[1] = 1;
			dx[2] = 2;
			dx[3] = 3;		/* Black is variable */
		}

		/* Itterate throught the CMY clut grid points */
		for (co[0] = 0; co[0] < gres; co[0]++) {
			for (co[1] = 0; co[1] < gres; co[1]++) {
				for (co[2] = 0; co[2] < gres; co[2]++) {
					int j, k, ck, nm;
					double dev[MGR][4];
					double pcs[MGR][3];
					double apcs[3], ss;

					/* Run up the variable axis */
					for (ck = 0; ck < gres; ck++) {

						dev[ck][dx[0]] = co[0]/(gres-1.0);
						dev[ck][dx[1]] = co[1]/(gres-1.0);
						dev[ck][dx[2]] = co[2]/(gres-1.0);
						dev[ck][dx[3]] = ck/(gres-1.0);

						/* Device to PCS */
						if ((rv = luluto->clut(luluto, pcs[ck], dev[ck])) > 1)
							error ("%d, %s",rd_icco->errc,rd_icco->err);

//						if (dovrml)
//							add_vertex(wrl, pcs[ck]);
					}

					/* Compute average vector direction */
					for (ss = 0.0, k = 0; k < 3; k++) {
						double tt;
						tt = pcs[gres-1][k] - pcs[0][k];
						ss += tt * tt;
						apcs[k] = tt;
					}
					for (k = 0; k < 3; k++)
						apcs[k] /= ss;

					/* Now compute the dot product for each vector, */
					/* and check for reversals. */
					j = 0;
//printf("Checking          CMYK %f %f %f %f Lab %f %f %f\n",
//       dev[j][0], dev[j][1], dev[j][2], dev[j][3],
//       pcs[j][0], pcs[j][1], pcs[j][2]);
					for (nm = 0, j = 1; j < gres; j++) {
						for (ss = 0.0, k = 0; k < 3; k++)	/* Dot product */
							ss += (pcs[j][k] - pcs[j-1][k]) * apcs[k];

//printf("Checking %f CMYK %f %f %f %f Lab %f %f %f\n",
//       ss, dev[j][0], dev[j][1], dev[j][2], dev[j][3],
//       pcs[j][0], pcs[j][1], pcs[j][2]);

						if (ss <= 0.0) {
							nm = 1;
							printf("NonMon %f at CMYK %f %f %f %f Lab %f %f %f\n",
							       ss, dev[j][0], dev[j][1], dev[j][2], dev[j][3],
							       pcs[j][0], pcs[j][1], pcs[j][2]);
						}
					}
//printf("\n");

					/* Display just the non mono threads */
					if (nm && dovrml) {
						for (j = 0; j < gres; j++)
							add_vertex(wrl, pcs[j]);
					}
					if (verb) {
						printf("."); fflush(stdout);
					}
				}
			}
		}
	}

	if (dovrml) {
		make_lines(wrl, gres);
		end_vrml(wrl);
	}

	/* Done with lookup object */
	luo->del(luo);

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


/* ------------------------------------------------ */
/* Basic printf type error() and warning() routines */

void
error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"icctest: Error - ");
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

	fprintf(stderr,"icctest: Warning - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}
