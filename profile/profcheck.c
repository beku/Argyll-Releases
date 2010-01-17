
/* 
 * Argyll Color Correction System
 * Color Device profile checker.
 *
 * Author: Graeme W. Gill
 * Date:   15/7/2001
 *
 * Copyright 2001 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes in the .ti3 scattered test chart
 * points, and checks them against an ICC profile.
 * forward ICC device profile.
 */

/*
 * TTBD:
 *		Switch to generic colorant read code rather than Grey/RGB/CMYK,
 *		and allow checking ICC profiles > 4 colors
 */

#undef DEBUG

#define IMP_MONO			/* Turn on development code */

#define verbo stdout

#include <stdio.h>
#include <string.h>
#include <math.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "cgats.h"
#include "xicc.h"
#include "insttypes.h"
#include "sort.h"

void
usage(void) {
	fprintf(stderr,"Check accuracy of ICC profile, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	fprintf(stderr,"usage: profcheck [-options] data.ti3 iccprofile.icm\n");
	fprintf(stderr," -v [level]      Verbosity level (default 1), 2 to print each DE\n");
	fprintf(stderr," -c              Show CIE94 delta E values\n");
	fprintf(stderr," -k              Show CIEDE2000 delta E values\n");
	fprintf(stderr," -w              create VRML visualisation (iccprofile.wrl)\n");
	fprintf(stderr," -x              Use VRML axes\n");
	fprintf(stderr," -m              Make VRML lines a minimum of 0.5\n");
	fprintf(stderr," -e              Color vectors acording to delta E\n");
	fprintf(stderr," -d devval1,deval2,devvalN\n");
	fprintf(stderr,"                 Specify a device value to sort against\n");
	fprintf(stderr," -p              Sort device value by PCS (Lab) target\n");
	fprintf(stderr," -i illum        Choose illuminant for print/transparency spectral data:\n");
	fprintf(stderr,"                 A, C, D50 (def.), D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ       Choose CIE Observer for spectral data:\n");
	fprintf(stderr,"                 1931_2 (def), 1964_10, S&B 1955_2, shaw, J&V 1978_2\n");
	fprintf(stderr," -f              Use Fluorescent Whitening Agent compensation\n");
	fprintf(stderr," data.ti3        Test data file\n");
	fprintf(stderr," iccprofile.icm  Profile to check against\n");
	exit(1);
	}

FILE *start_vrml(char *name, int doaxes);
void start_line_set(FILE *wrl);
void add_vertex(FILE *wrl, double pp[3]);
void make_lines(FILE *wrl, int ppset);
void make_de_lines(FILE *wrl);
void end_vrml(FILE *wrl);

/* Patch value type */
typedef struct {
	char sid[50];		/* sample id */
	double p[MAX_CHAN];	/* Device value */
	double v[3];		/* CIE value */
	double dp;			/* Delta from target value */
	double dv;			/* Delta from CIE value */
} pval;

int main(int argc, char *argv[])
{
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	int cie94 = 0;
	int cie2k = 0;
	int dovrml = 0;
	int dominl = 0;
	int doaxes = 0;
	int dodecol = 0;
	char ti3name[MAXNAMEL+1] = { 0 };	/* Input cgats file base name */
	cgats *icg;				/* input cgats structure */
	char iccname[MAXNAMEL+1] = { 0 };	/* Input icc file base name */
	icmFile *rd_fp;
	icc *rd_icco;
	icmLuBase *luo;
	char out_name[MAXNAMEL+1], *xl;		/* VRML name */
	FILE *wrl = NULL;

	int fwacomp = 0;			/* FWA compensation */
	int isdisp = 0;				/* nz if this is a display device, 0 if output */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType illum = icxIT_D50;	/* Spectral defaults */
	xspect cust_illum;			/* Custom illumination spectrum */
	icxObserverType observ = icxOT_CIE_1931_2;

	int ddevv = 0;				/* Do device value sort */
	double devval[MAX_CHAN];	/* device value to sort on */
	int sortbypcs = 0;			/* Sort by PCS */

	int npat;					/* Number of patches */
	pval *tpat;					/* Patch input values */
	int i, j, rv = 0;
	icColorSpaceSignature devspace = 0;	/* The device colorspace */
	int isAdditive = 0;			/* 0 if subtractive, 1 if additive colorspace */
	int isLab = 0;				/* 0 if input is XYZ, 1 if input is Lab */
	int devchan = 0;			/* Number of device chanels */

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	error_program = "profcheck";

	if (argc <= 1)
		usage();

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {

		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {	/* Look for any flags */
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
				if (na != NULL && isdigit(na[0])) {
					verb = atoi(na);
				}
			}

			/* VRML */
			else if (argv[fa][1] == 'w' || argv[fa][1] == 'W')
				dovrml = 1;

			/* Minimum line length */
			else if (argv[fa][1] == 'm' || argv[fa][1] == 'M')
				dominl = 1;

			/* Axes */
			else if (argv[fa][1] == 'x' || argv[fa][1] == 'X')
				doaxes = 1;

			/* Delta E coloring */
			else if (argv[fa][1] == 'e' || argv[fa][1] == 'E')
				dodecol = 1;

			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				cie94 = 1;
				cie2k = 0;
			}

			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				cie94 = 0;
				cie2k = 1;
			}

			/* Device sort value */
			else if (argv[fa][1] == 'd' || argv[fa][1] == 'D') {
				char *tp, buf[200];
				int ndv;
				fa = nfa;
				if (na == NULL) usage();

				ddevv = 1;
				strcpy(buf, na);
				
				/* Replace ',' with '\000' */
				for (ndv = 1,tp = buf; *tp != '\000'; tp++) {
					if (*tp == ',') {
						*tp = '\000';
						ndv++;
					}
				}
				if (ndv >= MAX_CHAN)
					ndv = MAX_CHAN;

				for (tp = buf, i = 0; i < ndv; i++, tp += strlen(tp) + 1) {
					devval[i] = atof(tp);
				}
			}

			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P')
				sortbypcs = 1;

			/* Spectral Illuminant type */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL) usage();
				if (strcmp(na, "A") == 0) {
					spec = 1;
					illum = icxIT_A;
				} else if (strcmp(na, "C") == 0) {
					spec = 1;
					illum = icxIT_C;
				} else if (strcmp(na, "D50") == 0) {
					spec = 1;
					illum = icxIT_D50;
				} else if (strcmp(na, "D65") == 0) {
					spec = 1;
					illum = icxIT_D65;
				} else if (strcmp(na, "F5") == 0) {
					spec = 1;
					illum = icxIT_F5;
				} else if (strcmp(na, "F8") == 0) {
					spec = 1;
					illum = icxIT_F8;
				} else if (strcmp(na, "F10") == 0) {
					spec = 1;
					illum = icxIT_F10;
				} else {	/* Assume it's a filename */
					spec = 1;
					illum = icxIT_custom;
					if (read_xspect(&cust_illum, na) != 0)
						usage();
				}
			}

			/* Spectral Observer type */
			else if (argv[fa][1] == 'o' || argv[fa][1] == 'O') {
				fa = nfa;
				if (na == NULL) usage();
				if (strcmp(na, "1931_2") == 0) {			/* Classic 2 degree */
					spec = 1;
					observ = icxOT_CIE_1931_2;
				} else if (strcmp(na, "1964_10") == 0) {	/* Classic 10 degree */
					spec = 1;
					observ = icxOT_CIE_1964_10;
				} else if (strcmp(na, "1955_2") == 0) {		/* Stiles and Burch 1955 2 degree */
					spec = 1;
					observ = icxOT_Stiles_Burch_2;
				} else if (strcmp(na, "1978_2") == 0) {		/* Judd and Voss 1978 2 degree */
					spec = 1;
					observ = icxOT_Judd_Voss_2;
				} else if (strcmp(na, "shaw") == 0) {		/* Shaw and Fairchilds 1997 2 degree */
					spec = 1;
					observ = icxOT_Shaw_Fairchild_2;
				} else
					usage();
			}

			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F')
				fwacomp = 1;

			else 
				usage();
		} else
			break;
	}

	/* Get the file name arguments */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(ti3name,argv[fa++],MAXNAMEL); ti3name[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(iccname,argv[fa++],MAXNAMEL); iccname[MAXNAMEL] = '\000';

	strncpy(out_name,iccname,MAXNAMEL-4); out_name[MAXNAMEL-4] = '\000';
	if ((xl = strrchr(out_name, '.')) == NULL)	/* Figure where extention is */
		xl = out_name + strlen(out_name);
	strcpy(xl,".wrl");

	if (fwacomp && spec == 0)
		error ("FWA compensation only works when viewer and/or illuminant selected");

	/* Open and look at the .ti3 profile patches file */
	icg = new_cgats();			/* Create a CGATS structure */
	icg->add_other(icg, "CTI3"); 	/* our special input type is Calibration Target Information 3 */

	if (icg->read_name(icg, ti3name))
		error("CGATS file read error on '%s': %s",ti3name,icg->err);

	if (icg->ntables == 0 || icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file '%s' isn't a CTI3 format file",ti3name);
	if (icg->ntables < 1)
		error ("Input file '%s' doesn't contain at least one table",ti3name);

	/* See if CIE is actually available - some sources of .TI3 don't provide it */
	if (!spec
	 && icg->find_field(icg, 0, "LAB_L") < 0
	 && icg->find_field(icg, 0, "XYZ_X") < 0) {

		if (icg->find_kword(icg, 0, "SPECTRAL_BANDS") < 0)
			error ("Neither CIE nor spectral data found in file '%s'",ti3name);

		/* Switch to using spectral information */
		if (verb)
			printf("No CIE data found, switching to spectral with standard observer & D50\n");
		spec = 1;
		illum = icxIT_D50;
		observ = icxOT_CIE_1931_2;
	}
	
	/* Figure out what sort of device it is */
	{
		int ti;

		if ((ti = icg->find_kword(icg, 0, "DEVICE_CLASS")) < 0)
			error ("Input file '%s' doesn't contain keyword DEVICE_CLASS",ti3name);

		if (strcmp(icg->t[0].kdata[ti],"DISPLAY") == 0) {
			isdisp = 1;
		}
	}

	/* Figure out what sort of device colorspace it is */
	{
		int ti;

		if ((ti = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
			error("Input file '%s' doesn't contain keyword COLOR_REPS",ti3name);

		if (strcmp(icg->t[0].kdata[ti],"CMYK_XYZ") == 0) {
			devspace = icSigCmykData;
			devchan = 4;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMYK_LAB") == 0) {
			devspace = icSigCmykData;
			devchan = 4;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMY_XYZ") == 0) {
			devspace = icSigCmyData;
			devchan = 3;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"CMY_LAB") == 0) {
			devspace = icSigCmyData;
			devchan = 3;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"RGB_XYZ") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"RGB_LAB") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 1;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"iRGB_XYZ") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"iRGB_LAB") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 1;
			isAdditive = 1;
		/* Scanner .ti3 files: */
		} else if (strcmp(icg->t[0].kdata[ti],"XYZ_RGB") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"LAB_RGB") == 0) {
			devspace = icSigRgbData;
			devchan = 3;
			isLab = 1;
			isAdditive = 1;
#ifdef IMP_MONO
		} else if (strcmp(icg->t[0].kdata[ti],"K_XYZ") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 0;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"K_LAB") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 1;
			isAdditive = 0;
		} else if (strcmp(icg->t[0].kdata[ti],"W_XYZ") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 0;
			isAdditive = 1;
		} else if (strcmp(icg->t[0].kdata[ti],"W_LAB") == 0) {
			devspace = icSigGrayData;
			devchan = 1;
			isLab = 1;
			isAdditive = 1;
#endif /* IMP_MONO */

		} else 
			error("Device input file '%s' has unhandled color representation '%s'",
			                                                     ti3name, icg->t[0].kdata[ti]);
	}

	if ((npat = icg->t[0].nsets) <= 0)
		error("Input file '%s' has no sets of data",ti3name);

	if (verb) {
		fprintf(verbo,"No of test patches = %d\n",npat);
	}

	/* Allocate arrays to hold test patch input and output values */
	if ((tpat = (pval *)malloc(sizeof(pval) * npat)) == NULL)
		error("Malloc failed - tpat[]");

	/* Read in the CGATs fields */
	{
		int sidx;					/* Sample ID index */
		int ti, ci, mi, yi, ki;
		int Xi, Yi, Zi;

		if ((sidx = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
			error("Input file '%s' doesn't contain field SAMPLE_ID",ti3name);
		if (icg->t[0].ftype[sidx] != nqcs_t)
			error("Input file '%s' field SAMPLE_ID is wrong type",ti3name);

		if (devspace == icSigGrayData) {
			if (isAdditive) {
				if ((ci = icg->find_field(icg, 0, "GRAY_W")) < 0)
					error("Input file doesn't contain field GRAY_W");
				if (icg->t[0].ftype[ci] != r_t)
					error("Field GRAY_W is wrong type - corrupted file ?");
			} else {
				if ((ci = icg->find_field(icg, 0, "GRAY_K")) < 0)
					error("Input file doesn't contain field GRAY_K");
				if (icg->t[0].ftype[ci] != r_t)
					error("Field GRAY_K is wrong type - corrupted file ?");
			}
			mi = yi = ki = ci;

		} else if (devspace == icSigRgbData) {
			if ((ci = icg->find_field(icg, 0, "RGB_R")) < 0)
				error("Input file '%s' doesn't contain field RGB_R",ti3name);
			if (icg->t[0].ftype[ci] != r_t)
				error("Input file '%s' field RGB_R is wrong type",ti3name);
			if ((mi = icg->find_field(icg, 0, "RGB_G")) < 0)
				error("Input file '%s' doesn't contain field RGB_G",ti3name);
			if (icg->t[0].ftype[mi] != r_t)
				error("Input file '%s' field RGB_G is wrong type",ti3name);
			if ((yi = icg->find_field(icg, 0, "RGB_B")) < 0)
				error("Input file '%s' doesn't contain field RGB_B",ti3name);
			if (icg->t[0].ftype[yi] != r_t)
				error("Input file '%s' field RGB_B is wrong type",ti3name);
			ki = yi;

		} else if (devspace == icSigCmyData) {

			if ((ci = icg->find_field(icg, 0, "CMY_C")) < 0)
				error("Input file '%s' doesn't contain field CMY_C",ti3name);
			if (icg->t[0].ftype[ci] != r_t)
				error("Input file '%s' field CMY_C is wrong type",ti3name);
			if ((mi = icg->find_field(icg, 0, "CMY_M")) < 0)
				error("Input file '%s' doesn't contain field CMY_M",ti3name);
			if (icg->t[0].ftype[mi] != r_t)
				error("Input file '%s' field CMY_M is wrong type",ti3name);
			if ((yi = icg->find_field(icg, 0, "CMY_Y")) < 0)
				error("Input file '%s' doesn't contain field CMY_Y",ti3name);
			ki = yi;
		} else {	/* Assume CMYK */

			if ((ci = icg->find_field(icg, 0, "CMYK_C")) < 0)
				error("Input file '%s' doesn't contain field CMYK_C",ti3name);
			if (icg->t[0].ftype[ci] != r_t)
				error("Input file '%s' field CMYK_C is wrong type",ti3name);
			if ((mi = icg->find_field(icg, 0, "CMYK_M")) < 0)
				error("Input file '%s' doesn't contain field CMYK_M",ti3name);
			if (icg->t[0].ftype[mi] != r_t)
				error("Input file '%s' field CMYK_M is wrong type",ti3name);
			if ((yi = icg->find_field(icg, 0, "CMYK_Y")) < 0)
				error("Input file '%s' doesn't contain field CMYK_Y",ti3name);
			if (icg->t[0].ftype[yi] != r_t)
				error("Input file '%s' field CMYK_Y is wrong type",ti3name);
			if ((ki = icg->find_field(icg, 0, "CMYK_K")) < 0)
				error("Input file '%s' doesn't contain field CMYK_K",ti3name);
			if (icg->t[0].ftype[ki] != r_t)
				error("Input file '%s' field CMYK_K is wrong type",ti3name);
		}

		if (spec == 0) { 		/* Using instrument tristimulous value */

			if (isLab) {		/* Expect Lab */
				if ((Xi = icg->find_field(icg, 0, "LAB_L")) < 0)
					error("Input file '%s' doesn't contain field LAB_L",ti3name);
				if (icg->t[0].ftype[Xi] != r_t)
					error("Input file '%s' field LAB_L is wrong type",ti3name);
				if ((Yi = icg->find_field(icg, 0, "LAB_A")) < 0)
					error("Input '%s' file doesn't contain field LAB_A",ti3name);
				if (icg->t[0].ftype[Yi] != r_t)
					error("Input file '%s' field LAB_A is wrong type",ti3name);
				if ((Zi = icg->find_field(icg, 0, "LAB_B")) < 0)
					error("Input file '%s' doesn't contain field LAB_B",ti3name);
				if (icg->t[0].ftype[Zi] != r_t)
					error("Input file '%s' field LAB_B is wrong type",ti3name);

			} else { 		/* Expect XYZ */
				if ((Xi = icg->find_field(icg, 0, "XYZ_X")) < 0)
					error("Input file '%s' doesn't contain field XYZ_X",ti3name);
				if (icg->t[0].ftype[Xi] != r_t)
					error("Input file '%s' field XYZ_X is wrong type",ti3name);
				if ((Yi = icg->find_field(icg, 0, "XYZ_Y")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Y",ti3name);
				if (icg->t[0].ftype[Yi] != r_t)
					error("Input file '%s' field XYZ_Y is wrong type",ti3name);
				if ((Zi = icg->find_field(icg, 0, "XYZ_Z")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Z",ti3name);
				if (icg->t[0].ftype[Zi] != r_t)
					error("Input file '%s' field XYZ_Z is wrong type",ti3name);
			}

			for (i = 0; i < npat; i++) {
				strcpy(tpat[i].sid, (char *)icg->t[0].fdata[i][sidx]);
				tpat[i].p[0] = *((double *)icg->t[0].fdata[i][ci]) / 100.0;
				tpat[i].p[1] = *((double *)icg->t[0].fdata[i][mi]) / 100.0;
				tpat[i].p[2] = *((double *)icg->t[0].fdata[i][yi]) / 100.0;
				tpat[i].p[3] = *((double *)icg->t[0].fdata[i][ki]) / 100.0;
				if (tpat[i].p[0] > 1.0
				 || tpat[i].p[1] > 1.0
				 || tpat[i].p[2] > 1.0
				 || tpat[i].p[3] > 1.0) {
					error("Input file '%s' device value field value exceeds 100.0 !",ti3name);
				}
				tpat[i].v[0] = *((double *)icg->t[0].fdata[i][Xi]);
				tpat[i].v[1] = *((double *)icg->t[0].fdata[i][Yi]);
				tpat[i].v[2] = *((double *)icg->t[0].fdata[i][Zi]);
				if (!isLab) {
					tpat[i].v[0] /= 100.0;		/* Normalise XYZ to range 0.0 - 1.0 */
					tpat[i].v[1] /= 100.0;
					tpat[i].v[2] /= 100.0;
				}
				if (!isLab) { /* Convert test patch result XYZ to PCS (D50 Lab) */
					icmXYZ2Lab(&icmD50, tpat[i].v, tpat[i].v);
				}
			}

		} else { 		/* Using spectral data */
			int ii;
			xspect sp;
			char buf[100];
			int  spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
			xsp2cie *sp2cie;	/* Spectral conversion object */

			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_BANDS")) < 0)
				error ("Input file '%s' doesn't contain keyword SPECTRAL_BANDS",ti3name);
			sp.spec_n = atoi(icg->t[0].kdata[ii]);
			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_START_NM")) < 0)
				error ("Input file '%s' doesn't contain keyword SPECTRAL_START_NM",ti3name);
			sp.spec_wl_short = atof(icg->t[0].kdata[ii]);
			if ((ii = icg->find_kword(icg, 0, "SPECTRAL_END_NM")) < 0)
				error ("Input file '%s; doesn't contain keyword SPECTRAL_END_NM",ti3name);
			sp.spec_wl_long = atof(icg->t[0].kdata[ii]);
			sp.norm = 100.0;

			/* Find the fields for spectral values */
			for (j = 0; j < sp.spec_n; j++) {
				int nm;
		
				/* Compute nearest integer wavelength */
				nm = (int)(sp.spec_wl_short + ((double)j/(sp.spec_n-1.0))
				            * (sp.spec_wl_long - sp.spec_wl_short) + 0.5);
				
				sprintf(buf,"SPEC_%03d",nm);

				if ((spi[j] = icg->find_field(icg, 0, buf)) < 0)
					error("Input file '%s' doesn't contain field %s",ti3name,buf);
			}

			if (isdisp) {
				illum = icxIT_none;		/* Displays are assumed to be self luminous */
			}

			/* Create a spectral conversion object */
			if ((sp2cie = new_xsp2cie(illum, illum == icxIT_none ? NULL : &cust_illum,
			                          observ, NULL, icSigLabData)) == NULL)
				error("Creation of spectral conversion object failed");

			if (fwacomp) {
				double nw = 0.0;		/* Number of media white patches */
				xspect mwsp;			/* Medium spectrum */
				instType itype;			/* Spectral instrument type */
				xspect insp;			/* Instrument illuminant */

				mwsp = sp;		/* Struct copy */

	 			if ((ti = icg->find_kword(icg, 0, "TARGET_INSTRUMENT")) < 0)
					error ("Input file '%s' can't find target instrument needed for FWA compensation",ti3name);

				if ((itype = inst_enum(icg->t[0].kdata[ti])) == instUnknown)
					error ("Input file '%s' unrecognised target instrument '%s'",ti3name, icg->t[0].kdata[ti]);

				if (inst_illuminant(&insp, itype) != 0)
					error ("Instrument doesn't have an FWA illuminent");

				/* Find the media white spectral reflectance */
				for (j = 0; j < mwsp.spec_n; j++)
					mwsp.spec[j] = 0.0;

				/* Compute the mean of all the media white patches */
				for (i = 0; i < npat; i++) {
					int use = 0;
	
					if (devspace == icSigGrayData) {
						if (isAdditive) {
							if (*((double *)icg->t[0].fdata[i][ci]) > (100.0 - 0.1))
								use = 1;
						} else {
							if (*((double *)icg->t[0].fdata[i][ci]) < 0.1)
								use = 1;
						}
					} else if (devspace == icSigRgbData) {
						if (*((double *)icg->t[0].fdata[i][ci]) > (100.0 - 0.1)
						 && *((double *)icg->t[0].fdata[i][mi]) > (100.0 - 0.1)
						 && *((double *)icg->t[0].fdata[i][yi]) > (100.0 - 0.1))
							use = 1;
					} else if (devspace == icSigCmyData) {
						if (*((double *)icg->t[0].fdata[i][ci]) < 0.1
						 && *((double *)icg->t[0].fdata[i][mi]) < 0.1
						 && *((double *)icg->t[0].fdata[i][yi]) < 0.1)
							use = 1;
					} else {	/* Assume CMYK */

						if (*((double *)icg->t[0].fdata[i][ci]) < 0.1
						 && *((double *)icg->t[0].fdata[i][mi]) < 0.1
						 && *((double *)icg->t[0].fdata[i][yi]) < 0.1
						 && *((double *)icg->t[0].fdata[i][ki]) < 0.1) {
							use = 1;
						}
					}

					if (use) {
						/* Read the spectral values for this patch */
						for (j = 0; j < mwsp.spec_n; j++) {
							mwsp.spec[j] += *((double *)icg->t[0].fdata[i][spi[j]]);
						}
						nw++;
					}
				}
				if (nw == 0.0) {
					warning("Input file '%s' can't find a media white patch to init FWA",ti3name);

					/* Track the maximum reflectance for any band to determine white. */
					/* This might give bogus results if there is no white patch... */
					for (i = 0; i < npat; i++) {
						for (j = 0; j < mwsp.spec_n; j++) {
							double rv = *((double *)icg->t[0].fdata[i][spi[j]]);
							if (rv > mwsp.spec[j])
								mwsp.spec[j] = rv;
						}
					}
					nw++;
				}

				for (j = 0; j < mwsp.spec_n; j++)
					mwsp.spec[j] /= nw;	/* Compute average */

				if (sp2cie->set_fwa(sp2cie, &insp, &mwsp)) 
					error ("Set FWA on sp2cie failed");

				if (verb) {
					double FWAc;
					sp2cie->get_fwa_info(sp2cie, &FWAc);
					fprintf(verbo,"FWA content = %f\n",FWAc);
				}
				
			}

			for (i = 0; i < npat; i++) {

				strcpy(tpat[i].sid, (char *)icg->t[0].fdata[i][sidx]);
				tpat[i].p[0] = *((double *)icg->t[0].fdata[i][ci]) / 100.0;
				tpat[i].p[1] = *((double *)icg->t[0].fdata[i][mi]) / 100.0;
				tpat[i].p[2] = *((double *)icg->t[0].fdata[i][yi]) / 100.0;
				tpat[i].p[3] = *((double *)icg->t[0].fdata[i][ki]) / 100.0;
				if (tpat[i].p[0] > 1.0
				 || tpat[i].p[1] > 1.0
				 || tpat[i].p[2] > 1.0
				 || tpat[i].p[3] > 1.0) {
					error("Input file '%s' device value field value exceeds 100.0 !",ti3name);
				}

				/* Read the spectral values for this patch */
				for (j = 0; j < sp.spec_n; j++) {
					sp.spec[j] = *((double *)icg->t[0].fdata[i][spi[j]]);
				}

				/* Convert it to CIE space */
				sp2cie->convert(sp2cie, tpat[i].v, &sp);
			}

			sp2cie->del(sp2cie);		/* Done with this */

		}

		icg->del(icg);		/* Clean up */
	}	/* End of reading in CGATs file */

	/* - - - - - - - - - - */
	/* Check the forward profile accuracy against the data points */
	{
		double merr = 0.0;		/* Max */
		double aerr = 0.0;		/* Avg */
		double rerr = 0.0;		/* RMS */
		double nsamps = 0.0;
		int inn, outn;			/* Chanells for input and output spaces */

		if (dovrml) {
			wrl = start_vrml(out_name, doaxes);
			start_line_set(wrl);
		}

		/* Open up the file for reading */
		if ((rd_fp = new_icmFileStd_name(iccname,"r")) == NULL)
			error("Write: Can't open file '%s'",iccname);

		if ((rd_icco = new_icc()) == NULL)
			error("Read: Creation of ICC object failed");

		/* Read the header and tag list */
		if ((rv = rd_icco->read(rd_icco,rd_fp,0)) != 0)
			error("Read: %d, %s",rv,rd_icco->err);

		/* Get the Fwd table, absolute with Lab override */
		if ((luo = rd_icco->get_luobj(rd_icco, icmFwd, icAbsoluteColorimetric,
		                              icSigLabData, icmLuOrdNorm)) == NULL) {
			error("%d, %s",rd_icco->errc, rd_icco->err);
		}

		/* Get details of conversion (Arguments may be NULL if info not needed) */
		luo->spaces(luo, NULL, &inn, NULL, &outn, NULL, NULL, NULL, NULL, NULL);

		for (i = 0; i < npat; i++) {
			double out[3];
			double mxd;

			/* Lookup the patch value in the profile */
			if (luo->lookup(luo, out, tpat[i].p) > 1)
				error("%d, %s",rd_icco->errc,rd_icco->err);

			if (verb > 1) {
				if (devspace == icSigCmykData) {
					printf("[%f] %s: %f %f %f %f -> %f %f %f should be %f %f %f\n",
					       cie2k ? icmCIE2K(tpat[i].v, out) : 
					               cie94 ? icmCIE94(tpat[i].v, out) : icmLabDE(tpat[i].v, out),
					       tpat[i].sid,
					       tpat[i].p[0],tpat[i].p[1],tpat[i].p[2],tpat[i].p[3],
					       out[0],out[1],out[2],
					       tpat[i].v[0],tpat[i].v[1],tpat[i].v[2]);
				} else {	/* Assume RGB/CMY */
					printf("[%f] %s: %f %f %f -> %f %f %f should be %f %f %f\n",
					       cie2k ? icmCIE2K(tpat[i].v, out) : 
					               cie94 ? icmCIE94(tpat[i].v, out) : icmLabDE(tpat[i].v, out),
					       tpat[i].sid,
					       tpat[i].p[0],tpat[i].p[1],tpat[i].p[2],
					       out[0],out[1],out[2],
					       tpat[i].v[0],tpat[i].v[1],tpat[i].v[2]);
				}
			}
			if (dovrml) {
				if (dominl && icmLabDE(tpat[i].v, out) < 0.5) {
					double cent[3], vec[3], vlen;
					double p1[3], p2[3];

					/* Compute center */
					icmAdd3(cent, tpat[i].v, out);
					icmScale3(cent, cent, 0.5);
					if ((vlen = icmLabDE(tpat[i].v, out)) < 1e-6) {
						vec[0] = 0.25; vec[1] = 0.0; vec[2] = 0.0;
					} else {
						icmSub3(vec, tpat[i].v, out);
						icmScale3(vec, vec, 0.25/vlen);
					}
					icmSub3(p1, cent, vec);
					icmAdd3(p2, cent, vec);
					add_vertex(wrl, p1);
					add_vertex(wrl, p2);
				} else {
					add_vertex(wrl, tpat[i].v);
					add_vertex(wrl, out);
				}
			}

			/* Check the result */
			if (cie2k)
				mxd = icmCIE2K(tpat[i].v, out);
			else if (cie94)
				mxd = icmCIE94(tpat[i].v, out);
			else
				mxd = icmLabDE(tpat[i].v, out);

			aerr += mxd;
			rerr += mxd * mxd;

			nsamps++;
			if (mxd > merr)
				merr = mxd;

		}
		if (dovrml) {
			if (dodecol)
				make_de_lines(wrl);
			else
				make_lines(wrl, 2);
			end_vrml(wrl);
		}
		printf("Profile check complete, errors%s: max. = %f, avg. = %f, RMS = %f\n",
            cie2k ? "(CIEDE2000)" : cie94 ? " (CIE94)" : "", merr, aerr/nsamps, sqrt(rerr/nsamps));

		/* ------------------------------- */
		/* If we want sort by target value */
		if (ddevv) {
			double cieval[3];

			/* Lookup the CIE value of the target */
			if (luo->lookup(luo, cieval, devval) > 1)
				error("%d, %s",rd_icco->errc,rd_icco->err);

			/* Compute deltas to target value. */
			for (i = 0; i < npat; i++) {
				if (cie2k)
					tpat[i].dv = icmCIE2K(tpat[i].v, cieval);
				else if (cie94)
					tpat[i].dv = icmCIE94(tpat[i].v, cieval);
				else
					tpat[i].dv = icmLabDE(tpat[i].v, cieval);

				tpat[i].dp = 0.0;
				for (j = 0; j < inn; j++) {
					double tt;
					tt = tpat[i].p[j] - devval[j];
					tpat[i].dp += tt * tt;
				}
				tpat[i].dp = sqrt(tpat[i].dp);
			}

			if (sortbypcs) {
				/* Sort by pcs delta */
#define HEAP_COMPARE(A,B) (A.dv < B.dv)
				HEAPSORT(pval, tpat, npat);
#undef HEAP_COMPARE
			} else {
				/* Sort by device delta */
#define HEAP_COMPARE(A,B) (A.dp < B.dp)
				HEAPSORT(pval, tpat, npat);
#undef HEAP_COMPARE
			}

			printf("Target point:\n");
			if (devspace == icSigCmykData) {
				printf("%f %f %f %f -> %f %f %f\n",
				       devval[0],devval[1],devval[2],devval[3],
				       cieval[0],cieval[1],cieval[2]);
			} else {	/* Assume RGB/CMY */
				printf("%f %f %f -> %f %f %f\n",
				       devval[0],devval[1],devval[2],
				       cieval[0],cieval[1],cieval[2]);
			}
			printf("\n");
	
			for (i = 0; i < npat; i++) {
				if (devspace == icSigCmykData) {
					printf("%s: %f %f %f %f [%f] -> %f %f %f [%f]\n",
					       tpat[i].sid,
					       tpat[i].p[0],tpat[i].p[1],tpat[i].p[2],tpat[i].p[3],
					       tpat[i].dp,
					       tpat[i].v[0],tpat[i].v[1],tpat[i].v[2],
					       tpat[i].dv);
				} else {	/* Assume RGB/CMY */
					printf("%s: %f %f %f [%f] -> %f %f %f [%f]\n",
					       tpat[i].sid,
					       tpat[i].p[0],tpat[i].p[1],tpat[i].p[2],
					       tpat[i].dp,
					       tpat[i].v[0],tpat[i].v[1],tpat[i].v[2],
					       tpat[i].dv);
				}
			}
		}

		/* Done with lookup object */
		luo->del(luo);

		/* Close the file */
		rd_icco->del(rd_icco);
		rd_fp->del(rd_fp);
	}

	return 0;
}


/* ------------------------------------------------ */
/* Some simple functions to do basix VRML work */
/* !!! Should change to plot/vrml lib !!! */

#define GAMUT_LCENT 50.0
static int npoints = 0;
static int paloc = 0;
static struct { double pp[3]; } *pary;

static void Lab2RGB(double *out, double *in);
static void DE2RGB(double *out, double in);

FILE *start_vrml(char *name, int doaxes) {
	FILE *wrl;

	/* Define the axis boxes */
	struct {
		double x, y, z;			/* Box center */
		double wx, wy, wz;		/* Box size */
		double r, g, b;			/* Box color */
	} axes[5] = {
		{ 0, 0,   50-GAMUT_LCENT, 2, 2, 100, .7, .7, .7 },	/* L axis */
		{ 50, 0,  0-GAMUT_LCENT,  100, 2, 2,  1,  0,  0 },	/* +a (red) axis */
		{ 0, -50, 0-GAMUT_LCENT,  2, 100, 2,  0,  0,  1 },	/* -b (blue) axis */
		{ -50, 0, 0-GAMUT_LCENT,  100, 2, 2,  0,  1,  0 },	/* -a (green) axis */
		{ 0,  50, 0-GAMUT_LCENT,  2, 100, 2,  1,  1,  0 },	/* +b (yellow) axis */
	};

	/* Define the labels */
	struct {
		double x, y, z;
		double size;
		char *string;
		double r, g, b;
	} labels[6] = {
		{ -2, 2, -GAMUT_LCENT + 100 + 10, 10, "+L*",  .7, .7, .7 },	/* Top of L axis */
		{ -2, 2, -GAMUT_LCENT - 10,      10, "0",    .7, .7, .7 },	/* Bottom of L axis */
		{ 100 + 5, -3,  0-GAMUT_LCENT,  10, "+a*",  1,  0,  0 },	/* +a (red) axis */
		{ -5, -100 - 10, 0-GAMUT_LCENT,  10, "-b*",  0,  0,  1 },	/* -b (blue) axis */
		{ -100 - 15, -3, 0-GAMUT_LCENT,  10, "-a*",  0,  0,  1 },	/* -a (green) axis */
		{ -5,  100 + 5, 0-GAMUT_LCENT,  10, "+b*",  1,  1,  0 },	/* +b (yellow) axis */
	};

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
		int n;
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
		fprintf(wrl,"    # Axes identification:\n");
		for (n = 0; n < 6; n++) {
			fprintf(wrl,"    Transform { translation %f %f %f\n", labels[n].x, labels[n].y, labels[n].z);
			fprintf(wrl,"      children [\n");
			fprintf(wrl,"        Shape{\n");
			fprintf(wrl,"          geometry Text { string [\"%s\"]\n",labels[n].string);
			fprintf(wrl,"            fontStyle FontStyle { family \"SANS\" style \"BOLD\" size %f }\n",
			                                  labels[n].size);
			fprintf(wrl,"                        }\n");
			fprintf(wrl,"          appearance Appearance { material Material ");
			fprintf(wrl,"{ diffuseColor %f %f %f} }\n", labels[n].r, labels[n].g, labels[n].b);
			fprintf(wrl,"        }\n");
			fprintf(wrl,"      ]\n");
			fprintf(wrl,"    }\n");
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

/* Assume 2 ppset, and make line color prop to length */
void make_de_lines(FILE *wrl) {
	int i, j;

	fprintf(wrl,"      ]\n");
	fprintf(wrl,"    }\n");
	fprintf(wrl,"  coordIndex [\n");

	for (i = 0; i < npoints;) {
		for (j = 0; j < 2; j++, i++) {
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
		double rgb[3], ss;
		for (ss = 0.0, j = 0; j < 3; j++) {
			double tt = (pary[i & ~1].pp[j] - pary[i | 1].pp[j]);
			ss += tt * tt;
		}
		ss = sqrt(ss);
		DE2RGB(rgb, ss);
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

/* Convert a delta E value into a signal color: */
static void
DE2RGB(double *out, double in) {
	struct {
		double de;
		double r, g, b;
	} range[6] = {
		{ 10.0, 1, 1, 0 },		/* yellow */
		{ 4.0,  1, 0, 0 },		/* red */
		{ 2.0, 1, 0, 1 },		/* magenta */
		{ 1.0, 0, 0, 1 },		/* blue */
		{ 0.5, 0, 1, 1 },		/* cyan */
		{ 0.0, 0, 1, 0 }		/* green */
	};
	int i;
	double bl;

	/* Locate the range we're in */
	if (in > range[0].de) {
		out[0] = range[0].r;
		out[1] = range[0].g;
		out[2] = range[0].b;
	} else {
		for (i = 0; i < 5; i++) {
			if (in <= range[i].de && in >= range[i+1].de)
				break;
		}
		bl = (in - range[i+1].de)/(range[i].de - range[i+1].de);
		out[0] = bl * range[i].r + (1.0 - bl) * range[i+1].r;
		out[1] = bl * range[i].g + (1.0 - bl) * range[i+1].g;
		out[2] = bl * range[i].b + (1.0 - bl) * range[i+1].b;
	}
}



