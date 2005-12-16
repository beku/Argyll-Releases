/* 
 * Argyll Color Correction System
 * Verify two sets of PCS values.
 *
 * Author: Graeme W. Gill
 * Date:   7/6/2005
 *
 * Copyright 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/*
 * This program takes in two CGATS files (probably .ti3 files) of PCS
 * values (either XYZ or L*a*b*), matches the values, and computes
 * overall errors. This is useful for verifying proofing systems.
 */

/*
 * TTBD:
 */

#define DEBUG

#define verbo stdout

#include <stdio.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include "copyright.h"
#include "config.h"
#include "numsup.h"
#include "cgats.h"
#include "xicc.h"
#include "insttypes.h"
#include "sort.h"

void
usage(void) {
	fprintf(stderr,"Verify CIE values, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	fprintf(stderr,"usage: verify [-options] target.ti3 measured.ti3\n");
	fprintf(stderr," -v              Verbose - print each patch value\n");
	fprintf(stderr," -c              Show CIE94 delta E values\n");
	fprintf(stderr," -k              Show CIEDE2000 delta E values\n");
	fprintf(stderr," -s              Sort patch values by error\n");
	fprintf(stderr," -w              create VRML vector visualisation (measured.wrl)\n");
	fprintf(stderr," -x              Use VRML axies\n");
	fprintf(stderr," -i illum        Choose illuminant for print/transparency spectral data:\n");
	fprintf(stderr,"                 A, D50 (def.), D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ       Choose CIE Observer for spectral data:\n");
	fprintf(stderr,"                 1931_2, 1964_10, S&B 1955_2, J&V 1978_2 (def.)\n");
	fprintf(stderr," -f              Use Fluorescent Whitening Agent compensation\n");
	fprintf(stderr," target.ti3      Target (reference) PCS or spectral values.\n");
	fprintf(stderr," measured.ti3    Measured (actual) PCS or spectral values\n");
	exit(1);
	}

/* VRML support */
FILE *start_vrml(char *name, int doaxes);
void start_line_set(FILE *wrl);
void add_vertex(FILE *wrl, double pp[3]);
void make_lines(FILE *wrl, int ppset);
void end_vrml(FILE *wrl);

/* Patch value type */
typedef struct {
	char sid[50];		/* sample id */
	double v[3];		/* CIE value */
	double de;			/* Delta E */
} pval;

main(int argc, char *argv[])
{
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	int cie94 = 0;
	int cie2k = 0;
	int dovrml = 0;
	int doaxes = 0;
	int dosort = 0;

	struct {
		char name[MAXNAMEL+1];	/* Patch filename  */
		int npat;				/* Number of patches */
		pval *pat;				/* patch values */
	} cg[2];					/* Target and current patch file information */

	int *match;					/* Array mapping first list indexes to corresponding second */
	int *sort;					/* Array of first list indexes in sorted order */
	int fwacomp = 0;			/* FWA compensation */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType illum = icxIT_D50;	/* Spectral defaults */
	xspect cust_illum;			/* Custom illumination spectrum */
	icxObserverType observ = icxOT_Judd_Voss_2;

	char out_name[MAXNAMEL+4+1]; /* VRML name */
	FILE *wrl;

	int i, j, n, rv = 0;

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

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

			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V')
				verb = 1;

			/* VRML */
			else if (argv[fa][1] == 'w' || argv[fa][1] == 'W')
				dovrml = 1;

			/* Axes */
			else if (argv[fa][1] == 'x' || argv[fa][1] == 'X')
				doaxes = 1;

			/* CIE94 delta E */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				cie94 = 1;
				cie2k = 0;
			}

			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				cie94 = 0;
				cie2k = 1;
			}

			/* Sort */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S')
				dosort = 1;

			/* Spectral Illuminant type */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				fa = nfa;
				if (na == NULL) usage();
				if (strcmp(na, "A") == 0) {
					spec = 1;
					illum = icxIT_A;
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

			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F') {
				spec = 1;
				fwacomp = 1;
			}

			else 
				usage();
		} else
			break;
	}

	/* Get the file name arguments */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(cg[0].name,argv[fa++],MAXNAMEL); cg[0].name[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(cg[1].name,argv[fa],MAXNAMEL); cg[1].name[MAXNAMEL] = '\000';

	/* Create VRML name */
	{
		char *xl;
		strcpy(out_name, cg[1].name);
		if ((xl = strrchr(out_name, '.')) == NULL)	/* Figure where extention is */
			xl = out_name + strlen(out_name);
		strcpy(xl,".wrl");
	}

	if (fwacomp && spec == 0)
		error ("FWA compensation only works when viewer and/or illuminant selected");

	/* Open up each file in turn, target then measured, */
	/* and read in the CIE values. */
	for (n = 0; n < 2; n++) {
		cgats *cgf = NULL;			/* cgats file data */
		int isLab = 0;				/* 0 if file CIE is XYZ, 1 if is Lab */
		int sidx;					/* Sample ID index */
		int xix, yix, zix;

		/* Open CIE target values */
		cgf = new_cgats();			/* Create a CGATS structure */
		cgf->add_other(cgf, ""); 	/* Allow any signature file */
	
		if (cgf->read_name(cgf, cg[n].name))
			error("CGATS file '%s' read error : %s",cg[n].name,cgf->err);
	
		if (cgf->ntables < 1)
			error ("Input file '%s' doesn't contain at least one table",cg[n].name);
	
		/* Check if the file is suitable */
		if (!spec
		 && cgf->find_field(cgf, 0, "LAB_L") < 0
		 && cgf->find_field(cgf, 0, "XYZ_X") < 0) {
	
			if (cgf->find_kword(cgf, 0, "SPECTRAL_BANDS") < 0)
				error ("Neither CIE nor spectral data found in file '%s'",cg[n].name);
	
			/* Switch to using spectral information */
			if (verb)
				printf("No CIE data found, switching to spectral with standard observer & D50 for file '%s'\n",cg[n].name);
			spec = 1;
			illum = icxIT_D50;
			observ = icxOT_CIE_1931_2;
		}
		if (spec && cgf->find_kword(cgf, 0, "SPECTRAL_BANDS") < 0)
			error ("No spectral data data found in file '%s' when spectral expected",cg[n].name);
	
		if (!spec && cgf->find_field(cgf, 0, "LAB_L") >= 0)
			isLab = 1;
		
		cg[n].npat = cgf->t[0].nsets;		/* Number of patches */
	
		/* Read all the target patches */
		if (cg[n].npat <= 0)
			error("No sets of data in file '%s'",cg[n].name);
	
		if (verb && n == 0) {
			fprintf(verbo,"No of test patches = %d\n",cg[n].npat);
		}
	
		/* Allocate arrays to hold test patch input and output values */
		if ((cg[n].pat = (pval *)malloc(sizeof(pval) * cg[n].npat)) == NULL)
			error("Malloc failed - pat[]");
	
		/* Read in the CGATs fields */
		if ((sidx = cgf->find_field(cgf, 0, "SampleName")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "Sample_Name")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "SAMPLE_NAME")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "SAMPLE_ID")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "SAMPLE_LOC")) < 0)
			error("Input file '%s' doesn't contain field SampleName, Sample_Name, SAMPLE_NAME, SAMPLE_ID or SAMPLE_LOC",cg[n].name);
		if (cgf->t[0].ftype[sidx] != nqcs_t)
			error("Sample ID field isn't non quoted character string");

		if (spec == 0) { 		/* Using instrument tristimulous value */

			if (isLab) {		/* Expect Lab */
				if ((xix = cgf->find_field(cgf, 0, "LAB_L")) < 0)
					error("Input file '%s' doesn't contain field LAB_L",cg[n].name);
				if (cgf->t[0].ftype[xix] != r_t)
					error("Field LAB_L is wrong type");
				if ((yix = cgf->find_field(cgf, 0, "LAB_A")) < 0)
					error("Input file '%s' doesn't contain field LAB_A",cg[n].name);
				if (cgf->t[0].ftype[yix] != r_t)
					error("Field LAB_A is wrong type");
				if ((zix = cgf->find_field(cgf, 0, "LAB_B")) < 0)
					error("Input file '%s' doesn't contain field LAB_B",cg[n].name);
				if (cgf->t[0].ftype[zix] != r_t)
					error("Field LAB_B is wrong type");

			} else { 		/* Expect XYZ */
				if ((xix = cgf->find_field(cgf, 0, "XYZ_X")) < 0)
					error("Input file '%s' doesn't contain field XYZ_X",cg[n].name);
				if (cgf->t[0].ftype[xix] != r_t)
					error("Field XYZ_X is wrong type");
				if ((yix = cgf->find_field(cgf, 0, "XYZ_Y")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Y",cg[n].name);
				if (cgf->t[0].ftype[yix] != r_t)
					error("Field XYZ_Y is wrong type");
				if ((zix = cgf->find_field(cgf, 0, "XYZ_Z")) < 0)
					error("Input file '%s' doesn't contain field XYZ_Z",cg[n].name);
				if (cgf->t[0].ftype[zix] != r_t)
					error("Field XYZ_Z is wrong type");
			}

			for (i = 0; i < cg[n].npat; i++) {
				strcpy(cg[n].pat[i].sid, (char *)cgf->t[0].fdata[i][sidx]);
				cg[n].pat[i].v[0] = *((double *)cgf->t[0].fdata[i][xix]);
				cg[n].pat[i].v[1] = *((double *)cgf->t[0].fdata[i][yix]);
				cg[n].pat[i].v[2] = *((double *)cgf->t[0].fdata[i][zix]);
				if (!isLab) {
					cg[n].pat[i].v[0] /= 100.0;		/* Normalise XYZ to range 0.0 - 1.0 */
					cg[n].pat[i].v[1] /= 100.0;
					cg[n].pat[i].v[2] /= 100.0;
				}
				if (!isLab) { /* Convert test patch result XYZ to PCS (D50 Lab) */
					icmXYZ2Lab(&icmD50, cg[n].pat[i].v, cg[n].pat[i].v);
				}
			}

		} else { 		/* Using spectral data */
			int ii;
			xspect sp;
			char buf[100];
			int  spi[XSPECT_MAX_BANDS];	/* CGATS indexes for each wavelength */
			xsp2cie *sp2cie;	/* Spectral conversion object */

			if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_BANDS")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_BANDS");
			sp.spec_n = atoi(cgf->t[0].kdata[ii]);
			if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_START_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_START_NM");
			sp.spec_wl_short = atof(cgf->t[0].kdata[ii]);
			if ((ii = cgf->find_kword(cgf, 0, "SPECTRAL_END_NM")) < 0)
				error ("Input file doesn't contain keyword SPECTRAL_END_NM");
			sp.spec_wl_long = atof(cgf->t[0].kdata[ii]);
			sp.norm = 100.0;

			/* Find the fields for spectral values */
			for (j = 0; j < sp.spec_n; j++) {
				int nm;
		
				/* Compute nearest integer wavelength */
				nm = (int)(sp.spec_wl_short + ((double)j/(sp.spec_n-1.0))
				            * (sp.spec_wl_long - sp.spec_wl_short) + 0.5);
				
				sprintf(buf,"SPEC_%03d",nm);

				if ((spi[j] = cgf->find_field(cgf, 0, buf)) < 0)
					error("Input file doesn't contain field %s",buf);
			}

			/* Figure out what sort of device it is */
			{
				int ti;
		
				if ((ti = cgf->find_kword(cgf, 0, "DEVICE_CLASS")) < 0)
					error ("Input file '%s' doesn't contain keyword DEVICE_CLASS",cg[n].name);
		
				if (strcmp(cgf->t[0].kdata[ti],"DISPLAY") == 0) {
					illum = icxIT_none;		/* Displays are assumed to be self luminous */
				}
			}

			/* Create a spectral conversion object */
			if ((sp2cie = new_xsp2cie(illum, illum == icxIT_none ? NULL : &cust_illum,
			                          observ, NULL, icSigLabData)) == NULL)
				error("Creation of spectral conversion object failed");

			if (fwacomp) {
				int ti;
				xspect mwsp;			/* Medium spectrum */
				instType itype;			/* Spectral instrument type */
				xspect *insp;			/* Instrument illuminant */

				mwsp = sp;		/* Struct copy */

	 			if ((ti = cgf->find_kword(cgf, 0, "TARGET_INSTRUMENT")) < 0)
					error ("Can't find target instrument in '%s' needed for FWA compensation",cg[n].name);

				if ((itype = inst_enum(cgf->t[0].kdata[ti])) == instUnknown)
					error ("Unrecognised target instrument '%s'", cgf->t[0].kdata[ti]);

				if ((insp = inst_illuminant(itype)) == NULL)
					error ("Instrument doesn't have an FWA illuminent");

				/* Determine a media white spectral reflectance */
				for (j = 0; j < mwsp.spec_n; j++)
					mwsp.spec[j] = 0.0;

				/* Track the maximum reflectance for any band to determine white. */
				/* This might silently fail, if there isn't white in the sampe set. */
				for (i = 0; i < cg[0].npat; i++) {
					for (j = 0; j < mwsp.spec_n; j++) {
						double rv = *((double *)cgf->t[0].fdata[i][spi[j]]);
						if (rv > mwsp.spec[j])
							mwsp.spec[j] = rv;
					}
				}
				if (sp2cie->set_fwa(sp2cie, insp, &mwsp)) 
					error ("Set FWA on sp2cie failed");
			}

			for (i = 0; i < cg[0].npat; i++) {

				strcpy(cg[n].pat[i].sid, (char *)cgf->t[0].fdata[i][sidx]);

				/* Read the spectral values for this patch */
				for (j = 0; j < sp.spec_n; j++) {
					sp.spec[j] = *((double *)cgf->t[0].fdata[i][spi[j]]);
				}

				/* Convert it to CIE space */
				sp2cie->convert(sp2cie, cg[n].pat[i].v, &sp);
			}

			sp2cie->del(sp2cie);		/* Done with this */

		}	/* End of reading in CGATs file */
		cgf->del(cgf);		/* Clean up */
	}

	/* Check that the number of test patches matches */
	if (cg[0].npat != cg[1].npat)
		error("Number of patches between '%s' and '%s' doesn't match",cg[0].name,cg[1].name);
	
	/* Create a list to map the second list of patches to the first */
	if ((match = (int *)malloc(sizeof(int) * cg[0].npat)) == NULL)
		error("Malloc failed - match[]");
	for (i = 0; i < cg[0].npat; i++) {
		for (j = 0; j < cg[1].npat; j++) {
			if (strcmp(cg[0].pat[i].sid, cg[1].pat[j].sid) == 0)
				break;			/* Found it */
		}
		if (j < cg[1].npat) {
			match[i] = j;
		} else {
			error("Failed to find matching patch to '%s'",cg[0].pat[i].sid);
		}
	}

	/* Compute the delta E's */
	for (i = 0; i < cg[0].npat; i++) {
		if (cie2k)
			cg[0].pat[i].de = icmCIE2K(cg[0].pat[i].v, cg[1].pat[match[i]].v);
		else if (cie94)
			cg[0].pat[i].de = icmCIE94(cg[0].pat[i].v, cg[1].pat[match[i]].v);
		else
			cg[0].pat[i].de = icmLabDE(cg[0].pat[i].v, cg[1].pat[match[i]].v);
	}

	/* Create sorted list, from worst to best. */
	if ((sort = (int *)malloc(sizeof(int) * cg[0].npat)) == NULL)
		error("Malloc failed - sort[]");
	for (i = 0; i < cg[0].npat; i++) 
		sort[i] = i;

#define HEAP_COMPARE(A,B) (cg[0].pat[A].de > cg[0].pat[B].de)
	HEAPSORT(int, sort, cg[0].npat);
#undef HEAP_COMPARE

	/* - - - - - - - - - - */
	/* Figure out the report */
	{
		double merr = 0.0, aerr = 0.0;
		int n90;
		double merr90 = 0.0, aerr90 = 0.0;
		int n10;
		double merr10 = 0.0, aerr10 = 0.0;

		if (dovrml) {
			wrl = start_vrml(out_name, doaxes);
			start_line_set(wrl);
		}

		/* Do overall results */
		for (i = 0; i < cg[0].npat; i++) {
			double de;
			if (dosort)
				j = sort[i];
			else
				j = i;

			de = cg[0].pat[j].de;
			aerr += de;

			if (verb) {
				printf("%s: %f %f %f <=> %f %f %f  de %f\n",
					cg[0].pat[j].sid,
					cg[0].pat[j].v[0], cg[0].pat[j].v[1], cg[0].pat[j].v[2],
					cg[1].pat[match[j]].v[0], cg[1].pat[match[j]].v[1], cg[1].pat[match[j]].v[2],
					de);
			}

			if (de > merr)
				merr = de;

			if (dovrml && de > 1e-6) {
				add_vertex(wrl, cg[0].pat[j].v);
				add_vertex(wrl, cg[1].pat[j].v);
			}

		}
		if (cg[0].npat > 0)
			aerr /= (double)cg[0].npat;

		if (dovrml) {
			make_lines(wrl, 2);
			end_vrml(wrl);
		}

		/* Do best 90% */
		n90 = (int)(cg[0].npat * 9.0/10.0 + 0.5);
		for (i = (cg[0].npat-n90); i < cg[0].npat; i++) {
			double de = cg[0].pat[sort[i]].de;
			aerr90 += de;
			if (de > merr90)
				merr90 = de;
		}
		if (n90 > 0)
			aerr90 /= (double)n90;

		/* Do worst 10% */
		n10 = (int)(cg[0].npat * 1.0/10.0 + 0.5);
		for (i = 0; i < n10; i++) {
			double de = cg[0].pat[sort[i]].de;
			aerr10 += de;
			if (de > merr10)
				merr10 = de;
		}
		if (n10 > 0)
			aerr10 /= (double)n10;

		if (verb) {
			fprintf(verbo,"No of test patches in worst 10%% are = %d\n",n10);
			fprintf(verbo,"No of test patches in best 90%% are = %d\n",n90);
		}
		printf("Verify results:\n");
		printf("  Total errors%s:     peak = %f, avg = %f\n", cie2k ? " (CIEDE2000)" : cie94 ? " (CIE94)" : "", merr, aerr);
		printf("  Worst 10%% errors%s: peak = %f, avg = %f\n", cie2k ? " (CIEDE2000)" : cie94 ? " (CIE94)" : "", merr10, aerr10);
		printf("  Best  90%% errors%s: peak = %f, avg = %f\n", cie2k ? " (CIEDE2000)" : cie94 ? " (CIE94)" : "", merr90, aerr90);
	
		free(sort);
		free(match);
		free(cg[0].pat);
		free(cg[1].pat);
	}

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




