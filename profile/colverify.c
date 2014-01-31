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
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes in two CGATS files (probably but not necesserily .ti3 files) of PCS
 * values (either XYZ, L*a*b* or spectral), matches the values, and computes
 * overall errors. This is useful for verifying proofing systems.
 */

/*
 * TTBD:
 */

#undef DEBUG

#define verbo stdout

#include <stdio.h>
#include <string.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include "copyright.h"
#include "aconfig.h"
#include "numlib.h"
#include "vrml.h"
#include "cgats.h"
#include "xicc.h"
#include "ccmx.h"
#include "insttypes.h"
#include "sort.h"

void
usage(void) {
	fprintf(stderr,"Verify CIE values, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the AGPL Version 3\n");
	fprintf(stderr,"usage: verify [-options] target.ti3 measured.ti3\n");
	fprintf(stderr," -v [n]          Verbose mode, n >= 2 print each value\n");
	fprintf(stderr," -n              Normalise each files reading to its white Y\n");
	fprintf(stderr," -N              Normalise each files reading to its white XYZ\n");
	fprintf(stderr," -m              Normalise each files reading to its white X+Y+Z\n");
	fprintf(stderr," -D              Use D50 100.0 as L*a*b* white reference\n");
	fprintf(stderr," -c              Show CIE94 delta E values\n");
	fprintf(stderr," -k              Show CIEDE2000 delta E values\n");
	fprintf(stderr," -s              Sort patch values by error\n");
	fprintf(stderr," -w              create VRML vector visualisation (measured.wrl)\n");
	fprintf(stderr," -W              create VRML marker visualisation (measured.wrl)\n");
	fprintf(stderr," -x              Use VRML axes\n");
	fprintf(stderr," -f [illum]      Use Fluorescent Whitening Agent compensation [opt. simulated inst. illum.:\n");
	fprintf(stderr,"                  M0, M1, M2, A, C, D50 (def.), D50M2, D65, F5, F8, F10 or file.sp]\n");
	fprintf(stderr," -i illum        Choose illuminant for computation of CIE XYZ from spectral data & FWA:\n");
	fprintf(stderr,"                  A, C, D50 (def.), D50M2, D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ       Choose CIE Observer for spectral data:\n");
	fprintf(stderr,"                 1931_2 (def), 1964_10, S&B 1955_2, shaw, J&V 1978_2\n");
	fprintf(stderr," -L profile.%s  Skip any first file out of profile gamut patches\n",ICC_FILE_EXT_ND);
	fprintf(stderr," -X file.ccmx    Apply Colorimeter Correction Matrix to second file\n");
	fprintf(stderr," target.ti3      Target (reference) PCS or spectral values.\n");
	fprintf(stderr," measured.ti3    Measured (actual) PCS or spectral values\n");
	exit(1);
	}

/* Patch value type */
typedef struct {
	char sid[50];		/* sample id */
	char loc[100];		/* sample location (empty if none) */
	int og;				/* Out of gamut flag */
	double xyz[3];		/* XYZ value */
	double v[3];		/* Lab value */
	double de;			/* Delta E */
	double ixde[3];		/* XYZ Component DE */
	double ide[3];		/* Lab Component DE */
} pval;

int main(int argc, char *argv[])
{
	int fa,nfa,mfa;			/* current argument we're looking at */
	int verb = 0;       	/* Verbose level */
	int norm = 0;			/* 1 = norm to White Y, 2 = norm to White XYZ */
							/* 3 = norm to White X+Y+Z */
	int usestdd50 = 0;		/* Use standard D50 instead of avg white as reference */
	int cie94 = 0;
	int cie2k = 0;
	int dovrml = 0;
	int doaxes = 0;
	int dosort = 0;
	char ccmxname[MAXNAMEL+1] = "\000";  /* Colorimeter Correction Matrix name */
	ccmx *cmx = NULL;					/* Colorimeter Correction Matrix */
	char gprofname[MAXNAMEL+1] = "\000";  /* Gamut limit profile name */
	icmFile *fp = NULL;
	icc *icco = NULL;
	xicc *xicco = NULL;
	icxLuBase *luo = NULL;

	struct {
		char name[MAXNAMEL+1];	/* Patch filename  */
		int isdisp;				/* nz if display */
		int isdnormed;      	/* Has display data been normalised to 100 ? */
		int npat;				/* Number of patches */
		int nig;				/* Number of patches in gamut */
		double w[3];			/* XYZ of "white" */
		double nw[3];			/* Normalised XYZ of "white" */
		pval *pat;				/* patch values */
	} cg[2];					/* Target and current patch file information */

	int *match;					/* Array mapping first list indexes to corresponding second */
	int *sort;					/* Array of first list indexes in sorted order */
	int fwacomp = 0;			/* FWA compensation */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType tillum = icxIT_none;	/* Target/simulated instrument illuminant */ 
	xspect cust_tillum, *tillump = NULL; /* Custom target/simulated illumination spectrum */
	icxIllumeType illum = icxIT_D50;	/* Spectral defaults */
	xspect cust_illum;			/* Custom illumination spectrum */
	icxObserverType observ = icxOT_CIE_1931_2;

	icmXYZNumber labw = icmD50;	/* The Lab white reference */

	char out_name[MAXNAMEL+4+1]; /* VRML name */
	vrml *wrl = NULL;

	int i, j, n;

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	if (argc <= 1)
		usage();

	/* Process the arguments */
	mfa = 2;        				/* Minimum final arguments */
	for(fa = 1;fa < argc;fa++) {

		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1+mfa) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage();

			/* Verbose */
			else if (argv[fa][1] == 'v') {
				verb = 1;

				if (na != NULL && na[0] >= '0' && na[0] <= '9') {
					verb = atoi(na);
					fa = nfa;
				}
			}

			/* normalize */
			else if (argv[fa][1] == 'n'
			      || argv[fa][1] == 'N') {
				norm = 1;
				if (argv[fa][1] == 'N')
					norm = 2;
			}

			else if (argv[fa][1] == 'm') {
				norm = 3;
			}

			else if (argv[fa][1] == 'D')
				usestdd50 = 1;

			/* VRML */
			else if (argv[fa][1] == 'w')
				dovrml = 1;
			else if (argv[fa][1] == 'W')
				dovrml = 2;

			/* Axes */
			else if (argv[fa][1] == 'x')
				doaxes = 1;

			/* CIE94 delta E */
			else if (argv[fa][1] == 'c') {
				cie94 = 1;
				cie2k = 0;
			}

			else if (argv[fa][1] == 'k') {
				cie94 = 0;
				cie2k = 1;
			}

			/* Sort */
			else if (argv[fa][1] == 's')
				dosort = 1;

			/* FWA compensation */
			else if (argv[fa][1] == 'f') {
				fwacomp = 1;

				if (na != NULL) {	/* Argument is present - target/simulated instr. illum. */
					fa = nfa;
					if (strcmp(na, "A") == 0
					 || strcmp(na, "M0") == 0) {
						spec = 1;
						tillum = icxIT_A;
					} else if (strcmp(na, "C") == 0) {
						spec = 1;
						tillum = icxIT_C;
					} else if (strcmp(na, "D50") == 0
					        || strcmp(na, "M1") == 0) {
						spec = 1;
						tillum = icxIT_D50;
					} else if (strcmp(na, "D50M2") == 0
					        || strcmp(na, "M2") == 0) {
						spec = 1;
						tillum = icxIT_D50M2;
					} else if (strcmp(na, "D65") == 0) {
						spec = 1;
						tillum = icxIT_D65;
					} else if (strcmp(na, "F5") == 0) {
						spec = 1;
						tillum = icxIT_F5;
					} else if (strcmp(na, "F8") == 0) {
						spec = 1;
						tillum = icxIT_F8;
					} else if (strcmp(na, "F10") == 0) {
						spec = 1;
						tillum = icxIT_F10;
					} else {	/* Assume it's a filename */
						spec = 1;
						tillum = icxIT_custom;
						if (read_xspect(&cust_tillum, na) != 0)
							usage();
					}
				}
			}

			/* Spectral to CIE Illuminant type */
			else if (argv[fa][1] == 'i') {
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
				} else if (strcmp(na, "D50M2") == 0) {
					spec = 1;
					illum = icxIT_D50M2;
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
			else if (argv[fa][1] == 'o') {
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

			/* Gamut limit profile for first file */
			else if (argv[fa][1] == 'L') {
				fa = nfa;
				if (na == NULL) usage();
				strncpy(gprofname,na,MAXNAMEL-1); gprofname[MAXNAMEL-1] = '\000';
			}

			/* Colorimeter Correction Matrix for second file */
			else if (argv[fa][1] == 'X') {
				fa = nfa;
				if (na == NULL) usage();
				strncpy(ccmxname,na,MAXNAMEL-1); ccmxname[MAXNAMEL-1] = '\000';

			} else 
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

	/* Gamut limit profile */
	if (gprofname[0] != '\000') {
		int rv;

		if ((fp = new_icmFileStd_name(gprofname,"r")) == NULL)
			error ("Can't open file '%s'",gprofname);
	
		if ((icco = new_icc()) == NULL)
			error ("Creation of ICC object failed");
	
		if ((rv = icco->read(icco,fp,0)) != 0)
			error("Reading profile '%s' failed failed with error %d:'%s'\n",
		     	       gprofname, icco->errc,  icco->err);

		if (icco->header->deviceClass != icSigInputClass
		 && icco->header->deviceClass != icSigDisplayClass
		 && icco->header->deviceClass != icSigOutputClass)
			error("Profile '%s' must be a device profile to filter by gamut",gprofname);

		/* Wrap with an expanded icc */
		if ((xicco = new_xicc(icco)) == NULL)
			error ("Creation of xicc failed");

		/* Get a expanded color conversion object */
		if ((luo = xicco->get_luobj(xicco, ICX_CLIP_NEAREST | ICX_FAST_SETUP,
		    icmFwd, icRelativeColorimetric, icSigXYZData, icmLuOrdNorm, NULL, NULL)) == NULL)
			error ("%d, %s",xicco->errc, xicco->err);
	}

	/* Colorimeter Correction Matrix */
	if (ccmxname[0] != '\000') {
		if ((cmx = new_ccmx()) == NULL)
			error("new_ccmx failed\n");
		if (cmx->read_ccmx(cmx,ccmxname))
			error("Reading Colorimeter Correction Matrix file '%s' failed with error %d:'%s'\n",
		     	       ccmxname, cmx->errc,  cmx->err);
	}


	/* Open up each file in turn, target then measured, */
	/* and read in the CIE values. */
	for (n = 0; n < 2; n++) {
		cgats *cgf = NULL;			/* cgats file data */
		int isLab = 0;				/* 0 if file CIE is XYZ, 1 if is Lab */
		int sidx;					/* Sample ID index */
		int sldx = -1;				/* Sample location index, < 0 if invalid */
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
		
		cg[n].nig = cg[n].npat = cgf->t[0].nsets;		/* Number of patches */
	
		/* Figure out what sort of device it is */
		{
			int ti;
	
			cg[n].isdisp = 0;
			cg[n].isdnormed = 0;
			cg[n].w[0] = cg[n].w[1] = cg[n].w[2] = 0.0;

			if ((ti = cgf->find_kword(cgf, 0, "DEVICE_CLASS")) < 0)
				error ("Input file '%s' doesn't contain keyword DEVICE_CLASS",cg[n].name);
	
			if (strcmp(cgf->t[0].kdata[ti],"DISPLAY") == 0) {
				cg[n].isdisp = 1;
				cg[n].isdnormed = 1;	/* Assume display type is normalised to 100 */
				illum = icxIT_none;		/* Displays are assumed to be self luminous */
				/* ?? What if two files are different ?? */
			}

			if (cg[n].isdisp) {

				if ((ti = cgf->find_kword(cgf, 0, "LUMINANCE_XYZ_CDM2")) >= 0) {
					if (sscanf(cgf->t[0].kdata[ti], " %lf %lf %lf ",&cg[n].w[0], &cg[n].w[1], &cg[n].w[2]) != 3)
						cg[n].w[0] = cg[n].w[1] = cg[n].w[2] = 0.0;
				}

				/* See if there is an explicit tag indicating data has been normalised to Y = 100 */
				if ((ti = cgf->find_kword(cgf, 0, "NORMALIZED_TO_Y_100")) >= 0) {
					if (strcmp(cgf->t[0].kdata[ti],"NO") == 0) {
						cg[n].isdnormed = 0;
					} else {
						cg[n].isdnormed = 1;
					}
				}
			}
		}

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
		if ((sidx = cgf->find_field(cgf, 0, "SAMPLE_ID")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "SampleName")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "Sample_Name")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "SAMPLE_NAME")) < 0
		 && (sidx = cgf->find_field(cgf, 0, "SAMPLE_LOC")) < 0)
			error("Input file '%s' doesn't contain field SAMPLE_ID, SampleName, Sample_Name, SAMPLE_NAME or SAMPLE_LOC",cg[n].name);
		if (cgf->t[0].ftype[sidx] != nqcs_t
		 && cgf->t[0].ftype[sidx] != cs_t)
			error("Sample ID/Name field isn't a quoted or non quoted character string");

		if ((sldx = cgf->find_field(cgf, 0, "SAMPLE_LOC")) < 0
		 || cgf->t[0].ftype[sldx] != cs_t)
			sldx = -1;

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
				if (sldx >= 0)
					strcpy(cg[n].pat[i].loc, (char *)cgf->t[0].fdata[i][sldx]);
				else
					cg[n].pat[i].loc[0] = '\000';
				cg[n].pat[i].og = 0;
				cg[n].pat[i].xyz[0] = *((double *)cgf->t[0].fdata[i][xix]);
				cg[n].pat[i].xyz[1] = *((double *)cgf->t[0].fdata[i][yix]);
				cg[n].pat[i].xyz[2] = *((double *)cgf->t[0].fdata[i][zix]);

				if (isLab) {	/* Convert to XYZ */
					icmLab2XYZ(&icmD50, cg[n].pat[i].xyz, cg[n].pat[i].xyz);
				}
//printf("~1 file %d patch %d = XYZ %f %f %f\n", n,i,cg[n].pat[i].xyz[0],cg[n].pat[i].xyz[1],cg[n].pat[i].xyz[2]);

				/* restore normalised display values to absolute */
				if (cg[n].isdnormed) {
					if (cg[n].w[1] > 0.0) {
						cg[n].pat[i].xyz[0] *= cg[n].w[1]/100.0;
						cg[n].pat[i].xyz[1] *= cg[n].w[1]/100.0;
						cg[n].pat[i].xyz[2] *= cg[n].w[1]/100.0;
					}

				} else if (!cg[n].isdisp) {
					/* If reflective or transmissive that are 0..100%, */
					/* scale back to 0.. 1 */
					cg[n].pat[i].xyz[0] /= 100.0;		/* scale back to XYZ 1.0 */
					cg[n].pat[i].xyz[1] /= 100.0;
					cg[n].pat[i].xyz[2] /= 100.0;
				}

				/* Apply ccmx */
				if (n == 1 && cmx != NULL) {
					cmx->xform(cmx, cg[n].pat[i].xyz, cg[n].pat[i].xyz);
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
			if (!cg[n].isdisp || cg[n].isdnormed != 0)
				sp.norm = 100.0;
			else
				sp.norm = 1.0;

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

			/* Create a spectral conversion object */
			if ((sp2cie = new_xsp2cie(illum, illum == icxIT_none ? NULL : &cust_illum,
			                          observ, NULL, icSigXYZData, icxClamp)) == NULL)
				error("Creation of spectral conversion object failed");

			if (fwacomp) {
				int ti;
				xspect mwsp;			/* Medium spectrum */
				instType itype;			/* Spectral instrument type */
				xspect insp;			/* Instrument illuminant */

				mwsp = sp;		/* Struct copy */

	 			if ((ti = cgf->find_kword(cgf, 0, "TARGET_INSTRUMENT")) < 0)
					error ("Can't find target instrument in '%s' needed for FWA compensation",cg[n].name);

				if ((itype = inst_enum(cgf->t[0].kdata[ti])) == instUnknown)
					error ("Unrecognised target instrument '%s'", cgf->t[0].kdata[ti]);

				if (inst_illuminant(&insp, itype) != 0)
					error ("Instrument doesn't have an FWA illuminent");

				/* Determine a media white spectral reflectance */
				for (j = 0; j < mwsp.spec_n; j++)
					mwsp.spec[j] = 0.0;

				/* Since we don't want to assume that there are any associated device */
				/* values present in each file, we can't use this as means of */
				/* determining the media color. Use an alternative approach here, */
				/* which may give slightly different results to profile. */

				/* Track the maximum reflectance for any band to determine white. */
				/* This might silently fail, if there isn't white in the sample set. */
				for (i = 0; i < cg[0].npat; i++) {
					for (j = 0; j < mwsp.spec_n; j++) {
						double rv = *((double *)cgf->t[0].fdata[i][spi[j]]);
						if (rv > mwsp.spec[j])
							mwsp.spec[j] = rv;
					}
				}

				/* If we are setting a specific simulated instrument illuminant */
				if (tillum != icxIT_none) {
					tillump = &cust_tillum;
					if (tillum != icxIT_custom) {
						if (standardIlluminant(tillump, tillum, 0.0)) {
							error("simulated inst. illum. not recognised");
						}
					}
				}

				if (sp2cie->set_fwa(sp2cie, &insp, tillump, &mwsp)) 
					error ("Set FWA on sp2cie failed");

				if (verb) {
					double FWAc;
					sp2cie->get_fwa_info(sp2cie, &FWAc);
					fprintf(verbo,"FWA content = %f\n",FWAc);
				}
			}

			for (i = 0; i < cg[0].npat; i++) {

				strcpy(cg[n].pat[i].sid, (char *)cgf->t[0].fdata[i][sidx]);
				if (sldx >= 0)
					strcpy(cg[n].pat[i].loc, (char *)cgf->t[0].fdata[i][sldx]);
				else
					cg[n].pat[i].loc[0] = '\000';
				cg[n].pat[i].og = 0;

				/* Read the spectral values for this patch */
				for (j = 0; j < sp.spec_n; j++) {
					sp.spec[j] = *((double *)cgf->t[0].fdata[i][spi[j]]);
				}

				/* Convert it to XYZ space */
				sp2cie->convert(sp2cie, cg[n].pat[i].xyz, &sp);

				/* restore normalised display values to absolute */
				if (cg[n].isdnormed) {
					if (cg[n].w[1] > 0.0) {
						cg[n].pat[i].xyz[0] *= cg[n].w[1];
						cg[n].pat[i].xyz[1] *= cg[n].w[1];
						cg[n].pat[i].xyz[2] *= cg[n].w[1];
					}

				}

				/* Apply ccmx */
				if (n == 1 && cmx != NULL) {
					cmx->xform(cmx, cg[n].pat[i].xyz, cg[n].pat[i].xyz);
				}
			}

			sp2cie->del(sp2cie);		/* Done with this */

		}	/* End of reading in CGATs file */


		/* Locate the patch with maximum Y, a possible white patch */
		if (norm) {
			int ii;

			if (cg[n].w[1] == 0.0) {	/* No white patch */

				/* Locate patch with biggest Y, assume it is white... */
				for (i = 0; i < cg[n].npat; i++) {
					if (cg[n].pat[i].xyz[1] > cg[n].w[1]) {
						icmCpy3(cg[n].w, cg[n].pat[i].xyz);
						ii = i;
					}
				}
				if (verb) printf("File %d Chose patch %d as white, xyz %f %f %f\n",
				                       n, ii+1,cg[n].w[0],cg[n].w[1],cg[n].w[2]);
			} else {
				if (verb) printf("File %d White is from display luminance ref. xyz %f %f %f\n",
				                       n, cg[n].w[0],cg[n].w[1],cg[n].w[2]);
			}
			icmCpy3(cg[n].nw, cg[n].w);
		}


		/* Normalise this file to white = 1.0 or D50 */
		if (norm) {
			int ii;

			double chmat[3][3];				/* Chromatic adapation matrix */

			if (norm == 2) {		/* Norm to white XYZ */ 
				icmXYZNumber s_wp;
				icmAry2XYZ(s_wp, cg[n].w);
				icmChromAdaptMatrix(ICM_CAM_BRADFORD, icmD50, s_wp, chmat);
			}

			for (i = 0; i < cg[n].npat; i++) {
				if (norm == 1) {
					cg[n].pat[i].xyz[0] *= 100.0 / cg[n].w[1];
					cg[n].pat[i].xyz[1] *= 100.0 / cg[n].w[1];
					cg[n].pat[i].xyz[2] *= 100.0 / cg[n].w[1];
				} else if (norm == 2) { 
					icmMulBy3x3(cg[n].pat[i].xyz, chmat, cg[n].pat[i].xyz);
				} else {
					cg[n].pat[i].xyz[0] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
					cg[n].pat[i].xyz[1] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
					cg[n].pat[i].xyz[2] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
				}
//printf("~1 file %d patch %d = norm XYZ %f %f %f\n", n,i,cg[n].pat[i].xyz[0],cg[n].pat[i].xyz[1],cg[n].pat[i].xyz[2]);
			}
			/* Compute normalised white too */
			if (norm == 1) {
				cg[n].nw[0] *= 100.0 / cg[n].w[1];
				cg[n].nw[1] *= 100.0 / cg[n].w[1];
				cg[n].nw[2] *= 100.0 / cg[n].w[1];
			} else if (norm == 2) { 
				icmMulBy3x3(cg[n].nw, chmat, cg[n].w);
			} else {
				cg[n].nw[0] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
				cg[n].nw[1] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
				cg[n].nw[2] *= 100.0 / (cg[n].w[0] + cg[n].w[1] + cg[n].w[2]);
			}
//printf("~1 file %d norm white XYZ %f %f %f\n", n,cg[n].nw[0], cg[n].nw[1], cg[n].nw[2]);
		}
		cgf->del(cgf);		/* Clean up */
	}
	if (cmx != NULL) 
		cmx->del(cmx);
	cmx = NULL;

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

	/* Figure out which patches to skip because they are out of gamut */
	if (luo != NULL) {
		double chmat[3][3];				/* Chromatic adapation matrix */
		double out[MAX_CHAN], in[3], check[3];
		icmXYZNumber s_wp;
		int rv;

		/* Convert sample PCS to relative */
		icmAry2XYZ(s_wp, cg[0].nw);
		icmChromAdaptMatrix(ICM_CAM_BRADFORD, icmD50, s_wp, chmat);

		for (i = 0; i < cg[0].npat; i++) {
			icmMulBy3x3(in, chmat, cg[0].pat[i].xyz);
//printf("~1 %d: xyz %f %f %f, rel %f %f %f\n", i+1, cg[0].pat[i].xyz[0], cg[0].pat[i].xyz[1], cg[0].pat[i].xyz[2], in[0], in[1], in[2]);

			if ((rv = luo->inv_lookup(luo, out, in)) > 0 || 1) {
				double de;

				luo->lookup(luo, check, out);
				de = icmXYZLabDE(&icmD50,check, in);
//printf("~1 %d: rv %d, de %f, check XYZ %f %f %f\n",i+1,rv, de, check[0],check[1],check[2]);

				if (de >= 0.01) {
					cg[0].pat[i].og = 1;
//printf("~1 Patch %d is out of gamut by DE %f RGB %f %f %f\n",i+1,de,out[0],out[1],out[2]);
					if (verb >= 3)
						printf("Patch %d is out of gamut by DE %f\n",i+1,de);
					cg[0].nig--;
				}
			}
		}
		if (verb)
			fprintf(verbo,"No of test patches in gamut = %d/%d\n",cg[0].npat - cg[0].nig,cg[0].npat);
	}

	/* Adjust the Lab reference white to be the mean of the white of the two files */
	if (norm != 0 && !usestdd50) {
		labw.X = 0.5 * (cg[0].nw[0] + cg[1].nw[0]);
		labw.Y = 0.5 * (cg[0].nw[1] + cg[1].nw[1]);
		labw.Z = 0.5 * (cg[0].nw[2] + cg[1].nw[2]);

		if (verb)
			printf("L*a*b* white reference = XYZ %f %f %f\n",labw.X,labw.Y,labw.Z);
	}

	/* Convert XYZ to Lab */
	for (n = 0; n < 2; n++) {
		for (i = 0; i < cg[n].npat; i++) {
			icmXYZ2Lab(&labw, cg[n].pat[i].v, cg[n].pat[i].xyz);
		}
	}

	/* Compute the delta E's */
	for (i = 0; i < cg[0].npat; i++) {

		if (cg[0].pat[i].og)		/* Skip out of gamut patches */
			continue;

		cg[0].pat[i].ixde[0] = fabs(cg[0].pat[i].xyz[0] - cg[1].pat[match[i]].xyz[0]);
		cg[0].pat[i].ixde[1] = fabs(cg[0].pat[i].xyz[1] - cg[1].pat[match[i]].xyz[1]);
		cg[0].pat[i].ixde[2] = fabs(cg[0].pat[i].xyz[2] - cg[1].pat[match[i]].xyz[2]);

		if (cie2k)
			cg[0].pat[i].de = icmCIE2K(cg[0].pat[i].v, cg[1].pat[match[i]].v);
		else if (cie94)
			cg[0].pat[i].de = icmCIE94(cg[0].pat[i].v, cg[1].pat[match[i]].v);
		else
			cg[0].pat[i].de = icmLabDE(cg[0].pat[i].v, cg[1].pat[match[i]].v);

		cg[0].pat[i].ide[0] = fabs(cg[0].pat[i].v[0] - cg[1].pat[match[i]].v[0]);
		cg[0].pat[i].ide[1] = fabs(cg[0].pat[i].v[1] - cg[1].pat[match[i]].v[1]);
		cg[0].pat[i].ide[2] = fabs(cg[0].pat[i].v[2] - cg[1].pat[match[i]].v[2]);
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
		double rad;
		double aierr[3] = { 0.0, 0.0, 0.0 };
		double aixerr[3] = { 0.0, 0.0, 0.0 };

		if (dovrml) {
			wrl = new_vrml(out_name, doaxes, 0);
			wrl->start_line_set(wrl, 0);

			/* Fudge sphere diameter */
			rad = 10.0/pow(cg[0].npat, 1.0/3.0);
		}

		/* Do overall results */
		for (i = 0; i < cg[0].npat; i++) {
			double de;

			if (cg[0].pat[i].og)		/* Skip out of gamut patches */
				continue;

			if (dosort)
				j = sort[i];
			else
				j = i;

			de = cg[0].pat[j].de;
			aerr += de;

			aierr[0] += cg[0].pat[j].ide[0];
			aierr[1] += cg[0].pat[j].ide[1];
			aierr[2] += cg[0].pat[j].ide[2];

			aixerr[0] += cg[0].pat[j].ixde[0];
			aixerr[1] += cg[0].pat[j].ixde[1];
			aixerr[2] += cg[0].pat[j].ixde[2];

			if (verb >= 2) {

				printf("%s%s%s: %f %f %f <=> %f %f %f  de %f\n",
					cg[0].pat[j].sid,
					cg[0].pat[j].loc[0] != '\000' ? " " : "",
					cg[0].pat[j].loc,
					cg[0].pat[j].v[0], cg[0].pat[j].v[1], cg[0].pat[j].v[2],
					cg[1].pat[match[j]].v[0], cg[1].pat[match[j]].v[1], cg[1].pat[match[j]].v[2],
					de);

#ifdef NEVER	/* Print XYZ as well */
				printf("  %f %f %f <=> %f %f %f\n",
					cg[0].pat[j].xyz[0], cg[0].pat[j].xyz[1], cg[0].pat[j].xyz[2],
					cg[1].pat[match[j]].xyz[0], cg[1].pat[match[j]].xyz[1], cg[1].pat[match[j]].xyz[2]);
#endif
			}

			if (de > merr)
				merr = de;

			if (dovrml) {
				if (de > 1e-6) {
					wrl->add_vertex(wrl, 0, cg[0].pat[j].v);
					wrl->add_vertex(wrl, 0, cg[1].pat[j].v);
				}
				if (dovrml == 2) {
					wrl->add_marker(wrl, cg[0].pat[j].v, NULL, rad);
					wrl->add_marker(wrl, cg[1].pat[j].v, NULL, rad);
				}
			}

		}
		if (cg[0].nig > 0) {
			aerr /= (double)cg[0].nig;
			aierr[0] /= (double)cg[0].nig;
			aierr[1] /= (double)cg[0].nig;
			aierr[2] /= (double)cg[0].nig;

			aixerr[0] /= (double)cg[0].nig;
			aixerr[1] /= (double)cg[0].nig;
			aixerr[2] /= (double)cg[0].nig;
		}

		if (dovrml) {
			wrl->make_lines(wrl, 0, 2);
			wrl->del(wrl);
			wrl = NULL;
		}

		/* Do best 90% */
		n90 = (int)(cg[0].nig * 9.0/10.0 + 0.5);
		for (i = j = 0; i < cg[0].npat; i++) {
			double de = cg[0].pat[sort[i]].de;
			if (cg[0].pat[i].og)		/* Skip out of gamut */
				continue;
			if (j >= (cg[0].nig-n90)) {	/* If in top 90% of in gamut patches */
				aerr90 += de;
				if (de > merr90)
					merr90 = de;
			}
			j++;						/* Index of within gamut patches */
		}
		if (n90 > 0)
			aerr90 /= (double)n90;

		/* Do worst 10% */
		n10 = (int)(cg[0].nig * 1.0/10.0 + 0.5);
		for (i = j = 0; i < cg[0].npat; i++) {
			double de = cg[0].pat[sort[i]].de;
			if (cg[0].pat[i].og)		/* Skip out of gamut */
				continue;
			if (j <= n10) {				/* If in worst 10% of in gamut patches */
				aerr10 += de;
				if (de > merr10)
					merr10 = de;
			}
			j++;
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
		printf("  avg err X  %f, Y  %f, Z  %f\n", aixerr[0], aixerr[1], aixerr[2]);
		printf("  avg err L* %f, a* %f, b* %f\n", aierr[0], aierr[1], aierr[2]);
	
		free(sort);
		free(match);
		free(cg[0].pat);
		free(cg[1].pat);
	}

	if (luo != NULL)
		luo->del(luo);
	if (xicco != NULL)
		xicco->del(xicco);		/* Expansion wrapper */
	if (icco != NULL)
		icco->del(icco);		/* Icc */
	if (fp != NULL)
		fp->del(fp);

	return 0;
}


