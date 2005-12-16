
/* 
 * Create/modify an abstract ICC transformation that will
 * correct measured color inacuracies of a proof reproduction.
 * Input is a list of target and actual CIE values.
 *
 * Author:  Graeme W. Gill
 * Date:    12/5/05
 * Version: 1.00
 *
 * Copyright 2005 Graeme W. Gill
 * Please refer to Licence.txt file for details.
 */

/* TTBD:
 *
 */

/* Basic idea:

	Given to .ti3 files (or equivalent CGATS files), one containing a
	spread of target patch values (XYZ, Lab or spectral), and the
	other containing the corresponding measured values, a PCS->PCS
	correction RSPL mapping is created to adjust for any innacuracy
	in the B2A table, which is then used to refine an existing abstract
	correction profile. (Any device values in the .ti3 tables are ignored.)

	The refined abstract profile can the be used to create an adjusted
	device profile, or an adjusted device link.

	Complications are:

		Measurements are absolute, while an abstract profile
		is relative.
		- use flag to mark this when abstract is used, or mark in profile header ?

		Corrections shouldn't go outside the target devices gamut, or
		they will lead to out of control regions on the gamut surface.
		- need output profile to clip changes to gamut surface.

		Refinement feedback could go unstable.
		- use a damping factor to improve stability.

	Interestingly, for CMYK the results are most stable (in simulation)
	when applied to the simple linked device link, and tends to
	be slightly unstable when applied to inverse A2B lookups
	that are used in profile and icclink -G. This could be a symtom
	of the black generation non-uniformity problem causing instability
	in the black inversion. Moving to optimised separation CMYK
	profile generation might overcome this problem.

	NOTE:- the current value for the rspl weak default weight seems OK
	for a reasonable number of points, but if refine was to be used for
	arbitrary tweaking, it should probably be made a tunable parameter
	that affects the radius of influence of each adjustment point.
*/


#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include "copyright.h"
#include "config.h"
#include "rspl.h"
#include "xicc.h"
#include "xicc.h"

#undef DEBUG		/* Print each value changed */

#define verbo stdout

#define DEF_DAMP1 0.95
#define DEF_DAMP2 0.50
#define DEF_CLUTRES 33
#define GAMRES 10.0

void usage(void) {
	fprintf(stderr,"Create abstract correction profile given table of absolute CIE correction values\n");
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	fprintf(stderr,"usage: refine [-options] cietarget ciecurrent [outdevicc] [inabs] outabs\n");
	fprintf(stderr," -v              Verbose\n");
	fprintf(stderr," -c              Create initial abstract correction profile\n");
	fprintf(stderr," -g              Don't impose output device gamut limit\n");
	fprintf(stderr," -r res          Set abstract profile clut resolution (default %d)\n",DEF_CLUTRES);
	fprintf(stderr," -d factor       Override default damping factor (default %f, then %f)\n",DEF_DAMP1,DEF_DAMP2);
	fprintf(stderr," -i illum        Choose illuminant for spectral data:\n");
	fprintf(stderr,"                 A, D50 (def.), D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ       Choose CIE Observer for spectral data:\n");
	fprintf(stderr,"                 1931_2, 1964_10, S&B 1955_2, J&V 1978_2 (def.)\n");
	fprintf(stderr," -f              Use Fluorescent Whitening Agent compensation on spectral data\n");
	fprintf(stderr," cietarget       Target CIE or spectral values, CGATS file (e.g. .ti3)\n");
	fprintf(stderr," ciecurrent      Actual CIE or spectral values, CGATS file (e.g. .ti3)\n");
	fprintf(stderr," [outdevicc]     Output device ICC profile to set gamut limit (not used if -g)\n");
	fprintf(stderr," [inabs]         Previous abstract correction ICC profile (not used if -c)\n");
	fprintf(stderr," outabs          Created/refined abstract correction ICC profile\n");
	exit(1);
}

/* ------------------------------------------- */
/* structure to support icc Lut initialisation/modification calbacks */

struct _callback {
	int verb;				/* Verbosity */
	int total, count, last;	/* Progress count information */
	rspl *r;				/* correction transform */
	icmLuBase *rd_luo;		/* Existing abstract profile (NULL if none) */
	gamut *dev_gam;			/* Gamut of output device (NULL if none) */
}; typedef struct _callback callback;


/* - - - - */
/*  clut  */

/* New CLUT table */
/* Correct for PCS errors */
void PCSp_PCSp(void *cntx, double *out, double *in) {
	double vv[3], temp[3];
	callback *p = (callback *)cntx;
	co pp;

#ifdef DEBUG
	printf("Got Lab in %f %f %f\n",in[0],in[1],in[2]);
#endif

//printf("~1 radial on %f %f %f = %f\n",in[0],in[1],in[2],p->dev_gam->nradial(p->dev_gam, temp, in));
	vv[0] = in[0];
	vv[1] = in[1];
	vv[2] = in[2];

	/* We compound new correction only if */
	/* this is the first one, 
	/* the source is in gamut, to prevent runaway corrections. */
	if (p->rd_luo == NULL
	 || p->dev_gam == NULL
	 || p->dev_gam->nradial(p->dev_gam, temp, in) <= 1.00) {
		pp.p[0] = vv[0];
		pp.p[1] = vv[1];
		pp.p[2] = vv[2];
		p->r->interp(p->r, &pp);				/* This correction */
		vv[0] = pp.v[0];
		vv[1] = pp.v[1];
		vv[2] = pp.v[2];
	}

	/* Compound with previous correction */
	if (p->rd_luo != NULL) {
		p->rd_luo->lookup(p->rd_luo, vv, vv);			/* Previous correction */
	}

	out[0] = vv[0];
	out[1] = vv[1];
	out[2] = vv[2];

#ifdef DEBUG
	printf("Got Lab out %f %f %f\n",out[0],out[1],out[2]);
	printf("\n");
#endif

	if (p->verb) {		/* Output percent intervals */
		int pc;
		p->count++;
		pc = p->count * 100.0/p->total + 0.5;
		if (pc != p->last) {
			printf("\r%2d%%",pc), fflush(stdout);
			p->last = pc;
		}
	}
}

/* ------------------------------------------- */

/* Patch value type */
typedef struct {
	char sid[50];		/* sample id */
	double v[3];		/* CIE value */
	double de;			/* Delta E */
} pval;

/* Weak default function */
static void wfunc(void *cbntx, double *out, double *in) {
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];
}

int
main(int argc, char *argv[]) {
	int fa,nfa;				/* argument we're looking at */
	struct {
		char name[MAXNAMEL+1];	/* Patch filename  */
		int npat;				/* Number of patches */
		pval *pat;				/* patch values */
	} cg[2];					/* Target and current patch file information */
	char dev_name[MAXNAMEL+1];	/* Output device ICC filename for gamut */
	char rd_name[MAXNAMEL+1];	/* Abstract profile ICC to modify */
	char wr_name[MAXNAMEL+1];	/* Modified/created abstract profile ICC */

	int *match;					/* Array mapping first list indexes to corresponding second */
	int fwacomp = 0;			/* FWA compensation on spectral ? */
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType illum = icxIT_D50;	/* Spectral defaults */
	xspect cust_illum;					/* Custom illumination spectrum */
	icxObserverType observ = icxOT_Judd_Voss_2;
	callback cb;				/* Callback support stucture for setting abstract profile */

	icmFile *dev_fp;			/* Output device profile for gamut limit */
	icc *dev_icc;
	xicc *dev_xicc;

	icmFile *rd_fp;				/* Existing abstract profile to modify */
	icc *rd_icc = NULL;

	icmFile *wr_fp;				/* Modified/created abstract profile to write */
	icc *wr_icc;

	int verb = 0;
	int nogamut = 0;					/* Don't impose a gamut limit */
	int docreate = 0;					/* Create an initial abstract correction profile */
	int clutres = DEF_CLUTRES;			/* Output abstract profile clut resolution */
	double damp = 0.0;					/* Damping factor */
	double smoothf = 1.0;				/* RSPL Smoothing factor */
	double avgdev = 0.005;				/* RSPL Average Deviation */
	double wweight = 1.0;				/* weak default function weight */
	double merr = 0.0, aerr = 0.0;		/* Stats on color change */
	int i, j, e, n, rv = 0;

	error_program = argv[0];

	if (argc < 6)
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
			/* Create initial abstract correction profile */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				docreate = 1;
			}
			/* Don't impose a gamut limit */
			else if (argv[fa][1] == 'g' || argv[fa][1] == 'G') {
				nogamut = 1;
			}
			/* Override the correction clut resolution */
			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R') {
				fa = nfa;
				if (na == NULL) usage();
				clutres = atoi(na);
			}
			/* Override the damping factor */
			else if (argv[fa][1] == 'd' || argv[fa][1] == 'D') {
				fa = nfa;
				if (na == NULL) usage();
				damp = atof(na);
			}


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

			/* FWA compensation */
			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F')
				fwacomp = 1;

			else 
				usage();
		} else
			break;
	}

	/* Grab all the filenames: */

	/* The two CIE value files */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(cg[0].name,argv[fa++],MAXNAMEL); cg[0].name[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(cg[1].name,argv[fa++],MAXNAMEL); cg[1].name[MAXNAMEL] = '\000';

	/* Optional output device name */
	if (nogamut == 0) {
		if (fa >= argc || argv[fa][0] == '-') usage();
		strncpy(dev_name,argv[fa++],MAXNAMEL); dev_name[MAXNAMEL] = '\000';
	}

	/* Optional input abstract profile name */
	if (docreate == 0) {
		if (fa >= argc || argv[fa][0] == '-') usage();
		strncpy(rd_name,argv[fa++],MAXNAMEL); rd_name[MAXNAMEL] = '\000';
	}

	/* Output abstract profile name */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(wr_name,argv[fa++],MAXNAMEL); wr_name[MAXNAMEL] = '\000';

	if (damp == 0.0) {		/* Use default damping */
		if (docreate)
			damp = DEF_DAMP1;
		else
			damp = DEF_DAMP2;
	}

	/* ======================= */
	/* Open up each CIE file in turn, target then measured, */
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

	/* Compute the delta E's just for information */
	for (i = 0; i < cg[0].npat; i++) {
		double de = icmLabDE(cg[0].pat[i].v, cg[1].pat[match[i]].v);
		cg[0].pat[i].de = de;
		if (de > merr)
			merr = de;
		aerr += de;
	}
	if (cg[0].npat > 0)
		aerr /= (double)cg[0].npat;

	if (verb) {
		fprintf(verbo,"No of correction patches = %d\n",cg[0].npat);
		fprintf(verbo,"Average dE = %f, Maximum dE = %f\n",aerr,merr);
	}

	/* ======================= */
	/* Create refining rspl */
	{
		cow *rp;		/* rspl setup points */
		int npnts = 0;	/* Total number of test points */
		int gres[MXDI];	/* rspl grid resolution */
		datai mn, mx;

		if ((rp = (cow *)malloc(sizeof(cow) * cg[0].npat)) == NULL)
			error("Malloc failed - rp[]");
		
		/* Create mapping points */
		for (i = 0; i < cg[0].npat; i++) {

			/* Input is target [0] */
			for (j = 0; j < 3; j++)
				rp[i].p[j] = cg[0].pat[i].v[j];

			/* Cull out of range points */
			if (rp[i].p[0] < 0.0 || rp[i].p[0] > 100.0
			 || rp[i].p[1] < -127.0 || rp[i].p[1] > 127.0
			 || rp[i].p[2] < -127.0 || rp[i].p[2] > 127.0)
				continue;
			
			/* Creat output as absolute Lab correction */
			/* [0] = target, [1] = measured */
			for (j = 0; j < 3; j++) {
				rp[i].v[j] = cg[0].pat[i].v[j] + damp * (cg[0].pat[i].v[j] - cg[1].pat[match[i]].v[j]);
			}

			/* Set weighting */
			rp[i].w = 1.0;
			npnts++;
		}

		/* Create refining rspl */
		mn[0] =   0.0, mn[1] = mn[2] = -128.0;			/* Allow for 16 bit grid range */
		mx[0] = 100.0, mx[1] = mx[2] =  (65535.0 * 255.0)/65280.0 - 128.0;
		cb.verb = verb;
		if ((cb.r = new_rspl(3, 3)) == NULL)
			error("new_rspl failed");

		for (e = 0; e < 3; e++)
			gres[e] = clutres;

		cb.r->fit_rspl_w_df(cb.r,
		           RSPL_EXTRAFIT		/* Extra fit flag */
		           | verb ? RSPL_VERBOSE : 0,
		           rp,					/* Test points */
		           npnts,				/* Number of test points */
		           mn, mx, gres,		/* Low, high, resolution of grid */
		           NULL, NULL,			/* Default data scale */
		           smoothf,				/* Smoothing */
		           avgdev,				/* Average Deviation */
                   wweight,				/* weak default function weight */
				   NULL,				/* No context */
		           wfunc				/* Weak function */
		);
		if (verb) printf("\n");
		free(rp);
	}

	/* ======================= */
	/* Possible limiting gamut */
	if (nogamut == 0) {
		icmFile *dev_fp;
		icc *dev_icc;
		xicc *dev_xicc;
		icxLuBase *dev_luo;
		icxInk ink;							/* Ink parameters */

		/* Open up the device ICC profile, so that we can create a gamut */
		/* and get an absolute PCS->device conversion */
		if ((dev_fp = new_icmFileStd_name(dev_name,"r")) == NULL)
			error ("Can't open file '%s'",dev_name);
	
		if ((dev_icc = new_icc()) == NULL)
			error ("Creation of ICC object failed");
	
		/* Read header etc. */
		if ((rv = dev_icc->read(dev_icc,dev_fp,0)) != 0)
			error ("%d, %s",rv,dev_icc->err);
	
		/* Check that the profile is appropriate */
		if (dev_icc->header->deviceClass != icSigInputClass
		 && dev_icc->header->deviceClass != icSigDisplayClass
		 && dev_icc->header->deviceClass != icSigOutputClass)
			error("Device Profile '%s' isn't a device profile",dev_name);

		ink.tlimit = -1.0;		/* No ink limit by default */
		ink.klimit = -1.0;

		/* Use a heuristic to guess the ink limit */
		if (icmCSSig2nchan(dev_icc->header->colorSpace) > 3) {
			int ttres = 17;
			int co[3];
			double in[3], out[MAX_CHAN];
			icmLuBase *luo;
			int nooch;
			double maxink = 0.0;

			nooch = icmCSSig2nchan(dev_icc->header->colorSpace);
			/* Lab to device lookup */
			if ((luo = dev_icc->get_luobj(dev_icc, icmBwd, icRelativeColorimetric, icSigLabData, icmLuOrdNorm)) == NULL)
				error ("%d, %s",dev_icc->errc, dev_icc->err);

			for (co[0] = 0; co[0] < ttres; co[0]++) {
				for (co[1] = 0; co[1] < ttres; co[1]++) {
					for (co[2] = 0; co[2] < ttres; co[2]++) {
						double sum;

						in[0] = 60.0 * co[0]/(ttres-1.0);
						in[1] = 200.0 * co[1]/(ttres-1.0) - 100.0;
						in[2] = 200.0 * co[2]/(ttres-1.0) - 100.0;

						/* PCS -> Device */
						if ((rv = luo->lookup(luo, out, in)) > 1)
							error ("%d, %s",rd_icc->errc,rd_icc->err);
		
						for (i = 0, sum = 0.0; i < nooch; i++)
							sum += out[i];
						if (sum > maxink)
							maxink = sum;
					}
				}
			}
			maxink += 0.05;		/* allow a slight margine */
			if (verb)
				printf("Guessed inklimit is %f%%\n",100.0 * maxink);
			ink.tlimit = maxink;
			luo->del(luo);
		}

		/* Wrap with an expanded icc */
		if ((dev_xicc = new_xicc(dev_icc)) == NULL)
			error ("Creation of xicc failed");

		/* Get a expanded color conversion object suitable for gamut */
		if ((dev_luo = dev_xicc->get_luobj(dev_xicc, 0, icmFwd, icAbsoluteColorimetric, icSigLabData, icmLuOrdNorm, NULL, &ink)) == NULL)
			error ("%d, %s",dev_xicc->errc, dev_xicc->err);
	
		/* Creat a gamut surface */
		if ((cb.dev_gam = dev_luo->get_gamut(dev_luo, GAMRES)) == NULL)
			error ("%d, %s",dev_xicc->errc, dev_xicc->err);

		dev_luo->del(dev_luo);
		dev_xicc->del(dev_xicc);
		dev_icc->del(dev_icc);
		dev_fp->del(dev_fp);
	} else {
		cb.dev_gam = NULL;
	}

	/* ======================= */
	/* Open up the existing abstract profile that is to be refined. */
	if (docreate == 0) {
		if ((rd_fp = new_icmFileStd_name(rd_name,"r")) == NULL)
			error ("Can't open file '%s'",rd_name);
	
		if ((rd_icc = new_icc()) == NULL)
			error ("Creation of ICC object failed");
	
		/* Read header etc. */
		if ((rv = rd_icc->read(rd_icc,rd_fp,0)) != 0)
			error ("%d, %s",rv,rd_icc->err);
	
		if (rd_icc->header->deviceClass != icSigAbstractClass)
			error("Input Profile '%s' isn't abstract type",rd_name);

		if ((cb.rd_luo = rd_icc->get_luobj(rd_icc, icmFwd, icAbsoluteColorimetric, icSigLabData, icmLuOrdNorm)) == NULL)
				error ("%d, %s",rd_icc->errc, rd_icc->err);
	} else {
		cb.rd_luo = NULL;
	}

	/* ======================= */
	/* Create new abstract ICC profile */
	if ((wr_fp = new_icmFileStd_name(wr_name,"w")) == NULL)
		error ("Can't open file '%s' for writing",wr_name);

	if ((wr_icc = new_icc()) == NULL)
		error ("Creation of write ICC object failed");

	/* Add all the tags required */

	/* The header: */
	{
		icmHeader *wh = wr_icc->header;

		/* Values that must be set before writing */
		wh->deviceClass     = icSigAbstractClass;
    	wh->colorSpace      = icSigLabData;
    	wh->pcs             = icSigLabData;
    	wh->renderingIntent = icAbsoluteColorimetric;	/* Instrument reading based */
	}
	/* Profile Description Tag: */
	{
		icmTextDescription *wo;
		char *dst = "Argyll refine output";
		if ((wo = (icmTextDescription *)wr_icc->add_tag(
		           wr_icc, icSigProfileDescriptionTag,	icSigTextDescriptionType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

		wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, dst);		/* Copy the string in */
	}
	/* Copyright Tag: */
	{
		icmText *wo;
		char *crt = "Copyright the user who created it";
		if ((wo = (icmText *)wr_icc->add_tag(
		           wr_icc, icSigCopyrightTag,	icSigTextType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

		wo->size = strlen(crt)+1; 	/* Allocated and used size of text, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->data, crt);		/* Copy the text in */
	}
	/* White Point Tag: */
	{
		icmXYZArray *wo;
		/* Note that tag types icSigXYZType and icSigXYZArrayType are identical */
		if ((wo = (icmXYZArray *)wr_icc->add_tag(
		           wr_icc, icSigMediaWhitePointTag, icSigXYZArrayType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

		wo->size = 1;
		wo->allocate((icmBase *)wo);	/* Allocate space */
		wo->data[0] = icmD50;			/* So absolute/relative rendering is the same */
	}
	/* 16 bit pcs -> pcs lut: */
	{
		icmLut *wo;

		/* Intent 0 = default/perceptual */
		if ((wo = (icmLut *)wr_icc->add_tag(
		           wr_icc, icSigAToB0Tag,	icSigLut16Type)) == NULL) 
			error("add_tag failed: %d, %s",wr_icc->errc,wr_icc->err);

		wo->inputChan = 3;
		wo->outputChan = 3;
    	wo->clutPoints = clutres;
    	wo->inputEnt = 256;				/* Not actually used */
    	wo->outputEnt = 256;
		wo->allocate((icmBase *)wo);/* Allocate space */

		/* The matrix is only applicable to XYZ input space, */
		/* so it is not used here. */

		/* Use helper function to do the hard work. */
		if (cb.verb) {
			for (cb.total = 1, i = 0; i < 3; i++, cb.total *= wo->clutPoints)
				; 
			cb.count = 0;
			cb.last = -1;
			printf(" 0%%"), fflush(stdout);
		}
		if (wo->set_tables(wo, &cb,
				icSigLabData, 			/* Input color space */
				icSigLabData, 			/* Output color space */
				NULL,					/* Linear input transform Lab->Lab' (NULL = default) */
				NULL, NULL,				/* Use default Maximum range of Lab' values */
				PCSp_PCSp,				/* Lab' -> Lab' transfer function */
				NULL, NULL,				/* Use default Maximum range of Lab' values */
				NULL					/* Linear output transform Lab'->Lab */
		) != 0)
			error("Setting 16 bit Lab->Lab Lut failed: %d, %s",wr_icc->errc,wr_icc->err);

		if (verb)
			printf("\nDone filling abstract table\n");
	}
	/* Write the file out */
	if ((rv = wr_icc->write(wr_icc,wr_fp,0)) != 0)
		error ("Write file: %d, %s",rv,wr_icc->err);
	
	/* ======================================= */
	
	/* Clean everything up */
	wr_icc->del(wr_icc);
	wr_fp->del(wr_fp);

	if (docreate == 0) {
		cb.rd_luo->del(cb.rd_luo);
		rd_icc->del(rd_icc);
		rd_fp->del(rd_fp);
	}

	if (nogamut == 0) {
		cb.dev_gam->del(cb.dev_gam);
	}

	cb.r->del(cb.r);

	free(match);
	free(cg[0].pat);
	free(cg[1].pat);

	return 0;
}

