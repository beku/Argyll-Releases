/* 
 * Argyll Color Correction System
 * film target chart reader
 *
 * Author: Neil Okamoto
 * Date:   3/1/2001
 *
 * Copyright 2001, DreamWorks LLC
 * Copyright 2006 - 2007, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE :-
 * see the License2.txt file for licencing details.
 */

/* This program reads a film target chart using a colorimeter. */

#undef DEBUG

#define COMPORT 1		/* Default com port 1..4 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "xspect.h"
#include "cgats.h"
#include "insttypes.h"
#include "icoms.h"
#include "../h/sort.h"
#include "inst.h"

#include <stdarg.h>

#if defined (NT)
#include <conio.h>
#endif

/* A color structure */
/* This can hold all representations simultaneously */
typedef struct {
	double gy;
	double r,g,b;
	double cmyk[4];
	double XYZ[4];		/* Colorimeter readings */
	double eXYZ[4];		/* Expected XYZ values */
	double spect[36];	/* Raw spectral data */
	char *id;			/* Id string */
	char *loc;			/* Location string */
	int loci;			/* Location integer = pass * 256 + step */
} col;

int
read_patches(col *cols, int npat, int comport, int verb, int qspect)
{
	inst *gt;
	int i,j, rv;

	gt = new_inst(comport, instSpectroScanT, 0, 0);

	/* Establish communications */
	if ((rv = gt->init_coms(gt, comport, baud_nc, -1.0)) != inst_ok) {
		error("couldn't init SpectroScan/T on port %d",comport);
	}

	/* Initialise the instrument */
	if ((rv = gt->init_inst(gt)) != inst_ok) {
		error("failed to init SpectroScan/T");
	}

	/* Configure the instrument mode */
	{
		inst_capability cap;
		inst_mode mode = 0;

		cap = gt->capabilities(gt);

		if ((cap & inst_trans_spot) == 0) {
			printf("Need spot transmission measurement capability,\n");
			printf("but instrument doesn't support it\n");
			return -1;
		}

		if (qspect && (cap & inst_spectral) == 0) {
			printf("Spectral information was requested,\n");
			printf("but instrument doesn't support it\n");
			qspect = 0;
		}

		mode = inst_mode_trans_spot;

		if (qspect)
			mode |= inst_mode_spectral;

		if ((rv = gt->set_mode(gt,mode)) != inst_ok) {
			error("couldn't set instrument  mode");
		}

		/* Do a calibration up front, so as not to get in the users way. */
		if (gt->needs_calibration(gt) & inst_calt_needs_cal_mask) {
			inst_code ev;
			ev = inst_handle_calibrate(gt, inst_calt_all, inst_calc_none, NULL);
			if (ev != inst_ok) {	/* Abort or fatal error */
				return -1;
			}
		}
	}

	/* iterate over all patches */
	for (i=0; i<npat; i++) {

		ipatch pat;

		printf("Reading patch id=%s loc=%s\n", cols[i].id, cols[i].loc);

		/* TODO: a clean way to save data and resume later... */
		for (;;) {
			if ((rv = gt->read_sample(gt,cols[i].id,&pat)) != inst_ok
			 && (rv & inst_mask) != inst_user_trig) {

				if ((rv & inst_mask) == inst_needs_cal) {
					inst_code ev;

					printf("Spectrolino needs to recalibrate. Please clear table.\n");

					printf("\nStrip read failed because instruments needs calibration\n");
					ev = inst_handle_calibrate(gt, inst_calt_all, inst_calc_none);
					if (ev != inst_ok) {	/* Abort or fatal error */
						error("Calibration failed or was aborted");
					}
				}
				else
					error("bad measurement... aborting");
			}
			else break;
		}

		/* copy and convert data */
		if (pat.XYZ_v)
			for (j=0; j<3; j++) cols[i].XYZ[j] = pat.XYZ[j];
		if (qspect && pat.sp.spec_n > 0)
			for (j=0; j<pat.sp.spec_n; j++) cols[i].spect[j] = pat.sp.spec[j]*100.0;

		printf("patch %s: RGB %f %f %f XYZ %f %f %f\n",
			   cols[i].id,
			   cols[i].r,cols[i].g,cols[i].b,
			   cols[i].XYZ[0],cols[i].XYZ[1],cols[i].XYZ[2]);

	}

	/* clean up */
	gt->del(gt);
	return 0;
}

void
usage(void) {
	icoms *icom;
	fprintf(stderr,"Read Film Target Test Chart colors, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Licensed under the GPL\n");
	fprintf(stderr,"usage: filmread [-v] [-s] [-c comport] outfile\n");
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -s              Gather spectral data also\n");
	fprintf(stderr," -c listno       Set communication port from the following list (default %d)\n",COMPORT);
	if ((icom = new_icoms()) != NULL) {
		icompath **paths;
		if ((paths = icom->get_paths(icom)) != NULL) {
			int i;
			for (i = 0; ; i++) {
				if (paths[i] == NULL)
					break;
				fprintf(stderr,"    %d: '%s'\n",i+1,paths[i]->path);
			}
		} else
			fprintf(stderr,"    ** No ports found **\n");
		icom->del(icom);
	}
	fprintf(stderr," outfile         Base name for input[.ti2]/output[.ti3] file\n");
	exit(1);
	}

int main(int argc, char *argv[])
{
	int i,j;
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	int spect = 0;
	int comport = COMPORT;	/* COM port used */
	static char inname[200] = { 0 };		/* Input cgats file base name */
	static char outname[200] = { 0 };		/* Output cgats file base name */
	cgats *icg;			/* input cgats structure */
	cgats *ocg;			/* output cgats structure */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	col *cols;			/* Internal storage of all the patch colors */
	int dim = 0;		/* Dimensionality - 1, 3, or 4 */
	int npat;			/* Number of patches */
	int si;				/* Sample id index */
	int li;				/* Location id index */
	int ti;				/* Temp index */
	int fi;				/* Colorspace index */
	int rstart = -1;	/* Random start index */

	set_exe_path(argv[0]);			/* Set global exe_path and error_program */

	if (argc <= 1)
		usage();

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++)
		{
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')		/* Look for any flags */
			{
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else
				{
				if ((fa+1) < argc)
					{
					if (argv[fa+1][0] != '-')
						{
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
						}
					}
				}

			if (argv[fa][1] == '?')
				usage();

			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V')
				verb = 1;

			else if (argv[fa][1] == 's' || argv[fa][1] == 'S')
				spect = 1;

			/* COM port  */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage();
				comport = atoi(na);
				if (comport < 1 || comport > 40) usage();
			}

			else 
				usage();
			}
		else
			break;
		}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(inname,argv[fa]);
	strcat(inname,".ti2");
	strcpy(outname,argv[fa]);
	strcat(outname,".ti3");

	icg = new_cgats();			/* Create a CGATS structure */
	icg->add_other(icg, "CTI2"); 	/* our special input type is Calibration Target Information 2 */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->ntables == 0 || icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI2 format file");
	if (icg->ntables != 1)
		error ("Input file doesn't contain exactly one table");

	if ((npat = icg->t[0].nsets) <= 0)
		error ("No sets of data");

	/* Setup output cgats file */
	ocg = new_cgats();	/* Create a CGATS structure */
	ocg->add_other(ocg, "CTI3"); 	/* our special type is Calibration Target Information 3 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll filmread", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	ocg->add_kword(ocg, 0, "DEVICE_CLASS","OUTPUT", NULL);	/* What sort of device this is */

	if ((ti = icg->find_kword(icg, 0, "APPROXIMATE_CHART")) >= 0)
		ocg->add_kword(ocg, 0, "APPROXIMATE_CHART",icg->t[0].kdata[ti], NULL);
	else {
		if ((ti = icg->find_kword(icg, 0, "SINGLE_DIM_STEPS")) >= 0)
			ocg->add_kword(ocg, 0, "SINGLE_DIM_STEPS",icg->t[0].kdata[ti], NULL);
	
		if ((ti = icg->find_kword(icg, 0, "MULTI_DIM_STEPS")) >= 0)
			ocg->add_kword(ocg, 0, "MULTI_DIM_STEPS",icg->t[0].kdata[ti], NULL);

		if ((ti = icg->find_kword(icg, 0, "FULL_SPREAD_PATCHES")) >= 0)
			ocg->add_kword(ocg, 0, "FULL_SPREAD_PATCHES",icg->t[0].kdata[ti], NULL);
	}
	

	if ((ti = icg->find_kword(icg, 0, "RANDOM_START")) >= 0)
		rstart = atoi(icg->t[0].kdata[ti]);

	if (verb) {
		printf("Total patches to read = %d\n", npat);
		if (rstart >=0)
			printf("Random starting index = %d\n",rstart);
	}

	/* Fields we want */
	ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);

	if ((si = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
		error ("Input file doesn't contain field SAMPLE_ID");
	if (icg->t[0].ftype[si] != nqcs_t)
		error ("Field SAMPLE_ID is wrong type");

	if ((li = icg->find_field(icg, 0, "SAMPLE_LOC")) < 0)
		error ("Input file doesn't contain field SAMPLE_LOC");
	if (icg->t[0].ftype[li] != cs_t)
		error ("Field SAMPLE_LOC is wrong type");

	if ((cols = (col *)malloc(sizeof(col) * npat)) == NULL)
		error("Malloc failed!");

	/* Figure out the color space */
	if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
		error ("Input file doesn't contain keyword COLOR_REPS");
	if (strcmp(icg->t[0].kdata[fi],"CMYK") == 0) {
		error ("We don't support CMYK yet");
	} else if (strcmp(icg->t[0].kdata[fi],"RGB") == 0) {
		int ri, gi, bi;
		dim = 3;
		if ((ri = icg->find_field(icg, 0, "RGB_R")) < 0)
			error ("Input file doesn't contain field RGB_R");
		if (icg->t[0].ftype[ri] != r_t)
			error ("Field RGB_R is wrong type");
		if ((gi = icg->find_field(icg, 0, "RGB_G")) < 0)
			error ("Input file doesn't contain field RGB_G");
		if (icg->t[0].ftype[gi] != r_t)
			error ("Field RGB_G is wrong type");
		if ((bi = icg->find_field(icg, 0, "RGB_B")) < 0)
			error ("Input file doesn't contain field RGB_B");
		if (icg->t[0].ftype[bi] != r_t)
			error ("Field RGB_B is wrong type");
		ocg->add_field(ocg, 0, "RGB_R", r_t);
		ocg->add_field(ocg, 0, "RGB_G", r_t);
		ocg->add_field(ocg, 0, "RGB_B", r_t);
		ocg->add_kword(ocg, 0, "COLOR_REP","RGB_XYZ", NULL);
		ocg->add_field(ocg, 0, "XYZ_X", r_t);
		ocg->add_field(ocg, 0, "XYZ_Y", r_t);
		ocg->add_field(ocg, 0, "XYZ_Z", r_t);

		if (spect) {
			/* transmission spectrum from 380nm to 730nm */
			ocg->add_kword(ocg, 0, "SPECTRAL_BANDS","36", NULL);
			ocg->add_kword(ocg, 0, "SPECTRAL_START_NM","380.0", NULL);
			ocg->add_kword(ocg, 0, "SPECTRAL_END_NM","730.0", NULL);

			for (i=0; i<36; i++) {
				char buf[20];
				sprintf(buf, "SPEC_%03d", 380+i*10);
				ocg->add_field(ocg, 0, buf, r_t);
			}
		}

		for (i = 0; i < npat; i++) {
			cols[i].id = ((char *)icg->t[0].fdata[i][si]);
			cols[i].loc = ((char *)icg->t[0].fdata[i][li]);
			cols[i].r = *((double *)icg->t[0].fdata[i][ri]) / 100.0;
			cols[i].g = *((double *)icg->t[0].fdata[i][gi]) / 100.0;
			cols[i].b = *((double *)icg->t[0].fdata[i][bi]) / 100.0;
			cols[i].XYZ[0] = cols[i].XYZ[1] = cols[i].XYZ[2] = -1.0;
			for (j=0; j<36; j++) cols[i].spect[j] = -1.0;
		}
	} else if (strcmp(icg->t[0].kdata[fi],"W") == 0) {
		error ("We don't support GRAY yet");
	} else
		error ("Input file keyword COLOR_REPS has unknown value");

	/* Read all of the patches in */
	read_patches(cols, npat, comport, verb, spect);

	/* Write out the patch info to the output CGATS file */
	for (i = 0; i < npat; i++) {
		if (dim == 4)
			ocg->add_set(ocg, 0, cols[i].id, 100.0 * cols[i].cmyk[0], 100.0 * cols[i].cmyk[1],
			                         100.0 * cols[i].cmyk[2], 100.0 * cols[i].cmyk[3],
			                         cols[i].XYZ[0], cols[i].XYZ[1], cols[i].XYZ[2]);
		else if (dim == 3)

			if (!spect)
				ocg->add_set(ocg, 0, cols[i].id, 100.0 * cols[i].r,
						 100.0 * cols[i].g, 100.0 * cols[i].b,
						 cols[i].XYZ[0], cols[i].XYZ[1], cols[i].XYZ[2]);
			else
				ocg->add_set(ocg, 0, cols[i].id, 100.0 * cols[i].r,
						 100.0 * cols[i].g, 100.0 * cols[i].b,
						 cols[i].XYZ[0], cols[i].XYZ[1], cols[i].XYZ[2],
						 cols[i].spect[0], cols[i].spect[1], cols[i].spect[2],
						 cols[i].spect[3], cols[i].spect[4], cols[i].spect[5],
						 cols[i].spect[6], cols[i].spect[7], cols[i].spect[8],
						 cols[i].spect[9], cols[i].spect[10], cols[i].spect[11],
						 cols[i].spect[12], cols[i].spect[13], cols[i].spect[14],
						 cols[i].spect[15], cols[i].spect[16], cols[i].spect[17],
						 cols[i].spect[18], cols[i].spect[19], cols[i].spect[20],
						 cols[i].spect[21], cols[i].spect[22], cols[i].spect[23],
						 cols[i].spect[24], cols[i].spect[25], cols[i].spect[26],
						 cols[i].spect[27], cols[i].spect[28], cols[i].spect[29],
						 cols[i].spect[30], cols[i].spect[31], cols[i].spect[32],
						 cols[i].spect[33], cols[i].spect[34], cols[i].spect[35]);


		else if (dim == 1)
			ocg->add_set(ocg, 0, cols[i].id, 100.0 * cols[i].gy,
			                         cols[i].XYZ[0], cols[i].XYZ[1], cols[i].XYZ[2]);
	}

	if (ocg->write_name(ocg, outname))
		error("Write error : %s",ocg->err);

	free(cols);
	ocg->del(ocg);		/* Clean up */
	icg->del(icg);		/* Clean up */

	return 0;
}


