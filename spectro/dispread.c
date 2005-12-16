/* 
 * Argyll Color Correction System
 * DTP92/Spectrolino display target reader
 *
 * Author: Graeme W. Gill
 * Date:   4/10/96
 *
 * Copyright 1996 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* This program displays test patches, and takes readings from a display device */

/* TTBD
 *
 */

#undef DEBUG
#undef DEBUG_OFFSET	/* Keep test window out of the way */

#define COMPORT 1		/* Default com port 1..4 */
#define VERBOUT stdout

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "xspect.h"
#include "cgats.h"
#include "serio.h"
#include "inst.h"
#include "dtp92.h"
#include "gretag.h"
#include "dispwin.h"
#include "dispsup.h"
#include "sort.h"
#if defined (NT)
#include <conio.h>
#endif

void
usage(void) {
	serio *sio;
	fprintf(stderr,"Read a Display, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	fprintf(stderr,"usage: dispread [-v] [-c comport] [-a] [-i inst] [-s]%s outfile\n",
#if defined(UNIX) && !defined(__APPLE__)
	" [-n]"
#else
	""
#endif
    );
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -d              Print debug diagnostics to stderr\n");
	fprintf(stderr," -c comport      Set COM port, 1..N (default %d)\n",COMPORT);
	if ((sio = new_serio()) != NULL) {
		char **paths;
		if ((paths = sio->get_paths(sio)) != NULL) {
			int i;
			for (i = 0; ; i++) {
				if (paths[i] == NULL)
					break;
				fprintf(stderr,"    %d = '%s'\n",i+1,paths[i]);
			}
		} else
			fprintf(stderr,"    ** No ports found **\n");
		sio->del(sio);
	}
	fprintf(stderr," -a              Run instrument calibration\n");
 	fprintf(stderr," -i 92|SO        Select target instrument (default DTP92)\n");
	fprintf(stderr,"                 92 = DTP92, SO = Spectrolino\n");
	fprintf(stderr," -k file.cal     Apply display calibration file while reading\n");
	fprintf(stderr," -s              Save spectral information (default don't save)\n");
#if defined(UNIX) && !defined(__APPLE__)
	fprintf(stderr," -n              Don't set override redirect on test window\n");
#endif
	fprintf(stderr," outfile         Base name for input[ti1]/output[ti3] file\n");
	exit(1);
	}

main(int argc, char *argv[])
{
	int i,j;
	int fa,nfa;							/* current argument we're looking at */
	double patsize = 80.0;				/* size of displayed color patch */
	double ho = 0.0, vo = 0.0;
	int verb = 0;
	int debug = 0;
	int override = 1;					/* Override redirect on X11 */
	int comport = COMPORT;				/* COM port used */
	instType itype = instDTP92;			/* Default target instrument */
	int doCalib = 0;					/* Do a calibration */
	int spectral = 0;					/* Don't save spectral information */
	char inname[MAXNAMEL+1] = "\000";	/* Input cgats file base name */
	char outname[MAXNAMEL+1] = "\000";	/* Output cgats file base name */
	char calname[MAXNAMEL+1] = "\000";	/* Calibration file name */
	double cal[3][256];					/* Display calibration */
	cgats *icg;			/* input cgats structure */
	cgats *ocg;			/* output cgats structure */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	col *cols;			/* Internal storage of all the patch colors */
	int dim;			/* Dimensionality - 1, 3, or 4 */
	int npat;			/* Number of patches/colors */
	int xpat = 0;		/* Set to number of extra patches */
	int wpat;			/* Set to index of white patch */
	int si;				/* Sample id index */
	int li;				/* Location id index */
	int ti;				/* Temp index */
	int fi;				/* Colorspace index */
	int nsetel = 0;
	cgats_set_elem *setel;	/* Array of set value elements */
	disprd *dr;			/* Display patch read object */
	int rv;

	error_program = "Dispread";
	if (argc <= 1)
		usage();

	/* Process the arguments */
	for (fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {		/* Look for any flags */
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

			else if (argv[fa][1] == 'd' || argv[fa][1] == 'D')
				debug = 1;

#if defined(UNIX) && !defined(__APPLE__)
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N')
				override = 0;
#endif /* UNIX */
			/* COM port  */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage();
				comport = atoi(na);
				if (comport < 1 || comport > 4) usage();
			}
			/* Target Instrument */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				char *p;
				int tt;
				fa = nfa;
				if (na == NULL) usage();

				if (strcmp("92", na) == 0) {
					itype = instDTP92;
					patsize = 80.0;
				}
				else if (strcmp("so", na) == 0
				      || strcmp("SO", na) == 0) {
					itype = instSpectrolino;
					patsize = 60.0;
				}
				else
					usage();
			}
			else if (argv[fa][1] == 'a' || argv[fa][1] == 'A') {
				doCalib = 1;

			}
			/* Calibration file */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				char *p;
				int tt;
				fa = nfa;
				if (na == NULL) usage();
				strncpy(calname,na,MAXNAMEL); calname[MAXNAMEL] = '\000';
			}
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				spectral = 1;

			} else 
				usage();
			}
		else
			break;
	}

	if (doCalib) {
		if ((rv = disprd_calibration(itype, comport, override, patsize, verb, debug)) != 0) {
			error("docalibration failed with return value %d\n",rv);
		}
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(inname,argv[fa++],MAXNAMEL-4); inname[MAXNAMEL-4] = '\000';
	strcpy(outname,inname);
	strcat(inname,".ti1");
	strcat(outname,".ti3");

	icg = new_cgats();			/* Create a CGATS structure */
	icg->add_other(icg, "CTI1"); 	/* our special input type is Calibration Target Information 1 */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI1 format file");
	if (icg->ntables != 2)		/* We don't use second table at the moment */
		error ("Input file doesn't contain exactly two tables");

	if ((npat = icg->t[0].nsets) <= 0)
		error ("No sets of data");

	/* Setup output cgats file */
	ocg = new_cgats();					/* Create a CGATS structure */
	ocg->add_other(ocg, "CTI3"); 		/* our special type is Calibration Target Information 3 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll dispread", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	ocg->add_kword(ocg, 0, "DEVICE_CLASS","DISPLAY", NULL);	/* What sort of device this is */

	if ((ti = icg->find_kword(icg, 0, "SINGLE_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "SINGLE_DIM_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "COMP_GREY_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "COMP_GREY_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "MULTI_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "MULTI_DIM_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "FULL_SPREAD_PATCHES")) >= 0)
		ocg->add_kword(ocg, 0, "FULL_SPREAD_PATCHES",icg->t[0].kdata[ti], NULL);
	
	if (verb) {
		printf("Number of patches = %d\n",npat);
	}

	/* Fields we want */
	ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);

	if ((si = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
		error ("Input file doesn't contain field SAMPLE_ID");
	if (icg->t[0].ftype[si] != nqcs_t)
		error ("Field SAMPLE_ID is wrong type");

	if ((cols = (col *)malloc(sizeof(col) * (npat+1))) == NULL)
		error("Malloc failed!");

	/* Figure out the color space */
	if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
		error ("Input file doesn't contain keyword COLOR_REPS");
	if (strcmp(icg->t[0].kdata[fi],"RGB") == 0) {
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
		for (i = 0; i < npat; i++) {
			cols[i].id = ((char *)icg->t[0].fdata[i][si]);
			cols[i].r = *((double *)icg->t[0].fdata[i][ri]) / 100.0;
			cols[i].g = *((double *)icg->t[0].fdata[i][gi]) / 100.0;
			cols[i].b = *((double *)icg->t[0].fdata[i][bi]) / 100.0;
			cols[i].XYZ[0] = cols[i].XYZ[1] = cols[i].XYZ[2] = -1.0;
		}
	} else
		error ("Input file keyword COLOR_REPS has illegal value");

	/* Check that there is a white patch, and if not, add one */
	for (wpat = 0; wpat < npat; wpat++) {
		if (cols[wpat].r > 0.9999999 && 
		    cols[wpat].g > 0.9999999 && 
		    cols[wpat].b > 0.9999999) {
			break;
		}
	}
	if (wpat >= npat) {	/* Create a white patch */
		xpat = 1;
		cols[wpat].r = cols[wpat].g = cols[wpat].b = 1.0;
	}

	/* Setup a display calibration set if we are given one */
	if (calname[0] != '\000') {
		cgats *ccg;			/* calibration cgats structure */
		int ncal;
		int ii, ri, gi, bi;
		
		ccg = new_cgats();			/* Create a CGATS structure */
		ccg->add_other(ccg, "CAL"); /* our special calibration type */
	
		if (ccg->read_name(ccg, calname))
			error("CGATS calibration file read error %s on file '%s'",ccg->err,calname);
	
		if (ccg->t[0].tt != tt_other || ccg->t[0].oi != 0)
			error ("Calibration file isn't a CAL format file");
		if (ccg->ntables < 1)
			error ("Calibration file '%s' doesn't contain at least one table",calname);
	
		if ((ncal = ccg->t[0].nsets) <= 0)
			error ("No data in set of file '%s'",calname);
	
		if (ncal != 256)
			error ("Expect 256 data sets in file '%s'",calname);
	
		if ((fi = ccg->find_kword(ccg, 0, "DEVICE_CLASS")) < 0)
			error ("Calibration file '%s' doesn't contain keyword COLOR_REPS",calname);
		if (strcmp(ccg->t[0].kdata[fi],"DISPLAY") != 0)
			error ("Calibration file '%s' doesn't have DEVICE_CLASS of DISPLAY",calname);

		if ((ii = ccg->find_field(ccg, 0, "RGB_I")) < 0)
			error ("Calibration file '%s' doesn't contain field RGB_I",calname);
		if (ccg->t[0].ftype[ii] != r_t)
			error ("Field RGB_R in file '%s' is wrong type",calname);
		if ((ri = ccg->find_field(ccg, 0, "RGB_R")) < 0)
			error ("Calibration file '%s' doesn't contain field RGB_R",calname);
		if (ccg->t[0].ftype[ri] != r_t)
			error ("Field RGB_R in file '%s' is wrong type",calname);
		if ((gi = ccg->find_field(ccg, 0, "RGB_G")) < 0)
			error ("Calibration file '%s' doesn't contain field RGB_G",calname);
		if (ccg->t[0].ftype[gi] != r_t)
			error ("Field RGB_G in file '%s' is wrong type",calname);
		if ((bi = ccg->find_field(ccg, 0, "RGB_B")) < 0)
			error ("Calibration file '%s' doesn't contain field RGB_B",calname);
		if (ccg->t[0].ftype[bi] != r_t)
			error ("Field RGB_B in file '%s' is wrong type",calname);
		for (i = 0; i < ncal; i++) {
			cal[0][i] = *((double *)ccg->t[0].fdata[i][ri]);
			cal[1][i] = *((double *)ccg->t[0].fdata[i][gi]);
			cal[2][i] = *((double *)ccg->t[0].fdata[i][bi]);
		}
		ccg->del(ccg);
	} else {
		cal[0][0] = -1.0;	/* Not used */
	}

#ifdef DEBUG_OFFSET
	ho = 0.8;
	vo = -0.8;
#endif

	if ((dr = new_disprd(itype, comport, 0, cal, override, patsize, ho, vo,
	                     spectral, verb, VERBOUT, debug)) == NULL)
		error("dispread failed to create a disprd object\n");

	/* Test the CRT with all of the test points */
	if ((rv = dr->read(dr, cols, npat + xpat, 0, 0)) != 0) {
		dr->del(dr);
		error("test_crt returned error code %d\n",rv);
	}
	dr->del(dr);

	/* And save the result: */

	/* Convert from absolute XYZ to relative XYZ if that's all we have */
	if (npat > 0 && cols[0].XYZ_v == 0) {
		double nn;

		if (npat > 0 && cols[0].aXYZ_v == 0)
			error("Neither relative or absolute XYZ is valid!");

		nn = 100.0 / cols[wpat].aXYZ[1];		/* Normalise Y of white to 100 */

		for (i = 0; i < (npat+ xpat); i++) {
			for (j = 0; j < 3; j++) {
				cols[i].XYZ[j] = nn * cols[i].aXYZ[j];
				cols[i].XYZ_v = 1;
			}
		}
	}

	nsetel += 1;		/* For id */
	nsetel += dim;		/* For device values */
	nsetel += 3;		/* For XYZ */

	/* If we have spectral information, output it too */
	if (npat > 0 && cols[0].spec_n > 0) {
		char buf[100];

		nsetel += cols[0].spec_n;		/* Spectral values */
		sprintf(buf,"%d", cols[0].spec_n);
		ocg->add_kword(ocg, 0, "SPECTRAL_BANDS",buf, NULL);
		sprintf(buf,"%f", cols[0].spec_wl_short);
		ocg->add_kword(ocg, 0, "SPECTRAL_START_NM",buf, NULL);
		sprintf(buf,"%f", cols[0].spec_wl_long);
		ocg->add_kword(ocg, 0, "SPECTRAL_END_NM",buf, NULL);

		/* Generate fields for spectral values */
		for (i = 0; i < cols[0].spec_n; i++) {
			int nm;
	
			/* Compute nearest integer wavelength */
			nm = (int)(cols[0].spec_wl_short + ((double)i/(cols[0].spec_n-1.0))
			            * (cols[0].spec_wl_long - cols[0].spec_wl_short) + 0.5);
			
			sprintf(buf,"SPEC_%03d",nm);
			ocg->add_field(ocg, 0, buf, r_t);
		}
	}

	if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * nsetel)) == NULL)
		error("Malloc failed!");

	/* Write out the patch info to the output CGATS file */
	for (i = 0; i < npat; i++) {
		int k = 0;

		setel[k++].c = cols[i].id;
		setel[k++].d = 100.0 * cols[i].r;
		setel[k++].d = 100.0 * cols[i].g;
		setel[k++].d = 100.0 * cols[i].b;

		setel[k++].d = cols[i].XYZ[0];
		setel[k++].d = cols[i].XYZ[1];
		setel[k++].d = cols[i].XYZ[2];

		for (j = 0; j < cols[i].spec_n; j++) {
			setel[k++].d = cols[i].spec[j];
		}

		ocg->add_setarr(ocg, 0, setel);
	}

	free(setel);

	/* If we have the absolute brightness of the display, record it */
	if (cols[wpat].aXYZ_v != 0) {
		char buf[100];

		sprintf(buf,"%f %f %f", cols[wpat].aXYZ[0], cols[wpat].aXYZ[1], cols[wpat].aXYZ[2]);
		ocg->add_kword(ocg, 0, "LUMINANCE_XYZ_CDM2",buf, NULL);
	}

	/* Write out the calibration if we have it */
	if (cal[0][0] >= 0.0) {
		int calres = 256;	/* 256 steps in calibration */

		ocg->add_other(ocg, "CAL"); 		/* our special type is Calibration file */
		ocg->add_table(ocg, tt_other, 1);	/* Add another table for RAMDAC values */
		ocg->add_kword(ocg, 1, "DESCRIPTOR", "Argyll Device Calibration State",NULL);
		ocg->add_kword(ocg, 1, "ORIGINATOR", "Argyll dispread", NULL);
		ocg->add_kword(ocg, 1, "CREATED",atm, NULL);

		ocg->add_kword(ocg, 1, "DEVICE_CLASS","DISPLAY", NULL);

		ocg->add_field(ocg, 1, "RGB_I", r_t);
		ocg->add_field(ocg, 1, "RGB_R", r_t);
		ocg->add_field(ocg, 1, "RGB_G", r_t);
		ocg->add_field(ocg, 1, "RGB_B", r_t);

		if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * 4)) == NULL)
			error("Malloc failed!");

		for (i = 0; i < calres; i++) {
			double vv = i/(calres-1.0);

			setel[0].d = vv;
			setel[1].d = cal[0][i];
			setel[2].d = cal[1][i];
			setel[3].d = cal[2][i];

			ocg->add_setarr(ocg, 1, setel);
		}

		free(setel);
	}
	if (ocg->write_name(ocg, outname))
		error("Write error : %s",ocg->err);

	free(cols);
	ocg->del(ocg);		/* Clean up */
	icg->del(icg);		/* Clean up */

	return 0;
}


