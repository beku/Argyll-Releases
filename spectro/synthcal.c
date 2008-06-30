
/* 
 * Argyll Color Correction System
 * Create a linear display calibration file.
 *
 * Author: Graeme W. Gill
 * Date:   15/11/2005
 *
 * Copyright 1996 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#undef DEBUG

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "cgats.h"

#include "sort.h"

#include <stdarg.h>

#if defined (NT)
#include <conio.h>
#endif

void
usage(void) {
	fprintf(stderr,"Create linear display calibration file, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	fprintf(stderr,"usage: displin outfile\n");
	fprintf(stderr," -s s1,s2,s3,    Set non-linear curve scale (default 1.0)\n");
	fprintf(stderr," -p p1,p2,p3,    Set non-linear curve powers (default 1.0)\n");
	fprintf(stderr," outfile         Base name for output .cal file\n");
	exit(1);
	}

int main(int argc, char *argv[])
{
	int fa,nfa;							/* current argument we're looking at */
	int verb = 0;
	static char outname[200] = { 0 };	/* Output cgats file base name */
	double max[3] = { 1.0, 1.0, 1.0 };	/* Output scale for RGB */
	double gam[3] = { 1.0, 1.0, 1.0 };	/* Gamm applied */

	error_program = "displin";
	if (argc <= 1)
		usage();

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')		/* Look for any flags */
			{
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1) < argc) {
					if (argv[fa+1][0] != '-')
						{
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage();

			else if (argv[fa][1] == 'v')
				verb = 1;

			/* curve scale */
			else if (argv[fa][1] == 's' || argv[fa][1] == 'S') {
				fa = nfa;
				if (na == NULL) usage();
				if (sscanf(na, " %lf,%lf,%lf ", &max[0], &max[1], &max[2]) != 3)
					usage();
			}

			/* curve power */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				fa = nfa;
				if (na == NULL) usage();
				if (sscanf(na, " %lf,%lf,%lf ", &gam[0], &gam[1], &gam[2]) != 3)
					usage();
			}

			else 
				usage();
		} else
			break;
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(outname,argv[fa]);
	strcat(outname,".cal");

	/* Write out the resulting calibration file */
	{
		int i, j, calres = 256;	/* 256 steps in calibration */
		cgats *ocg;			/* output cgats structure */
		time_t clk = time(0);
		struct tm *tsp = localtime(&clk);
		char *atm = asctime(tsp);	/* Ascii time */
		cgats_set_elem *setel;		/* Array of set value elements */

		ocg = new_cgats();				/* Create a CGATS structure */
		ocg->add_other(ocg, "CAL"); 	/* our special type is Calibration file */

		ocg->add_table(ocg, tt_other, 0);	/* Add a table for RAMDAC values */
		ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Device Calibration State",NULL);
		ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll dispcal", NULL);
		atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
		ocg->add_kword(ocg, 0, "CREATED",atm, NULL);

		ocg->add_kword(ocg, 0, "DEVICE_CLASS","DISPLAY", NULL);

		ocg->add_field(ocg, 0, "RGB_I", r_t);
		ocg->add_field(ocg, 0, "RGB_R", r_t);
		ocg->add_field(ocg, 0, "RGB_G", r_t);
		ocg->add_field(ocg, 0, "RGB_B", r_t);

		if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * 4)) == NULL)
			error("Malloc failed!");

		for (i = 0; i < calres; i++) {
			double vv = i/(calres-1.0);

			setel[0].d = vv;
			for (j = 0; j < 3; j++) {
				setel[j+1].d = max[j] * pow(vv, gam[j]);
				if (setel[j+1].d < 0.0)
					setel[j+1].d = 0.0;
				else if (setel[j+1].d > 1.0)
					setel[j+1].d = 1.0;
			}

			ocg->add_setarr(ocg, 0, setel);
		}

		free(setel);

		if (ocg->write_name(ocg, outname))
			error("Write error : %s",ocg->err);

		ocg->del(ocg);		/* Clean up */
	}
	return 0;
}


