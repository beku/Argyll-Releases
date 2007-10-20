/* 
 * Argyll Color Correction System
 * Split a .ti3 file into two parts.
 *
 * Author: Graeme W. Gill
 * Date:   14/12/2005
 *
 * Copyright 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * This program takes in a CGATS .ti3 file, and splits it into
 * two .ti3 files, spreading the readings between them. 
 * This is intended for use in verifying the profiler.
 */

/*
 * TTBD:
 */

#define DEBUG

#define verbo stdout

#include <stdio.h>
#include <string.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include <sys/types.h>
#include <time.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "cgats.h"
#include "xicc.h"
#include "insttypes.h"
#include "sort.h"

void
usage(void) {
	fprintf(stderr,"Split a .ti3 into two, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	fprintf(stderr,"usage: splitcgats [-options] input.ti3 output1.ti3 output2.ti3\n");
	fprintf(stderr," -v              Verbose - print each patch value\n");
	fprintf(stderr," -n no           Put no sets in first file, and balance in second file.\n");
	fprintf(stderr," -p percent      Put percent%% sets in first file, and balance in second file. (def. 50%%)\n");
	fprintf(stderr," -r seed         Use given random seed.\n");
	fprintf(stderr," input.ti3       File to be split up.\n");
	fprintf(stderr," output1.ti3     First output file\n");
	fprintf(stderr," output2.ti3     Second output file\n");
	exit(1);
}

int main(int argc, char *argv[]) {
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	int numb = -1;			/* Number to put in first */
	double prop = 0.5;		/* Proportion to put in first */
	int seed = 0x12345678;
	int doseed = 0;

	cgats *cgf = NULL;			/* cgats file data */
	char in_name[MAXNAMEL+1];	/* Patch filename  */

	cgats *cg1 = NULL;			/* cgats file data */
	char out_name1[MAXNAMEL+4+1]; /* VRML name */
	cgats *cg2 = NULL;			/* cgats file data */
	char out_name2[MAXNAMEL+4+1]; /* VRML name */

	cgats_set_elem *setel;		/* Array of set value elements */
	int *flags;					/* Point to destination of set */

	int i, n;

	error_program = "spitcgats";

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

			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;

			} 
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N') {
				fa = nfa;
				if (na == NULL) usage();
				numb = atoi(na);
				if (numb < 0) usage();
			}
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				fa = nfa;
				if (na == NULL) usage();
				numb = atoi(na);
				if (numb < 0) usage();
			}
			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R') {
				fa = nfa;
				if (na == NULL) usage();
				seed = atoi(na);
				doseed = 1;
			}

			else 
				usage();
		} else
			break;
	}

	if (doseed)
		rand32(seed);			/* Init seed deterministicaly */
	else
		rand32(time(NULL));		/* Init seed randomly */

	/* Get the file name arguments */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(in_name,argv[fa++],MAXNAMEL); in_name[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(out_name1,argv[fa++],MAXNAMEL); out_name1[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(out_name2,argv[fa++],MAXNAMEL); out_name2[MAXNAMEL] = '\000';

	if ((cgf = new_cgats()) == NULL)
		error("Failed to create cgats object");
	cgf->add_other(cgf, ""); 	/* Allow any signature file */
	
	if (cgf->read_name(cgf, in_name))
		error("CGATS file '%s' read error : %s",in_name,cgf->err);
	
	if (cgf->ntables < 1)
		error ("Input file '%s' doesn't contain at least one table",in_name);
	
	/* Create the two output files */
	if ((cg1 = new_cgats()) == NULL)
		error("Failed to create cgats object");
	if ((cg2 = new_cgats()) == NULL)
		error("Failed to create cgats object");

	/* Duplicate the type of the file */
	if (cgf->t[0].tt == cgats_X) {
		cg1->add_other(cg1, cgf->cgats_type);
		cg1->add_table(cg1, tt_other, 0);
		cg2->add_other(cg2, cgf->cgats_type);
		cg2->add_table(cg2, tt_other, 0);
	} else if (cgf->t[0].tt == tt_other) {
		cg1->add_other(cg1, cgf->others[cgf->t[0].oi]);
		cg1->add_table(cg1, tt_other, 0);
		cg2->add_other(cg2, cgf->others[cgf->t[0].oi]);
		cg2->add_table(cg2, tt_other, 0);
	} else {
		cg1->add_table(cg1, cgf->t[0].tt, 0);
		cg2->add_table(cg1, cgf->t[0].tt, 0);
	}

	/* Duplicate all the keywords */
	for (i = 0; i < cgf->t[0].nkwords; i++) {
		cg1->add_kword(cg1, 0, cgf->t[0].ksym[i], cgf->t[0].kdata[i], NULL);
		cg2->add_kword(cg2, 0, cgf->t[0].ksym[i], cgf->t[0].kdata[i], NULL);
	}

	/* Duplicate all of the fields */
	for (i = 0; i < cgf->t[0].nfields; i++) {
		cg1->add_field(cg1, 0, cgf->t[0].fsym[i], cgf->t[0].ftype[i]);
		cg2->add_field(cg2, 0, cgf->t[0].fsym[i], cgf->t[0].ftype[i]);
	}

	if ((setel = (cgats_set_elem *)malloc(
	     sizeof(cgats_set_elem) * cgf->t[0].nfields)) == NULL)
		error("Malloc failed!");

	if ((flags = (int *)calloc(cgf->t[0].nsets, sizeof(int))) == NULL)
		error("Malloc failed!");
	
	if (numb < 0) {	/* Use percentage */
		numb = (int)(cgf->t[0].nsets * prop + 0.5);
	}
	if (numb > cgf->t[0].nsets)
		numb = cgf->t[0].nsets;

	if (verb)
		printf("Putting %d sets in '%s' and %d in '%s'\n",numb,out_name1,cgf->t[0].nsets-numb,out_name2);

	/* Chose which of the sets go into file 1 */
	for (n = 0; n < numb;) {
		i = i_rand(0, cgf->t[0].nsets-1);
		if (flags[i] == 0) {
			flags[i] = 1;
			n++;
		}
	}

	/* Copy themn approproately */
	for (i = 0; i < cgf->t[0].nsets; i++) {
		cgf->get_setarr(cgf, 0, i, setel);
		if (flags[i]) {
			cg1->add_setarr(cg1, 0, setel);
		} else {
			cg2->add_setarr(cg2, 0, setel);
		}
	}

	/* Write out the files */
	if (cg1->write_name(cg1, out_name1))
		error("CGATS file '%s' write error : %s",out_name1,cg1->err);
	if (cg2->write_name(cg2, out_name2))
		error("CGATS file '%s' write error : %s",out_name2,cg2->err);
	

	free(flags);
	free(setel);
	
	return 0;
}





