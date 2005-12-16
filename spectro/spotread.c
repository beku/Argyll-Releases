/* 
 * Argyll Color Correction System
 * Spectrometer/Colorimeter color spot reader utility
 *
 * Author: Graeme W. Gill
 * Date:   3/10/2001
 *
 * Derived from printread.c
 * Was called printspot.
 *
 * Copyright 2001 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* This program reads a spot reflection/transmission/emission value using */
/* a spectrometer or colorimeter. */

/* TTBD
 *
 *     Complete support for spectrolino and spectroscan.
 *
 */

#define DEBUG

#define COMPORT 1		/* Default com port 1..4 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include "copyright.h"
#include "config.h"
#include "cgats.h"
#include "xicc.h"
#include "serio.h"
#include "inst.h"
#include "sort.h"

#include <stdarg.h>

#if defined (NT)
#include <conio.h>
#endif

/* Convert control chars to ^[A-Z] notation in a string */
static char *
fix_asciiz(char *s) {
	static char buf [200];
	char *d;
	for(d = buf; ;) {
		if (*s < ' ' && *s > '\000') {
			*d++ = '^';
			*d++ = *s++ + '@';
		} else
		*d++ = *s++;
		if (s[-1] == '\000')
			break;
	}
	return buf;
}

/* Deal with an instrument error. */
/* Return 0 to retry, 1 to abort */
static int ierror(inst *it, inst_code ic) {
	int c;
	empty_con_chars();
	printf("Got '%s' (%s) error.\nHit Esc to give up, any other key to retry:",
	       it->inst_interp_error(it, ic), it->interp_error(it, ic));
	fflush(stdout);
	c = next_con_char();
	printf("\n");
	if (c == 0x1b)	/* Escape */
		return 1;
	return 0;
}

/* A color structure */
/* This can hold all representations simultaniously */
typedef struct {
	double gy;
	double r,g,b;
	double cmyk[4];
	double XYZ[4];		/* Colorimeter readings */
	double eXYZ[4];		/* Expected XYZ values */

	int    spec_n;					/* Number of spectral bands, 0 if not valid */
	double spec_wl_short;			/* First reading wavelength in nm (shortest) */
	double spec_wl_long;			/* Last reading wavelength in nm (longest) */
	double spec[INSTR_MAX_BANDS];	/* Spectral reflectance %, shortest to longest */
	
	char *id;			/* Id string */
	char *loc;			/* Location string */
	int loci;			/* Location integer = pass * 256 + step */
} col;

void
usage(void) {
	serio *sio;
	fprintf(stderr,"Read Print Spot values, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	fprintf(stderr,"usage: spotread [-options]\n");
	fprintf(stderr," -v                   Verbose mode\n");
	fprintf(stderr," -d                   Print debug diagnostics to stderr\n");
	fprintf(stderr," -s                   Print spectrum for each reading\n");
	fprintf(stderr," -c comport           Set COM port, 1..N (default %d)\n",COMPORT);
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
	fprintf(stderr," -i 41 | 51 | 92 | SO | SS Select target instrument (no default)\n");
	fprintf(stderr," -t                   Use transmission measurement mode\n");
	fprintf(stderr," -e                   Use emissive  measurement mode (absolute results)\n");
	fprintf(stderr," -l illum             Choose illuminant for print/transparency spectral data:\n");
	fprintf(stderr,"                      A, D50 (def.), D65, F5, F8, F10 or file.sp\n");
	fprintf(stderr," -o observ            Choose CIE Observer for spectral data:\n");
	fprintf(stderr,"                      1931_2, 1964_10, S&B 1955_2, J&V 1978_2 (def.)\n");
	fprintf(stderr,"                      (Choose FWA during operation)\n");
	exit(1);
	}

main(int argc, char *argv[])
{
	int i,j;
	int fa,nfa;						/* current argument we're looking at */
	int verb = 0;
	int debug = 0;
	int pspec = 0;					/* Print out the spectrum for each reading */
	int trans = 0;					/* Use transmissioin mode */
	int emiss = 0;					/* Use emissive mode */
	int comport = COMPORT;			/* COM port used */
	instType itype = instUnknown;	/* No default target instrument */
	inst_capability cap = inst_unknown;	/* Instrument capabilities */
	double lx, ly;					/* Read location on xy table */
	baud_rate br = baud_38400;		/* Target baud rate */
	inst *it;						/* Instrument object */
	inst_code rv;
	int spec = 0;				/* Use spectral data flag */
	icxIllumeType illum = icxIT_D50;	/* Spectral defaults */
	xspect cust_illum;			/* Custom illumination spectrum */
	icxObserverType observ = icxOT_Judd_Voss_2;
	xspect sp;
	xsp2cie *sp2cie = NULL;			/* default conversion */
	xsp2cie *sp2cief[26];			/* FWA corrected conversions */
	double Lab[3] = { -10.0, 0, 0};	/* Last Lab */
	double rLab[3] = { -10.0, 0, 0};	/* Reference Lab */

	error_program = "SpotRead";

	for (i = 0; i < 26; i++)
		sp2cief[i] = NULL;

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

			else if (argv[fa][1] == 'd' || argv[fa][1] == 'D')
				debug = 1;

			else if (argv[fa][1] == 's' || argv[fa][1] == 'S')
				pspec = 1;

			/* COM port  */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage();
				comport = atoi(na);
				if (comport < 1 || comport > 4) usage();
			}
			/* Target Instrument size */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				char *p;
				int tt;
				fa = nfa;
				if (na == NULL) usage();

				if (strcmp("92", na) == 0)
					itype = instDTP92;
				else if (strcmp("51", na) == 0)
					itype = instDTP51;
				else if (strcmp("41", na) == 0)
					itype = instDTP41;
				else if (strcmp("SO", na) == 0 || strcmp("so", na) == 0)
					itype = instSpectrolino;
				else if (strcmp("SS", na) == 0 || strcmp("so", na) == 0)
					itype = instSpectroScan;
				else
					usage();
			}

			/* Spectral Illuminant type */
			else if (argv[fa][1] == 'l' || argv[fa][1] == 'L') {
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

			/* Request transmission measurement */
			else if (argv[fa][1] == 't' || argv[fa][1] == 'T') {
				trans = 1;
				emiss = 0;
			}

			/* Request emissive measurement */
			else if (argv[fa][1] == 'e' || argv[fa][1] == 'E') {
				emiss = 1;
				trans = 0;
			}

			else 
				usage();
		}
		else
			break;
	}


	/* - - - - - - - - - - - - - - - - - - -  */
	/* Setup the instrument ready to do reads */
	if ((it = new_inst(itype)) == NULL) {
		warning("Unknown or no instrument selected");
		usage();
	}
	if (debug)
		it->debug = debug;	/* Turn on debugging */
	if (verb)
		it->verb = verb;	/* Turn on verbosity */

	if (verb)
		printf("Connecting to the instrument ..\n");

#ifdef DEBUG
	printf("About to init the comms\n");
#endif

	/* Establish communications */
	if ((rv = it->init_coms(it, comport, br, 15.0)) != inst_ok) {
		printf("Failed to initialise communications with instrument\n"
		       "or wrong instrument or bad configuration!\n"
		       "('%s' + '%s')\n", it->inst_interp_error(it, rv), it->interp_error(it, rv));
		return -1;
	}

#ifdef DEBUG
	printf("Established comms\n");
#endif

#ifdef DEBUG
	printf("About to init the instrument\n");
#endif

	/* Initialise the instrument */
	if ((rv = it->init_inst(it)) != inst_ok) {
#ifdef DEBUG
		printf("init_inst returned '%s' (%s)\n",
		       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
		printf("Instrument initialisation failed!\n");
		return -1;
	}
	
	/* Configure the instrument mode */
	{
		inst_mode mode = 0;

		cap = it->capabilities(it);

		if (trans) {
			if ((cap & inst_trans_spot) == 0) {
				printf("Need transmission spot capability,\n");
				printf("and instrument doesn't support it\n");
				return -1;
			}
		} else if (emiss) {
			if ((cap & inst_emis_spot) == 0) {
				printf("Need emissive spot capability\n");
				printf("and instrument doesn't support it\n");
				return -1;
			}
		} else {
			if ((cap & inst_ref_spot) == 0) {
				printf("Need reflection spot reading capability,\n");
				printf("and instrument doesn't support it\n");
				return -1;
			}
		}

		if ((spec || pspec) && (cap & inst_spectral) == 0) {
			printf("Need spectral information for custom illuminant or observer\n");
			printf("and instrument doesn't support it\n");
			return -1;
		}

		/* Set it to the appropriate mode */

		/* Should look at instrument type & user spec ??? */
		if (trans)
			mode = inst_mode_trans_spot;
		else if (emiss)
			mode = inst_mode_emis_spot;
		else
			mode = inst_mode_ref_spot;

		if (spec || pspec)
			mode |= inst_mode_spectral;

		if ((rv = it->set_mode(it, mode)) != inst_ok) {
			printf("\nSetting instrument mode failed with error :'%s' (%s)\n",
		     	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
			return -1;
		}

	}

#ifdef DEBUG
	printf("About to enter read loop\n");
#endif

	if (verb)
		printf("Success !\n");

	if (spec) {
		/* Create a spectral conversion object */
		if ((sp2cie = new_xsp2cie(illum, &cust_illum, observ, NULL, icSigXYZData)) == NULL)
			error("Creation of spectral conversion object failed");
	}

	if (verb)
		printf("[Hit R to store reference value while taking a reading]\n");

	/* Hold table */
	if (cap & inst_xy_holdrel) {
		for (;;) {		/* retry loop */
			if ((rv = it->xy_sheet_hold(it)) == inst_ok)
				break;

			if (ierror(it, rv)) {
				it->xy_clear(it);
				return -1;
			}
		}
	}

	/* Read spots until the user quits */
	for (;;) {
		ipatch val;
		double XYZ[3], tXYZ[3];
		int nc;
		int dofwa = 0;		/* Setup for FWA compensation */
		int fidx = -1;		/* FWA compensated index */
	
		/* Do a calibration before the user places the instrument on a desired spot */
		if (it->needs_calibration(it) == inst_needs_cal) {
			for (;;) {		/* retry loop */
				if ((rv = it->calibrate(it, 0)) == inst_ok)
					break;
				if (ierror(it, rv) == 0)	/* Ignore */
					continue;
				break;						/* Abort */
			}
			if (rv != inst_ok)
				break;			/* Abort */
		}

		/* Allow user to locate measurement point, and move to it. */
		if ((cap & inst_xy_locate) && (cap & inst_xy_position)) {
			for (;;) {		/* retry loop */
				if ((rv = it->xy_sheet_hold(it)) == inst_ok)
					break;
				if (ierror(it, rv) == 0)	/* Ignore */
					continue;
				break;						/* Abort */
			}
			if (rv != inst_ok)
				break;			/* Abort */

			for (;;) {		/* retry loop */
				if ((rv = it->xy_locate_start(it)) == inst_ok)
					break;
				if (ierror(it, rv) == 0)	/* Ignore */
					continue;
				break;						/* Abort */
			}
			if (rv != inst_ok)
				break;			/* Abort */

			empty_con_chars();
			printf("\nUsing the XY table controls, locate the point to measure with the sight,\n");

		/* Purely manual instrument */
		} else {

			empty_con_chars();
			printf("Place instrument on spot to be measured,\n");
		}

		printf("and hit [A-Z] to setup FWA compensation (keyed to letter)\n");
		printf("[a-z] to make FWA compensated reading from keyed reference\n");
		printf("'r' to take previous reading as the reference\n");
		printf("Hit ESC or Q to exit, any other key to take a reading: "); fflush(stdout);

		nc = next_con_char();
		if (nc == 0x1b || nc == 0x03 || nc == 'q' || nc == 'Q') {	/* Or ^C */
			printf("\n");
			if ((cap & inst_xy_locate) && (cap & inst_xy_position))
				it->xy_locate_end(it);
			break;
		}
		printf("\n");
		if (nc == 'R' || nc == 'r') {		/* Make last reading the reference */
			if (Lab[0] >= 0.0) {
				rLab[0] = Lab[0];
				rLab[1] = Lab[1];
				rLab[2] = Lab[2];
				printf("Reference is now Lab: %f %f %f\n", rLab[0], rLab[1], rLab[2]);
			} else {
				printf("No previous reading to use as reference\n");
			}
			printf("\n");
			continue;
		}

		if (nc >= 'A' && nc <= 'Z') {
			printf("Measuring media to setup FWA compensation '%c'\n",nc);
			dofwa = 1;
			fidx = nc - 'A';
		}
		if (nc >= 'a' && nc <= 'z') {
			fidx = nc - 'a';
		}

		if ((cap & inst_xy_locate) && (cap & inst_xy_position)) {
			for (;;) {		/* retry loop */
				if ((rv = it->xy_get_location(it, &lx, &ly)) == inst_ok)
					break;
				if (ierror(it, rv) == 0)	/* Ignore */
					continue;
				break;						/* Abort */
			}
			if (rv != inst_ok)
				break;			/* Abort */

			for (;;) {		/* retry loop */
				if ((rv = it->xy_locate_end(it)) == inst_ok)
					break;
				if (ierror(it, rv) == 0)	/* Ignore */
					continue;
				break;						/* Abort */
			}
			if (rv != inst_ok)
				break;			/* Abort */

			for (;;) {		/* retry loop */
				if ((rv = it->xy_position(it, 1, lx, ly)) == inst_ok)
					break;
				if (ierror(it, rv) == 0)	/* Ignore */
					continue;
				break;						/* Abort */
			}
			if (rv != inst_ok)
				break;			/* Abort */
		}

		if ((rv = it->read_sample(it, "SPOT", &val)) != inst_ok) {
#ifdef DEBUG
			printf("read_sample returned '%s' (%s)\n",
		       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */

			/* Deal with a user abort */
			if ((rv & inst_mask) == inst_user_abort) {
				printf("\nSpot read stopped at user request!\n");
				printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
				nc = next_con_char();
				if (nc == 0x1b || nc == 0x03 || nc == 'q' || nc == 'Q') {
					printf("\n");
					break;
				}
				printf("\n");
				continue;
			/* Deal with a misread or needs calibration */
			} else if ((rv & inst_mask) == inst_misread || (rv & inst_mask) == inst_needs_cal) {
				empty_con_chars();
				if ((rv & inst_mask) == inst_needs_cal) {
					printf("\nSpot read failed because instruments needs calibration.\n");

					if (it->capabilities(it) & inst_calib) {	/* Host controlled */
						printf("Hit Esc to give up, or place instrument on calibration tile\n");
						printf("and hit any other key:"); fflush(stdout);
						nc = next_con_char();
						if (nc == 0x1b || nc == 0x03 || nc == 'q' || nc == 'Q') {
							printf("\n");
							break;
						}
						printf("\n");
						if ((rv = it->calibrate(it, 0)) != inst_ok) {
							printf("Calibrate failed with '%s' (%s)\n",
							       it->inst_interp_error(it, rv), it->interp_error(it, rv));
							break;
						} else {
							printf("calibration successful\n");
						}
						continue;
					} else {	/* User controlled */
						printf("Hit Esc to give up, or do calibration and then hit any\n");
						printf("other key to retry:"); fflush(stdout);
					}
				} else {
					printf("\nSpot read failed due to misread\n");
					printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
				}
				nc = next_con_char();
				if (nc == 0x1b || nc == 0x03 || nc == 'q' || nc == 'Q') {
					printf("\n");
					break;
				}
				printf("\n");
				continue;
			/* Deal with a communications error */
			} else if ((rv & inst_mask) == inst_coms_fail) {
				int tt = it->last_sioerr(it);
				if (tt & (SIO_BRK | SIO_FER | SIO_PER | SIO_OER)) {
					if (br == baud_57600) br = baud_38400;
					else if (br == baud_38400) br = baud_9600;
					else if (br == baud_9600) br = baud_4800;
					else if (br == baud_9600) br = baud_4800;
					else if (br == baud_2400) br = baud_1200;
					else br = baud_1200;
				}
				/* Communication problem, allow retrying at a lower baud rate */
				empty_con_chars();
				printf("\nSpot read failed due to communication problem.\n");
				printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
				nc = next_con_char();
				if (nc == 0x1b || nc == 0x03 || nc == 'q' || nc == 'Q') {
					printf("\n");
					break;
				}
				printf("\n");
				if ((rv = it->init_coms(it, comport, br, 15.0)) != inst_ok) {
#ifdef DEBUG
					printf("init_coms returned '%s' (%s)\n",
				       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
					break;
				}
				continue;

			/* Some other fatal error */
			} else {
				printf("Got fatal error '%s'\n", it->inst_interp_error(it, rv));
				break;
			}
		}

		if ((cap & inst_xy_locate) && (cap & inst_xy_position)) {
			for (;;) {		/* retry loop */
				if ((rv = it->xy_position(it, 0, lx, ly)) == inst_ok)
					break;
				if (ierror(it, rv) == 0)	/* Ignore */
					continue;
				break;						/* Abort */
			}
			if (rv != inst_ok)
				break;			/* Abort */
		}

#ifdef DEBUG
		printf("read_sample returned '%s' (%s)\n",
	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */


		/* We got a valid spot reading */

		if (pspec) {	/* Print out spectrum */

			if (val.spec_n <= 0)
				error("Instrument didn't return spectral data");

			printf("Spectrum from %f yp %f nm in %d steps\n",
		                val.spec_wl_short, val.spec_wl_long, val.spec_n);

			for (j = 0; j < val.spec_n; j++)
				printf("%s%8.3f",j > 0 ? ", " : "", val.spec[j]);
			printf("\n");
		}

		if (!spec) {
			if (val.XYZ_v == 0 && val.aXYZ_v == 0)
				error("Instrument didn't return XYZ value");
			
			if (val.XYZ_v != 0) {
				for (j = 0; j < 3; j++)
					XYZ[j] = val.XYZ[j];
			} else {
				// Hmm. Should we really allow an absolute emmisive and relative mode ? */
				for (j = 0; j < 3; j++)
					XYZ[j] = val.aXYZ[j];
			}

			/* Compute D50 Lab from XYZ */
			tXYZ[0] = XYZ[0]/100.0;			/* Scale to Y = 1.0 */
			tXYZ[1] = XYZ[1]/100.0;
			tXYZ[2] = XYZ[2]/100.0;
			icmXYZ2Lab(&icmD50, Lab, tXYZ);
	
		} else {
			/* Compute XYZ given spectral */

			if (val.spec_n <= 0) {
				error("Instrument didn't return spectral data");
			}

			/* Transfer from vals to sp structure */
			sp.spec_n = val.spec_n;
			sp.spec_wl_short = val.spec_wl_short;
			sp.spec_wl_long  = val.spec_wl_long;
			sp.norm  = 100.0;
			
			for (j = 0; j < val.spec_n; j++)
				sp.spec[j] = val.spec[j];

			if (dofwa) {	/* Setup FWA compensation */
				xspect *insp;			/* Instrument illuminant */

				if ((insp = inst_illuminant(itype)) == NULL)
					error ("Instrument doesn't have an FWA illuminent");

				/* Creat the base conversion object */
				if (sp2cief[fidx] == NULL) {
					if ((sp2cief[fidx] = new_xsp2cie(illum, &cust_illum, observ,
					                                 NULL, icSigXYZData)) == NULL)
						error("Creation of spectral conversion object failed");
				}

				if (sp2cief[fidx]->set_fwa(sp2cief[fidx], insp, &sp)) 
					error ("Set FWA on sp2cie failed");

			} else {

				if (fidx < 0 || sp2cief[fidx] == NULL) {
					/* Convert it to XYZ space using uncompensated */
					sp2cie->convert(sp2cie, XYZ, &sp);
				} else {
					/* Convert using compensated conversion */
					sp2cief[fidx]->convert(sp2cief[fidx], XYZ, &sp);
				}
			}

			/* Compute D50 Lab from XYZ */
			icmXYZ2Lab(&icmD50, Lab, XYZ);

			XYZ[0] *= 100.0;
			XYZ[1] *= 100.0;
			XYZ[2] *= 100.0;
		}

		if (dofwa == 0) {
			printf("Result is XYZ: %f %f %f, D50 Lab: %f %f %f\n",
			XYZ[0], XYZ[1], XYZ[2], Lab[0], Lab[1], Lab[2]);

			if (rLab[0] >= -1.0) {
				printf("Delta E to reference is %f %f %f (%f, CIE94 %f)\n",
				Lab[0] - rLab[0], Lab[1] - rLab[1], Lab[2] - rLab[2],
				icmLabDE(Lab, rLab), icmCIE94(Lab, rLab));
			}
		} 
		printf("\n");
	}

	/* Release paper */
	if (cap & inst_xy_holdrel) {
		it->xy_clear(it);
	}

#ifdef DEBUG
	printf("About to exit\n");
#endif

	if (sp2cie != NULL)
		sp2cie->del(sp2cie);
	for (i = 0; i < 26; i++)
		if (sp2cief[i] != NULL)
			sp2cief[i]->del(sp2cief[i]);

	it->del(it);

	return 0;
}



