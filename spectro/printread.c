
/* 
 * Argyll Color Correction System
 * Spectrometer/Colorimeter print target chart reader
 *
 * Author: Graeme W. Gill
 * Date:   4/10/96
 *
 * Copyright 1996 - 2005, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* This program reads a reflective or transmissive print target chart */
/* using a spectrometer or colorimeter. */

/*
 * Nomencalture:
 *
 *	Largely due to how the strip readers name things, the following terms
 *  are used for how patches are grouped:
 *
 *  Step:   One test patch in a pass, usually labelled with a number.
 *  Pass:	One row of patches in a strip. A pass is usually labeled
 *          with a unique alphabetic label. 
 *  Strip:  A group of passes that can be read by a strip reader.
 *          For an XY instrument, the strip is a complete sheet, and
 *          a each pass is one column. The rows of an XY chart are
 *          the step numbers within a pass.
 *  Sheet:  One sheet of paper, containing full and partial strips.
 *          For an XY instrument, there will be only one strip per sheet.
 *           
 */

/* TTBD
 *
 */
 
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
#include "cgats.h"
#include "serio.h"
#include "numlib.h"
#include "xspect.h"
#include "xcolorants.h"
#include "inst.h"
#include "alphix.h"
#include "sort.h"
#include "icc.h"

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

/* Return the normal Delta E, given two 100 scaled XYZ values */
static double xyzLabDE(double *in1, double *in2) {
	int i;
	double Lab1[3];
	double Lab2[3];

	for (i = 0; i < 3; i++) {
		Lab1[i] = in1[i]/100.0;
		Lab2[i] = in2[i]/100.0;
	}

	icmXYZ2Lab(&icmD50, Lab1, Lab1);
	icmXYZ2Lab(&icmD50, Lab2, Lab2);

	return icmLabDE(Lab1, Lab2);
}

/* A color structure */
/* This can hold all representations simultaniously */
typedef struct {
	char *id;					/* Id string (e.g. "1") */
	char *loc;					/* Location string (e.g. "A1") */
	int loci;					/* Location integer = pass * 256 + step */

	int n;						/* Number of colorants */
	double dev[ICX_MXINKS];		/* Value of colorants */
	double eXYZ[3];				/* Expected XYZ values (100.0 scale) */

	int rr;						/* nz if reading read (used for tracking manual reading) */
	double XYZ[3];				/* Colorimeter readings (100.0 scale) */

	int    spec_n;					/* Number of spectral bands, 0 if not valid */
	double spec_wl_short;			/* First reading wavelength in nm (shortest) */
	double spec_wl_long;			/* Last reading wavelength in nm (longest) */
	double spec[INSTR_MAX_BANDS];	/* Spectral reflectance %, shortest to longest */
	
} col;

/* Convert a base 62 character into a number */
/* (This is used for converting the PASSES_IN_STRIPS string */
/* (Could convert this to using an alphix("0-9A-Za-Z")) */
static int b62_int(char *p) {
	int cv, rv;

	cv = *p;
	if (cv == '\000')
		rv = 0;
	else if (cv <= '9')
		rv = cv - '0';
	else if (cv <= 'Z')
		rv = cv - 'A' + 10;
	else
		rv = cv - 'a' + 36;
	return rv;
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

/* Read all the strips, and return nonzero on abort/error */
int
read_strips(
instType itype,		/* Instrument type */
col **scols,		/* Location sorted pointers to cols (return values) */
int npat,			/* Total valid patches */
int totpa, 			/* Total passes */
int stipa,			/* Steps in Pass (Excluding DTP51 Max/Min) */
char *paists, 		/* Passs in each Strip string, base 62 */
alphix *paix,		/* Pass index generators */
alphix *saix,		/* Step index generators */
int ixord,			/* Index order, 0 = pass then step */
int hex,			/* Hexagon test patches */
int comport, 		/* COM port used */
double plen,		/* Patch length in mm (used by DTP41) */
double glen,		/* Gap length in mm  (used by DTP41) */
double tlen,		/* Trailer length in mm  (used by DTP41T) */
int trans,			/* Use transmissioin mode */
int emiss,			/* Use emissive mode */
int xtern,			/* Use external values */
int spectral,		/* Generate spectral info flag */
int verb,			/* Verbosity flag */
int debug			/* Debug flag */
) {
	inst *it = NULL;
	inst_capability cap;
	int n, i, j;
	int rmode = 0;		/* Read mode, 0 = spot, 1 = strip, 2 = xy */
	inst_code rv;
	char buf[100];
	baud_rate br = baud_38400;	/* Target baud rate */
	int skipp = 0;		/* Initial strip readings to skip */
	int	nextrap = 0;	/* Number of extra patches for max and min */

	if (xtern == 0) {

		if (itype == instDTP51) {
			skipp = 1;		/* First reading is the Max density patch */
			nextrap = 2;
		}

		it = new_inst(itype);
		if (debug)
			it->debug = debug;	/* Turn on debugging */
		if (verb)
			it->verb = verb;	/* Turn on verbosity */

		/* Establish communications */
		if ((rv = it->init_coms(it, comport, br, 15.0)) != inst_ok) {
			printf("Establishing communications with instrument failed with message '%s' (%s)\n",
			       it->inst_interp_error(it, rv), it->interp_error(it, rv));
			return -1;
		}

		/* Set it up the way we want */
		if ((rv = it->init_inst(it)) != inst_ok) {
			printf("Initialising instrument failed with message '%s' (%s)\n",
			       it->inst_interp_error(it, rv), it->interp_error(it, rv));
			return -1;
		}

		{
			inst_mode mode = 0;

			cap = it->capabilities(it);

			if (trans) {
				if ((cap & (inst_trans_strip | inst_trans_xy | inst_trans_spot)) == 0) {
					printf("Need transmission spot, strip or xy reading capability,\n");
					printf("and instrument doesn't support it\n");
					return -1;
				}
			} else if (emiss) {
				if ((cap & (inst_emis_spot)) == 0) {
					printf("Need emissive spot reading capability,\n");
					printf("and instrument doesn't support it\n");
					return -1;
				}
			} else {
				if ((cap & (inst_ref_strip | inst_ref_xy | inst_ref_spot)) == 0) {
					printf("Need reflection spot, strip or xy reading capability,\n");
					printf("and instrument doesn't support it\n");
					return -1;
				}
			}

			if (spectral && (cap & inst_spectral) == 0) {
				printf("Instrument isn't capable of spectral measurement\n");
				spectral = 0;
			}

			/* Set it to the appropriate mode */

			/* Should look at instrument type & user spec ??? */
			if (trans) {
				if (cap & inst_trans_xy) {
					mode = inst_mode_trans_xy;
					rmode = 2;
				} else if (cap & inst_trans_strip) {
					mode = inst_mode_trans_strip;
					rmode = 1;
				} else {
					mode = inst_mode_trans_spot;
					rmode = 0;
				}
			} else if (emiss) {
				mode = inst_mode_emis_spot;
				rmode = 0;
			} else {
				if (cap & inst_ref_xy) {
					mode = inst_mode_ref_xy;
					rmode = 2;
				} else if (cap & inst_ref_strip) {
					mode = inst_mode_ref_strip;
					rmode = 1;
				} else {
					mode = inst_mode_ref_spot;
					rmode = 0;
				}
			}
			if (spectral)
				mode |= inst_mode_spectral;

			if ((rv = it->set_mode(it, mode)) != inst_ok) {
				printf("\nSetting instrument mode failed with error :'%s' (%s)\n",
			     	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
				return -1;
			}
		}
	}

	/* -------------------------------------------------- */
	if (rmode == 2) { /* For xy mode, read each sheet */
		ipatch *vals;
		int nsheets, sheet;		/* Total sheets/current sheet (sheet == pass) */
		int rpat = npat;		/* Remaining total patches */
		int pai;				/* Overall pass index */
		int sti;				/* Overall step index */
		char *pn[3] = { NULL, NULL, NULL} ;	/* Location 1/2/3 Pass name (letters) */
		char *sn[3] = { NULL, NULL, NULL} ;	/* Location 1/2/3 Step name (numbers) */
		int k;
		
		{	/* Figure the maximum number sheets and of patches in a sheet, for allocation */
			int i, lpaist = 0;			/* Largest passes in strip */
			char *p;
			for (nsheets = 0, p = paists; *p != '\000'; p++, nsheets++) {
				i = b62_int(p);
				if (i > lpaist)
					lpaist = i;
			}

			if ((vals = (ipatch *)calloc(sizeof(ipatch), (lpaist * stipa))) == NULL)
				error("Malloc failed!");
		}

		/* Make sure we can access the instrument table */
		if (cap & inst_xy_holdrel) {
			it->xy_clear(it);
		}

		/* For each pass (==sheet) */
		for (sheet = 1, pai = sti = 0; *paists != '\000'; paists++, sheet++) {
			int paist;		/* Passes in current Strip (== columns in current sheet) */ 
			int rnpatch;	/* Rounded up (inc. padding) Patches in current pass (sheet) */
			int npatch;		/* Patches in pass (sheet), excluding padding */
			int fspipa;		/* First pass steps in pass */
			int nloc;		/* Number of fiducial locations needed */
			double ox = 0.0, oy = 0.0;		/* Origin */
			double ax = 1.0, ay = 0.0;		/* pass increment */
			double aax = 0.0, aay = 0.0;	/* pass offset for hex odd steps */
			double px = 0.0, py = 1.0;		/* step (==patch) increment */

			fspipa = stipa;
			paist = b62_int(paists);			/* columns (letters) in sheet (strip) */
			npatch = rnpatch = paist * stipa;	/* Total patches including padding */
			if (npatch > rpat) {				/* This is a non-full pass */
				if (paist == 1) {
					fspipa -= (npatch - rpat);/* Last patch in first strip */
					if (fspipa < 1)
						error ("Assert in read_strips, fspipa = %d",fspipa);
				}
				npatch = rpat;				/* Total patches excluding padding */
			}

			nloc = 3;
			if (paist == 1) {
				nloc = 2;					/* Only one strip, so only 2 locations needed */
				if (fspipa == 1)
					nloc = 1;				/* Only one strip, one patch, so one location */
			}
			for (k = 0; k < 3; k++) {
				if (pn[k] != NULL)
					free(pn[k]);
				if (sn[k] != NULL)
					free(sn[k]);
			}
			pn[0] = paix->aix(paix, pai);			/* First pass (letter) */
			sn[0] = saix->aix(saix, 0);				/* First step (patch) (number) */
			pn[1] = paix->aix(paix, pai);			/* First pass (letter) */
			sn[1] = saix->aix(saix, 0 + fspipa-1);	/* Last step (patch) (number) */
			pn[2] = paix->aix(paix, pai + paist-1);	/* Last pass (letter) */
			sn[2] = saix->aix(saix, 0);				/* First step (patch) (number) */

			empty_con_chars();
			if (sheet == 1) {
				printf("Please place sheet %d of %d on table, then\n",sheet, nsheets);
			} else
				printf("\nPlease remove previous sheet, then place sheet %d of %d on table, then\n",sheet, nsheets);
			printf("hit return to continue, Esc to give up"); fflush(stdout);
				if (next_con_char() == 0x1b) {
					printf("\n");
					for (k = 0; k < 3; k++) {
						if (pn[k] != NULL)
							free(pn[k]);
						if (sn[k] != NULL)
							free(sn[k]);
					}
					it->xy_clear(it);
					return -1;
				}
			printf("\n");

			if (cap & inst_xy_holdrel) {

				/* Hold table */
				for (;;) {		/* retry loop */
					if ((rv = it->xy_sheet_hold(it)) == inst_ok)
						break;
					if (ierror(it, rv) == 0)	/* Ignore */
						continue;
					break;						/* Abort */
				}
				if (rv != inst_ok) {
					for (k = 0; k < 3; k++) {
						if (pn[k] != NULL)
							free(pn[k]);
						if (sn[k] != NULL)
							free(sn[k]);
					}
					it->xy_clear(it);
					return -1;
				}
			}

			if (cap & inst_xy_locate) {
				int ll;
				double x[3], y[3];

				/* Allow user location of points */
				for (;;) {		/* retry loop */
					if ((rv = it->xy_locate_start(it)) == inst_ok)
						break;
					if (ierror(it, rv) == 0)	/* Ignore */
						continue;
					break;						/* Abort */
				}
				if (rv != inst_ok) {
					for (k = 0; k < 3; k++) {
						if (pn[k] != NULL)
							free(pn[k]);
						if (sn[k] != NULL)
							free(sn[k]);
					}
					it->xy_clear(it);
					return -1;
				}

				/* For each location point */
				for (ll = 0; ll < nloc; ll++) {

					empty_con_chars();
					printf("\nUsing the XY table controls, locate patch %s%s with the sight,\n",
					        pn[ll], sn[ll]);
					printf("then hit return to continue, Esc to give up"); fflush(stdout);
						if (next_con_char() == 0x1b) {
							printf("\n");
							for (k = 0; k < 3; k++) {
								if (pn[k] != NULL)
									free(pn[k]);
								if (sn[k] != NULL)
									free(sn[k]);
							}
							it->xy_clear(it);
							return -1;
						}
					printf("\n");

					for (;;) {		/* retry loop */
						if ((rv = it->xy_get_location(it, &x[ll], &y[ll])) == inst_ok)
							break;
						if (ierror(it, rv) == 0)	/* Ignore */
							continue;
						break;						/* Abort */
					}
					if (rv != inst_ok) {
						for (k = 0; k < 3; k++) {
							if (pn[k] != NULL)
								free(pn[k]);
							if (sn[k] != NULL)
								free(sn[k]);
						}
						it->xy_clear(it);
						return -1;
					}
				}

				/* We're done with user control */
				for (;;) {		/* retry loop */
					if ((rv = it->xy_locate_end(it)) == inst_ok)
						break;
					if (ierror(it, rv) == 0)	/* Ignore */
						continue;
					break;						/* Abort */
				}
				if (rv != inst_ok) {
					for (k = 0; k < 3; k++) {
						if (pn[k] != NULL)
							free(pn[k]);
						if (sn[k] != NULL)
							free(sn[k]);
					}
					it->xy_clear(it);
					return -1;
				}

				/* Convert point locations into navigation values */
				ox = x[0];
				oy = y[0];
				if (hex) {
					double kk = sqrt(1.0/3.0);
					double nn = fspipa - 1.0;
					if (((fspipa-1) & 1) == 0) {	/* [0] & [1] are lined up */
						if (nloc == 3) {
							px = (x[1] - x[0])/nn;
							py = (y[1] - y[0])/nn;
							ax = (x[2] - x[0])/(paist-1);
							ay = (y[2] - y[0])/(paist-1);
							aax = 0.5 * ax;
							aay = 0.5 * ay;
						} else if (nloc == 2) {
							px = (x[1] - x[0])/nn;
							py = (y[1] - y[0])/nn;
							aax = kk * py;			/* Scale and rotate */
							aay = kk * -px;
						}
					} else {						/* [0] & [1] are offset by aa[xy] */
						if (nloc == 3) {
							ax = (x[2] - x[0])/(paist-1);
							ay = (y[2] - y[0])/(paist-1);
							aax = 0.5 * ax;
							aay = 0.5 * ay;
							px = (x[1] - x[0] - aax)/nn;
							py = (y[1] - y[0] - aay)/nn;
						} else if (nloc == 2) {
							px =  (nn * (x[1] - x[0]) - kk * (y[1] - y[0]))/(kk * kk + nn * nn);
							py =  (nn * (y[1] - y[0]) + kk * (x[1] - x[0]))/(kk * kk + nn * nn);
							aax = kk * py;			/* Scale and rotate */
							aay = kk * -px;
						}
					}

				} else {	/* Rectangular patches */
					if (paist > 1) {
						ax = (x[2] - x[0])/(paist-1);
						ay = (y[2] - y[0])/(paist-1);
					}
					if (fspipa > 1) {
						px = (x[1] - x[0])/(fspipa-1);
						py = (y[1] - y[0])/(fspipa-1);
					}
					aax = aay = 0.0;
				}
			}

			/* Read the sheets patches */
			for (;;) {		/* retry loop */
				if ((rv = it->read_xy(it, paist, stipa, npatch, pn[0], sn[0],
					                  ox, oy, ax, ay, aax, aay, px, py, vals)) == inst_ok)
					break;
				if (ierror(it, rv) == 0)	/* Ignore */
					continue;
				break;						/* Abort */
			}
			if (rv != inst_ok) {
				for (k = 0; k < 3; k++) {
					if (pn[k] != NULL)
						free(pn[k]);
					if (sn[k] != NULL)
						free(sn[k]);
				}
				it->xy_clear(it);
				return -1;
			}

			printf("Sheet %d of %d read OK\n",sheet, nsheets);

			/* Transfer the values */
			/* We assume they are all in the right order */
			for (i = 0; i < npatch; i++, sti++) {
				/* Copy XYZ */
				if (vals[i].XYZ_v == 0)
					error("Instrument didn't return XYZ value for patch %d, loc %s",i,scols[sti]->loc);
				for (j = 0; j < 3; j++)
					scols[sti]->XYZ[j] = vals[i].XYZ[j];

				/* Copy spectral */
				if (vals[i].spec_n > 0) {
					scols[sti]->spec_n = vals[i].spec_n;
					scols[sti]->spec_wl_short = vals[i].spec_wl_short;
					scols[sti]->spec_wl_long  = vals[i].spec_wl_long;
					
					for (j = 0; j < vals[i].spec_n; j++)
						scols[sti]->spec[j] = vals[i].spec[j];
				}
				scols[sti]->rr = 1;		/* Has been read */
			}

			if (cap & inst_xy_holdrel) {

				/* Release table and reset head */
				for (;;) {		/* retry loop */
					if ((rv = it->xy_clear(it)) == inst_ok)
						break;
					if (ierror(it, rv) == 0)	/* Ignore */
						continue;
					break;						/* Abort */
				}
				if (rv != inst_ok) {
					for (k = 0; k < 3; k++) {
						if (pn[k] != NULL)
							free(pn[k]);
						if (sn[k] != NULL)
							free(sn[k]);
					}
					return -1;
				}
			}

			pai += paist;		/* Tracj next first pass in strip */
			rpat -= npatch;		/* Track remaining patches */
		}
		for (k = 0; k < 3; k++) {
			if (pn[k] != NULL)
				free(pn[k]);
			if (sn[k] != NULL)
				free(sn[k]);
		}
		free(vals);

		printf("\nPlease remove last sheet from table\n"); fflush(stdout);

	/* -------------------------------------------------- */
	} else if (rmode == 1) {	/* For strip mode, simply read each strip */
		ipatch *vals;		/* Values read for a strip pass */
		int opai;			/* Overall pass index */

		/* Allocate space for values from a pass */
		if ((vals = (ipatch *)calloc(sizeof(ipatch), (stipa+nextrap))) == NULL)
			error("Malloc failed!");

		/* For all passes */
		for (opai = 0; *paists != '\000'; paists++) {
			int paist;		/* Passes in current strip */
			int pai;		/* Current pass in current strip */

			/* For all strips in this pass */
			for (pai = 0, paist = b62_int(paists); pai < paist; pai++, opai++) {
				char *nn = NULL;	/* Pass name */
				int guide;
				col **scb;
				int boff = 0;			/* Best offset to read from (DTP51) */

				/* Convert overall pass number index into alpha label */
				if (nn != NULL)
					free(nn);
				nn = paix->aix(paix, opai);

				guide = (paist - pai) * 5;		/* Mechanical guide offset */

				for (;;) {	/* Until we give up retrying */
					/* Read a strip pass */
					printf("About to read strip pass %s : ",nn); fflush(stdout);
					if ((rv = it->read_strip(it, "STRIP", stipa+nextrap, nn, guide, plen, glen, tlen, vals)) != inst_ok) {
#ifdef DEBUG
						printf("read_strip returned '%s' (%s)\n",
					       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
						/* Deal with a user abort */
						if ((rv & inst_mask) == inst_user_abort) {
							empty_con_chars();
							printf("\nStrip read stopped at user request!\n");
							printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
							if (next_con_char() == 0x1b) {
								printf("\n");
								return -1;
							}
							printf("\n");
						continue;
						/* Deal with a misread or needs calibration */
						} else if ((rv & inst_mask) == inst_misread || (rv & inst_mask) == inst_needs_cal) {
							empty_con_chars();
							if ((rv & inst_mask) == inst_needs_cal) {
								printf("\nStrip read failed because instruments needs calibration\n");
								printf("Hit Esc to give up, or do calibration and then hit any\n");
								printf("other key to retry:"); fflush(stdout);
							} else {
								printf("\nStrip read failed due to misread\n");
								printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
							}
							if (next_con_char() == 0x1b) {
								printf("\n");
								return -1;
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
							printf("\nStrip read failed due to communication problem.\n");
							printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
							if (next_con_char() == 0x1b) {
								printf("\n");
								return -1;
							}
							printf("\n");
							if ((rv = it->init_coms(it, comport, br, 15.0)) != inst_ok) {
#ifdef DEBUG
								printf("init_coms returned '%s' (%s)\n",
							       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
								return -1;
							}
							continue;
						} else {
							/* Some other error. Treat it as fatal */
							printf("\nStrip read failed due unexpected error :'%s' (%s)\n",
					       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
							printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
							if (next_con_char() == 0x1b) {
								printf("\n");
								return -1;
							}
							printf("\n");
							continue;
						}

					/* Successfully read the strip */
					/* See which expected pass best correlates with the one we've read. */
					/* Figure out if there is an "off by one" error for a DTP51 */
					} else {
						int chopai;			/* Check overall pass index */
						double corr;		/* Correlation with expected value */
						double werror;		/* Worst case error in best correlation strip */
						double xwerror;		/* Expected strips worst case error in best correlation strip */
						double bcorr;		/* Best correlation value */
						int bopai;			/* Best overall pass index */
						int loff = 0, hoff = 0;	/* DTP51 offset test range */
						int toff;			/* Test offset */

						/* DTP51 has a nasty habit of misaligning test squares by +/- 1 */
						/* See if this might have happened */
						if (itype == instDTP51) {
							loff = -1;
							hoff = 1;
						}

						bcorr = 1e6;
						bopai = -1;
						for (chopai = 0; chopai < totpa; chopai++) {
							double pwerr;	/* This passes worst error */
							scb = &scols[chopai * stipa];

							/* Explore off by +/-1 error for DTP51 */
							for (toff = loff; toff <= hoff; toff++) {

								/* Compare just sample patches (not padding Max/Min) */
								for (pwerr = corr = 0.0, n = 0, i = 0; i < stipa; i++, n++) {
									double vcorr;
									if (vals[i+skipp+toff].XYZ_v == 0)
										error("Instrument didn't return XYZ value");
									vcorr = xyzLabDE(vals[i+skipp+toff].XYZ, scb[i]->eXYZ);
									corr += vcorr;
									if (vcorr > pwerr)
										pwerr = vcorr;
								}
								corr /= (double)n;
#ifdef DEBUG
								printf("  Strip %d offset %d correlation = %f\n",chopai,off,corr);

#endif
								if (chopai == opai)
									xwerror = werror;	/* Expected passes worst error */

								if (corr < bcorr) {
									bcorr = corr;
									bopai = chopai;
									boff = toff;
									werror = pwerr;
								}
							}
						}
						if (bopai != opai) {	/* Looks like the wrong strip */
							char *mm = NULL, ch;
							mm = paix->aix(paix, bopai);
#ifdef DEBUG
							printf("Strip pass %s (%d) seems to have a better correlation that strip %s (%d)\n",
							mm, bopai, nn, opai);
#endif
							empty_con_chars();
							printf("\nSeem to have read strip pass %s rather than %s!\n",mm,nn);
							printf("Hit Return to use it anyway, any other key to retry:"); fflush(stdout);
							free(mm);
							if ((ch = next_con_char()) == 0x1b) {	/* Esc */
								printf("\n");
								return -1;
							}
							if (ch != 0x0d && ch != 0x0a) { 			/* !Cr */
								printf("\n");
								continue;
							}
							printf("\n");
							werror = xwerror;
						}
						if (werror >= 80.0) { /* Arbitrary threshold. Good seems about 15-35, bad 95-130 */
							char ch;
#ifdef DEBUG
							printf("Patch error %f (>35 not good, >95 bad)\n",werror);
#endif
							empty_con_chars();
							printf("\nThere is at least one patch with an very unexpected response! (DeltaE %f)\n",werror);
							printf("Hit Return to use it anyway, any other key to retry:"); fflush(stdout);
							if ((ch = next_con_char()) == 0x1b) {	/* Esc */
								printf("\n");
								return -1;
							}
							if (ch != 0x0d && ch != 0x0a) { 			/* !Cr */
								printf("\n");
								continue;
							}
							printf("\n");
							break;
						}

						/* Must be OK - save the readings */
						if (boff != 0)
							printf(" (DTP51 offset fix of %d applied)",boff);
						printf(" Strip read OK\n");
						break;		/* Break out of retry loop */
					}
				}

				if (nn != NULL)		/* Finished with strip alpha index */
					free(nn);
				nn = NULL;

				/* Transfer the values (including DTP51 offset) */
				scb = &scols[opai * stipa];
				for (n = 0, i = 0; i < stipa; i++, n++) {
					/* Copy XYZ */
					if (vals[i+skipp+boff].XYZ_v == 0)
						error("Instrument didn't return XYZ value");
					for (j = 0; j < 3; j++)
						scb[i]->XYZ[j] = vals[i+skipp+boff].XYZ[j];

					/* Copy spectral */
					if (vals[i+skipp+boff].spec_n > 0) {
						scb[i]->spec_n        = vals[i+skipp+boff].spec_n;
						scb[i]->spec_wl_short = vals[i+skipp+boff].spec_wl_short;
						scb[i]->spec_wl_long  = vals[i+skipp+boff].spec_wl_long;
						
						for (j = 0; j < vals[i+skipp+boff].spec_n; j++)
							scb[i]->spec[j] = vals[i+skipp+boff].spec[j];
					}
					scb[i]->rr = 1;		/* Has been read */
				}
			}
		}
		free(vals);

	/* -------------------------------------------------- */
	/* Spot mode. This will be used if xtern != 0 */
	} else {
		int pix = -1;
		ipatch val;

		/* Until we're  done */
		for(;;) {
			char buf[200], *bp;
			char ch;

			/* Advance to the next patch */
			if (++pix >= npat)
				pix = 0;
			else if (pix < 0)
				pix = npat-1;

			/* See if there aren any unread patches */
			for (i = 0; i < npat; i++)
				if (scols[i]->rr == 0)
					break;

			if (xtern != 0) {
				printf("\nReady to read patch '%s'%s\n",scols[pix]->loc,
				       i >= npat ? "(All patches read!)" : scols[pix]->rr ? " (Already read)" : "");
				printf("Enter %s values or  'f' to move forward, 'b' to move back,\n"
				       " 'd' when done, 'x' to abort, then press <return>: ",
				       xtern == 1 ? "L*a*b*" : "XYZ");
				fflush(stdout);
	
				/* Read in the next line from stdin. */
				if (fgets(buf, 200, stdin) == NULL) {
					printf("Error - unrecognised input\n");
					pix--;		/* Retry */
					continue;
				}
				/* Skip whitespace */
				for (bp = buf; *bp != '\000' && isspace(*bp); bp++)
					;

				ch = *bp;
				if (ch == '\000') {
					printf("Error - unrecognised input\n");
					pix--;		/* Retry */
					continue;
				}
			} else {
				empty_con_chars();
				printf("\nReady to read patch '%s'%s\n",scols[pix]->loc,
				       i >= npat ? "(All patches read!)" : scols[pix]->rr ? " (Already read)" : "");
				printf("hit 'f' to move forward, 'b' to move back,\n"
				       " 'd' when done, <return> or <space> to read, <esc> to abort:"); fflush(stdout);
	
				ch = next_con_char();
				printf("\n");
			}		
	
			if (ch == 'x' || ch == 0x1b || ch == 0x03) {	/* x or Esc or ^C */
				empty_con_chars();
				printf("\nAbort ? - Are you sure ? [y/n]:"); fflush(stdout);
				if ((ch = next_con_char()) == 'y' || ch == 'Y') {
					printf("\n");
					return -1;
				}
				printf("\n");
				pix--;
				continue;
			} else if (ch == 'f') {
				continue;
			} else if (ch == 'F') {
				pix += 9;
				continue;
			} else if (ch == 'b') {
				pix -= 2;
				continue;
			} else if (ch == 'B') {
				pix -= 11;
				continue;
			} else if (ch == 'd' || ch == 'D' || ch == 'q' || ch == 'Q') {
				int i;
				for (i = 0; i < npat; i++)
					if (scols[i]->rr == 0)
						break;
				if (i >= npat)
					break;			/* None unread, so done */

				/* Not all patches have been read */
				empty_con_chars();
				printf("\nDone ? - At least one unread patch (%s), Are you sure [y/n]: ",
				       scols[i]->loc);
				fflush(stdout);
				if ((ch = next_con_char()) == 0x1b) {
					printf("\n");
					return -1;
				}
				printf("\n");
				if (ch == 'y' || ch == 'Y')
					break;
				pix = i-1;
				continue;

			/* Read the external sample */
	        } else if (xtern != 0 && (isdigit(*bp) || ch == '-' || ch == '+' || ch == '.')) {
				int i;
				double vals[3];

				/* For each input number */
				for (i = 0; *bp != '\000' && i < 3; i++) {
					char *tp, *nbp;

					/* Find the start of the number */
					while(*bp != '\000' && !isdigit(*bp)
					      && *bp != '-' && *bp != '+' && *bp != '.')
						bp++;
					if (!isdigit(*bp) && *bp != '-' && *bp != '+' && *bp != '.')
						break;

					/* Find the end of the number */
					for (tp = bp+1; isdigit(*tp) || *tp == 'e' || *tp == 'E'
					                || *tp == '-' || *tp == '+' || *tp == '.'; tp++)
						;
					if (*tp != '\000')
						nbp = tp+1;
					else
						nbp = tp;
					*tp = '\000';

					/* Read the number */
					scols[pix]->XYZ[i] = atof(bp);

					bp = nbp;
				}
				if (i < 3) {	/* Didn't find 3 numbers */
					printf("Error - unrecognised input\n");
					pix--;		/* Retry */
					continue;
				}
				if (xtern == 1)
					icmLab2XYZ(&icmD50, scols[pix]->XYZ,scols[pix]->XYZ);

				scols[pix]->rr = 1;		/* Has been read */
				/* Advance to next patch. */
				
			/* Read the spot sample */
			} else if (xtern == 0 && (ch == '\r' || ch == '\n' || ch == ' ')) {	
				if ((rv = it->read_sample(it, "SPOT", &val)) != inst_ok) {
#ifdef DEBUG
					printf("read_sample returned '%s' (%s)\n",
				       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
					/* Deal with a user abort */
					if ((rv & inst_mask) == inst_user_abort) {
						empty_con_chars();
						printf("\nSpot read stopped at user request!\n");
						printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
						if (next_con_char() == 0x1b) {
							printf("\n");
							return -1;
						}
						printf("\n");
						pix--;
						continue;
					/* Deal with a misread or needs calibration */
					} else if ((rv & inst_mask) == inst_misread || (rv & inst_mask) == inst_needs_cal) {
						empty_con_chars();
						if ((rv & inst_mask) == inst_needs_cal) {
							printf("\nStrip read failed because instruments needs calibration\n");
							printf("Hit Esc to give up, or do calibration and then hit any\n");
							printf("other key to retry:"); fflush(stdout);
						} else {
							printf("\nStrip read failed due to misread\n");
							printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
						}
						if (next_con_char() == 0x1b) {
							printf("\n");
							return -1;
						}
						printf("\n");
						pix--;
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
						printf("\nStrip read failed due to communication problem.\n");
						printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
						if (next_con_char() == 0x1b) {
							printf("\n");
							return -1;
						}
						printf("\n");
						if ((rv = it->init_coms(it, comport, br, 15.0)) != inst_ok) {
#ifdef DEBUG
							printf("init_coms returned '%s' (%s)\n",
						       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
							return -1;
						}
						pix--;
						continue;
					} else {
						/* Some other error. Treat it as fatal */
						printf("\nStrip read failed due unexpected error :'%s' (%s)\n",
				       	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
						printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
						if (next_con_char() == 0x1b) {
							printf("\n");
							return -1;
						}
						printf("\n");
						pix--;
						continue;
					}
				}

				/* Save the reading */
				if (val.XYZ_v == 0 && val.aXYZ_v == 0)
					error("Instrument didn't return XYZ value");
				
				if (val.XYZ_v != 0) {
					for (j = 0; j < 3; j++)
						scols[pix]->XYZ[j] = val.XYZ[j];
				} else {
					for (j = 0; j < 3; j++)
						scols[pix]->XYZ[j] = val.aXYZ[j];
				}

				/* Copy spectral */
				if (val.spec_n > 0) {
					scols[pix]->spec_n =        val.spec_n;
					scols[pix]->spec_wl_short = val.spec_wl_short;
					scols[pix]->spec_wl_long  = val.spec_wl_long;
					
					for (j = 0; j < val.spec_n; j++)
						scols[pix]->spec[j] = val.spec[j];
				}
				scols[pix]->rr = 1;		/* Has been read */
				printf(" Read OK\n");
				/* Advance to next patch. */

			} else {	/* Unrecognised response */
				pix--;
				continue;
			}
		}
	}
	/* clean up */
	if (it != NULL)
		it->del(it);
	return 0;
}

void
usage(void) {
	serio *sio;
	fprintf(stderr,"Read Print Target Test Chart colors, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	fprintf(stderr,"usage: printread [-options] outfile\n");
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
	fprintf(stderr," -i 41 | 51 | SO | SS Choose target instrument\n");
	fprintf(stderr,"                 41 = DTP41, 51 = DTP51, SS = Spectrolino or SpectroScan or SpectroScanT\n");
	fprintf(stderr," -t              Use transmission measurement mode\n");
	fprintf(stderr," -e              Use emissive  measurement mode (absolute results)\n");
	fprintf(stderr," -x [lx]         Take external values, either L*a*b* (-xl) or XYZ (-xx).\n");
	fprintf(stderr," -n              Don't save spectral information (default saves spectral)\n");
	fprintf(stderr," -l              Save CIE as D50 L*a*b* rather than XYZ\n");
	fprintf(stderr," outfile         Base name for input[ti2]/output[ti3] file\n");
	exit(1);
	}

main(int argc, char *argv[]) {
	int i,j;
	int fa,nfa;							/* current argument we're looking at */
	int verb = 0;
	int debug = 0;
	int comport = COMPORT;				/* COM port used */
	instType itype = instUnknown;		/* Default target instrument */
	int trans = 0;						/* Use transmission mode */
	int emiss = 0;						/* Use emissive mode */
	int xtern = 0;						/* Take external values, 1 = Lab, 2 = XYZ */
	int spectral = 1;					/* Save spectral information */
	int dolab = 0;						/* Save CIE as Lab */
	static char inname[200] = { 0 };	/* Input cgats file base name */
	static char outname[200] = { 0 };	/* Output cgats file base name */
	cgats *icg;					/* input cgats structure */
	cgats *ocg;					/* output cgats structure */
	int nmask = 0;				/* Device colorant mask */
	char *pixpat = "A-Z, A-Z";			/* Pass index pattern */		
	char *sixpat = "0-9,@-9,@-9;1-999";	/* Step index pattern */		
	alphix *paix, *saix;		/* Pass and Step index generators */
	int ixord = 0;				/* Index order, 0 = pass then step */
	int hex = 0;				/* Hexagon pattern layout */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	col *cols;					/* Internal storage of all the patch colors */
	col **scols;				/* Location sorted pointers to cols */
	int nchan;					/* Number of device chanels */
	int npat;					/* Number of overall patches */
	char *paists;				/* Passes in each Strip string, one char per pass, base 36 */
	int stipa;					/* Steps in each Pass */
	int totpa;					/* Total Passes Needed */
	int runpat;					/* Rounded Up to (totpa * stipa) Number of patches */
	int si;						/* Sample id index */
	int li;						/* Location id index */
	int ti;						/* Temp index */
	int fi;						/* Colorspace index */
	double plen = 0.0, glen = 0.0, tlen = 0.0;	/* Patch, gap and trailer length in mm */

	error_program = "Printread";
	if (argc <= 1)
		usage();

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
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

			/* COM port  */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage();
				comport = atoi(na);
				if (comport < 1 || comport > 4) usage();
			}

			/* Override Target Instrument type */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				char *p;
				int tt;
				fa = nfa;
				if (na == NULL) usage();

				if (strcmp("51", na) == 0)
					itype = instDTP51;
				else if (strcmp("41", na) == 0)
					itype = instDTP41;
				else if (strcmp("SO", na) == 0 || strcmp("so", na) == 0)
					itype = instSpectrolino;
				else if (strcmp("SS", na) == 0 || strcmp("ss", na) == 0) {
					itype = instSpectroScan;
				} else
					usage();
			}

			/* Request transmission measurement */
			else if (argv[fa][1] == 't' || argv[fa][1] == 'T') {
				emiss = 0;
				trans = 1;
			}

			/* Request emissive measurement */
			else if (argv[fa][1] == 'e' || argv[fa][1] == 'E') {
				emiss = 1;
				trans = 0;
			}

			/* Request external values */
			else if (argv[fa][1] == 'x' || argv[fa][1] == 'X') {
				char *p;
				int tt;
				fa = nfa;
				if (na == NULL) usage();

				if (na[0] == 'l' || na[0] == 'L')
					xtern = 1;
				else if (na[0] == 'x' || na[0] == 'X')
					xtern = 2;
				else
					usage();
			}

			/* Turn off spectral measurement */
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N')
				spectral = 0;

			/* Save as Lab */
			else if (argv[fa][1] == 'l' || argv[fa][1] == 'L')
				dolab = 1;

			else 
				usage();
		} else
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

	if (icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI2 format file");
	if (icg->ntables != 1)
		error ("Input file doesn't contain exactly one table");

	if ((npat = icg->t[0].nsets) <= 0)
		error ("No sets of data");

	/* Setup output cgats file */
	ocg = new_cgats();				/* Create a CGATS structure */
	ocg->add_other(ocg, "CTI3"); 		/* our special type is Calibration Target Information 3 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll printread", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	ocg->add_kword(ocg, 0, "DEVICE_CLASS","OUTPUT", NULL);	/* What sort of device this is */

	if (itype == instUnknown) {
		if ((ti = icg->find_kword(icg, 0, "TARGET_INSTRUMENT")) >= 0) {

				if ((itype = inst_enum(icg->t[0].kdata[ti])) == instUnknown)
					error ("Unrecognised target instrument '%s'", icg->t[0].kdata[ti]);
		} else {
			itype = instDTP41;		/* Default target instrument */
		}
	}

	/* Note what instrument the chart was read with */
	ocg->add_kword(ocg, 0, "TARGET_INSTRUMENT", inst_name(itype) , NULL);

	if ((ti = icg->find_kword(icg, 0, "SINGLE_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "SINGLE_DIM_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "COMP_GREY_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "COMP_GREY_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "MULTI_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "MULTI_DIM_STEPS",icg->t[0].kdata[ti], NULL);

	if ((ti = icg->find_kword(icg, 0, "FULL_SPREAD_PATCHES")) >= 0)
		ocg->add_kword(ocg, 0, "FULL_SPREAD_PATCHES",icg->t[0].kdata[ti], NULL);
	
	if ((ti = icg->find_kword(icg, 0, "STEPS_IN_PASS")) < 0)
		error ("Input file doesn't contain keyword STEPS_IN_PASS");
	stipa = atoi(icg->t[0].kdata[ti]);

	if ((ti = icg->find_kword(icg, 0, "PASSES_IN_STRIPS")) < 0)
		error ("Input file doesn't contain keyword PASSES_IN_STRIPS");
	paists = icg->t[0].kdata[ti];

	/* Get specified location indexing patterns */
	if ((ti = icg->find_kword(icg, 0, "STRIP_INDEX_PATTERN")) >= 0)
		pixpat = icg->t[0].kdata[ti];
	if ((ti = icg->find_kword(icg, 0, "PATCH_INDEX_PATTERN")) >= 0)
		sixpat = icg->t[0].kdata[ti];
	if ((ti = icg->find_kword(icg, 0, "INDEX_ORDER")) >= 0) {
		if (strcmp(icg->t[0].kdata[ti], "PATCH_THEN_STRIP") == 0)
			ixord = 1;
	}
	if ((ti = icg->find_kword(icg, 0, "HEXAGON_PATCHES")) >= 0)
		hex = 1;

	if ((paix = new_alphix(pixpat)) == NULL)
		error("Strip indexing pattern '%s' doesn't parse",pixpat);

	if ((saix = new_alphix(sixpat)) == NULL)
		error("Patch in strip indexing pattern '%s' doesn't parse",sixpat);

	if ((ti = icg->find_kword(icg, 0, "TOTAL_INK_LIMIT")) >= 0)
		ocg->add_kword(ocg, 0, "TOTAL_INK_LIMIT",icg->t[0].kdata[ti], NULL);

	if (itype == instDTP41) {	/* DTP41 specific info */
		if ((ti = icg->find_kword(icg, 0, "PATCH_LENGTH")) < 0)
			error ("Input file doesn't contain keyword PATCH_LENGTH");
		plen = atof(icg->t[0].kdata[ti]);
		if ((ti = icg->find_kword(icg, 0, "GAP_LENGTH")) < 0)
			error ("Input file doesn't contain keyword GAP_LENGTH");
		glen = atof(icg->t[0].kdata[ti]);
		if ((ti = icg->find_kword(icg, 0, "TRAILER_LENGTH")) < 0) {
			tlen = 18.0;	/* Default for backwards compatibility */
		} else 
			tlen = atof(icg->t[0].kdata[ti]);
	}

	if (verb) {
		printf("Steps in each Pass = %d\n",stipa);
		printf("Passes in each Strip = ",paists);
		for (i = 0; paists[i] != '\000'; i++) {
			printf("%s%d",i > 0 ? ", " : "", b62_int(&paists[i]));
		}
		printf("\n");
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

	totpa = (npat + stipa -1)/stipa;	/* Total passes for all strips */
	runpat = stipa * totpa;				/* Rounded up totao number of patches */
	if ((cols = (col *)malloc(sizeof(col) * runpat)) == NULL)
		error("Malloc failed!");
	if ((scols = (col **)calloc(sizeof(col *), runpat)) == NULL)
		error("Malloc failed!");

	/* Figure out the color space */
	if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
		error ("Input file doesn't contain keyword COLOR_REPS");

	if ((nmask = icx_char2inkmask(icg->t[0].kdata[fi])) != 0) {
		int i, j, ii;
		int chix[ICX_MXINKS];	/* Device chanel indexes */
		int xyzix[3];			/* XYZ/Lab chanel indexes */
		char *ident;
		char *xyzfname[3] = { "XYZ_X", "XYZ_Y", "XYZ_Z" };
		char *labfname[3] = { "LAB_L", "LAB_A", "LAB_B" };
		int gotexyz = 1;		/* Flag set if input file already has approx XYZ */ 
		icxColorantLu *clu;	/* Xcolorants model based Device -> CIE */

		nchan = icx_noofinks(nmask);
		ident = icx_inkmask2char(nmask); 

		/* Device channels */
		for (j = 0; j < nchan; j++) {
			int imask;
			char fname[100];

			imask = icx_index2ink(nmask, j);
			sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : ident,
			                      icx_ink2char(imask));

			if ((ii = icg->find_field(icg, 0, fname)) < 0)
				error ("Input file doesn't contain field %s",fname);
			if (icg->t[0].ftype[ii] != r_t)
				error ("Field %s is wrong type",fname);
	
			ocg->add_field(ocg, 0, fname, r_t);
			chix[j] = ii;
		}

		/* Approximate XYZ and real XYZ */
		for (j = 0; j < 3; j++) {
			if ((ii = icg->find_field(icg, 0, xyzfname[j])) >= 0) {
				if (icg->t[0].ftype[ii] != r_t)
					error ("Field %s is wrong type",xyzfname[j]);
				xyzix[j] = ii;
			} else {
				gotexyz = 0;
			}
	
			if (dolab)
				ocg->add_field(ocg, 0, labfname[j], r_t);
			else
				ocg->add_field(ocg, 0, xyzfname[j], r_t);
		}

		if ((clu = new_icxColorantLu(nmask)) == NULL)
			error ("Creation of xcolorant lu object failed");

		{
			char fname[100];
			if (dolab)
				sprintf(fname, "%s_LAB", ident);
			else
				sprintf(fname, "%s_XYZ", ident);
			ocg->add_kword(ocg, 0, "COLOR_REP", fname, NULL);
		}

		/* Read all the test patches in */
		for (i = 0; i < npat; i++) {
			cols[i].id = ((char *)icg->t[0].fdata[i][si]);
			cols[i].loc = ((char *)icg->t[0].fdata[i][li]);
			cols[i].n  = nchan;
			for (j = 0; j < nchan; j++)
				cols[i].dev[j] = *((double *)icg->t[0].fdata[i][chix[j]]) / 100.0;
			if (gotexyz) {
				for (j = 0; j < 3; j++)
					cols[i].eXYZ[j] = *((double *)icg->t[0].fdata[i][xyzix[j]]);
			} else {
				clu->dev_to_XYZ(clu, cols[i].eXYZ, cols[i].dev);
				for (j = 0; j < 3; j++)
					cols[i].eXYZ[j] *= 100.0;
			}
			cols[i].XYZ[0] = cols[i].XYZ[1] = cols[i].XYZ[2] = -1.0;
		}

		for (; i < runpat; i++) {
			cols[i].id = cols[i].loc = "-1";
			for (j = 0; j < nchan; j++)
				cols[i].dev[j] = 0.0;
			clu->dev_to_XYZ(clu, cols[i].eXYZ, cols[i].dev);
			for (j = 0; j < 3; j++)
				cols[i].eXYZ[j] *= 100.0;
			cols[i].XYZ[0] = cols[i].XYZ[1] = cols[i].XYZ[2] = -1.0;
		}

		clu->del(clu);
		free(ident);
	} else
		error ("Input file keyword COLOR_REPS has unknown value");

	/* Set up the sorted array of pointers */
	for (i = 0; i < npat; i++) {
		scols[i] = &cols[i];
		if ((cols[i].loci = patch_location_order(paix, saix, ixord, cols[i].loc)) < 0)
			error ("Bad location field value '%s' on patch %d", cols[i].loc, i);
	}
	for (; i < runpat; i++) {
		scols[i] = &cols[i];
		cols[i].loci = (totpa-1) * (256 - stipa) + i; 
/* printf("~~extra = %d, %d\n",cols[i].loci >> 8, cols[i].loci & 255); */
	}

	/* Reset 'read' flag and all data */
	for (i = 0; i < runpat; i++) {
		cols[i].rr = 0;
		cols[i].XYZ[0] = -1.0;
		cols[i].XYZ[1] = -1.0;
		cols[i].XYZ[2] = -1.0;
		cols[i].spec_n = 0;
	}

#define HEAP_COMPARE(A,B) (A->loci < B->loci)
	HEAPSORT(col *, scols, npat);

	/* Read all of the strips in */
	if (read_strips(itype, scols, npat, totpa, stipa, paists, paix, saix, ixord, hex, comport,
	                plen, glen, tlen, trans, emiss, xtern, spectral, verb, debug) == 0) {
		/* And save the result */

		int nsetel = 0;
		cgats_set_elem *setel;	/* Array of set value elements */

		nsetel += 1;		/* For id */
		nsetel += nchan;	/* For device values */
		nsetel += 3;		/* For XYZ or Lab */

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

			if (strcmp(cols[i].id, "0") == 0)	/* This is a padding patch */
				continue;			/* Skip it */

			setel[k++].c = cols[i].id;
			
			for (j = 0; j < nchan; j++)
				setel[k++].d = 100.0 * cols[i].dev[j];

			if (dolab) {
				double lab[3];
				double xyz[3];

				xyz[0] = cols[i].XYZ[0]/100.0;
				xyz[1] = cols[i].XYZ[1]/100.0;
				xyz[2] = cols[i].XYZ[2]/100.0;
				icmXYZ2Lab(&icmD50, lab, xyz);
				setel[k++].d = lab[0];
				setel[k++].d = lab[1];
				setel[k++].d = lab[2];
			} else {
				setel[k++].d = cols[i].XYZ[0];
				setel[k++].d = cols[i].XYZ[1];
				setel[k++].d = cols[i].XYZ[2];
			}

			for (j = 0; j < cols[i].spec_n; j++) {
				setel[k++].d = cols[i].spec[j];
			}

			ocg->add_setarr(ocg, 0, setel);
		}

		free(setel);

		if (ocg->write_name(ocg, outname))
			error("Write error : %s",ocg->err);
	}

	saix->del(saix);
	paix->del(paix);
	free(cols);
	ocg->del(ocg);		/* Clean up */
	icg->del(icg);		/* Clean up */

	return 0;
}


