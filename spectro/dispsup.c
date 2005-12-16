

/* 
 * Argyll Color Correction System
 * Common display patch reading support.
 *
 * Author: Graeme W. Gill
 * Date:   2/11/2005
 *
 * Copyright 1998 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
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
#include "xspect.h"
#include "serio.h"
#include "inst.h"
#include "dtp92.h"
#include "gretag.h"
#include "dispwin.h"
#include "dispsup.h"


/* User requested calibration of the instrument */
int disprd_calibration(
instType itype,		/* Instrument type */
int comport, 		/* COM port used */
int override,		/* Override_redirect on X11 */
double patsize,		/* Size of dispwin */
int verb,			/* Verbosity flag */
int debug			/* Debug flag */
) {
	inst *it;
	int n, i, j, c;
	inst_code rv;
	int patch;
	char buf[100];
	ipatch val;
	baud_rate br = baud_19200;
	dispwin *dw;

	if (verb)
		printf("Setting up the instrument\n");

	it = new_inst(itype);
	if (debug)
		it->debug = debug;	/* Turn on debugging */
	if (verb)
		it->verb = verb;	/* Turn on verbosity */

	/* Establish communications */
	if ((rv = it->init_coms(it, comport, br, 15.0)) != inst_ok) {
#ifdef DEBUG
		printf("init_coms returned '%s' (%s)\n",
		       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
		return -1;
	}

	/* Initialise the instrument */
	if ((rv = it->init_inst(it)) != inst_ok) {
#ifdef DEBUG
		printf("init_inst returned '%s' (%s)\n",
		       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
		return -1;
	}

	if (itype == instDTP92) {
		dispwin *dw;

		/* Offset calibration */
		empty_con_chars();
		printf("Offset calibration:\n");
		printf("Place instrument on a dark surface.\n");
		printf("Hit Esc to skip, any other key to start calibration:"); fflush(stdout);
		if ((c = next_con_char()) == 0x1b) {
			printf("\n");
		} else {
			printf("\n");

			if ((rv = it->calibrate(it, 0)) != inst_ok) {
#ifdef DEBUG
				printf("calibrate returned '%s' (%s)\n",
				       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
				printf("Calibrate failed with '%s' (%s)\n",
			       it->inst_interp_error(it, rv), it->interp_error(it, rv));
				return -1;
			} else {
				printf("Offset calibration successful\n");
			}

			/* The user may remove the instrument */
			printf("The instrument can be removed from the dark surface.\n");
		}

		/* Cell ratio calibration */
		/* Open display window */
		if ((dw = new_dispwin(patsize, patsize, 0.0, 0.0, 0, override)) == NULL) {
#ifdef DEBUG
			printf("Failed to creat a display window \n");
#endif /* DEBUG */
			return -1;
		}

		empty_con_chars();
		printf("Cell ratio calibration:\n");
		printf("Place instrument on test window.\n");
		printf("Hit Esc to skip, any other key to start calibration:"); fflush(stdout);
		if ((c = next_con_char()) == 0x1b) {
			printf("\n");
		} else {
			int i;
			double gy;

			printf("\n");

			for (i = 0, gy = 0.6; i < 3; i++) {

				dw->set_color(dw, gy, gy, gy);

				if ((rv = it->calibrate(it, 1)) != inst_ok) {
#ifdef DEBUG
					printf("calibrate returned '%s' (%s)\n",
					       it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
					if ((rv & 0xff) == DTP92_TOO_MUCH_LIGHT)
						gy *= 0.7;
					else if ((rv & 0xff) == DTP92_NOT_ENOUGH_LIGHT) {
						gy *= 1.4;
						if (gy > 1.0)
							gy = 1.0;
					} else  {
						printf("Calibrate failed with '%s' (%s)\n",
					       it->inst_interp_error(it, rv), it->interp_error(it, rv));
						dw->del(dw);
						return -1;
					}
				} else {
					printf("Cell ratio calibration successful\n");
					break;		/* Succeeded */
				}

				/* Try again */
			}

			if (i >= 3) {
				printf("Calibration failed - too many tries.\n");
				dw->del(dw);
				return -1;
			}

			/* Remove sample window */
			dw->del(dw);

			/* The user may remove the instrument */
			printf("The instrument can be removed from the screen.\n");

		}
	}
	else if (itype == instSpectrolino) {

		if ((rv = it->set_mode(it, inst_mode_emis_disp)) != inst_ok) {
#ifdef DEBUG
			printf("set_mode returned '%s' (%s)\n",
				   it->inst_interp_error(it, rv), it->interp_error(it, rv));
#endif /* DEBUG */
			printf("Set_mode failed with '%s' (%s)\n",
			       it->inst_interp_error(it, rv), it->interp_error(it, rv));
			return -1;
		}

		if ((rv = it->calibrate(it, 0)) != inst_ok) {
			printf("Calibrate failed with '%s' (%s)\n",
		       it->inst_interp_error(it, rv), it->interp_error(it, rv));
			return -1;
		}
	}

	/* clean up */
	it->del(it);

	return 0;
}

/* Take a series of readings from the display */
/* Return nz on fail/abort */
static int disprd_read(disprd *p,
	col *cols,		/* Array of patch colors to be tested */
	int npat, 		/* Number of patches to be tested */
	int spat,		/* Start patch index for "verb", 0 if not used */
	int tpat		/* Total patch index for "verb", 0 if not used */
) {
	int j, rv;
	int ttpat = tpat;
	int patch;
	ipatch val;		/* Return value */

	if (ttpat < npat)
		ttpat = npat;

	for (patch = 0; patch < npat; patch++) {
		col *scb = &cols[patch];

		scb->XYZ_v = 0;		/* No readings are valid */
		scb->aXYZ_v = 0;
		scb->spec_n = 0;

		if (p->verb) {
			fprintf(p->df,"\rpatch %d of %d",spat+patch+1,ttpat);
			fflush(p->df);
		}
#ifdef DEBUG
		printf("About to read patch %d\n",patch);
#endif /* DEBUG */

		p->dw->set_color(p->dw, scb->r, scb->g, scb->b);

		for (;;) {	/* Until we give up retrying */
			val.XYZ_v = 0;		/* No readings are valid */
			val.aXYZ_v = 0;
			val.spec_n = 0;
	
			if ((rv = p->it->read_sample(p->it, scb->id, &val)) != inst_ok) {
#ifdef DEBUG
				printf("read_sample returned '%s' (%s)\n",
			       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv));
#endif /* DEBUG */
				/* Deal with a user abort */
				if ((rv & inst_mask) == inst_user_abort) {
					empty_con_chars();
					printf("\nSample read stopped at user request!\n");
					printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
					if (next_con_char() == 0x1b) {
						printf("\n");
						return 1;
					}
					printf("\n");
				continue;
				/* Deal with a misread or needs calibration */
				} else if ((rv & inst_mask) == inst_misread || (rv & inst_mask) == inst_needs_cal) {
					empty_con_chars();
					if ((rv & inst_mask) == inst_needs_cal) {
						printf("\nSample read failed because instruments needs calibration\n");
						printf("Hit Esc to give up, or do calibration and then hit any\n");
						printf("other key to retry:"); fflush(stdout);
					} else {
						printf("\nSample read failed due to misread\n");
						printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
					}
					if (next_con_char() == 0x1b) {
						printf("\n");
						return 1;
					}
					printf("\n");
					continue;
				/* Deal with a communications error */
				} else if ((rv & inst_mask) == inst_coms_fail) {
					int tt = p->it->last_sioerr(p->it);
					if (tt & (SIO_BRK | SIO_FER | SIO_PER | SIO_OER)) {
						if (p->br == baud_19200) p->br = baud_9600;
						else if (p->br == baud_9600) p->br = baud_4800;
						else if (p->br == baud_2400) p->br = baud_1200;
						else p->br = baud_1200;
					}
					/* Communication problem, allow retrying at a lower baud rate */
					empty_con_chars();
					printf("\nSample read failed due to communication problem.\n");
					printf("Hit Esc to give up, any other key to retry:"); fflush(stdout);
					if (next_con_char() == 0x1b) {
						printf("\n");
						return 1;
					}
					printf("\n");
					if ((rv = p->it->init_coms(p->it, p->comport, p->br, 15.0)) != inst_ok) {
#ifdef DEBUG
						printf("init_coms returned '%s' (%s)\n",
					       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv));
#endif /* DEBUG */
						return 1;
					}
				continue;
				}
			} else {
				break;		/* Sucesful reading */
			}
		}
		/* We only fall through with a valid reading */
#ifdef DEBUG
		printf("got reading, transfering to col\n");
#endif /* DEBUG */

		/* Copy relative XYZ */
		if (val.XYZ_v != 0) {
			for (j = 0; j < 3; j++)
				scb->XYZ[j] = val.XYZ[j];
			scb->XYZ_v = 1;
		}

		/* Copy absolute XYZ */
		if (val.aXYZ_v != 0) {
			for (j = 0; j < 3; j++)
				scb->aXYZ[j] = val.aXYZ[j];
			scb->aXYZ_v = 1;
		}

		/* Copy spectral */
		if (p->spectral && val.spec_n > 0) {
			scb->spec_n = val.spec_n;
			scb->spec_wl_short = val.spec_wl_short;
			scb->spec_wl_long  = val.spec_wl_long;
			
			for (j = 0; j < val.spec_n; j++)
				scb->spec[j] = val.spec[j];
		}
#ifdef DEBUG
		printf("on to next reading\n");
#endif /* DEBUG */
	}
	if (p->verb && tpat == 0)
		fprintf(p->df,"\n");
	return 0;
}

static void disprd_del(disprd *p) {

	/* The user may remove the instrument */
	printf("The instrument can be removed from the screen.\n");

	if (p->it != NULL)
		p->it->del(p->it);
	if (p->dw != NULL) {
		if (p->or != NULL)
			p->dw->set_ramdac(p->dw,p->or);
		p->dw->del(p->dw);
	}
	free(p);
}

/* Create a display reading object. */
/* Return NULL if error */
disprd *new_disprd(
instType itype,		/* Instrument type */
int comport, 		/* COM port used */
int donat,			/* Use ramdac for native output, else run through current or set ramdac */
double cal[3][256],	/* Calibration set/return (cal[0][0] < 0.0 if can't/not to be used) */
int override,		/* Override_redirect on X11 */
double patsize,		/* Size of dispwin */
double ho,			/* Horizontal offset */
double vo,			/* Vertical offset */
int spectral,		/* Generate spectral info flag */
int verb,			/* Verbosity flag */
FILE *df,			/* Verbose output - NULL = stdout */
int debug			/* Debug flag */
) {
	disprd *p;
	int c;
	inst_code rv;
	
	/* Allocate a disprd */
	if ((p = (disprd *)calloc(sizeof(disprd), 1)) == NULL) {
		return NULL;
	}
	p->del = disprd_del;
	p->read = disprd_read;

	p->verb = verb;
	if (df)
		p->df = df;
	else
		p->df = stdout;
	p->comport = comport;
	p->br = baud_19200;

	if (verb)
		fprintf(p->df,"Setting up the instrument\n");

	if ((p->it = new_inst(itype)) == NULL) {
		p->del(p);
		return NULL;
	}

	if (debug)
		p->it->debug = debug;	/* Turn on debugging */

	/* Establish communications */
	if ((rv = p->it->init_coms(p->it, p->comport, p->br, 15.0)) != inst_ok) {
#ifdef DEBUG
		printf("init_coms returned '%s' (%s)\n",
		       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv));
#endif /* DEBUG */
		p->del(p);
		return NULL;
	}

	/* Initialise the instrument */
	if ((rv = p->it->init_inst(p->it)) != inst_ok) {
#ifdef DEBUG
		printf("init_inst returned '%s' (%s)\n",
		       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv));
#endif /* DEBUG */
		p->del(p);
		return NULL;
	}

	/* Configure the instrument mode */
	{
		inst_capability cap;
		inst_mode mode = 0;

		cap = p->it->capabilities(p->it);

		if ((cap & inst_emis_disp) == 0) {
			printf("Need display measurement capability,\n");
			printf("but instrument doesn't support it\n");
			p->del(p);
			return NULL;
		}

		if (spectral && (cap & inst_spectral) == 0) {
			printf("Spectral information was requested,\n");
			printf("but instrument doesn't support it\n");
			spectral = 0;
		}

		mode = inst_mode_emis_disp;

		if (spectral) {
			mode |= inst_mode_spectral;
			p->spectral = 1;
		} else {
			p->spectral = 0;
		}
	

		if ((rv = p->it->set_mode(p->it,mode)) != inst_ok) {
#ifdef DEBUG
			printf("set_mode returned '%s' (%s)\n",
			       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv));
#endif /* DEBUG */
			p->del(p);
			return NULL;
		}

		/* Do a calibration up front, so as not to get in the users way. */
		if (p->it->needs_calibration(p->it) == inst_needs_cal) {
			if ((rv = p->it->calibrate(p->it, 0)) != inst_ok) {
				printf("Calibrate failed with '%s' (%s)\n",
			       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv));
				p->del(p);
				return NULL;
			}
		}
	}

	/* Open display window */
	if ((p->dw = new_dispwin(patsize, patsize, ho, vo, donat, override)) == NULL) {
#ifdef DEBUG
		printf("Failed to creat a display window \n");
#endif /* DEBUG */
		p->del(p);
		return NULL;
	}

	/* Save current RAMDAC so that we can restore it */
	if ((p->or = p->dw->get_ramdac(p->dw)) == NULL) {
		if (p->verb)
			fprintf(p->df,"Unable to read or set display RAMDAC\n");
	}

	/* Set the given RAMDAC so we can characterise through it */
	if (cal != NULL && cal[0][0] >= 0.0 && p->or != NULL) {
		ramdac *r;
		
		r = p->or->clone(p->or);
		if (r->nent != 256) {
			if (p->verb)
				fprintf(p->df,"Can't currently use a RAMDAC which doesn't have 256 entries");
		} else {
			int j, i;
			/* Set the ramdac contents */
			for (j = 0; j < 3; j++) {
				for (i = 0; i < r->nent; i++) {
					r->v[j][i] = cal[j][i];
				}
			}
			if (p->dw->set_ramdac(p->dw, r)) {
				if (p->verb) {
					fprintf(p->df,"Failed to set RAMDAC to desired calibration.\n");
					fprintf(p->df,"Perhaps the operating system is being fussy ?\n");
				}
				r->del(r);
				p->del(p);
				return NULL;
			}
			r->del(r);
		}
	}

	/* Return the ramdac being used */
	if (p->or != NULL && cal != NULL) {
		ramdac *r;
		
		if ((r = p->dw->get_ramdac(p->dw)) == NULL) {
			if (p->verb)
				fprintf(p->df,"Failed to read current RAMDAC");
		} else {
			if (r->nent != 256) {
				if (p->verb)
					fprintf(p->df,"Can't currently read a RAMDAC which doesn't have 256 entries");

			} else {
				int j, i;

				/* Get the ramdac contents */
				for (j = 0; j < 3; j++) {
					for (i = 0; i < r->nent; i++) {
						cal[j][i] = r->v[j][i];
					}
				}
			}
			r->del(r);
		}
	}
	

	/* Ask user to put instrument on screen */
	empty_con_chars();
	printf("Place instrument on test window.\n");
	printf("Hit Esc to give up, any other key to continue:"); fflush(stdout);
	if ((c = next_con_char()) == 0x1b) {
		printf("\n");
		p->del(p);
		return NULL;
	}
	printf("\n");

	return p;
}

