

/* 
 * Argyll Color Correction System
 * Common display patch reading support.
 *
 * Author: Graeme W. Gill
 * Date:   2/11/2005
 *
 * Copyright 1998 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#ifdef __MINGW32__
# define WINVER 0x0500
#endif

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
#include "insttypes.h"
#include "icoms.h"
#include "inst.h"
#include "dispwin.h"
#include "dispsup.h"

#undef DEBUG

#ifdef DEBUG
#define DBG(xxx) printf xxx ;
#else
#define DBG(xxx) 
#endif

#define FAKE_NOISE 0.0		/* Add noise to fake devices XYZ */

/* -------------------------------------------------------- */

/* A default calibration user interaction handler using the console. */
/* This handles both normal and display based calibration interaction */
/* with the instrument. */
inst_code inst_handle_calibrate(
	inst *p,
	inst_cal_type calt,		/* type of calibration to do. inst_calt_all for all pending. */
	inst_cal_cond calc,		/* Current condition. inst_calc_none for not setup */
	disp_win_info *dwi      /* Optional Information to be able to open a display test patch */
) {
	inst_code rv = inst_ok, ev;
	int usermes = 0;		/* User was given a message */
	char id[200];			/* Condition identifier */
	int ch;
	dispwin *dw, *_dw = NULL;	/* Display window to display test patches on, NULL if none. */

	DBG(("inst_handle_calibrate called\n"))

	/* Untill we're done with the calibration */
	for (;;) {

		DBG(("About to call calibrate at top of loop\n"))
	    ev = p->calibrate(p, calt, &calc, id);

		/* We're done */
		if ((ev & inst_mask) == inst_ok) {
			if (calc == inst_calc_message)
				printf("%s\n",id);
			if (usermes)
				printf("Calibration complete\n");
			if (_dw != NULL)
				_dw->del(_dw);
			DBG(("inst_handle_calibrate done 0x%x\n",ev))
			return ev;
		}

		/* User aborted */
		if ((ev & inst_mask) == inst_user_abort
		 || (ev & inst_mask) == inst_user_term) {
			if (_dw != NULL)
				_dw->del(_dw);
			DBG(("inst_handle_calibrate user aborted 0x%x\n",ev))
			return ev;
		}

		/* Retry on an error */
		if ((ev & inst_mask) != inst_cal_setup) {
			if ((ev & inst_mask) == inst_unsupported) {
				DBG(("inst_handle_calibrate calibration not supported\n"))
				return inst_unsupported;
			}

			printf("Calibration failed with '%s' (%s)\n",
				       p->inst_interp_error(p, ev), p->interp_error(p, ev));
			printf("Hit any key to retry, or Esc, ^C or Q to abort:\n");

			empty_con_chars();
			ch = next_con_char();
			printf("\n");
			if (ch == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
				if (_dw != NULL)
					_dw->del(_dw);
				DBG(("inst_handle_calibrate user aborted 0x%x\n",inst_user_abort))
				return inst_user_abort;
			}

		/* Get user to do/setup calibration */
		} else {

			switch (calc) {
				case inst_calc_uop_ref_white:
					printf("Do a reflective white calibration,\n");
					printf(" and then hit any key to continue,\n"); 
					printf(" or hit Esc, ^C or Q to abort:"); 
					break;

				case inst_calc_uop_trans_white:
					printf("Do a transmissive white calibration,\n");
					printf(" and then hit any key to continue,\n"); 
					printf(" or hit Esc, ^C or Q to abort:"); 
					break;
			
				case inst_calc_uop_trans_dark:
					printf("Do a transmissive dark calibration,\n");
					printf(" and then hit any key to continue,\n"); 
					printf(" or hit Esc, ^C or Q to abort:"); 
					break;
			
				case inst_calc_man_ref_white:
					printf("Place the instrument on its reflective white reference %s,\n",id);
					printf(" and then hit any key to continue,\n"); 
					printf(" or hit Esc, ^C or Q to abort:"); 
					break;

				case inst_calc_man_ref_whitek:
					printf("Click the instrument on its reflective white reference %s,\n",id);
					printf(" or hit Esc, ^C or Q to abort:"); 
					break;

				case inst_calc_man_ref_dark:
					printf("Place the instrument in the dark, not in contact with any surface,\n");
					printf(" and then hit any key to continue,\n"); 
					printf(" or hit Esc, ^C or Q to abort:"); 
					break;

				case inst_calc_man_em_dark:
					printf("Place the instrument on a dark surface,\n");
					printf(" and then hit any key to continue,\n"); 
					printf(" or hit Esc, ^C or Q to abort:"); 
					break;

				case inst_calc_man_trans_white:
					printf("Place the instrument on its transmissive white source,\n");
					printf(" and then hit any key to continue,\n"); 
					printf(" or hit Esc, ^C or Q to abort:"); 
					break;

				case inst_calc_man_trans_dark:
					printf("Use the appropriate tramissive blocking to block the transmission path,\n");
					printf(" and then hit any key to continue,\n"); 
					printf(" or hit Esc, ^C or Q to abort:"); 
					break;

				case inst_calc_change_filter:
					printf("Change filter on instrument to %s,\n",id);
					printf(" and then hit any key to continue,\n"); 
					printf(" or hit Esc, ^C or Q to abort:"); 
					break;

				case inst_calc_message:
					printf("%s\n",id);
					printf(" Hit any key to continue,\n"); 
					printf(" or hit Esc, ^C or Q to abort:"); 
					break;

				case inst_calc_disp_white:
					if (dwi == NULL) {	/* No way of creating a test window */
						DBG(("inst_handle_calibrate no way of creating test window 0x%x\n",inst_internal_error))
						return inst_internal_error; 
					}
				
					if (dwi->dw == NULL) {
						if ((dw = new_dispwin(dwi->disp, dwi->patsize, dwi->patsize,
						                      dwi->ho, dwi->vo, 0, 0, dwi->blackbg,
						                      dwi->override)) == NULL) {
							DBG(("inst_handle_calibrate failed to create test window 0x%x\n",inst_other_error))
							return inst_other_error; 
						}
						_dw = dw;	/* window to delete on exit */
						printf("Frequency calibration, Place instrument on test window.\n");
						printf(" Hit any key to continue,\n"); 
						printf(" or hit Esc, ^C or Q to abort:"); 
					} else {
						dw = dwi->dw;
					}
					p->cal_gy_level = 1.0;
					dw->set_color(dw, 1.0, 1.0, 1.0);
					break;

				case inst_calc_disp_grey:
				case inst_calc_disp_grey_darker: 
				case inst_calc_disp_grey_ligher:
					if (dwi == NULL) {	/* No way of creating a test window */
						DBG(("inst_handle_calibrate no way of creating test window 0x%x\n",inst_internal_error))
						return inst_internal_error; 
					}

					if (dwi->dw == NULL) {
						if ((dw = new_dispwin(dwi->disp, dwi->patsize, dwi->patsize,
						                      dwi->ho, dwi->vo, 0, 0, dwi->blackbg,
						                      dwi->override)) == NULL) {
							DBG(("inst_handle_calibrate failed to create test window 0x%x\n",inst_other_error))
							return inst_other_error; 
						}
						_dw = dw;	/* window to delete on exit */
						printf("Cell ratio calibration, Place instrument on test window.\n");
						printf(" Hit any key to continue,\n"); 
						printf(" or hit Esc, ^C or Q to abort:"); 
					} else {
						dw = dwi->dw;
					}
					if (calc == inst_calc_disp_grey) {
						p->cal_gy_level = 0.6;
						p->cal_gy_count = 0;
					} else if (calc == inst_calc_disp_grey_darker) {
						p->cal_gy_level *= 0.7;
						p->cal_gy_count++;
					} else if (calc == inst_calc_disp_grey_ligher) {
						p->cal_gy_level *= 1.4;
						if (p->cal_gy_level > 1.0)
							p->cal_gy_level = 1.0;
						p->cal_gy_count++;
					}
					if (p->cal_gy_count > 4) {
						printf("Cell ratio calibration failed - too many tries at setting grey level.\n");
						if (_dw != NULL)
							_dw->del(_dw);
						DBG(("inst_handle_calibrate too many tries at setting grey level 0x%x\n",inst_internal_error))
						return inst_internal_error; 
					}
					dw->set_color(dw, p->cal_gy_level, p->cal_gy_level, p->cal_gy_level);
					break;

				default:
					/* Something isn't being handled */
					if (_dw != NULL)
						_dw->del(_dw);
					DBG(("inst_handle_calibrate unhandled calc case 0x%x, err 0x%x\n",calc,inst_internal_error))
					return inst_internal_error;
			}
			fflush(stdout);

			usermes = 1;

			if (calc != inst_calc_man_ref_whitek) {
				empty_con_chars();
				ch = next_con_char();
				printf("\n");
				if (ch == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
					if (_dw != NULL)
						_dw->del(_dw);
					DBG(("inst_handle_calibrate user aborted 0x%x\n",inst_user_abort))
					return inst_user_abort;
				}
			}
		}
	}

	if (_dw != NULL)
		_dw->del(_dw);
	DBG(("inst_handle_calibrate done 0x%x\n",ev))
	return inst_ok;
}

/* -------------------------------------------------------- */
/* User requested calibration of the display instrument */
/* (Not supporting inst_crt_freq_cal since this */
/*  is done automatically in dispsup.) */
/* Should add an argument to be able to select type of calibration, */
/* rather than guessing what the user wants ? */
int disprd_calibration(
instType itype,			/* Instrument type (usually instUnknown) */
int comport, 			/* COM port used */
flow_control fc,		/* Serial flow control */
int dtype,				/* Display type, 0 = unknown, 1 = CRT, 2 = LCD */
int nocal,				/* NZ to disable auto instrument calibration */
disppath *disp,			/* display to calibrate. */
int blackbg,			/* NZ if whole screen should be filled with black */
int override,			/* Override_redirect on X11 */
double patsize,			/* Size of dispwin */
double ho, double vo,	/* Position of dispwin */
int verb,				/* Verbosity flag */
int debug				/* Debug flag */
) {
	inst *it;
	int c;
	inst_code rv;
	baud_rate br = baud_19200;
	disp_win_info dwi;
	inst_cal_type calt = inst_calt_all;
	inst_capability  cap;
	inst2_capability cap2;

	dwi.disp = disp; 
	dwi.blackbg = blackbg;
	dwi.override = override;
	dwi.patsize = patsize;
	dwi.ho = ho;
	dwi.vo = vo;
	dwi.dw = NULL;

	if (verb)
		printf("Setting up the instrument\n");

	if ((it = new_inst(comport, itype, debug, verb)) == NULL) {
		DBG(("new_inst failed\n"))
		return -1;
	}

	/* Establish communications */
	if ((rv = it->init_coms(it, comport, br, fc, 15.0)) != inst_ok) {
		DBG(("init_coms returned '%s' (%s)\n",
		       it->inst_interp_error(it, rv), it->interp_error(it, rv)))
		it->del(it);
		return -1;
	}

	/* Initialise the instrument */
	if ((rv = it->init_inst(it)) != inst_ok) {
		DBG(("init_inst returned '%s' (%s)\n",
		       it->inst_interp_error(it, rv), it->interp_error(it, rv)))
		it->del(it);
		return -1;
	}

	itype = it->get_itype(it);			/* Actual type */
	cap  = it->capabilities(it);
	cap2 = it->capabilities2(it);

	/* Set to emission mode to read a display */
	if ((rv = it->set_mode(it, inst_mode_emis_disp)) != inst_ok) {
		DBG(("Set_mode failed with '%s' (%s)\n",
		       it->inst_interp_error(it, rv), it->interp_error(it, rv)))
		return -1;
	}

	/* Set CRT or LCD mode */
	if ((cap & (inst_emis_disp_crt | inst_emis_disp_lcd))
	 && (dtype == 1 || dtype == 2)) {
		inst_opt_mode om;

		if (dtype == 1)
			om = inst_opt_disp_crt;
		else
			om = inst_opt_disp_lcd;

		if ((rv = it->set_opt_mode(it,om)) != inst_ok) {
			DBG(("Setting display type failed failed with '%s' (%s)\n",
			       it->inst_interp_error(it, rv), it->interp_error(it, rv)))
			it->del(it);
			return -1;
		}
	} else if (cap & (inst_emis_disp_crt | inst_emis_disp_lcd)) {
		printf("Either CRT or LCD must be selected\n");
		it->del(it);
		return -1;
	}

	/* Disable autocalibration of machine if selected */
	if (nocal != 0){
		if ((rv = it->set_opt_mode(it,inst_opt_noautocalib)) != inst_ok) {
			DBG(("Setting no-autocalibrate failed failed with '%s' (%s)\n",
			       it->inst_interp_error(it, rv), it->interp_error(it, rv)))
			it->del(it);
			return -1;
		}
	}

	/* Guess a calibration type */
	if (cap2 & inst2_cal_ref_white) {
		/* White or emmission mode calibration (ie. Spectrolino) */
		calt = inst_calt_ref_white;
	} else if (cap2 & inst2_cal_disp_offset) {
		/* Offset calibration */
		calt = inst_calt_disp_offset;
	} else if (cap2 & inst2_cal_disp_ratio) {
		/* Cell ratio calibration */
		calt = inst_calt_disp_ratio;
	} else if (cap2 & inst2_cal_crt_freq) {
		/* Frequency calibration for CRT */
		calt = inst_calt_crt_freq;
	}
	
	/* Do the calibration */
	rv = inst_handle_calibrate(it, calt, inst_calc_none, &dwi);
	if (rv == inst_unsupported) {
		printf("No calibration available for instrument in this mode\n");
	} else if (rv != inst_ok) {	/* Abort or fatal error */
		printf("Calibrate failed with '%s' (%s)\n",
	       it->inst_interp_error(it, rv), it->interp_error(it, rv));
		it->del(it);
		return -1;
	}

	/* clean up */
	it->del(it);

	return 0;
}

/* Take a series of readings from the display */
/* Return nz on fail/abort */
/* 1 = user aborted */
/* 2 = instrument access failed */
/* 3 = window access failed */ 
/* 4 = user hit terminate key */
/* 5 = system error */
/* Use disprd_err() to interpret it */
static int disprd_read(
	disprd *p,
	col *cols,		/* Array of patch colors to be tested */
	int npat, 		/* Number of patches to be tested */
	int spat,		/* Start patch index for "verb", 0 if not used */
	int tpat,		/* Total patch index for "verb", 0 if not used */
	int acr,		/* If nz, do automatic final carriage return */
	int tc			/* If nz, termination key */
) {
	int j, rv;
	int patch;
	ipatch val;		/* Return value */
	int ch;			/* Character */
	int cal_type;
	char id[CALIDLEN];

	/* Setup user termination character */
	p->it->icom->set_uih(p->it->icom, tc, tc, ICOM_TERM);

	/* See if we should do a frequency calibration first */
	if ((p->it->capabilities2(p->it) & inst2_cal_crt_freq) != 0
	 && p->it->needs_calibration(p->it) == inst_calt_crt_freq
	 && npat > 0
	 && (cols[0].r != 1.0 || cols[1].r != 1.0 || cols[1].r != 1.0)) {
		col tc;
		inst_cal_cond calc = inst_calc_disp_white;

		if ((rv = p->dw->set_color(p->dw, 1.0, 1.0, 1.0)) != 0) {
			DBG(("set_color() returned %s\n",rv))
			return 3;
		}
		/* Do calibrate, but ignore return code. Press on regardless. */
		if ((rv = p->it->calibrate(p->it, inst_calt_crt_freq, &calc, id)) != inst_ok) {
			DBG(("warning, frequency calibrate failed with '%s' (%s)\n",
				      p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv)))
		}
	}

	for (patch = 0; patch < npat; patch++) {
		col *scb = &cols[patch];

		scb->XYZ_v = 0;		/* No readings are valid */
		scb->aXYZ_v = 0;
		scb->sp.spec_n = 0;

		if (p->verb && spat != 0 && tpat != 0) {
			fprintf(p->df,"\rpatch %d of %d",spat+patch,tpat);
			fflush(p->df);
		}
		DBG(("About to read patch %d\n",patch))

		if ((rv = p->dw->set_color(p->dw, scb->r, scb->g, scb->b)) != 0) {
			DBG(("set_color() returned %s\n",rv))
			return 3;
		}

		/* Until we give up retrying */
		for (;;) {
			val.XYZ_v = 0;		/* No readings are valid */
			val.aXYZ_v = 0;
			val.sp.spec_n = 0;
	
			if ((rv = p->it->read_sample(p->it, scb->id, &val)) != inst_ok
			     && (rv & inst_mask) != inst_user_trig) {
				DBG(("read_sample returned '%s' (%s)\n",
			       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv)))

				/* Deal with a user terminate */
				if ((rv & inst_mask) == inst_user_term) {
					return 4;

				/* Deal with a user abort */
				} else if ((rv & inst_mask) == inst_user_abort) {
					empty_con_chars();
					printf("\nSample read stopped at user request!\n");
					printf("Hit Esc, ^C or Q to give up, any other key to retry:"); fflush(stdout);
					if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
						printf("\n");
						return 1;
					}
					printf("\n");
					continue;

				/* Deal with a misread or needs calibration */
				} else if ((rv & inst_mask) == inst_needs_cal) {
					disp_win_info dwi;
					dwi.dw = p->dw;		/* Set window to use */
					printf("\nSample read failed because instruments needs calibration\n");
					rv = inst_handle_calibrate(p->it, inst_calt_all, inst_calc_none, &dwi);
					if (rv != inst_ok) {	/* Abort or fatal error */
						return 1;
					}
					continue;

				/* Deal with a misread */
				} else if ((rv & inst_mask) == inst_misread) {
					empty_con_chars();
					printf("\nSample read failed due to misread\n");
					printf("Hit Esc, ^C or Q to give up, any other key to retry:"); fflush(stdout);
					if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
						printf("\n");
						return 1;
					}
					printf("\n");
					continue;

				/* Deal with a communications error */
				} else if ((rv & inst_mask) == inst_coms_fail) {
					int tt = p->it->last_comerr(p->it);
					if (tt & (ICOM_BRK | ICOM_FER | ICOM_PER | ICOM_OER)) {
						if (p->br == baud_19200) p->br = baud_9600;
						else if (p->br == baud_9600) p->br = baud_4800;
						else if (p->br == baud_2400) p->br = baud_1200;
						else p->br = baud_1200;
					}
					/* Communication problem, allow retrying at a lower baud rate */
					empty_con_chars();
					printf("\nSample read failed due to communication problem.\n");
					printf("Hit Esc, ^C or Q to give up, any other key to retry:"); fflush(stdout);
					if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
						printf("\n");
						return 1;
					}
					printf("\n");
					if ((rv = p->it->init_coms(p->it, p->comport, p->br, p->fc, 15.0)) != inst_ok) {
						DBG(("init_coms returned '%s' (%s)\n",
					       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv)))
						return 2;
					}
				continue;
				}
			} else {
				break;		/* Sucesful reading */
			}
		}
		/* We only fall through with a valid reading */
		DBG(("got reading abs. %f %f %f, transfering to col\n",
		                val.aXYZ[0], val.aXYZ[1], val.aXYZ[2]))

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
		if (p->spectral && val.sp.spec_n > 0) {
			scb->sp = val.sp;
		}
		DBG(("on to next reading\n"))
	}
	/* Final return. */
	if (acr && p->verb && spat != 0 && tpat != 0 && (spat+patch-1) == tpat)
		fprintf(p->df,"\n");
	return 0;
}


/* Test without a spectrometer using a fake device */
static int disprd_fake_read(disprd *p,
	col *cols,		/* Array of patch colors to be tested */
	int npat, 		/* Number of patches to be tested */
	int spat,		/* Start patch index for "verb", 0 if not used */
	int tpat,		/* Total patch index for "verb", 0 if not used */
	int acr,		/* If nz, do automatic final carriage return */
	int tc			/* If nz, termination key */
) {
	icmXYZNumber white;		/* White point */
	icmXYZNumber red;		/* Red colorant */
	icmXYZNumber green;		/* Green colorant */
	icmXYZNumber blue;		/* Blue colorant */
	double doff[3];			/* device offsets */
	double mat[3][3];		/* Destination matrix */
	double ooff[3];			/* XYZ offsets */
	double rgb[3];
	double br = 35.4;		/* Overall brightness */
	int patch, j;
	int ttpat = tpat;

	if (ttpat < npat)
		ttpat = npat;

	/* Setup fake device */
	white.X = br * 0.955;		/* Somewhere between D50 and D65 */
	white.Y = br * 1.00;
	white.Z = br * 0.97;
	red.X = br * 0.41;
	red.Y = br * 0.21;
	red.Z = br * 0.02;
	green.X = br * 0.30;
	green.Y = br * 0.55;
	green.Z = br * 0.15;
	blue.X = br * 0.15;
	blue.Y = br * 0.10;
	blue.Z = br * 0.97;
	/* Input offset, equivalent to RGB offsets having various values */
	doff[0] = 0.10;
	doff[1] = 0.06;
	doff[2] = 0.08;
	/* Output offset - equivalent to flare */
	ooff[0] = 0.04;
	ooff[1] = 0.09;
	ooff[2] = 0.16;

	if (icmRGBprim2matrix(white, red, green, blue, mat))
		error("Fake read unexpectedly got singular matrix\n");

	for (patch = 0; patch < npat; patch++) {
		if (p->verb && spat != 0 && tpat != 0) {
			fprintf(p->df,"\rpatch %d of %d",spat+patch,tpat);
			fflush(p->df);
		}
		rgb[0] = cols[patch].r;
		rgb[1] = cols[patch].g;
		rgb[2] = cols[patch].b;

		/* Apply device offset */
		for (j = 0; j < 3; j++) {
			if (rgb[j] < 0.0)
				rgb[j] = 0.0;
			else if (rgb[j] > 1.0)
				rgb[j] = 1.0;
			rgb[j] = doff[j] + (1.0 - doff[j]) * rgb[j];
		}

		/* Apply gamma */
		for (j = 0; j < 3; j++) {
			if (rgb[j] >= 0.0)
				rgb[j] = pow(rgb[j], 2.5);
			else
				rgb[j] = 0.0;
		}

		/* Convert to XYZ */
		icmMulBy3x3(cols[patch].aXYZ, mat, rgb);

		/* Apply XYZ offset */
		for (j = 0; j < 3; j++)
			cols[patch].aXYZ[j] += ooff[j];

		cols[patch].aXYZ_v = 1;
	}
	if (acr && p->verb && spat != 0 && tpat != 0 && (spat+patch-1) == tpat)
		fprintf(p->df,"\n");
	return 0;
}

/* Test without a spectrometer using a fake ICC profile device */
static int disprd_fake_read_lu(disprd *p,
	col *cols,		/* Array of patch colors to be tested */
	int npat, 		/* Number of patches to be tested */
	int spat,		/* Start patch index for "verb", 0 if not used */
	int tpat,		/* Total patch index for "verb", 0 if not used */
	int acr,		/* If nz, do automatic final carriage return */
	int tc			/* If nz, termination key */
) {
	int patch, j;
	int ttpat = tpat;
	double br = 120.4;		/* Overall brightness */

	if (ttpat < npat)
		ttpat = npat;

	for (patch = 0; patch < npat; patch++) {
		double rgb[3];;
		if (p->verb && spat != 0 && tpat != 0) {
			fprintf(p->df,"\rpatch %d of %d",spat+patch,tpat);
			fflush(p->df);
		}
		rgb[0] = cols[patch].r;
		rgb[1] = cols[patch].g;
		rgb[2] = cols[patch].b;

		p->fake_lu->lookup(p->fake_lu, cols[patch].aXYZ, rgb); 
		for (j = 0; j < 3; j++) 
			cols[patch].aXYZ[j] *= br;
#ifdef FAKE_NOISE
		cols[patch].aXYZ[0] += 2.0 * FAKE_NOISE * d_rand(-1.0, 1.0);
		cols[patch].aXYZ[1] += FAKE_NOISE * d_rand(-1.0, 1.0);
		cols[patch].aXYZ[2] += 4.0 * FAKE_NOISE * d_rand(-1.0, 1.0);
		for (j = 0; j < 3; j++) { 
			if (cols[patch].aXYZ[j] < 0.0)
				cols[patch].aXYZ[j] = 0.0;
		}
#endif
		cols[patch].aXYZ_v = 1;
	}
	if (acr && p->verb && spat != 0 && tpat != 0 && (spat+patch-1) == tpat)
		fprintf(p->df,"\n");
	return 0;
}

/* Return a string describing the error code */
char *disprd_err(int en) {
	switch(en) {
		case 1:
			return "User Aborted";
		case 2:
			return "Instrument Access Failed";
		case 3:
			return "Window Access Failed";
		case 4:
			return "VideoLUT Access Failed";
		case 5:
			return "User Terminated";
		case 6:
			return "System Error";
		case 7:
			return "Either CRT or LCD must be selected";
	}
	return "Unknown";
}

static void disprd_del(disprd *p) {

	/* The user may remove the instrument */
	printf("The instrument can be removed from the screen.\n");

	if (p->fake_lu != NULL)
		p->fake_lu->del(p->fake_lu);
	if (p->fake_icc != NULL)
		p->fake_icc->del(p->fake_icc);
	if (p->fake_fp != NULL)
		p->fake_fp->del(p->fake_fp);

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
/* Set *errc to code: */
/* 1 = user aborted */
/* 2 = instrument access failed */
/* 3 = window access failed */ 
/* 4 = RAMDAC access failed */ 
/* 5 = user hit terminate key */
/* 6 = system error */
/* 7 = CRT or LCD must be selected */
/* Use disprd_err() to interpret *errc */
disprd *new_disprd(
int *errc,          /* Error code. May be NULL */
instType itype,		/* Nominal instrument type (usually instUnknown) */
int comport, 		/* COM port used. -99 == fake Display */
flow_control fc,	/* Flow control */
int dtype,			/* Display type, 0 = unknown, 1 = CRT, 2 = LCD */
int nocal,			/* No automatic instrument calibration */
int highres,		/* Use high res mode if available */
int donat,			/* Use ramdac for native output, else run through current or set ramdac */
double cal[3][256],	/* Calibration set/return (cal[0][0] < 0.0 if can't/not to be used) */
disppath *disp,		/* Display to calibrate. */
int blackbg,		/* NZ if whole screen should be filled with black */
int override,		/* Override_redirect on X11 */
char *callout,      /* Shell callout on set color */
double patsize,		/* Size of dispwin */
double ho,			/* Horizontal offset */
double vo,			/* Vertical offset */
int spectral,		/* Generate spectral info flag */
int verb,			/* Verbosity flag */
FILE *df,			/* Verbose output - NULL = stdout */
int debug,			/* Debug flag */
char *fake_name		/* Name of profile to use as a fake device */
) {
	disprd *p;
	int ch;
	inst_code rv;
	
	if (errc != NULL) *errc = 0;

	/* Allocate a disprd */
	if ((p = (disprd *)calloc(sizeof(disprd), 1)) == NULL) {
		if (errc != NULL) *errc = 6;
		return NULL;
	}
	p->del = disprd_del;
	p->read = disprd_read;
	p->fake_name = fake_name;

	p->verb = verb;
	if (df)
		p->df = df;
	else
		p->df = stdout;
	p->comport = comport;
	p->br = baud_19200;
	p->fc = fc;

	if (comport == -99) {
		p->fake = 1;

		/* See if there is a profile we should use as the fake device */
		p->fake_fp = NULL;
		p->fake_icc = NULL;
		p->fake_lu = NULL;
		if (p->fake_name != NULL && (p->fake_fp = new_icmFileStd_name(p->fake_name,"r")) != NULL) {
			if ((p->fake_icc = new_icc()) != NULL) {
				if (p->fake_icc->read(p->fake_icc,p->fake_fp,0) == 0) {
					icColorSpaceSignature ins;
					p->fake_lu = p->fake_icc->get_luobj(p->fake_icc, icmFwd, icAbsoluteColorimetric,
					                            icSigXYZData, icmLuOrdNorm);
					p->fake_lu->spaces(p->fake_lu, &ins, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
					if (ins != icSigRgbData) {
						p->fake_lu->del(p->fake_lu);
						p->fake_lu = NULL;
					}
				}
			}
		}

		if (p->fake_lu != NULL) {
//			if (verb)
				printf("Using profile '%s' rather than real device\n",p->fake_name);
			p->read = disprd_fake_read_lu;
		} else
			p->read = disprd_fake_read;
		return p;
	}

	if (verb)
		fprintf(p->df,"Setting up the instrument\n");

	if ((p->it = new_inst(comport, itype, debug, verb)) == NULL) {
		p->del(p);
		if (errc != NULL) *errc = 2;
		return NULL;
	}

	/* Establish communications */
	if ((rv = p->it->init_coms(p->it, p->comport, p->br, p->fc, 15.0)) != inst_ok) {
		DBG(("init_coms returned '%s' (%s)\n",
		       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv)))
		p->del(p);
		if (errc != NULL) *errc = 2;
		return NULL;
	}

	/* Initialise the instrument */
	if ((rv = p->it->init_inst(p->it)) != inst_ok) {
		DBG(("init_inst returned '%s' (%s)\n",
		       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv)))
		p->del(p);
		if (errc != NULL) *errc = 2;
		return NULL;
	}
	itype = p->it->get_itype(p->it);			/* Actual type */

	/* Configure the instrument mode */
	{
		inst_capability cap;
		inst2_capability cap2;
		inst_mode mode = 0;

		cap = p->it->capabilities(p->it);
		cap2 = p->it->capabilities2(p->it);

		if ((cap & inst_emis_disp) == 0) {
			printf("Need emissive measurement capability,\n");
			printf("but instrument doesn't support it\n");
			p->del(p);
			if (errc != NULL) *errc = 2;
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
	

		/* Set CRT or LCD mode */
		if ((cap & (inst_emis_disp_crt | inst_emis_disp_lcd))
		 && (dtype == 1 || dtype == 2)) {
			inst_opt_mode om;

			if (dtype == 1)
				om = inst_opt_disp_crt;
			else
				om = inst_opt_disp_lcd;

			if ((rv = p->it->set_opt_mode(p->it, om)) != inst_ok) {
				DBG(("Setting display type failed failed with '%s' (%s)\n",
			       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv)))
				p->del(p);
				if (errc != NULL) *errc = 2;
				return NULL;
			}
		} else if (cap & (inst_emis_disp_crt | inst_emis_disp_lcd)) {
			printf("Either CRT or LCD must be selected\n");
			p->del(p);
			if (errc != NULL) *errc = 7;
			return NULL;
		}

		/* Disable autocalibration of machine if selected */
		if (nocal != 0){
			if ((rv = p->it->set_opt_mode(p->it,inst_opt_noautocalib)) != inst_ok) {
				DBG(("Setting no-autocalibrate failed failed with '%s' (%s)\n",
			       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv)))
				p->del(p);
				if (errc != NULL) *errc = 2;
				return NULL;
			}
		}

		if (highres) {
			if (cap & inst_highres) {
				inst_code ev;
				if ((ev = p->it->set_opt_mode(p->it, inst_opt_highres)) != inst_ok) {
					printf("\nSetting high res mode failed with error :'%s' (%s)\n",
			       	       p->it->inst_interp_error(p->it, ev), p->it->interp_error(p->it, ev));
					if (errc != NULL) *errc = 2;
					return NULL;
				}
			} else if (p->verb) {
				printf("high resolution ignored - instrument doesn't support high res. mode\n");
			}
		}
		if ((rv = p->it->set_mode(p->it, mode)) != inst_ok) {
			DBG(("set_mode returned '%s' (%s)\n",
			       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv)))
			p->del(p);
			if (errc != NULL) *errc = 2;
			return NULL;
		}

		/* Set the trigger mode to program triggered */
		if ((rv = p->it->set_opt_mode(p->it,inst_opt_trig_prog)) != inst_ok) {
			DBG(("Setting program trigger mode failed failed with '%s' (%s)\n",
		       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv)))
			p->del(p);
			if (errc != NULL) *errc = 2;
			return NULL;
		}
	}

	/* Open display window */
	if ((p->dw = new_dispwin(disp, patsize, patsize, ho, vo, 0, donat, blackbg,
	                                                        override)) == NULL) {
		DBG(("Failed to creat a display window \n"))
		p->del(p);
		if (errc != NULL) *errc = 3;
		return NULL;
	}

	/* Set color change callout */
	if (callout) {
		p->dw->set_callout(p->dw, callout);
	}

	/* Save current RAMDAC so that we can restore it */
	if ((p->or = p->dw->get_ramdac(p->dw)) == NULL) {
		warning("Unable to read or set display RAMDAC");
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
				if (errc != NULL) *errc = 4;
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
	
	/* Do a calibration up front, so as not to get in the users way, */
	/* but ignore a CRT frequency calibration, since it will be done */
	/* automatically. */
	if ((p->it->needs_calibration(p->it) & inst_calt_needs_cal_mask) != 0
	 && p->it->needs_calibration(p->it) != inst_calt_crt_freq) {
		disp_win_info dwi;
		dwi.dw = p->dw;		/* Set window to use */

		rv = inst_handle_calibrate(p->it, inst_calt_all, inst_calc_none, &dwi);
		printf("\n");
		if (rv != inst_ok) {	/* Abort or fatal error */
			printf("Calibrate failed with '%s' (%s)\n",
			       p->it->inst_interp_error(p->it, rv), p->it->interp_error(p->it, rv));
			p->del(p);
			if (errc != NULL) *errc = 2;
			return NULL;
		}
	}

	/* Ask user to put instrument on screen */
	empty_con_chars();
	printf("Place instrument on test window.\n");
	printf("Hit Esc, ^C or Q to give up, any other key to continue:"); fflush(stdout);
	if ((ch = next_con_char()) == 0x1b || ch == 0x3 || ch == 'q' || ch == 'Q') {
		printf("\n");
		p->del(p);
		if (errc != NULL) *errc = 1;
		return NULL;
	}
	printf("\n");

	return p;
}

