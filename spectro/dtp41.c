/* 
 * Argyll Color Correction System
 *
 * Xrite DTP41 related functions
 *
 * Author: Graeme W. Gill
 * Date:   10/3/2001
 *
 * Copyright 1996 - 2001, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 *
 * Derived from DTP51.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include "xspect.h"
#include "serio.h"
#include "dtp41.h"

#undef DEBUG

#include <stdarg.h>
void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);

static inst_code interp_code(inst *pp, int ec);
static inst_code activate_mode(dtp41 *p);

#if defined (NT)
# define WIN32_LEAN_AND_MEAN
# include <windows.h>
#define sleep(secs) Sleep((secs) * 1000)
#endif

#define MAX_MES_SIZE 500		/* Maximum normal message reply size */
#define MAX_RD_SIZE 100000		/* Maximum reading messagle reply size */

static void
delay (double tt) {
	long etime;
	etime = clock() + (long)(CLOCKS_PER_SEC * tt + 0.5);
	while (clock() < etime);
}

/* Convert control chars to ^[A-Z] notation in a string */
static char *
fix(char *s) {
	static char buf [MAX_RD_SIZE];
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

/* Extract an error code from a reply string */
/* Return -1 if no error code can be found */
static int 
extract_ec(char *s) {
	char *p;
	char tt[3];
	int rv;
	p = s + strlen(s);
	/* Find the trailing '>' */
	for (p--; p >= s;p--) {
		if (*p == '>')
			break;
	}
	if (  (p-3) < s
	   || p[0] != '>'
	   || p[-3] != '<')
		return -1;
	tt[0] = p[-2];
	tt[1] = p[-1];
	tt[2] = '\000';
	if (sscanf(tt,"%x",&rv) != 1)
		return -1;
	return rv;
}

/* Do a full featured command/response echange with the dtp41 */
/* End on the specified number of characters, or expiry if */
/* the specified timeout. */
/* Assume standard error code if tc = '>' and ntc = 1 */
/* Return a DTP41 error code */
static int
dtp41_fcommand(
dtp41 *p,
char *in,			/* In string */
char *out,			/* Out string buffer */
int bsize,			/* Out buffer size */
char tc,			/* Terminating character */
int ntc,			/* Number of terminating characters */
double to) {		/* Timout in seconts */
	int rv, se;

	if ((se = p->sio->write_read(p->sio, in, out, bsize, tc, ntc, to)) != 0) {
#ifdef DEBUG
	printf("dtp41 fcommand: serial i/o failure on write_read '%s'\n",fix(in));
#endif
		if (se & SIO_USER)
			return DTP41_USER_ABORT;
		return DTP41_COMS_FAIL;
	}
	rv = DTP41_OK;
	if (tc == '>' && ntc == 1) {
		rv = extract_ec(out);
		if (rv > 0) {
			rv &= 0x7f;
			if (rv != DTP41_OK)	/* Clear the error */
				p->sio->write(p->sio, "CE\r", 0.5);
		}
	}
#ifdef DEBUG
	printf("command '%s'",fix(in));
	printf(" returned '%s', value 0x%x\n",fix(out),rv);
#endif
	return rv;
}

/* Do a standard command/response echange with the dtp41 */
/* Return the instrument error code */
static inst_code
dtp41_command(inst *pp, char *in, char *out, int bsize) {
	dtp41 *p = (dtp41 *)pp;
	int rv = dtp41_fcommand(p, in, out, bsize, '>', 1, 1.5);
	return interp_code(pp, rv);
}

/* Establish communications with a DTP41 */
/* Use the baud rate given, and timeout in to secs */
/* Return DTP_COMS_FAIL on failure to establish communications */
static inst_code
dtp41_init_coms(inst *pp, int port, baud_rate br, double tout) {
	dtp41 *p = (dtp41 *)pp;
	static char buf[MAX_MES_SIZE];
	baud_rate brt[9] = { baud_9600, baud_19200, baud_38400, baud_57600,
	                     baud_4800, baud_2400, baud_1200, baud_600, baud_300 };
	char *brc[9] =     { "9600BR\r", "19200BR\r", "38400BR\r", "57600BR\r",
	                     "4800BR\r", "2400BR\r", "1200BR\r", "600BR\r", "300BR\r" };
	int bi;
	long etime;
	int i, rv;
	inst_code ev = inst_ok;

	if (p->debug)
		p->sio->debug = p->debug;	/* Turn on debugging */

	/* Figure DTP41 baud rate being asked for */
	for (bi = 0; bi < 9; bi++) {
		if (brt[bi] == br)
			break;
	}
	if (bi >= 9)
		bi = 0;	

	/* The tick to give up on */
	etime = clock() + (long)(CLOCKS_PER_SEC * tout + 0.5);

	while (clock() < etime) {

		/* Until we time out, find the correct baud rate */
		for (i=0; clock() < etime;) {
			p->sio->set_port(p->sio, NULL, port, brt[i], parity_none, stop_1, length_8);
			if (((ev = interp_code(pp, dtp41_fcommand(p, "\r", buf, MAX_MES_SIZE, '>', 1, 0.5)))
			         & inst_mask) != inst_coms_fail)
				break;		/* We've got coms or user abort */
			if (++i >= 9)
				i = 0;
		}

		if ((ev & inst_mask) == inst_user_abort)
			return ev;

		break;		/* Got coms */
	}

	if (clock() >= etime) {		/* We haven't established comms */
		return inst_coms_fail;
	}

	/* set the protocol to RCI */
	if (((ev = p->command((inst *)p, "0012CF\r", buf, 500)) & inst_mask) != inst_ok)
		return ev;

	/* Disable handshaking */
	if (((ev = p->command((inst *)p, "0004CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	if (i != bi) {
		/* Change the baud rate to the rate we've been told */
		if ((rv = p->sio->write_read(p->sio, brc[bi], buf, MAX_MES_SIZE, '>', 1, 1.5)) != 0) {
			if (extract_ec(buf) != DTP41_OK)
				return inst_coms_fail;
		}
	
		/* Configure our baud rate as well */
		p->sio->set_port(p->sio, NULL, port, brt[bi], parity_none, stop_1, length_8);
	
		/* Loose a character - not sure why. */
		p->sio->write_read(p->sio, "\r", buf, MAX_MES_SIZE, '>', 1, 0.5);

		/* Check instrument is responding */
		if (((ev = p->command((inst *)p, "\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
			return inst_coms_fail;
	}

#ifdef DEBUG
	printf("Got communications\n");
#endif
	p->gotcoms = 1;
	return inst_ok;
}

/* Build a strip definition as a set of passes, including DS command */
static void
build_strip(
dtp41 *p,
char *tp,			/* pointer to string buffer */
char *name,			/* Strip name (7 chars) (not used) */
int npatch,			/* Number of patches in the pass */
char *pname,		/* Pass name (3 chars) (not used) */
int sguide,			/* Guide number (not used) */
double pwid,		/* Patch length in mm */
double gwid,		/* Gap length in mm */
double twid			/* Trailer length in mm (DTP41T only) */
) {
	int i;

	/* Number of patches in strip */
	sprintf(tp, "%03d",npatch);
	tp += 3;

	/* Patch width in mm, dd.dd */
	sprintf(tp, "%05.2f",pwid);
	tp[2] = tp[3];	/* Remove point */
	tp[3] = tp[4];
	tp += 4;

	/* Gap width in mm, dd.dd */
	sprintf(tp, "%05.2f",gwid);
	tp[2] = tp[3];	/* Remove point */
	tp[3] = tp[4];
	tp += 4;

	*tp++ = '0';	/* Normal strip */

	*tp++ = '8';	/* Auto type */

	if (p->mode & inst_mode_transmission) {

		if (twid >= 9999.5)
			error ("build_strip given trailer length > 9999 mm");
		/* Trailer length in mm, dddd */
		sprintf(tp, "%04.0f",twid);
		tp += 4;
	}

	*tp++ = 'D';				/* The DS command */
	*tp++ = 'S';
	*tp++ = '\r';				/* The CR */
	*tp++ = '\000';				/* The end */

}
		
/* Initialise the DTP41. */
/* return non-zero on an error, with instrument error code */
static inst_code
dtp41_init_inst(inst *pp) {
	dtp41 *p = (dtp41 *)pp;
	static char tbuf[100], buf[MAX_MES_SIZE];
	int i;
	inst_code rv = inst_ok;

	if (p->gotcoms == 0)
		return inst_internal_error;		/* Must establish coms before calling init */

	/* Resetting instrument resets the baud rate, so do manual reset. */

	/* Set emulation mode to DTP41 */
	if (((rv = p->command((inst *)p, "0010CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Turn echoing of characters off */
	if (((rv = p->command((inst *)p, "0009CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set Response delimeter to CR */
	if (((rv = p->command((inst *)p, "0008CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Get the model and version number */
	if (((rv = p->command((inst *)p, "SV\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Check that it is a DTP41 */
	if (   strlen(buf) < 12
	    || strncmp(buf,"X-Rite DTP41",11) != 0
	    || (buf[11] != '1' && buf[11] != '2'))
		return inst_unknown_model;

	/* Set Language to English */
	if (((rv = p->command((inst *)p, "0000CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set Beeper to medium */
	if (((rv = p->command((inst *)p, "0201CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set Automatic Transmit off */
	if (((rv = p->command((inst *)p, "0005CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set decimal point on */
	if (((rv = p->command((inst *)p, "0106CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set color data separator to TAB */
	if (((rv = p->command((inst *)p, "0207CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set 2 decimal digit resolution */
	if (((rv = p->command((inst *)p, "020ACF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set Min/Max mode off */
	if (((rv = p->command((inst *)p, "000CCF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set persistent errors off */
	if (((rv = p->command((inst *)p, "000DCF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set show data labels mode off */
	if (((rv = p->command((inst *)p, "000FCF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set drive motor calibration at power up to off */
	if (((rv = p->command((inst *)p, "0011CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set Reflection calibration timeout to 24 Hrs */
	if (((rv = p->command((inst *)p, "181ECF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set Trailer timout to 2 seconds */
	if (((rv = p->command((inst *)p, "021FCF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set Transmission calibration timeout to 24 Hrs */
	if (((rv = p->command((inst *)p, "1820CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* - - - - - - - - - - - - - - - - - - - - - - - - */
	/* Setup for the type of measurements we want to do */
	/* Enable the read microswitch */
	if (((rv = p->command((inst *)p, "01PB\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set dynamic measurement mode */
	if (((rv = p->command((inst *)p, "0113CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set instrument to reflectance mode */
	if (((rv = p->command((inst *)p, "0019CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set data format to Reflectance, so TS can select. */
	if (((rv = p->command((inst *)p, "0318CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set density format to spectral */
	if (((rv = p->command((inst *)p, "0417CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set Illuminant to D50_2 */
	if (((rv = p->command((inst *)p, "0416CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set static samples to 10 */
	if (((rv = p->command((inst *)p, "0A14CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* Set static readings to configured number (usually 5) */
	sprintf(tbuf, "%02x15CF\r", p->nstaticr);
	if (((rv = p->command((inst *)p, tbuf, buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return rv;

	/* We are configured in this mode now */
	p->mode = inst_mode_ref_strip;

	if (p->lastmode != p->mode)
		return activate_mode(p);

	if (rv == inst_ok)
		p->inited = 1;

	return inst_ok;
}

/* For an xy instrument, release the paper */
/* Return the inst error code */
static inst_code
dtp41_xy_sheet_release(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, hold the paper */
/* Return the inst error code */
static inst_code 
dtp41_xy_sheet_hold(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, allow the user to locate a point */
/* Return the inst error code */
static inst_code
dtp41_xy_locate_start(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, read back the location */
/* Return the inst error code */
static inst_code
dtp41_xy_get_location(
struct _inst *p,
double *x, double *y) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, ends allowing the user to locate a point */
/* Return the inst error code */
static inst_code
dtp41_xy_locate_end(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, try and clear the table after an abort */
/* Return the inst error code */
static inst_code
dtp41_xy_clear(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* Read a sheet full of patches using xy mode */
/* Return the inst error code */
static inst_code
dtp41_read_xy(
inst *p,
int pis,				/* Passes in strip (letters in sheet) */
int sip,				/* Steps in pass (numbers in sheet) */
int npatch,				/* Total patches in strip (skip in last pass) */
char *pname,			/* Starting pass name (' A' to 'ZZ') */
char *sname,			/* Starting step name (' 1' to '99') */
double ox, double oy,	/* Origin of first patch */
double ax, double ay,	/* pass increment */
double aax, double aay,	/* pass offset for odd patches */
double px, double py,	/* step/patch increment */
ipatch *vals) { 		/* Pointer to array of values */
	/* This is not supported by this device */
	return inst_unsupported;
}

/* Read a set of strips */
/* Return the instrument error code */
static inst_code
dtp41_read_strip(
inst *pp,
char *name,			/* Strip name (7 chars) */
int npatch,			/* Number of patches in the pass */
char *pname,		/* Pass name (3 chars) */
int sguide,			/* Guide number */
double pwid,		/* Patch length in mm (For DTP41) */
double gwid,		/* Gap length in mm (For DTP41) */
double twid,		/* Trailer length in mm (For DTP41T) */
ipatch *vals) {		/* Pointer to array of instrument patch values */
	dtp41 *p = (dtp41 *)pp;
	char tbuf[200], *tp;
	static char buf[MAX_RD_SIZE];
	int i, rv = 0;
	inst_code ev = inst_ok;

	/* Configure for dynamic mode */
	p->lastmode = (p->lastmode & ~inst_mode_sub_mask) | inst_mode_strip;
	activate_mode(p);

	build_strip(p, tbuf, name, npatch, pname, sguide, pwid, gwid, twid);
	
	/* Send strip definition */
	if (((ev = p->command((inst *)p, tbuf, buf, 500)) & inst_mask) != inst_ok)
		return ev;

#ifdef NEVER
	/* Do a strip read */
	if (((ev = p->command((inst *)p, "RM\r", buf, 500)) & inst_mask) != inst_ok)
		return ev;
#endif

	/* Wait for the Read status, or a user abort - allow 5 minuits */
	rv = dtp41_fcommand(p, "", buf, MAX_RD_SIZE, '>', 1, 5 * 60.0);

	if (rv != DTP41_OK) {
		return interp_code(pp, rv); 	/* Strip misread */
	}

	/* Gather the results in D50_2 XYZ */
	if ((rv = dtp41_fcommand(p, "0405TS\r", buf, MAX_RD_SIZE, '>', 1, 0.5 + npatch * 0.1)) != DTP41_OK)
		return interp_code(pp, rv); 	/* Strip misread */

	/* Parse the buffer */
	/* Replace '\r' with '\000' */
	for (tp = buf; *tp != '\000'; tp++) {
		if (*tp == '\r')
			*tp = '\000';
	}
	for (tp = buf, i = 0; i < npatch; i++) {
		if (*tp == '\000')
			return inst_protocol_error;
		if (sscanf(tp, " %lf %lf %lf ",
		           &vals[i].XYZ[0], &vals[i].XYZ[1], &vals[i].XYZ[2]) != 3) {
			if (sscanf(tp, " %lf %lf %lf ",
			           &vals[i].XYZ[0], &vals[i].XYZ[1], &vals[i].XYZ[2]) != 3) {
				return inst_protocol_error;
			}
		}
		vals[i].XYZ_v = 1;
		tp += strlen(tp) + 1;
	}

	if (p->mode & inst_mode_spectral) {

		/* Gather the results in Spectral reflectance */
		if ((rv = dtp41_fcommand(p, "0403TS\r", buf, MAX_RD_SIZE, '>', 1, 0.5 + npatch * 0.1)) != DTP41_OK)
			return interp_code(pp, rv); 	/* Strip misread */
	
		/* Parse the buffer */
		/* Replace '\r' with '\000' */
		for (tp = buf; *tp != '\000'; tp++) {
			if (*tp == '\r')
				*tp = '\000';
		}
		/* Get each patches spetra */
		for (tp = buf, i = 0; i < npatch; i++) {
			int cnp, j;
			char *tpp;
			if (strlen(tp) < (31 * 8 - 1)) {
				return inst_protocol_error;
			}

			/* Read the spectral value */
			for (tpp = tp, j = 0; j < 31; j++, tpp += 8) {
				char c;
				c = tpp[7];
				tpp[7] = '\000';
				vals[i].spec[j] = atof(tpp);
				tpp[7] = c;
			}

			vals[i].spec_n = 31;
			vals[i].spec_wl_short = 400.0;
			vals[i].spec_wl_long = 700.0;
			tp += strlen(tp) + 1;
		}
	}

	return inst_ok;
}


/* Read a single sample */
/* Return the instrument error code */
static inst_code
dtp41_read_sample(
inst *pp,
char *name,			/* Strip name (7 chars) */
ipatch *val) {		/* Pointer to instrument patch value */
	dtp41 *p = (dtp41 *)pp;
	char tbuf[200], *tp;
	static char buf[MAX_RD_SIZE];
	int i, rv = 0;
	inst_code ev = inst_ok;

	/* Configure for static mode */
	p->lastmode = (p->lastmode & ~inst_mode_sub_mask) | inst_mode_spot;
	activate_mode(p);

	/* Set static measurement mode */
	if (((ev = p->command((inst *)p, "0013CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Do a spot read */
	if ((rv = dtp41_fcommand(p, "RM\r", buf, MAX_RD_SIZE, '>', 1, 5 * 60.0)) != DTP41_OK) {
		return interp_code(pp, rv); 	/* Strip misread */
	}

#ifdef NEVER
	/* Wait for the Read status, or a user abort */
	rv = dtp41_fcommand(p, "", buf, MAX_RD_SIZE, '>', 1, 5 * 60.0);

	if (rv != DTP41_OK) {
		return interp_code(pp, rv); 	/* Strip misread */
	}
#endif

	/* Gather the results in D50_2 XYZ */
	if ((rv = dtp41_fcommand(p, "0405TS\r", buf, MAX_RD_SIZE, '>', 1, 0.5)) != DTP41_OK)
		return interp_code(pp, rv); 	/* Strip misread */

	/* Parse the buffer */
	/* Replace '\r' with '\000' */
	for (tp = buf; *tp != '\000'; tp++) {
		if (*tp == '\r')
			*tp = '\000';
	}
	
	val->XYZ[0] = val->XYZ[1] = val->XYZ[2] = 0.0;

	/* for all the readings taken */
	for (tp = buf, i = 0; i < p->nstaticr; i++) {
		double XYZ[3];

		if (*tp == '\000')
			return inst_protocol_error;

		if (sscanf(tp, " %lf %lf %lf ",
		           &XYZ[0], &XYZ[1], &XYZ[2]) != 3) {
			if (sscanf(tp, " %lf %lf %lf ",
			           &XYZ[0], &XYZ[1], &XYZ[2]) != 3) {
				return inst_protocol_error;
			}
		}
		val->XYZ[0] += XYZ[0];
		val->XYZ[1] += XYZ[1];
		val->XYZ[2] += XYZ[2];
		tp += strlen(tp) + 1;
	}

	/* Average */
	val->XYZ[0] /= (double)p->nstaticr;
	val->XYZ[1] /= (double)p->nstaticr;
	val->XYZ[2] /= (double)p->nstaticr;
	val->XYZ_v = 1;

	if (p->mode & inst_mode_spectral) {
		int j;

		/* Gather the results in Spectral reflectance */
		if ((rv = dtp41_fcommand(p, "0403TS\r", buf, MAX_RD_SIZE, '>', 1, 0.5)) != DTP41_OK)
			return interp_code(pp, rv); 	/* Strip misread */
	
		/* Parse the buffer */
		/* Replace '\r' with '\000' */
		for (tp = buf; *tp != '\000'; tp++) {
			if (*tp == '\r')
				*tp = '\000';
		}

		for (j = 0; j < 31; j++)
			val->spec[j] = 0.0;

		/* Get each readings spetra */
		for (tp = buf, i = 0; i < p->nstaticr; i++) {
			int cnp;
			char *tpp;
			if (strlen(tp) < (31 * 8 - 1)) {
				return inst_protocol_error;
			}

			/* Read the spectral value */
			for (tpp = tp, j = 0; j < 31; j++, tpp += 8) {
				char c;
				c = tpp[7];
				tpp[7] = '\000';
				val->spec[j] += atof(tpp);
				tpp[7] = c;
			}

			tp += strlen(tp) + 1;
		}

		/* Average the result */
		for (j = 0; j < 31; j++)
			val->spec[j] /= (double)p->nstaticr;

		val->spec_n = 31;
		val->spec_wl_short = 400.0;
		val->spec_wl_long = 700.0;
	}

	/* Set back to dynamic measurement mode */
	if (((ev = p->command((inst *)p, "0113CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	return inst_ok;
}


/* Determine if a calibration is needed. Returns inst_ok if not, */
/* inst_unsupported if it is unknown, or inst_needs_cal if needs calibration */
inst_code dtp41_needs_calibration(struct _inst *p) {
	/* This can't be known ? */
	return inst_unsupported;
}

/* Perform an instrument calibration */
inst_code dtp41_calibrate(
struct _inst *p,
int cal_no
) {
	/* This is not currently supported */
	return inst_unsupported;

}

/* Error codes interpretation */
static char *
dtp41_interp_error(inst *pp, int ec) {
	dtp41 *p = (dtp41 *)pp;
	ec &= 0x7f;
	switch (ec) {

		case DTP41_USER_ABORT:
			return "User hit Escape";
		case DTP41_INTERNAL_ERROR:
			return "Internal software error";
		case DTP41_COMS_FAIL:
			return "Communications failure";
		case DTP41_UNKNOWN_MODEL:
			return "Not a DTP41";
		case DTP41_DATA_PARSE_ERROR:
			return "Data from DTP41 didn't parse as expected";
		case DTP41_OK:
			return "No error";
		case DTP41_MEASUREMENT_STATUS:
			return "Measurement complete";
		case DTP41_CALIBRATION_STATUS:
			return "Calibration complete";
		case DTP41_KEYPRESS_STATUS:
			return "A key was pressed";
		case DTP41_DEFAULTS_LOADED_STATUS:
			return "Default configuration values have been loaded";
		case DTP41_BAD_COMMAND:
			return "Unrecognised command";
		case DTP41_BAD_PARAMETERS:
			return "Wrong number of parameters";
		case DTP41_PRM_RANGE_ERROR:
			return "One or more parameters are out of range";
		case DTP41_BUSY:
			return "Instrument is busy - command ignored";
		case DTP41_USER_ABORT_ERROR:
			return "User aborted process";
		case DTP41_MEASUREMENT_ERROR:
			return "General measurement error";
		case DTP41_TIMEOUT:
			return "Receive timeout";
		case DTP41_BAD_STRIP:
			return "Bad strip";
		case DTP41_BAD_COLOR:
			return "Bad color";
		case DTP41_BAD_STEP:
			return "Bad step";
		case DTP41_BAD_PASS:
			return "Bad pass";
		case DTP41_BAD_PATCHES:
			return "Bad patches";
		case DTP41_BAD_READING:
			return "Bad reading";
		case DTP41_NEEDS_CAL_ERROR:
			return "Instrument needs calibration";
		case DTP41_CAL_FAILURE_ERROR:
			return "Calibration failed";
		case DTP41_INSTRUMENT_ERROR:
			return "General instrument error";
		case DTP41_LAMP_ERROR:
			return "Reflectance lamp error";
		case DTP41_FILTER_ERROR:
			return "Filter error";
		case DTP41_FILTER_MOTOR_ERROR:
			return "Filter motor error";
		case DTP41_DRIVE_MOTOR_ERROR:
			return "Strip drive motor error";
		case DTP41_KEYPAD_ERROR:
			return "Keypad error";
		case DTP41_DISPLAY_ERROR:
			return "Display error";
		case DTP41_MEMORY_ERROR:
			return "Memory error";
		case DTP41_ADC_ERROR:
			return "ADC error";
		case DTP41_PROCESSOR_ERROR:
			return "Processor error";
		case DTP41_BATTERY_ERROR:
			return "Battery error";
		case DTP41_BATTERY_LOW_ERROR:
			return "Battery low error";
		case DTP41_INPUT_POWER_ERROR:
			return "Input power error";
		case DTP41_TEMPERATURE_ERROR:
			return "Temperature error";
		case DTP41_BATTERY_ABSENT_ERROR:
			return "Battery absent error";
		case DTP41_TRAN_LAMP_ERROR:
			return "Transmission lamp error";
		case DTP41_INVALID_COMMAND_ERROR:
			return "Invalid command";
		default:
			return "Unknown error code";
	}
}

/* Convert a machine specific error code into an abstract dtp code */
static inst_code 
interp_code(inst *pp, int ec) {
	dtp41 *p = (dtp41 *)pp;

	ec &= 0x7f;
	switch (ec) {

		case DTP41_OK:
			return inst_ok;

		case DTP41_INTERNAL_ERROR:
			return inst_internal_error | ec;

		case DTP41_COMS_FAIL:
			return inst_coms_fail | ec;

		case DTP41_UNKNOWN_MODEL:
			return inst_unknown_model | ec;

		case DTP41_DATA_PARSE_ERROR:
			return inst_protocol_error | ec;

		case DTP41_USER_ABORT:
			return inst_user_abort | ec;

		case DTP41_BUSY:
		case DTP41_TIMEOUT:
		case DTP41_BAD_READING:
		case DTP41_INSTRUMENT_ERROR:
		case DTP41_DRIVE_MOTOR_ERROR:
		case DTP41_ADC_ERROR:
		case DTP41_TRAN_LAMP_ERROR:
			return inst_misread | ec;

		case DTP41_NEEDS_CAL_ERROR:
			return inst_needs_cal | ec;
	}
	return inst_other_error | ec;
}

/* Return the last serial I/O error code */
static int
dtp41_last_sioerr(inst *pp) {
	dtp41 *p = (dtp41 *)pp;
	return p->sio->lerr;
}

/* Destroy ourselves */
static void
dtp41_del(inst *pp) {
	dtp41 *p = (dtp41 *)pp;
	if (p->sio != NULL)
		p->sio->del(p->sio);
	free (p);
}

/* Interogate the device to discover its capabilities */
static void	discover_capabilities(dtp41 *p) {
	static char buf[MAX_MES_SIZE];
	inst_code rv = inst_ok;

	p->cap = inst_ref_spot
	       | inst_ref_strip
	       | inst_colorimeter
	       | inst_spectral
	         ;

	/* Check whether we have transmission capability */
	if (((rv = p->command((inst *)p, "0119CF\r", buf, MAX_MES_SIZE)) & inst_mask) == inst_ok) {
		p->cap |= inst_trans_spot | inst_trans_strip;
	}
	rv = p->command((inst *)p, "0019CF\r", buf, MAX_MES_SIZE);

}

/* Return the instrument capabilities */
inst_capability dtp41_capabilities(inst *pp) {
	dtp41 *p = (dtp41 *)pp;

	if (p->cap == inst_unknown)
		discover_capabilities(p);
	return p->cap;
}

/* Activate the last set mode */
static inst_code
activate_mode(dtp41 *p)
{
	static char buf[MAX_MES_SIZE];
	inst_code rv;

	/* Setup for transmission or reflection */
	if ((p->lastmode & inst_mode_illum_mask) == inst_mode_reflection
	 && (p->mode     & inst_mode_illum_mask) != inst_mode_reflection) {
		if (((rv = p->command((inst *)p, "0019CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
			return rv;
	}
	if ((p->lastmode & inst_mode_illum_mask) == inst_mode_transmission
	 && (p->mode     & inst_mode_illum_mask) != inst_mode_transmission) {
		if (((rv = p->command((inst *)p, "0119CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
			return rv;
	}

	/* Setup for static or dynamic reading */
	if ((p->lastmode & inst_mode_sub_mask) == inst_mode_spot
	 && (p->mode     & inst_mode_sub_mask) != inst_mode_spot) {
		if (((rv = p->command((inst *)p, "0013CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
			return rv;
	}
	if ((p->lastmode & inst_mode_sub_mask) == inst_mode_strip
	 && (p->mode     & inst_mode_sub_mask) != inst_mode_strip) {
		if (((rv = p->command((inst *)p, "0113CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
			return rv;
	}

	p->mode = p->lastmode;

	return inst_ok;
}

/* 
 * set measurement mode
 * We assume that the instrument has been initialised.
 */
static inst_code
dtp41_set_mode(inst *pp, inst_mode m)
{
	dtp41 *p = (dtp41 *)pp;
	inst_capability cap = pp->capabilities(pp);
	inst_mode mm;		/* Measurement mode */

	/* The measurement mode portion of the mode */
	mm = m & inst_mode_measuremet_mask;

	/* General check mode against specific capabilities logic: */
	if (mm == inst_mode_ref_spot) {
 		if (!(cap & inst_ref_spot))
			return inst_unsupported;
	} else if (mm == inst_mode_ref_strip) {
		if (!(cap & inst_ref_strip))
			return inst_unsupported;
	} else if (mm == inst_mode_ref_xy) {
		if (!(cap & inst_ref_xy))
			return inst_unsupported;
	} else if (mm == inst_mode_trans_spot) {
		if (!(cap & inst_trans_spot))
			return inst_unsupported;
	} else if (mm == inst_mode_trans_strip) {
		if (!(cap & inst_trans_strip))
			return inst_unsupported;
	} else if (mm == inst_mode_trans_xy) {
		if (!(cap & inst_trans_xy))
			return inst_unsupported;
	} else if (mm == inst_mode_emis_disp) {
		if (!(cap & inst_emis_disp))
			return inst_unsupported;
	} else if (mm == inst_mode_emis_illum) {
		if (!(cap & inst_emis_illum))
			return inst_unsupported;
	} else {
		return inst_unsupported;
	}

	if (m & inst_mode_colorimeter)
		if (!(cap & inst_colorimeter))
			return inst_unsupported;
		
	if (m & inst_mode_spectral)
		if (!(cap & inst_spectral))
			return inst_unsupported;

	p->lastmode = m;

	if (p->lastmode != p->mode) {
		return activate_mode(p);
	}
	return inst_ok;
}

/* 
 * set or reset an optional mode
 * We assume that the instrument has been initialised.
 */
static inst_code
dtp41_set_opt_mode(inst *pp, inst_opt_mode m)
{
	dtp41 *p = (dtp41 *)pp;
	return inst_unsupported;
}

/* Constructor */
extern dtp41 *new_dtp41()
{
	dtp41 *p;
	if ((p = (dtp41 *)calloc(sizeof(dtp41),1)) == NULL)
		error("dtp41: malloc failed!");

	p->sio = new_serio();

	p->init_coms    = dtp41_init_coms;
	p->init_inst    = dtp41_init_inst;
	p->capabilities = dtp41_capabilities;
	p->set_mode     = dtp41_set_mode;
	p->set_opt_mode     = dtp41_set_opt_mode;
	p->xy_sheet_release = dtp41_xy_sheet_release;
	p->xy_sheet_hold    = dtp41_xy_sheet_hold;
	p->xy_locate_start  = dtp41_xy_locate_start;
	p->xy_get_location  = dtp41_xy_get_location;
	p->xy_locate_end  = dtp41_xy_locate_end;
	p->xy_clear     = dtp41_xy_clear;
	p->read_xy      = dtp41_read_xy;
	p->read_strip   = dtp41_read_strip;
	p->read_sample  = dtp41_read_sample;
	p->needs_calibration = dtp41_needs_calibration;
	p->calibrate    = dtp41_calibrate;
	p->command      = dtp41_command;
	p->interp_error = dtp41_interp_error;
	p->inst_interp_error = NULL;				/* virtual constructor will do this */
	p->last_sioerr  = dtp41_last_sioerr;
	p->del          = dtp41_del;

	p->itype = instDTP41;
	p->cap = inst_unknown;						/* Unknown until initialised */
	p->mode = inst_mode_unknown;				/* Not in a known mode yet */
	p->nstaticr = 5;							/* Number of static readings */

	return p;
}
