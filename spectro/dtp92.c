/* 
 * Argyll Color Correction System
 *
 * Xrite DTP92 related functions
 *
 * Author: Graeme W. Gill
 * Date:   5/10/96
 *
 * Copyright 1996, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include "xspect.h"
#include "serio.h"
#include "dtp92.h"

#undef DEBUG

#include <stdarg.h>
void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);

static inst_code intep_code(inst *pp, int ec);

#if defined (NT)
# define WIN32_LEAN_AND_MEAN
# include <windows.h>
#define sleep(secs) Sleep((secs) * 1000)
#endif

#define MAX_MES_SIZE 500		/* Maximum normal message reply size */
#define MAX_RD_SIZE 5000		/* Maximum reading messagle reply size */

static void
delay (double tt) {
	long etime;
	etime = clock() + (long)(CLOCKS_PER_SEC * tt + 0.5);
	while (clock() < etime);
}

/* Convert control chars to ^[A-Z] notation in a string */
static char *
fix(char *s) {
	static char buf [5000];
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

/* Do a full featured command/response echange with the dtp92 */
/* Return the dtp error code. End on the specified number */
/* of specified characters, or expiry if the specified timeout */
/* Assume standard error code if tc = '>' and ntc = 1 */
/* Return a DTP92 error code */
static int
dtp92_fcommand(
	struct _dtp92 *p,
	char *in,			/* In string */
	char *out,			/* Out string buffer */
	int bsize,			/* Out buffer size */
	char tc,			/* Terminating character */
	int ntc,			/* Number of terminating characters */
	double to) {		/* Timout in seconds */
	int rv, se;

	if ((se = p->sio->write_read(p->sio, in, out, bsize, tc, ntc, to)) != 0) {
#ifdef DEBUG
	printf("dtp92 fcommand: serial i/o failure on write_read '%s'\n",fix(in));
#endif
		if (se & SIO_USER)
			return DTP92_USER_ABORT;
		return DTP92_COMS_FAIL;
	}
	rv = DTP92_OK;
	if (tc == '>' && ntc == 1) {
		rv = extract_ec(out);
		if (rv > 0) {
			rv &= 0x7f;
			if (rv != DTP92_OK)	/* Clear the error */
				p->sio->write(p->sio, "CE\r", 0.5);
		}
	}
#ifdef DEBUG
	printf("command '%s'",fix(in));
	printf(" returned '%s', value 0x%x\n",fix(out),rv);
#endif
	return rv;
}

/* Do a standard command/response echange with the dtp92 */
/* Return the dtp error code */
static inst_code
dtp92_command(inst *pp, char *in, char *out, int bsize) {
	dtp92 *p = (dtp92 *)pp;
	int rv = dtp92_fcommand(p, in, out, bsize, '>', 1, 0.5);
	return intep_code(pp, rv);
}

/* Establish communications with a DTP92 */
/* Use the baud rate given, and timeout in to secs */
/* Return DTP_COMS_FAIL on failure to establish communications */
static inst_code
dtp92_init_coms(inst *pp, int port, baud_rate br, double tout) {
	dtp92 *p = (dtp92 *)pp;
	static char buf[MAX_MES_SIZE];
	baud_rate brt[5] = { baud_9600, baud_19200, baud_4800, baud_2400, baud_1200 };
	char *brc[5]     = { "30BR\r",  "60BR\r",   "18BR\r",  "0CBR\r",  "06BR\r" };
	int bi;
	long etime;
	int i, rv;
	inst_code ev = inst_ok;

	if (p->debug)
		p->sio->debug = p->debug;	/* Turn on debugging */

	/* Figure DTP92 baud rate being asked for */
	for (bi = 0; bi < 5; bi++) {
		if (brt[bi] == br)
			break;
	}
	if (bi >= 5)
		bi = 0;	

	/* The tick to give up on */
	etime = clock() + (long)(CLOCKS_PER_SEC * tout + 0.5);

	while (clock() < etime) {

		/* Until we time out, find the correct baud rate */
		for (i=0;clock() < etime;) {
			p->sio->set_port(p->sio, NULL, port, brt[i], parity_none, stop_1, length_8);
			if (((ev = p->command((inst *)p, "\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_coms_fail) {
				break;		/* We've got coms or user abort */
			}
			if (++i >= 5)
				i = 0;
		}

		if ((ev & inst_mask) == inst_user_abort)
			return ev;

		break;		/* Got coms */
	}

	if (clock() >= etime) {		/* We haven't established comms */
		return inst_coms_fail;
	}

	/* Disable handshaking */
	if (((ev = p->command((inst *)p, "0004CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok) {
#ifdef DEBUG
		printf("Disabling handshaking failed\n");
#endif
		return ev;
	}

	if (i != bi) {
		/* Change the baud rate to the rate we've been told */
		if ((rv = p->sio->write_read(p->sio, brc[bi], buf, MAX_MES_SIZE, '>', 1, 0.5)) == 0) {
			if (extract_ec(buf) != DTP92_OK) {
#ifdef DEBUG
				printf("Changing baud rate failed\n");
#endif
				return inst_coms_fail;
			}
		}

		/* Configure our baud rate as well */
		p->sio->set_port(p->sio, NULL, port, brt[bi], parity_none, stop_1, length_8);

		/* Loose a character - just in case. */
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

/* Initialise the DTP92 */
/* return non-zero on an error, with dtp error code */
static inst_code
dtp92_init_inst(inst *pp) {
	dtp92 *p = (dtp92 *)pp;
	static char tbuf[100], buf[MAX_MES_SIZE];
	int i;
	inst_code ev = inst_ok;

	if (p->gotcoms == 0)
		return inst_internal_error;		/* Must establish coms before calling init */

	/* Reset it */
	if (((ev = dtp92_fcommand(p, "0PR\r", buf, MAX_MES_SIZE, '>', 1, 4.0)) & inst_mask) != inst_ok)
		return ev;

	/* Turn echoing of characters off */
	if (((ev = p->command((inst *)p, "DEC\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Get the model and version number */
	if (((ev = p->command((inst *)p, "SV\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Check that it is a DTP92 or DTP92Q */
	if (   strlen(buf) < 12
	    || strncmp(buf,"X-Rite DTP92",12) != 0
	    || (strlen(buf) > 12 && buf[12] != ' ' && buf[12] != 'Q'))
		return inst_unknown_model;

	/* Set decimal point on */
	if (((ev = p->command((inst *)p, "0106CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set color data separator to TAB */
	if (((ev = p->command((inst *)p, "0207CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set delimeter to CR */
	if (((ev = p->command((inst *)p, "0008CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set extra digit resolution (X10) */
	if (((ev = p->command((inst *)p, "010ACF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set absolute (luminance) calibration */
	if (((ev = p->command((inst *)p, "0118CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set no black point subtraction */
	if (((ev = p->command((inst *)p, "0019CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set to factory calibration */
	if (((ev = p->command((inst *)p, "EFC\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set to transmit just the colorimetric data */
	if (((ev = p->command((inst *)p, "0120CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Transmit colorimetric data as XYZ */
	if (((ev = p->command((inst *)p, "0221CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Disable Luminance group (~1) */
	if (((ev = p->command((inst *)p, "0022CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Disable Frequency group */
	if (((ev = p->command((inst *)p, "0023CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Disable color temperature group */
	if (((ev = p->command((inst *)p, "0024CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Disable RGB group */
	if (((ev = p->command((inst *)p, "0025CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Disable pushbutton */
	if (((ev = p->command((inst *)p, "DPB\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set sample size to 16 (default is 16) */
	if (((ev = p->command((inst *)p, "10SS\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	if (ev == inst_ok)
		p->inited = 1;

	return ev;
}

/* For an xy instrument, release the paper */
/* Return the inst error code */
static inst_code
dtp92_xy_sheet_release(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, hold the paper */
/* Return the inst error code */
static inst_code 
dtp92_xy_sheet_hold(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, allow the user to locate a point */
/* Return the inst error code */
static inst_code
dtp92_xy_locate_start(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, read back the location */
/* Return the inst error code */
static inst_code
dtp92_xy_get_location(
struct _inst *p,
double *x, double *y) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, ends allowing the user to locate a point */
/* Return the inst error code */
static inst_code
dtp92_xy_locate_end(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, try and clear the table after an abort */
/* Return the inst error code */
static inst_code
dtp92_xy_clear(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* Read a sheet full of patches using xy mode */
/* Return the inst error code */
static inst_code
dtp92_read_xy(
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
/* Return the dtp error code */
static inst_code
dtp92_read_strip(
inst *pp,
char *name,			/* Strip name (7 chars) */
int npatch,			/* Number of patches in the pass */
char *pname,		/* Pass name (3 chars) */
int sguide,			/* Guide number */
double pwid,		/* Patch length in mm (DTP41) */
double gwid,		/* Gap length in mm (DTP41) */
double twid,		/* Trailer length in mm (DTP41T) */
ipatch *vals) {		/* Pointer to array of instrument patch values */
	dtp92 *p = (dtp92 *)pp;

	/* This is not supported by DTP92 */
	return inst_unsupported;
}

/* The DTP92 seems to have a bug whereby it adds a spurious */
/* digit after the 'Z' of the Z value. Try and discard this. */

/* Read a single sample */
/* Return the dtp error code */
static inst_code
dtp92_read_sample(
inst *pp,
char *name,			/* Strip name (7 chars) */
ipatch *val) {		/* Pointer to instrument patch value */
	dtp92 *p = (dtp92 *)pp;
	static char buf[MAX_RD_SIZE];
	int tries;
	int rv = inst_protocol_error;

	/* Could change SS to suite level expected. */
#ifdef NEVER
		/* Set sample size to 31 (default is 16) for low level readings */
		if (((rv = p->command((inst *)p, "1fSS\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
			return rv;
#endif


	/* Until we get a valid return */
	for (tries = 0; tries < 5; tries++) {
		/* Take a reading */
		if ((rv = dtp92_fcommand(p, "RM\r", buf, MAX_RD_SIZE, '>', 1, 3.5)) != DTP92_OK)
			return intep_code(pp, rv);

		if (sscanf(buf, " X%*c %lf\t Y%*c %lf\t Z%*c %lf ",
	           &val->aXYZ[0], &val->aXYZ[1], &val->aXYZ[2]) == 3) {

			val->aXYZ_v = 1;		/* These are absolute XYZ readings */
			rv = inst_ok;
			break;
		}
	}

#ifdef NEVER
		/* Set sample size back to 16 */
		if (((rv = p->command((inst *)p, "10SS\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
			return rv;
#endif

	return rv;
}

/* Determine if a calibration is needed. Returns inst_ok if not, */
/* inst_unsupported if it is unknown, or inst_needs_cal if needs calibration */
inst_code dtp92_needs_calibration(struct _inst *p) {
	/* This can't be known ? */
	return inst_unsupported;
}

/* Perform an instrument calibration */
inst_code dtp92_calibrate(
inst *pp,
int cal_no
) {
	dtp92 *p = (dtp92 *)pp;
	static char buf[MAX_RD_SIZE];
	int rv = 0;

	if (cal_no == 0) {

		/* Do offset calibration */
		if ((rv = dtp92_fcommand(p, "CO\r", buf, MAX_RD_SIZE, '>', 1, 3.5)) != DTP92_OK)
			return intep_code(pp, rv);

		return inst_ok;

	} else if (cal_no == 1) {

		/* Do ratio calibration */
		if ((rv = dtp92_fcommand(p, "CR\r", buf, MAX_RD_SIZE, '>', 1, 25.0)) != DTP92_OK)
			return intep_code(pp, rv);

		return inst_ok;
	}

	return inst_unsupported;
}

/* Error codes interpretation */
static char *
dtp92_interp_error(inst *pp, int ec) {
	dtp92 *p = (dtp92 *)pp;
	ec &= 0x7f;
	switch (ec) {
		case DTP92_USER_ABORT:
			return "User hit Escape";
		case DTP92_INTERNAL_ERROR:
			return "Internal software error";
		case DTP92_COMS_FAIL:
			return "Communications failure";
		case DTP92_UNKNOWN_MODEL:
			return "Not a DTP92 or DTP52";
		case DTP92_DATA_PARSE_ERROR:
			return "Data from DTP didn't parse as expected";
		case DTP92_OK:
			return "No error";
		case DTP92_BAD_COMMAND:
			return "Unrecognized command";

		case DTP92_PRM_RANGE:
			return "Command parameter out of range";
		case DTP92_MEMORY_OVERFLOW:
			return "Memory bounds error";
		case DTP92_INVALID_BAUD_RATE:
			return "Invalid baud rate";
		case DTP92_TIMEOUT:
			return "Receive timeout";
		case DTP92_SYNTAX_ERROR:
			return "Badly formed parameter";
		case DTP92_NO_DATA_AVAILABLE:
			return "No data available";
		case DTP92_MISSING_PARAMETER:
			return "Paramter is missing";
		case DTP92_CALIBRATION_DENIED:
			return "Invalid calibration enable code";
		case DTP92_NEEDS_OFFSET_CAL:
			return "Offset calibration checksum failed";
		case DTP92_NEEDS_RATIO_CAL:
			return "Ratio calibration checksum failed";
		case DTP92_NEEDS_LUMINANCE_CAL:
			return "Luminance calibration checksum failed";
		case DTP92_NEEDS_WHITE_POINT_CAL:
			return "White point calibration checksum failed";
		case DTP92_NEEDS_BLACK_POINT_CAL:
			return "Black point calibration checksum failed";
		case DTP92_INVALID_READING:
			return "Unable to take a reading";
		case DTP92_BAD_COMP_TABLE:
			return "Bad compensation table";
		case DTP92_TOO_MUCH_LIGHT:
			return "Too much light entering instrument";
		case DTP92_NOT_ENOUGH_LIGHT:
			return "Not enough light to complete operation";
		case DTP92_BAD_SERIAL_NUMBER:
			return "New serial number is invalid";
		case DTP92_NO_MODULATION:
			return "No refresh modulation detected";
		case DTP92_EEPROM_FAILURE:
			return "Eeprom write failure";
		case DTP92_FLASH_WRITE_FAILURE:
			return "Flash memory write failure";
		case DTP92_INST_INTERNAL_ERROR:
			return "Internal instrument error";
		default:
			return "Unknown error code";
	}
}


/* Convert a machine specific error code into an abstract dtp code */
static inst_code 
intep_code(inst *pp, int ec) {
	dtp92 *p = (dtp92 *)pp;

	ec &= 0x7f;
	switch (ec) {

		case DTP92_OK:
			return inst_ok;

		case DTP92_INTERNAL_ERROR:
			return inst_internal_error | ec;

		case DTP92_COMS_FAIL:
			return inst_coms_fail | ec;

		case DTP92_UNKNOWN_MODEL:
			return inst_unknown_model | ec;

		case DTP92_DATA_PARSE_ERROR:
			return inst_protocol_error | ec;

		case DTP92_USER_ABORT:
			return inst_user_abort | ec;

		case DTP92_NO_MODULATION:		/* ?? */
		case DTP92_TOO_MUCH_LIGHT:
		case DTP92_NOT_ENOUGH_LIGHT:
		case DTP92_INVALID_READING:
			return inst_misread | ec;

		case DTP92_NEEDS_OFFSET_CAL:
		case DTP92_NEEDS_RATIO_CAL:
		case DTP92_NEEDS_LUMINANCE_CAL:
		case DTP92_NEEDS_WHITE_POINT_CAL:
		case DTP92_NEEDS_BLACK_POINT_CAL:
			return inst_needs_cal | ec;
	}
	return inst_other_error | ec;
}

/* Return the last serial I/O error code */
static int
dtp92_last_sioerr(inst *pp) {
	dtp92 *p = (dtp92 *)pp;
	return p->sio->lerr;
}

/* Destroy ourselves */
static void
dtp92_del(inst *pp) {
	dtp92 *p = (dtp92 *)pp;
	if (p->sio != NULL)
		p->sio->del(p->sio);
	free(p);
}

/* Return the instrument capabilities */
inst_capability dtp92_capabilities(inst *pp) {
	dtp92 *p = (dtp92 *)pp;

	return 
	  inst_emis_spot
	| inst_emis_disp
	| inst_colorimeter
	| inst_calib
	  ;
}

/* Set device measurement mode */
inst_code dtp92_set_mode(inst *pp, inst_mode m)
{
	inst_mode mm;		/* Measurement mode */

	/* The measurement mode portion of the mode */
	mm = m & inst_mode_measuremet_mask;

	/* only display emission mode supported */
	if (mm != inst_mode_emis_disp
	 && mm != inst_mode_emis_spot) {
		return inst_unsupported;
	}

	/* Spectral mode is not supported */
	if (m & inst_mode_spectral)
		return inst_unsupported;

	return inst_ok;
}

/* 
 * set or reset an optional mode
 * We assume that the instrument has been initialised.
 */
static inst_code
dtp92_set_opt_mode(inst *pp, inst_opt_mode m)
{
	dtp92 *p = (dtp92 *)pp;
	return inst_unsupported;
}

/* Constructor */
extern dtp92 *new_dtp92()
{
	dtp92 *p;
	if ((p = (dtp92 *)calloc(sizeof(dtp92),1)) == NULL)
		error("dtp92: malloc failed!");

	p->sio = new_serio();

	p->init_coms    = dtp92_init_coms;
	p->init_inst    = dtp92_init_inst;
	p->capabilities = dtp92_capabilities;
	p->set_mode     = dtp92_set_mode;
	p->set_opt_mode     = dtp92_set_opt_mode;
	p->xy_sheet_release = dtp92_xy_sheet_release;
	p->xy_sheet_hold    = dtp92_xy_sheet_hold;
	p->xy_locate_start  = dtp92_xy_locate_start;
	p->xy_get_location  = dtp92_xy_get_location;
	p->xy_locate_end  = dtp92_xy_locate_end;
	p->xy_clear     = dtp92_xy_clear;
	p->read_xy      = dtp92_read_xy;
	p->read_strip   = dtp92_read_strip;
	p->read_sample  = dtp92_read_sample;
	p->needs_calibration = dtp92_needs_calibration;
	p->calibrate    = dtp92_calibrate;
	p->command      = dtp92_command;
	p->interp_error = dtp92_interp_error;
	p->inst_interp_error = NULL;				/* virtual constructor will do this */
	p->last_sioerr  = dtp92_last_sioerr;
	p->del          = dtp92_del;

	p->itype = instDTP92;

	return p;
}
