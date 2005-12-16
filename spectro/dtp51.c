/* 
 * Argyll Color Correction System
 *
 * Xrite DTP51 related functions
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
#include "dtp51.h"

#undef DEBUG

#include <stdarg.h>
void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);

static inst_code interp_code(inst *pp, int ec);

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

/* Do a full featured command/response echange with the dtp51 */
/* Return the dtp error code. End on the specified number */
/* of specified characters, or expiry if the specified timeout */
/* Assume standard error code if tc = '>' and ntc = 1 */
/* Return a DTP51 error code */
static int
dtp51_fcommand(
	struct _dtp51 *p,
	char *in,			/* In string */
	char *out,			/* Out string buffer */
	int bsize,			/* Out buffer size */
	char tc,			/* Terminating character */
	int ntc,			/* Number of terminating characters */
	double to) {		/* Timout in seconts */
	int rv, se;

	if ((se = p->sio->write_read(p->sio, in, out, bsize, tc, ntc, to)) != 0) {
#ifdef DEBUG
	printf("dtp51 fcommand: serial i/o failure on write_read '%s'\n",fix(in));
#endif
		if (se & SIO_USER)
			return DTP51_USER_ABORT;
		return DTP51_COMS_FAIL;
	}
	rv = DTP51_OK;
	if (tc == '>' && ntc == 1) {
		rv = extract_ec(out);
		if (rv > 0) {
			rv &= 0x7f;
			if (rv != DTP51_OK)	/* Clear the error */
				p->sio->write(p->sio, "CE\r", 0.5);
		}
	}
#ifdef DEBUG
	printf("command '%s'",fix(in));
	printf(" returned '%s', value 0x%x\n",fix(out),rv);
#endif
	return rv;
}

/* Do a standard command/response echange with the dtp51 */
/* Return the dtp error code */
static inst_code
dtp51_command(inst *pp, char *in, char *out, int bsize) {
	dtp51 *p = (dtp51 *)pp;
	int rv = dtp51_fcommand(p, in, out, bsize, '>', 1, 1.5);
	return interp_code(pp, rv);
}

/* Establish communications with a DTP51 */
/* Use the baud rate given, and timeout in to secs */
/* Return DTP_COMS_FAIL on failure to establish communications */
static inst_code
dtp51_init_coms(inst *pp, int port, baud_rate br, double tout) {
	dtp51 *p = (dtp51 *)pp;
	static char buf[MAX_MES_SIZE];
	baud_rate brt[5] = { baud_9600, baud_19200, baud_4800, baud_2400, baud_1200 };
	char *brc[5]     = { "30BR\r",  "60BR\r",   "18BR\r",  "0CBR\r",  "06BR\r" };
	int bi;
	long etime;
	int i, rv;
	inst_code ev = inst_ok;

	if (p->debug)
		p->sio->debug = p->debug;	/* Turn on debugging */

	/* Figure DTP51 baud rate being asked for */
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
			if (((ev = interp_code(pp, dtp51_fcommand(p, "\r", buf, MAX_MES_SIZE, '>', 1, 0.5)))
			         & inst_mask) != inst_coms_fail)
				break;		/* We've got coms or user abort */
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
	if (((ev = p->command((inst *)p, "0004CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	if (i != bi) {
		/* Change the baud rate to the rate we've been told */
		if ((rv = p->sio->write_read(p->sio, brc[bi], buf, MAX_MES_SIZE, '>', 1, 0.5)) == 0) {
			if (extract_ec(buf) != DTP51_OK)
				return inst_coms_fail;
		}
	
		/* Configure our baud rate as well */
		p->sio->set_port(p->sio, NULL, port, brt[bi], parity_none, stop_1, length_8);

		/* Loose a character - just in case */
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

/* Build a strip definition as a set of passes */
static void
build_strip(
char *tp,			/* pointer to string buffer */
char *name,			/* Strip name (7 chars) */
int npatch,			/* Number of patches in the pass */
char *pname,		/* Pass name (3 chars) */
int sguide) {		/* Guide number */
	int i;

	/* Per strip */
	for (i = 0; i < 7 && *name != '\000'; i++)
		*tp++ = *name++;	/* Strip name padded to 7 characters */
	for (; i < 7; i++)
		*tp++ = ' ';
	*tp++ = '1'; *tp++ = '0'; *tp++ = '0'; *tp++ = '0';	/* N factor */
	*tp++ = 'F';				/* Forward order */
	*tp++ = '1';				/* No min/max transmitted */
	*tp++ = '0';				/* No absolute data, no extra steps */
	*tp++ = '0'; *tp++ = '0'; *tp++ = '0'; *tp++ = '0'; *tp++ = '0';	/* Reserved */

	/* Per pass */
	for (i = 0; i < 3 && *pname != '\000'; i++)
		*tp++ = *pname++;		/* Pass name padded to 3 characters */
	for (; i < 3; i++)
		*tp++ = ' ';
	/* *tp++ = '4'; */			/* Lab data */
	*tp++ = '5';				/* XYZ data */
	*tp++ = '8';				/* Auto color */
	*tp++ = '0' + npatch/10;	/* Number of patches MS */
	*tp++ = '0' + npatch%10;	/* Number of patches LS */
	*tp++ = '0' + sguide/10;	/* Guide location MS */
	*tp++ = '0' + sguide%10;	/* Guide location LS */
	*tp++ = '0';				/* (Data output type) */
	*tp++ = '0';				/* Extra steps */
	*tp++ = '0';				/* Reserved */

	*tp++ = '\r';				/* The CR */
	*tp++ = '\000';				/* The end */
}
		
/* Initialise the DTP51. The spectral flag is ignored. */
/* return non-zero on an error, with dtp error code */
static inst_code
dtp51_init_inst(inst *pp) {
	dtp51 *p = (dtp51 *)pp;
	static char tbuf[100], buf[MAX_MES_SIZE];
	int i;
	inst_code ev = inst_ok;

	if (p->gotcoms == 0)
		return inst_internal_error;		/* Must establish coms before calling init */

	/* Reset it */
	if (((ev = p->command((inst *)p, "0PR\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;
	sleep(2);	/* Let it recover from reset */

	/* Turn echoing of characters off */
	if (((ev = p->command((inst *)p, "EC\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Get the model and version number */
	if (((ev = p->command((inst *)p, "SV\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Check that it is a DTP51 or 52 */
	if (   strlen(buf) < 12
	    || strncmp(buf,"X-Rite DTP5",11) != 0
	    || (buf[11] != '1' && buf[11] != '2'))
		return inst_unknown_model;

	/* Set the A/D Conversion rate to normal */
	if (((ev = p->command((inst *)p, "00AD\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Disable Bar code reading to save time */
	if (((ev = p->command((inst *)p, "0BC\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set the Calibration Check tolerance to default of 0.15D */
	if (((ev = p->command((inst *)p, "0FCC\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set Automatic Transmit off */
	if (((ev = p->command((inst *)p, "0005CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set Read Status on */
	if (((ev = p->command((inst *)p, "1RS\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set decimal pont on */
	if (((ev = p->command((inst *)p, "0106CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set color data separator to TAB */
	if (((ev = p->command((inst *)p, "0207CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set default strip format to off */
	if (((ev = p->command((inst *)p, "0009CF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set extra digit resolution */
	if (((ev = p->command((inst *)p, "010ACF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Set data after pass to off */
	if (((ev = p->command((inst *)p, "000BCF\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

#ifdef NEVER	/* Doesn't seem to work on DTP51 */
	/* Set the patch detection window to default HRW = 3, LRW 0.122% */
	if (((ev = p->command((inst *)p, "2CW\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;
#endif

	/* Enable the LCD display */
	if (((ev = p->command((inst *)p, "EL\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;

#ifdef NEVER
	/* Disable the read microswitch */
	if (((ev = p->command((inst *)p, "0SM\r", buf, MAX_MES_SIZE)) & inst_mask) != inst_ok)
		return ev;
#endif

	/* Set a strip length of 1, to ensure parsing is invalidated */
	build_strip(tbuf, "       ", 1, "   ", 30);
	
	if (((ev = dtp51_fcommand(p, "0105DS\r", buf, MAX_MES_SIZE, '*', 1, 0.5)) & inst_mask) != inst_ok)
		return ev;

	/* Expect '*' as response */
	if (buf[0] != '*' || buf[1] != '\000')
		return inst_coms_fail;

	/* Send strip definition */
	if (((ev = dtp51_fcommand(p, tbuf, buf, MAX_MES_SIZE, '>', 1, 4.0)) & inst_mask) != inst_ok)
		return ev;

	if (ev == inst_ok)
		p->inited = 1;

	return ev;
}

/* For an xy instrument, release the paper */
/* Return the inst error code */
static inst_code
dtp51_xy_sheet_release(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, hold the paper */
/* Return the inst error code */
static inst_code 
dtp51_xy_sheet_hold(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, allow the user to locate a point */
/* Return the inst error code */
static inst_code
dtp51_xy_locate_start(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, read back the location */
/* Return the inst error code */
static inst_code
dtp51_xy_get_location(
struct _inst *p,
double *x, double *y) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, ends allowing the user to locate a point */
/* Return the inst error code */
static inst_code
dtp51_xy_locate_end(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* For an xy instrument, try and clear the table after an abort */
/* Return the inst error code */
static inst_code
dtp51_xy_clear(
struct _inst *p) {
	/* This is not supported by this device */
	return inst_unsupported;
}

/* Read a sheet full of patches using xy mode */
/* Return the inst error code */
static inst_code
dtp51_read_xy(
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
dtp51_read_strip(
inst *pp,
char *name,			/* Strip name (up to first 7 chars used) */
int npatch,			/* Number of patches in the pass */
char *pname,		/* Pass name (up to first 3 chars used) */
int sguide,			/* Guide number */
double pwid,		/* Patch length in mm (For DTP41) */
double gwid,		/* Gap length in mm (For DTP41) */
double twid,		/* Trailer length in mm (For DTP41T) */
ipatch *vals) {		/* Pointer to array of instrument patch values */
	dtp51 *p = (dtp51 *)pp;
	char tbuf[200], *tp;
	static char buf[MAX_RD_SIZE];
	int i, rv = 0;
	inst_code ev = inst_ok;

	build_strip(tbuf, name, npatch, pname, sguide);
	
	if ((rv = dtp51_fcommand(p, "0105DS\r", buf, MAX_RD_SIZE, '*', 1, 0.5)) != DTP51_OK)
		return interp_code(pp, rv); 

	/* Expect '*' as response */
	if (buf[0] != '*' || buf[1] != '\000')
		return inst_coms_fail;

	/* Send strip definition */
	if ((rv = dtp51_fcommand(p, tbuf, buf, MAX_RD_SIZE, '>', 1, 4.0)) != DTP51_OK)
		return interp_code(pp, rv); 

	/* Set Read Status on */
	if (((ev = p->command((inst *)p, "1RS\r", buf, MAX_RD_SIZE)) & inst_mask) != inst_ok)
		return ev;

#ifdef NEVER
	/* Enable the read microswitch */
	if (((ev = p->command((inst *)p, "1SM\r", buf, MAX_RD_SIZE)) & inst_mask) != inst_ok)
		return ev;
#endif

	/* Do a strip read */
	if (((ev = p->command((inst *)p, "5SW\r", buf, MAX_RD_SIZE)) & inst_mask) != inst_ok)
		return ev;

	/* Wait for the Read status, or a user abort - allow 5 munutes. */
	rv = dtp51_fcommand(p, "", buf, MAX_RD_SIZE, '>', 1, 5 * 60.0);

#ifdef NEVER
	/* disable the read switch */
	if ((interp_code(pp, rv) & inst_mask) == inst_misread)
		p->command((inst *)p, "0SM\r", buf, MAX_RD_SIZE);		/* Disable the microswitch */
#endif

	/* Soft reset the unit */
	p->command((inst *)p, "0PR\r", buf, MAX_RD_SIZE); /* Soft Reset it */

	if (rv != DTP51_USER_ABORT)
		sleep(2);	/* Let it recover from reset */

	if (rv != DTP51_OK) {
		return interp_code(pp, rv); 	/* misread */
	}


	/* Gather the results */
	if ((rv = dtp51_fcommand(p, "TS\r", buf, MAX_RD_SIZE, '>', 1, 0.5 + npatch * 0.1)) != DTP51_OK)
		return interp_code(pp, rv);

	/* Parse the buffer */
	/* Replace '\r' with '\000' */
	for (tp = buf; *tp != '\000'; tp++) {
		if (*tp == '\r')
			*tp = '\000';
	}
	for (tp = buf, i = 0; i < npatch; i++) {
		if (*tp == '\000')
			return inst_protocol_error;
#ifdef NEVER	/* Lab */
		if (sscanf(tp, " L %lf a %lf b %lf ",
		           &vals[i].Lab[0], &vals[i].Lab[1], &vals[i].Lab[2]) != 3) {
			if (sscanf(tp, " l %lf a %lf b %lf ",
			           &vals[i].Lab[0], &vals[i].Lab[1], &vals[i].Lab[2]) != 3) {
				return inst_protocol_error;
			}
		}
		vals[i].Lab_v = 1;
#else /* XYZ */
		if (sscanf(tp, " X %lf Y %lf Z %lf ",
		           &vals[i].XYZ[0], &vals[i].XYZ[1], &vals[i].XYZ[2]) != 3) {
			if (sscanf(tp, " x %lf y %lf z %lf ",
			           &vals[i].XYZ[0], &vals[i].XYZ[1], &vals[i].XYZ[2]) != 3) {
				return inst_protocol_error;
			}
		}
		vals[i].XYZ_v = 1;
#endif
		tp += strlen(tp) + 1;
	}

#ifdef NEVER
	/* Disable the read microswitch */
	if (((ev = p->command((inst *)p, "0SM\r", buf, MAX_RD_SIZE)) & inst_mask) != inst_ok)
		return ev;
#endif

	return inst_ok;
}

/* Read a single sample */
/* Return the dtp error code */
static inst_code
dtp51_read_sample(
inst *pp,
char *name,			/* Strip name (7 chars) */
ipatch *val) {		/* Pointer to instrument patch value */

	/* This is not supported by DTP51 */
	return inst_unsupported;
}

/* Determine if a calibration is needed. Returns inst_ok if not, */
/* inst_unsupported if it is unknown, or inst_needs_cal if needs calibration */
inst_code dtp51_needs_calibration(struct _inst *p) {
	/* This can't be known ? */
	return inst_unsupported;
}

/* Perform an instrument calibration */
inst_code dtp51_calibrate(
struct _inst *p,
int cal_no
) {
	/* This is not currently supported */
	return inst_unsupported;
}

/* Error codes interpretation */
static char *
dtp51_interp_error(inst *pp, int ec) {
	dtp51 *p = (dtp51 *)pp;
	ec &= 0x7f;
	switch (ec) {
		case DTP51_USER_ABORT:
			return "User hit Escape";
		case DTP51_INTERNAL_ERROR:
			return "Internal software error";
		case DTP51_COMS_FAIL:
			return "Communications failure";
		case DTP51_UNKNOWN_MODEL:
			return "Not a DTP51 or DTP52";
		case DTP51_DATA_PARSE_ERROR:
			return "Data from DTP didn't parse as expected";
		case DTP51_OK:
			return "No error";
		case DTP51_BAD_COMMAND:
			return "Unrecognized command";
		case DTP51_PRM_RANGE:
			return "Command parameter out of range";
		case DTP51_DISPLAY_OVERFLOW:
			return "Display overflow";
		case DTP51_MEMORY_OVERFLOW:
			return "Memory bounds error";
		case DTP51_INVALID_BAUD_RATE:
			return "Invalid baud rate";
		case DTP51_TIMEOUT:
			return "Receive timeout";
		case DTP51_INVALID_PASS:
			return "Invalid pass";
		case DTP51_INVALID_STEP:
			return "Invalid step";
		case DTP51_NO_DATA_AVAILABLE:
			return "No data availble";
		case DTP51_LAMP_MARGINAL:
			return "Lamp marginal";
		case DTP51_LAMP_FAILURE:
			return "Lamp failure";
		case DTP51_STRIP_RESTRAINED:
			return "Strip was restrained";
		case DTP51_BAD_CAL_STRIP:
			return "Bad calibration strip";
		case DTP51_MOTOR_ERROR:
			return "Motor error";
		case DTP51_BAD_BARCODE:
			return "Bad barcode on cal strip";
		case DTP51_INVALID_READING:
			return "Invalid strip reading";
		case DTP51_WRONG_COLOR:
			return "Wrong color strip";
		case DTP51_BATTERY_TOO_LOW:
			return "Battery too low";
		case DTP51_NEEDS_CALIBRATION:
			return "Needs calibration";
		case DTP51_COMP_TABLE_MISMATCH:
			return "Compensation table mismatch";
		case DTP51_BAD_COMP_TABLE:
			return "Bad compensation table";
		case DTP51_NO_VALID_DATA:
			return "No valid data found";
		case DTP51_BAD_PATCH:
			return "Bad patch in strip";
		case DTP51_BAD_STRING_LENGTH:
			return "Bad strip def. string length";
		case DTP51_BAD_CHARACTER:
			return "Bad chareter";
		case DTP51_BAD_MEAS_TYPE:
			return "Bad measure type field";
		case DTP51_BAD_COLOR:
			return "Bad color field";
		case DTP51_BAD_STEPS:
			return "Bad step field";
		case DTP51_BAD_STOP_LOCATION:
			return "Bad guide stop field";
		case DTP51_BAD_OUTPUT_TYPE:
			return "Bad output type field";
		case DTP51_MEMORY_ERROR:
			return "Memory error (need AC supply)";
		case DTP51_BAD_N_FACTOR:
			return "Bad N factore";
		case DTP51_STRIP_DOESNT_EXIST:
			return "Strip doesn't exist";
		case DTP51_BAD_MIN_MAX_VALUE:
			return "Bad min/max field value";
		case DTP51_BAD_SERIAL_NUMBER:
			return "Bad serial number";
		default:
			return "Unknown error code";
	}
}


/* Convert a machine specific error code into an abstract dtp code */
static inst_code 
interp_code(inst *pp, int ec) {
	dtp51 *p = (dtp51 *)pp;

	ec &= 0x7f;
	switch (ec) {

		case DTP51_OK:
			return inst_ok;

		case DTP51_INTERNAL_ERROR:
			return inst_internal_error | ec;

		case DTP51_COMS_FAIL:
			return inst_coms_fail | ec;

		case DTP51_UNKNOWN_MODEL:
			return inst_unknown_model | ec;

		case DTP51_DATA_PARSE_ERROR:
			return inst_protocol_error | ec;

		case DTP51_USER_ABORT:
			return inst_user_abort | ec;

		case DTP51_NO_DATA_AVAILABLE:
		case DTP51_LAMP_MARGINAL:
		case DTP51_LAMP_FAILURE:
		case DTP51_STRIP_RESTRAINED:
		case DTP51_MOTOR_ERROR:
		case DTP51_INVALID_READING:
		case DTP51_WRONG_COLOR:
		case DTP51_BATTERY_TOO_LOW:
		case DTP51_COMP_TABLE_MISMATCH:
		case DTP51_BAD_COMP_TABLE:
		case DTP51_NO_VALID_DATA:
		case DTP51_BAD_PATCH:
			return inst_misread | ec;

		case DTP51_NEEDS_CALIBRATION:
			return inst_needs_cal | ec;
	}
	return inst_other_error | ec;
}

/* Return the last serial I/O error code */
static int
dtp51_last_sioerr(inst *pp) {
	dtp51 *p = (dtp51 *)pp;
	return p->sio->lerr;
}

/* Destroy ourselves */
static void
dtp51_del(inst *pp) {
	dtp51 *p = (dtp51 *)pp;
	if (p->sio != NULL)
		p->sio->del(p->sio);
	free(p);
}

/* Return the instrument capabilities */
inst_capability dtp51_capabilities(inst *pp) {
	dtp51 *p = (dtp51 *)pp;

	return 
	  inst_ref_strip
	| inst_colorimeter
	  ;
}

/* Set device measurement mode */
inst_code dtp51_set_mode(inst *pp, inst_mode m)
{
	inst_mode mm;		/* Measurement mode */

	/* The measurement mode portion of the mode */
	mm = m & inst_mode_measuremet_mask;

	/* only reflection strip measurement mode suported */
	if (mm != inst_mode_ref_strip) {
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
dtp51_set_opt_mode(inst *pp, inst_opt_mode m)
{
	dtp51 *p = (dtp51 *)pp;
	return inst_unsupported;
}

/* Constructor */
extern dtp51 *new_dtp51()
{
	dtp51 *p;
	if ((p = (dtp51 *)calloc(sizeof(dtp51),1)) == NULL)
		error("dtp51: malloc failed!");

	p->sio = new_serio();

	p->init_coms    	= dtp51_init_coms;
	p->init_inst    	= dtp51_init_inst;
	p->capabilities 	= dtp51_capabilities;
	p->set_mode     	= dtp51_set_mode;
	p->set_opt_mode     = dtp51_set_opt_mode;
	p->xy_sheet_release = dtp51_xy_sheet_release;
	p->xy_sheet_hold    = dtp51_xy_sheet_hold;
	p->xy_locate_start  = dtp51_xy_locate_start;
	p->xy_get_location  = dtp51_xy_get_location;
	p->xy_locate_end	= dtp51_xy_locate_end;
	p->xy_clear     	= dtp51_xy_clear;
	p->read_xy      	= dtp51_read_xy;
	p->read_strip   	= dtp51_read_strip;
	p->read_sample  	= dtp51_read_sample;
	p->needs_calibration = dtp51_needs_calibration;
	p->calibrate    	= dtp51_calibrate;
	p->command      	= dtp51_command;
	p->interp_error 	= dtp51_interp_error;
	p->inst_interp_error = NULL;				/* virtual constructor will do this */
	p->last_sioerr  	= dtp51_last_sioerr;
	p->del          	= dtp51_del;

	p->itype = instDTP51;

	return p;
}
