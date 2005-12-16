
/* 
 * Argyll Color Correction System
 *
 * Gretag Spectrolino and Spectroscan related
 * defines and declarations.
 *
 * Author: Graeme W. Gill
 * Date:   13/7/2005
 *
 * Copyright 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 *
 * Derived from DTP41.h
 *
 * This is an alternative driver to spm/gretag.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include "xspect.h"
#include "serio.h"
#include "ss.h"

#undef DEBUG

#include <stdarg.h>
void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);


/* Convert control chars to ^[A-Z] notation in a string */
static char *
fix(char *s) {
	static char buf[3 * SS_MAX_WR_SIZE];
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

/* Some tables to convert between emums and text descriptions */

/* Filter type */
char* filter_desc[] = {
	"Filter not defined",
	"No Filter (U)",
	"Polarizing Filter",
	"D65 Filter",
	"(Illegal)",
	"Custon Filter"
};

/* Track the number of measurements taken, so that a recalibration will */
/* be done at the right time. */
static void inc_calcount(ss *p) {
	p->calcount++;
	if (((p->mode & inst_mode_illum_mask) == inst_mode_reflection && p->calcount >= 50 ) ||
		((p->mode & inst_mode_illum_mask) == inst_mode_transmission
		 && (p->mode & inst_mode_sub_mask) != inst_mode_spot		/* ??? */
	     && p->calcount >= 10)) {
		p->calcount = 0;
		p->need_cal = 1;
	}
}

/* Establish communications with a Spectrolino/Spectroscan */
/* Use the baud rate given, and timeout in to secs */
/* Return DTP_COMS_FAIL on failure to establish communications */
static inst_code
ss_init_coms(inst *pp, int port, baud_rate br, double tout) {
	ss *p = (ss *)pp;
	/* We're a bit stuffed if the Specrolino/scan is set to 28800, since */
	/* this rate isn't universally supported by computer systems. */
	baud_rate brt[7] = { baud_9600,           baud_19200,         baud_57600,
	                     baud_2400,           baud_1200,          baud_600,
	                     baud_300 };
	ss_bt ssbrc[7]   = { ss_bt_9600,          ss_bt_19200,        ss_bt_57600,
	                     ss_bt_2400,          ss_bt_1200,         ss_bt_600,
	                     ss_bt_300 };
	ss_ctt sobrc[7]  = { ss_ctt_SetBaud9600, ss_ctt_SetBaud19200, ss_ctt_SetBaud57600,
	                     ss_ctt_SetBaud2400, ss_ctt_SetBaud1200,  ss_ctt_SetBaud600,
	                     ss_ctt_SetBaud300 };
	int bi;
	long etime;
	int i, rv;
	inst_code ev = inst_ok;

#ifdef DEBUG
	printf("Init comms called for ss\n");
#endif

	if (p->debug)
		p->sio->debug = p->debug;	/* Turn on debugging */

	/* Figure Spectrolino baud rate being asked for */
	for (bi = 0; bi < 7; bi++) {
		if (brt[bi] == br)
			break;
	}
	if (bi >= 7)
		bi = 0;

	/* The tick to give up on */
	etime = clock() + (long)(CLOCKS_PER_SEC * tout + 0.5);

	/* Until we establish comms or give up */
	while (clock() < etime) {

		/* Until we time out, find the correct baud rate */
		for (i = 0; clock() < etime;) {
#ifdef DEBUG
			printf("Trying baud rate %d\n",i);
#endif
			p->sio->set_port(p->sio, NULL, port, brt[i], parity_none, stop_1, length_8);

			/* Try a SpectroScan Output Status */
			ss_init_send(p);
			ss_add_ssreq(p, ss_OutputStatus);
			ss_command(p, SH_TMO);
			
			if (ss_sub_1(p) == ss_AnsPFX				/* Got comms */
			 || p->snerr == ss_et_UserAbort)			/* User aborted */
				break;

			/* Try a Spectrolino Parameter Request */
			ss_init_send(p);
			ss_add_soreq(p, ss_ParameterRequest);
			ss_command(p, SH_TMO);
			
			if (ss_sub_1(p) == ss_ParameterAnswer		/* Got comms */
			 || p->snerr == ss_et_UserAbort)			/* User aborted */
				break;

			if (++i >= 7)
				i = 0;
		}

		if (p->snerr == ss_et_UserAbort)
			return inst_user_abort;

		break;			/* Got coms */
	}

	if (clock() >= etime) {		/* We haven't established comms */
		return inst_coms_fail;
	}

#ifdef DEBUG
	printf("Got basic communications\n");
#endif

	p->itype = instUnknown;

	/* See if we have a Spectroscan or SpectroscanT */
	{
		char devn[19];

		ev = ss_do_OutputType(p, devn);

		if (ev == inst_ok) {
#ifdef DEBUG
			printf("Got device name '%s'\n",devn);
#endif
			if (strcmp(devn, "SpectroScan ") == 0) {
				p->itype = instSpectroScan;
			} else if (strcmp(devn, "SpectroScanT") == 0) {
				p->itype = instSpectroScanT;
			}
		}
	}
	/* Check if there's a Spectrolino */
	{
		char devn[19];
		int sn, sr, yp, mp, dp, hp, np;		/* Date & Time of production */
		int fswl, nosw, dsw;				/* Wavelengths sampled */
		ss_ttt tt;							/* Target Type */

		ev = so_do_TargetIdRequest(p, devn, &sn, &sr, &yp, &mp, &dp, &hp, &np,
		                           &tt, &fswl, &nosw, &dsw);

		if (ev == inst_ok) {

#ifdef DEBUG
			printf("Got device name '%s'\n",devn);
			printf("Got target type '%d'\n",tt);
			printf("Start wl %d, no wl %d, wl space %d\n",fswl, nosw, dsw);
#endif
			if (tt != ss_ttt_Spectrolino
			  || strcmp(devn, "Spectrolino") != 0)
				return inst_unknown_model;
	
			if (p->itype == instUnknown)	/* No SpectrScan */
			 	p->itype = instSpectrolino;
		}
	}
	
#ifdef DEBUG
	printf("Trying to change baud rate etc.\n");
#endif
	if (p->itype == instSpectrolino) {

		if ((ev = so_do_MeasControlDownload(p, ss_ctt_ProtokolWithoutXonXoff)) != inst_ok)
			return ev;

		if (bi != i) { 	/* Set the requested baud rate now. */
			/* Do baudrate change without checking results */
			so_do_MeasControlDownload(p, sobrc[bi]);
			p->sio->set_port(p->sio, NULL, port, brt[bi], parity_none, stop_1, length_8);
		}
	} else {	/* Spectroscan */

		ss_do_SetDeviceOnline(p);	/* Put the device online */

		/* Make sure other communication parameters are right */
		if ((ev = ss_do_ChangeHandshake(p, ss_hst_None)) != inst_ok)
			return ev;

		if (bi != i) { 	/* Set the requested baud rate now. */
			/* Do baudrate change without checking results */
			ss_do_ChangeBaudRate(p, ssbrc[bi]); 
			p->sio->set_port(p->sio, NULL, port, brt[bi], parity_none, stop_1, length_8);
		}

		/* Make sure the Spectrolino is talking to us. */
		if ((ev = ss_do_ScanSpectrolino(p)) != inst_ok) {
			return ev;
		}
	}

#ifdef DEBUG
	printf("Establish communications\n");
#endif
#ifdef EMSST
	printf("DEBUG: Emulating SpectroScanT with SpectroScan!\n");
#endif

	p->gotcoms = 1;
	return inst_ok;
}

/* Initialise the Spectrolino/SpectroScan. */
/* return non-zero on an error, with dtp error code */
static inst_code
ss_init_inst(inst *pp) {
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;

	if (p->gotcoms == 0)
		return inst_internal_error;		/* Must establish coms before calling init */

	/* Reset the instrument to a known state */
	if (p->itype != instSpectrolino) { 
		
		/* Initialise the device without resetting the baud rate */
		if (p->itype == instSpectroScanT) {
			if ((rv = ss_do_SetTableMode(p, ss_tmt_Reflectance)) != inst_ok)
				return rv; 
		}
		if ((rv = ss_do_SetDeviceOnline(p)) != inst_ok)
			return rv;
		if ((rv = ss_do_ResetKeyAcknowlge(p)) != inst_ok)
			return rv;
		if ((rv = ss_do_ReleasePaper(p)) != inst_ok)
			return rv;
		if ((rv = ss_do_InitMotorPosition(p)) != inst_ok)
			return rv;

		if (p->verb) {
			char buf[1000]; 
			char dn[19];		/* Device name */
			unsigned int sn;	/* Serial number */
			char pn[9];			/* Part number */
			int yp;				/* Year of production (e.g. 1996) */
			int mp;				/* Month of production (1-12) */
			int dp;				/* Day of production (1-31) */
			char sv[13];		/* Software version */

			if ((rv = ss_do_OutputType(p, dn)) != inst_ok)
				return rv;
			if ((rv = ss_do_OutputSerialNumber(p, &sn)) != inst_ok)
				return rv;
			if ((rv = ss_do_OutputArticleNumber(p, pn)) != inst_ok)
				return rv;
			if ((rv = ss_do_OutputProductionDate(p, &yp, &mp, &dp)) != inst_ok)
				return rv;
			if ((rv = ss_do_OutputSoftwareVersion(p, sv)) != inst_ok)
				return rv;

			sprintf(buf, "Device:     %s\n"
			             "Serial No:  %u\n"
			             "Part No:    %s\n"
			             "Prod Date:  %d/%d/%d\n"
			             "SW Version: %s\n", dn, sn, pn, dp, mp, yp, sv);
			p->user((inst *)p, inst_verb, 1, 0, buf);
		}
	}
	/* Do Spectrolino stuff */
	if ((rv = so_do_ResetStatusDownload(p, ss_smt_InitWithoutRemote)) != inst_ok)
		return rv;
	if ((rv = so_do_ExecWhiteRefToOrigDat(p)) != inst_ok)
		return rv;

	if (p->verb) {
		char buf[1000]; 
		char dn[19];		/* device name */
		ss_dnot dno;		/* device number */
		char pn[9];			/* part number */
		unsigned int sn;	/* serial number */
		char sv[13];		/* software release */
		int yp;				/* Year of production (e.g. 1996) */
		int mp;				/* Month of production (1-12) */
		int dp;				/* Day of production (1-31) */
		char devn[19];
		int sn2, sr, hp, np;	/* Date & Time of production */
		int fswl, nosw, dsw;	/* Wavelengths sampled */
		ss_ttt tt;				/* Target Type */

		if ((rv = so_do_DeviceDataRequest(p, dn, &dno, pn, &sn, sv)) != inst_ok)
			return rv;

		if ((rv = so_do_TargetIdRequest(p, devn, &sn2, &sr, &yp, &mp, &dp, &hp, &np,
	                           &tt, &fswl, &nosw, &dsw)) != inst_ok)
			return rv;

		sprintf(buf, "Device:     %s\n"
		             "Serial No:  %u\n"
		             "Part No:    %s\n"
		             "Prod Date:  %d/%d/%d\n"
		             "SW Version: %s\n", dn, sn, pn, dp, mp, yp, sv);
		p->user((inst *)p, inst_verb, 1, 0, buf);
	}

	/* Set the default colorimetric parameters */
	if ((rv = so_do_ParameterDownload(p, p->dstd, p->wbase, p->illum, p->obsv)) != inst_ok)
		return rv;

	/* Set the capabilities mask */
	p->cap = ( inst_ref_spot     |
	           inst_emis_spot    |
	           inst_emis_disp    |
	           inst_emis_illum   |
	           inst_colorimeter  |
	           inst_spectral     |
	           inst_calib 
	         );

	if (p->itype == instSpectrolino)
		p->cap |= inst_trans_spot;		/* Support this manually using a light table */

	if (p->itype == instSpectroScan
	 || p->itype == instSpectroScanT)	/* Only in reflective mode */
		p->cap |= (inst_ref_xy       |
	               inst_xy_holdrel   |
	               inst_xy_locate    |
	               inst_xy_position);

	if (p->itype == instSpectroScanT)
		p->cap |= inst_trans_spot;

	if (rv == inst_ok)
		p->inited = 1;

	return rv;
}

/* For an xy instrument, release the paper */
/* Return the inst error code */
static inst_code
ss_xy_sheet_release(
struct _inst *pp) {
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;

	rv = ss_do_ReleasePaper(p);
	return rv;
}

/* For an xy instrument, hold the paper */
/* Return the inst error code */
static inst_code 
ss_xy_sheet_hold(
struct _inst *pp) {
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;

	rv = ss_do_HoldPaper(p);
	return rv;
}

/* For an xy instrument, allow the user to locate a point */
/* Return the inst error code */
static inst_code
ss_xy_locate_start(
struct _inst *pp) {
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;

	rv = ss_do_SetDeviceOffline(p);
	return rv;
}

/* For an xy instrument, position the reeading point */
/* Return the inst error code */
static inst_code ss_xy_position(
struct _inst *pp,
int measure,			/* nz if position measure point, z if locate point */
double x, double y
) {
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;

	if ((rv = ss_do_MoveAbsolut(p, measure ? ss_rt_SensorRef : ss_rt_SightRef, x, y)) != inst_ok)
		return rv;

	return rv;
}

/* For an xy instrument, read back the location */
/* Return the inst error code */
static inst_code
ss_xy_get_location(
struct _inst *pp,
double *x, double *y) {
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;
	ss_rt rr;
	ss_zkt zk;

	if ((rv = ss_do_OutputActualPosition(p, ss_rt_SightRef, &rr, x, y, &zk)) != inst_ok)
		return rv;

	return rv;
}

/* For an xy instrument, ends allowing the user to locate a point */
/* Return the inst error code */
static inst_code
ss_xy_locate_end(
struct _inst *pp) {
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;

	rv = ss_do_SetDeviceOnline(p);
	return rv;
}

/* For an xy instrument, try and clear the table after an abort */
/* Return the inst error code */
static inst_code
ss_xy_clear(
struct _inst *pp) {
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;

	ss_do_SetDeviceOnline(p);	/* Put the device online */
	ss_do_MoveUp(p);			/* Raise the sensor */
	ss_do_ReleasePaper(p);		/* Release the paper */
	ss_do_MoveHome(p);			/* Move to the home position */

	return rv;
}

static inst_code ss_calibrate_imp(ss *p, int cal_no, int recal);

/* Read a sheet full of patches using xy mode */
/* Return the inst error code */
static inst_code
ss_read_xy(
inst *pp,
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
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;
	int pass, step, patch;
	inst_code err;
	unsigned char execerror;
	char buf[255];
	int tries, tc;			/* Total read tries */
	int fstep = 0;			/* NZ if step is fast & quiet direction */
	int cstep;				/* step value closest to calibration tile */

	if (p->itype != instSpectroScan
	 && p->itype != instSpectroScanT)
		return inst_unsupported;

	/* Move quickest in X direction to minimise noise, */
	/* and maximise speed */
	if (fabs(px) > fabs(ax)) {
		fstep = 1;
		if (px > 0.0)
			cstep = 0;
		else
			cstep = sip - 1;
	} else {
		if (ax > 0.0)
			cstep = 0;
		else
			cstep = pis - 1;
	}

	tries = sip * pis;

	/* Read all the patches */
	for (step = pass = tc = 0; tc < tries; tc++) {
		int astep, apass;	/* Actual step and pass to use */
		double ix, iy;

		/* Move in serpentine order */
		if (fstep) {
			astep = (pass & 1) ? sip - 1 - step : step;
			apass = pass;
		} else {
			astep = step;
			apass = (step & 1) ? pis - 1 - pass : pass;
		}

		patch = apass * sip + astep;

		if (patch < npatch) {	/* Over a valid patch */

			ix = ox + apass * ax + astep * px;
			iy = oy + apass * ay + astep * py;

			if ((step & 1) == 1) {	/* Offset for odd hex patches */
				ix += aax;
				iy += aay;
			}

			/* Do calibration if it is needed, */
			/* when closest to white tile. */
			if ((fstep ? astep == cstep : apass == cstep)
			    && p->need_cal && p->noutocalib == 0) {
				if ((rv = ss_calibrate_imp(p, 0, 1)) != inst_ok)
					return rv;
			}

			{
				ss_rvt refvalid;
				double col[3], spec[36];
				int i;

				vals[patch].XYZ_v = 0;
				vals[patch].aXYZ_v = 0;
				vals[patch].Lab_v = 0;
				vals[patch].spec_n = 0;
			
				/* move and measure gives us spectrum data anyway */
				if ((rv = ss_do_MoveAndMeasure(p, ix, iy, spec, &refvalid)) != inst_ok)
					return rv;

				vals[patch].spec_n = 36;
				vals[patch].spec_wl_short = 380;
				vals[patch].spec_wl_long = 730;
				for (i = 0; i < INSTR_MAX_BANDS; i++)
					if (i < vals[patch].spec_n)
						vals[patch].spec[i] = 100.0 * (double)spec[i];
					else
						vals[patch].spec[i] = -1.0;

				/* Get the XYZ */
				{
					ss_cst rct;
					ss_rvt rvf;
					ss_aft af;
					ss_wbt wb;
					ss_ilt it;
					ss_ot  ot;

					if ((rv = so_do_CParameterRequest(p, ss_cst_XYZ, &rct, col, &rvf,
					     &af, &wb, &it, &ot)) != inst_ok)
						return rv;
				}
				vals[patch].XYZ_v = 1;
				vals[patch].XYZ[0] = col[0];
				vals[patch].XYZ[1] = col[1];
				vals[patch].XYZ[2] = col[2];

				/* Track the need for a calibration */
				inc_calcount(p);
			}
		}

		/* Move on to the next patch */
		if (fstep) {
			if (++step >= sip) {
				step = 0;
				pass++;
			}
		} else {
			if (++pass >= pis) {
				pass = 0;
				step++;
			}
		}
	}

	return rv;
}

/* Read a set of strips */
/* Return the inst error code */
static inst_code
ss_read_strip(
inst *pp,
char *name,			/* Strip name (7 chars) */
int npatch,			/* Number of patches in the pass */
char *pname,		/* Pass name (3 chars) */
int sguide,			/* Guide number */
double pwid,		/* Patch length in mm (For DTP41) */
double gwid,		/* Gap length in mm (For DTP41) */
double twid,		/* Trailer length in mm (For DTP41T) */
ipatch *vals) {		/* Pointer to array of instrument patch values */
	ss *p = (ss *)pp;
	inst_code rv = inst_unsupported;
	return rv;
}


/* Observer weightings for Spectrolino spectrum */
/* 2 degree/10 degree, X, Y, Z */
double obsv[2][3][36] = {
	{
		{
			0.001446366705, 0.004669560650, 0.015038215600, 0.047557971000, 0.140987051500,
			0.276828785000, 0.342502310000, 0.334490535000, 0.287317520000, 0.195959900000,
			0.098192668500, 0.034533608000, 0.007030228050, 0.013030656650, 0.066802270500,
			0.166799475500, 0.291746645000, 0.434847990000, 0.594927790000, 0.761146250000,
			0.912975800000, 1.021235340000, 1.055716880000, 0.996420100000, 0.848984685000,
			0.644757190000, 0.449904640000, 0.287097725000, 0.167828905000, 0.090282337500,
			0.047660644500, 0.023668598500, 0.011722312400, 0.005977020750, 0.003003636650,
			0.001489205265
		},
		{
			0.000041511229, 0.000132336684, 0.000417218810, 0.001328120585, 0.004421233900,
			0.011872663850, 0.023193248000, 0.038522636000, 0.060559925000, 0.092293997500,
			0.140178090000, 0.211441485000, 0.328653775000, 0.505678185000, 0.704742665000,
			0.857418505000, 0.950276520000, 0.992144120000, 0.991708975000, 0.949110500000,
			0.867608170000, 0.756175865000, 0.630904045000, 0.503568860000, 0.381089945000,
			0.267281115000, 0.176511260000, 0.108674442500, 0.062212162000, 0.033089676500,
			0.017336597000, 0.008564816100, 0.004233960250, 0.002158411500, 0.001084668205,
			0.000537779300
		},
		{
			0.006818967200, 0.022079022300, 0.071358445500, 0.226964156000, 0.678939510000,
			1.354105715000, 1.721483290000, 1.766641960000, 1.649692095000, 1.285482130000,
			0.822264325000, 0.476534375000, 0.278198295000, 0.160407804000, 0.081826455500,
			0.042977953000, 0.021066016500, 0.009229350450, 0.004104619900, 0.002201453250,
			0.001622560100, 0.001153694700, 0.000797717335, 0.000383841335, 0.000179449332,
			0.000058813330, 0.000020030666, 0.000002926667, 0.000000000000, 0.000000000000,
			0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000,
			0.000000000000
		}
	},
	{
		{
			0.000272700780, 0.003297449000, 0.022707435000, 0.088927850000, 0.203483000000,
			0.312590900000, 0.377132950000, 0.367243000000, 0.300031250000, 0.194373900000,
			0.084108500000, 0.020172350000, 0.007532100000, 0.040903550000, 0.120323550000,
			0.237736250000, 0.377452300000, 0.532077950000, 0.705061800000, 0.875160750000,
			1.013217150000, 1.110392000000, 1.116849000000, 1.024365500000, 0.854241200000,
			0.646132200000, 0.436043700000, 0.271629400000, 0.155857600000, 0.083471760000,
			0.042258530000, 0.020691880000, 0.009953238000, 0.004740058000, 0.002262416000,
			0.001086550000
		},
		{
			0.000029496310, 0.000351341000, 0.002371945000, 0.009177750000, 0.021731100000,
			0.039180250000, 0.062144850000, 0.090069850000, 0.129007450000, 0.185882000000,
			0.256457250000, 0.343273100000, 0.462339900000, 0.607733100000, 0.757658700000,
			0.874167100000, 0.956879800000, 0.991016600000, 0.993553000000, 0.951687050000,
			0.869684200000, 0.774904950000, 0.657614200000, 0.527900150000, 0.399541600000,
			0.283781050000, 0.182358050000, 0.109413600000, 0.061649085000, 0.032699280000,
			0.016462845000, 0.008041705000, 0.003864443500, 0.001841166000, 0.000880087000,
			0.000423607700
		},
		{
			0.001204339800, 0.014693135000, 0.102659390000, 0.410999600000, 0.971019700000,
			1.545538000000, 1.936093000000, 1.976725000000, 1.734917000000, 1.303624000000,
			0.788122500000, 0.427392300000, 0.225584500000, 0.116998350000, 0.061871800000,
			0.031261400000, 0.014025850000, 0.004297150000, 0.000313350000, 0.000000000000,
			0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000,
			0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000,
			0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000, 0.000000000000,
			0.000000000000
		}
	}
};

/* Read a single sample */
/* Return the dtp error code */
static inst_code
ss_read_sample(
inst *pp,
char *name,			/* Strip name (7 chars) */
ipatch *val) {		/* Pointer to instrument patch value */
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;
	double col[3], spec[36];
	int i;

	if (!val)
		return inst_internal_error;

	val->XYZ_v = 0;
	val->aXYZ_v = 0;
	val->Lab_v = 0;
	val->spec_n = 0;

	/* Do calibration if it is needed */
	if (p->need_cal && p->noutocalib == 0) {
		if ((rv = ss_calibrate_imp(p, 0, 1)) != inst_ok)
			return rv;
	}

	/* For reflection spot mode on a SpectroScan, lower the head. */
	/* (A SpectroScanT in transmission will position automatically) */
	if (p->itype != instSpectrolino
	 && (p->mode & inst_mode_illum_mask) != inst_mode_transmission) {
		if ((rv = ss_do_MoveDown(p)) != inst_ok)
				return rv;
	}

	/* measure */
	if ((rv = so_do_ExecMeasurement(p)) != inst_ok)
		return rv;

	/* For reflection spot mode on a SpectroScan, raise the head. */
	/* (A SpectroScanT in transmission will position automatically) */
	if (p->itype != instSpectrolino
	 && (p->mode & inst_mode_illum_mask) != inst_mode_transmission) {
		if ((rv = ss_do_MoveUp(p)) != inst_ok)
				return rv;
	}

	/* Track the need for a calibration */
	inc_calcount(p);

	/* Get the XYZ: */

	/* Emulated spot transmission mode: */
	if ((p->mode & inst_mode_illum_mask) == inst_mode_transmission
	 && p->itype == instSpectrolino) {
		ss_st rst;		/* Return Spectrum Type (Reflectance/Density) */
		ss_rvt rvf;		/* Return Reference Valid Flag */
		ss_aft af;		/* Return filter being used (None/Pol/D65/UV/custom */
		ss_wbt wb;		/* Return white base (Paper/Absolute) */
		double norm;	/* Y normalisation factor */
		int j, tix = 0;	/* Default 2 degree */

		/* Get the spectrum */
		if ((rv = so_do_SpecParameterRequest(p, ss_st_LinearSpectrum,
		          &rst, spec, &rvf, &af, &wb)) != inst_ok)
			return rv;

		/* Divide by transmission white reference to get transmission level. */
		for (i = 0; i < 36; i++) {
			if (p->tref[i] >= 0.0001)
				spec[i] = spec[i]/p->tref[i];
			else
				spec[i] = 0.0;
		}

		if (p->mode & inst_mode_spectral) {
			val->spec_n = 36;
			val->spec_wl_short = 380;
			val->spec_wl_long = 730;
			for (i=0; i < INSTR_MAX_BANDS; i++)
				if (i < val->spec_n)
					val->spec[i] = 100.0 * (double)spec[i];
				else
					val->spec[i] = -1.0;
		}

		/* Convert to desired illuminant XYZ */
		val->XYZ[0] = 0.0;
		val->XYZ[1] = 0.0;
		val->XYZ[2] = 0.0;
		if (p->obsv == ss_ot_TenDeg)
			tix = 1;

		/* Compute normalisation factor */
		for (norm = 0.0, i = 0; i < 36; i++) {
			if (p->tref[i] >= 0.0001)
				norm += obsv[tix][1][i] * p->cill[i];
		}
		norm = 100.0/norm;
		
		/* Compute XYZ */
		for (i = 0; i < 36; i++) {
			if (p->tref[i] >= 0.0001) {
				for (j = 0; j < 3; j++)
					val->XYZ[j] += obsv[tix][j][i] * p->cill[i] * spec[i];
			}
		}
		for (j = 0; j < 3; j++)
			val->XYZ[j] *= norm;

		val->XYZ_v = 1;
	 

	/* Normal instrument values */
	} else {
		ss_cst rct;
		ss_rvt rvf;
		ss_aft af;
		ss_wbt wb;
		ss_ilt it;
		ss_ot  ot;

		if ((rv = so_do_CParameterRequest(p, ss_cst_XYZ, &rct, col, &rvf,
		     &af, &wb, &it, &ot)) != inst_ok)
			return rv;
		if ((p->mode & inst_mode_illum_mask) == inst_mode_emission) {
			val->aXYZ_v = 1;
			val->aXYZ[0] = col[0];
			val->aXYZ[1] = col[1];
			val->aXYZ[2] = col[2];
		} else {
			val->XYZ_v = 1;
			val->XYZ[0] = col[0];
			val->XYZ[1] = col[1];
			val->XYZ[2] = col[2];
		}

		/* spectrum data is returned only if requested */
		if (p->mode & inst_mode_spectral) {
			ss_st rst;		/* Return Spectrum Type (Reflectance/Density) */
			ss_rvt rvf;		/* Return Reference Valid Flag */
			ss_aft af;		/* Return filter being used (None/Pol/D65/UV/custom */
			ss_wbt wb;		/* Return white base (Paper/Absolute) */
	
			if ((rv = so_do_SpecParameterRequest(p, ss_st_LinearSpectrum,
			          &rst, spec, &rvf, &af, &wb)) != inst_ok)
				return rv;
	 
			val->spec_n = 36;
			val->spec_wl_short = 380;
			val->spec_wl_long = 730;
			for (i=0; i<INSTR_MAX_BANDS; i++)
				if (i < val->spec_n)
					val->spec[i] = 100.0 * (double)spec[i];
				else
					val->spec[i] = -1.0;
		}
	}
	return rv;
}


/* Determine if a calibration is needed. Returns inst_ok if not, */
/* inst_unsupported if it is unknown, or inst_needs_cal if needs calibration */
inst_code ss_needs_calibration(struct _inst *pp) {
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;

	if (p->need_cal)
		return inst_needs_cal;
	return rv;
}

/* Perform an instrument calibration (implementation). */
/* Switches to next mode. */
static inst_code ss_calibrate_imp(
ss *p,
int cal_no,				/* Type of calibration 0-n */
int recal				/* nz if recalibration rather than first calibration */
) {
	inst_code rv = inst_ok;
	ss_aft afilt;	/* Actual filter */
	double sp[36];	/* Spectral values of filter */
	ss_owrt owr;	/* Original white reference */

//printf("~1 Attempting calibration, recal = %d\n",recal);
	/* There are different procedures depending on the intended mode, */
	/* whether this is a Spectrolino or SpectrScan, and whether this is a */
	/* first calibration since changing a mode, or a recalibration. */
	p->mode = p->nextmode;
	if ((p->nextmode & inst_mode_illum_mask) == inst_mode_emission)
		p->filt = ss_aft_NoFilter;	/* Need no filter for emission */

	/* All first time calibrations do an initial reflective white calibration. */
	if (recal == 0
	 || (p->mode & inst_mode_illum_mask) == inst_mode_reflection) {
		char wrn[19];	/* White reference name */

		/* Set mode to reflection */
		if (p->itype == instSpectroScanT) {
			if ((rv = ss_do_SetTableMode(p, ss_tmt_Reflectance)) != inst_ok)
				return rv;
		} else {
			if ((rv = so_do_MeasControlDownload(p, ss_ctt_RemissionMeas)) != inst_ok)
				return rv;
		}

		/* Get the name of the expected white reference */
		if ((rv = so_do_WhiteReferenceRequest(p, p->filt, &afilt, sp, &owr, wrn)) != inst_ok)
			return rv;

		if (p->itype == instSpectrolino) {
			char buf1[100];
			sprintf(buf1, "Place Spectrolino on white reference %s and press:\n",wrn);
			if (p->user((inst *)p, inst_question, 3, 2, buf1, "Continue\n", "Abort\n") != 1)
				return inst_user_abort;
		} else if (recal == 0) { /* SpectroScan/T on first calibration */
			char buf1[100];
			sprintf(buf1, "Confirm white reference %s is in slot 1 and press:\n",wrn);
			if (p->user((inst *)p, inst_question, 3, 2, buf1, "Continue\n", "Abort\n") != 1)
				return inst_user_abort;
		}

		/* Do the white calibration */
		for (;;) {		/* Untill everything is OK */
			/* Set the desired colorimetric parameters + absolute white base */
			if ((rv = so_do_ParameterDownload(p, p->dstd, ss_wbt_Abs, p->illum, p->obsv))
			                                                                  != inst_ok)
				return rv;

			/* For SpectroScan, move to the white reference in slot 1 and lower */
			if (p->itype != instSpectrolino) {
				if ((rv = ss_do_MoveToWhiteRefPos(p, ss_wrpt_RefTile1)) != inst_ok)
					return rv;
				if ((rv = ss_do_MoveDown(p)) != inst_ok)
					return rv;
			}

			/* Calibrate */
			if ((rv = so_do_ExecRefMeasurement(p, ss_mmt_WhiteCalWithWarn))
			                          != (inst_notify | ss_et_WhiteMeasOK))
				return rv;
			rv = inst_ok;

			/* For SpectroScan, raise */
			if (p->itype != instSpectrolino) {
				if ((rv = ss_do_MoveUp(p)) != inst_ok)
					return rv;
			}

			/* Verify the filter. */
			{
				char buf[200];
				ss_dst ds;
				ss_wbt wb;
				ss_ilt it;
				ss_ot  ot;
				ss_aft af;

				if ((rv = so_do_ParameterRequest(p, &ds, &wb, &it, &ot, &af)) != inst_ok)
					return rv;
				if (af == p->filt)
					break;
				sprintf(buf, "Filter mismatch - replace filter with '%s' and press:\n",filter_desc[p->filt]);
				if (p->user((inst *)p, inst_question, 3, 2, buf, "Continue\n", "Abort\n") != 1)
					return inst_user_abort;
			}
		}
//printf("~1 reflection calibration and filter verify is complete\n");

		/* Emission or spot transmission mode */
		if ((p->mode & inst_mode_illum_mask) == inst_mode_emission
		 || ((p->mode & inst_mode_illum_mask) == inst_mode_transmission
		     && p->itype == instSpectrolino)) {
//printf("~1 emmission dark calibration:\n");
				/* Set emission mode */
				if ((rv = so_do_MeasControlDownload(p, ss_ctt_EmissionMeas)) != inst_ok)
					return rv;
				/* Do dark calibration (Assume we're still on white reference) */
				if ((rv = so_do_ExecRefMeasurement(p, ss_mmt_EmissionCal))
				                    != (inst_notify | ss_et_EmissionCalOK))
					return rv;
				rv = inst_ok;
//printf("~1 emmission dark calibration done\n");

		}
		/* SpectroScanT - Transmission mode */
		if ((p->mode & inst_mode_illum_mask) == inst_mode_transmission
		     && p->itype == instSpectroScanT) {

			/* Advise user to change aperture before switching to transmission mode. */
			if (p->user((inst *)p, inst_question, 3, 2,
			     "Ensure that desired transmission aperture is fitted and press:\n",
			     "Continue\n", "Abort\n") != 1)
				return inst_user_abort;

			if ((p->mode & inst_mode_illum_mask) == inst_mode_transmission) {
				if ((rv = ss_do_SetTableMode(p, ss_tmt_Transmission)) != inst_ok)
					return rv;

			}
		}
	}

//printf("~1 got out of first calibrate/filter loop\n");

	/* ??? If White Base Type is not Absolute, where is Paper type set, */
	/* ??? and how is it calibrated ????? */

	/* For non-reflective measurement, do the recalibration or 2nd part of calibration. */

	/* Transmission mode calibration: */
	if ((p->mode & inst_mode_illum_mask) == inst_mode_transmission) {

		/* Emulated spot transmission */
		if (p->itype == instSpectrolino) {
			if ((p->mode & inst_mode_illum_mask) == inst_mode_transmission) {
				ss_st rst;		/* Return Spectrum Type (Reflectance/Density) */
				ss_rvt rvf;		/* Return Reference Valid Flag */
				ss_aft af;		/* Return filter being used (None/Pol/D65/UV/custom */
				ss_wbt wb;		/* Return white base (Paper/Absolute) */
				ss_ilt it;		/* Return illuminant type */
				int i;

				if (p->user((inst *)p, inst_question, 3, 2,
				    "Place Spectrolino on transmission white reference and press:\n",
				    "Continue\n", "Abort\n") != 1)
					return inst_user_abort;
	
				/* Measure white reference spectrum */
				if ((rv = so_do_ExecMeasurement(p)) != inst_ok)
					return rv;
				if ((rv = so_do_SpecParameterRequest(p, ss_st_LinearSpectrum,
				          &rst, p->tref, &rvf, &af, &wb)) != inst_ok)
					return rv;

				/* See how good a source it is */
				for (i = 0; i < 36; i++)  {
					if (p->tref[i] < 0.0001)
						break;
				}
				if (i < 36) {
					p->user((inst *)p, inst_question, 1, 0,
					    "Warning: Transmission light source is low at some wavelengths!\n");
				}
				
				/* Get the instrument illuminant */
				if ((rv = so_do_IllumTabRequest(p, p->illum, &it, p->cill)) != inst_ok)
					return rv; 
			}

		/* SpectroScanT */
		} else {
			/* Presuming this is the right return code */
			if ((rv = so_do_ExecRefMeasurement(p, ss_mmt_WhiteCalWithWarn))
			                          != (inst_notify | ss_et_WhiteMeasOK))
				return rv;
			rv = inst_ok;
		}
	}
	
	p->need_cal = 0;

//printf("~1 calibration completed\n");
	return rv;
}

/* Perform an instrument calibration. */
/* Switches to next mode. */
inst_code ss_calibrate(
struct _inst *pp,
int cal_no				/* Type of calibration 0-n */
) {
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;

	/* Call the implementation code */
	rv =  ss_calibrate_imp(p, cal_no, 0);

	return rv;
}

/* Instrument specific error codes interpretation */
static char *
ss_interp_error(inst *pp, int ec) {
	ss *p = (ss *)pp;
	ec &= 0xff;

	switch (ec) {

		/* Device errors */
		case ss_et_NoError:
			return "No device error";
		case ss_et_MemoryFailure:
			return "Memory failure";
		case ss_et_PowerFailure:
			return "Power failure";
		case ss_et_LampFailure:
			return "Lamp failure";
		case ss_et_HardwareFailure:
			return "Hardware failure";
		case ss_et_FilterOutOfPos:
			return "Filter wheel out of position";
		case ss_et_SendTimeout:
			return "Data transmission timout";
		case ss_et_DriveError:
			return "Data drive defect";
		case ss_et_MeasDisabled:
			return "Measuring disabled";
		case ss_et_DensCalError:
			return "Incorrect input during densitometric calibration";
		case ss_et_EPROMFailure:
			return "Defective EPROM";
		case ss_et_RemOverFlow:
			return "Too much light or wrong white calibration";
		case ss_et_MemoryError:
			return "Checksum error in memory";
		case ss_et_FullMemory:
			return "Memory is full";
		case ss_et_WhiteMeasOK:
			return "White measurement is OK";
		case ss_et_NotReady:
			return "Instrument is not ready - please wait";
		case ss_et_WhiteMeasWarn:
			return "White measurement warning";
		case ss_et_ResetDone:
			return "Reset is done";
		case ss_et_EmissionCalOK:
			return "Emission calibration is OK";
		case ss_et_OnlyEmission:
			return "Only for emission (not reflection)";
		case ss_et_CheckSumWrong:
			return "Wrong checksum";
		case ss_et_NoValidMeas:
			return "No valid measurement (e.g. no white measurement)";
		case ss_et_BackupError:
			return "Error in backing up values";
		case ss_et_ProgramROMError:
			return "Errors in programming ROM";

		/* Incorporate remote error set codes thus: */
		case ss_et_NoValidDStd:
			return "No valid Density standard set";
		case ss_et_NoValidWhite:
			return "No valid White standard set";
		case ss_et_NoValidIllum:
			return "No valid Illumination set";
		case ss_et_NoValidObserver:
			return "No valid Observer set";
		case ss_et_NoValidMaxLambda:
			return "No valid maximum Lambda set";
		case ss_et_NoValidSpect:
			return "No valid spectrum";
		case ss_et_NoValidColSysOrIndex:
			return "No valid color system or index";
		case ss_et_NoValidChar:
			return "No valid character";
		case ss_et_DorlOutOfRange:
			return "Density is out of range";
		case ss_et_ReflectanceOutOfRange:
			return "Reflectance is out of range";
		case ss_et_Color1OutOfRange:
			return "Color 1 is out of range";
		case ss_et_Color2OutOfRange:
			return "Color 2 is out of range";
		case ss_et_Color3OutOfRange:
			return "Color 3 is out of range";
		case ss_et_NotAnSROrBoolean:
			return "Not an SR or Boolean";
		case ss_et_NoValidValOrRef:
			return "No valid value or reference";

		/* Translated scan error codes thus: */
		case ss_et_DeviceIsOffline:
			return "Device has been set offline";
		case ss_et_OutOfRange:
			return "A parameter of the command is out of range";
		case ss_et_ProgrammingError:
			return "Error writing to Flash-EPROM";
		case ss_et_NoUserAccess:
			return "No access to internal function";
		case ss_et_NoValidCommand:
			return "Unknown command sent";
		case ss_et_NoDeviceFound:
			return "Spectrolino can't be found";
		case ss_et_MeasurementError:
			return "Measurement error";
		case ss_et_NoTransmTable:
			return "SpectroScanT command when no tansmission table";
		case ss_et_NotInTransmMode:
			return "SpectroScanT transmission command in reflection mode";
		case ss_et_NotInReflectMode:
			return "SpectroScanT reflection command in transmission mode";

		/* Translated device communication errors */
		case ss_et_StopButNoStart:
			return "No start character received by instrument";
		case ss_et_IllegalCharInRec:
			return "Invalid character received by instrument";
		case ss_et_IncorrectRecLen:
			return "Record length received by instrument incorrect";
		case ss_et_IllegalRecType:
			return "Invalid message number receivec by instrument";
		case ss_et_NoTagField:
			return "No message number received by instrument";
		case ss_et_ConvError:
			return "Received data couldn't be converted by instrument";
		case ss_et_InvalidForEmission:
			return "Invalid message number for emission instrument";
		case ss_et_NoAccess:
			return "Failure in user identification by instrument";

		/* Our own communication errors here too. */
		case ss_et_SerialFail:
			return "Serial communications failure";
		case ss_et_UserAbort:
			return "User requested abort";
		case ss_et_SendBufferFull:
			return "Message send buffer is full";
		case ss_et_RecBufferEmpty:
			return "Message receive buffer is full";
		case ss_et_BadAnsFormat:
			return "Message received from instrument is badly formatted";
		case ss_et_BadHexEncoding:
			return "Message received from instrument has bad Hex encoding";
		case ss_et_RecBufferOverun:
			return "Message received from instrument would overflow recieve buffer";
		default:
			return "Unknown error code";
	}
}

/* Return the last serial I/O error code */
static int
ss_last_sioerr(inst *pp) {
	ss *p = (ss *)pp;
	return p->sio->lerr;
}

/* Return the instrument capabilities */
inst_capability ss_capabilities(inst *pp) {
	ss *p = (ss *)pp;

	return p->cap;
}

/* 
 * set measurement mode
 * We assume that the instrument has been initialised.
 * The measurement mode is not activated until it's actually needed.
 */
static inst_code
ss_set_mode(inst *pp, inst_mode m) {
	ss *p = (ss *)pp;
	inst_code rv = inst_ok;

	int cap = pp->capabilities(pp);
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
	} else if (mm == inst_mode_emis_spot) {
		if (!(cap & inst_emis_spot))
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

	p->nextmode = m;

	/* Now activate the next mode */
	if (((p->nextmode & inst_mode_illum_mask) == inst_mode_reflection
	  && (p->mode     & inst_mode_illum_mask) != inst_mode_reflection)
	 || ((p->nextmode & inst_mode_illum_mask) == inst_mode_emission
	  && (p->mode     & inst_mode_illum_mask) != inst_mode_emission)
	 || ((p->nextmode & inst_mode_illum_mask) == inst_mode_transmission
	  && (p->mode     & inst_mode_illum_mask) != inst_mode_transmission)) {

		/* Mode has changed */
		p->mode = p->nextmode;
		p->need_cal = 1;

		rv = ss_calibrate_imp(p, 0, 0);
	}
	return rv;
}

/* 
 * set or reset an optional mode
 * We assume that the instrument has been initialised.
 */
static inst_code
ss_set_opt_mode(inst *pp, inst_opt_mode m)
{
	ss *p = (ss *)pp;

	switch(m) {
		case inst_opt_noautocalib:
			p->noutocalib = 1;
			break;
	
		case inst_opt_autocalib:
			p->noutocalib = 0;
			break;

		default:
			return inst_unsupported;
	}

	return inst_ok;
}

/* Destroy ourselves */
static void
ss_del(inst *pp) {
	ss *p = (ss *)pp;

	if (p->inited) {
		ss_xy_clear(pp);
	}

	if (p->sio != NULL)
		p->sio->del(p->sio);
	free (p);
}

/* Constructor */
extern ss *new_gretag() {
	ss *p;
	if ((p = (ss *)calloc(sizeof(ss),1)) == NULL)
		error("ss: malloc failed!");

	p->sio = new_serio();

	/* Init public methods */
	p->init_coms    	= ss_init_coms;
	p->init_inst    	= ss_init_inst;
	p->capabilities 	= ss_capabilities;
	p->set_mode     	= ss_set_mode;
	p->set_opt_mode     = ss_set_opt_mode;
	p->xy_sheet_release = ss_xy_sheet_release;
	p->xy_sheet_hold    = ss_xy_sheet_hold;
	p->xy_locate_start  = ss_xy_locate_start;
	p->xy_get_location  = ss_xy_get_location;
	p->xy_locate_end	= ss_xy_locate_end;
	p->xy_position		= ss_xy_position;
	p->xy_clear     	= ss_xy_clear;
	p->read_xy      	= ss_read_xy;
	p->read_strip   	= ss_read_strip;
	p->read_sample  	= ss_read_sample;
	p->needs_calibration = ss_needs_calibration;
	p->calibrate    	= ss_calibrate;
	p->interp_error 	= ss_interp_error;
	p->inst_interp_error = NULL;				/* virtual constructor will do this */
	p->last_sioerr  	= ss_last_sioerr;
	p->del          	= ss_del;

	/* Init state */
	p->itype = instUnknown;						/* Unknown until initialised */
	p->cap = inst_unknown;						/* Unknown until initialised */
	p->mode = inst_mode_unknown;				/* Not in a known mode yet */
	p->nextmode = inst_mode_unknown;			/* Not in a known mode yet */

	/* Set default measurement configuration */
	p->filt = ss_aft_NoFilter;
	p->dstd = ss_dst_ANSIT;
	p->illum = ss_ilt_D50;
	p->obsv = ss_ot_TwoDeg;
	p->wbase = ss_wbt_Abs;
	p->phmode = ss_ctt_PhotometricAbsolute;
	p->phref = 1.0;

	/* Init serialisation stuff */
	p->sbuf  = &p->_sbuf[0];
	p->sbufe = &p->_sbuf[SS_MAX_WR_SIZE-2];		/* Allow one byte for nul */
	p->rbuf  = &p->_rbuf[0];
	p->rbufe = &p->_rbuf[0];					/* Initialy empty */
	p->snerr = ss_et_NoError;					/* COMs error */

#ifdef EMSST
	p->tmode = 0;			/* Reflection mode */
	p->sbr = ss_rt_SensorRef;
	p->sbx = 100.0;
	p->sby = 200.0;
#endif

	return p;
}
