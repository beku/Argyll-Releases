/* Gretag Spectrolino/Spectroscan command parser
 *
 * Author: Neil Okamoto
 * Date: 1/9/2001
 *
 * Copyright 2001, DreamWorks LLC
 * All Rights Reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include "xspect.h"
#include "gretag.h"
#include "serio.h"

#define SPM_BUFFER_SZ 1000

/* forward decls */
static inst_code scanerror_code(unsigned char);
static inst_code remoteerrorset_code(unsigned short);
static inst_code execerror_code(unsigned char);
static inst_code activate_mode(gretag *p);
void error(char *fmt, ...);

/* Should allow for XON/XOFF handshaking, because SS by default has this mode ? */
/* <Xon> = ASCII(17), <Xoff> = ASCII(19) */
/* default Gretag speed is 9600 */

/**********************************************************************
 * Public Methods
 **********************************************************************/

/* gt_init_coms()
 * setup serial communications with the device
 */
static inst_code
gt_init_coms(inst *pp, int port, baud_rate baud, double tout)
{
	gretag *p = (gretag*)pp;
	static baud_rate br[] = {baud_9600, baud_19200, baud_57600,
							 baud_4800, baud_2400,
							 baud_1200, baud_600, baud_300};
	static int bri[] = {9600, 19200, 57600,
						4800, 2400,
						1200, 600, 300};
	static int nbr = sizeof(br)/sizeof(baud_rate);
	inst_code rv = inst_unknown_model;
	unsigned char n1,n2,n3,n4,n5;
	char buf[255];
	int i;

	if (p->debug)
		p->sio->debug = p->debug;	/* Turn on debugging */

	/* set to user-specified baud rate */
	i = 0;
	if (baud != baud_nc) {			/* Set user specified baud rate */
		for (i=0; i<nbr; i++) {
			if (baud == br[i])
				break;
		}
	}

	for (i=0; i<nbr; i++) {
		/* Try SpectrScan status check */
		p->sio->set_port(p->sio, NULL, port, br[i], parity_none, stop_1, length_8);
		rv = p->spm_command(p, SPM_PROTO_OUTPUTSTATUS, &n1);
		if (rv == inst_ok || rv == inst_unexpected_reply) {
			break;
		}
		if (baud != baud_nc) {	/* Was user request, switch to trying all */
			baud = baud_nc;
			i = 0;
		}
	}

	if (rv == inst_unexpected_reply) {
		/* See if it is a Spectrolino */
		rv = p->spm_command(p, SPM_PROTO_PARAMETERREQUEST, &n1, &n2, &n3, &n4, &n5);
		if (rv == inst_ok)
			p->itype = instSpectrolino;

	} else {
		char name18[19];
		p->itype = instSpectroScan;

		/* Check if it is a SpectroScanT */
		rv = p->spm_command(p, SPM_PROTO_OUTPUTTYPE, name18);
		name18[18] = '\000';
		if (rv == inst_ok && strncmp("SpectroScanT", name18, 12) == 0) {
			p->itype = instSpectroScanT;
		}
	}

	if (rv != inst_ok) {
		sprintf(buf, "gt_init_coms: failed to establish coms on port %d\n", port);
		p->logwarn(buf);
	}
	else {
		sprintf(buf, "gt_init_coms: established coms on port %d baud %d, device %s\n",
		port, bri[i], inst_name(p->itype));
		p->logwarn(buf);
	}

	return rv;
}


/* gt_init_gretag()
 * reset/initialize the device
 */
static inst_code
gt_init_inst(inst *pp)
{
	gretag *p = (gretag*)pp;
	unsigned char devicename[19 + 20];
	unsigned char devicenumber;
	unsigned char articlenumber[9 + 20];
	unsigned short serialnumber[2];
	unsigned char softwarerelease[13 + 20];
	unsigned char reserve[17 + 20];
	unsigned char scanerror;
	char buf[255];
	int rv;

	p->initialized = false;

	/* initialize spectroscan table, if present */

	if (p->itype == instSpectroScan || p->itype == instSpectroScanT) {
		if ((rv = p->spm_command(p, SPM_PROTO_SCANSPECTROLINO, &scanerror)) != inst_ok) {
			p->logerr("gt_init_gretag: couldn't init");
			return rv;
		}
		if (scanerror != SPM_SCANERROR_NOERROR) {
			sprintf(buf, "gt_init_gretag: scanerror %x", scanerror);
			p->logerr(buf);
			return scanerror_code(scanerror);
		}
		if ((rv = p->spm_command(p, SPM_PROTO_OUTPUTTYPE, devicename)) != inst_ok) {
			p->logerr("gt_init_gretag: couldn't get device name");
			return rv;
		}
		devicename[18] = '\000';
		if ((rv = p->spm_command(p, SPM_PROTO_OUTPUTSERIALNUMBER, serialnumber)) != inst_ok) {
			p->logerr("gt_init_gretag: couldn't get serial number");
			return rv;
		}
		if ((rv = p->spm_command(p, SPM_PROTO_OUTPUTARTICLENUMBER, articlenumber)) != inst_ok) {
			p->logerr("gt_init_gretag: couldn't get article number");
			return rv;
		}
		articlenumber[8] = '\000';
		if ((rv = p->spm_command(p, SPM_PROTO_OUTPUTSOFTWAREVERSION, softwarerelease)) != inst_ok) {
			p->logerr("gt_init_gretag: couldn't get firmware version");
			return rv;
		}
		softwarerelease[12] = '\000';

		sprintf(buf, "Device: %s\nPart#: %s\nSerial#: %d\nFirmware: v%s\n",
				devicename, articlenumber,
				serialnumber[0]+serialnumber[1]*65536,
				softwarerelease);
		p->logwarn(buf);
	}

	/* initialize spectrolino */
	if ((rv = p->spm_command(p, SPM_PROTO_DEVICEDATAREQUEST,
						  devicename, &devicenumber,
						  articlenumber, serialnumber,
						  softwarerelease, reserve)) != inst_ok) {
		p->logerr("gt_init_gretag: couldn't init");
		return rv;
	}
	sprintf(buf, "Device: %s\nPart#: %s\nSerial#: %d\nFirmware: v%s\n",
			devicename, articlenumber,
			serialnumber[0]+serialnumber[1]*65536,
			softwarerelease);
	p->logwarn(buf);

	/* Set to standard XYZ readings */
 	/* Absolute white */
	/* 2 degree observer */
	/* D50 illuminant */
	/* Status T density */
	/* No filter fitted */

	if ((rv = p->set_whitebase(p, GT_WHITEBASE_ABSOLUTE)) != inst_ok)
		return rv;
	if((rv = p->set_observer(p, GT_OBSERVER_2)) != inst_ok)
		return rv;
	if((rv = p->set_illum(p, GT_ILLUM_D50)) != inst_ok)
		return rv;
	if((rv = p->set_dstd(p, GT_DSTD_ANSIT)) != inst_ok)
		return rv;
	if((rv = p->set_filter(p, GT_FILTER_NONE)) != inst_ok)
		return rv;

	p->initialized = true;

	/* We are configured in this mode now ?? */
	p->mode = inst_mode_unknown;

	return inst_ok;
}


/* gt_capabilities()
 * return capabilities flags
 */
static inst_capability
gt_capabilities(inst *pp)
{
	gretag *p = (gretag*)pp;
	int cap = ( inst_ref_spot     |
				inst_emis_disp    |
				inst_emis_illum   |		// ?????
				inst_colorimeter  |
				inst_spectral     |
				inst_calib        |
	            inst_xy_holdrel   |
	            inst_xy_locate
	          );

	if (p->itype == instSpectroScan)
		cap |= inst_ref_xy;

	if (p->itype == instSpectroScanT)
		cap |= (inst_ref_xy       |
				inst_trans_spot );

	return cap;
}

/* gt_set_mode
 * set the next measurement mode.
 * Don't activate this mode until we actually have to.
 */
static inst_code
gt_set_mode(inst *pp, inst_mode m) {
	gretag *p = (gretag*)pp;
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
	return inst_ok;
}

/* Activate the last set mode if it different to the current mode. */
/* The Spectrolino needs calibration every time the mode actually changes. */
static inst_code
activate_mode(gretag *p)
{
	int err;
	unsigned char scanerror;
	char buf[255];
	unsigned short remoteerrorset;
	unsigned char execerror;

	if (((p->lastmode & inst_mode_illum_mask) == inst_mode_reflection
	  && (p->mode     & inst_mode_illum_mask) != inst_mode_reflection)
	 || ((p->lastmode & inst_mode_illum_mask) == inst_mode_emission
	  && (p->mode     & inst_mode_illum_mask) != inst_mode_emission)
	 || ((p->lastmode & inst_mode_illum_mask) == inst_mode_transmission
	  && (p->mode     & inst_mode_illum_mask) != inst_mode_transmission)) {
		inst *pp = (inst *)p;
		int iact = 0;		/* Interactive calibration */

		p->mode = p->lastmode;
		if (p->itype == instSpectroScan ||
			p->itype == instSpectroScanT) {
			iact = 1;		/* Can calibrate without user interaction */
		}

		/* Signal that a calibration will be needed. */
		p->need_ref_measurement = true;

		/* Activate some of mode in case some instructions validity depend on it */
		/* (How much mode ?) */
		if (p->itype == instSpectroScanT) {

			if ((p->mode & inst_mode_illum_mask) == inst_mode_reflection) {

				/* Set table mode to reflectance */
				err = p->spm_command(p, SPM_PROTO_SETTABLEMODE,
								 SPM_TABLEMODE_REFLECTANCE,
								 &scanerror);
				if (err != inst_ok) return err;
				if (scanerror != SPM_SCANERROR_NOERROR) {
					sprintf(buf,"gt_ref_measurement: got scan error %x",scanerror);
					p->logerr(buf);
					return scanerror_code(scanerror);
				}
			}
			else if ((p->mode & inst_mode_illum_mask) == inst_mode_transmission) {

				/* set table mode to transmission */
				err = p->spm_command(p, SPM_PROTO_SETTABLEMODE,
								 SPM_TABLEMODE_TRANSMISSION,
								 &scanerror);
				if (err != inst_ok) return err;
				if (scanerror != SPM_SCANERROR_NOERROR) {
					sprintf(buf,"gt_ref_measurement: got scan error %x",scanerror);
					p->logerr(buf);
					return scanerror_code(scanerror);
				}
			}
		}
	}

	return inst_ok;
}

/* 
 * set or reset an optional mode
 * We assume that the instrument has been initialised.
 */
static inst_code
gt_set_opt_mode(inst *pp, inst_opt_mode m)
{
	ggretag *p = (gretag *)pp;
	return inst_unsupported;
}

/* For an xy instrument, release the paper */
/* Return the inst error code */
static inst_code
gt_xy_sheet_release(
struct _inst *pp) {
	gretag *p = (gretag*)pp;
	int err;
	unsigned char execerror;
	char buf[255];

	/* Raise the sensor */
	err = p->spm_command(p, SPM_PROTO_MOVEUP, &execerror);
	if (err != inst_ok) return err;
	if (execerror != SPM_SCANERROR_NOERROR) {
		sprintf(buf,"gt_xy_sheet_release, raise sensor: got scan error %x",execerror);
		p->logerr(buf);
		return execerror_code(execerror);
	}
	
	/* Release the paper */
	err = p->spm_command(p, SPM_PROTO_RELEASEPAPER, &execerror);
	if (err != inst_ok) return err;
	if (execerror != SPM_SCANERROR_NOERROR) {
		sprintf(buf,"gt_xy_sheet_release, release paper: got scan error %x",execerror);
		p->logerr(buf);
		return execerror_code(execerror);
	}
	
	/* Move to the home position */
	err = p->spm_command(p, SPM_PROTO_MOVEHOME, &execerror);
	if (err != inst_ok) return err;
	if (execerror != SPM_SCANERROR_NOERROR) {
		sprintf(buf,"gt_xy_sheet_release, move to home: got scan error %x",execerror);
		p->logerr(buf);
		return execerror_code(execerror);
	}
	
	return inst_ok;
}

/* For an xy instrument, hold the paper */
/* Return the inst error code */
static inst_code 
gt_xy_sheet_hold(
struct _inst *pp) {
	gretag *p = (gretag*)pp;
	int err;
	unsigned char execerror;
	char buf[255];

	/* Release the paper */
	err = p->spm_command(p, SPM_PROTO_HOLDPAPER, &execerror);
	if (err != inst_ok) return err;
	if (execerror != SPM_SCANERROR_NOERROR) {
		sprintf(buf,"gt_xy_sheet_hold, hold paper: got scan error %x",execerror);
		p->logerr(buf);
		return execerror_code(execerror);
	}
	
	return inst_ok;
}

/* For an xy instrument, allow the user to locate a point */
/* Return the inst error code */
static inst_code
gt_xy_locate_start(
struct _inst *pp) {
	gretag *p = (gretag*)pp;
	int err;
	unsigned char execerror;
	char buf[255];

	/* The SpectroscanT can't do this in transmission mode, */
	/* so activate the next selected mode now. */
	if ((err = activate_mode(p)) != inst_ok)
		return err;

#ifdef NEVER
	/* Lower the sensor */
	err = p->spm_command(p, SPM_PROTO_MOVEDOWN, &execerror);
	if (err != inst_ok) return err;
	if (execerror != SPM_SCANERROR_NOERROR) {
		sprintf(buf,"gt_xy_locate_start, lower sensor: got scan error %x",execerror);
		p->logerr(buf);
		return execerror_code(execerror);
	}
#endif
	
	/* Put the device offline */
	err = p->spm_command(p, SPM_PROTO_SETDEVICEOFFLINE, &execerror);
	if (err != inst_ok) return err;
	if (execerror != SPM_SCANERROR_NOERROR) {
		sprintf(buf,"gt_xy_locate_start, set offline: got scan error %x",execerror);
		p->logerr(buf);
		return execerror_code(execerror);
	}

	return inst_ok;
}

/* For an xy instrument, read back the location */
/* Return the inst error code */
static inst_code
gt_xy_get_location(
struct _inst *pp,
double *x, double *y) {
	gretag *p = (gretag*)pp;
	int err;
	unsigned char reftype, zcoord;
	unsigned short dumy, xcoord, ycoord;
	char buf[255];

	/* Get the current sight position */
	err = p->spm_command(p, SPM_PROTO_OUTPUTACTUALPOSITION, SPM_REFERENCE_SIGHT,
	                     &reftype, &dumy, &xcoord, &ycoord, &zcoord);
	if (err != inst_ok) return err;
	
	*x = (double)xcoord;
	*y = (double)ycoord;

	return inst_ok;
}

/* For an xy instrument, ends allowing the user to locate a point */
/* Return the inst error code */
static inst_code
gt_xy_locate_end(
struct _inst *pp) {
	gretag *p = (gretag*)pp;
	int err;
	unsigned char execerror;
	char buf[255];

	/* The SpectroscanT can't do this in transmission mode, */
	/* so activate the next selected mode now. */
	if ((err = activate_mode(p)) != inst_ok)
		return err;

	/* Put the device online */
	err = p->spm_command(p, SPM_PROTO_SETDEVICEONLINE, &execerror);
	if (err != inst_ok) return err;
	if (execerror != SPM_SCANERROR_NOERROR) {
		sprintf(buf,"gt_xy_locate_end, set online: got scan error %x",execerror);
		p->logerr(buf);
		return execerror_code(execerror);
	}

#ifdef NEVER
	/* Raise the sensor */
	err = p->spm_command(p, SPM_PROTO_MOVEUP, &execerror);
	if (err != inst_ok) return err;
	if (execerror != SPM_SCANERROR_NOERROR) {
		sprintf(buf,"gt_xy_locate_end, raise sensor: got scan error %x",execerror);
		p->logerr(buf);
		return execerror_code(execerror);
	}
#endif

	return inst_ok;
}

/* For an xy instrument, try and clear the table after an abort */
/* Ignore any errors and retun OK */
static inst_code
gt_xy_clear(
struct _inst *pp) {
	gretag *p = (gretag*)pp;
	unsigned char execerror;

	/* Put the device online */
	p->spm_command(p, SPM_PROTO_SETDEVICEONLINE, &execerror);

	/* Raise the sensor */
	p->spm_command(p, SPM_PROTO_MOVEUP, &execerror);

	/* Release the paper */
	p->spm_command(p, SPM_PROTO_RELEASEPAPER, &execerror);
	
	/* Move to the home position */
	p->spm_command(p, SPM_PROTO_MOVEHOME, &execerror);

	return inst_ok;
}

/* Read a sheet full of patches using xy mode */
/* Return the inst error code */
static inst_code
gt_read_xy(
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
ipatch *vals) { 		/* Pointer to array of return values */
	gretag *p = (gretag*)pp;
	int pass, step, patch;
	inst_code err;
	unsigned char execerror;
	char buf[255];
	int tries, tc;			/* Total read tries */
	int fstep = 0;		/* NZ if step is fast direction */

	if (p->itype != instSpectroScan)
		return inst_unsupported;

	/* Activate the next selected mode */
	if ((err = activate_mode(p)) != inst_ok)
		return err;

	/* Move quickest in X direction to minimise noise */
	if (fabs(px) > fabs(ax))
		fstep = 1;

	tries = sip * pis;

	/* Read all the patches */
	for (step = pass = tc = 0; tc < tries; tc++) {
		int astep, apass;	/* Actual step and pass to use */
		unsigned short ix, iy;

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

			ix = (int)((ox + apass * ax + astep * px) + 0.5);
			iy = (int)((oy + apass * ay + astep * py) + 0.5);

			if ((step & 1) == 1) {	/* Offset for odd hex patches */
				ix += aax;
				iy += aay;
			}

//			printf("~1 patch = %d out of %d, pass %d, step %d\n",patch, npatch, patch, npatch,apass,astep);

			/* Do calibration if it is needed */
			if (p->need_ref_measurement) {
				if ((err = pp->calibrate(pp, 1)) != inst_ok)
					return err;
			}

#ifdef NEVER	/* Separate move and measure */
			/* Move to the right XY */
			err = p->spm_command(p, SPM_PROTO_MOVEABSOLUT, SPM_REFERENCE_SENSOR, ix, iy, &execerror);
			if (err != inst_ok) return err;
			if (execerror != SPM_SCANERROR_NOERROR) {
				sprintf(buf,"gt_read_xy, move absolute: got scan error %x",execerror);
				p->logerr(buf);
				return execerror_code(execerror);
			}
			
			vals[patch].XYZ_v = 0;
			vals[patch].aXYZ_v = 0;
			vals[patch].Lab_v = 0;
			vals[patch].spec_n = 0;

			/* Do a measurement and get the values */
			err = p->read_sample(pp, pname, &vals[patch]);
			if (err != inst_ok) return err;

#else /* !NEVER */
			{
				unsigned char spectype, refvalid;
				unsigned short remoteerrorset;
				float col[3], spec[36];
				int i;

				vals[patch].XYZ_v = false;
				vals[patch].aXYZ_v = false;
				vals[patch].Lab_v = false;
				vals[patch].spec_n = 0;
			
				/* move and measure gives us spectrum data anyway */
				err = p->spm_command(p, SPM_PROTO_MOVEANDMEASURE, ix, iy,
								 &spectype,		/* Will be 0, spectrum not density */
								 spec,
								 &refvalid,
								 &remoteerrorset);
				if (err != inst_ok) return err;
				if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
					sprintf(buf, "gt_move_and_measure: got remote error %x", remoteerrorset);
					p->logerr(buf);
					return remoteerrorset_code(remoteerrorset);
				}

				vals[patch].spec_n = 36;
				vals[patch].spec_wl_short = 380;
				vals[patch].spec_wl_long = 730;
				for (i = 0; i < INSTR_MAX_BANDS; i++)
					if (i < vals[patch].spec_n)
						vals[patch].spec[i] = 100.0 * (double)spec[i];
					else
						vals[patch].spec[i] = -1.0;

				/* Get the XYZ */
				if ((err = p->query_color(p, GT_COLOR_XYZ, col)) != inst_ok)
					return err;
				vals[patch].XYZ_v = true;
				vals[patch].XYZ[0] = col[0];
				vals[patch].XYZ[1] = col[1];
				vals[patch].XYZ[2] = col[2];

				/* Track the need for a calibration */
				p->incr_mcounter(p);
			}
#endif /* !NEVER */
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

	return inst_ok;
}


/* gt_read_strip()
 * read a set of strips
 */
static inst_code
gt_read_strip(inst *pp, char* name, int npatch, char* pname,
			  int sguide, double pwid, double gwid, double twid, ipatch* vals)
{
	return inst_unsupported;
}


/* gt_read_sample()
 * do a measurement and 
 * read a single sample.
 * We assume that calibration has been taken care of before calling this.
 */
static inst_code
gt_read_sample(inst *pp, char* name, ipatch *val)
{
	gretag *p = (gretag*)pp;
	float col[3], spec[36];
	int i, err;

	if (!val) return inst_internal_error;

	/* Configure for spot reading */
	p->lastmode = (p->lastmode & ~inst_mode_sub_mask) | inst_mode_spot;
	if ((err = activate_mode(p)) != inst_ok)
		return err;

	val->XYZ_v = false;
	val->aXYZ_v = false;
	val->Lab_v = false;
	val->spec_n = 0;

	if ((err = p->measurement(p)) != inst_ok)
		return err;

	if ((err = p->query_color(p, GT_COLOR_XYZ, col)) != inst_ok)
		return err;
	if ((p->mode & inst_mode_illum_mask) == inst_mode_emission) {
		val->aXYZ_v = true;
		val->aXYZ[0] = col[0];
		val->aXYZ[1] = col[1];
		val->aXYZ[2] = col[2];
	} else {
		val->XYZ_v = true;
		val->XYZ[0] = col[0];
		val->XYZ[1] = col[1];
		val->XYZ[2] = col[2];
	}

	/* spectrum data is returned only if requested */
	if (p->mode & inst_mode_spectral) {
		if ((err = p->query_spectrum(p, GT_SPECTYPE_REMISSION, spec)) != inst_ok)
			return err;

		val->spec_n = 36;
		val->spec_wl_short = 380;
		val->spec_wl_long = 730;
		for (i=0; i<INSTR_MAX_BANDS; i++)
			if (i < val->spec_n)
				val->spec[i] = 100.0 * (double)spec[i];
			else
				val->spec[i] = -1.0;
	}

	return inst_ok;
}

/* Determine if a calibration is needed. Returns inst_ok if not, */
/* inst_unsupported if it is unknown, or inst_needs_cal if needs calibration */
inst_code gt_needs_calibration(struct _inst *pp) {
	gretag *p = (gretag*)pp;

	if (p->need_ref_measurement)
		return inst_needs_cal;
	return inst_ok;
}


/* gt_calibrate
 * perform a reference measurement based on current mode.
 * Note this is usually used to activate the current mode.
 * If cal_no != 0, don't ask any questions for spectroscan.
 */
static inst_code
gt_calibrate(inst *pp, int cal_no)
{
	gretag *p = (gretag*)pp;
	unsigned short remoteerrorset;
	unsigned char control, execerror, scanerror;
	unsigned char dstd, wbase, illum, observer, actfilter;
	unsigned char origwhiteref, referencename[19];
	float whitespectrum[36];
	char buf[256];
	int cap, err;

	p->mode = p->lastmode;		/* In case not called from activate */
	p->need_ref_measurement = true;

	if ((p->mode & inst_mode_illum_mask) == inst_mode_reflection) {

		/* if spectroscan, set table mode to reflectance */
		if (p->itype == instSpectroScanT) {
			err = p->spm_command(p, SPM_PROTO_SETTABLEMODE,
							 SPM_TABLEMODE_REFLECTANCE,
							 &scanerror);
			if (err != inst_ok) return err;
			if (scanerror != SPM_SCANERROR_NOERROR) {
				sprintf(buf,"gt_ref_measurement: got scan error %x",scanerror);
				p->logerr(buf);
				return scanerror_code(scanerror);
			}
		}

		/* (redundant?) set remission mode */
		err = p->spm_command(p, SPM_PROTO_MEASCONTROLDOWNLOAD,
						 SPM_CONTROL_REMISSIONMEAS, &remoteerrorset);
		if (err != inst_ok) return err;
		if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
			sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
			p->logerr(buf);
			return remoteerrorset_code(remoteerrorset);
		}

		/* verify white reference */
		err = p->spm_command(p, SPM_PROTO_WHITEREFERENCEREQUEST,
						 SPM_FILTER(p->filter),
						 &actfilter, whitespectrum, &origwhiteref, referencename);
		if (err != inst_ok) return err;
		sprintf(buf, "gt_ref_measurement: white reference == %s", referencename);
		p->logwarn(buf);
		if (p->itype == instSpectrolino) {
			sprintf(buf, "Place device on white reference %s", referencename);
			if (p->confirm(buf, "OK", "Cancel", NULL) != 1)
				return inst_user_abort;
		} else if (cal_no == 0) { /* SPECTROSCAN/T */
			sprintf(buf, "Please confirm white reference %s is in slot 1", referencename);
			if (p->confirm(buf, "OK", "Cancel", NULL) != 1)
				return inst_user_abort;
		}

		/* attempt white calibration until successful */
		for (;;) {
			err = p->spm_command(p, SPM_PROTO_PARAMETERDOWNLOAD,
							 SPM_DSTD(p->dstd),
							 SPM_WBASE_ABS,
							 SPM_ILLUM(p->illum),
							 SPM_OBSERVER(p->observer),
							 &remoteerrorset);
			if (err != inst_ok) return err;
			if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
				sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
				p->logerr(buf);
				return remoteerrorset_code(remoteerrorset);
			}

			/* if spectroscan, move to white ref */
			if (p->itype == instSpectroScan ||
				p->itype == instSpectroScanT) {
				err = p->spm_command(p, SPM_PROTO_MOVETOWHITEREFPOS, 0, &scanerror);
				if (err != inst_ok) return err;
				if (scanerror != SPM_SCANERROR_NOERROR) {
					sprintf(buf,"gt_ref_measurement: got scan error %x",scanerror);
					p->logerr(buf);
					return scanerror_code(scanerror);
				}
				err = p->spm_command(p, SPM_PROTO_MOVEDOWN, &scanerror);
				if (err != inst_ok) return err;
				if (scanerror != SPM_SCANERROR_NOERROR) {
					sprintf(buf,"gt_ref_measurement: got scan error %x",scanerror);
					p->logerr(buf);
					return scanerror_code(scanerror);
				}
			}

			/* calibrate */
			err = p->spm_command(p, SPM_PROTO_EXECREFMEASUREMENT,
							 SPM_MEASUREMENTMODE_WHITECALWITHWARN,
							 &execerror);
			if (err != inst_ok) return err;
			if (execerror !=  SPM_ERROR_WHITEMEASOK) {
				p->logerr("gt_ref_measurement: white reference measurement failed");
				return execerror_code(execerror);
			}

			/* if spectroscan, raise up */
			if (p->itype == instSpectroScan ||
				p->itype == instSpectroScanT) {
				err = p->spm_command(p, SPM_PROTO_MOVEUP, &scanerror);
				if (err != inst_ok) return err;
				if (scanerror != SPM_SCANERROR_NOERROR) {
					sprintf(buf,"gt_ref_measurement: got scan error %x",scanerror);
					p->logerr(buf);
					return scanerror_code(scanerror);
				}
			}

			/* verify filter */
			err = p->spm_command(p, SPM_PROTO_PARAMETERREQUEST,
							 &dstd, &wbase, &illum, &observer, &actfilter);
			if (err != inst_ok) return err;
			if (SPM_FILTER(p->filter) == actfilter)
				break;  /* success */

			sprintf(buf, "gt_ref_measurement: filter mismatch");
			p->logwarn(buf);
			sprintf(buf, "Replace filter with %s", TXT_FILTER(p->filter));
			if (p->confirm(buf, "OK", "Cancel", NULL) != 1)
				return inst_user_abort;
		}

		/* (re)set requested parameters */
		if (dstd != SPM_DSTD(p->dstd) ||
			wbase != SPM_WBASE(p->whitebase) ||
			illum != SPM_ILLUM(p->illum) ||
			observer != SPM_OBSERVER(p->observer))
		{
			err = p->spm_command(p, SPM_PROTO_PARAMETERDOWNLOAD,
							 SPM_DSTD(p->dstd),
							 SPM_WBASE(p->whitebase),
							 SPM_ILLUM(p->illum),
							 SPM_OBSERVER(p->observer),
							 &remoteerrorset);
			if (err != inst_ok) return err;
			if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
				sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
				p->logerr(buf);
				return remoteerrorset_code(remoteerrorset);
			}
		}

	} /* end remission */

	else if ((p->mode & inst_mode_illum_mask) == inst_mode_emission) {

		/* set remission mode for reference reading */
		err = p->spm_command(p, SPM_PROTO_MEASCONTROLDOWNLOAD,
						 SPM_CONTROL_REMISSIONMEAS, &remoteerrorset);
		if (err != inst_ok) return err;
		if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
			sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
			p->logerr(buf);
			return remoteerrorset_code(remoteerrorset);
		}

		/* verify white reference */
		err = p->spm_command(p, SPM_PROTO_WHITEREFERENCEREQUEST,
						 SPM_FILTER(p->filter),
						 &actfilter, whitespectrum, &origwhiteref, referencename);
		if (err != inst_ok) return err;
		sprintf(buf, "gt_ref_measurement: white reference == %s", referencename);
		p->logwarn(buf);
		sprintf(buf, "Install NoFilter and place device on white reference %s", referencename);
		if (p->confirm(buf, "OK", "Cancel", NULL) != 1)
			return inst_user_abort;

		/* REQUIRED for emission mode */
		p->filter = GT_FILTER_NONE;
		p->whitebase = GT_WHITEBASE_ABSOLUTE;

		/* perform white calibration */
		for (;;) {
			err = p->spm_command(p, SPM_PROTO_PARAMETERDOWNLOAD,
							 SPM_DSTD(p->dstd),
							 SPM_WBASE_ABS,
							 SPM_ILLUM(p->illum),
							 SPM_OBSERVER(p->observer),
							 &remoteerrorset);
			if (err != inst_ok) return err;
			if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
				sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
				p->logerr(buf);
				return remoteerrorset_code(remoteerrorset);
			}

			err = p->spm_command(p, SPM_PROTO_EXECREFMEASUREMENT,
							 SPM_MEASUREMENTMODE_WHITECALWITHWARN,
							 &execerror);
			if (err != inst_ok) return err;
			if (execerror !=  SPM_ERROR_WHITEMEASOK) {
				p->logerr("gt_ref_measurement: white reference measurement failed");
				return execerror_code(execerror);
			}

			/* verify filter */
			err = p->spm_command(p, SPM_PROTO_PARAMETERREQUEST,
							 &dstd, &wbase, &illum, &observer, &actfilter);
			if (err != inst_ok) return err;
			if (SPM_FILTER(p->filter) == actfilter)
				break;  /* success */

			sprintf(buf, "gt_ref_measurement: filter mismatch");
			p->logwarn(buf);
			sprintf(buf, "Replace filter with %s", TXT_FILTER(p->filter));
			if (p->confirm(buf, "OK", "Cancel", NULL) != 1)
				return inst_user_abort;
		}

		/* set emission mode */
		err = p->spm_command(p, SPM_PROTO_MEASCONTROLDOWNLOAD,
						 SPM_CONTROL_EMISSIONMEAS, &remoteerrorset);
		if (err != inst_ok) return err;
		if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
			sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
			p->logerr(buf);
			return remoteerrorset_code(remoteerrorset);
		}

		/* perform dark measurement */
		err = p->spm_command(p, SPM_PROTO_EXECREFMEASUREMENT,
						 SPM_MEASUREMENTMODE_EMISSIONCAL,
						 &execerror);
		if (err != inst_ok) return err;
		if (execerror !=  SPM_ERROR_EMISSIONCALOK) {
			p->logerr("gt_ref_measurement: dark reference measurement failed");
			return execerror_code(execerror);
		}

		/* set photometric mode (absolute or relative) */
		err = p->spm_command(p, SPM_PROTO_MEASCONTROLDOWNLOAD,
						 SPM_PHOTOMODE(p->photomode), &remoteerrorset);
		if (err != inst_ok) return err;
		if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
			sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
			p->logerr(buf);
			return remoteerrorset_code(remoteerrorset);
		}

		if (p->photomode == GT_PHOTOMODE_RELATIVE) {
			err = p->spm_command(p, SPM_PROTO_FLOATDOWNLOAD,
							 SPM_COMFLOAT_VPHOTOMETRICYREF,
							 p->photoref,
							 &remoteerrorset);
			if (err != inst_ok) return err;
			if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
				sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
				p->logerr(buf);
				return remoteerrorset_code(remoteerrorset);
			}
		}

	} /* end emission */

	else if ((p->mode & inst_mode_illum_mask) == inst_mode_transmission) {

		if (p->itype != instSpectroScanT) {
			p->logerr("gt_ref_measurement: transmission supported on spectroscan/t device only");
			return inst_unsupported;
		}

		/* REQUIRED for transmission mode */
		p->filter = GT_FILTER_NONE;
		p->whitebase = GT_WHITEBASE_ABSOLUTE;

		sprintf(buf, "Please clear the light table and replace filter with %s", TXT_FILTER(p->filter));
		if (p->confirm(buf, "OK", "Cancel", NULL) != 1)
			return inst_user_abort;

		/* set table mode to transmission */
		err = p->spm_command(p, SPM_PROTO_SETTABLEMODE,
						 SPM_TABLEMODE_TRANSMISSION,
						 &scanerror);
		if (err != inst_ok) return err;
		if (scanerror != SPM_SCANERROR_NOERROR) {
			sprintf(buf,"gt_ref_measurement: got scan error %x",scanerror);
			p->logerr(buf);
			return scanerror_code(scanerror);
		}

		/* attempt white calibration until successful */
		for (;;) {
			err = p->spm_command(p, SPM_PROTO_PARAMETERDOWNLOAD,
							 SPM_DSTD(p->dstd),
							 SPM_WBASE_ABS,
							 SPM_ILLUM(p->illum),
							 SPM_OBSERVER(p->observer),
							 &remoteerrorset);
			if (err != inst_ok) return err;
			if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
				sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
				p->logerr(buf);
				return remoteerrorset_code(remoteerrorset);
			}

			/* calibrate */
			err = p->spm_command(p, SPM_PROTO_EXECREFMEASUREMENT,
							 SPM_MEASUREMENTMODE_WHITECALWITHWARN,
							 &execerror);
			if (err != inst_ok) return err;
			if (execerror !=  SPM_ERROR_WHITEMEASOK) {
				p->logerr("gt_ref_measurement: white reference measurement failed");
				return execerror_code(execerror);
			}

			/* verify filter */
			err = p->spm_command(p, SPM_PROTO_PARAMETERREQUEST,
							 &dstd, &wbase, &illum, &observer, &actfilter);
			if (err != inst_ok) return err;
			if (SPM_FILTER(p->filter) == actfilter)
				break;  /* success */

			sprintf(buf, "gt_ref_measurement: filter mismatch");
			p->logwarn(buf);
			sprintf(buf, "Replace filter with %s", TXT_FILTER(p->filter));
			if (p->confirm(buf, "OK", "Cancel", NULL) != 1)
				return inst_user_abort;
		}

		/* (re)set requested parameters */
		if (dstd != SPM_DSTD(p->dstd) ||
			wbase != SPM_WBASE(p->whitebase) ||
			illum != SPM_ILLUM(p->illum) ||
			observer != SPM_OBSERVER(p->observer))
		{
			err = p->spm_command(p, SPM_PROTO_PARAMETERDOWNLOAD,
							 SPM_DSTD(p->dstd),
							 SPM_WBASE(p->whitebase),
							 SPM_ILLUM(p->illum),
							 SPM_OBSERVER(p->observer),
							 &remoteerrorset);
			if (err != inst_ok) return err;
			if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
				sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
				p->logerr(buf);
				return remoteerrorset_code(remoteerrorset);
			}
		}

	} /* end transmission */

	else {
		p->logerr("gt_ref_measurement: unknown measurement mode");
		return inst_unsupported;
	}

	if (p->itype == instSpectrolino || cal_no == 0) {
		if (p->confirm("Ready to measure.", "OK", NULL, NULL) != 1)
			return inst_user_abort;
	}

	p->mcounter = 0;
	p->need_ref_measurement = false;
	return inst_ok;
}


/* gt_command()
 * do a command/response with the instrument
 */
static inst_code
gt_command(inst *pp, char* in, char* out, int bsize)
{
	return inst_unsupported;
}


/* gt_interp_error()
 * instrument specific error codes interpretation
 */
static char*
gt_interp_error(inst *pp, int ec)
{
	return "interp_error not yet implemented";
}


/* gt_last_sioerr()
 * return the last serial I/O error code
 */
static int
gt_last_sioerr(inst *pp) {
	return pp->sio->lerr;
}


/* gt_del
 * free associated resources
 */
static void
gt_del(inst *pp)
{
	gretag *p = (gretag *)pp;
	if (p->sio != NULL)
		p->sio->del(p->sio);
	free(p);
}

/**********************************************************************
 * Gretag specific methods
 **********************************************************************/

/* gt_measurement
 * perform a measurement
 * Assumes mode has been activated. 
 */
static inst_code
gt_measurement(gretag *p)
{
	int err;
	unsigned char execerror;
	char buf[255];

	/* do we need to calibrate? */
	if (p->need_ref_measurement)
		return inst_needs_cal;

	/* if spectroscan and mode is emission, move down */
	if ((p->itype == instSpectroScan ||
		 p->itype == instSpectroScanT) &&
		(p->mode & inst_mode_illum_mask) == inst_mode_reflection)
	{
		err = p->spm_command(p, SPM_PROTO_MOVEDOWN, &execerror);
		if (err != inst_ok) return err;
		if (execerror != SPM_SCANERROR_NOERROR) {
			sprintf(buf,"gt_ref_measurement: got scan error %x",execerror);
			p->logerr(buf);
			return execerror_code(execerror);
		}
	}

	/* measure */
	err = p->spm_command(p, SPM_PROTO_EXECMEASUREMENT, &execerror);
	if (err != inst_ok) return err;
	if (execerror != SPM_ERROR_NOERROR) {
		sprintf(buf, "gt_measurement: failed measurement");
		p->logerr(buf);
		return execerror_code(execerror);
	}

	/* if spectroscan and mode is remission, raise up */
	if ((p->itype == instSpectroScan ||
		 p->itype == instSpectroScanT) &&
		(p->mode & inst_mode_illum_mask) == inst_mode_reflection)
	{
		err = p->spm_command(p, SPM_PROTO_MOVEUP, &execerror);
		if (err != inst_ok) return err;
		if (execerror != SPM_SCANERROR_NOERROR) {
			sprintf(buf,"gt_ref_measurement: got scan error %x",execerror);
			p->logerr(buf);
			return execerror_code(execerror);
		}
	}

	p->incr_mcounter(p);
	return inst_ok;
}


/* gt_query_color
 * query color values from last measurement
 */
static inst_code
gt_query_color(gretag *p, gt_color color, float *cvec)
{
	unsigned char ctype, refvalid, actfilter, wbase, illum, observer;
	unsigned short remoteerrorset;
	char buf[256];
	int err;

	err = p->spm_command(p, SPM_PROTO_CPARAMETERREQUEST,
					 SPM_COLOR(color),
					 &ctype,
					 cvec,
					 &refvalid,
					 &actfilter,
					 &wbase,
					 &illum,
					 &observer,
					 &remoteerrorset);
	if (err != inst_ok) return err;
	if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
		sprintf(buf, "gt_query_color: got remote error %x", remoteerrorset);
		p->logerr(buf);
		return remoteerrorset_code(remoteerrorset);
	}
	return inst_ok;
}


/* gt_query_spectrum
 * query spectrum values from last measurement
 */
static inst_code
gt_query_spectrum(gretag *p, gt_spectype spec, float *svec)
{
	unsigned char spectype, refvalid;
	unsigned short remoteerrorset;
	char buf[256];
	int err;

	err = p->spm_command(p, SPM_PROTO_SPECTRUMREQUEST,
					 SPM_SPECTYPE(spec),
					 &spectype,
					 svec,
					 &refvalid,
					 &remoteerrorset);
	if (err != inst_ok) return err;
	if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
		sprintf(buf, "gt_query_spectrum: got remote error %x", remoteerrorset);
		p->logerr(buf);
		return remoteerrorset_code(remoteerrorset);
	}
	return inst_ok;
}



/* gt_set_whitebase
 * set white base for subsequent measurements
 */
static inst_code
gt_set_whitebase(gretag *p, gt_whitebase wb)
{
	int err;
	unsigned short remoteerrorset;
	char buf[255];

	p->whitebase = wb;

	err = p->spm_command(p, SPM_PROTO_PARAMETERDOWNLOAD,
					 SPM_DSTD(p->dstd),
					 SPM_WBASE(p->whitebase),
					 SPM_ILLUM(p->illum),
					 SPM_OBSERVER(p->observer),
					 &remoteerrorset);
	if (err != inst_ok) return err;
	if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
		sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
		p->logerr(buf);
		return remoteerrorset_code(remoteerrorset);
	}

	return inst_ok;
}

static gt_whitebase
gt_get_whitebase(gretag *p)
{
	return p->whitebase;
}


/* gt_set_observer
 * set observer for subsequent measurements
 */
static inst_code
gt_set_observer(gretag *p, gt_observer ob)
{
	int err;
	unsigned short remoteerrorset;
	char buf[255];

	p->observer = ob;

	err = p->spm_command(p, SPM_PROTO_PARAMETERDOWNLOAD,
					 SPM_DSTD(p->dstd),
					 SPM_WBASE(p->whitebase),
					 SPM_ILLUM(p->illum),
					 SPM_OBSERVER(p->observer),
					 &remoteerrorset);
	if (err != inst_ok) return err;
	if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
		sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
		p->logerr(buf);
		return remoteerrorset_code(remoteerrorset);
	}

	return inst_ok;
}

static gt_observer
gt_get_observer(gretag *p)
{
	return p->observer;
}


/* gt_set_illum
 * set illum for subsequent measurements
 */
static inst_code
gt_set_illum(gretag *p, gt_illum il)
{
	int err;
	unsigned short remoteerrorset;
	char buf[255];

	p->illum = il;

	err = p->spm_command(p, SPM_PROTO_PARAMETERDOWNLOAD,
					 SPM_DSTD(p->dstd),
					 SPM_WBASE(p->whitebase),
					 SPM_ILLUM(p->illum),
					 SPM_OBSERVER(p->observer),
					 &remoteerrorset);
	if (err != inst_ok) return err;
	if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
		sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
		p->logerr(buf);
		return remoteerrorset_code(remoteerrorset);
	}

	return inst_ok;
}

static gt_illum
gt_get_illum(gretag *p)
{
	return p->illum;
}


/* gt_set_dstd
 * set dstd for subsequent measurements
 */
static inst_code
gt_set_dstd(gretag *p, gt_dstd d)
{
	int err;
	unsigned short remoteerrorset;
	char buf[255];

	p->dstd = d;

	err = p->spm_command(p, SPM_PROTO_PARAMETERDOWNLOAD,
					 SPM_DSTD(p->dstd),
					 SPM_WBASE(p->whitebase),
					 SPM_ILLUM(p->illum),
					 SPM_OBSERVER(p->observer),
					 &remoteerrorset);
	if (err != inst_ok) return err;
	if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
		sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
		p->logerr(buf);
		return remoteerrorset_code(remoteerrorset);
	}

	return inst_ok;
}

static gt_dstd
gt_get_dstd(gretag *p)
{
	return p->dstd;
}


/* gt_set_photomode (absolute or relative)
 * set photomode for subsequent measurements
 */
static inst_code
gt_set_photomode(gretag *p, gt_photomode pm)
{
	int err;
	unsigned short remoteerrorset;
	char buf[255];

	p->photomode = pm;

	err = p->spm_command(p, SPM_PROTO_MEASCONTROLDOWNLOAD,
					 SPM_PHOTOMODE(p->photomode), &remoteerrorset);
	if (err != inst_ok) return err;
	if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
		sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
		p->logerr(buf);
		return remoteerrorset_code(remoteerrorset);
	}

	return inst_ok;
}

static gt_photomode
gt_get_photomode(gretag *p)
{
	return p->photomode;
}


/* gt_set_photoref
 * set photoref for subsequent measurements
 */
static inst_code
gt_set_photoref(gretag *p, float pr)
{
	int err;
	unsigned short remoteerrorset;
	char buf[255];

	p->photoref = pr;

	err = p->spm_command(p, SPM_PROTO_FLOATDOWNLOAD,
					 SPM_COMFLOAT_VPHOTOMETRICYREF,
					 p->photoref,
					 &remoteerrorset);
	if (err != inst_ok) return err;
	if (remoteerrorset != SPM_REMOTEERROR_NOERROR) {
		sprintf(buf,"gt_ref_measurement: got remote error %x",remoteerrorset);
		p->logerr(buf);
		return remoteerrorset_code(remoteerrorset);
	}

	return inst_ok;
}

static float
gt_get_photoref(gretag *p)
{
	return p->photoref;
}


/* gt_set_filter
 * set desired filter for subsequent measurements, require new ref measurement
 */
static inst_code
gt_set_filter(gretag *p, gt_filter f)
{
	p->filter = f;
	p->need_ref_measurement = true;
	return inst_ok;
}

static gt_filter
gt_get_filter(gretag *p)
{
	return p->filter;
}


/* gt_set_logwarn
 * set warning callback
 */
static inst_code
gt_set_logwarn(gretag *p, int (*f)(char*))
{
	p->logwarn = f;
	return inst_ok;
}

static inst_code
gt_set_logerr(gretag *p, int (*f)(char*))
{
	p->logerr = f;
	return inst_ok;
}

static inst_code
gt_set_confirm(gretag *p, int (*f)(char*, char*, char*, char*))
{
	p->confirm = f;
	return inst_ok;
}


/**********************************************************************
 * Gretag private methods
 **********************************************************************/

/* scanerror_code()
 * convert scanerror to inst_code
 */
static inst_code
scanerror_code(unsigned char scanerror)
{
	switch (scanerror) {
	case SPM_SCANERROR_NOERROR:
		return inst_ok;
	case SPM_SCANERROR_DEVICEISOFFLINE:
	case SPM_SCANERROR_NODEVICEFOUND:
	case SPM_SCANERROR_NOUSERACCESS:
		return inst_coms_fail;
	case SPM_SCANERROR_NOVALIDCOMMAND:
	case SPM_SCANERROR_OUTOFRANGE:
		return inst_protocol_error;
	case SPM_SCANERROR_MEASUREMENTERROR:
		return inst_misread;
	case SPM_SCANERROR_NOTRANSMTABLE:
		return inst_unknown_model;
	case SPM_SCANERROR_PROGRAMMINGERROR:
	case SPM_SCANERROR_NOTINTRANSMMODE:
	case SPM_SCANERROR_NOTINREFLECTMODE:
	default:
		return inst_other_error;
	}
}


/* remoteerrorset_code()
 * convert remoteerrorset to inst_code
 */
static inst_code
remoteerrorset_code(unsigned short remoteerrorset)
{
	switch(remoteerrorset) {
	case SPM_REMOTEERROR_NOERROR:
		return inst_ok;
	case SPM_REMOTEERROR_SLOPEOUTOFRANGE:
	case SPM_REMOTEERROR_DORLOUTOFRANGE:
	case SPM_REMOTEERROR_REFLECTANCEOUTOFRANGE:
	case SPM_REMOTEERROR_COLOR1OUTOFRANGE:
	case SPM_REMOTEERROR_COLOR2OUTOFRANGE:
	case SPM_REMOTEERROR_COLOR3OUTOFRANGE:
	case SPM_REMOTEERROR_NOTANSRORBOOLEAN:
	case SPM_REMOTEERROR_NOVALIDVALORREF:
	case SPM_REMOTEERROR_NOVALIDDSTD:
	case SPM_REMOTEERROR_NOVALIDWHITE:
	case SPM_REMOTEERROR_NOVALIDILLUM:
	case SPM_REMOTEERROR_NOVALIDOBSERVER:
	case SPM_REMOTEERROR_NOVALIDMAXLAMBDA:
	case SPM_REMOTEERROR_NOVALIDSPECT:
	case SPM_REMOTEERROR_NOVALIDCOLSYSORINDEX:
	case SPM_REMOTEERROR_NOVALIDCHAR:
		return inst_protocol_error;
	default:
		return inst_other_error;
	}
}

/* execerror_code()
 * convert execerror to inst_code
 */
static inst_code
execerror_code(unsigned char execerror)
{
	switch (execerror) {
	case SPM_ERROR_NOERROR:
	case SPM_ERROR_WHITEMEASOK:
	case SPM_ERROR_EMISSIONCALOK:
	case SPM_ERROR_RESETDONE:
		return inst_ok;
	case SPM_ERROR_MEMORYFAILURE:
	case SPM_ERROR_POWERFAILURE:
	case SPM_ERROR_LAMPFAILURE:
	case SPM_ERROR_HARDWAREFAILURE:
	case SPM_ERROR_FILTEROUTOFPOS:
	case SPM_ERROR_DRIVEERROR:
	case SPM_ERROR_EPROMFAILURE:
	case SPM_ERROR_MEMORYERROR:
	case SPM_ERROR_FULLMEMORY:
	case SPM_ERROR_BACKUPERROR:
	case SPM_ERROR_PROGRAMROMERROR1:
	case SPM_ERROR_PROGRAMROMERROR2:
	case SPM_ERROR_PROGRAMROMERROR3:
	case SPM_ERROR_PROGRAMROMERROR4:
	case SPM_ERROR_PROGRAMROMERROR5:
	case SPM_ERROR_PROGRAMROMERROR6:
	case SPM_ERROR_PROGRAMROMERROR7:
	case SPM_ERROR_PROGRAMROMERROR8:
		return inst_hardware_fail;
	case SPM_ERROR_SENDTIMEOUT:
	case SPM_ERROR_CHECKSUMWRONG:
		return inst_coms_fail;
	case SPM_ERROR_MEASDISABLED:
	case SPM_ERROR_REMOVERFLOW:
	case SPM_ERROR_NOTREADY:
		return inst_misread;
	case SPM_ERROR_DENSCALERROR:
	case SPM_ERROR_WHITEMEASWARN:
	case SPM_ERROR_ONLYEMISSION:
	case SPM_ERROR_NOVALIDMEAS:
	default:
		return inst_other_error;
	}
}


/* gt_incr_mcounter()
 * increment internal measurement counter, and
 * decide if a new reference measurement is needed.
 * gretag policy is to recalibrate after every 50
 * remission measurements or every 10 transmission.
 */
static inst_code
gt_incr_mcounter(gretag *p)
{
	p->mcounter++;
	if ( ((p->mode & inst_mode_illum_mask) == inst_mode_reflection && p->mcounter >= 50) ||
		 /*	((p->mode & inst_mode_illum_mask) == inst_mode_emission && p->mcounter >= 50) ||   what's the actual policy ? It's inconvenient when calibrating a monitor in any case ! */
		 ((p->mode & inst_mode_illum_mask) == inst_mode_transmission && p->mcounter >= 10) )
	{
		p->need_ref_measurement = true;
	}
	return inst_ok;
}

/* gt_spm_command()
 * format a spm command and send to the device,
 * then receive and parse the reply
 * vararg list includes first the args for the command
 * followed by the args for the reply
 */
static inst_code
gt_spm_command(gretag *p, int cmd, ...)
{
	unsigned char wbuf[SPM_BUFFER_SZ], rbuf[SPM_BUFFER_SZ];
	int reply = spm_expect_reply(cmd);

	int err;
	va_list ap;

	p->unexrep = SPM_PROTO_NONE;
	p->unexarg = -1;

	va_start(ap, cmd);
	err = spm_commandv(wbuf, cmd, &ap);
	if (!err) {
		sprintf(wbuf, "gt_spm_command: problem formatting command %s", spm_cmdname(cmd));
		p->logerr(wbuf);
		va_end(ap);
		return inst_protocol_error;
	}
	if ( 0 != (err = p->sio->write_read(p->sio, wbuf, rbuf, SPM_BUFFER_SZ, '\n', 1, spm_timeout(cmd)))) {
		sprintf(wbuf, "gt_spm_command: serial i/o error %x in command/reply %s/%s\n",
				err, spm_cmdname(cmd), spm_cmdname(reply));
		p->logerr(wbuf);
		va_end(ap);
		return inst_coms_fail;
	}
	err = spm_replyv(rbuf, reply, &ap);
	if (err == 1) {

		/* Get the actual reply type */
		p->unexrep = spm_reply_type(rbuf);

		sprintf(wbuf, "gt_spm_command: Got unexpected reply, expected %s, got %s",
		        spm_cmdname(reply), spm_cmdname(p->unexrep));
		p->logerr(wbuf);

		/* If it's an error type, get the one argument */
		if (   p->unexrep == SPM_PROTO_EXECERROR
		    || p->unexrep == SPM_PROTO_EXECERROR
	        || p->unexrep == SPM_PROTO_COMERR
	        || p->unexrep == SPM_PROTO_DOWNLOADERROR
	        || p->unexrep == SPM_PROTO_ACTERRORANSWER
	        || p->unexrep == SPM_PROTO_ERRORANSWER
	        || p->unexrep == SPM_PROTO_SSCOMERR) {
			unsigned char arg_char;

			err = spm_reply(rbuf, p->unexrep, &arg_char);
			p->unexarg = arg_char;
		}

		va_end(ap);
		return inst_unexpected_reply;
	}

	if (err != 2) {
		sprintf(rbuf, "gt_spm_command: problem parsing reply %s", spm_cmdname(reply));
		p->logerr(rbuf);
		va_end(ap);
		return inst_protocol_error;
	}

	va_end(ap);
	return inst_ok;
}

/* gt_unexpected_reply()
 * Return the spm code of the last recieved, unexpected reply.
 * If there was no unexpected reply to the last command,
 * return SPM_PROTO_NONE.
 */
static int
gt_unexpected_reply(gretag *p)
{
	return p->unexrep;
}

/* gt_unexpected_argument()
 * Return the value of the argument to the unexpected
 * reply, if it was one of 
    SPM_PROTO_EXECERROR
    SPM_PROTO_COMERR
    SPM_PROTO_DOWNLOADERROR
    SPM_PROTO_ACTERRORANSWER
    SPM_PROTO_ERRORANSWER
    SPM_PROTO_SSCOMERR
 * otherwise return -1
 */
static int
gt_unexpected_argument(gretag *p)
{
	switch (p->unexrep) {
    	case SPM_PROTO_EXECERROR:
    	case SPM_PROTO_COMERR:
    	case SPM_PROTO_DOWNLOADERROR:
    	case SPM_PROTO_ACTERRORANSWER:
    	case SPM_PROTO_ERRORANSWER:
    	case SPM_PROTO_SSCOMERR:
			return p->unexarg;
			break;
	}
	return -1;
}


static int
gt_default_log(char *msg)
{
	printf("%s\n", msg);
	fflush(stdout);
	return inst_ok;
}

static int
gt_default_confirm(char *msg, char *b1, char *b2, char *b3)
{
	int val = 0;
	char c;

	while (val < 1 || val > 3) {
		printf("%s\n", msg);
		if (b1) printf("1)%s\n", b1);
		if (b2) printf("2)%s\n", b2);
		if (b3) printf("3)%s\n", b3);
		printf("? ");
		fflush(stdout);
		empty_con_chars();
		c = next_con_char();
		printf("\n");
		val = (c == '1') ? 1 : (c == '2') ? 2 : (c == '3') ? 3 : 0;
	}
	return val;
}

/* gt_new_gretag
 * create gretag object
 */
gretag* new_gretag(void)
{
	gretag *p;
	if ((p = (gretag *)calloc(sizeof(gretag),1)) == NULL)
		error("gretag: malloc failed!");

	p->sio = new_serio();

	/* inst methods */
	p->init_coms    = gt_init_coms;
	p->init_inst    = gt_init_inst;
	p->capabilities = gt_capabilities;
	p->set_mode     = gt_set_mode;
	p->set_opt_mode     = gt_set_opt_mode;
	p->xy_sheet_release = gt_xy_sheet_release;
	p->xy_sheet_hold    = gt_xy_sheet_hold;
	p->xy_locate_start  = gt_xy_locate_start;
	p->xy_get_location  = gt_xy_get_location;
	p->xy_locate_end  = gt_xy_locate_end;
	p->xy_clear     = gt_xy_clear;
	p->read_xy      = gt_read_xy;
	p->read_strip   = gt_read_strip;
	p->read_sample  = gt_read_sample;
	p->needs_calibration = gt_needs_calibration;
	p->calibrate    = gt_calibrate;
	p->command      = gt_command;
	p->interp_error = gt_interp_error;
	p->inst_interp_error = NULL;            /* virtual constructor will do this */
	p->last_sioerr  = gt_last_sioerr;
	p->del          = gt_del;

	/* gretag specific methods */
	p->measurement = gt_measurement;
	p->query_color = gt_query_color;
	p->query_spectrum = gt_query_spectrum;
	p->set_whitebase = gt_set_whitebase;
	p->get_whitebase = gt_get_whitebase;
	p->set_observer = gt_set_observer;
	p->get_observer = gt_get_observer;
	p->set_illum = gt_set_illum;
	p->get_illum = gt_get_illum;
	p->set_dstd = gt_set_dstd;
	p->get_dstd = gt_get_dstd;
	p->set_photomode = gt_set_photomode;
	p->get_photomode = gt_get_photomode;
	p->set_photoref = gt_set_photoref;
	p->get_photoref = gt_get_photoref;
	p->set_filter = gt_set_filter;
	p->get_filter = gt_get_filter;
	p->set_logwarn = gt_set_logwarn;
	p->set_logerr = gt_set_logerr;
	p->set_confirm = gt_set_confirm;
	p->incr_mcounter = gt_incr_mcounter;
	p->spm_command = gt_spm_command;
	p->unexpected_reply = gt_unexpected_reply;
	p->unexpected_argument = gt_unexpected_argument;
	p->logerr = gt_default_log;
	p->logwarn = gt_default_log;
	p->confirm = gt_default_confirm;

	/* default state */
	p->itype = instUnknown;						/* Unknown until initialised */
	p->cap = inst_unknown;						/* Unknown until initialised */
	p->mode = inst_mode_unknown;				/* Not in a known mode yet */

	p->lastmode = inst_mode_ref_spot;
	p->mode = inst_mode_ref_spot;
	p->whitebase = GT_WHITEBASE_ABSOLUTE;
	p->observer = GT_OBSERVER_2;
	p->illum = GT_ILLUM_D50;
	p->dstd = GT_DSTD_ANSIA;
	p->photomode = GT_PHOTOMODE_ABSOLUTE;
	p->photoref = 80.0;
	p->filter = GT_FILTER_NONE;
	p->mcounter = 0;
	p->unexrep = SPM_PROTO_NONE;
	p->unexarg = 0;
	p->need_ref_measurement = true;       	/* needs ref meas at startup */ 
	p->lastmode |= inst_mode_spectral;		/* Return spectral by default */
	p->mode |= inst_mode_spectral;			/* Return spectral by default */

	return p;
}


/**********************************************************************
 * data conversion
 **********************************************************************/

unsigned char gt_whitebase_spm[] = {
	SPM_WBASE_PAP,
	SPM_WBASE_ABS
};

unsigned char gt_observer_spm[] = {
	SPM_OBSERVER_TWODEG,
	SPM_OBSERVER_TENDEG
};

unsigned char gt_illum_spm[] = {
	SPM_ILLUM_ILLUMINANTD50,
	SPM_ILLUM_ILLUMINANTD65
};

unsigned char gt_dstd_spm[] = {
	SPM_DSTD_ANSIA,
	SPM_DSTD_ANSIT
};

unsigned char gt_photomode_spm[] = {
	SPM_CONTROL_PHOTOMETRICABSOLUTE,
	SPM_CONTROL_PHOTOMETRICRELATIVE
};

unsigned char gt_filter_spm[] = {
	SPM_ACTUALFILTER_NOFILTER,
	SPM_ACTUALFILTER_D65FILTER,
	SPM_ACTUALFILTER_POLFILTER,
	SPM_ACTUALFILTER_UVCUTFILTER
};

unsigned char* gt_filter_txt[] = {
	"NoFilter",
	"D65Filter",
	"PolFilter",
	"UVCutFilter"
};

unsigned char gt_spectype_spm[] = {
	SPM_SPECT_RS,
	SPM_SPECT_DS
};

unsigned char gt_color_spm[] = {
	SPM_C_XYY,
	SPM_C_XYZ,
	SPM_C_LAB,
	SPM_C_LUV
};




/**********************************************************************
 * test/debug
 **********************************************************************/

#ifdef STANDALONE_TEST

#include <stdio.h>
#include <sys/time.h>
#include <GL/glut.h>

typedef enum {MEASURE_ONE, MEASURE_ALL} testtype;
testtype mytest = MEASURE_ONE;

struct _delay {
	int t;
	char* name;
};
struct _delay delay[] = {
	{0, "stop"},
	{5, "5s"},
	{10, "10s"},
	{15, "15s"},
	{30, "30s"},
	{60, "1m"},
	{5*60, "5m"},
	{30*60, "30m"},
	{60*60, "1h"},
	{6*60*60, "6h"},
	{24*60*60, "24h"}
};
int ndelay = sizeof(delay)/sizeof(struct _delay);
int mydelay =0;

struct _color {
	int r,g,b;
	char *name;
};
struct _color color[] = {
	{255,0,0,     "red"},
	{0,255,0,     "green"},
	{0,0,255,     "blue"},
	{0,0,0,       "gray0"},
	{16,16,16,    "gray16"},
	{32,32,32,    "gray32"},
	{48,48,48,    "gray48"},
	{64,64,64,    "gray64"},
	{80,80,80,    "gray80"},
	{96,96,96,    "gray96"},
	{112,112,112, "gray112"},
	{128,128,128, "gray128"},
	{144,144,144, "gray144"},
	{160,160,160, "gray160"},
	{176,176,176, "gray176"},
	{192,192,192, "gray192"},
	{208,208,208, "gray208"},
	{224,224,224, "gray224"},
	{240,240,240, "gray240"},
	{255,255,255, "gray255"}
};
int ncolor = sizeof(color)/sizeof(struct _color);
int mycolor = 11;

inst *p;
FILE *fp;
void logprintf(char *fmt, ...);
struct timezone tz = {0,0};
struct timeval begin, base, now;
#define ELAPSED(x,y) (x.tv_sec - y.tv_sec)

void
reshape(int w, int h)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glutPostRedisplay();
}

void
display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
}

void
measure(void)
{
	int i;
	float c[3];
	ipatch pat;
	gretag *g = (gretag*) p;

	begin = now;

	switch (mytest) {
	case MEASURE_ONE:
		gettimeofday(&now, &tz);
		g->measurement(g);
		g->query_color(g, GT_COLOR_XYY, c);
		logprintf("time %d %s xyY %f %f %f\n", ELAPSED(now,base), color[mycolor].name, c[0], c[1], c[2]);

		p->read_sample(p, "name", &pat);
		logprintf("time %d %s XYZ %lf %lf %lf xyY %lf %lf %lf\n", ELAPSED(now,base), color[mycolor].name,
				  pat.aXYZ[0],pat.aXYZ[1],pat.aXYZ[2],
				  pat.aXYZ[0] / (pat.aXYZ[0]+pat.aXYZ[1]+pat.aXYZ[2]),
				  pat.aXYZ[1] / (pat.aXYZ[0]+pat.aXYZ[1]+pat.aXYZ[2]),
				  pat.aXYZ[1]);

		break;
	case MEASURE_ALL:
		for (i=0;i<ncolor;i++) {
			glClearColor(color[i].r/255., color[i].g/255.,color[i].b/255.,0);
			glClear(GL_COLOR_BUFFER_BIT);
			glutSwapBuffers();

			gettimeofday(&now, &tz);
			g->measurement(g);
			g->query_color(g, GT_COLOR_XYY, c);
			logprintf("time %d %s xyY %f %f %f\n", ELAPSED(now,base), color[i].name, c[0], c[1], c[2]);
		}
		break;
	}
}

void
idle(void)
{
	gettimeofday(&now,&tz);
	if (mydelay && (ELAPSED(now,begin) >= delay[mydelay].t))
		measure();
}

void
key( unsigned char key, int x, int y )
{
   	switch (key) {
	case ' ':
		measure();
		break;
	case 'a':
		if (mytest == MEASURE_ONE) {
			mytest = MEASURE_ALL;
			logprintf("measure all patches\n");
		} else {
			mytest = MEASURE_ONE;
			logprintf("measure one patch only\n");
		}
		break;
	case 'c':
		mydelay = 0;
		logprintf("recalibrating\n");
		p->calibrate(p, 0);
		break;
	case 27:
		exit(1);
	}
}

void
colmenu(int i)
{
	mycolor = i;
	glClearColor(color[i].r/255., color[i].g/255.,color[i].b/255.,0);
	glutPostRedisplay();
}

void
delaymenu(int i)
{
	mydelay = i;
	gettimeofday(&begin,&tz);
	base = begin;
	logprintf("delay set to %s\n", delay[mydelay].name);
	glutPostRedisplay();
}

int
main(int argc, char **argv)
{
	inst *p;
	int i;

	if (argc == 2) fp = fopen(argv[1],"w+");
	else fp = fopen("spectro.log", "w+");
	if (!fp) perror("can't open logfile");

	/* scan serial ports for gretag device */
	p = new_inst(instSpectrolino);
	for (i=1; i<5; i++)
		if (p->init_coms(p,i,baud_nc,1) == inst_ok) break;
	if (i == 5) return 0;

	p->init_inst(p);
	p->set_mode(p,inst_mode_emis_disp);

	/* setup window stuff */
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
    glutCreateWindow("gretag");
    glutReshapeWindow(512,512);

	glutCreateMenu(colmenu);
	for (i=0;i<ncolor;i++) {
		glutAddMenuEntry(color[i].name,i);
	}
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	glutCreateMenu(delaymenu);
	for (i=0; i<ndelay; i++)
		glutAddMenuEntry(delay[i].name, i);
	glutAttachMenu(GLUT_MIDDLE_BUTTON);

    glutReshapeFunc (reshape);
    glutKeyboardFunc (key);
    glutDisplayFunc(display);
    glutIdleFunc(idle);

	gettimeofday(&begin,&tz);
	base=begin;

	/* init color */
	glClearColor(color[i].r/255., color[i].g/255.,color[i].b/255.,0);
    glutMainLoop();
    return 0;
}


/******************************************************************/
/* Basic printf type error() and warning() routines */

void
logprintf(char *fmt, ...)
{
	va_list args;

	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fflush(stderr);
	va_start(args, fmt);
	vfprintf(fp, fmt, args);
	va_end(args);
	fflush(fp);
}


void
error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"gt_gretag: Error - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	fflush(stdout);
	exit (-1);
}

#endif	/* STANDALONE_TEST */

