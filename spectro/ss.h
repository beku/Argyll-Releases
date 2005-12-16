
#ifndef SS_H

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

#undef EMSST		/* Debug - emulate a SpectroScanT with a SpectroScan */

#include "inst.h"
#include "ss_imp.h"

#define SS_MAX_WR_SIZE 1000	/* Assumed maximum normal message query size */
#define SS_MAX_RD_SIZE 1000	/* Assumed maximum normal messagle answer size */

/* Gretag Spectrolino/Spectroscan communication object */
struct _ss {
	/* **** base instrument class **** */
	INST_OBJ_BASE

	/* *** Spectroscan/lino private data **** */
	inst_capability	cap;		/* Instrument capability */
	inst_mode	nextmode;		/* Next requested mode */
	inst_mode	mode;			/* Currently instrument mode */

	/* Desired measurement configuration */
	ss_aft     filt;			/* Filter type (None/UV/D65 etc.) */
	ss_dst     dstd;			/* Density standard (ANSI A/ANSI T/DIN etc.) */
	ss_ilt     illum;			/* Illuminant type (A/C/D50 etc.) */
	ss_ot      obsv;			/* Observer type (2deg/10deg) */
	ss_wbt     wbase;			/* White base type (Paper/Abs>) */
	ss_ctt     phmode;			/* Photometric mode (Absolute/Relative) */
	double     phref;			/* Photometric reference (cd/m^2) */

	int		calcount;			/* Calibration needed counter */
	int		need_cal;			/* Calibration needed flag */
	int     noutocalib;			/* Don't auto calibrate */

	/* Emulated manual transmission mode support */
	double tref[36];			/* White reference spectrum */
	double cill[36];			/* Colorimetry illuminant */

#ifdef EMSST
	int tmode;					/* Transmission mode */
	ss_rt sbr;					/* Standby reference */
	double sbx, sby;			/* Standby location */
#endif

	/* Serialisation support: */

	/* Send buffer */
	char _sbuf[SS_MAX_WR_SIZE];	/* Buffer allocation */
	char *sbufe;				/* Pointer to last valid location in _sbuf[] */
	char *sbuf;					/* Next character output */

	/* Receive buffer */
	char _rbuf[SS_MAX_RD_SIZE];	/* Buffer allocation */
	char *rbufe;				/* Pointer to last valid location in _rbuf[] */
	char *rbuf;					/* Next character output */

	/* Accumulated device error */
	ss_et snerr;

	}; typedef struct _ss ss;

/* Constructor */
extern ss *new_ss(void);

#define SS_H
#endif /* SS_H */
