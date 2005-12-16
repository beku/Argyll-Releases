
#ifndef DISPSUP_H

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

/* User requested calibration of the display instrument */
int disprd_calibration(
instType itype,		/* Instrument type */
int comport, 		/* COM port used */
int override,		/* Override_redirect on X11 */
double patsize,		/* Size of dispwin */
int verb,			/* Verbosity flag */
int debug			/* Debug flag */
);


/* A color structure to return values with. */
/* This can hold all representations simultaniously */
typedef struct {
	double r,g,b;
	char *id;			/* Id string */

	int    XYZ_v;
	double XYZ[3];		/* Colorimeter readings */

	int    aXYZ_v;
	double aXYZ[3];		/* Absolute colorimeter readings */

	int    spec_n;					/* Number of spectral bands, 0 if not valid */
	double spec_wl_short;			/* First reading wavelength in nm (shortest) */
	double spec_wl_long;			/* Last reading wavelength in nm (longest) */
	double spec[INSTR_MAX_BANDS];	/* Spectral reflectance %, shortest to longest */
} col;

/* Display reading context */
struct _disprd {

/* private: */
	int verb;			/* Verbosity flag */
	FILE *df;			/* Verbose output */
	int debug;			/* Debug flag */
	int comport; 		/* COM port used */
	baud_rate br;
	inst *it;			/* Instrument */
	dispwin *dw;		/* Window */
	ramdac *or;			/* Original ramdac if we set one */
	int spectral;		/* Spectral values requested */

/* public: */

	/* Destroy ourselves */
	void (*del)(struct _disprd *p);

	/* Take a series of readings from the display */
	/* return nz on fail/abort */
	int (*read)(struct _disprd *p,
		col *cols,		/* Array of patch colors to be tested */
		int npat, 		/* Number of patches to be tested */
		int spat,		/* Start patch index for "verb", 0 if not used */
		int tpat		/* Total patch index for "verb", 0 if not used */
	);
}; typedef struct _disprd disprd;

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
);

#define DISPSUP_H
#endif /* DISPSUP_H */

