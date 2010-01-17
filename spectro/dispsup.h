
#ifndef DISPSUP_H

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
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* A default calibration user interaction handler using the console. */
/* this will function only for non-display specific types of calibration. */
typedef struct {
	disppath *disp;			/* display to calibrate. */
	int blackbg;			/* NZ if whole screen should be filled with black */
	int override;			/* Override_redirect on X11 */
	double patsize;			/* Size of dispwin */
	double ho, vo;			/* Position of dispwin */
	dispwin *dw;			/* Display window if already open */
} disp_win_info;

inst_code inst_handle_calibrate(
	inst *p,
	inst_cal_type calt,		/* type of calibration to do. inst_calt_all for all */
	inst_cal_cond calc,		/* Current condition. inst_calc_none for not setup */
	disp_win_info *dwi		/* Information to be able to open a display test patch */
);

/* User requested calibration of the display instrument */
int disprd_calibration(
instType itype,		/* Instrument type */
int comport, 		/* COM port used */
flow_control fc,	/* Serial flow control */
int dtype,			/* Display type, 0 = unknown, 1 = CRT, 2 = LCD */
int proj,			/* NZ for projector mode */
int adaptive,		/* NZ for adaptive mode */
int nocal,			/* NZ to disable auto instrument calibration */
disppath *screen,	/* Screen to calibrate. */
int blackbg,		/* NZ if whole screen should be filled with black */
int override,		/* Override_redirect on X11 */
double patsize,		/* Size of dispwin */
double ho, double vo,	/* Position of dispwin */
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

	xspect sp;			/* Spectrum. sp.spec_n > 0 if valid */

	double duration;	/* Total duration in seconds (flash measurement) */
} col;


/* Maximum number of entries to setup for calibration */
#define MAX_CAL_ENT 4096

/* Display reading context */
struct _disprd {

/* private: */
	int verb;			/* Verbosity flag */
	FILE *df;			/* Verbose output */
	int fake;			/* Fake display/instrument flag */
	char *fake_name;	/* Fake profile name */
	icmFile *fake_fp;
	icc *fake_icc;		/* NZ if ICC profile is being used for fake */
	double cal[3][MAX_CAL_ENT];	/* Calibration being worked through (cal[0][0] < 0.0 if not used */
	int ncal;			/* Number of entries used in cal[] */
	icmLuBase *fake_lu;
	char *mcallout;		/* fake instrument shell callout */
	int debug;			/* Debug flag */
	int comport; 		/* COM port used */
	baud_rate br;
	flow_control fc;
	inst *it;			/* Instrument */
	instType itype;		/* Instrument type */
	int dtype;			/* Display type, 0 = unknown, 1 = CRT, 2 = LCD */
	int proj;			/* NZ for projector mode */
	int adaptive;		/* NZ for adaptive mode */
	int spectral;		/* Spectral values requested/used */
	int nocal;			/* No automatic instrument calibration */
	int highres;		/* Use high res mode if available */
	dispwin *dw;		/* Window */
	ramdac *or;			/* Original ramdac if we set one */

/* public: */

	/* Destroy ourselves */
	void (*del)(struct _disprd *p);

	/* Take a series of readings from the display */
	/* return nz on fail/abort */
	/* 1 = user aborted */
	/* 2 = instrument access failed */
	/* 3 = window access failed */ 
	/* 4 = user hit terminate key */
	/* 5 = system error */
	int (*read)(struct _disprd *p,
		col *cols,		/* Array of patch colors to be tested */
		int npat, 		/* Number of patches to be tested */
		int spat,		/* Start patch index for "verb", 0 if not used */
		int tpat,		/* Total patch index for "verb", 0 if not used */
		int acr,		/* If nz, do automatic final carriage return */
		int tc			/* If nz, termination key */
	);

	/* Take an ambient reading if the instrument has the capability. */
	/* return nz on fail/abort */
	/* 1 = user aborted */
	/* 2 = instrument access failed */
	/* 3 = no ambient capability */ 
	/* 4 = user hit terminate key */
	/* 5 = system error */
	int (*ambient)(struct _disprd *p,
		double *ambient,	/* return ambient in cd/m^2 */
		int tc				/* If nz, termination key */
	);

}; typedef struct _disprd disprd;

/* Create a display reading object. */
/* Return NULL if error */
/* Set *errc to code: */
/* 1 = user aborted */
/* 2 = instrument access failed */
/* 3 = window access failed */
/* 4 = user hit terminate key */
/* 5 = system error */
/* Use disprd_err() to interpret errc */
disprd *new_disprd(
int *errc,			/* Error code. May be NULL */ 
instType itype,		/* Instrument type */
int comport, 		/* COM port used, -99 for fake display */
flow_control fc,	/* Serial flow control */
int dtype,			/* Display type, 0 = unknown, 1 = CRT, 2 = LCD */
int proj,			/* NZ for projector mode */
int adaptive,		/* NZ for adaptive mode */
int nocal,			/* No automatic instrument calibration */
int highres,		/* Use high res mode if available */
int donat,			/* Use ramdac for native output, else run through current or set ramdac */
double cal[3][MAX_CAL_ENT],	/* Calibration set/return (cal[0][0] < 0.0 if can't/not to be used) */
int ncal,			/* number of entries use in cal */
disppath *screen,	/* Screen to calibrate. */
int blackbg,		/* NZ if whole screen should be filled with black */
int override,		/* Override_redirect on X11 */
char *ccallout,		/* Shell callout on set color */
char *mcallout,		/* Shell callout on measure color (forced fake) */
double patsize,		/* Size of dispwin */
double ho,			/* Horizontal offset */
double vo,			/* Vertical offset */
int spectral,		/* Generate spectral info flag */
int verb,			/* Verbosity flag */
FILE *df,			/* Verbose output - NULL = stdout */
int debug,			/* Debug flag */
char *fake_name		/* Name of profile to use as a fake device */
);

/* Return a string describing the error code */
char * disprd_err(int en);

#define DISPSUP_H
#endif /* DISPSUP_H */

