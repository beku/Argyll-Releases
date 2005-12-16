#ifndef INST_H

/* 
 * Argyll Color Correction System
 *
 * Abstract base class for common color instrument interface
 * and other common instrument stuff.
 *
 * Author: Graeme W. Gill
 * Date:   15/3/2001
 *
 * Copyright 2001 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 *
 */

#include "insttypes.h"

/* ------------------------------------------------- */
/* Structure for holding an instrument patch reading */

#define INSTR_MAX_BANDS 75

struct _ipatch {
	int XYZ_v;			/* XYZ valid */
	double XYZ[3];		/* XYZ values, 0.0 .. 100.0 */

	int aXYZ_v;			/* Absolute XYZ valid */
	double aXYZ[3];		/* XYZ values in cd/m^2 */

	int Lab_v;			/* Lab valid */
	double Lab[3];		/* Lab value */

	int    spec_n;					/* Number of spectral bands, 0 if not valid */
	double spec_wl_short;			/* First reading wavelength in nm (shortest) */
	double spec_wl_long;			/* Last reading wavelength in nm (longest) */
	double spec[INSTR_MAX_BANDS];	/* Spectral reflectance %, shortest to longest */

	}; typedef struct _ipatch ipatch;

/* ---------------------------------------- */
/* Instrument interface abstract base class */

/* Instrument capabilities */
typedef enum {
	inst_unknown           =  0x0000,	/* Capabilities can't be determined */
	inst_ref_spot          =  0x0001,	/* Capable of reflection spot measurement */
	inst_ref_strip         =  0x0002,	/* Capable of reflection strip measurement */
	inst_ref_xy            =  0x0004,	/* Capable of reflection X-Y measurement */
	inst_reflection        =  0x0007,	/* Capable of general reflection measurements */

	inst_trans_spot        =  0x0010,	/* Capable of transmission spot measurement */
	inst_trans_strip       =  0x0020,	/* Capable of transmission strip measurement */
	inst_trans_xy          =  0x0040,	/* Capable of transmission X-Y measurement */
	inst_transmission      =  0x0070,	/* Capable of general transmission measurements */

	inst_emis_spot         =  0x0100,	/* Capable of general emission spot measurement */
	inst_emis_disp         =  0x0200,	/* Capable of display emission measurement */
	inst_emis_illum        =  0x0400,	/* Capable of illuminant emission measurement */
	inst_emission          =  0x0700,	/* Capable of general emission measurements */

	inst_colorimeter       =  0x1000,	/* Colorimetric capability */
	inst_spectral          =  0x2000,	/* Spectral capability */
	inst_calib             =  0x4000,	/* Uses computer conrolled calibration */

	inst_xy_holdrel        = 0x10000,	/* Needs paper hold/release between each sheet */
	inst_xy_locate         = 0x20000,	/* Needs user to locate patch locations */
	inst_xy_position       = 0x40000	/* Can be positioned at a given location */

} inst_capability;

/* Instrument modes and sub-modes */
/* We assume that we only want to be in one measurement mode at a time */
typedef enum {
	inst_mode_unknown           = 0x0000,	/* Mode not specified */

	inst_mode_reflection        = 0x0001,	/* General reflection mode */
	inst_mode_transmission      = 0x0002,	/* General transmission mode */
	inst_mode_emission          = 0x0003,	/* General emission mode */
	inst_mode_illum_mask        = 0x000f,	/* Mask of sample illumination mode */

	inst_mode_spot              = 0x0010,	/* General spot measurement mode */
	inst_mode_strip             = 0x0020,	/* General strip measurement mode */
	inst_mode_xy                = 0x0030,	/* General X-Y measurement mode */
	inst_mode_disp              = 0x0040,	/* General Display measurement mode */
	inst_mode_illum             = 0x0050,	/* General illuminant measurement mode */
	inst_mode_sub_mask          = 0x00f0,	/* Mask of sub-mode */

	inst_mode_ref_spot          = 0x0011,	/* Reflection spot measurement mode */
	inst_mode_ref_strip         = 0x0021,	/* Reflection strip measurement mode */
	inst_mode_ref_xy            = 0x0031,	/* Reflection X-Y measurement mode */
	inst_mode_trans_spot        = 0x0012,	/* Transmission spot measurement mode */
	inst_mode_trans_strip       = 0x0022,	/* Transmission strip measurement mode */
	inst_mode_trans_xy          = 0x0032,	/* Transmission X-Y measurement mode */
	inst_mode_emis_spot         = 0x0013,	/* Spot emission measurement mode */
	inst_mode_emis_disp         = 0x0043,	/* Display emission measurement mode */
	inst_mode_emis_illum        = 0x0053,	/* Illuminant emission measurement mode */

	inst_mode_measuremet_mask   = 0x00ff,	/* Mask of exclusive measurement modes */

	inst_mode_colorimeter       = 0x1000,	/* Colorimetric mode */
	inst_mode_spectral          = 0x2000	/* Spectral mode */
} inst_mode;

/* Instrument option-modes */
typedef enum {
	inst_opt_unknown           = 0x0000,	/* Option not specified */

	inst_opt_noautocalib        = 0x0001,	/* Disable auto calibration */
	inst_opt_autocalib          = 0x0002	/* Re-enable auto calibration */
} inst_opt_mode;

/* Abstract return codes in ms byte. */
/* Machine dependant codes in ls byte. */
typedef enum {
	inst_ok                = 0x0000,
	inst_notify            = 0x0100,	/* A Notification */
	inst_warning           = 0x0200,	/* A Warning */
	inst_internal_error    = 0x0300,	/* Internal software error */
	inst_coms_fail         = 0x0400,	/* Communication failure */
	inst_unknown_model     = 0x0500,	/* Not the expected instrument */
	inst_protocol_error    = 0x0600, 	/* Read or Write protocol error */
	inst_user_abort        = 0x0700,	/* User hit escape */
	inst_misread           = 0x0800,	/* Strip misread */
	inst_needs_cal         = 0x0900,	/* Instrument needs calibration */
	inst_unsupported       = 0x0A00,	/* Unsupported function */
	inst_unexpected_reply  = 0x0B00,	/* Unexpected Reply */
	inst_wrong_config      = 0x0C00,    /* Configuration is wrong */
	inst_hardware_fail     = 0x0D00,    /* Hardware failure */
	inst_other_error       = 0x0E00,	/* Some other error */
	inst_mask              = 0xff00		/* Mask value */
} inst_code;

/* Type of user interaction */
typedef enum {
	inst_verb              = 0x00,	/* verbose message */
	inst_question          = 0x01	/* A question requiring answers */
} inst_uiact;

/* Color instrument interface base object */
#define INST_OBJ_BASE															\
																				\
	int debug;		/* debug flag */											\
	int verb;		/* Verbosity level */                                       \
	instType  itype;															\
	serio *sio;		/* Serial coms object */									\
	int gotcoms;	/* Coms established flag */                                 \
	int inited;		/* Initialised flag */                                      \
																				\
	/* Establish communications at the indicated baud rate. */					\
	/* Timout in to seconds, and return non-zero error code */					\
	inst_code (*init_coms)(														\
        struct _inst *p,														\
        int port,			/* Serial port number */							\
        baud_rate br,		/* Baud rate */										\
        double tout);		/* Timeout */										\
																				\
	/* Initialise or re-initialise the INST */									\
	/* return non-zero on an error, with inst error code */						\
	inst_code (*init_inst)(  													\
        struct _inst *p);														\
																				\
	/* Return the instrument capabilities */									\
	inst_capability (*capabilities)(struct _inst *p);							\
																				\
	/* Callback to allow instrument to interact with user. */                   \
	/* Default will be stdio. */                                                \
	/* Returns users choice 1..nch, 0 if no choices */                          \
	int (*user)(struct _inst *p,                                                \
	           inst_uiact it,		/* Interaction type */                      \
	           int nm,				/* Number of messages */                    \
	           int nch,				/* Number of choices to accept */           \
	           ...);                /* sequence of pointers to strings */       \
                                                                                \
    /* Set the device measurement mode */                                       \
    inst_code (*set_mode)(														\
        struct _inst *p,														\
        inst_mode m);		/* Requested mode */								\
																				\
    /* Set or reset an option */												\
    inst_code (*set_opt_mode)(													\
        struct _inst *p,														\
        inst_opt_mode m);	/* Requested mode */								\
																				\
	/* For an xy instrument, release the paper */								\
	/* (if cap has inst_xy_holdrel) */											\
	/* Return the inst error code */											\
	inst_code (*xy_sheet_release)(												\
		struct _inst *p);														\
																				\
	/* For an xy instrument, hold the paper */									\
	/* (if cap has inst_xy_holdrel) */											\
	/* Return the inst error code */											\
	inst_code (*xy_sheet_hold)(													\
		struct _inst *p);														\
																				\
	/* For an xy instrument, allow the user to locate a point */				\
	/* (if cap has inst_xy_locate) */											\
	/* Return the inst error code */											\
	inst_code (*xy_locate_start)(												\
		struct _inst *p);														\
																				\
	/* For an xy instrument, read back the location */							\
	/* (if cap has inst_xy_locate) */											\
	/* Return the inst error code */											\
	inst_code (*xy_get_location)(												\
		struct _inst *p,														\
		double *x, double *y);													\
																				\
	/* For an xy instrument, end allowing  the user to locate a point */		\
	/* (if cap has inst_xy_locate) */											\
	/* Return the inst error code */											\
	inst_code (*xy_locate_end)(													\
		struct _inst *p);														\
																				\
	/* For an xy instrument, move the measurement point */						\
	/* (if cap has inst_xy_position) */											\
	/* Return the inst error code */											\
	inst_code (*xy_position)(													\
		struct _inst *p,														\
		int measure,	/* nz for measure point, z for locate point */	        \
		double x, double y);													\
																				\
	/* For an xy instrument, try and clear the table after an abort */			\
	/* Return the inst error code */											\
	inst_code (*xy_clear)(														\
		struct _inst *p);														\
																				\
	/* Read a sheet full of patches using xy mode */							\
	/* Return the inst error code */											\
	inst_code (*read_xy)(														\
		struct _inst *p,														\
		int pis,			/* Passes in strip (letters in sheet) */			\
		int sip,			/* Steps in pass (numbers in sheet) */				\
		int npatch,			/* Total patches in strip (skip in last pass) */	\
		char *pname,		/* Starting pass name (' A' to 'ZZ') */             \
		char *sname,		/* Starting step name (' 1' to '99') */             \
		double ox, double oy,	/* Origin of first patch */						\
		double ax, double ay,	/* pass increment */							\
		double aax, double aay,	/* pass offset for odd patches */				\
		double px, double py,	/* step/patch increment */						\
		ipatch *vals);		/* Pointer to array of values */					\
																				\
	/* Read a set of strips (applicable to strip reader) */						\
	/* Return the inst error code */											\
	inst_code (*read_strip)(													\
		struct _inst *p,														\
		char *name,			/* Strip name (up to 7 chars) */					\
		int npatch,			/* Number of patches in each pass */				\
		char *pname,		/* Pass name (3 chars) */							\
		int sguide,			/* Guide number (decrements by 5) */				\
		double pwid,		/* Patch width in mm (For DTP41) */					\
		double gwid,		/* Gap width in mm  (For DTP41) */					\
		double twid,		/* Trailer width in mm  (For DTP41T) */				\
		ipatch *vals);		/* Pointer to array of values */					\
																				\
	/* Read a single sample (applicable to spot instruments */					\
	/* Return the inst error code */											\
	inst_code (*read_sample)(													\
		struct _inst *p,														\
		char *name,			/* Patch identifier (up to 7 chars) */				\
		ipatch *val);		/* Pointer to value to be returned */				\
																				\
	/* Determine if a calibration is needed. Returns inst_ok if not, */         \
	/* inst_unsupported if it is unknown, or inst_needs_cal if needs calibration */ \
	inst_code (*needs_calibration)(												\
		struct _inst *p);														\
																				\
	/* Request an instrument calibration. */									\
	/* (Some instruments will do so anyway when mode is change, or */           \
	/*  if it is needed. Use this if it is voluntary) */                        \
	inst_code (*calibrate)(														\
		struct _inst *p,														\
		int cal_no);		/* Calibration type, 0 - n */						\
																				\
	/* Generic inst error codes interpretation */								\
	char * (*inst_interp_error)(struct _inst *p, inst_code ec);					\
																				\
	/* Instrument specific error codes interpretation */						\
	char * (*interp_error)(struct _inst *p, int ec);							\
																				\
	/* Return the last serial i/o error code */									\
	int (*last_sioerr)(struct _inst *p);										\
																				\
	/* Destroy ourselves */														\
	void (*del)(struct _inst *p);												\

/* The base object type */
struct _inst {
	INST_OBJ_BASE
	}; typedef struct _inst inst;

/* Constructor */
extern inst *new_inst(instType itype);

#define INST_H
#endif /* INST_H */
