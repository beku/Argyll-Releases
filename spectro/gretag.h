#ifndef GRETAG_H
#define GRETAG_H

/* GretagMacbeth Spectrolino command protocol
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

#include "serio.h"
#include "stdarg.h"
#include "inst.h"
#include "spm.h"

typedef int bool;
#ifndef true
#define true 1
#endif
#ifndef false
#define false 0
#endif
#ifndef NULL
#define NULL 0
#endif


/* the following typedefs support only a subset of the entire
 * SPM protocol, i.e. just what's needed for our purposes.
 * we provide macros to convert these enums to their SPM
 * equivalents.  extending this ought to be trivial.
 */

typedef enum  {
	GT_WHITEBASE_PAPER,
	GT_WHITEBASE_ABSOLUTE
} gt_whitebase;

extern unsigned char gt_whitebase_spm[];
#define SPM_WBASE(x) gt_whitebase_spm[x]


typedef enum {
	GT_OBSERVER_2,
	GT_OBSERVER_10
} gt_observer;

extern unsigned char gt_observer_spm[];
#define SPM_OBSERVER(x) gt_observer_spm[x]


typedef enum {
	GT_ILLUM_D50,
	GT_ILLUM_D65
} gt_illum;

extern unsigned char gt_illum_spm[];
#define SPM_ILLUM(x) gt_illum_spm[x]


typedef enum {
	GT_DSTD_ANSIA,
	GT_DSTD_ANSIT
} gt_dstd;

extern unsigned char gt_dstd_spm[];
#define SPM_DSTD(x) gt_dstd_spm[x]


typedef enum {
	GT_PHOTOMODE_ABSOLUTE,
	GT_PHOTOMODE_RELATIVE
} gt_photomode;

extern unsigned char gt_photomode_spm[];
#define SPM_PHOTOMODE(x) gt_photomode_spm[x]


typedef enum {
	GT_FILTER_NONE,
	GT_FILTER_D65,
	GT_FILTER_POL,
	GT_FILTER_UVCUT
} gt_filter;

extern unsigned char gt_filter_spm[];
#define SPM_FILTER(x) gt_filter_spm[x]

extern unsigned char* gt_filter_txt[];
#define TXT_FILTER(x) gt_filter_txt[x]

typedef enum {
	GT_SPECTYPE_REMISSION,
	GT_SPECTYPE_DENSITY
} gt_spectype;

extern unsigned char gt_spectype_spm[];
#define SPM_SPECTYPE(x) gt_spectype_spm[x]


typedef enum {
	GT_COLOR_XYY,
	GT_COLOR_XYZ,
	GT_COLOR_LAB,
	GT_COLOR_LUV
} gt_color;

extern unsigned char gt_color_spm[];
#define SPM_COLOR(x) gt_color_spm[x]


/* Gretag Spectrolino/Spectroscan communication object */

struct _gretag {
/* base instrument class ************************************************ */
	INST_OBJ_BASE

/* gretag private data ************************************************** */

	bool         initialized;

	/* probably a bad idea caching device state like this */
	inst_mode    mode;			/* Currently instrument mode */
	inst_mode    lastmode;		/* Last requested mode */
	gt_whitebase whitebase;
	gt_observer  observer;
	gt_illum     illum;
	gt_dstd      dstd;
	gt_photomode photomode;
	float        photoref;
	gt_filter    filter;

	int          mcounter;
	bool         need_ref_measurement;
	bool		 read_returns_spectrum;

	int          unexrep;		/* last unexpected reply type */
	int          unexarg;		/* last unexpected reply first argument value if error */

/* gretag specific methods ************************************************ */

	/* Perform measurement */
	inst_code (*measurement)(struct _gretag *p);
	inst_code (*query_color)(struct _gretag *p, gt_color color, float *cvec);
	inst_code (*query_spectrum)(struct _gretag *p, gt_spectype spectype, float *svec);

	/* Manipulate device state */
	inst_code (*set_whitebase)(struct _gretag *p, gt_whitebase wbase);
	gt_whitebase (*get_whitebase)(struct _gretag *p);
	inst_code (*set_observer)(struct _gretag *p, gt_observer observer);
	gt_observer (*get_observer)(struct _gretag *p);
	inst_code (*set_illum)(struct _gretag *p, gt_illum illum);
	gt_illum (*get_illum)(struct _gretag *p);
	inst_code (*set_dstd)(struct _gretag *p, gt_dstd dstd);
	gt_dstd (*get_dstd)(struct _gretag *p);
	inst_code (*set_photomode)(struct _gretag *p, gt_photomode photomode);
	gt_photomode (*get_photomode)(struct _gretag *p);
	inst_code (*set_photoref)(struct _gretag *p, float photoref);
	float (*get_photoref)(struct _gretag *p);
	inst_code (*set_filter)(struct _gretag *p, gt_filter filter);
	gt_filter (*get_filter)(struct _gretag *p);

	/* Manage callbacks */
	inst_code (*set_logwarn)(struct _gretag *p, int (*f)(char*));
	inst_code (*set_logerr) (struct _gretag *p, int (*f)(char*));
	inst_code (*set_confirm)(struct _gretag *p, int (*f)(char*, char*, char*, char*));
	
/* gretag private methods *********************************************** */

	inst_code (*incr_mcounter)(struct _gretag *p);

	inst_code (*spm_command)(struct _gretag *p, int cmd, ...);

	int (*unexpected_reply)(struct _gretag *p);
	int (*unexpected_argument)(struct _gretag *p);

	/* callbacks for logging and user feedback */
	int (*logwarn)(char *msg);
	int (*logerr)(char *msg);
	int (*confirm)(char *msg, char *b1, char *b2, char *b3);

	}; typedef struct _gretag gretag;

/* Constructor */
extern gretag *new_gretag(void);


#endif
