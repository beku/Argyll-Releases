#ifndef XSPECT_H
#define XSPECT_H

/* 
 * International Color Consortium color transform expanded support
 *
 * Author:  Graeme W. Gill
 * Date:    21/6/01
 * Version: 1.00
 *
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the LICENCE.TXT file for licencing details.
 *
 */

/*
 * This class supports converting spectral samples
 * into CIE XYZ or D50 Lab tristimulous values.
 */

/*
 * TTBD:
 *
 */

#include "icc.h"		/* icclib ICC definitions */ 

/* ------------------------------------------------------------------------------ */

/* Structure for conveying spectral information */

#define XSPECT_MAX_BANDS 500		/* Enought for 1nm */

typedef struct {
	int    spec_n;					/* Number of spectral bands, 0 if not valid */
	double spec_wl_short;			/* First reading wavelength in nm (shortest) */
	double spec_wl_long;			/* Last reading wavelength in nm (longest) */
	double norm;					/* Normalising scale value */
	double spec[XSPECT_MAX_BANDS];	/* Spectral value, shortest to longest */
} xspect;

/* Spectrum utility functions */
int write_xspect(char *fname, xspect *s);
int read_xspect(xspect *sp, char *fname);

/* ------------------------------------------------------------------------------ */
/* Class for converting between spectral and CIE */

/* We build in some useful spectra */

/* Type of illumination */
typedef enum {
    icxIT_default	= 0,	/* Default illuminant (usually D50) */
    icxIT_none		= 1,	/* No illuminant - self luminous spectrum */
    icxIT_custom	= 2,	/* Custom illuminant spectrum */
    icxIT_A			= 3,	/* Standard Illuminant A */
    icxIT_D50		= 4,	/* Daylight 5000K */
    icxIT_D65		= 5,	/* Daylight 6500K */
    icxIT_F5		= 6,	/* Fluorescent, Standard, 6350K, CRI 72 */
    icxIT_F8		= 7,	/* Fluorescent, Broad Band 5000K, CRI 95 */
    icxIT_F10		= 8,	/* Fluorescent Narrow Band 5000K, CRI 81 */
	icxIT_Spectrocam = 9	/* Spectrocam Xenon Lamp */
} icxIllumeType;

/* Return a pointer to a standard illuminant spectrum */
/* return NULL if not matched */
xspect *standardIlluminant(icxIllumeType ilType);

/* Type of observer */
typedef enum {
    icxOT_default			= 0,	/* Default observer (usually CIE_1931_2) */
    icxOT_custom			= 1,	/* Custom observer type weighting */
    icxOT_CIE_1931_2		= 2,	/* Standard CIE 1931 2 degree */
    icxOT_Stiles_Burch_2	= 3,	/* Stiles & Burch 1955 2 degree */
    icxOT_Judd_Voss_2		= 4,	/* Judd & Voss 1978 2 degree */
    icxOT_CIE_1964_10		= 5,	/* Standard CIE 1964 10 degree */
    icxOT_CIE_1964_10c		= 6,	/* Standard CIE 1964 10 degree, 2 degree compatible */
    icxOT_Shaw_Fairchild_2	= 7		/* Shaw & Fairchild 1997 2 degree */
} icxObserverType;

/* The conversion object */
struct _xsp2cie {
	/* Private: */
	xspect illuminant;
	xspect observer[3];
	int doLab;					/* Return D50 Lab result */

	/* FWA compensation */
	double bw;		/* Intergration bandwidth */
	xspect emit;	/* Estimated FWA emmission spectrum */
	xspect media;	/* Estimated base media (ie. minus FWA) */
	xspect instr;	/* Normalised instrument illuminant spectrum */
	xspect illum;	/* Normalised target illuminant spectrum */
	double Sm;		/* FWA Stimulation level for emit contribution */
	double FWAc;	/* FWA content (informational) */

	/* Public: */
	void (*del)(struct _xsp2cie *p);

	/* Convert (and possibly fwa correct) reflectance spectrum */
	void (*convert) (struct _xsp2cie *p,	/* this */
	                 double *out,			/* Return XYZ or D50 Lab value */
	                 xspect *in				/* Spectrum to be converted, normalised by norm */
	                );

	/* Convert and also return (possibly corrected) reflectance spectrum */
	/* Spectrum will be same wlength range and readings as input spectrum */
	void (*sconvert) (struct _xsp2cie *p,	/* this */
	                 xspect *sout,			/* Return corrected reflectance spectrum */
	                 double *out,			/* Return XYZ or D50 Lab value */
	                 xspect *in				/* Spectrum to be converted, normalised by norm */
	                );

	/* Set Fluorescent Whitening Agent compensation */
	/* return NZ if error */
	int (*set_fwa) (struct _xsp2cie *p,	/* this */
					xspect *inst,		/* Spectrum of instrument illuminant */
	                xspect *white		/* Spectrum of plain media */
	                );

}; typedef struct _xsp2cie xsp2cie;

xsp2cie *new_xsp2cie(
	icxIllumeType ilType,			/* Illuminant */
	xspect        *custIllum,

	icxObserverType obType,			/* Observer */
	xspect        *custObserver[3],

	icColorSpaceSignature  rcs		/* Return color space, icSigXYZData or icSigLabData */
);

/* --------------------------- */
/* Density and other functions */

/* Given a reflectance or transmition spectral product, */
/* return status T CMY + V density values */
void xsp_Tdensity(double *out,			/* Return CMYV density */
                 xspect *in				/* Spectral product to be converted */
                );

/* Given a reflectance or transmission XYZ value, */
/* return approximate status T CMYV log10 density values */
void icx_XYZ2Tdens(
double *out,			/* Return aproximate CMYV log10 density */
double *in				/* Input XYZ values */
);

/* Given a reflectance or transmission XYZ value, */
/* return log10 XYZ density values */
void icx_XYZ2dens(
double *out,			/* Return log10 XYZ density */
double *in				/* Input XYZ values */
);

/* Given an XYZ value, */
/* return sRGB values */
void icx_XYZ2sRGB(
double *out,			/* Return sRGB value */
double *wp,				/* Input XYZ white point (may be NULL) */
double *in				/* Input XYZ values */
);


/* Given a daylight color temperature in degrees K, */
/* return the corresponding XYZ value */
void icx_DTEMP2XYZ(
double *out,			/* Return XYZ value with Y == 1 */
double in				/* Input temperature in degrees K */
);


#endif /* XSPECTFM_H */






































