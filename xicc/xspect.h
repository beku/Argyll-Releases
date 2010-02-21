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
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
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

#define XSPECT_MAX_BANDS 601		/* Enought for 1nm from 300 to 900 */

typedef struct {
	int    spec_n;					/* Number of spectral bands, 0 if not valid */
	double spec_wl_short;			/* First reading wavelength in nm (shortest) */
	double spec_wl_long;			/* Last reading wavelength in nm (longest) */
	double norm;					/* Normalising scale value */
	double spec[XSPECT_MAX_BANDS];	/* Spectral value, shortest to longest */
} xspect;

/* Some helpful macro's: */

/* Given an index and the sampling ranges, compute the sample wavelegth */
#define XSPECT_WL(SHORT, LONG, N, IX) \
((SHORT) + (double)(IX) * ((LONG) - (SHORT))/((N)-1.0))

/* Given the address of an xspect and an index, compute the sample wavelegth */
#define XSPECT_XWL(PXSP, IX) \
(((PXSP)->spec_wl_short) + (double)(IX) * (((PXSP)->spec_wl_long) - ((PXSP)->spec_wl_short))/(((PXSP)->spec_n)-1.0))

/* Given a wavelength and the sampling ranges, compute the nearest index */
#define XSPECT_IX(SHORT, LONG, N, WL) \
(int)(((N)-1.0) * ((WL) - (SHORT))/((LONG) - (SHORT)) + 0.5)

/* Spectrum utility functions. Return NZ if error */
int write_xspect(char *fname, xspect *s);
int read_xspect(xspect *sp, char *fname);

/* Get interpolated value at wavelenth (not normalised) */
double value_xspect(xspect *sp, double wl);

/* Convert from one xspect type to another */
void xspect2xspect(xspect *dst, xspect *targ, xspect *src);
/* ------------------------------------------------------------------------------ */
/* Class for converting between spectral and CIE */

/* We build in some useful spectra */

/* Type of illumination */
typedef enum {
    icxIT_default	 = 0,	/* Default illuminant (usually D50) */
    icxIT_none		 = 1,	/* No illuminant - self luminous spectrum */
    icxIT_custom	 = 2,	/* Custom illuminant spectrum */
    icxIT_A			 = 3,	/* Standard Illuminant A */
	icxIT_C          = 4,	/* Standard Illuminant C */
    icxIT_D50		 = 5,	/* Daylight 5000K */
    icxIT_D65		 = 6,	/* Daylight 6500K */
    icxIT_F5		 = 7,	/* Fluorescent, Standard, 6350K, CRI 72 */
    icxIT_F8		 = 8,	/* Fluorescent, Broad Band 5000K, CRI 95 */
    icxIT_F10		 = 9,	/* Fluorescent Narrow Band 5000K, CRI 81 */
	icxIT_Spectrocam = 10,	/* Spectrocam Xenon Lamp */
    icxIT_Dtemp		 = 11,	/* Daylight at specified temperature */
    icxIT_Ptemp		 = 12	/* Planckian at specified temperature */
} icxIllumeType;

/* Fill in an xpsect with a standard illuminant spectrum */
/* return 0 on sucecss, nz if not matched */
int standardIlluminant(
xspect *sp,					/* Xspect to fill in */
icxIllumeType ilType,		/* Type of illuminant */
double temp);				/* Optional temperature in degrees kelvin, for Dtemp and Ptemp */


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

/* Fill in three xpsects with a standard observer weighting curves */
/* return 0 on sucecss, nz if not matched */
int standardObserver(
xspect *sp0,
xspect *sp1,
xspect *sp2,				/* Xspects to fill in */
icxObserverType obType);	/* Type of observer */


/* The conversion object */
struct _xsp2cie {
	/* Private: */
	xspect illuminant;
	int isemis;					/* nz if we are doing an emission conversion */
	xspect observer[3];
	int doLab;					/* Return D50 Lab result */

	/* FWA compensation */
	double bw;		/* Integration bandwidth */
	xspect emits;	/* Estimated FWA emmission spectrum */
	xspect media;	/* Estimated base media (ie. minus FWA) */
	xspect instr;	/* Normalised instrument illuminant spectrum */
	xspect illum;	/* Normalised target illuminant spectrum */
	double Sm;		/* FWA Stimulation level for emits contribution */
	double FWAc;	/* FWA content (informational) */

	/* Public: */
	void (*del)(struct _xsp2cie *p);

	/* Convert (and possibly fwa correct) reflectance spectrum */
	/* Note that XYZ is 0..1 range */
	void (*convert) (struct _xsp2cie *p,	/* this */
	                 double *out,			/* Return XYZ or D50 Lab value */
	                 xspect *in				/* Spectrum to be converted, normalised by norm */
	                );

	/* Convert and also return (possibly corrected) reflectance spectrum */
	/* Spectrum will be same wlength range and readings as input spectrum */
	/* Note that XYZ is 0..1 range */
	void (*sconvert) (struct _xsp2cie *p,	/* this */
	                 xspect *sout,			/* Return corrected refl. spectrum (may be NULL) */
	                 double *out,			/* Return XYZ or D50 Lab value (may be NULL) */
	                 xspect *in				/* Spectrum to be converted, normalised by norm */
	                );

	/* Set Media White value */
	/* return NZ if error */
	int (*set_mw) (struct _xsp2cie *p,	/* this */
	                xspect *white		/* Spectrum of plain media */
	                );

	/* Set Fluorescent Whitening Agent compensation */
	/* return NZ if error */
	int (*set_fwa) (struct _xsp2cie *p,	/* this */
					xspect *inst,		/* Spectrum of instrument illuminant */
	                xspect *white		/* Spectrum of plain media */
	                );

	/* Get Fluorescent Whitening Agent compensation information */
	/* return NZ if error */
	void (*get_fwa_info) (struct _xsp2cie *p,	/* this */
					double *FWAc		/* FWA content as a ratio. */
	                );

	/* Extract the colorant reflectance value from the media. Takes FWA */
	/* into account if set. Media white or FWA must be set. */
	/* return NZ if error */
	int (*extract) (struct _xsp2cie *p,	/* this */
	                 xspect *out,			/* Extracted colorant refl. spectrum */
	                 xspect *in				/* Spectrum to be converted, normalised by norm */
	                );


	/* Apply the colorant reflectance value from the media. Takes FWA */
	/* into account if set. DOESN'T convert to FWA target illumination! */
	/* FWA must be set. */
	int (*apply) (struct _xsp2cie *p,	/* this */
	                 xspect *out,			/* Applied refl. spectrum */
	                 xspect *in				/* Colorant reflectance to be applied */
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



/* Given an illuminant definition and an observer model, return */
/* the normalised XYZ value for that spectrum. */
/* Return 0 on sucess, 1 on error */
int icx_ill_sp2XYZ(
double xyz[3],			/* Return XYZ value with Y == 1 */
icxObserverType obType,	/* Observer */
xspect *custObserver[3],/* Optional custom observer */
icxIllumeType ilType,	/* Type of illuminant */
double ct,				/* Input temperature in degrees K */
xspect *custIllum);		/* Optional custom illuminant */


/* Given a choice of temperature dependent illuminant (icxIT_Dtemp or icxIT_Ptemp), */
/* return the closest correlated color temperature to the given spectrum or XYZ. */
/* An observer type can be chosen for interpretting the spectrum of the input and */
/* the illuminant. */
/* Note we're using CICDE94, rather than the traditional L*u*v* 2/3 space for CCT */
/* Return -1 on erorr */
double icx_XYZ2ill_ct(
double txyz[3],			/* If not NULL, return the XYZ of the black body temperature */
icxIllumeType ilType,	/* Type of illuminant, icxIT_Dtemp or icxIT_Ptemp */
icxObserverType obType,	/* Observer */
xspect *custObserver[3],/* Optional custom observer */
double xyz[3],			/* Input XYZ value, NULL if spectrum intead */
xspect *insp0,			/* Input spectrum value, NULL if xyz[] instead */
int viscct);			/* nz to use visual CIEDE2000, 0 to use CCT CIE 1960 UCS. */

/* Compute the CIE1995 CRI: Ra */
/* Return < 0.0 on error */
/* If invalid is not NULL, set it to nz if CRI */
/* is invalid because the sample is not white enough. */
double icx_CIE1995_CRI(
int *invalid,			/* if not NULL, set to nz if invalid */
xspect *sample			/* Illuminant sample to compute CRI of */
);

#endif /* XSPECTFM_H */






































