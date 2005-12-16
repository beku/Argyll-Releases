
/* 
 * cam02
 *
 * Color Appearance Model, based on
 * CIECAM02, "The CIECAM02 Color Appearance Model"
 * by Nathan Moroney, Mark D. Fairchild, Robert W.G. Hunt, Changjun Li,
 * M. Ronnier Luo and Todd Newman, IS&T/SID Tenth Color Imaging
 * Conference, with the addition of the Viewing Flare
 * model described on page 487 of "Digital Color Management",
 * by Edward Giorgianni and Thomas Madden, and the
 * Helmholtz-Kohlraush effect, using the equation
 * the Bradford-Hunt 96C model as detailed in Mark Fairchilds
 * book "Color Appearance Models". 
 *
 * Author:  Graeme W. Gill
 * Date:    17/1/2004
 * Version: 1.00
 *
 * This file is based on cam97s3.h by Graeme Gill.
 *
 * Copyright 2004 Graeme W. Gill
 * Please refer to COPYRIGHT file for details.
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the LICENCE.TXT file for licencing details.
 */

/* Definitions assumed here:

  Stimulus field is the 2 degrees field of view of the sample
	being characterised.

  Viewing/Scene/Image field is the area of the whole image or
	display surface that the stimulus is part of.

  Background field is 10-12 degree field surrounding the stimulus field.
	This may be within, overlap or encompass the Viewing/Scene/Image field.

  Surround/Adapting field is the visual field minus the background field.

  Visual field is the 130 degree angular field that is seen by the eyes.

  Illuminating field is the field that illuminates the reflective
  Scene/Image. It may be the same as the Ambient field.

  Ambient field is the whole surrounding environmental light field.

  NOTE: In "Digital Color Management", Giorgianni and Madden use the term
  "Surround" to mean the same thing as "Background" in the CIECAM02 terminology.
  The ICC standard doesn't define what it means by "Surround illumination".

*/

/* Rules of Thumb: */
/* Ambient Luminance (Lat, cd/m^2) is Ambient Illuminance (Lux) divided by PI. */
/* i.e. Lat = Iat/PI */

/* The Adapting/Surround Luminance is often taken to be */
/* the 20% of the Ambient Luminance. (gray world) */
/* i.e. La = Lat/5  */

/* For a reflective print, the Viewing/Scene/Image luminance (Lv, cd/m^2), */
/* will be the Illuminating Luminance (Li, cd/m^2) reflected by the */
/* media white point (Yw) */

/* If there is no special illumination for a reflective print, */
/* then the Illuminating Luminance (Li) will be the Ambient Luminance (Lat) */

/* An emisive display will have an independently determined Lv. */

/* The classification of the type of surround is */
/* determined by comparing the Adapting/Surround luminance (La, cd/m^2) */
/* with the average luminance of the Viewing/Scene/Image field (Lv, cd/m^2) */

/* La/Lv == 0%,   dark */
/* La/Lv 0 - 20%, dim */
/* La/Lv > 20%,   average */
/* special,       cut sheet */

/* The Background relative luminance Yb is typically assumed to */
/* be 0.18 .. 0.2, and is assumed to be grey. */

/* The source of flare light depends on the type of display system. */
/* For a CRT, it will be the Ambient light reflecting off the glass surface. */
/* (This implies Yf = Lat * reflectance/Lv) */
/* For a reflection print, it will be the Illuminant reflecting from the media */
/* surface. (Yf = Li * reflectance) */
/* For a projected image, it will be stray projector light, scattered by the */
/* surround, screen and air particles. (Yf = Li * reflectance_and_scattering) */

/*

  Typical Field brightness
  (cd/m2) Condition 
    30	Subdued indoor lighting 
    60	Less than typical office light; sometimes recommended for
  		display-only workplaces 
   120	Typical office 
   240	Bright indoor office 
   480	Very bright; precision indoor tasks 
   960	Usual outdoors 
  1920	Bright afternoon 

*/

/* ---------------------------------- */

struct _cam02 {
/* Public: */
	void (*del)(struct _cam02 *s);	/* We're done with it */

	int (*set_view)(
		struct _cam02 *s,
		ViewingCondition Ev,	/* Enumerated Viewing Condition */
		double Wxyz[3],	/* Reference/Adapted White XYZ (Y scale 1.0) */
		double La,		/* Adapting/Surround Luminance cd/m^2 */
		double Yb,		/* Luminance of Background relative to reference white (range 0.0 .. 1.0) */
		double Lv,		/* Luminance of white in the Viewing/Scene/Image field (cd/m^2) */
						/* Ignored if Ev is set */
		double Yf,		/* Flare as a fraction of the reference white (range 0.0 .. 1.0) */
		double Fxyz[3],	/* The Flare white coordinates (typically the Ambient color) */
		int hk			/* Flag, NZ to use Helmholtz-Kohlraush effect */
	);

	/* Conversions */
	int (*XYZ_to_cam)(struct _cam02 *s, double *out, double *in);
	int (*cam_to_XYZ)(struct _cam02 *s, double *out, double *in);

/* Private: */
	/* Scene parameters */
	ViewingCondition Ev;	/* Enumerated Viewing Condition */
	double Wxyz[3];	/* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
	double Yb;		/* Relative Luminance of Background to reference white (Y range 0.0 .. 1.0) */
	double La;		/* Adapting/Surround Luminance cd/m^2 */
	double Yf;		/* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
	double Fxyz[3];	/* The Flare white coordinates (typically the Ambient color) */

	/* Internal parameters */
	double  C;		/* Surround Impact */
	double Nc;		/* Chromatic Induction */
	double  F;		/* Adaptation Degree */

	/* Pre-computed values */
	double Fsc;			/* Flare scale */
	double Fisc;		/* Inverse flare scale */
	double Fsxyz[3];	/* Scaled Flare color coordinates */
	double rgbW[3];		/* Sharpened cone response white values */
	double D;			/* Degree of chromatic adaption */
	double Drgb[3];		/* Chromatic transformation value */
	double rgbcW[3];	/* Chromatically transformed white value */
	double rgbpW[3];	/* Hunt-Pointer-Estevez cone response space white */
	double n;			/* Background induction factor */
	double nn;			/* Precomuted function of n */
	double Fl;			/* Lightness contrast factor ?? */
	double Nbb;			/* Background brightness induction factors */
	double Ncb;			/* Chromatic brightness induction factors */
	double z;			/* Base exponential nonlinearity */
	double rgbaW[3];	/* Post-adapted cone response of white */
	double Aw;			/* Achromatic response of white */

	/* Option flags */
	int hk;				/* Use Helmholtz-Kohlraush effect */

}; typedef struct _cam02 cam02;


/* Create a cam02 conversion class, with default viewing conditions */
cam02 *new_cam02(void);
