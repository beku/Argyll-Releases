#ifndef INSTTYPES_H

/* 
 * Argyll Color Correction System
 *
 * Instrument suported types utilities.
 *
 * Author: Graeme W. Gill
 * Date:   15/3/2001
 *
 * Copyright 2001 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 */

/* ----------------------------- */
/* Possible types of instruments */
typedef enum {
    instUnknown      = -1,		/* Undefined Instrument */
    instDTP20        = 0,		/* Xrite DTP20 (Pulse)  */
    instDTP22,					/* Xrite DTP22 (Digital Swatchbook)  */
    instDTP41,       			/* Xrite DTP41 */
    instDTP51, 					/* Xrite DTP51 */
    instDTP92, 					/* Xrite DTP92 */
    instDTP94, 					/* Xrite DTP94 (Optix) */
    instSpectrolino, 			/* GretagMacbeth Spectrolino */
    instSpectroScan, 			/* GretagMacbeth SpectroScan */
    instSpectroScanT, 			/* GretagMacbeth SpectroScanT */
	instSpectrocam,				/* Avantes Spectrocam */
	instI1Display,				/* GretagMacbeth i1 Display */
	instI1Monitor,				/* GretagMacbeth i1 Monitor */
	instI1Pro,					/* GretagMacbeth i1 Pro */
	instColorMunki,				/* X-Rite ColorMunki */
	instHCFR,					/* Colorimtre HCFR */
	instSpyder2,				/* Datacolor/ColorVision Spyder2 */
	instSpyder3,				/* Datacolor Spyder3 */
	instHuey,					/* GretagMacbeth Huey */
} instType;

/* Utility functions in libinsttypes */

/* Return the instrument identification name (static string), */
/* given its instrument type. */
extern char *inst_name(instType itype);


/* Given an instrument identification name, return the matching */
/* instType, or instUnknown if not matched */
extern instType inst_enum(char *name);


#ifdef ENABLE_USB
/* Given a USB vendor and product ID, */
/* return the matching instrument type, or */
/* instUnknown if none match. */
extern instType inst_usb_match(
unsigned short idVendor,
unsigned short idProduct);
#endif /* ENABLE_USB */


/* Should deprecate the following. It should be replaced with a */
/* method in the instrument class that returns its configured spectrum, */
/* and the spectrum should be embedded in the .ti3 file, not the instrument */
/* name. */

#ifdef NOXSPECT
#    undef NOXSPECT
#else

/* Fill in an instruments illuminant spectrum. */
/* Return 0 on sucess, 1 if not not applicable. */
extern int inst_illuminant(xspect *sp, instType itype);

#endif

#define INSTTYPES_H
#endif /* INSTTYPES_H */

