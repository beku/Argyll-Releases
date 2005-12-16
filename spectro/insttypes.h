#ifndef INSTTYPES_H

/* 
 * Argyll Color Correction System
 *
 * Instrument suported types.
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

/* ----------------------------- */
/* Possible types of instruments */
typedef enum {
    instUnknown      = -1,		/* Undefined Instrument */
    instDTP41        = 0,		/* Xrite DTP41 */
    instDTP51, 					/* Xrite DTP51 */
    instDTP92, 					/* Xrite DTP92 */
    instSpectrolino, 			/* GretagMacbeth Spectrolino */
    instSpectroScan, 			/* GretagMacbeth SpectroScan */
    instSpectroScanT, 			/* GretagMacbeth SpectroScanT */
	instSpectrocam				/* Avantes Spectrocam */
} instType;

/* Utility functions */

/* Return the instrument identification name (static string) */
extern char *inst_name(instType itype);

/* Given an instrument identification name, return the matching */
/* instType, or instUnknown if not matched */
extern instType inst_enum(char *name);

/* Return an instruments illuminant spectrum. Return NULL */
/* if not applicable. */
extern xspect *inst_illuminant(instType itype);

#define INSTTYPES_H
#endif /* INSTTYPES_H */
