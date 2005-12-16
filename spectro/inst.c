/* 
 * Argyll Color Correction System
 *
 * Abstract instrument class implemenation
 *
 * Author: Graeme W. Gill
 * Date:   10/3/2001
 *
 * Copyright 2001 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include "xspect.h"
#include "serio.h"
#include "inst.h"
#include "dtp41.h"
#include "dtp51.h"
#include "dtp92.h"
#include "gretag.h"
#include "ss.h"

/* inst error code interpretation */
static char *interp_code(struct _inst *p, inst_code ev) {

	switch(ev & inst_mask) {
		case inst_ok:
			return "No error";
		case inst_notify:
			return "Notification";
		case inst_warning:
			return "Warning";
		case inst_internal_error:
			return "Internal software error";
		case inst_coms_fail:
			return "Communications failure";
		case inst_unknown_model:
			return "Not expected instrument model";
		case inst_protocol_error:
			return "Communication protocol breakdown";
		case inst_user_abort:
			return "User hit Escape";
		case inst_misread:
			return "Measurement misread";
		case inst_needs_cal:
			return "Instrument needs calibration";
		case inst_unsupported:
			return "Unsupported function";
		case inst_unexpected_reply:
			return "Unexpected Reply";
		case inst_wrong_config:
			return "Wrong Configuration";
		case inst_hardware_fail:
			return "Hardware Failure";
		case inst_other_error:
			return "Non-specific error";
	}
	return "Unknown inst error code";

}

/* Default instrument user interaction callback. */
/* This will use stdio. */
static int inst_user(
	inst *p,
	inst_uiact it,		/* Interaction type */
	int nm,				/* Number of messages */
	int nch,			/* Number of choices to accept */
    ...					/* sequence of pointers to strings */
) {
	va_list args;
	int i, rv = 0;
	char c;

	if (nch > 35)
		nch = 35;

	if (nch > nm)
		nch = nm;

	for (;;) {
		va_start(args, nch);
		for (i = 0; i < nm; i++) {
			char *mes = va_arg(args, char *);
			if ((nm-i) <= nch) {
				int j = i-nm+nch+1;
				if (j < 10)
					c = '0' + j;
				else
					c = 'A' + j - 10; 
				printf(" %c: %s",c,mes);
			} else {
				printf("%s",mes);
			}
			
		}
		va_end(args);
		if (nch <= 0)
			return rv;

		fflush(stdout);
		empty_con_chars();
		c = next_con_char();
		printf("\n");
		if (c >= '0' && c <= '9')
			rv = c - '0';
		else if (c >= 'A' && c <= 'Z')
			rv = c - 'A' + 10;
		else if (c >= 'a' && c <= 'z')
			rv = c - 'a' + 10;
		else
			rv = 0;
		if (rv >= 1 && rv <= nch)
			return rv;
	}
}

/* Virtual constructor */
/* Return NULL for unknown instrument */
extern inst *new_inst(
instType itype) {
	inst *p;

	if (itype == instDTP41)
		p = (inst *)new_dtp41();
	else if (itype == instDTP51)
		p = (inst *)new_dtp51();
	else if (itype == instDTP92)
		p = (inst *)new_dtp92();
	else if (itype == instSpectrolino ||
			 itype == instSpectroScan ||
			 itype == instSpectroScanT)
		p = (inst *)new_gretag();
/* NYI
	else if (itype == instSpectrocam)
		p = (inst *)new_spc();
*/
	else
		return NULL;

	p->inst_interp_error = interp_code;

	/* Provide default user interaction callback */
	if (p->user == NULL)
		p->user = inst_user;

	return p;
}


/* Utility functions */

/* Return the instrument identification name (static string) */
char *inst_name(instType itype) {
	switch (itype) {
		case instDTP41:
			return "Xrite DTP41";
		case instDTP51:
			return "Xrite DTP51";
		case instDTP92:
			return "Xrite DTP92";
		case instSpectrolino:
			return "GretagMacbeth Spectrolino";
		case instSpectroScan:
			return "GretagMacbeth SpectroScan";
		case instSpectroScanT:
			return "GretagMacbeth SpectroScanT";
		case instSpectrocam:
			return "Spectrocam";
	}
	return "Unknown Instrument";
}

/* Given an instrument identification name, return the matching */
/* instType, or instUnknown if not matched */
extern instType inst_enum(char *name) {

	if (strcmp(name, "Xrite DTP41") == 0)
		return instDTP41;
	else if (strcmp(name, "Xrite DTP51") == 0)
		return instDTP51;
	else if (strcmp(name, "Xrite DTP92") == 0)
		return instDTP92;
	else if (strcmp(name, "GretagMacbeth Spectrolino") == 0)
		return instSpectrolino;
	else if (strcmp(name, "GretagMacbeth SpectroScan") == 0)
		return instSpectroScan;
	else if (strcmp(name, "GretagMacbeth SpectroScanT") == 0)
		return instSpectroScanT;
	else if (strcmp(name, "Spectrocam") == 0)
		return instSpectrocam;

	return instUnknown;
}

/* Return an instruments illuminant spectrum. Return NULL */
/* if not applicable. */
xspect *inst_illuminant(instType itype) {

	switch (itype) {
		case instDTP41:
			return standardIlluminant(icxIT_A);		/* 2850K */

		case instDTP51:
			return standardIlluminant(icxIT_A);		/* 2850K */

		case instDTP92:
			return NULL;							/* Not applicable */

		/* Strictly the Spectrolino could have the UV or D65 filter on, */
		/* but since we're not currently passing this inform through, we */
		/* are simply assuming the default A type illuminant. */

		case instSpectrolino:
			return standardIlluminant(icxIT_A);		/* Standard A type assumed */

		case instSpectroScan:
			return standardIlluminant(icxIT_A);		/* Standard A type assumed */

		case instSpectroScanT:
			return standardIlluminant(icxIT_A);		/* Standard A type assumed */

		case instSpectrocam:
			return standardIlluminant(icxIT_Spectrocam);    /* Spectrocam Xenon Lamp */
	}
	return NULL;
}


