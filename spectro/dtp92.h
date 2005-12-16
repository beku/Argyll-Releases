#ifndef DTP92_H

/* 
 * Argyll Color Correction System
 *
 * Xrite DTP92 related defines
 *
 * Author: Graeme W. Gill
 * Date:   5/10/96
 *
 * Copyright 2001, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include "inst.h"

/* Fake Error codes */
#define DTP92_INTERNAL_ERROR		0x61			/* Internal software error */
#define DTP92_COMS_FAIL				0x62			/* Communication failure */
#define DTP92_UNKNOWN_MODEL			0x63			/* Not a DPT92 */
#define DTP92_DATA_PARSE_ERROR  	0x64			/* Read data parsing error */
#define DTP92_USER_ABORT		    0x65			/* User hit escape */

/* Real error code */
#define DTP92_OK   					0x00

#define DTP92_BAD_COMMAND			0x01
#define DTP92_PRM_RANGE				0x02
#define DTP92_MEMORY_OVERFLOW		0x04
#define DTP92_INVALID_BAUD_RATE		0x05
#define DTP92_TIMEOUT				0x07
#define DTP92_SYNTAX_ERROR			0x08
#define DTP92_NO_DATA_AVAILABLE		0x0B
#define DTP92_MISSING_PARAMETER		0x0C
#define DTP92_CALIBRATION_DENIED	0x0D
#define DTP92_NEEDS_OFFSET_CAL		0x16
#define DTP92_NEEDS_RATIO_CAL		0x17
#define DTP92_NEEDS_LUMINANCE_CAL	0x18
#define DTP92_NEEDS_WHITE_POINT_CAL	0x19

#define DTP92_NEEDS_BLACK_POINT_CAL	0x2A
#define DTP92_INVALID_READING		0x20
#define DTP92_BAD_COMP_TABLE		0x25
#define DTP92_TOO_MUCH_LIGHT		0x28
#define DTP92_NOT_ENOUGH_LIGHT		0x29

#define DTP92_BAD_SERIAL_NUMBER		0x40

#define DTP92_NO_MODULATION			0x50

#define DTP92_EEPROM_FAILURE		0x70
#define DTP92_FLASH_WRITE_FAILURE	0x71
#define DTP92_INST_INTERNAL_ERROR	0x7F


/* DTP92 communication object */
struct _dtp92 {
	INST_OBJ_BASE

	/* *** DTP92 private methods *** */
	/* Do a command/response echange with the instrument */
	/* Return the inst error code */
	inst_code (*command)(struct _inst *p, char *in, char *out, int bsize);

	}; typedef struct _dtp92 dtp92;

/* Constructor */
extern dtp92 *new_dtp92(void);



#define DTP92_H
#endif /* DTP92_H */
