
/* 
 * Argyll Color Correction System
 *
 * ColorVision Spyder 2 related software.
 *
 * Author: Graeme W. Gill
 * Date:   17/9/2007
 *
 * Copyright 2006 - 2007, Graeme W. Gill
 * All rights reserved.
 *
 * (Based on i1disp.c)
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
	IMPORTANT NOTE:

    This instrument cannot function without the driver software
	having access to the vendor supplied PLD firmware pattern for it.
    This firmware is not provided with Argyll, since it is not available
    under a compatible license.

    The purchaser of a Spyder 2 instrument should have received a copy
    of this firmware along with their instrument, and should therefore be able to
    enable the Argyll driver for this instrument by using the spyd2en utility
	to create a spyd2PLD.bin file.

 */

/* 
   If you make use of the instrument driver code here, please note
   that it is the author(s) of the code who take responsibility
   for its operation. Any problems or queries regarding driving
   instruments with the Argyll drivers, should be directed to
   the Argyll's author(s), and not to any other party.

   If there is some instrument feature or function that you
   would like supported here, it is recommended that you
   contact Argyll's author(s) first, rather than attempt to
   modify the software yourself, if you don't have firm knowledge
   of the instrument communicate protocols. There is a chance
   that an instrument could be damaged by an incautious command
   sequence, and the instrument companies generally cannot and
   will not support developers that they have not qualified
   and agreed to support.
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <stdarg.h>
#include <math.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "xspect.h"
#include "insttypes.h"
#include "icoms.h"
#include "spyd2.h"

#undef DEBUG

#ifdef DEBUG
#define DBG(xxx) printf xxx ;
#else
#define DBG(xxx) 
#endif

#define CLKRATE 1000000			/* Clockrate the Spyder 2 hardware runs at */
#define MSECDIV (CLKRATE/1000)	/* Divider to turn clocks into msec */
#define DEFRRATE 50				/* Default display refresh rate */
#define RDTIME 2.0				/* Time to integrate reading over */ 
#define DO_ADAPTIVE				/* Adapt the integration time to the light level */
								/* This helps repeatability at low levels A LOT */
#define LEVEL2					/* Second level (nonliniarity) calibration */
#define RETRIES 4				/* usb_reads are unreliable - bug in spyder H/W ?*/

static inst_code spyd2_interp_code(inst *pp, int ec);

/* ------------------------------------------------------------------------ */
/* Implementation */

/* Interpret an icoms error into a SPYD2 error */
static int icoms2spyd2_err(int se) {
	if (se & ICOM_USERM) {
		se &= ICOM_USERM;
		if (se == ICOM_USER)
			return SPYD2_USER_ABORT;
		if (se == ICOM_TERM)
			return SPYD2_USER_TERM;
		if (se == ICOM_TRIG)
			return SPYD2_USER_TRIG;
		if (se == ICOM_CMND)
			return SPYD2_USER_CMND;
	}
	if (se != ICOM_OK)
		return SPYD2_COMS_FAIL;
	return SPYD2_OK;
}

/* ----------------------------------------------------- */
/* Big endian wire format conversion routines */

/* Take an int, and convert it into a byte buffer big endian */
static void int2buf(unsigned char *buf, int inv) {
	buf[0] = (inv >> 24) & 0xff;
	buf[1] = (inv >> 16) & 0xff;
	buf[2] = (inv >> 8) & 0xff;
	buf[3] = (inv >> 0) & 0xff;
}

/* Take a short, and convert it into a byte buffer big endian */
static void short2buf(unsigned char *buf, int inv) {
	buf[0] = (inv >> 8) & 0xff;
	buf[1] = (inv >> 0) & 0xff;
}

/* Take a word sized buffer, and convert it to an int */
static int buf2int(unsigned char *buf) {
	int val;
	val = buf[0];
	val = ((val << 8) + (0xff & buf[1]));
	val = ((val << 8) + (0xff & buf[2]));
	val = ((val << 8) + (0xff & buf[3]));
	return val;
}

/* Take a 3 byte buffer, and convert it to an unsigned int */
static unsigned int buf2uint24(unsigned char *buf) {
	unsigned int val;
	val = buf[0];
	val = ((val << 8) + (0xff & buf[1]));
	val = ((val << 8) + (0xff & buf[2]));
	return val;
}

/* Take a 24 bit unsigned sized buffer in little endian */
/* nibble swapped format, and return an int */
static unsigned int buf2uint24lens(unsigned char *buf) {
	unsigned int val;
	val =              (0xf &  buf[2]);
	val = (val << 4) + (0xf & (buf[2] >> 4));
	val = (val << 4) + (0xf &  buf[1]);
	val = (val << 4) + (0xf & (buf[1] >> 4));
	val = (val << 4) + (0xf &  buf[0]);
	val = (val << 4) + (0xf & (buf[0] >> 4));
	return val;
}

/* Take a short sized buffer, and convert it to an int */
static int buf2short(unsigned char *buf) {
	int val;
	val = buf[0];
	val = ((val << 8) + (0xff & buf[1]));
	return val;
}

/* Take a unsigned short sized buffer, and convert it to an int */
static int buf2ushort(unsigned char *buf) {
	int val;
	val = (0xff & buf[0]);
	val = ((val << 8) + (0xff & buf[1]));
	return val;
}

/* ============================================================ */
/* Low level commands */

/* USB Instrument commands */

/* Reset the instrument */
inst_code
spyd2_reset(
	spyd2 *p
) {
	int rwbytes;			/* Data bytes read or written */
	int se;
	int isdeb = 0;
	int retr;
	inst_code rv = inst_ok;

	/* Turn off low level debug messages, and sumarise them here */
	isdeb = p->icom->debug;
	p->icom->debug = 0;

	if (isdeb) fprintf(stderr,"\nspyd2: Instrument reset\n");

	for (retr = 0; ; retr++) {
		se = p->icom->usb_control(p->icom,
		               USB_ENDPOINT_OUT | USB_TYPE_VENDOR | USB_RECIP_DEVICE,
	                   0xC7, 0, 0, NULL, 0, 5.0);

		if (se == ICOM_OK) {
			if (isdeb) fprintf(stderr,"Reset complete, ICOM err 0x%x\n",se);
			break;
		}
		if ((se & ICOM_USERM) || retr >= RETRIES ) {
			if (isdeb) fprintf(stderr,"\nspyd2: Reset failed with  ICOM err 0x%x\n",se);
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
		}
		msec_sleep(500);
		if (isdeb) fprintf(stderr,"\nspyd2: Reset retry with  ICOM err 0x%x\n",se);
	}

	p->icom->debug = isdeb;
	return rv;
}

/* Get status */
/* return pointer may be NULL if not needed. */
inst_code
spyd2_getstatus(
	spyd2 *p,
	int *stat	/* Return the 1 byte status code */
) {
	int rwbytes;			/* Data bytes read or written */
	unsigned char pbuf[8];	/* status bytes read */
	int _stat;
	int se;
	int isdeb = 0;
	int retr;
	inst_code rv = inst_ok;

	/* Turn off low level debug messages, and sumarise them here */
	isdeb = p->icom->debug;
	p->icom->debug = 0;

	if (isdeb) fprintf(stderr,"\nspyd2: Get Status\n");

	for (retr = 0; ; retr++) {
		se = p->icom->usb_control(p->icom,
		               USB_ENDPOINT_IN | USB_TYPE_VENDOR | USB_RECIP_DEVICE,
	                   0xC6, 0, 0, pbuf, 8, 5.0);

		if (se == ICOM_OK)
			break;
		if ((se & ICOM_USERM) || retr >= RETRIES ) {
			if (isdeb) fprintf(stderr,"\nspyd2: Get Status failed with ICOM err 0x%x\n",se);
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
		}
		msec_sleep(500);
		if (isdeb) fprintf(stderr,"\nspyd2: Get Status retry with ICOM err 0x%x\n",se);
	}
	msec_sleep(100);		/* Limit rate status commands can be given */

	_stat    = pbuf[0];		/* Only the first byte is examined. */
							/* Other bytes have information, but SW ignores them */

	if (isdeb) fprintf(stderr,"Get Status returns %d ICOM err 0x%x\n", _stat, se);

	p->icom->debug = isdeb;

	if (stat != NULL) *stat = _stat;

	return rv;
}

/* Read Serial EEProm bytes (implementation) */
/* Can't read more than 256 in one go */
inst_code
spyd2_readEEProm_imp(
	spyd2 *p,
	unsigned char *buf,	/* Buffer to return bytes in */
	int addr,	/* Serial EEprom address, 0 - 511 */
	int size	/* Number of bytes to read, 0 - 128 (ie. max of  bank) */
) {
	int rwbytes;			/* Data bytes read or written */
	int se;
	int isdeb = 0;
	int retr;
	inst_code rv = inst_ok;

	/* Turn off low level debug messages, and sumarise them here */
	isdeb = p->icom->debug;
	p->icom->debug = 0;

	if (isdeb >= 2) fprintf(stderr,"\nspyd2: Read EEProm addr %d, bytes %d\n",addr,size);

	if (addr < 0 || (addr + size) > 512)
		return spyd2_interp_code((inst *)p, SPYD2_BAD_EE_ADDRESS);

	if (size >= 256)
		return spyd2_interp_code((inst *)p, SPYD2_BAD_EE_SIZE);

	for (retr = 0; ; retr++) {
		se = p->icom->usb_control(p->icom,
			               USB_ENDPOINT_IN | USB_TYPE_VENDOR | USB_RECIP_DEVICE,
		                   0xC4, addr, size, buf, size, 5.0);
		if (se == ICOM_OK)
			break;
		if ((se & ICOM_USERM) || retr >= RETRIES ) {
			if (isdeb) fprintf(stderr,"\nspyd2: Read bytes failed with ICOM err 0x%x\n",se);
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
		}
		msec_sleep(500);
		if (isdeb) fprintf(stderr,"\nspyd2: Read bytes retry with ICOM err 0x%x\n",se);
	}


	if (isdeb >= 2) fprintf(stderr,"Read EEProm ICOM err 0x%x\n", se);

	p->icom->debug = isdeb;

	return rv;
}

/* Read Serial EEProm bytes */
/* (Handles reads > 256 bytes) */
inst_code
spyd2_readEEProm(
	spyd2 *p,
	unsigned char *buf,	/* Buffer to return bytes in */
	int addr,	/* Serial EEprom address, 0 - 511 */
	int size	/* Number of bytes to read, 0 - 511 */
) {

	if (addr < 0 || (addr + size) > 512)
		return spyd2_interp_code((inst *)p, SPYD2_BAD_EE_ADDRESS);
	
	while (size > 255) {		/* Single read is too big */
		inst_code rv;
		if ((rv = spyd2_readEEProm_imp(p, buf, addr, 255)) != inst_ok)
			return rv;
		size -= 255;
		buf  += 255;
		addr += 255;
	}
	return spyd2_readEEProm_imp(p, buf, addr, size);
}

/* Download PLD pattern */
inst_code
spyd2_loadPLD(
	spyd2 *p,
	unsigned char *buf,		/* Bytes to download */
	int size				/* Number of bytes */
) {
	int rwbytes;			/* Data bytes read or written */
	int se;
	int isdeb = 0;
	int retr;
	inst_code rv = inst_ok;

	/* Turn off low level debug messages, and sumarise them here */
	isdeb = p->icom->debug;
	p->icom->debug = 0;

	if (isdeb >= 2) fprintf(stderr,"\nspyd2: Load PLD %d bytes\n",size);

	for (retr = 0; ; retr++) {
		se = p->icom->usb_control(p->icom,
		               USB_ENDPOINT_OUT | USB_TYPE_VENDOR | USB_RECIP_DEVICE,
	                   0xC0, 0, 0, buf, size, 5.0);

		if (se == ICOM_OK)
			break;
		if ((se & ICOM_USERM) || retr >= RETRIES ) {
			if (isdeb) fprintf(stderr,"\nspyd2: Load PLD failed with ICOM err 0x%x\n",se);
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
		}
		msec_sleep(500);
		if (isdeb) fprintf(stderr,"\nspyd2: Load PLD retry with ICOM err 0x%x\n",se);
	}

	if (isdeb >= 2) fprintf(stderr,"Load PLD returns ICOM err 0x%x\n", se);

	p->icom->debug = isdeb;

	return rv;
}

/* Get minmax command */
/* Figures out the current minimum and maximum frequency periods */
/* so as to be able to set a frame detect threshold. */
/* Note it returns 0,0 if there is not enough light. */ 
/* (The light to frequency output period size is inversly */
/*  related to the lightness level) */
spyd2_GetMinMax(
	spyd2 *p,
	int *clocks,		/* Number of clocks to use (may get limited) */
	int *min,			/* Return min and max light->frequency periods */
	int *max
) {
	int rwbytes;			/* Data bytes read or written */
	int se;
	int isdeb = 0;
	inst_code rv = inst_ok;
	int value;
	int index;
	int retr;
	unsigned char buf[8];	/* return bytes read */

	/* Turn off low level debug messages, and sumarise them here */
	isdeb = p->icom->debug;
	p->icom->debug = 0;

	if (isdeb >= 2) fprintf(stderr,"\nspyd2: Get Min/Max, %d clocks\n",*clocks);

	/* Issue the triggering command */
	if (*clocks > 0xffffff)
		*clocks = 0xffffff;		/* Maximum count hardware will take ? */
	value = *clocks >> 8;
	value = (value >> 8) | ((value << 8) & 0xff00);		/* Convert to big endian */
	index = (*clocks << 8) & 0xffff;
	index = (index >> 8) | ((index << 8) & 0xff00);		/* Convert to big endian */

	for (retr = 0; ; retr++) {
		/* Issue the trigger command */
		se = p->icom->usb_control(p->icom,
		               USB_ENDPOINT_OUT | USB_TYPE_VENDOR | USB_RECIP_DEVICE,
	                   0xC2, value, index, NULL, 0, 5.0);

		if ((se & ICOM_USERM) || (se != ICOM_OK && retr >= RETRIES)) {
			if (isdeb) fprintf(stderr,"\nspyd2: Get Min/Max Trig failed with ICOM err 0x%x\n",se);
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
		}
		if (se != ICOM_OK) {
			msec_sleep(500);
			if (isdeb) fprintf(stderr,"\nspyd2: Get Min/Max Trig retry with ICOM err 0x%x\n",se);
			continue;
		}

		if (isdeb >= 2) fprintf(stderr,"Trigger Min/Max returns ICOM err 0x%x\n", se);
	
		/* Allow some time for the instrument to respond */
		msec_sleep(*clocks/MSECDIV);
	
		/* Now read the bytes */
		se = p->icom->usb_read(p->icom, 0x81, buf, 8, &rwbytes, 5.0);
		if (se == ICOM_OK)
			break;
		if ((se & ICOM_USERM) || retr >= RETRIES ) {
			if (isdeb) fprintf(stderr,"\nspyd2: Get Min/Max failed with ICOM err 0x%x\n",se);
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
		}
		msec_sleep(500);
		if (isdeb) fprintf(stderr,"\nspyd2: Get Min/Max retry with ICOM err 0x%x\n",se);
	}

	if (rwbytes != 8) {
		if (isdeb) fprintf(stderr,"\nspyd2: Get Min/Max got short data read %d",rwbytes);
		p->icom->debug = isdeb;
		return spyd2_interp_code((inst *)p, SPYD2_BADREADSIZE);
	}

	*min = buf2ushort(&buf[0]);
	*max = buf2ushort(&buf[2]);

	if (isdeb >= 2) fprintf(stderr,"Get Min/Max got %d/%d returns ICOM err 0x%x\n", *min, *max, se);
	p->icom->debug = isdeb;

	return rv;
}

/* Get refresh rate (low level) command */
spyd2_GetRefRate_ll(
	spyd2 *p,
	int *clocks,		/* Maximum number of clocks to use */
	int nframes,		/* Number of frames to count */
	int thresh,			/* Frame detection threshold */
	int *minfclks,		/* Minimum number of clocks per frame */
	int *maxfclks,		/* Maximum number of clocks per frame */
	int *clkcnt			/* Return number of clocks for nframes frames */
) {
	int rwbytes;			/* Data bytes read or written */
	int se;
	int isdeb = 0;
	inst_code rv = inst_ok;
	int value;
	int index;
	int flag;
	int retr;
	unsigned char buf1[8];	/* send bytes */
	unsigned char buf2[8];	/* return bytes read */

	/* Turn off low level debug messages, and sumarise them here */
	isdeb = p->icom->debug;
	p->icom->debug = 0;

	if (isdeb >= 2) fprintf(stderr,"\nspyd2: Get Refresh Rate, %d clocks\n",*clocks);

	/* Setup the triggering parameters */
	if (*clocks > 0xffffff)		/* Enforce hardware limits */
		*clocks = 0xffffff;
	if (*minfclks > 0x7fff)
		*minfclks = 0x7fff;
	if (*maxfclks > 0x7fff)
		*maxfclks = 0x7fff;
	value = *clocks >> 8;
	value = (value >> 8) | ((value << 8) & 0xff00);		/* Convert to big endian */
	index = (*clocks << 8) & 0xffff;
	index = (index >> 8) | ((index << 8) & 0xff00);		/* Convert to big endian */

	/* Setup parameters in send buffer */
	short2buf(&buf1[0], thresh);
	short2buf(&buf1[2], nframes);
	short2buf(&buf1[4], *minfclks);
	short2buf(&buf1[6], *maxfclks);

	/* Issue the triggering command */
	for (retr = 0; ; retr++) {
		se = p->icom->usb_control(p->icom,
		               USB_ENDPOINT_OUT | USB_TYPE_VENDOR | USB_RECIP_DEVICE,
	                   0xC3, value, index, buf1, 8, 5.0);

		if ((se & ICOM_USERM) || (se != ICOM_OK && retr >= RETRIES)) {
			if (isdeb) fprintf(stderr,"\nspyd2: Get Refresh Rate Trig failed with ICOM err 0x%x\n",se);
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
		}
		if (se != ICOM_OK) {
			msec_sleep(500);
			if (isdeb) fprintf(stderr,"\nspyd2: Get Refresh Rate Trig retry with ICOM err 0x%x\n",se);
			continue;
		}

		if (isdeb >= 2) fprintf(stderr,"Trigger Get Refresh Rate returns ICOM err 0x%x\n", se);

		/* Allow some time for the instrument to respond */
		msec_sleep(*clocks/MSECDIV);
	
		/* Now read the bytes */
		se = p->icom->usb_read(p->icom, 0x81, buf2, 8, &rwbytes, 5.0);
		if (se == ICOM_OK)
			break;
		if ((se & ICOM_USERM) || retr >= RETRIES ) {
			if (isdeb) fprintf(stderr,"\nspyd2: Get Refresh Rate failed with ICOM err 0x%x\n",se);
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
		}
		msec_sleep(500);
		if (isdeb) fprintf(stderr,"\nspyd2: Get Refresh Rate retry with ICOM err 0x%x\n",se);
	}

	if (rwbytes != 8) {
		if (isdeb) fprintf(stderr,"\nspyd2: Get Refresh Rate got short data read %d",rwbytes);
		p->icom->debug = isdeb;
		return spyd2_interp_code((inst *)p, SPYD2_BADREADSIZE);
	}

	flag = buf2[0];
	*clkcnt = buf2uint24(&buf2[1]);

	if (flag == 1) {
		if (isdeb) fprintf(stderr,"\nspyd2: Get Refresh Rate got trigger timeout");
		p->icom->debug = isdeb;
		return spyd2_interp_code((inst *)p, SPYD2_TRIGTIMEOUT);
	}

	if (flag == 2) {
		if (isdeb) fprintf(stderr,"\nspyd2: Get Refresh Rate got overall timeout");
		p->icom->debug = isdeb;
		return spyd2_interp_code((inst *)p, SPYD2_OVERALLTIMEOUT);
	}

	if (isdeb >= 2) fprintf(stderr,"Get Refresh Rate got %d, returns ICOM err 0x%x\n", *clkcnt, se);
	p->icom->debug = isdeb;

	return rv;
}

/* Get a reading (low level) command */
spyd2_GetReading_ll(
	spyd2 *p,
	int *clocks,		/* Maximum number of clocks to use */
	int nframes,		/* Number of frames being measured */
	int thresh,			/* Frame detection threshold */
	int *minfclks,		/* Minimum number of clocks per frame */
	int *maxfclks,		/* Maximum number of clocks per frame */
	double *sensv		/* Return the 8 sensor readings (may be NULL) */ 
) {
	int rwbytes;			/* Data bytes read or written */
	int se;
	int isdeb = 0;
	inst_code rv = inst_ok;
	int value;
	int index;
	int flag;
	int retr;
	unsigned char buf1[8];		/* send bytes */
	unsigned char buf2[9 * 8];	/* return bytes read */
	int rvals[3][8];			/* Raw values */
	int i, j, k;

	/* Turn off low level debug messages, and sumarise them here */
	isdeb = p->icom->debug;
	p->icom->debug = 0;

	if (isdeb >= 2) fprintf(stderr,"\nspyd2: Get Reading, %d clocks\n",*clocks);

	/* Setup the triggering parameters */
	if (*clocks > 0xffffff)		/* Enforce hardware limits */
		*clocks = 0xffffff;
	if (*minfclks > 0x7fff)
		*minfclks = 0x7fff;
	if (*maxfclks > 0x7fff)
		*maxfclks = 0x7fff;
	value = *clocks >> 8;
	value = (value >> 8) | ((value << 8) & 0xff00);		/* Convert to big endian */
	index = (*clocks << 8) & 0xffff;
	index = (index >> 8) | ((index << 8) & 0xff00);		/* Convert to big endian */

	/* Setup parameters in send buffer */
	thresh *= 256;
	int2buf(&buf1[0], thresh);
	short2buf(&buf1[4], *minfclks);
	short2buf(&buf1[6], *maxfclks);

	/* Issue the triggering command */
	for (retr = 0; ; retr++) {
		se = p->icom->usb_control(p->icom,
		               USB_ENDPOINT_OUT | USB_TYPE_VENDOR | USB_RECIP_DEVICE,
	                   0xC1, value, index, buf1, 8, 5.0);

		if ((se & ICOM_USERM) || (se != ICOM_OK && retr >= RETRIES)) {
			if (isdeb) fprintf(stderr,"\nspyd2: Get Reading Trig failed with ICOM err 0x%x\n",se);
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
		}
		if (se != ICOM_OK) {
			msec_sleep(500);
			if (isdeb) fprintf(stderr,"\nspyd2: Get Reading Trig retry with ICOM err 0x%x\n",se);
			continue;
		}

		if (isdeb >= 2) fprintf(stderr,"Trigger Get Reading returns ICOM err 0x%x\n", se);

		/* Allow some time for the instrument to respond */
		msec_sleep(*clocks/MSECDIV);
	
		/* Now read the first 8 bytes */
		se = p->icom->usb_read(p->icom, 0x81, buf2, 8, &rwbytes, 5.0);
		if ((se & ICOM_USERM) || (se != ICOM_OK && retr >= RETRIES)) {
			if (isdeb) fprintf(stderr,"\nspyd2: Get Reading failed with ICOM err 0x%x\n",se);
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
		}

		if (se != ICOM_OK) {
			msec_sleep(500);
			if (isdeb >= 2) fprintf(stderr,"Get Reading retry with ICOM err 0x%x\n", se);
			continue;
		}

		if (rwbytes != 8) {
			if (isdeb) fprintf(stderr,"\nspyd2: Get Reading got short data read %d",rwbytes);
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, SPYD2_BADREADSIZE);
		}
	
		flag = buf2[0];

		if (flag == 1) {
			if (isdeb) fprintf(stderr,"\nspyd2: Get Reading got trigger timeout");
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, SPYD2_TRIGTIMEOUT);
		}

		if (flag == 2) {
			if (isdeb) fprintf(stderr,"\nspyd2: Get Reading got overall timeout");
			p->icom->debug = isdeb;
			return spyd2_interp_code((inst *)p, SPYD2_OVERALLTIMEOUT);
		}

		/* Now read the following 9 x 8 bytes of sensor data */
		for (i = 0; i < 9; i++) {

			se = p->icom->usb_read(p->icom, 0x81, buf2 + i * 8, 8, &rwbytes, 5.0);
			if ((se & ICOM_USERM) || (se != ICOM_OK && retr >= RETRIES)) {
				if (isdeb) fprintf(stderr,"\nspyd2: Get Reading failed with ICOM err 0x%x\n",se);
				p->icom->debug = isdeb;
				return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
			}
			if (se != ICOM_OK)
				break;
			if (rwbytes != 8) {
				if (isdeb) fprintf(stderr,"\nspyd2: Get Reading got short data read %d",rwbytes);
				p->icom->debug = isdeb;
				return spyd2_interp_code((inst *)p, SPYD2_BADREADSIZE);
			}
		}

		if (i >= 9)
			break;		/* We're done */

		msec_sleep(500);
		if (isdeb) fprintf(stderr,"\nspyd2: Get Reading retry with ICOM err 0x%x\n",se);
	}

	if (sensv == NULL) {
		p->icom->debug = isdeb;
		return rv;
	}

//printf("~1 got bytes:\n");
//for (i = 0; i < 9; i++) {
//	printf("~1 %d: 0x%x 0x%x 0x%x 0x%x 0x%x 0x%x 0x%x 0x%x\n",i * 8,
//	buf[i * 8 + 0], buf[i * 8 + 1], buf[i * 8 + 2], buf[i * 8 + 3],
//	buf[i * 8 + 4], buf[i * 8 + 5], buf[i * 8 + 6], buf[i * 8 + 7]);
//}

	/* Convert the raw buffer readings into 3 groups of 8 integers. */
	/* At the start of each reading, the HW starts counting master */
	/* (1MHz) clocks. When the first transition after the start of */
	/* the reading is received from a light->frequency sensor (TAOS TSL237), */
	/* the clock count is recorded, and returned as the second of the three */
	/* numbers returned for each sensor. */
	/* When the last transition before the end of the reading period is */
	/* received, the clock count is recorded, and returned as the first */
	/* of the three numbers. The integration period is therefore */
	/* the first number minus the second. */
	/* The third number is the number of transitions from the sensor */
	/* counted during the integration period. Since this 24 bit counter is */
	/* not reset between readings, the previous count is recorded (prevraw[]), */
	/* and subtracted (modulo 24 bits) from the current value. */
	/* The light level is directly proportional to the frequency, */
	/* hence the transitions counted. */
	/* In the case of a CRT, the total number of clocks is assumed to be */
	/* set to an integer number of refresh cycles, and the total transitions */
	/* over that period are counted. */
	for (i = j = 0; j < 3; j++) {
		for (k = 0; k < 8; k++, i += 3) {
			rvals[j][k] = buf2uint24lens(buf2 + i);
//printf("~1 got rvals[%d][%d] = 0x%x\n",j,k,rvals[j][k]);
		}
	}

	/* And convert 3 integers per sensor into sensor values */
	for (k = 0; k < 8; k++) {
		int transcnt, intclks;

//printf("~1 sensor %d:\n",k);

		/* Compute difference of third integer to previous value */
		/* read, modulo 24 bits */
		if (p->prevraw[k] <= rvals[2][k]) {
			transcnt = rvals[2][k] - p->prevraw[k];
//printf("~1 transcnt = %d - %d = %d\n",rvals[2][k], p->prevraw[k],transcnt);
		} else {
			transcnt = rvals[2][k] + 0x1000000 - p->prevraw[k];
//printf("~1 transcnt = %d + %d - %d = %d\n",rvals[2][k], 0x1000000, p->prevraw[k],transcnt);
		}

		p->prevraw[k] = rvals[2][k];	/* New previuos value */

		/* Compute difference of first integer to second */
		intclks = rvals[0][k] - rvals[1][k];
//printf("~1 intclks = %d - %d = %d\n",rvals[0][k], rvals[1][k]);

		if (transcnt == 0 || intclks == 0) {		/* It's too dark ... */
			sensv[k] = 0.0;
		} else {			/* Transitions within integration period */
							/* hence one is discarded ? */
			sensv[k] = ((double)transcnt - 1.0) * (double)CLKRATE/(double)intclks;
//printf("~1 %d: initial senv %f\n",k,sensv[k]);
		}

#ifdef NEVER	/* This seems to make repeatability worse ??? */
// ~~99999
		/* If CRT and bright enough */
		if (sensv[k] > 1.5 && p->lcd == 0) {
			sensv[k] = ((double)transcnt) * (double)p->rrate/(double)nframes;
//printf("~1 %d: corrected senv %f\n",k,sensv[k]);
		}
#endif
	}

	p->icom->debug = isdeb;
	return rv;
}

/* ============================================================ */
/* Medium level commands */

/* Read a 16 bit word from the EEProm */
static inst_code
spyd2_rd_ee_ushort(
	spyd2 *p,				/* Object */
	unsigned int *outp,		/* Where to write value */
	int addr				/* EEprom Address, 0 - 510 */
) {
	inst_code ev;
	unsigned char buf[2];
	int v, val;

	if ((ev = spyd2_readEEProm(p, buf, addr, 2)) != inst_ok)
		return ev;

	*outp = buf2ushort(buf);

	return inst_ok;
}

/* Read a 32 bit word from the EEProm */
static inst_code
spyd2_rd_ee_int(
	spyd2 *p,				/* Object */
	int *outp,				/* Where to write value */
	int addr				/* EEprom Address, 0 - 508 */
) {
	inst_code ev;
	unsigned char buf[4];
	int v, val;

	if ((ev = spyd2_readEEProm(p, buf, addr, 4)) != inst_ok)
		return ev;

	*outp = buf2int(buf);

	return inst_ok;
}

/* Read a float from the EEProm */
static inst_code
spyd2_rdreg_float(
	spyd2 *p,				/* Object */
	double *outp,			/* Where to write value */
	int addr				/* Register Address, 0 - 508 */
) {
	inst_code ev;
	int val;

	if ((ev = spyd2_rd_ee_int(p, &val, addr)) != inst_ok)
		return ev;

	*outp = IEEE754todouble((unsigned int)val);
	return inst_ok;
}

/* Special purpose float read, */
/* Read three 9 vectors of floats from the EEprom */
static inst_code
spyd2_rdreg_3x9xfloat(
	spyd2 *p,				/* Object */
	double *out0,			/* Where to write first 9 doubles */
	double *out1,			/* Where to write second 9 doubles */
	double *out2,			/* Where to write third 9 doubles */
	int addr				/* Register Address, 0 - 404 */
) {
	inst_code ev;
	unsigned char buf[3 * 9 * 4], *bp;
	int i;

	if ((ev = spyd2_readEEProm(p, buf, addr, 3 * 9 * 4)) != inst_ok)
		return ev;

	bp = buf;
	for (i = 0; i < 9; i++, bp +=4, out0++) {
		int val;
		val = buf2int(bp);
		*out0 = IEEE754todouble((unsigned int)val);
	}

	for (i = 0; i < 9; i++, bp +=4, out1++) {
		int val;
		val = buf2int(bp);
		*out1 = IEEE754todouble((unsigned int)val);
	}

	for (i = 0; i < 9; i++, bp +=4, out2++) {
		int val;
		val = buf2int(bp);
		*out2 = IEEE754todouble((unsigned int)val);
	}

	return inst_ok;
}

/* Get refresh rate command. Return 0.0 */
/* if no refresh rate can be established */
spyd2_GetRefRate(
	spyd2 *p,
	double *refrate		/* return the refresh rate */
) {
	inst_code ev;
	int clocks;			/* Clocks to run commands */
	int min, max;		/* min and max light intensity frequency periods */

	if (p->debug) fprintf(stderr,"spyd2: about to get the refresh rate\n");

	/* Establish the frame rate detect threshold level */
	clocks = (10 * CLKRATE)/DEFRRATE;

	if ((ev = spyd2_GetMinMax(p, &clocks, &min, &max)) != inst_ok)
		return ev; 

	if (min == 0 || max < (5 * min)) {
		if (p->debug) fprintf(stderr,"spyd2: no refresh rate detectable\n");
		*refrate = 0.0;
	} else {
		int frclocks;		/* notional clocks per frame */
		int nframes;		/* Number of frames to count */
		int thresh;			/* Frame detection threshold */
		int minfclks;		/* Minimum number of clocks per frame */
		int maxfclks;		/* Maximum number of clocks per frame */
		int clkcnt;		/* Return number of clocks for nframes frames */

		frclocks = CLKRATE/DEFRRATE;
		nframes = 50;
		thresh = (max - min)/5 + min;			/* Threshold is at 80% of max brightness */
		minfclks = frclocks/3;					/* Allow for 180 Hz */
		maxfclks = (frclocks * 5)/2;			/* Allow for 24 Hz */
		clocks = nframes * frclocks * 2;		/* Allow for 120 Hz */

		if ((ev = spyd2_GetRefRate_ll(p, &clocks, nframes, thresh, &minfclks, &maxfclks,
		                                 &clkcnt)) != inst_ok)
			return ev; 

		/* Compute the refresh rate */
		*refrate = ((double)nframes * (double)CLKRATE)/(double)clkcnt;
		if (p->debug) fprintf(stderr,"spyd2: refresh rate is %f Hz\n",*refrate);
	}
	return ev;
}

/* Do a reading */
spyd2_GetReading(
	spyd2 *p,
	double *XYZ		/* return the XYZ values */
) {
	inst_code ev;
	int clocks1, clocks2;	/* Clocks to run commands */
	int min, max;		/* min and max light intensity frequency periods */
	int frclocks;		/* notional clocks per frame */
	int nframes;		/* Number of frames to measure over */
	int thresh;			/* Frame detection threshold */
	int minfclks;		/* Minimum number of clocks per frame */
	int maxfclks;		/* Maximum number of clocks per frame */
	double sensv1[8];	/* The 8 preliminary sensor value readings */
	double sensv2[8];	/* The 8 final sensor value readings */
	double pows[9];		/* Power combinations of initial XYZ */
	double bxyz;		/* Biggest of XYZ */
	int j, k;

	if (p->debug) fprintf(stderr,"spyd2: about to get a reading\n");

	/* Compute number of frames for desired read time */
	nframes = (int)(RDTIME * p->rrate + 0.5);

	/* Establish the frame rate detect threshold level */
	clocks1 = (nframes * CLKRATE)/(10 * p->rrate);		/* Use 10% of measurement clocks */

	if ((ev = spyd2_GetMinMax(p, &clocks1, &min, &max)) != inst_ok)
		return ev; 

	/* Setup for measurement */
	thresh = (max - min)/5 + min;			/* Threshold is at 80% of max brightness */
	if (thresh == 0)
		thresh = 65535;						/* Set to max, otherwise reading will be 0 */
	frclocks = CLKRATE/p->rrate;
	minfclks = frclocks/3;					/* Allow for 180 Hz */
	maxfclks = (frclocks * 5)/2;			/* Allow for 24 Hz */

#ifdef DO_ADAPTIVE

	/* Do a reading at 10% of the standard intergration time */
	clocks1 = (nframes * CLKRATE)/(10 * p->rrate);		/* Use 10% of measurement clocks */

	if ((ev = spyd2_GetReading_ll(p, &clocks1, nframes/10, thresh, &minfclks, &maxfclks,
	                                 sensv1)) != inst_ok)
		return ev; 

	if (p->debug) {
		for (k = 1; k < 8; k++)
			printf("Sensor1 %d value = %f\n",k,sensv1[k]);
	}

	/* Convert sensor readings to XYZ value */
	bxyz = 0.0;
	for (j = 0; j < 3; j++) {
		XYZ[j] = p->cal_A[p->lcd][j][0];		/* First entry is a constant */
		for (k = 1; k < 8; k++)
			XYZ[j] += sensv1[k] * p->cal_A[p->lcd][j][k+1];
		if (XYZ[j] > bxyz)
			bxyz = XYZ[j];
	}
//printf("~1 Prelim Y = %f, bxyz = %f\n",XYZ[1],bxyz);

	if (bxyz < 5.0) {
		nframes *= 3;						/* Triple intergration time */
		if (p->debug) printf("Tripling integration time\n");
	} else if (bxyz < 10.0) {
		nframes *= 2;						/* Double intergration time */
		if (p->debug) printf("Doubling integration time\n");
	} else if (bxyz < 20.0) {
		nframes = (nframes * 3)/2;			/* 1.5 times intergration time */
		if (p->debug) printf("Extra 50%% integration time\n");
	}

#endif /* DO_ADAPTIVE */

	clocks2 = (nframes * CLKRATE)/p->rrate;		/* Normal intergration time */

	if ((ev = spyd2_GetReading_ll(p, &clocks2, nframes, thresh, &minfclks, &maxfclks,
	                                 sensv2)) != inst_ok)
		return ev; 

#ifdef DO_ADAPTIVE
	if (p->debug) {
		for (k = 1; k < 8; k++)
			printf("Sensor2 %d value = %f\n",k,sensv2[k]);
	}
#endif

	if (p->debug) {
		for (k = 1; k < 8; k++)
			printf("Sensor %d value = %f\n",k,sensv2[k]);
	}

	/* Convert sensor readings to XYZ values */
	for (j = 0; j < 3; j++) {
		XYZ[j] = p->cal_A[p->lcd][j][0];		/* First entry is a constant */

		for (k = 1; k < 8; k++) {		/* Skip first sensor reading since it isn't present */
			XYZ[j] += sensv2[k] * p->cal_A[p->lcd][j][k+1];
		}
	}
//printf("~1 real Y = %f\n",XYZ[1]);

#ifdef LEVEL2
	/* Add "level 2" correction factors */
	pows[0] = XYZ[0];
	pows[1] = XYZ[1];
	pows[2] = XYZ[2];
	pows[3] = XYZ[0] * XYZ[1];
	pows[4] = XYZ[0] * XYZ[2];
	pows[5] = XYZ[1] * XYZ[2];
	pows[6] = XYZ[0] * XYZ[0];
	pows[7] = XYZ[1] * XYZ[1];
	pows[8] = XYZ[2] * XYZ[2];

	if (p->debug) fprintf(stderr,"spyd2: got initial XYZ reading %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]);

	for (j = 0; j < 3; j++) {
		XYZ[j] = 0.0;

		for (k = 0; k < 9; k++) {
			XYZ[j] += pows[k] * p->cal_B[p->lcd][j][k];
		}
	}
#endif

	/* Guard against silliness */
	for (j = 0; j < 3; j++) {
		if (XYZ[j] < 0.0)
			XYZ[j] = 0.0;
	}

	if (p->debug) fprintf(stderr,"spyd2: got XYZ reading %f %f %f\n",XYZ[0], XYZ[1], XYZ[2]);

	return ev;
}

/* ------------------------------------------------------------ */

/* Read all the relevant register values */
static inst_code
spyd2_read_all_regs(
	spyd2 *p				/* Object */
) {
	inst_code ev;

	if (p->debug) fprintf(stderr,"spyd2: about to read all the EEProm values\n");

	/* HW version */
	if ((ev = spyd2_rd_ee_ushort(p, &p->hwver, 5)) != inst_ok)
		return ev;

	/* Serial number */
	if ((ev = spyd2_readEEProm(p, (unsigned char *)p->serno, 8, 8)) != inst_ok)
		return ev;
	p->serno[8] = '\000';

	/* CRT calibration values */
	if ((ev = spyd2_rdreg_3x9xfloat(p, p->cal_A[0][0], p->cal_A[0][1], p->cal_A[0][2], 16))
	                                                                            != inst_ok)
		return ev;
	if ((ev = spyd2_rdreg_3x9xfloat(p, p->cal_B[0][0], p->cal_B[0][1], p->cal_B[0][2], 128))
	                                                                            != inst_ok)
		return ev;

	/* LCD calibration values */
	if ((ev = spyd2_rdreg_3x9xfloat(p, p->cal_A[1][0], p->cal_A[1][1], p->cal_A[1][2], 256))
	                                                                            != inst_ok)
		return ev;
	if ((ev = spyd2_rdreg_3x9xfloat(p, p->cal_B[1][0], p->cal_B[1][1], p->cal_B[1][2], 384))
	                                                                            != inst_ok)
		return ev;

	/* Luminence only calibration values ??? */
	if ((ev = spyd2_rdreg_float(p, &p->cal_F[0], 240)) != inst_ok)
		return ev;
	if ((ev = spyd2_rdreg_float(p, &p->cal_F[1], 244)) != inst_ok)
		return ev;
	if ((ev = spyd2_rdreg_float(p, &p->cal_F[2], 248)) != inst_ok)
		return ev;
	if ((ev = spyd2_rdreg_float(p, &p->cal_F[3], 252)) != inst_ok)
		return ev;
	if ((ev = spyd2_rdreg_float(p, &p->cal_F[4], 364)) != inst_ok)
		return ev;
	if ((ev = spyd2_rdreg_float(p, &p->cal_F[5], 368)) != inst_ok)
		return ev;
	if ((ev = spyd2_rdreg_float(p, &p->cal_F[6], 372)) != inst_ok)
		return ev;

#ifdef NEVER
	{
		int i, j, k;

		
		DBG(("Cal_A:\n"))
		for (i = 0; i < 2;i++) {
			for (j = 0; j < 3; j++) {
				for (k = 0; k < 9; k++) {
					DBG(("Cal_A [%d][%d][%d] = %f\n",i,j,k,p->cal_A[i][j][k]))
				}
			}
		}
		DBG(("\nCal_B:\n"))
		for (i = 0; i < 2;i++) {
			for (j = 0; j < 3; j++) {
				for (k = 0; k < 9; k++) {
					DBG(("Cal_B [%d][%d][%d] = %f\n",i,j,k,p->cal_B[i][j][k]))
				}
			}
		}
		DBG(("\nCal_F:\n"))
		for (i = 0; i < 7;i++) {
			DBG(("Cal_F [%d] = %f\n",i,p->cal_F[i]))
		}
		DBG(("\n"))
	}
#endif/* NEVER */

	if (p->debug) fprintf(stderr,"spyd2: all EEProm read OK\n");

	return inst_ok;
}


/* Table to hold Spyder 2 Firmware, if it's installed */
unsigned int _spyder2_pld_size = 0;  /* Number of bytes to download */
unsigned int *spyder2_pld_size = &_spyder2_pld_size;
unsigned char *spyder2_pld_bytes = NULL;

/* Download the PLD if it is available, and check and check status */
static inst_code
spyd2_download_pld(
	spyd2 *p				/* Object */
) {
	inst_code ev;
	int stat;
	int i;

	if (p->debug) fprintf(stderr,"spyd2: about to download the PLD pattern\n");

	if (*spyder2_pld_size == 0 || *spyder2_pld_size == 0x11223344) {
		if (p->debug) fprintf(stderr,"spyd2: No PLD pattern available!\n");
		return spyd2_interp_code((inst *)p, SPYD2_NO_PLD_PATTERN) ;
	}
		
	for (i = 0; i < *spyder2_pld_size; i += 8) {
		if ((ev = spyd2_loadPLD(p, spyder2_pld_bytes + i, 8)) != inst_ok)
			return ev;
	}

	/* Check the status */
	if ((ev = spyd2_getstatus(p, &stat)) != inst_ok)
		return ev;

	if (stat != 0) {
		if (p->debug) fprintf(stderr,"spyd2: PLD download failed!\n");
		return spyd2_interp_code((inst *)p, SPYD2_PLDLOAD_FAILED);
	}

	/* Let the PLD initialize */
	msec_sleep(500);
		
	if (p->debug) fprintf(stderr,"spyd2: PLD pattern downloaded\n");

	return inst_ok;
}


/* ============================================================ */


/* Establish communications with a SPYD2 */
/* If it's a serial port, use the baud rate given, and timeout in to secs */
/* Return DTP_COMS_FAIL on failure to establish communications */
static inst_code
spyd2_init_coms(inst *pp, int port, baud_rate br, double tout) {
	spyd2 *p = (spyd2 *) pp;
	char buf[16];
	int rsize;
	long etime;
	int bi, i, rv;
	icomuflags usbflags = icomuf_none;
	inst_code ev = inst_ok;

	if (p->debug) {
		p->icom->debug = p->debug;	/* Turn on debugging */
		fprintf(stderr,"spyd2: About to init coms\n");
	}

	if (p->icom->is_usb_portno(p->icom, port) == instUnknown) {
		if (p->debug) fprintf(stderr,"spyd2: init_coms called to wrong device!\n");
			return spyd2_interp_code((inst *)p, SPYD2_UNKNOWN_MODEL);
	}

	if (p->debug) fprintf(stderr,"spyd2: About to init USB\n");

	/* On Linux the Spyder doesn't work reliably unless each */
	/* read is preceeded by a reset endpoint. */
#if defined(UNIX) && !defined(__APPLE__)		/* Linux*/
	usbflags |= icomuf_resetep_before_read;		/* The spyder USB is buggy ? */
#endif
#if defined(UNIX) && defined(__APPLE__)			/* OS X*/
	usbflags |= icomuf_reset_not_close;		/* The spyder USB is buggy ? */
#endif

	/* Set config, interface, write end point, read end point */
	/* ("serial" end points aren't used - the spyd2lay uses USB control messages) */
	p->icom->set_usb_port(p->icom, port, 1, 0x00, 0x00, usbflags); 

	if (p->debug) fprintf(stderr,"spyd2: init coms has suceeded\n");

	p->gotcoms = 1;
	return inst_ok;
}

/* Initialise the SPYD2 */
/* return non-zero on an error, with dtp error code */
static inst_code
spyd2_init_inst(inst *pp) {
	spyd2 *p = (spyd2 *)pp;
	inst_code ev = inst_ok;
	int stat;
	int i;

	if (p->debug) fprintf(stderr,"spyd2: About to init instrument\n");

	if (p->gotcoms == 0) /* Must establish coms before calling init */
		return spyd2_interp_code((inst *)p, SPYD2_NO_COMS);

	p->rrset = 0;
	p->rrate = DEFRRATE;
	for (i = 0; i < 8; i++)
		p->prevraw[i] = 0;			/* Internal counters will be reset */

	/* Reset the instrument */
	if ((ev = spyd2_reset(p)) != inst_ok)
		return ev;

	/* Fetch status until we get a status = 1 */
	for (i = 0; i < 50; i++) {
		if ((ev = spyd2_getstatus(p, &stat)) != inst_ok)
			return ev;

		if (stat == 1)
			break;
	}
	if (i >= 50)
		return spyd2_interp_code((inst *)p, SPYD2_BADSTATUS);

	/* Read the Serial EEProm contents */
	if ((ev = spyd2_read_all_regs(p)) != inst_ok)
		return ev;

	/* Download the PLD pattern and check the status */
	if ((ev = spyd2_download_pld(p)) != inst_ok)
		return ev;

	/* Do a dumy read */
	{
		int clocks = 500;
		int minfclks = 0;
		int maxfclks = 0;
		if ((ev = spyd2_GetReading_ll(p, &clocks, 100, 0, &minfclks, &maxfclks, NULL)) != inst_ok)
			return ev;
	}

	p->trig = inst_opt_trig_keyb;		/* default trigger mode */

	if (ev == inst_ok) {
		p->inited = 1;
		if (p->debug) fprintf(stderr,"spyd2: instrument inited OK\n");
	}

	if (p->verb) {
		printf("Instrument Type:   Spyder 2\n");
		printf("Serial Number:     %s\n",p->serno);
		printf("Hardwar version:   0x%04x\n",p->hwver);
	}

	return ev;
}

/* Read a single sample */
/* Return the dtp error code */
static inst_code
spyd2_read_sample(
inst *pp,
char *name,			/* Strip name (7 chars) */
ipatch *val) {		/* Pointer to instrument patch value */
	spyd2 *p = (spyd2 *)pp;
	int user_trig = 0;
	inst_code ev = inst_protocol_error;

	if (p->trig == inst_opt_trig_keyb) {
		int se;
		if ((se = icoms_poll_user(p->icom, 1)) != ICOM_TRIG) {
			/* Abort, term or command */
			return spyd2_interp_code((inst *)p, icoms2spyd2_err(se));
		}
		user_trig = 1;
		if (p->trig_return)
			printf("\n");
	}

	/* Attemp a CRT frame rate calibration if needed */
	if (p->lcd == 0 && p->rrset == 0) { 
		double refrate;

		if ((ev = spyd2_GetRefRate(p, &refrate)) != inst_ok)
			return ev; 
		if (refrate != 0.0) {
			p->rrate = refrate;
			p->rrset = 1;
		}
	}

	/* Read the XYZ value */
	if ((ev = spyd2_GetReading(p, val->aXYZ)) != inst_ok)
		return ev;

	val->XYZ_v = 0;
	val->aXYZ_v = 1;		/* These are absolute XYZ readings ? */
	val->Lab_v = 0;
	val->sp.spec_n = 0;

	if (user_trig)
		return inst_user_trig;
	return ev;
}

/* Determine if a calibration is needed. Returns inst_calt_none if not, */
/* inst_calt_unknown if it is unknown, or inst_calt_crt_freq */
/* if the display refresh rate has not bee determined, */
/* and we are in CRT mode */
inst_cal_type spyd2_needs_calibration(inst *pp) {
	spyd2 *p = (spyd2 *)pp;

	if (p->lcd == 0 && p->rrset == 0)
		return inst_calt_crt_freq;
	return inst_ok;
}

/* Request an instrument calibration. */
/* This is use if the user decides they want to do a calibration, */
/* in anticipation of a calibration (needs_calibration()) to avoid */
/* requiring one during measurement, or in response to measureing */
/* returning inst_needs_cal. Initially us an inst_cal_cond of inst_calc_none, */
/* and then be prepared to setup the right conditions, or ask the */
/* user to do so, each time the error inst_cal_setup is returned. */
inst_code spyd2_calibrate(
inst *pp,
inst_cal_type calt,		/* Calibration type. inst_calt_all for all neeeded */
inst_cal_cond *calc,	/* Current condition/desired condition */
char id[CALIDLEN]		/* Condition identifier (ie. white reference ID) */
) {
	spyd2 *p = (spyd2 *)pp;
	inst_code ev = inst_ok;

	id[0] = '\000';

	/* Translate default into what's needed or expected default */
	if (calt == inst_calt_all) {
		if (p->lcd == 0)
			calt = inst_calt_crt_freq;
	}

	if (calt == inst_calt_crt_freq && p->lcd == 0) {
		double refrate;

		if (*calc != inst_calc_disp_white) {
			*calc = inst_calc_disp_white;
			return inst_cal_setup;
		}

		/* Do CRT frame rate calibration */
		if ((ev = spyd2_GetRefRate(p, &refrate)) != inst_ok)
			return ev; 

		if (refrate == 0.0) {
			p->rrate = DEFRRATE;
		} else {
			p->rrate = refrate;
			p->rrset = 1;
		}

		return inst_ok;
	}

	return inst_unsupported;
}

/* Error codes interpretation */
static char *
spyd2_interp_error(inst *pp, int ec) {
//	spyd2 *p = (spyd2 *)pp;
	ec &= inst_imask;
	switch (ec) {
		case SPYD2_INTERNAL_ERROR:
			return "Non-specific software internal software error";
		case SPYD2_COMS_FAIL:
			return "Communications failure";
		case SPYD2_UNKNOWN_MODEL:
			return "Not a i1 Display";
		case SPYD2_DATA_PARSE_ERROR:
			return "Data from i1 Display didn't parse as expected";
		case SPYD2_USER_ABORT:
			return "User hit Abort key";
		case SPYD2_USER_TERM:
			return "User hit Terminate key";
		case SPYD2_USER_TRIG:
			return "User hit Trigger key";
		case SPYD2_USER_CMND:
			return "User hit a Command key";

		case SPYD2_OK:
			return "No error";

		case SPYD2_BADSTATUS:
			return "Too many retries waiting for status to come good";
		case SPYD2_PLDLOAD_FAILED:
			return "Wrong status after download of PLD";
		case SPYD2_BADREADSIZE:
			return "Didn't read expected amount of data";
		case SPYD2_TRIGTIMEOUT:
			return "Trigger timout";
		case SPYD2_OVERALLTIMEOUT:
			return "Overall timout";

		/* Internal errors */
		case SPYD2_BAD_EE_ADDRESS:
			return "Serial EEProm read is out of range 0 - 511";
		case SPYD2_BAD_EE_SIZE:
			return "Serial EEProm read size > 256";
		case SPYD2_NO_PLD_PATTERN:
			return "No PLD pattern is available";
		case SPYD2_NO_COMS:
			return "Communications hasn't been established";
		case SPYD2_NOT_INITED:
			return "Insrument hasn't been initialised";
		default:
			return "Unknown error code";
	}
}


/* Convert a machine specific error code into an abstract dtp code */
static inst_code 
spyd2_interp_code(inst *pp, int ec) {
//	spyd2 *p = (spyd2 *)pp;

	ec &= inst_imask;
	switch (ec) {

		case SPYD2_OK:
			return inst_ok;

		case SPYD2_INTERNAL_ERROR:
		case SPYD2_NO_COMS:
		case SPYD2_NOT_INITED:
		case SPYD2_BAD_EE_ADDRESS:
		case SPYD2_BAD_EE_SIZE:
		case SPYD2_NO_PLD_PATTERN:
			return inst_internal_error | ec;

		case SPYD2_COMS_FAIL:
		case SPYD2_BADREADSIZE:
		case SPYD2_TRIGTIMEOUT:
		case SPYD2_BADSTATUS:
		case SPYD2_OVERALLTIMEOUT:
			return inst_coms_fail | ec;

		case SPYD2_UNKNOWN_MODEL:
			return inst_unknown_model | ec;

//		return inst_protocol_error | ec;

		case SPYD2_USER_ABORT:
			return inst_user_abort | ec;
		case SPYD2_USER_TERM:
			return inst_user_term | ec;
		case SPYD2_USER_TRIG:
			return inst_user_trig | ec;
		case SPYD2_USER_CMND:
			return inst_user_cmnd | ec;

		case SPYD2_PLDLOAD_FAILED:
			return inst_hardware_fail | ec;
	}
	return inst_other_error | ec;
}

/* Destroy ourselves */
static void
spyd2_del(inst *pp) {
	spyd2 *p = (spyd2 *)pp;
	if (p->icom != NULL)
		p->icom->del(p->icom);
	free(p);
}

/* Return the instrument capabilities */
inst_capability spyd2_capabilities(inst *pp) {
	spyd2 *p = (spyd2 *)pp;
	inst_capability rv;

	rv = inst_emis_spot
	   | inst_emis_disp
	   | inst_emis_disp_crt
	   | inst_emis_disp_lcd
	   | inst_colorimeter
	     ;

	return rv;
}

/* Return the instrument capabilities 2 */
inst2_capability spyd2_capabilities2(inst *pp) {
	spyd2 *p = (spyd2 *)pp;
	inst2_capability rv;

	rv  = inst2_cal_crt_freq;
	rv |= inst2_prog_trig;
	rv |= inst2_keyb_trig;

	return rv;
}

/* Set device measurement mode */
inst_code spyd2_set_mode(inst *pp, inst_mode m)
{
	spyd2 *p = (spyd2 *)pp;
	inst_mode mm;		/* Measurement mode */

	/* The measurement mode portion of the mode */
	mm = m & inst_mode_measurement_mask;

	/* only display emission mode supported */
	if (mm != inst_mode_emis_spot
	 && mm != inst_mode_emis_disp) {
		return inst_unsupported;
	}

	/* Spectral mode is not supported */
	if (m & inst_mode_spectral)
		return inst_unsupported;

	p->mode = m;
	return inst_ok;
}

/* 
 * set or reset an optional mode
 * We assume that the instrument has been initialised.
 */
static inst_code
spyd2_set_opt_mode(inst *pp, inst_opt_mode m, ...)
{
	spyd2 *p = (spyd2 *)pp;
	inst_code ev = inst_ok;

	if (m == inst_opt_disp_crt) {
		if (p->lcd == 1)
			p->rrset = 0;		/* This is a hint we may have swapped displays */
		p->lcd = 0;
		return inst_ok;
	} else if (m == inst_opt_disp_lcd) {
		if (p->lcd != 1)
			p->rrset = 0;		/* This is a hint we may have swapped displays */
		p->lcd = 1;
		return inst_ok;

	}
	/* Record the trigger mode */
	if (m == inst_opt_trig_prog
	 || m == inst_opt_trig_keyb) {
		p->trig = m;
		return inst_ok;
	}
	if (m == inst_opt_trig_return) {
		p->trig_return = 1;
		return inst_ok;
	} else if (m == inst_opt_trig_no_return) {
		p->trig_return = 0;
		return inst_ok;
	}

	return inst_unsupported;
}

/* Constructor */
extern spyd2 *new_spyd2(icoms *icom, int debug, int verb)
{
	spyd2 *p;
	if ((p = (spyd2 *)calloc(sizeof(spyd2),1)) == NULL)
		error("spyd2: malloc failed!");

	if (icom == NULL)
		p->icom = new_icoms();
	else
		p->icom = icom;

	p->debug = debug;
	p->verb = verb;

	p->init_coms        = spyd2_init_coms;
	p->init_inst        = spyd2_init_inst;
	p->capabilities     = spyd2_capabilities;
	p->capabilities2    = spyd2_capabilities2;
	p->set_mode         = spyd2_set_mode;
	p->set_opt_mode     = spyd2_set_opt_mode;
	p->read_sample      = spyd2_read_sample;
	p->needs_calibration = spyd2_needs_calibration;
	p->calibrate        = spyd2_calibrate;
	p->interp_error     = spyd2_interp_error;
	p->del              = spyd2_del;

	p->itype = instUnknown;		/* Until initalisation */

	return p;
}

