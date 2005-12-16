/* 
 * Argyll Color Correction System
 *
 * Windows NT serial I/O class
 *
 * Author: Graeme W. Gill
 * Date:   28/9/97
 *
 * Copyright 1997 - 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <conio.h>
#include "serio.h"

#include <stdarg.h>
void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);

/* return the next character from the keyboard */
int next_con_char(void) {
	return _getch();
}

/* Suck all characters from the keyboard */
void empty_con_chars(void) {
	Sleep(50);					/* _kbhit seems to have a bug */
	while (_kbhit()) {
		if (_getch() == 0x3)	/* ^C Safety */
			break;
	}
}

/* Convert control chars to ^[A-Z] notation in a string */
static char *
fix(char *s) {
	static char buf [200];
	char *d;
	for(d = buf; ;) {
		if (*s < ' ' && *s > '\000') {
			*d++ = '^';
			*d++ = *s++ + '@';
		} else
		*d++ = *s++;
		if (s[-1] == '\000')
			break;
	}
	return buf;
}

/* Create and return a list of available serial ports for this system. */
/* We look at the registry key "HKEY_LOCAL_MACHINE\HARDWARE\DEVICEMAP\SERIALCOMM" */
/* to determine this. Return NULL if function fails */
static char **
serio_get_paths(
struct _serio *p 
) {
	int i;
	LONG stat;
	HKEY sch;		/* Serial coms handle */

	/* Free any old list */
	if (p->paths != NULL) {
		for (i = 0; i < p->npaths; i++)
			free(p->paths[i]);
		free(p->paths);
		p->npaths = 0;
		p->paths = NULL;
	}

	/* Look in the registry */
	if ((stat = RegOpenKeyEx(HKEY_LOCAL_MACHINE, "HARDWARE\\DEVICEMAP\\SERIALCOMM",
	                         0, KEY_READ, &sch)) != ERROR_SUCCESS) {
		warning("RegOpenKeyEx failed with %d",stat);
		return NULL;
	}

	/* Look at all the values in this key */
	for (i = 0; ; i++) {
		char valname[100];
		DWORD vnsize = 100;
		DWORD vtype;
		char value[100];
		DWORD vsize = 100;

		stat = RegEnumValue(
			sch,		/* handle to key to enumerate */
			i,			/* index of subkey to enumerate */
			valname,    /* address of buffer for value name */
			&vnsize,	/* address for size of value name buffer */
			NULL,		/* reserved */
			&vtype,		/* Address of value type */
			value,		/* Address of value buffer */
			&vsize		/* Address of value buffer size */
		);
		if (stat == ERROR_NO_MORE_ITEMS) {
			break;
		}
		if (stat != ERROR_SUCCESS) {
			warning("RegEnumValue failed with %d",stat);
			break;
		}
		valname[100-1] = '\000';
		value[100-1] = '\000';

		if (vtype != REG_SZ) {
			warning("RegEnumValue didn't return stringz type");
			continue;
		}

		/* Add the port to the list */
		if (p->paths == NULL) {
			if ((p->paths = (char **)calloc(sizeof(char *), 1 + 1)) == NULL)
				error("serio: calloc failed!");
		} else {
			if ((p->paths = (char **)realloc(p->paths,
			                     sizeof(char *) * (p->npaths + 2))) == NULL)
				error("serio: realloc failed!");
			p->paths[p->npaths+1] = NULL;
		}
		if ((p->paths[p->npaths] = strdup(value)) == NULL)
			error("serio: strdup failed!");
		p->npaths++;

	}
	if ((stat = RegCloseKey(sch)) != ERROR_SUCCESS)
		warning("RegCloseKey failed with %d\n",stat);

	return p->paths;
}

/* Set the port number and characteristics */
static void
serio_set_port(
struct _serio *p, 
char       *ppath,		/* Path to serial device, NULL if use port */
int 		port,		/* com port, 1 - 4 */
baud_rate	baud,
parity		parity,
stop_bits	stop,
word_length	word)
{
	DCB dcb;
	COMMTIMEOUTS tmo;
	long rwait = (int)(0.1 * 1000.0 + 0.5);		/* Timout every 100msec secs */

	if (p->debug) fprintf(stderr,"serio: About to set port characteristics\n");

	if (ppath != NULL || (port >= 1)) {
		if (p->phandle != NULL
		 && ((ppath != NULL && strcmp(ppath, p->ppath) != 0)
		 ||  (ppath == NULL && port != p->port))) {	/* If port number changes */
			CloseHandle(p->phandle);
			p->phandle = NULL;
		}
	}
	if (baud != baud_nc)
		p->baud = baud;
	if (parity != parity_nc)
		p->parity = parity;
	if (stop != stop_nc)
		p->stop_bits = stop;
	if (word != length_nc)
		p->word_length = word;

	/* Make sure the port is open */
	if (p->phandle == NULL) {
		if (p->ppath != NULL)
			free(p->ppath);

		if (p->ppath != NULL) {	/* Use given device path */
			p->ppath = strdup(ppath);
			p->port = -1;

		} else {		/* Choose from ports */
			if (p->paths == NULL)
				serio_get_paths(p);
	
			if (port <= 0 || port > p->npaths)
				error("serio - set_port: port number out of range!");

			p->ppath = strdup(p->paths[port-1]);
			p->port = port;
		}
		if (p->ppath == NULL)
			error("strdup() failed on com port path");

		if (p->debug) fprintf(stderr,"serio: About to open port '%s'\n",p->ppath);

		if ((p->phandle = CreateFile(p->ppath,
			GENERIC_READ|GENERIC_WRITE,
			0,				/* Exclusive access */
			NULL,			/* No security */
			OPEN_EXISTING,	/* Does not make sense to create */
			0,				/* No overlapped I/O */
			NULL)			/* NULL template */
		) == INVALID_HANDLE_VALUE)
			error("Opening COM port '%s' failed",p->ppath);
	}

	if (GetCommState(p->phandle, &dcb) == FALSE)
		error("Reading Comm State failed");

	/* Set misc stuff */
	dcb.fBinary = TRUE; 	/* We need binary mode */
	dcb.fOutxCtsFlow = FALSE;
	dcb.fOutxDsrFlow = FALSE;
	dcb.fDtrControl = DTR_CONTROL_ENABLE;
	dcb.fDsrSensitivity = FALSE;
	dcb.fTXContinueOnXoff = TRUE;
	dcb.fOutX = FALSE;
	dcb.fInX = FALSE;
	dcb.fErrorChar = FALSE;
	dcb.fNull = FALSE;
	dcb.fRtsControl = RTS_CONTROL_ENABLE;
	dcb.fAbortOnError = TRUE;

	switch (p->baud) {
		case baud_110:
			dcb.BaudRate = CBR_110;
			break;
		case baud_300:
			dcb.BaudRate = CBR_300;
			break;
		case baud_600:
			dcb.BaudRate = CBR_600;
			break;
		case baud_1200:
			dcb.BaudRate = CBR_1200;
			break;
		case baud_2400:
			dcb.BaudRate = CBR_2400;
			break;
		case baud_4800:
			dcb.BaudRate = CBR_4800;
			break;
		case baud_9600:
			dcb.BaudRate = CBR_9600;
			break;
		case baud_14400:
			dcb.BaudRate = CBR_14400;
			break;
		case baud_19200:
			dcb.BaudRate = CBR_19200;
			break;
		case baud_38400:
			dcb.BaudRate = CBR_38400;
			break;
		case baud_57600:
			dcb.BaudRate = CBR_57600;
			break;
		case baud_115200:
			dcb.BaudRate = CBR_115200;
			break;
		default:
			error("nt serio - set_port: illegal baud rate! (0x%x)",p->baud);
			break;
			
	}

	switch (p->parity) {
		case parity_nc:
			error("serio - set_port: illegal parity setting!");
			break;
		case parity_none:
			dcb.fParity = FALSE;
			dcb.Parity = NOPARITY;
			break;
		case parity_odd:
			dcb.fParity = TRUE;
			dcb.Parity = ODDPARITY;
			break;
		caseparity_even:
			dcb.fParity = TRUE;
			dcb.Parity = EVENPARITY;
			break;
	}

	switch (p->stop_bits) {
		case stop_nc:
			error("serio - set_port: illegal stop bits!");
			break;
		case stop_1:
			dcb.StopBits = ONESTOPBIT;
			break;
		case stop_2:
			dcb.StopBits = TWOSTOPBITS;
			break;
	}

	switch (p->word_length) {
		case length_nc:
			error("serio - set_port: illegal word length!");
		case length_5:
			dcb.ByteSize = 5;
			break;
		case length_6:
			dcb.ByteSize = 6;
			break;
		case length_7:
			dcb.ByteSize = 7;
			break;
		case length_8:
			dcb.ByteSize = 8;
			break;
	}

	if (!SetCommState(p->phandle, &dcb))
		error("SetCommState failed");

	/* Set the timout value */
	tmo.ReadIntervalTimeout = rwait;
	tmo.ReadTotalTimeoutMultiplier = 0;
	tmo.ReadTotalTimeoutConstant = rwait;
	tmo.WriteTotalTimeoutMultiplier = rwait;
	tmo.WriteTotalTimeoutConstant = 0;
	if (!SetCommTimeouts(p->phandle, &tmo))
		error("SetCommTimouts failed");

	if (p->debug) fprintf(stderr,"serio: port characteristics set ok\n");
}

/* Write the characters in the buffer out */
/* Data will be written up to the terminating nul */
/* Return relevant error status bits */
static int
serio_write(
struct _serio *p,
char *wbuf,
double tout)
{
	DWORD wbytes;
	int len;
	long i, rwait = (int)(tout * 10.0 + 0.5);	/* Number of 100msec intervals */

	if (rwait < 1)
		rwait = 1;
	if (p->debug) fprintf(stderr,"serio: About to write '%s'\n",fix(wbuf));
	if (p->phandle == NULL)
		error("serio_write: not initialised");

	p->lerr = 0;
	i = rwait;
	len = strlen(wbuf);

	/* Until timed out, aborted, or transmitted */
	while (i-- > 0 && len > 0) {
		if (!WriteFile(p->phandle, wbuf, len, &wbytes, NULL)) {
			DWORD errs;
			if (!ClearCommError(p->phandle,&errs,NULL))
				error("Write to COM port failed, and Clear error failed");
			if (errs & CE_BREAK)
				p->lerr |= SIO_BRK; 
			if (errs & CE_FRAME)
				p->lerr |= SIO_FER; 
			if (errs & CE_RXPARITY)
				p->lerr |= SIO_PER; 
			if (errs & CE_RXOVER)
				p->lerr |= SIO_OER; 
			break;
		}
		if (len > 0) { /* Account for bytes done */
			len -= wbytes;
			wbuf += len;
			i = rwait;
		}

		/* Check for user abort (This should be a callback!) */
		if (_kbhit() && _getch() == 0x1b) {	/* User hit escape */
			p->lerr |= SIO_USER;
			break;
		}
	}
	if (i <= 0) {		/* Timed out */
		p->lerr |= SIO_TO;
	}
	return p->lerr;
}

/* Read characters into the buffer */
/* Return string will be terminated with a nul */
static int
serio_read(struct _serio *p,
char *rbuf,			/* Buffer to store characters read */
int bsize,			/* Buffer size */
char tc,			/* Terminating characer */
int ntc,			/* Number of terminating characters */
double tout)		/* Time out in seconds */
{
	DWORD rbytes;
	int j = 0;
	long i, rwait = (int)(tout * 10.0 + 0.5);	/* number of 100msec intervals */
	char *rrbuf = rbuf;		/* Start of return buffer */

	if (rwait < 1)
		rwait = 1;

	if (p->phandle == NULL)
		error("serio_read: not initialised");

	if (bsize < 3)
		error("serio_read given too small a buffer");

	if (tc != p->tc) {	/* Set the termination char */
		DCB dcb;
		if (GetCommState(p->phandle, &dcb) == FALSE)
			error("Reading Comm State failed");

		dcb.EofChar = tc;

		if (!SetCommState(p->phandle, &dcb))
			error("SetCommState failed");
		p->tc = tc;
	}
	p->lerr = 0;
	i = rwait;
	bsize--;	/* Allow space for null */
	/* Until timed out, aborted, or we get what we want */
	while (i-- > 0 && bsize > 1 && j < ntc) {
		if (!ReadFile(p->phandle, rbuf, bsize, &rbytes, NULL)) {
			DWORD errs;
			if (!ClearCommError(p->phandle,&errs,NULL))
				error("Write to COM port failed, and Clear error failed");
			if (errs & CE_BREAK)
				p->lerr |= SIO_BRK; 
			if (errs & CE_FRAME)
				p->lerr |= SIO_FER; 
			if (errs & CE_RXPARITY)
				p->lerr |= SIO_PER; 
			if (errs & CE_RXOVER)
				p->lerr |= SIO_OER; 
			break;
		}

		/* Account for bytes done */
		if (rbytes > 0) { /* Account for bytes done */
			bsize -= rbytes;
			i = rwait;

			while(rbytes--) {	/* Count termination characters */
				if (*rbuf++ == tc)
					j++;
			}
		}
		/* Check for user abort (This should be a callback!) */
		if (_kbhit() && _getch() == 0x1b) {
			p->lerr |= SIO_USER;
			break;
		}
	}
	*rbuf = '\000';
	if (i <= 0) {			/* timed out */
		p->lerr |= SIO_TO;
	}
	if (p->debug) fprintf(stderr,"serio: About to return read '%s' lerr 0x%x\n",fix(rrbuf),p->lerr);
	return p->lerr;
}

/* write and read */
static int
serio_write_read(
struct _serio *p,
char *wbuf,			/* Write puffer */
char *rbuf,			/* Read buffer */
int bsize,			/* Buffer size */
char tc,			/* Terminating characer */
int ntc,			/* Number of terminating characters */
double tout)
{
	int debug;

	/* Flush any stray chars */
	for (p->lerr = 0;;) {
		int debug = p->debug; p->debug = 0;
		serio_read(p, rbuf, bsize, '\000', 100000, 0.1);
		p->debug = debug;
		if (p->lerr != 0)
			break;				/* Expect timeout with nothing to read */
	}
	p->lerr = 0;

	/* Write the write data */
	serio_write(p, wbuf, tout);
	if (p->lerr != 0)
		return p->lerr;

	/* Read response */
	serio_read(p, rbuf, bsize, tc, ntc, tout);
	return p->lerr;
}

/* Destroy ourselves */
static void
serio_del(serio *p) {
	if (p->paths != NULL) {
		int i;
		for (i = 0; i < p->npaths; i++)
			free(p->paths[i]);
		free(p->paths);
	}
	if (p->phandle != NULL)
		CloseHandle(p->phandle);
	if (p->ppath != NULL)
		free(p->ppath);
	free (p);
}

/* Constructor */
extern serio *new_serio()
{
	serio *p;
	if ((p = (serio *)calloc(sizeof(serio), 1)) == NULL)
		error("serio: calloc failed!");

	/* Init things to null values */
	p->pollcalib = 0;
	p->lerr = 0;
	p->phandle = NULL;
	p->ppath = NULL;
	p->port = -1;
	p->baud = baud_nc;
	p->parity = parity_nc;
	p->stop_bits = stop_nc;
	p->word_length = length_nc;
	p->tc = -1;
	p->debug = 0;
	
	p->get_paths = serio_get_paths;
	p->set_port = serio_set_port;
	p->write = serio_write;
	p->read = serio_read;
	p->write_read = serio_write_read;
	p->del = serio_del;

	return p;
}

