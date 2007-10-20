
/* 
 * Argyll Color Correction System
 *
 * Windows NT serial I/O class
 *
 * Author: Graeme W. Gill
 * Date:   28/9/97
 *
 * Copyright 1997 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <conio.h>
#include "copyright.h"
#include "config.h"
#include "xspect.h"
#include "insttypes.h"
#include "icoms.h"
#include "usbio.h"

#undef DEBUG

#ifdef DEBUG
# define errout stderr
# define DBG(xx)	fprintf(errout, xx )
# define DBGF(xx)	fprintf xx
#else
# define DBG(xx)
# define DBGF(xx)
#endif


void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);

#ifdef __BORLANDC__
#define _kbhit kbhit
#endif

/* wait for and then return the next character from the keyboard */
int next_con_char(void) {
	return _getch();
}

/* If here is one, return the next character from the keyboard, else return 0 */
int poll_con_char(void) {
	if (_kbhit() != 0) {
		int c = _getch();
		return c;
	}
	return 0; 
}

/* Suck all characters from the keyboard */
void empty_con_chars(void) {
	Sleep(50);					/* _kbhit seems to have a bug */
	while (_kbhit()) {
		if (_getch() == 0x3)	/* ^C Safety */
			break;
	}
}

/* Sleep for the given number of seconds */
void sleep(unsigned int secs) {
	Sleep(secs * 1000);
}

/* Sleep for the given number of msec */
void msec_sleep(unsigned int msec) {
	Sleep(msec);
}

/* Return the current time in msec since */
/* the process started. */
unsigned int msec_time() {
	return GetTickCount();
}

static athread *beep_thread = NULL;
static int beep_delay;
static int beep_freq;
static int beep_msec;

/* Delayed beep handler */
static int delayed_beep(void *pp) {
	msec_sleep(beep_delay);
	Beep(beep_freq, beep_msec);
	return 0;
}

/* Activate the system beeper */
void msec_beep(int delay, int freq, int msec) {
	if (delay > 0) {
		if (beep_thread != NULL)
			beep_thread->del(beep_thread);
		beep_delay = delay;
		beep_freq = freq;
		beep_msec = msec;
		if ((beep_thread = new_athread(delayed_beep, NULL)) == NULL)
			error("Delayed beep failed to create thread");
	} else {
		Beep(freq, msec);
	}
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef NEVER    /* Not currently needed, or effective */

/* Set the current threads priority */
int set_interactive_priority() {
	if (SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_HIGHEST) == 0)
		return 1;
	return 0;
}

int set_normal_priority() {
	if (SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_NORMAL) == 0)
		return 1;
	return 0;
}

#endif /* NEVER */

/* - - - - - - - - - - - - - - - - - - - - - - - - */

/* Destroy the thread */
static void athread_del(
athread *p
) {
	DBG("athread_del called\n");

	if (p == NULL)
		return;

	if (p->th != NULL) {
		CloseHandle(p->th);
	}

	free(p);
}

/* _beginthread doesn't leak memory, but */
/* needs to be linked to a different library */
#ifdef NEVER
/* Thread function */
static void __cdecl threadproc(
	void *lpParameter
) {
#else
DWORD WINAPI threadproc(
	LPVOID lpParameter
) {
#endif
	athread *p = (athread *)lpParameter;

	p->result = p->function(p->context);
#ifdef NEVER
#else
	return 0;
#endif
}
 

athread *new_athread(
	int (*function)(void *context),
	void *context
) {
	athread *p = NULL;

	DBG("new_athread called\n");
	if ((p = (athread *)calloc(sizeof(athread), 1)) == NULL)
		return NULL;

	p->function = function;
	p->context = context;
	p->del = athread_del;

	/* Create a thread */
#ifdef NEVER
	p->th = _beginthread(threadproc, 0, (void *)p);
	if (p->th == -1) {
#else
	p->th = CreateThread(NULL, 0, threadproc, (void *)p, 0, NULL);
	if (p->th == NULL) {
#endif
		DBG("Failed to create thread\n");
		athread_del(p);
		return NULL;
	}

	DBG("About to exit new_athread()\n");
	return p;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */

/* Delete a file */
void delete_file(char *fname) {
	_unlink(fname);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */

/* Create and return a list of available serial ports or USB instruments for this system. */
/* We look at the registry key "HKEY_LOCAL_MACHINE\HARDWARE\DEVICEMAP\SERIALCOMM" */
/* to determine serial ports, and use libusb to discover USB instruments. */
/* Return NULL if function fails */
icompath **
icoms_get_paths(
icoms *p 
) {
	int i, j;
	int usbend = 0;
	LONG stat;
	HKEY sch;		/* Serial coms handle */

	/* Free any old list */
	if (p->paths != NULL) {
		for (i = 0; i < p->npaths; i++) {
			if (p->paths[i]->path != NULL)
				free(p->paths[i]->path);
			if (p->paths[i]->hev != NULL)
				hid_del_hid_device(p->paths[i]->hev);
			free(p->paths[i]);
		}
		free(p->paths);
		p->npaths = 0;
		p->paths = NULL;
	}
	
	usb_get_paths(p);
	hid_get_paths(p);
	usbend = p->npaths;

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
			if ((p->paths = (icompath **)calloc(sizeof(icompath *), 1 + 1)) == NULL)
				error("icoms: calloc failed!");
		} else {
			if ((p->paths = (icompath **)realloc(p->paths,
			                     sizeof(icompath *) * (p->npaths + 2))) == NULL)
				error("icoms: realloc failed!");
			p->paths[p->npaths+1] = NULL;
		}
		if ((p->paths[p->npaths] = malloc(sizeof(icompath))) == NULL)
			error("icoms: malloc failed!");
		if ((p->paths[p->npaths]->path = strdup(value)) == NULL)
			error("icoms: strdup failed!");
#ifdef ENABLE_USB
		p->paths[p->npaths]->dev = NULL;
		p->paths[p->npaths]->hev = NULL;
#endif /* ENABLE_USB */
		p->npaths++;
		p->paths[p->npaths] = NULL;
	}
	if ((stat = RegCloseKey(sch)) != ERROR_SUCCESS)
		warning("RegCloseKey failed with %d\n",stat);

	/* Sort the COM keys so people don't get confused... */
	for (i = usbend; i < (p->npaths-1); i++) {
		for (j = i+1; j < p->npaths; j++) {
			if (strcmp(p->paths[i]->path, p->paths[j]->path) > 0) {
				icompath *tt = p->paths[i];
				p->paths[i] = p->paths[j];
				p->paths[j] = tt;
			}
		}
	}

	return p->paths;
}


/* Close the port */
void icoms_close_port(icoms *p) {
	if (p->is_open) {
		if (p->is_usb) {
			usb_close_port(p);
		} else if (p->is_hid) {
			hid_close_port(p);
		} else {
			if (p->phandle != NULL)
				CloseHandle(p->phandle);
			p->is_open = 0;
			if (p->ppath != NULL) {
				if (p->ppath->path != NULL)
					free(p->ppath->path);
				free(p->ppath);
				p->ppath = NULL;
			}
		}
	}
}

static int icoms_usb_ser_write(icoms *p, char *wbuf, double tout);
static int icoms_usb_ser_read(icoms *p, char *rbuf, int bsize, char tc, int ntc, double tout);
static int icoms_ser_write(icoms *p, char *wbuf, double tout);
static int icoms_ser_read(icoms *p, char *rbuf, int bsize, char tc, int ntc, double tout);

/* Set the serial port number and characteristics */
static void
icoms_set_ser_port(
icoms *p, 
int 		 port,		/* Serial com port, 1 - N, 0 for no change. */
flow_control fc,
baud_rate	 baud,
parity		 parity,
stop_bits	 stop,
word_length	 word)
{

	if (p->debug) fprintf(stderr,"icoms: About to set serial port characteristics\n");

	if (port >= 1) {
		if (p->is_open && port != p->port) {	/* If port number changes */
			icoms_close_port(p);
		}
	}

	if (p->is_usb_portno(p, port) == instUnknown) {
		DCB dcb;

		if (fc != fc_nc)
			p->fc = fc;
		if (baud != baud_nc)
			p->baud = baud;
		if (parity != parity_nc)
			p->parity = parity;
		if (stop != stop_nc)
			p->stop_bits = stop;
		if (word != length_nc)
			p->word_length = word;

		/* Make sure the port is open */
		if (!p->is_open) {
			if (p->ppath != NULL) {
				if (p->ppath->path != NULL)
					free(p->ppath->path);
				free(p->ppath);
				p->ppath = NULL;
			}

			if (p->paths == NULL)
				icoms_get_paths(p);

			if (port <= 0 || port > p->npaths)
				error("icoms - set_ser_port: port number out of range!");

			if (p->paths[port-1]->dev != NULL)
				error("icoms - set_ser_port: USB port sent to set_serial!");

			if ((p->ppath = malloc(sizeof(icompath))) == NULL)
				error("malloc() failed on com port path");
			*p->ppath = *p->paths[port-1];				/* Structure copy */
			if ((p->ppath->path = strdup(p->paths[port-1]->path)) == NULL)
				error("strdup() failed on com port path");
			p->port = port;

			if (p->debug) fprintf(stderr,"icoms: About to open port '%s'\n",p->ppath->path);

			if ((p->phandle = CreateFile(p->ppath->path,
				GENERIC_READ|GENERIC_WRITE,
				0,				/* Exclusive access */
				NULL,			/* No security */
				OPEN_EXISTING,	/* Does not make sense to create */
				0,				/* No overlapped I/O */
				NULL)			/* NULL template */
			) == INVALID_HANDLE_VALUE)
				error("Opening COM port '%s' failed",p->ppath->path);
			p->is_open = 1;
		}

		if (GetCommState(p->phandle, &dcb) == FALSE)
			error("Reading Comm State failed");

		/* Set misc stuff to default */
		dcb.fBinary = TRUE; 						/* Binary mode is only one that works */
		dcb.fOutxCtsFlow = FALSE;					/* Not using Cts flow control */
		dcb.fOutxDsrFlow = FALSE;					/* Not using Dsr Flow control */
		dcb.fDtrControl = DTR_CONTROL_ENABLE;		/* Enable DTR during connection */
		dcb.fDsrSensitivity = FALSE;				/* Not using Dsr Flow control */
		dcb.fTXContinueOnXoff = TRUE;				/* */
		dcb.fOutX = FALSE;							/* No default Xon/Xoff flow control */
		dcb.fInX = FALSE;							/* No default Xon/Xoff flow control */
		dcb.fErrorChar = FALSE;
		dcb.fNull = FALSE;
		dcb.fRtsControl = RTS_CONTROL_ENABLE;		/* Turn RTS on during connection */
		dcb.fAbortOnError = TRUE;

		switch (p->fc) {
			case fc_nc:
				error("icoms - set_ser_port: illegal flow control!");
				break;
			case fc_XonXOff:
				/* Use Xon/Xoff bi-directional flow control */
				dcb.fOutX = TRUE;
				dcb.fInX = TRUE;
				dcb.XonChar = 0x11;					/* ^Q */
				dcb.XoffChar = 0x13;				/* ^S */
				// dcb.fTXContinueOnXoff = FALSE;		/* Stop transmitting if input is full */
				break;
			case fc_Hardware:
				/* Use RTS/CTS bi-directional flow control */
				dcb.fOutxCtsFlow = TRUE;
				dcb.fRtsControl = RTS_CONTROL_HANDSHAKE; 
				break;
			default:
				break;
		}

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
				error("nt icoms - set_ser_port: illegal baud rate! (0x%x)",p->baud);
				break;
				
		}

		switch (p->parity) {
			case parity_nc:
				error("icoms - set_ser_port: illegal parity setting!");
				break;
			case parity_none:
				dcb.fParity = FALSE;
				dcb.Parity = NOPARITY;
				break;
			case parity_odd:
				dcb.fParity = TRUE;
				dcb.Parity = ODDPARITY;
				break;
			case parity_even:
				dcb.fParity = TRUE;
				dcb.Parity = EVENPARITY;
				break;
			default:
				break;
		}

		switch (p->stop_bits) {
			case stop_nc:
				error("icoms - set_ser_port: illegal stop bits!");
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
				error("icoms - set_ser_port: illegal word length!");
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

		PurgeComm(p->phandle, PURGE_TXCLEAR | PURGE_RXCLEAR);

		if (!SetCommState(p->phandle, &dcb))
			error("SetCommState failed with LastError %d", GetLastError());

		PurgeComm(p->phandle, PURGE_TXCLEAR | PURGE_RXCLEAR);

		p->write = icoms_ser_write;
		p->read = icoms_ser_read;

	}
	if (p->debug) fprintf(stderr,"icoms: serial port characteristics set ok\n");
}

/* ---------------------------------------------------------------------------------*/
/* Serial write/read */

/* Write the characters in the buffer out */
/* Data will be written up to the terminating nul */
/* Return relevant error status bits */
static int
icoms_ser_write(
icoms *p,
char *wbuf,
double tout)
{
	COMMTIMEOUTS tmo;
	DWORD wbytes;
	int c, len;
	long toc, i, top;		/* Timout count, counter, timeout period */

	if (p->debug) fprintf(stderr,"icoms: About to write '%s' ",icoms_fix(wbuf));

	if (p->phandle == NULL)
		error("icoms_write: not initialised");

	len = strlen(wbuf);
	tout *= 1000.0;		/* Timout in msec */
	p->lerr = 0;

	top = 100;						/* Timeout period in msecs */
	toc = (int)(tout/top + 0.5);	/* Number of timout periods in timeout */
	if (toc < 1)
		toc = 1;

	/* Set the timout value */
	tmo.ReadIntervalTimeout = top;
	tmo.ReadTotalTimeoutMultiplier = 0;
	tmo.ReadTotalTimeoutConstant = top;
	tmo.WriteTotalTimeoutMultiplier = top;
	tmo.WriteTotalTimeoutConstant = 0;
	if (!SetCommTimeouts(p->phandle, &tmo))
		error("SetCommTimouts failed");

	/* Until data is all written, we time out, or the user aborts */
	for(i = toc; i > 0 && len > 0;) {
		if (!WriteFile(p->phandle, wbuf, len, &wbytes, NULL)) {
			DWORD errs;
			if (!ClearCommError(p->phandle,&errs,NULL))
				error("Write to COM port failed, and Clear error failed");
			if (errs & CE_BREAK)
				p->lerr |= ICOM_BRK; 
			if (errs & CE_FRAME)
				p->lerr |= ICOM_FER; 
			if (errs & CE_RXPARITY)
				p->lerr |= ICOM_PER; 
			if (errs & CE_RXOVER)
				p->lerr |= ICOM_OER; 
			break;
		} else if (wbytes == 0) { 
			i--;			/* Timeout */
		} else if (wbytes > 0) { /* Account for bytes done */
			i = toc;
			len -= wbytes;
			wbuf += len;
		}

		/* Check for user abort/term/command */
		if (_kbhit() && p->uih[c = _getch()] != ICOM_OK) {
			p->cut = c;
			p->lerr = p->uih[c];
			if (p->uih[c] == ICOM_TERM
			 || p->uih[c] == ICOM_CMND)
				break;
		}
	}
	if (i <= 0) {		/* Timed out */
		p->lerr |= ICOM_TO;
	}
	if (p->debug) fprintf(stderr,"ICOM err 0x%x\n",p->lerr);
	return p->lerr;
}


/* Read characters into the buffer */
/* Return string will be terminated with a nul */
static int
icoms_ser_read(
icoms *p,
char *rbuf,			/* Buffer to store characters read */
int bsize,			/* Buffer size */
char tc,			/* Terminating characer */
int ntc,			/* Number of terminating characters */
double tout)		/* Time out in seconds */
{
	COMMTIMEOUTS tmo;
	DWORD rbytes;
	int c, j;
	long toc, i, top;		/* Timout count, counter, timeout period */
	char *rrbuf = rbuf;		/* Start of return buffer */
	DCB dcb;

	if (p->phandle == NULL)
		error("icoms_read: not initialised");

	if (bsize < 3)
		error("icoms_read given too small a buffer");

	if (GetCommState(p->phandle, &dcb) == FALSE)
		error("Reading Comm State failed");

	dcb.EofChar = tc;

	if (!SetCommState(p->phandle, &dcb))
		error("SetCommState failed");
	p->tc = tc;

	p->lerr = 0;
	tout *= 1000.0;		/* Timout in msec */
	bsize--;	/* Allow space for null */

	top = 100;						/* Timeout period in msecs */
	toc = (int)(tout/top + 0.5);	/* Number of timout periods in timeout */
	if (toc < 1)
		toc = 1;

	/* Set the timout value */
	tmo.ReadIntervalTimeout = top;
	tmo.ReadTotalTimeoutMultiplier = 0;
	tmo.ReadTotalTimeoutConstant = top;
	tmo.WriteTotalTimeoutMultiplier = top;
	tmo.WriteTotalTimeoutConstant = 0;
	if (!SetCommTimeouts(p->phandle, &tmo))
		error("SetCommTimouts failed");

	/* Until data is all read, we time out, or the user aborts */
	for (i = toc, j = 0; i > 0 && bsize > 1 && j < ntc ;) {
		if (!ReadFile(p->phandle, rbuf, bsize, &rbytes, NULL)) {
			DWORD errs;
			if (!ClearCommError(p->phandle,&errs,NULL))
				error("Read from COM port failed, and Clear error failed");
			if (errs & CE_BREAK)
				p->lerr |= ICOM_BRK; 
			if (errs & CE_FRAME)
				p->lerr |= ICOM_FER; 
			if (errs & CE_RXPARITY)
				p->lerr |= ICOM_PER; 
			if (errs & CE_RXOVER)
				p->lerr |= ICOM_OER; 
			break;
		} else if (rbytes == 0) { 
			i--;			/* Timeout */
		} else if (rbytes > 0) { /* Account for bytes done */
			i = toc;
			bsize -= rbytes;
			while(rbytes--) {	/* Count termination characters */
				if (*rbuf++ == tc)
					j++;
			}
		}
		/* Check for user abort/term/command */
		if (_kbhit() && p->uih[c = _getch()] != ICOM_OK) {
			p->cut = c;
			p->lerr = p->uih[c];
			if (p->uih[c] == ICOM_USER
			 || p->uih[c] == ICOM_TERM
			 || p->uih[c] == ICOM_TRIG
			 || p->uih[c] == ICOM_CMND)
				break;
		}
	}
	if (i <= 0) {			/* timed out */
		p->lerr |= ICOM_TO;
	}
	*rbuf = '\000';
	if (p->debug) fprintf(stderr,"icoms: About to return read '%s' ICOM err 0x%x\n",icoms_fix(rrbuf),p->lerr);
	return p->lerr;
}

/* ---------------------------------------------------------------------------------*/


/* Destroy ourselves */
static void
icoms_del(icoms *p) {
	if (p->debug) fprintf(stderr,"icoms: delete called\n");
	if (p->is_open) {
		if (p->debug) fprintf(stderr,"icoms: closing port\n");
		icoms_close_port(p);
	}
	if (p->paths != NULL) {
		int i;
		for (i = 0; i < p->npaths; i++) {
			if (p->paths[i]->path != NULL)
				free(p->paths[i]->path);
			if (p->paths[i]->hev != NULL)
				hid_del_hid_device(p->paths[i]->hev);
			free(p->paths[i]);
		}
		free(p->paths);
	}
	if (p->ppath != NULL) {
		if (p->ppath->path != NULL)
			free(p->ppath->path);
		free(p->ppath);
	}
	free (p);
}

/* Constructor */
extern icoms *new_icoms()
{
	icoms *p;
	if ((p = (icoms *)calloc(sizeof(icoms), 1)) == NULL)
		error("icoms: calloc failed!");

	/* Init things to null values */
	p->lerr = 0;
	p->phandle = NULL;
	p->ppath = NULL;
	p->port = -1;
	p->fc = fc_nc;
	p->baud = baud_nc;
	p->parity = parity_nc;
	p->stop_bits = stop_nc;
	p->word_length = length_nc;
	p->tc = -1;
	p->debug = 0;
	
	p->get_paths = icoms_get_paths;
	p->set_ser_port = icoms_set_ser_port;

	p->write = NULL;
	p->read = NULL;

	p->del = icoms_del;

	usb_set_usb_methods(p);
	hid_set_hid_methods(p);

	return p;
}

