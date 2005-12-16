#ifndef SERIO_H
/* 
 * Argyll Color Correction System
 *
 * An abstracted PC serial I/O class
 *
 * Author: Graeme W. Gill
 * Date:   5/10/96
 *
 * Copyright 1996, 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* !!!! Should add support for XON/XOFF handshaking ? */
/* !!!! At the moment handshaking is disabled for Windows & Unix */
/* !!!! We're assuming nothing blocks ? */
/* (!!!! default SS has this mode ?) */

/* baud rate available on all systems */
typedef enum {
	baud_nc      = 0,		/* Not Configured */
	baud_110     = 1,
	baud_300     = 2,
	baud_600     = 3,
	baud_1200    = 4,
	baud_2400    = 5,
	baud_4800    = 6,
	baud_9600    = 7,
	baud_14400   = 8,
	baud_19200   = 9,
	baud_38400   = 10,
	baud_57600   = 11,
	baud_115200  = 12
} baud_rate;

/* Possible parity */
typedef enum {
	parity_nc,
	parity_none,
	parity_odd,
	parity_even
} parity;

/* Possible stop bits */
typedef enum {
	stop_nc,
	stop_1,
	stop_2
} stop_bits;

/* Possible word length */
typedef enum {
	length_nc,
	length_5,
	length_6,
	length_7,
	length_8
} word_length;

/* Status bits/ return values */
#define SIO_TO	0x80		/* Timed out */
#define SIO_XRE	0x40		/* Xmit shift reg empty */
#define SIO_XHE	0x20		/* Xmit hold reg empty */
#define SIO_BRK	0x10		/* Break detected */
#define SIO_FER	0x08		/* Framing error */
#define SIO_PER	0x04		/* Parity error */
#define SIO_OER	0x02		/* Overun error */
#define SIO_DRY	0x01		/* Recv data ready */

#define SIO_USER 0x100		/* User abort */

#if defined (NT)
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

struct _serio {
  /* Private: */
	char *ppath;				/* path to port */
	int port;					/* com port number: 1, 2, 3, 4, -1 == defined by path only. */
#if defined (MSDOS)
	unsigned short porta;		/* Hardware port address */
#endif
#if defined (NT)
	HANDLE phandle;				/* NT handle */
#endif
#if defined (UNIX) || defined(__APPLE__)
	int fd;						/* Unix file descriptor */
#endif
	baud_rate	baud;
	parity		parity;
	stop_bits	stop_bits;
	word_length	word_length;
	int lerr;					/* Last error code */
	long pollcalib;				/* Number of loops for a seconds delay */
	int tc;						/* Current termination character (-1 if not set) */
	long timeout;				/* Current timout value in msec */
	int npaths;
	char **paths;				/* Paths if any */
	int debug;					/* Flag, nz to print serial coms info to stderr */

  /* Public: */

	/* Return a list of serial port paths */
	/* Return NULL on failure. */
	char ** (*get_paths)(			/* Return pointer to list of paths */
		struct _serio *p 			/* End is marked with NULL pointer, */
	);								/* Storage will be freed with serio object */

	/* Set the port number and characteristics */
	void (*set_port)(
		struct _serio *p, 
		char       *ppath,		/* Path to device if not NULL */
		int 		port,		/* Serial port number if ppath == NULL */
		baud_rate	baud,
		parity		parity,
		stop_bits	stop_bits,
		word_length	word_length);

	/* Write the characters in the buffer out */
	/* Data will be written up to the terminating nul. */
	/* Return relevant error status bits */
	int (*write)(
		struct _serio *p,
		char *buf,
		double tout);		/* Timeout in seconds */

	/* Read characters into the buffer */
	/* The returned data will be terminated by a nul. */
	int (*read)(
		struct _serio *p,
		char *buf,			/* Buffer to store characters read */
		int bsize,			/* Buffer size */
		char tc,			/* Terminating character */
		int ntc,			/* Number of terminating characters */
		double tout);		/* Timeout in seconds */

	/* write and read */
	int (*write_read)(
		struct _serio *p,
		char *wbuf,			/* Write puffer */
		char *rbuf,			/* Read buffer */
		int bsize,			/* Buffer size */
		char tc,			/* Terminating characer */
		int ntc,			/* Number of terminating characters */
		double tout);		/* Timeout in seconds */

	/* Destroy ourselves */
	void (*del)(struct _serio *p);

	}; typedef struct _serio serio;

/* Constructor */
extern serio *new_serio(void);

/* System dependent convenience functions */
int next_con_char(void); /* return the next character from the keyboard */
void empty_con_chars(void);

#define SERIO_H
#endif /* SERIO_H */
