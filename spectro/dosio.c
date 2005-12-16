/* 
 * Argyll Color Correction System
 *
 * Simple PC direct polled serial I/O class
 *
 * Author: Graeme W. Gill
 * Date:   5/10/96
 *
 * Copyright 1996, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include <stdio.h>
#include <pc.h>
#include <time.h>
#include "serio.h"

#include <stdarg.h>
void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);

#undef DEBUG

/* return the next character from the keyboard */
int next_con_char(void) {
	return getch();
}

/* Suck all characters from the keyboard */
void empty_con_chars(void) {
	while (kbhit())
		getch();
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

/* Comm port base addresses */
static unsigned short comporta[4] = {
	0x3f8, 0x2f8, 0x3e8, 0x2e8 };

/* Register offsets */
#define RBR 0		/* Reveive buffer register */
#define THR 0		/* Transmit hold register */
#define IER 1		/* Interrupt enable register */
#define IIR 2		/* Interrupt ident. register */
#define LCR 3		/* Line control register */
#define MCR 4		/* Modem control register */
#define LSR 5		/* Line status register */
#define MSR 6		/* Modem status register */
#define SCR 7		/* Scratch register */
#define DLL 0		/* Divisor latch low */
#define DLH 1		/* Divisor latch high */

/* Set the port number and characteristics */
static void
serio_set_port(
struct _serio *p, 
int 		port,		/* com port, 1 - 4 */
baud_rate	baud,
parity		parity,
stop_bits	stop,
word_length	word)
{
	unsigned short div;		/* Baud rate divisor from 1.8432 Hz */
	unsigned char bt;

	if (port >= 1 && port <= 4) {
		p->port = port;
		p->porta = comporta[port-1];
	}
	if (baud != baud_nc)
		p->baud = baud;
	if (parity != parity_nc)
		p->parity = parity;
	if (stop != stop_nc)
		p->stop_bits = stop;
	if (word != length_nc)
		p->word_length = word;

	if (p->port <= 0)
		error("serio - set_port: illegal port!");

	switch (p->baud) {
		case baud_110:
			div = 1047;
			break;
		case baud_300:
			div = 384;
			break;
		case baud_600:
			div = 192;
			break;
		case baud_1200:
			div = 96;
			break;
		case baud_2400:
			div = 48;
			break;
		case baud_4800:
			div = 24;
			break;
		case baud_9600:
			div = 12;
			break;
		case baud_19200:
			div = 6;
			break;
		case baud_38400:
			div = 3;
			break;
		case baud_57600:
			div = 2;
			break;
		default:
			error("serio - set_port: illegal baud rate!");
			break;
	}

	/* figure the line control reg value */
	bt = 0x00;
	switch (p->parity) {
		case parity_nc:
			error("serio - set_port: illegal parity setting!");
			break;
		case parity_none:
			break;
		case parity_odd:
			bt |= 0x04;
			break;
		caseparity_even:
			bt |= 0x14;
			break;
	}

	switch (p->stop_bits) {
		case stop_nc:
			error("serio - set_port: illegal stop bits!");
			break;
		case stop_1:
			break;
		case stop_2:
			bt |= 0x02;
			break;
	}

	switch (p->word_length) {
		case length_nc:
			error("serio - set_port: illegal word length!");
		case length_5:
			break;
		case length_6:
			bt |= 0x01;
			break;
		case length_7:
			bt |= 0x02;
			break;
		case length_8:
			bt |= 0x03;
			break;
	}

	/* Set the chip up */
	outportb(p->porta + IER, 0x00);	/* Turn interrupts off */
	outportb(p->porta + LCR, 0x80);	/* Set DLAB bit */
	outportb(p->porta + DLL, div & 0xff);			/* Set low byte of divisor */
	outportb(p->porta + DLH, (div >> 8) & 0xff);	/* Set high byte of divisor */
	outportb(p->porta + LCR, bt);	/* Set line control and reset DLAB bit */
	outportb(p->porta + MCR, 0x00);	/* Reset all modem control signals and loopback */

	/* Figure the polling timeout callibration */
	if (p->pollcalib == 0) {
		int i,j,k, bt;
		long stime, ttime;
		stime = clock();
		for (j=k=i=0;i < 500000;i++) {
			bt = inportb(comporta[0] + LSR);
			if (bt & SIO_XHE)
				j++;
			else
				k++;
			if (kbhit())
				k++;
		}
		ttime = clock() - stime;
		p->pollcalib = (long)(((500000.0 * CLOCKS_PER_SEC)/ttime) + 0.5);
#ifdef DEBUG
		printf("serio: Polling rate constant = %d\n",p->pollcalib);
#endif
		if (p->pollcalib < 300)
			error("Polling rate from serial chip seems to low! (%d)",p->pollcalib);
	}
}

/* Write the characters in the buffer out */
/* Return relevant error status bits */
static int
serio_write(
struct _serio *p,
char *wbuf,
double tout)
{
	int rv;
	unsigned char bt;
	char c;
	long i, rwait = (int)(tout * p->pollcalib + 0.5);

	if (p->pollcalib == 0)
		error("serio_write: not initialised");
	p->lerr = 0;
	i = rwait;
	c = *wbuf++;
	while (i-- > 0 && c != '\000') {
		bt = inportb(p->porta + LSR);
		if (bt & SIO_XHE) {
			outportb(p->porta + THR, c); /* Send the character */
			c = *wbuf++;
			i = rwait;
		}
		if (kbhit() && getch() == 0x1b) {	/* User hit escape */
			p->lerr = SIO_USER;
#ifdef DEBUG
			printf("serio read user escape\n");
#endif
			break;
		}
	}
	if (i == 0) {
#ifdef DEBUG
		printf("serio write timed out\n");
#endif
		p->lerr = SIO_TO;
	}
	return p->lerr;
}

/* Read characters into the buffer */
static int
serio_read(struct _serio *p,
char *rbuf,			/* Buffer to store characters read */
int bsize,			/* Buffer size */
char tc,			/* Terminating characer */
int ntc,			/* Number of terminating characters */
double tout)		/* Time out in seconds */
{
	int rv = 0;
	unsigned char bt;
	int j = 0;
	char c;
	long i, rwait = (int)(tout * p->pollcalib + 0.5);

	if (p->pollcalib == 0)
		error("serio_read: not initialised");

	if (bsize < 2)
		error("serio_read given too small a buffer");

	p->lerr = 0;
	i = rwait;
	bt = 0;
	/* Poll till we get a char, or we time out */
	while (   i-- > 0
	       && !(bt & (SIO_BRK | SIO_FER | SIO_PER | SIO_OER))
	       && bsize > 1
	       && j < ntc) {
		bt = inportb(p->porta + LSR);
		if (bt & SIO_DRY) {
			*rbuf++ = c = inportb(p->porta + RBR);
	        if (c == tc)
				j++;
			bsize--;	
			i = rwait;
		}
		if (kbhit() && getch() == 0x1b) {	/* User hit escape */
			p->lerr = SIO_USER;
#ifdef DEBUG
			printf("serio read user escape\n");
#endif
			break;
		}
	}
	*rbuf = '\000';
	if (i <= 0) {			/* timed out */
#ifdef DEBUG
		printf("serio read timed out\n");
#endif
		p->lerr = SIO_TO;
	} else if ((rv = bt & (SIO_BRK | SIO_FER | SIO_PER | SIO_OER)) != 0) {
#ifdef DEBUG
		printf("serio got error 0x%x\n",rv);
#endif
		p->lerr = rv;
	}
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
	int rv = 0;
	unsigned char bt;
	char c;
	int j = 0;
	long i, rwait;

	if (bsize < 2)
		error("serio_read given too small a buffer");

	if (p->pollcalib == 0)
		error("serio_write_read: not initialised");

	p->lerr = 0;

	/* Flush any read characters */
	rwait = (int)(0.1 * p->pollcalib + 0.5);
	for (i=0;i < rwait;) {
		bt = inportb(p->porta + LSR);
		if (bt & SIO_DRY)
			bt = inportb(p->porta + RBR);
		else
			i++;
		if (kbhit() && getch() == 0x1b) {	/* User hit escape */
			p->lerr = SIO_USER;
#ifdef DEBUG
			printf("serio read user escape\n");
#endif
			break;
		}
	}

	if (p->lerr != 0)
		return p->lerr;

	/* Write the write data */
	rwait = (int)(tout * p->pollcalib + 0.5);
	i = rwait;
	c = *wbuf++;
	while (i-- > 0 && c != '\000' && j < ntc) {
		bt = inportb(p->porta + LSR);
		if ((bt & SIO_DRY) && bsize > 1) {
			char tt;
			if ((tt = inportb(p->porta + RBR)) == tc)
				j++;
			*rbuf++ = tt;
			bsize--;
		}
		if (bt & SIO_XHE) {
			outportb(p->porta + THR, c); /* Send the character */
			c = *wbuf++;
			i = rwait;
		}
		if (kbhit() && getch() == 0x1b) {	/* User hit escape */
			p->lerr = SIO_USER;
#ifdef DEBUG
			printf("serio read user escape\n");
#endif
			break;
		}
	}

	if (p->lerr != 0)
		return p->lerr;

	if (i <= 0) {
#ifdef DEBUG
		printf("serio write (read) timed out\n");
#endif
		p->lerr = SIO_TO;
		return p->lerr;
	}
	if (j >= ntc)	/* Got number of tc's */
		return p->lerr;

	/* Read response */
	i = rwait;
	bt = 0;
	/* Poll till we get a char, or we time out */
	while (   i-- > 0
	       && !(bt & (SIO_BRK | SIO_FER | SIO_PER | SIO_OER))
	       && bsize > 1
	       && j < ntc) {
		bt = inportb(p->porta + LSR);
		if (bt & SIO_DRY) {
			*rbuf++ = c = inportb(p->porta + RBR);
			if (c == tc)
				j++;
			bsize--;	
			i = rwait;
		}
		if (kbhit() && getch() == 0x1b) {	/* User hit escape */
			p->lerr = SIO_USER;
#ifdef DEBUG
			printf("serio read user escape\n");
#endif
			break;
		}
	}
	*rbuf = '\000';
	if (i <= 0) {			/* timed out */
#ifdef DEBUG
		printf("serio read (write) timed out\n");
#endif
		p->lerr = SIO_TO;
	} else if ((rv = bt & (SIO_BRK | SIO_FER | SIO_PER | SIO_OER)) != 0) {
#ifdef DEBUG
		printf("serio got error 0x%x\n",rv);
#endif
		p->lerr |= rv;
	}
	return p->lerr;
}

/* Destroy ourselves */
static void
serio_free(serio *p) {
	free (p);
}

/* Constructor */
extern serio *new_serio()
{
	serio *p;
	if ((p = (serio *)malloc(sizeof(serio))) == NULL)
		error("serio: malloc failed!");

	/* Init things to null values */
	p->pollcalib = 0;
	p->lerr = 0;
	p->port = -1;
	p->baud = baud_nc;
	p->parity = parity_nc;
	p->stop_bits = stop_nc;
	p->word_length = length_nc;
	
	p->set_port = serio_set_port;
	p->write = serio_write;
	p->read = serio_read;
	p->write_read = serio_write_read;
	p->free = serio_free;

	return p;
}

