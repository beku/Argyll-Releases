
/* 
 * Argyll Color Correction System
 *
 * Linux serial I/O class
 *
 * Author: Graeme W. Gill
 * Date:   18/11/2000
 *
 * Copyright 1997 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include <sys/types.h>		/* Include sys/select.h ? */
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>
#include <errno.h>
#include <string.h>
#include "serio.h"	/* Class we are implementing */

/* select() defined, but not poll(), so emulate poll() */
#if defined(FD_CLR) && !defined(POLLIN)
#include "pollem.h"
#define poll_x pollem
#else
#include <sys/poll.h>	/* Else assume poll() is native */
#define poll_x poll
#endif

#ifdef __APPLE__
#include <sys/param.h>
#include <CoreFoundation/CoreFoundation.h>
#include <IOKit/IOKitLib.h>
#include <IOKit/serial/IOSerialKeys.h>
#include <IOKit/IOBSD.h>
#endif /* __APPLE__ */

void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);

/* return the next character from the keyboard */
/* NOTE: should fix this so it doesn't echo !! */
int next_con_char(void) {
	struct pollfd pa[1];		/* Poll array to monitor stdin */
	struct termios origs, news;
	char buf[3];

	/* Configure stdin to be ready with just one character */
	if (tcgetattr(STDIN_FILENO, &origs) < 0)
		error("ycgetattr failed with '%s'", strerror(errno));
	news = origs;
	news.c_lflag &= (~ICANON);
	news.c_cc[VTIME] = 0;
	news.c_cc[VMIN] = 1;
	if (tcsetattr(STDIN_FILENO,TCSANOW, &news) < 0)
		error("tcsetattr failed with '%s'", strerror(errno));

	/* Wait for stdin to have a character */
	pa[0].fd = STDIN_FILENO;
	pa[0].events = POLLIN | POLLPRI;
	pa[0].revents = 0;

	if (poll_x(pa, 1, -1) > 0
	 && (pa[0].revents == POLLIN
		 || pa[0].revents == POLLPRI)) {
		char tb[3];
		if (read(STDIN_FILENO, tb, 1) > 0) {	/* User hit a key */
			/* Restore stdin */
			if (tcsetattr(STDIN_FILENO, TCSANOW, &origs) < 0)
				error("tcsetattr failed with '%s'", strerror(errno));
			return tb[0] ;
		}
	} else {
		error("poll on stdin returned unexpected value 0x%x",pa[0].revents);
	}
}

/* Suck all characters from the keyboard */
void empty_con_chars(void) {
	tcflush(STDIN_FILENO, TCIFLUSH);
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

/* Create and return a list of available serial ports for this system */
static char **
serio_get_paths(
struct _serio *p 
) {
	int i;

	/* Free any old list */
	if (p->paths != NULL) {
		for (i = 0; i < p->npaths; i++)
			free(p->paths[i]);
		free(p->paths);
		p->npaths = 0;
		p->paths = NULL;
	}

#ifdef __APPLE__
	/* Search the OSX registry for serial ports */
	{
	    kern_return_t kstat; 
	    mach_port_t mp;						/* Master IO port */
	    CFMutableDictionaryRef sdict;		/* Serial Port  dictionary */
		io_iterator_t mit;					/* Matching itterator */
		io_object_t ioob;					/* Serial object found */

		/* Open a master port */
		if ((kstat = IOMasterPort(MACH_PORT_NULL, &mp)) != KERN_SUCCESS) {
			warning("Opening IO master port failed with %d",kstat);
			return NULL;
		}

		/* Get dictionary of serial ports */
    	if ((sdict = IOServiceMatching(kIOSerialBSDServiceValue)) == NULL) {
        	warning("IOServiceMatching returned a NULL dictionary");
		}

		/* Set value to match to RS232 type serial */
        CFDictionarySetValue(sdict, CFSTR(kIOSerialBSDTypeKey), CFSTR(kIOSerialBSDRS232Type));

		/* Init itterator to find matching types */
		if ((kstat = IOServiceGetMatchingServices(mp, sdict, &mit)) != KERN_SUCCESS) 
        	error("IOServiceGetMatchingServices returned %d\n", kstat);

		/* Find all the matching serial ports */
		for (;;) {
			char pname[100];
			
	        CFTypeRef dfp;		/* Device file path */

		    if ((ioob = IOIteratorNext(mit)) == NULL)
				break;

		    /* Get the callout device's path (/dev/cu.xxxxx). */
			if ((dfp = IORegistryEntryCreateCFProperty(ioob, CFSTR(kIOCalloutDeviceKey),
			                                      kCFAllocatorDefault, 0)) == NULL)
				goto continue1;

			/* Convert from CF string to C string */
			if (!CFStringGetCString(dfp, pname, 100, kCFStringEncodingASCII))
				goto continue2;

			/* Ignore infra red port or Bluetooth */
			if (strstr(pname, "IrDA") != NULL
			 || strstr(pname, "Bluetooth") != NULL)
				goto continue2;

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
			if ((p->paths[p->npaths] = strdup(pname)) == NULL)
				error("serio: strdup failed!");
			p->npaths++;

		continue2:
            CFRelease(dfp);
		continue1:
		    IOObjectRelease(ioob);		/* Release found object */
		}
	    IOObjectRelease(mit);			/* Release the itterator */
	}
#else
	/* Other UNIX like systems */
	/* Many are crude and list every available device name, whether */
	/* it's usable or not. Do any UNIX systems have a mechanism for listing */
	/* serial ports ?? */
	/* For the moment, just create four ttyS* paths */

	p->npaths = 4;
	if ((p->paths = (char **)calloc(sizeof(char *), p->npaths + 1)) == NULL)
		error("serio: calloc failed!");

	for (i = 0; i < p->npaths; i++) {
		char pname[100];

		sprintf(pname,"/dev/ttyS%d",i);
		if ((p->paths[i] = strdup(pname)) == NULL)
			error("serio: strdup failed!");
	}
#endif

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
word_length	word
) {
	struct termios tio;
	speed_t speed;
	int i;
	long rwait = (int)(0.1 * 1000.0 + 0.5);		/* Timout every 100msec secs */

	if (p->debug) fprintf(stderr,"serio: About to set port characteristics\n");
	if (ppath != NULL || (port >= 1)) {
		if (p->fd != -1
		 && ((ppath != NULL && strcmp(ppath, p->ppath) != 0)
		 ||  (ppath == NULL && port != p->port))) {	/* If port number changes */
			close(p->fd);
			p->fd = -1;
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
	if (p->fd == -1) {
		if (p->ppath != NULL)
			free(p->ppath);

		if (p->ppath != NULL) {	/* Use given device path */
			p->ppath = strdup(ppath);
			p->port = -1;

		} else {		/* Synthesise com number */
			char pname[20];

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

		if ((p->fd = open(p->ppath, O_RDWR | O_NOCTTY )) < 0)
			error("Opening COM port '%s' failed with '%s'",p->ppath, strerror(errno));
		/* O_NONBLOCK O_SYNC */
		if (p->debug) fprintf(stderr,"serio: Opened port OK, fd = %d\n",p->fd);
	}

	if (tcgetattr(p->fd, &tio) < 0) {
		error("tcgetattr failed with '%s'", strerror(errno));
	}

	/* Set misc stuff */
	/* Input : */
	tio.c_iflag &= ~(
				  BRKINT		/* Don't interupt on Break */
				| IGNPAR		/* Don't ignore framing and parity errors */
				| PARMRK		/* Don't prefix framing or parity error characters with \377 and \0 */
				| ISTRIP		/* Don't strip off eighth bit */
				| INLCR			/* Don't translate NL to CR */
				| IGNCR			/* Don't ignore CR */
				| ICRNL			/* Don't translate CR to NL */
#ifndef __APPLE__
				| IUCLC			/* Don't map uppercase to lower case */
#endif
				| IXON			/* Don't enable XON/XOFF flow control on output */
				| IXANY			/* Don't enable any character to restart input */
				| IXOFF			/* Don't enable XON/XOFF flow control on input */
				| IMAXBEL		/* Don't ring the bell when the queue is full */
					);
	tio.c_iflag |= (
				  IGNBRK		/* Ignore Break */
				   );

	/* Output : */
	tio.c_oflag &= ~(
				  OPOST			/* Don't enable implementation defined output processing */
#ifndef __APPLE__
				| OLCUC			/* Don't map lower case to upper case */
				| OCRNL			/* Don't map CR to NL */
				| ONOCR			/* Don't not output CR at column 0 */
				| ONLRET		/* Don't not output CR */
				| OFILL			/* Don't send fill character for delay */
				| OFDEL			/* Don't use DEL as delay fill */
				| NLDLY			/* No NL delay */
				| CRDLY			/* No CR delay */
				| TABDLY		/* No Tab delay */
				| BSDLY			/* No BS delay */
				| VTDLY			/* No VT delay */
				| FFDLY			/* No FF delay */
#endif
				| ONLCR			/* Don't map NL to CR NL */
					);
	tio.c_oflag |= ( 0 );

	/* Control : zero everything out, then add it in again */
	tio.c_cflag &= ~(
				  CSIZE		/* Bits/char mask */
				| CSTOPB	/* Default 1 stop bits */
				| PARENB	/* No Parity checking */
				| PARODD	/* Default even parity */
				| HUPCL		/* No hang up */
#ifndef __APPLE__
				| CIBAUD	/* Input baud rate mask (not used) */
#endif
				| CRTSCTS	/* No Hardware flow control */
					);
	tio.c_cflag |= (
				  CREAD		/* Enable the receiver */
				| CLOCAL	/* Ignore modem control lines */
					);

	/* Local : zero everything out, then add it in again */
	tio.c_lflag &= ~(
				  ISIG		/* No signal on interrupt characters */
#ifndef __APPLE__
				| XCASE		/* No upper case mode */
#endif
				| ECHO		/* No echo of input characters */
				| ECHOE		/* No erase handling */
				| ECHOK		/* No kill handling */
				| ECHONL	/* No NL echoing */
				| ECHOCTL	/* No cntrl character echoing */
				| ECHOPRT	/* Don't print erased characters */
				| ECHOKE	/* No KILL erase handling */
				| FLUSHO	/* Output is not being flushed */
				| NOFLSH	/* Don't flush the input and output on signals */
				| TOSTOP	/* Don't sent SIGTTOU signal */
				| PENDIN	/* Don't repring input queue characters when next is read. */
				| IEXTEN	/* Disable implementation defined processing */
					);
	tio.c_lflag |= (
				  ICANON	/* Enable canonical, line buffered mode */
					);

	switch (p->parity) {
		case parity_nc:
			error("serio - set_port: illegal parity setting!");
			break;
		case parity_none:
			tio.c_iflag &= ~INPCK;		/* Disable input parity checking */
			break;
		case parity_odd:
			tio.c_iflag |= INPCK;		/* Enable input parity checking */
			tio.c_cflag |= PARENB;		/* Enable input and output parity checking */
			tio.c_cflag |= PARODD;		/* Input and output parity is odd */
			break;
		caseparity_even:
			tio.c_iflag |= INPCK;		/* Enable input parity checking */
			tio.c_cflag |= PARENB;		/* Enable input and output parity checking */
			break;
	}

	switch (p->stop_bits) {
		case stop_nc:
			error("serio - set_port: illegal stop bits!");
			break;
		case stop_1:
			break;		/* defaults to 1 */
		case stop_2:
			tio.c_cflag |= CSTOPB;
			break;
	}

	switch (p->word_length) {
		case length_nc:
			error("serio - set_port: illegal word length!");
		case length_5:
			tio.c_cflag |= CS5;
			break;
		case length_6:
			tio.c_cflag |= CS6;
			break;
		case length_7:
			tio.c_cflag |= CS7;
			break;
		case length_8:
			tio.c_cflag |= CS8;
			break;
	}

	/* Set the default canonical control codes */
	for (i = 0; i < NCCS; i++)
		tio.c_cc[i] = '\0';

	/* Set the baud rate */
	switch (p->baud) {
		case baud_110:
			speed = B110;
			break;
		case baud_300:
			speed = B300;
			break;
		case baud_600:
			speed = B600;
			break;
		case baud_1200:
			speed = B1200;
			break;
		case baud_2400:
			speed = B2400;
			break;
		case baud_4800:
			speed = B4800;
			break;
		case baud_9600:
			speed = B9600;
			break;
		case baud_19200:
			speed = B19200;
			break;
		case baud_38400:
			speed = B38400;
			break;
		case baud_57600:
			speed = B57600;
			break;
		case baud_115200:
			speed = B115200;
			break;
		default:
			error("serio - set_port: illegal baud rate!");
			break;
	}

	if (cfsetispeed(&tio,  speed) < 0)
		error("cfsetispeed failed with '%s'", strerror(errno));
	if (cfsetospeed(&tio,  speed) < 0)
		error("cfsetospeed failed with '%s'", strerror(errno));

	/* Make change immediately */
	tcflush(p->fd, TCIFLUSH);
	if (tcsetattr(p->fd, TCSANOW, &tio) < 0)
		error("tcsetattr failed with '%s'", strerror(errno));

	/* Set the timout value */
	p->timeout = rwait;
	if (p->debug) fprintf(stderr,"serio: port characteristics set ok\n");
}

/* Write the characters in the buffer out */
/* Data will be written up to the terminating nul */
/* Return relevant error status bits */
static int
serio_write(
struct _serio *p,
char *wbuf,
double tout
) {
	int wbytes;
	int len;
	long i, rwait = (int)(tout / 0.1 + 0.5);	/* Number of 100msec intervals to loop */
	struct pollfd pa[2];		/* Poll array to monitor serial write and stdin */
	struct termios origs, news;

	if (p->debug) fprintf(stderr,"About to write '%s'\n",fix(wbuf));
	if (p->fd == -1)
		error("serio_write: not initialised");

	/* Configure stdin to be ready with just one character */
	if (tcgetattr(STDIN_FILENO, &origs) < 0)
		error("ycgetattr failed with '%s'", strerror(errno));
	news = origs;
	news.c_lflag &= (~ICANON);
	news.c_cc[VTIME] = 0;
	news.c_cc[VMIN] = 1;
	if (tcsetattr(STDIN_FILENO,TCSANOW, &news) < 0)
		error("tcsetattr failed with '%s'", strerror(errno));

	/* Wait for serial output not block */
	pa[0].fd = p->fd;
	pa[0].events = POLLOUT;
	pa[0].revents = 0;

	/* Wait for stdin to have a character */
	pa[1].fd = STDIN_FILENO;
	pa[1].events = POLLIN | POLLPRI;
	pa[1].revents = 0;

	/* Until timed out, aborted, or transmitted */
	len = strlen(wbuf);
	p->lerr = 0;
	i = rwait;

	while (i-- > 0 && len > 0) {
		if (poll_x(pa, 2, p->timeout) > 0) {
			if (pa[0].revents != 0) {
				if (pa[0].revents == POLLOUT) {
					/* We can write it without blocking */
					if (write(p->fd, wbuf, len) != len) {
						error("write to serial failed with '%s'", strerror(errno));
					}
					break;
				} else {
					error("poll on serial out returned unexpected value 0x%x",pa[0].revents);
				}
			} else if (pa[1].revents != 0) {
				if (pa[1].revents == POLLIN
				 || pa[1].revents == POLLPRI) {
					char tb[10];
					/* Check for user abort (This should be a callback!) */
					if (read(STDIN_FILENO, tb, 10) > 0) {	/* User hit escape */
						if (tb[0] == 0x1b) {	/* User hit escape */
							p->lerr |= SIO_USER;
							break;
						}
						/* Else ignore it */
					}
				} else {
					error("poll on stdin  returned unexpected value 0x%x",pa[1].revents);
				}
			}
		}
		/* We timed out on both - loop */
	}
	if (i <= 0) {		/* Timed out */
		p->lerr |= SIO_TO;
	}

	/* Restore stdin */
	if (tcsetattr(STDIN_FILENO, TCSANOW, &origs) < 0)
		error("tcsetattr failed with '%s'", strerror(errno));

	if (p->debug) fprintf(stderr,"serio: Write returning with 0x%x\n",p->lerr);
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
double tout			/* Time out in seconds */
) {
	int rbytes;
	int j = 0;
	long i, rwait = (int)(tout * 10.0 + 0.5);	/* number of 100msec intervals to loop */
	struct pollfd pa[2];		/* Poll array to monitor serial read and stdin */
	struct termios origs, news;
	char *rrbuf = rbuf;		/* Start of return buffer */

	if (p->debug) fprintf(stderr,"serio: Read called\n");
	if (p->fd == -1)
		error("serio_read: not initialised");

	if (bsize < 3)
		error("serio_read given too small a buffer");

	if (tc != p->tc) {	/* Set the termination char */
		struct termios tio;

		if (tcgetattr(p->fd, &tio) < 0)
			error("ycgetattr failed with '%s'", strerror(errno));

		tio.c_cc[VEOL] = tc;

		/* Make change immediately */
		tcflush(p->fd, TCIFLUSH);
		if (tcsetattr(p->fd, TCSANOW, &tio) < 0)
			error("tcsetattr failed with '%s'", strerror(errno));

		p->tc = tc;
	}

	/* Configure stdin to be ready with just one character */
	if (tcgetattr(STDIN_FILENO, &origs) < 0)
		error("ycgetattr failed with '%s'", strerror(errno));
	news = origs;
	news.c_lflag &= (~ICANON);
	news.c_cc[VTIME] = 0;
	news.c_cc[VMIN] = 1;
	if (tcsetattr(STDIN_FILENO,TCSANOW, &news) < 0)
		error("tcsetattr failed with '%s'", strerror(errno));

	/* Wait for serial input to have data */
	pa[0].fd = p->fd;
	pa[0].events = POLLIN | POLLPRI;
	pa[0].revents = 0;

	/* Wait for stdin to have a character */
	pa[1].fd = STDIN_FILENO;
	pa[1].events = POLLIN | POLLPRI;
	pa[1].revents = 0;

	p->lerr = 0;
	i = rwait;
	bsize--;	/* Allow space for null */
	/* Until timed out, aborted, or we get what we want */
	while (i-- > 0 && bsize > 1 && j < ntc) {

		rbytes = 0;
		if (poll_x(pa, 2, p->timeout) > 0) {
			if (pa[0].revents != 0) {
				if (pa[0].revents == POLLIN
				 || pa[0].revents == POLLPRI) {
					/* We have data to read from input */
					if ((rbytes = read(p->fd, rbuf, bsize)) < 0) {
						error("read from serial failed with '%s'", strerror(errno));
					}
				} else {
					error("poll on serial in returned unexpected value 0x%x",pa[0].revents);
				}
			} else if (pa[1].revents != 0) {
				if (pa[1].revents == POLLIN
				 || pa[1].revents == POLLPRI) {
					char tb[10];
					/* Check for user abort (This should be a callback!) */
					if (read(STDIN_FILENO, tb, 10) > 0) {	/* User hit escape */
						if (tb[0] == 0x1b) {	/* User hit escape */
							p->lerr |= SIO_USER;
							break;
						}
						/* Else ignore it */
					}
				} else {
					error("poll on stdin  returned unexpected value 0x%x",pa[1].revents);
				}
			}
		}

		/* Account for bytes done (if any) */
		if (rbytes > 0) { /* Account for bytes done */
			bsize -= rbytes;
			i = rwait;		/* Reset time */

			while(rbytes--) {	/* Count termination characters */
				if (*rbuf++ == tc)
					j++;
			}
		}
	}

	*rbuf = '\000';
	if (i <= 0) {			/* timed out */
		p->lerr |= SIO_TO;
	}
	if (p->debug) fprintf(stderr,"serio: About to return read '%s' lerr 0x%x\n",fix(rrbuf),p->lerr);

	/* Restore stdin */
	if (tcsetattr(STDIN_FILENO, TCSANOW, &origs) < 0)
		error("tcsetattr failed with '%s'", strerror(errno));

	if (p->debug) fprintf(stderr,"serio: Read returning with 0x%x\n",p->lerr);

	return p->lerr;
}

#ifdef NEVER
/* ~~ haven't added serial error handling ~~~~~ !!!  */
		if (!ReadFile(p->fd, rbuf, bsize, &rbytes, NULL)) {
			DWORD errs;
			if (!ClearCommError(p->fd,&errs,NULL))
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
#endif

/* write and read */
static int
serio_write_read(
struct _serio *p,
char *wbuf,			/* Write puffer */
char *rbuf,			/* Read buffer */
int bsize,			/* Buffer size */
char tc,			/* Terminating characer */
int ntc,			/* Number of terminating characters */
double tout
) {
	if (p->debug) fprintf(stderr,"\nserio: Write_Read called with '%s'\n",fix(wbuf));

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
	if (p->lerr != 0) {
		if (p->debug) fprintf(stderr,"serio: Write_Read Write failed - returning 0x%x\n",p->lerr);
		return p->lerr;
	}

	/* Read response */
	serio_read(p, rbuf, bsize, tc, ntc, tout);
	if (p->debug) {
		if (p->lerr != 0)
			fprintf(stderr,"serio: Write_Read Write failed - returning 0x%x\n",p->lerr);
		else
			fprintf(stderr,"serio: Write_Read Write_Read success, returning '%s'\n",fix(rbuf));
	}
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
	if (p->fd != -1)
		close(p->fd);
	if (p->ppath != NULL)
		free(p->ppath);
	free (p);
}

/* Constructor */
extern serio *new_serio() {
	serio *p;
	if ((p = (serio *)calloc(sizeof(serio), 1)) == NULL)
		error("serio: malloc failed!");

	/* Init things to null values */
	p->fd = -1;
	p->pollcalib = 0;
	p->lerr = 0;
	p->ppath = NULL;
	p->port = -1;
	p->baud = baud_nc;
	p->parity = parity_nc;
	p->stop_bits = stop_nc;
	p->word_length = length_nc;
	p->debug = 0;
	
	p->get_paths = serio_get_paths;
	p->set_port = serio_set_port;
	p->write = serio_write;
	p->read = serio_read;
	p->write_read = serio_write_read;
	p->del = serio_del;

	return p;
}

