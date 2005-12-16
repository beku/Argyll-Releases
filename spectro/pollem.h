#ifndef POLLEN_H

/* 
 * Argyll Color Correction System
 *
 * Unix serial I/O class poll() emulation.
 *
 * Author: Graeme W. Gill
 * Date:   12/9/2004
 *
 * Copyright 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* Fake up poll() support on systems that only support select() */

/* Fake poll array structure */
struct pollfd {
	int fd;			/* File descriptor */
	short events;	/* Requested events */
	short revents;	/* Returned events */
};

/* Fake Poll flag values supported */
#define POLLIN	0x01
#define	POLLPRI	0x02
#define	POLLOUT	0x04

int pollem(struct pollfd fds[], unsigned long nfds, int timeout);

#define POLLEN_H
#endif /* POLLEN_H */
