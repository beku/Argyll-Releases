
#ifndef CONV_H

/* 
 * Argyll Color Correction System
 *
 * Some system dependent comvenience functions.
 * Implemented in unixio.c and ntio.c
 *
 * Author: Graeme W. Gill
 * Date:   2008/2/9
 *
 * Copyright 1996 - 2008 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 * 
 * Derived from icoms.h
 */

#if defined (NT)
#if !defined(_WIN32_WINNT) || _WIN32_WINNT < 0x0501
# define _WIN32_WINNT 0x0501
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#if defined (UNIX) || defined(__APPLE__)
#include <pthread.h>
#endif

/* - - - - - - - - - - - - - - - - - - -- */
/* System dependent convenience functions */

/* wait for and then return the next character from the keyboard */
int next_con_char(void);

/* If there is one, return the next character from the keyboard, else return 0 */
int poll_con_char(void);

/* Empty the console of any pending characters */
void empty_con_chars(void);

/* Sleep for the given number of msec */
void msec_sleep(unsigned int msec);

/* Return the current time in msec since */
/* the process started. */
unsigned int msec_time();

/* Activate the system beeper after a delay */
/* (Note frequancy and duration may not be honoured on all systems) */
void msec_beep(int delay, int freq, int msec);

void normal_beep(); /* Emit a "normal" beep */
void good_beep(); /* Emit a "good" beep */
void bad_beep(); /* Emit a "bad" double beep */

/* - - - - - - - - - - - - - - - - - - -- */

#ifdef NEVER	/* Not currently needed, or effective */

/* Set the current threads priority */
/* return nz if this fails */
int set_interactive_priority();

int set_normal_priority();

#endif /* NEVER */

/* - - - - - - - - - - - - - - - - - - -- */

/* An Argyll thread. */
struct _athread {
#if defined (NT)
	HANDLE th;				/* Thread */
	DWORD thid;				/* Thread ID */
#endif
#if defined (UNIX) || defined(__APPLE__)
	pthread_t thid;			/* Thread ID */
#endif
	int result;				/* Return code from thread function */

	/* Thread function to call */
	int (*function)(void *context);

	/* And the context to call it with */
	void *context;

    /* Kill the thread and delete the object */
	/* (Killing it may have side effects, so this is a last */
	/*  resort if the thread hasn't exited) */
    void (*del)(struct _athread *p);

}; typedef struct _athread athread;

/* Create and start a thread */
/* Thread function should only return on completion or error. */
/* It should return 0 on completion or exit, nz on error. */
athread *new_athread(int (*function)(void *context), void *context);

#ifdef NEVER

/* Ideas for worker variant on thread: */

	/* Create a new worker thread, and put it to sleep */
	athread *new_aworker();
	/* Give the worker a job to do */
   	         ->start_work(int (*function)(void *context), void *context);
	/* See if the worker is finished its job */
	result = ->poll_work();
	/* Wait until the worker has finished its job */
	result = ->wait_work();

#endif

/* - - - - - - - - - - - - - - - - - - -- */

/* Delete a file */
void delete_file(char *fname);

/* - - - - - - - - - - - - - - - - - - -- */

#ifdef __APPLE__

/* Kill a particular named process. */
/* return < 0 if this fails. */
/* return 0 if there is no such process */
/* return 1 if a process was killed */
int kill_nprocess(char *pname);

#endif /* __APPLE__ */

#define CONV_H
#endif /* CONV_H */
