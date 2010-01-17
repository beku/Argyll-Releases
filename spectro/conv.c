
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
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>

#ifdef NT
#include <conio.h>
#endif

#if defined(UNIX)
#include <sys/types.h>		/* Include sys/select.h ? */
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>
#include <errno.h>
#include <string.h>
#endif

#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "xspect.h"
#include "insttypes.h"
#include "icoms.h"
#include "conv.h"
#include "usbio.h"
/* select() defined, but not poll(), so emulate poll() */
#if defined(FD_CLR) && !defined(POLLIN)
#include "pollem.h"
#define poll_x pollem
#else
#include <sys/poll.h>	/* Else assume poll() is native */
#define poll_x poll
#endif

#ifdef __APPLE__
//#include <stdbool.h>
#include <sys/sysctl.h>
#include <sys/param.h>
#include <CoreFoundation/CoreFoundation.h>
#include <IOKit/IOKitLib.h>
#include <IOKit/serial/IOSerialKeys.h>
#include <IOKit/IOBSD.h>
#include <mach/mach_init.h>
#include <mach/task_policy.h>
#endif /* __APPLE__ */

#undef DEBUG

#ifdef DEBUG
# define errout stderr
# define DBG(xx)	fprintf(errout, xx )
# define DBGF(xx)	fprintf xx
#else
# define DBG(xx)
# define DBGF(xx)
#endif

#ifdef __BORLANDC__
#define _kbhit kbhit
#endif

/* ============================================================= */
#ifdef NT

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

#undef USE_BEGINTHREAD

/* Destroy the thread */
static void athread_del(
athread *p
) {
	DBG("athread_del called\n");

	if (p == NULL)
		return;

	if (p->th != NULL) {		/* Oops. this isn't good. */
		DBG("athread_del calling TerminateThread() because thread hasn't finished\n");
		TerminateThread(p->th, -1);		/* But it is worse to leave it hanging around */
		CloseHandle(p->th);
	}

	free(p);
}

/* _beginthread doesn't leak memory, but */
/* needs to be linked to a different library */
#ifdef USE_BEGINTHREAD
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
	CloseHandle(p->th);
	p->th = NULL;
#ifdef USE_BEGINTHREAD
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
#ifdef USE_BEGINTHREAD
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

#endif /* NT */
/* ============================================================= */

/* ============================================================= */
#if defined(UNIX)

/* wait for and return the next character from the keyboard */
int next_con_char(void) {
	struct pollfd pa[1];		/* Poll array to monitor stdin */
	struct termios origs, news;
	char rv = 0;

	/* Configure stdin to be ready with just one character */
	if (tcgetattr(STDIN_FILENO, &origs) < 0)
		error("tcgetattr failed with '%s' on stdin", strerror(errno));
	news = origs;
	news.c_lflag &= ~(ICANON | ECHO);
	news.c_cc[VTIME] = 0;
	news.c_cc[VMIN] = 1;
	if (tcsetattr(STDIN_FILENO,TCSANOW, &news) < 0)
		error("next_con_char: tcsetattr failed with '%s' on stdin", strerror(errno));

	/* Wait for stdin to have a character */
	pa[0].fd = STDIN_FILENO;
	pa[0].events = POLLIN | POLLPRI;
	pa[0].revents = 0;

	if (poll_x(pa, 1, -1) > 0
	 && (pa[0].revents == POLLIN
		 || pa[0].revents == POLLPRI)) {
		char tb[3];
		if (read(STDIN_FILENO, tb, 1) > 0) {	/* User hit a key */
			rv = tb[0] ;
		}
	} else {
		tcsetattr(STDIN_FILENO, TCSANOW, &origs);
		error("poll on stdin returned unexpected value 0x%x",pa[0].revents);
	}

	/* Restore stdin */
	if (tcsetattr(STDIN_FILENO, TCSANOW, &origs) < 0)
		error("tcsetattr failed with '%s' on stdin", strerror(errno));

	return rv;
}

/* If here is one, return the next character from the keyboard, else return 0 */
int poll_con_char(void) {
	struct pollfd pa[1];		/* Poll array to monitor stdin */
	struct termios origs, news;
	char rv = 0;

	/* Configure stdin to be ready with just one character */
	if (tcgetattr(STDIN_FILENO, &origs) < 0)
		error("tcgetattr failed with '%s' on stdin", strerror(errno));
	news = origs;
	news.c_lflag &= ~(ICANON | ECHO);
	news.c_cc[VTIME] = 0;
	news.c_cc[VMIN] = 1;
	if (tcsetattr(STDIN_FILENO,TCSANOW, &news) < 0)
		error("next_con_char: tcsetattr failed with '%s' on stdin", strerror(errno));

	/* Wait for stdin to have a character */
	pa[0].fd = STDIN_FILENO;
	pa[0].events = POLLIN | POLLPRI;
	pa[0].revents = 0;

	if (poll_x(pa, 1, 0) > 0
	 && (pa[0].revents == POLLIN
		 || pa[0].revents == POLLPRI)) {
		char tb[3];
		if (read(STDIN_FILENO, tb, 1) > 0) {	/* User hit a key */
			rv = tb[0] ;
		}
	}

	/* Restore stdin */
	if (tcsetattr(STDIN_FILENO, TCSANOW, &origs) < 0)
		error("tcsetattr failed with '%s' on stdin", strerror(errno));

	return rv;
}

/* Suck all characters from the keyboard */
void empty_con_chars(void) {
	tcflush(STDIN_FILENO, TCIFLUSH);
}

/* Sleep for the given number of msec */
void msec_sleep(unsigned int msec) {
#ifdef NEVER
	if (msec > 1000) {
		unsigned int secs;
		secs = msec / 1000;
		msec = msec % 1000;
		sleep(secs);
	}
	usleep(msec * 1000);
#else
	struct timespec ts;

	ts.tv_sec = msec / 1000;
	ts.tv_nsec = (msec % 1000) * 1000000;
	nanosleep(&ts, NULL);
#endif
}

/* Return the current time in msec. This is not related to any particular epoch */
unsigned int msec_time() {
	unsigned int rv;
	static struct timeval startup = { 0, 0 };
	struct timeval cv;

	gettimeofday(&cv, NULL);

	/* Set time to 0 on first invocation */
	if (startup.tv_sec == 0 && startup.tv_usec == 0)
		startup = cv;

	/* Subtract, taking care of carry */
	cv.tv_sec -= startup.tv_sec;
	if (startup.tv_usec > cv.tv_usec) {
		cv.tv_sec--;
		cv.tv_usec += 1000000;
	}
	cv.tv_usec -= startup.tv_usec;

	/* Convert usec to msec */
	rv = cv.tv_sec * 1000 + cv.tv_usec / 1000;

	return rv;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */

#ifdef NEVER    /* Not currently needed, or effective */

/* Set the current threads priority */
int set_interactive_priority() {
#ifdef __APPLE__
#ifdef NEVER
    int rv = 0;
    struct task_category_policy tcatpolicy;

    tcatpolicy.role = TASK_FOREGROUND_APPLICATION;

    if (task_policy_set(mach_task_self(),
        TASK_CATEGORY_POLICY, (thread_policy_t)&tcatpolicy,
        TASK_CATEGORY_POLICY_COUNT) != KERN_SUCCESS)
		rv = 1;
printf("~1 set to forground got %d\n",rv);
	return rv;
#else
    int rv = 0;
    struct thread_precedence_policy tppolicy;

    tppolicy.importance = 500;

    if (thread_policy_set(mach_thread_self(),
        THREAD_PRECEDENCE_POLICY, (thread_policy_t)&tppolicy,
        THREAD_PRECEDENCE_POLICY_COUNT) != KERN_SUCCESS)
		rv = 1;
printf("~1 set to important got %d\n",rv);
	return rv;
#endif /* NEVER */
#else /* !APPLE */
	int rv;
	struct sched_param param;
	param.sched_priority = 32;

	/* This doesn't work unless we're running as su :-( */
	rv = pthread_setschedparam(pthread_self(), SCHED_FIFO, &param);
//printf("Set got %d\n",rv);
	return rv;
#endif /* !APPLE */
}

int set_normal_priority() {
#ifdef __APPLE__
#ifdef NEVER
    int rv = 0;
    struct task_category_policy tcatpolicy;

    tcatpolicy.role = TASK_UNSPECIFIED;

    if (task_policy_set(mach_task_self(),
        TASK_CATEGORY_POLICY, (thread_policy_t)&tcatpolicy,
        TASK_CATEGORY_POLICY_COUNT) != KERN_SUCCESS)
		rev = 1;
printf("~1 set to normal got %d\n",rv);
#else
    int rv = 0;
    struct thread_precedence_policy tppolicy;

    tppolicy.importance = 1;

    if (thread_policy_set(mach_thread_self(),
        THREAD_STANDARD_POLICY, (thread_policy_t)&tppolicy,
        THREAD_STANDARD_POLICY_COUNT) != KERN_SUCCESS)
		rv = 1;
printf("~1 set to standard got %d\n",rv);
	return rv;
#endif /* NEVER */
#else /* !APPLE */
	struct sched_param param;
	param.sched_priority = 0;
	int rv;

	rv = pthread_setschedparam(pthread_self(), SCHED_OTHER, &param);
//printf("Reset got %d\n",rv);
	return rv;
#endif /* !APPLE */
}

#endif /* NEVER */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
static athread *beep_thread = NULL;
static int beep_delay;
static int beep_freq;
static int beep_msec;

/* Delayed beep handler */
static int delayed_beep(void *pp) {
	msec_sleep(beep_delay);
#ifdef __APPLE__
	SysBeep((beep_msec * 60)/1000);
#else
	fprintf(stdout, "\a"); fflush(stdout);
#endif 
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
#ifdef __APPLE__
		SysBeep((msec * 60)/1000);
#else
		/* Linux is pretty lame in this regard... */
		fprintf(stdout, "\a"); fflush(stdout);
#endif
	}
}


#ifdef NEVER
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* If we're UNIX and running on X11, we could do this */
/* sort of thing (from xset), IF we were linking with X11: */

static void
set_bell_vol(Display *dpy, int percent)
{
XKeyboardControl values;
XKeyboardState kbstate;
values.bell_percent = percent;
if (percent == DEFAULT_ON)
  values.bell_percent = SERVER_DEFAULT;
XChangeKeyboardControl(dpy, KBBellPercent, &values);
if (percent == DEFAULT_ON) {
  XGetKeyboardControl(dpy, &kbstate);
  if (!kbstate.bell_percent) {
    values.bell_percent = -percent;
    XChangeKeyboardControl(dpy, KBBellPercent, &values);
  }
}
return;
}

static void
set_bell_pitch(Display *dpy, int pitch)
{
XKeyboardControl values;
values.bell_pitch = pitch;
XChangeKeyboardControl(dpy, KBBellPitch, &values);
return;
}

static void
set_bell_dur(Display *dpy, int duration)
{
XKeyboardControl values;
values.bell_duration = duration;
XChangeKeyboardControl(dpy, KBBellDuration, &values);
return;
}

XBell(..);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#endif /* NEVER */

/* - - - - - - - - - - - - - - - - - - - - - - - - */

/* Destroy the thread */
static void athread_del(
athread *p
) {
	DBG("athread_del called\n");

	if (p == NULL)
		return;

	pthread_cancel(p->thid);
	pthread_join(p->thid, NULL);

	free(p);
}

static void *threadproc(
	void *param
) {
	athread *p = (athread *)param;

	p->result = p->function(p->context);
	return 0;
}
 

athread *new_athread(
	int (*function)(void *context),
	void *context
) {
	int rv;
	athread *p = NULL;

	DBG("new_athread called\n");
	if ((p = (athread *)calloc(sizeof(athread), 1)) == NULL)
		return NULL;

	p->function = function;
	p->context = context;
	p->del = athread_del;

	/* Create a thread */
	rv = pthread_create(&p->thid, NULL, threadproc, (void *)p);
	if (rv != 0) {
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
	unlink(fname);
}

/* - - - - - - - - - - - - - - - - - - - - - - - - */
#ifdef __APPLE__

/* Kill a particular named process. */
/* return < 0 if this fails. */
/* return 0 if there is no such process */
/* return 1 if a process was killed */
int kill_nprocess(char *pname) {
	struct kinfo_proc *procList = NULL;
	size_t procCount = 0;
	static int name[] = { CTL_KERN, KERN_PROC, KERN_PROC_ALL, 0 };
	size_t length;
	int rv = 0;
	int i;

	procList = NULL;
	for (;;) {
		int err;

		/* Establish the amount of memory needed */
		length = 0;
		if (sysctl(name, (sizeof(name) / sizeof(*name)) - 1,
                      NULL, &length,
                      NULL, 0) == -1) {
			DBGF((errout,"sysctl #1 failed with %d\n", errno));
			return -1;
		}
	
		/* Add some more entries in case the number of processors changed */
		length += 10 * sizeof(struct kinfo_proc);
		if ((procList = malloc(length)) == NULL) {
			DBGF((errout,"malloc failed for %d bytes\n", length));
			return -1;
		}

		/* Call again with memory */
		if ((err = sysctl( name, (sizeof(name) / sizeof(*name)) - 1,
                          procList, &length,
                          NULL, 0)) == -1) {
			DBGF((errout,"sysctl #1 failed with %d\n", errno));
			free(procList);
			return -1;
		}
		if (err == 0) {
			break;
		} else if (err == ENOMEM) {
			free(procList);
			procList = NULL;
		}
	}

	procCount = length / sizeof(struct kinfo_proc);

	/* Locate the ninjad */
	for (i = 0; i < procCount; i++) {
		if (strcmp(procList[i].kp_proc.p_comm,pname) == 0) {
			DBGF((errout,"killing process '%s' pid %d\n",pname,procList[i].kp_proc.p_pid));
			if (kill(procList[i].kp_proc.p_pid, SIGTERM) != 0) {
				fprintf(stderr, "kill failed with %d\n", errno);
				DBGF((errout,"kill process '%s' failed with %d\n",pname,errno));
				free(procList);
				return -1;
			}
			rv = 1;
		}
	}
	free(procList);
	return rv;
}

#endif /* __APPLE__ */

#endif /* defined(UNIX) */

/* ============================================================= */

