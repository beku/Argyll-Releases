
#ifndef DISPWIN_H

/* 
 * Argyll Color Correction System
 * Display target patch window
 *
 * Author: Graeme W. Gill
 * Date:   4/10/96
 *
 * Copyright 1998 - 2007 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */


int do_plot(double *x, double *y1, double *y2, double *y3, int n);

#ifdef NT
#define OEMRESOURCE
#include <windows.h>

#if(WINVER < 0x0500)
#error Require WINVER >= 0x500 to compile (multi-monitor API needed)
#endif

#ifndef COLORMGMTCAPS	/* In case SDK is out of date */

#define COLORMGMTCAPS   121

#define CM_NONE             0x00000000
#define CM_DEVICE_ICM       0x00000001
#define CM_GAMMA_RAMP       0x00000002
#define CM_CMYK_COLOR       0x00000004

#endif	/* !COLORMGMTCAPS */
#endif /* NT */

#ifdef __APPLE__	/* Assume OSX Carbon */
#include <Carbon/Carbon.h>
#include <CoreServices/CoreServices.h>
#include <IOKit/Graphics/IOGraphicsLib.h>
#endif /* __APPLE__ */

#if defined(UNIX) && !defined(__APPLE__)
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/extensions/xf86vmode.h>
#include <X11/extensions/dpms.h>
#include <X11/extensions/Xinerama.h>
#include <X11/extensions/scrnsaver.h>
#endif /* UNIX */

/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Enumerate and list all the available displays */

/* Structure to store infomation about possible displays */
typedef struct {
	char *description;	/* Description of display */
	int sx,sy;			/* Displays offset in pixels */
	int sw,sh;			/* Displays width and height in pixels*/
#ifdef NT
	char name[CCHDEVICENAME];	/* Display path */
#endif /* NT */
#ifdef __APPLE__
	CGDirectDisplayID ddid;
#endif /* __APPLE__ */
#if defined(UNIX) && !defined(__APPLE__)
	char *name;				/* Display name */
	int screen;				/* Screen to select */
	int uscreen;			/* Underlying screen */
	int rscreen;			/* Underlying RAMDAC screen */
	char *icc_atom_name;	/* Name of ICC profile atom for this display */
#endif /* UNIX */
} disppath;

/* Return pointer to list of disppath. Last will be NULL. */
/* Return NULL on failure. */
/* Call free_disppaths() to free up allocation */
disppath **get_displays();

void free_disppaths(disppath **paths);

/* Return the given display given its index 0..n-1 */
disppath *get_a_display(int ix);

void free_a_disppath(disppath *path);


/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Structure to handle RAMDAC values */
struct _ramdac {
	int pdepth;		/* Plane depth, usually 8 */
	int nent;		/* Number of entries, =  2^pdepth */
	double *v[3];	/* 2^pdepth entries for RGB, values 0.0 - 1.0 */

	/* Clone ourselves */
	struct _ramdac *(*clone)(struct _ramdac *p);

	/* Set the curves to linear */
	void (*setlin)(struct _ramdac *p);

	/* Destroy ourselves */
	void (*del)(struct _ramdac *p);
}; typedef struct _ramdac ramdac;


/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Dispwin object */
struct _dispwin {

/* private: */
	/* Plot instance information */
	int sx,sy;			/* Screen offset in pixels */
	int sw,sh;			/* Screen width and height in pixels*/
	int ww,wh;			/* Window width and height */
	int tx,ty;			/* Test area within window offset in pixels */
	int tw,th;			/* Test area width and height in pixels */

	double rgb[3];		/* Current color (full resolution) */
	double r_rgb[3];	/* Current color (raster value) */
	int nowin;			/* Don't create a test window */
	int donat;			/* Do native output */
	ramdac *or;			/* Original ramdac contents */
	ramdac *r;			/* Ramdac in use for native mode */
	int blackbg;		/* NZ if black full screen background */

	char *callout;		/* if not NULL - set color Shell callout routine */

#ifdef NT
	HDC hdc;			/* Handle to display */
	char *AppName;
	HWND hwnd;			/* Window handle */
	HCURSOR curs;		/* Invisible cursor */
	
	MSG msg;
	ATOM arv;
#endif /* NT */

#ifdef __APPLE__
	CGDirectDisplayID ddid;
	WindowRef mywindow;
	CGrafPtr port;
	CGContextRef mygc;
	int firstdraw;			/* Flag, nz if black background has been drawn */
	int btf;				/* Flag, nz if window has been brought to the front once */
	int winclose;			/* Flag, set to nz if window was closed */
#endif /* __APPLE__ */

#if defined(UNIX) && !defined(__APPLE__)
	Display *mydisplay;
	int myscreen;			/* Usual or virtual screen with Xinerama */
	int myuscreen;			/* Underlying screen */
	int myrscreen;			/* Underlying RAMDAC screen */
	Window mywindow;
	GC mygc;

    /* Screensaver state */
	int sssuspend;				/* Was able to suspend the screensaver */

	int ssvalid;				/* Was able to save & disable screensaver */
	int timeout, interval;
	int prefer_blanking;
	int allow_exposures;

	int dpmsdisabled;			/* monitor disabled by DPM */
	
#endif /* UNIX */

/* public: */
	int pdepth;				/* Plane depth of display */

	/* Get RAMDAC values. ->del() when finished. */
	/* Return NULL if not possible */
	ramdac *(*get_ramdac)(struct _dispwin *p);

	/* Set the RAMDAC values. */
	/* Return nz if not possible */
	int (*set_ramdac)(struct _dispwin *p, ramdac *r);

	/* Set a color (values 0.0 - 1.0) */
	/* Return nz on error */
	int (*set_color)(struct _dispwin *p, double r, double g, double b);

	/* Set a shell set color callout command line */
	void (*set_callout)(struct _dispwin *p, char *callout);

	/* Destroy ourselves */
	void (*del)(struct _dispwin *p);

}; typedef struct _dispwin dispwin;

/* Create a RAMDAC access and display test window, default white */
dispwin *new_dispwin(
	disppath *screen,				/* Screen to calibrate. */
	double width, double height,	/* Width and height in mm */
	double hoff, double voff,		/* Offset from c. in fraction of screen, range -1.0 .. 1.0 */
	int nowin,						/* NZ if no window should be created - RAMDAC access only */
	int native,						/* NZ if ramdac should be bypassed instead of being used */
	int blackbg,					/* NZ if whole screen should be filled with black */
	int override					/* NZ if override_redirect is to be used on X11 */
);


#define DISPWIN_H
#endif /* DISPWIN_H */

