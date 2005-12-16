
#ifndef DISPWIN_H

/* 
 * Argyll Color Correction System
 * Display target patch window
 *
 * Author: Graeme W. Gill
 * Date:   4/10/96
 *
 * Copyright 1998 - 2004 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */


int do_plot(double *x, double *y1, double *y2, double *y3, int n);

#ifdef NT
#include <windows.h>
//#include <winuser.h>

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
#include <IOKit/Graphics/IOGraphicsLib.h>
#endif /* __APPLE__ */

#if defined(UNIX) && !defined(__APPLE__)
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/extensions/xf86vmode.h>
#include <X11/extensions/dpms.h>
#endif /* UNIX */


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


/* Dispwin object */
struct _dispwin {

/* private: */
	/* Plot instance information */
	int sx,sy;			/* Screen offset in pixels */
	int sw,sh;			/* Screen width and height in pixels*/

	int donat;			/* Do native output */
	ramdac *or;			/* Original ramdac contents */
	ramdac *r;			/* Ramdac in use for native mode */

#ifdef NT
	char *AppName;
	HWND hwnd;
	MSG msg;
	ATOM arv;
#endif /* NT */

#ifdef __APPLE__
	CGDirectDisplayID ddid;
	WindowRef mywindow;
	CGrafPtr port;
	CGContextRef mygc;
	double rgb[3];			/* Current color */
	int btf;				/* Flag, nz if window has been brought to the front once */
	int winclose;			/* Flag, set to nz if window was closed */
#endif /* __APPLE__ */

#if defined(UNIX) && !defined(__APPLE__)
	Display *mydisplay;
	int myscreen;
	Window mywindow;
	GC mygc;

    /* Screensaver state */
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
	int (*set_color)(struct _dispwin *p, double r, double g, double b);

	/* Destroy ourselves */
	void (*del)(struct _dispwin *p);

}; typedef struct _dispwin dispwin;

/* Create a display window, default white */
dispwin *new_dispwin(
	double width, double height,	/* Width and height in mm */
	double hoff, double voff,		/* Offset from c. in fraction of screen, range -1.0 .. 1.0 */
	int native,						/* NZ if ramdac should be bypassed */
	int override					/* NZ if override_redirect is to be used on X11 */
);

#define DISPWIN_H
#endif /* DISPWIN_H */

