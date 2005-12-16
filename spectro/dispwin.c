
/* 
 * Argyll Color Correction System
 * Display target patch window
 *
 * Author: Graeme W. Gill
 * Date:   4/10/96
 *
 * Copyright 1998 - 2004, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* This program displays test patches on a WinNT, MAC OSX or X11 windowing system. */

/*
   We are rather rough with how we handle window messages. We should
   really start another thread/process to handle the messages, rather
   that only servicing messages when we feel like it.
*/

/* TTBD
 * Add interface for listing and selecting the display to profile
 * (similar to selecting a serial port).
 *
 * Should probably check the display attributes (like visual depth)
 * and complain if we aren't using 24 bit color or better. 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include "copyright.h"
#include "config.h"
#include "cgats.h"
#include "dispwin.h"

#undef DEBUG
//#define STANDALONE_TEST	//~~~1
#define DONT_MESS_WITH_SCREENSAVER

#ifdef DEBUG
# define debug(xx)	fprintf(stderr, xx )
#else
# define debug(xx)
#endif

/* ----------------------------------------------- */

static ramdac *dispwin_clone_ramdac(ramdac *r);
static void dispwin_setlin_ramdac(ramdac *r);
static void dispwin_del_ramdac(ramdac *r);

/* For RAMDAC use, we assume that the number of entries in the RAMDAC */
/* meshes perfectly with the display raster depth, so that we can */
/* figure out how to apportion device values. We fail if they don't */
/* seem to mesh. */

/* Get RAMDAC values. ->del() when finished. */
/* Return NULL if not possible */
static ramdac *dispwin_get_ramdac(dispwin *p) {
	ramdac *r = NULL;
	int i, j;

#ifdef NT
	HDC hdc;
	WORD vals[3][256];		/* 16 bit elements */

	debug("dispwin_get_ramdac called\n");

	if ((hdc = GetDC(p->hwnd)) == NULL)
		return NULL;
	
	if ((GetDeviceCaps(hdc, COLORMGMTCAPS) & CM_GAMMA_RAMP) == 0) {
		DeleteDC(hdc);
		return NULL;
	}

	/* Allocate a ramdac */
	if ((r = (ramdac *)calloc(sizeof(ramdac), 1)) == NULL) {
		DeleteDC(hdc);
		return NULL;
	}
	r->pdepth = p->pdepth;
	r->nent = (1 << p->pdepth);
	r->clone = dispwin_clone_ramdac;
	r->setlin = dispwin_setlin_ramdac;
	r->del = dispwin_del_ramdac;
	for (j = 0; j < 3; j++) {

		if ((r->v[j] = (double *)calloc(sizeof(double), r->nent)) == NULL) {
			for (j--; j >= 0; j--)
				free(r->v[j]);
			free(r);
			DeleteDC(hdc);
			return NULL;
		}
	}

	if (GetDeviceGammaRamp(hdc, vals) == 0) {
		DeleteDC(hdc);
		return NULL;
	}
	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			r->v[j][i] = vals[j][i]/65535.0;
		}
	}

	DeleteDC(hdc);
#endif	/* NT */

#ifdef __APPLE__
	int nent;
	CGGammaValue vals[3][16385];

	debug("dispwin_get_ramdac called\n");

	if (CGGetDisplayTransferByTable(p->ddid, 163845, vals[0], vals[1], vals[2], &nent) != 0) {
		debug("CGGetDisplayTransferByTable failed\n");
		return NULL;
	}

	if (nent == 16385) {	/* oops - we didn't provide enought space! */
		debug("CGGetDisplayTransferByTable has more entries than we can handle\n");
		return NULL;
	}

	if (nent != (1 << p->pdepth)) {
		debug("CGGetDisplayTransferByTable number of entries mismatches screen depth\n");
		return NULL;
	}

	/* Allocate a ramdac */
	if ((r = (ramdac *)calloc(sizeof(ramdac), 1)) == NULL)
		return NULL;

	r->pdepth = p->pdepth;
	r->nent = (1 << p->pdepth);
	r->clone = dispwin_clone_ramdac;
	r->setlin = dispwin_setlin_ramdac;
	r->del = dispwin_del_ramdac;
	for (j = 0; j < 3; j++) {

		if ((r->v[j] = (double *)calloc(sizeof(double), r->nent)) == NULL) {
			for (j--; j >= 0; j--)
				free(r->v[j]);
			free(r);
			return NULL;
		}
	}

	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			r->v[j][i] = vals[j][i];
		}
	}
#endif /* __APPLE__ */

#if defined(UNIX) && !defined(__APPLE__)
	unsigned short vals[3][16384];
	int nent = 0;
	int evb = 0, erb = 0;

	debug("dispwin_get_ramdac called\n");

	if (XF86VidModeQueryExtension(p->mydisplay, &evb, &erb) == 0) {
		debug("XF86VidModeQueryExtension failed\n");
		return NULL;
	}

	if (XF86VidModeGetGammaRampSize(p->mydisplay, p->myscreen, &nent) == 0) {
		debug("XF86VidModeGetGammaRampSize failed\n");
		return NULL;
	}
	if (nent == 0) {
		debug("XF86VidModeGetGammaRampSize returned 0 size\n");
		return NULL;
	}

	if (nent > 16384) {
		debug("XF86VidModeGetGammaRampSize has more entries than we can handle\n");
		return NULL;
	}

	if (XF86VidModeGetGammaRamp(p->mydisplay, p->myscreen, nent,  vals[0], vals[1], vals[2]) == 0) {
		debug("XF86VidModeGetGammaRamp failed\n");
		return NULL;
	}

	if (nent != (1 << p->pdepth)) {
		debug("CGGetDisplayTransferByTable number of entries mismatches screen depth\n");
		return NULL;
	}

//for (i = 0; i < nent; i++) printf("red %d = %d\n",i,vals[0][i]);

	/* Allocate a ramdac */
	if ((r = (ramdac *)calloc(sizeof(ramdac), 1)) == NULL)
		return NULL;

	r->pdepth = p->pdepth;
	r->nent = (1 << p->pdepth);
	r->clone = dispwin_clone_ramdac;
	r->setlin = dispwin_setlin_ramdac;
	r->del = dispwin_del_ramdac;
	for (j = 0; j < 3; j++) {

		if ((r->v[j] = (double *)calloc(sizeof(double), r->nent)) == NULL) {
			for (j--; j >= 0; j--)
				free(r->v[j]);
			free(r);
			return NULL;
		}
	}

	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			r->v[j][i] = vals[j][i]/65535.0;
		}
	}
#endif /* UNIX && !__APPLE__ */
	return r;
}

/* Set the RAMDAC values. */
/* Return nz if not possible */
static int dispwin_set_ramdac(dispwin *p, ramdac *r) {
	int i, j;

#ifdef NT
	HDC hdc;
	WORD vals[3][256];		/* 16 bit elements */

	debug("dispwin_set_ramdac called\n");

	if ((hdc = GetDC(p->hwnd)) == NULL)
		return 1;
	
	if ((GetDeviceCaps(hdc, COLORMGMTCAPS) & CM_GAMMA_RAMP) == 0) {
		DeleteDC(hdc);
		return 1;
	}

	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			vals[j][i] = (int)(65535.0 * r->v[j][i] + 0.5);
		}
	}

	if (SetDeviceGammaRamp(hdc, vals) == 0) {
		DeleteDC(hdc);
		return 1;
	}

	DeleteDC(hdc);
#endif	/* NT */

#ifdef __APPLE__
	CGGammaValue vals[3][16384];

	debug("dispwin_set_ramdac called\n");

	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			vals[j][i] = r->v[j][i];
		}
	}

	if (CGSetDisplayTransferByTable(p->ddid, r->nent, vals[0], vals[1], vals[2]) != 0) {
		debug("CGSetDisplayTransferByTable failed\n");
		return 1;
	}
#endif /* __APPLE__ */

#if defined(UNIX) && !defined(__APPLE__)
	unsigned short vals[3][16384];

	debug("dispwin_set_ramdac called\n");

	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			vals[j][i] = (int)(r->v[j][i] * 65535.0 + 0.5);
		}
	}

	if (XF86VidModeSetGammaRamp(p->mydisplay, p->myscreen, r->nent, vals[0], vals[1], vals[2]) == 0) {
		debug("XF86VidModeSetGammaRamp failed\n");
		return 1;
	}
#endif /* UNIX && !__APPLE__ */

	return 0;
}


/* Clone ourselves */
static ramdac *dispwin_clone_ramdac(ramdac *r) {
	ramdac *nr;
	int i, j;

	debug("dispwin_clone_ramdac called\n");

	/* Allocate a ramdac */
	if ((nr = (ramdac *)calloc(sizeof(ramdac), 1)) == NULL) {
		return NULL;
	}

	*nr = *r;		/* Structrure copy */

	for (j = 0; j < 3; j++) {

		if ((nr->v[j] = (double *)calloc(sizeof(double), r->nent)) == NULL) {
			for (j--; j >= 0; j--)
				free(nr->v[j]);
			free(nr);
			return NULL;
		}
	}

	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			nr->v[j][i] = r->v[j][i];
		}
	}

	debug("clone is done\n");
	return nr;
}

/* Set the ramdac values to linear */
static void dispwin_setlin_ramdac(ramdac *r) {
	int i, j;

	debug("dispwin_setlin_ramdac called\n");

	for (i = 0; i < r->nent; i++) {
		double val = i/(r->nent - 1.0);
		for (j = 0; j < 3; j++) {
			r->v[j][i] = val;
		}
	}
}

/* We're done with a ramdac structure */
static void dispwin_del_ramdac(ramdac *r) {
	int j;

	debug("dispwin_del_ramdac called\n");

	for (j = 0; j < 3; j++)
		free(r->v[j]);

	free(r);
}

/* ----------------------------------------------- */
/* Change the window color. */
/* Return 1 on error, 2 on window being closed */
static int dispwin_set_color(
dispwin *p,
double r, double g, double b	/* Color values 0.0 - 1.0 */
) {
	int j;
	double val[3];

	debug("dispwin_set_color called\n");

	val[0] = r;
	val[1] = g;
	val[2] = b;

	for (j = 0; j < 3; j++) {
		if (val[j] < 0.0)
			val[j] = 0.0;
		else if (val[j] > 1.0)
			val[j] = 1.0;
	}

	if (p->donat) {		/* Output high precision native using RAMDAC */
		double prange = p->r->nent - 1.0;

		for (j = 0; j < 3; j++) {
			int tt;

//printf("~1 %d: in %f, ",j,val[j]);
			tt = (int)(val[j] * prange);
			p->r->v[j][tt] = val[j];
			val[j] = (double)tt/prange;
//printf(" cell[%d], val %f\n",tt, val[j]);
		}
		if (p->set_ramdac(p,p->r)) {
			debug("set_ramdac() failed\n");
			return 1;
		}
	}

	/* - - - - - - - - - - - - - - */
#ifdef NT
	{
		HDC hdc;
		PAINTSTRUCT ps;
		RECT rect;
		HBRUSH hbr;
		MSG msg;
		int vali[3];

		/* Force a repaint with the new data */
		if (!InvalidateRgn(p->hwnd,NULL,TRUE)) {
#ifdef DEBUG
			printf("InvalidateRgn failed, lasterr = %d\n",GetLastError());
#endif
			return 1;
		}

		/* Convert to 8 bit color */
		for (j = 0; j < 3; j++)
			vali[j] = (int)(255.0 * val[j] + 0.5);

		hdc = BeginPaint(p->hwnd, &ps);
		SaveDC(hdc);
//	hdc = CreateDC("DISPLAY", NULL, NULL, NULL);

		/* Try and turn ICM off */
#ifdef ICM_DONE_OUTSIDEDC
		if (!SetICMMode(hdc, ICM_DONE_OUTSIDEDC)) {
			/* This seems to fail with "invalid handle" under NT4 */
			/* Does it work under Win98 or Win2K ? */
			printf("SetICMMode failed, lasterr = %d\n",GetLastError());
		}
#endif

		hbr = CreateSolidBrush(RGB(vali[0],vali[1],vali[2]));
		SelectObject(hdc,hbr);

		GetClientRect(p->hwnd, &rect);
		FillRect(hdc, &rect, hbr);

		RestoreDC(hdc,-1);
		EndPaint(p->hwnd, &ps);
		
		UpdateWindow(p->hwnd);		/* Flush the paint */
 
		/* Process any pending messages */
		while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
//			GetMessage(&msg, NULL, 0, 0);
//printf("~1 got message %d\n",msg.message);
			if (msg.message == WM_SYSCOMMAND
			 && (msg.wParam == SC_SCREENSAVE
			  || msg.wParam == SC_MONITORPOWER)) {
//printf("~1 was SCREENSAVE/MONITORPOWER\n");
			} else {
				TranslateMessage(&msg);
				DispatchMessage(&msg);
			}
		}

		DeleteDC(hdc);

#ifdef NEVER	/* Only Win98 & 2K */
		/* Reset the system display idle timer */
		SetThreadExecutionState(ES_DISPLAY_REQUIRED);
#endif /* NEVER */
	}
#endif /* NT */

	/* - - - - - - - - - - - - - - */

#ifdef __APPLE__

	/* Set the color */
	p->rgb[0] = val[0];
	p->rgb[1] = val[1];
	p->rgb[2] = val[2];

	/* Cause window repaint with the new data */
	{
		OSStatus stat;
		Rect wrect;
		GetPortBounds(p->port, &wrect);
		if ((stat = InvalWindowRect(p->mywindow, &wrect)) != noErr) {
#ifdef DEBUG
			printf("InvalWindowRect failed with %d\n",stat);
#endif
			return 1;
		}
	}

	/* Make sure our window is brought to the front at least once, */
	/* but not every time, in case the user wants to kill the application. */
	if (p->btf == 0){
		OSStatus stat;
		ProcessSerialNumber cpsn;
		if ((stat = GetCurrentProcess(&cpsn)) != noErr) {
#ifdef DEBUG
			printf("GetCurrentProcess returned error %d\n",stat);
#endif
		} else {
		if ((stat = SetFrontProcess(&cpsn)) != noErr) {
#ifdef DEBUG
				printf("SetFrontProcess returned error %d\n",stat);
#endif
			}
		}
		p->btf = 1;
	}

	/* Run the event loop until refreshed */
	RunApplicationEventLoop();

	if (p->winclose) {
		return 2;
	}

#endif /* __APPLE__ */

	/* - - - - - - - - - - - - - - */

#if defined(UNIX) && !defined(__APPLE__)
	{
		XWindowAttributes mywattributes;
		Colormap mycmap;
		XColor col;
		int vali[3];

		/* Convert to 16 bit color */
		for (j = 0; j < 3; j++)
			vali[j] = (int)(65535.0 * val[j] + 0.5);

		mycmap = DefaultColormap(p->mydisplay, 0);
		col.red = vali[0];
		col.green = vali[1];
		col.blue = vali[2];
		XAllocColor(p->mydisplay, mycmap, &col);
		XSetForeground(p->mydisplay, p->mygc, col.pixel);

		XGetWindowAttributes(p->mydisplay, p->mywindow, &mywattributes);

		XFillRectangle(p->mydisplay, p->mywindow, p->mygc,
		               0, 0, mywattributes.width, mywattributes.height);

		XFlush(p->mydisplay);		/* Make sure it happens */
	}
#endif /* UNIX */
	return 0;
}

/* ----------------------------------------------- */
/* Destroy ourselves */
static void dispwin_del(
dispwin *p
) {

	debug("dispwin_del called\n");

	/* Restore original RAMDAC if we were in native mode */
	if (p->donat && p->or != NULL) {
		p->set_ramdac(p, p->or);
		p->or->del(p->or);
		p->r->del(p->r);
		debug("Restored original ramdac\n");
	}

#ifdef NT
	if (!DestroyWindow(p->hwnd)) {
#ifdef DEBUG
		printf("DestroyWindow failed, lasterr = %d\n",GetLastError());
#endif
	}

	if (!UnregisterClass(p->AppName, NULL)) {
#ifdef DEBUG
		printf("UnregisterClass failed, lasterr = %d\n",GetLastError());
#endif
	}
#endif /* NT */

#ifdef __APPLE__

	if (p->mygc != NULL)
		QDEndCGContext(p->port, &p->mygc);		/* Ignore any errors */
	if (p->mywindow != NULL)
		DisposeWindow(p->mywindow);	

	CGDisplayShowCursor(p->ddid);
//	CGDisplayRelease(p->ddid);

#endif /* __APPLE__ */

#if defined(UNIX) && !defined(__APPLE__)
	debug("About to close display\n");

#ifndef DONT_MESS_WITH_SCREENSAVER
	/* Restore the screensaver state */
	XSetScreenSaver(p->mydisplay, p->timeout, p->interval,
	                p->prefer_blanking, p->allow_exposures);
#endif /* DONT_MESS_WITH_SCREENSAVER */
	
	XFreeGC(p->mydisplay, p->mygc);
	XDestroyWindow(p->mydisplay, p->mywindow);
	XCloseDisplay(p->mydisplay);
	debug("finished\n");
#endif /* UNIX */

	free(p);
}

/* ----------------------------------------------- */
/* Event handler callbacks */

#ifdef NT
	/* None */
#endif /* NT */

#ifdef __APPLE__

/* The OSX event handler */
pascal OSStatus HandleEvent(
EventHandlerCallRef nextHandler,
EventRef theEvent,
void* userData
) {
	dispwin *p = (dispwin *)userData;
	OSStatus result = eventNotHandledErr;
	
	switch (GetEventClass(theEvent)) {
		case kEventClassWindow: {
			switch (GetEventKind(theEvent)) {
				case kEventWindowClose:
					debug("Event: Close Window, exiting event loop\n");
					p->winclose = 1;			/* Mark the window closed */
					QuitApplicationEventLoop();
					result = noErr;
					break;
				
				case kEventWindowBoundsChanged: {
					OSStatus stat;
					Rect wrect;
					debug("Event: Bounds Changed\n");
					GetPortBounds(p->port, &wrect);
					if ((stat = InvalWindowRect(p->mywindow, &wrect)) != noErr) {
#ifdef DEBUG
						printf("InvalWindowRect failed with %d\n",stat);
#endif
					}
					break;
				}
				case kEventWindowDrawContent: {
					OSStatus stat;
					Rect rect;
					CGRect frect;

					debug("Event: Draw Content\n");
					GetPortBounds(p->port, &rect);		/* Bounds is inclusive, global coords */
					frect = CGRectMake(0.0, 0.0,
					  (float)(1.0 + rect.right - rect.left), (float)(1.0 + rect.bottom - rect.top));
					CGContextSetRGBFillColor(p->mygc, p->rgb[0], p->rgb[1], p->rgb[2], 1.0);
					CGContextFillRect(p->mygc, frect);
					CGContextFlush(p->mygc);		/* Force draw to display */
					QuitApplicationEventLoop();		/* Break out of event loop after draw */
					SelectWindow(p->mywindow);
					result = noErr;
					break;
				}
			}
		}
	}

	/* Call next hander if not handled here */
	if (result == eventNotHandledErr) {
		result =  CallNextEventHandler(nextHandler, theEvent);	/* Close window etc */
	}

	return result;
}

#endif /* __APPLE__ */

#if defined(UNIX) && !defined(__APPLE__)
	/* None */
#endif /* UNIX */

/* ----------------------------------------------- */
/* Create a display window, default color of white */
dispwin *new_dispwin(
double width, double height,	/* Width and height in mm */
double hoff, double voff,		/* Offset from center in fraction of screen, range -1.0 .. 1.0 */
int native,						/* NZ if ramdac should be bypassed rather than used. */
int override					/* NZ if override_redirect is to be used on X11 */
) {
	dispwin *p = NULL;

	debug("new_dispwin called\n");

	if ((p = (dispwin *)calloc(sizeof(dispwin), 1)) == NULL)
		return NULL;

	p->donat = native;
	p->get_ramdac   = dispwin_get_ramdac;
	p->set_ramdac   = dispwin_set_ramdac;
	p->set_color    = dispwin_set_color;
	p->del          = dispwin_del;

	/* Basic object is initialised, so create a window */

#ifdef NT
	{
		MSG msg;
		WNDCLASS wc;
		HDC hdc;
		DWORD ctid;					/* Current thread ID */
		int disp_hsz, disp_vsz;		/* Display horizontal/vertical size in mm */
		int disp_hrz, disp_vrz;		/* Display horizontal/vertical resolution in pixels */
		int wi, he;					/* Width and height of window in pixels */
		int ho, vo;					/* Horizontal and vertical offset from center in pixels */
		int bpp;
		
		p->AppName = "Argyll Window";
		
		/* Fill in window class structure with parameters that describe the */
		/* main window. */

		wc.style         = CS_SAVEBITS;	/* Class style(s). */
		wc.lpfnWndProc   = DefWindowProc; /* Function to retrieve messages for class windows */
		wc.cbClsExtra    = 0;			/* No per-class extra data. */
		wc.cbWndExtra    = 0;			/* No per-window extra data. */
		wc.hInstance     = NULL;		/* Application that owns the class. */
		wc.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
		wc.hCursor       = LoadCursor(NULL, IDC_ARROW);
		wc.hbrBackground = GetStockObject(WHITE_BRUSH);
		wc.lpszMenuName  = NULL;
		wc.lpszClassName = p->AppName;

		p->arv = RegisterClass(&wc);

		if (!p->arv) {
#ifdef DEBUG
			printf("RegisterClass failed, lasterr = %d\n",GetLastError());
#endif
			free(p);
			return NULL;
		}

		/* Check out EnumDesktops() for selection between multiple monitors ?? */
		/* as well as EnumDisplayDevices() */

		/* Get device context to main display */
		hdc = CreateDC("DISPLAY", NULL, NULL, NULL);
		disp_hsz = GetDeviceCaps(hdc, HORZSIZE);	/* mm */
		disp_vsz = GetDeviceCaps(hdc, VERTSIZE);	/* pixels */
		disp_hrz = GetDeviceCaps(hdc, HORZRES);
		disp_vrz = GetDeviceCaps(hdc, VERTRES);
		wi = (int)(width * (double)disp_hrz/(double)disp_hsz + 0.5);
		he = (int)(height * (double)disp_vrz/(double)disp_vsz + 0.5);
		ho = (int)(hoff * 0.5 * (disp_hrz - wi) + 0.5);
		vo = (int)(voff * 0.5 * (disp_vrz - he) + 0.5);

		/* It's a bit difficult to know how windows defines the display */
		/* depth. Microsofts doco is fuzzy, and typical values */
		/* for BITSPIXEL and PLANES are confusing (What does "32" and "1" */
		/* mean ?) The doco for COLORRES is also fuzzy, but it returns */
		/* a meaningful number on my box (24) */
		bpp = GetDeviceCaps(hdc, COLORRES);
		if (bpp <= 0)
			p->pdepth = 8;		/* Assume this is so */
		else
			p->pdepth = bpp/3;

		DeleteDC(hdc);

		p->hwnd = CreateWindowEx(
			WS_EX_PALETTEWINDOW,
			p->AppName,
			"Argyll Display Calibration Window",
			WS_DISABLED | WS_CAPTION | WS_THICKFRAME, 
			(disp_hrz - wi)/2 + ho,	/* Position of top left */
			(disp_vrz - he)/2 + vo,
			wi, he,					/* Size */
			NULL,					/* Handle to parent or owner */
			NULL,					/* Handle to menu or child window */
			NULL, 					/* hInstance Handle to appication instance */
			NULL);					/* pointer to window creation data */

		if (!p->hwnd) {
#ifdef DEBUG
			printf("CreateWindow failed, lasterr = %d\n",GetLastError());
#endif
			/* Should UnregisterClass(p->arv, ??) ?? */
			free(p);
			return NULL;
		}
		
		ShowWindow(p->hwnd, SW_SHOWNA);

		/* Process any pending messages */
		while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
			GetMessage(&msg, NULL, 0, 0);
//printf("~1 got message %d\n",msg.message);
			if (msg.message == WM_SYSCOMMAND
			 && (msg.wParam == SC_SCREENSAVE
			  || msg.wParam == SC_MONITORPOWER)) {
//printf("~1 was SCREENSAVE/MONITORPOWER\n");
			} else {
				TranslateMessage(&msg);
				DispatchMessage(&msg);
			}
		}
	}
#endif /* NT */

#ifdef __APPLE__
	{
		OSStatus stat;
		EventHandlerRef fHandler;		/* Event handler reference */
		const EventTypeSpec	kEvents[] =
		{ 
			{ kEventClassWindow, kEventWindowClose },
			{ kEventClassWindow, kEventWindowBoundsChanged },
			{ kEventClassWindow, kEventWindowDrawContent }
		};

		int wi, he;				/* Width and height in pixels */
		int ho, vo;				/* Horizontal and vertical offset from center in pixels */
		int xo, yo;				/* Window location */
	    Rect	wRect;
		WindowClass wclass = 0;
		WindowAttributes attr = kWindowNoAttributes;

		/* Choose the windows class and attributes */
//		wclass = kAlertWindowClass;				/* Above everything else */
//		wclass = kModalWindowClass;
//		wclass = kFloatingWindowClass;
		wclass = kUtilityWindowClass;
		attr |= kWindowDoesNotCycleAttribute;
		attr |= kWindowNoActivatesAttribute;
//		attr |= kWindowStandardFloatingAttributes;
		attr |= kWindowStandardHandlerAttribute;
								/* This window has the standard Carbon Handler */
		attr |= kWindowNoShadowAttribute;

		/* Using Quartz services to discover displays */
		/* Could capture the display with CGDisplayCapture(), and then get a */
		/* GC using CGDisplayGetDrawingContext(), and hide the cursor with */
		/* CGDisplayHideCursor(), CGDisplayShowCursor() to ensure nothing */
		/* interferes with the profiling ?? - doesn't work quite that simply, */
		/* and how would the user interact with the application ? */
		{
			int target = 0;			/* Display device we are targeting */
			int i;
			CGDisplayErr dstat;
			CGDisplayCount dcount;		/* Number of display IDs */
			CGDirectDisplayID *dids;	/* Array of display IDs */
			io_connect_t dport;
			CFDictionaryRef dr;
			CFNumberRef value;			/* Temp value */
			CGRect dbound;				/* Bounding rectangle of chosen display */
			int wmm, hmm;				/* Width and height in mm */

			if ((dstat = CGGetActiveDisplayList(0, NULL, &dcount)) != kCGErrorSuccess) {
#ifdef DEBUG
				printf("CGGetActiveDisplayList #1 returned error %d\n",dstat);
#endif
				dispwin_del(p);
				return NULL;
			}
			if ((dids = (CGDirectDisplayID *)malloc(dcount * sizeof(CGDirectDisplayID))) == NULL) {
#ifdef DEBUG
				printf("malloc of %d CGDirectDisplayID's failed\n",dcount);
#endif
				dispwin_del(p);
				return NULL;
			}
			if ((dstat = CGGetActiveDisplayList(dcount, dids, &dcount)) != kCGErrorSuccess) {
#ifdef DEBUG
				printf("CGGetActiveDisplayList #2 returned error %d\n",dstat);
#endif
				free(dids);
				dispwin_del(p);
				return NULL;
			}

#ifdef DEBUG
			printf("Found %d displays, choosing index %d\n",dcount,target);
#endif

#ifdef NEVER
			/* Got displays, now have a look through them */
			for (i = 0; i < dcount; i++) {
				GDHandle gdh;
				GDPtr gdp;

				/* Dump display mode dictionary */
				CFIndex nde, j;
				CFDictionaryRef dr;
				void **keys, **values;
				
				dr = CGDisplayCurrentMode(dids[i]);
				nde = CFDictionaryGetCount(dr);

				printf("Dict contains %d entries \n", nde);
				if ((keys = (void **)malloc(nde * sizeof(void *))) == NULL) {
#ifdef DEBUG
					printf("malloc failed for disp mode keys\n");
#endif
					free(dids);
					dispwin_del(p);
					return NULL;
				}
				if ((values = (void **)malloc(nde * sizeof(void *))) == NULL) {
#ifdef DEBUG
					printf("malloc failed for disp mode values\n");
#endif
					free(keys);
					free(dids);
					dispwin_del(p);
					return NULL;
				}
				CFDictionaryGetKeysAndValues(dr, (const void **)keys, (const void **)values);
				for (j = 0; j < nde; j++) {
					printf("Entry %d key = %s\n", j, CFStringGetCStringPtr(keys[j], kCFStringEncodingMacRoman));
				}
				free(values);
				free(keys);
			}
#endif /* NEVER */

			/* Choose the display we want to use */
			if (target > (dcount-1))
				target = dcount-1;
			p->ddid = dids[target];
			free(dids);

			/* Get enough information about the display to position */
			/* the test window in the middle */
			dbound = CGDisplayBounds(p->ddid);
			
			/* We need to figure the dpi of the display to compute a window size. */
			if ((dport = CGDisplayIOServicePort(p->ddid)) == MACH_PORT_NULL) {
#ifdef DEBUG
				printf("CGDisplayIOServicePort returned NULL\n");
#endif
				dispwin_del(p);
				return NULL;
			}
			if ((dr = IODisplayCreateInfoDictionary(dport, 0)) == NULL) {
#ifdef DEBUG
				printf("IOCreateDisplayInfoDictionary returned NULL\n");
#endif
				dispwin_del(p);
				return NULL;
			}
			if ((value = CFDictionaryGetValue(dr, CFSTR(kDisplayHorizontalImageSize))) == NULL
			  || CFGetTypeID(value) != CFNumberGetTypeID()) {
#ifdef DEBUG
				printf("Get DisplayHorizontalImageSize failed\n");
#endif
				CFRelease(dr);
				dispwin_del(p);
				return NULL;
			}
			CFNumberGetValue(value, kCFNumberIntType, &wmm);

			if ((value = CFDictionaryGetValue(dr, CFSTR(kDisplayVerticalImageSize))) == NULL
			  || CFGetTypeID(value) != CFNumberGetTypeID()) {
#ifdef DEBUG
				printf("Get DisplayVerticalImageSize failed\n");
#endif
				CFRelease(dr);
				dispwin_del(p);
				return NULL;
			}
			CFNumberGetValue(value, kCFNumberIntType, &hmm);
			CFRelease(dr);

//printf("~1 Display size = %d x %d mm\n",wmm,hmm);
			p->pdepth = CGDisplayBitsPerSample(p->ddid);

			wi = (int)(width * dbound.size.width/wmm + 0.5);
			he = (int)(height * dbound.size.height/hmm + 0.5);
			xo = (int)(0.5 * (dbound.size.width - wi) + 0.5);
			yo = (int)(0.5 * (dbound.size.height - he) + 0.5);
			ho = (int)(hoff * 0.5 * (dbound.size.width - wi) + 0.5);
			vo = (int)(voff * 0.5 * (dbound.size.height - he) + 0.5);
			xo += ho;
			yo += vo;
		}

//printf("~1 Got width size %d x %d at %d %d from %fmm x %fmm\n",wi,he,xo,yo,width,height);
	    SetRect(&wRect,xo,yo,xo+wi,yo+he); /* left, top, right, bottom */

		/* Create invisible new window of given class, attributes and size */
	    stat = CreateNewWindow(wclass, attr, &wRect, &p->mywindow);
		if (stat != noErr || p->mywindow == nil) {
#ifdef DEBUG
			printf("CreateNewWindow failed with code %d\n",stat);
#endif
			dispwin_del(p);
			return NULL;
		}

		/* Set a title */
		SetWindowTitleWithCFString(p->mywindow, CFSTR("Argyll Window"));

		/* Install the event handler */
	    stat = InstallEventHandler(
			GetWindowEventTarget(p->mywindow),	/* Objects events to handle */
			NewEventHandlerUPP(HandleEvent),	/* Handler to call */
			sizeof(kEvents)/sizeof(EventTypeSpec),	/* Number in type list */
			kEvents,						/* Table of events to handle */
			(void *)p,						/* User Data is our object */
			&fHandler						/* Event handler reference return value */
		);
		if (stat != noErr) {
#ifdef DEBUG
			printf("InstallEventHandler failed with code %d\n",stat);
#endif
			dispwin_del(p);
			return NULL;
		}

		/* Create a GC */
		p->port = GetWindowPort(p->mywindow);

		if ((stat = QDBeginCGContext(p->port, &p->mygc)) != noErr) {
#ifdef DEBUG
			printf("QDBeginCGContext returned error %d\n",stat);
#endif
			dispwin_del(p);
			return NULL;
		}
		p->winclose = 0;
		p->rgb[0] = p->rgb[0] = p->rgb[0] = 1.0;	/* Set White */

		/* Activate the window */

#ifdef NEVER
		/* Make window pop to the top */
		{
			ProcessSerialNumber cpsn;
			if ((stat = GetCurrentProcess(&cpsn)) != noErr) {
#ifdef DEBUG
				printf("GetCurrentProcess returned error %d\n",stat);
#endif
			} else {
			if ((stat = SetFrontProcess(&cpsn)) != noErr) {
#ifdef DEBUG
					printf("SetFrontProcess returned error %d\n",stat);
#endif
				}
			}
		}
#endif /* NEVER */
// ~~99
//		BringToFront(p->mywindow);
		ShowWindow(p->mywindow);	/* Makes visible and triggers update event */
//		ActivateWindow(p->mywindow, false);
		SendBehind(p->mywindow, NULL);
//		SelectWindow(p->mywindow);

//printf("~1 created window, about to enter event loop\n");
		/* Run the event loop ?? */
		/* Seems to hang because draw happens before event loop gets run ???/ */
		/* RunApplicationEventLoop(); */

//		CGDisplayHideCursor(p->ddid);
//		CGDisplayCapture(p->ddid);

//printf("~1 done create window\n");
		if (p->winclose) {
			dispwin_del(p);
			return NULL;
		}
	}
#endif /* __APPLE__ */


#if defined(UNIX) && !defined(__APPLE__)
	{
		/* NOTE: That we're not doing much to detect if the display/window
		   we open is unsuitable for high quality color (ie. at least
		   24 bit etc.
		 */

		/* stuff for X windows */
		Window rootwindow;
		char *appname = "TestWin";
		XSetWindowAttributes myattr;
		XWindowAttributes mywattributes;
		XEvent      myevent;
		XTextProperty myappname;
		XSizeHints  mysizehints;
		XWMHints    mywmhints;
		unsigned long myforeground,mybackground;
		int disp_hsz, disp_vsz;		/* Display horizontal/vertical size in mm */
		int disp_hrz, disp_vrz;		/* Display horizontal/vertical resolution in pixels */
		int wi, he;				/* Width and height of window in pixels */
		int ho, vo;				/* Horizontal and vertical offset from center in pixels */
		int rv;
	
		/* open the display */
		p->mydisplay = XOpenDisplay("");
		if(!p->mydisplay) {
			debug("Unable to open display\n");
			dispwin_del(p);
			return NULL;
		}
		debug("Opened display OK\n");

		p->myscreen = DefaultScreen(p->mydisplay);
		rootwindow = DefaultRootWindow(p->mydisplay);

		mybackground = WhitePixel(p->mydisplay,p->myscreen);
		myforeground = BlackPixel(p->mydisplay,p->myscreen);
	
		/* Get device context to main display */
		disp_hsz = DisplayWidthMM(p->mydisplay, p->myscreen);
		disp_vsz = DisplayHeightMM(p->mydisplay, p->myscreen);
		disp_hrz = DisplayWidth(p->mydisplay, p->myscreen);
		disp_vrz = DisplayHeight(p->mydisplay, p->myscreen);

		wi = (int)(width * (double)disp_hrz/(double)disp_hsz + 0.5);
		he = (int)(height * (double)disp_vrz/(double)disp_vsz + 0.5);
		ho = (int)(hoff * 0.5 * (disp_hrz - wi) + 0.5);
		vo = (int)(voff * 0.5 * (disp_vrz - he) + 0.5);

		/* Setup Size Hints */
		mysizehints.flags = PPosition | USSize;
		mysizehints.x = (disp_hrz - wi)/2 + ho;
		mysizehints.y = (disp_vrz - he)/2 + vo;
		mysizehints.width = wi;
		mysizehints.height = he;
	
		/* Setup Window Manager Hints */
		mywmhints.flags = InputHint | StateHint;
		mywmhints.input = 0;
		mywmhints.initial_state = NormalState;

		/* Setup Window Attributes */
		myattr.background_pixel = mybackground;
		myattr.bit_gravity = CenterGravity;
		myattr.win_gravity = CenterGravity;
		myattr.backing_store = WhenMapped;		/* Since we aren't listning to events */
		if (override)
			myattr.override_redirect = True;		/* Takes the WM out of the picture */
		else
			myattr.override_redirect = False;

		debug("Opening window\n");
		p->mywindow = XCreateWindow(
			p->mydisplay, rootwindow,
			mysizehints.x,mysizehints.y,mysizehints.width,mysizehints.height,
			5,							/* Border width */
			CopyFromParent,				/* Depth */
			InputOutput,				/* Class */
			CopyFromParent,				/* Visual */
			CWBackPixel | CWBitGravity	/* Attributes Valumask */
			| CWWinGravity | CWBackingStore | CWOverrideRedirect,
			&myattr						/* Attribute details */
		);

		/* Get the windows attributes */
		if (XGetWindowAttributes(
			p->mydisplay, p->mywindow,
			&mywattributes) == 0) {
			debug("XGetWindowAttributes failed\n");
			dispwin_del(p);
			return NULL;
		}

		p->pdepth = mywattributes.depth/3;

		/* Setup TextProperty */
		XStringListToTextProperty(&appname, 1, &myappname);

		XSetWMProperties(
			p->mydisplay, p->mywindow,
			&myappname,				/* Window name */
			&myappname,				/* Icon name */
			NULL, 0,				/* argv, argc */
			&mysizehints,
			&mywmhints,
			NULL);					/* No class hints */
	
		/* Set aditional properties */
		XChangeProperty(
			p->mydisplay, p->mywindow,
			XA_WM_TRANSIENT_FOR,		/* Property */
			XA_WINDOW,					/* Type */
			sizeof(Window) * 8,			/* format */
			PropModeReplace,			/* Change mode */
			(char *)(&rootwindow),		/* Data is Root Window XID */
			1							/* Number of elements of data */
		);

		p->mygc = XCreateGC(p->mydisplay,p->mywindow,0,0);
		XSetBackground(p->mydisplay,p->mygc,mybackground);
		XSetForeground(p->mydisplay,p->mygc,myforeground);
		
		XSelectInput(p->mydisplay,p->mywindow, ExposureMask);
	
		XMapRaised(p->mydisplay,p->mywindow);
		debug("Raised window\n");
	
#ifndef DONT_MESS_WITH_SCREENSAVER
		/* Save the screensaver state, and then disable it */
		XGetScreenSaver(p->mydisplay, &p->timeout, &p->interval,
		                &p->prefer_blanking, &p->allow_exposures);
		XSetScreenSaver(p->mydisplay, 0, 0, DefaultBlanking, DefaultExposures);
#endif
	
		/* Deal with any pending events */
		debug("About to enter main loop\n");
		while(XCheckIfEvent(p->mydisplay, &myevent, NULL, NULL)) {
			switch(myevent.type) {
				case Expose:
					if(myevent.xexpose.count == 0) {	/* Repare the exposed region */
						XWindowAttributes mywattributes;
						XGetWindowAttributes(p->mydisplay, p->mywindow, &mywattributes);

						XFillRectangle(p->mydisplay, p->mywindow, p->mygc,
						               0, 0, mywattributes.width, mywattributes.height);
					}
					break;
			}
		}
	}

#endif /* UNIX */

	/* Setup for native mode */
	if (p->donat) {
		debug("About to setup native mode\n");
		if ((p->or = p->get_ramdac(p)) == NULL
		 || (p->r = p->or->clone(p->or)) == NULL) {
			debug("Native mode can't work, no RAMDAC support\n");
			dispwin_del(p);
			return NULL;
		}
		p->r->setlin(p->r);
		debug("Saved original RAMDAC\n");
	}

	return p;
}

/* ---------------------------------------------------------------- */
#ifdef STANDALONE_TEST
/* test code */

#include "numlib.h"

#if defined (NT)
#define sleep(secs) Sleep((secs) * 1000)
#endif

static void
usage(void) {
	fprintf(stderr,"Test display patch test, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	fprintf(stderr,"usage: dispwin [options]\n");
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -i              Run forever with random values\n");
	fprintf(stderr," -f              Test grey ramp fade\n");
	fprintf(stderr," -r              Test just Ramdac\n");
	fprintf(stderr," -n              Test native output\n");
	fprintf(stderr," -k file.cal     Load display calibration and exit\n");
	exit(1);
}

/* 32 bit pseudo random sequencer based on XOR feedback */
/* generates number between 1 and 4294967295 */
#define PSRAND32(S) (((S) & 0x80000000) ? (((S) << 1) ^ 0xa398655d) : ((S) << 1))

int 
main(int argc, char *argv[]) {
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;			/* Verbose flag */
	int ramd = 0;			/* Just test ramdac */
	int fade = 0;			/* Test greyramp fade */
	int donat = 0;			/* Test native output */
	int inf = 0;			/* Infnite patches flag */
	char calname[MAXNAMEL+1] = "\000";	/* Calibration file name */
	dispwin *dw;
	unsigned int seed = 0x56781234;
	int i, j;
	ramdac *or, *r;

	error_program = "Dispwin";
	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {

		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage();

			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V')
				verb = 1;

			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I')
				inf = 1;

			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F')
				fade = 1;

			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R')
				ramd = 1;

			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N')
				donat = 1;

			/* Calibration file */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				char *p;
				int tt;
				fa = nfa;
				if (na == NULL) usage();
				strncpy(calname,na,MAXNAMEL); calname[MAXNAMEL] = '\000';
			}
			else
				usage();
		}
		else
			break;
	}


	if (verb)
		printf("About to open window on the display\n");
	if ((dw = new_dispwin(40.0, 40.0, 0.0, 0.0, donat, 1)) == NULL) {
		printf("Error - new_dispwin failed!\n");
		return -1;
	}

	/* Setup a display calibration set if we are given one */
	if (calname[0] != '\000') {
		if ((r = dw->get_ramdac(dw)) != NULL) {
			cgats *ccg;			/* calibration cgats structure */
			int ncal;
			int ii, fi, ri, gi, bi;
			double cal[3][256];					/* Display calibration */
			
			ccg = new_cgats();			/* Create a CGATS structure */
			ccg->add_other(ccg, "CAL"); /* our special calibration type */
		
			if (ccg->read_name(ccg, calname))
				error("CGATS calibration file read error %s on file '%s'",ccg->err,calname);
		
			if (ccg->t[0].tt != tt_other || ccg->t[0].oi != 0)
				error ("Calibration file isn't a CAL format file");
			if (ccg->ntables < 1)
				error ("Calibration file '%s' doesn't contain at least one table",calname);
		
			if ((ncal = ccg->t[0].nsets) <= 0)
				error ("No data in set of file '%s'",calname);
		
			if (ncal != 256)
				error ("Expect 256 data sets in file '%s'",calname);
		
			if ((fi = ccg->find_kword(ccg, 0, "DEVICE_CLASS")) < 0)
				error ("Calibration file '%s' doesn't contain keyword COLOR_REPS",calname);
			if (strcmp(ccg->t[0].kdata[fi],"DISPLAY") != 0)
				error ("Calibration file '%s' doesn't have DEVICE_CLASS of DISPLAY",calname);

			if ((ii = ccg->find_field(ccg, 0, "RGB_I")) < 0)
				error ("Calibration file '%s' doesn't contain field RGB_I",calname);
			if (ccg->t[0].ftype[ii] != r_t)
				error ("Field RGB_R in file '%s' is wrong type",calname);
			if ((ri = ccg->find_field(ccg, 0, "RGB_R")) < 0)
				error ("Calibration file '%s' doesn't contain field RGB_R",calname);
			if (ccg->t[0].ftype[ri] != r_t)
				error ("Field RGB_R in file '%s' is wrong type",calname);
			if ((gi = ccg->find_field(ccg, 0, "RGB_G")) < 0)
				error ("Calibration file '%s' doesn't contain field RGB_G",calname);
			if (ccg->t[0].ftype[gi] != r_t)
				error ("Field RGB_G in file '%s' is wrong type",calname);
			if ((bi = ccg->find_field(ccg, 0, "RGB_B")) < 0)
				error ("Calibration file '%s' doesn't contain field RGB_B",calname);
			if (ccg->t[0].ftype[bi] != r_t)
				error ("Field RGB_B in file '%s' is wrong type",calname);
			for (i = 0; i < ncal; i++) {
				r->v[0][i] = *((double *)ccg->t[0].fdata[i][ri]);
				r->v[1][i] = *((double *)ccg->t[0].fdata[i][gi]);
				r->v[2][i] = *((double *)ccg->t[0].fdata[i][bi]);
			}

			if (verb)
				printf("About to set given calibration\n");
			if (dw->set_ramdac(dw,r)) {
				error("Failed to set ramdac");
			}

			r->del(r);
			ccg->del(ccg);
		} else {
			printf("We don't have access to the RAMDAC\n");
		}
	}

	if (calname[0] == '\000' && ramd == 0) {

		if (fade) {
			int i;
			int steps = 256;
			for (i = 0; i < steps; i++) {
				double tt;
				tt = i/(steps - 1.0);
				dw->set_color(dw, tt, tt, tt);
				/* Need a 20msec delay here */
				printf("Val = %f\n",tt);
			}
		} else {

			printf("Setting Red\n");
			dw->set_color(dw, 1.0, 0.0, 0.0);	/* Red */

			sleep(2);
			printf("Setting Green\n");
			dw->set_color(dw, 0.0, 1.0,  0.0);	/* Green */

			sleep(2);
			printf("Setting Blue\n");
			dw->set_color(dw, 0.0, 0.0, 1.0);	/* Blue */

			sleep(2);
			printf("Setting Cyan\n");
			dw->set_color(dw, 0.0, 1.0, 1.0);	/* Cyan */

			sleep(2);
			printf("Setting Magenta\n");
			dw->set_color(dw, 1.0, 0.0,  1.0);	/* Magenta */

			sleep(2);
			printf("Setting Yellow\n");
			dw->set_color(dw, 1.0, 1.0, 0.0);	/* Yellow */

			sleep(2);

			for (;inf != 0;) {
				double col[3];

				for (i = 0; i < 3; i++) {
					seed = PSRAND32(seed);
					col[i] = seed/4294967295.0;
				}

				printf("Setting %f %f %f\n",col[0],col[1],col[2]);
				dw->set_color(dw, col[0],col[1],col[2]);
				sleep(2);
			}
		}
	}

	/* Test out the RAMDAC access */
	if (calname[0] == '\000') {
		if ((or = dw->get_ramdac(dw)) != NULL) {
			
			r = or->clone(or);

			/* Try darkening it */
			for (j = 0; j < 3; j++) {
				for (i = 0; i < r->nent; i++) {
					r->v[j][i] = 0.5 * or->v[j][i];
				}
			}
			printf("Darkening screen\n");
			if (dw->set_ramdac(dw,r)) {
				error("Failed to set ramdac");
			}
			sleep(1);

			/* Try lightening it */
			for (j = 0; j < 3; j++) {
				for (i = 0; i < r->nent; i++) {
					double tt;
					tt = 2.0 * or->v[j][i];
					if (tt > 1.0)
						tt = 1.0;
					r->v[j][i] = tt;
				}
			}
			printf("Lightening screen\n");
			if (dw->set_ramdac(dw,r)) {
				error("Failed to set ramdac");
			}
			sleep(1);

			/* restor it */
			printf("Restoring screen\n");
			if (dw->set_ramdac(dw,or)) {
				error("Failed to set ramdac");
			}

			r->del(r);
			or->del(or);

		} else {
			printf("We don't have access to the RAMDAC\n");
		}
	}
	
	if (verb)
		printf("About to destroy window\n");

	dw->del(dw);

	return 0;
}

#endif /* STANDALONE_TEST */
/* ---------------------------------------------------------------- */

