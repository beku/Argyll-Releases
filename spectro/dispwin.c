
/* 
 * Argyll Color Correction System
 * Display target patch window
 *
 * Author: Graeme W. Gill
 * Date:   4/10/96
 *
 * Copyright 1998 - 2007, Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/* This program displays test patches on a WinNT, MAC OSX or X11 windowing system. */

/*
   We are rather rough with how we handle window messages. We should
   really start another thread/process to handle the messages, rather
   that only servicing messages when we feel like it. It seems to
   work OK for test patches. 
*/

/* TTBD
 *
 *	We could improve the robustnes against LCD warm up
 *	effects by reading the white every N readings,
 *	and then linearly interpolating the white readings,
 *	and scaling them back to a reference white.
 *
 * Should probably check the display attributes (like visual depth)
 * and complain if we aren't using 24 bit color or better. 
 *
 */

#ifdef __MINGW32__
# define WINVER 0x0500
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifndef NT
#include <unistd.h>
#endif
#include <time.h>
#include "copyright.h"
#include "config.h"
#include "icc.h"
#include "numsup.h"
#include "cgats.h"
#include "dispwin.h"

#define VERIFY_TOL (1.0/255.0)

#undef DEBUG
//#define STANDALONE_TEST

#ifdef DEBUG
# define errout stderr
# define debug(xx)	fprintf(stderr, xx )
# define debug2(xx)	fprintf xx
#else
# define debug(xx)
# define debug2(xx)
#endif

/* ----------------------------------------------- */
/* Dealing with locating displays */

#ifdef NT

#define sleep(secs) Sleep((secs) * 1000)

static BOOL CALLBACK MonitorEnumProc(
  HMONITOR hMonitor,  /* handle to display monitor */
  HDC hdcMonitor,     /* NULL, because EnumDisplayMonitors hdc is NULL */
  LPRECT lprcMonitor, /* Virtual screen coordinates of this monitor */
  LPARAM dwData       /* Context data */
) {
	disppath ***pdisps = (disppath ***)dwData;
	disppath **disps = *pdisps;
	MONITORINFOEX pmi;
	int ndisps = 0;
	
	/* Add the display to the list */
	if (disps == NULL) {
		if ((disps = (disppath **)calloc(sizeof(disppath *), 1 + 1)) == NULL) {
			debug("get_displays failed on malloc\n");
			return FALSE;
		}
	} else {
		/* Count current number on list */
		for (ndisps = 0; disps[ndisps] != NULL; ndisps++)
			;
		if ((disps = (disppath **)realloc(disps,
		                     sizeof(disppath *) * (ndisps + 2))) == NULL) {
			debug("get_displays failed on malloc\n");
			return FALSE;
		}
		disps[ndisps+1] = NULL;	/* End marker */
	}

	if ((disps[ndisps] = calloc(sizeof(disppath),1)) == NULL) {
		debug("get_displays failed on malloc\n");
		return FALSE;
	}

	pmi.cbSize = sizeof(MONITORINFOEX);
	if (GetMonitorInfo(hMonitor, (MONITORINFO *)&pmi) == 0) {
		debug("get_displays failed GetMonitorInfo\n");
		return FALSE;
	}

	strcpy(disps[ndisps]->name, pmi.szDevice);
	disps[ndisps]->description = NULL;

	disps[ndisps]->sx = lprcMonitor->left;
	disps[ndisps]->sy = lprcMonitor->top;
	disps[ndisps]->sw = lprcMonitor->right - lprcMonitor->left;
	disps[ndisps]->sh = lprcMonitor->bottom - lprcMonitor->top;

	*pdisps = disps;
	return TRUE;
}

#endif /* NT */


#if defined(UNIX) && !defined(__APPLE__)
/* A noop X11 error handler */
int null_error_handler(Display *disp, XErrorEvent *ev) {
	return 0;
}
#endif	/* X11 */

/* Return pointer to list of disppath. Last will be NULL. */
/* Return NULL on failure. Call free_disppaths() to free up allocation */
disppath **get_displays() {
	disppath **disps = NULL;

#ifdef NT
	BOOL (WINAPI* pEnumDisplayDevices)(PVOID,DWORD,PVOID,DWORD);
	DISPLAY_DEVICE dd;
	char buf[200];
	int i, j;

	/* EnumDisplayDevicesA was left out of lib32.lib on earlier SDK's ... */
//	(FARPROC)pEnumDisplayDevices = GetProcAddress(LoadLibrary("USER32"), "EnumDisplayDevicesA");
	pEnumDisplayDevices = (BOOL (WINAPI*)(PVOID,DWORD,PVOID,DWORD)) GetProcAddress(LoadLibrary("USER32"), "EnumDisplayDevicesA");

	if (EnumDisplayMonitors(NULL, NULL, MonitorEnumProc, (LPARAM)&disps) == 0) {
		free_disppaths(disps);
		debug("EnumDisplayMonitors failed\n");
		return NULL;
	}

	/* Now locate the displays */
	for (i = 0; ; i++) {
		if (disps[i] == NULL)
			break;

		/* Find the matching device, and get further information */
		dd.cb = sizeof(DISPLAY_DEVICE);
		for (j = 0; ; j++) {
			if ((*pEnumDisplayDevices)(NULL, j, &dd, 0) == 0)
				break;
			/* (Could add dd.DeviceString, which is the graphics card name) */
			if (strcmp(disps[i]->name, dd.DeviceName) == 0) {
				sprintf(buf,"%s, at %d, %d, width %d, height %d%s",dd.DeviceName+4,
				        disps[i]->sx, disps[i]->sy, disps[i]->sw, disps[i]->sh,
				        dd.StateFlags & DISPLAY_DEVICE_PRIMARY_DEVICE ? " (Primary Display)" : "");

				if ((disps[i]->description = strdup(buf)) == NULL) {
					debug("get_displays failed on malloc\n");
					free_disppaths(disps);
					return NULL;
				}
			}
		}
		if (disps[i]->description == NULL) {
			if ((disps[i]->description = strdup("Unknown")) == NULL) {
				free_disppaths(disps);
				debug("get_displays malloc failed\n");
				return NULL;
			}
		}
	}
#endif /* NT */

#ifdef __APPLE__
	int i;
	CGDisplayErr dstat;
	CGDisplayCount dcount;		/* Number of display IDs */
	CGDirectDisplayID *dids;	/* Array of display IDs */

	if ((dstat = CGGetActiveDisplayList(0, NULL, &dcount)) != kCGErrorSuccess || dcount < 1) {
		debug("CGGetActiveDisplayList #1 returned error\n");
		return NULL;
	}
	if ((dids = (CGDirectDisplayID *)malloc(dcount * sizeof(CGDirectDisplayID))) == NULL) {
		debug("malloc of CGDirectDisplayID's failed\n");
		return NULL;
	}
	if ((dstat = CGGetActiveDisplayList(dcount, dids, &dcount)) != kCGErrorSuccess) {
		debug("CGGetActiveDisplayList #2 returned error\n");
		free(dids);
		return NULL;
	}

	/* Found dcount displays */
	debug2((errout,"Found %d screens\n",dcount));

	/* Allocate our list */
	if ((disps = (disppath **)calloc(sizeof(disppath *), dcount + 1)) == NULL) {
		debug("get_displays failed on malloc\n");
		free(dids);
		return NULL;
	}
	for (i = 0; i < dcount; i++) {
		if ((disps[i] = calloc(sizeof(disppath), 1)) == NULL) {
			debug("get_displays failed on malloc\n");
			free_disppaths(disps);
			free(dids);
			return NULL;
		}
		disps[i]->ddid = dids[i];
	}

	/* Got displays, now figure out a description for each one */
	for (i = 0; i < dcount; i++) {
		CGRect dbound;				/* Bounding rectangle of chosen display */
		io_service_t dport;
		CFDictionaryRef ddr, pndr;
		CFIndex dcount;
		char *dp = NULL, desc[50];
		char buf[200];

		dbound = CGDisplayBounds(dids[i]);
		disps[i]->sx = dbound.origin.x;
		disps[i]->sy = dbound.origin.y;
		disps[i]->sw = dbound.size.width;
		disps[i]->sh = dbound.size.height;
			
		/* Try and get some information about the display */
		if ((dport = CGDisplayIOServicePort(dids[i])) == MACH_PORT_NULL) {
			debug("CGDisplayIOServicePort returned error\n");
			free_disppaths(disps);
			free(dids);
			return NULL;
		}

#ifdef NEVER
		{
			io_name_t name;
			if (IORegistryEntryGetName(dport, name) != KERN_SUCCESS) {
				debug("IORegistryEntryGetName returned error\n");
				free_disppaths(disps);
				free(dids);
				return NULL;
			}
			printf("Driver %d name = '%s'\n",i,name);
		}
#endif
		if ((ddr = IODisplayCreateInfoDictionary(dport, 0)) == NULL) {
			debug("IODisplayCreateInfoDictionary returned NULL\n");
			free_disppaths(disps);
			free(dids);
			return NULL;
		}
		if ((pndr = CFDictionaryGetValue(ddr, CFSTR(kDisplayProductName))) == NULL) {
			debug("CFDictionaryGetValue returned NULL\n");
			CFRelease(ddr);
			free_disppaths(disps);
			free(dids);
			return NULL;
		}
		if ((dcount = CFDictionaryGetCount(pndr)) > 0) {
			const void **keys;
			const void **values;
			int j;

			keys = (const void **)calloc(sizeof(void *), dcount);
			values = (const void **)calloc(sizeof(void *), dcount);
			if (keys == NULL || values == NULL) {
				if (keys != NULL)
					free(keys);
				if (values != NULL)
					free(values);
				debug("malloc failed\n");
				CFRelease(ddr);
				free_disppaths(disps);
				free(dids);
				return NULL;
			}
			CFDictionaryGetKeysAndValues(pndr, keys, values);
			for (j = 0; j < dcount; j++) {
				const char *k, *v;
				char kbuf[50], vbuf[50];
				k = CFStringGetCStringPtr(keys[j], kCFStringEncodingMacRoman);
				if (k == NULL) {
					if (CFStringGetCString(keys[j], kbuf, 50, kCFStringEncodingMacRoman))
						k = kbuf;
				}
				v = CFStringGetCStringPtr(values[j], kCFStringEncodingMacRoman);
				if (v == NULL) {
					if (CFStringGetCString(values[j], vbuf, 50, kCFStringEncodingMacRoman))
						v = vbuf;
				}
//printf("~1 got key %s and value %s\n",k,v);
				/* We're only grabing the english description... */
				if (k != NULL && v != NULL && strcmp(k, "en_US") == 0) {
					strncpy(desc, v, 49);
					desc[49] = '\000';
					dp = desc;
				}
			}
			free(keys);
			free(values);
		}
		CFRelease(ddr);

		if (dp == NULL) {
			strcpy(desc, "(unknown)");
			dp = desc;
		}
		sprintf(buf,"%s, at %d, %d, width %d, height %d%s",dp,
	        disps[i]->sx, disps[i]->sy, disps[i]->sw, disps[i]->sh,
	        CGDisplayIsMain(dids[i]) ? " (Primary Display)" : "");

		if ((disps[i]->description = strdup(buf)) == NULL) {
			debug("get_displays failed on malloc\n");
			free_disppaths(disps);
			free(dids);
			return NULL;
		}
	}

	free(dids);
#endif /* __APPLE__ */

#if defined(UNIX) && !defined(__APPLE__)
	int i, j;
	int dcount;		/* Number of screens */
	char *dname;
	char buf[100];
	int evb = 0, erb = 0;
	Display *mydisplay;
	disppath *tdispp;
	XineramaScreenInfo *xai = NULL;

	/* There seems to be no way of getting the available displays */
	/* on an X11 system. Attempting to open them in sequence */
	/* takes too long. We just rely on the user supplying the */
	/* right display. We can enumerate screens though. */
	if ((dname = getenv("DISPLAY")) != NULL) {
		strncpy(buf,dname,99); buf[99] = '\000';
	} else
		strcpy(buf,":0.0");

	if ((mydisplay = XOpenDisplay(buf)) == NULL) {
		debug2((errout, "failed to open display '%s'\n",buf));
		return NULL;
	}

	if (XineramaQueryExtension(mydisplay, &evb, &erb) != 0
	 && XineramaIsActive(mydisplay)) {

		xai = XineramaQueryScreens(mydisplay, &dcount);

		if (xai == NULL || dcount == 0) {
			debug("XineramaQueryScreens failed\n");
			XCloseDisplay(mydisplay);
			return NULL;
		}
		j = 0;
	} else {
		dcount = ScreenCount(mydisplay);
		j = DefaultScreen(mydisplay);
	}

	/* Allocate our list */
	if ((disps = (disppath **)calloc(sizeof(disppath *), dcount + 1)) == NULL) {
		debug("get_displays failed on malloc\n");
		XCloseDisplay(mydisplay);
		return NULL;
	}
	for (i = 0; i < dcount; i++) {
		if ((disps[i] = calloc(sizeof(buf), 1)) == NULL) {
			debug("get_displays failed on malloc\n");
			free_disppaths(disps);
			XCloseDisplay(mydisplay);
			return NULL;
		}
	}

	/* Create a description for each screen */
	for (i = 0; i < dcount; i++) {
		char desc1[100], desc2[200];
	    XF86VidModeMonitor monitor;
		int evb = 0, erb = 0;

		if ((disps[i]->name = strdup(buf)) == NULL) {
			debug("get_displays failed on malloc\n");
			free_disppaths(disps);
			XCloseDisplay(mydisplay);
			return NULL;
		}
		if (xai != NULL) {					/* Xinerama */
			disps[i]->screen = 0;
			disps[i]->uscreen = i;			/* We are assuming xinerama lists screens in the same order */
			disps[i]->rscreen = i;
			disps[i]->sx = xai[i].x_org;
			disps[i]->sy = xai[i].y_org;
			disps[i]->sw = xai[i].width;
			disps[i]->sh = xai[i].height;
		} else {
			disps[i]->screen = i;
			disps[i]->uscreen = i;
			disps[i]->rscreen = i;
			disps[i]->sx = 0;			/* Must be 0 */
			disps[i]->sy = 0;
			disps[i]->sw = DisplayWidth(mydisplay, disps[i]->screen);
			disps[i]->sh = DisplayHeight(mydisplay, disps[i]->screen);
		}

		if ((disps[i]->icc_atom_name = (char *)malloc(30)) == NULL) {
			debug("get_displays failed on malloc\n");
			free_disppaths(disps);
			XCloseDisplay(mydisplay);
			return NULL;
		}

		/* Create the X11 root atom of the default screen */
		/* that may contain the associated ICC profile */
		if (disps[i]->uscreen == 0)
			strcpy(disps[i]->icc_atom_name, "_ICC_PROFILE");
		else
			sprintf(disps[i]->icc_atom_name, "_ICC_PROFILE_%d",disps[i]->uscreen);

		if (XF86VidModeQueryExtension(mydisplay, &evb, &erb) != 0) {
			/* Some propietary multi-screen drivers (ie. TwinView & MergeFB) */
			/* don't implement the XVidMode extension properly. */
			if (XSetErrorHandler(null_error_handler) == 0) {
				debug("get_displays failed on XSetErrorHandler\n");
				if (xai != NULL)
					XFree(xai);
				free_disppaths(disps);
				XCloseDisplay(mydisplay);
				return NULL;
			}
			monitor.model = NULL;
			if (XF86VidModeGetMonitor(mydisplay, disps[i]->uscreen, &monitor) != 0
			 && monitor.model != NULL && monitor.model[0] != '\000')
				sprintf(desc1, "%s",monitor.model);
			else
				sprintf(desc1,"Screen %d",i+1);
			XSetErrorHandler(NULL);
		} else
			sprintf(desc1,"Screen %d",i+1);

		sprintf(desc2,"%s at %d, %d, width %d, height %d",desc1,
	        disps[i]->sx, disps[i]->sy, disps[i]->sw, disps[i]->sh);
		if ((disps[i]->description = strdup(desc2)) == NULL) {
			debug("get_displays failed on malloc\n");
			free_disppaths(disps);
			XCloseDisplay(mydisplay);
			return NULL;
		}
	}

	/* Put the screen given by the display name at the top */
	tdispp = disps[j];
	disps[j] = disps[0];
	disps[0] = tdispp;

	if (xai != NULL)
		XFree(xai);

	XCloseDisplay(mydisplay);

#endif /* UNIX X11 */

	return disps;
}

void free_disppaths(disppath **disps) {
	if (disps != NULL) {
		int i;
		for (i = 0; ; i++) {
			if (disps[i] == NULL)
				break;

			if (disps[i]->description != NULL)
				free(disps[i]->description);
#if defined(UNIX) && !defined(__APPLE__)
			if (disps[i]->name != NULL)
				free(disps[i]->name);
			if (disps[i]->icc_atom_name != NULL)
				free(disps[i]->icc_atom_name);
#endif	/* UNIX X11 */
			free(disps[i]);
		}
		free(disps);
	}
}
	
/* ----------------------------------------------- */
/* Deal with selecting a display */

/* Return the given display given its index 0..n-1 */
disppath *get_a_display(int ix) {
	disppath **paths, *rv = NULL;
	int i;

	if ((paths = get_displays()) == NULL)
		return NULL;

	for (i = 0; ;i++) {
		if (paths[i] == NULL) {
			free_disppaths(paths);
			return NULL;
		}
		if (i == ix)
			break;
	}
	if ((rv = malloc(sizeof(disppath))) == NULL) {
		debug("get_a_display failed malloc\n");
		free_disppaths(paths);
		return NULL;
	}
	*rv = *paths[i];		/* Structure copy */
	if ((rv->description = strdup(paths[i]->description)) == NULL) {
		debug("get_displays failed on malloc\n");
		free(rv);
		free_disppaths(paths);
		return NULL;
	}
#if defined(UNIX) && !defined(__APPLE__)
	if ((rv->name = strdup(paths[i]->name)) == NULL) {
		debug("get_displays failed on malloc\n");
		free(rv->description);
		free(rv);
		free_disppaths(paths);
		return NULL;
	}
	if ((rv->icc_atom_name = strdup(paths[i]->icc_atom_name)) == NULL) {
		debug("get_displays failed on malloc\n");
		free(rv->description);
		free(rv->name);
		free(rv);
		free_disppaths(paths);
		return NULL;
	}
#endif	/* UNXI X11 */
	free_disppaths(paths);
	return rv;
}

void free_a_disppath(disppath *path) {
	if (path != NULL) {
		if (path->description != NULL)
			free(path->description);
#if defined(UNIX) && !defined(__APPLE__)
			if (path->name != NULL)
				free(path->name);
			if (path->icc_atom_name != NULL)
				free(path->icc_atom_name);
#endif	/* UNXI X11 */
		free(path);
	}
}

/* ----------------------------------------------- */

static ramdac *dispwin_clone_ramdac(ramdac *r);
static void dispwin_setlin_ramdac(ramdac *r);
static void dispwin_del_ramdac(ramdac *r);

/* For VideoLUT/RAMDAC use, we assume that the number of entries in the RAMDAC */
/* meshes perfectly with the display raster depth, so that we can */
/* figure out how to apportion device values. We fail if they don't */
/* seem to mesh. */

/* Get RAMDAC values. ->del() when finished. */
/* Return NULL if not possible */
static ramdac *dispwin_get_ramdac(dispwin *p) {
	ramdac *r = NULL;
	int i, j;

#ifdef NT
	WORD vals[3][256];		/* 16 bit elements */

	debug("dispwin_get_ramdac called\n");

#ifdef NEVER	/* GetDeviceCaps(COLORMGMTCAPS) doesn't work for hdc from CreateDC() ! */
	if ((GetDeviceCaps(p->hdc, COLORMGMTCAPS) & CM_GAMMA_RAMP) == 0) {
		debug("dispwin_get_ramdac failed on GetDeviceCaps(CM_GAMMA_RAMP)\n");
		return NULL;
	}
#endif

	/* Allocate a ramdac */
	if ((r = (ramdac *)calloc(sizeof(ramdac), 1)) == NULL) {
		debug("dispwin_get_ramdac failed on malloc()\n");
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
			debug("dispwin_get_ramdac failed on malloc()\n");
			return NULL;
		}
	}

	/* GetDeviceGammaRamp() is hard coded for 3 x 256 entries */
	if (r->nent != 256) {
		debug("GetDeviceGammaRamp() is hard coded for nent == 256\n");
		return NULL;
	}

	if (GetDeviceGammaRamp(p->hdc, vals) == 0) {
		debug("dispwin_get_ramdac failed on GetDeviceGammaRamp()\n");
		return NULL;
	}
	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			r->v[j][i] = vals[j][i]/65535.0;
		}
	}
#endif	/* NT */

#ifdef __APPLE__
	unsigned int nent;
	CGGammaValue vals[3][16385];

	debug("dispwin_get_ramdac called\n");

	/* Could try and use the ColorSync CMSetGammaByAVID() and */
	/* CMGetGammaByAVID() routines instead - they may be persistent */
	/* after the application exits ? */
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
	if ((r = (ramdac *)calloc(sizeof(ramdac), 1)) == NULL) {
		debug("dispwin_get_ramdac failed on malloc()\n");
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
			debug("dispwin_get_ramdac failed on malloc()\n");
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
	/* Some propietary multi-screen drivers (ie. TwinView & MergedFB) */
	/* don't implement the XVidMode extenstion properly. */
	if (XSetErrorHandler(null_error_handler) == 0) {
		debug("get_displays failed on XSetErrorHandler\n");
		return NULL;
	}
	nent = -1;
	if (XF86VidModeGetGammaRampSize(p->mydisplay, p->myrscreen, &nent) == 0
	 || nent == -1) {
		XSetErrorHandler(NULL);
		debug("XF86VidModeGetGammaRampSize failed\n");
		return NULL;
	}
	XSetErrorHandler(NULL);		/* Restore handler */
	if (nent == 0) {
		debug("XF86VidModeGetGammaRampSize returned 0 size\n");
		return NULL;
	}

	if (nent > 16384) {
		debug("XF86VidModeGetGammaRampSize has more entries than we can handle\n");
		return NULL;
	}

	if (XF86VidModeGetGammaRamp(p->mydisplay, p->myrscreen, nent,  vals[0], vals[1], vals[2]) == 0) {
		debug("XF86VidModeGetGammaRamp failed\n");
		return NULL;
	}

	if (nent != (1 << p->pdepth)) {
		debug("CGGetDisplayTransferByTable number of entries mismatches screen depth\n");
		return NULL;
	}

	/* Allocate a ramdac */
	if ((r = (ramdac *)calloc(sizeof(ramdac), 1)) == NULL) {
		debug("dispwin_get_ramdac failed on malloc()\n");
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
			debug("dispwin_get_ramdac failed on malloc()\n");
			return NULL;
		}
	}

	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			r->v[j][i] = vals[j][i]/65535.0;
		}
	}
#endif	/* UNXI X11 */
	return r;
}

/* Set the RAMDAC values. */
/* Return nz if not possible */
static int dispwin_set_ramdac(dispwin *p, ramdac *r) {
	int i, j;

#ifdef NT
	WORD vals[3][256];		/* 16 bit elements */

	debug("dispwin_set_ramdac called\n");

#ifdef NEVER	/* GetDeviceCaps(COLORMGMTCAPS) doesn't work for hdc from CreateDC() ! */
	if ((GetDeviceCaps(p->hdc, COLORMGMTCAPS) & CM_GAMMA_RAMP) == 0) {
		debug("dispwin_set_ramdac failed on GetDeviceCaps(CM_GAMMA_RAMP)\n");
		return 1;
	}
#endif

	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			double vv = r->v[j][i];
			if (vv < 0.0)
				vv = 0.0;
			else if (vv > 1.0)
				vv = 1.0;
			vals[j][i] = (int)(65535.0 * vv + 0.5);
		}
	}

	if (SetDeviceGammaRamp(p->hdc, vals) == 0) {
		debug("dispwin_set_ramdac failed on SetDeviceGammaRamp()\n");
		return 1;
	}
#endif	/* NT */

#ifdef __APPLE__
	CGGammaValue vals[3][16384];

	debug("dispwin_set_ramdac called\n");

	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			double vv = r->v[j][i];
			if (vv < 0.0)
				vv = 0.0;
			else if (vv > 1.0)
				vv = 1.0;
			vals[j][i] = vv;
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
			double vv = r->v[j][i];
			if (vv < 0.0)
				vv = 0.0;
			else if (vv > 1.0)
				vv = 1.0;
			vals[j][i] = (int)(vv * 65535.0 + 0.5);
		}
	}
	/* Some propietary multi-screen drivers (ie. TwinView & MergedFB) */
	/* don't implement the XVidMode extenstion properly. */
	if (XSetErrorHandler(null_error_handler) == 0) {
		debug("get_displays failed on XSetErrorHandler\n");
		return 1;
	}
	if (XF86VidModeSetGammaRamp(p->mydisplay, p->myrscreen, r->nent, vals[0], vals[1], vals[2]) == 0) {
		XSetErrorHandler(NULL);
		debug("XF86VidModeSetGammaRamp failed\n");
		return 1;
	}
	XSetErrorHandler(NULL);
#endif	/* UNXI X11 */

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

	debug("dispwin_set_color called\n");

	if (p->nowin)
		return 1;

	p->rgb[0] = r;
	p->rgb[1] = g;
	p->rgb[2] = b;

	for (j = 0; j < 3; j++) {
		if (p->rgb[j] < 0.0)
			p->rgb[j] = 0.0;
		else if (p->rgb[j] > 1.0)
			p->rgb[j] = 1.0;
		p->r_rgb[j] = p->rgb[j];
	}

	if (p->donat) {		/* Output high precision native using RAMDAC */
		double prange = p->r->nent - 1.0;

		for (j = 0; j < 3; j++) {
			int tt;

//printf("~1 %d: in %f, ",j,p->rgb[j]);
			tt = (int)(p->rgb[j] * prange);
			p->r->v[j][tt] = p->rgb[j];
			p->r_rgb[j] = (double)tt/prange;	/* Quantized value */
//printf(" cell[%d], val %f, rast val %f\n",tt, p->rgb[j], p->r_rgb[j]);
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
		HRGN regn;
		HBRUSH hbr;
		MSG msg;
		int vali[3];

		/* Stop the system going to sleep */
		SetThreadExecutionState(ES_DISPLAY_REQUIRED);

		if ((regn = CreateRectRgn(p->tx, p->ty, p->tx + p->tw, p->ty + p->th)) == NULL) {
			debug2((errout,"CreateRectRgn failed, lasterr = %d\n",GetLastError()));
			return 1;
		}

		/* Force a repaint with the new data */
		if (!InvalidateRgn(p->hwnd,regn,TRUE)) {
			debug2((errout,"InvalidateRgn failed, lasterr = %d\n",GetLastError()));
			return 1;
		}

		/* Convert to 8 bit color */
		for (j = 0; j < 3; j++) {
			vali[j] = (int)(255.0 * p->r_rgb[j] + 0.5);
		}

		hdc = BeginPaint(p->hwnd, &ps);
		SaveDC(hdc);

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

		SetRect(&rect, p->tx, p->ty, p->tx + p->tw, p->ty + p->th);
		FillRect(hdc, &rect, hbr);

		RestoreDC(hdc,-1);
		EndPaint(p->hwnd, &ps);
		
		UpdateWindow(p->hwnd);		/* Flush the paint */
 
		/* Process any pending messages */
		while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
//printf("~1 got message %d\n",msg.message);
			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}
		DeleteDC(hdc);
		DeleteObject(regn);
	}
#endif /* NT */

	/* - - - - - - - - - - - - - - */

#ifdef __APPLE__

	/* Stop the system going to sleep */
    UpdateSystemActivity(OverallAct);

	/* Cause window repaint with the new data */
	{
		OSStatus stat;
		Rect wRect;
	    SetRect(&wRect,p->tx,p->ty,p->tw,p->th); /* left, top, right, bottom */
		if ((stat = InvalWindowRect(p->mywindow, &wRect)) != noErr) {
			debug2((errout,"InvalWindowRect failed with %d\n",stat));
			return 1;
		}
	}

	/* Make sure our window is brought to the front at least once, */
	/* but not every time, in case the user wants to kill the application. */
	if (p->btf == 0){
		OSStatus stat;
		ProcessSerialNumber cpsn;
		if ((stat = GetCurrentProcess(&cpsn)) != noErr) {
			debug2((errout,"GetCurrentProcess returned error %d\n",stat));
		} else {
#if MAC_OS_X_VERSION_MAX_ALLOWED >= 1030
			if ((stat = TransformProcessType(&cpsn, kProcessTransformToForegroundApplication)) != noErr) {
				debug2((errout,"TransformProcessType returned error %d\n",stat));
			}
#endif /* OS X 10.3 */
			if ((stat = SetFrontProcess(&cpsn)) != noErr) {
				debug2((errout,"SetFrontProcess returned error %d\n",stat));
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
        static time_t ltime = 0;
        time_t ttime;
		Colormap mycmap;
		XColor col;
		int vali[3];

        /* Because this is really slow, we do it every 60 seconds */
        if (p->sssuspend == 0 && (ttime = time(NULL)) - ltime > 60) {
			if (fork() == 0) {
				if (fork() == 0) {
					/* Try and stop xscreensaver messing things up */
					/* It's a pitty xscreensaver isn't more cooperative with X11... */
					freopen("/dev/null", "r", stdin);
					freopen("/dev/null", "a", stdout);		/* Hide output */
					freopen("/dev/null", "a", stderr);
					execlp("xscreensaver-command", "xscreensaver-command", "-deactivate", NULL);
					/* exec won't normally return */
					_exit(0);

				} else {
					sleep(5);		/* Stagger this a little */

					/* Try and stop GnomeScreensaver messing things up */
					/* It's a pitty GnomeScreensaver isn't more cooperative with X11 too... */
					if (fork() == 0) {
						freopen("/dev/null", "r", stdin);
						freopen("/dev/null", "a", stdout);		/* Hide output */
						freopen("/dev/null", "a", stderr);
						execlp("dbus-send", "dbus-send", "--dest=org.gnome.ScreenSaver", "--type=method_call", "/org/gnome/ScreenSaver", "org.gnome.ScreenSaver.SimulateUserActivity", NULL);
						/* exec won't normally return */
						_exit(0);
					}
				}
				_exit(0);
			}
            ltime = ttime;
		}

		/* Convert to 16 bit color */
		for (j = 0; j < 3; j++)
			vali[j] = (int)(65535.0 * p->r_rgb[j] + 0.5);

		mycmap = DefaultColormap(p->mydisplay, p->myscreen);
		col.red = vali[0];
		col.green = vali[1];
		col.blue = vali[2];
		XAllocColor(p->mydisplay, mycmap, &col);
		XSetForeground(p->mydisplay, p->mygc, col.pixel);

		XFillRectangle(p->mydisplay, p->mywindow, p->mygc,
		               p->tx, p->ty, p->tw, p->th);

		XSync(p->mydisplay, False);		/* Make sure it happens */
	}
#endif	/* UNXI X11 */

	if (p->callout != NULL) {
		int rv;
		char *cmd;

		if ((cmd = malloc(strlen(p->callout) + 200)) == NULL)
			error("Malloc of command string failed");

		sprintf(cmd, "%s %d %d %d %f %f %f",p->callout,
			        (int)(r * 255.0 + 0.5),(int)(g * 255.0 + 0.5),(int)(b * 255.0 + 0.5), r, g, b);
		if ((rv = system(cmd)) != 0)
			warning("System command '%s' failed with %d",cmd,rv); 
		free(cmd);
	}

	/* Allow some time for the display to update before */
	/* a measurement can take place. This allows for CRT */
	/* refresh, or LCD processing/update time. */
	msec_sleep(60);

	return 0;
}

/* ----------------------------------------------- */
/* Set the shell set color callout */
void dispwin_set_callout(
dispwin *p,
char *callout
) {
	debug2((errout,"dispwin_set_callout called with '%s'\n",callout));

	p->callout = strdup(callout);
}

/* ----------------------------------------------- */
/* Destroy ourselves */
static void dispwin_del(
dispwin *p
) {

	debug("dispwin_del called\n");

	if (p == NULL)
		return;

	if (p->callout != NULL)
		free(p->callout);

	/* Restore original RAMDAC if we were in native mode */
	if (p->donat && p->or != NULL) {
		p->set_ramdac(p, p->or);
		p->or->del(p->or);
		p->r->del(p->r);
		debug("Restored original ramdac\n");
	}

	/* -------------------------------------------------- */
#ifdef NT
	if (p->hwnd != NULL) {
		if (!DestroyWindow(p->hwnd)) {
			debug2((errout, "DestroyWindow failed, lasterr = %d\n",GetLastError()));
		}	
//		DestroyCursor(p->curs); 
	}

	UnregisterClass(p->AppName, NULL);

	if (p->hdc != NULL)
		DeleteDC(p->hdc);

#endif /* NT */
	/* -------------------------------------------------- */

	/* -------------------------------------------------- */
#ifdef __APPLE__

	if (p->nowin == 0) {	/* We have a window up */
		QDEndCGContext(p->port, &p->mygc);		/* Ignore any errors */
		DisposeWindow(p->mywindow);	
	}

	CGDisplayShowCursor(p->ddid);
//	CGDisplayRelease(p->ddid);

#endif /* __APPLE__ */
	/* -------------------------------------------------- */

	/* -------------------------------------------------- */
#if defined(UNIX) && !defined(__APPLE__)
	debug("About to close display\n");

	if (p->nowin == 0) {	/* We have a window up */

		/* Restore the screensaver state */
#if ScreenSaverMajorVersion >= 1 && ScreenSaverMinorVersion >= 1	/* X11R7.1 */
		if (p->sssuspend) {
			XScreenSaverSuspend(p->mydisplay, False);
			p->sssuspend = 0;
		}
#endif
		if (p->ssvalid) {
			XSetScreenSaver(p->mydisplay, p->timeout, p->interval,
		                p->prefer_blanking, p->allow_exposures);
		}
	
		XFreeGC(p->mydisplay, p->mygc);
		XDestroyWindow(p->mydisplay, p->mywindow);
	}

	XCloseDisplay(p->mydisplay);
	debug("finished\n");
#endif	/* UNXI X11 */
	/* -------------------------------------------------- */

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
					Rect wRect;
					debug("Event: Bounds Changed\n");
					GetPortBounds(p->port, &wRect);
					if ((stat = InvalWindowRect(p->mywindow, &wRect)) != noErr) {
						debug2((errout,"InvalWindowRect failed with %d\n",stat));
					}
					break;
				}
				case kEventWindowDrawContent: {
					CGRect frect;

					debug("Event: Draw Content\n");
					/* If we're using an overlay window, paint it all black */
					if (p->blackbg && p->firstdraw == 0) {
						frect = CGRectMake(0.0, 0.0, (float)p->ww, (float)p->wh);
						CGContextSetRGBFillColor(p->mygc, 0.0, 0.0, 0.0, 1.0);
						CGContextFillRect(p->mygc, frect);
						p->firstdraw = 1;
					}
					frect = CGRectMake((float)p->tx, (float)(p->wh - p->ty - p->th - 1),
					  (float)(1.0 + p->tw), (float)(1.0 + p->th ));
					CGContextSetRGBFillColor(p->mygc, p->r_rgb[0], p->r_rgb[1], p->r_rgb[2], 1.0);
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
#endif	/* UNXI X11 */

/* ----------------------------------------------- */
/* Create a RAMDAC access and display test window, default white */
dispwin *new_dispwin(
disppath *disp,					/* Display to calibrate. */
double width, double height,	/* Width and height in mm */
double hoff, double voff,		/* Offset from center in fraction of screen, range -1.0 .. 1.0 */
int nowin,						/* NZ if no window should be created - RAMDAC access only */
int native,						/* NZ if ramdac should be bypassed rather than used. */
int blackbg,					/* NZ if whole screen should be filled with black */
int override					/* NZ if override_redirect is to be used on X11 */
) {
	dispwin *p = NULL;

	debug("new_dispwin called\n");

	if ((p = (dispwin *)calloc(sizeof(dispwin), 1)) == NULL)
		return NULL;

	p->nowin = nowin;
	p->donat = native;
	p->blackbg = blackbg;
	p->get_ramdac   = dispwin_get_ramdac;
	p->set_ramdac   = dispwin_set_ramdac;
	p->set_color    = dispwin_set_color;
	p->set_callout  = dispwin_set_callout;
	p->del          = dispwin_del;

	/* Basic object is initialised, so create a window */

	/* -------------------------------------------------- */
#ifdef NT
	{
		MSG msg;
		WNDCLASS wc;
		int disp_hsz, disp_vsz;		/* Display horizontal/vertical size in mm */
		int disp_hrz, disp_vrz;		/* Display horizontal/vertical resolution in pixels */
		int wi, he;					/* Width and height of window in pixels */
		int xo, yo;					/* Window location in pixels */
		int bpp;
		
		p->AppName = "Argyll Test Window";
		
		debug2((errout, "About to open display '%s'\n",disp->name));

		/* Get device context to main display */
		if ((p->hdc = CreateDC(disp->name, disp->name, NULL, NULL)) == NULL) {
			debug2((errout, "CreateDC failed, lasterr = %d\n",GetLastError()));
			dispwin_del(p);
			return NULL;
		}

		disp_hsz = GetDeviceCaps(p->hdc, HORZSIZE);	/* mm */
		disp_vsz = GetDeviceCaps(p->hdc, VERTSIZE);
		disp_hrz = GetDeviceCaps(p->hdc, HORZRES);	/* pixels */
		disp_vrz = GetDeviceCaps(p->hdc, VERTRES);

		wi = (int)(width * (double)disp_hrz/(double)disp_hsz + 0.5);
		if (wi > disp_hrz)
			wi = disp_hrz;
		he = (int)(height * (double)disp_vrz/(double)disp_vsz + 0.5);
		if (he > disp_vrz)
			he = disp_vrz;

		if (p->blackbg) {	/* Window fills the screen, test area is within it */
			p->tx = (int)((hoff * 0.5 + 0.5) * (disp->sw - wi) + 0.5);
			p->ty = (int)((voff * 0.5 + 0.5) * (disp->sh - he) + 0.5);
			p->tw = wi;
			p->th = he;
			wi = disp->sw;
			he = disp->sh;
			xo = disp->sx;
			yo = disp->sy;
		} else {			/* Test area completely fills the window */
			p->tx = 0;
			p->ty = 0;
			p->tw = wi;
			p->th = he;
			xo = disp->sx + (int)((hoff * 0.5 + 0.5) * (disp->sw - wi) + 0.5);
			yo = disp->sy + (int)((voff * 0.5 + 0.5) * (disp->sh - he) + 0.5);
		}
		p->ww = wi;
		p->wh = he;

		/* It's a bit difficult to know how windows defines the display */
		/* depth. Microsofts doco is fuzzy, and typical values */
		/* for BITSPIXEL and PLANES are confusing (What does "32" and "1" */
		/* mean ?) The doco for COLORRES is also fuzzy, but it returns */
		/* a meaningful number on my box (24) */
		bpp = GetDeviceCaps(p->hdc, COLORRES);
		if (bpp <= 0)
			p->pdepth = 8;		/* Assume this is so */
		else
			p->pdepth = bpp/3;

		if (nowin == 0) {

			/* Fill in window class structure with parameters that describe the */
			/* main window. */
			wc.style         = CS_SAVEBITS ;	/* Class style(s). */
			wc.lpfnWndProc   = DefWindowProc; /* Function to retrieve messages for class windows */
			wc.cbClsExtra    = 0;			/* No per-class extra data. */
			wc.cbWndExtra    = 0;			/* No per-window extra data. */
			wc.hInstance     = NULL;		/* Application that owns the class. */
			wc.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
			wc.hCursor       = LoadCursor(NULL, IDC_CROSS);
			wc.hbrBackground = GetStockObject(BLACK_BRUSH);
			wc.lpszMenuName  = NULL;
			wc.lpszClassName = p->AppName;

			/* Make the cursor disapear over our window */
			/* (How does it know our window ??) */
			ShowCursor(FALSE);

			if ((p->arv = RegisterClass(&wc)) == 0) {
				debug2((errout, "RegisterClass failed, lasterr = %d\n",GetLastError()));
				dispwin_del(p);
				return NULL;
			}

			p->rgb[0] = p->rgb[1] = p->rgb[2] = 1.0;	/* Set White */

			p->hwnd = CreateWindowEx(
				WS_EX_TOPMOST,
				p->AppName,
				"Argyll Display Calibration Window",
				WS_DISABLED | WS_POPUP, 
				xo, yo,					/* Location */
				wi, he,					/* Size */
				NULL,					/* Handle to parent or owner */
				NULL,					/* Handle to menu or child window */
				NULL, 					/* hInstance Handle to appication instance */
				NULL);					/* pointer to window creation data */

			if (!p->hwnd) {
				debug2((errout, "CreateWindow failed, lasterr = %d\n",GetLastError()));
				dispwin_del(p);
				return NULL;
			}
			
			/* Make cursor dissapear in our window. */
			/* (Probably have to use something more complicated for */
			/*  a true GUI application - event loop & window class cursor ?) */
			ShowWindow(p->hwnd, SW_SHOWNA);

			/*
				Should we call BOOL SystemParametersInfo()
				to disable high contrast, powertimout and screensaver timeout ?

			 */

			/* Process any pending messages */
			while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
				TranslateMessage(&msg);
				DispatchMessage(&msg);
			}
		}
	}
#endif /* NT */
	/* -------------------------------------------------- */

	/* -------------------------------------------------- */
#ifdef __APPLE__

	p->ddid = disp->ddid;		/* Display we're working on */
	p->pdepth = CGDisplayBitsPerSample(p->ddid);

	if (nowin == 0) {			/* Create a window */
		OSStatus stat;
		EventHandlerRef fHandler;		/* Event handler reference */
		const EventTypeSpec	kEvents[] =
		{ 
			{ kEventClassWindow, kEventWindowClose },
			{ kEventClassWindow, kEventWindowBoundsChanged },
			{ kEventClassWindow, kEventWindowDrawContent }
		};

		CGSize sz;				/* Display size in mm */
		int wi, he;				/* Width and height in pixels */
		int xo, yo;				/* Window location */
	    Rect	wRect;
		WindowClass wclass = 0;
		WindowAttributes attr = kWindowNoAttributes;

		/* Choose the windows class and attributes */
//		wclass = kAlertWindowClass;				/* Above everything else */
//		wclass = kModalWindowClass;
//		wclass = kFloatingWindowClass;
//		wclass = kUtilityWindowClass; /* The usefule one */
//		wclass = kHelpWindowClass;
		wclass = kOverlayWindowClass;
//		wclass = kSimpleWindowClass;

		attr |= kWindowDoesNotCycleAttribute;
		attr |= kWindowNoActivatesAttribute;
//		attr |= kWindowStandardFloatingAttributes;

		attr |= kWindowStandardHandlerAttribute; /* This window has the standard Carbon Handler */
		attr |= kWindowNoShadowAttribute;		/* usual */
//		attr |= kWindowNoTitleBarAttribute;		/* usual - but not with Overlay Class */
		attr |= kWindowIgnoreClicksAttribute;		/* usual */

		sz = CGDisplayScreenSize(p->ddid);
		debug2((errout," Display size = %f x %f mm\n",sz.width,sz.height));

		wi = (int)(width * disp->sw/sz.width + 0.5);
		if (wi > disp->sw)
			wi = disp->sw;
		he = (int)(height * disp->sh/sz.height + 0.5);
		if (he > disp->sh)
			he = disp->sh;

		if (p->blackbg) {	/* Window fills the screen, test area is within it */
			p->tx = (int)((hoff * 0.5 + 0.5) * (disp->sw - wi) + 0.5);
			p->ty = (int)((voff * 0.5 + 0.5) * (disp->sh - he) + 0.5);
			p->tw = wi;
			p->th = he;
			wi = disp->sw;
			he = disp->sh;
			xo = disp->sx;
			yo = disp->sy;
		} else {			/* Test area completely fills the window */
			p->tx = 0;
			p->ty = 0;
			p->tw = wi;
			p->th = he;
			xo = disp->sx + (int)((hoff * 0.5 + 0.5) * (disp->sw - wi) + 0.5);
			yo = disp->sy + (int)((voff * 0.5 + 0.5) * (disp->sh - he) + 0.5);
		}
		p->ww = wi;
		p->wh = he;

//printf("~1 Got size %d x %d at %d %d from %fmm x %fmm\n",wi,he,xo,yo,width,height);
	    SetRect(&wRect,xo,yo,xo+wi,yo+he); /* left, top, right, bottom */

		/* Create invisible new window of given class, attributes and size */
	    stat = CreateNewWindow(wclass, attr, &wRect, &p->mywindow);
		if (stat != noErr || p->mywindow == nil) {
			debug2((errout,"CreateNewWindow failed with code %d\n",stat));
			dispwin_del(p);
			return NULL;
		}

		/* Set a title */
		SetWindowTitleWithCFString(p->mywindow, CFSTR("Argyll Window"));

		/* Set the window blackground color to black */
		{
			RGBColor col = { 0,0,0 };
			SetWindowContentColor(p->mywindow, &col);
		}

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
			debug2((errout,"InstallEventHandler failed with code %d\n",stat));
			dispwin_del(p);
			return NULL;
		}

		/* Create a GC */
		p->port = GetWindowPort(p->mywindow);

		if ((stat = QDBeginCGContext(p->port, &p->mygc)) != noErr) {
			debug2((errout,"QDBeginCGContext returned error %d\n",stat));
			dispwin_del(p);
			return NULL;
		}
		p->winclose = 0;
		p->rgb[0] = p->rgb[1] = p->rgb[2] = 1.0;	/* Set White */

		/* Activate the window */

#ifdef NEVER
		/* Make window pop to the top */
		{
			ProcessSerialNumber cpsn;
			if ((stat = GetCurrentProcess(&cpsn)) != noErr) {
				debug2((errout,"GetCurrentProcess returned error %d\n",stat));
			} else {
				if ((stat = SetFrontProcess(&cpsn)) != noErr) {
					debug2((errout,"SetFrontProcess returned error %d\n",stat));
				}
			}
		}
#endif /* NEVER */
//		BringToFront(p->mywindow);
		ShowWindow(p->mywindow);	/* Makes visible and triggers update event */
//		ActivateWindow(p->mywindow, false);
		SendBehind(p->mywindow, NULL);
//		SelectWindow(p->mywindow);

//printf("~1 created window, about to enter event loop\n");
		/* Run the event loop ?? */
		/* Seems to hang because draw happens before event loop gets run ???/ */
		/* RunApplicationEventLoop(); */

		CGDisplayHideCursor(p->ddid);
//		CGDisplayCapture(p->ddid);

#ifdef NEVER
		/* Carbon doesn't support setting a custom cursor yet :-( */
		/* The QuickDraw way of making the cursor vanish. */
		/* unfortuneately, it vanished for the whole screen, */
		/* and comes back when you click the mouse, just like. */
		/* CGDisplayHideCursor() */
		{
			InitCursor();
			Cursor myc = {
				{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},		/* data */
				{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},		/* Mask */
				{ 0, 0 }								/* hotSpot */
			};
			
			SetCursor(&myc);
		}
#endif
/*
		To disable screensaver,
		could try exec("defaults ....."); ??
		see "man defaults"
		Maybe the universal access, "Enhance Contrast" can
		be zero'd this way too ?? 

		To disable Universal Access contrast, could try
		exec("com.apple.universalaccess contrast 0");
		but this doesn't activate the change.

		Applescript may be able to do this in an ugly way,
		something like:

			tell application "System Preferences"
				activate
			end tell

			tell application "System Events"
				get properties
				tell process "System Preferences"
					click menu item "Universal Access" of menu "View" of menu bar 1
					tell window "Energy Saver"
						tell group 1
							tell tab group 1
								set the value of slider 2 to 600
							end tell
						end tell
					end tell
				end tell
			end tell

 */

//printf("~1 done create window\n");
		if (p->winclose) {
			dispwin_del(p);
			return NULL;
		}
	}
#endif /* __APPLE__ */
	/* -------------------------------------------------- */

	/* -------------------------------------------------- */
#if defined(UNIX) && !defined(__APPLE__)
	{
		/* NOTE: That we're not doing much to detect if the display/window
		   we open is unsuitable for high quality color (ie. at least
		   24 bit etc.
		 */

		/* stuff for X windows */
		Window rootwindow;
		char *appname = "TestWin";
		Visual *myvisual;
		XSetWindowAttributes myattr;
		XEvent      myevent;
		XTextProperty myappname;
		XSizeHints  mysizehints;
		XWMHints    mywmhints;
		int evb = 0, erb = 0;		/* Extension version */
		unsigned long myforeground,mybackground;
		int disp_hsz, disp_vsz;		/* Display horizontal/vertical size in mm */
		int disp_hrz, disp_vrz;		/* Display horizontal/vertical resolution in pixels (virtual screen) */
		int wi, he;				/* Width and height of window in pixels */
		int xo, yo;				/* Window location in pixels */
	
		/* open the display */
		p->mydisplay = XOpenDisplay(disp->name);
		if(!p->mydisplay) {
			debug2((errout,"Unable to open display '%s'\n",disp->name));
			dispwin_del(p);
			return NULL;
		}
		debug("Opened display OK\n");

		p->myscreen = disp->screen;
		p->myuscreen = disp->uscreen;
		p->myrscreen = disp->rscreen;

		//p->pdepth = DefaultDepth(p->mydisplay, p->myscreen)/3;
		myvisual = DefaultVisual(p->mydisplay, p->myscreen);
		p->pdepth = myvisual->bits_per_rgb;

		if (nowin == 0) {			/* Create a window */
			rootwindow = RootWindow(p->mydisplay, p->myscreen);

			myforeground = WhitePixel(p->mydisplay, p->myscreen);
			mybackground = BlackPixel(p->mydisplay, p->myscreen);
		
			/* Get device context to main display */
			disp_hsz = DisplayWidthMM(p->mydisplay, p->myscreen);
			disp_vsz = DisplayHeightMM(p->mydisplay, p->myscreen);
			disp_hrz = DisplayWidth(p->mydisplay, p->myscreen);
			disp_vrz = DisplayHeight(p->mydisplay, p->myscreen);

			/* Compute width and offset from overal display in case Xinerama is active */
			wi = (int)(width * (double)disp_hrz/(double)disp_hsz + 0.5);
			if (wi > disp_hrz)
				wi = disp_hrz;
			he = (int)(height * (double)disp_vrz/(double)disp_vsz + 0.5);
			if (he > disp_vrz)
				he = disp_vrz;

			if (p->blackbg) {	/* Window fills the screen, test area is within it */
				p->tx = (int)((hoff * 0.5 + 0.5) * (disp->sw - wi) + 0.5);
				p->ty = (int)((voff * 0.5 + 0.5) * (disp->sh - he) + 0.5);
				p->tw = wi;
				p->th = he;
				wi = disp->sw;
				he = disp->sh;
				xo = disp->sx;
				yo = disp->sy;
			} else {			/* Test area completely fills the window */
				p->tx = 0;
				p->ty = 0;
				p->tw = wi;
				p->th = he;
				xo = disp->sx + (int)((hoff * 0.5 + 0.5) * (disp->sw - wi) + 0.5);
				yo = disp->sy + (int)((voff * 0.5 + 0.5) * (disp->sh - he) + 0.5);
			}
			p->ww = wi;
			p->wh = he;

			/* Setup Size Hints */
			mysizehints.flags = PPosition | USSize;
			mysizehints.x = xo;
			mysizehints.y = yo;
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
				0,							/* Border width */
				CopyFromParent,				/* Depth */
				InputOutput,				/* Class */
				CopyFromParent,				/* Visual */
				CWBackPixel | CWBitGravity	/* Attributes Valumask */
				| CWWinGravity | CWBackingStore | CWOverrideRedirect,
				&myattr						/* Attribute details */
			);

#ifdef NEVER
			XWindowAttributes mywattributes;

			/* Get the windows attributes */
			if (XGetWindowAttributes(
				p->mydisplay, p->mywindow,
				&mywattributes) == 0) {
				debug("XGetWindowAttributes failed\n");
				dispwin_del(p);
				return NULL;
			}
			p->pdepth = mywattributes.depth/3;
#endif

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
			{
				unsigned int xid = (unsigned int)rootwindow;	/* Hope this is 32 bit */
				XChangeProperty(
					p->mydisplay, p->mywindow,
					XA_WM_TRANSIENT_FOR,		/* Property */
					XA_WINDOW,					/* Type */
					32,							/* format = bits in type of unsigned int */
					PropModeReplace,			/* Change mode */
					(char *)(&xid),				/* Data is Root Window XID */
					1							/* Number of elements of data */
				);
			}

			p->mygc = XCreateGC(p->mydisplay,p->mywindow,0,0);
			XSetBackground(p->mydisplay,p->mygc,mybackground);
			XSetForeground(p->mydisplay,p->mygc,myforeground);
			
			/* Create an invisible cursor over our window */
			{
				Cursor mycursor;
				Pixmap mypixmap;
				Colormap mycmap;
				XColor col;
				char pmdata[1] = { 0 };

				col.red = col.green = col.blue = 0;

				mycmap = DefaultColormap(p->mydisplay, p->myscreen);
				XAllocColor(p->mydisplay, mycmap, &col);
				mypixmap = XCreatePixmapFromBitmapData(p->mydisplay, p->mywindow, pmdata, 1, 1, 0, 0, 1); 
				mycursor = XCreatePixmapCursor(p->mydisplay, mypixmap, mypixmap, &col, &col, 0,0);
				XDefineCursor(p->mydisplay, p->mywindow, mycursor);
			}

			XSelectInput(p->mydisplay,p->mywindow, ExposureMask);
		
			XMapRaised(p->mydisplay,p->mywindow);
			debug("Raised window\n");
		
			/* Suspend any screensavers if we can */

#if ScreenSaverMajorVersion >= 1 && ScreenSaverMinorVersion >= 1	/* X11R7.1 ??? */

			if (XScreenSaverQueryExtension (p->mydisplay, &evb, &erb) != 0) {
				int majv, minv;
				XScreenSaverSuspend(p->mydisplay, True);
					p->sssuspend = 1;

				/* Else we'd have to register as a screensaver to */
				/* prevent another one activating ?? */
			}
#endif	/* X11R7.1 screensaver extension */

			/* ~~~ Look into running xdg-screensaver suspend, resume ~~~ */
			if (p->sssuspend == 0) {
				/* Save the screensaver state, and then disable it */
				XGetScreenSaver(p->mydisplay, &p->timeout, &p->interval,
				                &p->prefer_blanking, &p->allow_exposures);
				XSetScreenSaver(p->mydisplay, 0, 0, DefaultBlanking, DefaultExposures);
				p->ssvalid = 1;
			}
		
			/* Deal with any pending events */
			debug("About to enter main loop\n");
			while(XPending(p->mydisplay) > 0) {
				XNextEvent(p->mydisplay, &myevent);
				switch(myevent.type) {
					case Expose:
						if(myevent.xexpose.count == 0) {	/* Repare the exposed region */
							debug("Servicing final expose\n");
							XFillRectangle(p->mydisplay, p->mywindow, p->mygc,
							               p->tx, p->ty, p->tw, p->th);
							debug("Finished expose\n");
						}
						break;
				}
			}
		}
	}
#endif	/* UNXI X11 */
	/* -------------------------------------------------- */

	/* Setup for native mode */
	if (p->donat) {
		debug("About to setup native mode\n");
		if ((p->or = p->get_ramdac(p)) == NULL
		 || (p->r = p->or->clone(p->or)) == NULL) {
			warning("Native mode can't work, no VideoLUT support");
			dispwin_del(p);
			return NULL;
		}
		p->r->setlin(p->r);
		debug("Saved original VideoLUT\n");
	}

	debug("About to exit new_dispwin()\n");
	return p;
}

/* ---------------------------------------------------------------- */
#ifdef STANDALONE_TEST
/* test code */

#include "numlib.h"

static void
usage(void) {
	disppath **dp;
	fprintf(stderr,"Test display patch window & test display LUT access, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	fprintf(stderr,"usage: dispwin [options] [calfile] \n");
	fprintf(stderr," -v                   Verbose mode\n");
#if defined(UNIX) && !defined(__APPLE__)
	fprintf(stderr," -display displayname Choose X11 display name\n");
	fprintf(stderr," -d n[,m]             Choose the display n from the following list (default 1)\n");
	fprintf(stderr,"                      Optionally choose different display m for Video LUT access\n"); 
#else
	fprintf(stderr," -d n                 Choose the display from the following list (default 1)\n");
#endif
	dp = get_displays();
	if (dp == NULL || dp[0] == NULL) {
		fprintf(stderr,"    ** No displays found **\n");
	} else {
		int i;
		for (i = 0; ; i++) {
			if (dp[i] == NULL)
				break;
			fprintf(stderr,"    %d = '%s'\n",i+1,dp[i]->description);
		}
	}
	free_disppaths(dp);
	fprintf(stderr," -p ho,vo,ss          Position test window and scale it\n");
	fprintf(stderr," -B                   Fill whole screen with black background\n");
	fprintf(stderr," -i                   Run forever with random values\n");
	fprintf(stderr," -m                   Manually cycle through initial values\n");
	fprintf(stderr," -f                   Test grey ramp fade\n");
	fprintf(stderr," -r                   Test just Video LUT loading\n");
	fprintf(stderr," -n                   Test native output (rather than through Video LUT)\n");
	fprintf(stderr," -c                   Load a linear display calibration\n");
	fprintf(stderr," -x                   Don't exit after loading a display calibration\n");
	fprintf(stderr," -V                   Verify that calfile is currently loaded\n");
#if defined(UNIX) && !defined(__APPLE__)
	fprintf(stderr," -S                   Set X11 ICC_PROFILE property to profile\n");
	fprintf(stderr," -L                   Load X11 ICC_PROFILE property profile into LUT\n");
#endif	/* UNXI X11 */
	fprintf(stderr," calfile              Load display calibration (.cal or .icm) into LUT, and exit\n");
	exit(1);
}

/* 32 bit pseudo random sequencer based on XOR feedback */
/* generates number between 1 and 4294967295 */
#define PSRAND32(S) (((S) & 0x80000000) ? (((S) << 1) ^ 0xa398655d) : ((S) << 1))

int 
main(int argc, char *argv[]) {
	int fa, nfa, mfa;			/* current argument we're looking at */
	int verb = 0;				/* Verbose flag */
	disppath *disp = NULL;		/* Display being used */
	double patscale = 1.0;		/* scale factor for test patch size */
	double ho = 0.0, vo = 0.0;	/* Test window offsets, -1.0 to 1.0 */
	int blackbg = 0;       		/* NZ if whole screen should be filled with black */
	int nowin = 0;				/* Don't create test window */
	int ramd = 0;				/* Just test ramdac */
	int fade = 0;				/* Test greyramp fade */
	int donat = 0;				/* Test native output */
	int inf = 0;				/* Infnite patches flag */
	int clear = 0;				/* Clear any display calibration (any calname is ignored) */
	int noexit = 0;				/* Don't exit after loading a calibration */
	int verify = 0;				/* Verify that calname is currently loaded */
	int setatom = 0;			/* Set X11 ICC_PROFILE atom to profile */
	int loadatom = 0;			/* Load X11 ICC_PROFILE atom profile into LUT */
	int loadfile = 0;			/* Load given profile into LUT */
	unsigned char *atomv = NULL;	/* Profile loaded from/to atom */
	char calname[MAXNAMEL+1] = "\000";	/* Calibration file name */
	dispwin *dw;
	unsigned int seed = 0x56781234;
	int i, j;
	ramdac *or, *r;

	error_program = "Dispwin";

	/* Process the arguments */
	mfa = 0;        /* Minimum final arguments */
	for(fa = 1;fa < argc;fa++) {

		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-') {	/* Look for any flags */
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1+mfa) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage();

			else if (argv[fa][1] == 'v')
				verb = 1;

			/* Display number */
			else if (argv[fa][1] == 'd') {
#if defined(UNIX) && !defined(__APPLE__)
				int ix, iv;

				if (strcmp(&argv[fa][2], "isplay") == 0 || strcmp(&argv[fa][2], "ISPLAY") == 0) {
					if (++fa >= argc || argv[fa][0] == '-') usage();
					setenv("DISPLAY", argv[fa], 1);
				} else {
					if (na == NULL) usage();
					fa = nfa;
					if (sscanf(na, "%d,%d",&ix,&iv) != 2) {
						ix = atoi(na);
						iv = 0;
					}
					if ((disp = get_a_display(ix-1)) == NULL)
						usage();
					if (iv > 0)
						disp->rscreen = iv-1;
				}
#else
				int ix;
				if (na == NULL) usage();
				fa = nfa;
				ix = atoi(na);
				if ((disp = get_a_display(ix-1)) == NULL)
					usage();
#endif
			}

			/* Test patch offset and size */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				fa = nfa;
				if (na == NULL) usage();
				if (sscanf(na, " %lf,%lf,%lf ", &ho, &vo, &patscale) != 3)
					usage();
				if (ho < 0.0 || ho > 1.0
				 || vo < 0.0 || vo > 1.0
				 || patscale <= 0.0 || patscale > 50.0)
					usage();
				ho = 2.0 * ho - 1.0;
				vo = 2.0 * vo - 1.0;

			/* Black background */
			} else if (argv[fa][1] == 'B') {
				blackbg = 1;

			} else if (argv[fa][1] == 'i' || argv[fa][1] == 'I')
				inf = 1;

			else if (argv[fa][1] == 'm' || argv[fa][1] == 'M')
				inf = 2;

			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F')
				fade = 1;

			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R')
				ramd = 1;

			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N')
				donat = 1;

			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C')
				clear = 1;

			else if (argv[fa][1] == 'x' || argv[fa][1] == 'X')
				noexit = 1;

			else if (argv[fa][1] == 'V')
				verify = 1;

#if defined(UNIX) && !defined(__APPLE__)
			else if (argv[fa][1] == 'S')
				setatom = 1;

			else if (argv[fa][1] == 'L')
				loadatom = 1;
#endif  /* UNXI X11 */
			else
				usage();
		}
		else
			break;
	}

	if (disp == NULL && (disp = get_a_display(0)) == NULL) {
		error("Unable to open the default display");
	}

	/* See if there's a calibration file */
	if (fa < argc) {
		strncpy(calname,argv[fa++],MAXNAMEL); calname[MAXNAMEL] = '\000';
		if (setatom == 0 && loadatom == 0 && verify == 0)
			loadfile = 1;
	}

	/* Bomb on bad combinations (not all are being detected) */
	if (setatom && calname[0] == '\000')
		error("Can't set X11 atom without a profile argument");

	if (verify && calname[0] == '\000' && loadatom == 0)
		error("No calibration/profile provided to verify against");

	if (loadatom && setatom)
		error("Can't load from X11 atom and set atom at the same time");

	if (verify && setatom)
		error("Can't verify and set X11 atom at the same time");

	/* Don't create a window if it won't be used */
	if (ramd != 0 || clear != 0 || verify != 0 || loadfile != 0 || setatom != 0 || loadatom != 0)
		nowin = 1;

	if (verb)
		printf("About to open dispwin object on the display\n");

	if ((dw = new_dispwin(disp, 100.0 * patscale, 100.0 * patscale, ho, vo, nowin, donat, blackbg, 1)) == NULL) {
		printf("Error - new_dispwin failed!\n");
		return -1;
	}

	/* Clear the display calibration curve */
	if (clear != 0) {
		
		if ((r = dw->get_ramdac(dw)) == NULL) {
			error("We don't have access to the VideoLUT");
		}

		for (i = 0; i < r->nent; i++) {
			double iv = i/(r->nent-1.0);
			r->v[0][i] = iv;
			r->v[1][i] = iv;
			r->v[2][i] = iv;
		}
		if (verb)
			printf("About to clear the calibration\n");
		if (dw->set_ramdac(dw,r)) {
			error("Failed to set ramdac");
		}

		if (noexit) {
			for (;;) {
				sleep(1000);
			}
		}
		r->del(r);

#if defined(UNIX) && !defined(__APPLE__)
	/* Read in the ICC profile, then set the X11 atom value */
	} else if (setatom != 0) {
		icmFile *rd_fp = NULL;
		icc *icco = NULL;
		FILE *fp;
		unsigned long psize, bread;
		Display *mydisplay;
		Atom icc_atom;
		
		/* Open up the profile for reading */
		if ((rd_fp = new_icmFileStd_name(calname,"r")) == NULL)
			error("Can't open file '%s'",calname);

		if ((icco = new_icc()) == NULL)
			error("Creation of ICC object failed");

		/* Read header etc. */
		if (icco->read(icco, rd_fp,0) != 0) 		/* Read ICC OK */
			error("File '%s' doesn't seem to be an ICC profile!",calname);

		icco->del(icco);
		rd_fp->del(rd_fp);

#if defined(O_BINARY) || defined(_O_BINARY)
		if ((fp = fopen(calname,"rb")) == NULL)
#else
		if ((fp = fopen(calname,"r")) == NULL)
#endif
			error ("Can't open file '%s'",calname);

		/* Figure out how big it is */
		if (fseek(fp, 0, SEEK_END))
			error ("Seek '%s' to EOF failed",calname);
		psize = (unsigned long)ftell(fp);
	
		if (fseek(fp, 0, SEEK_SET))
			error ("Seek '%s' to SOF failed",calname);
	
		if ((atomv = (unsigned char *)malloc(psize)) == NULL)
			error("Failed to allocate buffer for profile '%s'",calname);
	
		if ((bread = fread(atomv, 1, psize, fp)) != psize)
			error("Failed to read profile '%s' into buffer",calname);
		
		fclose(fp);

		if ((mydisplay = XOpenDisplay(disp->name)) == NULL)
			error("Unable to open display '%s'",disp->name);

		if ((icc_atom = XInternAtom(mydisplay, disp->icc_atom_name, False)) == None)
			error("Unable to intern atom '%s'",disp->icc_atom_name);

		XChangeProperty(mydisplay, RootWindow(mydisplay, 0), icc_atom,
		                XA_CARDINAL, 8, PropModeReplace, atomv, psize);

		XCloseDisplay(mydisplay);

		free(atomv);
		atomv = NULL;

#endif  /* UNXI X11 */

	/* Setup or verify a display calibration curve set if we are given one */
	/* Note this won't be permanent on OSX. To fix this, a special */
	/* case would have to be made, to load the given profile using */
	/* the appropriate OSX call. */
	} else if (loadfile != 0 || verify != 0 || loadatom != 0) {
		icmFile *rd_fp = NULL;
		icc *icco = NULL;
		cgats *ccg = NULL;			/* calibration cgats structure */
		
		if ((r = dw->get_ramdac(dw)) == NULL) {
			error("We don't have access to the VideoLUT");
		}

#if defined(UNIX) && !defined(__APPLE__)
		if (loadatom) {
			Display *mydisplay;
			Atom icc_atom, ret_type;
			int ret_format;
			long ret_len, ret_togo;

			if ((mydisplay = XOpenDisplay(disp->name)) == NULL)
				error("Unable to open display '%s'",disp->name);

			if ((icc_atom = XInternAtom(mydisplay, disp->icc_atom_name, False)) == None)
				error("Unable to intern atom '%s'",disp->icc_atom_name);

			/* Get the ICC profile property */
			if (XGetWindowProperty(mydisplay, RootWindow(mydisplay, 0), icc_atom,
			            0, 0x7ffffff, False, XA_CARDINAL, 
                        &ret_type, &ret_format, &ret_len, &ret_togo, &atomv) != Success || ret_len == 0)
				error("Getting property '%s' failed",disp->icc_atom_name); 

			XCloseDisplay(mydisplay);
	
			if ((rd_fp = new_icmFileMem((void *)atomv, ret_len)) == NULL)
				error("Creating memory file from X11 atom failed");

			strcpy(calname, disp->icc_atom_name);

		} else
#endif	  /* UNXI X11 */

		{
			/* Open up the profile for reading */
			if ((rd_fp = new_icmFileStd_name(calname,"r")) == NULL)
				error("Can't open file '%s'",calname);
		} 

		if ((icco = new_icc()) == NULL)
			error("Creation of ICC object failed");

		/* Read header etc. */
		if (icco->read(icco, rd_fp,0) == 0) {		/* Read ICC OK */
			icmVideoCardGamma *wo;
			double iv;

			if ((wo = (icmVideoCardGamma *)icco->read_tag(icco, icSigVideoCardGammaTag)) == NULL)
				error("Profile '%s' has no vcgt tag",calname);

			if (wo->u.table.channels == 3) {
				for (i = 0; i < r->nent; i++) {
					iv = i/(r->nent-1.0);
					r->v[0][i] = wo->lookup(wo, 0, iv);
					r->v[1][i] = wo->lookup(wo, 1, iv);
					r->v[2][i] = wo->lookup(wo, 2, iv);
				}
			} else if (wo->u.table.channels == 1) {
				for (i = 0; i < r->nent; i++) {
					iv = i/(r->nent-1.0);
					r->v[0][i] = 
					r->v[1][i] = 
					r->v[2][i] = wo->lookup(wo, 0, iv);
				}
			} else {
				error("Profile '%s' vcgt tag doesn't have 1 or 3 channels",calname);
			}
		} else {	/* See if it's a .cal file */
			int ncal;
			int ii, fi, ri, gi, bi;
			
			icco->del(icco);			/* Don't need these now */
			icco = NULL;
			rd_fp->del(rd_fp);
			rd_fp = NULL;

			ccg = new_cgats();			/* Create a CGATS structure */
			ccg->add_other(ccg, "CAL"); /* our special calibration type */
		
			if (ccg->read_name(ccg, calname)) {
				ccg->del(ccg);
				error("File '%s' is not a valid ICC profile or Argyll .cal file",calname);
			}
		
			if (ccg->ntables == 0 || ccg->t[0].tt != tt_other || ccg->t[0].oi != 0)
				error("Calibration file isn't a CAL format file");
			if (ccg->ntables < 1)
				error("Calibration file '%s' doesn't contain at least one table",calname);
		
			if ((ncal = ccg->t[0].nsets) <= 0)
				error("No data in set of file '%s'",calname);
		
			if (ncal != 256)
				error("Expect 256 data sets in file '%s'",calname);
		
			if ((fi = ccg->find_kword(ccg, 0, "DEVICE_CLASS")) < 0)
				error("Calibration file '%s' doesn't contain keyword COLOR_REPS",calname);
			if (strcmp(ccg->t[0].kdata[fi],"DISPLAY") != 0)
				error("Calibration file '%s' doesn't have DEVICE_CLASS of DISPLAY",calname);

			if ((ii = ccg->find_field(ccg, 0, "RGB_I")) < 0)
				error("Calibration file '%s' doesn't contain field RGB_I",calname);
			if (ccg->t[0].ftype[ii] != r_t)
				error("Field RGB_R in file '%s' is wrong type",calname);
			if ((ri = ccg->find_field(ccg, 0, "RGB_R")) < 0)
				error("Calibration file '%s' doesn't contain field RGB_R",calname);
			if (ccg->t[0].ftype[ri] != r_t)
				error("Field RGB_R in file '%s' is wrong type",calname);
			if ((gi = ccg->find_field(ccg, 0, "RGB_G")) < 0)
				error("Calibration file '%s' doesn't contain field RGB_G",calname);
			if (ccg->t[0].ftype[gi] != r_t)
				error("Field RGB_G in file '%s' is wrong type",calname);
			if ((bi = ccg->find_field(ccg, 0, "RGB_B")) < 0)
				error("Calibration file '%s' doesn't contain field RGB_B",calname);
			if (ccg->t[0].ftype[bi] != r_t)
				error("Field RGB_B in file '%s' is wrong type",calname);
			for (i = 0; i < ncal; i++) {
				r->v[0][i] = *((double *)ccg->t[0].fdata[i][ri]);
				r->v[1][i] = *((double *)ccg->t[0].fdata[i][gi]);
				r->v[2][i] = *((double *)ccg->t[0].fdata[i][bi]);
			}
		}

		/* We've loaded r with the contents of calname */
		if (verb)
			printf("About to set given calibration\n");
		if (verify) {
			int ver = 1;
			double berr = 0.0;
			if ((or = dw->get_ramdac(dw)) == NULL)
				error("Unable to get current VideoLUT for verify");
		
			for (j = 0; j < 3; j++) {
				for (i = 0; i < r->nent; i++) {
					double err;
					err = fabs(r->v[j][i] - or->v[j][i]);
					if (err > berr)
						berr = err;
					if (err > VERIFY_TOL) {
						ver = 0;
					}
				}
			}
			if (ver)
				printf("Verify: '%s' IS loaded (discrepancy %.1f%%)\n", calname, berr * 100);
			else
				printf("Verify: '%s' is NOT loaded (discrepancy %.1f%%)\n", calname, berr * 100);
			or->del(or);
		} else {
			if (dw->set_ramdac(dw,r))
				error("Failed to set VideoLUT");
		}
		r->del(r);

		if (ccg != NULL)
			ccg->del(ccg);
		if (icco != NULL)
			icco->del(icco);
		if (rd_fp != NULL)
			rd_fp->del(rd_fp);

#if defined(UNIX) && !defined(__APPLE__)
		if (loadatom && atomv != NULL)
			XFree(atomv);
#endif	  /* UNXI X11 */

		if (verb)
			printf("Calibration set\n");

		if (noexit) {
			for (;;) {
				sleep(1000);
			}
		}

	/* Window or VideoLUT test */
	} else {

		if (ramd == 0) {

			if (fade) {
				int i;
				int steps = 256;
				for (i = 0; i < steps; i++) {
					double tt;
					tt = i/(steps - 1.0);
					dw->set_color(dw, tt, tt, tt);
					msec_sleep(20);
					printf("Val = %f\n",tt);
				}
			} else {

				if (inf == 2)
					printf("\nHit return to advance each color\n");

				printf("Setting White\n");
				dw->set_color(dw, 1.0, 1.0, 1.0);	/* White */

				if (inf == 2)
					getchar();				
				else
					sleep(2);

				printf("Setting 75%% Grey\n");
				dw->set_color(dw, 0.75, 0.75, 0.75);	/* Grey */

				if (inf == 2)
					getchar();				
				else
					sleep(2);

				printf("Setting 50%% Grey\n");
				dw->set_color(dw, 0.5, 0.5, 0.5);	/* Grey */

				if (inf == 2)
					getchar();				
				else
					sleep(2);

				printf("Setting 25%% Grey\n");
				dw->set_color(dw, 0.25, 0.25, 0.25);	/* Grey */

				if (inf == 2)
					getchar();				
				else
					sleep(2);

				printf("Setting 12.5%% Grey\n");
				dw->set_color(dw, 0.125, 0.125, 0.125);	/* Grey */

				if (inf == 2)
					getchar();				
				else
					sleep(2);

				printf("Setting Black\n");
				dw->set_color(dw, 0.0, 0.0, 0.0);

				if (inf == 2)
					getchar();				
				else
					sleep(2);

				printf("Setting Red\n");
				dw->set_color(dw, 1.0, 0.0, 0.0);	/* Red */

				if (inf == 2)
					getchar();				
				else
					sleep(2);

				printf("Setting Green\n");
				dw->set_color(dw, 0.0, 1.0,  0.0);	/* Green */

				if (inf == 2)
					getchar();				
				else
					sleep(2);

				printf("Setting Blue\n");
				dw->set_color(dw, 0.0, 0.0, 1.0);	/* Blue */

				if (inf == 2)
					getchar();				
				else
					sleep(2);

				printf("Setting Cyan\n");
				dw->set_color(dw, 0.0, 1.0, 1.0);	/* Cyan */

				if (inf == 2)
					getchar();				
				else
					sleep(2);

				printf("Setting Magenta\n");
				dw->set_color(dw, 1.0, 0.0,  1.0);	/* Magenta */

				if (inf == 2)
					getchar();				
				else
					sleep(2);

				printf("Setting Yellow\n");
				dw->set_color(dw, 1.0, 1.0, 0.0);	/* Yellow */

				if (inf == 2)
					getchar();				
				else
					sleep(2);


				if (inf == 1) {
					for (;inf != 0;) {
						double col[3];
	
						for (i = 0; i < 3; i++) {
							seed = PSRAND32(seed);
							col[i] = seed/4294967295.0;
						}
	
						printf("Setting %f %f %f\n",col[0],col[1],col[2]);
						dw->set_color(dw, col[0],col[1],col[2]);
	
						if (inf == 2)
							getchar();
						else
							sleep(2);
	
					}
				}
			}
		}

		if (inf != 2) {
			/* Test out the VideoLUT access */
			if ((or = dw->get_ramdac(dw)) != NULL) {
				
				r = or->clone(or);
	
				/* Try darkening it */
				for (j = 0; j < 3; j++) {
					for (i = 0; i < r->nent; i++) {
						r->v[j][i] = pow(or->v[j][i], 2.0);
					}
				}
				printf("Darkening screen\n");
				if (dw->set_ramdac(dw,r)) {
					dw->set_ramdac(dw,or);
					error("Failed to set ramdac");
				}
				sleep(1);
	
				/* Try lightening it */
				for (j = 0; j < 3; j++) {
					for (i = 0; i < r->nent; i++) {
						r->v[j][i] = pow(or->v[j][i], 0.5);
					}
				}
				printf("Lightening screen\n");
				if (dw->set_ramdac(dw,r)) {
					dw->set_ramdac(dw,or);
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
				printf("We don't have access to the VideoLUT\n");
			}
	
			/* Test out the beeps */
			printf("Normal beep\n");
			normal_beep();
	
			sleep(1);
		
			printf("Good beep\n");
			good_beep();
	
			sleep(1);
		
			printf("Bad double beep\n");
			bad_beep();
	
			sleep(2);		/* Allow beep to complete */
		}
	}
	
	free_a_disppath(disp);

	if (verb)
		printf("About to destroy dispwin object\n");

	dw->del(dw);

	return 0;
}

#endif /* STANDALONE_TEST */

/* ---------------------------------------------------------------- */
/* Unused code */

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
			debug("malloc failed for disp mode keys\n");
			free_disppaths(disps);
			free(dids);
			return NULL;
		}
		if ((values = (void **)malloc(nde * sizeof(void *))) == NULL) {
			debug("malloc failed for disp mode values\n");
			free(keys);
			free_disppaths(disps);
			free(dids);
			return NULL;
		}
		CFDictionaryGetKeysAndValues(dr, (const void **)keys, (const void **)values);
		for (j = 0; j < nde; j++) {
			printf("Entry %d key = %s\n", j, CFStringGetCStringPtr(keys[j], kCFStringEncodingMacRoman));
		}
		free(values);
		free(keys);
	}
#endif
