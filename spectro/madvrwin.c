

/* 
 * Argyll Color Correction System
 * Web Display target patch window
 *
 * Author: Graeme W. Gill
 * Date:   3/4/12
 *
 * Copyright 2013 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#include <stdio.h>
#include <string.h>
#ifdef NT
# include <winsock2.h>
#endif
#ifdef UNIX
# include <sys/types.h>
# include <ifaddrs.h>
# include <netinet/in.h> 
# include <arpa/inet.h>
# ifdef __FreeBSD__
#  include <sys/socket.h>
# endif /* __FreeBSD__ */
#endif
#include "copyright.h"
#include "aconfig.h"
#include "icc.h"
#include "numsup.h"
#include "cgats.h"
#include "conv.h"
#include "dispwin.h"
#include "madvrwin.h"
#include "conv.h"
#include "mongoose.h"

#define ENABLE_RAMDAC

#undef DEBUG
//#define STANDALONE_TEST

#ifdef DEBUG
#define errout stderr
# define debug(xx)	fprintf(errout, xx )
# define debug2(xx)	fprintf xx
# define debugr(xx)	fprintf(errout, xx )
# define debugr2(xx)	fprintf xx
# define debugrr(xx)	fprintf(errout, xx )
# define debugrr2(xx)	fprintf xx
# define debugrr2l(lev, xx)	fprintf xx
#else
#define errout stderr
# define debug(xx) 
# define debug2(xx)
# define debugr(xx) if (p->ddebug) fprintf(errout, xx ) 
# define debugr2(xx) if (p->ddebug) fprintf xx
# define debugrr(xx) if (callback_ddebug) fprintf(errout, xx ) 
# define debugrr2(xx) if (callback_ddebug) fprintf xx
# define debugrr2l(lev, xx) if (callback_ddebug >= lev) fprintf xx
#endif

/* ----------------------------------------------- */
/* MadVR functions */

#ifndef KEY_WOW64_32KEY
# define KEY_WOW64_32KEY (0x0200)
#endif

HMODULE HcNetDll = NULL;
BOOL (WINAPI *madVR_BlindConnect)(BOOL searchLan, DWORD timeOut) = NULL;
BOOL (WINAPI *madVR_SetOsdText)(LPCWSTR text);
BOOL (WINAPI *madVR_Disable3dlut)() = NULL;
BOOL (WINAPI *madVR_GetDeviceGammaRamp)(LPVOID ramp) = NULL;
BOOL (WINAPI *madVR_SetDeviceGammaRamp)(LPVOID ramp) = NULL;
BOOL (WINAPI *madVR_SetBackground)(int patternAreaInPercent, COLORREF backgroundColor) = NULL;
BOOL (WINAPI *madVR_ShowRGB)(double r, double g, double b) = NULL;
BOOL (WINAPI *madVR_SetProgressBarPos)(int currentPos, int maxPos);
BOOL (WINAPI *madVR_Disconnect)() = NULL;

/* Locate and populate the MadVR functions */
/* Return NZ on failure */
static int initMadVR(dispwin *p) {

	if ((HcNetDll = LoadLibraryW(L"madHcNet32.dll")) == NULL) {
		HKEY hk1;
    	if (RegOpenKeyExW(HKEY_CLASSES_ROOT, L"CLSID\\{E1A8B82A-32CE-4B0D-BE0D-AA68C772E423}\\InprocServer32", 0, KEY_QUERY_VALUE | KEY_WOW64_32KEY, &hk1) == ERROR_SUCCESS) {
			DWORD size;
			LPWSTR us1;
			int i1;
			size = MAX_PATH * 2 + 2;
		    us1 = (LPWSTR) LocalAlloc(LPTR, size + 20);
		    i1 = RegQueryValueExW(hk1, NULL, NULL, NULL, (LPBYTE) us1, &size);
			if (i1 == ERROR_MORE_DATA) {
				LocalFree(us1);
				us1 = (LPWSTR) LocalAlloc(LPTR, size + 20);
				i1 = RegQueryValueExW(hk1, NULL, NULL, NULL, (LPBYTE) us1, &size);
			}
			if (i1 == ERROR_SUCCESS) {
				for (i1 = lstrlenW(us1) - 2; i1 > 0; i1--)
					if (us1[i1] == L'\\') {
						us1[i1 + 1] = 0;
						break;
					}
				wcscat(us1, L"madHcNet32.dll");
				HcNetDll = LoadLibraryW(us1);
			}
			LocalFree(us1);
			RegCloseKey(hk1);
		}
	}
	if (HcNetDll != NULL) {
		*(FARPROC*)&madVR_BlindConnect       = GetProcAddress(HcNetDll, "madVR_BlindConnect");
		*(FARPROC*)&madVR_SetOsdText         = GetProcAddress(HcNetDll, "madVR_SetOsdText");
		*(FARPROC*)&madVR_Disable3dlut       = GetProcAddress(HcNetDll, "madVR_Disable3dlut");
		*(FARPROC*)&madVR_GetDeviceGammaRamp = GetProcAddress(HcNetDll, "madVR_GetDeviceGammaRamp");
		*(FARPROC*)&madVR_SetDeviceGammaRamp = GetProcAddress(HcNetDll, "madVR_SetDeviceGammaRamp");
		*(FARPROC*)&madVR_SetBackground      = GetProcAddress(HcNetDll, "madVR_SetBackground");
		*(FARPROC*)&madVR_ShowRGB            = GetProcAddress(HcNetDll, "madVR_ShowRGB");
		*(FARPROC*)&madVR_SetProgressBarPos  = GetProcAddress(HcNetDll, "madVR_SetProgressBarPos");
		*(FARPROC*)&madVR_Disconnect         = GetProcAddress(HcNetDll, "madVR_Disconnect");
	
		if (madVR_BlindConnect
		 && madVR_SetOsdText
		 && madVR_Disable3dlut
		 && madVR_GetDeviceGammaRamp
		 && madVR_SetDeviceGammaRamp
		 && madVR_SetBackground
		 && madVR_ShowRGB
		 && madVR_SetProgressBarPos
		 && madVR_Disconnect) {
			return 0;
		}
		debugr2((errout,"Failed to locate function in madHcNet32.dll\n"));
		FreeLibrary(HcNetDll);
	} else {
		debugr2((errout,"Failed to load madHcNet32.dll\n"));
	}
	return 1;
}

static void deinitMadVR(dispwin *p) {
	if (HcNetDll != NULL) {
		FreeLibrary(HcNetDll);
		HcNetDll = NULL;
	}
}

/* ----------------------------------------------- */

/* Get RAMDAC values. ->del() when finished. */
/* Return NULL if not possible */
static ramdac *madvrwin_get_ramdac(dispwin *p) {
	ramdac *r = NULL;

#ifdef ENABLE_RAMDAC
	WORD vals[3][256];		/* 256 x 16 bit elements (Quantize) */
	int i, j;

	debugr("madvrwin_get_ramdac called\n");

	/* Allocate a ramdac */
	if ((r = (ramdac *)calloc(sizeof(ramdac), 1)) == NULL) {
		debugr("madvrwin_get_ramdac failed on malloc()\n");
		return NULL;
	}
	r->pdepth = p->pdepth;
	r->nent = (1 << p->pdepth);
	r->clone =  dispwin_clone_ramdac;
	r->setlin = dispwin_setlin_ramdac;
	r->del =    dispwin_del_ramdac;

	for (j = 0; j < 3; j++) {

		if ((r->v[j] = (double *)calloc(sizeof(double), r->nent)) == NULL) {
			for (j--; j >= 0; j--)
				free(r->v[j]);
			free(r);
			debugr("madvrwin_get_ramdac failed on malloc()\n");
			return NULL;
		}
	}

	/* GetDeviceGammaRamp() is hard coded for 3 x 256 entries (Quantize) */
	if (r->nent != 256) {
		free(r);
		debugr2((errout,"GetDeviceGammaRamp() is hard coded for nent == 256, and we've got nent = %d!\n",r->nent));
		return NULL;
	}

	if (madVR_GetDeviceGammaRamp(vals) == 0) {
		free(r);
		debugr2((errout,"madvrwin_get_ramdac failed on madVR_GetDeviceGammaRamp()\n"));
		return NULL;
	}

	for (j = 0; j < 3; j++) {
		for (i = 0; i < r->nent; i++) {
			r->v[j][i] = vals[j][i]/65535.0;
		}
	}
#endif // ENABLE_RAMDAC
	debugr2((errout,"madvrwin_get_ramdac returning 0x%x\n",r));

	return r;
}

/* Set the RAMDAC values. */
/* Return nz if not possible */
static int madvrwin_set_ramdac(dispwin *p, ramdac *r, int persist) {
	int rv = 1;

#ifdef ENABLE_RAMDAC
	int i, j;
	WORD vals[3][256];		/* 16 bit elements */

	debugr("madvrwin_set_ramdac called\n");

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

	if (madVR_SetDeviceGammaRamp(vals) == 0) {
		debugr2((errout,"madvrwin_set_ramdac failed on madVR_SetDeviceGammaRamp()\n"));
		rv = 1;
	} else {
		debugr("madvrwin_set_ramdac set\n");
		rv = 0;
	}
#endif // ENABLE_RAMDAC
	return rv;
}

/* ----------------------------------------------- */
/* Install a display profile and make */
/* it the default for this display. */
/* Return nz if failed */
int madvrwin_install_profile(dispwin *p, char *fname, ramdac *r, p_scope scope) {
	debugr("madVRdisp doesn't support installing profiles\n"); 
	return 1;
}

/* Un-Install a display profile */
/* Return nz if failed, */
int madvrwin_uninstall_profile(dispwin *p, char *fname, p_scope scope) {
	debugr("madVRdisp doesn't support uninstalling profiles\n"); 
	return 1;
}

/* Get the currently installed display profile. */
/* Return NULL if failed. */
icmFile *madvrwin_get_profile(dispwin *p, char *name, int mxlen) {
	debugr("madVRdisp doesn't support getting the current profile\n"); 
	return NULL;
}

/* ----------------------------------------------- */

/* Change the window color. */
/* Return 1 on error, 2 on window being closed */
static int madvrwin_set_color(
dispwin *p,
double r, double g, double b	/* Color values 0.0 - 1.0 */
) {
	int j;
	double orgb[3];		/* Previous RGB value */
	double kr, kf;
	int update_delay = p->update_delay; 
	double xdelay = 0.0;		/* Extra delay for response time */

	debugr("madvrwin_set_color called\n");

	if (p->nowin)
		return 1;

	orgb[0] = p->rgb[0]; p->rgb[0] = r;
	orgb[1] = p->rgb[1]; p->rgb[1] = g;
	orgb[2] = p->rgb[2]; p->rgb[2] = b;

	if (!madVR_ShowRGB(p->rgb[0], p->rgb[1], p->rgb[2])) {
		debugr2((errout,"madVR_ShowRGB failed\n"));
		return 1;
	}

	/* Don't want extra delay if we're measuring update delay */
	if (update_delay != 0 && p->do_resp_time_del) {
		/* Compute am expected response time for the change in level */
		kr = DISPLAY_RISE_TIME/log(1 - 0.9);	/* Exponent constant */
		kf = DISPLAY_FALL_TIME/log(1 - 0.9);	/* Exponent constant */
//printf("~1 k2 = %f\n",k2);
		for (j = 0; j < 3; j++) {
			double el, dl, n, t;
	
			el = pow(p->rgb[j], 2.2);
			dl = el - pow(orgb[j], 2.2);	/* Change in level */
			if (fabs(dl) > 0.01) {		/* More than 1% change in level */
				n = DISPLAY_SETTLE_AIM * el;
				if (n < DISPLAY_ABS_AIM)
					n = DISPLAY_ABS_AIM;
//printf("~1 sl %f, el %f, log (%f / %f)\n",sl,el,n,fabs(sl - el));
				if (dl > 0.0)
					t = kr * log(n/dl);
				else
					t = kf * log(n/-dl);
	
				if (t > xdelay)
					xdelay = t;
			}
		}
//printf("~1 xdelay = %f secs\n",xdelay);
		xdelay *= 1000.0;		/* To msec */
		/* This is kind of a fudge since update delay is after latency, */
		/* but displays with long delay (ie. CRT) have short latency, and visa versa */
		if ((int)xdelay > update_delay)
			update_delay = (int)xdelay;
	}

	/* Allow some time for the display to update before */
	/* a measurement can take place. This allows for CRT */
	/* refresh, or LCD processing/update time, + */
	/* display settling time (quite long for smaller LCD changes). */
	msec_sleep(update_delay);

	return 0;
}

/* ----------------------------------------------- */
/* set patch info */
/* Return nz on error */
static int madvrwin_set_pinfo(dispwin *p, int pno, int tno) {

	if (!madVR_SetProgressBarPos(pno, tno))
		return 1;

	return 0;
}

/* ----------------------------------------------- */
/* Set an update delay, and return the previous value */
/* Value can be set to zero, but othewise will be forced */
/* to be >= min_update_delay */
static int madvrwin_set_update_delay(
dispwin *p,
int update_delay) {
	int cval = p->update_delay;
	p->update_delay = update_delay;
	if (update_delay != 0 && p->update_delay < p->min_update_delay)
		p->update_delay = p->min_update_delay;
	return cval;
}

/* ----------------------------------------------- */
/* Set the shell set color callout */
void madvrwin_set_callout(
dispwin *p,
char *callout
) {
	debugr2((errout,"madvrwin_set_callout called with '%s'\n",callout));

	p->callout = strdup(callout);
}

/* ----------------------------------------------- */
/* Destroy ourselves */
static void madvrwin_del(
dispwin *p
) {

	debugr("madvrwin_del called\n");

	if (p == NULL)
		return;

	deinitMadVR(p);

	if (p->name != NULL)
		free(p->name);
	if (p->description != NULL)
		free(p->description);
	if (p->callout != NULL)
		free(p->callout);

	/* Since we don't restore the display, delete these here */
	if (p->oor != NULL) {
		p->oor->del(p->oor);
		p->oor = NULL;
	}
	if (p->or != NULL) {
		p->or->del(p->or);
		p->or = NULL;
	}
	if (p->r != NULL) {
		p->r->del(p->r);
		p->r = NULL;
	}

	free(p);
}

/* ----------------------------------------------- */

/* Create a web display test window, default grey */
dispwin *new_madvrwin(
double width, double height,	/* Width and height in mm - turned into % of screen */
double hoff, double voff,		/* Offset from center in fraction of screen - ignored */
int nowin,						/* NZ if no window should be created - RAMDAC access only */
int native,						/* X0 = use current per channel calibration curve */
								/* X1 = set native linear output and use ramdac high precn. */
								/* 0X = use current color management cLut (MadVR) */
								/* 1X = disable color management cLUT (MadVR) */
int *noramdac,					/* Return nz if no ramdac access. native is set to X0 */
int *nocm,						/* Return nz if no CM cLUT access. native is set to 0X */
int out_tvenc,					/* 1 = use RGB Video Level encoding */
int blackbg,					/* NZ if whole screen should be filled with black */
int verb,						/* NZ for verbose prompts */
int ddebug						/* >0 to print debug statements to stderr */
) {
	dispwin *p = NULL;
	char *cp;
	struct mg_context *mg;
	const char *options[3];
	char port[50];

	debug("new_madvrwin called with native = %d\n");

	if (out_tvenc)
		return NULL;

	if ((p = (dispwin *)calloc(sizeof(dispwin), 1)) == NULL) {
		if (ddebug) fprintf(stderr,"new_madvrwin failed because malloc failed\n");
		return NULL;
	}

	/* !!!! Make changes in dispwin.c & webwin.c as well !!!! */
	p->name = strdup("Web Window");
	p->nowin = nowin;
	p->native = native;
	p->out_tvenc = 0;
	p->blackbg = blackbg;
	p->ddebug = ddebug;
	p->get_ramdac        = madvrwin_get_ramdac;
	p->set_ramdac        = madvrwin_set_ramdac;
	p->install_profile   = madvrwin_install_profile;
	p->uninstall_profile = madvrwin_uninstall_profile;
	p->get_profile       = madvrwin_get_profile;
	p->set_color         = madvrwin_set_color;
	p->set_pinfo         = madvrwin_set_pinfo;
	p->set_update_delay  = madvrwin_set_update_delay;
	p->set_callout       = madvrwin_set_callout;
	p->del               = madvrwin_del;

	debugr2((errout, "new_madvrwin got native = %d\n",native));

#ifndef ENABLE_RAMDAC
	if (noramdac != NULL)
		*noramdac = 1;
	p->native &= ~1;
#endif

	p->rgb[0] = p->rgb[1] = p->rgb[2] = 0.5;	/* Set Grey as the initial test color */

	p->min_update_delay = 20;

	if ((cp = getenv("ARGYLL_MIN_DISPLAY_UPDATE_DELAY_MS")) != NULL) {
		p->min_update_delay = atoi(cp);
		if (p->min_update_delay < 20)
			p->min_update_delay = 20;
		if (p->min_update_delay > 60000)
			p->min_update_delay = 60000;
		debugr2((errout, "new_madvrwin: Minimum display update delay set to %d msec\n",p->min_update_delay));
	}

	p->update_delay = DISPLAY_UPDATE_DELAY;		/* Default update delay */
	if (p->update_delay < p->min_update_delay)
		p->update_delay = p->min_update_delay;

	p->pdepth = 8;		/* Assume this */
	p->edepth = 16;

	if (initMadVR(p)) {
		free(p);
		return NULL;
	}

	if (!madVR_BlindConnect(0, 1000)) {
		debugr2((errout,"Failed to locate MadVR\n"));
		free(p);
		return NULL;
	}

	if (p->native & 2) {
		debugr2((errout,"new_madvrwin: disbling 3dLuts\n"));
		madVR_Disable3dlut();
	}

	if (blackbg) {
		int perc;

		perc = (int)(width * height * 100 + 0.5);
		madVR_SetBackground(perc, 0x0000);
	}

	/* Create a suitable description */
	{
		char buf[1000];

		sprintf(buf,"ArgyllCMS Patches");
		p->description = strdup(buf);

		if (verb)
			printf("Created MadVR window\n");

		madVR_SetOsdText(L"ArgyllCMS Patches");
	}

#ifdef ENABLE_RAMDAC

	/* Save the original ramdac, which gets restored on exit */
	if ((p->or = p->get_ramdac(p)) != NULL) {

		debugr("Saved original VideoLUT\n");

		if (noramdac != NULL)
			*noramdac = 0;

		/* Copy original ramdac that never gets altered */
		if ((p->oor = p->or->clone(p->or)) == NULL) {
			madvrwin_del(p);
			debugr("ramdac clone failed - memory ?\n");
			return NULL;
		}

		/* Create a working ramdac for native or other use */
		if ((p->r = p->or->clone(p->or)) == NULL) {
			madvrwin_del(p);
			debugr("ramdac clone failed - memory ?\n");
			return NULL;
		}

		/* Setup for native mode == linear RAMDAC */
		if ((p->native & 1) == 1) {
int ii = 0;
			debug("About to setup native mode\n");

			/* Since we MadVR does dithering, we don't need to use the VideoLUTs */
			/* to emulate higher precision, so simply set it to linear here. */
			if ((ii = madVR_SetDeviceGammaRamp(NULL)) == 0) {
				madvrwin_del(p);
				debugr("Clear gamma ramp failed\n");
				return NULL;
			}
		}

	} else {
		debugr2((errout,"Unable to access VideoLUT\n"));
		if (noramdac != NULL)
			*noramdac = 1;
		p->oor = p->or = p->r = NULL;
	}

	if (!p->nowin) {

		/* Make sure initial test color is displayed */
		madvrwin_set_color(p, p->rgb[0], p->rgb[1], p->rgb[2]);
	}
#endif

	debugr("new_madvrwin: return sucessfully\n");

	return p;
}

