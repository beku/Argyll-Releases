

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
#include "webwin.h"
#include "conv.h"
#include "mongoose.h"

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

// A handler for the /ajax/get_messages endpoint.
// Return a list of messages with ID greater than requested.
static void ajax_get_messages(struct mg_connection *conn,
                              const struct mg_request_info *request_info) {
char src_addr[20];

//	dispwin *p = (dispwin *)(request_info->user_data);
	dispwin *p = (dispwin *)mg_get_user_data(conn);

//sockaddr_to_string(src_addr, sizeof(src_addr), &conn->client.rsa);

//	printf("ajax_messages query_string '%s'\n",request_info->query_string);

	p->ccix++;

	while(p->ncix == p->ccix && p->mg_stop == 0) {
		msec_sleep(50);
	}

	mg_printf(conn,
			"\r\n#%02X%02X%02X",
			(int)(p->r_rgb[0] * 255.0 + 0.5),
			(int)(p->r_rgb[1] * 255.0 + 0.5),
			(int)(p->r_rgb[2] * 255.0 + 0.5));
}

/* Event handler */
static void *webwin_ehandler(enum mg_event event,
                           struct mg_connection *conn) {
//                           const struct mg_request_info *request_info) {
	const struct mg_request_info *request_info = mg_get_request_info(conn);

	if (event != MG_NEW_REQUEST) {
		return NULL;
	}
#ifdef DEBUG
	printf("Got event with uri = '%s'\n",request_info->uri);
#endif
	if (strcmp(request_info->uri, "/ajax/messages") == 0) {
		ajax_get_messages(conn, request_info);
	} else if (strcmp(request_info->uri, "/webdisp.js") == 0) {
#ifndef NEVER
		char *webdisp_js = 
	"\r\n"
	"if (typeof XMLHttpRequest == \"undefined\") {\r\n"
	"	XMLHttpRequest = function () {\r\n"
	"		try { return new ActiveXObject(\"Msxml2.XMLHTTP.6.0\"); }\r\n"
	"			catch (e) {}\r\n"
	"		try { return new ActiveXObject(\"Msxml2.XMLHTTP.3.0\"); }\r\n"
	"			catch (e) {}\r\n"
	"		try { return new ActiveXObject(\"Microsoft.XMLHTTP\"); }\r\n"
	"			catch (e) {}\r\n"
	"		throw new Error(\"This browser does not support XMLHttpRequest.\");\r\n"
	"	};\r\n"
	"}\r\n"
	"\r\n"
	"var ccolor = \"\";\r\n"
	"var oXHR;\r\n"
	"\r\n"
	"function XHR_response() {\r\n"
	"	if (oXHR.readyState != 4)\r\n"
	"		return;\r\n"
	"\r\n"
	"	if (oXHR.status != 200) {\r\n"
	"		return;\r\n"
	"	}\r\n"
	"	var rt = oXHR.responseText;\r\n"
	"	if (rt.charAt(0) == '\\r' && rt.charAt(1) == '\\n')\r\n"
	"		rt = rt.slice(2);\r\n"
	"	if (ccolor != rt) {\r\n"
	"		ccolor = rt;\r\n"
	"		document.body.style.background = ccolor;\r\n"
	"	}\r\n"
	"	oXHR.open(\"GET\", \"/ajax/messages?\" + document.body.style.background + \" \" + Math.random(), true);\r\n"
	"	oXHR.onreadystatechange = XHR_response;\r\n"
	"	oXHR.send();\r\n"
	"}\r\n"
	"\r\n"
	"window.onload = function() {\r\n"
	"	ccolor = \"#808080\";\r\n"
	"	document.body.style.background = ccolor;\r\n"
	"\r\n"
	"	oXHR = new XMLHttpRequest();\r\n"
	"	oXHR.open(\"GET\", \"/ajax/messages?\" + document.body.style.background, true);\r\n"
	"	oXHR.onreadystatechange = XHR_response;\r\n"
	"	oXHR.send();\r\n"
	"};\r\n";
	    mg_write(conn, webdisp_js, strlen(webdisp_js));
#else
		return NULL;	/* Read webdisp.js */
#endif
	} else {
	    mg_printf(conn, "HTTP/1.1 200 OK\r\n"
		"Cache-Control: no-cache\r\n"
		"Content-Type: text/html\r\n\r\n"
		"<html>\r\n"
		"<head>\r\n"
		"<title>ArgyllCMS Web Display</title>\r\n"
		"<script src=\"webdisp.js\"></script>\r\n"
		"</head>\r\n"
		"</html>\r\n"
		);
	}

//		"<script type=\"text/javascript\"src=\"webdisp.js\"></script>"
	return "yes";
}

/* ----------------------------------------------- */

/* Get RAMDAC values. ->del() when finished. */
/* Return NULL if not possible */
static ramdac *webwin_get_ramdac(dispwin *p) {
	debugr("webdisp doesn't have a RAMDAC\n"); 
	return NULL;
}

/* Set the RAMDAC values. */
/* Return nz if not possible */
static int webwin_set_ramdac(dispwin *p, ramdac *r, int persist) {
	debugr("webdisp doesn't have a RAMDAC\n"); 
	return 1;
}

/* ----------------------------------------------- */
/* Install a display profile and make */
/* it the default for this display. */
/* Return nz if failed */
int webwin_install_profile(dispwin *p, char *fname, ramdac *r, p_scope scope) {
	debugr("webdisp doesn't support installing profiles\n"); 
	return 1;
}

/* Un-Install a display profile */
/* Return nz if failed, */
int webwin_uninstall_profile(dispwin *p, char *fname, p_scope scope) {
	debugr("webdisp doesn't support uninstalling profiles\n"); 
	return 1;
}

/* Get the currently installed display profile. */
/* Return NULL if failed. */
icmFile *webwin_get_profile(dispwin *p, char *name, int mxlen) {
	debugr("webdisp doesn't support getting the current profile\n"); 
	return NULL;
}

/* ----------------------------------------------- */

/* Change the window color. */
/* Return 1 on error, 2 on window being closed */
static int webwin_set_color(
dispwin *p,
double r, double g, double b	/* Color values 0.0 - 1.0 */
) {
	int j;
	double orgb[3];		/* Previous RGB value */
	double kr, kf;
	int update_delay = p->update_delay; 
	double xdelay = 0.0;		/* Extra delay for response time */

	debugr("webwin_set_color called\n");

	if (p->nowin)
		return 1;

	orgb[0] = p->rgb[0]; p->rgb[0] = r;
	orgb[1] = p->rgb[1]; p->rgb[1] = g;
	orgb[2] = p->rgb[2]; p->rgb[2] = b;

	for (j = 0; j < 3; j++) {
		if (p->rgb[j] < 0.0)
			p->rgb[j] = 0.0;
		else if (p->rgb[j] > 1.0)
			p->rgb[j] = 1.0;
		p->r_rgb[j] = p->s_rgb[j] = p->rgb[j];
		if (p->out_tvenc) {
			p->r_rgb[j] = p->s_rgb[j] = ((235.0 - 16.0) * p->s_rgb[j] + 16.0)/255.0;

			/* For video encoding the extra bits of precision are created by bit shifting */
			/* rather than scaling, so we need to scale the fp value to account for this. */
			if (p->pdepth > 8)
				p->r_rgb[j] = (p->s_rgb[j] * 255 * (1 << (p->pdepth - 8)))
				            /((1 << p->pdepth) - 1.0); 	
		}
	}

	/* This is probably not actually thread safe... */
	p->ncix++;

	while(p->ncix != p->ccix) {
		msec_sleep(50);
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
/* Set an update delay, and return the previous value */
/* Value can be set to zero, but othewise will be forced */
/* to be >= min_update_delay */
static int webwin_set_update_delay(
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
void webwin_set_callout(
dispwin *p,
char *callout
) {
	debugr2((errout,"webwin_set_callout called with '%s'\n",callout));

	p->callout = strdup(callout);
}

/* ----------------------------------------------- */
/* Destroy ourselves */
static void webwin_del(
dispwin *p
) {

	debugr("webwin_del called\n");

	if (p == NULL)
		return;

	p->mg_stop = 1;
	mg_stop((struct mg_context *)p->pcntx);

	if (p->name != NULL)
		free(p->name);
	if (p->description != NULL)
		free(p->description);
	if (p->callout != NULL)
		free(p->callout);

	free(p);
}

/* ----------------------------------------------- */

/* Create a web display test window, default grey */
dispwin *new_webwin(
int webdisp,					/* Port number */
double width, double height,	/* Width and height in mm */
double hoff, double voff,		/* Offset from center in fraction of screen, range -1.0 .. 1.0 */
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

	debug("new_webwin called\n");

	if ((p = (dispwin *)calloc(sizeof(dispwin), 1)) == NULL) {
		if (ddebug) fprintf(stderr,"new_webwin failed because malloc failed\n");
		return NULL;
	}

	/* !!!! Make changes in dispwin.c & madvrwin.c as well !!!! */
	p->name = strdup("Web Window");
	p->nowin = nowin;
	p->native = native;
	p->out_tvenc = out_tvenc;
	p->blackbg = blackbg;
	p->ddebug = ddebug;
	p->get_ramdac        = webwin_get_ramdac;
	p->set_ramdac        = webwin_set_ramdac;
	p->install_profile   = webwin_install_profile;
	p->uninstall_profile = webwin_uninstall_profile;
	p->get_profile       = webwin_get_profile;
	p->set_color         = webwin_set_color;
	p->set_update_delay  = webwin_set_update_delay;
	p->set_callout       = webwin_set_callout;
	p->del               = webwin_del;

	if (noramdac != NULL)
		*noramdac = 1;
	p->native &= ~1;

	if (nocm != NULL)
		*nocm = 1;
	p->native &= ~2;

	p->rgb[0] = p->rgb[1] = p->rgb[2] = 0.5;	/* Set Grey as the initial test color */

	p->min_update_delay = 20;

	if ((cp = getenv("ARGYLL_MIN_DISPLAY_UPDATE_DELAY_MS")) != NULL) {
		p->min_update_delay = atoi(cp);
		if (p->min_update_delay < 20)
			p->min_update_delay = 20;
		if (p->min_update_delay > 60000)
			p->min_update_delay = 60000;
		debugr2((errout, "new_webwin: Minimum display update delay set to %d msec\n",p->min_update_delay));
	}

	p->update_delay = DISPLAY_UPDATE_DELAY;		/* Default update delay */
	if (p->update_delay < p->min_update_delay)
		p->update_delay = p->min_update_delay;

	p->ncix = 1;

	p->pdepth = 8;		/* Assume this by API */
	p->edepth = 8;

	/* Basic object is initialised, so create a web server */

	options[0] = "listening_ports";
	sprintf(port,"%d", webdisp);
	options[1] = port;
	options[2] = NULL;

	mg = mg_start(&webwin_ehandler, (void *)p, options);
	p->pcntx = (void *)mg;

//printf("Domain = %s'\n",mg_get_option(mg, "authentication_domain"));

	/* Create a suitable description */
#if NT
	{
		char szHostName[255];
		struct hostent *host_entry;
		char *localIP;
		char buf[1000];

		/* We assume WinSock has been started by mongoose */

		// Get the local hostname
		gethostname(szHostName, 255);
		host_entry=gethostbyname(szHostName);
		/* Get first entry */
		localIP = inet_ntoa(*(struct in_addr *)*host_entry->h_addr_list);

		sprintf(buf,"Web Window at http://%s:%d",localIP,webdisp);
		p->description = strdup(buf);

		if (verb)
			printf("Created web server at 'http://%s:%d', now waiting for browser to connect\n",localIP,webdisp);
	}
#else
	{
		struct ifaddrs * ifAddrStruct=NULL;
		struct ifaddrs * ifa=NULL;
		void *tmpAddrPtr=NULL;
		char abuf[INET_ADDRSTRLEN] = "";
		char abuf6[INET6_ADDRSTRLEN] = "";
		char *addr = abuf;
		char buf[1000];
	
		getifaddrs(&ifAddrStruct);
	
		/* Stop at the first non local adderss */
		for (ifa = ifAddrStruct; ifa != NULL; ifa = ifa->ifa_next) {
#ifdef AF_INET6
			if (ifa->ifa_addr->sa_family==AF_INET) { /* IP4 ? */
#endif
				if (strncmp(ifa->ifa_name, "lo",2) == 0 || abuf[0] != '\000')
					continue;
				tmpAddrPtr=&((struct sockaddr_in *)ifa->ifa_addr)->sin_addr;
				inet_ntop(AF_INET, tmpAddrPtr, abuf, INET_ADDRSTRLEN);
//				printf("%s IP Address %s\n", ifa->ifa_name, addressBuffer); 
#ifdef AF_INET6
			} else if (ifa->ifa_addr->sa_family==AF_INET6) { /* IP6 ? */
				if (strncmp(ifa->ifa_name, "lo",2) == 0 || abuf6[0] != '\000')
					continue;
				tmpAddrPtr=&((struct sockaddr_in6 *)ifa->ifa_addr)->sin6_addr;
				inet_ntop(AF_INET6, tmpAddrPtr, abuf6, INET6_ADDRSTRLEN);
//				printf("%s IP Address %s\n", ifa->ifa_name, addressBuffer); 
			} 
#endif
		}
		if (ifAddrStruct!=NULL)
			freeifaddrs(ifAddrStruct);
		if (addr[0] == '\000')
			addr = abuf6;
		if (addr[0] == '\000')
			addr = "Unknown";

		sprintf(buf,"Web Window at http://%s:%d",addr,webdisp);
		p->description = strdup(buf);

		if (verb)
			printf("Created web server at 'http://%s:%d', now waiting for browser to connect\n",addr,webdisp);
	}
#endif

	/* Wait for the web server to connect */
	debugr("new_webwin: waiting for web browser to connect\n");
	while(p->ccix == 0) {
		msec_sleep(50);
	}

	debugr("new_webwin: return sucessfully\n");

	return p;
}

