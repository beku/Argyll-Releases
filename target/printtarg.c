
/* 
 * Argyll Color Correction System
 * PostScript print chart generator module.
 *
 * Author: Graeme W. Gill
 * Date:   28/9/96
 *
 * Copyright 1996 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* This program generates a PostScript print target file. */
/* containing color test patches, given the .ti1 file specifying */
/* what the colors are. */

/* The output is designed to suite a general XY spectrometer (such as
/* the Gretag SpectrScan), a handheld, manual instrument, or */
/* an Xrite DTP51, or an Xrite DTP41 strip spectrometer. */

/* Description:

	This program simply generates a PostScript file containing
	the patches layed out for an Xrite DTP51/DTP41/SpectroScan. It
    allows them to be layed out on a choice of paper sizes,
	with the appropriate contrasting color spacers between
	each patch for the strip reading instruments. Unlike other
	charts, Argyll charts are generated as required, rather
	that being fixed. Also unlike most other strip reading charts,
	the spacers may colored, so that the density contrast ratio is
	guaranteed, even when two patches are about 50% density.

	Another feature is the pseudo random patch layout. This has
	three purposes. One is to try and average out any variation
	in the device response in relationship to the location of
	the patch on the paper. Color copiers and printing presses
	(for instance), are notorious in having side to side density
	variations.

	Another purpose of the random patch layout, is that it gives
	the reading program a good mechanism for detecting user error.
	It can guess the expected values, compare them to the readings,
	and complain if it seems that the strip is probably the wrong
	one. 

	The final purpose of the random patch layout is to optimse the
	contrast between patches in a strip, to improve the robustness
	of the strip reading. Using this, small charts may be even be
	generated without any gaps between the test patches.

 */

/*
 * Nomencalture:
 *
 *	Largely due to how the strip readers name things, the following terms
 *  are used for how patches are grouped:
 *
 *  Pass:	One row of patches in a strip. A pass is usually labeled
 *          with a unique alphabetic label. 
 *  Strip:  A group of passes that can be read by a strip reader.
 *          For an XY instrument, the strip is a complete sheet, and
 *          a each pass is one column. The rows of an XY chart are
 *          the step numbers within a pass.
 *  Step:   One test patch in a pass.
 *  Sheet:  One sheet of paper, containing full and partial strips.
 *          For an XY instrument, there will be only one strip per sheet.
 *           
 */

/* TTBD:
 *
 *    Add pixmap output support (ie. TIFF)
 *
 *    Improve EPS support to add a preview to each eps file.
 */

#undef DEBUG
#undef FORCEN			/* For testing, force DeviceN */
#define DEN_COMPRESS	/* Compress density estimates > 1.0 */
						/* - this biases it towards white spacers */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <ctype.h>
#include "copyright.h"
#include "config.h"
#include "cgats.h"
#include "icc.h"
#include "xicc.h"
#include "insttypes.h"
#include "randix.h"
#include "alphix.h"
#include "rspl.h"
#include "sort.h"
#include "numlib.h"

#include <stdarg.h>
void error(char *fmt, ...), warning(char *fmt, ...), verbose(int level, char *fmt, ...);

/* Convert inches into mm */
#define inch2mm(xx) ((xx) * 25.4)

/* Convert mm into points */
#define mm2pnt(xx) ((xx) * 72.0/25.4)

/* A color structure */
struct _col {
	int nmask;		/* colorant mask */
	int pgreyt;		/* printer grey representation type 0..6 */
	int i;			/* cols list index */
	int ix;			/* random list index */
	char *id;		/* Id string */
	char loc[10];	/* Location ID string */
	int t;			/* Tag */
#define T_XYZ 0x1
#define T_DEN 0x2
#define T_RGB 0x4
#define T_N   0x8
#define T_NFB 0x2000			/* DeviceN fallback enabled */
#define T_PRESET 0x4000			/* A preset color rather than a test patch */
#define T_PAD 0x8000			/* A padding color patch */
	double XYZ[3];				/* Aproximate XYZ */
	double den[4];				/* Approx statusT density + visual density */
	double rgb[3];				/* Aproximate sRGB */
	int n;						/* Number of colorants */
	double dev[ICX_MXINKS];		/* Value of colorants */
	
	struct _col *nc[2];			/* Neigborhood colors */
	double       wnd;			/* Worst neigborhood contrast density */
}; typedef struct _col col;

#define min2(a,b) ((a) < (b) ? (a) : (b))
#define min3(a,b,c)  (min2((a), min2((b),(c))))
#define max2(a,b) ((a) > (b) ? (a) : (b))
#define max3(a,b,c)  (max2((a), max2((b),(c))))

/* Declare edge tracking functions */
void et_height(double height);
void et_media(double *rgb);
void et_color(double *rgb);
void et_edge(int isx, int negh, double mj, double mi0, double mi1);
void et_patch(char *id, double xo, double yo, double w, double h);
void et_write(char *fname, col *cols, int *rix, int si, int ei);
void et_clear(void);

/* Generate PS file prolog */
void
gen_prolog(
FILE *of,			/* Output stream */
int npages,			/* Number of pages needed */
int nmask,			/* Non zero if we are doing a DeviceN chart */
double pw, double ph,	/* Page width and height in mm */
int eps,			/* EPS flag */
int rstart			/* Random start number */
) {
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	int ipw, iph;

	ipw = ceil(mm2pnt(pw));
	iph = ceil(mm2pnt(ph));
	if (eps)
		fprintf(of,"%%!PS-Adobe-3.0 EPSF-3.0\n");
	else
		fprintf(of,"%%!PS-Adobe-3.0\n");
	fprintf(of,"%%%%Title: Argyll Color Calibration Target\n");
#ifdef FORCEN
	if (1) {
#else
	if (nmask != ICX_W		/* If not a Gray, RGB or CMYK device space */
	 && nmask != ICX_K
	 && nmask != ICX_RGB
	 && nmask != ICX_CMYK) {
#endif
		fprintf(of,"%%%%LanguageLevel: 3\n");
	} else {
		fprintf(of,"%%%%LanguageLevel: 1\n");
		if (nmask == ICX_CMYK)
			fprintf(of,"%%%%Extensions: CMYK\n");
	}

	fprintf(of,"%%%%Creator: Argyll target chart generator\n");
	fprintf(of,"%%%%For: The user who wants accurate color\n");
//	fprintf(of,"%%%%Version: %s\n",VERSION);
	fprintf(of,"%%%%CreationDate: %s",asctime(tsp));
	fprintf(of,"%%%%DocumentData: Clean7Bit\n");
	if (eps)
		fprintf(of,"%%%%Pages: %d\n",1);
	else
		fprintf(of,"%%%%Pages: %d\n",npages);
	fprintf(of,"%%%%PageOrder: Ascend\n");
	fprintf(of,"%%%%BoundingBox: %d %d %d %d\n",0,0,ipw-1,iph-1);
	fprintf(of,"%%%%Orientation: Portrait\n");		/* Rows are always virtical */
	fprintf(of,"%%%%EndComments\n");
	fprintf(of,"\n");
	fprintf(of,"<< /PageSize [%d %d] >> setpagedevice\n",ipw, iph);
	fprintf(of,"\n");

	fprintf(of,"%%%%BeginProlog\n\n");
#ifdef NEVER
	fprintf(of,"%% Duplicate nth element of stack\n");
	fprintf(of,"%% arguments: n, the offset from the element bellow the n\n");
	fprintf(of,"/dupn { 2 add dup -1 roll dup 3 -1 roll 1 roll } bind def\n");
	fprintf(of,"\n");
#endif
	fprintf(of,"%% arbitrary rectangle\n");
	fprintf(of,"%% arguments: w h x y\n");
	fprintf(of,"/rect { gsave \n");
	fprintf(of,"newpath\n");
	fprintf(of,"moveto\n");
	fprintf(of,"dup 0.0 exch rlineto\n");
	fprintf(of,"exch 0.0  rlineto\n");
	fprintf(of,"0.0 exch neg rlineto\n");
	fprintf(of,"closepath\n");
	fprintf(of,"fill\n");
	fprintf(of,"grestore } bind def\n");
	fprintf(of,"\n");
	fprintf(of,"%% hexagon with bottom left origin\n");
	fprintf(of,"%% arguments: w h x y\n");
	fprintf(of,"/hex { gsave \n");
	fprintf(of,"newpath\n");
	fprintf(of,"moveto\n");
	fprintf(of,"0 1 index rlineto\n");
	fprintf(of,"1 index 2 div 1 index 2 div rlineto\n");
	fprintf(of,"1 index 2 div 1 index 2 div neg rlineto\n");
	fprintf(of,"0 1 index neg rlineto\n");
	fprintf(of,"1 index 2 div neg 1 index 2 div neg rlineto\n");
	fprintf(of,"pop pop\n");
	fprintf(of,"closepath\n");
	fprintf(of,"fill\n");
	fprintf(of,"grestore } bind def\n");
	fprintf(of,"\n");
	fprintf(of,"%% set times-roman font\n");
	fprintf(of,"%% argument: scale\n");
	fprintf(of,"/scaleTimes {\n");
	fprintf(of,"/Times-Roman findfont\n");
	fprintf(of,"exch scalefont\n");
	fprintf(of,"setfont } bind def\n");
	fprintf(of,"\n");
	fprintf(of,"%% Print a centered string\n");
	fprintf(of,"%% argument: string, x, y\n");
	fprintf(of,"/centerShow {\n");
	fprintf(of,"gsave translate\n");
	fprintf(of,"newpath 0.0 0.0 moveto dup true charpath pathbbox\n");
	fprintf(of,"3 -1 roll sub exch 3 -1 roll sub\n");
	fprintf(of,"-0.5 mul exch -0.5 mul\n");
	fprintf(of,"moveto show grestore} bind def\n");
	fprintf(of,"\n");
	fprintf(of,"%% Print a vertically centered string\n");
	fprintf(of,"%% argument: string, x, y\n");
	fprintf(of,"/vcenterShow {\n");
	fprintf(of,"gsave translate 90.0 rotate\n");
	fprintf(of,"newpath 0.0 0.0 moveto dup true charpath pathbbox\n");
	fprintf(of,"3 -1 roll sub exch 3 -1 roll sub\n");
	fprintf(of,"-0.5 mul exch -0.5 mul\n");
	fprintf(of,"moveto show grestore} bind def\n");

	fprintf(of,"%%%%EndProlog\n");
	fprintf(of,"\n");

	if (rstart >= 0) {
		fprintf(of,"%% RandomStart %d\n",rstart);
		fprintf(of,"\n");
	}
}

/* Generate PS file epilog */
void
gen_epilog(
FILE *of)			/* Output stream */
{
	fprintf(of,"\n");
	fprintf(of,"%%%%EOF\n");
}


/* Generate PS file prolog */
void
gen_startpage(
FILE *of,			/* Output stream */
int pagen)			/* Pages number */
{
	fprintf(of,"%%%%Page: (Page %d) %d\n",pagen,pagen);
}

/* Generate PS file end of page */
void
gen_endpage(
FILE *of)			/* Output stream */
	{
	fprintf(of,"showpage\n");
	fprintf(of,"\n");
	}


/* Convert XYZ represention into XYZ density and RGB */
void
col_convert(col *cp, double *wp) {

	if ((cp->t & T_XYZ) == 0 || (cp->t & T_N) == 0)
		error("gen_color needs XYZ and device colors set !");

	if ((cp->t & T_DEN) == 0) {
		icx_XYZ2Tdens(cp->den, cp->XYZ);

#ifdef DEN_COMPRESS
		/* Compress densities > 1.0, 2:1 */
		{
			int e;
			for (e = 0; e < 3; e++) {
				if (cp->den[e] > 1.0)
					cp->den[e] = 1.0 + (cp->den[e] - 1.0) * 0.5;
			}
		}
#endif /* DEN_COMPRESS */
		cp->t |= T_DEN;
	}

	if ((cp->t & T_RGB) == 0) {
		icx_XYZ2sRGB(cp->rgb, wp, cp->XYZ);
		cp->t |= T_RGB;
	}
}

/* Set a device N color with fallback */
void
gen_ncolor(
FILE *of,			/* Output stream */
col *c) {
	int i;

	/* define the colorspace */
	fprintf(of,"[ /DeviceN [ ");
	for (i = 0; i < c->n; i++) {
		int imask = icx_index2ink(c->nmask, i);
		fprintf(of,"/%s ", icx_ink2psstring(imask));
	}

	if (c->t & T_NFB) {		/* Use color fallback */
		fprintf(of,"] /DeviceRGB ");	/* Fallback to RGB */
		fprintf(of,"{ ");
		for (i = 0; i < c->n; i++)		/* Remove N values */
			fprintf(of,"pop ");
		for (i = 0; i < 3; i++)			/* Set RGB values */
			fprintf(of,"%f ",c->rgb[i]);
	} else {
		fprintf(of,"] /DeviceGray ");	/* Fallback to Gray */
		fprintf(of,"{ ");
		for (i = 0; i < c->n; i++)		/* Remove N values */
			fprintf(of,"pop ");
		fprintf(of,"%f ",(c->rgb[0] + c->rgb[1] + c->rgb[2])/3.0); /* Set Gray value */
	}

	fprintf(of," } ] setcolorspace\n");

	/* Set the color */
	for (i = 0; i < c->n; i++)
		fprintf(of,"%f ",c->dev[i]);
	fprintf(of,"setcolor\n");
}


/* Set a device color */
/* Set it by the rep with most components */
void
gen_color(
FILE *of,			/* Output stream */
col *c) {

	if ((c->t & T_N) == 0)
		error("gen_color given color with no device values set");
	
#ifndef FORCEN
	if (c->nmask == ICX_W) {
		if ((c->t & T_PRESET) == 0)
			fprintf(of,"%% Ref %s %s %f %f %f %f\n",c->id, c->loc, 100.0 * c->dev[0]);

		if (c->pgreyt == 0) {	/* DeviceGray */
			fprintf(of,"%f setgray\n",c->dev[0]);
		} else if (c->pgreyt == 4) {	/* DeviceRGB */
			fprintf(of,"%f %f %f setrgbcolor\n",c->dev[0],c->dev[0],c->dev[0]);
		} else if (c->pgreyt == 5) {	/* Separation */
			fprintf(of,"[ /Separation (White) /DeviceGray { pop %f } ] setcolorspace\n",
			           c->dev[0]);
			fprintf(of,"%f setcolor\n",c->dev[0]);
		} else if (c->pgreyt == 6) {	/* DeviceN */
			gen_ncolor(of, c);
		} else {
			error("Device white encoding not approproate!");
		}

	} else if (c->nmask == ICX_K) {
		if ((c->t & T_PRESET) == 0)
			fprintf(of,"%% Ref %s %s %f %f %f %f\n",c->id, c->loc, 100.0 * c->dev[0]);
		if (c->pgreyt == 0) {	/* DeviceGray */
			fprintf(of,"%f setgray\n",1.0 - c->dev[0]);
		} else if (c->pgreyt == 1) {	/* DeviceCMYK */
			fprintf(of,"0.0 0.0 0.0 %f setcmykcolor\n",c->dev[0]);
		} else if (c->pgreyt == 2) {	/* Separation */
			fprintf(of,"[ /Separation (Black) /DeviceGray { pop %f } ] setcolorspace\n",
			           1.0 - c->dev[0]);
			fprintf(of,"%f setcolor\n",c->dev[0]);
		} else if (c->pgreyt == 3) {	/* DeviceN */
			gen_ncolor(of, c);
		} else {
			error("Device black encoding not approproate!");
		}

	} else if (c->nmask == ICX_RGB) {
		if ((c->t & T_PRESET) == 0)
			fprintf(of,"%% Ref %s %s %f %f %f %f\n",c->id, c->loc,
			       100.0 * c->dev[0], 100.0 *c->dev[1], 100.0 *c->dev[2]);
		fprintf(of,"%f %f %f setrgbcolor\n",c->dev[0],c->dev[1],c->dev[2]);

	} else if (c->nmask == ICX_CMYK) {
		if ((c->t & T_PRESET) == 0)
			fprintf(of,"%% Ref %s %s %f %f %f %f\n", c->id, c->loc,
			        100.0 * c->dev[0], 100.0 * c->dev[1], 100.0 * c->dev[2], 100.0 * c->dev[3]);
		fprintf(of,"%f %f %f %f setcmykcolor\n",c->dev[0],c->dev[1],c->dev[2],c->dev[3]);

	} else
#endif /* !FORCEN */
	       {	/* Device N */
		int i;
		if ((c->t & T_PRESET) == 0) {
			fprintf(of,"%% Ref %s %s",c->id, c->loc);
			for (i = 0; i < c->n; i++)
				fprintf(of,"%f ", 100.0 * c->dev[i]);
			fprintf(of,"\n");
		}
		gen_ncolor(of, c);
	}

	/* Remember edge tracking color */
	et_color(c->rgb);
}

/* Generate one testpad square */
/* Note the page coordinate origin is bottom left. */
void
gen_rect(
FILE *of,				/* Output stream */
double x, double y,		/* Top left corner of rectangle in mm from origin */
double w, double h,		/* Width and height */
char *id				/* Patch id, NULL if not a test patch */
) {
	if (w < 1e-6 || h < 1e-6)
		return;			/* Skip zero sized rectangle */
	y -= h;				/* Convert to bottom left corner */
	x = mm2pnt(x);
	y = mm2pnt(y);
	w = mm2pnt(w);
	h = mm2pnt(h);
	fprintf(of,"%f %f %f %f rect\n",w,h,x,y);

	et_patch(id, x, y, w, h);

	et_edge(1, 0, x, y, y + h);
	et_edge(1, 1, x + w, y, y + h);
	et_edge(0, 0, y, x, x + w); 
	et_edge(0, 1, y + h, x, x + w);
}

/* Generate one testpad hexagon. */
/* Note the page coordinate origin is bottom left. */
/* The hex always has left/right sides */
/* and peaks at the top and the bottom. */
void
gen_hex(
FILE *of,				/* Output stream */
double x, double y,		/* Top left vertex of hex in mm from origin */
double w, double h,		/* Width and height */
int step,				/* Step number from 0 to figure odd/even */
char *id				/* Patch id, NULL if not a test patch */
) {
	if (w < 1e-6 || h < 1e-6)
		return;			/* Skip zero sized rectangle */
	if ((step & 1) == 0) /* Even so left side of stagger */
		x -= 0.25 * w;
	else 				 /* Odd so right side of stagger */
		x += 0.25 * w;
	y = y - 5.0/6.0 * h;
	h *= 2.0/3.0;		/* Convert to hex side length */
	x = mm2pnt(x);
	y = mm2pnt(y);
	w = mm2pnt(w);
	h = mm2pnt(h);
	fprintf(of,"%f %f %f %f hex\n",w,h,x,y);
}

/* Generate a centered string */
void
gen_string(
FILE *of,				/* Output stream */
double x, double y,		/* Bot Left Corner of rectangle in mm from origin */
double w, double h,		/* Width and height */
char *s)				/* String */
{
	if (fabs(w) < 1e-6 || fabs(h) < 1e-6)
		return;			/* Skip zero sized string */
	x = mm2pnt(x);
	y = mm2pnt(y);
	w = mm2pnt(w);
	h = mm2pnt(h);
	fprintf(of,"%f scaleTimes\n",h * 0.75);
	fprintf(of,"(%s) %f %f centerShow\n",s,x+w/2.0,y+h/2.0);
}

/* Generate a vertically centered string */
void
gen_vstring(
FILE *of,				/* Output stream */
double x, double y,		/* Bot Right Corner of rectangle in mm from origin */
double w, double h,		/* Width and height */
char *s)				/* String */
{
	if (fabs(w) < 1e-6 || fabs(h) < 1e-6)
		return;			/* Skip zero sized string */
	x = mm2pnt(x);
	y = mm2pnt(y);
	w = mm2pnt(w);
	h = mm2pnt(h);
	fprintf(of,"%f scaleTimes\n",w * 0.75);
	fprintf(of,"(%s) %f %f vcenterShow\n",s,x-w/2.0,y+h/2.0);
}

void
gen_dotted_line(
FILE *of,
double x,
double y1,
double y2,
double w) {
	if (fabs(w) < 1e-6 || fabs(y2 - y1) < 1e-6)
		return;			/* Skip zero sized line */
	x = mm2pnt(x);
	y1 = mm2pnt(y1);
	y2 = mm2pnt(y2);
	w = mm2pnt(w);
	fprintf(of,"[%f %f] %f setdash\n",mm2pnt(1.0),mm2pnt(2.0),mm2pnt(0.0));
	fprintf(of,"%f setlinewidth\n",w);
	fprintf(of,"newpath %f %f moveto %f %f lineto stroke\n",x,y1,x,y2);
}

/* return the middle of 3 values */
double mid3(double a, double b, double c) {

	if ((a < b && a > c) || (a < c && a > b))
		return a;
	if ((b < a && b > c) || (b < c && b > a))
		return b;
	return c;
}

/* return the vector difference of 3 values */
double vec3(double a, double b, double c) {

	return sqrt(a * a + b * b + c * c);
}

/* Type of rule to use with color patches/spacers. */
/* This really depends on the algorithm used by the strip */
/* reading instrument. For an Xrite DTP41, it appears that */
/* using max3 is the best. This agrees with their documentation, */
/* that it is looking for the largest change in one of the three */
/* channels. Another make of instrument might use a different */
/* algorithm. */
/* Choices are: */
/* max3 - aim for largest change to be greatest */
/* mid3 - aim for middle change to be greatest */
/* min3 - aim for minimum change to be greatest */
/* vec3 - aim for the vector change to be greatest */

#define RULE3 max3

/* Setup a suitable spacer color, and */
/* Return the worst case density difference. */
double setup_spacer(
col **psc,		/* Return pointer to spacer color */
col *pp,		/* Previous patch color */
col *cp,		/* Current patch color */
col *pcol,		/* 8 pre-defined spacer colors */
int sptype		/* Spacer type code, -1 = none, */
				/* 0 = No spacer, 1 = b&W spacer, */
				/* 2 = colored     */
) {
	col *sc = NULL;	/* Spacer chosen */

//printf("~1\n");
//printf("~1 setting spacer between %s (%s) and %s (%s)\n",pp->loc, pp->id, cp->loc, cp->id);

	if (sptype <= 0) {			/* No spacer */

		/* return the density contrast between the patches */
		if (pp->nmask == ICX_W
		 || pp->nmask == ICX_K) {	/* If only capable of single density */
			double dd;

			dd = fabs(pp->den[3] - cp->den[3]);
//printf("~1 computed B&W diff of %f\n",dd);

			return dd;
		} else {
			double dd;

			dd = RULE3(fabs(pp->den[0] - cp->den[0]),
			           fabs(pp->den[1] - cp->den[1]),
			           fabs(pp->den[2] - cp->den[2]));

//printf("~1 computed color diff of %f\n",dd);
			return dd;
		}
	}

	if (pp->nmask == ICX_W
	 || pp->nmask == ICX_K) {	/* If only capable of single density */
		sptype = 1;				/* Force to B&W spacer */
	}

	if (sptype == 1) {			/* B&W spacer */
		double d1, d2;
		double gg;

		/* Choose whether space should be white or black */
		/* Shoose color that has greatest worst contrast */
		d1 = min2(fabs(pcol[0].den[3] - pp->den[3]), fabs(pcol[0].den[3] - cp->den[3]));
		d2 = min2(fabs(pcol[7].den[3] - pp->den[3]), fabs(pcol[7].den[3] - cp->den[3]));

//printf("~1 worst difference to white = %f\n", d1);
//printf("~1 worst difference to black = %f\n", d2);

		if (d1 > d2) {
//printf("~1 chosen white\n");
			sc = &pcol[0];
		} else {
//printf("~1 chosen black\n");
			sc = &pcol[7];
		}

		d1 = fabs(sc->den[3] - pp->den[3]);		/* diference between previous and spacer */
		d2 = fabs(sc->den[3] - cp->den[3]);		/* diference between current and spacer */
//printf("~1 diff between prev & space %f, curr & space %f\n",d1,d2);
//printf("~1 returning %f\n",min2(d1,d2));

		if (psc != NULL)
			*psc = sc;			/* Return pointer to chosen spacer */

		return min2(d1, d2);	/* Return worst contrast */

	} else {				/* else colored spacer */
		int ii, i;
		double p1, p2;		/* Best contrast, best other edge contrast */
		double ww;			/* Worst contrast */

		/* Check out all the possible space values for the one that gives the best */
		/* and second best contrast to each edge */

		/* for all possible spacer colors */
		p1 = p2 = ww = -1.0;
		for (i = 0; i < 8; i++) {
			double b1, b2, bb;

			b1 = RULE3(fabs(pcol[i].den[0] - pp->den[0]),
			          fabs(pcol[i].den[1] - pp->den[1]),
			          fabs(pcol[i].den[2] - pp->den[2]));

			b2 = RULE3(fabs(pcol[i].den[0] - cp->den[0]),
			          fabs(pcol[i].den[1] - cp->den[1]),
			          fabs(pcol[i].den[2] - cp->den[2]));

			bb = min2(b1, b2);	/* Worst of two edges */

			/* Worst of two edges best is better than any previous */
			if (bb > ww) {
				p1 = b1;
				p2 = b2;
				ww = bb;		/* Worst color of worst edge */
				ii = i;
				sc = &pcol[i];
			}
		}

//printf("~1 best diff between prev & space %f, curr & space %f\n",p1,p2);
//printf("~1 returning %f\n",ww);
		if (psc != NULL)
			*psc = sc;			/* Return pointer to chosen spacer */

		return ww;		/* Return worst contrast */
	}

	return 0.0;
}


/* Setup the randomised index */
void setup_randix(
int *rix,			/* Index lookup array to fill in */
int npat,			/* Number of test targets needed */
int rand,			/* Randomise flag */
int rstart,			/* Starting index for random */
int verb,			/* Verbose flag */
col *cols,			/* Array of colors to be put on target chart */
col *pcol,			/* 8 spacer colors */
int pprow,			/* Patches per row including max/min patches */
int spacer,			/* Spacer code, 0 = None, 1 = b&w, 2 = colored */
int needpc,			/* Need patch to patch contrast in a row */
int domaxmin,		/* Top and tail strip with max and min values */
col *media			/* Alias for media color */
) {
	int i, ix;
	randix *r;	/* Random index order object */

	if (rand)
		r = new_randix(npat, rstart);

	for (i = 0; i < npat; i++) {
		if (rand) {
			rix[i] = r->next(r);
		} else {
			rix[i] = i;
		}
		cols[rix[i]].ix = i;	 /* This colors random index */
	}
	rix[i] = 0;		/* npat+1 may be read */

	if (rand)
		r->del(r);

	/* Setup initial contrast check */

	if (domaxmin)
		pprow -= 2;		/* Number of test patches per row */

	{
		col *pp, *cp, *np;	/* Previous, current and next patch */
		col *maxd;	/* Alias for minimum density  */
		col *mind;	/* Alias for maximum density */
		double wc = 10.0;	/* Current worst case contrast */
		col **slist;		/* Array of pointers to colors, */
							/* sorted by worst case contrast */
		double temp, trate;	/* Annealing temperature & rate */
		double tstart, tend;/* Annealing chedule range */
		int mis = 0;		/* Number of sucessive misses */

		mind = &pcol[0];		/* White */
		maxd = &pcol[7];		/* Black */

		for (i = 0; i < npat; i++) {
			int j;		/* Row index */
			double tt;
	
			j = (i % pprow);

			if (j == 0) {			/* First in row */
				if (domaxmin)
					pp = maxd;
				else
					pp = media;
			} else {
				pp = &cols[rix[i-1]];
			}

			cp = &cols[rix[i]];

			if (j == (pprow-1) || i == (npat-1)) { /* Last in row or last patch */
				if (j == (pprow-1) && domaxmin)
					np = mind;
				else
					np = media;
			} else {
				np = &cols[rix[i+1]];
			}

			/* Setup pointers and worst case contrast */
			cp->nc[0] = pp;
			cp->nc[1] = np;
			cp->wnd = setup_spacer(NULL, pp, cp, pcol, spacer);
			tt      = setup_spacer(NULL, cp, np, pcol, spacer);
			if (tt < cp->wnd)
				cp->wnd = tt;
		}

		if ((slist = (col **)malloc(sizeof(col *) * npat)) == NULL)
			error("Malloc failed!");

		for (i = 0; i < npat; i++) {
			slist[i] = &cols[i];
		}

#define HEAP_COMPARE(A,B) (A->wnd < B->wnd)
		HEAPSORT(col *, slist, npat);
#undef HEAP_COMPARE

		if (verb)
			printf("Worst case density contrast = %f\n", slist[0]->wnd);

		if (needpc == 0 || !rand)
			return;		/* Current order is sufficient */

		if (verb) {
			printf("Optimising layout for strip reader:\n");
			printf(" 0%%"); fflush(stdout);
		}

		if (spacer == 2) {	/* Colored spacer, don't optimise too hard */
			tstart = 0.1;
			tend   = 0.002;
			trate  = 0.7;
		} else {			/* No spacer or B&W spacer, do more optimisation */
			tstart = 0.2;
			tend   = 0.00002;
			trate  = 0.9;
		}
		for (temp = tstart; temp > tend; temp *= 0.9) {

			int ii, itlim;	/* Maximum passes at a temperature */
			int nsuc = 0;	/* Number that succeed */
			int suclim;		/* Number of successful changes before continueing */

			if (spacer == 2) {	/* Colored spacer, don't optimise too hard */
				itlim = npat * 2;
				suclim = npat/2;
			} else {
				itlim = npat * 6;
				suclim = npat;
			}

			if (verb) {		/* Output percent intervals */
				double pc;
	
				pc = (log(temp) - log(tstart))/(log(tend) - log(tstart));
				printf("\r%2d%%",(int)(100.0 * pc+0.5)); fflush(stdout);
			}

			/* Improve the ordering */
			for (ii = 0; ii < itlim ; ii++) { 
				col *p1, *p2;
				double tt, de;

				/* Chose another patch to try swapping worst with */
				tt = d_rand(0.0, 1.0);
				tt *= tt;		/* Skew odds to patches with worse contrast */
				p1 = slist[0];	/* Worst */
				p2 = slist[(int)(tt * (npat - 2.0) + 1.5)];	/* Swap candidate */

				/* Check p1 in p2's place */
				de = setup_spacer(NULL, p2->nc[0], p1, pcol, spacer);
				tt = setup_spacer(NULL, p1, p2->nc[1], pcol, spacer);
				if (tt < de)
					de = tt;

				/* Check p2 in p1's place */
				tt = setup_spacer(NULL, p1->nc[0], p2, pcol, spacer);
				if (tt < de)
					de = tt;
				tt = setup_spacer(NULL, p2, p1->nc[1], pcol, spacer);
				if (tt < de)
					de = tt;

				de = de - p1->wnd;		/* Increase in worst difference */

				/* If this swap will improve things, or temp is high enough */
				if (de > 0.0
				   || d_rand(0.0, 1.0) < exp(de/temp)) {
					int t;
					col *tp0, *tp1;

					nsuc++;

					/* Swap them in random index list */
					rix[p1->ix] = p2->i; 
					rix[p2->ix] = p1->i; 
					t = p1->ix; 
					p1->ix = p2->ix; 
					p2->ix = t;

					/* Swap their neighbors, taking care */
					/* of the situation if they are neigbors */
					tp0 = p1->nc[0];
					tp1 = p2->nc[0];
					if (tp0 == p1)
						tp0 = p2;
					else if (tp0 == p2)
						tp0 = p1;
					if (tp1 == p1)
						tp1 = p2;
					else if (tp1 == p2)
						tp1 = p1;
					p2->nc[0] = tp0;
					p1->nc[0] = tp1;

					tp0 = p1->nc[1];
					tp1 = p2->nc[1];
					if (tp0 == p1)
						tp0 = p2;
					else if (tp0 == p2)
						tp0 = p1;
					if (tp1 == p1)
						tp1 = p2;
					else if (tp1 == p2)
						tp1 = p1;
					p2->nc[1] = tp0;
					p1->nc[1] = tp1;

					/* Reset backwards references */
					p1->nc[0]->nc[1] = p1;
					p1->nc[1]->nc[0] = p1;
					p2->nc[0]->nc[1] = p2;
					p2->nc[1]->nc[0] = p2;

					/* re-compute contrast to neighbors */
					p1->wnd = setup_spacer(NULL, p1->nc[0], p1, pcol, spacer);
					tt      = setup_spacer(NULL, p1, p1->nc[1], pcol, spacer);
					if (tt < p1->wnd)
						p1->wnd = tt;
					p2->wnd = setup_spacer(NULL, p2->nc[0], p2, pcol, spacer);
					tt      = setup_spacer(NULL, p2, p2->nc[1], pcol, spacer);
					if (tt < p2->wnd)
						p2->wnd = tt;

					/* Re-sort the list */
#define HEAP_COMPARE(A,B) (A->wnd < B->wnd)
					HEAPSORT(col *, slist, npat);
#undef HEAP_COMPARE

//printf("~1 worst case contrast = %1.20f\n",slist[0]->wnd);

					if (nsuc > suclim)
						break;
				}
			}
		}

		if (verb)
			printf("\r100%%\nAfter optimisation, density contrast = %f\n", slist[0]->wnd);

		free(slist);
	}

}

#define MAXPPROW 400	/* Absolute maximum patches per row permitted */

void
generate_ps(
instType itype,		/* Target instrument type */
char *bname,		/* Output file basename */
col *cols,			/* Array of colors to be put on target chart */
int npat,			/* Number of test targets needed */
char *label,		/* Per set label */
double pw,			/* Page width */
double ph,			/* Page height */
double bord, 		/* Border in mm */
int rand,			/* Randomise flag */
int rstart,			/* Starting index for random */
alphix *saix,		/* Strip alpha index object */
alphix *paix,		/* Patch alpha index object */
int ixord,			/* Index order, 0 = strip then patch */
double scale,		/* Test patch scale factor */
int hex,			/* Hexagon patches flag */
int verb,			/* Verbose flag */
int scanc,			/* Scan compatible flag, 1 = .cht gen, 2 = wide first row */
int eps,			/* EPS output flag */
int spacer,			/* Spacer code, -1 = default, 0 = None, 1 = b&w, 2 = colored */
int nmask,			/* DeviceN mask */
col *pcol,			/* 8 spacer colors */
double *wp,			/* Approximate white XYZ point */
int *ppprow,		/* Return steps per pass */
unsigned char **pprps,	/* Return malloced array holding Passes per strip */
double *p_patchlen,	/* Return patch length in mm */
double *p_gaplen,	/* Return gap length in mm */
double *p_taplen,	/* Return trailer length in mm */
int *p_npat			/* Return number of patches including padding */
) {
	char psname[200];		/* Name of output file */
	FILE *of = NULL;		/* File to write the PS to */
	double x1, y1, x2, y2;	/* Bounding box in mm */
	double iw, ih;			/* Imagable areas width and height in mm */
	double aph;				/* Available patch height */

	/* Chart definition variables */
	int domaxmin;	/* Top and tail strip with max and min values */
	int nextrap;	/* Number of extra patches for max and min */
	int needpc;		/* Need patch to patch contrast in a row */
	int dorspace;	/* Do a rrsp from center of last patch to cut line */
	int padlrow;	/* flag - Pad the last row with white */
	double lspa;	/* Leader space before first patch */
	double lcar;	/* Leading clear area before first patch */
	double plen;	/* Patch min length */
	double pspa;	/* Patch to patch spacer */
	double tspa;	/* Clear space after last patch */
	double txhi;	/* Text Height */
	int docutmarks;	/* Generate set cut marks */
	double clwi;	/* Cut line width */

	int dorowlabel;	/* Generate a set of row labels */
	double rlwi;	/* Row label test width */

	double hxew;	/* Hexagon chart extra width padding around patches */
	double hxeh;	/* Hexagon chart extra height padding around patches */

	double pwid;	/* Patch min width */
	double rrsp;	/* Row to row spacing */
	double pwex;	/* Patch width expansion between rows of a set */

	int mxpprow;	/* Maximum patches per row permitted */
	int mxrowl;		/* Maximum row length */
	int rpset;		/* Rows per set */

	double mints, minbs;	/* Minimum top & bottom space from paper edges */
	double amints, aminbs;	/* Actual mints & minbs, allowing for unused space */
	double swid;	/* Whole set width */
	int pprow;		/* patches per row (inc. max/min patches) */
	int sppage;		/* whole & partial sets per page */
	int rppset;		/* rows per partial set on whole page */
	int npages;		/* Number whole & partial pages */
	int lsppage;	/* Last page whole & partial sets per page */
	int lrpset;		/* Last sets whole & partial rows per set */
	int lpprow;		/* last row patches per row (inc. max/min patches) */
	int ppset;		/* Real patches per whole set */
	int pppage;		/* Real patches per whole page */
	int rem;		/* temporary */
	double sxwi = 0.0;	/* Scan compatible first row extra width */

	int *rix;	/* Random index lookup (Logical -> patch index) */
	int i;		/* Logical patch in target list */
	int l_si;	/* Last start of page value of i */
	int ix;		/* Patch index in target list */
	int pir;	/* Patch in row */
	int ris;	/* Row in set */
	int sip;	/* Set in page (including partial set) */
	int pif;	/* Page in file */
	double x,y;	/* Current position */
	col *pp, *cp = NULL;	/* Previous and current printed patch colors */
	int slix;	/* Strip label index */
	char *slab = NULL;	/* Strip label string */
	unsigned char *rpsp;	/* Rows per set, pointer */
	col *mark;	/* Alias for mark color */
	col *media;	/* Alias for media color */
	col *maxd;	/* Alias for minimum density  */
	col *mind;	/* Alias for maximum density */
	col *sc;	/* Alias for current spacer color */

	/* We assume that since this is intended for a printer, */
	/* the media is always white. This may not be the case */
	/* on other output media. */
	mark  = &pcol[7];		/* Black */
	media = &pcol[0];		/* White */
	mind = &pcol[0];		/* White */
	maxd = &pcol[7];		/* Black */

	/* Setup .cht edge tracking information */
	et_height(mm2pnt(ph));
	et_media(media->rgb);

	/* Set Instrument specific parameters */
	if (itype == instDTP51) { 	/* Xrite DTP51 */
		hex = 0;				/* No hex for strip instruments */
		hxew = hxeh = 0.0;		/* No extra padding because no hex */
		domaxmin = 1;			/* Print max and min patches */
		nextrap = 2;			/* Number of extra patches for max and min */
		dorspace = 1;			/* Do a rrsp from center of last patch to cut line */
		padlrow = 1;			/* Pad the last row with white */
		if (spacer < 0)
			spacer = 2;			/* Colored Spacer */
		needpc = 1;				/* Need patch to patch contrast in a row */
		lspa  = inch2mm(1.2);	/* Leader space before first patch */
		lcar  = inch2mm(0.25);	/* Leading clear area before first patch */
		plen  = scale * inch2mm(0.4);	/* Patch min length */
		if (spacer > 0)
			pspa  = inch2mm(0.07);	/* Inbetween Patch spacer */
		else
			pspa  = 0.0;		/* No spacer */
		tspa  = inch2mm(0.0);	/* Clear space after last patch */
		pwid  = inch2mm(0.4);	/* Patch min width */
		rrsp  = inch2mm(0.5);	/* Row center to row center spacing */
		pwex  = (rrsp - pwid)/2.0;	/* Patch width expansion between rows of a set */
		mxpprow = 72;			/* Maximum patches per row permitted */
		mxrowl = inch2mm(40.0);	/* Maximum row length */
		rpset = 6;				/* Rows per set */
		txhi  = 5.0;			/* Text Height */
		docutmarks = 1;			/* Generate set cut marks */
		clwi  = 0.3;			/* Cut line width */
		dorowlabel = 0;			/* Don't generate row labels */
		rlwi = 0.0;				/* Row label width */

	} else if (itype == instDTP41) {	/* Xrite DTP41 */
		hex = 0;				/* No hex for strip instruments */
		hxew = hxeh = 0.0;		/* No extra padding because no hex */
		domaxmin = 0;			/* Don't print max and min patches */
		nextrap = 0;			/* Number of extra patches for max and min */
		dorspace = 0;			/* Maximise number of rows */
		padlrow = 1;			/* Pad the last row with white */
		if (spacer < 0)
			spacer = 2;			/* Colored Spacer */
		needpc = 1;				/* Need patch to patch contrast in a row */
		lspa  = inch2mm(1.5);	/* Leader space before first patch */
		lcar  = inch2mm(0.5);	/* Leading clear area before first patch */
		plen  = scale * inch2mm(0.29);	/* Patch min length (should be 7.0 mm min.) */
		if (spacer > 0)
			pspa  = inch2mm(0.08);	/* Inbetween Patch spacer (should be 2.0 mm min.)*/
		else
			pspa  = 0.0;		/* No spacer */
		tspa  = 2 * (plen + pspa);	/* Clear space after last patch */
		pwid  = inch2mm(0.5);	/* Patch min width */
		rrsp  = inch2mm(0.5);	/* Row center to row center spacing */
		pwex  = (rrsp - pwid)/2.0;	/* Patch width expansion between rows of a set */
		mxpprow = 100;			/* Maximum patches per row permitted */
		mxrowl = inch2mm(55.0);	/* Maximum row length */
		rpset = 8;				/* Rows per set */
		txhi  = 5.0;			/* Text Height */
		docutmarks = 1;			/* Generate set cut marks */
		clwi  = 0.3;			/* Cut line width */
		dorowlabel = 0;			/* Don't generate row labels */
		rlwi = 0.0;				/* Row label width */

	} else if (itype == instSpectroScan ) {	/* GretagMacbeth SpectroScan */
		domaxmin = 0;			/* Don't print max and min patches */
		nextrap = 0;			/* Number of extra patches for max and min */
		dorspace = 0;			/* Maximise number of rows */
		padlrow = 0;			/* Pad the last row with white */
		spacer = 0;				/* No spacer */
		needpc = 0;				/* Don't need patch to patch contrast in a row */
		lspa  = 7.0 + bord;		/* Leader space before first patch - allow space for text */
		lcar  = 0.0;			/* Leading clear area before first patch */
		if (hex) {
			plen = scale * sqrt(0.75) * 7.0;	/* Patch min length */
			hxeh = 1.0/6.0 * plen;				/* Extra border for hex tops & bottoms */
			hxew = scale * 0.25 * 7.0;			/* Extra border for hex sides */
		} else {
			plen = scale * 7.0;	/* Patch min length */
			hxew = hxeh = 0.0;		/* No extra padding because no hex */
		}
		pspa  = 0.0;			/* Inbetween Patch spacer */
		tspa  = 0.0;			/* Clear space after last patch */
		pwid  = scale * 7.0;	/* Patch min width */
		rrsp  = scale * 7.0;	/* Row center to row center spacing */
		pwex  = 0.0;			/* Patch width expansion between rows of a set */
		mxpprow = MAXPPROW;		/* Maximum patches per row permitted */
		mxrowl = 1000.0;		/* Maximum row length */
		rpset = 999;			/* Rows per set */
		txhi  = 5.0;			/* Text Height */
		docutmarks = 0;			/* Don't generate set cut marks */
		clwi  = 0.0;			/* Cut line width */
		dorowlabel = 1;			/* Generate row labels */
		rlwi = 10.0;			/* Row label width */

	} else {
		error("Unsupported intrument type");
	}
	
	/* Compute page limits */
	x1 = bord;			/* Bounding box in mm */
	y1 = bord;
	x2 = pw - bord;
	y2 = ph - bord;
	iw = x2 - x1;		/* Imagable areas width and height in mm */
	ih = y2 - y1;

	*p_patchlen = plen;		/* Return patch lenth in mm */
	*p_gaplen = pspa;		/* Return gap lenth in mm */
	*p_taplen = tspa;		/* Return trailer lenth in mm */

	if (scanc & 2) 			/* Scan compatiblity */
		sxwi = pwid/2.0;	/* First row patches extra width */

	/* Compute limits for this page size */
	/* Figure the available space for patches */
	mints = bord + txhi + lcar;	/* Minimum top space due to border, text and clear area */
	if (mints < lspa)
		mints = lspa;				/* Minimum top space due to leader */
	minbs = bord;					/* Minimum botom space due to border */
	if (minbs < tspa)
		minbs = tspa;				/* Minimum botom space due to trailer */
	aph = ph - mints - minbs - 2.0 * hxeh;	/* Available space for printing test patches */
	if (aph > mxrowl)
		aph = mxrowl;				/* Limit maximum row length */

	/* We are assuming that every patch may be surounded by a spacer */
	/* (ie. there are always pprow+1 gaps spacers are used) */
	pprow = (aph - pspa)/(plen + pspa);	/* Raw Patches per row */
	if (pprow > mxpprow)			/* Limit to maximum */
		pprow = mxpprow;
	
	if (pprow < (1+nextrap))
		error("Paper size not long enought for a single patch per row!");

	*ppprow = pprow - nextrap;	/* Real patches per row */

	if ((*pprps = (unsigned char *)malloc(sizeof(unsigned char) * (2 + (npat/pprow)))) == NULL)
		error("Malloc failed!");
	rpsp = *pprps;
	
	/* Compute actual lowest coordinate used */
	aminbs = ph - mints - pspa - pprow * (plen + pspa);
	amints = mints + 0.5 * (aminbs - minbs); 	/* Distribute extra space */
	aminbs = minbs + 0.5 * (aminbs - minbs);

	/* Compute whole set width */
	if (dorspace)
		swid = rpset * rrsp + pwid/2.0;			/* set gutter is rrsp - pwid/2 wide */
	else
		swid = (rpset-1) * rrsp + pwid + clwi;	/* set gutter is 0, but allow for cut line */
		
	/* Compute sets per page */
	sppage = (iw - rlwi - sxwi - 2.0 * hxew)/swid + 1; /* Number of whole sets + partial sets */

	/* Compute rows per partial set on whole page */
	if (dorspace)
		rppset = (iw - rlwi - sxwi - 2.0 * hxew - swid * (sppage-1) - pwid/2.0)/rrsp;
	else
		rppset = (iw - rlwi - sxwi - 2.0 * hxew - swid * (sppage-1) - pwid + rrsp)/rrsp;
	if (rppset < 0)
		rppset = 0;
	if (rppset == 0) {	/* Make last partial set a full set */
		sppage--;
		rppset = rpset;
	}
	
	if (sppage <= 0)
		error("Not enough width for even one row!");

	/* The number of pages needed */
	pppage = (pprow-nextrap) * ((sppage-1) * rpset + rppset);	/* Real patches per page */
	npages = (npat + pppage -1)/pppage;						/* whole & partial pages */
	ppset = (pprow-nextrap) * rpset;						/* Real patches per full set */

	rem = npat;							/* Total test patches to print */
	rem -= (npages-1) * pppage;			/* Remaining patches to be printed on last page */

	lsppage = (rem + ppset -1)/ppset;	/* Last pages whole & partial sets per page */
	rem -= (lsppage - 1) * ppset;		/* remaining patches to be printed in last set */

	lrpset = (rem + (pprow-nextrap) - 1)/(pprow-nextrap);
										/* Last sets whole & partial rows per set */

	rem -= (lrpset - 1) * (pprow-nextrap); /* remaining patches to be printed in last row */

	lpprow = rem + nextrap;				/* Patches in last row of last set of last page */

	if (verb) {
		fprintf(stderr,"Patches = %d\n",npat);
		fprintf(stderr,"Test patches per row = %d\n",(pprow-nextrap));
		fprintf(stderr,"Sets per page = %d, rows per partial set = %d\n",sppage, rppset);
		fprintf(stderr,"Rows in last set = %d, patches in last row = %d\n", lrpset, lpprow-nextrap);
		fprintf(stderr,"Total pages needed = %d\n",npages);
	}

	if (padlrow) {		/* Add in extra padding patches */
		int i;
//printf("~1 adding %d padding patches, going from %d to %d patches\n",
//pprow - lpprow, npat, npat + pprow - lpprow);
		for (i = 0; lpprow < pprow; lpprow++, npat++, i = (i + 1) & 7) {
			if (needpc && rand)
				cols[npat] = pcol[i];		/* structure copy */
			else
				cols[npat] = *media;		/* structure copy */
			cols[npat].i = npat;			/* Now has an index */
			cols[npat].t &= ~T_PRESET;
			cols[npat].t |= T_PAD; 
			cols[npat].id = "0";			/* Padding identification */
		}
	}
	*p_npat = npat;			/* Return number of patches including padding */
//printf("~1 padded no of patches = %d\n", npat);

	/* Setup logical to random patch mapping */
	if ((rix = (int *)malloc(sizeof(int) * (npat + 1))) == NULL)
		error("Malloc failed!");

	setup_randix(rix, npat, rand, rstart, verb, cols, pcol,
	             pprow, spacer, needpc, domaxmin, media);
	rix[npat] = -1;		/* Shouldn't use this */

	/* Init everything */
	l_si = i = 0;	/* Physical test patch index */
	ix = rix[i];	/* First index */

	pir = 1;		/* Starting patch in row (includes max/min patches) */
	ris = 1;		/* Starting row in set */
	sip = 1;		/* Starting set in page */
	pif = 1;		/* Starting page in file */

	slix = -1;		/* Start at one befor first strip index */

	/* Until there are no more patches to do */
	for (;;) {
		char *sp = NULL;	/* String pointer - label */
		double w;			/* Width of next patch */
		int flags;			/* flags for current patch */
#define IS_FPIR   0x00001	/* Is first patch in row (possibly max density patch) */
#define IS_LPIR   0x00002	/* Is last patch in row (possibly min density patch) */
#define IS_SLPIR  0x00004	/* Is second last patch in row */
#define IS_FRIS   0x00008	/* Is first row in set */
#define IS_LRIS   0x00010	/* Is last row in set */
#define IS_FSIP   0x00020	/* Is first set in page */
#define IS_LSIP   0x00040	/* Is last set in page */
#define IS_FPIF   0x00080	/* Is first page in file */
#define IS_LPIF   0x00100	/* Is last page in file */
#define IS_PAD    0x04000	/* Is fake padding patch, to round out very last row */
#define IS_SPAD   0x08000	/* Is second to fake padding patch */

		/* Init flags for this patch */
		flags = 0;
		if (pir == 1)
			flags |= IS_FPIR;
		if (pir == pprow-1)
			flags |= IS_SLPIR;
		if (pir == pprow)
			flags |= IS_LPIR;

		if (ris == 1)
			flags |= IS_FRIS;
		if (ris == rpset)
			flags |= IS_LRIS;

		if (sip == 1)
			flags |= IS_FSIP;
		if (sip == sppage) {
			flags |= IS_LSIP;
			if (ris == rppset)
				flags |= IS_LRIS;
		}
		if (pif == 1)
			flags |= IS_FPIF;
		if (pif == npages) {
			flags |= IS_LPIF;
			if (sip == lsppage) {
				flags |= IS_LSIP;
				if (ris == lrpset) {
					flags |= IS_LRIS;

					if (padlrow) {				/* Last row in chart may be a runt */
						if (pir >= lpprow)
							flags |= IS_SPAD;
						if (pir > lpprow)
							flags |= IS_PAD;
					} else {
						if (pir == (lpprow-1))
							flags |= IS_SLPIR;
						if (pir == lpprow)
							flags |= IS_LPIR;
					}
				}
			}
		}

//printf("~1 pir %d, ris %d, sip %d, pif %d\n", pir, ris, sip, pif);

		/* Set initial patch width */
		w = pwid;
		if (!(flags & IS_FRIS))
			w += pwex;		/* Extend patch into previous row to row gap */
		if (!(flags & IS_LRIS))
			w += pwex;		/* Extend patch into next row to row gap */

		if ((flags & IS_FRIS) && (flags & IS_FSIP))	/* First row on page */
			w += sxwi;		/* Make first row fatter for scan compatiblity */

		if (flags & IS_FPIR) {			/* Start of row */
			y = ph - amints;			/* Start ready for spacer or patch */ 

			if (flags & IS_FRIS) {							/* Start of set */
				if (flags & IS_FSIP) {						/* Start of page */
					x = x1;									/* Start at leftmost position */
					if (eps) {	/* EPS */
						if (npages > 1)
							sprintf(psname,"%s%d.eps",bname,pif);
						else
							sprintf(psname,"%s.eps",bname);
						if ((of = fopen(psname,"w")) == NULL)
							error ("Unable to open output file '%s'",psname);
						gen_prolog(of,npages,nmask,pw,ph,eps,rstart);
					} else {	/* PS */

						if (flags & IS_FPIF) {					/* First page */
							sprintf(psname,"%s.ps",bname);
							if ((of = fopen(psname,"w")) == NULL)
								error ("Unable to open output file '%s'",psname);
							gen_prolog(of,npages,nmask,pw,ph,eps,rstart);
						}
					}
					gen_startpage(of,pif);

					/* Print all the row labels */
					/* (we are assuming no min/max density patches !!) */
					if (dorowlabel) {
						double ty = y;		/* Temporary y coord */
						int tpir;

						for (tpir = 0; tpir < pprow; tpir++) {
							char *rlabl;
							if ((rlabl = paix->aix(paix, tpir)) == NULL)
								error ("Patch in row label %d out of range",tpir);
							gen_color(of, mark);
							gen_string(of, x, ty-plen, rlwi, plen, rlabl);
							free(rlabl);
						
							ty -= plen + pspa;
						}

						x += rlwi;
					}

					x += hxew;		/* Allow space for extra bits of hexagons */

					/* Clear edge list tracking */
					et_clear();
				}
			}

			/* Increment strip label */
			slix++;

			if (slab != NULL)
				free(slab);
			if ((slab = saix->aix(saix, slix)) == NULL)
				error("strip index %d out of range",slix);

			if ((lspa - lcar - bord) >= txhi) {	/* There is room for label */
				/* Print strip label */
				gen_color(of, mark);
				gen_string(of,x,y2-txhi,w,txhi,slab);
			}

			/* Start with background = media */
			pp = media;

		}

		/* Figure the current patch color */
		if ((flags & IS_FPIR) && domaxmin) {		/* Max at start of row */
			sp = NULL;								/* Not a test patch (no label) */
			cp = maxd;								/* Maximum density patch at start */

		} else if ((flags & IS_LPIR) && domaxmin) { /* Min at end of rows */
			sp = NULL;								/* Not a test patch (no label) */
			cp = mind;

		} else if (flags & IS_PAD) {				/* Fake blank patch */
			sp = NULL;								/* Not a test patch (no label) */
			cp = media;

		} else {	/* set test patch location and color */
			int apir = domaxmin ? pir - 1 : pir;	/* Adjusted pir for max/min patches */
			if (sp != NULL)
				free(sp);
			if ((sp = patch_location(saix, paix, ixord, slix, apir-1)) == NULL)
				error ("Patch location out of range, strip %d patch %d",slix,apir-1);
			strcpy(cols[ix].loc, sp);		/* Record location */
			cp = &cols[ix];					/* Get color for this patch */

			i++;							/* Consumed a test patch */
			if (i > npat)
				error("Internal - ran out of test patches !");
				 
			ix = rix[i];					/* Next patch index */
		}

		/* Print a spacer in front of patch if requested */
		if (spacer > 0) {
			setup_spacer(&sc, pp, cp, pcol, spacer);
			gen_color(of, sc);
			gen_rect(of, x, y, w, pspa, NULL);
			y -= pspa;
		}

		/* Print patch */
		gen_color(of, cp);
		if (hex) {
			int apir = domaxmin ? pir - 1 : pir;	/* Adjusted pir for max/min patches */
			gen_hex(of, x, y, w, plen, apir-1, sp);
		} else
			gen_rect(of, x, y, w, plen, sp);
		y -= plen;
		if (sp != NULL) {		/* Done with sp for the moment */
			free(sp);
			sp = NULL;
		}

		/* Advance the patch count */
		pir++;
		pp = cp;	/* Current color becomes last color */

		/* If this is the last patch in the row, */
		/* print a possible last spacer. */
		if (flags & IS_LPIR) {								/* End of row */
			cp = media;
			if (spacer > 0) {
				setup_spacer(&sc, pp, cp, pcol, spacer);
				gen_color(of, sc);
				gen_rect(of, x, y, w, pspa, NULL);
				y -= pspa;
			}
		}

		if (flags & IS_LPIR) {								/* Last patch in row */
			pir = 1;
			ris++;

			/* First step to the middle of the patch */
			if (flags & IS_FRIS)
				x += pwid/2.0;
			else
				x += pwid/2.0 + pwex;

			/* Then step to the start of the next patch */
			if (flags & IS_LRIS) {
				if (dorspace)
					x += rrsp;				/* row to row space, making gutter */
				else
					x += pwid/2.0 + clwi;	/* no gutter, but room for cut line */
			} else
				x += pwid/2.0 + pwex;

			if ((flags & IS_FRIS) && (flags & IS_FSIP))	/* First row on page */
				x += sxwi;					/* Allow for scan compatible fatter first row */

			if (flags & IS_LRIS) {							/* End of set */
				*rpsp++ = ris-1;	/* Record how many rows in this set */
				ris = 1;
				sip++;
				
				/* Print end of set crop line */
				gen_color(of, mark);
				if (docutmarks)			/* Generate set cut marks */
					gen_dotted_line(of,x-0.3/2.0,y1,y2,0.3);	/* 0.3 wide dotted line */

				/* Print end of set identification if there is space */
				if (dorspace)
					gen_vstring(of,x,y1,rrsp-pwid/2.0-pwex,y2-y1,label);

				if (flags & IS_LSIP) {						/* End of page */
					sip = 1;

					gen_endpage(of);
					if (eps) {
						gen_epilog(of);
						if (fclose(of))
							error ("Unable to close output file '%s'",psname);
					}
					if (flags & IS_LPIF) {					/* Last page in file */
						if (!eps) {
							gen_epilog(of);
							if (fclose(of))
								error ("Unable to close output file '%s'",psname);
						}
					}

					/* If we are anticipating scanner input, create the */
					/* scanner recognition file for this page. */
					if (scanc & 1) {
						char chtname[200];	/* Name of .cht file */

						if (npages > 1)
							sprintf(chtname,"%s%d.cht",bname,pif);
						else
							sprintf(chtname,"%s.cht",bname);
						et_write(chtname, cols, rix, l_si, i);
					}
					l_si = i;		/* New last start i */

					if (flags & IS_LPIF) {					/* Last page in file */
						break;								/* Done */
					}
					pif++;
				}
			}
		}
	}
	if (slab != NULL)
		free(slab);
	free(rix);

	*rpsp++ = 0;	/* End of rows per set stuff */

	et_clear();		/* Cleanup edge list structures */
}

/* A paper size structure */
typedef struct {
	char *name;			/* User name (lower case) */
	double w,h;			/* Width and height in mm */
	int def;			/* Non-zero if default */
} paper;

static paper psizes[] = {
	{ "A4", 	 210.0,	297.0, 0 },
	{ "A4R", 	 297.0,	210.0, 0 },
	{ "A3", 	 297.0,	420.0, 1 },
	{ "A2", 	 420.0,	594.0, 0 },
	{ "Letter",	 215.9,	279.4, 0 },
	{ "LetterR", 279.4,	215.9, 0 },
	{ "Legal",	 215.9,	355.6, 0 },
	{ "11x17",	 279.4,	431.8, 0 },
	{ NULL, 0.0, 0.0, 0 }
};

/* Case independent string compare */
int
cistrcmp(char *s1, char *s2) {
	for (;;s1++, s2++) {
		if (tolower(*s1) != tolower(*s2))
			return 1;
		if (*s1 == '\000')
			return 0;
	}
}

#define DEF_SIXPAT	"A-Z, A-Z"				/* Default strip index pattern */		
#define DEF_PIXPAT	"0-9,@-9,@-9;1-999"		/* Default patch index pattern */		

void usage(char *diag, ...) {
	paper *pp;
	fprintf(stderr,"Generate Target PostScrip file, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	if (diag != NULL) {
		va_list args;
		fprintf(stderr,"Diagnostic: ");
		va_start(args, diag);
		vfprintf(stderr, diag, args);
		va_end(args);
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"usage: printtarg [-v] [-i instr] [-r] [-s] [-p size] basename\n");
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -i 41 | 51 | SS Select target instrument (default DTP41)\n");
	fprintf(stderr,"                 41 = DTP41, 51 = DTP51, SS = SpectroScan\n");
	fprintf(stderr," -h              Use hexagon patches for SS\n");
	fprintf(stderr," -a scale        Scale patch size by factor (e.g. 0.857 or 1.5 etc.)\n");
	fprintf(stderr," -r              Don't randomize patch location\n");
	fprintf(stderr," -s              Create a scan image recognition (.cht) file\n");
	fprintf(stderr," -S              Same as -s, but don't generate wide orientation strip.\n");
	fprintf(stderr," -c              Force colored spacers\n");
	fprintf(stderr," -b              Force B&W spacers\n");
	fprintf(stderr," -n              Force no spacers\n");
	fprintf(stderr," -f              Create DeviceN Color fallback\n");
	fprintf(stderr," -w g|r|s|n      White test encoding DeviceGray (def), DeviceRGB, Separation or DeviceN\n");            
	fprintf(stderr," -k g|c|s|n      Black test encoding DeviceGray (def), DeviceCMYK, Separation or DeviceN\n");            
	fprintf(stderr," -e              Output EPS compatible file\n");
	fprintf(stderr," -t rsnum        Use given random start number\n");
	fprintf(stderr," -x pattern      Use given strip indexing pattern (Default = \"%s\")\n",DEF_SIXPAT);
	fprintf(stderr," -y pattern      Use given patch indexing pattern (Default = \"%s\")\n",DEF_PIXPAT);
	fprintf(stderr," -p size         Select page size from:\n");
	for (pp = psizes; pp->name != NULL; pp++)
		fprintf(stderr,"                 %s	[%.1f x %.1f mm]%s\n", pp->name, pp->w, pp->h,
		pp->def ? " (default)" : "");
	fprintf(stderr," -p WWWxHHH      Custom size, WWW mm wide by HHH mm high\n");
	fprintf(stderr," basname         Base name for input(.ti1)/output(.ti2)\n");
	exit(1);
	}

int
main(argc,argv)
int argc;
char *argv[];
{
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	int hex = 0;			/* Hexagon patches */
	double scale = 1.0;		/* Patch size scale */
	int rand = 1;
	int eps = 0;
	int spacer = -1;		/* -1 = default for instrument */
							/* 0 = forse no spacer, 1 = Force B&W spacers */
							/* 2 = Force colored spacer */
	int rstart = -1;		/* Random sequence start value */
	char *sixpat = DEF_SIXPAT;	/* Strip index pattern */		
	char *pixpat = DEF_PIXPAT;	/* Patch index pattern */		
	alphix *saix, *paix;	/* Strip and Patch index generators */
	int ixord = 0;			/* Index order, 0 = strip then patch */
	int scanc = 0;			/* Scan compatible flag, 1 = .cht, 2 = wide first row */
	int devnfb = 0;			/* Add device N fallback colors */
	int pgreyt = 0;			/* Device K/W color type 0..6 */
	static char inname[200] = { 0 };		/* Input cgats file base name */
	static char outname[200] = { 0 };		/* Output cgats file base name */
	static char psname[200] = { 0 };		/* Output postscrip file base name */
	static char sname[200] = { 0 };			/* Output scanner .cht file name */
	cgats *icg;				/* input cgats structure */
	cgats *ocg;				/* output cgats structure */
	instType itype = instDTP41;		/* Default target instrument */
	int nmask = 0;			/* Device colorant mask */
	int nchan = 0;			/* Number of device chanels */
	int i;
	int si, fi, wi;			/* sample id index, field index, keyWord index */
	unsigned long ss = 0x65de4523;
	char *label = "Argyll Color Management System - Test chart";
	double marg = 6.0;		/* Margin from paper edge in mm */
	paper *pap = NULL;		/* Paper size pointer, NULL if custom */
	double cwidth, cheight;	/* Custom paper width and height in mm */
	col *cols;				/* test patch colors */
	int npat;				/* Number of patches */
	int nppat;				/* Number of patches including padding */
	col pcol[8];			/* pre-defined colors */
	double wp[3];			/* Approximate XYZ white point */
	time_t clk = time(0);
	struct tm *tsp = localtime(&clk);
	char *atm = asctime(tsp); /* Ascii time */
	int sip;				/* Steps in Pass */
	unsigned char *pis;		/* Passes in strip array */
	double plen, glen, tlen;/* Patch, gap and trailer length in mm */
	char buf[100];			/* general sprintf buffer */

#if defined(__IBMC__)
	_control87(EM_UNDERFLOW, EM_UNDERFLOW);
	_control87(EM_OVERFLOW, EM_OVERFLOW);
#endif

	if (argc <= 1)
		usage("Not enough arguments");

#ifdef DEBUG
	printf("Debug is set\n");
#endif

	/* Find the default paper size */
	for (pap = psizes; pap->name != NULL; pap++) {
		if (pap->def != 0)
			break;
	}
	if (pap->name == NULL)
		error ("Internal - can't find default paper size");

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')	{	/* Look for any flags */
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
				usage("Requested usage");

			/* Verbosity */
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V')
				verb = 1;

			/* hexagon patches */
			else if (argv[fa][1] == 'h' || argv[fa][1] == 'H')
				hex = 1;

			/* Patch scaling */
			else if (argv[fa][1] == 'a' || argv[fa][1] == 'A') {
				fa = nfa;
				if (na == NULL) usage("Expected scale factor to -a");
				scale = atof(na);
				if (scale < 0.1 || scale > 4.0)
					usage("Scale factor %f is outside expected range 0.1 - 4.0",scale);
			}

			/* Scan compatible */
			else if (argv[fa][1] == 's')
				scanc = 3;

			else if (argv[fa][1] == 'S')
				scanc = 1;

			/* Force colored spacer */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C')
				spacer = 2;

			/* Force B&W spacer */
			else if (argv[fa][1] == 'b' || argv[fa][1] == 'B')
				spacer = 1;

			/* No spacer */
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N')
				spacer = 0;

			/* Randomisation */
			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R')
				rand = 0;

			/* Enable DeviceN color fallback */
			else if (argv[fa][1] == 'f' || argv[fa][1] == 'F')
				devnfb = 1;

			/* Select the printer W color representation */
			else if (argv[fa][1] == 'w' || argv[fa][1] == 'W') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -w");
				switch(na[0]) {
					case 'g':
					case 'G':
						pgreyt = 0;
						break;
					case 'r':
					case 'R':
						pgreyt = 4;
						break;
					case 's':
					case 'S':
						pgreyt = 5;
						break;
					case 'n':
					case 'N':
						pgreyt = 6;
						break;
					default:
						usage("Unexpected argument to -w");
				}
			}

			/* Select the printer K color representation */
			else if (argv[fa][1] == 'k' || argv[fa][1] == 'K') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -k");
				switch(na[0]) {
					case 'g':
					case 'G':
						pgreyt = 0;
						break;
					case 'c':
					case 'C':
						pgreyt = 1;
						break;
					case 's':
					case 'S':
						pgreyt = 2;
						break;
					case 'n':
					case 'N':
						pgreyt = 3;
						break;
					default:
						usage("Unexpected argument to -k");
				}
			}

			/* EPS */
			else if (argv[fa][1] == 'e' || argv[fa][1] == 'E')
				eps = 1;

			/* Specify random seed */
			else if (argv[fa][1] == 't' || argv[fa][1] == 'T') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -t");
				rstart = atoi(na);
				if (rstart < 0)
					usage("Argument to -t must be positive");
			}

			/* Specify strip index pattern */
			else if (argv[fa][1] == 'x' || argv[fa][1] == 'X') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -x");
				sixpat = na;
			}

			/* Specify patch index pattern */
			else if (argv[fa][1] == 'y' || argv[fa][1] == 'Y') {
				fa = nfa;
				if (na == NULL) usage("Expected argument to -y");
				pixpat = na;
			}

			/* Page size */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				char *p;
				int tt;
				fa = nfa;
				if (na == NULL) usage("Expected an argument to -p");
				for (pap = psizes; pap->name != NULL; pap++) {
					if (cistrcmp(na, pap->name) == 0)
						break;
				}
				
				if (pap->name == NULL) {	/* See if it matches a custom size */
					if (sscanf(na,"%lfx%lf",&cwidth, &cheight) == 2) {
						pap = NULL;			/* Indicate custom */
						if (cwidth < 1.0 || cwidth > 4000.0
						 || cheight < 1.0 || cheight > 4000.0)
							usage("Argument to -p was of unexpected size");		/* Sanity check */
					} else {
						usage("Failed to recognise argument to -p");
					}
				}
			}
			/* Target Instrument type */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				char *p;
				int tt;
				fa = nfa;
				if (na == NULL) usage("Expected an argument to -i");

				if (strcmp("51", na) == 0)
					itype = instDTP51;
				else if (strcmp("41", na) == 0)
					itype = instDTP41;
				else if (strcmp("SS", na) == 0 || strcmp("ss", na) == 0)
					itype = instSpectroScan;
				else
					usage("Argument to -i wasn't recognised");
			} else 
				usage("Unknown flag");
		}
		else
			break;
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage("Expecting basename argument");
	strcpy(inname,argv[fa]);
	strcat(inname,".ti1");
	strcpy(outname,argv[fa]);
	strcat(outname,".ti2");
	strcpy(psname,argv[fa]);
	strcpy(sname,argv[fa]);
	strcat(sname,".cht");

	if (hex && itype != instSpectroScan) {
		if (verb)
			printf("Can only select hexagonal patches for SpectrScan - ignored!\n");
		hex = 0;
	}

	if (hex && scanc) {
		if (verb)
			printf("Can only select hexagonal patches if no scan recognition is needed - ignored!\n");
		hex = 0;
	}

	if ((saix = new_alphix(sixpat)) == NULL)
		error("Strip indexing pattern '%s' doesn't parse",sixpat);

	if ((paix = new_alphix(pixpat)) == NULL)
		error("Patch in strip indexing pattern '%s' doesn't parse",pixpat);

	icg = new_cgats();	/* Create a CGATS structure */
	icg->add_other(icg, "CTI1"); 	/* our special input type is Calibration Target Information 1 */

	if (icg->read_name(icg, inname))
		error("CGATS file read error : %s",icg->err);

	if (icg->t[0].tt != tt_other || icg->t[0].oi != 0)
		error ("Input file isn't a CTI1 format file");
	if (icg->ntables != 2)
		error ("Input file doesn't contain exactly two tables");

	if ((npat = icg->t[0].nsets) <= 0)
		error ("No sets of data");

	/* Allocate room for test patches and maximum padding patches */
	if ((cols = (col *)malloc(sizeof(col) * (npat + MAXPPROW))) == NULL)
		error("Malloc failed!");

	/* Setup output cgats file */
	ocg = new_cgats();	/* Create a CGATS structure */
	ocg->add_other(ocg, "CTI2"); 	/* our special type is Calibration Target Information 2 */
	ocg->add_table(ocg, tt_other, 0);	/* Start the first table */

	ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 2",NULL);
	ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll printtarg", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	ocg->add_kword(ocg, 0, "CREATED",atm, NULL);

	/* Note what instrument the chart is setup for */
	ocg->add_kword(ocg, 0, "TARGET_INSTRUMENT", inst_name(itype) , NULL);

	/* Copy various parameters through */
	if ((wi = icg->find_kword(icg, 0, "SINGLE_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "SINGLE_DIM_STEPS",icg->t[0].kdata[wi], NULL);
	
	if ((wi = icg->find_kword(icg, 0, "COMP_GREY_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "COMP_GREY_STEPS",icg->t[0].kdata[wi], NULL);
	
	if ((wi = icg->find_kword(icg, 0, "MULTI_DIM_STEPS")) >= 0)
		ocg->add_kword(ocg, 0, "MULTI_DIM_STEPS",icg->t[0].kdata[wi], NULL);

	if ((wi = icg->find_kword(icg, 0, "FULL_SPREAD_PATCHES")) >= 0)
		ocg->add_kword(ocg, 0, "FULL_SPREAD_PATCHES",icg->t[0].kdata[wi], NULL);

	/* Fields we want */
	ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);
	ocg->add_field(ocg, 0, "SAMPLE_LOC", cs_t);

	if ((si = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
		error ("Input file doesn't contain field SAMPLE_ID");
	if (icg->t[0].ftype[si] != nqcs_t)
		error ("Field SAMPLE_ID is wrong type");

	/* Read the approximate white point */
	if ((fi = icg->find_kword(icg, 0, "APPROX_WHITE_POINT")) < 0)
		error ("Input file doesn't contain keyword APPROX_WHITE_POINT");
	if (sscanf(icg->t[0].kdata[fi], "%lf %lf %lf", &wp[0], &wp[1], &wp[2]) != 3)
		error ("Couldn't parse the white point data correctly");
	wp[0] /= 100.0; wp[1] /= 100.0; wp[2] /= 100.0;
	ocg->add_kword(ocg, 0, "APPROX_WHITE_POINT",icg->t[0].kdata[fi], NULL);

//printf("~1 got approx white point of %f %f %f\n",wp[0],wp[1],wp[2]);

	/* Figure out the color space */
	if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
		error ("Input file doesn't contain keyword COLOR_REPS");

	if ((nmask = icx_char2inkmask(icg->t[0].kdata[fi])) != 0) {
		int i, j, ii;
		int chix[ICX_MXINKS];	/* Device chanel indexes */
		int xyzix[3];			/* XYZ chanel indexes */
		char *ident;
		char *xyzfname[3] = { "XYZ_X", "XYZ_Y", "XYZ_Z" };

		if ((ii = icg->find_kword(icg, 0, "TOTAL_INK_LIMIT")) >= 0)
			ocg->add_kword(ocg, 0, "TOTAL_INK_LIMIT",icg->t[0].kdata[ii], NULL);

		nchan = icx_noofinks(nmask);
		ident = icx_inkmask2char(nmask); 

		for (j = 0; j < nchan; j++) {
			int imask;
			char fname[100];

			imask = icx_index2ink(nmask, j);
			sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : ident,
			                      icx_ink2char(imask));

			if ((ii = icg->find_field(icg, 0, fname)) < 0)
				error ("Input file doesn't contain field %s",fname);
			if (icg->t[0].ftype[ii] != r_t)
				error ("Field %s is wrong type",fname);
	
			ocg->add_field(ocg, 0, fname, r_t);
			chix[j] = ii;
		}

		for (j = 0; j < 3; j++) {
			if ((ii = icg->find_field(icg, 0, xyzfname[j])) < 0)
				error ("Input file doesn't contain field %s",xyzfname[j]);
			if (icg->t[0].ftype[ii] != r_t)
				error ("Field %s is wrong type",xyzfname[j]);
	
			ocg->add_field(ocg, 0, xyzfname[j], r_t);
			xyzix[j] = ii;
		}

		ocg->add_kword(ocg, 0, "COLOR_REP", ident, NULL);

		/* Read all the test patches in */
		for (i = 0; i < npat; i++) {
			cols[i].i = i;
			cols[i].t = T_N | T_XYZ;
			if (devnfb)
				cols[i].t |= T_NFB;
			cols[i].nmask = nmask;
			cols[i].pgreyt = pgreyt;
			cols[i].n  = nchan;
			cols[i].id = ((char *)icg->t[0].fdata[i][si]);
			sprintf(cols[i].loc, "???");
			for (j = 0; j < nchan; j++)
				cols[i].dev[j] = *((double *)icg->t[0].fdata[i][chix[j]]) / 100.0;
			for (j = 0; j < 3; j++)
				cols[i].XYZ[j] = *((double *)icg->t[0].fdata[i][xyzix[j]]) / 100.0;
			col_convert(&cols[i], wp);	/* Ensure other representations */
		}

		free(ident);
	} else
		error ("Input file keyword COLOR_REPS has unknown value");

	/* Load up the pre-defined spacer colors */
	{
		int i, j, ii;
		int nsp;
		int chix[ICX_MXINKS];	/* Device chanel indexes */
		int xyzix[3];			/* XYZ chanel indexes */
		char *ident;
		char *xyzfname[3] = { "XYZ_X", "XYZ_Y", "XYZ_Z" };

		nchan = icx_noofinks(nmask);
		ident = icx_inkmask2char(nmask); 

		if ((nsp = icg->t[1].nsets) <= 0)
			error ("No sets of data in second table");

		for (j = 0; j < nchan; j++) {
			int imask;
			char fname[100];

			imask = icx_index2ink(nmask, j);
			sprintf(fname,"%s_%s",nmask == ICX_W || nmask == ICX_K ? "GRAY" : ident,
			                      icx_ink2char(imask));

			if ((ii = icg->find_field(icg, 1, fname)) < 0)
				error ("Input file doesn't contain field %s",fname);
			if (icg->t[1].ftype[ii] != r_t)
				error ("Field %s is wrong type",fname);
			chix[j] = ii;
		}

		for (j = 0; j < 3; j++) {
			if ((ii = icg->find_field(icg, 1, xyzfname[j])) < 0)
				error ("Input file doesn't contain field %s",xyzfname[j]);
			if (icg->t[1].ftype[ii] != r_t)
				error ("Field %s is wrong type",xyzfname[j]);
			xyzix[j] = ii;
		}

		if (nsp != 8)
			error ("Expect second set of data to have 8 sets, found %d",nsp);

		/* Read all the spacer patches in */
		for (i = 0; i < nsp; i++) {
			pcol[i].i = -1;
			pcol[i].t = T_N | T_XYZ | T_PRESET;
			if (devnfb)
				pcol[i].t |= T_NFB;
			pcol[i].nmask = nmask;
			pcol[i].pgreyt = pgreyt;
			pcol[i].n  = nchan;
			pcol[i].id = "";
			sprintf(cols[i].loc, "???");
			for (j = 0; j < nchan; j++)
				pcol[i].dev[j] = *((double *)icg->t[1].fdata[i][chix[j]]) / 100.0;
			for (j = 0; j < 3; j++)
				pcol[i].XYZ[j] = *((double *)icg->t[1].fdata[i][xyzix[j]]) / 100.0;
			col_convert(&pcol[i], wp);	/* Ensure other representations */
		}

		free(ident);
	}

	if (verb) {
		if (pap != NULL)
			printf("Paper chosen is %s	[%.1f x %.1f mm]\n", pap->name, pap->w, pap->h);
		else
			printf("Paper chosen is custom %.1f x %.1f mm\n", cwidth, cheight);
	}

	if (rand) {
		if (rstart == -1) {
			rstart = clk % npat;
		} else {
			rstart = rstart % npat;
		}
		sprintf(buf,"%d",rstart);
		ocg->add_kword(ocg, 0, "RANDOM_START", buf, NULL);
	}

	if (hex) {
		ocg->add_kword(ocg, 0, "HEXAGON_PATCHES", "True", NULL);
	}

	generate_ps(itype, psname, cols, npat, label,
	            pap != NULL ? pap->w : cwidth, pap != NULL ? pap->h : cheight,
	            marg, rand, rstart, saix, paix,	ixord,
	            scale, hex, verb, scanc, eps, spacer, nmask, pcol, wp,
	            &sip, &pis, &plen, &glen, &tlen, &nppat);

	if (itype == instDTP41) {	/* DTP41 needs this */
		sprintf(buf,"%f",plen);
		ocg->add_kword(ocg, 0, "PATCH_LENGTH", buf, NULL);
		sprintf(buf,"%f",glen);
		ocg->add_kword(ocg, 0, "GAP_LENGTH", buf, NULL);
		sprintf(buf,"%f",tlen);
		ocg->add_kword(ocg, 0, "TRAILER_LENGTH", buf, NULL);
	}

	sprintf(buf,"%d",sip);
	ocg->add_kword(ocg, 0, "STEPS_IN_PASS", buf, NULL);

	/* Convert pass in strips count to base 62 */
	for (i = 0; ;i++) {
		if (pis[i] == 0)
			break;
		if (pis[i] <= 9)
			pis[i] += '0';
		else if (pis[i] <= 35)
			pis[i] += 'A' - 10;
		else if (pis[i] <= 62)
			pis[i] += 'a' - 36;
		else
			error("Too many passes in strip to be encode as base 62 (%d)",pis[i]);
	}
	ocg->add_kword(ocg, 0, "PASSES_IN_STRIPS", pis, NULL);

	/* Output the default Argyll style strip and patch numbering */
	ocg->add_kword(ocg, 0, "STRIP_INDEX_PATTERN", sixpat, NULL);
	ocg->add_kword(ocg, 0, "PATCH_INDEX_PATTERN", pixpat, NULL);
	ocg->add_kword(ocg, 0, "INDEX_ORDER", ixord ? "PATCH_THEN_STRIP" : "STRIP_THEN_PATCH", NULL);

	/* Write out the patch info to the output CGATS file */
	for (i = 0; i < nppat; i++) {
		cgats_set_elem ary[2 + ICX_MXINKS + 3];
		int j;

		if (strcmp(cols[i].loc, "???") == 0)
			warning ("Internal, patch %s (%d) wasn't given a valid location string",cols[i].id,i+1);
		ary[0].c = cols[i].id;
		ary[1].c = cols[i].loc;
		for (j = 0; j < nchan; j++)
			ary[2 + j].d = 100.0 * cols[i].dev[j];
		for (j = 0; j < 3; j++)
			ary[2 + nchan + j].d = 100.0 * cols[i].XYZ[j];
		ocg->add_setarr(ocg, 0, ary);
	}

	if (ocg->write_name(ocg, outname))
		error("Write error : %s",ocg->err);

	paix->del(paix);
	saix->del(saix);
	free(pis);
	free(cols);
	ocg->del(ocg);		/* Clean up */
	icg->del(icg);		/* Clean up */

	return 0;
}

/******************************************************************/
/* Edge tracking support, for generating the scanner image        */
/* recognition reference chart file. */

/* Establish width and height to convert between topleft and */
/* bottom left origin ??? ~~9 */

/*
	Basic algorithm strategy:

	First we simply accumulate the raw recognition and patch
	identification information. Once the chart is generated, we:
		sort into horizontal and vertical half edges
		sort into +ve and -ve edges
		match +ve and -ve edges
		for each match, generate a delta edge segment
		Assume any non-matched edges are against the media.
		Coalesce delta edges into X & Y edge lists
		Compute normalised strength.
		Compute crossing count.
		Figure average box size, and compute shrink.
		
*/

/* A half edge structure */
/* coordinate origin is top left */
struct _hedge {
	double rgb[3];	/* Color this half edge transitions to */
	int negh;		/* 1 if this is a -ve major coordinate side half edge, 0 otherwise */
	double mj;		/* Major coordinate offset (ie. X coord for vertical edge) */
	double mi0;		/* Minor coordinate smaller value (ie. Y for vertical edge) */
	double mi1;		/* Minor coordinate larger value (ie. Y for vertical edge) */
	struct _hedge *next;	/* Next in linked list */
}; typedef struct _hedge hedge;

/* A patch identifier */
/* coordinate origin is top left */
struct _patch {
	char id[20];	/* ID string, Zeri length if a diagnostic rectangle */
	double xo;		/* Location of the rectangle origin (bottom left ???) */
	double yo;
	double w;		/* Size of the patch */
	double h;
	struct _patch *next;	/* Next in linked list */
}; typedef struct _patch patch;


/* Structure to one edge */
struct _edge {
	double p1,p2;	/* Start and end of line in orthogonal direction */
	struct _edge *next;	/* next in the linked list */
}; typedef struct _edge edge;

/* Structure of an edge list */
struct _elist {
	double pos;		/* Position of edges along major axis */
	double len;		/* Total length of edges atthis position */
	double cc;		/* Crossing count */
	int ne;			/* Count of edges */
	edge *e;		/* Head of linked list of edges at this position */
	struct _elist *next;	/* Next in linked list */
}; typedef struct _elist elist;

/* - - - - - - - - - - - - - - - - - - - */
/* Structure to track recognition edges */
struct {
	double height;	/* Height of the page */
	double mrgb[3];	/* Media RGB */	
	double rgb[3];	/* Currently set RGB */	

	/* Raw half edge lists, [vertical, horizontal] */
	int nhe[2];
	hedge *he[2];

	/* Patch identity information */
	int npatches;
	patch *patches;

	/* Processed information */
	hedge **she[2];	/* Sorted half edges */

	int nel[2];		/* Number of edge positions */
	elist *el[2];	/* Head of edge linked list */
	elist **nelp;	/* Next edge to append to */

} et = {
	0.0,				/* Height */
	1.0, 1.0, 1.0,		/* media RGB */
	1.0, 1.0, 1.0,		/* Current RGB */
	{ 0, 0},
	{ NULL, NULL }, 
	0, NULL,
	{ NULL, NULL },
	{ 0, 0 },			/* No edge positions */
	{ NULL, NULL },		/* Head of edge linked list */
	NULL
};


/* Tell et of the height, so the Y coordinate can be flipped */
void et_height(double height) {
	et.height = height;
//printf("~1 media height set to %f\n",height);
}

/* Tell et of the media color */
void et_media(double *rgb) {
	int e;
	for (e = 0; e < 3; e++)
		et.mrgb[e] = rgb[e];
}

/* Track the current GC color */
void et_color(
double *rgb		/* New RGB values */
) {
	int e;
	for (e = 0; e < 3; e++)
		et.rgb[e] = rgb[e];
}

/* Track a drawn object half edge */
/* We assume that no object is written over any other object, */
/* and that each half edge has a perfect opposite edge (ie. same */
/* mi0 and mi1), or no matching half edge (it is over the media) - */
/* ie. no partialy overlapping half edges. */
/* The arguments origin is assumed to be bottom left */
void et_edge(
int isx,		/* NZ if this is a vertical edge */
int negh,		/* NZ if this is a -ve major coordinate side half edge */
double mj,		/* Major coordinate offset (ie. X coord for vertical edge) */
double mi0,		/* Minor coordinate smaller value (ie. Y for vertical edge) */
double mi1		/* Minor coordinate larger value (ie. Y for vertical edge) */
) {
	int e, h;
	hedge *he;

	if (mi1 < mi0)
		error ("et_edge, minor coords wern't ordered");

	if ((he = (hedge *)calloc(sizeof(hedge), 1)) == NULL)
		error("Malloc of half edge structure failed");

	for (e = 0; e < 3; e++)
		he->rgb[e] = et.rgb[e];

	/* Flip the Y coordinate */
	if (isx) {
		double tmi0, tmi1;
		tmi0 = et.height - mi1;		/* swap to keep smallest small */
		tmi1 = et.height - mi0;
		mi0 = tmi0;
		mi1 = tmi1;
	} else {
		mj = et.height - mj;
	}

	he->negh = negh ? 1 : 0;
	he->mj   = mj;
	he->mi0  = mi0;
	he->mi1  = mi1;

	/* Add half edges to the list */
	h = isx ? 0 : 1;
	et.nhe[h]++;
	he->next = et.he[h];
	et.he[h] = he;
}

/* Track a patch identity */
/* The arguments origin is assumed to be bottom left */
void et_patch(
char  *id,		/* ID string, NULL if a diagnostic rectangle */
double xo,		/* Bottom left of the rectangle */
double yo,
double w,		/* Size of the patch */
double h
) {
	patch *p;

//printf("~1 got patch at %f %f, w %f, h %f\n", xo, yo, w, h);

	if ((p = (patch *)calloc(sizeof(patch), 1)) == NULL)
		error("Malloc of patch structure failed");

	/* Flip Y */
	yo = et.height - (yo + h);

	if (id != NULL) {
		strncpy(p->id, id, 19);
		p->id[19] = '\000';
	} else {
		p->id[0] = '\000';
	}
	p->xo = xo;
	p->yo = yo;
	p->h  = h;
	p->w  = w;

	/* Add patch to list */
	et.npatches++;
	p->next = et.patches;
	et.patches = p;
}

/* Compute the image recognition information, and write the */
/* .cht file. */
void et_write(char *fname, col *cols, int *rix, int si, int ei) {
	FILE *of;
	hedge *ep;
	int i, h;

//printf("~1 et has %d vertical and %d horizontal half edges\n", et.nhe[0], et.nhe[1]);
//printf("~1 et has %d patches\n", et.npatches);

	for (h = 0; h < 2; h++) {
		/* Create sorted list of vertical half edges */
		if ((et.she[h] = (hedge **)malloc(sizeof(patch*) * et.nhe[h])) == NULL)
			error("Malloc of array of vertical halfedge pointers failed");
	
		for (ep = et.he[h], i = 0; i < et.nhe[h]; i++, ep = ep->next)
			et.she[h][i] = ep;
	
		/* Sort helf edges by their X location, then their Y0 location */
#define HEAP_COMPARE(A,B) (fabs(A->mj - B->mj) < 1e-6 ? A->mi0 < B->mi0 : A->mj < B->mj)
		HEAPSORT(hedge *, et.she[h], et.nhe[h]);
#undef HEAP_COMPARE

#ifdef NEVER
for (i = 0; i < et.nhe[h]; i++) {
printf("%s %d at %c = %f from %c = %f to %f\n",
h == 0 ? "Vert" : "Horiz", i,
h == 0 ? 'X' : 'Y', et.she[h][i]->mj,
h == 0 ? 'Y' : 'X', et.she[h][i]->mi0, et.she[h][i]->mi1);
}
#endif /* NEVER */

		et.nel[h] = 0;
		et.nelp = &et.el[h];		/* Append next edge list here */
		*et.nelp = NULL;			/* No edge list at this position yet */

		/* Create the edge list information */
		for (i = 0; i < et.nhe[h];) {
			int j, ii, nj;
			double *rgb;	/* Contrast RGB */
			elist *el;		/* Current elist */

			el = *et.nelp;

			/* Locate the end of the half edges at the same position */
			for (ii = i; ii < et.nhe[h]; ii++) {
				if (fabs(et.she[h][i]->mj - et.she[h][ii]->mj) > 1e-6)
					break;
			}

//printf("~1 doing group from %d to %d\n",i, ii);
			/* Find half edge pairs */
			/* Note that we assume that the half edges match perfectly, */
			/* or not at all. This will be normaly be the case with targets */
			/* generated by printtarg. */
			for (j = i; j < ii; j = nj, j++) {
				int e, k = j+1;
				double vv;

				if (k < ii
				 && et.she[h][j]->negh != et.she[h][k]->negh
				 && fabs(et.she[h][j]->mi0 - et.she[h][k]->mi0) < 1e-5
				 && fabs(et.she[h][j]->mi1 - et.she[h][k]->mi1) < 1e-5) {
					/* Found a matching half edge */

					nj = k;
					rgb = et.she[h][k]->rgb;

				} else if (k < ii	/* Assert */
				 && (   (et.she[h][j]->mi0+1e-6) < et.she[h][k]->mi1
				     && et.she[h][j]->mi1 > (et.she[h][k]->mi0+1e-6))) {

					/* Found an overlapping non-matching edge */
					nj = k;

#ifdef NEVER
fprintf(stderr,"%s %d at %c = %f from %c = %f to %f, half %s\n",
h == 0 ? "Vert" : "Horiz", i,
h == 0 ? 'X' : 'Y', et.she[h][j]->mj,
h == 0 ? 'Y' : 'X', et.she[h][j]->mi0, et.she[h][j]->mi1,
et.she[h][j]->negh ? "Neg" : "Pos");
fprintf(stderr,"%s %d at %c = %f from %c = %f to %f, half %s\n",
h == 0 ? "Vert" : "Horiz", i,
h == 0 ? 'X' : 'Y', et.she[h][k]->mj,
h == 0 ? 'Y' : 'X', et.she[h][k]->mi0, et.she[h][k]->mi1,
et.she[h][k]->negh ? "Neg" : "Pos");
#endif /* NEVER */
						error("Internal - half edges don't match");

				} else {
					/* Must be a non-matching edge */
					nj = j;
					rgb = et.mrgb;	/* Edge must be against media */
				}

				/* Compute vector delta in rgb */
				/* Add entry to edge list */
				for (e = 0, vv = 0.0; e < 3; e++) {
					double tt = rgb[e] - et.she[h][j]->rgb[e];
					vv += tt * tt;
				}

//printf("~1 h %d, mj %f, mi %f, vv = %f\n",h, et.she[h][j]->mj, et.she[h][j]->mi0, vv);
				/* If edge is of sufficient magnitude */
				// ~~99
				if (vv > 0.2) {
					edge *ep;

					/* See if we need to add a first elist for this position */
					if (el == NULL) {		/* We do */
						if ((el = (elist *)calloc(sizeof(elist), 1)) == NULL)
							error("Malloc of elist structure failed");
						*et.nelp = el;
						el->pos = et.she[h][j]->mj;
						et.nel[h]++;
					}

					/* Add another edge entry */
					if ((ep = (edge *)calloc(sizeof(edge), 1)) == NULL)
						error("Malloc of edge structure failed");

					ep->next = el->e;
					ep->p1 = et.she[h][j]->mi0;
					ep->p2 = et.she[h][j]->mi1;

					el->e = ep;		/* Add to edge list */
					el->ne++;
				}
			}

			if (el != NULL) {
				/* We've done that position, so get ready for next */
				et.nelp = &el->next;		/* Append any more positions here */
				*et.nelp = NULL;
			}
			i = ii;		/* Start search for next group here */
		}
	}

	/* Figure the crossing count */
	for (h = 0; h < 2; h++) {
		int oh = 1 - h;	/* The other list */
		elist *el, *pl;		/* Current, previous elist */

		for (pl = NULL, el = et.el[h]; el != NULL; pl = el, el = el->next) {
			edge *ep;
			double pp, np;	/* Window in pos direction for crossing */
			double ppos = el->pos;

			if (pl != NULL)
				pp = (pl->pos + el->pos)/2.0;	/* Half distance to next line */
			else
				pp = -1e6;

			if (el->next != NULL)
				np = (el->next->pos + el->pos)/2.0;	/* Half distance to next line */
			else
				np = 1e6;

			/* For each edge on this edge position */
			for (ep = el->e; ep != NULL; ep = ep->next) {
				elist *oel;		/* Other edge list pointer */

				/* For each edge in other list */
				for (oel = et.el[oh]; oel != NULL; oel = oel->next) {
					edge *oep;

					if (oel->pos < ep->p1 || oel->pos > ep->p2)
						continue;	/* Other edge doesn't intersect this segment */

					for (oep = oel->e; oep != NULL; oep = oep->next) {

						/* If crosses on this line within +-0.5 of line each side */
						if (oep->p1 <= np && oep->p2 >= pp) {
							el->cc++;	/* Count crossing */
						}

					}
				}
			}
		}
	}

	/* Compute and normalise the length (strength) & crossing count of each edge */
	for (h = 0; h < 2; h++) {
		elist *el;		/* Current elist */
		double maxlen = 0.0;  
		double maxcc = 0.0;  

		for (el = et.el[h]; el != NULL; el = el->next) {
			edge *ep;
			double tlen;

			for (tlen = 0.0, ep = el->e; ep != NULL; ep = ep->next) {
				tlen += ep->p2 - ep->p1;
			}
			el->len = tlen;
			if (maxlen < tlen)
				maxlen = tlen;
			if (maxcc < el->cc)
				maxcc = el->cc;
		}

		/* Normalise */
		for (el = et.el[h]; el != NULL; el = el->next) {
			el->len /= maxlen;
			el->cc  /= maxcc;
		}
	}

	/* Output the .cht file */
	if ((of = fopen(fname,"w")) == NULL)
			error ("Unable to open output file '%s' for writing",fname);
	
	fprintf(of,"\n\n");
	fprintf(of, "BOXES %d\n",et.npatches);

	{
		int fidc = 0;
		patch *pp;
		double mins;		/* Minimum sample box size */

		for (pp = et.patches; pp != NULL; pp = pp->next) {

			if (pp->id[0] == '\000') {
				fprintf(of, "  D FID%d FID%d _ _ %f %f %f %f 0 0\n",
				        fidc, fidc, pp->w, pp->h, pp->xo, pp->yo);
				fidc++;
			}
		}
		
		mins = 1e6;
		for (pp = et.patches; pp != NULL; pp = pp->next) {

			if (pp->id[0] != '\000') {
				fprintf(of, "  X %s %s _ _ %f %f %f %f 0 0\n",
				        pp->id, pp->id, pp->w, pp->h, pp->xo, pp->yo);
				if (mins > pp->w)
					mins = pp->w;
				if (mins > pp->h)
					mins = pp->h;
			}
		}
		fprintf(of,"\n");
		
		/* Use a 15% box shrink */
		fprintf(of, "BOX_SHRINK %f\n", mins * 0.15);
		fprintf(of,"\n");

	}

	fprintf(of,"REF_ROTATION %f\n", 0.0);
	fprintf(of,"\n");

	{
		elist *el;		/* Current elist */

		fprintf(of,"XLIST %d\n",et.nel[0]);
		for (el = et.el[0]; el != NULL; el = el->next)
			fprintf(of,"  %f %f %f\n",el->pos, el->len, el->cc);
		fprintf(of,"\n");
	
		fprintf(of,"YLIST %d\n",et.nel[1]);
		for (el = et.el[1]; el != NULL; el = el->next)
			fprintf(of,"  %f %f %f\n",el->pos, el->len, el->cc);
		fprintf(of,"\n");
	
		fprintf(of,"\n");
	}


	fprintf(of, "EXPECTED XYZ %d\n",ei - si);

	for (i = si; i < ei; i++) {
		int ix = rix[i];
		fprintf(of, "  %s %f %f %f\n", cols[ix].loc, 100.0 * cols[ix].XYZ[0], 100.0 * cols[ix].XYZ[1], 100.0 * cols[ix].XYZ[2]);
	}
	fprintf(of,"\n");

	if (fclose(of))
		error ("Unable to close output file '%s'",fname);
}


/* Cleanup any allocation */
void et_clear(void) {
	int h;
	patch *p;

	for (h = 0; h < 2; h++) {
		hedge *he;
		elist *el;

		/* Free up half edges */
		he = et.he[h];
		while (he != NULL) {
			hedge *ne = he->next;
			free(he);
			he = ne;
		}
		et.nhe[h] = 0;
		et.he[h] = NULL;
	
		/* Free up sorted half edge lists */
		if (et.she[h] != NULL)
			free(et.she[h]);
		et.she[h] = NULL;

		/* Free up edge lists and edges */
		el = et.el[h];
		while (el != NULL) {
			elist *nel;
			edge  *ep;

			ep = el->e;
			while (ep != NULL) {
				edge *nep;

				nep = ep->next;
				free(ep);
				ep = nep;
			}
			el->ne = 0;
			el->e = NULL;

			nel = el->next;
			free(el);
			el = nel;
		}
		et.nel[h] = 0;
		et.el[h] = NULL;
	}
	et.nelp = NULL;

	p = et.patches;
	while (p != NULL) {
		patch *np = p->next;
		free(p);
		p = np;
	}
	et.patches = NULL;
	et.npatches = 0;
}

/******************************************************************/
/* Error/debug output routines */
/******************************************************************/


/* Basic printf type error() and warning() routines */

void
error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"printtarg: Error - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	exit (-1);
}

void
warning(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"printtarg: Warning - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}







