
/* 
 * Argyll Color Correction System
 * DTP92/Spectrolino display callibrator.
 *
 * Author: Graeme W. Gill
 * Date:   14/10/2005
 *
 * Copyright 1996 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* This program displays test patches, and takes readings from a display device */
/* in order to create a RAMDAC calibration curve (usually stored in the ICC vcgt tag) */

/* This is the third version of the program. */

/* TTBD
 */

#undef DEBUG
#undef DEBUG_OFFSET		/* Keep test window out of the way */
#undef DEBUG_PLOT		/* Plot curve each time around */
#undef FAKE_DEVICE		/* Test without a spectrometer using a fake device */
#define MATRIX			/* Use matrix for aprox inversion (seems to be slightly better) */

#define COMPORT 1		/* Default com port 1..4 */
#define REFINE_GAIN 0.5	/* Refinement correction gain */

#define VERBOUT stdout

#if defined(__APPLE__)
#define DEF_GAMMA 1.8
#else
#define DEF_GAMMA 2.2
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <string.h>
#include "copyright.h"
#include "config.h"
#include "numlib.h"
#include "xspect.h"
#include "xcolorants.h"
#include "cgats.h"
#include "serio.h"
#include "inst.h"
#include "dtp92.h"
#include "gretag.h"
#include "dispwin.h"
#include "dispsup.h"
#include "rspl.h"
#include "moncurve.h"
#include "targen.h"
#include "ofps.h"

#include "sort.h"

#include <stdarg.h>

#if defined (NT)
#include <conio.h>
#endif

#define VLENGTH(aa) (sqrt(aa[0] * aa[0] + aa[1] * aa[1] + aa[2] * aa[2]))

/* =================================================================== */

/* Make aXYZ black relative */
static void krel(double dst[3], double src[3], double bk[3]) {
	int j;
	for (j = 0; j < 3; j++) {
		dst[j] = src[j] - bk[j];
		if (dst[j] < 0.0)
			dst[j] = 0.0;
	}
}

/* Make aXYZ black relative */
static void XYZNkrel(icmXYZNumber *dst, double src[3], double bk[3]) {
	int j;
	double dd[3];
	for (j = 0; j < 3; j++) {
		dd[j] = src[j] - bk[j];
		if (dd[j] < 0.0)
			dd[j] = 0.0;
	}
	dst->X = dd[0];
	dst->Y = dd[1];
	dst->Z = dd[2];
}

/* - - - - - - - - - - - - - - - - - - */
/* Adaptive sampling support. */
/* Simple, non-optimised implementation... */

typedef struct {
	double v;			/* Input value */
	double xyz[3];		/* Reading */
	double lab[3];		/* Reading */
} sxyz;

typedef struct {
	int no;				/* Number of samples */
	int _no;			/* Allocation */
	sxyz *s;			/* List of samples */
} asamp;

static void init_asamp(asamp *p) {
	p->no = 0;
	p->_no = 0;
	p->s = NULL;
}

static void free_alloc_asamp(asamp *p) {
	if (p->s != NULL)
		free(p->s);
	p->s = NULL;
}

/* Add a new sample. Samples are kept sorted. */
static void asamp_add(asamp *p, double v, double *lab, double *xyz) {
	if (p->no >= p->_no) {
		if (p->_no == 0) {
			p->_no = 10;
			if ((p->s = (sxyz *)malloc(p->_no * sizeof(sxyz))) == NULL)
				error("asamp malloc failed");
		} else {
			p->_no *= 2;
			if ((p->s = (sxyz *)realloc(p->s, p->_no * sizeof(sxyz))) == NULL)
				error("asamp realloc failed");
		}
	}

	p->s[p->no].v = v;
	p->s[p->no].xyz[0] = xyz[0];
	p->s[p->no].xyz[1] = xyz[1];
	p->s[p->no].xyz[2] = xyz[2];
	p->s[p->no].lab[0] = lab[0];
	p->s[p->no].lab[1] = lab[1];
	p->s[p->no].lab[2] = lab[2];
	p->no++;

	/* Sort the current samples */
#define HEAP_COMPARE(A,B) (A.v < B.v)
	HEAPSORT(sxyz,p->s,p->no)
#undef HEAP_COMPARE

}

/* Return the next sample location */
static double asamp_next(asamp *p) {
	int i, j;
	double bdist;
	int bi;

	if (p->no < 2)
		error("Calling asamp_next with less than two existing samples");

	/* Locate the pair of sample with the largest reading + input difference */
	bdist = -1.0;
	for (i = 0; i < (p->no-1); i++) {
		double dd, tt;

		/* 50% Lab difference */
		for (dd = 0.0, j = 0; j < 3; j++) {
			tt = p->s[i].lab[j] - p->s[i+1].lab[j]; 
			tt /= 100.0;
			dd += tt * tt;
		}
		tt = p->s[i].v - p->s[i+1].v; 				/* Input difference */

		/* and 50% input difference */	
		dd += 2.0 * tt * tt;
		
		if (dd > bdist) {
			bdist = dd;
			bi = i;
		}
	}
//printf("Biggest was %d - %d\n",bi,bi+1);

	return (0.5 * (p->s[bi].v + p->s[bi+1].v)); 
}

/* Diagnostic: */
/* Return a linear interpolation */
static void asamp_interp(asamp *p, double xyz[3], double v) {
	int i, j;
	double b;

	if (p->no < 2)
		error("Calling asamp_interp with less than two existing samples");

	/* Locate the pair surrounding our input value */
	for (i = 0; i < (p->no-1); i++) {
		if (v >= p->s[i].v && v <= p->s[i+1].v)
			break;
	}
	if (i >= (p->no-1))
		error("asamp_interp out of range");

	
	b = (v - p->s[i].v)/(p->s[i+1].v - p->s[i].v);

	for (j = 0; j < 3; j++) {
		xyz[j] = b * p->s[i+1].xyz[j] + (1.0 - b) * p->s[i].xyz[j];
	}
}

/* ------------------------------------------------------------------- */

/* device RGB inverse solution code, using powell() */

/* Context for calibration solution */
typedef struct {
	double wh[3];		/* White absolute XYZ value (after subtracting bk) */
	double bk[3];		/* Black absolute XYZ value */
	double twh[3];		/* Target white absolute XYZ value */
	icmXYZNumber twN;	/* Same as above as XYZNumber */
	double fm[3][3];	/* Forward, aprox. linear RGB -> XYZ */
	double bm[3][3];	/* Backwards, aprox. XYZ -> linear RGB */
#ifdef MATRIX
	mcv *cvs[3];		/* Device RGB channel to sum XYZ curves */
#else
	mcv *cvs[3][3];		/* Device RGB channel XYZ curves */
#endif	/* !MATRIX */
	mcv *rdac[3];		/* Current RGB to RGB ramdac curves */

	double xyz[3];		/* Target xyz value */
} cctx;

#ifdef MATRIX		/* Version using matrix */

/* Return the xyz that is expected to be generated */
/* by the given device RGB. */
static void fwddev(cctx *x, double xyz[3], double rgb[3]) {
	double lrgb[3];
	int j;

	/* Convert device RGB into linear light RGB via curves */
	for (j = 0; j < 3; j++)
		lrgb[j] = x->cvs[j]->interp(x->cvs[j], rgb[j])
		        / x->cvs[j]->interp(x->cvs[j], 1.0);

	/* Convert linear light RGB into XYZ via the matrix */
	icmMulBy3x3(xyz, x->fm, lrgb);
}

/* Return the closest device RGB that is expected */
/* to generate the given xyz. */
/* Return 0 on exact match, 1 on nearest match */
static int invdev(cctx *x, double rgb[3], double xyz[3]) {
	double lrgb[3];
	int j;

	/* Convert XYZ to linear light RGB via the inverse matrix */
	icmMulBy3x3(lrgb, x->bm, xyz);

	/* Convert linear light RGB to device RGB via inverse curves */
	for (j = 0; j < 3; j++)
		rgb[j] = x->cvs[j]->inv_interp(x->cvs[j], lrgb[j] * x->cvs[j]->interp(x->cvs[j], 1.0));

	return 0;
}

#else /* !MATRIX */			/* Version using powell and sum of curve XYZs */

/* Return the xyz that is expected to be generated */
/* by the given device RGB. */
static void fwddev(cctx *x, double xyz[3], double rgb[3]) {
	int i, j;

	for (j = 0; j < 3; j++)
		xyz[j] = 0.0;

	/* Sum RGB's XYZs */
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			xyz[j] += x->cvs[i][j]->interp(x->cvs[i][j], rgb[i]);
		}
	}
}


/* Function being solved */
static double rev_fcn(	/* Return < 0 on abort */
	void *fdata,		/* Opaque data pointer */
	double *rgb			/* Multivariate input values */
) {
	double dd, rv = 0.0;
	int i, j;
	double xyz[3];		/* XYZ from current rgb */
	cctx *x = (cctx *)fdata;

//printf("\n~1 got rgb %f %f %f\n",rgb[0],rgb[1],rgb[2]);
	/* Discourage out of gamut */
	for(j = 0; j < 3; j++) {
		if (rgb[j] < 0.0) {
			rv = -rgb[j];
		} else if (rgb[j] > 1.0) {
			rv = rgb[j] - 1.0;
		}
	}
	rv *= 20000.0;
//printf("~1 rv after clip = %f\n",rv);

	/* Encourage highest device value, */
	/* to skip any "dead" band from 0 */
	for(j = 0; j < 3; j++) {
		rv += 0.5 * (1.0 - rgb[j]);
	}
//printf("~1 rv enc high = %f\n",rv);

	fwddev(x, xyz, rgb);

//printf("~1 got xyz %f %f %f for target %f %f %f\n", xyz[0], xyz[1], xyz[2], x->xyz[0], x->xyz[1], x->xyz[2]);
	/* Difference to XYZ target */
	for(dd = 0.0, j = 0; j < 3; j++) {
		double tt;
		tt = xyz[j] - x->xyz[j];
		dd += tt * tt;
	}
	rv += dd;

//printf("~1 rev_fnc returning %f (err %f)\n", rv);
	return rv;
}

/* Return the closest device RGB that is expected */
/* to generate the given xyz. */
/* Assume rgb[] is set to starting values. */
/* Return 0 on exact match, 1 on nearest match */
static int invdev(cctx *x, double rgb[3], double xyz[3]) {
	int i, j;
	int rv = 0;
	double sa[3] = { 0.2, 0.2, 0.2 };	/* Search area */
	double txyz[3];				/* Function values at solution */

	icmAry2Ary(x->xyz, xyz);	/* Targe XYZ */

	if (powell(3, rgb, sa, 1e-7, 2000, rev_fcn, x) < 0.0)
		error("Powell failed");  

	/* Clip to gamut */
	for (j = 0; j < 3; j++) {
		if (rgb[j] < 0.0) {
			rgb[j] = 0.0;
			rv = 1;
		} else if (rgb[j] > 1.0) {
			rgb[j] = 1.0;
			rv = 1;
		}
	}

	return rv;
}

#endif /* !MATRIX */

/* - - - - - - - - - - - - - - - - - - - */

#ifdef NEVER
/* Do a linear interp of the ramdac */
static void interp_ramdac(double cal[256][3], double drgb[3], double rgb[3]) {
	int i, j;
	int gres = 256;
	double w;

	/* For r,g & b */
	for (j = 0; j < 3; j++) {
		int mi, gres_1 = gres-1;
		double t, vv = rgb[j];
		t = gres * vv;
		mi = (int)floor(t);			/* Grid coordinate */
		if (mi < 0)					/* Limit to valid cube base index range */
			mi = 0;
		else if (mi >= gres_1)
			mi = gres_1-1;
		w = t - (double)mi;	 		/* 1.0 - weight */

		drgb[j] = (1.0 - w) * cal[mi][j] + w * cal[mi+1][j];
	}
}
#endif	/* NEVER */

#ifdef FAKE_DEVICE
/* Test without a spectrometer using a fake device */
static int fake_read(
	col *cols,		/* Array of patch colors to be tested */
	int npat, 		/* Number of patches to be tested */
	int spat,		/* Start patch index for "verb", 0 if not used */
	int tpat		/* Total patch index for "verb", 0 if not used */
) {
	icmXYZNumber white;		/* White point */
	icmXYZNumber red;		/* Red colorant */
	icmXYZNumber green;		/* Green colorant */
	icmXYZNumber blue;		/* Blue colorant */
	double mat[3][3];		/* Destination matrix */
	double rgb[3];
	double br = 35.4;		/* Overall brightness */
	int patch, j;
	int ttpat = tpat;

	if (ttpat < npat)
		ttpat = npat;

	/* Setup fake device */
	white.X = br * 0.955;		/* Somewhere between D50 and D65 */
	white.Y = br * 1.00;
	white.Z = br * 0.97;
	red.X = br * 0.41;
	red.Y = br * 0.21;
	red.Z = br * 0.02;
	green.X = br * 0.30;
	green.Y = br * 0.55;
	green.Z = br * 0.15;
	blue.X = br * 0.15;
	blue.Y = br * 0.10;
	blue.Z = br * 0.97;

	if (icmRGBprim2matrix(white, red, green, blue, mat))
		error("Fake read unexpectedly got singular matrix\n");

	for (patch = 0; patch < npat; patch++) {
		printf("\rpatch %d of %d",spat+patch+1,ttpat); fflush(stdout);

		rgb[0] = cols[patch].r;
		rgb[1] = cols[patch].g;
		rgb[2] = cols[patch].b;

		/* Apply gamma */
		for (j = 0; j < 3; j++)
			rgb[j] = pow(rgb[j], 2.5);

		/* Convert to XYZ */
		icmMulBy3x3(cols[patch].aXYZ, mat, rgb);

		cols[patch].aXYZ_v = 1;
	}

	return 0;
}




#endif

/* =================================================================== */

void
usage(void) {
	serio *sio;
	fprintf(stderr,"Calibrate a Display, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL\n");
	fprintf(stderr,"usage: dispcal [-v] [-c comport] [-a] [-i inst] [-s]%s outfile\n",
#if defined(UNIX) && !defined(__APPLE__)
	" [-n]"
#else
	""
#endif
    );
	fprintf(stderr," -v              Verbose mode\n");
	fprintf(stderr," -d              Print debug diagnostics to stderr\n");
	fprintf(stderr," -c comport      Set COM port, 1..N (default %d)\n",COMPORT);
	if ((sio = new_serio()) != NULL) {
		char **paths;
		if ((paths = sio->get_paths(sio)) != NULL) {
			int i;
			for (i = 0; ; i++) {
				if (paths[i] == NULL)
					break;
				fprintf(stderr,"    %d = '%s'\n",i+1,paths[i]);
			}
		} else
			fprintf(stderr,"    ** No ports found **\n");
		sio->del(sio);
	}
	fprintf(stderr," -a              Run instrument calibration\n");
 	fprintf(stderr," -i 92|SO        Select target instrument (default DTP92)\n");
	fprintf(stderr,"                 92 = DTP92, SO = Spectrolino\n");

	fprintf(stderr," -q [lmhu]       Quality - Low, Medium (def), High, Ultra\n");
	fprintf(stderr," -t temp         Set the target white point daylight color temperature in deg. K\n");
	fprintf(stderr," -w x,y        	 Set the target white point as chromaticity coordinates\n");
	fprintf(stderr," -b bright       Set the target brightness in cd/m^2\n");
	fprintf(stderr," -g gamma        Set the target response curve gamma (Def. %3.1f)\n",DEF_GAMMA);
	fprintf(stderr,"                 Use \"-gl\" for L*a*b* curve\n");
	fprintf(stderr,"                 Use \"-gs\" for sRGB curve\n");
	fprintf(stderr," -V              Run verify pass on final curves\n");

#if defined(UNIX) && !defined(__APPLE__)
	fprintf(stderr," -n              Don't set override redirect on test window\n");
#endif
	fprintf(stderr," outfile         Base name for output .cal file\n");
	exit(1);
	}

main(int argc, char *argv[])
{
	int i, j, k;
	int fa,nfa;							/* current argument we're looking at */
	double patsize = 80.0;				/* size of displayed color patch */
	int verb = 0;
	int debug = 0;
	int override = 1;					/* Override redirect on X11 */
	int comport = COMPORT;				/* COM port used */
	instType itype = instDTP92;			/* Default target instrument */
	int doCalib = 0;					/* Do an instrument calibration */
	double dtemp = 0.0;					/* Daylight color temperature (native) */
	double wpx = 0.0, wpy = 0.0;		/* White point xy (native) */
	double tbright = 0.0;				/* Target brightness (max)  */
	int gammat = 0;						/* Gamma type, 0 = number, 1 = Lab, 2 = sRGB */
	double gamma = DEF_GAMMA;			/* Gamma number */
	int isteps = 22;					/* Initial measurement steps/3 (medium) */
	int rsteps = 64;					/* Refinement measurement steps (medium) */
	int mxits = 3;						/* maximum itterations (medium) */
	int verify = 0;						/* Do a verify after last refinement */
	static char outname[200] = { 0 };	/* Output cgats file base name */
	int spectral = 0;					/* Want spectral data from instrument */
	double ho = 0.0, vo = 0.0;
	disprd *dr = NULL;					/* Display patch read object */
	cctx x;
	int rv;

	error_program = "dispcal";
	if (argc <= 1)
		usage();

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')		/* Look for any flags */
			{
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else {
				if ((fa+1) < argc) {
					if (argv[fa+1][0] != '-')
						{
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
					}
				}
			}

			if (argv[fa][1] == '?')
				usage();

			else if (argv[fa][1] == 'v')
				verb = 1;

			else if (argv[fa][1] == 'd' || argv[fa][1] == 'D')
				debug = 1;

			else if (argv[fa][1] == 'V')
				verify = 1;

#if defined(UNIX) && !defined(__APPLE__)
			else if (argv[fa][1] == 'n' || argv[fa][1] == 'N')
				override = 0;
#endif /* UNIX */
			/* COM port  */
			else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				fa = nfa;
				if (na == NULL) usage();
				comport = atoi(na);
				if (comport < 1 || comport > 4) usage();
			}
			/* Target Instrument */
			else if (argv[fa][1] == 'i' || argv[fa][1] == 'I') {
				char *p;
				int tt;
				fa = nfa;
				if (na == NULL) usage();

				if (strcmp("92", na) == 0) {
					itype = instDTP92;
					patsize = 80.0;
				}
				else if (strcmp("so", na) == 0
				      || strcmp("SO", na) == 0) {
					itype = instSpectrolino;
					patsize = 60.0;
				}
				else
					usage();
			}
			else if (argv[fa][1] == 'a' || argv[fa][1] == 'A') {
				doCalib = 1;

			/* Quality */
			} else if (argv[fa][1] == 'q' || argv[fa][1] == 'Q') {
				fa = nfa;
				if (na == NULL) usage();
    			switch (na[0]) {
					case 'L':
						isteps = 3;
						rsteps = 9;
						mxits = 1;
						break;
					case 'l':
						isteps = 11;
						rsteps = 32;
						mxits = 2;
						break;
					case 'm':
					case 'M':
						isteps = 22;
						rsteps = 64;
						mxits = 3;
						break;
					case 'h':
					case 'H':
						isteps = 32;
						rsteps = 96;
						mxits = 4;
						break;
					case 'u':
					case 'U':
						isteps = 43;
						rsteps = 128;
						mxits = 5;
						break;
					default:
						usage();
				}

			/* Daylight color temperature */
			} else if (argv[fa][1] == 't' || argv[fa][1] == 'T') {
				fa = nfa;
				if (na == NULL) usage();
				dtemp = atof(na);
				if (dtemp < 1000.0 || dtemp > 15000.0) usage();

			/* White point as x, y */
			} else if (argv[fa][1] == 'w' || argv[fa][1] == 'W') {
				fa = nfa;
				if (na == NULL) usage();
				if (sscanf(na, " %lf,%lf ", &wpx, &wpy) != 2)
					usage();

			/* Brightness */
			} else if (argv[fa][1] == 'b' || argv[fa][1] == 'B') {
				fa = nfa;
				if (na == NULL) usage();
				tbright = atof(na);
				if (tbright <= 0.0 || tbright > 100000.0) usage();

			/* Gamma */
			} else if (argv[fa][1] == 'g' || argv[fa][1] == 'G') {
				fa = nfa;
				if (na == NULL) usage();
				if ((na[0] == 'l' || na[0] == 'L') && na[1] == '\000')
					gammat = 1;
				else if ((na[0] == 's' || na[0] == 'S') && na[1] == '\000')
					gammat = 2;
				else {
					gamma = atof(na);
					if (gamma <= 0.0 || gamma > 10.0) usage();
					gammat = 0;
				}
			} else 
				usage();
		} else
			break;
	}

	if (doCalib) {
		if ((rv = disprd_calibration(itype, comport, override, patsize, verb, debug)) != 0) {
			error("docalibration failed with return value %d\n",rv);
		}
	}

	/* Get the file name argument */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(outname,argv[fa]);
	strcat(outname,".cal");

	if (verb) {
		if (wpx > 0.0 || wpy > 0.0)
			printf("Target white = xy %f %f\n",wpx,wpy);
		else if (dtemp > 0.0)
			printf("Target white = %f degrees kelvin\n",dtemp);
		else
			printf("Target white = native white\n");

		if (tbright > 0.0)
			printf("Target brightness = %f cd/m^2\n",tbright);
		else
			printf("Target brightness = native maximum\n",tbright);

		if (gammat == 0)
			printf("Target gamma = %f\n",gamma);
		else if (gammat == 1)
			printf("Target gamma = L* curve\n");
		else if (gammat == 2)
			printf("Target gamma = sRGB curve\n");
	}

	/* Get ready to do some readings */
#ifdef DEBUG_OFFSET
	ho = 0.8;
	vo = -0.8;
#endif
#ifndef FAKE_DEVICE
	if ((dr = new_disprd(itype, comport, 1, NULL, override, patsize, ho, vo,
	                     spectral, verb, VERBOUT, debug)) == NULL)
		error("dispread failed to create a disprd object\n");
#endif
	{
		icmXYZNumber mrd;		/* Number for matrix */
		icmXYZNumber mgn;
		icmXYZNumber mbl;
		icmXYZNumber mwh;
		double *mat[3];			/* Forward matrix in numlib format */
		double tt[3], xx[3];	/* RHS and solution */
		double k;				/* Smallest k */

		col base[6] = {	/* Base set of test colors */
			{ 0.0, 0.0, 0.0 },		/* 0 - Black */
			{ 1.0, 0.0, 0.0 },		/* 1 - Red */
			{ 0.0, 1.0, 0.0 },		/* 2 - Green */
			{ 0.0, 0.0, 1.0 },		/* 3 - Blue */
			{ 1.0, 1.0, 1.0 },		/* 4 - White */
			{ 0.0, 0.0, 0.0 }		/* 5 - Black */
		};

		/* Read the base test set */
#ifdef FAKE_DEVICE
		fake_read(base, 6, 0, 0);
#else
		if ((rv = dr->read(dr, base, 6, 0, 0)) != 0) {
			error("display read failed with value %d\n",rv);
		} 
#endif

		if (verb)
			printf("\n");
		if (base[0].aXYZ_v == 0)
			error("Failed to get an absolute XYZ value from the instrument!\n");

		/* Average black relative from 2 readings */
		x.bk[0] = 0.5 * (base[0].aXYZ[0] + base[5].aXYZ[0]);
		x.bk[1] = 0.5 * (base[0].aXYZ[1] + base[5].aXYZ[1]);
		x.bk[2] = 0.5 * (base[0].aXYZ[2] + base[5].aXYZ[2]);

		/* Make other readings black relative */
		krel(x.wh, base[1].aXYZ, x.bk); icmAry2XYZ(mrd, x.wh);
		krel(x.wh, base[2].aXYZ, x.bk); icmAry2XYZ(mgn, x.wh);
		krel(x.wh, base[3].aXYZ, x.bk); icmAry2XYZ(mbl, x.wh);
		krel(x.wh, base[4].aXYZ, x.bk); icmAry2XYZ(mwh, x.wh);

		if (verb) {
			printf("Black = XYZ %f %f %f\n",x.bk[0],x.bk[1],x.bk[2]);
			printf("Red   = XYZ %f %f %f\n",base[1].aXYZ[0], base[1].aXYZ[1], base[1].aXYZ[2]);
			printf("Green = XYZ %f %f %f\n",base[2].aXYZ[0], base[2].aXYZ[1], base[2].aXYZ[2]);
			printf("Blue  = XYZ %f %f %f\n",base[3].aXYZ[0], base[3].aXYZ[1], base[3].aXYZ[2]);
			printf("White = XYZ %f %f %f\n",base[4].aXYZ[0], base[4].aXYZ[1], base[4].aXYZ[2]);
		}

		/* Figure out the target white point */
		if (wpx > 0.0 || wpy > 0.0) {	/* xy coordinates */
			double Yxy[3];
			Yxy[0] = 1.0;
			Yxy[1] = wpx;
			Yxy[2] = wpy;
			icmYxy2XYZ(x.twh, Yxy);

		} else if (dtemp > 0.0) {		/* Daylight color temperature */
			icx_DTEMP2XYZ(x.twh, dtemp);

		} else {						/* Native white */
			x.twh[0] = x.wh[0]/x.wh[1];
			x.twh[1] = x.wh[1]/x.wh[1];
			x.twh[2] = x.wh[2]/x.wh[1];
		}

		/* Convert it to absolute white target */
		if (tbright > 0.0) {			/* Given brightness */
			x.twh[0] *= tbright;
			x.twh[1] *= tbright;
			x.twh[2] *= tbright;
		} else {						/* Maximum brightness */
			x.twh[0] *= x.wh[1];
			x.twh[1] *= x.wh[1];
			x.twh[2] *= x.wh[1];
		}

		/* Now see if target white will fit in gamut. */
		/* Compute addive fwd matrix */
		if (icmRGBprim2matrix(mwh, mrd, mgn, mbl, x.fm))
			error("Aprox. fwd matrix unexpectedly singular\n");

		if (verb) {
			printf("Forward matrix is:\n");
			printf("%f %f %f\n", x.fm[0][0], x.fm[0][1], x.fm[0][2]);
			printf("%f %f %f\n", x.fm[1][0], x.fm[1][1], x.fm[1][2]);
			printf("%f %f %f\n", x.fm[2][0], x.fm[2][1], x.fm[2][2]);
		}

		/* Comute bwd matrix */
		if (icmInverse3x3(x.bm, x.fm))
			error("Inverting aprox. fwd matrix failed");

		/* Setup for intersection solution */

		k = 1e6;
		mat[0] = x.fm[0];
		mat[1] = x.fm[1];
		mat[2] = x.fm[2];
		for (i = 0; i < 3; i++) {	/* Intersect with upper 3 faces of cube */
			/* Save column we're working on and set RHS of equation */
			tt[0] = mat[0][i];
			tt[1] = mat[1][i];
			tt[2] = mat[2][i];
			xx[0] = -mat[0][i];
			xx[1] = -mat[1][i];
			xx[2] = -mat[2][i];

			/* Set column to white point vector */
			mat[i][0] = -x.twh[0];
			mat[i][1] = -x.twh[1];
			mat[i][2] = -x.twh[2];

			if (verb) {
				printf("A.X = B to solve is:\n");
				printf("%f %f %f . X = %f\n", x.fm[0][0], x.fm[0][1], x.fm[0][2], xx[0]);
				printf("%f %f %f . X = %f\n", x.fm[1][0], x.fm[1][1], x.fm[1][2], xx[1]);
				printf("%f %f %f . X = %f\n", x.fm[2][0], x.fm[2][1], x.fm[2][2], xx[2]);
			}

			/* Solve for scale factor and other 2 device values */
			if (solve_se(mat, xx, 3)) {
				warning ("White point vector intesection failure");
			} else {
printf("~1 got k[%d] = %f\n",i,xx[i]);
				if (xx[i] > 0.0 && xx[i] < k)
					k = xx[i];
			}

			/* Restor matrix column */
			mat[0] = x.fm[0];		/* Rows might have been flipped */
			mat[1] = x.fm[1];
			mat[2] = x.fm[2];
			mat[0][i] = tt[0];
			mat[1][i] = tt[1];
			mat[2][i] = tt[2];
		}
		if (k < 1.0) {
			x.twh[0] *= k;		/* Scale brightness to fit */
			x.twh[1] *= k;
			x.twh[2] *= k;
			if (verb)
				printf("Had to scale brightness to %f to fit within gamut\n",x.twh[1]);
		}
		if (verb)
			printf("Target white value is XYZ %f %f %f\n",x.twh[0],x.twh[1],x.twh[2]);

		icmAry2XYZ(x.twN, x.twh);		/* Need this for Lab conversions */
	}

	{
		asamp asrgb[3];			/* Adaptive sampling for r, g & b */
		asamp asgrey;			/* Adaptive sampling for r=g=b */
		col set[3];
		int it;					/* verify & refine itteration */
		
#ifdef MATRIX
		for (k = 0; k < 3; k++) {
			if ((x.cvs[k] = new_mcv()) == NULL)
				error("new_mcv x.cvs[%d] failed",k);
		}
#else /* !MATRIX */
		for (k = 0; k < 3; k++) {
			for (j = 0; j < 3; j++) {
				if ((x.cvs[k][j] = new_mcv()) == NULL)
					error("new_mcv x.cvs[%d][%d] failed",k,j);
			}
		}
#endif /* MATRIX */

		for (j = 0; j < 3; j++)
			init_asamp(&asrgb[j]);

		/* Read bracket values */
		for (j = 0; j < 3; j++) {
			set[j].r = j == 0 ? 1.0 : 0.0;
			set[j].g = j == 1 ? 1.0 : 0.0;
			set[j].b = j == 2 ? 1.0 : 0.0;
			set[j].id = NULL;
		}
#ifdef FAKE_DEVICE
		fake_read(set, 3, 0, (isteps-1) * 3);
#else
		if ((rv = dr->read(dr, set, 3, 0, (isteps-1) * 3)) != 0) {
			error("display read failed with value %d\n",rv);
		} 
#endif

		/* setup the brackets */
		for (j = 0; j < 3; j++) {
			double xyz[3], lab[3];
			xyz[0] = xyz[1] = xyz[2] = 0.0;
			icmXYZ2Lab(&x.twN, lab, xyz);
			asamp_add(&asrgb[j], 0.0, lab, xyz);
			krel(set[j].aXYZ, set[j].aXYZ, x.bk);		/* Convert to bk relative */
			icmXYZ2Lab(&x.twN, lab, set[j].aXYZ);
			asamp_add(&asrgb[j], 1.0, lab, set[j].aXYZ);
		}

		for (i = 1; i < (isteps-1); i++) {
			double v[3];
			for (j = 0; j < 3; j++) {
// ~~99
				v[j] = asamp_next(&asrgb[j]);
//				v[j] = i/(isteps-1.0);
//				v[j] = pow(v[j], 2.2);
				if (j == 0)
					set[j].r = v[j];
				else if (j == 1)
					set[j].g = v[j];
				else
					set[j].b = v[j];
			}
//printf("~1 R = %f, G = %f, B = %f\n", v[0], v[1], v[2]);

#ifdef FAKE_DEVICE
			fake_read(set, 3, i * 3, (isteps-1) * 3);
#else
			if ((rv = dr->read(dr, set, 3, i * 3, (isteps-1) * 3)) != 0) {
				error("display read failed with value %d\n",rv);
			} 
#endif
			for (j = 0; j < 3; j++) {
				double lab[3];
				krel(set[j].aXYZ, set[j].aXYZ, x.bk);		/* Convert to bk relative */
				icmXYZ2Lab(&x.twN, lab, set[j].aXYZ);
				asamp_add(&asrgb[j], v[j], lab, set[j].aXYZ);
			}
		}
		if (verb)
			printf("\n");

		/* Convert RGB channel samples to curves */
		{
			mcvco *sdv;				/* Scattered data for mcv */
	
			if ((sdv = malloc(sizeof(mcvco) * isteps)) == NULL)
				error ("Malloc of scattered data points failed");
#ifdef MATRIX
			for (k = 0; k < 3; k++) {		/* RGB */
				for (i = 0; i < isteps; i++) {
					sdv[i].p = asrgb[k].s[i].v;
					sdv[i].v = VLENGTH(asrgb[k].s[i].xyz);
				}
				x.cvs[k]->fit(x.cvs[k], 0, 20, sdv, isteps, 1.0);
			}
#else /* !MATRIX */
			for (k = 0; k < 3; k++) {		/* RGB */
				for (j = 0; j < 3; j++) {	/* XYZ */
					for (i = 0; i < isteps; i++) {
						sdv[i].p = asrgb[k].s[i].v;
						sdv[i].v = asrgb[k].s[i].xyz[j];
					}
					x.cvs[k][j]->fit(x.cvs[k][j], 0, 20, sdv, isteps, 1.0);
				}
			}
#endif /* !MATRIX */
			free (sdv);
		}

#ifdef DEBUG_PLOT
		/* Plot the current calc curves */
#ifdef MATRIX
		{
			#define	XRES 256
			double xx[XRES];
			double yy[3][XRES];
			double xyz[3];
			for (i = 0; i < XRES; i++) {
				xx[i] = i/(XRES-1.0);
				for (j = 0; j < 3; j++)
					yy[j][i] = x.cvs[j]->interp(x.cvs[j], xx[i]);
			}
			printf("Channel curves\n");
			do_plot(xx,yy[0],yy[1],yy[2],XRES);
			#undef XRES
		}
#else /* !MATRIX */
		for (k = 0; k < 3; k++) {
			#define	XRES 256
			double xx[XRES];
			double yy[6][XRES];
			double xyz[3];
			for (i = 0; i < XRES; i++) {
				xx[i] = i/(XRES-1.0);
				for (j = 0; j < 3; j++)
					yy[j][i] = x.cvs[k][j]->interp(x.cvs[k][j], xx[i]);
				/* Plot linear for sample points too */
				asamp_interp(&asrgb[k], xyz, xx[i]);
				for (j = 0; j < 3; j++)
					yy[3+j][i] = xyz[j];
			}
			printf("Channel %d curves\n",k);
//			do_plot6(xx,yy[0],yy[1],yy[2],yy[3],yy[4],yy[5],XRES);
			do_plot(xx,yy[0],yy[1],yy[2],XRES);
			#undef XRES
		}
#endif /* !MATRIX */
#endif

		/* We're done with asrgb[] */
		for (j = 0; j < 3; j++)
			free_alloc_asamp(&asrgb[j]);

		/* Calculate the current calibration curve values */
		{
			int nsamp = 128;
			mcvco *sdv[3];				/* Scattered data for mcv */

			for (j = 0; j < 3; j++) {
				if ((x.rdac[j] = new_mcv()) == NULL)
					error("new_mcv x.rdac[%d] failed",j);
			}

			for (j = 0; j < 3; j++) {
				if ((sdv[j] = malloc(sizeof(mcvco) * nsamp)) == NULL)
					error ("Malloc of scattered data points failed");
			}

			/* Generate the sample points */
			for (i = 0; i < nsamp; i++) {
				double vv = i/(nsamp - 1.0);
				double Lab[3], xyz[3], rgb[3];

				if (gammat == 0 || gammat == 2) {
					double y;
					if (gammat == 0) {
						y = pow(vv, gamma);
					} else {
						if (vv <= 0.03928)
							y = vv/12.92;
						else
							y = pow((0.055 + vv)/1.055, 2.4);
					}

					/* Convert Y to L* */
					if (y > 0.008856451586)
						y = pow(y,1.0/3.0);
					else
						y = 7.787036979 * y + 16.0/116.0;
					Lab[0] = 116.0 * y - 16.0;

				} else {	/* L* curve */
					Lab[0] = 100.0 * vv;
				}
				Lab[1] = Lab[2] = 0.0;	/* Target is neutral */

				/* Got Lab aim value */
				icmLab2XYZ(&x.twN, xyz, Lab);		/* XYZ Value to look for */
//printf("~1 %d: target %f %f %f, ",i,xyz[0],xyz[1],xyz[2]);

				/* Lookup device RGB for that target */
				rgb[0] = rgb[1] = rgb[2] = vv;
				invdev(&x, rgb, xyz);
//printf("rgb %f %f %f\n",rgb[0],rgb[1],rgb[2]);

				for (j = 0; j < 3; j++) {
					sdv[j][i].p = vv;
					sdv[j][i].v = rgb[j];
				}
			}
			for (j = 0; j < 3; j++)
				x.rdac[j]->fit(x.rdac[j], 0, 20, sdv[j], nsamp, 1.0);

			/* Make sure that if we are using native brightness and white point, */
			/* that the curves go to a perfect 1.0 ... */
			if (wpx == 0.0 && wpy == 0.0  && dtemp == 0.0 && tbright == 0.0) {
				for (j = 0; j < 3; j++)
					x.rdac[j]->force_scale(x.rdac[j], 1.0);
			}

			for (j = 0; j < 3; j++)
				free (sdv[j]);
		}

#ifdef DEBUG_PLOT
		/* Plot the current curves */
		{
			#define	XRES 255
			double xx[XRES];
			double y1[XRES];
			double y2[XRES];
			double y3[XRES];
			double rgb[3];
			for (i = 0; i < XRES; i++) {
				double drgb[3], rgb[3];
				xx[i] = i/(XRES-1.0);
				rgb[0] = rgb[1] = rgb[2] = xx[i];
				for (j = 0; j < 3; j++)
					drgb[j] = x.rdac[j]->interp(x.rdac[j], rgb[j]);
				y1[i] = drgb[0];
				y2[i] = drgb[1];
				y3[i] = drgb[2];
			}
			printf("Current ramdac curves\n");
			do_plot(xx,y1,y2,y3,XRES);
			#undef XRES
		}
#endif

		/* Now we go into the verify & refine loop */
		for (it = 0; it < mxits || verify != 0; it++) {
			double drgb[3], rgb[3];
			double xyz[3], lab[3];

			/* Verify exit */
			if (it >= mxits)
				rsteps = 100;		/* Fixed verification resolution */

			/* Test out the current curves */
			init_asamp(&asgrey);

			rgb[0] = rgb[1] = rgb[2] = 1.0;
			for (j = 0; j < 3; j++)
				drgb[j] = x.rdac[j]->interp(x.rdac[j], rgb[j]);
			set[0].r = drgb[0];
			set[0].g = drgb[1];
			set[0].b = drgb[2];
			set[0].id = NULL;

#ifdef FAKE_DEVICE
			fake_read(set, 1, 0, rsteps-1);
#else
			if ((rv = dr->read(dr, set, 1, 0, rsteps-1)) != 0) {
				error("display read failed with value %d\n",rv);
			} 
#endif

			/* setup the brackets */
			xyz[0] = xyz[1] = xyz[2] = 0.0;
			icmXYZ2Lab(&x.twN, lab, xyz);
			asamp_add(&asgrey, 0.0, lab, xyz);			/* Black */

			krel(set[0].aXYZ, set[0].aXYZ, x.bk);		/* Convert to bk relative */
			icmXYZ2Lab(&x.twN, lab, set[0].aXYZ);
			asamp_add(&asgrey, 1.0, lab, set[0].aXYZ);	/* White */

			/* Continue with the sampling */
			for (i = 1; i < (rsteps-1); i++) {
				double v;

				if (it >= mxits) {		/* Verify, fixed spacing samples */
					v = i/(rsteps - 1.0);
					v = pow(v, 2.2);

				} else {				/* Normal, adaptive etc. */
// ~~999
					v = asamp_next(&asgrey);
//					v = i/(rsteps - 1.0);
//					v = pow(v, 2.2);
				}

				rgb[0] = rgb[1] = rgb[2] = v;
				for (j = 0; j < 3; j++)
					drgb[j] = x.rdac[j]->interp(x.rdac[j], rgb[j]);
				set[0].r = drgb[0];
				set[0].g = drgb[1];
				set[0].b = drgb[2];
				set[0].id = NULL;
#ifdef FAKE_DEVICE
				fake_read(set, 1, i, rsteps-1);
#else
				if ((rv = dr->read(dr, set, 1, i, rsteps-1)) != 0) {
					error("display read failed with value %d\n",rv);
				} 
#endif
				krel(set[0].aXYZ, set[0].aXYZ, x.bk);		/* Convert to bk relative */
				icmXYZ2Lab(&x.twN, lab, set[0].aXYZ);
				asamp_add(&asgrey, v, lab, set[0].aXYZ);
			}
			if (verb)
				printf("\n");

#ifdef DEBUG_PLOT
			/* Plot the measured response XYZ */
			{
				#define	XRES 256
				double xx[XRES];
				double yy[3][XRES];
				double xyz[3];
				for (i = 0; i < XRES; i++) {
					xx[i] = i/(XRES-1.0);
					asamp_interp(&asgrey, xyz, xx[i]);
					for (j = 0; j < 3; j++)
						yy[j][i] = xyz[j];
				}
				printf("Measured neutral axis XYZ\n",k);
				do_plot(xx,yy[0],yy[1],yy[2],XRES);
				#undef XRES
			}
#endif
			/* Check out the results: */
			{
				double ctwh[3];		/* Current target white */
				icmXYZNumber ctwN;	/* Same as above as XYZNumber */
				double brerr;		/* Brightness error */
				double cterr;		/* Color temperature delta E */
				double mnerr;		/* Maximum neutral error */
				double mnv;			/* Value where maximum error is */
				double anerr;		/* Average neutral error */
				double lab1[3], lab2[3];
				
				/* Brightness */
				brerr = asgrey.s[asgrey.no-1].xyz[1] - x.twh[1];
			
				/* Compensate for brightness error */
				for (j = 0; j < 3; j++)
					ctwh[j] = x.twh[j] * asgrey.s[asgrey.no-1].xyz[1]/x.twh[1];
				icmAry2XYZ(ctwN, ctwh);		/* Need this for Lab conversions */
				
				/* Color temperature error */
				icmXYZ2Lab(&ctwN, lab1, ctwh);		/* Should be 100,0,0 */
				icmXYZ2Lab(&ctwN, lab2, asgrey.s[asgrey.no-1].xyz);
				cterr = icmLabDE(lab1, lab2);

				/* Compensate for white point error */
				icmAry2Ary(ctwh, asgrey.s[asgrey.no-1].xyz);
				icmAry2XYZ(ctwN, ctwh);		/* Need this for Lab conversions */
				
				/* check delta E of all the other sample points */
				mnerr = anerr = 0.0;
				for (i = asgrey.no-2; i >= 0; i--) {
					double err;
					double vv = asgrey.s[i].v;

					/* Compute aim Lab */
					if (gammat == 0 || gammat == 2) {
						double y;
						if (gammat == 0) {
							y = pow(vv, gamma);
						} else {
							if (vv <= 0.03928)
								y = vv/12.92;
							else
								y = pow((0.055 + vv)/1.055, 2.4);
						}

						/* Convert Y to L* */
						if (y > 0.008856451586)
							y = pow(y,1.0/3.0);
						else
							y = 7.787036979 * y + 16.0/116.0;
						lab1[0] = 116.0 * y - 16.0;

					} else {	/* L* curve */
						lab1[0] = 100.0 * vv;
					}
					lab1[1] = lab1[2] = 0.0;	/* Target is neutral */
					icmXYZ2Lab(&ctwN, lab2, asgrey.s[i].xyz);
					err = icmLabDE(lab1, lab2);
					if (err > mnerr) {
						mnerr = err;
						mnv = vv;
					}
					anerr += err;
				}
				anerr /= (asgrey.no-1.0);		/* Don't count white */
	
				if (verb || verify) {
					printf("Brightness error = %f cd/m^2\n",brerr);
					printf("White point error = %f deltaE\n",cterr);
					printf("Maximum neutral error (@ %f) = %f deltaE\n",mnv, mnerr);
					printf("Average neutral error = %f deltaE\n",anerr);
				}
			}

			/* Verify exit */
			if (it >= mxits) {
				free_alloc_asamp(&asgrey);
				break;
			}

			/* Compute a set of corrected device values */
			/* for each of the measured points, and then */
			/* compute new ramdac curves */
			{
				double lab[3], xyz[3];		/* Target values */
				double rgb[3], txyz[3];		/* Theoretical output xyz */
				double trgb[3], crgb[3];	/* Theoretical target rgb, theoretical current rgb */
				double fxyz[3];
				double frgb1[3] /*, frgb2[3] */;	/* Fixed rgb's by two methods */
				mcvco *sdv[3];				/* Scattered data for mcv */

				for (j = 0; j < 3; j++) {
					if ((sdv[j] = malloc(sizeof(mcvco) * asgrey.no)) == NULL)
						error ("Malloc of scattered data points failed");
				}

				for (i = 0; i < asgrey.no; i++) {
					double vv;

					rgb[0] = rgb[1] = rgb[2] = vv = asgrey.s[i].v;
					/* Compute aim Lab and XYZ */
					if (gammat == 0 || gammat == 2) {
						double y;
						if (gammat == 0) {
							y = pow(vv, gamma);
						} else {
							if (vv <= 0.03928)
								y = vv/12.92;
							else
								y = pow((0.055 + vv)/1.055, 2.4);
						}

						/* Convert Y to L* */
						if (y > 0.008856451586)
							y = pow(y,1.0/3.0);
						else
							y = 7.787036979 * y + 16.0/116.0;
						lab[0] = 116.0 * y - 16.0;

					} else {	/* L* curve */
						lab[0] = 100.0 * vv;
					}
					lab[1] = lab[2] = 0.0;	/* Target is neutral */
					icmLab2XYZ(&x.twN, xyz, lab);		/* XYZ Value we're aiming for */

					/* Lookup theoretical device RGB for target xyz */
					trgb[0] = trgb[1] = trgb[2] = vv;
					invdev(&x, trgb, xyz);

					/* Lookup theoretical device RGB for measured xyz */ 
					crgb[0] = crgb[1] = crgb[2] = vv;
					invdev(&x, crgb, asgrey.s[i].xyz);
				
					/* Comute fixed rgb's using method 1, with damping. */
					for (j = 0; j < 3; j++)
						frgb1[j] = x.rdac[j]->interp(x.rdac[j], vv)
					    + REFINE_GAIN * (trgb[j] - crgb[j]);

#ifdef NEVER
					/* Lookup theoretical xyz for target RGB */
					fwddev(&x, txyz, rgb);

					/* Comute fixed rgb's using method 2 */
					for (j = 0; j < 3; j++)
						fxyz[j] = txyz[j] + xyz[j] - asgrey.s[i].xyz[j];
					frgb2[0] = frgb2[1] = frgb2[2] = vv;
					invdev(&x, frgb2, fxyz);
#endif	/* NEVER */

					/* Use fixed rgb's */
					for (j = 0; j < 3; j++) {
						sdv[j][i].p = vv;
						sdv[j][i].v = frgb1[j];
					}
				}
				for (j = 0; j < 3; j++)
					x.rdac[j]->fit(x.rdac[j], 0, 15, sdv[j], asgrey.no, 1.0);

				/* Make sure that if we are using native brightness and white point, */
				/* that the curves go to a perfect 1.0 ... */
				if (wpx == 0.0 && wpy == 0.0  && dtemp == 0.0 && tbright == 0.0) {
					for (j = 0; j < 3; j++)
						x.rdac[j]->force_scale(x.rdac[j], 1.0);
				}

				for (j = 0; j < 3; j++)
					free (sdv[j]);
			}

#ifdef DEBUG_PLOT
			/* Plot the current curves */
			{
				#define	XRES 255
				double xx[XRES];
				double y1[XRES];
				double y2[XRES];
				double y3[XRES];
				double rgb[3];
				for (i = 0; i < XRES; i++) {
					double drgb[3], rgb[3];
					xx[i] = i/(XRES-1.0);
					rgb[0] = rgb[1] = rgb[2] = xx[i];
					for (j = 0; j < 3; j++)
						drgb[j] = x.rdac[j]->interp(x.rdac[j], rgb[j]);
					y1[i] = drgb[0];
					y2[i] = drgb[1];
					y3[i] = drgb[2];
				}
				printf("Current ramdac curves\n");
				do_plot(xx,y1,y2,y3,XRES);
				#undef XRES
			}
#endif
			free_alloc_asamp(&asgrey);
		}	/* Next refine loop */

		/* Write out the resulting calibration file */
		{
			int calres = 256;	/* 256 steps in calibration */
			cgats *ocg;			/* output cgats structure */
			time_t clk = time(0);
			struct tm *tsp = localtime(&clk);
			char *atm = asctime(tsp);	/* Ascii time */
			cgats_set_elem *setel;		/* Array of set value elements */

			ocg = new_cgats();				/* Create a CGATS structure */
			ocg->add_other(ocg, "CAL"); 	/* our special type is Calibration file */

			ocg->add_table(ocg, tt_other, 0);	/* Add a table for RAMDAC values */
			ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Device Calibration State",NULL);
			ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll dispcal", NULL);
			atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
			ocg->add_kword(ocg, 0, "CREATED",atm, NULL);

			ocg->add_kword(ocg, 0, "DEVICE_CLASS","DISPLAY", NULL);

			ocg->add_field(ocg, 0, "RGB_I", r_t);
			ocg->add_field(ocg, 0, "RGB_R", r_t);
			ocg->add_field(ocg, 0, "RGB_G", r_t);
			ocg->add_field(ocg, 0, "RGB_B", r_t);

			if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * 4)) == NULL)
				error("Malloc failed!");

			for (i = 0; i < calres; i++) {
				double vv = i/(calres-1.0);
				double rgb[3];
	
				for (j = 0; j < 3; j++) {
					double cc;
					cc = x.rdac[j]->interp(x.rdac[j], vv);
					if (cc < 0.0)
						cc = 0.0;
					else if (cc > 1.0)
						cc = 1.0;
					rgb[j] = cc;
				}

				setel[0].d = vv;
				setel[1].d = rgb[0];
				setel[2].d = rgb[1];
				setel[3].d = rgb[2];
	
				ocg->add_setarr(ocg, 0, setel);
			}
	
			free(setel);

			if (ocg->write_name(ocg, outname))
				error("Write error : %s",ocg->err);

			ocg->del(ocg);		/* Clean up */
		}
	}

	for (j = 0; j < 3; j++)
		x.rdac[j]->del(x.rdac[j]);

#ifdef MATRIX
	for (k = 0; k < 3; k++) {
		x.cvs[k]->del(x.cvs[k]);
	}
#else /* !MATRIX */
	for (k = 0; k < 3; k++) {
		for (j = 0; j < 3; j++) {
			x.cvs[k][j]->del(x.cvs[k][j]);
		}
	}
#endif /* !MATRIX */

	/* Now we're done */
#ifndef FAKE_DEVICE
	dr->del(dr);
#endif

	return 0;
}


