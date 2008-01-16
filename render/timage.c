
/* 
 * Argyll Color Correction System
 *
 * RGB gamut boundary test image generator
 *
 * Author: Graeme W. Gill
 * Date:   28/12/2005
 *
 * Copyright 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

/*
 * Generate TIFF image with two RGB cube surface hexagons,
 * plus a rectangular grey wedges between them, on a grey
 * background, or a rectangular gamut surface test image.
 */

/*
 * TTBD:
 */

#undef DEBUG

#define verbo stdout

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "copyright.h"
#include "config.h"
#include "numsup.h"
#include "render.h"

#define DEF_DPI 200

void
usage(void) {
	fprintf(stderr,"Create test images, default hex RGB surface and wedge, Version %s\n",ARGYLL_VERSION_STR);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	fprintf(stderr,"usage: timage [-options] outfile.tif\n");
//	fprintf(stderr," -v             Verbose\n");
	fprintf(stderr," -t             Generate rectangular gamut boundary test chart\n");
	fprintf(stderr," -p steps       Generate a colorspace step chart with L* steps^2\n");
	fprintf(stderr," -r res         Resolution in DPI (default %d)\n",DEF_DPI);
	fprintf(stderr," -s             Smooth blend\n");
	fprintf(stderr," -x             16 bit output\n");
	fprintf(stderr," -g prop        Percentage towards grey (default 0%%)\n");
	fprintf(stderr," outfile.tif    Profile to check against\n");
	exit(1);
	}

int main(int argc, char *argv[])
{
	int fa,nfa;				/* current argument we're looking at */
	int verb = 0;
	int rchart = 0;			/* Rectangular chart */
	int schart = 0;			/* Step chart with steps^2 */
	int smooth = 0;			/* Use smooth blending */
	double res = DEF_DPI;
	depth2d depth = bpc8_2d;
	char outname[MAXNAMEL+1] = { 0 };	/* Output TIFF name */
	render2d *r;
	color2d c;
	double vv[4][2];
	color2d cc[4];
	double gbf = 1.0;		/* Grey blend factor */
	double w, h;			/* Size of page in mm */
	int i, j;

	error_program = "timage";

	if (argc <= 1)
		usage();

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

			/* Rectangular chart */
			else if (argv[fa][1] == 't' || argv[fa][1] == 'T') {
				rchart = 1;
				schart = 0;

			/* Smooth blending */
			} else if (argv[fa][1] == 's' || argv[fa][1] == 'S')
				smooth = 1;

			/* 16 bit depth */
			else if (argv[fa][1] == 'x' || argv[fa][1] == 'X')
				depth = bpc16_2d;

			/* step chart */
			else if (argv[fa][1] == 'p' || argv[fa][1] == 'P') {
				fa = nfa;
				if (na == NULL) usage();
				schart = atoi(na);
				if (schart <= 0) usage();
				rchart = 0;
			}

			/* resolution */
			else if (argv[fa][1] == 'r' || argv[fa][1] == 'R') {
				fa = nfa;
				if (na == NULL) usage();
				res = atof(na);
				if (res <= 0.0) usage();
			}

			/* grey blend */
			else if (argv[fa][1] == 'g' || argv[fa][1] == 'G') {
				fa = nfa;
				if (na == NULL) usage();
				gbf = 1.0 - 0.01 * atof(na);
				if (gbf < 0.0 || gbf > 1.0) usage();
			}
			else 
				usage();
		} else
			break;
	}

	/* Get the file name arguments */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strncpy(outname,argv[fa++],MAXNAMEL); outname[MAXNAMEL] = '\000';

	res /= 25.4;				/* Convert to DPmm */ 

	/* RGB Hexagon chart */
	if (rchart == 0 && schart == 0) {
		double r3o2;			/* 0.866025 */
		double bb = 0.07;		/* Border proportion */
		double hh = 40.0;		/* Height of hexagon in mm */

		r3o2 = sqrt(3.0)/2.0;		/* Width to heigh of hexagon */
		h = (1.0 + 2.0 * bb) * hh;
		w = (4.0 * bb + 0.25 + 2.0 * r3o2) * hh;
	
		if ((r = new_render2d(w, h, res, res, rgb_2d, depth)) == NULL) {
			error("new_render2d() failed");
		}
	
		/* Set the default color */
		c[0] = 0.5;
		c[1] = 0.5;
		c[2] = 0.5;
		r->set_defc(r, c);
	
		/* Left hand hex */
		vv[0][0] = hh * bb + r3o2 * 0.5 * hh;
		vv[0][1] = hh * bb + hh/2.0;
		cc[0][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[0][1] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[0][2] = 1.0 * gbf + (1.0 - gbf) * 0.5;
	
		vv[1][0] = hh * bb + r3o2 * 0.5 * hh;
		vv[1][1] = hh * bb;
		cc[1][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
	
		vv[2][0] = hh * bb;
		vv[2][1] = hh * bb + 0.25 * hh;
		cc[2][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][2] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
		vv[1][0] = hh * bb;
		vv[1][1] = hh * bb + 0.75 * hh;
		cc[1][0] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][2] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
		vv[2][0] = hh * bb + r3o2 * 0.5 * hh;
		vv[2][1] = hh * bb + 1.0 * hh;
		cc[2][0] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][1] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][2] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
		vv[1][0] = hh * bb + r3o2 * 1.0 * hh;;
		vv[1][1] = hh * bb + 0.75 * hh;
		cc[1][0] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][1] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
		vv[2][0] = hh * bb + r3o2 * 1.0 * hh;
		vv[2][1] = hh * bb + 0.25 * hh;
		cc[2][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][1] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
		vv[1][0] = hh * bb + r3o2 * 0.5 * hh;;
		vv[1][1] = hh * bb;
		cc[1][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
	
		/* Right hand hex */
		vv[0][0] = hh * (3.0 * bb + 0.25 + r3o2) + r3o2 * 0.5 * hh;
		vv[0][1] = hh * bb + hh/2.0;
		cc[0][0] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[0][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[0][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
	
		vv[1][0] = hh * (3.0 * bb + 0.25 + r3o2) + r3o2 * 0.5 * hh;
		vv[1][1] = hh * bb;
		cc[1][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
	
		vv[2][0] = hh * (3.0 * bb + 0.25 + r3o2);
		vv[2][1] = hh * bb + 0.25 * hh;
		cc[2][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][2] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
		vv[1][0] = hh * (3.0 * bb + 0.25 + r3o2);
		vv[1][1] = hh * bb + 0.75 * hh;
		cc[1][0] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][2] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
		vv[2][0] = hh * (3.0 * bb + 0.25 + r3o2) + r3o2 * 0.5 * hh;
		vv[2][1] = hh * bb + 1.0 * hh;
		cc[2][0] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][1] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][2] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
		vv[1][0] = hh * (3.0 * bb + 0.25 + r3o2) + r3o2 * 1.0 * hh;;
		vv[1][1] = hh * bb + 0.75 * hh;
		cc[1][0] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][1] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
		vv[2][0] = hh * (3.0 * bb + 0.25 + r3o2) + r3o2 * 1.0 * hh;
		vv[2][1] = hh * bb + 0.25 * hh;
		cc[2][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][1] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
		vv[1][0] = hh * (3.0 * bb + 0.25 + r3o2) + r3o2 * 0.5 * hh;;
		vv[1][1] = hh * bb;
		cc[1][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_trivs2d(vv, cc));
	
		/* Center wedge */
		cc[0][0] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[0][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[0][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][0] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[1][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][1] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[2][2] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[3][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[3][1] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		cc[3][2] = 1.0 * gbf + (1.0 - gbf) * 0.5;
		r->add(r, new_rectvs2d((2.0 * bb + r3o2) * hh, bb * hh, 0.25 * hh, hh, cc));

	/* RGB Rectangular chart */
	} else if (schart == 0) {
		double bb = 0.07;		/* Border proportion */
		double hh = 50.0;		/* Height of hexagon in mm */
		double sc[6][3] = {		/* Saturated color sequence */
			{ 1, 0, 0 },
			{ 1, 0, 1 },
			{ 0, 0, 1 },
			{ 0, 1, 1 },
			{ 0, 1, 0 },
			{ 1, 1, 0 }
		};

		h = (1.0 + 2.0 * bb) * hh;
		w = (2.0 * bb + 0.20 * 7.0) * hh;
	
		if ((r = new_render2d(w, h, res, res, rgb_2d, depth)) == NULL) {
			error("new_render2d() failed");
		}
	
		/* Set the default color */
		c[0] = 0.5;
		c[1] = 0.5;
		c[2] = 0.5;
		r->set_defc(r, c);
	
		for (i = 0; i < 7; i++) {
			prim2d *p;

			/* Top rectangle */
			cc[0][0] = sc[i % 6][0] * gbf + (1.0 - gbf) * 0.5;
			cc[0][1] = sc[i % 6][1] * gbf + (1.0 - gbf) * 0.5;
			cc[0][2] = sc[i % 6][2] * gbf + (1.0 - gbf) * 0.5;
			cc[1][0] = sc[(i+1) % 6][0] * gbf + (1.0 - gbf) * 0.5;
			cc[1][1] = sc[(i+1) % 6][1] * gbf + (1.0 - gbf) * 0.5;
			cc[1][2] = sc[(i+1) % 6][2] * gbf + (1.0 - gbf) * 0.5;
			cc[2][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
			cc[2][1] = 1.0 * gbf + (1.0 - gbf) * 0.5;
			cc[2][2] = 1.0 * gbf + (1.0 - gbf) * 0.5;
			cc[3][0] = 1.0 * gbf + (1.0 - gbf) * 0.5;
			cc[3][1] = 1.0 * gbf + (1.0 - gbf) * 0.5;
			cc[3][2] = 1.0 * gbf + (1.0 - gbf) * 0.5;
			p = new_rectvs2d((bb + i * 0.2) * hh, (bb + 0.5) * hh, 0.2 * hh, 0.5 * hh, cc);
			if (smooth) {
				rectvs2d *pp = (rectvs2d *)p;
				pp->x_blend = 2;
				pp->y_blend = 3;
			}
			r->add(r, p);

			/* Bottom rectangle */
			cc[0][0] = 0.0 * gbf + (1.0 - gbf) * 0.5;
			cc[0][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
			cc[0][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
			cc[1][0] = 0.0 * gbf + (1.0 - gbf) * 0.5;
			cc[1][1] = 0.0 * gbf + (1.0 - gbf) * 0.5;
			cc[1][2] = 0.0 * gbf + (1.0 - gbf) * 0.5;
			cc[2][0] = sc[i % 6][0] * gbf + (1.0 - gbf) * 0.5;
			cc[2][1] = sc[i % 6][1] * gbf + (1.0 - gbf) * 0.5;
			cc[2][2] = sc[i % 6][2] * gbf + (1.0 - gbf) * 0.5;
			cc[3][0] = sc[(i+1) % 6][0] * gbf + (1.0 - gbf) * 0.5;
			cc[3][1] = sc[(i+1) % 6][1] * gbf + (1.0 - gbf) * 0.5;
			cc[3][2] = sc[(i+1) % 6][2] * gbf + (1.0 - gbf) * 0.5;
			p = new_rectvs2d((bb + i * 0.2) * hh, bb * hh, 0.2 * hh, 0.5 * hh, cc);
			if (smooth) {
				rectvs2d *pp = (rectvs2d *)p;
				pp->x_blend = 2;
				pp->y_blend = 2;
			}
			r->add(r, p);
		}

	} else {	/* Lab step chart */
		double hh = 50.0;		/* Height of hexagon in mm */
		double bb = 0.05;		/* Border proportion */
		double ss, bs;			/* Step size, border size */

		h = hh;
		w = hh;

		bs = (bb * hh)/(schart + 1.0);
		ss = hh * (1.0 - bb)/schart;
	
		if ((r = new_render2d(w, h, res, res, lab_2d, depth)) == NULL) {
			error("new_render2d() failed");
		}
	
		/* Set the default color */
		c[0] = 0.0;
		c[1] = 0.0;
		c[2] = 0.0;
		r->set_defc(r, c);
	
		for (i = 0; i < schart; i++) {
			for (j = 0; j < schart; j++) {
				double lv;
				
				lv = (double)(j * schart + i)/(schart * schart - 1.0) * 100.0;

				cc[0][0] = lv;
				cc[0][1] = -127.0;
				cc[0][2] = -127.0;
				cc[1][0] = lv;
				cc[1][1] =  127.0;
				cc[1][2] = -127.0;
				cc[2][0] = lv;
				cc[2][1] = -127.0;
				cc[2][2] =  127.0;
				cc[3][0] = lv;
				cc[3][1] =  127.0;
				cc[3][2] =  127.0;
				r->add(r, new_rectvs2d(bs + i * (bs + ss),
				                       bs + j * (bs + ss),
				                       ss, ss, cc));
			}
		}
	}

	r->write(r, outname);
	r->del(r);

	return 0;
}


