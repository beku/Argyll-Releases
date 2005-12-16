
/* 
 * International Color Consortium color transform expanded support
 *
 * Author:  Graeme W. Gill
 * Date:    11/11/2004
 * Version: 1.00
 *
 * Copyright 2000-2004 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the LICENCE.TXT file for licencing details.
 *
 */

/*
 * This is some test code to test the CIECAM02 functionality. 
 * This test creates a .tiff file containing values that
 * have been converted through CIECAM02 space. 
 */


#include <stdio.h>
#include <math.h>
#include "xcam.h"
#include "cam02.h"
#include "tiffio.h"

#define USE_HK 0	/* Use Helmholtz-Kohlraush */

#ifndef _isnan
#define _isnan(x) ((x) != (x))
#define _finite(x) ((x) == (x))
#endif

static void
Lab2XYZ(double *out, double *in) {
	double L = in[0], a = in[1], b = in[2];
	double x,y,z,fx,fy,fz;

	if (L > 8.0) {
		fy = (L + 16.0)/116.0;
		y = pow(fy,3.0);
	} else {
		y = L/903.2963058;
		fy = 7.787036979 * y + 16.0/116.0;
	}

	fx = a/500.0 + fy;
	if (fx > 24.0/116.0)
		x = pow(fx,3.0);
	else
		x = (fx - 16.0/116.0)/7.787036979;

	fz = fy - b/200.0;
	if (fz > 24.0/116.0)
		z = pow(fz,3.0);
	else
		z = (fz - 16.0/116.0)/7.787036979;

	out[0] = x * 0.9642;
	out[1] = y * 1.0000;
	out[2] = z * 0.8249;
}

/* CIE XYZ to perceptual Lab */
static void
XYZ2Lab(double *out, double *in) {
	double X = in[0], Y = in[1], Z = in[2];
	double x,y,z,fx,fy,fz;
	double L;

	x = X/0.9642;
	y = Y/1.0000;
	z = Z/0.8249;

	if (x > 0.008856451586)
		fx = pow(x,1.0/3.0);
	else
		fx = 7.787036979 * x + 16.0/116.0;

	if (y > 0.008856451586)
		fy = pow(y,1.0/3.0);
	else
		fy = 7.787036979 * y + 16.0/116.0;

	if (z > 0.008856451586)
		fz = pow(z,1.0/3.0);
	else
		fz = 7.787036979 * z + 16.0/116.0;

	out[0] = 116.0 * fy - 16.0;
	out[1] = 500.0 * (fx - fy);
	out[2] = 200.0 * (fy - fz);
}

/* Convert from XYZ to sRGB */
void XYZ2sRGB(double *out, double *in) {
	double tmp[3];
	int i;

	/* Now convert to sRGB assuming D50 (??) white point */
	tmp[0] = in[0] * 3.2410  + in[1] * -1.5374 + in[2] * -0.4986;
	tmp[1] = in[0] * -0.9692 + in[1] * 1.8760  + in[2] * 0.0416;
	tmp[2] = in[0] * 0.0556  + in[1] * -0.2040 + in[2] * 1.0570;

	/* Apply gamma and clip */
	for(i = 0; i < 3; i++) {
		if (tmp[i] < 0.00304) {
			out[i] = tmp[i] * 12.92;
		} else {
			out[i] = 1.055 * pow(tmp[i], 1.0/2.4) - 0.055;
		}
		if (out[i] < 0.0)
			out[i] = 0.0;
		else if (out[i] > 1.0)
			out[i] = 1.0;
	}
}

int
main(void) {
	char *tiffname = "cam02plot.tif";
	TIFF *wh = NULL;
	uint16 depth, bps;
	uint16 pconfig, photometric;
	uint16 resunits = RESUNIT_INCH;
	float resx = 75.0, resy = 75.0;
	tdata_t *obuf;
	unsigned char *ob;
	
	double white[6][3] = {
		{ 0.9505, 1.000, 1.0888 },
		{ 0.9505, 1.000, 1.0888 },
		{ 1.0985, 1.000, 0.3558 },
		{ 1.0985, 1.000, 0.3558 },
		{ 0.9505, 1.0000, 1.0890 }, /* D65 */
		{ 0.9642, 1.000, 0.8249 }	/* D50 for inversion tests */
	};
	double La[6] = { 318.31, 31.83, 318.31, 31.83, 318.31, 150.0 };

	int i, j, k, m;
	int e;
	double xyz[3], Jab[3], rgb[3];
	cam02 *cam1, *cam2;

	int res = 32;	/* Resolution of scan through space */
	int ares;		/* Array of 2d sub-rasters */
	int w, h;		/* Width and height of raster */
	int x, y;
	int sp = 5;			/* Spacing beween resxres sub-images */

	/* Setup cam to convert to Jab */
	cam1 = new_cam02();
	cam1->set_view(
		cam1,
		vc_average,	/* Enumerated Viewing Condition */
		white[4],	/* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
		0.20,		/* Relative Luminance of Background to reference white */
		1000.0,		/* Adapting/Surround Luminance cd/m^2 */
		0.0,		/* Luminance of white in image - not used */
		0.00,		/* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
		white[4],	/* The Flare color coordinates (typically the Ambient color) */
		USE_HK		/* use Helmholtz-Kohlraush flag */ 
	);
	
	/* Setup cam to convert from Jab */
	cam2 = new_cam02();
	cam2->set_view(
		cam2,
		vc_average,	/* Enumerated Viewing Condition */
		white[5],	/* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
		0.20,		/* Relative Luminance of Background to reference white */
		100.0,		/* Adapting/Surround Luminance cd/m^2 */
		0.0,		/* Luminance of white in image - not used */
		0.00,		/* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
		white[5],	/* The Flare color coordinates (typically the Ambient color) */
		USE_HK		/* use Helmholtz-Kohlraush flag */ 
	);
	
	/* Figure out the size of the raster */
	ares = (int)ceil(sqrt((double)res));
	printf("~1 for res %d, ares = %d\n",res,ares);
	
	w = ares * res + sp * (ares+1);
	h = ares * res + sp * (ares+1);

	printf("~1 raster width = %d, height = %d\n",w,h);
	/* Setup the tiff file */
	if ((wh = TIFFOpen(tiffname, "w")) == NULL)
		error("Can\'t create TIFF file '%s'!",tiffname);
	
	TIFFSetField(wh, TIFFTAG_IMAGEWIDTH,  w);
	TIFFSetField(wh, TIFFTAG_IMAGELENGTH, h);
	TIFFSetField(wh, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
	TIFFSetField(wh, TIFFTAG_SAMPLESPERPIXEL, 3);
	TIFFSetField(wh, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(wh, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(wh, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	TIFFSetField(wh, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(wh, TIFFTAG_RESOLUTIONUNIT, resunits);
	TIFFSetField(wh, TIFFTAG_XRESOLUTION, resx);
	TIFFSetField(wh, TIFFTAG_YRESOLUTION, resy);
	TIFFSetField(wh, TIFFTAG_IMAGEDESCRIPTION, "cam02plot");

	obuf = _TIFFmalloc(TIFFScanlineSize(wh));
	ob = (unsigned char *)obuf;

	y = 0;

	/* Fill sp lines with black */
	for (x = 0, m = 0; m < w; m++) {
		ob[x] = ob[x+1] = ob[x+2] = 0;
		x += 3;
	}
	for (j = 0; j < sp; j++, y++) {
//printf("~1 writing first blank lines at %d\n",y);
		TIFFWriteScanline(wh, ob, y, 0);
	}

	for (i = 0; i < ares; i++) {				/* Vertical blocks */
		for (j = 0; j < res; j++, y++) {		/* Vertical in block */
			xyz[1] = j/(res-1.0);
			x = 0;

			/* Fill sp pixels with black */
			for (m = 0; m < sp; m++) {
				ob[x] = ob[x+1] = ob[x+2] = 0;
				x += 3;
			}
			for (k = 0; k < ares; k++) {		/* Horizontal blocks */
				int zv = i * ares + k;
				if (zv >= res) {
					/* Fill res pixels with black */
					for (m = 0; m < res; m++) {
						ob[x] = ob[x+1] = ob[x+2] = 0;
						x += 3;
					}
				} else {
					xyz[2] = zv/(res-1.0);
					for (m = 0; m < res; m++) {	/* Horizontal in block */
						double tmp[3];
						xyz[0] = m/(res-1.0);

#ifndef NEVER	/* Straight XYZ */
						tmp[0] = xyz[0];
						tmp[1] = xyz[1];
						tmp[2] = xyz[2];
#else
						tmp[0] = 100.0 * xyz[0];
						tmp[1] = 256.0 * xyz[1] - 128.0;
						tmp[2] = 256.0 * xyz[2] - 128.0;
						Lab2XYZ(tmp, tmp);
#endif
						/* Convert XYZ through cam and back, then to RGB */
						cam1->XYZ_to_cam(cam1, Jab, tmp);
						cam2->cam_to_XYZ(cam2, rgb, Jab);
						XYZ2sRGB(rgb, rgb);
						/* Fill with pixel value */
						ob[x+0] = (int)(rgb[0] * 255.0 + 0.5);
						ob[x+1] = (int)(rgb[1] * 255.0 + 0.5);
						ob[x+2] = (int)(rgb[2] * 255.0 + 0.5);
						x += 3;
					}
				}
				/* Fill sp pixels with black */
				for (m = 0; m < sp; m++) {
					ob[x] = ob[x+1] = ob[x+2] = 0;
					x += 3;
				}
			}
//printf("~1 writing scan line at %d\n",y);
			TIFFWriteScanline(wh, ob, y, 0);
		}
		/* Fill sp lines with black */
		for (x = m = 0; m < w; m++) {
			ob[x] = ob[x+1] = ob[x+2] = 0;
			x += 3;
		}
		for (j = 0; j < sp; j++, y++) {
//printf("~1 writing blank scan line at %d\n",y);
			TIFFWriteScanline(wh, ob, y, 0);
		}
	}

	/* Write TIFF file */
	TIFFClose(wh);

	return 0;
}

