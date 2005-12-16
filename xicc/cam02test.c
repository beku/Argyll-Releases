
/* 
 * International Color Consortium color transform expanded support
 *
 * Author:  Graeme W. Gill
 * Date:    18/1/2004
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
 */


#include <stdio.h>
#include <math.h>
#include "xcam.h"
#include "cam02.h"

#define DIAG		/* Print internal value diagnostics for each spot test conversion */
					/* and print diagnostics for excessive errors, nans etc. */
#undef VERBOSE		/* Print diagnostic values for every conversion */
#define CHECKBIG	/* Check for bug numbers */
#undef SPOTTEST		/* Test known spot colors */
#define INVTEST		/* Lab cube to XYZ to Jab to XYZ */
#undef INVTEST1		/* Single value */
#define TESTINV		/* Jab cube to XYZ to Jab */
#undef TESTINV1		/* Single value */
#undef TESTINV2		/* J= 0 test values */
#define TRES 31		/* Grid resolution */
#define USE_HK 1	/* Use Helmholtz-Kohlraush */

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


/* Return maximum difference */
double maxdiff(double in1[3], double in2[3]) {
	int i;
	double tt, rv = 0.0;

	/* Deal with the possibility that we have nan's */
	for (i = 0; i < 3; i++) {
		tt = fabs(in1[i] - in2[i]);
		if (_isnan(tt))
			return tt;
		if (tt > rv)
			rv = tt;
	}
	return rv;
}

/* Return absolute difference */
double absdiff(double in1[3], double in2[3]) {
	double tt, rv = 0.0;
	tt = in1[0] - in2[0];
	rv += tt * tt;
	tt = in1[1] - in2[1];
	rv += tt * tt;
	tt = in1[2] - in2[2];
	rv += tt * tt;
	return sqrt(rv);
}

/* Return maximum Lab difference of XYZ */
double maxxyzdiff(double i1[3], double i2[3]) {
	int i;
	double tt, rv = 0.0;
	double in1[3], in2[3];

	XYZ2Lab(in1, i1);
	XYZ2Lab(in2, i2);

	/* Deal with the possibility that we have nan's */
	for (i = 0; i < 3; i++) {
		tt = fabs(in1[i] - in2[i]);
		if (_isnan(tt))
			return tt;
		if (tt > rv)
			rv = tt;
	}
	return rv;
}

int
main(void) {
	int ok = 1;
	double white[6][3] = {
		{ 0.9505, 1.000, 1.0888 },
		{ 0.9505, 1.000, 1.0888 },
		{ 1.0985, 1.000, 0.3558 },
		{ 1.0985, 1.000, 0.3558 },
		{ 0.9505, 1.000, 1.0888 },
		{ 0.9642, 1.000, 0.8249 }	/* D50 for inversion tests */
	};
	double La[6] = { 318.31, 31.83, 318.31, 31.83, 318.31, 150.0 };
	double sample[5][3] = {
		{ 0.1901, 0.2000, 0.2178 },
		{ 0.5706, 0.4306, 0.3196 },
		{ 0.0353, 0.0656, 0.0214 },
		{ 0.1901, 0.2000, 0.2178 },
		{ 0.9505, 1.000, 1.0888 }		/* Check white */
	};

#ifdef CIECAM02_HK
	double correct[5][4] = {			/* Hue range is last two */
		{ 41.75,  0.10,    219.0, 219.0 },
		{ 69.13,  48.57,   19.56, 19.56 },
		{ 30.22,  46.94,   177.1, 177.1 },
		{ 52.31,  51.92,   248.9, 248.9 },
		{ 100.00, 0.14,    0.0,   360.0 }
	};
#else
	double correct[5][4] = {			/* Hue range is last two */
		{ 41.73,  0.10,   219.0, 219.0 },
		{ 65.96,  48.57,  19.6,  19.6 },
		{ 21.79,  46.94,  177.1, 177.1 },
		{ 42.53,  51.92,  248.9, 248.9 },
		{ 100.0,  0.14,   0.0,   360.0 }		/* Check white */
	};
#endif

	double Jab[3];
	double res[3];
	cam02 *cam;
	int c;

	cam = new_cam02();

#ifdef SPOTTEST
	for (c = 0; c < 5; c++) {

#ifdef DIAG
		printf("Case %d:\n",c);
#endif /* DIAG */
		cam->set_view(
			cam,
			vc_average,	/* Enumerated Viewing Condition */
			white[c],	/* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) */
			0.20,		/* Relative Luminance of Background to reference white */
			La[c],		/* Adapting/Surround Luminance cd/m^2 */
			0.0,		/* Luminance of white in image - not used */
			0.00,		/* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
			white[c],	/* The Flare color coordinates (typically the Ambient color) */
			USE_HK		/* use Helmholtz-Kohlraush flag */ 
		);
#ifdef DIAG
		printf("\n");
#endif /* DIAG */
	
		{
			double JCh[3], target[3];
			cam->XYZ_to_cam(cam, Jab, sample[c]);
	
			/* Convert to JCh for checking */
			JCh[0] = Jab[0];

			/* Compute hue angle */
    		JCh[2] = (180.0/3.14159265359) * atan2(Jab[2], Jab[1]);
			JCh[2] = (JCh[2] < 0.0) ? JCh[2] + 360.0 : JCh[2];
	
			/* Compute chroma value */
			JCh[1] = sqrt(Jab[1] * Jab[1] + Jab[2] * Jab[2]);

			target[0] = correct[c][0];
			target[1] = correct[c][1];
			if (JCh[2] >= correct[c][2] && JCh[2] <= correct[c][3]) {
				target[2] = JCh[2];
			} else {
				if (fabs(JCh[2] - correct[c][2]) < fabs(JCh[2] - correct[c][3]))
					target[2] = correct[c][2];
				else
					target[2] = correct[c][3];
			}

			cam->cam_to_XYZ(cam, res, Jab);

			if (maxdiff(JCh, target) > 0.05) {
				printf("Spottest: Excesive error in conversion to CAM %f\n",maxdiff(JCh, target));
				ok = 0;
			}
			if (maxdiff(sample[c], res) > 1e-5) {
				printf("Spottest: Excessive error in inversion %f\n",maxdiff(sample[c], res));
				ok = 0;
			}
#ifdef DIAG
			printf("Jab is %f %f %f, Jch is %f %f %f\n",
			        Jab[0], Jab[1], Jab[2],
			        JCh[0], JCh[1], JCh[2]);

			printf("Error to expected value = %f\n",maxdiff(JCh, target));
		
			printf("XYZ is %f %f %f\n",res[0], res[1], res[2]);
			printf("Inversion error = %f\n",maxdiff(sample[c], res));
#endif /* DIAG */
		}
	}
#endif /* SPOTTEST */

#if defined(INVTEST) || defined(INVTEST1)
	for (c = 5; c < 6; c++) {
		/* Get the range of Jab values */
		double min[3] = { 1e38, 1e38, 1e38 };
		double max[3] = { -1e38, -1e38, -1e38 };

		cam->set_view(
			cam,
			vc_average,	/* Enumerated Viewing Condition */
			white[c],	/* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) = D50 */
			0.20,		/* Relative Luminance of Background to reference white */
			34.0,		/* Adapting/Surround Luminance cd/m^2 */
			0.0,		/* Luminance of white in image - not used */
			0.01,		/* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
			white[c],	/* The Flare color coordinates (typically the Ambient color) */
			USE_HK		/* use Helmholtz-Kohlraush flag */ 
		);
	
		{
			int co0, co1, co2;		/* (using co[3] triggers compiler bug) */
			double merr = 0.0;
			double xyz[3], Lab[3], Jab[3], in[3], checkxyz[3];
			double mxd;

#ifdef INVTEST1
			/* Test case */
			Lab[0] = 0.0;
			Lab[1] = -128.0;
			Lab[2] = 36.571429;
			Lab2XYZ(xyz, Lab);
			cam->XYZ_to_cam(cam, Jab, xyz);
			cam->cam_to_XYZ(cam, checkxyz, Jab);

			/* Check the result */
			mxd = maxxyzdiff(checkxyz, xyz);
			if (_finite(merr) && (_isnan(mxd) || mxd > merr))
				merr = mxd;
#ifdef DIAG
#ifndef VERBOSE
			if (_isnan(mxd) || mxd > 0.1)
#endif
			{
				printf("\n");
				printf("#### Lab = %f %f %f\n",Lab[0], Lab[1], Lab[2]);
				printf("%f %f %f -> %f %f %f -> %f %f %f [%f]\n",
				xyz[0],xyz[1],xyz[2],Jab[0],Jab[1],Jab[2],
				checkxyz[0],checkxyz[1],checkxyz[2],
				maxxyzdiff(checkxyz, xyz));
			}
#endif /* DIAG */

#else /* !INVTEST1 */
			/* itterate through Lab space */
			for (co0 = 0; co0 < TRES; co0++) {
				Lab[0] = co0/(TRES-1.0);
				Lab[0] = Lab[0] * 100.0;
				for (co1 = 0; co1 < TRES; co1++) {
					Lab[1] = co1/(TRES-1.0);
					Lab[1] = (Lab[1] - 0.5) * 256.0;
					for (co2 = 0; co2 < TRES; co2++) {
						int i;
						Lab[2] = co2/(TRES-1.0);
						Lab[2] = (Lab[2] - 0.5) * 256.0;
		
						Lab2XYZ(xyz, Lab);
						cam->XYZ_to_cam(cam, Jab, xyz);
						for (i = 0; i < 3; i++) {
							if (Jab[i] < min[i])
								min[i] = Jab[i];
							if (Jab[i] > max[i])
								max[i] = Jab[i];
						}
						cam->cam_to_XYZ(cam, checkxyz, Jab);

						/* Check the result */
						mxd = maxxyzdiff(checkxyz, xyz);
						if (_finite(merr) && (_isnan(mxd) || mxd > merr))
							merr = mxd;
#ifdef DIAG
#ifndef VERBOSE
						if (_isnan(mxd) || mxd > 0.1)
#endif
						{
							double oLab[3];

							XYZ2Lab(oLab, checkxyz);
							printf("\n");
							printf("#### Lab = %f %f %f -> %f %f %f\n",
							Lab[0], Lab[1], Lab[2], oLab[0], oLab[1], oLab[2]);
							printf("%f %f %f -> %f %f %f -> %f %f %f [%f]\n",
							xyz[0],xyz[1],xyz[2],Jab[0],Jab[1],Jab[2],
							checkxyz[0],checkxyz[1],checkxyz[2],
							maxxyzdiff(checkxyz, xyz));
						}
#endif /* DIAG */
					}
				}
			}
			if (_isnan(merr) || merr > 0.1) {
				printf("Excessive error in inversion check %f\n",merr);
				ok = 0;
			}
#ifdef DIAG
			printf("\n",merr);
			printf("Inversion check complete, peak error = %f\n",merr);
#endif /* !INVTEST1 */

#endif /* DIAG */
			printf("\n");
			printf("Range of Jab values was:\n");
			printf("J:  %f -> %f\n", min[0], max[0]);
			printf("a:  %f -> %f\n", min[1], max[1]);
			printf("b:  %f -> %f\n", min[2], max[2]);
		}
	}
#endif /* INVTEST */

#if defined(TESTINV) || defined(TESTINV1) || defined(TESTINV2)
	for (c = 5; c < 6; c++) {
		/* Get the range of XYZ values */
		double min[3] = { 1e38, 1e38, 1e38 };
		double max[3] = { -1e38, -1e38, -1e38 };

		cam->set_view(
			cam,
			vc_average,	/* Enumerated Viewing Condition */
			white[c],	/* Reference/Adapted White XYZ (Y range 0.0 .. 1.0) = D50 */
			0.20,		/* Relative Luminance of Background to reference white */
			34.0,		/* Adapting/Surround Luminance cd/m^2 */
			0.0,		/* Luminance of white in image - not used */
			0.01,		/* Flare as a fraction of the reference white (Y range 0.0 .. 1.0) */
			white[c],	/* The Flare color coordinates (typically the Ambient color) */
			USE_HK		/* use Helmholtz-Kohlraush flag */ 
		);
	
		{
			int co0, co1, co2;		/* (using co[3] triggers compiler bug) */
			double merr = 0.0;
			double xyz[3], Jab[3], in[3], checkJab[3];
			double mxd;

#if defined(TESTINV1) || defined(TESTINV2)

#ifdef TESTINV1
			/* SIngle sample test case */
			Jab[0] = 44.500000;
			Jab[1] = 32.914286;
			Jab[2] = -117.028571;
			cam->cam_to_XYZ(cam, xyz, Jab);
			for (i = 0; i < 3; i++) {
				if (xyz[i] < min[i])
					min[i] = xyz[i];
				if (xyz[i] > max[i])
					max[i] = xyz[i];
			}
			cam->XYZ_to_cam(cam, checkJab, xyz);

			/* Check the result */
			mxd = maxdiff(checkJab, Jab);
			if (_finite(merr) && (_isnan(mxd) || mxd > merr))
				merr = mxd;
#ifdef DIAG
			if (_isnan(mxd) || mxd > 0.1)
#ifndef VERBOSE
			{
#endif
				printf("\n");
				printf("#### Jab = %f %f %f\n",Jab[0], Jab[1], Jab[2]);
				printf("%f %f %f -> %f %f %f -> %f %f %f [%f]\n",
				Jab[0],Jab[1],Jab[2], xyz[0],xyz[1],xyz[2],
				checkJab[0],checkJab[1],checkJab[2],
				maxdiff(checkJab, Jab));
			}
#endif /* DIAG */
#else	/* !TESTINV1 */
			{	/* J = 0 test cases */
				int i;
				for (i = 1; i < 128; i++) {
					Jab[0] = 0.0;
					Jab[1] = (double)i/2;
					Jab[2] = (double)i;
					cam->cam_to_XYZ(cam, xyz, Jab);

					Jab[0] = 0.0;
					Jab[1] = (double)-i/2;
					Jab[2] = (double)i;
					cam->cam_to_XYZ(cam, xyz, Jab);

					Jab[0] = 0.0;
					Jab[1] = (double)i/2;
					Jab[2] = (double)-i;
					cam->cam_to_XYZ(cam, xyz, Jab);

					Jab[0] = 0.0;
					Jab[1] = (double)-i/2;
					Jab[2] = (double)-i;
					cam->cam_to_XYZ(cam, xyz, Jab);

					Jab[0] = 0.0;
					Jab[1] = (double)i;
					Jab[2] = (double)i/2;
					cam->cam_to_XYZ(cam, xyz, Jab);

					Jab[0] = 0.0;
					Jab[1] = (double)-i;
					Jab[2] = (double)i/2;
					cam->cam_to_XYZ(cam, xyz, Jab);

					Jab[0] = 0.0;
					Jab[1] = (double)i;
					Jab[2] = (double)-i/2;
					cam->cam_to_XYZ(cam, xyz, Jab);

					Jab[0] = 0.0;
					Jab[1] = (double)-i;
					Jab[2] = (double)-i/2;
					cam->cam_to_XYZ(cam, xyz, Jab);
				}
			}
#endif	/* !TESTINV1 */

#else /* !TESTINV1 && !TESTINV2 */

			/* itterate through Jab space */
			for (co0 = 0; co0 < TRES; co0++) {
				Jab[0] = co0/(TRES-1.0);
				Jab[0] = -5.0 + Jab[0] * 105.0;	
				for (co1 = 0; co1 < TRES; co1++) {
					Jab[1] = co1/(TRES-1.0);
					Jab[1] = (Jab[1] - 0.5) * 256.0;
					for (co2 = 0; co2 < TRES; co2++) {
						int i;
						Jab[2] = co2/(TRES-1.0);
						Jab[2] = (Jab[2] - 0.5) * 256.0;
		
						cam->cam_to_XYZ(cam, xyz, Jab);
						for (i = 0; i < 3; i++) {
							if (xyz[i] < min[i])
								min[i] = xyz[i];
							if (xyz[i] > max[i])
								max[i] = xyz[i];
						}
						cam->XYZ_to_cam(cam, checkJab, xyz);

						/* Check the result */
						mxd = maxdiff(checkJab, Jab);
						if (_finite(merr) && (_isnan(mxd) || mxd > merr))
							merr = mxd;
#ifdef DIAG
#ifndef VERBOSE
						if (_isnan(mxd) || mxd > 0.1)
#endif
						{
							printf("\n");
							printf("#### Jab = %f %f %f\n",Jab[0], Jab[1], Jab[2]);
							printf("%f %f %f -> %f %f %f -> %f %f %f [%f]\n",
							Jab[0],Jab[1],Jab[2], xyz[0],xyz[1],xyz[2],
							checkJab[0],checkJab[1],checkJab[2],
							maxdiff(checkJab, Jab));
						}
#endif /* DIAG */
					}
				}
			}
			if (_isnan(merr) || merr > 1.0) {
				printf("Excessive error in inversion check 2 %f\n",merr);
				ok = 0;
			}
#ifdef DIAG
			printf("Inversion check 2 complete, peak error = %f\n",merr);
#endif /* DIAG */
			printf("\n");
			printf("Range of XYZ values was:\n");
			printf("X:  %f -> %f\n", min[0], max[0]);
			printf("Y:  %f -> %f\n", min[1], max[1]);
			printf("Z:  %f -> %f\n", min[2], max[2]);
#endif /* !TESTINV1 && !TESTINV2 */
		}
	}
#endif /* TESTINV */

	printf("\n");
	if (ok == 0) {
		printf("Cam testing FAILED\n");
	} else {
		printf("Cam testing OK\n");
	}

	return 0;
}

