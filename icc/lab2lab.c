
/* 
 * International Color Consortium Format Library (icclib)
 * Creat a unity Lab to Lab input profile.
 *
 * Author:  Graeme W. Gill
 * Date:    1999/11/29
 * Version: 2.05
 *
 * Copyright 1999 - 2005 Graeme W. Gill
 *
 * This material is licensed with a free use license:-
 * see the License.txt file in this directory for licensing details.
 */

/* TTBD:
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <time.h>
#if defined(__IBMC__) && defined(_M_IX86)
#include <float.h>
#endif
#include "icc.h"

void error(char *fmt, ...), warning(char *fmt, ...);

/* - - - - - - - - - - - - - */
/* Some support functions */

/* Clip value to range 0.0 to 1.0 */
void clip(double inout[3]) {
	int i;
	for (i = 0; i < 3; i++) {
		if (inout[i] < 0.0)
			inout[i] = 0.0;
		else if (inout[i] > 1.0)
			inout[i] = 1.0;
	}
}

/* Protected power function */
double ppow(double num, double p) {
	if (num < 0.0)
		return -pow(-num, p);
	else
		return pow(num, p);
}

/* - - - - - - - - - - - - - */
/* The overall model is: */
/* Lab -> Lab' -> Lab' -> Lab */

/* Lab -> Lab' */
void Lab_Labp(void *cntx, double out[3], double in[3]) {
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];
}

/* Lab' -> Lab */
/* (We are using linear) */
void Labp_Labp(void *cntx, double out[3], double in[3]) {
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];
}

/* Lab' -> Lab */
void Labp_Lab(void *cntx, double out[3], double in[3]) {
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in[2];
}

/* - - - - - - - - - - - - - */

int
main(
	int argc,
	char *argv[]
) {
	char *file_name = "lab2lab.icm";
	FILE *wr_fp;
	icc *wr_icco;		/* Keep object separate */
	int rv = 0;

	/* ---------------------------------------- */
	/* Create a Lut16 based Lab profile   */
	/* ---------------------------------------- */

	/* Open up the file for writing */
	if ((wr_fp = fopen(file_name,"w")) == NULL)
		error ("Write: Can't open file '%s'",file_name);

#if defined(O_BINARY)
	wr_fp = freopen(file_name,"wb",wr_fp);
#endif
	
	if ((wr_icco = new_icc()) == NULL)
		error ("Write: Creation of ICC object failed");

	/* Add all the tags required */

	/* The header: */
	{
		icmHeader *wh = wr_icco->header;

		/* Values that must be set before writing */
		wh->deviceClass     = icSigInputClass;
    	wh->colorSpace      = icSigLabData;
    	wh->pcs             = icSigLabData;
    	wh->renderingIntent = icAbsoluteColorimetric;	/* For want of something */

		/* Values that should be set before writing */
		wh->manufacturer = str2tag("argl");
    	wh->model        = str2tag("    ");
	}
	/* Profile Description Tag: */
	{
		icmTextDescription *wo;
		char *dst = "A unity Lab to Lab transform";
		if ((wo = (icmTextDescription *)wr_icco->add_tag(
		           wr_icco, icSigProfileDescriptionTag,	icSigTextDescriptionType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, dst);		/* Copy the string in */
	}
	/* Copyright Tag: */
	{
		icmText *wo;
		char *crt = "Not Copyright";
		if ((wo = (icmText *)wr_icco->add_tag(
		           wr_icco, icSigCopyrightTag,	icSigTextType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(crt)+1; 	/* Allocated and used size of text, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->data, crt);		/* Copy the text in */
	}
	/* White Point Tag: */
	{
		icmXYZArray *wo;
		/* Note that tag types icSigXYZType and icSigXYZArrayType are identical */
		if ((wo = (icmXYZArray *)wr_icco->add_tag(
		           wr_icco, icSigMediaWhitePointTag, icSigXYZArrayType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = 1;
		wo->allocate((icmBase *)wo);	/* Allocate space */
		wo->data[0].X = 0.9642;	/* D50 white point */
		wo->data[0].Y = 1.0000;
		wo->data[0].Z = 0.8249;
	}
	/* 16 bit dev -> pcs lut: */
	{
		double Labmin[3];
		double Labmax[3];
		icmLut *wo;

		/* Intent 1 = perceptual colorimetric */
		if ((wo = (icmLut *)wr_icco->add_tag(
		           wr_icco, icSigAToB0Tag,	icSigLut16Type)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->inputChan = 3;
		wo->outputChan = 3;
    	wo->clutPoints = 2;
    	wo->inputEnt = 2;
    	wo->outputEnt = 2;
		wo->allocate((icmBase *)wo);/* Allocate space */

		/* Don't try for symetry etc. */
		Labmin[0] = 0.0 * (100.0 * 65535.0)/65280.0;			/* L */
		Labmin[1] = (0.0 * (255.0 * 65535.0)/65280) - 128.0;	/* a */
		Labmin[2] = (0.0 * (255.0 * 65535.0)/65280) - 128.0;	/* b */
		Labmax[0] = 1.0 * (100.0 * 65535.0)/65280.0;			/* L */
		Labmax[1] = (1.0 * (255.0 * 65535.0)/65280) - 128.0;	/* a */
		Labmax[2] = (1.0 * (255.0 * 65535.0)/65280) - 128.0;	/* b */

		/* The matrix is only applicable to XYZ input space, */
		/* so it is not used here. */

		/* Use helper function to do the hard work. */
		if (wo->set_tables(wo, NULL,
				icSigLabData, 				/* Input color space */
				icSigLabData, 				/* Output color space */
				Lab_Labp,					/* Input transfer function, Lab->Lab' */
				NULL, NULL,				/* Use default Maximum range of Lab' values */
				Labp_Labp,					/* Lab' -> Lab' transfer function */
				NULL, NULL,				/* Use default Maximum range of Lab' values */
				Labp_Lab					/* Linear output transform Lab'->Lab */
		) != 0)
			error("Setting 16 bit Lab->Lab Lut failed: %d, %s",wr_icco->errc,wr_icco->err);
	}

	/* Write the file out */
	if ((rv = wr_icco->write(wr_icco,wr_fp,0)) != 0)
		error ("Write file: %d, %s",rv,wr_icco->err);
	
	wr_icco->free(wr_icco);
	fclose(wr_fp);

	printf("Profile creation completed OK\n");
	return 0;
}

/* ------------------------------------------------ */
/* Basic printf type error() and warning() routines */

void
error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"lutest: Error - ");
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

	fprintf(stderr,"lutest: Warning - ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}
