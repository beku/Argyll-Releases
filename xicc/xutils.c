
/* 
 * xicc standalone utilities
 *
 * Author:  Graeme W. Gill
 * Date:    2/7/00
 * Version: 1.00
 *
 * Copyright 2000 - 2006 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU AFFERO GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 *
 */

/*
 * This module provides expanded capabilities,
 * but is independent of other modules.
 */

#include <sys/types.h>
#include <string.h>
#include <ctype.h>
#ifdef __sun
#include <unistd.h>
#endif
#if defined(__IBMC__) && defined(_M_IX86)
#include <float.h>
#endif
#include "copyright.h"
#include "config.h"
#include "icc.h"
#include "tiffio.h"
#include "xutils.h"		/* definitions for this library */


#undef DEBUG

#ifdef DEBUG
# define errout stderr
# define debug(xx)	fprintf(errout, xx )
# define debug2(xx)	fprintf xx
#else
# define debug(xx)
# define debug2(xx)
#endif

/* ------------------------------------------------------ */
/* Common clut table code */

/* Default table of clut resolutions */
/* See discussion in imdi/imdi_gen.c for ideal numbers */
static int lut_resolutions[9][4] = {
	/* low, med, high, vhigh */
	{ 0,     0,    0,    0 },		/* 0 */
	{ 256, 772, 4370, 4370 },		/* 1 */
	{  86, 256,  256,  256 },		/* 2 */
	{   9,  17,   33,   52 },		/* 3 */
	{   6,   9,   18,   33 },		/* 4 */
	{   6,   9,   16,   18 },		/* 5 */
	{   6,   6,    9,   12 },		/* 6 */
	{   6,   7,    7,    9 },		/* 7 */
	{   3,   5,    5,    7 }		/* 8 */
};


/* return a lut resolution given the input dimesion and quality */
/* Input dimension [0-8], quality: low 0, medium 1, high 2, very high 3 . */
/* A returned value of 0 indicates illegal.  */
int dim_to_clutres(int dim, int quality) {
	if (dim < 0)
		dim = 0;
	else if (dim > 8)
		dim = 8;
	if (quality < 0)
		quality = 0;
	if (quality > 3)
		quality = 3;
	return lut_resolutions[dim][quality];
}

/* ------------------------------------------------------ */

/* Open an ICC file or an TIFF file with an embeded ICC profile for reading. */
/* Return NULL on error */
icc *read_embeded_icc(char *file_name) {
	TIFF *rh = NULL;
	int  size;
	void *tag, *buf;
	icmAlloc *al;
	icmFile *fp;
	icc *icco;
	TIFFErrorHandler oldhandler;
	int rv;

	/* First see if the file can be opened as an ICC profile */
	if ((fp = new_icmFileStd_name(file_name,"r")) == NULL) {
		debug2((errout,"Can't open file '%s'\n",file_name));
		return NULL;
	}

	if ((icco = new_icc()) == NULL) {
		debug("Creation of ICC object failed\n");
		fp->del(fp);
		return NULL;
	}

	if ((rv = icco->read_x(icco,fp,0,1)) == 0) {
		debug2((errout,"Opened '%s' as an icc profile\n",file_name));
		return icco;
	}

	debug2((errout,"icc read failed with %d, %s\n",rv,icco->err));
	icco->del(icco);		/* icc wil fp->del() */

	/* Not an ICC profile, see if it's a TIFF file */
	oldhandler = TIFFSetWarningHandler(NULL);

	if ((rh = TIFFOpen(file_name, "r")) == NULL) {
		debug2((errout,"TIFFOpen failed for '%s'\n",file_name));
		TIFFSetWarningHandler(oldhandler);
		return NULL;
	}

	if (TIFFGetField(rh, TIFFTAG_ICCPROFILE, &size, &tag) == 0 || size == 0) {
		TIFFClose(rh);
		TIFFSetWarningHandler(oldhandler);
		return NULL;
	}

	/* Make a copy of the profile to a memory buffer */
	if ((al = new_icmAllocStd()) == NULL) {
		debug("new_icmAllocStd failed\n");
		TIFFClose(rh);
		TIFFSetWarningHandler(oldhandler);
	    return NULL;
	}
	if ((buf = al->malloc(al, size)) == NULL) {
		debug("malloc of profile buffer failed\n");
		al->del(al);
		TIFFClose(rh);
		TIFFSetWarningHandler(oldhandler);
	    return NULL;
	}

	memcpy(buf, tag, size);
	TIFFClose(rh);
	TIFFSetWarningHandler(oldhandler);

	/* Memory File fp that will free the buffer when deleted: */
	if ((fp = new_icmFileMem_ad(buf, size, al)) == NULL) {
		debug("Creating memory file from CMProfileLocation failed");
		al->free(al, buf);
		al->del(al);
		return NULL;
	}

	if ((icco = new_icc()) == NULL) {
		debug("Creation of ICC object failed\n");
		fp->del(fp);	/* fp will delete al */
		return NULL;
	}

	if ((rv = icco->read_x(icco,fp,0,1)) == 0) {
		debug2((errout,"Opened '%s' embeded icc profile\n",file_name));
		return icco;
	}

	debug2((errout,"Failed to read '%s' embeded icc profile\n",file_name));
	icco->del(icco);	/* icco will delete fp and al */
	return NULL;
}

/* ------------------------------------------------------ */
























