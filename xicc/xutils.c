
/* 
 * xicc standalone utilities
 *
 * Author:  Graeme W. Gill
 * Date:    2/7/00
 * Version: 1.00
 *
 * Copyright 2000 - 2006 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
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
#include "xutils.h"		/* definitions for this library */


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
























