
/* 
 * Argyll Gamut Mapping Library
 *
 * Author:  Graeme W. Gill
 * Date:    1/10/2000
 * Version: 2.00
 *
 * Copyright 2000 - 2005 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/* Gamut mapping object */
struct _gammap {

/* Private: */
	rspl *grey;		/* Incoming L map, NULL if none */
	rspl *map;		/* Lab -> Lab gamut map */

/* Public: */

	/* Methods */
	void (*del)(struct _gammap *s);			/* Free ourselves */
	void (*domap)(struct _gammap *s, double *out, double *in);	/* Do the mapping */

}; typedef struct _gammap gammap;

/* Creator */
gammap *new_gammap(
	int verb,			/* Verbose flag */
	gamut *sc_gam,		/* Source colorspace gamut */
	gamut *s_gam,		/* Source image gamut (NULL if none) */
	gamut *d_gam,		/* Destination colorspace gamut */
	double greymf,		/* Grey axis hue matching factor, 0.0 - 1.0 */
	double glumwcpf,	/* Grey axis White luminance compression factor, 0.0 - 1.0 */
	double glumwexf,	/* Grey axis White luminance expansion factor,   0.0 - 1.0 */
	double glumbcpf,	/* Grey axis Black luminance compression factor, 0.0 - 1.0 */
	double glumbexf,	/* Grey axis Black luminance expansion factor,   0.0 - 1.0 */
	double glumknf,		/* Grey axis luminance knee factor, 0.0 - 1.0 */
	double gamcpf,		/* Gamut compression factor, 0.0 - 1.0 */
	double gamexf,		/* Gamut expansion factor, 0.0 - 1.0 */
	double gamknf,		/* Gamut knee factor, 0.0 - 1.0 */
	double gampwf,		/* Gamut Perceptual Map weighting factor, 0.0 - 1.0 */
	double gamswf,		/* Gamut Saturation Map weighting factor, 0.0 - 1.0 */
	double satenh,		/* Saturation enhancement value, 0.0 - Inf */
	int    mapres,		/* Gamut map resolution, typically 9 - 33 */
	double *mn,			/* If not NULL, set minimum mapping input range */
	double *mx			/* for rspl grid */
);


	
