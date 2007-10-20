
/* 
 * Argyll Color Correction System
 *
 * Scanin: Input the scan of a test chart, and output cgats data
 *         Uses scanrd to do the hard work.
 *
 * Author: Graeme W. Gill
 * Date:   29/1/97
 *
 * Copyright 1995 - 2002 Graeme W. Gill
 * All rights reserved.
 * This material is licenced under the GNU GENERAL PUBLIC LICENSE Version 3 :-
 * see the License.txt file for licencing details.
 */

#define VERSION "1.5"

#include <stdio.h>
#include <fcntl.h>		/* In case DOS binary stuff is needed */
#include <ctype.h>
#include <string.h>
#include <math.h>

#include <time.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "cgats.h"
#include "scanrd.h"
#include <numsup.h>
#include <tiffio.h>
#include <icc.h>
#include <xcolorants.h>

void fix_it8(char *o, char *i);

#define TXBUF (256*1024L)

/* NOTE: We aren't handling the libtiff error/warning messages !! */

/* Read a line of data from the input Grey, RGB or CMYK tiff file */
/* return non-zero on error */
int read_line(
void *fdata,
int y,
char *dst
) {
	if (TIFFReadScanline((TIFF *)fdata, (tdata_t)dst, y, 0) < 0)
		return 1;
	return 0;
}

/* Write a line of data to the diagnostic RGB tiff file */
/* return non-zero on error */
static int
write_line(
void *ddata,
int y,
char *src
) {
	if (TIFFWriteScanline((TIFF *)ddata, (tdata_t)src, y, 0) < 0)
		return 1;
	return 0;
}

void
usage(void) {
	fprintf(stderr,"Scanin, Version %s\n",VERSION);
	fprintf(stderr,"Author: Graeme W. Gill, licensed under the GPL Version 3\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"usage: scanin [options] input.tif recogin.cht valin.cie [diag.tif]\n");
	fprintf(stderr,"   :- inputs 'input.tif' and outputs scanner 'input.ti3', or\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"usage: scanin -g [options] input.tif recogout.cht [diag.tif]\n");
	fprintf(stderr,"   :- outputs file 'recogout.cht', or\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"usage: scanin -o [options] input.tif recogin.cht [diag.tif]\n");
	fprintf(stderr,"   :- outputs file 'input.val', or\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"usage: scanin -c [options] input.tif recogin.cht scanprofile.icc pbase [diag.tif]\n");
	fprintf(stderr,"   :- inputs pbase.ti2 and outputs printer pbase.ti3, or\n");
	fprintf(stderr,"\n");
	fprintf(stderr,"usage: scanin -r [options] input.tif recogin.cht pbase [diag.tif]\n");
	fprintf(stderr,"   :- inputs pbase.ti2+.ti3 and outputs pbase.ti3\n");
	fprintf(stderr,"\n");
	fprintf(stderr," -g                   Generate a chart reference (.cht) file\n");
	fprintf(stderr," -o                   Output patch values in .val file\n");
	fprintf(stderr," -c                   Use scanner as colorimeter to\n");
	fprintf(stderr,"                       convert printer .ti2 to .ti3\n");
	fprintf(stderr," -ca                  Same as -c, but accumulates more values to .ti3\n");
	fprintf(stderr,"                       from subsequent pages\n");
	fprintf(stderr," -r                   Replace device values in .ti2/.ti3\n");
	fprintf(stderr,"                      Default is to create a scanner .ti3 file\n");
	fprintf(stderr," -F x1,y1,x2,y2,x3,y3 Don't auto recognize, locate using three fiducual marks\n");
	fprintf(stderr," -a                   Recognise chart in normal orientation only (-A fallback as is)\n");
	fprintf(stderr,"                      Default is to recognise all possible chart angles\n");
	fprintf(stderr," -m                   Return true mean (default is robust mean)\n");
	fprintf(stderr," -G gamma             Approximate gamma encoding of image\n");
	fprintf(stderr," -v [n]               Verbosity level 0-9\n");
	fprintf(stderr," -d [ihvglLrsonap]    Generate diagnostic output (try -dipn)\n");
	fprintf(stderr,"     i                diag - B&W of input image\n");
	fprintf(stderr,"     h                diag - Horizontal edge/tick detection\n");
	fprintf(stderr,"     v                diag - Vertical edge/tick detection\n");
	fprintf(stderr,"     g                diag - Groups detected\n");
	fprintf(stderr,"     l                diag - Lines detected\n");
	fprintf(stderr,"     L                diag - All lines detected\n");
	fprintf(stderr,"     r                diag - lines rotated\n");
	fprintf(stderr,"     s                diag - diagnostic sample boxes rotated\n");
	fprintf(stderr,"     o                diag - sample box outlines\n");
	fprintf(stderr,"     n                diag - sample box names\n");
	fprintf(stderr,"     a                diag - sample box areas\n");
	fprintf(stderr,"     p                diag - pixel areas sampled\n");
	exit(1);
	}


int main(int argc, char *argv[])
{
	int fa,nfa;					/* current argument we're looking at */
	static char tiffin_name[200] = { 0 };	/* TIFF Input file name (.tif) */
	static char datin_name[200] = { 0 };	/* Data input name (.cie/.q60) */
	static char datout_name[200] = { 0 };	/* Data output name (.ti3/.val) */
	static char recog_name[200] = { 0 };	/* Reference chart name (.cht) */
	static char prof_name[200] = { 0 };		/* scanner profile name (.cht) */
	static char diag_name[200] = { 0 };		/* Diagnostic Output (.tif) name, if used */
	int verb = 1;
	int tmean = 0;		/* Return true mean, rather than robust mean */
	int repl = 0;		/* Replace .ti3 device values from raster file */
	int outo = 0;		/* Output the values read, rather than creating scanner .ti3 */
	int colm = 0;		/* Use scan values as colorimter for print profile. > 1 == append */
	int flags = SI_GENERAL_ROT;	/* Default allow all rotations */

	TIFF *rh = NULL, *wh = NULL;
	uint16 depth, bps;
	uint16 pconfig, photometric;
	uint16 resunits;
	float resx, resy;

	icColorSpaceSignature tiffs = 0;	/* Type of tiff color space */

	int i, j;
	double gamma = 0.0;		/* default */
	double _sfid[6], *sfid = NULL;		/* Specified fiducials */
	int width, height;		/* x and y size */

	scanrd *sr;				/* Scanrd object */
	int err;	
	char *errm;

	if (argc <= 1)
		usage();

	error_program = argv[0];

	/* Process the arguments */
	for(fa = 1;fa < argc;fa++) {
		nfa = fa;					/* skip to nfa if next argument is used */
		if (argv[fa][0] == '-')		/* Look for any flags */
			{
			char *na = NULL;		/* next argument after flag, null if none */

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		/* next is directly after flag */
			else
				{
				if ((fa+1) < argc)
					{
					if (argv[fa+1][0] != '-')
						{
						nfa = fa + 1;
						na = argv[nfa];		/* next is seperate non-flag argument */
						}
					}
				}

			if (argv[fa][1] == '?') {
				usage();
			} else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 2;
				if (na != NULL && isdigit(na[0])) {
					verb = atoi(na);
				}
			} else if (argv[fa][1] == 'm' || argv[fa][1] == 'M') {
				tmean = 1;

			} else if (argv[fa][1] == 'g') {
				flags |= SI_BUILD_REF;
				repl = 0;
				outo = 0;
				colm = 0;

			} else if (argv[fa][1] == 'r' || argv[fa][1] == 'R') {
				repl = 1;
				outo = 0;
				colm = 0;

			} else if (argv[fa][1] == 'o' || argv[fa][1] == 'O') {
				repl = 0;
				outo = 1;
				colm = 0;

			} else if (argv[fa][1] == 'c' || argv[fa][1] == 'C') {
				repl = 0;
				outo = 0;
				colm = 1;
				if (na != NULL && (*na == 'a' || *na == 'A'))
					colm = 2;

			/* Approximate gamma encoding of image */
			} else if (argv[fa][1] == 'G') {
				fa = nfa;
				if (na == NULL) usage();
				gamma = atof(na);
				if (gamma < 0.0 || gamma > 5.0)
					usage();

			/* Use specified fiducials instead of auto recognition */
			} else if (argv[fa][1] == 'F') {
				fa = nfa;
				if (na == NULL) usage();
				if (sscanf(na, " %lf,%lf,%lf,%lf,%lf,%lf ", &_sfid[0], &_sfid[1], &_sfid[2], &_sfid[3], &_sfid[4], &_sfid[5]) != 6) {
					usage();
				}

				sfid = _sfid;

			/* Don't recognise rotations */
			} else if (argv[fa][1] == 'a') {
				flags &= ~SI_GENERAL_ROT;

			/* Don't recognise rotations, and read patches */
			/* anyway "as is", if everything else failes */
			} else if (argv[fa][1] == 'A') {
				flags &= ~SI_GENERAL_ROT;
				flags |= SI_ASISIFFAIL;

			} else if (argv[fa][1] == 'd' || argv[fa][1] == 'D') {
				while (na != NULL && *na != '\000') {
					switch (*na) {
						case 'i':
							flags |= SI_SHOW_IMAGE;
							break;
						case 'h':
							flags |= SI_SHOW_DIFFSH;
							break;
						case 'v':
							flags |= SI_SHOW_DIFFSV;
							break;
						case 'g':
							flags |= SI_SHOW_GROUPS;
							break;
						case 'l':
							flags |= SI_SHOW_LINES;
							break;
						case 'L':
							flags |= SI_SHOW_ALL_LINES;
							break;
						case 'r':
							flags |= SI_SHOW_ROT;
							break;
						case 's':
							flags |= SI_SHOW_SBOX;
							break;
						case 'o':
							flags |= SI_SHOW_SBOX_OUTLINES;
							break;
						case 'n':
							flags |= SI_SHOW_SBOX_NAMES;
							break;
						case 'a':
							flags |= SI_SHOW_SBOX_AREAS;
							break;
						case 'p':
							flags |= SI_SHOW_SAMPLED_AREA;
							break;
						default:
							usage();
					}
					na++;
				}
			} else 
				usage();
		} else
			break;
	}

	/* TIFF Raster input file name */
	if (fa >= argc || argv[fa][0] == '-') usage();
	strcpy(tiffin_name,argv[fa]);

	if ((flags & SI_BUILD_REF) == 0
	     && repl == 0 && colm == 0) {	/* Not generate ref or replacing .ti3 dev */
		char *xl;
		strcpy(datout_name,argv[fa]);
		if ((xl = strrchr(datout_name, '.')) == NULL)	/* Figure where extention is */
			xl = datout_name + strlen(datout_name);
		if (outo == 0)	/* Creating scan calib data */
			strcpy(xl,".ti3");
		else			/* Just outputing values for some other purpose */
			strcpy(xl,".val");
	}

	/* .cht Reference file in or out */
	if (++fa >= argc || argv[fa][0] == '-') usage();
	strcpy(recog_name,argv[fa]);

	if (colm > 0) {
		if (++fa >= argc || argv[fa][0] == '-') usage();
		strcpy(prof_name,argv[fa]);
	}

	/* CGATS Data file input/output */
	if ((flags & SI_BUILD_REF) == 0 && outo == 0) {	/* Not generate ref or just outputing */
		if (++fa >= argc || argv[fa][0] == '-') usage();
		if (outo == 0) {	/* Creating scan calib data */
			/* Data file */
			strcpy(datin_name,argv[fa]);
		}
		if (repl != 0 || colm > 0) {	/* Colorimter emulation or replacing .ti3 device data */
			strcpy(datin_name,argv[fa]);
			strcat(datin_name,".ti2");
			strcpy(datout_name,argv[fa]);
			strcat(datout_name,".ti3");
		}
	}

	/* optional diagnostic file */
	if (++fa < argc) {
		if (argv[fa][0] == '-')
			usage();
		strcpy(diag_name,argv[fa]);
	} else {	/* Provide a default name */
		strcpy(diag_name,"diag.tif");
	}

	/* ----------------------------------------- */
	/* Open up input tiff file ready for reading */
	/* Got arguments, so setup to process the file */
	if ((rh = TIFFOpen(tiffin_name, "r")) == NULL)
		error("error opening read file '%s'",tiffin_name);

	TIFFGetField(rh, TIFFTAG_IMAGEWIDTH,  &width);
	TIFFGetField(rh, TIFFTAG_IMAGELENGTH, &height);

	TIFFGetField(rh, TIFFTAG_BITSPERSAMPLE, &bps);
	if (bps != 8 && bps != 16)
		error("TIFF Input file must be 8 or 16 bits/channel");

	TIFFGetField(rh, TIFFTAG_SAMPLESPERPIXEL, &depth);
	if (depth != 1 && depth != 3 && depth != 4)
		error("Input must be a Grey, RGB or CMYK tiff file");

	TIFFGetField(rh, TIFFTAG_PHOTOMETRIC, &photometric);
	if (depth == 1 && photometric != PHOTOMETRIC_MINISBLACK
	               && photometric != PHOTOMETRIC_MINISWHITE)
		error("1 chanel input must be a Grey tiff file");
	else if (depth == 3 && photometric != PHOTOMETRIC_RGB)
		error("3 chanel input must be an RGB tiff file");
	else if (depth == 4 && photometric != PHOTOMETRIC_SEPARATED)
		error("4 chanel input must be a CMYK tiff file");

	if (depth == 1)
		tiffs = icSigGrayData;
	else if (depth == 3) 
		tiffs = icSigRgbData;
	else if (depth == 4)
		tiffs = icSigCmykData;

	TIFFGetField(rh, TIFFTAG_PLANARCONFIG, &pconfig);
	if (pconfig != PLANARCONFIG_CONTIG)
		error ("TIFF Input file must be planar");

	TIFFGetField(rh, TIFFTAG_RESOLUTIONUNIT, &resunits);
	TIFFGetField(rh, TIFFTAG_XRESOLUTION, &resx);
	TIFFGetField(rh, TIFFTAG_YRESOLUTION, &resy);

	/* -------------------------- */
	/* setup the diag output file */
	if (flags & SI_SHOW_FLAGS) {
		if ((wh = TIFFOpen(diag_name, "w")) == NULL)
			error("Can\'t create TIFF file '%s'!",diag_name);
	
		TIFFSetField(wh, TIFFTAG_IMAGEWIDTH,  width);
		TIFFSetField(wh, TIFFTAG_IMAGELENGTH, height);
		TIFFSetField(wh, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
		TIFFSetField(wh, TIFFTAG_SAMPLESPERPIXEL, 3);
		TIFFSetField(wh, TIFFTAG_BITSPERSAMPLE, 8);
		TIFFSetField(wh, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
		TIFFSetField(wh, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
		TIFFSetField(wh, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
		TIFFSetField(wh, TIFFTAG_RESOLUTIONUNIT, resunits);
		TIFFSetField(wh, TIFFTAG_XRESOLUTION, resx);
		TIFFSetField(wh, TIFFTAG_YRESOLUTION, resy);
		TIFFSetField(wh, TIFFTAG_IMAGEDESCRIPTION, "Scanin diagnosis output");
	}

	/* -------------------------- */
	if (verb >= 2) {
		printf("Input file '%s': w=%d, h=%d, d = %d, bpp = %d\n",
		        tiffin_name, width, height, depth, bps);
		if (flags & SI_BUILD_REF)
			printf("Build Scan Chart reference file '%s'\n",recog_name);
		else {
			printf("Data input file '%s'\n",datin_name);
			printf("Data output file '%s'\n",datout_name);
			printf("Chart reference file '%s'\n",recog_name);
		}
		if (flags & SI_SHOW_FLAGS)
			printf("Creating diagnostic tiff file '%s'\n",diag_name);
	}

	/* -------------------------- */
	/* Do the operation */

	if ((sr = do_scanrd(
		flags,			/* option flags */
		verb,			/* verbosity level */

		gamma,
		sfid,			/* Specified fiducuals, if any */
		width, height, depth, bps,	/* Width, Height and Depth of input in pixels */
		read_line,		/* Read line function */
		(void *)rh,		/* Opaque data for read_line */

		recog_name,		/* reference file name */

		write_line,		/* Write line function */
		(void *)wh		/* Opaque data for write_line */
	)) == NULL) {
		if (flags & SI_SHOW_FLAGS)
			TIFFClose(wh);
		error ("Unable to allocate scanrd object");
	}

	if ((err = sr->error(sr, &errm)) != 0) {
		if ((flags & SI_SHOW_FLAGS) && err != SI_DIAG_WRITE_ERR)
			TIFFClose(wh);		/* Close diagnostic file */
		error ("Code 0x%x, %s",err,errm);
	}

	/* Read an output the values */
	if ((flags & SI_BUILD_REF) == 0) {	/* Not generate ref */

		/* -------------------------------------------------- */
		if (outo != 0) {		/* Just output the values */
								/* Note value range is raw 0..255, */
								/* while all others output formats are out of 100 */
			cgats *ocg;			/* output cgats structure */
			time_t clk = time(0);
			struct tm *tsp = localtime(&clk);
			char *atm = asctime(tsp); /* Ascii time */
	
			/* Setup output cgats file */
			ocg = new_cgats();	/* Create a CGATS structure */
			ocg->add_other(ocg, "VALS"); 	/* Dummy type */
			ocg->add_table(ocg, tt_other, 0);	/* Start the first table */
	
			ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration raster values",NULL);
			ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll scanin", NULL);
			atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
			ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	
			ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);
			if (depth == 1) {
				ocg->add_field(ocg, 0, "GREY", r_t);
			} else if (depth == 3) {
				ocg->add_field(ocg, 0, "RGB_R", r_t);
				ocg->add_field(ocg, 0, "RGB_G", r_t);
				ocg->add_field(ocg, 0, "RGB_B", r_t);
			} else if (depth == 4) {
				ocg->add_field(ocg, 0, "CMYK_C", r_t);
				ocg->add_field(ocg, 0, "CMYK_M", r_t);
				ocg->add_field(ocg, 0, "CMYK_Y", r_t);
				ocg->add_field(ocg, 0, "CMYK_K", r_t);
			}
	
			/* Initialise, ready to read out all the values */
			for (j = 0; ; j++) {
				char id[100];		/* Input patch id */
				double P[4];		/* Robust/true mean values */

				if (tmean) {
					if (sr->read(sr, id, NULL, P, NULL, NULL) != 0)
						break;
				} else {
					if (sr->read(sr, id, P, NULL, NULL, NULL) != 0)
						break;
				}
		
				if (depth == 1) {
					ocg->add_set( ocg, 0, id, P[0]);
				} else if (depth == 3) {
					ocg->add_set( ocg, 0, id, P[0], P[1], P[2]);
				} else if (depth == 4) {
					ocg->add_set( ocg, 0, id, P[0], P[1], P[2], P[3]);
				}
			}
	
			if (ocg->write_name(ocg, datout_name))
				error("Write error : %s",ocg->err);
	
			ocg->del(ocg);		/* Clean up */

		/* -------------------------------------------------- */
		} else if (repl != 0) {	/* Replace .ti3 device values */
			cgats *icg;			/* input .ti2 cgats structure */
			cgats *ocg;			/* input/output .ti3 cgats structure */
			int npat;			/* Number of test patches */
			int dim = 0;		/* Dimenstionality of device space */
			int fi;				/* Field index */
			int isi, ili;		/* Input file sample and location indexes */
			char *dfnames[5][4] = {	/* Device colorspace names */
				{ "" },
				{ "GRAY_W" },
				{ "" },
				{ "RGB_R", "RGB_G", "RGB_B" },
				{ "CMYK_C", "CMYK_M", "CMYK_Y", "CMYK_K" }
			};
			int odim = 0;		/* Output file device dimensionality */
			int dfi[5][4];		/* Output file device colorspace indexes */
			int osi;			/* Output file sample id index */
	
			/* Setup input .ti2 file */
			icg = new_cgats();			/* Create a CGATS structure */
			icg->add_other(icg, "CTI2"); 	/* Calibration Target Information 2 */
			if (icg->read_name(icg, datin_name))
				error("CGATS file read error : %s",icg->err);
	
			if (icg->t[0].tt != tt_other || icg->t[0].oi != 0)
				error ("Input file '%s' isn't a CTI2 format file",datin_name);

			if (icg->ntables != 1)
				error ("Input file '%s' doesn't contain exactly one table",datin_name);
	
			if ((npat = icg->t[0].nsets) <= 0)
				error ("Input file '%s' doesn't contain any data sets",datin_name);
	
			/* Figure out the color space */
			if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
				error ("Input file '%s' doesn't contain keyword COLOR_REP",datin_name);

			if (strcmp(icg->t[0].kdata[fi],"CMYK") == 0) {
				dim = 4;
			} else if (strcmp(icg->t[0].kdata[fi],"RGB") == 0) {
				dim = 3;
			} else if (strcmp(icg->t[0].kdata[fi],"W") == 0) {
				dim = 1;
			} else
				error ("Input file '%s' keyword COLOR_REP has unknown value",datin_name);

			/* Find fields we want in the input file */
			if ((isi = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
				error ("Input file '%s' doesn't contain field SAMPLE_ID", datin_name);
			if (icg->t[0].ftype[isi] != nqcs_t)
				error ("Field SAMPLE_ID is wrong type");
		
			if ((ili = icg->find_field(icg, 0, "SAMPLE_LOC")) < 0)
				error ("Input file '%s' doesn't contain field SAMPLE_LOC", datin_name);
			if (icg->t[0].ftype[ili] != cs_t
			 && icg->t[0].ftype[ili] != nqcs_t)
				error ("Field SAMPLE_LOC is wrong type");

			/* Setup input/output .ti3 file */
			ocg = new_cgats();			/* Create a CGATS structure */
			ocg->add_other(ocg, "CTI3"); 	/* Calibration Target Information 3 */
			if (ocg->read_name(ocg, datout_name))
				error("CGATS file read error : %s",ocg->err);
	
			if (ocg->t[0].tt != tt_other || ocg->t[0].oi != 0)
				error ("Input file '%s' isn't a CTI3 format file",datout_name);

			if (ocg->ntables != 1)
				error ("Input file '%s' doesn't contain exactly one table",datout_name);
	
			if (npat != ocg->t[0].nsets)
				error ("Input file '%s' doesn't contain same number of data sets",datout_name);
	

			/* Find the fields we want in the output file */

			/* Figure out the color space */
			if ((fi = ocg->find_kword(ocg, 0, "COLOR_REP")) < 0)
				error ("Input file '%s' doesn't contain keyword COLOR_REP",datin_name);

			if (strncmp(ocg->t[0].kdata[fi],"CMYK",4) == 0) {
				odim = 4;
			} else if (strncmp(ocg->t[0].kdata[fi],"RGB",3) == 0) {
				odim = 3;
			} else if (strncmp(ocg->t[0].kdata[fi],"W",1) == 0) {
				odim = 1;
			} else
				error ("Input file '%s' keyword COLOR_REP has unknown value",datout_name);

			if (odim != dim)
				error ("File '%s' has different device space to '%s'",datin_name, datout_name);

			if ((osi = ocg->find_field(ocg, 0, "SAMPLE_ID")) < 0)
				error ("Input file doesn't contain field SAMPLE_ID");
			if (ocg->t[0].ftype[osi] != nqcs_t)
				error ("Field SAMPLE_ID is wrong type");
	
			for (i = 0; i < dim; i++) {
				if ((dfi[dim][i] = ocg->find_field(ocg, 0, dfnames[dim][i])) < 0)
					error ("Input file doesn't contain field %s", datout_name, dfnames[dim][i]);
				if (ocg->t[0].ftype[dfi[dim][i]] != r_t)
					error ("Field %s is wrong type",dfnames[dim][i]);
			}
	
			/* Initialise, ready to read out all the values */
			for (i = sr->reset(sr); i > 0; i--) {	/* For all samples in .tiff file */
				char loc[100];		/* Target patch location */
				double P[4];		/* Robust/raw mean values */
				int k, e;

				if (tmean)
					sr->read(sr, loc, NULL, P, NULL, NULL);
				else
					sr->read(sr, loc, P, NULL, NULL, NULL);
		
				/* Search for this location in the .ti2 file */
				for (j = 0; j < npat; j++) {
					if (strcmp(loc, (char *)icg->t[0].fdata[j][ili]) == 0) {
						char *sidp = (char *)icg->t[0].fdata[j][isi];

						/* Search for this sample id in .ti3 file */
						for (k = 0; k < npat; k++) {
							if (strcmp(sidp, (char *)ocg->t[0].fdata[k][osi]) == 0) {
								/* Update the device values */
								for (e = 0; e < dim; e++) {
									double vv = 100.0 * P[e]/255.0;
									*((double *)ocg->t[0].fdata[k][dfi[dim][e]]) = vv;
								}
								break;
							}	
						}
						if (k >= npat && verb >= 1)
							printf("Warning: Couldn't find sample '%s' in '%s'\n",sidp,datout_name);
						break;
					}
				}
				if (j >= npat && verb >= 1)
					printf("Warning: Couldn't find location '%s' in '%s'\n",loc,datin_name);
			}
	
			/* Flush our changes */
			if (ocg->write_name(ocg, datout_name))
				error("Write error : %s",ocg->err);
	
			ocg->del(ocg);		/* Clean up */
			icg->del(icg);		/* Clean up */

		/* ---------------------------------------------------------- */
		} else if (colm > 0) {	/* Using the scanner as a colorimeter */
			/* All this needs to track the code in spectro/printread.c */
			cgats *icg;			/* input .ti2 cgats structure */
			cgats *ocg;			/* input/output .ti3 cgats structure */
			icmFile *rd_fp;		/* Scanner to CIE lookup */
			icc *rd_icco;
			icmLuBase *luo;
			time_t clk = time(0);
			struct tm *tsp = localtime(&clk);
			char *atm = asctime(tsp); /* Ascii time */
			int nmask = 0;		/* Device colorant mask */
			int nchan = 0;		/* Number of device chanels */
			int npat;			/* Number of input patches (inc. padding) */
			int nopat = 0;		/* Number of output patches */
			int si;				/* Sample id index */
			int li;				/* Location id index */
			int ti;				/* Temp index */
			int fi;				/* Colorspace index */
	
			icg = new_cgats();			/* Create a CGATS structure */
			icg->add_other(icg, "CTI2"); 	/* special type Calibration Target Information 2 */
		
			if (icg->read_name(icg, datin_name))
				error("CGATS file read error : %s",icg->err);
		
			if (icg->t[0].tt != tt_other || icg->t[0].oi != 0)
				error ("Input file isn't a CTI2 format file");
			if (icg->ntables != 1)
				error ("Input file doesn't contain exactly one table");
		
			if ((npat = icg->t[0].nsets) <= 0)
				error ("No sets of data");
		
			/* Setup output cgats file */
			ocg = new_cgats();			/* Create a CGATS structure */
			ocg->add_other(ocg, "CTI3");	/* special type Calibration Target Information 3 */

			if (colm > 1) {		/* Appending information to .ti3 */

				if (ocg->read_name(ocg, datout_name))
					error("CGATS file read error on '%s': %s",datout_name, ocg->err);
		
				if (ocg->t[0].tt != tt_other || ocg->t[0].oi != 0)
					error ("Input file isn't a CTI3 format file");
				if (ocg->ntables != 1)
					error ("Input .TI3 file doesn't contain exactly one table");
				if ((nopat = ocg->t[0].nsets) <= 0)
					error ("No existing sets of .TI3 data");

			} else {			/* Creating .ti3 */

				ocg->add_table(ocg, tt_other, 0);	/* Start the first table */
		
				ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
				ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll printread", NULL);
				atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
				ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
				ocg->add_kword(ocg, 0, "DEVICE_CLASS","OUTPUT", NULL);	/* What sort of device this is */
				if ((ti = icg->find_kword(icg, 0, "SINGLE_DIM_STEPS")) >= 0)
					ocg->add_kword(ocg, 0, "SINGLE_DIM_STEPS",icg->t[0].kdata[ti], NULL);
			
				if ((ti = icg->find_kword(icg, 0, "COMP_GREY_STEPS")) >= 0)
					ocg->add_kword(ocg, 0, "COMP_GREY_STEPS",icg->t[0].kdata[ti], NULL);
		
				if ((ti = icg->find_kword(icg, 0, "MULTI_DIM_STEPS")) >= 0)
					ocg->add_kword(ocg, 0, "MULTI_DIM_STEPS",icg->t[0].kdata[ti], NULL);
		
				if ((ti = icg->find_kword(icg, 0, "FULL_SPREAD_PATCHES")) >= 0)
					ocg->add_kword(ocg, 0, "FULL_SPREAD_PATCHES",icg->t[0].kdata[ti], NULL);

				if ((ti = icg->find_kword(icg, 0, "TOTAL_INK_LIMIT")) >= 0)
					ocg->add_kword(ocg, 0, "TOTAL_INK_LIMIT",icg->t[0].kdata[ti], NULL);
			}

			if ((si = icg->find_field(icg, 0, "SAMPLE_ID")) < 0)
				error ("Input file doesn't contain field SAMPLE_ID");
			if (icg->t[0].ftype[si] != nqcs_t)
				error ("Field SAMPLE_ID is wrong type");
		
			/* Fields we want */
			if (colm > 1) {		/* Appending information to .ti3 */
				if ((ti = ocg->find_field(ocg, 0, "SAMPLE_ID")) != 0)
					error ("Field SAMPLE_ID (%d) not in expected location (%d) in '%s'",
					       ti, 0, datout_name);
			} else {
				ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);
			}
		
			if ((li = icg->find_field(icg, 0, "SAMPLE_LOC")) < 0)
				error ("Input file doesn't contain field SAMPLE_LOC");
			if (icg->t[0].ftype[li] != cs_t)
				error ("Field SAMPLE_LOC is wrong type");
		
			/* Figure out the color space */
			if ((fi = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
				error ("Input file doesn't contain keyword COLOR_REPS");
		
			if ((nmask = icx_char2inkmask(icg->t[0].kdata[fi])) != 0) {
				int i, j, ii;
				int chix[ICX_MXINKS];	/* Device chanel indexes */
				int xyzix[3];			/* XYZ chanel indexes */
				char *ident;
				char *xyzfname[3] = { "XYZ_X", "XYZ_Y", "XYZ_Z" };
		
				nchan = icx_noofinks(nmask);
				ident = icx_inkmask2char(nmask); 
		
				/* Device channels */
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
			
					if (colm > 1) {		/* Appending information to .ti3 */
						if (ocg->find_field(ocg, 0, fname) != 1 + j)
							error ("Field %s not in expected location in '%s'",fname, datout_name);
					} else {
						ocg->add_field(ocg, 0, fname, r_t);
					}
					chix[j] = ii;
				}
		
				/* Approximate XYZ and real XYZ */
				for (j = 0; j < 3; j++) {
					if ((ii = icg->find_field(icg, 0, xyzfname[j])) >= 0) {

						if (icg->t[0].ftype[ii] != r_t)
							error ("Field %s is wrong type",xyzfname[j]);
					}
			
					if (colm > 1) {		/* Appending information to .ti3 */
						if (ocg->find_field(ocg, 0, xyzfname[j]) != 1 + nchan + j)
							error ("Field %s not in expected location in '%s'"
							                            ,xyzfname[j], datout_name);
					} else {
						ocg->add_field(ocg, 0, xyzfname[j], r_t);
					}
					xyzix[j] = ii;
				}
		
				if (colm <= 1) {		/* Creating .ti3 */
					char fname[100];
					sprintf(fname, "%s_XYZ", ident);
					ocg->add_kword(ocg, 0, "COLOR_REP", fname, NULL);
				}
		
				if (colm > 1) {		/* Appending .ti3 data */

					/* Check that all the patches match */
					for (ii = i = 0; i < npat; i++) {

						if (strcmp(((char *)icg->t[0].fdata[i][si]), "0") == 0)
							continue;			/* Padding, so skip it */
		
						/* Id's */
						if (strcmp (((char *)icg->t[0].fdata[i][si]),
								    ((char *)ocg->t[0].fdata[ii][si])) != 0)
							error (".TI2 and .Ti3 field id's don't match at patch %d\n",i+1);

						/* device values */
						for (j = 0; j < nchan; j++) {
							double ival, oval;
							ival = *((double *)icg->t[0].fdata[i][chix[j]]);
							oval = *((double *)ocg->t[0].fdata[ii][1 + j]);
							if (fabs(ival - oval) > 0.001)
								error (".TI2 and .Ti3 device values (%f %f) don't match at patch %d %d\n",ival, oval, i+1, ii+1);
						}
						ii++;
					}
					if (ii != nopat)
						error("Different number of patches in .ti3 (%d) to expected(%d)",nopat,ii);

				} else { /* Read all the test patches in, and create output slots */
					cgats_set_elem *setel;	/* Array of set value elements */

					if ((setel = (cgats_set_elem *)malloc(sizeof(cgats_set_elem) * (1 + nchan + 3))) == NULL)
						error("Malloc failed!");


					for (ii = i = 0; i < npat; i++) {
						int k = 0;

						if (strcmp(((char *)icg->t[0].fdata[i][si]), "0") == 0)
							continue;			/* Padding, so skip it */
		
						/* Id */
						setel[k++].c = ((char *)icg->t[0].fdata[i][si]);
					
						/* device values */
						for (j = 0; j < nchan; j++) {
							setel[k++].d = *((double *)icg->t[0].fdata[i][chix[j]]);
						}

						/* Unset XYZ values */
						setel[k++].d = -1.0;
						setel[k++].d = -1.0;
						setel[k++].d = -1.0;

						ocg->add_setarr(ocg, 0, setel);

						ii++;
					}
					nopat = ii;
					free(setel);
				}
				free(ident);

			} else
				error ("Input file keyword COLOR_REPS has unknown value");
		
			/* Setup scanner RGB to XYZ conversion */
			{
				int inn, outn;			/* Chanels for input and output spaces */
				icColorSpaceSignature ins, outs;	/* Type of input and output spaces */
				int rv;

				/* Open up the file for reading */
				if ((rd_fp = new_icmFileStd_name(prof_name,"r")) == NULL)
					error("Write: Can't open file '%s'",prof_name);
		
				if ((rd_icco = new_icc()) == NULL)
					error("Read: Creation of ICC object failed");
		
				/* Read the header and tag list */
				if ((rv = rd_icco->read(rd_icco,rd_fp,0)) != 0)
					error("Read: %d, %s",rv,rd_icco->err);
		
				/* Check that this is an input profile */
				if (rd_icco->header->deviceClass != icSigInputClass)
					error ("ICC profile is expected to be an input profile");

				/* Get the Fwd table, absolute with XYZ override */
				if ((luo = rd_icco->get_luobj(rd_icco, icmFwd, icAbsoluteColorimetric,
				                              icSigXYZData, icmLuOrdNorm)) == NULL) {
					error("%d, %s",rd_icco->errc, rd_icco->err);
				}

				/* Get details of conversion */
				luo->spaces(luo, &ins, &inn, &outs, &outn, NULL, NULL, NULL, NULL);

				/* Check that it matches what we expect */
				if (inn != depth || tiffs != ins)
					error ("ICC profile doesn't match TIFF file type");
			}

			/* Initialise, ready to read out all the values */
			for (i = sr->reset(sr); i > 0; i--) {	/* For all samples in .tiff file */
				char loc[100];		/* Target patch location */
				double P[ICX_MXINKS];	/* Robust/true mean values */
				double xyz[3];			/* profile XYZ value */
				int k, e;

				if (tmean)
					sr->read(sr, loc, NULL, P, NULL, NULL);
				else
					sr->read(sr, loc, P, NULL, NULL, NULL);
		
				/* Search for this location in the .ti2 file */
				for (j = 0; j < npat; j++) {
					if (strcmp(loc, (char *)icg->t[0].fdata[j][li]) == 0) {	/* Got location */
						char *sidp = (char *)icg->t[0].fdata[j][si];	/* Get id */

						if (strcmp(sidp, "0") == 0)
							continue;			/* Padding, so skip it */

						/* Search for this sample id in .ti3 file */
						for (k = 0; k < nopat; k++) {
							if (strcmp(sidp, (char *)ocg->t[0].fdata[k][si]) == 0) {
								
//printf("Loc %s, ID %s got RGB value %f %f %f\n",sidp, loc, P[0], P[1], P[2]);

								/* Convert RGB to XYZ */
								for (e = 0; e < depth; e++)
									P[e] /= 255.0;			/* Convert to 0.0 .. 1.0 range */

								/* Convert to XYZ */
								luo->lookup(luo, xyz, P);

								/* Sanity check XYZ ? */
								// ~~~99

								/* Update the XYZ values */
								for (e = 0; e < 3; e++) {
									double ev = *((double *)ocg->t[0].fdata[k][1 + nchan + e]);

									if (ev != -1.0)
										error("Found an existing value in .ti3 file (%f)",ev);

									*((double *)ocg->t[0].fdata[k][1 + nchan + e]) = 100.0 * xyz[e];
								}
								break;
							}	
						}
						if (k >= nopat)
							error("Couldn't find sample '%s' in '%s'\n",sidp,datout_name);
						break;
					}
				}
				if (j >= npat && verb >= 1)
					error("Couldn't find location '%s' in '%s'\n",loc, datin_name);
			}
	
			/* Warn if not all patch values have been filled */
			if (verb) {
				int e, k;
				for (k = 0; k < nopat; k++) {
					for (e = 0; e < 3; e++) {
						double ev = *((double *)ocg->t[0].fdata[k][1 + nchan + e]);

						if (ev == -1.0)
							break;
					}
					if (e < 3)
						break;
				}
				if (k < nopat)
					printf("Not all sample values have been filled\n");
				else
					printf("All sample values have been filled\n");
			}

			if (ocg->write_name(ocg, datout_name))
				error("Write error : %s",ocg->err);
		
			luo->del(luo);
			rd_icco->del(rd_icco);
			rd_fp->del(rd_fp);

			ocg->del(ocg);		/* Clean up */
			icg->del(icg);		/* Clean up */
		
		/* ----------------------------------- */
		} else {	/* Normal scan calibration */
			cgats *icg;			/* input cgats structure */
			cgats *ocg;			/* output cgats structure */
			time_t clk = time(0);
			struct tm *tsp = localtime(&clk);
			char *atm = asctime(tsp); /* Ascii time */
			int sx;				/* Sample id index */
			int Xx, Yx, Zx;		/* XYZ_X, XYZ_Y, XYZ_Z index */
			int npat;			/* Number of test patches in it8 chart */
	
			icg = new_cgats();			/* Create a CGATS structure */
			icg->add_other(icg, ""); 	/* Accept any type */
			if (icg->read_name(icg, datin_name))
				error("CGATS file read error : %s",icg->err);
	
			/* ~~ should accept ti2 file and convert RGB to XYZ using    */
			/*    device cal., to make W/RGB/CMYK ->XYZ reading chart ~~ */
			if (icg->ntables < 1)
				error ("Input file doesn't contain at least one table");
	
			if ((npat = icg->t[0].nsets) <= 0)
				error ("No sets of data in first table");
	
			/* Setup output cgats file */
			ocg = new_cgats();	/* Create a CGATS structure */
			ocg->add_other(ocg, "CTI3"); 	/* our special type is Calibration Target Information 3 */
			ocg->add_table(ocg, tt_other, 0);	/* Start the first table */
	
			ocg->add_kword(ocg, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
			ocg->add_kword(ocg, 0, "ORIGINATOR", "Argyll target", NULL);
			atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
			ocg->add_kword(ocg, 0, "CREATED",atm, NULL);
	
			ocg->add_kword(ocg, 0, "DEVICE_CLASS","INPUT", NULL);	/* What sort of device this is */
			ocg->add_kword(ocg, 0, "COLOR_REP","XYZ_RGB", NULL);
	
			/* Fields we want from input chart reference file */
			if ((sx = icg->find_field(icg, 0, "Sample_Name")) < 0) {
				if ((sx = icg->find_field(icg, 0, "SAMPLE_NAME")) < 0) {
					if ((sx = icg->find_field(icg, 0, "SAMPLE_LOC")) < 0) {
						if ((sx = icg->find_field(icg, 0, "SAMPLE_ID")) < 0) {
							error ("Input file doesn't contain field SAMPLE_ID, Sample_Name or SAMPLE_NAME");
						}
					}
				}
			}
			if (icg->t[0].ftype[sx] != nqcs_t && icg->t[0].ftype[sx] != cs_t)
				error ("Field %s is wrong type", icg->t[0].fsym[sx]);
			if ((Xx = icg->find_field(icg, 0, "XYZ_X")) < 0)
				error ("Input file doesn't contain field XYZ_X");
			if (icg->t[0].ftype[Xx] != r_t)
				error ("Field XYZ_X is wrong type");
			if ((Yx = icg->find_field(icg, 0, "XYZ_Y")) < 0)
				error ("Input file doesn't contain field XYZ_Y");
			if (icg->t[0].ftype[Yx] != r_t)
				error ("Field XYZ_Y is wrong type");
			if ((Zx = icg->find_field(icg, 0, "XYZ_Z")) < 0)
				error ("Input file doesn't contain field XYZ_Z");
			if (icg->t[0].ftype[Zx] != r_t)
				error ("Field XYZ_Z is wrong type");
	
			ocg->add_field(ocg, 0, "SAMPLE_ID", nqcs_t);
			ocg->add_field(ocg, 0, "XYZ_X", r_t);
			ocg->add_field(ocg, 0, "XYZ_Y", r_t);
			ocg->add_field(ocg, 0, "XYZ_Z", r_t);
			if (depth == 1) {
				ocg->add_field(ocg, 0, "GREY", r_t);
				ocg->add_field(ocg, 0, "STDEV_GREY", r_t);
			} else if (depth == 3) {
				ocg->add_field(ocg, 0, "RGB_R", r_t);
				ocg->add_field(ocg, 0, "RGB_G", r_t);
				ocg->add_field(ocg, 0, "RGB_B", r_t);
				ocg->add_field(ocg, 0, "STDEV_R", r_t);
				ocg->add_field(ocg, 0, "STDEV_G", r_t);
				ocg->add_field(ocg, 0, "STDEV_B", r_t);
			} else if (depth == 4) {
				ocg->add_field(ocg, 0, "CMYK_C", r_t);
				ocg->add_field(ocg, 0, "CMYK_M", r_t);
				ocg->add_field(ocg, 0, "CMYK_Y", r_t);
				ocg->add_field(ocg, 0, "CMYK_K", r_t);
				ocg->add_field(ocg, 0, "STDEV_C", r_t);
				ocg->add_field(ocg, 0, "STDEV_M", r_t);
				ocg->add_field(ocg, 0, "STDEV_Y", r_t);
				ocg->add_field(ocg, 0, "STDEV_K", r_t);
			}
	
			/* Initialise, ready to read out all the values */
			for (j = 0; j < npat; j++) {
				char id[100];				/* Input patch id */
	
				/* Normalise labels */
				fix_it8(id,((char *)icg->t[0].fdata[j][sx]));	/* Copy and fix */
	
				/* Search for matching id */
				for (i = sr->reset(sr); i > 0; i--) {
					char tod[100];		/* Output patch id */
					char od[100];		/* Output patch id */
					double P[4];		/* Robust/true mean values */
					double sdP[4];		/* Standard deviation */

					if (tmean)
						sr->read(sr, tod, NULL, P, sdP, NULL);
					else
						sr->read(sr, tod, P, NULL, sdP, NULL);
					fix_it8(od,tod);
		
					if (strcmp(id,od) == 0) {
						if (depth == 1) {
							ocg->add_set( ocg, 0, id, 
						        *((double *)icg->t[0].fdata[j][Xx]),
						        *((double *)icg->t[0].fdata[j][Yx]),
						        *((double *)icg->t[0].fdata[j][Zx]),
								P[0] * 100.0/255.0, sdP[0] * 100.0/255.0);
						} else if (depth == 3) {
							ocg->add_set( ocg, 0, id, 
						        *((double *)icg->t[0].fdata[j][Xx]),
						        *((double *)icg->t[0].fdata[j][Yx]),
						        *((double *)icg->t[0].fdata[j][Zx]),
								P[0] * 100.0/255.0, P[1] * 100.0/255.0, P[2] * 100.0/255.0,
								sdP[0] * 100.0/255.0, sdP[1] * 100.0/255.0, sdP[2] * 100.0/255.0);
						} else if (depth == 4) {
							ocg->add_set( ocg, 0, id, 
						        *((double *)icg->t[0].fdata[j][Xx]),
						        *((double *)icg->t[0].fdata[j][Yx]),
						        *((double *)icg->t[0].fdata[j][Zx]),
								P[0] * 100.0/255.0, P[1] * 100.0/255.0,
							    P[2] * 100.0/255.0, P[3] * 100.0/255.0,
								sdP[0] * 100.0/255.0, sdP[1] * 100.0/255.0,
							    sdP[2] * 100.0/255.0, sdP[3] * 100.0/255.0);
						}
						break;
					}
				}
				if (i <= 0 && verb >= 1)
					printf("Warning: Couldn't match field '%s'\n",id);
			}
	
			if (ocg->write_name(ocg, datout_name))
				error("Write error : %s",ocg->err);
	
			ocg->del(ocg);		/* Clean up */
			icg->del(icg);		/* Clean up */

		}
	}

	/* Clean up */
	sr->free(sr);

	TIFFClose(rh);

	if (flags & SI_SHOW_FLAGS)
		TIFFClose(wh);

	return 0;
}

/* Fix IT8 chart labels */
void fix_it8(char *o, char *i) {
	if (strcmp(i,"Dmin")==0) {
		strcpy(o,"GS00");
		return;
	}
	if (strcmp(i,"Dmax")==0) {
		strcpy(o,"GS23");
		return;
	}
	while (!isdigit(*i) && *i != '\000') 	/* Skip non-numbers */
		*o++ = *i++; 
	if (i[0] != '\000' && i[1] == '\000')	/* Single last digit */
		*o++ = '0';				/* Add leading zero */
	strcpy(o, i);				/* Copy remainder */
}

/********************************************************************************/
#ifdef NEVER
/* Basic printf type error() and warning() routines */

#ifdef	__STDC__
void
error(char *fmt, ...)
#else
void
error(va_alist) 
va_dcl
#endif
{
	va_list args;
#ifndef	__STDC__
	char *fmt;
#endif

	fprintf(stderr,"scanin: Error - ");
#ifdef	__STDC__
	va_start(args, fmt);
#else
	va_start(args);
	fmt = va_arg(args, char *);
#endif
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	exit (-1);
}

#ifdef	__STDC__
void
warning(char *fmt, ...)
#else
void
warning(va_alist) 
va_dcl
#endif
{
	va_list args;
#ifndef	__STDC__
	char *fmt;
#endif

	fprintf(stderr,"scanin: Warning - ");
#ifdef	__STDC__
	va_start(args, fmt);
#else
	va_start(args);
	fmt = va_arg(args, char *);
#endif
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
}

#endif /* NEVER */
