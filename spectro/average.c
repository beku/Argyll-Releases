
/*
 * Argyll Color Correction System
 * Average measurement data
 *
 * Copyright 2008 Jordi Nodal
 * All rights reserved.
 *
 * Copyright 2006, 2007 Graeme Gill.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 * See the License3.txt file for details.
 *
 * This program performs an average/merge with two measurement files. 
 * All files have to be in CTI3 format (use Logo2Cgats or Spec2CIE).
 * Method: First it fills key field depending on device space and looks for this key 
 * in second file. If it exists, average is applied, else values are added directly increasing
 * patches number.
 */

/*
	NOTE :- this should really be re-written to use xcolorants rather than
	having hard codes RGB, CMYK and CMY !!!

	It also doesn't pass the calibration information though !!!

 */

#include <stdio.h>
#if defined(__IBMC__)
#include <float.h>
#endif
#include "config.h"
#include "numlib.h"
#include "cgats.h"
#include "xicc.h"
#include "insttypes.h"

#define Max(a,b) ((a > b) ? a : b)
#define Min(a,b) ((a < b) ? a : b)

void
usage (void)
{
	fprintf (stderr, "Average / Merge two measurement data files, Version %s\n", ARGYLL_VERSION_STR);
	fprintf (stderr, "Author: Jordi Nodal, licensed under the GPL Version 3\n");
	fprintf (stderr, "\n");
	fprintf (stderr, "Usage: average [options] input1.ti3 input2.ti3 output.ti3\n");
	fprintf (stderr, " -v          Verbose mode\n");
	fprintf (stderr, " input1.ti3   First measurement file\n");
	fprintf (stderr, " input2.ti3   Second measurement file\n");
	fprintf (stderr, " output.ti3   Average measurement file\n");
	exit (1);
}

/////////////////////////////////////////////////////////////
// Utils
////////////////////////////////////////////////////////////

int _findDouble(cgats *cg, int index, char* tag, double *dRet)
{
	int ti = 0;
	int ret = 0;

	if(cg){
		ti = cg->find_field(cg, 0, tag);
		if(ti >= 0){
			*dRet = *((double*)cg->t[0].fdata[index][ti]);
			ret = 1;
			}
		}

	return ret;
}

int _findString(cgats *cg, int index, char *tag, char *out_value, int out_size)
{
	char value[25];
	int ti = 0;
	int ret = 0;

	if(cg){
		ti = cg->find_field(cg, 0, tag);
		if(ti >= 0){
			strcpy(value, cg->t[0].fdata[index][ti]);
			ret = 1;
			}
		}

	if(strlen(value) <= out_size)
		strcpy(out_value, value);

	return ret;
}

int _getColorByDeviceSpace(icColorSpaceSignature space, cgats *cg, int index, double *out)
{
	int err = 0;

	if(!cg) 
		return -1;

	switch(space)
		{
		case icSigCmykData:
			_findDouble(cg, index, "CMYK_C", out);
			_findDouble(cg, index, "CMYK_M", out+1);
			_findDouble(cg, index, "CMYK_Y", out+2);
			_findDouble(cg, index, "CMYK_K", out+3);
			break;
		case icSigCmyData:
			_findDouble(cg, index, "CMY_C", out);
			_findDouble(cg, index, "CMY_M", out+1);
			_findDouble(cg, index, "CMY_Y", out+2);
			out[3] = 0.0;
			break;
		case icSigRgbData:
			_findDouble(cg, index, "RGB_R", out);
			_findDouble(cg, index, "RGB_G", out+1);
			_findDouble(cg, index, "RGB_B", out+2);
			out[3] = 0.0;
			break;
		default:
			err = -1;
			break;
		}	

	return err;
}

/////////////////////////////////////////////////////////////
// Main
////////////////////////////////////////////////////////////

int
main(int argc, char *argv[])
{
	long err;
	time_t clk;
	cgats *cgf;				/* output cgats file data */
	cgats *cg[2];				/* input cgats files*/
	cgats_set_elem *elems;
	struct tm *tsp;
	char *atm;
	long idProc;
	long *ptrIdProc;

	double dKey_f1[4];
	double dKey_f2[4];
	char dstBuff[10];

	int ii;
	long id; 

	/* The device colorspace
	*/
	icColorSpaceSignature devspace[2] = {icmSigDefaultData, icmSigDefaultData};	
	int isinv[2] = { 0, 0 };		/* Printer RGB ? */
	icColorSpaceSignature space = {icmSigDefaultData};

	int fa,nfa;				/* current argument we're looking at */
	int i, j;
	int verb = 0;

	char out_name[MAXNAMEL+1];		/* Output patch filename  */
	char in_name1[MAXNAMEL+4+1];	/* Input patch filename  */
	char in_name2[MAXNAMEL+4+1];	/* Input patch filename  */

	xspect sp1;		/* spectral file data */
	xspect sp2;		/* spectral file data */

	long spNum;
	double spLow, spHigh;

	int spi1[XSPECT_MAX_BANDS];		/* CGATS indexes for each wavelength */
	int spi2[XSPECT_MAX_BANDS];		/* CGATS indexes for each wavelength */
	
	double XYZ[3];
	double Lab[3];
	char buf[100];

	int isRgb[2];
	int isLab[2];
	int isXYZ[2];
	int isSpectral[2];

	int patches;
	int fieldNum;
	int nf;			// fields number
	int ti;
	int ciedata;
	int bFound;

	/* Init
	*/
	cgf = NULL;			/* output cgats file data */
	cg[0] = NULL;
	cg[1] = NULL;
	elems = NULL;
	nf = 0;
	ti = 0;	
	ciedata = 0;
	bFound = 0;

	clk = time(0);
	tsp = localtime(&clk);
	atm = asctime(tsp); /* Ascii time */

	isRgb[0]      = isRgb[1] = 0;
	isLab[0]	  = isLab[1] = 0;
	isXYZ[0]	  = isXYZ[1] = 0;
	isSpectral[0] = isSpectral[1] = 0;
 
	error_program = "Average";

	if (argc <= 1){
		usage ();
	}

	/* Process the arguments
	*/
	for(fa = 1;fa < argc;fa++) {

		nfa = fa;					
		if (argv[fa][0] == '-') {	
			char *na = NULL;		

			if (argv[fa][2] != '\000')
				na = &argv[fa][2];		
			else {
				if ((fa+1) < argc) {
					if (argv[fa+1][0] != '-') {
						nfa = fa + 1;
						na = argv[nfa];	
					}
				}
			}
			if (argv[fa][1] == '?'){
				usage ();
			}
			else if (argv[fa][1] == 'v' || argv[fa][1] == 'V') {
				verb = 1;
			} 
			else {
				usage ();
			}
		} else
			break;
	}

	/* Get the file name arguments */
	
	if (fa >= argc || argv[fa][0] == '-') error("Error");
 	strncpy(in_name1,argv[fa++],MAXNAMEL); in_name1[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') error("Error");
	strncpy(in_name2,argv[fa++],MAXNAMEL); in_name2[MAXNAMEL] = '\000';

	if (fa >= argc || argv[fa][0] == '-') error("Error");
	strncpy(out_name,argv[fa++],MAXNAMEL); out_name[MAXNAMEL] = '\000';

	if ((cg[0] = new_cgats()) == NULL){
		error("Failed to create cgats object");
	}
	cg[0]->add_other(cg[0], "CTI3");		/* Calibration Target Information 3 */

	if ((cg[1] = new_cgats()) == NULL){
		error("Failed to create cgats object");
	}
	cg[1]->add_other(cg[1], "CTI3");		/* Calibration Target Information 3 */

	if (cg[0]->read_name(cg[0], in_name1)){
		error("CGATS file '%s' read error : %s", in_name1,cg[0]->err);
	}
		
	if (cg[1]->read_name(cg[1], in_name2)){
		error("CGATS file '%s' read error : %s", in_name2,cg[1]->err);
	}
	
	/* at least one information table
	*/ 
	if (cg[0]->ntables < 1){
		error ("Input files '%s' don't contain at least one table", in_name1);
	}
	if (cg[1]->ntables < 1){
		error ("Input files '%s' don't contain at least one table", in_name2);
	}

	/* device color space
	*/
	if ((ti = cg[0]->find_kword(cg[0], 0, "COLOR_REP")) < 0){
		error("Input file doesn't contain keyword COLOR_REPS");
		}
	if (strncmp(cg[0]->t[0].kdata[ti],"CMYK_",5) == 0)
		devspace[0] = icSigCmykData;
	else if (strncmp(cg[0]->t[0].kdata[ti],"CMY_",4) == 0)
		devspace[0] = icSigCmyData;
	else if (strncmp(cg[0]->t[0].kdata[ti],"RGB_",4) == 0)
		devspace[0] = icSigRgbData;
	else if (strncmp(cg[0]->t[0].kdata[ti],"iRGB_",5) == 0) {
		devspace[0] = icSigRgbData;
		isinv[0] = 1;
	} else if (strncmp(cg[0]->t[0].kdata[ti],"K_",2) == 0) {
		devspace[0] = icSigGrayData;
	} else if (strncmp(cg[0]->t[0].kdata[ti],"W_",2) == 0) {
		devspace[0] = icSigGrayData;
	} else {
		error("Device has unhandled color representation '%s'",cg[0]->t[0].kdata[ti]);
	}

	if (devspace[0] == icSigRgbData) {
		if (cg[0]->find_field(cg[0], 0, "RGB_R") < 0 ||
			cg[0]->find_field(cg[0], 0, "RGB_G") < 0 ||
			cg[0]->find_field(cg[0], 0, "RGB_B") < 0){
			error("Input file doesn't contain correct fields");
			}
		}
	else{
		if(devspace[0] == icSigCmyData) {
			if (cg[0]->find_field(cg[0], 0, "CMY_C") < 0 ||
				cg[0]->find_field(cg[0], 0, "CMY_M") < 0 ||
				cg[0]->find_field(cg[0], 0, "CMY_Y") < 0){
				error("Input file doesn't contain correct fields");
				}
			}
		else{
			/* assume CMYK
			*/
			if (cg[0]->find_field(cg[0], 0, "CMYK_C") < 0 ||
				cg[0]->find_field(cg[0], 0, "CMYK_M") < 0 ||
				cg[0]->find_field(cg[0], 0, "CMYK_Y") < 0 ||
				cg[0]->find_field(cg[0], 0, "CMYK_K") < 0){
				error("Input file doesn't contain correct fields");
				}
			}
		}

	/* check spectral data in file 1
	*/
	if ((cg[0]->find_kword(cg[0], 0, "SPECTRAL_BANDS") >= 0 &&
		cg[0]->find_kword(cg[0], 0, "SPECTRAL_START_NM") >= 0	&&
		cg[0]->find_kword(cg[0], 0, "SPECTRAL_END_NM") >= 0))
		isSpectral[0] = 1;

	/* check CIE data in file 1
	*/
	if (cg[0]->find_field(cg[0], 0, "XYZ_X") >= 0) {
		if(cg[0]->find_field(cg[0], 0, "XYZ_Y") >= 0	&&
		   cg[0]->find_field(cg[0], 0, "XYZ_Z") >= 0)
		   isXYZ[0] = 1;
		}
	if (cg[0]->find_field(cg[0], 0, "LAB_L") >= 0) {
		if(cg[0]->find_field(cg[0], 0, "LAB_A") >= 0	&&
		   cg[0]->find_field(cg[0], 0, "LAB_B") >= 0)
		   isLab[0] = 1;
		}

	if ((ti = cg[1]->find_kword(cg[1], 0, "COLOR_REP")) < 0){
		error("Input file doesn't contain keyword COLOR_REPS");
		}
	if (strncmp(cg[1]->t[0].kdata[ti],"CMYK_",5) == 0)
		devspace[1] = icSigCmykData;
	else if (strncmp(cg[1]->t[0].kdata[ti],"CMY_",4) == 0)
		devspace[1] = icSigCmyData;
	else if (strncmp(cg[1]->t[0].kdata[ti],"RGB_",4) == 0)
		devspace[1] = icSigRgbData;
	else if (strncmp(cg[1]->t[0].kdata[ti],"iRGB_",5) == 0) {
		devspace[1] = icSigRgbData;
		isinv[1] = 1;
	} else if (strncmp(cg[1]->t[0].kdata[ti],"K_",2) == 0) {
		devspace[1] = icSigGrayData;
	} else if (strncmp(cg[1]->t[0].kdata[ti],"W_",2) == 0) {
		devspace[1] = icSigGrayData;
	} else {
		error("Device has unhandled color representation '%s'",cg[1]->t[0].kdata[ti]);
	}

	if (devspace[1] == icSigRgbData) {
		if (cg[1]->find_field(cg[1], 0, "RGB_R") < 0 ||
			cg[1]->find_field(cg[1], 0, "RGB_G") < 0 ||
			cg[1]->find_field(cg[1], 0, "RGB_B") < 0){
			error("Input file doesn't contain correct fields");
			}
		}
	else{
		if(devspace[1] == icSigCmyData) {
			if (cg[1]->find_field(cg[1], 0, "CMY_C") < 0 ||
				cg[1]->find_field(cg[1], 0, "CMY_M") < 0 ||
				cg[1]->find_field(cg[1], 0, "CMY_Y") < 0){
				error("Input file doesn't contain correct fields");
				}
			}
		else{
			/* assume CMYK
			*/
			if (cg[1]->find_field(cg[1], 0, "CMYK_C") < 0 ||
				cg[1]->find_field(cg[1], 0, "CMYK_M") < 0 ||
				cg[1]->find_field(cg[1], 0, "CMYK_Y") < 0 ||
				cg[1]->find_field(cg[1], 0, "CMYK_K") < 0){
				error("Input file doesn't contain correct fields");
				}
			}
		}

	/* check spectral data in file 2
	*/
	if ((cg[1]->find_kword(cg[1], 0, "SPECTRAL_BANDS") >= 0 &&
		cg[1]->find_kword(cg[1], 0, "SPECTRAL_START_NM") >= 0	&&
		cg[1]->find_kword(cg[1], 0, "SPECTRAL_END_NM") >= 0))
		isSpectral[1] = 1;

	/* check CIE data in file 2
	*/
	if (cg[1]->find_field(cg[1], 0, "XYZ_X") >= 0) {
		if(cg[1]->find_field(cg[1], 0, "XYZ_Y") >= 0	&&
		   cg[1]->find_field(cg[1], 0, "XYZ_Z") >= 0)
		   isXYZ[1] = 1;
		}
	if (cg[1]->find_field(cg[1], 0, "LAB_L") >= 0) {
		if(cg[1]->find_field(cg[1], 0, "LAB_A") >= 0	&&
		   cg[1]->find_field(cg[1], 0, "LAB_B") >= 0)
		   isLab[1] = 1;
		}

	/* check devices color spaces
	*/
	if(devspace[0] != devspace[1]){
		error("Can't perform average. Files are in different device space");
		}

	/* working space
	*/
	space = devspace[0];

	if (cg[0]->t[0].nsets < 0){
		error ("No patches in %s", in_name1);
	}
	if (cg[1]->t[0].nsets < 0){
		error ("No patches in %s", in_name2);
	}

	ciedata = 0;
	if(isXYZ[0] == 1 || isXYZ[1] == 1 || isLab[0] == 1 || isLab[1] == 1){
		ciedata = 1;
		if(verb)
			printf("Doing average with CIE data.\n");
		}
	else
		{
		if(verb) 
			printf("Doing average with spectral data.\n");
		}

	/* field counter
	*/
	nf = 0;

	/* Alloc
	*/
	idProc = 0;
	ptrIdProc = (long *) malloc(sizeof(long) * Max(cg[0]->t[0].nsets, cg[1]->t[0].nsets));

	/* Setup output cgats file */

	cgf = new_cgats();	/* Create a CGATS structure */
	cgf->add_other(cgf, "CTI3"); 	/* our special type is Calibration Target Information 3 */
	cgf->add_table(cgf, tt_other, 0);	/* Start the first table */
	cgf->add_kword(cgf, 0, "DESCRIPTOR", "Argyll Calibration Target chart information 3",NULL);
	cgf->add_kword(cgf, 0, "ORIGINATOR", "Argyll target", NULL);
	atm[strlen(atm)-1] = '\000';	/* Remove \n from end */
	cgf->add_kword(cgf, 0, "CREATED",atm, NULL);
	cgf->add_kword(cgf, 0, "DEVICE_CLASS","OUTPUT", NULL);
	cgf->add_kword(cgf, 0, "TARGET_INSTRUMENT", inst_name(instSpectrolino) , NULL);

	/* Fields we want */
	cgf->add_field(cgf, 0, "SAMPLE_ID", nqcs_t);nf++;

	switch(space)
		{
		case icSigCmykData:
			cgf->add_field(cgf, 0, "CMYK_C", r_t);nf++;
			cgf->add_field(cgf, 0, "CMYK_M", r_t);nf++;
			cgf->add_field(cgf, 0, "CMYK_Y", r_t);nf++;
			cgf->add_field(cgf, 0, "CMYK_K", r_t);nf++;
			if (isLab[0] == 1 || isLab[1] == 1)
				cgf->add_kword(cgf, 0, "COLOR_REP","CMYK_LAB", NULL);
			else
				cgf->add_kword(cgf, 0, "COLOR_REP","CMYK_XYZ", NULL);
			break;
		case icSigCmyData:
			cgf->add_field(cgf, 0, "CMY_C", r_t);nf++;
			cgf->add_field(cgf, 0, "CMY_M", r_t);nf++;
			cgf->add_field(cgf, 0, "CMY_Y", r_t);nf++;
			if (isLab[0] == 1 || isLab[1] == 1)
				cgf->add_kword(cgf, 0, "COLOR_REP","CMY_LAB", NULL);
			else
				cgf->add_kword(cgf, 0, "COLOR_REP","CMY_XYZ", NULL);
			break;
		case icSigRgbData:
			cgf->add_field(cgf, 0, "RGB_R", r_t);nf++;
			cgf->add_field(cgf, 0, "RGB_G", r_t);nf++;
			cgf->add_field(cgf, 0, "RGB_B", r_t);nf++;
			if (isLab[0] == 1 || isLab[1] == 1) {
				if (isinv[0] || isinv[1])
					cgf->add_kword(cgf, 0, "COLOR_REP","iRGB_LAB", NULL);
				else
					cgf->add_kword(cgf, 0, "COLOR_REP","RGB_LAB", NULL);
			} else {
				if (isinv[0] || isinv[1])
					cgf->add_kword(cgf, 0, "COLOR_REP","iRGB_XYZ", NULL);
				else
					cgf->add_kword(cgf, 0, "COLOR_REP","RGB_XYZ", NULL);
			}
			break;
		case icSigGrayData:
			error("Not supported space");
			break;
		}

	if(ciedata){
		/* perform average directly with CIE data values
		*/
		cgf->add_field(cgf, 0, "LAB_L", r_t);nf++;
		cgf->add_field(cgf, 0, "LAB_A", r_t);nf++;
		cgf->add_field(cgf, 0, "LAB_B", r_t);nf++;
		cgf->add_field(cgf, 0, "XYZ_X", r_t);nf++;
		cgf->add_field(cgf, 0, "XYZ_Y", r_t);nf++;
		cgf->add_field(cgf, 0, "XYZ_Z", r_t);nf++;
		}
	else{
		/* CIE data doesn't exist. Is there any spectral data?
		*/
		if(isSpectral[0] == 1 && isSpectral[1] == 1){
			i = cg[0]->find_kword (cg[0], 0, "SPECTRAL_BANDS");
			sp1.spec_n = atoi(cg[0]->t[0].kdata[i]);
			i = cg[0]->find_kword (cg[0], 0, "SPECTRAL_START_NM");
			sp1.spec_wl_short = atof(cg[0]->t[0].kdata[i]);
			i = cg[0]->find_kword (cg[0], 0, "SPECTRAL_END_NM");
			sp1.spec_wl_long = atof(cg[0]->t[0].kdata[i]);
			sp1.norm = 100.0;
		
			i = cg[0]->find_kword (cg[0], 0, "SPECTRAL_BANDS");
			sp2.spec_n = atoi (cg[0]->t[0].kdata[i]);
			i = cg[0]->find_kword (cg[0], 0, "SPECTRAL_START_NM");
			sp2.spec_wl_short = atof (cg[0]->t[0].kdata[i]);
			i = cg[0]->find_kword (cg[0], 0, "SPECTRAL_END_NM");
			sp2.spec_wl_long = atof (cg[0]->t[0].kdata[i]);
			sp2.norm = 100.0;

			/* spectral field counter
			*/
			spNum = Min(sp1.spec_n, sp2.spec_n);

			spLow = sp1.spec_wl_short;
			spLow = (spLow > sp2.spec_wl_short)?sp2.spec_wl_short:spLow;
			spHigh = sp1.spec_wl_long;
			spHigh = (spHigh < sp2.spec_wl_long)?sp2.spec_wl_long:spHigh;

			/* add spectral fields
			*/
			sprintf(buf,"%d", spNum);
			cgf->add_kword(cgf, 0, "SPECTRAL_BANDS",buf, NULL);
			sprintf(buf,"%d", spLow);
			cgf->add_kword(cgf, 0, "SPECTRAL_START_NM",buf, NULL);
			sprintf(buf,"%d", spHigh);
			cgf->add_kword(cgf, 0, "SPECTRAL_END_NM",buf, NULL);

			/* add spectral fields number to fields counter
			*/
			nf+=spNum;

			/* search field with spectral data
			*/
			for (j = 0; j < spNum; j++) 
				{
				int nm;
				nm = (int) (spLow +
							((double) j / (spNum - 1.0))
							* (spHigh - spLow) + 0.5);

				sprintf(buf, "SPEC_%03d", nm);
				cgf->add_field(cgf, 0, buf, r_t);
				if ((spi1[j] = cg[0]->find_field (cg[0], 0, buf)) < 0)
					error ("Input file doesn't contain field %s", buf);
				if ((spi2[j] = cg[1]->find_field (cg[1], 0, buf)) < 0)
					error ("Input file doesn't contain field %s", buf);
				}
			}
		}

	/* alloc for elements
	*/
	if ((elems = (cgats_set_elem *)
		calloc(nf, sizeof(cgats_set_elem))) == NULL){
		error("Out of memory");
		}
	
	for(i = 0;i < nf;i++){
		elems[i].c = NULL;	
		elems[i].d = 0.0; 	
		elems[i].i = 0;	
		}
	ii = 0;

	id = 1;
	for (i = 0; i < cg[0]->t[0].nsets; i++, id++) {
		fieldNum = 0;	/* reset */
		
		sprintf(dstBuff, "%d", id);
		elems[fieldNum++].c = dstBuff;

		if((err = _getColorByDeviceSpace(space, cg[0], i, &dKey_f1[0]))){
			error("Field error");
			}

		bFound = 0;

		/* search in file 2
		*/
		for(j = 0; j < cg[1]->t[0].nsets; j++){
			/*
			*/
			if((err = _getColorByDeviceSpace(space, cg[1], j, &dKey_f2[0]))){
				error("Field error");
				}

			/* key field compare
			*/
			if(dKey_f2[0] == dKey_f1[0] &&
				dKey_f2[1] == dKey_f1[1] &&
				dKey_f2[2] == dKey_f1[2] &&
				dKey_f2[3] == dKey_f1[3]){
					char Id[5];
					/* 
					*/
					bFound = 1;

					/* save processed id
					*/
					_findString(cg[1], j, "SAMPLE_ID", Id, sizeof(Id));
					ptrIdProc[idProc++] = atoi(Id);

					/* are identical, copy key field values and apply average with rest of fields
					*/
					switch(space){
						case icSigCmykData:
							ii = cgf->find_field(cgf, 0, "CMYK_C");
							elems[ii].d = dKey_f1[0];
							ii = cgf->find_field(cgf, 0, "CMYK_M");
							elems[ii].d = dKey_f1[1];
							ii = cgf->find_field(cgf, 0, "CMYK_Y");
							elems[ii].d = dKey_f1[2];
							ii = cgf->find_field(cgf, 0, "CMYK_K");
							elems[ii].d = dKey_f1[3];
							break;
						case icSigCmyData:
							ii = cgf->find_field(cgf, 0, "CMY_C");
							elems[ii].d = dKey_f1[0];
							ii = cgf->find_field(cgf, 0, "CMY_M");
							elems[ii].d = dKey_f1[1];
							ii = cgf->find_field(cgf, 0, "CMY_Y");
							elems[ii].d = dKey_f1[2];
							break;
						case icSigRgbData:
							ii = cgf->find_field(cgf, 0, "RGB_R");
							elems[ii].d = dKey_f1[0];
							ii = cgf->find_field(cgf, 0, "RGB_G");
							elems[ii].d = dKey_f1[1];
							ii = cgf->find_field(cgf, 0, "RGB_B");
							elems[ii].d = dKey_f1[2];
							break;
						default:
							error("Error");
						}

					if(ciedata){
						/* searching LAB_ and/or XYZ_ only
						*/
						if(isLab[0] > 0 && isLab[1] > 0){
							double labA[3];
							double labB[3];

							/* found Lab data in both files
							*/
							ii = cgf->find_field(cgf, 0, "LAB_L");
							_findDouble(cg[0], i, "LAB_L", &labA[0]);
							_findDouble(cg[1], j, "LAB_L", &labB[0]);
							elems[ii].d = (labA[0] + labB[0]) * .5;
							ii = cgf->find_field(cgf, 0, "LAB_A");
							_findDouble(cg[0], i, "LAB_A", &labA[1]);
							_findDouble(cg[1], j, "LAB_A", &labB[1]);
							elems[ii].d = (labA[1] + labB[1]) * .5;
							ii = cgf->find_field(cgf, 0, "LAB_B");
							_findDouble(cg[0], i, "LAB_B", &labA[2]);
							_findDouble(cg[1], j, "LAB_B", &labB[2]);
							elems[ii].d = (labA[2] + labB[2]) * .5;
							}
						if(isXYZ[0] > 0 && isXYZ[1] > 0){
							double xyzA[3];
							double xyzB[3];
							/* found XYZ data in both files
							*/
							ii = cgf->find_field(cgf, 0, "XYZ_X");
							_findDouble(cg[0], i, "XYZ_X", &xyzA[0]);
							_findDouble(cg[1], j, "XYZ_X", &xyzB[0]);
							elems[ii].d = (xyzA[0] + xyzB[0]) * .5;
							ii = cgf->find_field(cgf, 0, "XYZ_Y");
							_findDouble(cg[0], i, "XYZ_Y", &xyzA[1]);
							_findDouble(cg[1], j, "XYZ_Y", &xyzB[1]);
							elems[ii].d = (xyzA[1] + xyzB[1]) * .5;
							ii = cgf->find_field(cgf, 0, "XYZ_Z");
							_findDouble(cg[0], i, "XYZ_Z", &xyzA[2]);
							_findDouble(cg[1], j, "XYZ_Z", &xyzB[2]);
							elems[ii].d = (xyzA[2] + xyzB[2]) * .5;
							}
						if(isLab[0] > 0 && isLab[1] == 0){
							double labA[3];
							double labB[3];
							double XYZa[3];
							double XYZb[3];
							/* check if file 2 has XYZ data to convert to Lab 
							   NOTE: not sure if it's correct
							*/
							_findDouble(cg[0], i, "LAB_L", &labA[0]);
							_findDouble(cg[0], i, "LAB_A", &labA[1]);
							_findDouble(cg[0], i, "LAB_B", &labA[2]);

							if(isXYZ[0] > 0){
								_findDouble(cg[0], i, "XYZ_X", &XYZa[0]);
								_findDouble(cg[0], i, "XYZ_Y", &XYZa[1]);
								_findDouble(cg[0], i, "XYZ_Z", &XYZa[2]);
								}
							else{
								icmLab2XYZ(&icmD50, XYZa, labA);
								XYZa[0]*=100;
								XYZa[1]*=100;
								XYZa[2]*=100;
								}

							if(isXYZ[1] > 0){
								_findDouble(cg[1], j, "XYZ_X", &XYZb[0]);
								_findDouble(cg[1], j, "XYZ_Y", &XYZb[1]);
								_findDouble(cg[1], j, "XYZ_Z", &XYZb[2]);
								XYZb[0]/=100;
								XYZb[1]/=100;
								XYZb[2]/=100;
								/* Lab conversion
								*/
								icmXYZ2Lab(&icmD50, labB, XYZb);
								}
							else
								error("Error");

							ii = cgf->find_field(cgf, 0, "LAB_L");
							elems[ii].d = (labA[0] + labB[0]) * .5;
							ii = cgf->find_field(cgf, 0, "LAB_A");
							elems[ii].d = (labA[1] + labB[1]) * .5;
							ii = cgf->find_field(cgf, 0, "LAB_B");
							elems[ii].d = (labA[2] + labB[2]) * .5;
							ii = cgf->find_field(cgf, 0, "XYZ_X");
							elems[ii].d = (XYZa[0] + XYZb[0]) * .5;
							ii = cgf->find_field(cgf, 0, "XYZ_Y");
							elems[ii].d = (XYZa[1] + XYZb[1]) * .5;
							ii = cgf->find_field(cgf, 0, "XYZ_Z");
							elems[ii].d = (XYZa[2] + XYZb[2]) * .5;
							}
						else{
							if(isLab[0] == 0 && isLab[1] > 0){
								double labA[3];
								double labB[3];
								double XYZa[3];
								double XYZb[3];
								/* check if file 1 has XYZ data to convert to Lab 
								   NOTE: not sure if it's correct
								*/
								_findDouble(cg[1], j, "LAB_L", &labB[0]);
								_findDouble(cg[1], j, "LAB_A", &labB[1]);
								_findDouble(cg[1], j, "LAB_B", &labB[2]);
								if(isXYZ[1] > 0){
									_findDouble(cg[1], j, "XYZ_X", &XYZb[0]);
									_findDouble(cg[1], j, "XYZ_Y", &XYZb[1]);
									_findDouble(cg[1], j, "XYZ_Z", &XYZb[2]);
									}
								else{
									icmLab2XYZ(&icmD50, XYZb, labB);
									XYZb[0]*=100;
									XYZb[1]*=100;
									XYZb[2]*=100;
									}

								if(isXYZ[0] > 0){
									_findDouble(cg[0], i, "XYZ_X", &XYZa[0]);
									_findDouble(cg[0], i, "XYZ_Y", &XYZa[1]);
									_findDouble(cg[0], i, "XYZ_Z", &XYZa[2]);
									XYZa[0]/=100;
									XYZa[1]/=100;
									XYZa[2]/=100;
									/* Lab conversion
									*/
									icmXYZ2Lab(&icmD50, labA, XYZa);
									}
								else
									error("Error");

								ii = cgf->find_field(cgf, 0, "LAB_L");
								elems[ii].d = (labB[0] + labA[0]) * .5;
								ii = cgf->find_field(cgf, 0, "LAB_A");
								elems[ii].d = (labB[1] + labA[1]) * .5;
								ii = cgf->find_field(cgf, 0, "LAB_B");
								elems[ii].d = (labB[2] + labA[2]) * .5;
								ii = cgf->find_field(cgf, 0, "XYZ_X");
								elems[ii].d = (XYZb[0] + XYZa[0]) * .5;
								ii = cgf->find_field(cgf, 0, "XYZ_Y");
								elems[ii].d = (XYZb[1] + XYZa[1]) * .5;
								ii = cgf->find_field(cgf, 0, "XYZ_Z");
								elems[ii].d = (XYZb[2] + XYZa[2]) * .5;
								}
							}	
						}
					else
						{
						/* exist spectral data 
						*/ 
						long n;
						for ( n = 0; n < spNum; n++) 
							{
							int nm;
							nm = (int) (spLow +
										((double) n / (spNum - 1.0))
										* (spHigh - spLow) + 0.5);
							sprintf (buf, "SPEC_%03d", nm);

							ii = cgf->find_field (cgf, 0, buf);
							elems[ii].d  = (*((double *) cg[0]->t[0].fdata[i][cg[0]->find_field (cg[0], 0, buf)]) + 
											*((double *) cg[1]->t[0].fdata[j][cg[1]->find_field (cg[1], 0, buf)])) * .5;

							}
						}
					break;
					}
				else
					continue;
			}

		if(!bFound)
			{
			/* key field not found, add it
			*/
			if((err = _getColorByDeviceSpace(space, cg[0], i, &dKey_f1[0])))
				error("Field error");

			switch(space){
				case icSigCmykData:
					ii = cgf->find_field(cgf, 0, "CMYK_C");
					elems[ii].d = dKey_f1[0];
					ii = cgf->find_field(cgf, 0, "CMYK_M");
					elems[ii].d = dKey_f1[1];
					ii = cgf->find_field(cgf, 0, "CMYK_Y");
					elems[ii].d = dKey_f1[2];
					ii = cgf->find_field(cgf, 0, "CMYK_K");
					elems[ii].d = dKey_f1[3];
					break;
				case icSigCmyData:
					ii = cgf->find_field(cgf, 0, "CMY_C");
					elems[ii].d = dKey_f1[0];
					ii = cgf->find_field(cgf, 0, "CMY_M");
					elems[ii].d = dKey_f1[1];
					ii = cgf->find_field(cgf, 0, "CMY_Y");
					elems[ii].d = dKey_f1[2];
					break;
				case icSigRgbData:
					ii = cgf->find_field(cgf, 0, "RGB_R");
					elems[ii].d = dKey_f1[0];
					ii = cgf->find_field(cgf, 0, "RGB_G");
					elems[ii].d = dKey_f1[1];
					ii = cgf->find_field(cgf, 0, "RGB_B");
					elems[ii].d = dKey_f1[2];
					break;
				default:
					error("Error");
				}

			if(ciedata){
				ii = cgf->find_field(cgf, 0, "LAB_L");elems[ii].d = 0.0;
				ii = cgf->find_field(cgf, 0, "LAB_A");elems[ii].d = 0.0;
				ii = cgf->find_field(cgf, 0, "LAB_B");elems[ii].d = 0.0;
				ii = cgf->find_field(cgf, 0, "XYZ_X");elems[ii].d = 0.0;
				ii = cgf->find_field(cgf, 0, "XYZ_Y");elems[ii].d = 0.0;
				ii = cgf->find_field(cgf, 0, "XYZ_Z");elems[ii].d = 0.0;

				if(isLab[0] > 0){
					double lab[3];
					ii = cgf->find_field(cgf, 0, "LAB_L");
					_findDouble(cg[0], i, "LAB_L", &lab[0]);
					elems[ii].d = lab[0];
					ii = cgf->find_field(cgf, 0, "LAB_A");
					_findDouble(cg[0], i, "LAB_A", &lab[1]);
					elems[ii].d = lab[1];
					ii = cgf->find_field(cgf, 0, "LAB_B");
					_findDouble(cg[0], i, "LAB_B", &lab[2]);
					elems[ii].d = lab[2];
					}

				if(isXYZ[0] > 0){
					double XYZ[3];
					ii = cgf->find_field(cgf, 0, "XYZ_X");
					_findDouble(cg[0], i, "XYZ_X", &XYZ[0]);
					elems[ii].d = XYZ[0];
					ii = cgf->find_field(cgf, 0, "XYZ_Y");
					_findDouble(cg[0], i, "XYZ_Y", &XYZ[1]);
					elems[ii].d = XYZ[1];
					ii = cgf->find_field(cgf, 0, "XYZ_Z");
					_findDouble(cg[0], i, "XYZ_Z", &XYZ[2]);
					elems[ii].d = XYZ[2];
					}
				}	
			else
				{
				if(isSpectral[0] > 0){
					long n;
					for ( n = 0; n < spNum; n++) 
						{
						int nm;
						nm = (int) (spLow +
									((double) n / (spNum - 1.0))
									* (spHigh - spLow) + 0.5);
						sprintf(buf, "SPEC_%03d", nm);

						ii = cgf->find_field (cgf, 0, buf);
						elems[ii].d  = (*((double *) cg[0]->t[0].fdata[i][cg[0]->find_field (cg[0], 0, buf)]));
						}
					}
				}
				
			}
		cgf->add_setarr(cgf, 0, elems);			// salvamos array de datos
	}

	/* add fields not processeds from second file
	*/
	for (i = 0; i < cg[1]->t[0].nsets; i++, id++) {
		char si[5];
		int iid;

		fieldNum = 0;	// reset
		bFound = 0;
		
		sprintf(dstBuff, "%d", id);
		elems[fieldNum++].c = (dstBuff);

		_findString(cg[1], i, "SAMPLE_ID", si, sizeof(si));
		iid = atoi(si);

		/* search 'id' in processed items list. Add it when not found.
		*/
		for(j = 0;j < idProc;j++){
			if(ptrIdProc[j] == iid){
				bFound = 1;
				break;
				}
			}
		if(!bFound){
			if((err = _getColorByDeviceSpace(space, cg[1], i, &dKey_f2[0]))){
				error("Field error");
				}

			switch(space){
				case icSigCmykData:
					ii = cgf->find_field(cgf, 0, "CMYK_C");
					elems[ii].d = dKey_f2[0];
					ii = cgf->find_field(cgf, 0, "CMYK_M");
					elems[ii].d = dKey_f2[1];
					ii = cgf->find_field(cgf, 0, "CMYK_Y");
					elems[ii].d = dKey_f2[2];
					ii = cgf->find_field(cgf, 0, "CMYK_K");
					elems[ii].d = dKey_f2[3];
					break;
				case icSigCmyData:
					ii = cgf->find_field(cgf, 0, "CMY_C");
					elems[ii].d = dKey_f2[0];
					ii = cgf->find_field(cgf, 0, "CMY_M");
					elems[ii].d = dKey_f2[1];
					ii = cgf->find_field(cgf, 0, "CMY_Y");
					elems[ii].d = dKey_f2[2];
					break;
				case icSigRgbData:
					ii = cgf->find_field(cgf, 0, "RGB_R");
					elems[ii].d = dKey_f2[0];
					ii = cgf->find_field(cgf, 0, "RGB_G");
					elems[ii].d = dKey_f2[1];
					ii = cgf->find_field(cgf, 0, "RGB_B");
					elems[ii].d = dKey_f2[2];
					break;
				default:
					error("Error");
				}

			if(ciedata){
				ii = cgf->find_field(cgf, 0, "LAB_L");elems[ii].d = 0.0;
				ii = cgf->find_field(cgf, 0, "LAB_A");elems[ii].d = 0.0;
				ii = cgf->find_field(cgf, 0, "LAB_B");elems[ii].d = 0.0;
				ii = cgf->find_field(cgf, 0, "XYZ_X");elems[ii].d = 0.0;
				ii = cgf->find_field(cgf, 0, "XYZ_Y");elems[ii].d = 0.0;
				ii = cgf->find_field(cgf, 0, "XYZ_Z");elems[ii].d = 0.0;

				if(isLab[1] > 0){
					double lab[3];
					ii = cgf->find_field(cgf, 0, "LAB_L");
					_findDouble(cg[1], i, "LAB_L", &lab[0]);elems[ii].d = lab[0];
					ii = cgf->find_field(cgf, 0, "LAB_A");
					_findDouble(cg[1], i, "LAB_A", &lab[1]);elems[ii].d = lab[1];
					ii = cgf->find_field(cgf, 0, "LAB_B");
					_findDouble(cg[1], i, "LAB_B", &lab[2]);elems[ii].d = lab[2];
					}

				if(isXYZ[1] > 0){
					double XYZ[3];
					ii = cgf->find_field(cgf, 0, "XYZ_X");
					_findDouble(cg[1], i, "XYZ_X", &XYZ[0]);elems[ii].d = XYZ[0];
					ii = cgf->find_field(cgf, 0, "XYZ_Y");
					_findDouble(cg[1], i, "XYZ_Y", &XYZ[1]);elems[ii].d = XYZ[1];
					ii = cgf->find_field(cgf, 0, "XYZ_Z");
					_findDouble(cg[1], i, "XYZ_Z", &XYZ[2]);elems[ii].d = XYZ[2];
					}
				}
			else
				{
				if(isSpectral[1] > 0){
					long n;
					for ( n = 0; n < spNum; n++) 
						{
						int nm;

						/* Compute nearest integer wavelength 
						*/
						nm = (int) (spLow +
									((double) n / (spNum - 1.0))
									* (spHigh - spLow) + 0.5);
						sprintf(buf, "SPEC_%03d", nm);

						ii = cgf->find_field (cgf, 0, buf);
						elems[ii].d  = (*((double *) cg[1]->t[0].fdata[i][cg[1]->find_field (cg[1], 0, buf)]));
						}
					}
				}
			/* save data
			*/
			cgf->add_setarr(cgf, 0, elems);			
			}
		}

	/* Free
	*/
	if(ptrIdProc)
		free (ptrIdProc);

	/* save
	*/
	if (cgf->write_name(cgf, out_name)){
		error ("Write error: %s", cgf->err);
		}

	if(verb)
		printf("Resultant file %s contains %d patches\n", out_name, cgf->t[0].nsets);

	/* Clean up 
	*/
	if(cgf != NULL)
		cgf->del (cgf);					
	if(cg[0] != NULL)
		cg[0]->del (cg[0]);
	if(cg[1] != NULL)
		cg[1]->del (cg[1]);
	if(elems != NULL)
		free (elems);
	
	return err;
}

