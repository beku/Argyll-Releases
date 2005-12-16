
/* 
 * Argyll Color Correction System
 * Input device profile creator.
 *
 * Author: Graeme W. Gill
 * Date:   11/10/00
 *
 * Copyright 2000 Graeme W. Gill
 * All rights reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

/*
 * This program takes in the scattered test chart
 * points, and interpolates them into a gridded
 * forward ICC device profile.
 * It also creates (at the moment) a limited backward
 * profile based on the forward grid.
 *
 */

/*
 * TTBD:
 *      Need to make this more of a library:
 *  Fix error handling
 *  fix verbose output
 *  hand icc object back rather than writing file ?
 */

#undef DEBUG

#define verbo stdout

#include <stdio.h>
#include "cgats.h"
#include "icc.h"
#include "xicc.h"
#include "numlib.h"
#include "rspl.h"
#include "prof.h"

/*
   Basic algorithm outline:

 Scanner:

   Figure out the input curves to give
   the flattest grid.

   Figure out the grid values.

   Use them to generate the A2B table.

   Do all the calculations in Lab space,
   but represent the profile in XYZ space, so that
   the white/black point normalisation doesn't cause
   the clut values to be clipped.

   This leads to a poorer accuracy as an XYZ profile,
   but can then be compensated for using the ICX_MERGE_CLUT flag
   together with a PCS override.

*/


/* ---------------------------------------- */

/* Make an input device profile, where we create an A2B lut */
/* directly from the scattered input data. */
void
make_input_icc(
	prof_atype ptype,		/* Profile algorithm type */
	int verb,
	int iquality,			/* A2B table quality, 0..3 */
	int verify,
	int nsabs,				/* nz for non-standard absolute output */
	char *file_name,		/* output icc name */
	cgats *icg,				/* input cgats structure */
	double smooth,			/* RSPL smoothing factor, -ve if raw */
	double avgdev,			/* reading Average Deviation as a proportion of the input range */
	profxinf *xpi			/* Optional Profile creation extra data */
) {
	icmFile *wr_fp;
	icc *wr_icco;
	int npat;				/* Number of patches */
	co *tpat;				/* Patch input values */
	int i, rv = 0;
	int isLab = 0;			/* 0 if input is XYZ, 1 if input is Lab */
	int wantLab = 0;		/* 0 if want is XYZ, 1 want is Lab */
	int isLut = 0;			/* 0 if shaper+ matrix, 1 if lut type */
	int isShTRC = 0;		/* 0 if separate gamma/shaper TRC, 1 if shared */

	if (ptype == prof_clutLab) {		/* Lab lut */
		wantLab = 1;
		isLut = 1;
	} else if (ptype == prof_clutXYZ) {	/* XYZ lut */
		wantLab = 0;
		isLut = 1;
	} else {
		wantLab = 0;			/* gamma/shaper + matrix profile must be XYZ */
		isLut = 0;

		if (ptype == prof_gam1mat	
		 || ptype == prof_sha1mat) {
			isShTRC = 1;		/* Single curve */
		}
	}

	/* Open up the file for writing */
	if ((wr_fp = new_icmFileStd_name(file_name,"w")) == NULL)
		error ("Write: Can't open file '%s'",file_name);

	if ((wr_icco = new_icc()) == NULL)
		error ("Write: Creation of ICC object failed");

	/* Add all the tags required */

	/* The header: */
	{
		icmHeader *wh = wr_icco->header;

		/* Values that must be set before writing */
		wh->deviceClass     = icSigInputClass;
    	wh->colorSpace      = icSigRgbData;				/* It's an RGB profile */
		if (wantLab)
	    	wh->pcs         = icSigLabData;
		else
	    	wh->pcs         = icSigXYZData;
    	wh->renderingIntent = icRelativeColorimetric;	/* For want of something */

		/* Values that should be set before writing */
		if (xpi != NULL && xpi->manufacturer != 0L)
			wh->manufacturer = xpi->manufacturer;
		else
			wh->manufacturer = str2tag("????");

		if (xpi != NULL && xpi->model != 0L)
			wh->model = xpi->model;
		else
	    	wh->model = str2tag("????");

		/* Values that may be set before writing */
		if (xpi != NULL && xpi->creator != 0L)
			wh->creator = xpi->creator;
	}
	/* Profile Description Tag: */
	{
		icmTextDescription *wo;
		char *dst;			/* description */

		if (xpi != NULL && xpi->profDesc != NULL)
			dst = xpi->profDesc;
		else {
			dst = "This is a Lut style RGB - XYZ Input Profile";
		}

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
		char *crt;

		if (xpi != NULL && xpi->copyright != NULL)
			crt = xpi->copyright;
		else
			crt = "Copyright, the creator of this profile";

		if ((wo = (icmText *)wr_icco->add_tag(
		           wr_icco, icSigCopyrightTag,	icSigTextType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(crt)+1; 	/* Allocated and used size of text, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->data, crt);		/* Copy the text in */
	}
	/* Device Manufacturers Description Tag: */
	if (xpi != NULL && xpi->deviceMfgDesc != NULL) {
		icmTextDescription *wo;
		char *dst = xpi->deviceMfgDesc;

		if ((wo = (icmTextDescription *)wr_icco->add_tag(
		           wr_icco, icSigDeviceMfgDescTag,	icSigTextDescriptionType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, dst);		/* Copy the string in */
	}
	/* Model Description Tag: */
	if (xpi != NULL && xpi->modelDesc != NULL) {
		icmTextDescription *wo;
		char *dst = xpi->modelDesc;

		if ((wo = (icmTextDescription *)wr_icco->add_tag(
		           wr_icco, icSigDeviceModelDescTag,	icSigTextDescriptionType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = strlen(dst)+1; 	/* Allocated and used size of desc, inc null */
		wo->allocate((icmBase *)wo);/* Allocate space */
		strcpy(wo->desc, dst);		/* Copy the string in */
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
		wo->data[0].X = 0.9642;		/* Set a default value - D50 */
		wo->data[0].Y = 1.0000;
		wo->data[0].Z = 0.8249;
	}
	/* Black Point Tag: */
	{
		icmXYZArray *wo;
		if ((wo = (icmXYZArray *)wr_icco->add_tag(
		           wr_icco, icSigMediaBlackPointTag, icSigXYZArrayType)) == NULL) 
			error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

		wo->size = 1;
		wo->allocate((icmBase *)wo);	/* Allocate space */
		wo->data[0].X = 0.00;		/* Set a default value */
		wo->data[0].Y = 0.00;
		wo->data[0].Z = 0.00;
	}

	if (isLut) {		/* Lut type profile */

		/* 16 bit dev -> pcs lut: */
		{
			icmLut *wo;

			if ((wo = (icmLut *)wr_icco->add_tag(
			           wr_icco, icSigAToB0Tag,	icSigLut16Type)) == NULL) 
				error("add_tag failed: %d, %s",wr_icco->errc,wr_icco->err);

			wo->inputChan = 3;
			wo->outputChan = 3;
			if (iquality >= 3) {
		    	wo->clutPoints = 45;
		    	wo->inputEnt = 2048;
		    	wo->outputEnt = 2048;
			} else if (iquality == 2) {
		    	wo->clutPoints = 33;
		    	wo->inputEnt = 2048;
		    	wo->outputEnt = 2048;
			} else if (iquality == 1) {
		    	wo->clutPoints = 17;
		    	wo->inputEnt = 1024;
		    	wo->outputEnt = 1024;
			} else {
		    	wo->clutPoints = 9;
		    	wo->inputEnt = 512;
		    	wo->outputEnt = 512;
			}

			wo->allocate((icmBase *)wo);/* Allocate space */

			/* icxLuLut will set tables values */
		}

	} else {	/* shaper + matrix type */

		/* Red, Green and Blue Colorant Tags: */
		{
			icmXYZArray *wor, *wog, *wob;
			if ((wor = (icmXYZArray *)wr_icco->add_tag(
			           wr_icco, icSigRedColorantTag, icSigXYZArrayType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);
			if ((wog = (icmXYZArray *)wr_icco->add_tag(
			           wr_icco, icSigGreenColorantTag, icSigXYZArrayType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);
			if ((wob = (icmXYZArray *)wr_icco->add_tag(
			           wr_icco, icSigBlueColorantTag, icSigXYZArrayType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);

			wor->size = wog->size = wob->size = 1;
			wor->allocate((icmBase *)wor);	/* Allocate space */
			wog->allocate((icmBase *)wog);
			wob->allocate((icmBase *)wob);

			/* Setup some sane dummy values */
			/* icxMatrix will override these later */
			wor->data[0].X = 1.0; wor->data[0].Y = 0.0; wor->data[0].Z = 0.0;
			wog->data[0].X = 0.0; wog->data[0].Y = 1.0; wog->data[0].Z = 0.0;
			wob->data[0].X = 0.0; wob->data[0].Y = 0.0; wob->data[0].Z = 1.0;
		}

		/* Red, Green and Blue Tone Reproduction Curve Tags: */
		{
			icmCurve *wor, *wog, *wob;
			int i;
			if ((wor = (icmCurve *)wr_icco->add_tag(
			           wr_icco, icSigRedTRCTag, icSigCurveType)) == NULL) 
				error("add_tag failed: %d, %s",rv,wr_icco->err);

			if (isShTRC) {	/* Make all TRCs shared */
				if ((wog = (icmCurve *)wr_icco->link_tag(
				           wr_icco, icSigGreenTRCTag, icSigRedTRCTag)) == NULL) 
					error("link_tag failed: %d, %s",rv,wr_icco->err);
				if ((wob = (icmCurve *)wr_icco->link_tag(
				           wr_icco, icSigBlueTRCTag, icSigRedTRCTag)) == NULL) 
					error("link_tag failed: %d, %s",rv,wr_icco->err);

			} else {		/* Else individual */
				if ((wog = (icmCurve *)wr_icco->add_tag(
				           wr_icco, icSigGreenTRCTag, icSigCurveType)) == NULL) 
					error("add_tag failed: %d, %s",rv,wr_icco->err);
				if ((wob = (icmCurve *)wr_icco->add_tag(
				           wr_icco, icSigBlueTRCTag, icSigCurveType)) == NULL) 
					error("add_tag failed: %d, %s",rv,wr_icco->err);
			}
	
			if (ptype == prof_shamat || ptype == prof_sha1mat) {	/* Shaper */
				wor->flag = wog->flag = wob->flag = icmCurveSpec; 
				wor->size = wog->size = wob->size = 256;			/* Number of entries */
			} else {						/* Gamma */
				wor->flag = wog->flag = wob->flag = icmCurveGamma;
				wor->size = wog->size = wob->size = 1;				/* Must be 1 for gamma */
			}
			wor->allocate((icmBase *)wor);	/* Allocate space */
			wog->allocate((icmBase *)wog);
			wob->allocate((icmBase *)wob);

			/* icxMatrix will set curve values */
		}
	}

	if ((npat = icg->t[0].nsets) <= 0)
		error ("No sets of data");

	if (verb) {
		fprintf(verbo,"No of test patches = %d\n",npat);
	}

	/* Allocate arrays to hold test patch input and output values */
	if ((tpat = (co *)malloc(sizeof(co) * npat)) == NULL)
		error("Malloc failed - tpat[]");

	/* Read in the CGATs fields */
	{
		int ti, ii;
		int Xi, Yi, Zi;
		int ri, gi, bi;

		/* Check that we handle the color space */
		if ((ti = icg->find_kword(icg, 0, "COLOR_REP")) < 0)
			error ("Input file doesn't contain keyword COLOR_REPS");
		if (strcmp(icg->t[0].kdata[ti],"LAB_RGB") == 0) {
			isLab = 1;
		} else {
			if (strcmp(icg->t[0].kdata[ti],"XYZ_RGB") == 0) {
				isLab = 0;
			} else {
				error ("Input device input file has unhandled color representation");
			}
		}

		if (isLab) {
			if ((Xi = icg->find_field(icg, 0, "LAB_L")) < 0)
				error ("Input file doesn't contain field LAB_L");
			if (icg->t[0].ftype[Xi] != r_t)
				error ("Field LAB_L is wrong type - corrupted file ?");
			if ((Yi = icg->find_field(icg, 0, "LAB_A")) < 0)
				error ("Input file doesn't contain field LAB_A");
			if (icg->t[0].ftype[Yi] != r_t)
				error ("Field LAB_A is wrong type - corrupted file ?");
			if ((Zi = icg->find_field(icg, 0, "LAB_B")) < 0)
				error ("Input file doesn't contain field LAB_B");
			if (icg->t[0].ftype[Zi] != r_t)
				error ("Field LAB_B is wrong type - corrupted file ?");
		} else {
			if ((Xi = icg->find_field(icg, 0, "XYZ_X")) < 0)
				error ("Input file doesn't contain field XYZ_X");
			if (icg->t[0].ftype[Xi] != r_t)
				error ("Field XYZ_X is wrong type - corrupted file ?");
			if ((Yi = icg->find_field(icg, 0, "XYZ_Y")) < 0)
				error ("Input file doesn't contain field XYZ_Y");
			if (icg->t[0].ftype[Yi] != r_t)
				error ("Field XYZ_Y is wrong type - corrupted file ?");
			if ((Zi = icg->find_field(icg, 0, "XYZ_Z")) < 0)
				error ("Input file doesn't contain field XYZ_Z");
			if (icg->t[0].ftype[Zi] != r_t)
				error ("Field XYZ_Z is wrong type - corrupted file ?");
		}

		if ((ri = icg->find_field(icg, 0, "RGB_R")) < 0)
			error ("Input file doesn't contain field RGB_R");
		if (icg->t[0].ftype[ri] != r_t)
			error ("Field CMYK_C is wrong type - corrupted file ?");
		if ((gi = icg->find_field(icg, 0, "RGB_G")) < 0)
			error ("Input file doesn't contain field RGB_G");
		if (icg->t[0].ftype[gi] != r_t)
			error ("Field CMYK_M is wrong type - corrupted file ?");
		if ((bi = icg->find_field(icg, 0, "RGB_B")) < 0)
			error ("Input file doesn't contain field RGB_B");
		if (icg->t[0].ftype[bi] != r_t)
			error ("Field CMYK_Y is wrong type - corrupted file ?");

		for (i = 0; i < npat; i++) {
			tpat[i].p[0] = *((double *)icg->t[0].fdata[i][ri]) / 100.0;
			tpat[i].p[1] = *((double *)icg->t[0].fdata[i][gi]) / 100.0;
			tpat[i].p[2] = *((double *)icg->t[0].fdata[i][bi]) / 100.0;
			if (tpat[i].p[0] > 1.0
			 || tpat[i].p[1] > 1.0
			 || tpat[i].p[2] > 1.0) {
				error("Device value field exceeds 100.0!");
			}
			tpat[i].v[0] = *((double *)icg->t[0].fdata[i][Xi]);
			tpat[i].v[1] = *((double *)icg->t[0].fdata[i][Yi]);
			tpat[i].v[2] = *((double *)icg->t[0].fdata[i][Zi]);
			if (!isLab) {
				tpat[i].v[0] /= 100.0;		/* Normalise XYZ to range 0.0 - 1.0 */
				tpat[i].v[1] /= 100.0;
				tpat[i].v[2] /= 100.0;
			}
			if (!isLab && wantLab) { /* Convert test patch result XYZ to PCS (D50 Lab) */
				icmXYZ2Lab(&icmD50, tpat[i].v, tpat[i].v);
			} else if (isLab && !wantLab) {
				icmLab2XYZ(&icmD50, tpat[i].v, tpat[i].v);
			}
		}
	}	/* End of reading in CGATs file */

	if (isLut) {
		int flags = 0;

		xicc *wr_xicc;			/* extention object */
		icxLuBase *AtoB;		/* AtoB ixcLu */

		/* Wrap with an expanded icc */
		if ((wr_xicc = new_xicc(wr_icco)) == NULL)
			error ("Creation of xicc failed");

		if (verb)
			flags |= ICX_VERBOSE;

		if (nsabs == 0)
	        flags |= ICX_SET_WHITE | ICX_SET_BLACK;		/* Set white and black */

		/* Setup RGB -> Lab conversion object from scattered data. */
		/* Note that we've layered it on a native XYZ icc profile. */
		if ((AtoB = wr_xicc->set_luobj(
		               wr_xicc, icmFwd, icmDefaultIntent,
		               icmLuOrdNorm,
		               flags, 		/* Flags */
		               npat, tpat, smooth, avgdev, NULL, NULL, iquality)) == NULL)
			error ("%d, %s",wr_xicc->errc, wr_xicc->err);

		/* Free up xicc stuff */
		AtoB->del(AtoB);
		wr_xicc->del(wr_xicc);

	} else {		/* Gamma/Shaper + matrix profile */

		xicc *wr_xicc;			/* extention object */
		icxLuBase *xluo;		/* Forward ixcLu */
		int flags = 0;

		/* Wrap with an expanded icc */
		if ((wr_xicc = new_xicc(wr_icco)) == NULL)
			error("Creation of xicc failed");
		
		if (verb)
			flags |= ICX_VERBOSE;

		/* Setup Device -> XYZ conversion (Fwd) object from scattered data. */
		if ((xluo = wr_xicc->set_luobj(
//		               wr_xicc, icmFwd, icRelativeColorimetric,
		               wr_xicc, icmFwd, icmDefaultIntent,
		               icmLuOrdNorm,
		               flags | ICX_SET_WHITE | ICX_SET_BLACK, 		/* Flags */
		               npat, tpat, smooth, avgdev, NULL, NULL, iquality)) == NULL)
			error("%d, %s",wr_xicc->errc, wr_xicc->err);

		/* Free up xicc stuff */
		xluo->del(xluo);
		wr_xicc->del(wr_xicc);
	}

	/* Write the file (including all tags) out */
	if ((rv = wr_icco->write(wr_icco,wr_fp,0)) != 0)
		error ("Write file: %d, %s",rv,wr_icco->err);

	/* Close the file */
	wr_icco->del(wr_icco);
	wr_fp->del(wr_fp);

	/* Check the profile accuracy against the data points */
	if (verb || verify) {
		icmFile *rd_fp;
		icc *rd_icco;
		icmLuBase *luo;
		double merr = 0.0;
		double aerr = 0.0;
		double nsamps = 0.0;

		/* Open up the file for reading */
		if ((rd_fp = new_icmFileStd_name(file_name,"r")) == NULL)
			error ("Write: Can't open file '%s'",file_name);

		if ((rd_icco = new_icc()) == NULL)
			error ("Write: Creation of ICC object failed");

		/* Read the header and tag list */
		if ((rv = rd_icco->read(rd_icco,rd_fp,0)) != 0)
			error ("Read: %d, %s",rv,rd_icco->err);

		/* ~~ should use an xluobj with merge output ~~~ */
		/* Get the A2B table */
		if ((luo = rd_icco->get_luobj(rd_icco, icmFwd,
                           icAbsoluteColorimetric, icSigLabData, icmLuOrdNorm)) == NULL) {
			error ("%d, %s",rd_icco->errc, rd_icco->err);
		}

		for (i = 0; i < npat; i++) {
			double out[3], ref[3];
			double mxd;

			if (luo->lookup(luo, out, tpat[i].p) > 1)
				error ("%d, %s",rd_icco->errc,rd_icco->err);
		
			/* Our tpat data might be in XYZ, so generate an Lab ref value */
			if (!wantLab) { /* Convert test patch result XYZ to PCS (D50 Lab) */
				icmXYZ2Lab(&icmD50, ref, tpat[i].v);

			} else {
				ref[0] = tpat[i].v[0];
				ref[1] = tpat[i].v[1];
				ref[2] = tpat[i].v[2];
			}

			if (verb && verify) {
				printf("[%f] %f %f %f -> %f %f %f should be %f %f %f\n",
				       icmLabDE(ref, out),
				       tpat[i].p[0],tpat[i].p[1],tpat[i].p[2],
				       out[0],out[1],out[2],
				       ref[0],ref[1],ref[2]);
			}

			/* Check the result */
			mxd = icmLabDE(ref, out);
			if (mxd > merr)
				merr = mxd;

			aerr += mxd;
			nsamps++;
		}
		printf("profile check complete, peak err = %f, avg err = %f\n",merr,aerr/nsamps);

		/* Done with lookup object */
		luo->del(luo);

		/* Close the file */
		rd_icco->del(rd_icco);
		rd_fp->del(rd_fp);
	}
}


