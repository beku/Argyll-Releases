/* Gretag Spectrolino/Spectroscan command parser
 *
 * Author: Neil Okamoto
 * Date: 1/9/2001
 *
 * Copyright 2001, DreamWorks LLC
 * All Rights Reserved.
 *
 * This material is licenced under the GNU GENERAL PUBLIC LICENCE :-
 * see the Licence.txt file for licencing details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "spm.h"

void spm_error(char *fmt, ...);

#undef DEBUG

/***********************************************************************
 * message protocol definition table
 ***********************************************************************/

typedef struct {
	int   cmd;					/* command code */
	int   type;					/* command or reply? */
	int   timeout;				/* recommended timeout */
	int   reply;				/* if cmd, expect this reply */
	char *fmt;					/* argument format */
	char *desc;					/* text description */
} spm_message;

spm_message spm_protocol[] = {

	/* Spectrolino */

	{SPM_PROTO_EXECMEASUREMENT,        SPM_CMD, 6,  SPM_PROTO_EXECERROR,              NULL,                                  "ExecMeasurement"},
	{SPM_PROTO_EXECREFMEASUREMENT,     SPM_CMD, 6,  SPM_PROTO_EXECERROR,              "0x09 %c",                             "ExecRefMeasurement"},
	{SPM_PROTO_EXECWHITEMEASUREMENT,   SPM_CMD, 6,  SPM_PROTO_EXECERROR,              NULL,                                  "ExecWhiteMeasurement"},
	{SPM_PROTO_SETMEASUREMENTOUTPUT,   SPM_CMD, 6,  SPM_PROTO_DOWNLOADERROR,          "%c %c",                               "SetMeasurementOutput"},
	{SPM_PROTO_PARAMETERREQUEST,       SPM_CMD, 6,  SPM_PROTO_PARAMETERANSWER,        NULL,                                  "ParameterRequest"},
	{SPM_PROTO_PARAMETERDOWNLOAD,      SPM_CMD, 6,  SPM_PROTO_DOWNLOADERROR,          "%c %c %c %c",                         "ParameterDownload"},
	{SPM_PROTO_SPECPARAMETERREQUEST,   SPM_CMD, 6,  SPM_PROTO_SPECPARAMETERANSWER,    "0x09 %c",                             "SpecParameterRequest"},
	{SPM_PROTO_SPECTRUMREQUEST,        SPM_CMD, 6,  SPM_PROTO_SPECTRUMANSWER,         "0x09 %c",                             "SpectrumRequest"},
	{SPM_PROTO_CPARAMETERREQUEST,      SPM_CMD, 6,  SPM_PROTO_CPARAMETERANSWER,       "0x09 %c",                             "CParameterRequest"},
	{SPM_PROTO_CREQUEST,               SPM_CMD, 6,  SPM_PROTO_CANSWER,                "0x09 %c",                             "CRequest"},
	{SPM_PROTO_DENSITYPARAMETERREQUEST,SPM_CMD, 6,  SPM_PROTO_DENSITYPARAMETERANSWER, "0x09",                                "DensityParameterRequest"},
	{SPM_PROTO_DENSITYREQUEST,         SPM_CMD, 6,  SPM_PROTO_DENSITYANSWER,          "0x09",                                "DensityRequest"},
	{SPM_PROTO_DMAXREQUEST,            SPM_CMD, 6,  SPM_PROTO_DMAXANSWER,             "0x09",                                "DmaxRequest"},
	{SPM_PROTO_WHITEREFERENCEREQUEST,  SPM_CMD, 6,  SPM_PROTO_WHITEREFERENCEANSWER,   "%c",                                  "WhiteReferenceRequest"},
	{SPM_PROTO_WHITEREFERENCEDOWNLOAD, SPM_CMD, 6,  SPM_PROTO_DOWNLOADERROR,          "%c %f*36 %c*18",                      "WhiteReferenceDownload"},
	{SPM_PROTO_EXECWHITEREFTOORIGDAT,  SPM_CMD, 6,  SPM_PROTO_EXECERROR,              NULL,                                  "WhiteRefToOrigDat"},
	{SPM_PROTO_MEASCONTROLDOWNLOAD,    SPM_CMD, 6,  SPM_PROTO_DOWNLOADERROR,          "%c",                                  "MeasControlDownload"},
	{SPM_PROTO_MEASCONTROLREQUEST,     SPM_CMD, 6,  SPM_PROTO_MEASCONTROLANSWER,      "%c",                                  "MeasControlRequest"},
	{SPM_PROTO_FLOATREQUEST,           SPM_CMD, 6,  SPM_PROTO_FLOATANSWER,            "%c",                                  "FloatRequest"},
	{SPM_PROTO_FLOATDOWNLOAD,          SPM_CMD, 6,  SPM_PROTO_DOWNLOADERROR,          "%c %f",                               "FloatDownload"},
	{SPM_PROTO_ILLUMTABREQUEST,        SPM_CMD, 6,  SPM_PROTO_ILLUMTABANSWER,         "0x00 %c",                             "IllumTabRequest"},
	{SPM_PROTO_ILLUMTABDOWNLOAD,       SPM_CMD, 6,  SPM_PROTO_DOWNLOADERROR,          "0x08 %f*36",                          "IllumTabDownload"},
	{SPM_PROTO_GETVALNR,               SPM_CMD, 6,  SPM_PROTO_VALNRANSWER,            "0x60"                                 "GetValNr"},
	{SPM_PROTO_SETVALNR,               SPM_CMD, 6,  SPM_PROTO_DOWNLOADERROR,          "0x60 %d",                             "SetValNr"},
	{SPM_PROTO_SLOPEREQUEST,           SPM_CMD, 6,  SPM_PROTO_SLOPEANSWER,            NULL,                                  "SlopeRequest"},
	{SPM_PROTO_SLOPEDOWNLOAD,          SPM_CMD, 6,  SPM_PROTO_DOWNLOADERROR,          "%f*4",                                "SlopeDownload"},
	{SPM_PROTO_DENSTABREQUEST,         SPM_CMD, 6,  SPM_PROTO_DENSTABANSWER,          "0x00 %c",                             "DensTabRequest"},
	{SPM_PROTO_DENSTABDOWNLOAD,        SPM_CMD, 6,  SPM_PROTO_DOWNLOADERROR,          "%f*36 %f*36 %f*36 %f*36",             "DensTabDownload"},
	{SPM_PROTO_DEVICEDATAREQUEST,      SPM_CMD, 6,  SPM_PROTO_DEVICEDATAANSWER,       NULL,                                  "DeviceDataRequest"},
	{SPM_PROTO_TARGETIDREQUEST,        SPM_CMD, 6,  SPM_PROTO_TARGETIDANSWER,         NULL,                                  "TargetIdRequest"},
	{SPM_PROTO_NEWMEASUREREQUEST,      SPM_CMD, 6,  SPM_PROTO_NEWMEASUREANSWER,       NULL,                                  "NewMeasureRequest"},
	{SPM_PROTO_NEWKEYREQUEST,          SPM_CMD, 6,  SPM_PROTO_NEWKEYANSWER,           NULL,                                  "NewKeyRequest"},
	{SPM_PROTO_TARGETONOFFSTDOWNLOAD,  SPM_CMD, 6,  SPM_PROTO_DOWNLOADERROR,          "%c %c %c",                            "TargetOnOffStDownload"},
	{SPM_PROTO_RESETSTATUSDOWNLOAD,    SPM_CMD, 6,  SPM_PROTO_DOWNLOADERROR,          "0x01 0x04 %c",                        "ResetStatusDownload"},
	{SPM_PROTO_ACTERRORREQUEST,        SPM_CMD, 6,  SPM_PROTO_ACTERRORANSWER,         NULL,                                  "ActErrorRequest"},
	{SPM_PROTO_EXECERROR,              SPM_REP, 6,  SPM_PROTO_NONE,                   "%c",                                  "ExecError"},
	{SPM_PROTO_DOWNLOADERROR,          SPM_REP, 6,  SPM_PROTO_NONE,                   "%d",                                  "DownloadError"},
	{SPM_PROTO_PARAMETERANSWER,        SPM_REP, 6,  SPM_PROTO_NONE,                   "%c %c %c %c %c",                      "ParameterAnswer"},
	{SPM_PROTO_SPECPARAMETERANSWER,    SPM_REP, 6,  SPM_PROTO_NONE,                   "0x09 %c %f*36 %c %c %c 0x02 %d",      "SpecParameterAnswer"},
	{SPM_PROTO_SPECTRUMANSWER,         SPM_REP, 6,  SPM_PROTO_NONE,                   "0x09 %c %f*36 %c %d"                  "SpectrumAnswer"},
	{SPM_PROTO_CPARAMETERANSWER,       SPM_REP, 6,  SPM_PROTO_NONE,                   "0x09 %c %f*3 %c %c %c 0x02 %c %c %d", "CParameterAnswer"},
	{SPM_PROTO_CANSWER,                SPM_REP, 6,  SPM_PROTO_NONE,                   "0x09 %c %f*3 %c %d",                  "CAnswer"},
	{SPM_PROTO_DENSITYPARAMETERANSWER, SPM_REP, 6,  SPM_PROTO_NONE,                   "0x09 %f*4 %c %c %c %c 0x02 %c %d",    "DensityParameterAnswer"},
	{SPM_PROTO_DENSITYANSWER,          SPM_REP, 6,  SPM_PROTO_NONE,                   "0x09 %f*4 %c %c %d",                  "DensityAnswer"},
	{SPM_PROTO_DMAXANSWER,             SPM_REP, 6,  SPM_PROTO_NONE,                   "0x09 %c %c %c %c %d",                 "DmaxAnswer"},
	{SPM_PROTO_WHITEREFERENCEANSWER,   SPM_REP, 6,  SPM_PROTO_NONE,                   "%c %f*36 %c %c*18",                   "WhiteReferenceAnswer"},
	{SPM_PROTO_MEASCONTROLANSWER,      SPM_REP, 6,  SPM_PROTO_NONE,                   "%c %c %d",                            "MeasControlAnswer"},
	{SPM_PROTO_FLOATANSWER,            SPM_REP, 6,  SPM_PROTO_NONE,                   "%c %f %d",                            "FloatAnswer"},
	{SPM_PROTO_ILLUMTABANSWER,         SPM_REP, 6,  SPM_PROTO_NONE,                   "0x00 %c %f*36 %d",                    "IllumTabAnswer"},
	{SPM_PROTO_VALNRANSWER,            SPM_REP, 6,  SPM_PROTO_NONE,                   "0x60 %d %d",                          "ValNrAnswer"},
	{SPM_PROTO_SLOPEANSWER,            SPM_REP, 6,  SPM_PROTO_NONE,                   "%f*4",                                "SlopeAnswer"},
	{SPM_PROTO_DENSTABANSWER,          SPM_REP, 6,  SPM_PROTO_NONE,                   "0x00 %c %f*36 %f*36 %f*36 %f*36 %d",  "DensTabAnswer"},
	{SPM_PROTO_DEVICEDATAANSWER,       SPM_REP, 6,  SPM_PROTO_NONE,                   "%c*18 %c %c*8 %d*2 %c*12 %c*16",      "DeviceDataAnswer"},
	{SPM_PROTO_TARGETIDANSWER,         SPM_REP, 6,  SPM_PROTO_NONE,                   "%c*18 %d*7 %c %d*3",                  "TargetIdAnswer"},
	{SPM_PROTO_NEWMEASUREANSWER,       SPM_REP, 6,  SPM_PROTO_NONE,                   "%c 0x09",                             "NewMeasureAnswer"},
	{SPM_PROTO_NEWKEYANSWER,           SPM_REP, 6,  SPM_PROTO_NONE,                   "%c %d",                               "NewKeyAnswer"},
	{SPM_PROTO_COMERR,                 SPM_REP, 6,  SPM_PROTO_NONE,                   "%c",                                  "COMErr"},
	{SPM_PROTO_ACTERRORANSWER,         SPM_REP, 6,  SPM_PROTO_NONE,                   "%c 0x00 0x00",                        "ActErrorAnswer"},

	/* SpectroScan */

	{SPM_PROTO_MOVEABSOLUT,            SPM_CMD, 10, SPM_PROTO_ERRORANSWER,            "%c %d %d",                            "MoveAbsolut"},
	{SPM_PROTO_MOVERELATIVE,           SPM_CMD, 10, SPM_PROTO_ERRORANSWER,            "%d %d",                               "MoveRelative"},
	{SPM_PROTO_MOVEHOME,               SPM_CMD, 10, SPM_PROTO_ERRORANSWER,            NULL,                                  "MoveHome"},
	{SPM_PROTO_MOVEUP,                 SPM_CMD, 10, SPM_PROTO_ERRORANSWER,            NULL,                                  "MoveUp"},
	{SPM_PROTO_MOVEDOWN,               SPM_CMD, 10, SPM_PROTO_ERRORANSWER,            NULL,                                  "MoveDown"},
	{SPM_PROTO_MOVETOWHITEREFPOS,      SPM_CMD, 10, SPM_PROTO_ERRORANSWER,            "%c",                                  "MoveToWhiteRefPos"},
	{SPM_PROTO_MOVEANDMEASURE,         SPM_CMD, 10, SPM_PROTO_SPECTRUMANSWER,         "%d %d",                               "MoveAndMeasure"},
/* + SPM_PROTO_ERRORANSER. How is this coped with ???? */

	{SPM_PROTO_OUTPUTACTUALPOSITION,   SPM_CMD, 6,  SPM_PROTO_POSITIONANSWER,         "%c",                                  "OutputActualPosition"},
	{SPM_PROTO_INITIALIZEDEVICE,       SPM_CMD, 18, SPM_PROTO_ERRORANSWER,            NULL,                                  "InitializeDevice"},
	{SPM_PROTO_SCANSPECTROLINO,        SPM_CMD, 18, SPM_PROTO_ERRORANSWER,            NULL,                                  "ScanSpectrolino"},
	{SPM_PROTO_INITMOTORPOSITION,      SPM_CMD, 18, SPM_PROTO_ERRORANSWER,            NULL,                                  "InitMotorPosition"},
	{SPM_PROTO_SETDEVICEONLINE,        SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            NULL,                                  "SetDeviceOnline"},
	{SPM_PROTO_SETDEVICEOFFLINE,       SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            NULL,                                  "SetDeviceOffline"},
	{SPM_PROTO_HOLDPAPER,              SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            NULL,                                  "HoldPaper"},
	{SPM_PROTO_RELEASEPAPER,           SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            NULL,                                  "ReleasePaper"},
	{SPM_PROTO_SETTABLEMODE,           SPM_CMD, 18, SPM_PROTO_ERRORANSWER,            "%c",                                  "SetTableNode"},
	{SPM_PROTO_SETLIGHTLEVEL,          SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            "%c",                                  "SetLightLevel"},
	{SPM_PROTO_SETTRANSMSTANDBYPOS,    SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            "%c %d %d",                            "SetTransmStandbyPos"},
	{SPM_PROTO_SETDIGITIZINGMODE,      SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            NULL,                                  "SetDigitizingMode"},
	{SPM_PROTO_OUTPUTDIGITIZEDVALUES,  SPM_CMD, 6,  SPM_PROTO_POSITIONANSWER,         "%c",                                  "OutpuDigitizedValues"},
	{SPM_PROTO_CHANGEBAUDRATE,         SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            "%c",                                  "ChangeBaudRate"},
	{SPM_PROTO_CHANGEHANDSHAKE,        SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            "%c",                                  "ChangeHandshake"},
	{SPM_PROTO_OUTPUTACTUALKEY,        SPM_CMD, 6,  SPM_PROTO_KEYANSWER,              NULL,                                  "OutputActualKey"},
	{SPM_PROTO_OUTPUTLASTKEY,          SPM_CMD, 6,  SPM_PROTO_KEYANSWER,              NULL,                                  "OutputLastKey"},
	{SPM_PROTO_OUTPUTSTATUS,           SPM_CMD, 6,  SPM_PROTO_STATUSANSWER,           NULL,                                  "OutputStatus"},
	{SPM_PROTO_CLEARSTATUS,            SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            "%c",                                  "ClearStatus"},
	{SPM_PROTO_SETKEYACKNOWLEDGEMODE,  SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            NULL,                                  "SetKeyAcknoledgeMode"},
	{SPM_PROTO_RESETKEYACKNOWLEDGEMODE,SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            NULL,                                  "ResetKeyAcknoledgeMode"},
	{SPM_PROTO_SETSPECIALSTATUS,       SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            "%c",                                  "SetSpecialStatus"},
	{SPM_PROTO_CLEARSPECIALSTATUS,     SPM_CMD, 6,  SPM_PROTO_ERRORANSWER,            "%c",                                  "ClearSpecialStatus"},
	{SPM_PROTO_OUTPUTSPECIALSTATUS,    SPM_CMD, 6,  SPM_PROTO_STATUSANSWER,           NULL,                                  "OutputSpecialStatus"},
	{SPM_PROTO_OUTPUTTYPE,             SPM_CMD, 6,  SPM_PROTO_TYPEANSWER,             NULL,                                  "OutputType"},
	{SPM_PROTO_OUTPUTSERIALNUMBER,     SPM_CMD, 6,  SPM_PROTO_SERIALNUMBERANSWER,     NULL,                                  "OutputSerialNumber"},
	{SPM_PROTO_OUTPUTARTICLENUMBER,    SPM_CMD, 6,  SPM_PROTO_ARTICLENUMBERANSWER,    NULL,                                  "ArticleNumberAnswer"},
	{SPM_PROTO_OUTPUTPRODUCTIONDATE,   SPM_CMD, 6,  SPM_PROTO_PRODUCTIONDATEANSWER,   NULL,                                  "OutputProductionDate"},
	{SPM_PROTO_OUTPUTSOFTWAREVERSION,  SPM_CMD, 6,  SPM_PROTO_SOFTWAREVERSIONANSWER,  NULL,                                  "OutputSoftwareVersion"},
	{SPM_PROTO_ERRORANSWER,            SPM_REP, 6,  SPM_PROTO_NONE,                   "%c",                                  "ErrorAnswer"},
	{SPM_PROTO_POSITIONANSWER,         SPM_REP, 6,  SPM_PROTO_NONE,                   "%c %d %d %d %c",                      "PositionAnswer"},
	{SPM_PROTO_KEYANSWER,              SPM_REP, 6,  SPM_PROTO_NONE,                   "%c %c",                               "KeyAnswer"},
	{SPM_PROTO_STATUSANSWER,           SPM_REP, 6,  SPM_PROTO_NONE,                   "%c",                                  "StatusAnswer"},
	{SPM_PROTO_TYPEANSWER,             SPM_REP, 6,  SPM_PROTO_NONE,                   "%c*18",                               "TypeAnswer"},
	{SPM_PROTO_SERIALNUMBERANSWER,     SPM_REP, 6,  SPM_PROTO_NONE,                   "%d*2",                                "SerialNumberAnswer"},
	{SPM_PROTO_ARTICLENUMBERANSWER,    SPM_REP, 6,  SPM_PROTO_NONE,                   "%c*8",                                "ArticleNumberAnswer"},
	{SPM_PROTO_PRODUCTIONDATEANSWER,   SPM_REP, 6,  SPM_PROTO_NONE,                   "%d %d %d",                            "ProductionDate"},
	{SPM_PROTO_SOFTWAREVERSIONANSWER,  SPM_REP, 6,  SPM_PROTO_NONE,                   "%c*12",                               "SoftwareVersion"},
	{SPM_PROTO_SSCOMERR,               SPM_REP, 6,  SPM_PROTO_NONE,                   "%c",                                  "SpectroScanCOMError"},

	{SPM_PROTO_NONE, SPM_PROTO_NONE, SPM_PROTO_NONE, SPM_PROTO_NONE, NULL, NULL}
};

spm_message* spm_lookup(int cmd)
{
	spm_message* msg = spm_protocol;

	while (msg->cmd != SPM_PROTO_NONE) {
		if (msg->cmd == cmd) return msg;
		msg++;
	}

	return NULL;
}

char *spm_cmdname(int cmd)
{
	spm_message* msg = spm_lookup(cmd);
	if (!msg) return "<unknown>";
	else return msg->desc;
}

int spm_expect_reply(int cmd)
{
	spm_message *msg = spm_lookup(cmd);
	if (!msg) return SPM_PROTO_NONE;
	return msg->reply;
}

int spm_timeout(int cmd)
{
	spm_message *msg = spm_lookup(cmd);
	if (!msg) return SPM_PROTO_NONE;
	return msg->timeout;
}

/***********************************************************************
 * string tokenizer utils
 ***********************************************************************/

typedef struct {
	unsigned char  *buf;
	unsigned char **tokens;
	int    next;
	int    last;
} tokenlist;

tokenlist *tok_start(unsigned char *str)
{
	unsigned char *p;
	tokenlist *tok;

	if (!str) return NULL;

	tok = (tokenlist*) malloc(sizeof(tokenlist));
	tok->buf = strdup(str);
	tok->last = 0;

	/* first pass, count tokens */
	p = strtok(tok->buf, " ");
	while (p) {
		tok->last ++;
		p = strtok(NULL, " ");
	}
	free(tok->buf);

	/* second pass, save tokens */
	if (tok->last) {
		tok->next = 0;
		tok->tokens = (unsigned char **) malloc(tok->last * sizeof(unsigned char*));
		tok->buf = strdup(str);
		p = strtok(tok->buf, " ");
		while (p) {
			tok->tokens[tok->next] = p;
			tok->next ++;
			p = strtok(NULL, " ");
		}
	}
	
	tok->next = 0;

	return tok;
}

unsigned char *tok_next(tokenlist *tok)
{
	if (tok->next < tok->last)
		return tok->tokens[tok->next++];
	else
		return NULL;
}

void tok_end(tokenlist *tok)
{
	free(tok->buf);
	free(tok);
}


/***********************************************************************
 * basic formatting for output
 ***********************************************************************/

int spm_output_begin(unsigned char *buf)
{
	strcpy(buf, ";");
	return 1;
}

int spm_output_end(unsigned char *buf)
{
	strcat(buf, " \r\n");
	return 1;
}

int spm_output_byte(unsigned char *buf, unsigned char b)
{
	unsigned char tmp[80];
	sprintf(tmp, " %d", b);
	strcat(buf, tmp);
	return 1;
}

int spm_output_short(unsigned char *buf, unsigned short s)
{
	unsigned char tmp[80];
	sprintf(tmp, " %d", s);
	strcat(buf, tmp);
	return 1;
}

int spm_output_float(unsigned char *buf, float f)
{
	unsigned char tmp[80];
	sprintf(tmp, " %8.4f", f);
	strcat(buf, tmp);
	return 1;
}

int spm_output_string(unsigned char *buf, unsigned char *s)
{
	strcat(buf, s);
	return 1;
}


/***********************************************************************
 * basic formatting for input
 ***********************************************************************/

int spm_input_begin(unsigned char *buf)
{
	/* gobble leading characters until ":" is found */
	while (buf)
		if (*buf++ == ':') break;
	return (buf && *buf == 0);
}

int spm_input_end(unsigned char *buf)
{
	return (buf == NULL);
}

int spm_input_byte(unsigned char *buf, unsigned char *b)
{
	*b = (unsigned char)atoi(buf);
	return 1;
}

int spm_input_short(unsigned char *buf, unsigned short *s)
{
	*s = (unsigned short)atoi(buf);
	return 1;
}

int spm_input_float(unsigned char *buf, float *f)
{
	*f = (float)atof(buf);
	return 1;
}

int spm_input_string(unsigned char *buf, unsigned char **s)
{
	*s = strdup(buf);
	return 1;
}


/***********************************************************************
 * SPM public API
 ***********************************************************************/

int spm_command(char *buf, int cmd, ...)
{
	int err;
	va_list args;
	va_start(args, cmd);
	err = spm_command(buf, cmd, &args);
	va_end(args);
	return err;
}

int spm_commandv(char *buf, int cmd, va_list *ap)
{
	spm_message *msg;
	tokenlist *tlist;
	char *token;
	int i, toklen, veclen;
	float *arg_floatp, arg_float;
	unsigned char *arg_charp, arg_char;
	unsigned short *arg_shortp, arg_short;

	msg = spm_lookup(cmd);
	if (!msg) {
		spm_error("spm_command: bad command %x\n", cmd);
		return 0;
	}

	spm_output_begin(buf);

	/* output message type */
	if (msg->cmd <= 255)
		spm_output_byte(buf, (unsigned char)msg->cmd);					/* 8 bit only */
	else {
		spm_output_byte(buf, (unsigned char)((msg->cmd & 0xff00) >> 8)); /* high */
		spm_output_byte(buf, (unsigned char)(msg->cmd & 0x00ff));        /* low */
	}

	/* output remaining args */
	if (msg->fmt) {

		tlist = tok_start(msg->fmt);
		while ((token=tok_next(tlist)) != NULL) {

			/* take action based on token:
			 * 0x00 = output literal
			 * %c   = output char
			 * %d   = output int
			 * %f   = output float
			 * %s   = output string
			 * %c*n = output vector of n chars
			 * %d*n = output vector of n ints
			 * %f*n = output vector of n floats
			 */

			toklen = strlen(token);

			if (token[0] == '0') { /* literal */
				arg_char = strtol(token, NULL, 0);
				spm_output_byte(buf,arg_char);
			}
			else if (token[0] == '%' && toklen >= 2) { /* variable */

				if (toklen > 3 && token[2] == '*') { /* vector variable */
					veclen = atoi(&token[3]);
					switch (token[1]) {
					case 'c':
						arg_charp = va_arg(*ap, unsigned char *);
						for (i=0; i<veclen; i++)
							spm_output_byte(buf,arg_charp[i]);
						break;
					case 'd':
						arg_shortp = va_arg(*ap, unsigned short *);
						for (i=0; i<veclen; i++)
							spm_output_short(buf,arg_shortp[i]);
						break;
					case 'f':
						arg_floatp = va_arg(*ap, float *);
						for (i=0; i<veclen; i++)
							spm_output_float(buf,arg_floatp[i]);
						break;
					default:
						spm_error("spm_command: bad vector format type %c\n", token[1]);
						return 0;
					}
				}
				else {
					switch (token[1]) {	/* single variable */
					case 'c':
						arg_char = va_arg(*ap, /* unsigned char */ int);
						spm_output_byte(buf,arg_char);
						break;
					case 'd':
						arg_short = va_arg(*ap, /* unsigned short*/ int);
						spm_output_short(buf,arg_short);
						break;
					case 'f':
						arg_float = va_arg(*ap, double); /* why doesn't float work??? */
						spm_output_float(buf,arg_float);
						break;
					case 's':
						arg_charp = va_arg(*ap, unsigned char *);
						spm_output_string(buf,arg_charp);
						break;
					default:
						spm_error("spm_command: bad format type %c\n", token[1]);
						return 0;
					}
				}
			}
		}
		tok_end(tlist);
	}

	spm_output_end(buf);

#ifdef DEBUG
	fprintf(stderr, "spm_command [%s] %s", msg->desc, buf);
#endif

	return 1;
}

/* Process buffer for expected reply. */
/* return 0 on error */
/* return 1 on not expected reply */
/* return 2 on success */
int spm_reply(char *buf, int reply, ...)
{
	int err;
	va_list args;
	va_start(args, reply);
	err = spm_replyv(buf, reply, &args);
	va_end(args);
	return err;
}

int spm_replyv(char *buf, int reply, va_list *ap)
{
	spm_message *msg;
	tokenlist *tlfmt, *tlbuf;
	char *token;
	int i, toklen, veclen;
	float *arg_floatp, arg_float;
	unsigned char *arg_charp, arg_char, arg_char2;
	unsigned short *arg_shortp, arg_short;

	msg = spm_lookup(reply);
	if (!msg) {
		spm_error("spm_reply: bad reply %x\n", reply);
		return 0;
	}

#ifdef DEBUG
	fprintf(stderr, "spm_reply [%s] %s\n", msg->desc, buf);
#endif

	tlbuf = tok_start(buf);
	if (!spm_input_begin(tok_next(tlbuf))) {
		spm_error("spm_reply: bad input start char\n");
		tok_end(tlbuf);
		return 0;
	}

	/* input message type */
	if (reply <= 255) {
		spm_input_byte(tok_next(tlbuf), &arg_char);
		if (arg_char != reply) {
			tok_end(tlbuf);
			return 1;
		}
	} else {
		spm_input_byte(tok_next(tlbuf), &arg_char);
		spm_input_byte(tok_next(tlbuf), &arg_char2);
		if ( (arg_char != (reply & 0xff00)) && (arg_char2 != (reply & 0x00ff)) ) {
			tok_end(tlbuf);
			return 1;
		}
	}

	/* input remaining args */
	if (msg->fmt) {

		/* walk through format tokens, use them to parse input buf */
		tlfmt = tok_start(msg->fmt);
		while ((token = tok_next(tlfmt)) != NULL) {

			/* take action based on token:
			 * 0x00 = input literal
			 * %c   = input char
			 * %d   = input int
			 * %f   = input float
			 * %s   = input string
			 * %c*n = input vector of n chars
			 * %d*n = input vector of n ints
			 * %f*n = input vector of n floats
			 */

			toklen = strlen(token);

			if (token[0] == '0') { /* literal */
				unsigned char literal = strtol(token, NULL, 0);
				if (!spm_input_byte(tok_next(tlbuf), &arg_char) || arg_char != literal) {
					spm_error("spm_reply: unexpected literal value mismatch\n");
					tok_end(tlbuf);
					tok_end(tlfmt);
					return 0;
				}
			}
			else if (token[0] == '%' && toklen >= 2) { /* variable */

				if (toklen > 3 && token[2] == '*') { /* vector variable */
					veclen = atoi(&token[3]);
					switch (token[1]) {
					case 'c':
						arg_charp = va_arg(*ap, unsigned char *);
						for (i=0; i<veclen; i++)
							spm_input_byte(tok_next(tlbuf),&arg_charp[i]);
						break;
					case 'd':
						arg_shortp = va_arg(*ap, unsigned short *);
						for (i=0; i<veclen; i++)
							spm_input_short(tok_next(tlbuf),&arg_shortp[i]);
						break;
					case 'f':
						arg_floatp = va_arg(*ap, float *);
						for (i=0; i<veclen; i++)
							spm_input_float(tok_next(tlbuf),&arg_floatp[i]);
						break;
					default:
						spm_error("spm_reply: bad vector format type %c\n", token[1]);
						return 0;
					}
				}
				else {
					switch (token[1]) {	/* single variable */
					case 'c':
						arg_charp = va_arg(*ap, unsigned char *);
						spm_input_byte(tok_next(tlbuf),arg_charp);
						break;
					case 'd':
						arg_shortp = va_arg(*ap, unsigned short *);
						spm_input_short(tok_next(tlbuf),arg_shortp);
						break;
					case 'f':
						arg_floatp = va_arg(*ap, float *);
						spm_input_float(tok_next(tlbuf),arg_floatp);
						break;
					case 's':
						arg_charp = va_arg(*ap, unsigned char *);
						spm_input_string(tok_next(tlbuf),&arg_charp);
						break;
					default:
						spm_error("spm_reply: bad format type %c\n", token[1]);
						return 0;
					}
				}
			}
		}
		tok_end(tlfmt);
	}

	if (!spm_input_end(tok_next(tlbuf))) {
		spm_error("spm_reply: bad end-of-input char\n");
		tok_end(tlbuf);
		return 0;
	}

	tok_end(tlbuf);

	return 2;
}

/* Process buffer to return the reply message type */
int spm_reply_type(char *buf)
{
	int  rv;
	tokenlist *tlbuf;
	unsigned char arg_char;

	tlbuf = tok_start(buf);
	if (!spm_input_begin(tok_next(tlbuf))) {
		spm_error("spm_reply: bad input start char\n");
		tok_end(tlbuf);
		return 0;
	}

	spm_input_byte(tok_next(tlbuf), &arg_char);
	rv = arg_char;
	if (arg_char == 0xd0 || arg_char == 0xd1) {
		spm_input_byte(tok_next(tlbuf), &arg_char);
		rv = (rv << 8) | arg_char;
	}
		
	tok_end(tlbuf);

#ifdef DEBUG
	fprintf(stderr, "spm_reply_type %s returning 0x%x\n", buf, rv);
#endif

	return rv;
}


void
spm_error(char *fmt, ...)
{
	va_list args;

	fprintf(stderr,"spm_error: ");
	va_start(args, fmt);
	vfprintf(stderr, fmt, args);
	va_end(args);
	fprintf(stderr, "\n");
	fflush(stdout);
	exit (-1);
}



#ifdef MYDEBUG

main()
{
	char buf[256], *p;
	float fv[] = {2.0, 3.0, 4.0, 5.0};
	float spec[36];
	unsigned char a, b, c, d;
	unsigned short e;
	tokenlist *tok;

	spm_command(buf, SPM_PROTO_EXECMEASUREMENT);
	printf("%s",buf);
	spm_command(buf, SPM_PROTO_EXECREFMEASUREMENT, SPM_MEASUREMENTMODE_WHITECALIBRATION);
	printf("%s",buf);
	spm_command(buf, SPM_PROTO_PARAMETERDOWNLOAD,
				SPM_DSTD_ANSIA,
				SPM_WBASE_ABS,
				SPM_ILLUM_ILLUMINANTA,
				SPM_OBSERVER_TWODEG);
	printf("%s",buf);
	spm_command(buf, SPM_PROTO_RESETSTATUSDOWNLOAD, SPM_STATUSMODE_INITWITHOUTREMOTE);
	printf("%s",buf);
	spm_command(buf, SPM_PROTO_SLOPEDOWNLOAD, fv);
	printf("%s",buf);

	tok = tok_start("thisbuffer.");
	while ((p = tok_next(tok)) != NULL) {
		printf("token: %s\n", p);
	}
	tok_end(tok);

	spm_reply("; 185 9 25 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 21.0 22.0 23.0 24.0 25.0 26.0 27.0 28.0 29.0 30.0 31.0 32.0 33.0 34.0 35.0 36.0 255 254 253 2 12345",
			  SPM_PROTO_SPECPARAMETERANSWER,
			  &a, spec, &b, &c, &d, &e);
	printf("spm_reply = %d %d %d %d %d\n",a,b,c,d,e);

}

#endif

