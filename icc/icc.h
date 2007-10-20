#ifndef ICC_H
#define ICC_H

/* 
 * International Color Consortium Format Library (icclib)
 *
 * Author:  Graeme W. Gill
 * Date:    99/11/29
 * Version: 2.05
 *
 * Copyright 1997 - 2005 Graeme W. Gill
 *
 * This material is licensed with a free use license:-
 * see the License.txt file in this directory for licensing details.
 */

/* We can get some subtle errors if certain headers aren't included */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>

#ifdef __cplusplus
	extern "C" {
#endif

/* Version of icclib release */

#define ICCLIB_VERSION 0x020010
#define ICCLIB_VERSION_STR "2.10"

#undef ENABLE_V4		/* V4 is not fully implemented */

/*
 *  Note XYZ scaling to 1.0, not 100.0
 */

#undef ICC_DEBUG_MALLOC		/* Turns on partial support for filename and linenumber capture */

/* Make allowance for shared library use */
#ifdef ICCLIB_SHARED		/* Compiling or Using shared library version */
# ifdef ICCLIB_EXPORTS		/* Compiling shared library */
#  ifdef NT
#   define ICCLIB_API __declspec(dllexport)
#  endif /* NT */
# else						/* Using shared library */
#  ifdef NT
#   define ICCLIB_API __declspec(dllimport)
#   ifdef ICCLIB_DEBUG
#    pragma comment (lib, "icclibd.lib")
#   else
#    pragma comment (lib, "icclib.lib")
#   endif	/* DEBUG */
#  endif /* NT */
# endif
#else						/* Using static library */
# define ICCLIB_API	/* empty */
#endif


/* ---------------------------------------------- */
/* Platform specific defines */
/* Luckily, the following are largely always true on modern */
/* systems: */
#ifndef INR8
#define INR8   signed char		/* 8 bit signed */
#endif
#ifndef INR16
#define INR16  signed short		/* 16 bit signed */
#endif
#ifndef INR32
#define INR32  signed int		/* 32 bit signed */
#endif
#ifndef ORD8
#define ORD8   unsigned char	/* 8 bit unsigned */
#endif
#ifndef ORD16
#define ORD16  unsigned short	/* 16 bit unsigned */
#endif
#ifndef ORD32
#define ORD32  unsigned int		/* 32 bit unsigned */
#endif

#include "iccV42.h"	/* ICC Version 4.2 definitions. */ 

/* Note that the prefix icm is used for the native Machine */
/* equivalents of the ICC binary file structures, and other icclib */
/* specific definitions. */

/* ---------------------------------------------- */
/* System interface objects. The defaults can be replaced */
/* for adaption to different system environments */

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Heap allocator class interface definition */
#ifdef ICC_DEBUG_MALLOC

#define ICM_ALLOC_BASE																		\
	/* Public: */																			\
																							\
	void *(*dmalloc) (struct _icmAlloc *p, size_t size, char *file, int line);				\
	void *(*dcalloc) (struct _icmAlloc *p, size_t num, size_t size, char *file, int line);	\
	void *(*drealloc)(struct _icmAlloc *p, void *ptr, size_t size, char *file, int line);	\
	void  (*dfree)   (struct _icmAlloc *p, void *ptr, char *file, int line);				\
																							\
	/* we're done with the allocator object */												\
	void (*del)(struct _icmAlloc *p);														\

#ifdef _ICC_C_		/* only inside icc.c */
#define malloc( p, size )	    dmalloc( p, size, __FILE__, __LINE__ )
#define calloc( p, num, size )	dcalloc( p, num, size, __FILE__, __LINE__ )
#define realloc( p, ptr, size )	drealloc( p, ptr, size, __FILE__, __LINE__ )
#define free( p, ptr )	        dfree( p, ptr , __FILE__, __LINE__ )
#endif /* _ICC_C */

#else /* !ICC_DEBUG_MALLOC */

/* Heap allocator class interface definition */
#define ICM_ALLOC_BASE																		\
	/* Public: */																			\
																							\
	void *(*malloc) (struct _icmAlloc *p, size_t size);										\
	void *(*calloc) (struct _icmAlloc *p, size_t num, size_t size);							\
	void *(*realloc)(struct _icmAlloc *p, void *ptr, size_t size);							\
	void  (*free)   (struct _icmAlloc *p, void *ptr);										\
																							\
	/* we're done with the allocator object */												\
	void (*del)(struct _icmAlloc *p);														\

#endif /* !ICC_DEBUG_MALLOC */

/* Common heap allocator interface class */
struct _icmAlloc {
	ICM_ALLOC_BASE
}; typedef struct _icmAlloc icmAlloc;

/* - - - - - - - - - - - - - - - - - - - - -  */

/* Implementation of heap class based on standard system malloc */
struct _icmAllocStd {
	ICM_ALLOC_BASE
}; typedef struct _icmAllocStd icmAllocStd;

/* Create a standard alloc object */
icmAlloc *new_icmAllocStd(void);

/* File access class interface definition */
#define ICM_FILE_BASE																		\
	/* Public: */																			\
																							\
	/* Set current position to offset. Return 0 on success, nz on failure. */				\
	int    (*seek) (struct _icmFile *p, long int offset);									\
																							\
	/* Read count items of size length. Return number of items successfully read. */ 		\
	size_t (*read) (struct _icmFile *p, void *buffer, size_t size, size_t count);			\
																							\
	/* write count items of size length. Return number of items successfully written. */ 	\
	size_t (*write)(struct _icmFile *p, void *buffer, size_t size, size_t count);			\
																							\
	/* printf to the file */																\
	int (*printf)(struct _icmFile *p, const char *format, ...);								\
																							\
	/* flush all write data out to secondary storage. Return nz on failure. */				\
	int (*flush)(struct _icmFile *p);														\
																							\
	/* we're done with the file object, return nz on failure */								\
	int (*del)(struct _icmFile *p);															\


/* Common file interface class */
struct _icmFile {
	ICM_FILE_BASE
}; typedef struct _icmFile icmFile;


/* - - - - - - - - - - - - - - - - - - - - -  */

/* These are avalailable if SEPARATE_STD is not defined: */

/* Implementation of file access class based on standard file I/O */
struct _icmFileStd {
	ICM_FILE_BASE

	/* Private: */
	icmAlloc *al;		/* Heap allocator */
	int      del_al;	/* NZ if heap allocator should be deleted */
	FILE     *fp;
	int   doclose;		/* nz if free should close */
}; typedef struct _icmFileStd icmFileStd;

/* Create given a file name */
icmFile *new_icmFileStd_name(char *name, char *mode);

/* Create given a (binary) FILE* */
icmFile *new_icmFileStd_fp(FILE *fp);

/* Create given a file name and allocator */
icmFile *new_icmFileStd_name_a(char *name, char *mode, icmAlloc *al);

/* Create given a (binary) FILE* and allocator */
icmFile *new_icmFileStd_fp_a(FILE *fp, icmAlloc *al);


/* - - - - - - - - - - - - - - - - - - - - -  */
/* Implementation of file access class based on a memory image */
struct _icmFileMem {
	ICM_FILE_BASE

	/* Private: */
	icmAlloc *al;		/* Heap allocator */
	int      del_al;	/* NZ if heap allocator should be deleted */
	unsigned char *start, *cur, *end;

}; typedef struct _icmFileMem icmFileMem;

/* Create a memory image file access class with allocator */
icmFile *new_icmFileMem_a(void *base, size_t length, icmAlloc *al);

/* This is avalailable if SEPARATE_STD is not defined: */

/* Create a memory image file access class */
icmFile *new_icmFileMem(void *base, size_t length);


/* --------------------------------- */
/* Assumed constants                 */

#define MAX_CHAN 15		/* Maximum number of color channels */

/* --------------------------------- */
/* tag and other compound structures */

typedef int icmSig;	/* Otherwise un-enumerated 4 byte signature */

/* Non-standard Color Space Signatures - will be incompatible outside icclib! */

/* A monochrome CIE L* space */
#define icmSigLData ((icColorSpaceSignature) icmMakeTag('L',' ',' ',' '))

/* A monochrome CIE Y space */
#define icmSigYData ((icColorSpaceSignature) icmMakeTag('Y',' ',' ',' '))


/* Pseudo Color Space Signatures - just used within icclib */

/* Pseudo PCS colospace of profile */
#define icmSigPCSData ((icColorSpaceSignature) icmMakeTag('P','C','S',' '))

/* Pseudo PCS colospace to signal 8 bit Lab */
#define icmSigLab8Data ((icColorSpaceSignature) icmMakeTag('L','a','b','8'))

/* Pseudo PCS colospace to signal V2 16 bit Lab */
#define icmSigLabV2Data ((icColorSpaceSignature) icmMakeTag('L','a','b','2'))

/* Pseudo PCS colospace to signal V4 16 bit Lab */
#define icmSigLabV4Data ((icColorSpaceSignature) icmMakeTag('L','a','b','4'))

/* Pseudo PCS colospace to signal 8 bit L */
#define icmSigL8Data ((icColorSpaceSignature) icmMakeTag('L',' ',' ','8'))

/* Pseudo PCS colospace to signal V2 16 bit L */
#define icmSigLV2Data ((icColorSpaceSignature) icmMakeTag('L',' ',' ','2'))

/* Pseudo PCS colospace to signal V4 16 bit L */
#define icmSigLV4Data ((icColorSpaceSignature) icmMakeTag('L',' ',' ','4'))

/* Non-standard Platform Signature */
#define icmSig_nix ((icPlatformSignature) icmMakeTag('*','n','i','x'))


/* Internal Tag Type signature, used to handle any tag that */
/* isn't handled with a specific type object. */
#define icmSigUnknownType ((icTagTypeSignature)icmMakeTag('?','?','?','?'))

typedef struct {
	ORD32 l;			/* High and low components of signed 64 bit */
	INR32 h;
} icmInt64;

typedef struct {
	ORD32 l,h;			/* High and low components of unsigned 64 bit */
} icmUint64;

/* XYZ Number */
typedef struct {
    double  X;
    double  Y;
    double  Z;
} icmXYZNumber;

/* Response 16 number */
typedef struct {
	double	       deviceValue;	/* The device value in range 0.0 - 1.0 */
	double	       measurement;	/* The reading value */
} icmResponse16Number;

/*
 *  read and write method error codes:
 *  0 = sucess
 *  1 = file format/logistical error
 *  2 = system error
 */

#define ICM_BASE_MEMBERS																\
	/* Private: */																		\
	icTagTypeSignature  ttype;		/* The tag type signature */						\
	struct _icc    *icp;			/* Pointer to ICC we're a part of */				\
	int	           touched;			/* Flag for write bookeeping */						\
    int            refcount;		/* Reference count for sharing */					\
	unsigned int   (*get_size)(struct _icmBase *p);										\
	int            (*read)(struct _icmBase *p, unsigned long len, unsigned long of);	\
	int            (*write)(struct _icmBase *p, unsigned long of);						\
	void           (*del)(struct _icmBase *p);											\
																						\
	/* Public: */																		\
	void           (*dump)(struct _icmBase *p, icmFile *op, int verb);					\
	int            (*allocate)(struct _icmBase *p);									

/* Base tag element data object */
struct _icmBase {
	ICM_BASE_MEMBERS
}; typedef struct _icmBase icmBase;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Tag type to hold an unknown tag type, */
/* so that it can be read and copied. */
struct _icmUnknown {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int      _size;	/* size of data currently allocated */

	/* Public: */
	icTagTypeSignature uttype;	/* The unknown tag type signature */
	unsigned int	  size;		/* Allocated and used size of the array */
    unsigned char	  *data;  	/* tag data */
}; typedef struct _icmUnknown icmUnknown;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* UInt8 Array */
struct _icmUInt8Array {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _size;		/* Size currently allocated */

	/* Public: */
	unsigned long	size;		/* Allocated and used size of the array */
    unsigned int   *data;		/* Pointer to array of data */ 
}; typedef struct _icmUInt8Array icmUInt8Array;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* uInt16 Array */
struct _icmUInt16Array {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _size;		/* Size currently allocated */

	/* Public: */
	unsigned long	size;		/* Allocated and used size of the array */
    unsigned int	*data;		/* Pointer to array of data */ 
}; typedef struct _icmUInt16Array icmUInt16Array;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* uInt32 Array */
struct _icmUInt32Array {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _size;		/* Size currently allocated */

	/* Public: */
	unsigned long	size;		/* Allocated and used size of the array */
    unsigned int	*data;		/* Pointer to array of data */ 
}; typedef struct _icmUInt32Array icmUInt32Array;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* UInt64 Array */
struct _icmUInt64Array {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _size;		/* Size currently allocated */

	/* Public: */
	unsigned long	size;		/* Allocated and used size of the array */
    icmUint64		*data;		/* Pointer to array of hight data */ 
}; typedef struct _icmUInt64Array icmUInt64Array;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* u16Fixed16 Array */
struct _icmU16Fixed16Array {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _size;		/* Size currently allocated */

	/* Public: */
	unsigned long	size;		/* Allocated and used size of the array */
    double			*data;		/* Pointer to array of hight data */ 
}; typedef struct _icmU16Fixed16Array icmU16Fixed16Array;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* s15Fixed16 Array */
struct _icmS15Fixed16Array {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _size;		/* Size currently allocated */

	/* Public: */
	unsigned long	size;		/* Allocated and used size of the array */
    double			*data;		/* Pointer to array of hight data */ 
}; typedef struct _icmS15Fixed16Array icmS15Fixed16Array;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* XYZ Array */
struct _icmXYZArray {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _size;		/* Size currently allocated */

	/* Public: */
	unsigned long	size;		/* Allocated and used size of the array */
    icmXYZNumber	*data;		/* Pointer to array of data */ 
}; typedef struct _icmXYZArray icmXYZArray;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Curve */
typedef enum {
    icmCurveUndef           = -1, /* Undefined curve */
    icmCurveLin             = 0,  /* Linear transfer curve */
    icmCurveGamma           = 1,  /* Gamma power transfer curve */
    icmCurveSpec            = 2   /* Specified curve */
} icmCurveStyle;

/* Curve reverse lookup information */
typedef struct {
	int inited;				/* Flag */
	double rmin, rmax;		/* Range of reverse grid */
	double qscale;			/* Quantising scale factor */
	long rsize;				/* Number of reverse lists */
	unsigned int **rlists;	/* Array of list of fwd values that may contain output value */
							/* Offset 0 = allocated size */
							/* Offset 1 = next free index */
							/* Offset 2 = first fwd index */
	unsigned long size;		/* Copy of forward table size */
	double       *data;		/* Copy of forward table data */
} icmRevTable;

struct _icmCurve {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _size;		/* Size currently allocated */
	icmRevTable  rt;			/* Reverse table information */

	/* Public: */
    icmCurveStyle   flag;		/* Style of curve */
	unsigned long	size;		/* Allocated and used size of the array */
    double         *data;  		/* Curve data scaled to range 0.0 - 1.0 */
								/* or data[0] = gamma value */
	/* Translate a value through the curve, return warning flags */
	int (*lookup_fwd) (struct _icmCurve *p, double *out, double *in);	/* Forwards */
	int (*lookup_bwd) (struct _icmCurve *p, double *out, double *in);	/* Backwards */

}; typedef struct _icmCurve icmCurve;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Data */
typedef enum {
    icmDataUndef           = -1, /* Undefined data */
    icmDataASCII           = 0,  /* ASCII data */
    icmDataBin             = 1   /* Binary data */
} icmDataStyle;

struct _icmData {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _size;		/* Size currently allocated */

	/* Public: */
    icmDataStyle	flag;		/* Style of data */
	unsigned long	size;		/* Allocated and used size of the array (inc ascii null) */
    unsigned char	*data;  	/* data or string, NULL if size == 0 */
}; typedef struct _icmData icmData;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* text */
struct _icmText {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _size;		/* Size currently allocated */

	/* Public: */
	unsigned long	 size;		/* Allocated and used size of data, inc null */
	char             *data;		/* ascii string (null terminated), NULL if size==0 */

}; typedef struct _icmText icmText;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* The base date time number */
struct _icmDateTimeNumber {
	ICM_BASE_MEMBERS

	/* Public: */
    unsigned int      year;
    unsigned int      month;
    unsigned int      day;
    unsigned int      hours;
    unsigned int      minutes;
    unsigned int      seconds;
}; typedef struct _icmDateTimeNumber icmDateTimeNumber;

#ifdef NEW
/* - - - - - - - - - - - - - - - - - - - - -  */
/* DeviceSettings */

/*
   I think this all works like this:

Valid setting = (   (platform == platform1 and platform1.valid)
                 or (platform == platform2 and platform2.valid)
                 or ...
                )

where
	platformN.valid = (   platformN.combination1.valid
	                   or platformN.combination2.valid
	                   or ...
	                  )

where
	platformN.combinationM.valid = (    platformN.combinationM.settingstruct1.valid
	                                and platformN.combinationM.settingstruct2.valid
	                                and ...
	                               )

where
	platformN.combinationM.settingstructP.valid = (   platformN.combinationM.settingstructP.setting1.valid
	                                               or platformN.combinationM.settingstructP.setting2.valid
	                                               or ...
	                                              )

 */

/* The Settings Structure holds an array of settings of a particular type */
struct _icmSettingStruct {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _num;				/* Size currently allocated */

	/* Public: */
	icSettingsSig       settingSig;		/* Setting identification */
	unsigned long       numSettings; 	/* number of setting values */
	union {								/* Setting values - type depends on Sig */
		icUInt64Number      *resolution;
		icDeviceMedia       *media;
		icDeviceDither      *halftone;
	}
}; typedef struct _icmSettingStruct icmSettingStruct;

/* A Setting Combination holds all arrays of different setting types */
struct _icmSettingComb {
	/* Private: */
	unsigned int   _num;			/* number currently allocated */

	/* Public: */
	unsigned long       numStructs;   /* num of setting structures */
	icmSettingStruct    *data;
}; typedef struct _icmSettingComb icmSettingComb;

/* A Platform Entry holds all setting combinations */
struct _icmPlatformEntry {
	/* Private: */
	unsigned int   _num;			/* number currently allocated */

	/* Public: */
	icPlatformSignature platform;
	unsigned long       numCombinations;    /* num of settings and allocated array size */
	icmSettingComb      *data; 
}; typedef struct _icmPlatformEntry icmPlatformEntry;

/* The Device Settings holds all platform settings */
struct _icmDeviceSettings {
	/* Private: */
	unsigned int   _num;			/* number currently allocated */

	/* Public: */
	unsigned long       numPlatforms;	/* num of platforms and allocated array size */
	icmPlatformEntry    *data;			/* Array of pointers to platform entry data */
}; typedef struct _icmDeviceSettings icmDeviceSettings;

#endif /* NEW */

/* - - - - - - - - - - - - - - - - - - - - -  */
/* optimised simplex Lut lookup axis flipping structure */
typedef struct {
	double fth;		/* Flip threshold */
	char bthff;		/* Below threshold flip flag value */
	char athff;		/* Above threshold flip flag value */
} sx_flip_info;
 
/* Set method flags */
#define ICM_CLUT_SET_EXACT 0x0000	/* Set clut node values exactly from callback */
#define ICM_CLUT_SET_APXLS 0x0001	/* Set clut node values to aproximate least squares fit */

/* lut */
struct _icmLut {
	ICM_BASE_MEMBERS

	/* Private: */
	/* Cache appropriate normalization routines */
	int dinc[MAX_CHAN];				/* Dimensional increment through clut (in doubles) */
	int dcube[1 << MAX_CHAN];		/* Hyper cube offsets (in doubles) */
	icmRevTable rit[MAX_CHAN];		/* Reverse input table information */
	icmRevTable rot[MAX_CHAN];		/* Reverse output table information */
	sx_flip_info finfo[MAX_CHAN];	/* Optimised simplex flip information */

	unsigned int inputTable_size;	/* size allocated to input table */
	unsigned int clutTable_size;	/* size allocated to clut table */
	unsigned int outputTable_size;	/* size allocated to output table */

	/* Optimised simplex orientation information. oso_ffa is NZ if valid. */
	/* Only valid if inputChan > 1 && clutPoints > 1 */
									/* oso_ff[01] will be NULL if not valid */
	unsigned short *oso_ffa;		/* Flip flags for inputChan-1 dimensions, organised */
									/* [inputChan 0, 0..cp-2]..[inputChan ic-2, 0..cp-2] */
	unsigned short *oso_ffb;		/* Flip flags for dimemension inputChan-1, organised */
									/* [0..cp-2] */
	int odinc[MAX_CHAN];			/* Dimensional increment through oso_ffa */
	
	/* return the minimum and maximum values of the given channel in the clut */
	void (*min_max) (struct _icmLut *pp, double *minv, double *maxv, int chan);

	/* Translate color values through 3x3 matrix, input tables only, multi-dimensional lut, */
	/* or output tables, */
	int (*lookup_matrix)  (struct _icmLut *pp, double *out, double *in);
	int (*lookup_input)   (struct _icmLut *pp, double *out, double *in);
	int (*lookup_clut_nl) (struct _icmLut *pp, double *out, double *in);
	int (*lookup_clut_sx) (struct _icmLut *pp, double *out, double *in);
	int (*lookup_output)  (struct _icmLut *pp, double *out, double *in);

	/* Public: */

	/* return non zero if matrix is non-unity */
	int (*nu_matrix) (struct _icmLut *pp);

    unsigned int	inputChan;      /* Num of input channels */
    unsigned int	outputChan;     /* Num of output channels */
    unsigned int	clutPoints;     /* Num of grid points */
    unsigned int	inputEnt;       /* Num of in-table entries (must be 256 for Lut8) */
    unsigned int	outputEnt;      /* Num of out-table entries (must be 256 for Lut8) */
    double			e[3][3];		/* 3 * 3 array */
	double	        *inputTable;	/* The in-table: [inputChan * inputEnt] */
	double	        *clutTable;		/* The clut: [(clutPoints ^ inputChan) * outputChan] */
	double	        *outputTable;	/* The out-table: [outputChan * outputEnt] */
	/* inputTable  is organized [inputChan 0..ic-1][inputEnt 0..ie-1] */
	/* clutTable   is organized [inputChan 0, 0..cp-1]..[inputChan ic-1, 0..cp-1]
	                                                                [outputChan 0..oc-1] */
	/* outputTable is organized [outputChan 0..oc-1][outputEnt 0..oe-1] */

	/* Helper function to setup a Lut tables contents */
	int (*set_tables) (
		struct _icmLut *p,						/* Pointer to Lut object */
		int     flags,							/* Setting flags */
		void   *cbctx,							/* Opaque callback context pointer value */
		icColorSpaceSignature insig, 			/* Input color space */
		icColorSpaceSignature outsig, 			/* Output color space */
		void (*infunc)(void *cbctx, double *out, double *in),
								/* Input transfer function, inspace->inspace' (NULL = default) */
		double *inmin, double *inmax,			/* Maximum range of inspace' values */
												/* (NULL = default) */
		void (*clutfunc)(void *cbntx, double *out, double *in),
								/* inspace' -> outspace' transfer function */
		double *clutmin, double *clutmax,		/* Maximum range of outspace' values */
												/* (NULL = default) */
		void (*outfunc)(void *cbntx, double *out, double *in));
								/* Output transfer function, outspace'->outspace (NULL = deflt) */
		
}; typedef struct _icmLut icmLut;

/* Helper function to set multiple Lut tables simultaneously */
/* Note that these tables all have to be compatible in */
/* having the same configuration and resolutions. */
/* Set errc and return error number in underlying icc */
int icmSetMultiLutTables (
	int ntables,							/* Number of tables to be set, 1..3 */
	struct _icmLut *p[3],					/* Pointer to Lut object */
	int     flags,							/* Setting flags */
	void   *cbctx,							/* Opaque callback context pointer value */
	icColorSpaceSignature insig, 			/* Input color space */
	icColorSpaceSignature outsig, 			/* Output color space */
	void (*infunc)(void *cbctx, double *out, double *in),
							/* Input transfer function, inspace->inspace' (NULL = default) */
							/* Will be called ntables times each input grid value */
	double *inmin, double *inmax,			/* Maximum range of inspace' values */
											/* (NULL = default) */
	void (*clutfunc)(void *cbntx, double *out, double *in),
							/* inspace' -> outspace[ntables]' transfer function */
							/* will be called once for each input' grid value, and */
							/* ntables output values should be written consecutively */
							/* to out[]. */
	double *clutmin, double *clutmax,		/* Maximum range of outspace' values */
											/* (NULL = default) */
	void (*outfunc)(void *cbntx, double *out, double *in)
								/* Output transfer function, outspace'->outspace (NULL = deflt) */
								/* Will be called ntables times on each output value */
);
		

/* Helper function to set multiple Lut tables simultaneously */
/* Note that these tables all have to be compatible in */
/* having the same configuration and resolutions. */
/* Set errc and return error number in underlying icc */
/* Version that assumes independence of input and output curves */
int icmSetMultiLutTables2 (
	int ntables,							/* Number of tables to be set, 1..n */
	struct _icmLut *p[3],					/* Pointer to Lut objects */
	int     flags,							/* Setting flags */
	void   *cbctx,							/* Opaque callback context pointer value */
	icColorSpaceSignature insig, 			/* Input color space */
	icColorSpaceSignature outsig, 			/* Output color space */
	void (*infunc)(void *cbctx, double *out, double *in, int tab),
							/* Input transfer function, inspace->inspace' (NULL = default) */
	double *inmin[MAX_CHAN], double *inmax[MAX_CHAN],
							/* Maximum range of inspace' values (NULL = default) */
	void (*clutfunc)(void *cbntx, double *out, double *in, int tab),
							/* inspace' -> outspace[ntables]' transfer function */
							/* will be called once for each input' grid value */
	double *clutmin[MAX_CHAN], double *clutmax[MAX_CHAN],	
							/* Maximum range of outspace' values (NULL = default) */
	void (*outfunc)(void *cbntx, double *out, double *in, int tab)
								/* Output transfer function, outspace'->outspace (NULL = deflt) */
);
		

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Measurement Data */
struct _icmMeasurement {
	ICM_BASE_MEMBERS

	/* Public: */
    icStandardObserver           observer;       /* Standard observer */
    icmXYZNumber                 backing;        /* XYZ for backing */
    icMeasurementGeometry        geometry;       /* Meas. geometry */
    double                       flare;          /* Measurement flare */
    icIlluminant                 illuminant;     /* Illuminant */
}; typedef struct _icmMeasurement icmMeasurement;

#ifdef NEW
/* - - - - - - - - - - - - - - - - - - - - -  */
/* MultiLocalizedUnicode */

struct _icmMultiLocalizedUnicode {
	ICM_BASE_MEMBERS

	/* Private: */
	~~~999

	/* Public: */
	~~~999

}; typedef struct _icmMultiLocalizedUnicode icmMultiLocalizedUnicode;

#endif /* NEW */

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Named color */

/* Structure that holds each named color data */
typedef struct {
	struct _icc      *icp;				/* Pointer to ICC we're a part of */
	char              root[32];			/* Root name for color */
	double            pcsCoords[3];		/* icmNC2: PCS coords of color */
	double            deviceCoords[MAX_CHAN];	/* Dev coords of color */
} icmNamedColorVal;

struct _icmNamedColor {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int      _count;			/* Count currently allocated */

	/* Public: */
    unsigned int      vendorFlag;		/* Bottom 16 bits for IC use */
    unsigned int      count;			/* Count of named colors */
    unsigned int      nDeviceCoords;	/* Num of device coordinates */
    char              prefix[32];		/* Prefix for each color name (null terminated) */
    char              suffix[32];		/* Suffix for each color name (null terminated) */
    icmNamedColorVal  *data;			/* Array of [count] color values */
}; typedef struct _icmNamedColor icmNamedColor;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Colorant table */
/* (Contribution from Piet Vandenborre, derived from NamedColor) */

/* Structure that holds colorant table data */
typedef struct {
	struct _icc         *icp;			/* Pointer to ICC we're a part of */
	char                 name[32];		/* Name for colorant */
	double               pcsCoords[3];	/* PCS coords of colorant */
} icmColorantTableVal;

struct _icmColorantTable {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int         _count;		/* Count currently allocated */

	/* Public: */
    unsigned int         count;			/* Count of colorants */
    icmColorantTableVal *data;			/* Array of [count] colorants */
}; typedef struct _icmColorantTable icmColorantTable;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* textDescription */
struct _icmTextDescription {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned long  _size;			/* Size currently allocated */
	unsigned long  uc_size;			/* uc Size currently allocated */
	int            (*core_read)(struct _icmTextDescription *p, char **bpp, char *end);
	int            (*core_write)(struct _icmTextDescription *p, char **bpp);

	/* Public: */
	unsigned long	  size;			/* Allocated and used size of desc, inc null */
	char              *desc;		/* ascii string (null terminated) */

	unsigned int      ucLangCode;	/* UniCode language code */
	unsigned long	  ucSize;		/* Allocated and used size of ucDesc in wchars, inc null */
	ORD16             *ucDesc;		/* The UniCode description (null terminated) */

	ORD16             scCode;		/* ScriptCode code */
	unsigned long	  scSize;		/* Used size of scDesc in bytes, inc null */
	ORD8              scDesc[67];	/* ScriptCode Description (null terminated, max 67) */
}; typedef struct _icmTextDescription icmTextDescription;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Profile sequence structure */
struct _icmDescStruct {
	/* Private: */
	struct _icc      *icp;				/* Pointer to ICC we're a part of */

	/* Public: */
	int             (*allocate)(struct _icmDescStruct *p);	/* Allocate method */
    icmSig            deviceMfg;		/* Dev Manufacturer */
    unsigned int      deviceModel;		/* Dev Model */
    icmUint64         attributes;		/* Dev attributes */
    icTechnologySignature technology;	/* Technology sig */
	icmTextDescription device;			/* Manufacturer text (sub structure) */
	icmTextDescription model;			/* Model text (sub structure) */
}; typedef struct _icmDescStruct icmDescStruct;

/* Profile sequence description */
struct _icmProfileSequenceDesc {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int	 _count;			/* number currently allocated */

	/* Public: */
    unsigned int      count;			/* Number of descriptions */
	icmDescStruct     *data;			/* array of [count] descriptions */
}; typedef struct _icmProfileSequenceDesc icmProfileSequenceDesc;

#ifdef NEW
/* - - - - - - - - - - - - - - - - - - - - -  */
/* ResponseCurveSet16 */

struct _icmResponseCurveSet16 {
	ICM_BASE_MEMBERS

	/* Private: */
	~~~999

	/* Public: */
	~~~999

}; typedef struct _icmResponseCurveSet16 icmResponseCurveSet16;

#endif /* NEW */

/* - - - - - - - - - - - - - - - - - - - - -  */
/* signature (only ever used for technology ??) */
struct _icmSignature {
	ICM_BASE_MEMBERS

	/* Public: */
    icTechnologySignature sig;	/* Signature */
}; typedef struct _icmSignature icmSignature;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Per channel Screening Data */
typedef struct {
	/* Public: */
    double            frequency;		/* Frequency */
    double            angle;			/* Screen angle */
    icSpotShape       spotShape;		/* Spot Shape encodings below */
} icmScreeningData;

struct _icmScreening {
	ICM_BASE_MEMBERS

	/* Private: */
	unsigned int   _channels;			/* number currently allocated */

	/* Public: */
    unsigned int      screeningFlag;	/* Screening flag */
    unsigned int      channels;			/* Number of channels */
    icmScreeningData  *data;			/* Array of screening data */
}; typedef struct _icmScreening icmScreening;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Under color removal, black generation */
struct _icmUcrBg {
	ICM_BASE_MEMBERS

	/* Private: */
    unsigned int      UCR_count;		/* Currently allocated UCR count */
    unsigned int      BG_count;			/* Currently allocated BG count */
	unsigned long	  _size;			/* Currently allocated string size */

	/* Public: */
    unsigned int      UCRcount;			/* Undercolor Removal Curve length */
    double           *UCRcurve;		    /* The array of UCR curve values, 0.0 - 1.0 */
										/* or 0.0 - 100 % if count = 1 */
    unsigned int      BGcount;			/* Black generation Curve length */
    double           *BGcurve;			/* The array of BG curve values, 0.0 - 1.0 */
										/* or 0.0 - 100 % if count = 1 */
	unsigned long	  size;				/* Allocated and used size of desc, inc null */
	char              *string;			/* UcrBg description (null terminated) */
}; typedef struct _icmUcrBg icmUcrBg;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* viewingConditionsType */
struct _icmViewingConditions {
	ICM_BASE_MEMBERS

	/* Public: */
    icmXYZNumber    illuminant;		/* In candelas per sq. meter */
    icmXYZNumber    surround;		/* In candelas per sq. meter */
    icIlluminant    stdIlluminant;	/* See icIlluminant defines */
}; typedef struct _icmViewingConditions icmViewingConditions;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Postscript Color Rendering Dictionary names type */
struct _icmCrdInfo {
	ICM_BASE_MEMBERS
	/* Private: */
    unsigned long    _ppsize;		/* Currently allocated size */
	unsigned long    _crdsize[4];	/* Currently allocated sizes */

	/* Public: */
    unsigned long    ppsize;		/* Postscript product name size (including null) */
    char            *ppname;		/* Postscript product name (null terminated) */
	unsigned long    crdsize[4];	/* Rendering intent 0-3 CRD names sizes (icluding null) */
	char            *crdname[4];	/* Rendering intent 0-3 CRD names (null terminated) */
}; typedef struct _icmCrdInfo icmCrdInfo;

/* - - - - - - - - - - - - - - - - - - - - -  */
/* Apple ColorSync 2.5 video card gamma type (vcgt) */
struct _icmVideoCardGammaTable {
	unsigned short   channels;		/* No of gamma channels (1 or 3) */
	unsigned short   entryCount; 	/* Number of entries per channel */
	unsigned short   entrySize;		/* Size in bytes of each entry */
	void            *data;			/* Variable size data */
}; typedef struct _icmVideoCardGammaTable icmVideoCardGammaTable;

struct _icmVideoCardGammaFormula {
	unsigned short   channels;		/* No of gamma channels (always 3) */
	double           redGamma;		/* must be > 0.0 */
	double           redMin;		/* must be > 0.0 and < 1.0 */
	double           redMax;		/* must be > 0.0 and < 1.0 */
	double           greenGamma;   	/* must be > 0.0 */
	double           greenMin;		/* must be > 0.0 and < 1.0 */
	double           greenMax;		/* must be > 0.0 and < 1.0 */
	double           blueGamma;		/* must be > 0.0 */
	double           blueMin;		/* must be > 0.0 and < 1.0 */
	double           blueMax;		/* must be > 0.0 and < 1.0 */
}; typedef struct _icmVideoCardGammaFormula icmVideoCardGammaFormula;

typedef enum {
	icmVideoCardGammaTableType = 0,
	icmVideoCardGammaFormulaType = 1
} icmVideoCardGammaTagType;

struct _icmVideoCardGamma {
	ICM_BASE_MEMBERS
	icmVideoCardGammaTagType tagType;	/* eg. table or formula, use above enum */
	union {
		icmVideoCardGammaTable   table;
		icmVideoCardGammaFormula formula;
	} u;

	double (*lookup)(struct _icmVideoCardGamma *p, int chan, double iv); /* Read a value */

}; typedef struct _icmVideoCardGamma icmVideoCardGamma;

/* ------------------------------------------------- */
/* The Profile header */
struct _icmHeader {

  /* Private: */
	unsigned int           (*get_size)(struct _icmHeader *p);
	int                    (*read)(struct _icmHeader *p, unsigned long len, unsigned long of);
	int                    (*write)(struct _icmHeader *p, unsigned long of, int doid);
	void                   (*del)(struct _icmHeader *p);
	struct _icc            *icp;			/* Pointer to ICC we're a part of */
    unsigned int            size;			/* Profile size in bytes */

  /* public: */
	void                   (*dump)(struct _icmHeader *p, icmFile *op, int verb);

	/* Values that must be set before writing */
    icProfileClassSignature deviceClass;	/* Type of profile */
    icColorSpaceSignature   colorSpace;		/* Clr space of data (Cannonical input space) */
    icColorSpaceSignature   pcs;			/* PCS: XYZ or Lab (Cannonical output space) */
    icRenderingIntent       renderingIntent;/* Rendering intent */

	/* Values that should be set before writing */
    icmSig                  manufacturer;	/* Dev manufacturer */
    icmSig		            model;			/* Dev model */
    icmUint64               attributes;		/* Device attributes */
    unsigned int            flags;			/* Various bits */

	/* Values that may optionally be set before writing */
    /* icmUint64            attributes;		   Device attributes.h (see above) */
    icmSig                  creator;		/* Profile creator */

	/* Values that are not normally set, since they have defaults */
    icmSig                  cmmId;			/* CMM for profile */
    int            			majv, minv, bfv;/* Format version - major, minor, bug fix */
											/* Default is majv = 2, alternate value is 4 */
    icmDateTimeNumber       date;			/* Creation Date */
    icPlatformSignature     platform;		/* Primary Platform */
    icmXYZNumber            illuminant;		/* Profile illuminant */

	/* Values that are created automatically */
    unsigned char           id[16];			/* MD5 fingerprint value, lsB to msB  <V4.0+> */

}; typedef struct _icmHeader icmHeader;

/* ---------------------------------------------------------- */
/* Objects for accessing lookup functions. */
/* Note that the "effective" PCS colorspace is the one specified by the */
/* PCS override, and visible in the overal profile conversion. */
/* The "native" PCS colorspace is the one that the underlying tags actually */
/* represent, and the PCS that the individual stages within each profile type handle. */
/* The conversion between the native and effective PCS is done in the to/from_abs() */
/* conversions. */

/* Public: Parameter to get_luobj function */
typedef enum {
    icmFwd           = 0,  /* Device to PCS, or Device 1 to Last Device */
    icmBwd           = 1,  /* PCS to Device, or Last Device to Device */
    icmGamut         = 2,  /* PCS Gamut check */
    icmPreview       = 3   /* PCS to PCS preview */
} icmLookupFunc;

/* Public: Parameter to get_luobj function */
typedef enum {
    icmLuOrdNorm     = 0,  /* Normal profile preference: Lut, matrix, monochrome */
    icmLuOrdRev      = 1   /* Reverse profile preference: monochrome, matrix, monochrome */
} icmLookupOrder;

/* Public: Lookup algorithm object type */
typedef enum {
    icmMonoFwdType       = 0,	/* Monochrome, Forward */
    icmMonoBwdType       = 1,	/* Monochrome, Backward */
    icmMatrixFwdType     = 2,	/* Matrix, Forward */
    icmMatrixBwdType     = 3,	/* Matrix, Backward */
    icmLutType           = 4,	/* Multi-dimensional Lookup Table */
    icmNamedType         = 5	/* Named color data */
} icmLuAlgType;

/* Lookup class members common to named and non-named color types */
#define LU_ICM_BASE_MEMBERS																\
	/* Private: */																		\
	icmLuAlgType   ttype;		    	/* The object tag */							\
	struct _icc    *icp;				/* Pointer to ICC we're a part of */			\
	icRenderingIntent intent;			/* Effective (externaly visible) intent */		\
	icmLookupFunc function;				/* Functionality being used */					\
	icmXYZNumber pcswht, whitePoint, blackPoint;	/* White and black point info */	\
	double toAbs[3][3];					/* Matrix to convert from relative to absolute */ \
	double fromAbs[3][3];				/* Matrix to convert from absolute to relative */ \
    icColorSpaceSignature inSpace;		/* Native Clr space of input */					\
    icColorSpaceSignature outSpace;		/* Native Clr space of output */				\
	icColorSpaceSignature pcs;			/* Native PCS */								\
    icColorSpaceSignature e_inSpace;	/* Effective Clr space of input */				\
    icColorSpaceSignature e_outSpace;	/* Effective Clr space of output */				\
	icColorSpaceSignature e_pcs;		/* Effective PCS */								\
																						\
	/* Public: */																		\
	void           (*del)(struct _icmLuBase *p);										\
																						\
	               /* Internal native colorspaces: */	     							\
	void           (*lutspaces) (struct _icmLuBase *p, icColorSpaceSignature *ins, int *inn,	\
	                                                   icColorSpaceSignature *outs, int *outn,	\
	                                                   icColorSpaceSignature *pcs);			\
	                                                                                        \
	               /* External effecive colorspaces */										\
	void           (*spaces) (struct _icmLuBase *p, icColorSpaceSignature *ins, int *inn,	\
	                                 icColorSpaceSignature *outs, int *outn,				\
	                                 icmLuAlgType *alg, icRenderingIntent *intt, 			\
	                                 icmLookupFunc *fnc, icColorSpaceSignature *pcs); 		\
																							\
	               /* Relative to Absolute for this WP in XYZ */   							\
	void           (*XYZ_Rel2Abs)(struct _icmLuBase *p, double *xyzout, double *xyzin);		\
																							\
	               /* Absolute to Relative for this WP in XYZ */   							\
	void           (*XYZ_Abs2Rel)(struct _icmLuBase *p, double *xyzout, double *xyzin);		\


/* Non-algorithm specific lookup class. Used as base class of algorithm specific class. */
#define LU_ICM_NN_BASE_MEMBERS															\
    LU_ICM_BASE_MEMBERS                                                                 \
																						\
	/* Public: */																		\
																							\
	/* Get the native input space and output space ranges */								\
	/* This is an accurate number of what the underlying profile can hold. */				\
	void (*get_lutranges) (struct _icmLuBase *p,											\
		double *inmin, double *inmax,		/* Range of inspace values */					\
		double *outmin, double *outmax);	/* Range of outspace values */					\
																							\
	/* Get the effective input space and output space ranges */								\
	/* This may not be accurate when the effective type is different to the native type */	\
	void (*get_ranges) (struct _icmLuBase *p,												\
		double *inmin, double *inmax,		/* Range of inspace values */					\
		double *outmin, double *outmax);	/* Range of outspace values */					\
																							\
	/* Initialise the white and black points from the ICC tags */							\
	/* and the corresponding absolute<->relative conversion matrices */						\
	/* Return nz on error */																\
	int (*init_wh_bk)(struct _icmLuBase *p);												\
																							\
	/* Get the LU white and black points in XYZ space. */									\
	void (*wh_bk_points)(struct _icmLuBase *p, icmXYZNumber *wht, icmXYZNumber *blk);		\
																							\
	/* Translate color values through profile in effective in and out colorspaces, */		\
	/* return values: */																	\
	/* 0 = success, 1 = warning: clipping occured, 2 = fatal: other error */				\
	/* Note that clipping is not a reliable means of detecting out of gamut */				\
	/* in the lookup(bwd) call for clut based profiles. */									\
	int (*lookup) (struct _icmLuBase *p, double *out, double *in);							\
																							\
																							\
	/* Alternate to above, splits color conversion into three steps. */						\
	/* Colorspace of _in and _out and _core are the effective in and out */					\
	/* colorspaces. Note that setting colorspace overrides will probably. */				\
	/* force _in and/or _out to be dumy (unity) transforms. */								\
																							\
	/* Per channel input lookup (may be unity): */											\
	int (*lookup_in) (struct _icmLuBase *p, double *out, double *in);						\
																							\
	/* Intra channel conversion: */															\
	int (*lookup_core) (struct _icmLuBase *p, double *out, double *in);						\
																							\
	/* Per channel output lookup (may be unity): */											\
	int (*lookup_out) (struct _icmLuBase *p, double *out, double *in);						\
																							\
	/* Inverse versions of above - partially implemented */									\
	/* Inverse per channel input lookup (may be unity): */									\
	int (*lookup_inv_in) (struct _icmLuBase *p, double *out, double *in);					\
																							\


/* Base lookup object */
struct _icmLuBase {
	LU_ICM_NN_BASE_MEMBERS
}; typedef struct _icmLuBase icmLuBase;


/* Algorithm specific lookup classes */

/* Monochrome  Fwd & Bwd type object */
struct _icmLuMono {
	LU_ICM_NN_BASE_MEMBERS
	icmCurve    *grayCurve;

	/* Overall lookups */
	int (*fwd_lookup) (struct _icmLuBase *p, double *out, double *in);
	int (*bwd_lookup) (struct _icmLuBase *p, double *out, double *in);

	/* Components of lookup */
	int (*fwd_curve) (struct _icmLuMono *p, double *out, double *in);
	int (*fwd_map)   (struct _icmLuMono *p, double *out, double *in);
	int (*fwd_abs)   (struct _icmLuMono *p, double *out, double *in);
	int (*bwd_abs)   (struct _icmLuMono *p, double *out, double *in);
	int (*bwd_map)   (struct _icmLuMono *p, double *out, double *in);
	int (*bwd_curve) (struct _icmLuMono *p, double *out, double *in);

}; typedef struct _icmLuMono icmLuMono;

/* 3D Matrix Fwd & Bwd type object */
struct _icmLuMatrix {
	LU_ICM_NN_BASE_MEMBERS
	icmCurve    *redCurve, *greenCurve, *blueCurve;
	icmXYZArray *redColrnt, *greenColrnt, *blueColrnt;
    double		mx[3][3];	/* 3 * 3 conversion matrix */
    double		bmx[3][3];	/* 3 * 3 backwards conversion matrix */

	/* Overall lookups */
	int (*fwd_lookup) (struct _icmLuBase *p, double *out, double *in);
	int (*bwd_lookup) (struct _icmLuBase *p, double *out, double *in);

	/* Components of lookup */
	int (*fwd_curve)  (struct _icmLuMatrix *p, double *out, double *in);
	int (*fwd_matrix) (struct _icmLuMatrix *p, double *out, double *in);
	int (*fwd_abs)    (struct _icmLuMatrix *p, double *out, double *in);
	int (*bwd_abs)    (struct _icmLuMatrix *p, double *out, double *in);
	int (*bwd_matrix) (struct _icmLuMatrix *p, double *out, double *in);
	int (*bwd_curve)  (struct _icmLuMatrix *p, double *out, double *in);

}; typedef struct _icmLuMatrix icmLuMatrix;

/* Multi-D. Lut type object */
struct _icmLuLut {
	LU_ICM_NN_BASE_MEMBERS

	/* private: */
	icmLut *lut;								/* Lut to use */
	int    usematrix;							/* non-zero if matrix should be used */
    double imx[3][3];							/* 3 * 3 inverse conversion matrix */
	int    imx_valid;							/* Inverse matrix is valid */
	void (*in_normf)(double *out, double *in);	/* Lut input data normalizing function */
	void (*in_denormf)(double *out, double *in);/* Lut input data de-normalizing function */
	void (*out_normf)(double *out, double *in);	/* Lut output data normalizing function */
	void (*out_denormf)(double *out, double *in);/* Lut output de-normalizing function */
	void (*e_in_denormf)(double *out, double *in);/* Effective input de-normalizing function */
	void (*e_out_denormf)(double *out, double *in);/* Effecive output de-normalizing function */
	/* function chosen out of lut->lookup_clut_sx and lut->lookup_clut_nl to imp. clut() */
	int (*lookup_clut) (struct _icmLut *pp, double *out, double *in);	/* clut function */

	/* public: */

	/* Components of lookup */
	int (*in_abs)  (struct _icmLuLut *p, double *out, double *in);	/* Should be in icmLut ? */
	int (*matrix)  (struct _icmLuLut *p, double *out, double *in);
	int (*input)   (struct _icmLuLut *p, double *out, double *in);
	int (*clut)    (struct _icmLuLut *p, double *out, double *in);
	int (*output)  (struct _icmLuLut *p, double *out, double *in);
	int (*out_abs) (struct _icmLuLut *p, double *out, double *in);	/* Should be in icmLut ? */

	/* Some inverse components, in reverse order */
	/* Should be in icmLut ??? */
	int (*inv_out_abs) (struct _icmLuLut *p, double *out, double *in);
	int (*inv_output)  (struct _icmLuLut *p, double *out, double *in);
	/* inv_clut is beyond scope of icclib. See argyll for solution! */
	int (*inv_input)   (struct _icmLuLut *p, double *out, double *in);
	int (*inv_matrix)  (struct _icmLuLut *p, double *out, double *in);
	int (*inv_in_abs)  (struct _icmLuLut *p, double *out, double *in);

	/* Get various types of information about the LuLut */
	/* (Note that there's no reason why white/black point info */
	/*  isn't being made available for other Lu types) */
	void (*get_info) (struct _icmLuLut *p, icmLut **lutp,
	                 icmXYZNumber *pcswhtp, icmXYZNumber *whitep,
	                 icmXYZNumber *blackp);

	/* Get the matrix contents */
	void (*get_matrix) (struct _icmLuLut *p, double m[3][3]);

}; typedef struct _icmLuLut icmLuLut;

/* Named colors lookup object */
struct _icmLuNamed {
	LU_ICM_BASE_MEMBERS

  /* Private: */

  /* Public: */

	/* Get various types of information about the Named lookup */
	void (*get_info) (struct _icmLuLut *p, 
	                 icmXYZNumber *pcswhtp, icmXYZNumber *whitep,
	                 icmXYZNumber *blackp,
	                 int *maxnamesize,
	                 int *num_colors, char *prefix, char *suffix);


	/* Lookup a name that includes prefix and suffix */
	int (*fullname_lookup) (struct _icmLuNamed *p, double *pcs, double *dev, char *in);

	/* Lookup a name that doesn't include prefix and suffix */
	int (*name_lookup) (struct _icmLuNamed *p, double *pcs, double *dev, char *in);

	/* Fill in a list with the nout closest named colors to the pcs target */
	int (*pcs_lookup) (struct _icmLuNamed *p, char **out, int nout, double *in);

	/* Fill in a list with the nout closest named colors to the device target */
	int (*dev_lookup) (struct _icmLuNamed *p, char **out, int nout, double *in);

}; typedef struct _icmLuNamed icmLuNamed;

/* ---------------------------------------------------------- */
/* A tag */
typedef struct {
    icTagSignature      sig;			/* The tag signature */
	icTagTypeSignature  ttype;			/* The tag type signature */
    unsigned int        offset;			/* File offset to start header */
    unsigned int        size;			/* Size in bytes (not including padding) */
    unsigned int        pad;			/* Padding in bytes */
	icmBase            *objp;			/* In memory data structure */
} icmTag;

/* Pseudo enumerations valid as parameter to get_luobj(): */

/* To be specified where an intent is not appropriate */
#define icmDefaultIntent ((icRenderingIntent)98)

/* Pseudo PCS colospace used to indicate the native PCS */
#define icmSigDefaultData ((icColorSpaceSignature) 0x0)

/* Pseudo colospace used to indicate a named colorspace */
#define icmSigNamedData ((icColorSpaceSignature) 0x1)


/* Arguments to set a non-default creation version */
typedef enum {
    icmVersionDefault       = 0,	/* Version 2.2.0 is default */
    icmVersion2_3           = 1,	/* Version 2.3.0 - Chromaticity Tag */
    icmVersion2_4           = 2,	/* Version 2.4.0 - Display etc. have intents */
    icmVersion4_1           = 3,	/* Version 4.1.0 - General V4 features */
} icmICCVersion;

/* The ICC object */
struct _icc {
  /* Public: */
	int          (*set_version)(struct _icc *p, icmICCVersion ver);
	                                                       /* For creation, use ICC V4 etc. */
	unsigned int (*get_size)(struct _icc *p);				/* Return total size needed, 0 = err. */
	int          (*read)(struct _icc *p, icmFile *fp, unsigned long of);	/* Returns error code */
	int          (*write)(struct _icc *p, icmFile *fp, unsigned long of);/* Returns error code */
	void         (*dump)(struct _icc *p, icmFile *op, int verb);	/* Dump whole icc */
	void         (*del)(struct _icc *p);						/* Free whole icc */
	int          (*find_tag)(struct _icc *p, icTagSignature sig);
							/* Returns 0 if found, 1 if found but not readable, 2 of not found */
	icmBase *    (*read_tag)(struct _icc *p, icTagSignature sig);
															/* Returns pointer to object */
	icmBase *    (*add_tag)(struct _icc *p, icTagSignature sig, icTagTypeSignature ttype);
															/* Returns pointer to object */
	int          (*rename_tag)(struct _icc *p, icTagSignature sig, icTagSignature sigNew);
															/* Rename and existing tag */
	icmBase *    (*link_tag)(struct _icc *p, icTagSignature sig, icTagSignature ex_sig);
															/* Returns pointer to object */
	int          (*unread_tag)(struct _icc *p, icTagSignature sig);
														/* Unread a tag (free on refcount == 0 */
	int          (*read_all_tags)(struct _icc *p); /* Read all the tags, non-zero on error. */

	int          (*delete_tag)(struct _icc *p, icTagSignature sig);
															/* Returns 0 if deleted OK */
	int          (*check_id)(struct _icc *p, ORD8 *id); /* Returns 0 if ID is OK, 1 if not present etc. */
	double       (*get_tac)(struct _icc *p, double *chmax); /* Returns total ink limit and channel maximums */
	icmLuBase *  (*get_luobj) (struct _icc *p,
                                     icmLookupFunc func,			/* Functionality */
	                                 icRenderingIntent intent,		/* Intent */
	                                 icColorSpaceSignature pcsor,	/* PCS override (0 = def) */
	                                 icmLookupOrder order);			/* Search Order */
	                                 /* Return appropriate lookup object */
	                                 /* NULL on error, check errc+err for reason */

	
    icmHeader       *header;			/* The header */
	char             err[512];			/* Error message */
	int              errc;				/* Error code */

  /* Private: ? */
	icmAlloc        *al;				/* Heap allocator */
	int              del_al;			/* NZ if heap allocator should be deleted */
	icmFile         *fp;				/* File associated with object */
	unsigned long    of;				/* Offset of the profile within the file */
    unsigned int     count;				/* Num tags in the profile */
    icmTag          *data;    			/* The tagTable and tagData */
	icmICCVersion    ver;				/* Version class, see icmICCVersion enum */

	}; typedef struct _icc icc;

/* ========================================================== */
/* Utility structures and declarations */

/* Type of primitives that can be read or written */
typedef enum {
	icmUInt8Number,
	icmUInt16Number,
	icmUInt32Number,
	icmUInt64Number,
	icmU8Fixed8Number,
	icmU16Fixed16Number,
	icmSInt8Number,
	icmSInt16Number,
	icmSInt32Number,
	icmSInt64Number,
	icmS15Fixed16Number,
	icmDCS8Number,
	icmDCS16Number,
	icmPCSNumber,			/* 16 bit PCS of profile */
	icmPCSXYZNumber,		/* 16 bit XYZ */
	icmPCSLab8Number,		/* 8 bit Lab */
	icmPCSLabNumber,		/* 16 bit Lab of profile */
	icmPCSLabV2Number,		/* 16 bit Lab Version 2 */
	icmPCSLabV4Number		/* 16 bit Lab Version 4 */
} icmPrimType;

/* Structure to hold pseudo-hilbert counter info */
struct _psh {
	int      di;	/* Dimensionality */
	unsigned int res;	/* Resolution per coordinate */
	unsigned int bits;	/* Bits per coordinate */
	unsigned int ix;	/* Current binary index */
	unsigned int tmask;	/* Total 2^n count mask */
	unsigned int count;	/* Usable count */
}; typedef struct _psh psh;

/* Type of encoding to be returned as a string */
typedef enum {
    icmScreenEncodings,
    icmDeviceAttributes,
	icmProfileHeaderFlags,
	icmAsciiOrBinaryData,
	icmTagSignature,
	icmTechnologySignature,
	icmTypeSignature,
	icmColorSpaceSignature,
	icmProfileClassSignature,
	icmPlatformSignature,
	icmMeasurementFlare,
	icmMeasurementGeometry,
	icmRenderingIntent,
	icmTransformLookupFunc,
	icmSpotShape,
	icmStandardObserver,
	icmIlluminant,
	icmLuAlg
} icmEnumType;


/* A helper object that computes MD5 checksums */
struct _icmMD5 {
  /* Private: */
//	icc *icp;				/* ICC we're part of */
	icmAlloc *al;			/* Heap allocator */
	int fin;				/* Flag, nz if final has been called */
	ORD32 sum[4];			/* Current/final checksum */
	unsigned int tlen;		/* Total length added in bytes */
	ORD8 buf[64];			/* Partial buffer */

  /* Public: */
	void (*reset)(struct _icmMD5 *p);	/* Reset the checksum */
	void (*add)(struct _icmMD5 *p, ORD8 *buf, unsigned int len);	/* Add some bytes */
	void (*get)(struct _icmMD5 *p, ORD8 chsum[16]);		/* Finalise and get the checksum */
	void (*del)(struct _icmMD5 *p);		/* We're done with the object */

}; typedef struct _icmMD5 icmMD5;

/* Create a new MD5 checksumming object, with a reset checksum value */
/* Return it or NULL if there is an error */
icmMD5 *new_icmMD5(icmAlloc *al);


/* Implementation of file access class to compute an MD5 checksum */
struct _icmFileMD5 {
	ICM_FILE_BASE
	int (*get_errc)(struct _icmFile *p);		/* Check if there was an error */

	/* Private: */
	icmAlloc *al;		/* Heap allocator */
	icmMD5 *md5;		/* MD5 object */
	unsigned long of;	/* Current offset */
	int errc;			/* Error code, 0 for OK */

}; typedef struct _icmFileMD5 icmFileMD5;

/* Create a dumy file access class with allocator */
icmFile *new_icmFileMD5_a(icmMD5 *md5, icmAlloc *al);


/* ========================================================== */
/* Public function declarations */
/* Create an empty object. Return null on error */
extern ICCLIB_API icc *new_icc_a(icmAlloc *al);		/* With allocator class */

/* If SEPARATE_STD not defined: */
extern ICCLIB_API icc *new_icc(void);				/* Default allocator */

/* - - - - - - - - - - - - - */
/* Some useful utilities: */

/* Read a given primitive type. Return non-zero on error */
extern ICCLIB_API int read_Primitive(icc *icpp, icmPrimType ptype, void *prim, char *p);

/* Write a given primitive type. Return non-zero on error */
extern ICCLIB_API int write_Primitive(icc *icp, icmPrimType ptype, char *p, void *prim);

/* Return a string that represents a tag */
extern ICCLIB_API char *tag2str(int tag);

/* Return a tag created from a string */
extern ICCLIB_API int str2tag(const char *str);

/* Utility to define a non-standard tag, arguments are 4 character constants */
#define icmMakeTag(A,B,C,D)	\
	(   (((ORD32)(ORD8)(A)) << 24)	\
	  | (((ORD32)(ORD8)(B)) << 16)	\
	  | (((ORD32)(ORD8)(C)) << 8)	\
	  |  ((ORD32)(ORD8)(D)) )

/* Return a string description of the given enumeration value */
extern ICCLIB_API const char *icm2str(icmEnumType etype, int enumval);

/* Return the number of channels for the given color space. Return 0 if unknown. */
extern ICCLIB_API unsigned int icmCSSig2nchan(icColorSpaceSignature sig);

/* Return the individual channel names and number of channels give a colorspace signature. */
/* Return 0 if it is not a colorspace that itself defines particular channels, */ 
/* 1 if it is a colorant based colorspace, and 2 if it is not a colorant based space */
extern ICCLIB_API unsigned int icmCSSig2chanNames( icColorSpaceSignature sig, char *cvals[]);


/* Simple macro to transfer an array to an XYZ number */
#define icmAry2XYZ(xyz, ary) ((xyz).X = (ary)[0], (xyz).Y = (ary)[1], (xyz).Z = (ary)[2])

/* And the reverse */
#define icmXYZ2Ary(ary, xyz) ((ary)[0] = (xyz).X, (ary)[1] = (xyz).Y, (ary)[2] = (xyz).Z)

/* Simple macro to transfer an XYZ number to an XYZ number */
#define icmXYZ2XYZ(d_xyz, s_xyz) ((d_xyz).X = (s_xyz).X, (d_xyz).Y = (s_xyz).Y, \
                                  (d_xyz).Z = (s_xyz).Z)

/* Simple macro to transfer an 3array to 3array */
#define icmAry2Ary(d_ary, s_ary) ((d_ary)[0] = (s_ary)[0], (d_ary)[1] = (s_ary)[1], \
                                  (d_ary)[2] = (s_ary)[2])

/* CIE XYZ to perceptual Lab */
extern ICCLIB_API void icmXYZ2Lab(icmXYZNumber *w, double *out, double *in);

/* Perceptual Lab to CIE XYZ */
extern ICCLIB_API void icmLab2XYZ(icmXYZNumber *w, double *out, double *in);

/* LCh to Lab */
extern ICCLIB_API void icmLCh2Lab(double *out, double *in);

/* Lab to LCh */
extern ICCLIB_API void icmLab2LCh(double *out, double *in);

/* XYZ to Yxy */
extern ICCLIB_API void icmXYZ2Yxy(double *out, double *in);

/* Yxy to XYZ */
extern ICCLIB_API void icmYxy2XYZ(double *out, double *in);

/* CIE XYZ to perceptual Luv */
extern ICCLIB_API void icmXYZ2Luv(icmXYZNumber *w, double *out, double *in);

/* Perceptual Luv to CIE XYZ */
extern ICCLIB_API void icmLuv(icmXYZNumber *w, double *out, double *in);


/* CIE XYZ to perceptual CIE 1960 UCS */
extern ICCLIB_API void icmXYZ2UCS(icmXYZNumber *w, double *out, double *in);

/* Perceptual CIE 1960 UCS to CIE XYZ */
extern ICCLIB_API void icmUCS2XYZ(icmXYZNumber *w, double *out, double *in);


/* The standard D50 illuminant value */
extern ICCLIB_API icmXYZNumber icmD50;
extern ICCLIB_API icmXYZNumber icmD50_100;		/* Scaled to 100 */

/* The standard D65 illuminant value */
extern ICCLIB_API icmXYZNumber icmD65;
extern ICCLIB_API icmXYZNumber icmD65_100;		/* Scaled to 100 */

/* The default black value */
extern ICCLIB_API icmXYZNumber icmBlack;


/* Initialise a pseudo-hilbert grid counter, return total usable count. */
extern ICCLIB_API unsigned psh_init(psh *p, int di, unsigned int res, int co[]);

/* Reset the counter */
extern ICCLIB_API void psh_reset(psh *p);

/* Increment pseudo-hilbert coordinates */
/* Return non-zero if count rolls over to 0 */
extern ICCLIB_API int psh_inc(psh *p, int co[]);


/* RGB primaries to device to RGB->XYZ transform matrix */
/* Return non-zero if matrix would be singular */
int icmRGBprim2matrix(
	icmXYZNumber white,		/* White point */
	icmXYZNumber red,		/* Red colorant */
	icmXYZNumber green,		/* Green colorant */
	icmXYZNumber blue,		/* Blue colorant */
	double mat[3][3]		/* Destination matrix */
);

/* Chromatic Adaption transform utility */
/* Return a 3x3 chromatic adaption matrix */

#define ICM_CAM_BRADFORD	0x0001	/* Use Bradford sharpened response space */
#define ICM_CAM_MULMATRIX	0x0002	/* Transform the given matrix */

void icmChromAdaptMatrix(
	int flags,				/* Flags as defined below */
	icmXYZNumber d_wp,		/* Destination white point */
	icmXYZNumber s_wp,		/* Source white point */
	double mat[3][3]		/* Destination matrix */
);

/* - - - - - - - - - - - - - - */

/* Note icmAry2Ary() above */

/* Copy a 3x3 transform matrix */
void icmCpy3x3(double out[3][3], double mat[3][3]);

/* Scale each element of a 3x3 transform matrix */
void icmScale3x3(double dst[3][3], double src[3][3], double scale);

/* Multiply 3 vector by 3x3 transform matrix */
/* Organization is mat[out][in] */
void icmMulBy3x3(double out[3], double mat[3][3], double in[3]);

/* Add one 3x3 to another */
/* dst = src1 + src2 */
void icmAdd3x3(double dst[3][3], double src1[3][3], double src2[3][3]);

/* Tensor product. Multiply two 3 vectors to form a 3x3 matrix */
/* src1[] forms the colums, and src2[] forms the rows in the result */
void icmTensMul3(double dst[3][3], double src1[3], double src2[3]);

/* Multiply one 3x3 with another */
/* dst = src * dst */
void icmMul3x3(double dst[3][3], double src[3][3]);

/* Multiply one 3x3 with another #2 */
/* dst = src1 * src2 */
void icmMul3x3_2(double dst[3][3], double src1[3][3], double src2[3][3]);

/* Add two 3 vectors */
void icmAdd3(double out[3], double in1[3], double in2[3]);

#define ICMADD3(o, i, j) ((o)[0] = (i)[0] + (j)[0], (o)[1] = (i)[1] + (j)[1], (o)[2] = (i)[2] + (j)[2])

/* Subtract two 3 vectors, out = in1 - in2 */
void icmSub3(double out[3], double in1[3], double in2[3]);

#define ICMSUB3(o, i, j) ((o)[0] = (i)[0] - (j)[0], (o)[1] = (i)[1] - (j)[1], (o)[2] = (i)[2] - (j)[2])

/* Compute the dot product of two 3 vectors */
double icmDot3(double in1[3], double in2[3]);

#define ICMDOT3(o, i, j) ((o) = (i)[0] * (j)[0] + (i)[1] * (j)[1] + (i)[2] * (j)[2])

/* Compute the cross product of two 3D vectors, out = in1 x in2 */
void icmCross3(double out[3], double in1[3], double in2[3]);

/* Compute the norm squared (length squared) of a 3 vector */
double icmNorm3sq(double in[3]);

#define ICMNORM3SQ(i) ((i)[0] * (i)[0] + (i)[1] * (i)[1] + (i)[2] * (i)[2])

/* Compute the norm (length) of a 3 vector */
double icmNorm3(double in[3]);

#define ICMNORM3(i) sqrt((i)[0] * (i)[0] + (i)[1] * (i)[1] + (i)[2] * (i)[2])

/* Scale a 3 vector by the given ratio */
void icmScale3(double out[3], double in[3], double rat);

#define ICMSCALE3(o, i, j) ((o)[0] = (i)[0] * (j), (o)[1] = (i)[1] * (j), (o)[2] = (i)[2] * (j))

/* Normalise a 3 vector to the given length. Return nz if not normalisable */
int icmNormalize3(double out[3], double in[3], double len);

/* Compute the norm squared (length squared) of two point vector */
double icmNorm33sq(double in1[3], double in0[3]);

/* Compute the norm (length) of two point vector */
double icmNorm33(double in1[3], double in0[3]);

/* Scale a two point vector by the given ratio */
void icmScale33(double out[3], double in1[3], double in0[3], double rat);

/* Normalise a two point vector to the given length. */
/* The new location of in1[] is returned in out[], in0[] is the origin. */
/* Return nz if not normalisable */
int icmNormalize33(double out[3], double in1[3], double in0[3], double len);

/* Compute the determinant of a 3x3 matrix */
double icmDet3x3(double in[3][3]);

/* Invert a 3x3 transform matrix. Return 1 if error. */
int icmInverse3x3(double out[3][3], double in[3][3]);

/* Given two 3D points, create a matrix that rotates */
/* and scales one onto the other about the origin 0,0,0. */
/* Use icmMulBy3x3() to apply this to other points */
void icmRotMat(double m[3][3], double s[3], double t[3]);

/* Multiply 3 array by 3x4 transform matrix */
void icmMul3By3x4(double out[3], double mat[3][4], double in[3]);

/* Given two 3D vectors, create a matrix that translates, */
/* rotates and scales one onto the other. */
/* Use icmMul3By3x4 to apply this to other points */
void icmVecRotMat(double m[3][4], double s1[3], double s0[3], double t1[3], double t0[3]);

/* - - - - - - - - - - - - - - - - - - - - - - - */

/* Return the normal Delta E given two Lab values */
extern ICCLIB_API double icmLabDE(double *in0, double *in1);

/* Return the normal Delta E squared, given two Lab values */
extern ICCLIB_API double icmLabDEsq(double *in0, double *in1);

/* Return the normal Delta E given two XYZ values */
extern ICCLIB_API double icmXYZLabDE(icmXYZNumber *w, double *in0, double *in1);


/* Return the CIE94 Delta E color difference measure for two Lab values */
extern ICCLIB_API double icmCIE94(double *in0, double *in1);

/* Return the CIE94 Delta E color difference measure squared, for two Lab values */
extern ICCLIB_API double icmCIE94sq(double *in0, double *in1);

/* Return the CIE94 Delta E color difference measure for two XYZ values */
extern ICCLIB_API double icmXYZCIE94(icmXYZNumber *w, double *in0, double *in1);


/* Return the CIEDE2000 Delta E color difference measure for two Lab values */
extern ICCLIB_API double icmCIE2K(double *in0, double *in1);

/* Return the CIEDE2000 Delta E color difference measure squared, for two Lab values */
extern ICCLIB_API double icmCIE2Ksq(double *in0, double *in1);

/* Return the CIEDE2000 Delta E color difference measure for two XYZ values */
extern ICCLIB_API double icmXYZCIE2K(icmXYZNumber *w, double *in0, double *in1);

/* - - - - - - - - - - - - - - - - - - - - - - - */
/* Clip Lab, while maintaining hue angle. */
/* Return nz if clipping occured */
int icmClipLab(double out[3], double in[3]);

/* Clip XYZ, while maintaining hue angle */
/* Return nz if clipping occured */
int icmClipXYZ(double out[3], double in[3]);

/* ---------------------------------------------------------- */

#ifdef __cplusplus
	}
#endif

#endif /* ICC_H */

