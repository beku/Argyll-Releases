# Boilerplate Makefile for compiling tifflib under Windows

# Microsoft C++, WinNT setup
#include Makefile.WNT

# IBM C++, WinNT setup
include Makefile.IBMNT

# Generic UNIX setup
#include Makefile.UNIX

###############################

#Set optimisation on
#CCFLAGS = $(CCFLAGSDEF) $(CCOPTFLAG)
#Set debugging on
CCFLAGS = $(CCFLAGSDEF) $(CCDEBUGFLAG)
STDHDRS = $(STDHDRSDEF)
LINKFLAGS = $(LINKFLAGSDEF) $(LINKDEBUGFLAG)

all:: libtiff$(SUFLIB)

# support library for generation utilities

PORTOBJ = port/getopt$(SUFOBJ) port/strcasecmp$(SUFOBJ) port/strtoul$(SUFOBJ)

libport$(SUFLIB): $(PORTOBJ)
	$(LIBU) $(LIBOF)$@ $?

# generation utilities
mkg3states$(SUFEXE) : libtiff/mkg3states$(SUFOBJ) libport$(SUFLIB)
	$(LINK) $(LINKOF)$@ libtiff/mkg3states$(SUFOBJ) libport$(SUFLIB)

libtiff/tif_fax3sm.c : mkg3states$(SUFEXE)
	mkg3states -c const libtiff\tif_fax3sm.c

# version file
mkversion$(SUFEXE) : libtiff/mkversion$(SUFOBJ)
	$(LINK) $(LINKOF)$@ libtiff/mkversion$(SUFOBJ)

libtiff/version.h : mkversion$(SUFEXE)
	mkversion -v VERSION -a dist/tiff.alpha -r RELEASE-DATE libtiff\version.h


# tiff library itlsef
#
#IPATH on libtiff =	../jpeg ;

LIBOBJS = libtiff/tif_fax3$(SUFOBJ) libtiff/tif_aux$(SUFOBJ) libtiff/tif_close$(SUFOBJ) \
          libtiff/tif_codec$(SUFOBJ) libtiff/tif_fax3sm$(SUFOBJ) libtiff/tif_predict$(SUFOBJ) \
          libtiff/tif_compress$(SUFOBJ) libtiff/tif_dir$(SUFOBJ) libtiff/tif_dirinfo$(SUFOBJ) \
          libtiff/tif_dirread$(SUFOBJ) libtiff/tif_dirwrite$(SUFOBJ) \
          libtiff/tif_dumpmode$(SUFOBJ) libtiff/tif_error$(SUFOBJ) \
          libtiff/tif_getimage$(SUFOBJ) libtiff/tif_flush$(SUFOBJ) libtiff/tif_luv$(SUFOBJ) \
          libtiff/tif_lzw$(SUFOBJ) libtiff/tif_next$(SUFOBJ) libtiff/tif_open$(SUFOBJ) \
          libtiff/tif_packbits$(SUFOBJ) libtiff/tif_print$(SUFOBJ) libtiff/tif_read$(SUFOBJ) \
          libtiff/tif_swab$(SUFOBJ) libtiff/tif_strip$(SUFOBJ) libtiff/tif_thunder$(SUFOBJ) \
          libtiff/tif_tile$(SUFOBJ) libtiff/tif_version$(SUFOBJ) libtiff/tif_warning$(SUFOBJ) \
          libtiff/tif_write$(SUFOBJ) libtiff/tif_msdos$(SUFOBJ)

libtiff$(SUFLIB): $(LIBOBJS) libtiff/version.h
	$(LIBU) $(LIBOF)$@  $(LIBOBJS)

#ObjectHdrs libtiff.lib : $(IPATH) ;
#Library libtiff.lib : $(LIBSRCS) ;
#
## copy header file to lib directory
#File tiff.h : libtiff/tiff.h ;
#File tiffio.h : libtiff/tiffio.h ;
#
