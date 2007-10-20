#!/bin/sh
echo "Simple batch file to invoke Jam in all the subdirectories"
#
# Note that running xicc twice is a hack to workaround a circular dependency
# with regard to spectro/libinst.lib
#
# If you're compiling under VC++ on a multiprocessor box, with debugging turned on,
# you may have to comment out the multithread jam, and use the single thread one
# just below.
#
# Many people don't seem to have . in their path, so fix this up:
PATH=$PATH:.
if [ X$OS != "XWindows_NT" ] ; then
	# Fixup issues with the .zip format
	chmod +x *.ksh
	chmod +x tiff/configure
	chmod +x libusb/configure
fi
for i in `cat blddirs`
do
(echo ------------; echo $i; echo ----------; cd $i; jam -f../Jambase -j${NUMBER_OF_PROCESSORS:-1}; cd ..; echo)
#(echo ------------; echo $i; echo ----------; cd $i; jam -f../Jambase ; cd ..; echo)
done
