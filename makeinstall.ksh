#!/bin/sh
echo "Simple batch file to invoke Jam in all the subdirectories and install"
echo "binaries in the ./bin directory."
echo "(Hope you've unzip'd tiff/tiff.zip and if unix/OSX run tiff/configure/make ?)"
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
for i in numlib tiff plot icc rspl imdi cgats gamut xicc spectro xicc target scanin profile link tweak
do
(echo ------------; echo $i; echo ----------; cd $i; jam -j${NUMBER_OF_PROCESSORS:-1} install; cd ..; echo)
#(echo ------------; echo $i; echo ----------; cd $i; jam; cd ..; echo)
done
