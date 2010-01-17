#!/bin/sh
echo "Script to invoke Jam from the top"

if [ X$OS != "XWindows_NT" ] ; then
	# Fixup issues with the .zip format
	chmod +x *.sh
	chmod +x tiff/configure
	chmod +x libusb/configure
fi

# Make sure that some environment variable are visible to Jam:
export OSTYPE MACHTYPE HOSTTYPE

jam -q -fJambase -j${NUMBER_OF_PROCESSORS:-2}
