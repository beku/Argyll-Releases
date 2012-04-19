#!/bin/sh
echo "Script to invoke Jam and then package the binary release."

#	Typical environment variables:
#
#	Platform						$OSTYPE		MACHTYPE				HOSTTYPE
#
#	Win2K [CMD.EXE]					(none)		(none)					(none)		
#
#	Cygwin Win2K [bash]				cygwin		i686-pc-cygwin			i686
#
#	OS X PPC 10.3 [zsh]				darwin7.0	powerpc					(none)
#
#	OS X i386 10.4 [bash]			darwin8.0	i386-apple-darwin8.0	i386
#
#	OS X i386 10.5 [bash]			darwin9.0	i386-apple-darwin9.0	i386
#
#	OS X i386 10.6 [bash]			darwin10.0	i386-apple-darwin8.0	i386
#
#	Linux RH 4.0 [bash]				linux-gnu	i686-redhat-linux-gnu	i686
#
#	Linux Fedora 7.1 [bash]			linux-gnu	i386-redhat-linux-gnu	i386
#	Linux Ubuntu  ??7				linux-gnu	i486-pc-linux-gnu		i686
#
#	Linux Fedora 7.1 64 bit [bash]	linux-gnu	x86_64-redhat-linux-gnu	x86_64
#

# Set the environment string VERSION from the #define, ie 1.0.0
VERSION=`grep ARGYLL_VERSION_STR h/aconfig.h | sed 's/#define ARGYLL_VERSION_STR //' | sed 's/"//g'`

echo "About to make Argyll binary distribution $PACKAGE"

TOPDIR=Argyll_V$VERSION

if [ X$OS != "XWindows_NT" ] ; then
	# Fixup issues with the .zip format
	chmod +x *.sh
	chmod +x tiff/configure
#	chmod +x libusb/configure
fi

# Make sure that some environment variable are visible to Jam:
export OSTYPE MACHTYPE HOSTTYPE

# .sp come from profile, .cht from scanin and .ti3 from spectro
rm -f bin/*.exe bin/*.dll
rm -f ref/*.sp ref/*.cht ref/*.ti2

# Make sure it's built and installed
if ! jam -q -fJambase -j${NUMBER_OF_PROCESSORS:-2} -sBUILTIN_TIFF=true -sBUILTIN_JPEG=true install ; then
	echo "Build failed!"
	exit 1
fi 

# Maybe we could get Jam to do the following ?

if [ X$OS = "XWindows_NT" ] ; then
	echo "We're on MSWindows!"
	PACKAGE=Argyll_V${VERSION}_win32_exe.zip
	USBDIRS="libusb1"
	USBBINFILES="binfiles.msw"
	unset USETAR
else if [ X$OSTYPE = "Xdarwin7.0" ] ; then
	echo "We're on OSX 10.3 PPC!"
	PACKAGE=Argyll_V${VERSION}_osx10.3_ppc_bin.tgz
	USBDIRS="libusb1"
	USBBINFILES="binfiles.osx"
	USETAR=true
else if [ X$OSTYPE = "Xdarwin8.0" ] ; then
	if [ X$MACHTYPE = "Xi386-apple-darwin8.0" ] ; then
		echo "We're on OSX 10.4 i386!"
		PACKAGE=Argyll_V${VERSION}_osx10.4_i86_bin.tgz
	else if [ X$MACHTYPE = "Xpowerpc-apple-darwin8.0" ] ; then
		echo "We're on OSX 10.4 PPC!"
		PACKAGE=Argyll_V${VERSION}_osx10.4_ppc_bin.tgz
	fi
	fi
	USBDIRS="libusb1"
	USBBINFILES="binfiles.osx"
	USETAR=true
else if [ X$OSTYPE = "Xdarwin9.0" ] ; then
	if [ X$MACHTYPE = "Xi386-apple-darwin9.0" ] ; then
		echo "We're on OSX 10.5 i386!"
		PACKAGE=Argyll_V${VERSION}_osx10.5_i86_bin.tgz
	fi
	USBDIRS="libusb1"
	USBBINFILES="binfiles.osx"
	USETAR=true
else if [ X$OSTYPE = "Xlinux-gnu" ] ; then
	if [ X$MACHTYPE = "Xi686-redhat-linux-gnu" \
      -o X$MACHTYPE = "Xi386-redhat-linux-gnu" \
      -o X$MACHTYPE = "Xi486-pc-linux-gnu" ] ; then
		echo "We're on Linux x86!"
		PACKAGE=Argyll_V${VERSION}_linux_x86_bin.tgz
	else if [ X$MACHTYPE = "Xx86_64-redhat-linux-gnu" ] ; then
		echo "We're on Linux x86_64!"
		PACKAGE=Argyll_V${VERSION}_linux_x86_64_bin.tgz
	fi
	fi
	USBDIRS="libusb1"
	USBBINFILES="binfiles.unix"
	USETAR=true
fi
fi
fi
fi
fi

echo "Making GNU Argyll binary distribution $PACKAGE for Version $VERSION"

rm -rf $TOPDIR
mkdir $TOPDIR

# Collect the names of all the files that we're going to package
unset topfiles; for i in `cat binfiles`; do topfiles="$topfiles ${i}"; done
unset docfiles; for i in `cat doc/afiles`; do docfiles="$docfiles doc/${i}"; done
unset usbfiles;
for j in ${USBDIRS}; do
	if [ ${j} ]; then
		for i in `cat ${j}/${USBBINFILES}`; do usbfiles="$usbfiles ${j}/${i}"; done
	fi
done

allfiles="${topfiles} bin/* ref/* ${docfiles} ${usbfiles}"

# Copy all the files to the package top directory
for i in ${allfiles}; do
	path=${i%/*}		# extract path without filename
	file=${i##*/}		# extract filename
	if [ $path = $i ] ; then
		path=
	fi
	if [ X$path != "X" ] ; then
		mkdir -p $TOPDIR/${path}
	fi
	cp $i $TOPDIR/$i
done

# Create the package
rm -f $PACKAGE
if [ X$USETAR = "Xtrue" ] ; then
	tar -czvf $PACKAGE $TOPDIR
	# tar -xzf to extract
	# tar -tzf to list
else
	zip -9 -r $PACKAGE $TOPDIR
	# unzip to extract
	# unzip -l to list
fi
rm -rf $TOPDIR
echo "Done GNU Argyll binary distribution $PACKAGE"

exit 0

