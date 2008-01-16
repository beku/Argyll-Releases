#!/bin/sh
echo "Simple batch file to invoke Jam in all the subdirectories and then"
echo "to package the binary release."

#	Typical environment variables
#
#	Platform				$OSTYPE		MACHTYPE				HOSTTYPE
#
#	Win2K CMD.EXE			(none)		(none)					(none)		
#
#	Cygwin Win2K bash		cygwin		i686-pc-cygwin			i686
#
#	OS X PPC 10.3 zsh		darwin7.0	powerpc					(none)
#
#	OS X i386 10.4 bash		darwin8.0	i386-apple-darwin8.0	i386
#
#	Linux RH 4.0 bash		linux-gnu	i686-redhat-linux-gnu	i686
#

VERSION=0.70Beta8

PATH=$PATH:.

if [ X$OS != "XWindows_NT" ] ; then
	# Fixup issues with the .zip format
	chmod +x *.ksh
	chmod +x tiff/configure
	chmod +x libusb/configure
fi

rm -f bin/*.exe bin/*.dll
rm -f ref/*.sp ref/*.cht ref/*.ti2

for i in `cat blddirs`
do
(echo ------------; echo $i; echo ----------; cd $i; jam -f../Jambase ; jam -f../Jambase install ; cd ..; echo)
done

if [ X$OS = "XWindows_NT" ] ; then
	echo "We're on MSWindows!"
	PACKAGE=argyllV${VERSION}_win32_exe.zip
	USBDIR=libusbw
else if [ X$OSTYPE = "Xdarwin7.0" ] ; then
	echo "We're on OSX 10.3 PPC!"
	PACKAGE=argyllV${VERSION}_osx10.3_ppc_bin.zip
	unset USBDIR
else if [ X$OSTYPE = "Xdarwin8.0" ] ; then
	if [ X$MACHTYPE = "Xi386-apple-darwin8.0" ] ; then
		echo "We're on OSX 10.4 i386!"
		PACKAGE=argyllV${VERSION}_osx10.4_i86_bin.zip
	else if [ X$MACHTYPE = "Xpowerpc-apple-darwin8.0" ] ; then
		echo "We're on OSX 10.4 PPC!"
		PACKAGE=argyllV${VERSION}_osx10.4_ppc_bin.zip
	fi
	fi
	unset USBDIR
else if [ X$OSTYPE = "Xlinux-gnu" ] ; then
	echo "We're on Linux!"
	PACKAGE=argyllV${VERSION}_linux_x86_bin.zip
	USBDIR=libusb
fi
fi
fi
fi

echo "Package = " $PACKAGE


# Create zip archive of documents and exectutables
unset docfiles; for i in `cat doc/afiles`; do docfiles="$docfiles doc/${i}"; done
unset usbfiles;
if [ $USBDIR ] ; then
	for i in `cat $USBDIR/binfiles`; do usbfiles="$usbfiles $USBDIR/${i}"; done
fi
rm -f $PACKAGE
zip -9 -r $PACKAGE bin ref ${docfiles} ${usbfiles}


