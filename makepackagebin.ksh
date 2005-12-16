#!/bin/sh
echo "Simple batch file to invoke Jam in all the subdirectories and then"
echo "to package the binary release."

VERSION=0.53

PATH=$PATH:.

rm -f bin/*.exe
rm -f ref/*.sp ref/*.cht ref/*.ti2

for i in numlib tiff plot icc rspl imdi cgats gamut xicc spectro target scanin profile link tweak
do
(echo ------------; echo $i; echo ----------; cd $i; jam ; jam install ; cd ..; echo)
done

if [ X$OS = "XWindows_NT" ] ; then
	echo "We're on MSWindows!"
	PACKAGE=argyllV${VERSION}_win32_exe.zip
else if [ X$OSTYPE = "Xdarwin7.0" ] ; then
	echo "We're on OSX!"
	PACKAGE=argyllV${VERSION}_osx10.3_bin.zip
else if [ X$OSTYPE = "Xlinux-gnu" ] ; then
	echo "We're on Linux!"
	PACKAGE=argyllV${VERSION}_linux_x86_bin.zip
fi
fi
fi

echo "Package = " $PACKAGE


# Create zip archive of documents and exectutables
zip -9 -r $PACKAGE bin ref `cat doc/afiles`


