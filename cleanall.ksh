#!/bin/sh
echo "Simple batch file to clean all the subdirectories"
for i in numlib tiff plot icc rspl imdi cgats gamut xicc spectro target scanin profile link tweak
do
(echo ------------; echo $i; echo ----------; cd $i; rm *.obj *.lib *.exe *.o *.a ; cd ..; echo)
done
