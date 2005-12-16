#!/bin/sh
echo "Making Argyll distribution zip archive argyll.zip... "

rm -f argyll.zip

#do all
for i in . numlib plot tiff icc rspl imdi cgats target spectro scanin gamut xicc profile link tweak lib h ref doc
do
(echo
echo "#### Doing $i ####"
dos2unix `cat $i/afiles`
zip -9 argyll.zip `cat $i/afiles`
)
done

