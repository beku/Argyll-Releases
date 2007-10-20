#!/bin/sh
## Do a simple sanity check of the major toolchain

if [ X$OS = "XWindows_NT" ] ; then
    echo "We're on MSWindows!"
    ICM=icm
else if [ X$OSTYPE = "Xdarwin7.0" ] ; then
    echo "We're on OSX!"
    ICM=icc
else if [ X$OSTYPE = "Xdarwin8.0" ] ; then
    echo "We're on OSX!"
    ICM=icc
else if [ X$OSTYPE = "Xlinux-gnu" ] ; then
    echo "We're on Linux!"
    ICM=icc
fi
fi
fi
fi

echo "############### Checking dispcal ##################"
spectro/dispcal -v -m -dfake test

echo "############### Checking targen ##################"
target/targen -v -d4 -f500 -l280 test
target/targen -v -d3 -r -f500 rgbtest

echo "############### Checking printtarg ##################"
target/printtarg -v -iSS -S test

#echo "############### Checking printread ##################"
#spectro/printread -v test

echo "############### Checking dispread ##################"
spectro/dispread -v -dfake rgbtest

echo "############### Checking scanin ##################"
scanin/scanin -v -dipn scanin/it8_1150.tif scanin/it8.cht scanin/It8_1150.q60 diag.tif
cp scanin/it8_1150.ti3 .

echo "############### Checking cb2cgats ##################"
profile/cb2cgats -v profile/cbtest cbtest 

echo "############### Checking kodak2cgats ##################"
profile/kodak2cgats -v profile/Ktest Ktest

echo "############### Checking logo2cgats ##################"
profile/logo2cgats -v profile/LogoTest.txt Ltest

echo "############### Checking fakeread ##################"
spectro/fakeread -v profile/3dap.icm test

echo "############### Checking profile ##################"
profile/profile -v -ql -b -kr it8_1150
profile/profile -v -ql -b -kr cbtest
profile/profile -v -ql -b -kr Ktest
profile/profile -v -ql -b -kr Ltest
profile/profile -v -ql -as rgbtest
profile/profile -v -ql -b -kr test

echo "############### Checking fakecmy ##################"
xicc/fakecmy test.$ICM cmytest.ti3
profile/profile -v -ql -b cmytest

echo "############### Checking mppprof ##################"
profile/mpprof -v test

echo "############### Checking revfix ##################"
xicc/revfix -v -1 -r13 -kr test.$ICM testo.$ICM

echo "############### Checking icclink ##################"
link/icclink -v -s -ip -op profile/srgb.icm test.$ICM rgb2test.$ICM
link/icclink -v -g -ip -cmt -dpp profile/srgb.icm test.$ICM rgb2test.$ICM
link/icclink -v -G -ip -cmt -dpp -kr -l280 profile/srgb.icm test.$ICM rgb2test.$ICM

echo "############### Checking cctiff ##################"
imdi/cctiff -v rgb2test.$ICM imdi/rgbtest.tif cmykout.tif

echo "############### Checking icclu ##################"
echo 0.5 0.5 0.5 > temp.txt
icc/icclu -ff -ir -pl profile/srgb.icm < temp.txt

echo "############### Checking xicclu ##################"
echo 0.3 0.3 0.3 0.3 > temp.txt
xicc/xicclu -ff -ir -pl test.$ICM < temp.txt
echo 50.0 0.0 0.0 0.5 > temp.txt
xicc/xicclu -fif -ir -pl test.$ICM < temp.txt

echo "############### Checking mpplu ##################"
echo 0.3 0.3 0.3 0.3 > temp.txt
xicc/mpplu -ff -pl test.mpp < temp.txt

echo "############### Checking greytiff ##################"
imdi/greytiff -v profile/srgb.icm imdi/rgbtest.tif greyout.tif

echo "############### Checking refine ##################"
cd tweak
sanity.ksh
cd ..
 
echo "############### Checking iccgamut ##################"
xicc/iccgamut -ff -ia -w profile/srgb.icm
 
echo "############### Checking tiffgamut ##################"
xicc/tiffgamut -ia -w profile/srgb.icm imdi/rgbtest.tif

echo "############### Checking viewgam ##################"
gamut/viewgam profile/srgb.gam imdi/rgbtest.gam test.wrl

echo "############### Checking iccdump ##################"
icc/iccdump profile/srgb.icm

echo "############### Checking profcheck ##################"
profile/profcheck test.ti3 test.$ICM
 
echo "############### Checking invprofcheck ##################"
profile/invprofcheck test.$ICM

echo "############### Checking splitcgats ##################"
profile/splitcgats test.ti3 test1.ti3 test2.ti3

# render/timage ?

echo "############### Checking mppcheck ##################"
profile/mppcheck test.ti3 test.mpp

echo "############### Checking verify ##################"
cp test.ti1 test2.ti1
spectro/fakeread -v test.$ICM test2
profile/verify test.ti3 test2.ti3
 
echo "############### Checking dispwin ##################"
spectro/dispwin

# spec2cie ?

# Look at:
#   diag.tif
#   cmykout.tif
