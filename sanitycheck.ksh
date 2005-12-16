#!/bin/sh
## Do a simple sanity check of the major toolchain
target/targen -v -d4 -f500 -l280 test
#spectro/printread -v test
spectro/fakeread -v profile/3dap.icm test
profile/profile -v -ql -b -kr test
link/icclink -v -G -i4 -c2 -d0 -kr -l280 profile/srgb.icm test.icm rgb2test.icm
imdi/cctiff -v rgb2test.icm imdi/rgbtest.tif cmykout.tif

 
