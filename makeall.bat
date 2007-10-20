echo Simple batch file to invoke Jam in all the subdirectories
echo (Hope you've unzip'd tiff/tiff.zip ?)
rem Note that running xicc twice is a hack to workaround a circular dependency
rem with regard to spectro/libinst.lib
for %%i in (numlib tiff libusb libusbw plot icc cgats rspl gamut xicc spectro xicc imdi target scanin profile link tweak render spectro) do cd %%i& jam -f..\Jambase -j%NUMBER_OF_PROCESSORS%& cd ..
