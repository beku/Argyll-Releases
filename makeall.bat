echo Simple batch file to invoke Jam in all the subdirectories
echo (Hope you've unzip'd tiff/tiff.zip ?)
rem Note that running xicc twice is a hack to workaround a circular dependency
rem with regard to spectro/libinst.lib
for %%i in (numlib tiff plot icc rspl imdi cgats gamut xicc spectro xicc target scanin profile link tweak) do cd %%i& jam -j%NUMBER_OF_PROCESSORS%& cd ..
