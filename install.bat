echo Simple batch file to invoke Jam in all the subdirectories to install binaries

del /Q bin\*.*
del /Q ref\*.*

for %%i in (numlib tiff libusb libusbw plot icc cgats rspl gamut xicc spectro xicc imdi target scanin profile link tweak render) do cd %%i& jam install& cd ..

rem Create zip archive of documents and exectutables
zip -9 -r bin.zip bin ref -@ < doc\afiles


