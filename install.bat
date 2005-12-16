echo Simple batch file to invoke Jam in all the subdirectories to install binaries

del /Q bin\*.*
del /Q ref\*.*

for %%i in (numlib tiff plot icc rspl imdi cgats gamut xicc spectro target scanin profile link tweak) do cd %%i& jam install& cd ..

rem Create zip archive of documents and exectutables
zip -9 -r bin.zip bin ref -@ < doc\afiles


