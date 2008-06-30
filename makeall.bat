echo Simple batch file to invoke Jam from the top
echo (Hope you've unzip'd tiff/tiff.zip ?)
jam -q -fJambase -j%NUMBER_OF_PROCESSORS%
rem If you have trouble with the parallel build, try the
rem version with only one thread.
rem jam -q -fJambase 
