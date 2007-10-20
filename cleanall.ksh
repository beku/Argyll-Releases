#!/bin/sh
echo "Simple batch file to clean all the subdirectories"
for i in `cat blddirs`
do
(echo ------------; echo $i; echo ----------; cd $i; rm *.obj *.lib *.exe *.o *.a ; cd ..; echo)
done
