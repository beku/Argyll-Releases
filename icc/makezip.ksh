# Create icclib source distribution
#jam
#zip nt_iccdump.zip iccdump.exe
#zip nt_icclu.zip icclu.exe
rm icclib.zip
zip -9 -ll icclib.zip Readme.txt Licence.txt todo.txt log.txt Jamfile Makefile Makefile.WNT Makefile.IBMNT Makefile.UNIX Makefile.OSX icc.c iccstd.c icc.h iccV42.h iccdump.c icclu.c iccrw.c icctest.c lutest.c
