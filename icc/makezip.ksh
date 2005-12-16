# Create icclib source distribution
#jam
#zip nt_iccdump.zip iccdump.exe
#zip nt_icclu.zip icclu.exe
rm icclib.zip
zip -9 -ll icclib.zip Readme.txt Licence.txt icc_readme.html icc_problems.html todo.txt Jamfile Makefile Makefile.WNT Makefile.IBMNT Makefile.UNIX Makefile.OSX fbtest.c icc.c iccstd.c icc.h iccV42.h iccdump.c icclu.c iccrw.c icctest.c lutest.c
