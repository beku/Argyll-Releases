# default Makefile to involke Jam

all::
	jam -q -fJambase -j 3

clean::
	jam clean

install:;
	jam -q -fJambase -j 3 install
