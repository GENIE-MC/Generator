#########################################################
# Makefile for Neutrino Generator Validator
#
# Costas Andreopoulos, 2003-08-20
#########################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

all:
	cd XmlParser;          make; cd ..;    \
	cd DBUtils;            make; cd ..;    \
	cd Facades;            make; cd ..;    \
	cd NuVldGUI;           make; cd ..;    \
	cd app;                make; cd ..; 

libs:
	cd XmlParser;          make; cd ..;    \
	cd DBUtils;            make; cd ..;    \
	cd Facades;            make; cd ..;    \
	cd NuVldGUI;           make; cd ..;    

applications:
	cd app;                make; cd ..; 


clean: FORCE
	cd XmlParser;          make clean; cd ..;    \
	cd DBUtils;            make clean; cd ..;    \
	cd Facades;	       make clean; cd ..;    \
	cd NuVldGUI;           make clean; cd ..;    \
	cd app;                make clean; cd ..; 

clean-libs: FORCE
	cd XmlParser;          make clean; cd ..;    \
	cd DBUtils;            make clean; cd ..;    \
	cd Facades;            make clean; cd ..;    \
	cd NuVldGUI;           make clean; cd ..;    

clean-applications: FORCE
	cd app;                make clean; cd ..; 


FORCE:

