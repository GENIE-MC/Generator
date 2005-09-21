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

all-no-neugen:
	cd XmlParser;          make; cd ..;    \
	cd DBUtils;            make; cd ..;    \
	cd NuVldGUI;           make; cd ..;    \
	cd app;                make; cd ..; 

libs:
	cd XmlParser;          make; cd ..;    \
	cd DBUtils;            make; cd ..;    \
	cd Facades;            make; cd ..;    \
	cd NuVldGUI;           make; cd ..;    

libs-no-neugen:
	cd XmlParser;          make; cd ..;    \
	cd DBUtils;            make; cd ..;    \
	cd NuVldGUI;           make; cd ..;    

applications:
	cd app;                make; cd ..; 


purge: FORCE
	cd XmlParser;          make purge; cd ..;    \
	cd DBUtils;            make purge; cd ..;    \
	cd Facades;	       make purge; cd ..;    \
	cd NuVldGUI;           make purge; cd ..; 

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

