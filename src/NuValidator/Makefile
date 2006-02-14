#########################################################
# Makefile for Neutrino Generator Validator
#
# Costas Andreopoulos, 2003-08-20
#########################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/Make.include

ifndef GOPT_NEUGEN
  BUILD_TARGETS := dummy-neugen
else
  BUILD_TARGETS := neugen
endif
BUILD_TARGETS := $(BUILD_TARGETS) std

all: $(BUILD_TARGETS)

all:
	cd XmlParser;          make; cd ..;    \
	cd DBUtils;            make; cd ..;    \
	cd Facades;            make; cd ..;    \
	cd NuVldGUI;           make; cd ..;    \
	cd app;                make; cd ..; 

std:
	cd XmlParser;          make; cd ..;    \
	cd DBUtils;            make; cd ..;    \
	cd NuVldGUI;           make; cd ..;    \
	cd app;                make; cd ..; 

dummy-neugen:
	cd Facades/dummy/;     make; cd ../..;  

neugen:
	cd Facades;            make; cd ..;    \

libs:
	cd XmlParser;          make; cd ..;    \
	cd DBUtils;            make; cd ..;    \
	cd Facades;            make; cd ..;    \
	cd NuVldGUI;           make; cd ..;    

exe:
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
	cd Facades/dummy;      make clean; cd ../..; \
	cd app;                make clean; cd ..; 

clean-libs: FORCE
	cd XmlParser;          make clean; cd ..;    \
	cd DBUtils;            make clean; cd ..;    \
	cd Facades;            make clean; cd ..;    \
	cd NuVldGUI;           make clean; cd ..;    \
	cd Facades/dummy;      make clean; cd ../..; 

clean-exe: FORCE
	cd app;                make clean; cd ..; 


FORCE:

