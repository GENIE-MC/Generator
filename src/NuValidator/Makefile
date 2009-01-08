#############################################################
# Makefile for GENIE Neutrino Generator Validator
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#############################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/make/Make.include

BUILD_TARGETS = std-libs dummy-neugen neugen exe

all: $(BUILD_TARGETS)

std-libs: FORCE
	@echo "** Building standard nuvalidator libraries"
	cd XmlParser; make; cd ..; 
	cd DBUtils;   make; cd ..; 
	cd NuVldGUI;  make; cd ..;   

dummy-neugen: FORCE
ifeq ($(strip $(GOPT_ENABLE_NEUGEN)),NO)  
	cd Facades/dummy/; make; cd ../..;  
else
endif

neugen: FORCE
	@echo "** Building nuvalidator facades libraries"
	cd Facades;  make; cd ..;    

libs: FORCE
	cd XmlParser; make; cd ..;   
	cd DBUtils;   make; cd ..;   
	cd Facades;   make; cd ..;   
	cd NuVldGUI;  make; cd ..;    

exe: FORCE
	@echo "** Building standard nuvalidator executables"
	cd app;       make; cd ..; 

purge: FORCE
	cd XmlParser; make purge; cd ..;    
	cd DBUtils;   make purge; cd ..;    
	cd Facades;   make purge; cd ..;    
	cd NuVldGUI;  make purge; cd ..; 

clean: FORCE
	cd XmlParser;     make clean; cd ..;    
	cd DBUtils;       make clean; cd ..;    
	cd Facades;	  make clean; cd ..;    
	cd NuVldGUI;      make clean; cd ..;    
	cd Facades/dummy; make clean; cd ../..; 
	cd app;           make clean; cd ..; 

distclean: FORCE
	cd XmlParser;     make distclean; cd ..;    
	cd DBUtils;       make distclean; cd ..;    
	cd Facades;	  make distclean; cd ..;    
	cd NuVldGUI;      make distclean; cd ..;    
	cd app;           make distclean; cd ..; 
#	cd Facades/dummy; make distclean; cd ../..; 

clean-libs: FORCE
	cd XmlParser;     make clean; cd ..;    
	cd DBUtils;       make clean; cd ..;    
	cd Facades;       make clean; cd ..;    
	cd NuVldGUI;      make clean; cd ..;    
	cd Facades/dummy; make clean; cd ../..; 

clean-exe: FORCE
	cd app;           make clean; cd ..; 

FORCE:

