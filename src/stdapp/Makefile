#########################################################
# Makefile for Neutrino Generator Test Programs
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# May 02, 2004
#########################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
#NUVLD_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --nuvldlibs)
#LIBRARIES  := $(GENIE_LIBS) $(NUVLD_LIBS) $(LIBRARIES)  $(NEUGEN_LIBRARIES) $(CERN_LIBRARIES)

LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT =    gEvGen		  \
	 ghAEvGen         \
	 gEvGenGUI	  \
	 gMakeSplines	  \
	 gSplineAdd   	  \
	 gSplinePlotter	  \
	 gMaxPathLengths  \
	 gNtpConv	  \
 	 gNC2ElectronScat 

all: $(TGT)

# gevgen event generation driver
#
gEvGen: FORCE
	$(CXX) $(CXXFLAGS) -c gEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen 


# ghAgen program generating hadron + nucleus interactions using GENIE's INTRANUKE
#
ghAEvGen: FORCE
	$(CXX) $(CXXFLAGS) -c ghAEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) ghAEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/ghAevgen

# A gevgen GUI (in devel)
#
gEvGenGUI: FORCE
ifeq ($(strip $(GOPT_ENABLE_VIEWER)),YES)
	$(CXX) $(CXXFLAGS) -c gEvGenGUI.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvGenGUI.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_gui
else
	@echo "You need to enable the Viewer package to build the gevgen_gui utility"
endif

# gmkspl x-section spline building utility
#
gMakeSplines: FORCE
	$(CXX) $(CXXFLAGS) -c gMakeSplines.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gMakeSplines.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmkspl 

# gspladd utility merging x-section spline xml files
#
gSplineAdd: FORCE
	$(CXX) $(CXXFLAGS) -c gSplineAdd.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gSplineAdd.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gspladd 

# gsplt spline plotter utility
#
gSplinePlotter: FORCE
	$(CXX) $(CXXFLAGS) -c gSplinePlotter.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gSplinePlotter.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gsplt

# gmxpl utility computing maximum path lengths for a given root geometry
#
gMaxPathLengths: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c gMaxPathLengths.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gMaxPathLengths.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmxpl
else
	@echo "You need to enable the flux drivers to build the gmxpl utility"
endif

# gntpc ntuple conversion utility
#
gNtpConv: FORCE
	$(CXX) $(CXXFLAGS) -c gNtpConv.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gNtpConv.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gntpc 

# gnc2eN program converting NC events to electron scattering events
#
gNC2ElectronScat: FORCE
	$(CXX) $(CXXFLAGS) -c gNC2ElectronScat.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gNC2ElectronScat.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gnc2eN


#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(GENIE_BIN_PATH)/gevgen		
	$(RM) $(GENIE_BIN_PATH)/ghAevgen
	$(RM) $(GENIE_BIN_PATH)/gevgen_gui
	$(RM) $(GENIE_BIN_PATH)/gmkspl 	
	$(RM) $(GENIE_BIN_PATH)/gspladd 	
	$(RM) $(GENIE_BIN_PATH)/gsplt		
	$(RM) $(GENIE_BIN_PATH)/gmxpl		
	$(RM) $(GENIE_BIN_PATH)/gntpc		
	$(RM) $(GENIE_BIN_PATH)/gnc2eN  

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/ghAevgen	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_gui
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gmkspl 	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gspladd 	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gsplt		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gmxpl		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gntpc		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gnc2eN		

clean--gEvGen: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevgen

clean--ghAEvGen: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/ghAevgen

clean--gEvGenGUI: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevgen_gui

clean--gMakeSplines: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gmkspl

clean--gSplineAdd: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gspladd

clean--gSplinePlotter: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gsplt

clean--gMaxPathLengths: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gmxpl

clean--gNtpConv: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gntpc

clean-gNC2ElectronScat: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gnc2eN


FORCE:

# DO NOT DELETE
