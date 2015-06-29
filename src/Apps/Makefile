#
# Makefile for GENIE standard apps
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# May 02, 2004
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)

LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT =    gEvGen		  	\
	 gEvDump          	\
	 gEvGenHadronNucleus    \
	 gEvPick	        \
	 gMakeSplines	  	\
	 gSplineAdd   	  	\
	 gSplineXml2Root  	\
	 gMaxPathLengths  	\
	 gNtpConv	  

all: $(TGT)

# gevgen event generation driver
#
gEvGen: FORCE
	$(CXX) $(CXXFLAGS) -c gEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen 

# gevdump event print-out utility
#
gEvDump: FORCE
	$(CXX) $(CXXFLAGS) -c gEvDump.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvDump.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevdump


# gevgen_hadron program generating hadron + nucleus interactions using GENIE's INTRANUKE
#
gEvGenHadronNucleus: FORCE
	$(CXX) $(CXXFLAGS) -c gEvGenHadronNucleus.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvGenHadronNucleus.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_hadron

# gevpick event topology cherry-picking program
#
gEvPick: FORCE
	$(CXX) $(CXXFLAGS) -c gEvPick.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvPick.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevpick

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

# gspl2root utility for converting XML splines into ROOT
#
gSplineXml2Root: FORCE
	$(CXX) $(CXXFLAGS) -c gSplineXml2Root.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gSplineXml2Root.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gspl2root

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


#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(GENIE_BIN_PATH)/gevgen		
	$(RM) $(GENIE_BIN_PATH)/gevdump		
	$(RM) $(GENIE_BIN_PATH)/gevgen_hadron
	$(RM) $(GENIE_BIN_PATH)/gevpick
	$(RM) $(GENIE_BIN_PATH)/gmkspl 	
	$(RM) $(GENIE_BIN_PATH)/gspladd 	
	$(RM) $(GENIE_BIN_PATH)/gspl2root		
	$(RM) $(GENIE_BIN_PATH)/gmxpl		
	$(RM) $(GENIE_BIN_PATH)/gntpc		

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevdump		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_hadron
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevpick
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gmkspl 	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gspladd 	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gspl2root		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gmxpl		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gntpc		

FORCE:

# DO NOT DELETE
