####################################################################
# Makefile for the event generation driver customized for LBNE
#   Atmospheric neutrino studies  
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# Hugh Gallagher <hugh.gallagher \at tufts.edu>
####################################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = gLBNEAtmosEvGen

all: $(TGT)

gLBNEAtmosEvGen: FORCE 
	$(CXX) $(CXXFLAGS) -c gLBNEAtmosEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gLBNEAtmosEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gLBNEAtmosEvGen

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gLBNEAtmosEvGen

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gLBNEAtmosEvGen

FORCE:

# DO NOT DELETE
