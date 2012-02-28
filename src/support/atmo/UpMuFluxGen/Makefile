####################################################################
#
# Makefile for the GENIE upgoing muon flux generator
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#
####################################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = gUpMuFluxGen

all: $(TGT)

gUpMuFluxGen: FORCE 
	$(CXX) $(CXXFLAGS) -c gUpMuFluxGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gUpMuFluxGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_upmu

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevgen_upmu

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_upmu

FORCE:

# DO NOT DELETE
