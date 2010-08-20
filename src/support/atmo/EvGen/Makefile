####################################################################
#
# Makefile for the GENIE atmospheric event generation tools
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

TGT = gAtmoEvGen

all: $(TGT)

gAtmoEvGen: FORCE 
	$(CXX) $(CXXFLAGS) -c gAtmoEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gAtmoEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_atmo

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevgen_atmo

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_atmo

FORCE:

# DO NOT DELETE
