####################################################################
# Makefile for the event generation driver customized for NuMI
# beam-line experiments
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
####################################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = gNuMIevgen

all: $(TGT)

gNuMIevgen: FORCE
	$(CXX) $(CXXFLAGS) -c gNuMIExptEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gNuMIExptEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gNuMIevgen

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gNuMIevgen $(GENIE_BIN_PATH)/gNuMIevgen

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gNuMIevgen $(GENIE_BIN_INSTALLATION_PATH)/gNuMIevgen

FORCE:

# DO NOT DELETE
