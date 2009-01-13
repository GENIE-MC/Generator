#
# GENIE basic validation tools makefile
#
# Author: Costas Andreopoulos
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/make/Make.include

# Add extra libraries needed for these executables
GENIE_LIBS   = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES   := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

all:    gvld_sample_comp 

#
# Generated sample comparison utility
#
gVldSampleComp: FORCE
        $(CXX) $(CXXFLAGS) -c gVldSampleComp.cxx $(INCLUDES)
        $(LD) $(LDFLAGS) gVldSampleComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_sample_comp

#
# Cross section testing / comparison utility
#
gVldXSecComp: FORCE
        $(CXX) $(CXXFLAGS) -c gVldXSecComp.cxx $(INCLUDES)
        $(LD) $(LDFLAGS) gVldXSecComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_xsec_comp

#
# Utility checking conservation laws
#
gVldConservationTest: FORCE
        $(CXX) $(CXXFLAGS) -c gVldConservationTest.cxx $(INCLUDES)
        $(LD) $(LDFLAGS) gVldConservationTest.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_conservation_test



clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_sample_comp
	$(RM) $(GENIE_BIN_PATH)/gvld_conservation_test

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_sample_comp
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_conservation_test

FORCE:

# DO NOT DELETE
