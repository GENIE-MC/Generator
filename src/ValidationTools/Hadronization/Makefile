#
# Makefile for the GENIE hadronization validation package
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/make/Make.include

GENIE_LIBS   = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES   := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

FREE_NUC_HADRONZ_TEST_SOURCES   := HadPlots.cxx HadPlotter.cxx gVldFreeNucleonHadronzTest.cxx
FREE_NUC_HADRONZ_TEST_OBJECTS    = $(addsuffix .$(ObjSuf), $(basename $(FREE_NUC_HADRONZ_TEST_SOURCES)))

MEDIUM_EFF_HADRONZ_TEST_SOURCES := gVldMediumEffectsHadronzTest.cxx
MEDIUM_EFF_HADRONZ_TEST_OBJECTS  = $(addsuffix .$(ObjSuf), $(basename $(MEDIUM_EFF_HADRONZ_TEST_SOURCES)))

all: gVldFreeNucleonHadronzTest

gVldFreeNucleonHadronzTest: FORCE
	$(CXX) $(CXXFLAGS) -c $(FREE_NUC_HADRONZ_TEST_SOURCES) $(INCLUDES)
	$(LD) $(LDFLAGS) $(FREE_NUC_HADRONZ_TEST_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_free_nucleon_hadronz_test

gVldMediumEffectsHadronzTest: FORCE
	$(CXX) $(CXXFLAGS) -c $(MEDIUM_EFF_HADRONZ_TEST_SOURCES) $(INCLUDES)
	$(LD) $(LDFLAGS) $(MEDIUM_EFF_HADRONZ_TEST_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_medium_effects_hadronz_test

purge: FORCE
	$(RM) *.o *~ core

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(GENIE_BIN_PATH)/gvld_free_nucleon_hadronz_test
	$(RM) $(GENIE_BIN_PATH)/gvld_medium_effects_hadronz_test

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_free_nucleon_hadronz_test
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_medium_effects_hadronz_test

FORCE:

# DO NOT DELETE


