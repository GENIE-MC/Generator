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

HADRONZ_TEST_SOURCES   := HadPlots.cxx HadPlotter.cxx gVldHadronzTest.cxx
HADRONZ_TEST_OBJECTS    = $(addsuffix .$(ObjSuf), $(basename $(HADRONZ_TEST_SOURCES)))

all:    gVldHadronzTest 

gVldHadronzTest: FORCE
	$(CXX) $(CXXFLAGS) -c $(HADRONZ_TEST_SOURCES) $(INCLUDES)
	$(LD) $(LDFLAGS) $(HADRONZ_TEST_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_hadronz_test

purge: FORCE
	$(RM) *.o *~ core

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(GENIE_BIN_PATH)/gvld_hadronz_test

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_hadronz_test

FORCE:

# DO NOT DELETE


