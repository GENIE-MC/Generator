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

HADRONIZATION_TEST_SOURCES := $(wildcard *.cxx)
HADRONIZATION_TEST_OBJECTS  = $(addsuffix .$(ObjSuf), $(basename $(HADRONIZATION_TEST_SOURCES)))

all: gVldHadronizationTest

gVldHadronizationTest: FORCE
	$(CXX) $(CXXFLAGS) -c $(HADRONIZATION_TEST_SOURCES) $(INCLUDES)
	$(LD) $(LDFLAGS) $(HADRONIZATION_TEST_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_hadronization_test

purge: FORCE
	$(RM) *.o *~ core

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(GENIE_BIN_PATH)/gvld_hadronization_test

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_hadronization_test

FORCE:

# DO NOT DELETE


