# GENIE Merenyi test makefile
#
# Author: Pauli Kehayias
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/make/Make.include

GENIE_LIBS   = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES   := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

MERENYI_TEST_SOURCES := $(wildcard *.cxx)
MERENYI_TEST_OBJECTS  = $(addsuffix .$(ObjSuf), $(basename $(MERENYI_TEST_SOURCES)))

all: gVldMerenyiTest

gVldMerenyiTest: FORCE
	$(CXX) $(CXXFLAGS) -c $(MERENYI_TEST_SOURCES) $(INCLUDES)
	$(LD) $(LDFLAGS) $(MERENYI_TEST_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_merenyi_test

purge: FORCE
	$(RM) *.o *~ core

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(GENIE_BIN_PATH)/gvld_merenyi_test

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_merenyi_test

FORCE:

# DO NOT DELETE


