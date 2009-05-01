#
# GENIE basic validation tools makefile
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

# Add extra libraries needed for these executables
#
GENIE_LIBS   = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LDTOOL_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --vldlibs)
LIBRARIES   := $(GENIE_LIBS) $(LDTOOL_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

all: 	gCohTest

# 
#
gCohTest: FORCE
	$(CXX) $(CXXFLAGS) -c gCohTest.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gCohTest.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_cohtest

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_cohtest

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_cohtest

FORCE:

# DO NOT DELETE
