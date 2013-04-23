#
# GENIE Makefile for event sample and logfile anomaly scanners
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

all:    gVldSampleScan gVldRepeatabilityTest

gVldSampleScan: FORCE
	$(CXX) $(CXXFLAGS) -c gVldSampleScan.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldSampleScan.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_sample_scan

gVldRepeatabilityTest: FORCE
	$(CXX) $(CXXFLAGS) -c gVldRepeatabilityTest.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldRepeatabilityTest.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_repeatability_test

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_sample_scan
	$(RM) $(GENIE_BIN_PATH)/gvld_repeatability_test

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_sample_scan
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_repeatability_test

FORCE:

# DO NOT DELETE
