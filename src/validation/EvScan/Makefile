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

all:    gVldSampleScan 

#
# Utility scanning the input file and checking for conservation law violations,
# strange/undecayed particles etc
#
gVldSampleScan: FORCE
	$(CXX) $(CXXFLAGS) -c gVldSampleScan.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldSampleScan.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_sample_scan

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_sample_scan

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_sample_scan

FORCE:

# DO NOT DELETE
