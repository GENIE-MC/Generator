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
VLDTOOL_LIBS = $(shell $(GENIE)/src/scripts/setup/genie-config --vldlibs)
LIBRARIES   := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

all:    gVldSampleComp \
	gVldSampleScan \
	gVldXSecComp 

#
# Neutrino event sample comparison utility
#
gVldSampleComp: FORCE
	$(CXX) $(CXXFLAGS) -c gVldSampleComp.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldSampleComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_sample_comp

#
# Utility scanning the input file and checking for conservation law violations,
# strange/undecayed particles etc
#
gVldSampleScan: FORCE
	$(CXX) $(CXXFLAGS) -c gVldSampleScan.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldSampleScan.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_sample_scan

#
# Neutrino cross section comparison utility
#
gVldXSecComp: FORCE
	$(CXX) $(CXXFLAGS) -c gVldXSecComp.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldXSecComp.o $(LIBRARIES) $(VLDTOOL_LIBS) -o $(GENIE_BIN_PATH)/gvld_xsec_comp


purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_sample_comp
	$(RM) $(GENIE_BIN_PATH)/gvld_sample_scan
	$(RM) $(GENIE_BIN_PATH)/gvld_xsec_comp

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_sample_comp
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_sample_scan
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_xsec_comp

FORCE:

# DO NOT DELETE
