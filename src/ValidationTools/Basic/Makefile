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

all:    gVldNuEvComp \
	gVldNuEvScan \
	gVldNuXsComp 

#
# Neutrino event sample comparison utility
#
gVldNuEvComp: FORCE
	$(CXX) $(CXXFLAGS) -c gVldNuEvComp.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldNuEvComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_nuev_comp

#
# Utility scanning the input file and checking for conservation law violations,
# strange/undecayed particles etc
#
gVldNuEvScan: FORCE
	$(CXX) $(CXXFLAGS) -c gVldNuEvScan.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldNuEvScan.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_nuev_scan

#
# Neutrino cross section comparison utility
#
gVldNuXsComp: FORCE
	$(CXX) $(CXXFLAGS) -c gVldNuXsComp.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldNuXsComp.o $(LIBRARIES) $(VLDTOOL_LIBS) -o $(GENIE_BIN_PATH)/gvld_nuxs_comp


purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_nuev_comp
	$(RM) $(GENIE_BIN_PATH)/gvld_nuev_scan
	$(RM) $(GENIE_BIN_PATH)/gvld_nuxs_comp

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_nuev_comp
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_nuev_scan
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_nuxs_comp

FORCE:

# DO NOT DELETE
