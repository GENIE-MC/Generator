#
# GENIE Makefile for ValidationTools/NeutrinoCrossSections
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
#VLDTOOL_LIBS = $(shell $(GENIE)/src/scripts/setup/genie-config --vldlibs)
#LIBRARIES   := $(GENIE_LIBS) $(VLDTOOL_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

all: gVldNuXSecVsWorldData 

#
# Program comparing GENIE neutrino cross-sections against the world data
#
gVldNuXSecVsWorldData: FORCE
	$(CXX) $(CXXFLAGS) -c gVldNuXSecVsWorldData.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldNuXSecVsWorldData.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_nuxsec_vs_world_data

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_nuxsec_vs_world_data

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_nuxsec_vs_world_data

FORCE:

# DO NOT DELETE
