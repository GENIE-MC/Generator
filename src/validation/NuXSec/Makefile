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

VLD_NUXSEC_SOURCES   := NuXSecData.cxx NuXSecFunc.cxx NuXSecComparison.cxx gVldNuXSecVsWorldData.cxx
VLD_NUXSEC_OBJECTS    = $(addsuffix .$(ObjSuf), $(basename $(VLD_NUXSEC_SOURCES)))

all: gVldNuXSecVsWorldData 

#
# Program comparing GENIE neutrino cross-sections against the world data
#
gVldNuXSecVsWorldData: FORCE
	$(CXX) $(CXXFLAGS) -c $(VLD_NUXSEC_SOURCES) $(INCLUDES)
	$(LD) $(LDFLAGS) $(VLD_NUXSEC_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_nuxsec_vs_world_data

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_nuxsec_vs_world_data

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_nuxsec_vs_world_data

FORCE:

# DO NOT DELETE
