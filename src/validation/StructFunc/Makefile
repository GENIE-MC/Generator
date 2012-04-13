#
# GENIE makefile for ValidationTools/StructFunc package
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

#all: gVldStructFuncVsWorldData \
#     gVldStructFuncSumRuleTest

all: 

#
# Program comparing GENIE structure functions against the world data
#
gVldStructFuncVsWorldData: FORCE
	$(CXX) $(CXXFLAGS) -c gVldStructFuncVsWorldData.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldStructFuncVsWorldData.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_sf_vs_world_data

#
# Sum rule checks
#
gVldStructFuncSumRuleTest: FORCE
	$(CXX) $(CXXFLAGS) -c gVldStructFuncSumRuleTest.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldStructFuncSumRuleTest.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_sf_sum_rule_test

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_sf_vs_world_data
	$(RM) $(GENIE_BIN_PATH)/gvld_sf_sum_rule_test

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_sf_vs_world_data
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_sf_sum_rule_test

FORCE:

# DO NOT DELETE
