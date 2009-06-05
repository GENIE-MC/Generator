#
# GENIE SF validation package makefile
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

all:    gVldStructFuncComp 

#
# Structure function comparison utility
#
gVldStructFuncComp: FORCE
	$(CXX) $(CXXFLAGS) -c gVldStructFuncComp.cxx ExtractStructFunc.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldStructFuncComp.o $(LIBRARIES) $(VLDTOOL_LIBS) -o $(GENIE_BIN_PATH)/gvld_structfunc_comp

#
# Structure function comparison utility
#
gVldSumRules: FORCE
	$(CXX) $(CXXFLAGS) -c gVldSumRules.cxx ExtractStructFunc.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldSumRules.o $(LIBRARIES) $(VLDTOOL_LIBS) -o $(GENIE_BIN_PATH)/gvld_sumrules


purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_structfunc_comp

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_structfunc_comp

FORCE:

# DO NOT DELETE
