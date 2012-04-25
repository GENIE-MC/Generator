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

VLD_SF_SOURCES   := StructFunc.cxx gVldStructFunc.cxx
VLD_SF_OBJECTS    = $(addsuffix .$(ObjSuf), $(basename $(VLD_SF_SOURCES)))

VLD_SUMRULE_SOURCES   := StructFunc.cxx gVldStructFuncSumRuleTest.cxx
VLD_SUMRULE_OBJECTS    = $(addsuffix .$(ObjSuf), $(basename $(VLD_SUMRULE_SOURCES)))

all: gVldStructFunc \
     gVldStructFuncSumRuleTest

#
# Program comparing GENIE structure functions against the world data
#
gVldStructFunc: FORCE
	$(CXX) $(CXXFLAGS) -c $(VLD_SF_SOURCES) $(INCLUDES)
	$(LD) $(LDFLAGS) $(VLD_SF_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_sf

#
# Sum rule checks
#
gVldStructFuncSumRuleTest: FORCE
	$(CXX) $(CXXFLAGS) -c $(VLD_SUMRULE_SOURCES) $(INCLUDES)
	$(LD) $(LDFLAGS) $(VLD_SUMRULE_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_sf_sum_rule_test

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_sf
	$(RM) $(GENIE_BIN_PATH)/gvld_sf_sum_rule_test

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_sf
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_sf_sum_rule_test

FORCE:

# DO NOT DELETE
