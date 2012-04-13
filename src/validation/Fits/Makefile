#
# GENIE basic validation tools makefile
#
# Author: Costas Andreopoulos
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

all: 	gRijkFit \
	gMaQELFit \
	gDISNormFit

# 
#
gRijkFit: FORCE
	$(CXX) $(CXXFLAGS) -c gRijkFit.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gRijkFit.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtune_rijk

#
#
gMaQELFit: FORCE
	$(CXX) $(CXXFLAGS) -c gMaQELFit.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gMaQELFit.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtune_maqel

#
#
gDISNormFit: FORCE
	$(CXX) $(CXXFLAGS) -c gDISNormFit.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gDISNormFit.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtune_disnorm

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gtune_rijk
	$(RM) $(GENIE_BIN_PATH)/gtune_maqel
	$(RM) $(GENIE_BIN_PATH)/gtune_disnorm

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtune_rijk
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtune_maqel
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtune_disnorm

FORCE:

# DO NOT DELETE
