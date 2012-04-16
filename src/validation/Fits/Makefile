#
# GENIE Makefile for fitting apps
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
LIBRARIES   := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

all: 	gTuneRijk \
	gMaQELFit \
	gTuneDISNorm

# 
#
gTuneRijk: FORCE
	$(CXX) $(CXXFLAGS) -c gTuneRijk.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gTuneRijk.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtune_rijk

#
#
gMaQELFit: FORCE
	$(CXX) $(CXXFLAGS) -c gMaQELFit.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gMaQELFit.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtune_maqel

#
#
gTuneDISNorm: FORCE
	$(CXX) $(CXXFLAGS) -c gTuneDISNorm.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gTuneDISNorm.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtune_dis_norm

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gtune_rijk
	$(RM) $(GENIE_BIN_PATH)/gtune_maqel
	$(RM) $(GENIE_BIN_PATH)/gtune_dis_norm

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtune_rijk
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtune_maqel
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtune_dis_norm

FORCE:

# DO NOT DELETE
