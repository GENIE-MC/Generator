#
# GENIE Makefile for validation/NuXSec/
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

VLD_NUXSEC_SOURCES   := NuXSecData.cxx NuXSecFunc.cxx NuXSecComparison.cxx gVldNuXSec.cxx
VLD_NUXSEC_OBJECTS    = $(addsuffix .$(ObjSuf), $(basename $(VLD_NUXSEC_SOURCES)))

all: gVldNuXSec

#
# Program comparing GENIE neutrino cross-sections against the world data
#
gVldNuXSec: FORCE
	$(CXX) $(CXXFLAGS) -c $(VLD_NUXSEC_SOURCES) $(INCLUDES)
	$(LD) $(LDFLAGS) $(VLD_NUXSEC_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_nu_xsec

#
# Program tuning the GENIE non-resonance backgrounds
#
gTuneRijk: FORCE
	$(CXX) $(CXXFLAGS) -c gTuneRijk.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gTuneRijk.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtune_rijk

#
# Program tuning the GENIE DIS normalization factor
#
gTuneDISNorm: FORCE
	$(CXX) $(CXXFLAGS) -c gTuneDISNorm.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gTuneDISNorm.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtune_dis_norm

#
# Program fitting Ma_{CCQE}
#
gTuneMaQE: FORCE
	$(CXX) $(CXXFLAGS) -c gTuneMaQE.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gTuneMaQE.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtune_ma_qe


purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_nu_xsec
	$(RM) $(GENIE_BIN_PATH)/gtune_rijk
	$(RM) $(GENIE_BIN_PATH)/gtune_dis_norm
	$(RM) $(GENIE_BIN_PATH)/gtune_ma_qe

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_nu_xsec
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtune_rijk
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtune_dis_norm
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtune_ma_qe

FORCE:

# DO NOT DELETE
