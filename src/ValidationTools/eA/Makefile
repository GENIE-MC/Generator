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
#VLDTOOL_LIBS = $(shell $(GENIE)/src/scripts/setup/genie-config --vldlibs)
LIBRARIES   := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

all: 	gVldeQELXSec \
	gVldeRESXSec \
	gVldeHadroAtten

#
# Program comparing GENIE quasi-elastic electron cross sections 
# against data from D.Day's online dbase.
#
gVldeQELXSec: FORCE
	$(CXX) $(CXXFLAGS) -c gVldeQELXSec.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldeQELXSec.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_e_qel_xsec

#
# Program comparing GENIE electron cross sections in the resonance region 
# against data from JLab Hall C online dbase.
#
gVldeRESXSec: FORCE
	$(CXX) $(CXXFLAGS) -c gVldeRESXSec.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldeRESXSec.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_e_res_xsec

#
# Program comparing GENIE electron-nucleus samples with data on hadron attenuation
# from HERMES and JUPITER
#
gVldeHadroAtten: FORCE
	$(CXX) $(CXXFLAGS) -c gVldeHadroAtten.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldeHadroAtten.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_e_hadro_atten


purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_e_qel_xsec
	$(RM) $(GENIE_BIN_PATH)/gvld_e_res_xsec
	$(RM) $(GENIE_BIN_PATH)/gvld_e_hadro_atten

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_e_qel_xsec
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_e_res_xsec
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_e_hadro_atten

FORCE:

# DO NOT DELETE
