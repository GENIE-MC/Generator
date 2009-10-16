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
LIBRARIES   := $(GENIE_LIBS) $(VLDTOOL_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

all: gVldeRESXSec

#
# Program comparing GENIE neutrino cross sections against the world data
#
gVldeRESXSec: FORCE
	$(CXX) $(CXXFLAGS) -c gVldeRESXSec.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gVldeRESXSec.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gvld_e_res_xsec


purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gvld_e_res_xsec

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gvld_e_res_xsec

FORCE:

# DO NOT DELETE
