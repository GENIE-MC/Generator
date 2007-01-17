#########################################################
# Makefile for Neutrino Generator Validator binaries
#
# Costas Andreopoulos, 2003-08-20
#########################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/Make.include

# Add extra libraries needed for these executables
LIBRARIES := $(GENIE_LIBRARIES) $(NUVLD_LIBRARIES) \
             $(LIBRARIES) $(NEUGEN_LIBRARIES) $(CERN_LIBRARIES) 

EXE = $(GENIE_BIN_PATH)/parsexml \
      $(GENIE_BIN_PATH)/filldbase \
      $(GENIE_BIN_PATH)/nuvld \
      $(GENIE_BIN_PATH)/ngevgen 

all: parsexml filldbase nuvld neugen_evgen

parsexml:
	$(CXX) $(CXXFLAGS) -c parsexml.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) parsexml.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/parsexml

filldbase:
	$(CXX) $(CXXFLAGS) -c filldbase.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) filldbase.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/filldbase

nuvld:
	$(CXX) $(CXXFLAGS) -c nuvld_gui.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) nuvld_gui.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/nuvld_gui

neugen_evgen:
	$(CXX) $(CXXFLAGS) -c neugen_evgen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) neugen_evgen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/ngevgen

clean: FORCE
	$(RM) *.o *~ core $(EXE)

FORCE:

# DO NOT DELETE
