#
# Makefile for Neutrino Generator Validator binaries
#
# Author: Costas Andreopoulos
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/make/Make.include

# Add extra libraries needed for these executables
LIBRARIES := $(GENIE_LIBRARIES) $(NUVLD_LIBRARIES) \
             $(LIBRARIES) $(NEUGEN_LIBRARIES) $(CERN_LIBRARIES) 

all: parse_xml upload_dbase gui 

parse_xml:
	$(CXX) $(CXXFLAGS) -c parsexml.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) parsexml.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gnuvld-parsexml

upload_dbase:
	$(CXX) $(CXXFLAGS) -c filldbase.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) filldbase.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gnuvld-dbupload

gui:
	$(CXX) $(CXXFLAGS) -c nuvld_gui.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) nuvld_gui.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gnuvld-gui

neugen_evgen:
	$(CXX) $(CXXFLAGS) -c neugen_evgen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) neugen_evgen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gneugen-evgen

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gnuvld-parsexml 
	$(RM) $(GENIE_BIN_PATH)/gnuvld-dbupload 
	$(RM) $(GENIE_BIN_PATH)/gnuvld-gui 
	$(RM) $(GENIE_BIN_PATH)/gneugen-evgen

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gnuvld-parsexml 
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gnuvld-dbupload 
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gnuvld-gui 
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gneugen-evgen

FORCE:

# DO NOT DELETE
