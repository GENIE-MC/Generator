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
GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
NUVLD_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --nuvldlibs)
LIBRARIES  := $(GENIE_LIBS) $(NUVLD_LIBS) $(LIBRARIES)  $(NEUGEN_LIBRARIES) $(CERN_LIBRARIES)

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

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gnuvld-parsexml 
	$(RM) $(GENIE_BIN_PATH)/gnuvld-dbupload 
	$(RM) $(GENIE_BIN_PATH)/gnuvld-gui 

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gnuvld-parsexml 
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gnuvld-dbupload 
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gnuvld-gui 

FORCE:

# DO NOT DELETE
