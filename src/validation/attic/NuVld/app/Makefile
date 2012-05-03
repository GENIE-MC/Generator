#
# GENIE makefile fr nuvalidator executables
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

all:    gNuVldXmlRead \
	gNuVldDBUpload \
	gNuVldGUI

gNuVldXmlRead:
	$(CXX) $(CXXFLAGS) -c gNuVldXmlRead.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gNuVldXmlRead.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gnuvld_xmlread

gNuVldDBUpload:
	$(CXX) $(CXXFLAGS) -c gNuVldDBUpload.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gNuVldDBUpload.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gnuvld_dbupload

gNuVldGUI:
	$(CXX) $(CXXFLAGS) -c gNuVldGUI.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gNuVldGUI.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gnuvld_gui


purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gnuvld_xmlread 
	$(RM) $(GENIE_BIN_PATH)/gnuvld_dbupload 
	$(RM) $(GENIE_BIN_PATH)/gnuvld_gui 

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gnuvld_xmlread
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gnuvld_dbupload
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gnuvld_gui

FORCE:

# DO NOT DELETE
