#
# Makefile for GENIE test programs
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# May 02, 2004
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)


TGT = gtestINukeHadroXSec

all: $(TGT)

gtestINukeHadroXSec: FORCE
	$(CXX) $(CXXFLAGS) -c gtestINukeHadroXSec.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestINukeHadroXSec.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestINukeHadroXSec

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gtestINukeHadroXSec 	

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestINukeHadroXSec 	

FORCE:

# DO NOT DELETE
