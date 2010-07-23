#########################################################
# Makefile for the GENIE Event Server
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#########################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT =    gevserv

all: $(TGT)

gevserv: FORCE
	$(CXX) $(CXXFLAGS) -c gEvServ.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvServ.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevserv

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevserv 

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevserv 

FORCE:

# DO NOT DELETE
