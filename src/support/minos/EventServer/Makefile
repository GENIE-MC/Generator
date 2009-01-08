#########################################################
# Makefile for the MINOS GENIE Event Server
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

TGT =    gevserv_minos

all: $(TGT)

gevserv_minos: FORCE
	$(CXX) $(CXXFLAGS) -c gEvServMinos.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvServMinos.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevserv_minos

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevserv_minos 

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevserv_minos 

FORCE:

# DO NOT DELETE
