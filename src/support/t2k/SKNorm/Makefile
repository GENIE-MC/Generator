####################################################################
# Makefile for the gSKXSecTable T2K utility
# Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
####################################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = gSKXSecTable

all: $(TGT)

gSKXSecTable: FORCE
	$(CXX) $(CXXFLAGS) -c gSKXSecTable.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gSKXSecTable.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gSKXSecTable

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gSKXSecTable

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gSKXSecTable

FORCE:

# DO NOT DELETE
