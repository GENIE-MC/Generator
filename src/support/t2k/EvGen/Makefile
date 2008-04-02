####################################################################
# Makefile for the event generation driver customized for T2K
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

TGT = gT2Kevgen

all: $(TGT)

gT2Kevgen: FORCE
	$(CXX) $(CXXFLAGS) -c gT2KEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gT2KEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gT2Kevgen

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gT2Kevgen $(GENIE_BIN_PATH)/gT2Kgevgen

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gT2Kevgen $(GENIE_BIN_INSTALLATION_PATH)/gT2Kgevgen

FORCE:

# DO NOT DELETE
