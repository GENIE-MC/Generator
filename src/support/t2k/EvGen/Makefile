####################################################################
# Makefile for the event generation driver customized for T2K
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
####################################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = gevgen_t2k

all: $(TGT)

gevgen_t2k: FORCE
	$(CXX) $(CXXFLAGS) -c gT2KEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gT2KEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_t2k
	@echo "For historical purposes, add symbolic link to the T2K event generation app"
	ln -sf $(GENIE_BIN_PATH)/gevgen_t2k $(GENIE_BIN_PATH)/gT2Kevgen 

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gT2Kevgen $(GENIE_BIN_PATH)/gevgen_t2k

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_t2k

FORCE:

# DO NOT DELETE
