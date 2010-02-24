####################################################################
# Makefile for the event generation driver customized for INO
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

TGT = gINOevgen

all: $(TGT)

gINOevgen: FORCE
	$(CXX) $(CXXFLAGS) -c gINOEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gINOEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gINOevgen

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gINOevgen 

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gINOevgen 

FORCE:

# DO NOT DELETE
