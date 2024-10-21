#
# Makefile for GENIE test programs
#
# Costas Andreopoulos <c.andreopoulos \at cern.ch>
# May 02, 2004
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
GENIE_LIBS += -lGReWeight
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)


TGT = gRWIOExample1

all: $(TGT)

gRWIOExample1: FORCE
	$(CXX) $(CXXFLAGS) -c gRWIOExample1.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gRWIOExample1.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gRWIOExample1

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gRWIOExample1

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gRWIOExample1 	

FORCE:

# DO NOT DELETE
