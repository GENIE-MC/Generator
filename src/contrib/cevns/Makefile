#
# Makefile for GENIE Apps
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)

# Mac OS X -bind_at_load means we must list the low level libraries
# first, and each stage must fully resolved before the next library
LIBRARIES  := $(LIBRARIES) $(CERN_LIBRARIES) $(GENIE_LIBS)

all: gtestXSecCEvNS	

gtestXSecCEvNS: gtestXSecCEvNS.o $(call find_libs,gevgen)
	@echo "** Building gtestXSecCEvNS"
	$(LD) $(LDFLAGS) gtestXSecCEvNS.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestXSecCEvNS

# CLEANING-UP

purge: FORCE
	$(RM) *.o *~ *.d core

clean: FORCE
	$(RM) *.o *~ *.d core
	$(RM) $(TGT)

FORCE:

# DO NOT DELETE

# this should be last because it includes new targets...
-include $(wildcard *.d)

