#
# Makefile for the GENIE event reweighting applications 
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
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = gRwght1Scan 

all: $(TGT)

gRwght1Scan: FORCE
	$(CXX) $(CXXFLAGS) -c gRwght1Scan.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gRwght1Scan.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/grwght1scan

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/grwght1scan

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/grwght1scan

FORCE:

# DO NOT DELETE
