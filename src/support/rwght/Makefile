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

TGT = $(GENIE_BIN_PATH)/grwght1scan

all: $(TGT)

gRwght1Scan.o: gRwght1Scan.cxx
	$(CXX) $(CXXFLAGS) -c $^ $(INCLUDES)

$(GENIE_BIN_PATH)/grwght1scan: gRwght1Scan.o
	$(LD) $(LDFLAGS) $^ $(LIBRARIES) -o $@

purge:
	$(RM) *.o *~ core

.PHONY: clean
clean:
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/grwght1scan

.PHONY: distclean
distclean:
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/grwght1scan

FORCE:

# DO NOT DELETE
