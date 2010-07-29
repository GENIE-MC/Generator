#
# Makefile for the GENIE utility programs
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

TGT = gSKMergeData

all: $(TGT)

gSKMergeData: FORCE
	$(CXX) $(CXXFLAGS) -c gSKMergeData.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gSKMergeData.o $(LIBRARIES) -o gSKmergedata
#	$(LD) $(LDFLAGS) gSKMergeData.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gSKmergedata

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core gSKmergedata
#	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gSKmergedata

distclean: FORCE
#	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gSKmergedata

FORCE:

# DO NOT DELETE
