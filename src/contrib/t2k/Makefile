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

TGT =   gSKMergeData \
	gT2KRew1SystematicScan

all: $(TGT)

gSKMergeData: FORCE
	$(CXX) $(CXXFLAGS) -c gSKMergeData.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gSKMergeData.o $(LIBRARIES) -o gSKmergedata
#	$(LD) $(LDFLAGS) gSKMergeData.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gSKmergedata

gT2KRew1SystematicScan: FORCE
	$(CXX) $(CXXFLAGS) -c gT2KRew1SystematicScan.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gT2KRew1SystematicScan.o $(LIBRARIES) -o gT2Krew1systscan
#	$(LD) $(LDFLAGS) gT2KRew1SystematicScan.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gT2Krew1systscan

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core gSKmergedata
	$(RM) *.o *~ core gT2Krew1systscan
#	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gSKmergedata
#	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gT2Krew1systscan

distclean: FORCE
#	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gSKmergedata
#	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gT2Krew1systscan

FORCE:

# DO NOT DELETE
