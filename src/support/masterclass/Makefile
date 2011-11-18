####################################################################
# Makefile for the neutrino master class app `gnumc'
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

TGT = gnumc

MASTERCLASS_APP_SOURCES := $(wildcard *.$(SrcSuf)) _ROOT_DICT.cc
MASTERCLASS_APP_OBJECTS  = $(addsuffix .$(ObjSuf), $(basename $(MASTERCLASS_APP_SOURCES)))
DICTIONARY_HEADERS       = GNuMcMainFrame.h FastSimCherenkov.h FastSimScintCalo.h MCTruthDisplay.h LinkDef.h

all: $(TGT)

gnumc: _ROOT_DICT.$(ObjSuf)
	$(CXX) $(CXXFLAGS) -c $(MASTERCLASS_APP_SOURCES) $(INCLUDES)
	$(LD) $(LDFLAGS) $(MASTERCLASS_APP_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gnumc

purge: FORCE
	$(RM) *.o *~ core _ROOT_DICT.*

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gnumc _ROOT_DICT.*

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gnumc

_ROOT_DICT.cc: $(DICTIONARY_HEADERS) 
	@echo "Generating dictionary $@..."
	rootcint -f $@ -c -I$(ROOTSYS)/include -I$(GENIE)/src $^

FORCE:

# DO NOT DELETE
