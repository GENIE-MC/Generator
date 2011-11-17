####################################################################
# Makefile for the nucleon-decay event generation driver 
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

TGT = gevgen_ndcy

NUCLEON_DECAY_CODE_SOURCES := $(wildcard *.cxx)
NUCLEON_DECAY_CODE_OBJECTS  = $(addsuffix .$(ObjSuf), $(basename $(NUCLEON_DECAY_CODE_SOURCES)))

all: $(TGT)

gevgen_ndcy: FORCE
	$(CXX) $(CXXFLAGS) -c $(NUCLEON_DECAY_CODE_SOURCES) $(INCLUDES)
	$(LD) $(LDFLAGS) $(NUCLEON_DECAY_CODE_OBJECTS) $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_ndcy

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevgen_ndcy 

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_ndcy

FORCE:

# DO NOT DELETE
