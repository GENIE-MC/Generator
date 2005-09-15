#########################################################
# Makefile for GENIE Neutrino Generator
#
# Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
# Last update: 2003-12-23
#########################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/Make.include

# For building this package add GEANT and CLHEP includes to the include path
INCLUDES := $(INCLUDES) $(GEANT4_INCLUDES) $(CLHEP_INCLUDES)

GENIES	  := $(wildcard *.cxx)
CINTLIST  := $(addsuffix .h, $(basename $(GENIES)))
CINTLISTO := $(addsuffix .$(ObjSuf), $(basename $(CINTLIST)))
#GENIES   := $(GENIES) _ROOT_DICT_GeantInterface.cc $(wildcard *.cpp)
GENIEO 	   = $(addsuffix .$(ObjSuf), $(basename $(GENIES)))
GENIEbS	  := $(wildcard *.cxx)
GENIEbO    = $(addsuffix .$(ObjSuf), $(basename $(GENIEbS)))
GENIESO   := $(GENIE_GEANT_INTERFACE_SHARED_LIB)

#all:	rootcint lib
all:	lib

rootcint: _ROOT_DICT_GeantInterface.o

_ROOT_DICT_GeantInterface.cc: LinkDef.h $(CINTLISTO)
	$(RM) _ROOT_DICT_GeantInterface.*
	@echo "Generating ROOT dictionary ..."
	rootcint -f _ROOT_DICT_GeantInterface.cc -c				\
		-I$(ROOTSYS)/include -I$(GENIE)/src $(CINTLIST) LinkDef.h

lib: FORCE
	@echo "Generating shared library ..."
	@$(MAKE) -f $(MAKEFILE) $(GENIESO)

$(GENIESO): $(GENIEO)
	$(SOCMD) $(SOFLAGS) $(SOMINF)				\
		$(GENIEO) $(OutPutOpt) $(GENIESO)

purge: FORCE
	$(RM) $(GENIEbO) *~ core a.out Makefile.temp

clean: FORCE
	$(RM) *.o *~ core a.out $(GENIESO)		        \
		_ROOT_DICT_GeantInterface.* Makefile.temp

FORCE:

