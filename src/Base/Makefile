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

GENIES	  := $(wildcard *.cxx)
CINTLIST  := $(addsuffix .h, $(basename $(GENIES)))
CINTLISTO := $(addsuffix .$(ObjSuf), $(basename $(CINTLIST)))
GENIES 	  := $(GENIES) _ROOT_DICT_Base.cc $(wildcard *.cpp)
GENIEO 	   = $(addsuffix .$(ObjSuf), $(basename $(GENIES)))
GENIESO   := $(GENIE_BASE_SHARED_LIB)

all:	rootcint lib

rootcint: _ROOT_DICT_Base.o

_ROOT_DICT_Base.cc: LinkDef.h $(CINTLISTO)
	$(RM) _ROOT_DICT_Base.*
	@echo "Generating ROOT dictionary ..."
	rootcint -f _ROOT_DICT_Base.cc -c	\
		-I$(ROOTSYS)/include -I$(GENIE)/src $(CINTLIST) LinkDef.h

lib: FORCE
	@echo "Generating shared library ..."
	@$(MAKE) -f $(MAKEFILE) $(GENIESO)

$(GENIESO): $(GENIEO)
	$(SOCMD) $(SOFLAGS) $(SOMINF)				\
		$(GENIEO) $(OutPutOpt) $(GENIESO)

clean: FORCE
	$(RM) *.o *~ core a.out $(GENIESO)		        \
		_ROOT_DICT_Base.* Makefile.temp

FORCE:

