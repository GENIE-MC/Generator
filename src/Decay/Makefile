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
GENIES 	  := $(GENIES) _ROOT_DICT_Decay.cc $(wildcard *.cpp)
GENIEO 	   = $(addsuffix .$(ObjSuf), $(basename $(GENIES)))
GENIEbS   := $(wildcard *.cxx)
GENIEbO    = $(addsuffix .$(ObjSuf), $(basename $(GENIEbS)))

all: rootcint lib link

rootcint: _ROOT_DICT_Decay.o

_ROOT_DICT_Decay.cc: LinkDef.h $(CINTLISTO)
	$(RM) _ROOT_DICT_Decay.*
	@echo "Generating ROOT dictionary ..."
	rootcint -f _ROOT_DICT_Decay.cc -c	\
		-I$(ROOTSYS)/include -I$(GENIE)/src $(CINTLIST) LinkDef.h

lib: $(GENIEO) 
	@echo "Generating shared library ..."
	$(SOCMD) $(SOFLAGS) $(SOMINF) \
		$(GENIEO) $(OutPutOpt) $(GENIE_LIB_PATH)/$(GENIE_DECAY_LIBNAME).$(DllSuf)

link: FORCE
ifeq ($(strip $(ARCH)),macosx)
	@echo "Creating symbolic link to shared library ..."
	ln -sf $(GENIE_DECAY_LIBNAME).$(DllSuf) $(GENIE_LIB_PATH)/$(GENIE_DECAY_LIBNAME).$(DllLinkSuf)
else
endif

purge: FORCE
	$(RM) $(GENIEbO) *~ core a.out Makefile.temp

clean: FORCE
	$(RM) *.o *~ core a.out $(GENIE_LIB_PATH)/$(GENIE_DECAY_LIBNAME).* \
		_ROOT_DICT_Decay.* Makefile.temp

FORCE:

