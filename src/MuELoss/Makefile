#
# Makefile for GENIE Neutrino Generator
#
# Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/Make.include

GENIES	  := $(wildcard *.cxx)
CINTLIST  := $(addsuffix .h, $(basename $(GENIES)))
CINTLISTO := $(addsuffix .$(ObjSuf), $(basename $(CINTLIST)))
GENIES 	  := $(GENIES) _ROOT_DICT_MuELoss.cc $(wildcard *.cpp)
GENIEO 	   = $(addsuffix .$(ObjSuf), $(basename $(GENIES)))
GENIEbS   := $(wildcard *.cxx)
GENIEbO    = $(addsuffix .$(ObjSuf), $(basename $(GENIEbS)))

all: rootcint lib link

rootcint: _ROOT_DICT_MuELoss.o

_ROOT_DICT_MuELoss.cc: LinkDef.h $(CINTLISTO)
	$(RM) _ROOT_DICT_MuELoss.*
	@echo "Generating ROOT dictionary ..."
	rootcint -f _ROOT_DICT_MuELoss.cc -c	\
		-I$(ROOTSYS)/include -I$(GENIE)/src $(CINTLIST) LinkDef.h

lib: $(GENIEO) 
	@echo "Generating shared library ..."
	$(SOCMD) $(SOFLAGS) $(SOMINF) \
		$(GENIEO) $(OutPutOpt) $(GENIE_LIB_PATH)/$(GENIE_MUELOSS_LIBNAME)$(GLIBVERSION).$(DllSuf)

link: FORCE
	@echo "Creating symbolic link to shared library ..."
ifeq ($(strip $(GOPT_ENABLE_DYLIBVERSION)),YES)
	ln -sf $(GENIE_MUELOSS_LIBNAME)$(GLIBVERSION).$(DllSuf) $(GENIE_LIB_PATH)/$(GENIE_MUELOSS_LIBNAME).$(DllSuf)
endif
ifeq ($(strip $(ARCH)),macosx)
	ln -sf $(GENIE_MUELOSS_LIBNAME)$(GLIBVERSION).$(DllSuf) $(GENIE_LIB_PATH)/$(GENIE_MUELOSS_LIBNAME).$(DllLinkSuf)
else
endif

purge: FORCE
	$(RM) $(GENIEbO) *~ core a.out Makefile.temp

clean: FORCE
	$(RM) *.o *~ core a.out $(GENIE_LIB_PATH)/$(GENIE_MUELOSS_LIBNAME)* \
		_ROOT_DICT_MuELoss.* Makefile.temp

distclean: FORCE
	$(RM) $(GENIE_LIB_INSTALLATION_PATH)/$(GENIE_MUELOSS_LIBNAME)*

FORCE:

