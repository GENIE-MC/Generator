#########################################################
# Makefile for the Neutrino Generator Validator
#
# Costas Andreopoulos
#########################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/Make.include

#CSOURCES   = $(wildcard *.c) $(wildcard *.cxx)
#DEPSOURCES = $(CSOURCES) $(wildcard *.C)

PROJECT_S  := $(wildcard *.$(SrcSuf))
CINTLIST_H := $(addsuffix .h, $(basename $(PROJECT_S)))
CINTLIST_O := $(addsuffix .$(ObjSuf), $(basename $(CINTLIST_H)))
PROJECT_S  := $(PROJECT_S) _ROOT_DICT_NuVldGUI.cc
PROJECT_O   = $(addsuffix .$(ObjSuf), $(basename $(PROJECT_S)))

all:	rootcint lib

rootcint: _ROOT_DICT_NuVldGUI.o

_ROOT_DICT_NuVldGUI.cc: LinkDef.h $(CINTLIST_O)
	$(RM) _ROOT_DICT_NuVldGUI.*
	@echo "Generating CINT dictionary ..."
	rootcint -f _ROOT_DICT_NuVldGUI.cc -c \
	-I$(ROOTSYS)/include $(NUVLD_INC) $(GENIE_INC) $(CINTLIST_H) LinkDef.h

lib: FORCE
	@echo "Generating shared library ..."
	@$(MAKE) -f $(MAKEFILE) $(VLD_GUI_WIDGETS_SHARED_LIB)

$(VLD_GUI_WIDGETS_SHARED_LIB): $(PROJECT_O)
	$(SOCMD) $(SOFLAGS) $(SOMINF) $(PROJECT_O) $(OutPutOpt) $(VLD_GUI_WIDGETS_SHARED_LIB)

clean: FORCE
	$(RM) *.o *~ core a.out $(VLD_GUI_WIDGETS_SHARED_LIB) _ROOT_DICT_NuVldGUI.* 

FORCE:

# DO NOT DELETE
