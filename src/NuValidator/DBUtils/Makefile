#########################################################
# Makefile for Neutrino Generator Validator
#
# Costas Andreopoulos, 2003-08-20
#########################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/Make.include

CSOURCES   = $(wildcard *.c) $(wildcard *.cxx)
DEPSOURCES = $(CSOURCES) $(wildcard *.C)

PROJECT_S  := $(wildcard *.$(SrcSuf))
CINTLIST_H := $(addsuffix .h, $(basename $(PROJECT_S)))
CINTLIST_O := $(addsuffix .$(ObjSuf), $(basename $(CINTLIST_H)))
PROJECT_S  := $(PROJECT_S) _ROOT_DICT_DBUtils.cc
PROJECT_O   = $(addsuffix .$(ObjSuf), $(basename $(PROJECT_S)))
PROJECT_bS := $(wildcard *.$(SrcSuf))
PROJECT_bO  = $(addsuffix .$(ObjSuf), $(basename $(PROJECT_bS)))

all: rootcint lib link

rootcint: _ROOT_DICT_DBUtils.o

_ROOT_DICT_DBUtils.cc: LinkDef.h $(CINTLIST_O)
	$(RM) _ROOT_DICT_DBUtils.*
	@echo "Generating CINT dictionary ..."
	rootcint -f _ROOT_DICT_DBUtils.cc -c \
	    $(ROOT_INCLUDES) $(NUVLD_INCLUDES) $(CINTLIST_H) LinkDef.h

lib: $(PROJECT_O) FORCE
	@echo "Generating shared library ..."
	$(SOCMD) $(SOFLAGS) $(SOMINF) $(PROJECT_O) $(OutPutOpt) $(GENIE_LIB_PATH)/$(VLD_DBUTIL_LIBNAME).$(DllSuf)

link: FORCE
ifeq ($(strip $(ARCH)),macosx)
	@echo "Creating symbolic link to shared library ..."
	ln -sf $(VLD_DBUTIL_LIBNAME).$(DllSuf) $(GENIE_LIB_PATH)/$(VLD_DBUTIL_LIBNAME).$(DllLinkSuf)
else
endif

purge: FORCE
	$(RM) *.o *~ core a.out _ROOT_DICT_DBUtils.*

clean: FORCE
	$(RM) *.o *~ core a.out $(GENIE_LIB_PATH)/$(VLD_DBUTIL_LIBNAME).* _ROOT_DICT_DBUtils.* 

FORCE:

# DO NOT DELETE
