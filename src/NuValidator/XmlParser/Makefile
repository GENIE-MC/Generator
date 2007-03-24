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
PROJECT_O   = $(addsuffix .$(ObjSuf), $(basename $(PROJECT_S)))

all: lib

lib: $(PROJECT_O)
	@echo "Generating shared library ..."
	$(SOCMD) $(SOFLAGS) $(SOMINF) \
		$(PROJECT_O) $(OutPutOpt) $(GENIE_LIB_PATH)/$(VLD_PARSER_LIBNAME).$(DllSuf)

purge: FORCE
	$(RM) *.o *~ core a.out 

clean: FORCE
	$(RM) *.o *~ core a.out $(GENIE_LIB_PATH)/$(VLD_PARSER_LIBNAME).*

distclean: FORCE
	$(RM) $(GENIE_INSTALLATION_PATH)/lib/$(VLD_PARSER_LIBNAME).*

FORCE:

# DO NOT DELETE
