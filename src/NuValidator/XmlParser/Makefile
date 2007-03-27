#
# Makefile for Neutrino Generator Validator
#
# Author: Costas Andreopoulos
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/Make.include

CSOURCES   = $(wildcard *.c) $(wildcard *.cxx)
DEPSOURCES = $(CSOURCES) $(wildcard *.C)
PROJECT_S  := $(wildcard *.$(SrcSuf))
PROJECT_O   = $(addsuffix .$(ObjSuf), $(basename $(PROJECT_S)))

all: lib link

lib: $(PROJECT_O)
	@echo "Generating shared library ..."
	$(SOCMD) $(SOFLAGS) $(SOMINF) \
		$(PROJECT_O) $(OutPutOpt) $(GENIE_LIB_PATH)/$(VLD_PARSER_LIBNAME)$(GLIBVERSION).$(DllSuf)

link: FORCE
	@echo "Creating symbolic link to shared library ..."
ifeq ($(strip $(GOPT_ENABLE_DYLIBVERSION)),YES)
	ln -sf $(VLD_PARSER_LIBNAME)$(GLIBVERSION).$(DllSuf) $(GENIE_LIB_PATH)/$(VLD_PARSER_LIBNAME).$(DllSuf)
endif
ifeq ($(strip $(ARCH)),macosx)
	ln -sf $(VLD_PARSER_LIBNAME)$(GLIBVERSION).$(DllSuf) $(GENIE_LIB_PATH)/$(VLD_PARSER_LIBNAME).$(DllLinkSuf)
else
endif

purge: FORCE
	$(RM) *.o *~ core a.out 

clean: FORCE
	$(RM) *.o *~ core a.out $(GENIE_LIB_PATH)/$(VLD_PARSER_LIBNAME)*

distclean: FORCE
	$(RM) $(GENIE_LIB_INSTALLATION_PATH)/$(VLD_PARSER_LIBNAME)*

FORCE:

# DO NOT DELETE
