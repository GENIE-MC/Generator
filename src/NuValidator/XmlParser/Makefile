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

lib: FORCE
	@echo $(XML_INC_DIR)
	@echo "Generating shared library ..."
	@$(MAKE) -f $(MAKEFILE) $(VLD_PARSER_SHARED_LIB)

$(VLD_PARSER_SHARED_LIB): $(PROJECT_O)
	$(SOCMD) $(SOFLAGS) $(SOMINF)				\
		$(PROJECT_O) $(OutPutOpt) $(VLD_PARSER_SHARED_LIB)

clean: FORCE
	$(RM) *.o *~ core a.out $(VLD_PARSER_SHARED_LIB)

FORCE:

# DO NOT DELETE
