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

all: neugen3_dummy_lib

neugen3_dummy_lib: FORCE
	@echo "Generating NeuGEN dummy static library ..."
	$(FORT) $(FORTOPTS) -o $(GDUMMY_NEUGEN_STATIC_LIB) NeuGenDUMMY.f

clean: FORCE
	$(RM) *.o *~ core a.out $(GDUMMY_NEUGEN_STATIC_LIB) 

FORCE:

# DO NOT DELETE
