############################################################
# Makefile for dummy neugen3 library
# Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
############################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/make/Make.include

all: neugen3_dummy_lib

neugen3_dummy_lib: FORCE
	@echo "Generating NeuGEN dummy static library ..."
	$(FORT) $(FORTOPTS) -o $(GENIE_LIB_PATH)/$(GDUMMY_NEUGEN_LIBNAME).$(StaticLibSuf) NeuGenDUMMY.f 

clean: FORCE
	$(RM) *.o *~ core a.out $(GENIE_LIB_PATH)/$(GDUMMY_NEUGEN_LIBNAME).* 

FORCE:

# DO NOT DELETE
