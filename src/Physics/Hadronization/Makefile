#
# Makefile for GENIE Neutrino Generator
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

PACKAGE           = Physics/Hadronization
PACKAGE_ABBREV    = PhHadnz
DICTIONARY        = _ROOT_DICT_$(PACKAGE_ABBREV)
LIBNAME           = libG$(PACKAGE_ABBREV)
EXTRA_EXT_LIBS    =

all     : rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(GENIE)/src/make/Make.std-package-targets

ifeq ($(strip $(GOPT_ENABLE_PYTHIA8)),YES)
  CPP_INCLUDES           += $(PYTHIA8_INCLUDES)
  ROOT_DICT_GEN_INCLUDES += $(PYTHIA8_INCLUDES) 
  EXTRA_EXT_LIBS         += $(PYTHIA8_LIBRARIES)
else
  $(info $(PACKAGE) not build against Pythia8 )
endif

FORCE:
