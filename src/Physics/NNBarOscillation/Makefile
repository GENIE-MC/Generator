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

# ------------------------------------------
# This needs to disappear from here... (CA)
#
# Require C++11 for NeutronOsc module.
GCC_GTEQ_470 := $(shell expr `gcc -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40700)
ifeq "$(GCC_GTEQ_470)" "1"
	CXXFLAGS += -std=c++11
else
	CXXFLAGS += -std=c++0x
endif
# ------------------------------------------

PACKAGE           = Physics/NNBarOscillation
PACKAGE_ABBREV    = PhNNBarOsc
DICTIONARY        = _ROOT_DICT_$(PACKAGE_ABBREV)
LIBNAME           = libG$(PACKAGE_ABBREV)

all     : rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(GENIE)/src/make/Make.std-package-targets

FORCE:
