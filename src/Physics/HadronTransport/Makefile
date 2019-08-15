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

PACKAGE           = Physics/HadronTransport
PACKAGE_ABBREV    = PhHadTransp
DICTIONARY        = _ROOT_DICT_$(PACKAGE_ABBREV)
LIBNAME           = libG$(PACKAGE_ABBREV)
EXTRA_EXT_LIBS    =

all     : rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(GENIE)/src/make/Make.std-package-targets

# Only need these here ... not in every package

ifeq ($(strip $(GOPT_ENABLE_INCL)),YES)
  # extra flags, include paths, and libraries to link to
  CXXFLAGS               += $(INCL_FLAGS)
  CPP_INCLUDES           += $(INCL_INCLUDES)
  ROOT_DICT_GEN_INCLUDES += $(INCL_INCLUDES) ${INCL_FLAGS}
  EXTRA_EXT_LIBS         += $(INCL_LIBRARIES)
else
  $(info $(PACKAGE) not build against INCL++)
  # remove INCL++ related files from things to build
### ! can't do it this way because of the way Make.std-package-targets works
### ! has already used PACKAGE_SOURCE/OBJECTS to define how library is built
##  # VAR := $(filter-out A B C,$(VAR))
##  INCL_SOURCES = HINCLCascade.cxx INCLCascade.cxx INCLInterfacetoGENIE.cxx
##  PACKAGE_SOURCES := $(filter-out $(INCL_SOURCES), $(PACKAGE_SOURCES))
##  PACKAGE_OBJECTS  = $(addsuffix .$(ObjSuf), $(basename $(PACKAGE_SOURCES)))
endif

FORCE: