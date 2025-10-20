#
# Makefile for GENIE Neutrino Generator
#
# Costas Andreopoulos <c.andreopoulos \at cern.ch>
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
EXTRA_EXT_LIBS    = -L$(GENIE_LIB_PATH) -lGPhNuclSt

all     : rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(GENIE)/src/make/Make.std-package-targets

# Only need these here ... not in every package

ifeq ($(strip $(GOPT_ENABLE_INCL)),YES)
  # extra flags, include paths, and libraries to link to
  CXXFLAGS               += $(INCL_DFLAGS)
  CPP_INCLUDES           += $(INCL_INCLUDES)
  ROOT_DICT_GEN_INCLUDES += $(INCL_INCLUDES) $(INCL_DFLAGS)
  EXTRA_EXT_LIBS         += $(INCL_LIBRARIES)
else
  $(info $(PACKAGE) not built against INCL++)
endif

ifeq ($(strip $(GOPT_ENABLE_GEANT4_INTERFACE)),YES)
  # extra flags, include paths, and libraries to link to
  CXXFLAGS               += $(GEANT4_FLAGS)
  CPP_INCLUDES           += $(GEANT4_INCLUDES)
  ROOT_DICT_GEN_INCLUDES += $(GEANT4_INCLUDES)
  EXTRA_EXT_LIBS         += $(GEANT4_LIBRARIES)
else
  $(info $(PACKAGE) not built against Geant4)
endif

FORCE:
