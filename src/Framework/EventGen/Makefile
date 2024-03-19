#
# Makefile for GENIE Neutrino Generator
#
# Costas Andreopoulos <constantinos.andreopoulos \at cern.ch>
#

SHELL    = /bin/sh
NAME     = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/make/Make.include

PACKAGE           = Framework/EventGen
PACKAGE_ABBREV    = FwEG
DICTIONARY        = _ROOT_DICT_$(PACKAGE_ABBREV)
LIBNAME           = libG$(PACKAGE_ABBREV)
EXTRA_EXT_LIBS    =

all     : rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(GENIE)/src/make/Make.std-package-targets

# Only need these if we're building with HepMC3 support
ifeq ($(strip $(GOPT_ENABLE_HEPMC3_INTERFACE)),YES)
  # extra flags, include paths, and libraries to link to
  CXXFLAGS               += $(HEPMC3_FLAGS)
  CPP_INCLUDES           += $(HEPMC3_INCLUDES)
  ROOT_DICT_GEN_INCLUDES += $(HEPMC3_INCLUDES)
  EXTRA_EXT_LIBS         += $(HEPMC3_LIBRARIES)
else
  $(info $(PACKAGE) not built against HepMC3)
endif

FORCE:

