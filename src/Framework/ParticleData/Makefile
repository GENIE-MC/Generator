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

# Include INCL in here; INCL model include hyper-nuclei
# So, it is need to put hyper-nuclei into event record to
# avoid crashing.
ifeq ($(strip $(GOPT_ENABLE_INCL)),YES)
  # extra flags, include paths, and libraries to link to
  CXXFLAGS               += $(INCL_DFLAGS)
  CPP_INCLUDES           += $(INCL_INCLUDES)
  ROOT_DICT_GEN_INCLUDES += $(INCL_INCLUDES) $(INCL_DFLAGS)
  EXTRA_EXT_LIBS         += $(INCL_LIBRARIES)
else
  $(info $(PACKAGE) not built against INCL++)
endif

PACKAGE           = Framework/ParticleData
PACKAGE_ABBREV    = FwParDat
DICTIONARY        = _ROOT_DICT_$(PACKAGE_ABBREV)
LIBNAME           = libG$(PACKAGE_ABBREV)
EXTRA_EXT_LIBS    = -L$(GENIE_LIB_PATH) -lGFwMsg

all     : rootcint lib lib-link
install : install-inc install-lib

# Include standard package makefile targets
#
include $(GENIE)/src/make/Make.std-package-targets

FORCE:
