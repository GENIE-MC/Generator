#
# Makefile for GENIE Apps
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)

# Mac OS X -bind_at_load means we must list the low level libraries
# first, and each stage must fully resolved before the next library
LIBRARIES  := $(LIBRARIES) $(GENIE_LIBS)
#LIBRARIES  := $(LIBRARIES) $(CERN_LIBRARIES) $(GENIE_LIBS)

TGT_BASE =  gevgen          \
            gevgen_hadron   \
            gevdump         \
            gevpick         \
            gevscan         \
            gevcomp         \
            gxscomp         \
            gmkspl          \
            gspladd         \
            gspl2root       \
            gntpc           \
            gpdfcomp        \
            gsfcomp

ifeq ($(strip $(GOPT_ENABLE_FNAL)),YES)
TGT_BASE += gevgen_fnal
endif
ifeq ($(strip $(GOPT_ENABLE_T2K)),YES)
TGT_BASE += gevgen_t2k
endif
ifeq ($(strip $(GOPT_ENABLE_ATMO)),YES)
TGT_BASE += gevgen_atmo
TGT_BASE += gevgen_upmu
endif
ifeq ($(strip $(GOPT_ENABLE_NUCLEON_DECAY)),YES)
TGT_BASE += gevgen_ndcy
endif
ifeq ($(strip $(GOPT_ENABLE_NNBAR_OSCILLATION)),YES)
TGT_BASE += gevgen_nnbarosc
# ------------------------------------------
# This should dissapear from here (CA)
#
# Require C++11 for NNBarOsc module.
#GCC_GTEQ_470 := $(shell expr `gcc -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40700)
#ifeq "$(GCC_GTEQ_470)" "1"
#        CXXFLAGS += -std=c++11
#else
#        CXXFLAGS += -std=c++0x
#endif
# ------------------------------------------
endif
ifeq ($(strip $(GOPT_ENABLE_BOOSTED_DARK_MATTER)),YES)
TGT_BASE += gevgen_dm
TGT_BASE += gmkspl_dm
TGT_BASE += gevgen_lardm
endif
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
TGT_BASE += gmxpl
endif
ifeq ($(strip $(GOPT_ENABLE_MASTERCLASS)),YES)
TGT_BASE += gmstcl
endif

TGT = $(addprefix $(GENIE_BIN_PATH)/,$(TGT_BASE))

# this is a nasty-looking construction which runs ldd on an existing library
# to figure out what other libraries it depends on.
# (sort of like g++ -M)
find_libs = ${shell ldd $(GENIE_BIN_PATH)/$(1) 2>&1 | sed -n 's/.*=> \(.*\) (.*/\1/pg' | tr '\n' ' '}

all: $(TGT)

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -MMD -MP -c $(CPP_INCLUDES) $< -o $@


# a simple/generic gevgen event generation app
#
$(GENIE_BIN_PATH)/gevgen: gEvGen.o $(call find_libs,gevgen)
	@echo "** Building gevgen"
	$(LD) $(LDFLAGS) gEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen

# neutrino event generation app for Fermilab experiments (DUNE, or experiments in NuMI and Booster beamlines)
#
$(GENIE_BIN_PATH)/gevgen_fnal: gFNALExptEvGen.o $(call find_libs,gevgen_fnal)
	@echo "** Building gevgen_fnal"
	$(LD) $(LDFLAGS) gFNALExptEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_fnal
#	@echo "For historical purposes, add a NuMI event generation app sym link"
#	ln -sf $(GENIE_BIN_PATH)/gevgen_fnal $(GENIE_BIN_PATH)/gevgen_numi

# T2K neutrino event generation app
#
$(GENIE_BIN_PATH)/gevgen_t2k: gT2KEvGen.o $(call find_libs,gevgen_t2k)
	@echo "** Building gevgen_t2k"
	$(LD) $(LDFLAGS) gT2KEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_t2k
#	@echo "For historical purposes, add symbolic link to the T2K event generation app"
#	ln -sf $(GENIE_BIN_PATH)/gevgen_t2k $(GENIE_BIN_PATH)/gT2Kevgen

# atmospheric neutrino event generation app
#
$(GENIE_BIN_PATH)/gevgen_atmo: gAtmoEvGen.o $(call find_libs,gevgen_atmo)
	@echo "** Building gevgen_atmo"
	$(LD) $(LDFLAGS) gAtmoEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_atmo

# (atmospheric) neutrino-incuded upgoing-muon flux generation app
#
$(GENIE_BIN_PATH)/gevgen_upmu: gUpMuFluxGen.o $(call find_libs,gevgen_upmu)
	@echo "** Building gevgen_upmu"
	$(LD) $(LDFLAGS) gUpMuFluxGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_upmu

# hadron + nucleus event generation app
#
$(GENIE_BIN_PATH)/gevgen_hadron: gEvGenHadronNucleus.o $(call find_libs,gevgen_hadron)
	@echo "** Building gevgen_hadron"
	$(LD) $(LDFLAGS) gEvGenHadronNucleus.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_hadron

# nucleon decay event generation app
#
$(GENIE_BIN_PATH)/gevgen_ndcy: gNucleonDecayEvGen.o $(call find_libs,gevgen_ndcy)
	@echo "** Building gevgen_ndcy"
	$(LD) $(LDFLAGS) gNucleonDecayEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_ndcy

# neutron-antineutron event generation app
#
$(GENIE_BIN_PATH)/gevgen_nnbarosc: gNNBarOscEvGen.o $(call find_libs,gevgen_nnbarosc)
	@echo "** Building gevgen_nnbarosc"
	$(LD) $(LDFLAGS) gNNBarOscEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_nnbarosc

# a simple/generic gevgen_dm event generation app for dark matter scattering
#
$(GENIE_BIN_PATH)/gevgen_dm: gEvGenDM.o $(call find_libs,gevgen_dm)
	@echo "** Building gevgen_dm"
	$(LD) $(LDFLAGS) gEvGenDM.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_dm

# gevdump event print-out utility
#
$(GENIE_BIN_PATH)/gevdump: gEvDump.o $(call find_libs,gevdump)
	@echo "** Building gevdump"
	$(LD) $(LDFLAGS) gEvDump.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevdump

# gevgen event generation utility for use with simple fluxes
#
$(GENIE_BIN_PATH)/gevgen_lardm: gEvGenLArDM.o $(call find_libs,gevgen_lardm)
	@echo "** Building gevgen_lardm"
	$(LD) $(LDFLAGS) gEvGenLArDM.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_lardm

# event topology cherry-picking program
#
$(GENIE_BIN_PATH)/gevpick: gEvPick.o $(call find_libs,gevpick)
	@echo "** Building gevpick"
	$(LD) $(LDFLAGS) gEvPick.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevpick

# utility performing sanity checks on event samples
#
$(GENIE_BIN_PATH)/gevscan: gEvScan.o $(call find_libs,gevscan)
	@echo "** Building gevscan"
	$(LD) $(LDFLAGS) gEvScan.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevscan

# utility performing comparisons between two event samples
#
$(GENIE_BIN_PATH)/gevcomp: gEvComp.o $(call find_libs,gevcomp)
	@echo "** Building gevcomp"
	$(LD) $(LDFLAGS) gEvComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevcomp

# utility performing comparisons between two sets of pre-computed x-section splines
#
$(GENIE_BIN_PATH)/gxscomp: gXSecComp.o $(call find_libs,gxscomp)
	@echo "** Building gxscomp"
	$(LD) $(LDFLAGS) gXSecComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gxscomp

# x-section spline building utility
#
$(GENIE_BIN_PATH)/gmkspl: gMakeSplines.o $(call find_libs,gmkspl)
	@echo "** Building gmkspl"
	$(LD) $(LDFLAGS) gMakeSplines.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmkspl

# x-section spline building utility for dark matter scattering
#
$(GENIE_BIN_PATH)/gmkspl_dm: gMakeSplinesDM.o $(call find_libs,gmkspl_dm)
	@echo "** Building gmkspl_dm"
	$(LD) $(LDFLAGS) gMakeSplinesDM.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmkspl_dm

# utility merging x-section spline xml files
#
$(GENIE_BIN_PATH)/gspladd: gSplineAdd.o $(call find_libs,gspladd)
	@echo "** Building gspladd"
	$(LD) $(LDFLAGS) gSplineAdd.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gspladd

# utility for converting XML splines into ROOT format
#
$(GENIE_BIN_PATH)/gspl2root: gSplineXml2Root.o $(call find_libs,gspl2root)
	@echo "** Building gspl2root"
	$(LD) $(LDFLAGS) gSplineXml2Root.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gspl2root

# utility computing maximum path lengths for a given root geometry
#
$(GENIE_BIN_PATH)/gmxpl: gMaxPathLengths.o $(call find_libs,gmxpl)
	@echo "** Building gmxpl"
	$(LD) $(LDFLAGS) gMaxPathLengths.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmxpl

# ntuple conversion utility
#
$(GENIE_BIN_PATH)/gntpc: gNtpConv.o $(call find_libs,gntpc)
	@echo "** Building gntpc"
	$(LD) $(LDFLAGS) gNtpConv.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gntpc

# Masterclass app
#
$(GENIE_BIN_PATH)/gmstcl: gMasterclass.o $(call find_libs,gmstcl)
	@echo "** Building gmstcl";
	$(LD) $(LDFLAGS) gMasterclass.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmstcl

# App to compare PDF
#
$(GENIE_BIN_PATH)/gpdfcomp: gPDFComp.o $(call find_libs,gpdfcomp)
	@echo "** Building gpdfcomp"
	$(LD) $(LDFLAGS) gPDFComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gpdfcomp

# App to compare structure function models (see GENIE/Comparisons for an app to compare with data)
#
$(GENIE_BIN_PATH)/gsfcomp: gSFComp.o $(call find_libs,gsfcomp)
	@echo "** Building gsfcomp"
	$(LD) $(LDFLAGS) gSFComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gsfcomp

# CLEANING-UP

purge: FORCE
	$(RM) *.o *~ *.d core

clean: FORCE
	$(RM) *.o *~ *.d core
	$(RM) $(TGT)

distclean: FORCE
	$(RM) $(patsubst $(GENIE_BIN_DIR),$(GENIE_BIN_INSTALLATION_PATH),$(TGT))

FORCE:

# DO NOT DELETE

# this should be last because it includes new targets...
-include $(wildcard *.d)

