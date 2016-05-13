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

LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT =    gevgen		  	\
	 gevgen_fnal            \
	 gevgen_t2k             \
	 gevgen_atmo            \
	 gevgen_hadron          \
         gevgen_ndcy            \
	 gevgen_upmu            \
	 gevdump          	\
	 gevpick	        \
	 gevscan	        \
	 gevcomp	        \
	 gxscomp	        \
	 gmkspl  	  	\
	 gspladd   	  	\
	 gspl2root	  	\
	 gmxpl		  	\
	 grwght1p		\
	 grwghtnp		\
	 gntpc	  		\
	 gmstcl			\
	 gpdfcomp		\
	 gsfcomp

all: $(TGT)

# a simple/generic gevgen event generation app
#
gevgen: FORCE
	@echo "** Building gevgen"
	$(CXX) $(CXXFLAGS) -c gEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen 

# neutrino event generation app for Fermilab experiments (DUNE, or experiments in NuMI and Booster beamlines)
#
gevgen_fnal: FORCE
	@echo "** Building gevgen_fnal"
ifeq ($(strip $(GOPT_ENABLE_FNAL)),YES)
	$(CXX) $(CXXFLAGS) -c gFNALExptEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gFNALExptEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_fnal
	@echo "For historical purposes, add a NuMI event generation app sym link"
	ln -sf $(GENIE_BIN_PATH)/gevgen_fnal $(GENIE_BIN_PATH)/gevgen_numi
else
	@echo "FNAL support not enabled!"
endif

# T2K neutrino event generation app 
#
gevgen_t2k: FORCE
	@echo "** Building gevgen_t2k"
ifeq ($(strip $(GOPT_ENABLE_T2K)),YES)
	$(CXX) $(CXXFLAGS) -c gT2KEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gT2KEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_t2k
	@echo "For historical purposes, add symbolic link to the T2K event generation app"
	ln -sf $(GENIE_BIN_PATH)/gevgen_t2k $(GENIE_BIN_PATH)/gT2Kevgen
else
	@echo "T2K support not enabled!"
endif

# atmospheric neutrino event generation app
#
gevgen_atmo: FORCE
	@echo "** Building gevgen_atmo"
ifeq ($(strip $(GOPT_ENABLE_ATMO)),YES)
	$(CXX) $(CXXFLAGS) -c gAtmoEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gAtmoEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_atmo
else
	@echo "Atmospheric neutrino support not enabled!"
endif

# (atmospheric) neutrino-incuded upgoing-muon flux generation app
#
gevgen_upmu: FORCE
	@echo "** Building gevgen_upmu"
ifeq ($(strip $(GOPT_ENABLE_ATMO)),YES)
	$(CXX) $(CXXFLAGS) -c gUpMuFluxGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gUpMuFluxGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_upmu
else
	@echo "Atmospheric neutrino support not enabled!"
endif

# hadron + nucleus event generation app
#
gevgen_hadron: FORCE
	@echo "** Building gevgen_hadron"
	$(CXX) $(CXXFLAGS) -c gEvGenHadronNucleus.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvGenHadronNucleus.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_hadron

# nucleon decay event generation app
#
gevgen_ndcy: FORCE
	@echo "** Building gevgen_ndcy"
ifeq ($(strip $(GOPT_ENABLE_NUCLEON_DECAY)),YES)
	$(CXX) $(CXXFLAGS) -c gNucleonDecayEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gNucleonDecayEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_ndcy
else
	@echo "Nucleon decay not enabled!"
endif

# gevdump event print-out utility
#
gevdump: FORCE
	@echo "** Building gevdump"
	$(CXX) $(CXXFLAGS) -c gEvDump.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvDump.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevdump

# event topology cherry-picking program
#
gevpick: FORCE
	@echo "** Building gevpick"
	$(CXX) $(CXXFLAGS) -c gEvPick.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvPick.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevpick

# utility performing sanity checks on event samples
#
gevscan: FORCE
	@echo "** Building gevscan"
	$(CXX) $(CXXFLAGS) -c gEvScan.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvScan.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevscan

# utility performing comparisons between two event samples
#
gevcomp: FORCE
	@echo "** Building gevcomp"
	$(CXX) $(CXXFLAGS) -c gEvComp.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevcomp

# utility performing comparisons between two sets of pre-computed x-section splines
#
gxscomp: FORCE
	@echo "** Building gxscomp"
	$(CXX) $(CXXFLAGS) -c gXSecComp.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gXSecComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gxscomp

# x-section spline building utility
#
gmkspl: FORCE
	@echo "** Building gmkspl"
	$(CXX) $(CXXFLAGS) -c gMakeSplines.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gMakeSplines.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmkspl 

# utility merging x-section spline xml files
#
gspladd: FORCE
	@echo "** Building gspladd"
	$(CXX) $(CXXFLAGS) -c gSplineAdd.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gSplineAdd.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gspladd 

# utility for converting XML splines into ROOT format
#
gspl2root: FORCE
	@echo "** Building gspl2root"
	$(CXX) $(CXXFLAGS) -c gSplineXml2Root.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gSplineXml2Root.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gspl2root

# utility computing maximum path lengths for a given root geometry
#
gmxpl: FORCE
	@echo "** Building gmxpl"
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c gMaxPathLengths.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gMaxPathLengths.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmxpl
else
	@echo "You need to enable the flux drivers to build the gmxpl utility"
endif

# simple utility for calculating weights in response to tweaking a single systematic param
#
grwght1p: FORCE
	@echo "** Building grwght1p"
ifeq ($(strip $(GOPT_ENABLE_RWGHT)),YES)
	$(CXX) $(CXXFLAGS) -c gRwght1Param.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gRwght1Param.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/grwght1p
else
	@echo "Event reweighting not enabled!"
endif

# simple utility for calculating weights in response to tweaking a multiple (possibly correlated) systematic param
#
grwghtnp: FORCE
	@echo "** Building grwghtnp"
ifeq ($(strip $(GOPT_ENABLE_RWGHT)),YES)
	$(CXX) $(CXXFLAGS) -c gRwghtNCorrelatedParams.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gRwghtNCorrelatedParams.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/grwghtnp
else
	@echo "Event reweighting not enabled!"
endif


# ntuple conversion utility
#
gntpc: FORCE
	@echo "** Building gntpc"
	$(CXX) $(CXXFLAGS) -c gNtpConv.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gNtpConv.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gntpc 

# Masterclass app
#
gmstcl: FORCE
	@echo "** Building gmstcl";
ifeq ($(strip $(GOPT_ENABLE_MASTERCLASS)),YES)
	$(CXX) $(CXXFLAGS) -c gMasterclass.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gMasterclass.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmstcl 
else
	@echo "Masterclass not enabled!"
endif

# App to compare PDF 
#
gpdfcomp: FORCE
	@echo "** Building gpdfcomp"
	$(CXX) $(CXXFLAGS) -c gPDFComp.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gPDFComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gpdfcomp

# App to compare structure function models (see GENIE/Comparisons for an app to compare with data)
#
gsfcomp: FORCE
	@echo "** Building gsfcomp"
	$(CXX) $(CXXFLAGS) -c gSFComp.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gSFComp.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gsfcomp

# CLEANING-UP

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) $(GENIE_BIN_PATH)/gevgen		
	$(RM) $(GENIE_BIN_PATH)/gevgen_fnal		
	$(RM) $(GENIE_BIN_PATH)/gevgen_t2k		
	$(RM) $(GENIE_BIN_PATH)/gevgen_atmo		
	$(RM) $(GENIE_BIN_PATH)/gevgen_upmu		
	$(RM) $(GENIE_BIN_PATH)/gevgen_hadron
	$(RM) $(GENIE_BIN_PATH)/gevgen_ndcy
	$(RM) $(GENIE_BIN_PATH)/gevdump		
	$(RM) $(GENIE_BIN_PATH)/gevpick
	$(RM) $(GENIE_BIN_PATH)/gevscan
	$(RM) $(GENIE_BIN_PATH)/gevcomp
	$(RM) $(GENIE_BIN_PATH)/gxscomp
	$(RM) $(GENIE_BIN_PATH)/gmkspl 	
	$(RM) $(GENIE_BIN_PATH)/gspladd 	
	$(RM) $(GENIE_BIN_PATH)/gspl2root		
	$(RM) $(GENIE_BIN_PATH)/gmxpl		
	$(RM) $(GENIE_BIN_PATH)/grwght1p
	$(RM) $(GENIE_BIN_PATH)/grwghtnp
	$(RM) $(GENIE_BIN_PATH)/gntpc		
	$(RM) $(GENIE_BIN_PATH)/gmstcl		
	$(RM) $(GENIE_BIN_PATH)/gpdfcomp

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_fnal		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_t2k		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_atmo		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_upmu		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_hadron
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_ndcy
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevdump		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevpick
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevscan
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevcomp
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gxscomp
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gmkspl 	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gspladd 	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gspl2root		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gmxpl		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/grwght1p
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/grwghtnp
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gntpc		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gmstcl		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gpdfcomp

FORCE:

# DO NOT DELETE
