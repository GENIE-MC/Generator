#
# Makefile for GENIE test programs
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# May 02, 2004
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
#LIBRARIES  := $(LHAPDF_LIBRARIES) $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)
LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)


TGT =	gtestAlgorithms 	 \
	gtestBLI2DUnifGrid       \
	gtestCmdLnArg		 \
 	gtestConfigPool		 \
 	gtestDecay		 \
 	gtestDISSF		 \
 	gtestElFormFactors	 \
 	gtestEventLoop 		 \
 	gtestFluxAstro 		 \
 	gtestFluxAtmo 		 \
 	gtestFluxSimple 	 \
	gtestFGPauliBlockSuppr   \
        gtestGiBUUData           \
	gtestHadronization	 \
	gtestINukeHadroData      \
	gtestMessenger		 \
	gtestNaturalIsotopes	 \
	gtestPDFLIB		 \
	gtestPREM		 \
	gtestROOTGeometry	 \
	gtestFermiP		 \
	gtestRewght		 \
	gtestRegistry		 \
	gtestInteraction	 \
	gtestResonances		 \
	gtestRwMarginalization	 \
	gtestKPhaseSpace	

all: $(TGT)

gtestAlgorithms: FORCE
	$(CXX) $(CXXFLAGS) -c gtestAlgorithms.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestAlgorithms.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestAlgorithms

gtestBLI2DUnifGrid: FORCE
	$(CXX) $(CXXFLAGS) -c gtestBLI2DUnifGrid.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestBLI2DUnifGrid.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestBLI2DUnifGrid

gtestCmdLnArg: FORCE
	$(CXX) $(CXXFLAGS) -c gtestCmdLnArg.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestCmdLnArg.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestCmdLnArg

gtestConfigPool: FORCE
	$(CXX) $(CXXFLAGS) -c gtestConfigPool.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestConfigPool.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestConfigPool

gtestDecay: FORCE
	$(CXX) $(CXXFLAGS) -c gtestDecay.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestDecay.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestDecay

gtestDISSF: FORCE
	$(CXX) $(CXXFLAGS) -c gtestDISSF.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestDISSF.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestDISSF

gtestElFormFactors: FORCE
	$(CXX) $(CXXFLAGS) -c gtestElFormFactors.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestElFormFactors.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestElFormFactors

gtestEventLoop: FORCE
	$(CXX) $(CXXFLAGS) -c gtestEventLoop.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestEventLoop.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestEventLoop

gtestFluxAstro: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c gtestFluxAstro.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestFluxAstro.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFluxAstro
else
	@echo "You need to enable the flux drivers to build the gtestFluxAstro program"
endif

gtestFluxAtmo: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c gtestFluxAtmo.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestFluxAtmo.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFluxAtmo
else
	@echo "You need to enable the flux drivers to build the gtestFluxAtmo program"
endif

gtestFluxSimple: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c gtestFluxSimple.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestFluxSimple.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFluxSimple
else
	@echo "You need to enable the flux drivers to build the gtestFluxSimple program"
endif

gtestFGPauliBlockSuppr:
	$(CXX) $(CXXFLAGS) -c gtestFGPauliBlockSuppr.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestFGPauliBlockSuppr.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFGPauliBlockSuppr

gtestGiBUUData: FORCE
	$(CXX) $(CXXFLAGS) -c gtestGiBUUData.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestGiBUUData.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestGiBUUData

gtestHadronization: FORCE
	$(CXX) $(CXXFLAGS) -c gtestHadronization.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestHadronization.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestHadronization

gtestINukeHadroData: FORCE
	$(CXX) $(CXXFLAGS) -c gtestINukeHadroData.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestINukeHadroData.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestINukeHadroData

gtestXSec: FORCE
	$(CXX) $(CXXFLAGS) -c gtestXSec.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestXSec.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestXSec

gtestPDFLIB: FORCE
	$(CXX) $(CXXFLAGS) -c gtestPDFLIB.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestPDFLIB.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestPDFLIB

gtestPREM: FORCE
	$(CXX) $(CXXFLAGS) -c gtestPREM.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestPREM.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestPREM

gtestFermiP: FORCE
	$(CXX) $(CXXFLAGS) -c gtestFermiP.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestFermiP.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFermiP

gtestMessenger: FORCE
	$(CXX) $(CXXFLAGS) -c gtestMessenger.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestMessenger.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestMessenger

gtestMuELoss: FORCE
ifeq ($(strip $(GOPT_ENABLE_MUELOSS)),YES)
	$(CXX) $(CXXFLAGS) -c gtestMuELoss.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestMuELoss.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestMuELoss
else
	@echo "You need to enable the MuELoss package to build the gtestMuELoss program"
endif

gtestNaturalIsotopes: FORCE
	$(CXX) $(CXXFLAGS) -c gtestNaturalIsotopes.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestNaturalIsotopes.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestNaturalIsotopes

gtestRewght: FORCE
ifeq ($(strip $(GOPT_ENABLE_RWGHT)),YES)
	$(CXX) $(CXXFLAGS) -c gtestRewght.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestRewght.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestRewght
else
	@echo "You need to enable event reweighting first!"
endif
                   
gtestRegistry: FORCE
	$(CXX) $(CXXFLAGS) -c gtestRegistry.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestRegistry.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestRegistry

gtestInteraction: FORCE
	$(CXX) $(CXXFLAGS) -c gtestInteraction.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestInteraction.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestInteraction

gtestResonances: FORCE
	$(CXX) $(CXXFLAGS) -c gtestResonances.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestResonances.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestResonances

gtestRwMarginalization: FORCE
ifeq ($(strip $(GOPT_ENABLE_RWGHT)),YES)
	$(CXX) $(CXXFLAGS) -c gtestRwMarginalization.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestRwMarginalization.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestRwMarginalization
else
	@echo "You need to enable event reweighting first!"
endif

gtestKPhaseSpace: FORCE
	$(CXX) $(CXXFLAGS) -c gtestKPhaseSpace.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestKPhaseSpace.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestKPhaseSpace

gtestROOTGeometry: FORCE
ifeq ($(strip $(GOPT_ENABLE_GEOM_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c gtestROOTGeometry.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gtestROOTGeometry.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestROOTGeometry
else
	@echo "You need to enable the geometry drivers to build the gtestROOTGeometry program"
endif

#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gtestAlgorithms 	
	$(RM) $(GENIE_BIN_PATH)/gtestBLI2DUnifGrid	
	$(RM) $(GENIE_BIN_PATH)/gtestCmdLnArg		
	$(RM) $(GENIE_BIN_PATH)/gtestConfigPool		
	$(RM) $(GENIE_BIN_PATH)/gtestDecay		
	$(RM) $(GENIE_BIN_PATH)/gtestDISSF		
	$(RM) $(GENIE_BIN_PATH)/gtestElFormFactors
	$(RM) $(GENIE_BIN_PATH)/gtestEventLoop
	$(RM) $(GENIE_BIN_PATH)/gtestFluxAstro
	$(RM) $(GENIE_BIN_PATH)/gtestFluxAtmo
	$(RM) $(GENIE_BIN_PATH)/gtestFluxSimple
	$(RM) $(GENIE_BIN_PATH)/gtestFGPauliBlockSuppr
	$(RM) $(GENIE_BIN_PATH)/gtestGiBUUData
	$(RM) $(GENIE_BIN_PATH)/gtestHadronization	
	$(RM) $(GENIE_BIN_PATH)/gtestINukeHadroData	
	$(RM) $(GENIE_BIN_PATH)/gtestMessenger		
	$(RM) $(GENIE_BIN_PATH)/gtestNaturalIsotopes	
	$(RM) $(GENIE_BIN_PATH)/gtestPDFLIB		
	$(RM) $(GENIE_BIN_PATH)/gtestPREM		
	$(RM) $(GENIE_BIN_PATH)/gtestFermiP		
	$(RM) $(GENIE_BIN_PATH)/gtestRewght		
	$(RM) $(GENIE_BIN_PATH)/gtestRegistry		
	$(RM) $(GENIE_BIN_PATH)/gtestInteraction	
	$(RM) $(GENIE_BIN_PATH)/gtestResonances		
	$(RM) $(GENIE_BIN_PATH)/gtestRwMarginalization	
	$(RM) $(GENIE_BIN_PATH)/gtestKPhaseSpace	
ifeq ($(strip $(GOPT_ENABLE_MUELOSS)),YES)
	$(RM) $(GENIE_BIN_PATH)/gtestMuELoss		
endif
ifeq ($(strip $(GOPT_ENABLE_GEOM_DRIVERS)),YES)
	$(RM) $(GENIE_BIN_PATH)/gtestROOTGeometry		
endif

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestAlgorithms 	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestBLI2DUnifGrid	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestCmdLnArg		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestConfigPool		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestDecay		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestDISSF		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestElFormFactors
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestEventLoop
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestFluxAstro
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestFluxAtmo
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestFluxSimple
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestFGPauliBlockSuppr
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestGiBUUData
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestHadronization	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestINukeHadroData
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestMessenger		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestNaturalIsotopes		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestPDFLIB		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestPREM		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestFermiP		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestRewght		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestRegistry		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestInteraction	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestResonances		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestRwMarginalization		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestKPhaseSpace	
ifeq ($(strip $(GOPT_ENABLE_MUELOSS)),YES)
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestMuELoss		
endif
ifeq ($(strip $(GOPT_ENABLE_GEOM_DRIVERS)),YES)
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestROOTGeometry		
endif


FORCE:

# DO NOT DELETE
