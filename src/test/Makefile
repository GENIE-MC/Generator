#
# Makefile for Neutrino Generator Test Programs
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
 	gtestBaryonResDataSet	 \
	gtestBLI2DUnifGrid       \
 	gtestConfigPool		 \
 	gtestDecay		 \
 	gtestDISSF		 \
 	gtestElFormFactors	 \
 	gtestEventLoop 		 \
 	gtestFluxAstro 		 \
 	gtestFluxAtmo 		 \
 	gtestFluxSimple 	 \
        gtestGiBUUData           \
	gtestHadronization	 \
	gtestINukeHadroData      \
	gtestMessenger		 \
	gtestNumerical		 \
	gtestNaturalIsotopes	 \
	gtestPDFLIB		 \
	gtestQENuclSuppression   \
	gtestPREM		 \
	gtestROOTGeometry	 \
	gtestFermiP		 \
	gtestRegistry		 \
	gtestRecoINukeHadroFates \
	gtestInteraction	 \
	gtestResonances		 \
	gtestRwMarginalization	 \
	gtestRwNuXSecHelperP2P	 \
	gtestRwNuXSecHelperLoop	 \
	gtestKPhaseSpace	
# 	gtestIntranuke		 

all: $(TGT)

gtestAlgorithms: FORCE
	$(CXX) $(CXXFLAGS) -c testAlgorithms.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testAlgorithms.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestAlgorithms

gtestBaryonResDataSet: FORCE
	$(CXX) $(CXXFLAGS) -c testBaryonResDataSet.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testBaryonResDataSet.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestBaryonResDataSet

gtestBLI2DUnifGrid: FORCE
	$(CXX) $(CXXFLAGS) -c testBLI2DUnifGrid.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testBLI2DUnifGrid.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestBLI2DUnifGrid

gtestConfigPool: FORCE
	$(CXX) $(CXXFLAGS) -c testConfigPool.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testConfigPool.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestConfigPool

gtestDecay: FORCE
	$(CXX) $(CXXFLAGS) -c testDecay.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testDecay.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestDecay

gtestDISSF: FORCE
	$(CXX) $(CXXFLAGS) -c testDISSF.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testDISSF.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestDISSF

gtestElFormFactors: FORCE
	$(CXX) $(CXXFLAGS) -c testElFormFactors.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testElFormFactors.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestElFormFactors

gtestEventLoop: FORCE
	$(CXX) $(CXXFLAGS) -c testEventLoop.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testEventLoop.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestEventLoop

gtestFluxAstro: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c testFluxAstro.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testFluxAstro.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFluxAstro
else
	@echo "You need to enable the flux drivers to build the gtestFluxAstro program"
endif

gtestFluxAtmo: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c testFluxAtmo.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testFluxAtmo.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFluxAtmo
else
	@echo "You need to enable the flux drivers to build the gtestFluxAtmo program"
endif

gtestFluxSimple: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c testFluxSimple.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testFluxSimple.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFluxSimple
else
	@echo "You need to enable the flux drivers to build the gtestFluxSimple program"
endif

gtestGiBUUData: FORCE
	$(CXX) $(CXXFLAGS) -c testGiBUUData.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testGiBUUData.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestGiBUUData

gtestHadronization: FORCE
	$(CXX) $(CXXFLAGS) -c testHadronization.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testHadronization.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestHadronization

gtestIntranuke: FORCE
	$(CXX) $(CXXFLAGS) -c testIntranuke.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testIntranuke.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestIntranuke

gtestINukeHadroData: FORCE
	$(CXX) $(CXXFLAGS) -c testINukeHadroData.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testINukeHadroData.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestINukeHadroData

gtestXSec: FORCE
	$(CXX) $(CXXFLAGS) -c testXSec.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testXSec.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestXSec

gtestPDFLIB: FORCE
	$(CXX) $(CXXFLAGS) -c testPDFLIB.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testPDFLIB.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestPDFLIB

gtestQENuclSuppression: FORCE
	$(CXX) $(CXXFLAGS) -c testQENuclSuppression.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testQENuclSuppression.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestQENuclSuppression

gtestPREM: FORCE
	$(CXX) $(CXXFLAGS) -c testPREM.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testPREM.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestPREM

gtestFermiP: FORCE
	$(CXX) $(CXXFLAGS) -c testFermiP.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testFermiP.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFermiP

gtestMessenger: FORCE
	$(CXX) $(CXXFLAGS) -c testMessenger.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testMessenger.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestMessenger

gtestMuELoss: FORCE
ifeq ($(strip $(GOPT_ENABLE_MUELOSS)),YES)
	$(CXX) $(CXXFLAGS) -c testMuELoss.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testMuELoss.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestMuELoss
else
	@echo "You need to enable the MuELoss package to build the gtestMuELoss program"
endif

gtestNumerical: FORCE
	$(CXX) $(CXXFLAGS) -c testNumerical.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testNumerical.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestNumerical

gtestNaturalIsotopes: FORCE
	$(CXX) $(CXXFLAGS) -c testNaturalIsotopes.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testNaturalIsotopes.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestNaturalIsotopes

gtestRegistry: FORCE
	$(CXX) $(CXXFLAGS) -c testRegistry.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testRegistry.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestRegistry

gtestInteraction: FORCE
	$(CXX) $(CXXFLAGS) -c testInteraction.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testInteraction.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestInteraction

gtestResonances: FORCE
	$(CXX) $(CXXFLAGS) -c testResonances.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testResonances.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestResonances

gtestRwMarginalization: FORCE
	$(CXX) $(CXXFLAGS) -c testRwMarginalization.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testRwMarginalization.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestRwMarginalization

gtestRwNuXSecHelperP2P: FORCE
	$(CXX) $(CXXFLAGS) -c testRwNuXSecHelperP2P.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testRwNuXSecHelperP2P.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestRwNuXSecHelperP2P

gtestRwNuXSecHelperLoop: FORCE
	$(CXX) $(CXXFLAGS) -c testRwNuXSecHelperLoop.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testRwNuXSecHelperLoop.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestRwNuXSecHelperLoop

gtestRecoINukeHadroFates: FORCE
	$(CXX) $(CXXFLAGS) -c testRecoINukeHadroFates.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testRecoINukeHadroFates.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestRecoINukeHadroFates

gtestKPhaseSpace: FORCE
	$(CXX) $(CXXFLAGS) -c testKPhaseSpace.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testKPhaseSpace.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestKPhaseSpace

gtestROOTGeometry: FORCE
ifeq ($(strip $(GOPT_ENABLE_GEOM_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c testROOTGeometry.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testROOTGeometry.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestROOTGeometry
else
	@echo "You need to enable the geometry drivers to build the gtestROOTGeometry program"
endif

#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core 
	$(RM) $(GENIE_BIN_PATH)/gtestAlgorithms 	
	$(RM) $(GENIE_BIN_PATH)/gtestBaryonResDataSet	
	$(RM) $(GENIE_BIN_PATH)/gtestBLI2DUnifGrid	
	$(RM) $(GENIE_BIN_PATH)/gtestConfigPool		
	$(RM) $(GENIE_BIN_PATH)/gtestDecay		
	$(RM) $(GENIE_BIN_PATH)/gtestDISSF		
	$(RM) $(GENIE_BIN_PATH)/gtestElFormFactors
	$(RM) $(GENIE_BIN_PATH)/gtestEventLoop
	$(RM) $(GENIE_BIN_PATH)/gtestFluxAstro
	$(RM) $(GENIE_BIN_PATH)/gtestFluxAtmo
	$(RM) $(GENIE_BIN_PATH)/gtestFluxSimple
	$(RM) $(GENIE_BIN_PATH)/gtestGiBUUData
	$(RM) $(GENIE_BIN_PATH)/gtestHadronization	
	$(RM) $(GENIE_BIN_PATH)/gtestIntranuke	
	$(RM) $(GENIE_BIN_PATH)/gtestINukeHadroData	
	$(RM) $(GENIE_BIN_PATH)/gtestMessenger		
	$(RM) $(GENIE_BIN_PATH)/gtestNumerical		
	$(RM) $(GENIE_BIN_PATH)/gtestNaturalIsotopes	
	$(RM) $(GENIE_BIN_PATH)/gtestPDFLIB		
	$(RM) $(GENIE_BIN_PATH)/gtestQENuclSuppression
	$(RM) $(GENIE_BIN_PATH)/gtestPREM		
	$(RM) $(GENIE_BIN_PATH)/gtestFermiP		
	$(RM) $(GENIE_BIN_PATH)/gtestRegistry		
	$(RM) $(GENIE_BIN_PATH)/gtestInteraction	
	$(RM) $(GENIE_BIN_PATH)/gtestResonances		
	$(RM) $(GENIE_BIN_PATH)/gtestRwMarginalization	
	$(RM) $(GENIE_BIN_PATH)/gtestRwNuXSecHelperP2P	
	$(RM) $(GENIE_BIN_PATH)/gtestRwNuXSecHelperLoop		
	$(RM) $(GENIE_BIN_PATH)/gtestRecoINukeHadroFates		
	$(RM) $(GENIE_BIN_PATH)/gtestKPhaseSpace	
ifeq ($(strip $(GOPT_ENABLE_MUELOSS)),YES)
	$(RM) $(GENIE_BIN_PATH)/gtestMuELoss		
endif
ifeq ($(strip $(GOPT_ENABLE_GEOM_DRIVERS)),YES)
	$(RM) $(GENIE_BIN_PATH)/gtestROOTGeometry		
endif

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestAlgorithms 	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestBaryonResDataSet	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestBLI2DUnifGrid	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestConfigPool		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestDecay		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestDISSF		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestElFormFactors
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestEventLoop
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestFluxAstro
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestFluxAtmo
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestFluxSimple
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestGiBUUData
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestHadronization	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestIntranuke		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestINukeHadroData
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestMessenger		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestNumerical		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestNaturalIsotopes		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestPDFLIB		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestQENuclSuppression
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestPREM		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestFermiP		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestRegistry		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestInteraction	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestResonances		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestRwMarginalization		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestRwNuXSecHelperP2P		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestRwNuXSecHelperLoop		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestRecoINukeHadroFates		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestKPhaseSpace	
ifeq ($(strip $(GOPT_ENABLE_MUELOSS)),YES)
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestMuELoss		
endif
ifeq ($(strip $(GOPT_ENABLE_GEOM_DRIVERS)),YES)
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestROOTGeometry		
endif


FORCE:

# DO NOT DELETE
