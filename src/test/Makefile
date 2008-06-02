#########################################################
# Makefile for Neutrino Generator Test Programs
#
# Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
# May 02, 2004
#########################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

#GENIE_DICTS = $(shell $(GENIE)/src/scripts/setup/genie-config --dicts)
GENIE_DICTS =
GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(LHAPDF_LIBRARIES) $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)


TGT =	gtestAlgorithms 	\
 	gtestBaryonResDataSet	\
 	gtestConfigPool		\
 	gtestDecay		\
 	gtestDISSF		\
 	gtestElFormFactors	\
 	gtestHadronization	\
 	gtestIntranuke		\
	gtestINukeHadroData     \
	gtestMessenger		\
	gtestNumerical		\
	gtestNaturalIsotopes	\
	gtestPDFLIB		\
	gtestLHAPDF		\
	gtestFermiP		\
	gtestRegistry		\
	gtestInteraction	\
	gtestResonances		\
	gtestReWeightP2P	\
	gtestReWeightLoop	\
	gtestKPhaseSpace	\
	gtestReadGNtuple		

all: $(TGT)

gtestAlgorithms: FORCE
	$(CXX) $(CXXFLAGS) -c testAlgorithms.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testAlgorithms.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestAlgorithms

gtestBaryonResDataSet: FORCE
	$(CXX) $(CXXFLAGS) -c testBaryonResDataSet.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testBaryonResDataSet.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestBaryonResDataSet

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

gtestLHAPDF: FORCE
	$(CXX) $(CXXFLAGS) -c testLHAPDF.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testLHAPDF.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestLHAPDF

gtestFermiP: FORCE
	$(CXX) $(CXXFLAGS) -c testFermiP.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testFermiP.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFermiP

gtestFluxDriver: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c testFluxDriver.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testFluxDriver.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFluxDriver
else
	@echo "You need to enable the flux drivers to build the gtestFluxDriver program"
endif

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

gtestReWeightP2P: FORCE
	$(CXX) $(CXXFLAGS) -c testReWeightP2P.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testReWeightP2P.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestReWeightP2P

gtestReWeightLoop: FORCE
	$(CXX) $(CXXFLAGS) -c testReWeightLoop.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testReWeightLoop.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestReWeightLoop

gtestKPhaseSpace: FORCE
	$(CXX) $(CXXFLAGS) -c testKPhaseSpace.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testKPhaseSpace.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestKPhaseSpace

gtestReadGNtuple: FORCE
	$(CXX) $(CXXFLAGS) -c testReadGNtuple.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testReadGNtuple.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestReadGNtuple

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
	$(RM) $(GENIE_BIN_PATH)/gtestConfigPool		
	$(RM) $(GENIE_BIN_PATH)/gtestDecay		
	$(RM) $(GENIE_BIN_PATH)/gtestDISSF		
	$(RM) $(GENIE_BIN_PATH)/gtestElFormFactors
	$(RM) $(GENIE_BIN_PATH)/gtestHadronization	
	$(RM) $(GENIE_BIN_PATH)/gtestIntranuke	
	$(RM) $(GENIE_BIN_PATH)/gtestINukeHadroData	
	$(RM) $(GENIE_BIN_PATH)/gtestMessenger		
	$(RM) $(GENIE_BIN_PATH)/gtestNumerical		
	$(RM) $(GENIE_BIN_PATH)/gtestNaturalIsotopes	
	$(RM) $(GENIE_BIN_PATH)/gtestPDFLIB		
	$(RM) $(GENIE_BIN_PATH)/gtestLHAPDF		
	$(RM) $(GENIE_BIN_PATH)/gtestFermiP		
	$(RM) $(GENIE_BIN_PATH)/gtestRegistry		
	$(RM) $(GENIE_BIN_PATH)/gtestInteraction	
	$(RM) $(GENIE_BIN_PATH)/gtestResonances		
	$(RM) $(GENIE_BIN_PATH)/gtestReWeightP2P	
	$(RM) $(GENIE_BIN_PATH)/gtestReWeightLoop		
	$(RM) $(GENIE_BIN_PATH)/gtestKPhaseSpace	
	$(RM) $(GENIE_BIN_PATH)/gtestReadGNtuple		

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestAlgorithms 	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestBaryonResDataSet	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestConfigPool		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestDecay		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestDISSF		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestElFormFactors
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestHadronization	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestIntranuke		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestINukeHadroData
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestMessenger		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestNumerical		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestNaturalIsotopes		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestPDFLIB		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestLHAPDF		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestFermiP		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestRegistry		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestInteraction	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestResonances		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestReWeightP2P		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestReWeightLoop		
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestKPhaseSpace	
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gtestReadGNtuple		

clean--gtestAlgorithms: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestAlgorithms

clean--gtestBaryonResDataSet: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestBaryonResDataSet

clean--gtestConfigPool: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestConfigPool

clean--gtestDecay: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestDecay

clean--gtestDISSF: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestDISSF

clean--gtestElFormFactors: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestElFormFactors

clean--gtestHadronization: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestHadronization

clean--gtestIntranuke: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestIntranuke

clean--gtestINukeHadroData: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestINukeHadroData

clean--gtestXSec: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestXSec

clean--gtestMessenger: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestMessenger

clean--gtestNumerical: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestNumerical

clean--gtestNaturalIsotopes: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestNaturalIsotopes

clean--gtestPDFLIB: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestPDFLIB

clean--gtestLHAPDF: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestLHAPDF

clean--gtestFermiP: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestFermiP

clean--gtestFluxDriver: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestFluxDriver

clean--gtestRegistry: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestRegistry

clean--gtestInteraction: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestInteraction

clean--gtestResonances: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestInteraction

clean--gtestReWeightP2P: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestReWeightP2P

clean--gtestReWeightLoop: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestReWeightLoop

clean--gtestKPhaseSpace: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestKPhaseSpace

clean--gtestReadGNtuple: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestReadGNtuple

clean--gtestROOTGeometry: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestROOTGeometry

FORCE:

# DO NOT DELETE
