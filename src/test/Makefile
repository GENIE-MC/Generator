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
include $(GENIE)/src/Make.include

#GENIE_DICTS = $(shell $(GENIE)/src/scripts/setup/genie-config --dicts)
#GENIE_DICTS =
GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)
#LIBRARIES := $(LIBRARIES) -L /usr/local/lib -ltcmalloc
#LIBS_AND_HDRS := $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

EXE =	$(GENIE_BIN_PATH)/gtestAlgorithms	\
 	$(GENIE_BIN_PATH)/gtestBaryonResDataSet	\
 	$(GENIE_BIN_PATH)/gtestConfigPool	\
 	$(GENIE_BIN_PATH)/gtestDecay		\
 	$(GENIE_BIN_PATH)/gtestDISSF		\
 	$(GENIE_BIN_PATH)/gtestHadronization	\
 	$(GENIE_BIN_PATH)/gtestIntranuke	\
	$(GENIE_BIN_PATH)/gtestMessenger	\
	$(GENIE_BIN_PATH)/gtestMuELoss		\
	$(GENIE_BIN_PATH)/gtestNumerical	\
	$(GENIE_BIN_PATH)/gtestPDFLIB		\
	$(GENIE_BIN_PATH)/gtestFermiP		\
	$(GENIE_BIN_PATH)/gtestMCJobDriver	\
	$(GENIE_BIN_PATH)/gtestFluxDriver	\
	$(GENIE_BIN_PATH)/gtestRegistry		\
	$(GENIE_BIN_PATH)/gtestInteraction	\
	$(GENIE_BIN_PATH)/gtestResonances	\
	$(GENIE_BIN_PATH)/gtestReWeight		\
	$(GENIE_BIN_PATH)/gtestKPhaseSpace	\
	$(GENIE_BIN_PATH)/gtestReadGNtuple	\
	$(GENIE_BIN_PATH)/gtestROOTGeometry     \
	$(GENIE_BIN_PATH)/gtestINukeHadroData

TGT =	gtestAlgorithms 	\
 	gtestBaryonResDataSet	\
 	gtestConfigPool		\
 	gtestDecay		\
 	gtestDISSF		\
 	gtestHadronization	\
 	gtestIntranuke		\
	gtestMessenger		\
	gtestNumerical		\
	gtestPDFLIB		\
	gtestFermiP		\
	gtestMCJobDriver	\
	gtestRegistry		\
	gtestInteraction	\
	gtestResonances		\
	gtestReWeight		\
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

gtestHadronization: FORCE
	$(CXX) $(CXXFLAGS) -c testHadronization.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testHadronization.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestHadronization

gtestIntranuke: FORCE
	$(CXX) $(CXXFLAGS) -c testIntranuke.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testIntranuke.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestIntranuke

gtestXSec: FORCE
	$(CXX) $(CXXFLAGS) -c testXSec.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testXSec.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestXSec

gtestPDFLIB: FORCE
	$(CXX) $(CXXFLAGS) -c testPDFLIB.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testPDFLIB.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestPDFLIB

gtestFermiP: FORCE
	$(CXX) $(CXXFLAGS) -c testFermiP.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testFermiP.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestFermiP

gtestMCJobDriver: FORCE
	$(CXX) $(CXXFLAGS) -c testMCJobDriver.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testMCJobDriver.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestMCJobDriver

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

gtestRegistry: FORCE
	$(CXX) $(CXXFLAGS) -c testRegistry.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testRegistry.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestRegistry

gtestInteraction: FORCE
	$(CXX) $(CXXFLAGS) -c testInteraction.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testInteraction.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestInteraction

gtestResonances: FORCE
	$(CXX) $(CXXFLAGS) -c testResonances.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testResonances.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestResonances

gtestReWeight: FORCE
	$(CXX) $(CXXFLAGS) -c testReWeight.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testReWeight.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestReWeight

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

gtestINukeHadroData: FORCE
	$(CXX) $(CXXFLAGS) -c testINukeHadroData.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) testINukeHadroData.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gtestINukeHadroData

#################### CLEANING

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(EXE)

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

clean--gtestHadronization: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestHadronization

clean--gtestHadronization: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestIntranuke

clean--gtestXSec: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestXSec

clean--gtestMessenger: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestMessenger

clean--gtestNumerical: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestNumerical

clean--gtestPDFLIB: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestPDFLIB

clean--gtestFermiP: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestFermiP

clean--gtestMCJobDriver: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestMCJobDriver

clean--gtestFluxDriver: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestFluxDriver

clean--gtestRegistry: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestRegistry

clean--gtestInteraction: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestInteraction

clean--gtestResonances: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestInteraction

clean--gtestReWeight: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestReWeight

clean--gtestKPhaseSpace: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestKPhaseSpace

clean--gtestReadGNtuple: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestReadGNtuple

clean--gtestROOTGeometry: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestROOTGeometry

clean--gtestINukeHadroData: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gtestSAIDHadronXSec

FORCE:

# DO NOT DELETE
