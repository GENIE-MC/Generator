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

#GENIE_DICTS = $(shell $(GENIE)/src/scripts/genie-config --dicts)
GENIE_DICTS =
GENIE_LIBS  = $(shell $(GENIE)/src/scripts/genie-config --libs)

LIBRARIES := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

###LIBRARIES := $(LIBRARIES) -L /usr/local/lib -ltcmalloc

LIBS_AND_HDRS := $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

EXE =	$(GENIE_BIN_PATH)/testAlgorithms	\
 	$(GENIE_BIN_PATH)/testBaryonResDataSet	\
 	$(GENIE_BIN_PATH)/testConfigPool	\
 	$(GENIE_BIN_PATH)/testDecay		\
 	$(GENIE_BIN_PATH)/testDISSF		\
 	$(GENIE_BIN_PATH)/testHadronization	\
	$(GENIE_BIN_PATH)/testMessenger		\
	$(GENIE_BIN_PATH)/testMuELoss		\
	$(GENIE_BIN_PATH)/testNumerical		\
	$(GENIE_BIN_PATH)/testPDFLIB		\
	$(GENIE_BIN_PATH)/testFermiP		\
	$(GENIE_BIN_PATH)/testMCJobDriver	\
	$(GENIE_BIN_PATH)/testFluxDriver	\
	$(GENIE_BIN_PATH)/testRegistry		\
	$(GENIE_BIN_PATH)/testInteraction	\
	$(GENIE_BIN_PATH)/testResonances	\
	$(GENIE_BIN_PATH)/testReWeight		\
	$(GENIE_BIN_PATH)/testKPhaseSpace	\
	$(GENIE_BIN_PATH)/testReadGNtuple	\
	$(GENIE_BIN_PATH)/testROOTGeometry      \
	$(GENIE_BIN_PATH)/testINukeHadroData

TGT =	testAlgorithms 		\
 	testBaryonResDataSet	\
 	testConfigPool		\
 	testDecay		\
 	testDISSF		\
 	testHadronization	\
	testMessenger		\
	testNumerical		\
	testPDFLIB		\
	testFermiP		\
	testMCJobDriver		\
	testRegistry		\
	testInteraction		\
	testResonances		\
	testReWeight		\
	testKPhaseSpace		\
	testReadGNtuple		

ifdef GOPT_MUELOSS
  TGT := $(TGT) testMuELoss
endif

ifdef GOPT_GEOM_DRIVERS
  TGT := $(TGT) testROOTGeometry
endif
ifdef GOPT_FLUX_DRIVERS
  TGT := $(TGT) testFluxDriver
endif

all: $(TGT)

testAlgorithms: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testAlgorithms testAlgorithms.cxx $(LIBS_AND_HDRS)

testBaryonResDataSet: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testBaryonResDataSet testBaryonResDataSet.cxx $(LIBS_AND_HDRS) 

testConfigPool: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testConfigPool testConfigPool.cxx $(LIBS_AND_HDRS) 

testDecay: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testDecay testDecay.cxx $(LIBS_AND_HDRS) 

testDISSF: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testDISSF testDISSF.cxx $(LIBS_AND_HDRS) 

testHadronization: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testHadronization testHadronization.cxx $(LIBS_AND_HDRS) 

testXSec: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testXSec testXSec.cxx $(LIBS_AND_HDRS) 

testPDFLIB: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testPDFLIB testPDFLIB.cxx $(LIBS_AND_HDRS) 

testFermiP: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testFermiP testFermiP.cxx $(LIBS_AND_HDRS)

testMCJobDriver: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testMCJobDriver testMCJobDriver.cxx $(LIBS_AND_HDRS)

testFluxDriver: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testFluxDriver testFluxDriver.cxx $(LIBS_AND_HDRS) 

testMessenger: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testMessenger testMessenger.cxx $(LIBS_AND_HDRS) 

testMuELoss: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testMuELoss testMuELoss.cxx $(LIBS_AND_HDRS) 

testNumerical: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testNumerical testNumerical.cxx $(LIBS_AND_HDRS) 

testRegistry: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testRegistry testRegistry.cxx $(LIBS_AND_HDRS)

testInteraction: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testInteraction testInteraction.cxx $(LIBS_AND_HDRS) 

testResonances: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testResonances testResonances.cxx $(LIBS_AND_HDRS) 

testReWeight: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testReWeight testReWeight.cxx $(LIBS_AND_HDRS) 

testKPhaseSpace: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testKPhaseSpace testKPhaseSpace.cxx $(LIBS_AND_HDRS) 

testReadGNtuple: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testReadGNtuple testReadGNtuple.cxx $(LIBS_AND_HDRS)

testROOTGeometry: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testROOTGeometry testROOTGeometry.cxx $(LIBS_AND_HDRS)

testINukeHadroData: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testINukeHadroData testINukeHadroData.cxx $(LIBS_AND_HDRS) 

#################### CLEANING

clean: FORCE
	$(RM) *.o *~ core $(EXE)

clean--testAlgorithms: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testAlgorithms

clean--testBaryonResDataSet: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testBaryonResDataSet

clean--testConfigPool: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testConfigPool

clean--testDecay: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testDecay

clean--testDISSF: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testDISSF

clean--testHadronization: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testHadronization

clean--testXSec: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testXSec

clean--testMessenger: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testMessenger

clean--testNumerical: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testNumerical

clean--testPDFLIB: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testPDFLIB

clean--testFermiP: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testFermiP

clean--testMCJobDriver: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testMCJobDriver

clean--testFluxDriver: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testFluxDriver

clean--testRegistry: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testRegistry

clean--testInteraction: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testInteraction

clean--testResonances: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testInteraction

clean--testReWeight: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testReWeight

clean--testKPhaseSpace: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testKPhaseSpace

clean--testReadGNtuple: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testReadGNtuple

clean--testROOTGeometry: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testROOTGeometry

clean--testINukeHadroData: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testSAIDHadronXSec

FORCE:

# DO NOT DELETE
