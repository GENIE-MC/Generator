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

ifdef GOPT_GEANT_INTERFACE
  LIBRARIES := $(LIBRARIES) $(GEANT4_LIBRARIES) $(CLHEP_LIBRARIES)
endif

###LIBRARIES := $(LIBRARIES) -L /usr/local/lib -ltcmalloc


EXE =	$(GENIE_BIN_PATH)/testBaryonResDataSet	\
 	$(GENIE_BIN_PATH)/testConfigPool	\
 	$(GENIE_BIN_PATH)/testDecay		\
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
	$(GENIE_BIN_PATH)/testReadGNtuple	\
	$(GENIE_BIN_PATH)/testXSec		\
	$(GENIE_BIN_PATH)/testROOTGeometry      \
	$(GENIE_BIN_PATH)/testSAIDHadronXSec

TGT =	testBaryonResDataSet	\
 	testConfigPool		\
 	testDecay		\
 	testHadronization	\
 	testXSec		\
	testMessenger		\
	testNumerical		\
	testPDFLIB		\
	testFermiP		\
	testMCJobDriver		\
	testRegistry		\
	testInteraction		\
	testResonances		\
	testReWeight		\
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

testBaryonResDataSet: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testBaryonResDataSet \
			testBaryonResDataSet.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testConfigPool: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testConfigPool \
			testConfigPool.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testDecay: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testDecay \
			testDecay.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testHadronization: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testHadronization \
		testHadronization.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testXSec: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testXSec \
   		testXSec.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testPDFLIB: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testPDFLIB \
			testPDFLIB.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testFermiP: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testFermiP \
			testFermiP.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testMCJobDriver: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testMCJobDriver \
			testMCJobDriver.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testFluxDriver: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testFluxDriver \
			testFluxDriver.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testMessenger: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testMessenger \
			testMessenger.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testMuELoss: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testMuELoss \
			testMuELoss.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testNumerical: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testNumerical \
			testNumerical.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testRegistry: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testRegistry \
			testRegistry.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testInteraction: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testInteraction \
			testInteraction.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testResonances: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testResonances \
			testResonances.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testReWeight: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testReWeight \
			testReWeight.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testReadGNtuple: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testReadGNtuple \
			testReadGNtuple.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testROOTGeometry: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testROOTGeometry \
			testROOTGeometry.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testSAIDHadronXSec: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testSAIDHadronXSec \
	              testSAIDHadronXSec.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

#################### CLEANING

clean: FORCE
	$(RM) *.o *~ core $(EXE)

clean--testBaryonResDataSet: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testBaryonResDataSet

clean--testConfigPool: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testConfigPool

clean--testDecay: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testDecay

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

clean--testReadGNtuple: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testReadGNtuple

clean--testROOTGeometry: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testROOTGeometry

clean--testSAIDHadronXSec: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testSAIDHadronXSec

FORCE:

# DO NOT DELETE
