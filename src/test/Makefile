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

GENIE_DICTS = $(shell $(GENIE)/src/scripts/genie-config --dicts)
WITH_DICT   = $(GENIE_DICTS)

LIBRARIES := $(GENIE_LIBRARIES) $(LIBRARIES) $(CERN_LIBRARIES)

ifdef GOPT_GEANT_INTERFACE
  LIBRARIES := $(LIBRARIES) $(GEANT4_LIBRARIES) $(CLHEP_LIBRARIES)
endif

TGT =	testBaryonResDataSet	\
 	testConfigPool		\
 	testDecay		\
	testMessenger		\
	testNumerical		\
	testPDFLIB		\
	testMCJobDriver		\
	testRegistry		\
	testReWeight		\
	testReadGNtuple		\
	testROOTGeometry

EXE =	$(GENIE_BIN_PATH)/testBaryonResDataSet	\
 	$(GENIE_BIN_PATH)/testConfigPool	\
 	$(GENIE_BIN_PATH)/testDecay		\
	$(GENIE_BIN_PATH)/testMessenger		\
	$(GENIE_BIN_PATH)/testNumerical		\
	$(GENIE_BIN_PATH)/testPDFLIB		\
	$(GENIE_BIN_PATH)/testMCJobDriver	\
	$(GENIE_BIN_PATH)/testRegistry		\
	$(GENIE_BIN_PATH)/testReWeight		\
	$(GENIE_BIN_PATH)/testReadGNtuple	\
	$(GENIE_BIN_PATH)/testROOTGeometry  

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

testPDFLIB: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testPDFLIB \
			testPDFLIB.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testMCJobDriver: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testMCJobDriver \
			testMCJobDriver.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testMessenger: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testMessenger \
			testMessenger.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testNumerical: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testNumerical \
			testNumerical.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testRegistry: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testRegistry \
			testRegistry.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testReWeight: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testReWeight \
			testReWeight.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testReadGNtuple: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testReadGNtuple \
			testReadGNtuple.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

testROOTGeometry: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testROOTGeometry \
			testROOTGeometry.cxx $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

#################### CLEANING

clean: FORCE
	$(RM) *.o *~ core $(EXE)

clean--testBaryonResDataSet: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testBaryonResDataSet

clean--testConfigPool: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testConfigPool

clean--testDecay: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testDecay

clean--testMessenger: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testMessenger

clean--testNumerical: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testNumerical

clean--testPDFLIB: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testPDFLIB

clean--testMCJobDriver: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testPDFLIB

clean--testRegistry: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testRegistry

clean--testReWeight: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testReWeight

clean--testReadGNtuple: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testReadGNtuple

clean--testROOTGeometry: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testROOTGeometry

FORCE:

# DO NOT DELETE
