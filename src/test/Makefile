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

WITH_LIB   = -L $(ROOTSYS)/lib $(ROOT_LIBS) $(ROOT_GLIBS) \
		$(PYTHIA_LIBS) -lEG -lEGPythia6 \
		$(XML_LIBS)   \
		$(GEN_LIBS) \
		$(CERN_LIBS)  \
		$(LOG_LIBS) \
		$(GPT_LIB)

WITH_INC   = $(LOG_INC) $(XML_INC) $(ROOT_INC) $(GENIE_INC)

WITH_DICT  = $(GENIE_DICTS)

TGT =	testBaryonResDataSet	\
 	testConfigPool		\
 	testDecay		\
	testMessenger		\
	testNumerical		\
	testPDFLIB		\
	testRegistry		\
	testScatteringParams    \
	testROOTGeometry

EXE =	$(GENIE_BIN_PATH)/testBaryonResDataSet	\
 	$(GENIE_BIN_PATH)/testConfigPool	\
 	$(GENIE_BIN_PATH)/testDecay		\
	$(GENIE_BIN_PATH)/testMessenger		\
	$(GENIE_BIN_PATH)/testNumerical		\
	$(GENIE_BIN_PATH)/testPDFLIB		\
	$(GENIE_BIN_PATH)/testRegistry		\
	$(GENIE_BIN_PATH)/testScatteringParams  \
	$(GENIE_BIN_PATH)/testROOTGeometry  

all: $(TGT)

testBaryonResDataSet: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testBaryonResDataSet \
			testBaryonResDataSet.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

testConfigPool: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testConfigPool \
			testConfigPool.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

testDecay: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testDecay \
			testDecay.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

testPDFLIB: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testPDFLIB \
			testPDFLIB.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

testMessenger: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testMessenger \
			testMessenger.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

testNumerical: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testNumerical \
			testNumerical.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

testRegistry: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testRegistry \
			testRegistry.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

testScatteringParams: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testScatteringParams \
			testScatteringParams.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

testROOTGeometry: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testROOTGeometry \
			testROOTGeometry.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

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

clean--testRegistry: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testRegistry

clean--testScatteringParams: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testScatteringParams

clean--testROOTGeometry: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testROOTGeometry

FORCE:

# DO NOT DELETE
