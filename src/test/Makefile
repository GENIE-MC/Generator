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

GENIE_DICTS =   $(GENIE)/src/Algorithm/_ROOT_DICT_Algorithm.o \
		$(GENIE)/src/BaryonResonance/_ROOT_DICT_BaryonResonance.o \
		$(GENIE)/src/Base/_ROOT_DICT_Base.o \
		$(GENIE)/src/BodekYang/_ROOT_DICT_BodekYang.o \
		$(GENIE)/src/Charm/_ROOT_DICT_Charm.o \
		$(GENIE)/src/CrossSections/_ROOT_DICT_CrossSections.o \
		$(GENIE)/src/Decay/_ROOT_DICT_Decay.o \
		$(GENIE)/src/Elastic/_ROOT_DICT_Elastic.o \
		$(GENIE)/src/EventGeneration/_ROOT_DICT_EventGeneration.o \
		$(GENIE)/src/Fragmentation/_ROOT_DICT_Fragmentation.o \
		$(GENIE)/src/Interaction/_ROOT_DICT_Interaction.o \
		$(GENIE)/src/InverseMuonDecay/_ROOT_DICT_InverseMuDecay.o \
		$(GENIE)/src/LlewellynSmith/_ROOT_DICT_LlewellynSmith.o \
		$(GENIE)/src/Ntuple/_ROOT_DICT_Ntuple.o \
		$(GENIE)/src/Nuclear/_ROOT_DICT_Nuclear.o \
		$(GENIE)/src/Numerical/_ROOT_DICT_Numerical.o \
		$(GENIE)/src/PartonModel/_ROOT_DICT_PartonModel.o \
		$(GENIE)/src/Paschos/_ROOT_DICT_Paschos.o \
		$(GENIE)/src/PDG/_ROOT_DICT_PDG.o \
		$(GENIE)/src/PDF/_ROOT_DICT_PDF.o \
		$(GENIE)/src/Radiative/_ROOT_DICT_Radiative.o \
		$(GENIE)/src/Registry/_ROOT_DICT_Registry.o \
		$(GENIE)/src/ReinSeghal/_ROOT_DICT_ReinSeghal.o \
		$(GENIE)/src/Viewer/_ROOT_DICT_Viewer.o \
		$(GENIE)/src/Utils/_ROOT_DICT_Utils.o \
		$(GENIE)/src/ValidityDoc/_ROOT_DICT_ValidityDoc.o 

WITH_LIB   = -L $(ROOTSYS)/lib $(ROOT_LIBS) $(ROOT_GLIBS) \
		$(PYTHIA_LIBS) -lEG -lEGPythia6 \
		$(XML_LIBS)   \
		$(GEN_LIBS) \
		$(CERN_LIBS)  \
		$(LOG_LIBS) \
		$(GPT_LIB)

WITH_INC   = $(LOG_INC) $(XML_INC) $(ROOT_INC) $(GENIE_INC)

WITH_DICT  = $(GENIE_DICTS)

TGT =	testAlgorithms		\
 	testBaryonResDataSet	\
 	testConfigPool		\
 	testDecay		\
 	testEvGen       	\
	testMessenger		\
	testNumerical		\
	testPDFLIB		\
	testRegistry		\
	testScatteringParams    \
	makeValidityDoc		\
	gviewer

EXE =	$(GENIE_BIN_PATH)/testAlgorithms	\
 	$(GENIE_BIN_PATH)/testBaryonResDataSet	\
 	$(GENIE_BIN_PATH)/testConfigPool	\
 	$(GENIE_BIN_PATH)/testDecay		\
 	$(GENIE_BIN_PATH)/testEvGen       	\
	$(GENIE_BIN_PATH)/testMessenger		\
	$(GENIE_BIN_PATH)/testNumerical		\
	$(GENIE_BIN_PATH)/testPDFLIB		\
	$(GENIE_BIN_PATH)/testRegistry		\
	$(GENIE_BIN_PATH)/testScatteringParams  \
	$(GENIE_BIN_PATH)/makeValidityDoc	\
	$(GENIE_BIN_PATH)/gviewer

all: $(TGT)

all-no-test: makeValidityDoc gviewer testEvGen

testBaryonResDataSet: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testBaryonResDataSet \
			testBaryonResDataSet.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

testConfigPool: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testConfigPool \
			testConfigPool.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

testDecay: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testDecay \
			testDecay.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

testEvGen: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/testEvGen \
			testEvGen.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

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

makeValidityDoc: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/makeValidityDoc \
			 makeValidityDoc.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

gviewer: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gviewer \
			gviewer.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)


#################### CLEANING

clean: FORCE
	$(RM) *.o *~ core $(EXE)

clean--testBaryonResDataSet: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testBaryonResDataSet

clean--testConfigPool: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testConfigPool

clean--testDecay: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testDecay

clean--testEvGen: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/testEvGen

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

clean--makeValidityDoc: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/makeValidityDoc

clean--gviewer: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gviewer


FORCE:

# DO NOT DELETE
