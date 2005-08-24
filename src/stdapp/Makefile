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
		$(GEN_LIBS)   \
		$(CERN_LIBS)  \
		$(LOG_LIBS)   \
		$(GPT_LIB)

WITH_INC   = $(LOG_INC) $(XML_INC) $(ROOT_INC) $(GENIE_INC)

WITH_DICT  = $(GENIE_DICTS)

TGT =	gEvGen		\
	gVldDoc		\
	gViewer		

EXE =	$(GENIE_BIN_PATH)/gevgen  \
	$(GENIE_BIN_PATH)/gvlddoc \
	$(GENIE_BIN_PATH)/gviewer

all: $(TGT)

gEvGen: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gevgen \
			gEvGen.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

gVldDoc: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gvlddoc \
			 gVldDoc.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

gViewer: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gviewer \
			gViewer.cxx $(WITH_DICT) $(WITH_LIB) $(WITH_INC)

#################### CLEANING

clean: FORCE
	$(RM) *.o *~ core $(EXE)

clean--gEvGen: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevgen

clean--gVldDoc: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gvlddoc

clean--gviewer: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gviewer

FORCE:

# DO NOT DELETE
