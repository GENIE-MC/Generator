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

LIBRARIES := $(GENIE_LIBRARIES) $(LIBRARIES) $(CERN_LIBRARIES)

ifdef GOPT_GEANT_INTERFACE
  LIBRARIES := $(LIBRARIES) $(GEANT4_LIBRARIES) $(CLHEP_LIBRARIES)
endif

ALLTGT = gEvGen		\
	 gVldDoc	\
	 gMakeSplines	\
	 gMaxPathLengths\
	 g2Nuance	\
	 gSplinePlotter	\
	 gMCSampleTest	\
	 gViewer		

STDTGT = gEvGen		\
	 gVldDoc	\
	 gMakeSplines	\
	 gMaxPathLengths\
	 gSplinePlotter	

EXE =	$(GENIE_BIN_PATH)/gevgen   \
	$(GENIE_BIN_PATH)/gvlddoc  \
	$(GENIE_BIN_PATH)/gmkspl   \
	$(GENIE_BIN_PATH)/gmxpl    \
	$(GENIE_BIN_PATH)/g2nuance \
	$(GENIE_BIN_PATH)/gsplt	   \
	$(GENIE_BIN_PATH)/gmctest  \
	$(GENIE_BIN_PATH)/gviewer

all: $(STDTGT)

gEvGen: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gevgen gEvGen.cxx \
			$(GENIE_DICTS) $(LIBRARIES) $(INCLUDES) 

gVldDoc: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gvlddoc gVldDoc.cxx \
			$(GENIE_DICTS) $(LIBRARIES) $(INCLUDES) 

gMakeSplines: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gmkspl gMakeSplines.cxx \
			$(GENIE_DICTS) $(LIBRARIES) $(INCLUDES) 

gMaxPathLengths: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gmxpl gMaxPathLengths.cxx \
			$(GENIE_DICTS) $(LIBRARIES) $(INCLUDES) 

g2Nuance: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/g2nuance g2Nuance.cxx \
			$(GENIE_DICTS) $(LIBRARIES) $(INCLUDES) 

gSplinePlotter: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gsplt gSplinePlotter.cxx \
			$(GENIE_DICTS) $(LIBRARIES) $(INCLUDES) 

gMCSampleTest: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gmctest gMCSampleTest.cxx \
			$(GENIE_DICTS) $(LIBRARIES) $(INCLUDES) 

gViewer: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gviewer gViewer.cxx \
			$(GENIE_DICTS) $(LIBRARIES) $(INCLUDES) 


#################### CLEANING

clean: FORCE
	$(RM) *.o *~ core $(EXE)

clean--gEvGen: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevgen

clean--gVldDoc: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gvlddoc

clean--gMakeSplines: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gmkspl

clean--gMaxPathLengths: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gmxpl

clean--g2Nuance: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/g2nuance

clean--gSplinePlotter: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gsplt

clean--gMCSampleTest: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gmctest

clean--gViewer: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gviewer

FORCE:

# DO NOT DELETE
