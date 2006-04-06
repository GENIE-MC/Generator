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
GENIE_LIBS  = $(shell $(GENIE)/src/scripts/genie-config --libs)

LIBRARIES := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

ifdef GOPT_GEANT_INTERFACE
  LIBRARIES := $(LIBRARIES) $(GEANT4_LIBRARIES) $(CLHEP_LIBRARIES)
endif

EXE =	$(GENIE_BIN_PATH)/gevgen   \
	$(GENIE_BIN_PATH)/gvlddoc  \
	$(GENIE_BIN_PATH)/gmkspl   \
	$(GENIE_BIN_PATH)/gmxpl    \
	$(GENIE_BIN_PATH)/gntpc    \
	$(GENIE_BIN_PATH)/gsplt	   \
	$(GENIE_BIN_PATH)/gmctest  \
	$(GENIE_BIN_PATH)/gviewer

ALLTGT = gEvGen		\
	 gVldDoc	\
	 gMakeSplines	\
	 gMaxPathLengths\
	 gNtpConv	\
	 gSplinePlotter	\
	 gMCSampleTest	\
	 gViewer		

STDTGT = gEvGen		\
	 gMakeSplines	\
	 gNtpConv	\
	 gMaxPathLengths\
	 gSplinePlotter	\
	 gMCSampleTest	

ifdef GOPT_VALIDITY_DOC
  STDTGT := $(STDTGT) gVldDoc
endif
ifdef GOPT_VIEWER
  STDTGT := $(STDTGT) gViewer
endif

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

gNtpConv: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gntpc gNtpConv.cxx \
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

clean--gNtpConv: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gntpc

clean--gSplinePlotter: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gsplt

clean--gMCSampleTest: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gmctest

clean--gViewer: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gviewer

FORCE:

# DO NOT DELETE
