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

LIBS_AND_HDRS := $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

EXE =	$(GENIE_BIN_PATH)/gevgen   \
	$(GENIE_BIN_PATH)/gmkspl   \
	$(GENIE_BIN_PATH)/gmxpl    \
	$(GENIE_BIN_PATH)/gntpc    \
	$(GENIE_BIN_PATH)/gsplt	   \
	$(GENIE_BIN_PATH)/gmctest  \
	$(GENIE_BIN_PATH)/gviewer

ALLTGT = gEvGen		\
	 gMakeSplines	\
	 gMaxPathLengths\
	 gNtpConv	\
	 gSplinePlotter	\
	 gMCSampleTest	\
	 gViewer		

STDTGT = gEvGen		\
	 gNtpConv	\
	 gMakeSplines   \
	 gSplinePlotter	\
	 gMCSampleTest	

ifdef GOPT_FLUX_DRIVERS
  STDTGT := $(STDTGT) gMaxPathLengths
endif
ifdef GOPT_VIEWER
  STDTGT := $(STDTGT) gViewer
endif

all: $(STDTGT)

gEvGen: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gevgen gEvGen.cxx $(LIBS_AND_HDRS)

gMakeSplines: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gmkspl gMakeSplines.cxx $(LIBS_AND_HDRS)

gMaxPathLengths: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gmxpl gMaxPathLengths.cxx $(LIBS_AND_HDRS)

gNtpConv: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gntpc gNtpConv.cxx $(LIBS_AND_HDRS)

gSplinePlotter: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gsplt gSplinePlotter.cxx $(LIBS_AND_HDRS)

gMCSampleTest: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gmctest gMCSampleTest.cxx $(LIBS_AND_HDRS)

gViewer: FORCE
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/gviewer gViewer.cxx $(LIBS_AND_HDRS)


#################### CLEANING

clean: FORCE
	$(RM) *.o *~ core $(EXE)

clean--gEvGen: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevgen

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
