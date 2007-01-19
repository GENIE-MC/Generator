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
GENIE_DICTS = 
GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)

#LIBRARIES := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)
LIBRARIES := $(GENIE_LIBS) $(LIBRARIES) 

LIBS_AND_HDRS := $(GENIE_DICTS) $(LIBRARIES) $(INCLUDES)

EXE =	$(GENIE_BIN_PATH)/gevgen   \
	$(GENIE_BIN_PATH)/gmkspl   \
	$(GENIE_BIN_PATH)/gmxpl    \
	$(GENIE_BIN_PATH)/gntpc    \
	$(GENIE_BIN_PATH)/gsplt	   \
	$(GENIE_BIN_PATH)/gmctest  \
	$(GENIE_BIN_PATH)/gviewer

TGT =    gEvGen		\
	 gMakeSplines	\
	 gMaxPathLengths\
	 gNtpConv	\
	 gSplinePlotter	\
	 gMCSampleTest	\
	 gViewer	\
	 gMaxPathLengths\
	 gViewer

all: $(TGT)

# gevgen event generation driver
#
gEvGen: FORCE
	$(CXX) $(CXXFLAGS) -c gEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen 

# gmkspl x-section spline building utility
#
gMakeSplines: FORCE
	$(CXX) $(CXXFLAGS) -c gMakeSplines.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gMakeSplines.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmkspl 

# gmxpl utility computing maximum path lengths for a given root geometry
#
gMaxPathLengths: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c gMaxPathLengths.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gMaxPathLengths.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmxpl
else
	@echo "You need to enable the flux drivers to build the gmxpl utility"
endif


# gntpc ntuple conversion utility
#
gNtpConv: FORCE
	$(CXX) $(CXXFLAGS) -c gNtpConv.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gNtpConv.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gntpc 

# gsplt spline plotter utility
#
gSplinePlotter: FORCE
	$(CXX) $(CXXFLAGS) -c gSplinePlotter.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gSplinePlotter.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gsplt

# gmctest MC sample testing utility
#
gMCSampleTest: FORCE
	$(CXX) $(CXXFLAGS) -c gMCSampleTest.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gMCSampleTest.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmctest

# gviewer MC event viewing utility
#
gViewer: FORCE
ifeq ($(strip $(GOPT_ENABLE_VIEWER)),YES)
	$(CXX) $(CXXFLAGS) -c gViewer.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gViewer.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gviewer
else
	@echo "You need to enable the viewer to build the gviewer utility"
endif


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
