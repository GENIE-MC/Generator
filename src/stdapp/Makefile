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

LIBRARIES := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

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

DEPLTGT = MACOSX_DEPLOYMENT_TARGET=10.4 c++

# gevgen event generation driver
#
gEvGen: FORCE
ifeq ($(strip $(ARCH)),macosx)
	$(CXX) $(CXXFLAGS) -c gEvGen.cxx $(INCLUDES)
	$(DEPLTGT) $(LDFLAGS) gEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen 
else
	$(CXX) $(CXXFLAGS) -o $(GENIE_BIN_PATH)/gmkspl gMakeSplines.cxx $(LIBS_AND_HDRS)
endif

# gmkspl x-section spline building utility
#
gMakeSplines: FORCE
ifeq ($(strip $(ARCH)),macosx)
	$(CXX) $(CXXFLAGS) -c gMakeSplines.cxx $(INCLUDES)
	$(DEPLTGT) $(LDFLAGS) gMakeSplines.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmkspl 
else
	$(CXX) $(CXXFLAGS) -o $(GENIE_BIN_PATH)/gmkspl gMakeSplines.cxx $(LIBS_AND_HDRS)
endif


# gmxpl utility computing maximum path lengths for a given root geometry
#
gMaxPathLengths: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
ifeq ($(strip $(ARCH)),macosx)
	$(CXX) $(CXXFLAGS) -c gMaxPathLengths.cxx $(INCLUDES)
	$(DEPLTGT) $(LDFLAGS) gMaxPathLengths.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmxpl
else
	$(CXX) $(CXXFLAGS) -o $(GENIE_BIN_PATH)/gmxpl gMaxPathLengths.cxx $(LIBS_AND_HDRS)
endif
else
	@echo "You need to enable the flux drivers to build the gmxpl utility"
endif


# gntpc ntuple conversion utility
#
gNtpConv: FORCE
ifeq ($(strip $(ARCH)),macosx)
	$(CXX) $(CXXFLAGS) -c gNtpConv.cxx $(INCLUDES)
	$(DEPLTGT) $(LDFLAGS) gNtpConv.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gntpc 
else
	$(CXX) $(CXXFLAGS) -o $(GENIE_BIN_PATH)/gntpc gNtpConv.cxx $(LIBS_AND_HDRS)
endif

# gsplt spline plotter utility
#
gSplinePlotter: FORCE
ifeq ($(strip $(ARCH)),macosx)
	$(CXX) $(CXXFLAGS) -c gSplinePlotter.cxx $(INCLUDES)
	$(DEPLTGT) $(LDFLAGS) gSplinePlotter.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gsplt
else
	$(CXX) $(CXXFLAGS) -o $(GENIE_BIN_PATH)/gsplt gSplinePlotter.cxx $(LIBS_AND_HDRS)
endif

# gmctest MC sample testing utility
#
gMCSampleTest: FORCE
ifeq ($(strip $(ARCH)),macosx)
	$(CXX) $(CXXFLAGS) -c gMCSampleTest.cxx $(INCLUDES)
	$(DEPLTGT) $(LDFLAGS) gMCSampleTest.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gmctest
else
	$(CXX) $(CXXFLAGS) -o $(GENIE_BIN_PATH)/gmctest gMCSampleTest.cxx $(LIBS_AND_HDRS)
endif

# gviewer MC event viewing utility
#
gViewer: FORCE
ifeq ($(strip $(GOPT_ENABLE_VIEWER)),YES)
ifeq ($(strip $(ARCH)),macosx)
	$(CXX) $(CXXFLAGS) -c gViewer.cxx $(INCLUDES)
	$(DEPLTGT) $(LDFLAGS) gViewer.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gviewer
else
	$(CXX) $(CXXFLAGS) -o $(GENIE_BIN_PATH)/gviewer gViewer.cxx $(LIBS_AND_HDRS)
endif
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
