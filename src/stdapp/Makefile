#
# Makefile for GENIE standard apps
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# May 02, 2004
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)

LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = $(GENIE_BIN_PATH)/gevgen \
	$(GENIE_BIN_PATH)/gevdump \
	$(GENIE_BIN_PATH)/gevpick \
	$(GENIE_BIN_PATH)/gmkspl \
	$(GENIE_BIN_PATH)/gspladd \
	$(GENIE_BIN_PATH)/gspl2root \
	$(GENIE_BIN_PATH)/gntpc \
	$(GENIE_BIN_PATH)/gevgen_hadron

GMXPL = 	$(GENIE_BIN_PATH)/gmxpl
all: $(TGT) $(GMXPL)

# Sadly the names of the source files and final executable names aren't the same so we have to map them
$(GENIE_BIN_PATH)/gevgen : gEvGen.o
$(GENIE_BIN_PATH)/gevdump: gEvDump.o
$(GENIE_BIN_PATH)/gevgen_hadron : gEvGenHadronNucleus.o
$(GENIE_BIN_PATH)/gevpick : gEvPick.o
$(GENIE_BIN_PATH)/gmkspl : gMakeSplines.o
$(GENIE_BIN_PATH)/gspladd : gSplineAdd.o
$(GENIE_BIN_PATH)/gspl2root : gSplineXml2Root.o
$(GENIE_BIN_PATH)/gmxpl : gMaxPathLengths.o
$(GENIE_BIN_PATH)/gntpc : gNtpConv.o

$(TGT):
	$(LD) $(LDFLAGS) $^ $(LIBRARIES) -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $^ $(INCLUDES)

gMaxPathLengths.o: gMaxPathLengths.cxx
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(CXX) $(CXXFLAGS) -c $^ $(INCLUDES)
else
	@echo "Not creating $@ because flux drivers are disabled"
endif

$(GENIE_BIN_PATH)/gmxpl:
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	$(LD) $(LDFLAGS) $^ $(LIBRARIES) -o $@
else
	@echo "Not creating $@ because flux drivers are disabled"
endif


#################### CLEANING

.PHONY: purge clean distclean

purge:
	$(RM) *.o *~ core

clean:
	$(RM) *.o *~ core
	$(foreach exe, $(TGT), $(RM) $(exe) ;)

distclean:
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevdump
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_hadron
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevpick
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gmkspl
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gspladd
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gspl2root
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gmxpl
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gntpc

# DO NOT DELETE
