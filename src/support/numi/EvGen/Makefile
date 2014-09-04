####################################################################
# Makefile for the event generation driver customized for NuMI
# beam-line experiments
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
####################################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

# set extra stuff if $DK2NU is defined
ifneq ($(shell ls -1 ${DK2NU}/lib/libdk2nuGenie* 2>/dev/null ),)
    $(warning DK2NU/lib/libdk2nuGenie library exists, build against it)
    DK2NULIBS  := -L $(DK2NU)/lib -ldk2nuTree -ldk2nuGenie
    DK2NUFLAGS := -D__DK2NU_FLUX_DRIVER_AVAILABLE__ -I${DK2NU} -I${DK2NU}/include
    CXXFLAGS   += $(DK2NUFLAGS)
    # $(warning using $(DK2NULIBS))
    # $(warning using DK2NUFLAGS=$(DK2NUFLAGS))
#else
#    $(warning building without Dk2Nu support)
endif

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(DK2NULIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = gevgen_numi

all: $(TGT)

gevgen_numi: FORCE
	$(CXX) $(CXXFLAGS) -c gNuMIExptEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gNuMIExptEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_numi
	@echo "For historical purposes, add symbolic link to the NuMI event generation app"
	ln -sf $(GENIE_BIN_PATH)/gevgen_numi $(GENIE_BIN_PATH)/gNuMIevgen

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevgen_numi $(GENIE_BIN_PATH)/gNuMIevgen

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_numi

FORCE:

# DO NOT DELETE
