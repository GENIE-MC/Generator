####################################################################
# Makefile for the event generation driver customized for FNAL
# neutrino experiments
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
####################################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

## set extra stuff if $DK2NU is defined
#ifneq ($(shell ls -1 ${DK2NU}/lib/libdk2nuGenie* 2>/dev/null ),)
#    $(warning DK2NU/lib/libdk2nuGenie library exists, build against it)
#    DK2NULIBS  := -L $(DK2NU)/lib -ldk2nuTree -ldk2nuGenie
#    DK2NUFLAGS := -D__DK2NU_FLUX_DRIVER_AVAILABLE__ -I${DK2NU} -I${DK2NU}/include
#    CXXFLAGS   += $(DK2NUFLAGS)
#    # $(warning using $(DK2NULIBS))
#    # $(warning using DK2NUFLAGS=$(DK2NUFLAGS))
#else
#    $(warning building without Dk2Nu support)
#endif

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)
LIBRARIES  := $(GENIE_LIBS) $(DK2NULIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT = gevgen_fnal

all: $(TGT)

gevgen_fnal: FORCE
	$(CXX) $(CXXFLAGS) -c gFNALExptEvGen.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gFNALExptEvGen.o $(LIBRARIES) -o $(GENIE_BIN_PATH)/gevgen_fnal
	@echo "For historical purposes, add a NuMI event generation app sym link"
	ln -sf $(GENIE_BIN_PATH)/gevgen_fnal $(GENIE_BIN_PATH)/gevgen_numi

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core $(GENIE_BIN_PATH)/gevgen_fnal $(GENIE_BIN_PATH)/gevgen_numi

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/gevgen_fnal

FORCE:

# DO NOT DELETE
