#
# Makefile for GENIE Apps
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

GENIE_LIBS  = $(shell $(GENIE)/src/scripts/setup/genie-config --libs)

LIBRARIES  := $(GENIE_LIBS) $(LIBRARIES) $(CERN_LIBRARIES)

TGT =    grwghtzexpaxff		\
	 grwghtzexpdirect	


# reweighting for Zexpansion of axial form factor
#
grwghtzexpaxff: FORCE
	@echo "** Building grwghtzexpaxff"
ifeq ($(strip $(GOPT_ENABLE_RWGHT)),YES)
	$(CXX) $(CXXFLAGS) -c gRwghtZExpAxFF.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gRwghtZExpAxFF.o $(LIBRARIES) -o grwghtzexpaxff
else
	@echo "Event reweighting not enabled!"
endif

# reweighting for Zexpansion of axial form factor
#
grwghtzexpdirect: FORCE
	@echo "** Building grwghtzexpdirect"
ifeq ($(strip $(GOPT_ENABLE_RWGHT)),YES)
	$(CXX) $(CXXFLAGS) -c gRwghtZExpDirect.cxx $(INCLUDES)
	$(LD) $(LDFLAGS) gRwghtZExpDirect.o $(LIBRARIES) -o grwghtzexpdirect
else
	@echo "Event reweighting not enabled!"
endif

# CLEANING-UP

purge: FORCE
	$(RM) *.o *~ core 

clean: FORCE
	$(RM) *.o *~ core
	$(RM) grwghtzexpaxff		
	$(RM) grwghtzexpdirect

distclean: FORCE


FORCE:

# DO NOT DELETE
