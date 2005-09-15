#########################################################
# Makefile for Neutrino Generator Validator binaries
#
# Costas Andreopoulos, 2003-08-20
#########################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/Make.include

# Add extra libraries needed for these executables
LIBRARIES := $(GENIE_LIBRARIES) $(NUVLD_LIBRARIES) $(LIBRARIES) $(CERN_LIBRARIES)

EXE = $(GENIE_BIN_PATH)/parsexml \
      $(GENIE_BIN_PATH)/filldbase \
      $(GENIE_BIN_PATH)/nuvld \
      $(GENIE_BIN_PATH)/ngevgen 

all: parsexml filldbase nuvld neugen_evgen

parsexml:
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/parsexml parsexml.cxx $(LIBRARIES) $(INCLUDES) 

filldbase:
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/filldbase filldbase.cxx $(LIBRARIES) $(INCLUDES) 

nuvld:
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/nuvld nuvld_gui.cxx $(LIBRARIES) $(INCLUDES) 

neugen_evgen:
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/ngevgen neugen_evgen.cxx $(LIBRARIES) $(INCLUDES) 

clean: FORCE
	$(RM) *.o *~ core $(EXE)

FORCE:

# DO NOT DELETE
