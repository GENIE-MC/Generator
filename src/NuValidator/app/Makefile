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

APP_LIB   = -L $(ROOTSYS)/lib $(ROOT_LIBS) $(ROOT_GLIBS) \
                $(PYTHIA_LIBS) -lEG -lEGPythia6 \
                $(XML_LIBS)   \
                $(GEN_LIBS) \
                $(VLD_LIBS) \
                $(CERN_LIBS)  \
                $(LOG_LIBS)

APP_INC   = $(XML_INC) $(LOG_INC) $(ROOT_INC) $(NUVLD_INC) $(GENIE_INC)

EXE = $(GENIE_BIN_PATH)/parsexml \
      $(GENIE_BIN_PATH)/filldbase \
      $(GENIE_BIN_PATH)/nuvld \
      $(GENIE_BIN_PATH)/neugen_evgen 

all: parsexml filldbase nuvld neugen_evgen

parsexml:
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/parsexml parsexml.cxx $(APP_LIB) $(APP_INC) 

filldbase:
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/filldbase filldbase.cxx $(APP_LIB) $(APP_INC) 

nuvld:
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/nuvld nuvld_gui.cxx $(APP_LIB) $(APP_INC)

neugen_evgen:
	g++ -Wall -O2 -o $(GENIE_BIN_PATH)/ngevgen neugen_evgen.cxx $(APP_LIB) $(APP_INC)

clean: FORCE
	$(RM) *.o *~ core $(EXE)

FORCE:

# DO NOT DELETE
