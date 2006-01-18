#########################################################
# Makefile for GENIE Neutrino Generator
#
# Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
# Last update: 2003-12-23
#########################################################

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/Make.include

doxygen:
	$(DOXYGEN) genie_doxygen;

install:
	cp db-bootstrap $(GENIE_BIN_PATH); \
	cp genie-config $(GENIE_BIN_PATH); \
	chmod 755 $(GENIE_BIN_PATH)/db-bootstrap; \
	chmod 755 $(GENIE_BIN_PATH)/genie-config; 

clean: FORCE
	$(RM) $(GENIE_BIN_PATH)/db-bootstrap $(GENIE_BIN_PATH)/genie-config

FORCE:

