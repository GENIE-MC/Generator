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

doxygen: FORCE
	$(DOXYGEN) dox/genie-doxygen;

install: FORCE
	cp setup/gnuvld-dbbootstrap $(GENIE_BIN_PATH); \
	cp setup/genie-config $(GENIE_BIN_PATH); \
	cp gcint/genie $(GENIE_BIN_PATH); \
	chmod 755 $(GENIE_BIN_PATH)/gnuvld-dbbootstrap; \
	chmod 755 $(GENIE_BIN_PATH)/genie-config; \
	chmod 755 $(GENIE_BIN_PATH)/genie; 

clean: FORCE
	$(RM) \
	$(GENIE_BIN_PATH)/gnuvld-dbbootstrap \
	$(GENIE_BIN_PATH)/genie-config \
	$(GENIE_BIN_PATH)/genie

distclean: FORCE
	$(RM) \
	$(GENIE_BIN_INSTALLATION_PATH)/gnuvld-dbbootstrap \
	$(GENIE_BIN_INSTALLATION_PATH)/genie-config \
	$(GENIE_BIN_INSTALLATION_PATH)/genie

FORCE:

