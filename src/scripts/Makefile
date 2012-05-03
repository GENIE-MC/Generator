# Makefile for GENIE Neutrino Generator
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

doxygen: FORCE
	$(DOXYGEN) dox/genie-doxygen;

install: FORCE
	cp setup/genie-config $(GENIE_BIN_PATH); \
	cp gcint/genie $(GENIE_BIN_PATH); \
	chmod 755 $(GENIE_BIN_PATH)/genie-config; \
	chmod 755 $(GENIE_BIN_PATH)/genie; 

#
#install_db: FORCE
#	cp setup/gnuvld_dbbootstrap $(GENIE_BIN_PATH); \
#	chmod 755 $(GENIE_BIN_PATH)/gnuvld_dbbootstrap 
#

clean: FORCE
	$(RM) $(GENIE_BIN_PATH)/genie-config 
	$(RM) $(GENIE_BIN_PATH)/genie

#
#clean_db: FORCE
#	$(RM) $(GENIE_BIN_PATH)/gnuvld_dbbootstrap 
#

distclean: FORCE
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/genie-config 
	$(RM) $(GENIE_BIN_INSTALLATION_PATH)/genie

FORCE:

