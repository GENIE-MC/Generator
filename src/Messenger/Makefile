#
# Makefile for GENIE Neutrino Generator
#
# Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Include machine specific flags and locations (inc. files & libs)
include $(GENIE)/src/Make.include

GENIES	  := $(wildcard *.cxx)
CINTLIST  := $(addsuffix .h, $(basename $(GENIES)))
CINTLISTO := $(addsuffix .$(ObjSuf), $(basename $(CINTLIST)))
GENIEO 	   = $(addsuffix .$(ObjSuf), $(basename $(GENIES)))
GENIEbS   := $(wildcard *.cxx)
GENIEbO    = $(addsuffix .$(ObjSuf), $(basename $(GENIEbS)))
GENIESO   := $(GENIE_MSG_SHARED_LIB)

all: lib link

lib: $(GENIEO) FORCE
	$(SOCMD) $(SOFLAGS) $(SOMINF) \
		$(GENIEO) $(OutPutOpt) $(GENIE_LIB_PATH)/$(GENIE_MSG_LIBNAME)$(GLIBVERSION).$(DllSuf)

link: FORCE
	@echo "Creating symbolic link to shared library ..."
ifeq ($(strip $(GOPT_ENABLE_DYLIBVERSION)),YES)
	ln -sf $(GENIE_MSG_LIBNAME)$(GLIBVERSION).$(DllSuf) $(GENIE_LIB_PATH)/$(GENIE_MSG_LIBNAME).$(DllSuf)
endif
ifeq ($(strip $(ARCH)),macosx)
	ln -sf $(GENIE_MSG_LIBNAME)$(GLIBVERSION).$(DllSuf) $(GENIE_LIB_PATH)/$(GENIE_MSG_LIBNAME).$(DllLinkSuf)
else
endif

purge: FORCE
	$(RM) $(GENIEbO) *~ core a.out Makefile.temp

clean: FORCE
	$(RM) *.o *~ core a.out $(GENIE_LIB_PATH)/$(GENIE_MSG_LIBNAME)* \
		Cint.* Makefile.temp

distclean: FORCE
	$(RM) $(GENIE_LIB_INSTALLATION_PATH)/$(GENIE_MSG_LIBNAME)*

FORCE:

