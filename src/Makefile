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

# std build targets
BUILD_TARGETS = save-build-options generator-std-libs 

# build targets to be enabled by the user (nu-validator, 
# genie/geant interface, flux/geometry drivers doxygen documentation)

ifdef GOPT_FLUX_DRIVERS
  BUILD_TARGETS := $(BUILD_TARGETS) flux-drivers
endif
ifdef GOPT_GEOM_DRIVERS
  BUILD_TARGETS := $(BUILD_TARGETS) geom-drivers
endif
ifndef GOPT_NEUGEN
  BUILD_TARGETS := $(BUILD_TARGETS) dummy-neugen
endif
ifdef GOPT_NUVALIDATOR
  BUILD_TARGETS := $(BUILD_TARGETS) nuvld-std
endif
ifdef GOPT_GEANT_INTERFACE
  BUILD_TARGETS := $(BUILD_TARGETS) geanti
endif
ifdef GOPT_DOXYGEN_DOC
  BUILD_TARGETS := $(BUILD_TARGETS) doxygen-doc
endif

# after all libraries have been build, add some standard 
# build targets (standalone apps, etc)

BUILD_TARGETS := $(BUILD_TARGETS) generator-std-exe scripts

# define build targets

all: $(BUILD_TARGETS)

generator-std-libs:
	cd Algorithm;          make; cd ..; \
	cd BaryonResonance;    make; cd ..; \
	cd Base;               make; cd ..; \
	cd BodekYang;          make; cd ..; \
	cd Charm;              make; cd ..; \
	cd CrossSections;      make; cd ..; \
	cd Decay; 	       make; cd ..; \
	cd Elastic; 	       make; cd ..; \
	cd EVGCore;            make; cd ..; \
	cd EVGModules;         make; cd ..; \
	cd EVGDrivers;         make; cd ..; \
	cd Fragmentation;      make; cd ..; \
	cd GHEP;               make; cd ..; \
	cd ValidityDoc;        make; cd ..; \
	cd Interaction;        make; cd ..; \
	cd InverseMuonDecay;   make; cd ..; \
	cd LlewellynSmith;     make; cd ..; \
	cd Messenger;	       make; cd ..; \
	cd Ntuple;             make; cd ..; \
	cd Nuclear;            make; cd ..; \
	cd Numerical;          make; cd ..; \
	cd PartonModel;        make; cd ..; \
	cd Paschos; 	       make; cd ..; \
	cd PDF;                make; cd ..; \
	cd PDG;                make; cd ..; \
	cd Radiative;          make; cd ..; \
	cd Registry;           make; cd ..; \
	cd ReinSeghal;         make; cd ..; \
	cd Utils;              make; cd ..; \
	cd Viewer;             make; cd ..; \

flux-drivers:
	cd FluxDrivers;        make; cd ..;

geom-drivers:
	cd Geo;                make; cd ..; 

dummy-neugen: 
	cd NuValidator; make dummy-neugen; cd ..; 

nuvld-std: 
	cd NuValidator; make std; cd ..; 

nuvld-libs:
	cd NuValidator; make libs; cd ..; 

nuvld-exe:
	cd NuValidator; make applications; cd ..; 

geanti:
	cd GeantInterface; make; cd ..; \

doxygen-doc:
	cd ../doc/doxygen/; $(DOXYGEN) genie_doxygen; cd ../../src/

generator-std-exe:
	cd stdapp; make all; cd ..; \

generator-test-exe:
	cd test; make all; cd ..; \

scripts:
	cd scripts; make install; cd ..;

all-libs: generator-std-libs flux-drivers geom-drivers geanti nuvld-libs

save-build-options:
	perl scripts/save-build-env;

purge: FORCE
	cd Algorithm;          make purge; cd ..; \
	cd BaryonResonance;    make purge; cd ..; \
	cd Base;               make purge; cd ..; \
	cd BodekYang;          make purge; cd ..; \
	cd Charm;              make purge; cd ..; \
	cd CrossSections;      make purge; cd ..; \
	cd Decay; 	       make purge; cd ..; \
	cd Elastic; 	       make purge; cd ..; \
	cd EVGCore;            make purge; cd ..; \
	cd EVGModules;         make purge; cd ..; \
	cd EVGDrivers;         make purge; cd ..; \
	cd FluxDrivers;        make purge; cd ..; \
	cd Fragmentation;      make purge; cd ..; \
	cd GHEP;               make purge; cd ..; \
	cd GeantInterface;     make purge; cd ..; \
	cd Geo;                make purge; cd ..; \
	cd ValidityDoc;        make purge; cd ..; \
	cd Interaction;        make purge; cd ..; \
	cd InverseMuonDecay;   make purge; cd ..; \
	cd LlewellynSmith;     make purge; cd ..; \
	cd Messenger;	       make purge; cd ..; \
	cd Nuclear;            make purge; cd ..; \
	cd Ntuple;             make purge; cd ..; \
	cd Numerical;          make purge; cd ..; \
	cd NuValidator;        make purge; cd ..; \
	cd PartonModel;        make purge; cd ..; \
	cd Paschos;            make purge; cd ..; \
	cd PDF;                make purge; cd ..; \
	cd PDG;                make purge; cd ..; \
	cd Radiative;          make purge; cd ..; \
	cd Registry;           make purge; cd ..; \
	cd ReinSeghal;         make purge; cd ..; \
	cd Utils;              make purge; cd ..; \
	cd Viewer;             make purge; cd ..; \

clean: FORCE
	rm build.options; \
	cd Algorithm;          make clean; cd ..; \
	cd BaryonResonance;    make clean; cd ..; \
	cd Base;               make clean; cd ..; \
	cd BodekYang;          make clean; cd ..; \
	cd Charm;              make clean; cd ..; \
	cd CrossSections;      make clean; cd ..; \
	cd Decay; 	       make clean; cd ..; \
	cd Elastic; 	       make clean; cd ..; \
	cd EVGCore;            make clean; cd ..; \
	cd EVGModules;         make clean; cd ..; \
	cd EVGDrivers;         make clean; cd ..; \
	cd FluxDrivers;        make clean; cd ..; \
	cd Fragmentation;      make clean; cd ..; \
	cd GHEP;               make clean; cd ..; \
	cd GeantInterface;     make clean; cd ..; \
	cd Geo;                make clean; cd ..; \
	cd ValidityDoc;        make clean; cd ..; \
	cd Interaction;        make clean; cd ..; \
	cd InverseMuonDecay;   make clean; cd ..; \
	cd LlewellynSmith;     make clean; cd ..; \
	cd Messenger;	       make clean; cd ..; \
	cd Nuclear;            make clean; cd ..; \
	cd Ntuple;             make clean; cd ..; \
	cd Numerical;          make clean; cd ..; \
	cd NuValidator;        make clean; cd ..; \
	cd PartonModel;        make clean; cd ..; \
	cd Paschos;            make clean; cd ..; \
	cd PDF;                make clean; cd ..; \
	cd PDG;                make clean; cd ..; \
	cd Radiative;          make clean; cd ..; \
	cd Registry;           make clean; cd ..; \
	cd ReinSeghal;         make clean; cd ..; \
	cd Utils;              make clean; cd ..; \
	cd Viewer;             make clean; cd ..; \
	cd stdapp;             make clean; cd ..; \
	cd test;               make clean; cd ..; \
	cd scripts;	       make clean; cd ..;

FORCE:

