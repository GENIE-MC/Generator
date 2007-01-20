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
BUILD_TARGETS = make-bin-lib-dir \
		save-build-options \
		base-framework \
		utils \
		physics-models-modules \
		event-generation-modules \
		flux-drivers \
		geom-drivers \
		viewer \
		mueloss \
		dummy-neugen \
		nuvld-libs \
		neugen \
		nuvld-exe \
		doxygen-doc \
		generator-test-exe \
		generator-std-exe \
		install-scripts

# define build targets

all: $(BUILD_TARGETS)

base-framework: FORCE
	@echo " "
	@echo "** Building base-framework..."
	cd ${GENIE}/src;\
	cd Algorithm;  make; cd ..; \
	cd Messenger;  make; cd ..; \
	cd Registry;   make; cd ..; 

utils: FORCE
	@echo " "
	@echo "** Building utility libraries..."
	cd Numerical;  make; cd ..; 
	cd Utils;      make; cd ..; 

physics-models-modules: FORCE
	@echo " "
	@echo "** Building physics-models-modules..."
	cd ${GENIE}/src;\
	cd BaryonResonance;    make; cd ..; \
	cd Base;               make; cd ..; \
	cd BodekYang;          make; cd ..; \
	cd Charm;              make; cd ..; \
	cd CrossSections;      make; cd ..; \
	cd Decay; 	       make; cd ..; \
	cd Elastic; 	       make; cd ..; \
	cd Fragmentation;      make; cd ..; \
	cd Interaction;        make; cd ..; \
	cd InverseMuonDecay;   make; cd ..; \
	cd LlewellynSmith;     make; cd ..; \
	cd Nuclear;            make; cd ..; \
	cd PartonModel;        make; cd ..; \
	cd Paschos; 	       make; cd ..; \
	cd PDF;                make; cd ..; \
	cd PDG;                make; cd ..; \
	cd Radiative;          make; cd ..; \
	cd ReinSeghal;         make; cd ..; \
	cd VHE;                make; cd ..; 

event-generation-modules: FORCE
	@echo " "
	@echo "** Building event-generation-modules..."
	cd ${GENIE}/src;\
	cd EVGCore;            make; cd ..; \
	cd EVGModules;         make; cd ..; \
	cd EVGDrivers;         make; cd ..; \
	cd HadronTransport;    make; cd ..; \
	cd GHEP;               make; cd ..; \
	cd Ntuple;             make; cd ..; \
	cd ReWeight;           make; cd ..; 

flux-drivers: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	@echo " "
	@echo "** Building flux-drivers..."
	cd ${GENIE}/src;\
	cd FluxDrivers; make; cd ..;
else
	@echo " "
	@echo "** Building flux-drivers was not enabled. Skipping..."
endif

geom-drivers: FORCE
ifeq ($(strip $(GOPT_ENABLE_GEOM_DRIVERS)),YES)
	@echo " "
	@echo "** Building geometry-drivers..."
	cd ${GENIE}/src;\
	cd Geo; make; cd ..; 
else
	@echo " "
	@echo "** Building geometry-drivers was not enabled. Skipping..."
endif

viewer: FORCE
ifeq ($(strip $(GOPT_ENABLE_VIEWER)),YES)
	@echo " "
	@echo "** Buidling viewer..."
	cd ${GENIE}/src;\
	cd Viewer; make; cd ..; 
else
	@echo " "
	@echo "** Viewer was not enabled. Skipping..."
endif

mueloss: FORCE
ifeq ($(strip $(GOPT_ENABLE_MUELOSS)),YES)
	@echo " "
	@echo "** Building mueloss utility package..."
	cd ${GENIE}/src;\
	cd MuELoss; make; cd ..; 
else
	@echo "** Mueloss was not enabled. Skipping..."
endif

dummy-neugen: FORCE
ifeq ($(strip $(GOPT_ENABLE_NEUGEN)),NO)
	@echo " "
	@echo "** Neugen was not enabled. Building dummy-neugen..."
	cd ${GENIE}/src;\
	cd NuValidator; make dummy-neugen; cd ..; 
else
	@echo " "
	@echo "** Neugen was enabled. Will not build dummy version..."
endif

neugen: FORCE
ifeq ($(strip $(GOPT_ENABLE_NUVALIDATOR)),YES)
	@echo " "
	@echo "** Building nuvalidator's neugen interface..."
	cd ${GENIE}/src;\
	cd NuValidator; make neugen; cd ..; 
else
endif

nuvld-libs: FORCE
ifeq ($(strip $(GOPT_ENABLE_NUVALIDATOR)),YES)
	@echo " "
	@echo "** Building nuvalidator libraries..."
	cd ${GENIE}/src;\
	cd NuValidator; make libs; cd ..; 
else
	@echo " "
	@echo "** Nuvalidator was not enabled. Skipping..."
endif

nuvld-exe: FORCE
ifeq ($(strip $(GOPT_ENABLE_NUVALIDATOR)),YES)
	@echo " "
	@echo "** Building nuvalidator executables..."
	cd ${GENIE}/src;\
	cd NuValidator; make exe; cd ..; 
else
endif

doxygen-doc: FORCE
ifeq ($(strip $(GOPT_ENABLE_DOXYGEN_DOC)),YES)
	@echo " "
	@echo "** Building doxygen documentation..."
	cd ${GENIE}/src;\
	cd ./scripts; make doxygen; cd ..;
else
endif

generator-std-exe: FORCE
	@echo " "
	@echo "** Building GENIE applications..."
	cd ${GENIE}/src;\
	cd stdapp; make all; cd ..; 

generator-test-exe: FORCE
ifeq ($(strip $(GOPT_ENABLE_TEST)),YES)
	@echo " "
	@echo "** Building test applications..."
	cd ${GENIE}/src;\
	cd test; make all; cd ..; 
else
endif

install-scripts: FORCE
	@echo " "
	@echo "** Installing scripts..."
	cd ${GENIE}/src;\
	cd ./scripts; make install; cd ..;

all-libs: base-framework \
	  utils \
	  physics-models-modules \
	  event-generation-modules \
	  flux-drivers \
	  geom-drivers \
	  viewer \
	  mueloss \
	  nuvld-libs

save-build-options: FORCE
	@echo " "
	@echo "** Logging configuration options..."
	cd ${GENIE}/src;\
	perl scripts/setup/genie-save-build-env;

make-bin-lib-dir: FORCE
	@echo " "
	@echo "** Creating GENIE lib and bin directories..."
	cd $(GENIE); \
	mkdir bin; chmod 755 bin;\
	mkdir lib; chmod 755 lib;\
	cd ${GENIE}/src;

purge: FORCE
	@echo " "
	@echo "** Purging..."
	cd ${GENIE}/src;\
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
	cd Geo;                make purge; cd ..; \
	cd HadronTransport;    make purge; cd ..; \
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
	cd ReWeight;           make purge; cd ..; \
	cd Utils;              make purge; cd ..; \
	cd Viewer;             make purge; cd ..; \

clean: clean-files clean-dir

clean-files: FORCE
	@echo " "
	@echo "** Cleaning..."
	cd ${GENIE}/src;\
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
	cd Geo;                make clean; cd ..; \
	cd HadronTransport;    make clean; cd ..; \
	cd Interaction;        make clean; cd ..; \
	cd InverseMuonDecay;   make clean; cd ..; \
	cd LlewellynSmith;     make clean; cd ..; \
	cd Messenger;	       make clean; cd ..; \
	cd MuELoss;	       make clean; cd ..; \
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
	cd ReWeight;           make clean; cd ..; \
	cd Utils;              make clean; cd ..; \
	cd Viewer;             make clean; cd ..; \
	cd VHE;                make clean; cd ..; \
	cd stdapp;             make clean; cd ..; \
	cd test;               make clean; cd ..; \
	cd scripts;	       make clean; cd ..;

clean-dir: FORCE
	@echo "Deleting GENIE lib and bin directories...";\
	cd $(GENIE);\
	rmdir ./bin;\
	rmdir ./lib 

FORCE:

