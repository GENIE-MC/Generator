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

all:
	cd Algorithm;          make; cd ..; \
	cd BaryonResonance;    make; cd ..; \
	cd Base;               make; cd ..; \
	cd BodekYang;          make; cd ..; \
	cd Charm;              make; cd ..; \
	cd CrossSections;      make; cd ..; \
	cd Decay; 	       make; cd ..; \
	cd Elastic; 	       make; cd ..; \
	cd EventGeneration;    make; cd ..; \
	cd FKR;                make; cd ..; \
	cd Fragmentation;      make; cd ..; \
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
	cd NuValidator;        make; cd ..; \
	cd test;               make all-no-test; cd ..; \
	cd scripts;	       make install; cd ..;

all-no-nuvld:
	cd Algorithm;          make; cd ..; \
	cd BaryonResonance;    make; cd ..; \
	cd Base;               make; cd ..; \
	cd BodekYang;          make; cd ..; \
	cd Charm;              make; cd ..; \
	cd CrossSections;      make; cd ..; \
	cd Decay; 	       make; cd ..; \
	cd Elastic; 	       make; cd ..; \
	cd EventGeneration;    make; cd ..; \
	cd FKR;                make; cd ..; \
	cd Fragmentation;      make; cd ..; \
	cd ValidityDoc;        make; cd ..; \
	cd Interaction;        make; cd ..; \
	cd InverseMuonDecay;   make; cd ..; \
	cd LlewellynSmith;     make; cd ..; \
	cd Messenger;	       make; cd ..; \
	cd Ntuple;             make; cd ..; \
	cd Nuclear;            make; cd ..; \
	cd Numerical;          make; cd ..; \
	cd PartonModel;        make; cd ..; \
	cd Paschos;	       make; cd ..; \
	cd PDF;                make; cd ..; \
	cd PDG;                make; cd ..; \
	cd Radiative;          make; cd ..; \
	cd Registry;           make; cd ..; \
	cd ReinSeghal;         make; cd ..; \
	cd Utils;              make; cd ..; \
	cd Viewer;             make; cd ..; \
	cd test;               make all-no-test; cd ..; \
	cd scripts;	       make install; cd ..;

nuvld-only:
	cd NuValidator;        make; cd ..; 

doc:
	cd ../doc/doxygen/; doxygen genie_doxygen; cd ../../src/

clean: FORCE
	cd Algorithm;          make clean; cd ..; \
	cd BaryonResonance;    make clean; cd ..; \
	cd Base;               make clean; cd ..; \
	cd BodekYang;          make clean; cd ..; \
	cd Charm;              make clean; cd ..; \
	cd CrossSections;      make clean; cd ..; \
	cd Decay; 	       make clean; cd ..; \
	cd Elastic; 	       make clean; cd ..; \
	cd EventGeneration;    make clean; cd ..; \
	cd FKR;                make clean; cd ..; \
	cd Fragmentation;      make clean; cd ..; \
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
	cd test;               make clean; cd ..; \
	cd scripts;	       make clean; cd ..;

FORCE:

