#
# This is the main GENIE project Makefile.
# Each package has its own Makefile which is invoked by this one.
# Machine specific flags and locations are read from 'make/Make.include'.
# Configuration options are read from 'make/Make.config' generated by the 'configure' script.
#
# Costas Andreopoulos <c.andreopoulos \at cern.ch>
#

SHELL = /bin/sh
NAME = all
MAKEFILE = Makefile

# Add empty variable to add flags over command line
# CDBG +=
CFLAGS += -w


# Include machine specific flags and locations (inc. files & libs)
#
include $(GENIE)/src/make/Make.include

# define composite targets
#
INITIAL_BUILD_TARGETS = print-make-info \
		   make-bin-lib-dir \
		   save-build-env \
		   autogenerated-headers \
		   framework \
		   physics-utilities \
		   physics-nuclear-environment \
		   physics-hadronic-simulations \
		   physics-neutrino-scattering-modes \
		   physics-nucleon-decay \
		   physics-nnbar-oscillation \
		   physics-boosted-dark-matter \
		   physics-heavy-neutral-lepton \
		   physics-dark-neutrino \
		   tools-flux-drivers \
		   tools-geometry-drivers \
		   tools-evtlib \
		   tools-masterclass \
			 lib-professor2-build
FINAL_BUILD_TARGETS = doxygen-doc \
		   apps \
		   install-scripts
INSTALL_TARGETS =  print-makeinstall-info \
		   check-previous-installation \
		   make-install-dirs \
		   copy-install-files \
			 lib-professor2-install

# define targets

all:     $(FINAL_BUILD_TARGETS)
install: $(INSTALL_TARGETS)

print-make-info: FORCE
	@echo " "
	@echo " "
	@echo "***** Building GENIE from source tree at: $(GENIE)"
	@echo "***** The source tree corresponds to GENIE version $(GVERSION)"
	@echo " "


print-makeinstall-info: FORCE
	@echo " "
	@echo " "
	@echo "***** Installing GENIE version $(GVERSION) at $(GENIE_INSTALLATION_PATH)"
	@echo " "


framework: FORCE
	@echo " "
	@echo "** Building GENIE framework..."
	cd ${GENIE}/src/Framework && \
	cd Messenger       &&  $(MAKE) && cd .. && \
	cd Algorithm       &&  $(MAKE) && cd .. && \
	cd EventGen        &&  $(MAKE) && cd .. && \
	cd GHEP            &&  $(MAKE) && cd .. && \
	cd Interaction     &&  $(MAKE) && cd .. && \
	cd Ntuple          &&  $(MAKE) && cd .. && \
	cd Numerical       &&  $(MAKE) && cd .. && \
	cd ParticleData    &&  $(MAKE) && cd .. && \
	cd Registry        &&  $(MAKE) && cd .. && \
	cd Utils           &&  $(MAKE) && \
	cd ${GENIE}


physics-neutrino-scattering-modes: FORCE
	@echo " "
	@echo "** Building simulation modules for neutrino scattering modes..."
	cd ${GENIE}/src/Physics/AnomalyMediatedNuGamma/XSection  &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/AnomalyMediatedNuGamma/EventGen  &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/Charm/XSection                   &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/Coherent/XSection                &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/Coherent/EventGen                &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/DeepInelastic/XSection           &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/DeepInelastic/EventGen           &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/Diffractive/XSection             &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/Diffractive/EventGen             &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/HELepton/XSection                &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/HELepton/EventGen                &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/InverseBetaDecay/XSection        &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/InverseBetaDecay/EventGen        &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/Multinucleon/XSection            &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/Multinucleon/EventGen            &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/NuElectron/XSection              &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/NuElectron/EventGen              &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/QuasiElastic/XSection            &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/QuasiElastic/EventGen            &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/Resonance/XSection               &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/Resonance/EventGen               &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/Strange/XSection                 &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/Strange/EventGen                 &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/HEDIS/XSection                   &&  $(MAKE) &&   \
	cd ${GENIE}/src/Physics/HEDIS/EventGen                   &&  $(MAKE) &&   \
	cd ${GENIE}

physics-nucleon-decay:
	@echo " "
	@echo "** Building nucleon decay library..."
ifeq ($(strip $(GOPT_ENABLE_NUCLEON_DECAY)),YES)
	cd ${GENIE}/src/Physics && \
	cd NucleonDecay && \
	$(MAKE) && \
	cd ${GENIE}
else
	@echo " "
	@echo "** Nucleon decay library was not enabled. Skipping..."
endif


physics-nnbar-oscillation:
	@echo " "
	@echo "** Building n-nbar oscillation library..."
ifeq ($(strip $(GOPT_ENABLE_NNBAR_OSCILLATION)),YES)
	cd ${GENIE}/src/Physics && \
	cd NNBarOscillation && \
	$(MAKE) && \
	cd ${GENIE}
else
	@echo " "
	@echo "** N-Nbar oscillation library was not enabled. Skipping..."
endif


physics-boosted-dark-matter:
	@echo " "
	@echo "** Building boosted dark matter library..."
ifeq ($(strip $(GOPT_ENABLE_BOOSTED_DARK_MATTER)),YES)
	cd ${GENIE}/src/Physics/BoostedDarkMatter/EventGen && $(MAKE) && \
	cd ${GENIE}/src/Physics/BoostedDarkMatter/XSection && $(MAKE) && \
	cd ${GENIE}
else
	@echo " "
	@echo "** Boosted dark matter library was not enabled. Skipping..."
endif


physics-heavy-neutral-lepton:
	@echo " "
	@echo "** Building neutral heavy lepton library..."
ifeq ($(strip $(GOPT_ENABLE_HEAVY_NEUTRAL_LEPTON)),YES)
	cd ${GENIE}/src/Physics/BeamHNL && $(MAKE) && \
	cd ${GENIE}
else
	@echo " "
	@echo "** Neutral heavy lepton library was not enabled. Skipping..."
endif


physics-dark-neutrino:
	@echo " "
	@echo "** Building dark neutrino library..."
ifeq ($(strip $(GOPT_ENABLE_DARK_NEUTRINO)),YES)
	cd ${GENIE}/src/Physics/DarkNeutrino/XSection && $(MAKE) && \
	cd ${GENIE}/src/Physics/DarkNeutrino/EventGen && $(MAKE) && \
	cd ${GENIE}
else
	@echo " "
	@echo "** Dark neutrino library was not enabled. Skipping..."
endif


physics-utilities: FORCE
	@echo " "
	@echo "** Building misc physics utility libraries..."
	cd ${GENIE}/src/Physics && \
	cd Common              &&    $(MAKE) && cd .. && \
	cd Decay               &&    $(MAKE) && cd .. && \
	cd HadronTensors       &&    $(MAKE) && cd .. && \
	cd MuonEnergyLoss      &&    $(MAKE) && cd .. && \
	cd PartonDistributions &&    $(MAKE) && cd .. && \
	cd XSectionIntegration &&    $(MAKE) && \
	cd ${GENIE}


physics-nuclear-environment: FORCE
	@echo " "
	@echo "** Building libraries for simulation of nuclear environment..."
	cd ${GENIE}/src/Physics && \
	cd NuclearState         &&    $(MAKE) && cd .. && \
	cd NuclearDeExcitation  &&    $(MAKE) && \
	cd ${GENIE}


physics-hadronic-simulations: FORCE
	@echo " "
	@echo "** Building libraries for hadronic simulations..."
	cd ${GENIE}/src/Physics && \
	cd Hadronization     &&    $(MAKE) && cd .. && \
	cd HadronTransport   &&    $(MAKE) && \
	cd ${GENIE}

#

tools-flux-drivers: FORCE
ifeq ($(strip $(GOPT_ENABLE_FLUX_DRIVERS)),YES)
	@echo " "
	@echo "** Building flux-drivers..."
	cd ${GENIE}/src/Tools/Flux && \
	$(MAKE) && \
	cd ${GENIE}
else
	@echo " "
	@echo "** Building flux-drivers was not enabled. Skipping..."
endif


tools-geometry-drivers: FORCE
ifeq ($(strip $(GOPT_ENABLE_GEOM_DRIVERS)),YES)
	@echo " "
	@echo "** Building geometry-drivers..."
	cd ${GENIE}/src/Tools/Geometry && \
	$(MAKE) && \
	cd ${GENIE}
else
	@echo " "
	@echo "** Building geometry-drivers was not enabled. Skipping..."
endif

lib-professor2-build: FORCE
ifeq ($(strip $(GOPT_ENABLE_PROFESSOR2)),YES)
	@echo " "
	@echo "** Building Professor2..."
	cd ${GENIE}/src/ExternalLibs/professor && \
	$(MAKE) lib/libProfessor2.so && \
	ln -sf ${GENIE}/src/ExternalLibs/professor/lib/libProfessor2.so ${GENIE}/lib/libProfessor2.so && \
	cd ${GENIE}
else
	@echo " "
	@echo "** Building Professor2 was not enabled. Skipping..."
endif

lib-professor2-install: FORCE
ifeq ($(strip $(GOPT_ENABLE_PROFESSOR2)),YES)
	@echo " "
	@echo "** Building Professor2..."
	cp ${GENIE}/src/ExternalLibs/professor/lib/libProfessor2.so ${GENIE_LIB_INSTALLATION_PATH}
	cp -r ${GENIE}/src/ExternalLibs/professor/include/Professor ${GENIE_INC_INSTALLATION_PATH}
else
	@echo " "
	@echo "** Building Professor2 was not enabled. Skipping..."
endif

tools-evtlib: FORCE
ifeq ($(strip $(GOPT_ENABLE_EVTLIB)),YES)
	@echo " "
	@echo "** Building EvtLib..."
	cd ${GENIE}/src/Tools/EvtLib && \
	$(MAKE) && \
	cd ${GENIE}
else
	@echo " "
	@echo "** Building EvtLib was not enabled. Skipping..."
endif


tools-masterclass: FORCE
ifeq ($(strip $(GOPT_ENABLE_MASTERCLASS)),YES)
	@echo " "
	@echo "** Building Masterclass package..."
	cd ${GENIE}/src/Tools/Masterclass && \
	$(MAKE) && \
	cd ${GENIE}
else
	@echo " "
	@echo "** Masterclass was not enabled. Skipping..."
endif

# This target is used for generating the doxygen documentation
# during the genie build.
# It only does so if the option has been enabled explicitly by the user.
doxygen-doc: FORCE
ifeq ($(strip $(GOPT_ENABLE_DOXYGEN_DOC)),YES)
	@echo " "
	@echo "** Building doxygen documentation..."
	cd ${GENIE}/src/scripts && \
	$(MAKE) doxygen && \
	cd ${GENIE};
else
endif


# Use this target  to generate the doxygen documentation at any
# point, independently of your local genie build.
doxygen: FORCE
	@echo " "
	@echo "** Building doxygen documentation..."
	cd ${GENIE}/src/scripts \
	$(MAKE) doxygen && \
	cd ${GENIE}


apps: $(INITIAL_BUILD_TARGETS) FORCE
	@echo " "
	@echo "** Building GENIE applications..."
	cd ${GENIE}/src/Apps && \
	$(MAKE) all && \
	cd ${GENIE}


install-scripts: FORCE
	@echo " "
	@echo "** Installing scripts..."
	cd ${GENIE}/src/scripts && \
	$(MAKE) install && \
	cd ${GENIE}


all-libs: base-framework \
	  utils \
	  physics-models-modules \
	  event-generation-modules \
	  flux-drivers \
	  geom-drivers \
	  mueloss


save-build-env: FORCE
	@echo " "
	@echo "** Taking a snapshot of the build environment..."
	perl ${GENIE}/src/scripts/setup/genie-build-env-snapshot


autogenerated-headers: FORCE
	@echo " "
	@echo "** Adding automatically generated code..."
	perl ${GENIE}/src/scripts/setup/genie-write-gbuild
	perl ${GENIE}/src/scripts/setup/genie-write-gversion


make-bin-lib-dir: FORCE
	@echo " "
	@echo "** Creating GENIE lib and bin directories..."
	cd ${GENIE} && \
	[ -d bin ] || mkdir bin && chmod 755 bin && \
	[ -d lib ] || mkdir lib && chmod 755 lib;


check-previous-installation: FORCE
	@echo " "
	@echo "** Testing for existing GENIE installation at specified installation location..."
ifeq ($(strip $(GENIE_PREVIOUS_INSTALLATION)),YES)
	$(error Previous installation exists at your specified installation path: $(GENIE_INSTALLATION_PATH). Try '$(MAKE) distclean' first)
endif


make-install-dirs: FORCE
	@echo " "
	@echo "** Creating directory structure for GENIE installation..."
ifeq (${GENIE_INSTALLATION_PATH},/usr)
	cd ${GENIE}
	$(warning Installation in /usr or /usr/local is discouraged)
endif
ifeq (${GENIE_INSTALLATION_PATH},/usr/local)
	cd ${GENIE}
	$(warning Installation in /usr or /usr/local is discouraged)
endif
	[ -d ${GENIE_INSTALLATION_PATH} ] || mkdir ${GENIE_INSTALLATION_PATH}
	cd ${GENIE_INSTALLATION_PATH}
	[ -d ${GENIE_BIN_INSTALLATION_PATH}     ] || mkdir ${GENIE_BIN_INSTALLATION_PATH}
	[ -d ${GENIE_LIB_INSTALLATION_PATH}     ] || mkdir ${GENIE_LIB_INSTALLATION_PATH}
	[ -d ${GENIE_INCBASE_INSTALLATION_PATH} ] || mkdir ${GENIE_INCBASE_INSTALLATION_PATH}
	mkdir ${GENIE_INC_INSTALLATION_PATH}
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework/Algorithm
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework/Conventions
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework/GHEP
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework/Interaction
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework/Messenger
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework/Ntuple
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework/Numerical
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework/ParticleData
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework/Registry
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Framework/Utils
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/AnomalyMediatedNuGamma
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/AnomalyMediatedNuGamma/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/AnomalyMediatedNuGamma/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/BoostedDarkMatter
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/BoostedDarkMatter/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/BoostedDarkMatter/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Charm
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Charm/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Coherent
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Coherent/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Coherent/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Common
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/DarkNeutrino
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/DarkNeutrino/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/DarkNeutrino/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Decay
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/DeepInelastic
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/DeepInelastic/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/DeepInelastic/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Diffractive
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Diffractive/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Diffractive/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/HELepton
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/HELepton/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/HELepton/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Hadronization
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/HadronTensors
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/HadronTransport
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/InverseBetaDecay
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/InverseBetaDecay/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/InverseBetaDecay/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Multinucleon
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Multinucleon/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Multinucleon/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/MuonEnergyLoss
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/NNBarOscillation
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/BeamHNL
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/NuclearState
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/NuclearDeExcitation
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/NucleonDecay
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/NuElectron
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/NuElectron/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/NuElectron/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/PartonDistributions
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/QuasiElastic
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/QuasiElastic/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/QuasiElastic/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Resonance
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Resonance/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Resonance/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Strange
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Strange/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/Strange/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/HEDIS
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/HEDIS/XSection
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/HEDIS/EventGen
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Physics/XSectionIntegration
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Tools
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Tools/Flux
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Tools/Geometry
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Tools/Masterclass
	mkdir ${GENIE_INC_INSTALLATION_PATH}/Tools/EvtLib


copy-install-files: FORCE
	@echo " "
	@echo "** Copying libraries/binaries/headers to installation location..."
	cp ${GENIE_BIN_PATH}/* ${GENIE_BIN_INSTALLATION_PATH} && \
	cd ${GENIE}/src/Framework/Algorithm                      &&  $(MAKE) install && \
	cd ${GENIE}/src/Framework/Conventions                    &&  $(MAKE) install && \
	cd ${GENIE}/src/Framework/EventGen                       &&  $(MAKE) install && \
	cd ${GENIE}/src/Framework/GHEP                           &&  $(MAKE) install && \
	cd ${GENIE}/src/Framework/Interaction                    &&  $(MAKE) install && \
	cd ${GENIE}/src/Framework/Messenger                      &&  $(MAKE) install && \
	cd ${GENIE}/src/Framework/Ntuple                         &&  $(MAKE) install && \
	cd ${GENIE}/src/Framework/Numerical                      &&  $(MAKE) install && \
	cd ${GENIE}/src/Framework/ParticleData                   &&  $(MAKE) install && \
	cd ${GENIE}/src/Framework/Registry                       &&  $(MAKE) install && \
	cd ${GENIE}/src/Framework/Utils                          &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/AnomalyMediatedNuGamma/XSection  &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/AnomalyMediatedNuGamma/EventGen  &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/BoostedDarkMatter/XSection       &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/BoostedDarkMatter/EventGen       &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Charm/XSection                   &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Coherent/XSection                &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Coherent/EventGen                &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Common                           &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/DarkNeutrino/XSection            &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/DarkNeutrino/EventGen            &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Decay                            &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/DeepInelastic/XSection           &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/DeepInelastic/EventGen           &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Diffractive/XSection             &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Diffractive/EventGen             &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/HELepton/XSection                &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/HELepton/EventGen                &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Hadronization                    &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/HadronTransport                  &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/HadronTensors                    &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/InverseBetaDecay/XSection        &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/InverseBetaDecay/EventGen        &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Multinucleon/XSection            &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Multinucleon/EventGen            &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/MuonEnergyLoss                   &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/NNBarOscillation                 &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/BeamHNL               &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/NuclearState                     &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/NuclearDeExcitation              &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/NucleonDecay                     &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/NuElectron/XSection              &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/NuElectron/EventGen              &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/PartonDistributions              &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/QuasiElastic/XSection            &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/QuasiElastic/EventGen            &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Resonance/XSection               &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Resonance/EventGen               &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Strange/XSection                 &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/Strange/EventGen                 &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/HEDIS/XSection                   &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/HEDIS/EventGen                   &&  $(MAKE) install && \
	cd ${GENIE}/src/Physics/XSectionIntegration              &&  $(MAKE) install && \
	cd ${GENIE}/src/Tools/Flux                               &&  $(MAKE) install && \
	cd ${GENIE}/src/Tools/Geometry                           &&  $(MAKE) install && \
	cd ${GENIE}/src/Tools/Masterclass                        &&  $(MAKE) install && \
	cd ${GENIE}/src/Tools/EvtLib                             &&  $(MAKE) install && \
	cd ${GENIE}


purge: FORCE
	@echo " "
	@echo "** Purging..."
	cd ${GENIE}/src/Framework/Algorithm                      &&  $(MAKE) purge && \
	cd ${GENIE}/src/Framework/EventGen                       &&  $(MAKE) purge && \
	cd ${GENIE}/src/Framework/GHEP                           &&  $(MAKE) purge && \
	cd ${GENIE}/src/Framework/Interaction                    &&  $(MAKE) purge && \
	cd ${GENIE}/src/Framework/Messenger                      &&  $(MAKE) purge && \
	cd ${GENIE}/src/Framework/Ntuple                         &&  $(MAKE) purge && \
	cd ${GENIE}/src/Framework/Numerical                      &&  $(MAKE) purge && \
	cd ${GENIE}/src/Framework/ParticleData                   &&  $(MAKE) purge && \
	cd ${GENIE}/src/Framework/Registry                       &&  $(MAKE) purge && \
	cd ${GENIE}/src/Framework/Utils                          &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/AnomalyMediatedNuGamma/XSection  &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/AnomalyMediatedNuGamma/EventGen  &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/BoostedDarkMatter/XSection       &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/BoostedDarkMatter/EventGen       &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Charm/XSection                   &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Coherent/XSection                &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Coherent/EventGen                &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Common                           &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/DarkNeutrino/XSection            &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/DarkNeutrino/EventGen            &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Decay                            &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/DeepInelastic/XSection           &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/DeepInelastic/EventGen           &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Diffractive/XSection             &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Diffractive/EventGen             &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/HELepton/XSection                &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/HELepton/EventGen                &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Hadronization                    &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/HadronTensors                    &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/HadronTransport                  &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/InverseBetaDecay/XSection        &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/InverseBetaDecay/EventGen        &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Multinucleon/XSection            &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Multinucleon/EventGen            &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/MuonEnergyLoss                   &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/NNBarOscillation                 &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/BeamHNL               &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/NuclearState                     &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/NuclearDeExcitation              &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/NucleonDecay                     &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/NuElectron/XSection              &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/NuElectron/EventGen              &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/PartonDistributions              &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/QuasiElastic/XSection            &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/QuasiElastic/EventGen            &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Resonance/XSection               &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Resonance/EventGen               &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Strange/XSection                 &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/Strange/EventGen                 &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/HEDIS/XSection                   &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/HEDIS/EventGen                   &&  $(MAKE) purge && \
	cd ${GENIE}/src/Physics/XSectionIntegration              &&  $(MAKE) purge && \
	cd ${GENIE}/src/Tools/Flux                               &&  $(MAKE) purge && \
	cd ${GENIE}/src/Tools/Geometry                           &&  $(MAKE) purge && \
	cd ${GENIE}/src/Tools/Masterclass                        &&  $(MAKE) purge && \
	cd ${GENIE}/src/Tools/EvtLib                             &&  $(MAKE) purge && \
	cd ${GENIE}

clean: clean-files clean-dir clean-etc

clean-files: FORCE
	@echo " "
	@echo "** Cleaning..."
	cd ${GENIE}/src/Framework/Algorithm                      &&  $(MAKE) clean && \
	cd ${GENIE}/src/Framework/EventGen                       &&  $(MAKE) clean && \
	cd ${GENIE}/src/Framework/GHEP                           &&  $(MAKE) clean && \
	cd ${GENIE}/src/Framework/Interaction                    &&  $(MAKE) clean && \
	cd ${GENIE}/src/Framework/Messenger                      &&  $(MAKE) clean && \
	cd ${GENIE}/src/Framework/Ntuple                         &&  $(MAKE) clean && \
	cd ${GENIE}/src/Framework/Numerical                      &&  $(MAKE) clean && \
	cd ${GENIE}/src/Framework/ParticleData                   &&  $(MAKE) clean && \
	cd ${GENIE}/src/Framework/Registry                       &&  $(MAKE) clean && \
	cd ${GENIE}/src/Framework/Utils                          &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/AnomalyMediatedNuGamma/XSection  &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/AnomalyMediatedNuGamma/EventGen  &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/BoostedDarkMatter/XSection       &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/BoostedDarkMatter/EventGen       &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Charm/XSection                   &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Coherent/XSection                &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Coherent/EventGen                &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Common                           &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/DarkNeutrino/XSection            &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/DarkNeutrino/EventGen            &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Decay                            &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/DeepInelastic/XSection           &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/DeepInelastic/EventGen           &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Diffractive/XSection             &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Diffractive/EventGen             &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/HELepton/XSection                &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/HELepton/EventGen                &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Hadronization                    &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/HadronTensors                    &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/HadronTransport                  &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/InverseBetaDecay/XSection        &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/InverseBetaDecay/EventGen        &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Multinucleon/XSection            &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Multinucleon/EventGen            &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/MuonEnergyLoss                   &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/NNBarOscillation                 &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/BeamHNL               &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/NuclearState                     &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/NuclearDeExcitation              &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/NucleonDecay                     &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/NuElectron/XSection              &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/NuElectron/EventGen              &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/PartonDistributions              &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/QuasiElastic/XSection            &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/QuasiElastic/EventGen            &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Resonance/XSection               &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Resonance/EventGen               &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Strange/XSection                 &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/Strange/EventGen                 &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/HEDIS/XSection                   &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/HEDIS/EventGen                   &&  $(MAKE) clean && \
	cd ${GENIE}/src/Physics/XSectionIntegration              &&  $(MAKE) clean && \
	cd ${GENIE}/src/Tools/Flux                               &&  $(MAKE) clean && \
	cd ${GENIE}/src/Tools/Geometry                           &&  $(MAKE) clean && \
	cd ${GENIE}/src/Tools/Masterclass                        &&  $(MAKE) clean && \
	cd ${GENIE}/src/Tools/EvtLib                             &&  $(MAKE) clean && \
	cd ${GENIE}/src/Apps                                     &&  $(MAKE) clean && \
	cd ${GENIE}/src/scripts                                  &&  $(MAKE) clean && \
	cd ${GENIE}

clean-dir: FORCE
	@echo "Deleting GENIE lib and bin directories..." && \
	cd $(GENIE) && \
	[ ! -d ./bin ] || rmdir ./bin && \
	[ ! -d ./lib ] || rmdir ./lib

clean-etc: FORCE
	cd $(GENIE) && \
	rm -f ./*log && \
	cd ${GENIE}

distclean: FORCE
	@echo " "
	@echo "** Cleaning GENIE installation... "
ifeq (${GENIE_INSTALLATION_PATH},/usr)
	cd ${GENIE}
	$(warning Installation in /usr or /usr/local is discouraged)
	$(warning   as "distclean" is too aggressive and removes more than just GENIE libraries)
	$(error Installation directory ${GENIE_INSTALLATION_PATH} is too precious.  Must be removed by hand)
endif
ifeq (${GENIE_INSTALLATION_PATH},/usr/local)
	cd ${GENIE}
	$(warning Installation in /usr or /usr/local is discouraged)
	$(warning   as "distclean" is too aggressive and removes more than just GENIE libraries)
	$(error Installation directory ${GENIE_INSTALLATION_PATH} is too precious.  Must be removed by hand)
endif
	[ ! -d ${GENIE_INSTALLATION_PATH}/include/GENIE ] || rm -rf ${GENIE_INSTALLATION_PATH}/include/GENIE/
	cd ${GENIE}/src/Framework/Algorithm                      &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Framework/Conventions                    &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Framework/EventGen                       &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Framework/GHEP                           &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Framework/Interaction                    &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Framework/Messenger                      &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Framework/Ntuple                         &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Framework/Numerical                      &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Framework/ParticleData                   &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Framework/Registry                       &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Framework/Utils                          &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/AnomalyMediatedNuGamma/XSection  &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/AnomalyMediatedNuGamma/EventGen  &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/BoostedDarkMatter/XSection       &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/BoostedDarkMatter/EventGen       &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Charm/XSection                   &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Coherent/XSection                &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Coherent/EventGen                &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Common                           &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/DarkNeutrino/XSection            &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/DarkNeutrino/EventGen            &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Decay                            &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/DeepInelastic/XSection           &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/DeepInelastic/EventGen           &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Diffractive/XSection             &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Diffractive/EventGen             &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/HELepton/XSection                &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/HELepton/EventGen                &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Hadronization                    &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/HadronTransport                  &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/HadronTensors                    &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/InverseBetaDecay/XSection        &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/InverseBetaDecay/EventGen        &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Multinucleon/XSection            &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Multinucleon/EventGen            &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/MuonEnergyLoss                   &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/NNBarOscillation                 &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/BeamHNL                          &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/NuclearState                     &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/NuclearDeExcitation              &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/NucleonDecay                     &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/NuElectron/XSection              &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/NuElectron/EventGen              &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/PartonDistributions              &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/QuasiElastic/XSection            &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/QuasiElastic/EventGen            &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Resonance/XSection               &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Resonance/EventGen               &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Strange/XSection                 &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/Strange/EventGen                 &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/HEDIS/XSection                   &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/HEDIS/EventGen                   &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Physics/XSectionIntegration              &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Tools/Flux                               &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Tools/Geometry                           &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Tools/Masterclass                        &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Tools/EvtLib                             &&  $(MAKE) distclean && \
	cd ${GENIE}/src/Apps                                     &&  $(MAKE) distclean && \
	cd ${GENIE}/src/scripts                                  &&  $(MAKE) distclean && \
	cd ${GENIE}


FORCE:
