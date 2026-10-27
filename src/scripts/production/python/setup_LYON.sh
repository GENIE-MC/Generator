#!/bin/bash

# Environment for the CC-IN2P3 (Lyon) worker nodes.
# The GENIE dependencies (ROOT, log4cpp, LHAPDF6, Pythia6, GSL, libxml2, Geant4, HepMC3) are taken from an LCG view.
# Edit the copy of this file in your jobs top dir to use a different LCG release or platform.
# The GENIE build given with --genie-topdir must be compiled with this same environment.
export LCG_VIEW=/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt
source $LCG_VIEW/setup.sh

export LHAPDF_DATA_PATH=/cvmfs/sft.cern.ch/lcg/external/lhapdfsets/current:$LHAPDF_DATA_PATH
