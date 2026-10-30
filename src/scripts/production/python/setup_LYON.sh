#!/bin/bash

# Environment for GENIE at CC-IN2P3 (Lyon). It is used to compile GENIE (compile_at_Lyon.sh) and by the grid jobs:
# GENIE must run with the same environment it was compiled with.
# ROOT, LHAPDF6 and log4cpp are taken from the JUNO externals (gcc11) in CVMFS.
source /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J24.2.0/ExternalLibs/ROOT/6.30.08/bashrc
source /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J24.1.2/ExternalLibs/LHAPDF/6.3.0/bashrc
source /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J24.1.2/ExternalLibs/log4cpp/1.1.3/bashrc

# Pythia8 is not part of these externals: PYTHIA8PATH is your Pythia8 installation.
# compile_at_Lyon.sh writes it in the copy of this file saved with the GENIE build.
export PYTHIA8PATH=${PYTHIA8PATH:-}
if [ -n "$PYTHIA8PATH" ] ; then
    export PYTHIA8DATA=$PYTHIA8PATH/share/Pythia8/xmldoc
    export LD_LIBRARY_PATH=$PYTHIA8PATH/lib:$LD_LIBRARY_PATH
else
    echo "Warning: PYTHIA8PATH is not set"
fi
