#!/bin/bash

# GENIE setup for the CC-IN2P3 (Lyon) worker nodes. Requires setup_LYON.sh to be sourced first.
# Arguments as in setup_GENIE.sh: git location, git branch, INCL, Geant4, config dir, HepMC3

GITHUB_LOCATION=$1
GENIE_VERSION=master
CONFIGURE_INCL=$3
CONFIGURE_G4=$4
GENIE_CONFIG_DIR=$5
CONFIGURE_HEPMC3=$6

if [ -z "$LCG_VIEW" ] ; then
    echo "LCG_VIEW is not defined. Source setup_LYON.sh first"
    return 1
fi

if [ "$CONFIGURE_INCL" = "true" ] ; then
    echo "INCL is not available in the LCG views used at Lyon"
    return 1
fi

GENIEBASE=$(pwd)
if [ -n "$LYON_GENIE_PREBUILT" ] ; then
    echo Using prebuilt GENIE from $LYON_GENIE_PREBUILT
    export GENIE=$LYON_GENIE_PREBUILT
else
    git clone $GITHUB_LOCATION Generator || return 1
    if [ ! -z "$2" ] ; then
        GENIE_VERSION=$2
        git -C Generator checkout "$GENIE_VERSION" || return 1
    fi
    echo Requested Genie version $GENIE_VERSION
    export GENIE=$GENIEBASE/Generator
fi

export GENIE_BIN=$GENIE/bin
export GENIE_LIB=$GENIE/lib
export XSECSPLINEDIR=$GENIEBASE/data
export PATH=$GENIE/bin:$PATH
export LD_LIBRARY_PATH=$GENIE_LIB:$LD_LIBRARY_PATH
unset GENIEBASE

if [ ! -z "$GENIE_CONFIG_DIR" ] ; then
    echo Requested Genie config $GENIE_CONFIG_DIR
    if [ -d "$GENIE_CONFIG_DIR" ]; then
        export GALGCONF=$GENIE_CONFIG_DIR
    elif [ -d "$GENIE/$GENIE_CONFIG_DIR" ]; then
        export GALGCONF=$GENIE/$GENIE_CONFIG_DIR
    fi
    echo Genie config directory $GALGCONF
else
    unset GALGCONF
fi

if [ -n "$LYON_GENIE_PREBUILT" ] ; then
    return 0
fi

cd $GENIE

conf_string="--enable-gsl --enable-rwght --with-optimiz-level=O2 "
conf_string+=" --with-log4cpp-inc=${LCG_VIEW}/include --with-log4cpp-lib=${LCG_VIEW}/lib "
conf_string+=" --with-libxml2-inc=${LCG_VIEW}/include/libxml2 --with-libxml2-lib=${LCG_VIEW}/lib "
conf_string+=" --with-lhapdf6-inc=${LCG_VIEW}/include --with-lhapdf6-lib=${LCG_VIEW}/lib --disable-lhapdf5 --enable-lhapdf6 "
conf_string+=" --with-pythia6-lib=${LCG_VIEW}/lib --enable-pythia6 --disable-pythia8 "

if [ "$CONFIGURE_HEPMC3" == "true" ] ; then
    echo Requested HepMC3 library
    conf_string+=" --enable-hepmc3 "
fi

if [ "$CONFIGURE_G4" = "true" ] ; then
    conf_string+=" --enable-geant4 "
fi

./configure $conf_string || return 1

make -j${SLURM_CPUS_PER_TASK:-1} || return 1
