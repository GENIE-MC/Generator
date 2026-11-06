#!/bin/bash

# Compiles GENIE at CC-IN2P3 (Lyon), to be used by the Lyon grid jobs (--grid-system LYON --genie-topdir <genie dir>).
# This is a standalone script: it only needs to be run once, if you do not have a GENIE build yet.
#
# The FSI model is taken from the tune, Gxx_yyX_zz_www, where X is: a (hA) or b (hN), which need no extra
# libraries, c (INCL++) or d (Geant4). Use --fsi to set it when the tune does not follow this convention.
#
# The environment (setup_LYON.sh by default) is sourced before compiling and a copy of it, including the Pythia8
# path, is saved as <genie dir>/setup_LYON.sh. The grid submitters use that copy for the jobs, so that GENIE
# runs with the same environment it was compiled with.
#
# Author:
#       Julia Tena Vidal <jtenavidal \st tauex.tau.ac.il>
#       Tel Aviv University
# Copyright:
#    Copyright (c) 2003-2026, The GENIE Collaboration
#    For the full text of the license visit http://copyright.genie-mc.org

SCRIPTS_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

GIT_LOCATION=https://github.com/GENIE-MC/Generator
VERSION=master
TUNE=G18_02a_02_11b
FSI=""
TOPDIR=$(pwd)
GENIE_DIR=""
ENV_SCRIPT=$SCRIPTS_DIR/setup_LYON.sh
PYTHIA8_DIR=""
INCL_DIR=""
BOOST_DIR=""
GEANT4_SETUP=""
PROFESSOR2=false
EXTRA_OPTS=""
NJOBS=4
REBUILD=false
DRY_RUN=false

usage() {
    cat <<EOF
Usage: $(basename $0) [options]

  --git-location <url>      GENIE Generator repository. Default: $GIT_LOCATION
  --version <name>          Branch or tag to compile (i.e. master, R-3_06_02). Default: $VERSION
  --tune <tune>             Tune to be used with this build. Sets the FSI model. Default: $TUNE
  --fsi <hA|hN|INCL|G4>     FSI model, if it cannot be taken from the tune name
  --topdir <dir>            Directory where GENIE is cloned. Default: current directory
  --genie-dir <dir>         GENIE directory. Default: <topdir>/Generator_<version>[_INCL|_G4]
  --env-script <file>       Environment to compile and run GENIE. Default: $ENV_SCRIPT
  --pythia8-dir <dir>       Pythia8 installation. Default: \$PYTHIA8PATH
  --incl-dir <dir>          INCL++ installation (include/inclxx, lib). Needed for INCL tunes
  --boost-dir <dir>         Boost installation (include, lib). Needed for INCL tunes
  --geant4-setup <file>     Script setting up Geant4 (geant4-config). Needed for Geant4 tunes
  --enable-professor2       Also build Professor2 (needs cython)
  --configure-options <s>   Extra options passed to GENIE configure
  -j, --jobs <n>            Number of parallel make jobs. Default: $NJOBS
  --rebuild                 Reconfigure and recompile an existing GENIE directory
  --dry-run                 Print the commands without running them
  -h, --help                Print this message
EOF
}

die() {
    echo "*** Error *** $1"
    exit 1
}

# Errors on the environment are only warnings in dry-run mode
check() {
    if $DRY_RUN ; then
        echo "*** Warning *** $1"
    else
        die "$1"
    fi
}

run() {
    echo "> $*"
    if ! $DRY_RUN ; then
        "$@" || die "Command failed: $*"
    fi
}

while [ $# -gt 0 ] ; do
    case $1 in
        --enable-professor2) PROFESSOR2=true ; shift ; continue ;;
        --rebuild)           REBUILD=true ; shift ; continue ;;
        --dry-run)           DRY_RUN=true ; shift ; continue ;;
        -h|--help)           usage ; exit 0 ;;
    esac
    [ $# -ge 2 ] || { usage ; die "Missing value for option $1" ; }
    case $1 in
        --git-location)      GIT_LOCATION=$2 ;;
        --version)           VERSION=$2 ;;
        --tune)              TUNE=$2 ;;
        --fsi)               FSI=$2 ;;
        --topdir)            TOPDIR=$2 ;;
        --genie-dir)         GENIE_DIR=$2 ;;
        --env-script)        ENV_SCRIPT=$2 ;;
        --pythia8-dir)       PYTHIA8_DIR=$2 ;;
        --incl-dir)          INCL_DIR=$2 ;;
        --boost-dir)         BOOST_DIR=$2 ;;
        --geant4-setup)      GEANT4_SETUP=$2 ;;
        --configure-options) EXTRA_OPTS=$2 ;;
        -j|--jobs)           NJOBS=$2 ;;
        *) usage ; die "Unknown option $1" ;;
    esac
    shift 2
done

# FSI model from the tune name
if [ -z "$FSI" ] ; then
    TUNE_CMC=$(echo $TUNE | cut -d_ -f2)
    case ${TUNE_CMC: -1} in
        a) FSI=hA ;;
        b) FSI=hN ;;
        c) FSI=INCL ;;
        d) FSI=G4 ;;
        *) die "Cannot get the FSI model from tune $TUNE. Specify it with --fsi" ;;
    esac
fi
case $FSI in
    hA|hN) BUILD_SUFFIX="" ;;
    INCL)  BUILD_SUFFIX="_INCL" ;;
    G4)    BUILD_SUFFIX="_G4" ;;
    *)     die "Unknown FSI model $FSI. Use hA, hN, INCL or G4" ;;
esac

# hA and hN only differ at run time, so they share the same build
if [ -z "$GENIE_DIR" ] ; then
    GENIE_DIR=$TOPDIR/Generator_${VERSION//\//_}$BUILD_SUFFIX
fi
case $GENIE_DIR in
    /*) ;;
    *) GENIE_DIR=$(pwd)/$GENIE_DIR ;;
esac

echo "GENIE version  : $VERSION from $GIT_LOCATION"
echo "Tune           : $TUNE (FSI: $FSI)"
echo "GENIE directory: $GENIE_DIR"
echo "Environment    : $ENV_SCRIPT"
echo ""

if [ -x "$GENIE_DIR/bin/gevgen" ] && ! $REBUILD ; then
    echo "GENIE is already compiled in $GENIE_DIR. Use --rebuild to compile it again."
    echo "Use it with: --grid-system LYON --genie-topdir $GENIE_DIR --tune $TUNE"
    exit 0
fi

# Environment
[ -f "$ENV_SCRIPT" ] || die "Environment script $ENV_SCRIPT does not exist"
if [ -n "$PYTHIA8_DIR" ] ; then
    export PYTHIA8PATH=$PYTHIA8_DIR
fi
source $ENV_SCRIPT
PYTHIA8_DIR=$PYTHIA8PATH

[ -n "$PYTHIA8_DIR" ] || check "No Pythia8 installation given. Use --pythia8-dir or set PYTHIA8PATH"
[ -f "$PYTHIA8_DIR/include/Pythia8/Pythia.h" ] || check "Pythia8 headers not found in $PYTHIA8_DIR/include"
ls $PYTHIA8_DIR/lib/libpythia8.* > /dev/null 2>&1 || check "Pythia8 library not found in $PYTHIA8_DIR/lib"
for var in LOG4CPP_INC LOG4CPP_LIB LHAPDF_INC LHAPDF_LIB ; do
    [ -d "${!var}" ] || check "$var is not set by $ENV_SCRIPT, or it does not exist: '${!var}'"
done
which root-config > /dev/null 2>&1 || check "root-config not found. ROOT is not set up by $ENV_SCRIPT"
LIBXML2_INC=${LIBXML2_INC:-/usr/include/libxml2}
LIBXML2_LIB=${LIBXML2_LIB:-/usr/lib64}

CONF_OPTS=( --enable-gsl --enable-rwght --with-optimiz-level=O2
            --with-log4cpp-inc=$LOG4CPP_INC --with-log4cpp-lib=$LOG4CPP_LIB
            --with-libxml2-inc=$LIBXML2_INC --with-libxml2-lib=$LIBXML2_LIB
            --enable-lhapdf6 --disable-lhapdf5 --with-lhapdf6-inc=$LHAPDF_INC --with-lhapdf6-lib=$LHAPDF_LIB
            --disable-pythia6 --enable-pythia8 --with-pythia8-inc=$PYTHIA8_DIR/include --with-pythia8-lib=$PYTHIA8_DIR/lib
            --enable-flux-drivers --enable-geom-drivers --disable-doxygen --enable-test --enable-mueloss
            --enable-dylibversion --enable-t2k --enable-fnal --enable-atmo --enable-nucleon-decay
            --disable-masterclass --disable-debug --disable-profiler --enable-validation-tools --enable-gfortran )

if $PROFESSOR2 ; then
    which cython > /dev/null 2>&1 || check "cython not found, needed by Professor2"
    CONF_OPTS+=( --enable-professor2 )
fi

if [ "$FSI" = "INCL" ] ; then
    [ -d "$INCL_DIR/include/inclxx" ] || check "INCL tunes need INCL++: give its installation with --incl-dir"
    [ -d "$BOOST_DIR/include/boost" ] || check "INCL tunes need Boost: give its installation with --boost-dir"
    CONF_OPTS+=( --enable-incl --with-incl-inc=$INCL_DIR/include/inclxx --with-incl-lib=$INCL_DIR/lib
                 --with-boost-inc=$BOOST_DIR/include --with-boost-lib=$BOOST_DIR/lib )
elif [ "$FSI" = "G4" ] ; then
    if [ -n "$GEANT4_SETUP" ] ; then
        [ -f "$GEANT4_SETUP" ] || die "Geant4 setup $GEANT4_SETUP does not exist"
        source $GEANT4_SETUP
    fi
    which geant4-config > /dev/null 2>&1 || check "Geant4 tunes need Geant4: give its setup script with --geant4-setup"
    CONF_OPTS+=( --enable-geant4 )
fi

CONF_OPTS+=( $EXTRA_OPTS )

# Get the code
if [ -d "$GENIE_DIR" ] ; then
    echo "Using the existing code in $GENIE_DIR"
else
    run git clone $GIT_LOCATION $GENIE_DIR
    run git -C $GENIE_DIR checkout $VERSION
fi

# Compile
export GENIE=$GENIE_DIR
run cd $GENIE_DIR
if $REBUILD && [ -f Makefile ] ; then
    echo "> make distclean"
    $DRY_RUN || make distclean
fi
run ./configure "${CONF_OPTS[@]}"
run make -j$NJOBS

if ! $DRY_RUN ; then
    for app in gevgen gmkspl ; do
        [ -x "$GENIE_DIR/bin/$app" ] || die "$app was not built. Check the compilation output"
    done
fi

# Save the environment with the build, so that the jobs use the same one
echo "> writing $GENIE_DIR/setup_LYON.sh"
if ! $DRY_RUN ; then
    sed "s|^export PYTHIA8PATH=.*|export PYTHIA8PATH=$PYTHIA8_DIR|" $ENV_SCRIPT > $GENIE_DIR/setup_LYON.sh
    if [ "$FSI" = "G4" ] && [ -n "$GEANT4_SETUP" ] ; then
        echo "source $GEANT4_SETUP" >> $GENIE_DIR/setup_LYON.sh
    fi
fi

echo ""
echo "GENIE compiled in $GENIE_DIR"
echo "Use it with: --grid-system LYON --genie-topdir $GENIE_DIR --tune $TUNE"
