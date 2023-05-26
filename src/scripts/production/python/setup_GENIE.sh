#!/bin/bash 

setup git v2_15_1 

GITHUB_LOCATION=$1
git clone $GITHUB_LOCATION Generator 

export GENIEBASE=$(pwd)
export GENIE=$GENIEBASE/Generator 
export GENIE_BIN=$GENIE/bin
export GENIE_LIB=$GENIE/lib
export PYTHIA6=$PYTHIA_FQ_DIR/lib 
export LHAPDF5_INC=$LHAPDF_INC 
export LHAPDF5_LIB=$LHAPDF_LIB 
export XSECSPLINEDIR=$GENIEBASE/data 
export GENIE_REWEIGHT=$GENIEBASE/Reweight 
export PATH=$GENIE/bin:$GENIE_REWEIGHT/bin:$PATH 
export LD_LIBRARY_PATH=$GENIE_LIB:$GENIE_REWEIGHT/lib:$LD_LIBRARY_PATH 
unset GENIEBASE 

cd $GENIE 

GENIE_VERSION=master
if [ ! -z "$2" ] ; then 
  GENIE_VERSION=$2
  git checkout "$GENIE_VERSION" 
fi
echo Requested Genie version $GENIE_VERSION

if [ ! -z "$3" ] ; then
  GENIE_CONFIG_DIR=$3
  echo Requested Genie config $GENIE_CONFIG_DIR 
  if [ -d "$GENIE/$GENIE_CONFIG_DIR" ]; then
    GALGCONF=$GENIE/$GENIE_CONFIG_DIR
    export GALGCONF
    echo Genie config directory $GALGCONF
  else 
    if [ -d "$GENIE_CONFIG_DIR" ]; then 
      GALGCONF=$GENIE_CONFIG_DIR
      export GALGCONF
      echo Genie config directory $GALGCONF
    fi
  fi
else
  unset GALGCONF
fi

./configure \
--enable-rwght \
--enable-flux-drivers \
--enable-geom-drivers \
--disable-doxygen \
--enable-test \
--enable-mueloss \
--enable-dylibversion \
--enable-t2k \
--enable-fnal \
--enable-atmo \
--enable-nucleon-decay \
--enable-rwght \
--disable-masterclass \
--disable-debug \
--with-optimiz-level=O2 \
--disable-debug \
--disable-profiler  \
--enable-validation-tools \
--disable-cernlib \
--enable-lhapdf5 \
--with-log4cpp-inc=${LOG4CPP_INC} \
--with-log4cpp-lib=${LOG4CPP_LIB} \
--with-libxml2-inc=${LIBXML2_INC} \
--with-libxml2-lib=${LIBXML2_LIB} \
--with-lhapdf5-lib=${LHAPDF_LIB} \
--with-lhapdf5-inc=${LHAPDF_INC} \
--with-pythia6-lib=${PYTHIA_LIB}
#--enable-incl \
#--with-incl-inc=${INCL_FQ_DIR}/include/inclxx \
#--with-incl-lib=${INCL_FQ_DIR}/lib  \
#--with-boost-inc=${BOOST_FQ_DIR}/include \
#--with-boost-lib=${BOOST_FQ_DIR}/lib

make -j4
