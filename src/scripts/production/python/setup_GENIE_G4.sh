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

setup geant4 v4_10_3_p03e -q debug:e17

cd $GENIE 

GENIE_VERSION=master
if [ ! -z "$2" ] ; then 
  GENIE_VERSION=$2
  git checkout "$GENIE_VERSION" 
fi
echo Requested Genie version $GENIE_VERSION

echo $3
if [ ! -z "$3" ] ; then
  GENIE_CONFIG_DIR=$3
  echo Requested Genie config $GENIE_CONFIG_DIR 
  if [ -d "$GENIE_CONFIG_DIR" ]; then 
      GALGCONF=$GENIE_CONFIG_DIR
      export GALGCONF
      echo Genie config directory $GALGCONF
  else
      if [ -d "$GENIE/$GENIE_CONFIG_DIR" ]; then
	  GALGCONF=$GENIE/$GENIE_CONFIG_DIR
	  export GALGCONF
	  echo Genie config directory $GALGCONF
      fi
  fi
else
  unset GALGCONF
fi

./configure \
--enable-gsl \
--disable-debug \
--with-optimiz-level=O2 \
--with-log4cpp-inc=${LOG4CPP_INC} \
--with-log4cpp-lib=${LOG4CPP_LIB} \
--with-libxml2-inc=${LIBXML2_INC} \
--with-libxml2-lib=${LIBXML2_LIB} \
--with-lhapdf5-lib=${LHAPDF_LIB} \
--with-lhapdf5-inc=${LHAPDF_INC} \
--with-pythia6-lib=${PYTHIA_LIB} \
--enable-geant4

make -j4
