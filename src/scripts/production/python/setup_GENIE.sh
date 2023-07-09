#!/bin/bash 

setup git v2_15_1 

GITHUB_LOCATION=$1
GENIE_VERSION=master
if [ ! -z "$2" ] ; then 
  GENIE_VERSION=$2
  git checkout "$GENIE_VERSION" 
fi
echo Requested Genie version $GENIE_VERSION
CONFIGURE_INCL=$3
CONFIGURE_G4=$4
GENIE_CONFIG_DIR=$5

git clone $GITHUB_LOCATION Generator 

if [ "$CONFIGURE_INCL" = "true" ] ; then
    setup gcc v7_3_0
    setup cmake v3_13_1
    setup boost v1_66_0a -q e17:prof

    export BASE=$(pwd)
    export INCL_TOP=${BASE}/inclxx
    export INCL_PLATFORM=Linux64bit+2.6-2.12
    export INCL_VERSION=inclxx-v5.2.9.5-6aca7e6
    export INCL_DIR=${INCL_TOP}/${INCL_VERSION}
    export INCL_SRC_DIR=${INCL_DIR}/source/${INCL_VERSION}  # needed to find de-excitation data files
    export INCL_FQ_DIR=${INCL_DIR}/${INCL_PLATFORM}

    mkdir ${INCL_TOP}; mkdir -p ${INCL_DIR} ; mkdir ${INCL_DIR}/source/ ;
    cp /genie/app/rhatcher/genie_inclxx/${INCL_VERSION}.tar.gz ${INCL_DIR}/source ## THIS WONT WORK
    cd ${INCL_DIR}; cd source; tar xvzf ${INCL_VERSION}.tar.gz
    #CVMFS permissions
    find . -type f -exec chmod +r {} \;
    find . -type d -exec chmod +rx {} \;
    cd ${INCL_DIR} ; mkdir build-temp-${INCL_PLATFORM}; cd build-temp-${INCL_PLATFORM}
    cmake -DUSE_FPIC=ON -DCMAKE_INSTALL_PREFIX:PATH=${INCL_DIR}/${INCL_PLATFORM} ${INCL_SRC_DIR}
    make -j4; make install
    cd ${BASE}
elif [ "$CONFIGURE_G4" = "true" ] ; then
    source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
    setup geant4 v4_10_3_p03e -q debug:e17
fi

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

if [ ! -z "$5" ] ; then
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

if [ "$CONFIGURE_INCL" = "true" ] ; then
    echo Configuring INCL... 
    ./configure \
	--enable-gsl \
	--with-optimiz-level=O2 \
	--with-log4cpp-inc=${LOG4CPP_INC} \
	--with-log4cpp-lib=${LOG4CPP_LIB} \
	--with-libxml2-inc=${LIBXML2_INC} \
	--with-libxml2-lib=${LIBXML2_LIB} \
	--with-lhapdf5-lib=${LHAPDF_LIB} \
	--with-lhapdf5-inc=${LHAPDF_INC} \
	--with-pythia6-lib=${PYTHIA_LIB} \
	--enable-incl \
	--with-incl-inc=${INCL_FQ_DIR}/include/inclxx \
	--with-incl-lib=${INCL_FQ_DIR}/lib  \
	--with-boost-inc=${BOOST_FQ_DIR}/include \
	--with-boost-lib=${BOOST_FQ_DIR}/lib
elif [ "$CONFIGURE_G4" = "true" ] ; then
    echo Configuring G4...
    ./configure \
	--enable-gsl \
	--with-optimiz-level=O2 \
	--with-log4cpp-inc=${LOG4CPP_INC} \
	--with-log4cpp-lib=${LOG4CPP_LIB} \
	--with-libxml2-inc=${LIBXML2_INC} \
	--with-libxml2-lib=${LIBXML2_LIB} \
	--with-lhapdf5-lib=${LHAPDF_LIB} \
	--with-lhapdf5-inc=${LHAPDF_INC} \
	--with-pythia6-lib=${PYTHIA_LIB} \
	--enable-geant4
else 
    ./configure \
	--enable-gsl \
	--enable-rwght \
	--with-optimiz-level=O2 \
	--with-log4cpp-inc=${LOG4CPP_INC} \
	--with-log4cpp-lib=${LOG4CPP_LIB} \
	--with-libxml2-inc=${LIBXML2_INC} \
	--with-libxml2-lib=${LIBXML2_LIB} \
	--with-lhapdf5-lib=${LHAPDF_LIB} \
	--with-lhapdf5-inc=${LHAPDF_INC} \
	--with-pythia6-lib=${PYTHIA_LIB}
fi

#make -j4
