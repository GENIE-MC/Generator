#!/bin/bash

source ../genie_env.sh

./configure \
  --enable-gsl \
  --enable-geant4 \
  --enable-rwght \
  --with-lhapdf6-lib=${LHAPDF_LIB} \
  --with-lhapdf6-inc=${LHAPDF_INC} \
  --disable-lhapdf5 \
  --enable-lhapdf6 \
  --with-optimiz-level=O2 \
  --with-log4cpp-inc=${LOG4CPP_INC} \
  --with-log4cpp-lib=${LOG4CPP_LIB} \
  --with-libxml2-inc=${LIBXML2_INC} \
  --with-libxml2-lib=${LIBXML2_LIB} \
  --disable-pythia8 \
  --enable-pythia6 \
  --enable-hepmc3 \
  --enable-incl \
  --with-incl-inc=${INCLXX_FQ_DIR}/include/inclxx \
  --with-incl-lib=${INCLXX_FQ_DIR}/lib \
  --with-boost-inc=${BOOST_FQ_DIR}/include \
  --with-boost-lib=${BOOST_FQ_DIR}/lib
  --with-pythia6-lib=${PYTHIA_LIB}
