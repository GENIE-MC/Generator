#!/bin/bash
source ../genie_env.sh
./configure \
  --enable-gsl \
  --enable-rwght \
  --enable-nnbar-oscillation \
  --enable-nucleon-decay \
  --enable-atmo \
  --disable-lhapdf5 \
  --enable-lhapdf6 \
  --with-optimiz-level=O2 \
  --with-pythia6-lib=${PYTHIA6} \
  --with-log4cpp-lib=${LOG4CPP_LIB} \
  --with-log4cpp-inc=${LOG4CPP_INC} \
  --compiler=gcc
