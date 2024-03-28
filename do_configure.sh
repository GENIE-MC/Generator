#!/bin/bash
source ../genie_env.sh
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
