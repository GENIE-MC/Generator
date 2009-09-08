#!/usr/bin/env bash

#
# A script to build the external GENIE dependencies
#
# Usage: ./build_ext.sh 
#              [--log4cpp-version=...] 
#              [--gsl-version=...] 
#              [--pythia6-version=...] 
#              [--lhapdf-version=...] 
#              [--root-version=...] 
#              [--topdir=...] 
#              [--arch=linux|linuxx8664gcc|macosx|...]
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

log4cpp_version=1.0
gsl_version=1.8
pythia6_version=6.4.12
lhapdf_version=5.7.0
root_version=5.24.00
arch=linuxx8664gcc
topdir=`pwd`

while [ $# -gt 0 ] ; do
  case $1 in
  *log4cpp*)
     log4cpp_version=$1 ;;
  *gsl*)
     gsl_version=$1 ;;
  *pythia6*)
     pythia6_version=$1 ;;
  *lhapdf*)
     lhapdf_version=$1 ;;
  *root*)
     root_version=$1 ;;
  *arch*)
     arch=$1 ;;
  *topdir*)
     topdir=$1 ;;
  esac
  shift
done

geniewww="http://hepunx.rl.ac.uk/~candreop/generators/GENIE/etc/installation"

echo topdir $topdir

if [ ! -d ${topdir} ] ; then
  mkdir $topdir
fi
cd $topdir

for subdir in pythia6 lhapdf gsl root log4cpp stage; do
  if [ ! -d ${subdir} ] ; then
    mkdir ${subdir}
  fi
done

cd $topdir/pythia6
curl -O $geniewww/build_pythia6.sh
source build_pythia6.sh $pythia6_version

cd $topdir/lhapdf
curl -O $geniewww/build_lhapdf.sh
source build_lhapdf.sh $lhapdf_version

cd $topdir/gsl
curl -O $geniewww/build_gsl.sh
source build_gsl.sh $gsl_version

cd $topdir/root
curl -O $geniewww/build_root.sh
source build_root.sh $root_version --arch=$arch

cd $topdir/log4cpp
curl -O $geniewww/build_log4cpp.sh
source build_log4cpp.sh $log4cpp_version

