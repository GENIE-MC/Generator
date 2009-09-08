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
#

script_location="http://hepunx.rl.ac.uk/~candreop/generators/GENIE/devel/src/scripts/build/ext/"

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

echo topdir $topdir
if [ ! -d ${topdir} ] ; then
  mkdir $topdir
fi
cd $topdir

for subdir in pythia6 lhapdf gsl root log4cpp; do
  if [ ! -d ${subdir} ] ; then
    mkdir ${subdir}
  fi
done

#
# wget or curl for retreiving remote files?
# (OS X doesn't generally have wget on it, so fall back on curl in that case)
# 
whichfetchit=`which wget | grep -v "no wget in"`
if [ ! -z "${whichfetchit}" ] ; then
  echo use \"wget\" for fetching files
  fetchit='wget '
else
  whichfetchit=`which curl | grep -v "no curl in"`
  if [ ! -z "${whichfetchit}" ] ; then
    # -f = fail without creating dummy, -O output local named like remoteza
    echo use \"curl -f -O\" for fetching files
    fetchit='curl -f -O '
  else
    echo "Neither wget nor curl available -- can't download files"
    exit 1
  fi
fi

#
# build pythia6
#
cd $topdir/pythia6
$fetchit $script_location/build_pythia6.sh
source build_pythia6.sh $pythia6_version

#
# build lhapdf
#
cd $topdir/lhapdf
$fetchit $script_location/build_lhapdf.sh
source build_lhapdf.sh $lhapdf_version

#
# build GSL
#
cd $topdir/gsl
$fetchit $script_location/build_gsl.sh
source build_gsl.sh $gsl_version

#
# build ROOT
#
cd $topdir/root
$fetchit $script_location/build_root.sh
source build_root.sh $root_version --arch=$arch

#
# build log4cpp
#
cd $topdir/log4cpp
$fetchit $script_location/build_log4cpp.sh
source build_log4cpp.sh $log4cpp_version

