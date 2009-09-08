#!/usr/bin/env bash

#
# A script to build ROOT (adapted from R.Hatcher's build_pythia6.sh)
#
# Usage:  ./build_root.sh [version] 
#                         [--arch=linux|linuxx8664gcc|macosx|...] 
#                         [--refetch]
#                         [--with-pythia6-version=...]
#                         [--with-gsl-version=...]
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

version=5.24.00
arch=linux
doclean=0
refetch=0
while [ $# -gt 0 ] ; do
  case $1 in
  *clean*)
     doclean=1 ;;
  *refetch*)
     refetch=1 ;;
  *arch*)
     refetch=$1 ;;
  *pythia6*)
     pythia6_version=$1 ;;
  *gsl*)
     gsl_version=$1 ;;
  *)
     version=$1 ;;
  esac
  shift
done

major=`awk -v str=$version 'BEGIN {split(str, tk, "."); print tk[1]}'`
minor=`awk -v str=$version 'BEGIN {split(str, tk, "."); print tk[2]}'`
revis=`awk -v str=$version 'BEGIN {split(str, tk, "."); print tk[3]}'`

topdir=`pwd`/v${major}_${minor}_${revis}

echo version $version major $major minor $minor revision $revis
echo topdir $topdir

if [ ! -d ${topdir} ] ; then
  mkdir $topdir
fi
cd $topdir
for subdir in include lib download ; do
  if [ ! -d ${subdir} ] ; then
    mkdir ${subdir}
  fi
done

cd ${topdir}/download

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

echo "$fetchit http://root.cern.ch/root/root_v${major}.${minor}.${revis}.source.tar.gz"
$fetchit http://root.cern.ch/root/root_v${major}.${minor}.${revis}.source.tar.gz

tar xzvf http://root.cern.ch/root/root_v${major}.${minor}.${revis}.source.tar.gz
mv root ${topdir}/src

cd ${topdir}/src

export ROOTSYS=`pwd`
export PATH=$ROOTSYS/bin:$PATH
export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH

./configure $arch --enable-pythia6 --enable-mathmore --with-gsl-incdir= --with-gsl-libdir= --with-pythia6-libdir=
make


