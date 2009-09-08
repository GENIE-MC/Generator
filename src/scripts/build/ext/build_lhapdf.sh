#!/usr/bin/env bash

#
# A script to build LHAPDF (adapted from R.Hatcher's build_pythia6.sh)
#
# Usage:  ./build_lhapdf.sh [version] [--refetch]
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

version=5.7.0
doclean=0
refetch=0
while [ $# -gt 0 ] ; do
  case $1 in
  *clean*)
     doclean=1 ;;
  *refetch*)
     refetch=1 ;;
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
for subdir in stage download ; do
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

echo "$fetchit http://www.hepforge.org/archive/lhapdf/lhapdf-${major}.${minor}.${revis}.tar.gz"
$fetchit http://www.hepforge.org/archive/lhapdf/lhapdf-${major}.${minor}.${revis}.tar.gz

tar xzvf lhapdf-${major}.${minor}.${revis}.tar.gz
mv lhapdf-${major}.${minor}.${revis} ${topdir}/src
cd ${topdir}/src
./configure --prefix=${topdir}/stage/
make
make install
