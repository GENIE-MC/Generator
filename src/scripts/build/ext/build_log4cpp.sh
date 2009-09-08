#!/usr/bin/env bash

#
# A script to build log4cpp (adapted from R.Hatcher's build_pythia6.sh)
#
# Usage:  ./build_log4cpp.sh [version] [--refetch]
#
# Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

version=1.0
doclean=0
if [ $# > 0 ] ; then
  if [ -n "$1" ] ; then
    version=$1
  fi
  if [ "$2" == "clean" ] ; then
    doclean=1
  fi
fi
major=`awk -v str=$version 'BEGIN {split(str, tk, "."); print tk[1]}'`
minor=`awk -v str=$version 'BEGIN {split(str, tk, "."); print tk[2]}'`

topdir=`pwd`/v${major}_${minor}

echo version $version major $major minor $minor 
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

echo "$fetchit http://downloads.sourceforge.net/project/log4cpp/log4cpp-${major}.${minor}.x%20%28current%29/log4cpp-${major}.${minor}/log4cpp-${major}.${minor}.tar.gz"
$fetchit http://downloads.sourceforge.net/project/log4cpp/log4cpp-${major}.${minor}.x%20%28current%29/log4cpp-${major}.${minor}/log4cpp-${major}.${minor}.tar.gz

tar xzvf log4cpp-${major}.${minor}.tar.gz
mv log4cpp-${major}.${minor} ${topdir}/src

cd ${topdir}/src
./configure --prefix=${topdir}/stage/
make
make install


