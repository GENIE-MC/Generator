#!/usr/bin/env bash

#
# A script to build LHAPDF (adapted from R.Hatcher's build_pythia6.sh)
#
# Usage:  ./build_lhapdf.sh [version] [--refetch]
#
# Costas Andreopoulos <c.andreopoulos \at cern.ch>
#

version=6.5.6
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
  oopt="--output-document"
else
  whichfetchit=`which curl | grep -v "no curl in"`
  if [ ! -z "${whichfetchit}" ] ; then
    oopt="--output"
    # -f = fail without creating dummy
    # --location to follow redirection
    # use explicit --output as new hepforge URL is funky
    ####, -O output local named like remoteza
    echo use \"curl -f --location ${oopt}\" for fetching files
    fetchit='curl -f --location '
  else
    echo "Neither wget nor curl available -- can't download files"
    exit 1
  fi
fi

### example links
# https://lhapdf.hepforge.org/downloads/?f=old/lhapdf-5.7.0.tar.gz
# https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.5.6.tar.gz

if [ ${major} -eq 5 ]; then
   basepath=lhapdf-${major}.${minor}.${revis}
   tarfile=${basepath}.tar.gz
   fetchpath="?f=old/${tarfile}"
else
   basepath=LHAPDF-${major}.${minor}.${revis}
   tarfile=${basepath}.tar.gz
   fetchpath="?f=${tarfile}"
fi

echo "$fetchit ${oopt} ${tarfile} https://lhapdf.hepforge.org/downloads/${fetchpath}"
$fetchit ${oopt} ${tarfile}  "https://lhapdf.hepforge.org/downloads/${fetchpath}"

tar xzvf ${tarfile}
mv ${basepath} ${topdir}/src
cd ${topdir}/src
./configure --prefix=${topdir}/stage/
make
make install
