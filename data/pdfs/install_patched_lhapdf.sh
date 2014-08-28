#!/bin/sh

# Calculate the LHAPDF directory to store PDFs
LHAPDFBASE=""
if [ -n $LHAPATH ]; then
  echo "LHAPATH is $LHAPATH"
  echo " Using this as the PDF install directory."
  LHAPDFBASE=$LHAPATH
else 
  echo "LHAPATH is not set, checking for LIB or INC directories..."
  if [ -n $LHAPDF_LIB ]; then
    echo "LHAPDF_LIB is $LHAPDF_LIB"
    LHAPDFBASE=`dirname $LHAPDF_LIB`
  fi
  if [ -n $LHAPDF_INC ]; then
    echo "LHAPDF_INC is $LHAPDF_INC"
    LHAPDFTEMP=`dirname $LHAPDF_INC`
    if [ $LHAPDFBASE != $LHAPDFTEMP ]; then
      echo "The LHAPDF installation is confusing! Cannot install the PDF."
      exit 1
    fi
  fi
fi

if [ -z $LHAPDFBASE ]; then
  echo "Couldn't figure out where to put the PDF. Exiting..."
  exit 1
else
  echo "Installing the PDF set to $LHAPDFBASE"
fi
cp -v -b GRV98lo_pdflib.LHgrid $LHAPDFBASE/GRV98lo.LHgrid
