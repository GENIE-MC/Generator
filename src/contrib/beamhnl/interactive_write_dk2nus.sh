#!/bin/bash

# \author  John Plows <komninos-john.plows \at physics.ox.ac.uk>
#          University of Oxford
#
# \cpright Copyright (c) 2003-2025, The GENIE Collaboration
#          For the full text of the license visit http://copyright.genie-mc.org

INDIR="./" # Make sure sample_dk2nu.root lives in here.

OUTDIR="./dk2nuFluxes/" # Make sure this exists before running.
TAG=""

export INDIR=${INDIR}
export OUTDIR=${OUTDIR}

playlist="minervame"
#"minervamebar"
#"minervame"
beamconf="me000z200i"
#"me000z-200i"
#"me000z200i"

export GDK2NUFLUXXML="${PWD}/GNuMIFlux-dummy.xml"
echo "GDK2NUFLUXXML=${GDK2NUFLUXXML}"

XROOTDPREF="root://fndca1.fnal.gov:1095///pnfs/fnal.gov/usr"
#"root://fndca1.fnal.gov:1094///pnfs/fnal.gov/usr"

LSTR="${XROOTDPREF}/minerva/persistent/DataPreservation/flux/g4numiv6_dk2nu_${playlist}_${beamconf}"

FNUM="0000"
#"REPLACE_INDEX"
#"0001"
#"0000"
SUFFIX="0006.root"

FILE0="${LSTR}_${FNUM}_${SUFFIX}"

INSTRING="${FILE0}"
#"${LSTR}_${FNUM}_${SUFFIX}"

#echo ${INSTRING}
#exit 1

macro="write_dk2nus_reduced.C+(false,false,0,\"$INDIR\",\"$INSTRING\",\"$OUTDIR\",\"$TAG\",\"$playlist\",\"$beamconf\")" # don't keep all this info!
# macro="write_dk2nus.C+(false,false,0,\"$INDIR\",\"$INSTRING\",\"$OUTDIR\",\"$TAG\",\"$playlist\",\"$beamconf\")" # read and write all branches from dk2nu

echo "Calling genie"

genie -b -q $DK2NU/scripts/load_dk2nu.C\(true\) ${macro}
