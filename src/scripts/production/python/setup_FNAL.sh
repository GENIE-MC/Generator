#!/bin/bash

source /cvmfs/fermilab.opensciencegrid.org/products/genie/bootstrap_genie_ups.sh 
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups 
setup ifdhc v2_6_6
export IFDH_CP_MAXRETRIES=0

setup pdfsets v5_9_1b 
setup gdb v8_1 
