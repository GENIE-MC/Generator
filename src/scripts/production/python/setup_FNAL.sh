#!/bin/bash

source /cvmfs/fermilab.opensciencegrid.org/products/genie/bootstrap_genie_ups.sh 
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups 
setup ifdhc v2_6_6
export IFDH_CP_MAXRETRIES=0

setup root v6_12_06a -q e17:debug 
setup lhapdf v5_9_1k -q e17:debug 
setup log4cpp v1_1_3a -q e17:debug 
setup pdfsets v5_9_1b 
setup gdb v8_1 
