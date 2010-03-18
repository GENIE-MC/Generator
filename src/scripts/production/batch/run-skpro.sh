#
# Submit all jobs for a standard SuperK MC production.
#
# Syntax:
#   ./run-skpro.sh [version] [production] [cycle]
#
# Notes:
# - See production/rutherford_tier2/submit-evg_sk_fhst.pl for details about each actual MC job.
# - See misc/generate_sk_flux_histograms.C for the derivation of SK flux from the JNUBEAM flux ntuples.
# - Using standard SK statistics: 500 numu, 250 numubar, 500 nue, 250 nuesig files (2k events each, see submit-evg_sk_fhst.pl)
# - Remember to switch on the 'decay' flag for all charmed hadrons in UserPhysicsOptions.xml
#
# C.Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

#!/bin/bash

VRS=$1
PRO=$2
CYC=$3

perl bulk-submit.pl --cmd "perl submit-evg_sk_fhst.pl --neutrino numu    --version $VRS --production $PRO --cycle $CYC" --arg "--run" --min 1 --max 500
perl bulk-submit.pl --cmd "perl submit-evg_sk_fhst.pl --neutrino numubar --version $VRS --production $PRO --cycle $CYC" --arg "--run" --min 1 --max 250
perl bulk-submit.pl --cmd "perl submit-evg_sk_fhst.pl --neutrino nue     --version $VRS --production $PRO --cycle $CYC" --arg "--run" --min 1 --max 500
perl bulk-submit.pl --cmd "perl submit-evg_sk_fhst.pl --neutrino nuesig  --version $VRS --production $PRO --cycle $CYC" --arg "--run" --min 1 --max 250
