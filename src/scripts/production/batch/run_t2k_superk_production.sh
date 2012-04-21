#
# Submit all jobs for a standard SuperK MC production.
#
# Syntax:
#   ./run_t2k_superk_production.sh [version] [production] [cycle]
#
# Example:
#   ./run_t2k_superk_production.sh v2.6.2 2011a 01
#
# Notes:
# - See production/rutherford_tier2/submit-evg_sk_fhst.pl for details about each actual MC job.
# - See misc/generate_sk_flux_histograms.C for the derivation of SK flux from the JNUBEAM flux ntuples.
# - Using standard SK statistics: 500 numu, 250 numubar, 500 nue, 250 nuebar, 250 nuesig files (2k events each, see submit_t2k_superk_fhst_mc_jobs.pl)
# - Remember to switch on the 'decay' flag for all charmed hadrons in UserPhysicsOptions.xml
#
# C.Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

#!/bin/bash

VRS=$1
PRO=$2
CYC=$3

perl bulk_submit.pl --cmd "perl submit_t2k_superk_fhst_mc_jobs.pl --neutrino numu    --version $VRS --production $PRO --cycle $CYC" --arg "--run" --min 1 --max 500
perl bulk_submit.pl --cmd "perl submit_t2k_superk_fhst_mc_jobs.pl --neutrino numubar --version $VRS --production $PRO --cycle $CYC" --arg "--run" --min 1 --max 250
perl bulk_submit.pl --cmd "perl submit_t2k_superk_fhst_mc_jobs.pl --neutrino nue     --version $VRS --production $PRO --cycle $CYC" --arg "--run" --min 1 --max 500
perl bulk_submit.pl --cmd "perl submit_t2k_superk_fhst_mc_jobs.pl --neutrino nuebar  --version $VRS --production $PRO --cycle $CYC" --arg "--run" --min 1 --max 250
perl bulk_submit.pl --cmd "perl submit_t2k_superk_fhst_mc_jobs.pl --neutrino nuesig  --version $VRS --production $PRO --cycle $CYC" --arg "--run" --min 1 --max 250

