#-------------------------------------------------------------------------------------------
# A script to bulk submit other scripts...
#
# Syntax:
#    % perl submit_bulk.pl <options>
#
# Options:
#    --cmd : Actual script to execute *including* all fixed arguments
#    --arg : "Running" argument
#    --min : "Running" argument : minimum value
#    --max : "Running" argument : maximum value
#
# Example 1:
#
#   To submit 100 T2K/ND280 MC jobs (processing flux files 100-200) using GENIE v3.06.0, 
#   type:
#
#   % perl submit_bulk.pl \
#        --cmd 'perl submit_mc_jobs_prod_t2k_nd280.pl --gen-version v3.06.00' \
#        --arg '--flux-run' \
#        --min 100 \
#        --max 200
#   See submit_mc_jobs_prod_t2k_nd280.pl for more ND280 MC job options.
#
# Example 2:
#
#   To submit 100 standard `\nu_e signal' T2K/SuperK MC jobs using GENIE v3.06.00, 
#   type:
#   % perl bulk-submit.pl \
#        --cmd 'perl submit_mc_jobs_prod_t2k_superk_fhst.pl --neutrino nuesig --gen-version v3.06.00' \
#        --arg '--run' \
#        --min 1 \
#        --max 100
#   See submit_mc_jobs_prod_t2k_superk_fhst.pl for more SuperK MC job options.
#
# Author:
#   Costas Andreopoulos <c.andreopoulos \st cern.ch>
#   University of Liverpool
#
# Copyright:
#   Copyright (c) 2003-2025, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#-------------------------------------------------------------------------------------------

#!/usr/bin/perl

$syntax = "%perl bulk-submit.pl --cmd <command> --arg <variable arg> --min <min> --max <max>";

$iarg=0;
foreach (@ARGV) {
   if($_ eq '--cmd')   { $cmd = $ARGV[$iarg+1]; }
   if($_ eq '--arg')   { $arg = $ARGV[$iarg+1]; }
   if($_ eq '--min')   { $min = $ARGV[$iarg+1]; }
   if($_ eq '--max')   { $max = $ARGV[$iarg+1]; }
   $iarg++;
}
die("** Syntax: $syntax \n") unless defined $cmd;
die("** Syntax: $syntax \n") unless defined $arg;
die("** Syntax: $syntax \n") unless defined $min;
die("** Syntax: $syntax \n") unless defined $max;

$i=0;
for($i=$min; $i<=$max; $i++) {
  $full_cmd = "$cmd $arg $i";
  print "Running: $full_cmd \n";
  `$full_cmd`;
}
