#-------------------------------------------------------------------------------------------
# A script to bulk submit other scripts...
#
# Syntax:
#    % perl bulk-submit.pl <options>
#
# Options:
#    --cmd   : Actual script to execute *including* all fixed arguments
#    --arg   : "Running" argument
#    --min   : "Running" argument : minimum value
#    --max   : "Running" argument : maximum value
#
# Example 1:
#   To submit 100 nd280 MC jobs (processing flux files 100-200) using GENIE v2.6.0, type:
#   % perl bulk-submit.pl \
#        --cmd 'perl submit-evg_nd280.pl --version v2.6.0' \
#        --arg '--flux-run' \
#        --min 100 \
#        --max 200
#   See submit-evg_nd280.pl for more nd280 MC job options.
#
# Example 2:
#   To submit 100 standard `\nu_e signal' SuperK MC jobs using GENIE v2.6.0, type:
#   % perl bulk-submit.pl \
#        --cmd 'perl submit-evg_sk_fhst.pl --neutrino nuesig --version v2.6.0' \
#        --arg '--run' \
#        --min 1 \
#        --max 100
#   See submit-evg_sk_fhst.pl for more SuperK MC job options.
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
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
