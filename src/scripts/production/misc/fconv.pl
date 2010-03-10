#-------------------------------------------------------------------------------------------
# File converter
#
# Syntax:
#    % perl fconv.pl <options>
#
# Options:
#    --rmin    : Minimum run number
#    --rmax    : Maximum run number
#    --prefix  : File prefix
#    --suffix  : File suffix
#    --dirname : File location [default: ./]
#    --convcmd : Conversion command
#
# Example 1:
#   Input: JNUBEAM HBOOK files nu.sk_horn250ka.RUNNU.hbk (where RUNNU=1,...,100)
#   Located at: /opt/data/t2k_flux/10a/sk/
#   Task: Convert to ROOT format
#   Type:
#     % perl fconv.pl --rmin 1 --rmax 100 --prefix nu.sk_horn250ka. --suffix .hbk \
#                     --dirname /opt/data/t2k_flux/10a/sk/ --convcmd 'h2root'
#
# Example 2:
#   Input: GENIE GHEP event files gntp.ghep.RUNNU.root (where RUNNU=1000,...,1010) 
#   Located at: /opt/ghepfiles/
#   Task: Convert first 10k events of each file in the GENIE summary ntuple format (GST)
#   Type:
#     % perl fconv.pl --rmin 100 --rmax 1010 --prefix gntp. --suffix .ghep.root \
#                     --dirname /opt/ghepfiles/ --convcmd 'gntpc -f gst -n 10000 -i '
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#-------------------------------------------------------------------------------------------

#!/usr/bin/perl

$dirname="./";

$iarg=0;
foreach (@ARGV) {
   if($_ eq '--prefix')   { $file_prfx  = $ARGV[$iarg+1]; }
   if($_ eq '--suffix')   { $file_sufx  = $ARGV[$iarg+1]; }
   if($_ eq '--rmin')     { $run_min    = $ARGV[$iarg+1]; }
   if($_ eq '--rmax')     { $run_max    = $ARGV[$iarg+1]; }
   if($_ eq '--convcmd')  { $convcmd    = $ARGV[$iarg+1]; }
   if($_ eq '--dirname')  { $dirname    = $ARGV[$iarg+1]; }
   $iarg++;
}
die("** Aborting [You need to specify a minimum run number using the --rmin option]")
unless defined $run_min;
die("** Aborting [You need to specify a maximum run number using the --rmax option]")
unless defined $run_max;
die("** Aborting [You need to specify a file prefix using the --prefix option]")
unless defined $file_prfx;
die("** Aborting [You need to specify a file suffix using the --suffix option]")
unless defined $file_sufx;
die("** Aborting [You need to specify the conversion command using the --convcmd option]")
unless defined $convcmd;

$i=0;
for($i=$run_min; $i<=$run_max; $i++) {
  $filename = $dirname."/".$file_prfx.$i.$file_sufx;
  if(-e $filename) {
    $cmd = $convcmd." ".$filename;
    print "Running: $cmd \n";
    `$cmd`;
#   system("$cmd");
  }
}
