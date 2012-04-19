#-------------------------------------------------------------------------------------------
# Create XML file list expected by GENIE validation programs
#
# Syntax:
#    % perl mkgsimfile.pl <options>
#
# Options:
#    --dir    : file math
#    --model  : physics model tag
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#-------------------------------------------------------------------------------------------

#!/usr/bin/perl

$iarg=0;
foreach (@ARGV) {
   if($_ eq '--dir')     { $directory   = $ARGV[$iarg+1]; }
   if($_ eq '--model')   { $model_name  = $ARGV[$iarg+1]; }
   $iarg++;
}

die("** Aborting [Can not read the directory (if any) specified via the --dir option]")
unless -d $directory;
die("** Aborting [You need to specify a data file `tag' using the --model option]")
unless defined $model_name;

@files = <$directory/*.root>;

print "\n\n";
print "\t <model name=\"$model_name\"> \n";

foreach $file (@files) {
  if($file=~/.ghep.root$/) {
     print "\t\t <evt_file format=\"ghep\"> $file </evt_file> \n";
  }
  elsif($file=~/.gst.root$/) {
     print "\t\t <evt_file format=\"gst\"> $file </evt_file> \n";
  }
  else{
     print "\t\t <xsec_file> $file </xsec_file> \n";
  }
}

print "\t </model> \n\n";

