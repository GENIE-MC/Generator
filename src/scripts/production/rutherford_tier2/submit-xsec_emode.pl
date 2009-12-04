#!/usr/bin/perl

#----------------------------------------------------------------------
# Submit jobs for calculating GENIE eA cross section splines to be used 
# with GENIE's validation programs comparing GENIE against electron 
# scattering data.
#
# For use at the RAL/PPD Tier2 PBS batch farm.
#
# Syntax:
#   shell% perl submit-xsec_emode.pl <options>
#
# Options:
#    --version       : GENIE version number
#   [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]   : default: exspl_<version>
#   [--cycle]        : default: 01
#   [--use-valgrind] : default: off
#   [--queue]        : default: prod
#   [--softw-topdir] : default: /opt/ppd/t2k/GENIE
#
# Notes:
#   * Use GENIE gspladd utility to merge the job outputs
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#----------------------------------------------------------------------

use File::Path;

# inputs
#  
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir  = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$use_valgrind   = 0                       unless defined $use_valgrind;
$arch           = "SL5_64bit"             unless defined $arch;
$production     = "exspl\_$genie_version" unless defined $production;
$cycle          = "01"                    unless defined $cycle;
$queue          = "prod"                  unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/GENIE"    unless defined $softw_topdir;
$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$jobs_dir       = "$softw_topdir/scratch/$production\_$cycle/";

$nkots     = 200;
$emax      =  35;
$probes    = "11";
%targets = (
	'H1'    =>  '1000010020',
	'H2'    =>  '1000010020',
	'He4'   =>  '1000020040',
	'C12'   =>  '1000060120',
	'N14'   =>  '1000070140', 
	'O16'   =>  '1000080160', 
	'Ne20'  =>  '1000100200', 
	'Fe56'  =>  '1000260560',
	'Kr83'  =>  '1000360830',
	'Xe131' =>  '1000541310'  );

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

#
# loop over nuclear targets & submit jobs
#
while( my ($tgt_name, $tgt_code) = each %targets ) {

	$BATCH_SCRIPT = "$jobs_dir/job_$tgt_name.pbs";
	open(PBS, ">$BATCH_SCRIPT") 
	or die("Can not create the PBS batch script");

	$cmd = "gmkspl -p $probes -t $tgt_code -n $nkots -e $emax -o gxspl_emode_$tgt_name.xml &> job_$tgt_name.log";

	print PBS "#!/bin/bash \n";
        print PBS "#PBS -N $tgt_name \n";
        print PBS "#PBS -o $jobs_dir/job-$tgt_name.pbs_o \n";
        print PBS "#PBS -e $jobs_dir/job-$tgt_name.pbs_e \n";
	print PBS "source $genie_setup \n";
	print PBS "cd $jobs_dir \n";
	print PBS "export GEVGL=EM \n";
	print PBS "$cmd \n";

	print "EXEC: $cmd \n";

	# submit job
	`qsub -q $queue $BATCH_SCRIPT`;
}

