#!/usr/bin/perl

#----------------------------------------------------------------------
# Submit jobs to generate samples for the Wroclaw '09 Neutrino Winter 
# School generator comparisons.
# Script prepared for the RAL/PPD Tier2 batch farm.
#
# Syntax:
#   perl submit-wro09_winter_school_jobs.pl <options>
#
# Options:
#  --run           : Comma separated list of run numbers
#  --version       : GENIE version number
# [--production]   : production name
# [--cycle]        : cycle name
# [--use-valgrind] :
#
# Examples:
#  perl submit-wro09_winter_school_jobs.pl --production wroclaw09 --cycle 01 --version v2.5.1 --run 1001
#  perl submit-wro09_winter_school_jobs.pl --production wroclaw09 --cycle 01 --version v2.5.1 --run all
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#----------------------------------------------------------------------
#
# SAMPLES:
#......................................................................
#  run   | nev  |  init state      | energy  | processes
#  nu.   |      |                  | (GeV)   | enabled
#......................................................................
# 1000   | 500k | numu    + O16    | 0.6     | all
# 1001   | 500k | numu    + O16    | 1.0     | all
# 1002   | 500k | numu    + O16    | 5.0     | all
# 1100   | 500k | numu    + Ar40   | 0.6     | all
# 1101   | 500k | numu    + Ar40   | 1.0     | all
# 1102   | 500k | numu    + Ar40   | 5.0     | all
# 1200   | 500k | numu    + Fe56   | 0.6     | all
# 1201   | 500k | numu    + Fe56   | 1.0     | all
# 1202   | 500k | numu    + Fe56   | 5.0     | all
# 2000   | 500k | numubar + O16    | 0.6     | all
# 2001   | 500k | numubar + O16    | 1.0     | all
# 2002   | 500k | numubar + O16    | 5.0     | all
# 2100   | 500k | numubar + Ar40   | 0.6     | all
# 2101   | 500k | numubar + Ar40   | 1.0     | all
# 2102   | 500k | numubar + Ar40   | 5.0     | all
# 2200   | 500k | numubar + Fe56   | 0.6     | all
# 2201   | 500k | numubar + Fe56   | 1.0     | all
# 2202   | 500k | numubar + Fe56   | 5.0     | all
# 
#......................................................................
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--run')           { $runnu         = $ARGV[$iarg+1]; }
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined benchmark runs #. Use the --run option]")
unless defined $runnu;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$use_valgrind   = 0                         unless defined $use_valgrind;
$production     = "wroclaw\_$genie_version" unless defined $production;
$cycle          = "01"                      unless defined $cycle;

$queue          = "prod";
$time_limit     = "30:00:00";
$topdir         = "/opt/ppd/t2k";
$genie_inst_dir = "$topdir/GENIE/";
$genie_setup    = "$genie_inst_dir/$genie_version-setup";
$jobs_dir       = "$topdir/GENIE/scratch/$production\_$cycle";
$xspl_file      = "$topdir/GENIE/data/job_inputs/xspl/gxspl-t2k-$genie_version.xml";
$mcseed         = 210921029;

%nevents_hash = ( 
  '1000' =>  '500000',
  '1001' =>  '500000',
  '1002' =>  '100000',
  '1100' =>  '500000',
  '1101' =>  '500000',
  '1102' =>  '100000',
  '1200' =>  '500000',
  '1201' =>  '500000',
  '1202' =>  '100000',
  '2000' =>  '500000',
  '2001' =>  '500000',
  '2002' =>  '100000',
  '2100' =>  '500000',
  '2101' =>  '500000',
  '2102' =>  '100000',
  '2200' =>  '500000',
  '2201' =>  '500000',
  '2202' =>  '100000'
);

%nupdg_hash = ( 
  '1000' =>  '14',
  '1001' =>  '14',
  '1002' =>  '14',
  '1100' =>  '14',
  '1101' =>  '14',
  '1102' =>  '14',
  '1200' =>  '14',
  '1201' =>  '14',
  '1202' =>  '14',
  '2000' =>  '-14',
  '2001' =>  '-14',
  '2002' =>  '-14',
  '2100' =>  '-14',
  '2101' =>  '-14',
  '2102' =>  '-14',
  '2200' =>  '-14',
  '2201' =>  '-14',
  '2202' =>  '-14'
);

%tgtpdg_hash = ( 
  '1000' =>  '1000080160',
  '1001' =>  '1000080160',
  '1002' =>  '1000080160',
  '1100' =>  '1000180400',
  '1101' =>  '1000180400',
  '1102' =>  '1000180400',
  '1200' =>  '1000260560',
  '1201' =>  '1000260560',
  '1202' =>  '1000260560',
  '2000' =>  '1000080160',
  '2001' =>  '1000080160',
  '2002' =>  '1000080160',
  '2100' =>  '1000180400',
  '2101' =>  '1000180400',
  '2102' =>  '1000180400',
  '2200' =>  '1000260560',
  '2201' =>  '1000260560',
  '2202' =>  '1000260560'
);

%energy_hash = ( 
  '1000' =>  '0.6',
  '1001' =>  '1.0',
  '1002' =>  '5.0',
  '1100' =>  '0.6',
  '1101' =>  '1.0',
  '1102' =>  '5.0',
  '1200' =>  '0.6',
  '1201' =>  '1.0',
  '1202' =>  '5.0',
  '2000' =>  '0.6',
  '2001' =>  '1.0',
  '2002' =>  '5.0',
  '2100' =>  '0.6',
  '2101' =>  '1.0',
  '2102' =>  '5.0',
  '2200' =>  '0.6',
  '2201' =>  '1.0',
  '2202' =>  '5.0'
);

%gevgl_hash = ( 
  '1000' =>  'Default',
  '1001' =>  'Default',
  '1002' =>  'Default',
  '1100' =>  'Default',
  '1101' =>  'Default',
  '1102' =>  'Default',
  '1200' =>  'Default',
  '1201' =>  'Default',
  '1202' =>  'Default',
  '2000' =>  'Default',
  '2001' =>  'Default',
  '2002' =>  'Default',
  '2100' =>  'Default',
  '2101' =>  'Default',
  '2102' =>  'Default',
  '2200' =>  'Default',
  '2201' =>  'Default',
  '2202' =>  'Default'
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

print "Input runs: $runnu \n";

for my $curr_runnu (keys %gevgl_hash)  {
  print "Checking benchmark run: ...... $curr_runnu \n";

  if($runnu=~m/$curr_runnu/ || $runnu eq "all") {
    print "** matched -> submitting job \n";

    #
    # get runnu-dependent info
    #
    $nev   = $nevents_hash {$curr_runnu};
    $nu    = $nupdg_hash   {$curr_runnu};
    $tgt   = $tgtpdg_hash  {$curr_runnu};
    $en    = $energy_hash  {$curr_runnu};
    $gevgl = $gevgl_hash   {$curr_runnu};

    $batch_script  = "$jobs_dir/wro09job-$curr_runnu.pbs";
    $logfile_evgen = "$jobs_dir/wro09job-$curr_runnu.evgen.log";
    $logfile_conv  = "$jobs_dir/wro09job-$curr_runnu.conv.log";
    $logfile_pbse  = "$jobs_dir/wro09job-$curr_runnu.pbs_e.log";
    $logfile_pbso  = "$jobs_dir/wro09job-$curr_runnu.pbs_o.log";

    $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
    $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
    $evgen_cmd     = "gevgen -n $nev -s -e $en -p $nu -t $tgt -r $curr_runnu | grep_pipe &> $logfile_evgen";
    $conv_cmd      = "gntpc -f gst -i gntp.$curr_runnu.ghep.root | grep -B 100 -A 30 -i \"warn\\|error\\|fatal\" &> $logfile_conv";

    # create the PBS script
    #
    open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
    print PBS "#!/bin/bash \n";
    print PBS "#PBS -l cput=$time_limit \n";
    print PBS "#PBS -o $logfile_pbso \n";
    print PBS "#PBS -e $logfile_pbse \n";
    print PBS "source $genie_setup \n"; 
    print PBS "cd $jobs_dir \n";
    print PBS "export GSPLOAD=$xspl_file \n";
    print PBS "export GEVGL=$gevgl \n";
    print PBS "export GSEED=$mcseed  \n";
    print PBS "$evgen_cmd \n";
    print PBS "$conv_cmd \n";

    print "EXEC: $evgen_cmd \n";

    # submit job
    #
    `qsub -q $queue $batch_script`;
  }
}
