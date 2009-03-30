#!/usr/bin/perl

#----------------------------------------------------------------------
# Submit jobs to generate sample sfor the NuINT09 generator comparisons
#
# Syntax:
#   perl submit-nuint09_jobs.pl <options>
#
# Options:
#  --run           : Comma separated list of run numbers
#  --version       : GENIE version number
# [--production]   :
# [--cycle]        :
# [--use-valgrind] :
#
# Examples:
#  perl submit-nuint09_jobs.pl --production nuint09 --cycle 01 --version v2.5.1 --run 1001
#  perl submit-nuint09_jobs.pl --production nuint09 --cycle 01 --version v2.5.1 --run all
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
# 1001   | 500k | numu    + C12    | 0.5     | COH-CC
# 1002   | 500k | numu    + C12    | 1.0     | COH-CC
# 1003   | 500k | numu    + C12    | 1.5     | COH-CC
# 1011   | 500k | numu    + C12    | 0.5     | COH-NC
# 1012   | 500k | numu    + C12    | 1.0     | COH-NC
# 1013   | 500k | numu    + C12    | 1.5     | COH-NC
# 2001   | 500k | numu    + C12    | 0.5     | QEL-CC
# 2001   | 500k | numu    + C12    | 1.0     | QEL-CC
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

$use_valgrind = 0                            unless defined $use_valgrind;
$production   = "nuint09\_$genie_version"    unless defined $production;
$cycle        = "01"                         unless defined $cycle;

$queue          = "prod";
$time_limit     = "30:00:00";
$genie_inst_dir = "/opt/ppd/t2k/GENIE/";
$genie_setup    = "$genie_inst_dir/$genie_version-setup";
$jobs_dir       = "/opt/ppd/t2k/GENIE/scratch/vA-$production\_$cycle";
$xspl_file      = "/opt/ppd/t2k/GENIE/data/job_inputs/xspl/gxspl-t2k-$genie_version.xml";
$mcseed         = 210921029;

%nevents_hash = ( 
  '1001' =>  '500000',
  '1002' =>  '500000',
  '1003' =>  '100000',
  '1011' =>  '500000',
  '1012' =>  '500000',
  '1013' =>  '100000',
  '2001' =>  '500000',
  '2002' =>  '100000'
);

%nupdg_hash = ( 
  '1001' =>  '14',
  '1002' =>  '14',
  '1003' =>  '14',
  '1011' =>  '14',
  '1012' =>  '14',
  '1013' =>  '14',
  '2001' =>  '14',
  '2002' =>  '14'
);

%tgtpdg_hash = ( 
  '1001' =>  '1000060120',
  '1002' =>  '1000060120',
  '1003' =>  '1000060120',
  '1011' =>  '1000060120',
  '1012' =>  '1000060120',
  '1013' =>  '1000060120',
  '2001' =>  '1000060120',
  '2002' =>  '1000060120'
);

%energy_hash = ( 
  '1001' =>  '0.5',
  '1002' =>  '1.0',
  '1003' =>  '1.5',
  '1011' =>  '0.5',
  '1012' =>  '1.0',
  '1013' =>  '1.5',
  '2001' =>  '0.5',
  '2002' =>  '1.0'
);

%gevgl_hash = ( 
  '1001' =>  'COH-CC',
  '1002' =>  'COH-CC',
  '1003' =>  'COH-CC',
  '1011' =>  'COH-NC',
  '1012' =>  'COH-NC',
  '1013' =>  'COH-NC',
  '2001' =>  'QEL-CC',
  '2002' =>  'QEL-CC'
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

    $batch_script  = "$jobs_dir/nuint09job-$curr_runnu.pbs";
    $logfile_evgen = "$jobs_dir/nuint09job-$curr_runnu.evgen.log";
    $logfile_conv  = "$jobs_dir/nuint09job-$curr_runnu.conv.log";
    $logfile_pbse  = "$jobs_dir/nuint09job-$curr_runnu.pbs_e.log";
    $logfile_pbso  = "$jobs_dir/nuint09job-$curr_runnu.pbs_o.log";

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
