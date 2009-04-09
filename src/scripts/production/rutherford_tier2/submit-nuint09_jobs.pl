#!/usr/bin/perl

#-----------------------------------------------------------------------
# Submit jobs to generate requested samples for the NuINT09 `Confronting 
# theory, models & data' session organized by S.Dytman
#
# Syntax:
#   perl submit-nuint09_jobs.pl <options>
#
# Options:
#  --run           : Comma separated list of run numbers
#  --version       : GENIE version number
# [--nsubruns]     : number of subruns per run
# [--production]   :
# [--cycle]        :
# [--use-valgrind] :
#
# Examples:
#  perl submit-nuint09_jobs.pl --production nuint09 --cycle 01 --version v2.5.1 --run 1001 --nsubruns 5
#  perl submit-nuint09_jobs.pl --production nuint09 --cycle 01 --version v2.5.1 --run all  --nsubruns 5
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#----------------------------------------------------------------------
#
# SAMPLES:
#......................................................................
# run number      |  init state      | energy  | processes
#                 |                  | (GeV)   | enabled
# (xx=subrun_id)  |
#......................................................................
#
# 1001xx          | numu    + C12    | 0.5     | COH-CC
# 1002xx          | numu    + C12    | 1.0     | COH-CC
# 1003xx          | numu    + C12    | 1.5     | COH-CC
# 1011xx          | numu    + C12    | 0.5     | COH-NC
# 1012xx          | numu    + C12    | 1.0     | COH-NC
# 1013xx          | numu    + C12    | 1.5     | COH-NC
# 2001xx          | numu    + C12    | 0.5     | QEL-CC
# 2002xx          | numu    + C12    | 1.0     | QEL-CC
# 3002xx          | numu    + C12    | 1.0     | RES-CC
# 3003xx          | numu    + C12    | 1.5     | RES-CC 
# 9002xx          | numu    + C12    | 1.0     | CC (all)
# 9003xx          | numu    + C12    | 1.5     | CC (all)
#
# 1101xx          | numu    + O16    | 0.5     | COH-CC
# 1102xx          | numu    + O16    | 1.0     | COH-CC
# 1103xx          | numu    + O16    | 1.5     | COH-CC
# 1111xx          | numu    + O16    | 0.5     | COH-NC
# 1112xx          | numu    + O16    | 1.0     | COH-NC
# 1113xx          | numu    + O16    | 1.5     | COH-NC
# 2101xx          | numu    + O16    | 0.5     | QEL-CC
# 2102xx          | numu    + O16    | 1.0     | QEL-CC
# 3102xx          | numu    + O16    | 1.0     | RES-CC
# 3103xx          | numu    + O16    | 1.5     | RES-CC
# 9102xx          | numu    + O16    | 1.0     | CC (all)
# 9103xx          | numu    + O16    | 1.5     | CC (all)
#
# 1201xx          | numu    + Fe56   | 0.5     | COH-CC
# 1202xx          | numu    + Fe56   | 1.0     | COH-CC
# 1203xx          | numu    + Fe56   | 1.5     | COH-CC
# 1211xx          | numu    + Fe56   | 0.5     | COH-NC
# 1212xx          | numu    + Fe56   | 1.0     | COH-NC
# 1213xx          | numu    + Fe56   | 1.5     | COH-NC
# 2201xx          | numu    + Fe56   | 0.5     | QEL-CC
# 2202xx          | numu    + Fe56   | 1.0     | QEL-CC
# 3202xx          | numu    + Fe56   | 1.0     | RES-CC
# 3203xx          | numu    + Fe56   | 1.5     | RES-CC
# 9202xx          | numu    + Fe56   | 1.0     | CC (all)
# 9203xx          | numu    + Fe56   | 1.5     | CC (all)
#......................................................................
# run number key: IJKKxx
# I  : Enabled processes (1:COHCC, 2:QELCC, 3:RESCC, 9:ALLCC)
# J  : Initial state     (0:nu_mu+C12, 1:nu_mu+O16, 2:nu_mu+Fe56)
# KK : Energy            (01: 0.5 GeV, 02: 1.0 GeV, 03: 1.5 GeV)
# xx : Run id            (01-99, 100k events each)
#......................................................................
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--run')           { $runnu         = $ARGV[$iarg+1]; }
  if($_ eq '--nsubruns')      { $nsubruns      = $ARGV[$iarg+1]; }
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

$nsubruns       = 1                         unless defined $nsubruns;
$use_valgrind   = 0                         unless defined $use_valgrind;
$production     = "nuint09\_$genie_version" unless defined $production;
$cycle          = "01"                      unless defined $cycle;

$queue          = "prod";
$time_limit     = "60:00:00";
$topdir         = "/opt/ppd/t2k";
$genie_inst_dir = "$topdir/GENIE/";
$genie_setup    = "$genie_inst_dir/$genie_version-setup";
$jobs_dir       = "$topdir/GENIE/scratch/$production\_$cycle";
$xspl_file      = "$topdir/GENIE/data/job_inputs/xspl/gxspl-t2k-$genie_version.xml";
$mcseed         = 210921029;
$nev_per_subrun = 100000;

%nupdg_hash = ( 
  '1001' =>  '14',
  '1002' =>  '14',
  '1003' =>  '14',
  '1011' =>  '14',
  '1012' =>  '14',
  '1013' =>  '14',
  '2001' =>  '14',
  '2002' =>  '14',
  '3002' =>  '14',
  '3003' =>  '14',
  '9002' =>  '14',
  '9003' =>  '14',
  '1101' =>  '14',
  '1102' =>  '14',
  '1103' =>  '14',
  '1111' =>  '14',
  '1112' =>  '14',
  '1113' =>  '14',
  '2101' =>  '14',
  '2102' =>  '14',
  '3102' =>  '14',
  '3103' =>  '14',
  '9102' =>  '14',
  '9103' =>  '14',
  '1201' =>  '14',
  '1202' =>  '14',
  '1203' =>  '14',
  '1211' =>  '14',
  '1212' =>  '14',
  '1213' =>  '14',
  '2201' =>  '14',
  '2202' =>  '14',
  '3202' =>  '14',
  '3203' =>  '14',
  '9202' =>  '14',
  '9203' =>  '14'
);

%tgtpdg_hash = ( 
  '1001' =>  '1000060120',
  '1002' =>  '1000060120',
  '1003' =>  '1000060120',
  '1011' =>  '1000060120',
  '1012' =>  '1000060120',
  '1013' =>  '1000060120',
  '2001' =>  '1000060120',
  '2002' =>  '1000060120',
  '3002' =>  '1000060120',
  '3003' =>  '1000060120',
  '9002' =>  '1000060120',
  '9003' =>  '1000060120',
  '1101' =>  '1000080160',
  '1102' =>  '1000080160',
  '1103' =>  '1000080160',
  '1111' =>  '1000080160',
  '1112' =>  '1000080160',
  '1113' =>  '1000080160',
  '2101' =>  '1000080160',
  '2102' =>  '1000080160',
  '3102' =>  '1000080160',
  '3103' =>  '1000080160',
  '9102' =>  '1000080160',
  '9103' =>  '1000080160',
  '1201' =>  '1000260560',
  '1202' =>  '1000260560',
  '1203' =>  '1000260560',
  '1211' =>  '1000260560',
  '1212' =>  '1000260560',
  '1213' =>  '1000260560',
  '2201' =>  '1000260560',
  '2202' =>  '1000260560',
  '3202' =>  '1000260560',
  '3203' =>  '1000260560',
  '9202' =>  '1000260560',
  '9203' =>  '1000260560'
);

%energy_hash = ( 
  '1001' =>  '0.5',
  '1002' =>  '1.0',
  '1003' =>  '1.5',
  '1011' =>  '0.5',
  '1012' =>  '1.0',
  '1013' =>  '1.5',
  '2001' =>  '0.5',
  '2002' =>  '1.0',
  '3002' =>  '1.0',
  '3003' =>  '1.5',
  '9002' =>  '1.0',
  '9003' =>  '1.5',
  '1101' =>  '0.5',
  '1102' =>  '1.0',
  '1103' =>  '1.5',
  '1111' =>  '0.5',
  '1112' =>  '1.0',
  '1113' =>  '1.5',
  '2101' =>  '0.5',
  '2102' =>  '1.0',
  '3102' =>  '1.0',
  '3103' =>  '1.5',
  '9102' =>  '1.0',
  '9103' =>  '1.5',
  '1201' =>  '0.5',
  '1202' =>  '1.0',
  '1203' =>  '1.5',
  '1211' =>  '0.5',
  '1212' =>  '1.0',
  '1213' =>  '1.5',
  '2201' =>  '0.5',
  '2202' =>  '1.0',
  '3202' =>  '1.0',
  '3203' =>  '1.5',
  '9202' =>  '1.0',
  '9203' =>  '1.5'
);

%gevgl_hash = ( 
  '1001' =>  'COH-CC',
  '1002' =>  'COH-CC',
  '1003' =>  'COH-CC',
  '1011' =>  'COH-NC',
  '1012' =>  'COH-NC',
  '1013' =>  'COH-NC',
  '2001' =>  'QEL-CC',
  '2002' =>  'QEL-CC',
  '3002' =>  'RES-CC',
  '3003' =>  'RES-CC',
  '9002' =>  'CC',
  '9003' =>  'CC',
  '1101' =>  'COH-CC',
  '1102' =>  'COH-CC',
  '1103' =>  'COH-CC',
  '1111' =>  'COH-NC',
  '1112' =>  'COH-NC',
  '1113' =>  'COH-NC',
  '2101' =>  'QEL-CC',
  '2102' =>  'QEL-CC',
  '3102' =>  'RES-CC',
  '3103' =>  'RES-CC',
  '9102' =>  'CC',
  '9103' =>  'CC',
  '1201' =>  'COH-CC',
  '1202' =>  'COH-CC',
  '1203' =>  'COH-CC',
  '1211' =>  'COH-NC',
  '1212' =>  'COH-NC',
  '1213' =>  'COH-NC',
  '2201' =>  'QEL-CC',
  '2202' =>  'QEL-CC',
  '3202' =>  'RES-CC',
  '3203' =>  'RES-CC',
  '9202' =>  'CC',
  '9203' =>  'CC'
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

print "Input runs: $runnu \n";

for my $curr_runnu (keys %gevgl_hash)  {
  #print "Checking benchmark run: ...... $curr_runnu \n";

  if($runnu=~m/$curr_runnu/ || $runnu eq "all") {
    print "** submitting run: $curr_runnu \n";

    #
    # get runnu-dependent info
    #
    $nu    = $nupdg_hash   {$curr_runnu};
    $tgt   = $tgtpdg_hash  {$curr_runnu};
    $en    = $energy_hash  {$curr_runnu};
    $gevgl = $gevgl_hash   {$curr_runnu};

    #
    # submit subruns
    #
    for($isubrun = 0; $isubrun < $nsubruns; $isubrun++) {

       $curr_subrunnu = 100 * $curr_runnu + $isubrun;

       $batch_script  = "$jobs_dir/nuint09job-$curr_subrunnu.pbs";
       $logfile_evgen = "$jobs_dir/nuint09job-$curr_subrunnu.evgen.log";
       $logfile_conv  = "$jobs_dir/nuint09job-$curr_subrunnu.conv.log";
       $logfile_pbse  = "$jobs_dir/nuint09job-$curr_subrunnu.pbs_e.log";
       $logfile_pbso  = "$jobs_dir/nuint09job-$curr_subrunnu.pbs_o.log";

       $curr_seed     = $mcseed + $isubrun;
       $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $evgen_cmd     = "gevgen -n $nev_per_subrun -s -e $en -p $nu -t $tgt -r $curr_subrunnu | grep_pipe &> $logfile_evgen";
       $conv_cmd      = "gntpc -f gst -i gntp.$curr_subrunnu.ghep.root | grep -B 100 -A 30 -i \"warn\\|error\\|fatal\" &> $logfile_conv";

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
       print PBS "export GSEED=$curr_seed \n";
       print PBS "$evgen_cmd \n";
       print PBS "$conv_cmd \n";

       print "EXEC: $evgen_cmd \n";

       #
       # submit job
       #
       `qsub -q $queue $batch_script`;

    } #subruns

  }
}
