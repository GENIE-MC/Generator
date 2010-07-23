#!/usr/bin/perl

#-----------------------------------------------------------------------------------------------------------
# Submit jobs to generate data needed for validating GENIE's hadron transport model
#
# Syntax:
#   perl submit-vld_inuke.pl <options>
#
# Options:
#    --version       : GENIE version number
#    --run           : runs to submit (eg --run 1060680 / --run 1060680,1062000 / -run all)
#   [--model-enum]   : physics model enumeration, default: 0
#   [--nsubruns]     : number of subruns per run, default: 1
#   [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]   : production name, default: <model>_<version>
#   [--cycle]        : cycle in current production, default: 01
#   [--use-valgrind] : default: off
#   [--batch-system] : <PBS, LSF, none>, default: PBS
#   [--queue]        : default: prod
#   [--softw-topdir] : default: /opt/ppd/t2k/GENIE
#
# Tested at the RAL/PPD Tier2 PBS batch farm.
#
# Adapted from Steve Dytman's script
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#
#-----------------------------------------------------------------------------------------------------------
#
# EVENT SAMPLES:
#
# Run number key: PPTTTEEEEEMMxx
#
# PP    : probe  (10:pi+, 11:pi0, 12:pi-, 20:K+, 21:K0, 22:K-, 30:gamma, 40:p, 41:n)
# TTT   : target (006:C12, 008:O16, 013:Al27, 020:Ca40, 026:Fe56: 028:Ni58, 082:Pb208)
# EEEEE : kinetic energy (MeV)
# MM    : physics model enumeration, default: 00
# xx    : sub-run ID, 00-99, 100k events each
#
#.......................................................................................
# run number       |  init state      | kin energy   | req. statistics    | group of   | 
#                  |                  | (GeV)        | (num events)       | runs       |     
#.......................................................................................
# 4000600597MMxx   | p     + C12      |   0.597      | 1.0M               | amian      | 
# 1002600870MMxx   | pi+   + F56      |   0.870      | 0.5M               | iwamoto    | 
# 1202600870MMxx   | pi-   + F56      |   0.870      | 1.0M               | iwamoto    | 
# 1008200870MMxx   | pi+   + Pb208    |   0.870      | 1.0M               | iwamoto    | 
# 1002602100MMxx   | pi+   + F56      |   2.100      | 1.0M               | iwamoto    | 
# 1002800220MMxx   | pi+   + Ni58     |   0.220      | 1.0M               | mckeown    | 
# 1002600220MMxx   | pi+   + F56      |   0.220      | 1.0M               | mckeown    | 
# 4001300256MMxx   | p     + Al27     |   0.256      | 1.0M               | stamer     | 
# 4008200256MMxx   | p     + Pb208    |   0.256      | 1.0M               | stamer     | 
# 1000800114MMxx   | pi+   + O16      |   0.114      | 1.0M               | ingram     | 
# 1000800240MMxx   | pi+   + O16      |   0.240      | 1.0M               | ingram     | 
# 4000600800MMxx   | p     + C12      |   0.800      | 0.3M               | mcgill     | 
# 4002000800MMxx   | p     + Ca40     |   0.800      | 1.0M               | mcgill     | 
# 1008200220MMxx   | pi+   + Pb208    |   0.220      | 1.0M               | levenson   | 
# 1000600220MMxx   | pi+   + C12      |   0.220      | 1.0M               | levenson   | 
# 1002800220MMxx   | pi+   + Ni58     |   0.220      | 1.0M               | levenson   | 
#.......................................................................................
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--nsubruns')      { $nsubruns      = $ARGV[$iarg+1]; }
  if($_ eq '--run')           { $runnu         = $ARGV[$iarg+1]; }
  if($_ eq '--model-enum')    { $model_enum    = $ARGV[$iarg+1]; }
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')  { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir  = $ARGV[$iarg+1]; }  
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;
die("** Aborting [You need to specify which runs to submit. Use the --run option]")
unless defined $runnu;

$model_enum     = "0"                                     unless defined $model_enum;
$nsubruns       = 1                                       unless defined $nsubruns;
$use_valgrind   = 0                                       unless defined $use_valgrind;
$arch           = "SL5_64bit"                             unless defined $arch;
$production     = "$model_enum\_$genie_version"           unless defined $production;
$cycle          = "01"                                    unless defined $cycle;
$batch_system   = "PBS"                                   unless defined $batch_system;
$queue          = "prod"                                  unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/GENIE"                    unless defined $softw_topdir;
$time_limit     = "60:00:00";
$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$jobs_dir       = "$softw_topdir/scratch/vld\_nuclmod-$production\_$cycle";
$xspl_file      = "$softw_topdir/data/job_inputs/xspl/gxspl-emode-$genie_version.xml";
$mcseed         = 210921029;
$nev_per_subrun = 100000;

# inputs for event generation jobs
%evg_probepdg_hash = ( 
  '4000600597' => '2212',
  '1002600870' => '211',
  '1202600870' => '-211',
  '1008200870' => '211',
  '1002602100' => '211',
  '1002800220' => '211',
  '1002600220' => '211',
  '4001300256' => '2212',
  '4008200256' => '2212',
  '1000800114' => '211',
  '1000800240' => '211',
  '4000600800' => '2212',
  '4002000800' => '2212',
  '1008200220' => '211',
  '1000600220' => '211',
  '1002800220' => '211'
);
%evg_tgtpdg_hash = ( 
  '4000600597' => '1000060120',
  '1002600870' => '1000260560',
  '1202600870' => '1000260560',
  '1008200870' => '1000822080',
  '1002602100' => '1000260560',
  '1002800220' => '1000280580',
  '1002600220' => '1000260560',
  '4001300256' => '1000130270',
  '4008200256' => '1000822080',
  '1000800114' => '1000080160',
  '1000800240' => '1000080160',
  '4000600800' => '1000060120',
  '4002000800' => '1000200400',
  '1008200220' => '1000822080',
  '1000600220' => '1000060120',
  '1002800220' => '1000280580'
);
%evg_kinetic_energy_hash = ( 
  '4000600597' => '0.597',
  '1002600870' => '0.870',
  '1202600870' => '0.870',
  '1008200870' => '0.870',
  '1002602100' => '2.100',
  '1002800220' => '0.220',
  '1002600220' => '0.220',
  '4001300256' => '0.256',
  '4008200256' => '0.256',
  '1000800114' => '0.114',
  '1000800240' => '0.240',
  '4000600800' => '0.800',
  '4002000800' => '0.800',
  '1008200220' => '0.220',
  '1000600220' => '0.220',
  '1002800220' => '0.220',
);
%group_hash = ( 
  '4000600597' => 'amian',
  '1002600870' => 'iwamoto',
  '1202600870' => 'iwamoto',
  '1008200870' => 'iwamoto',
  '1002602100' => 'iwamoto',
  '1002800220' => 'mckeown',
  '1002600220' => 'mckeown',
  '4001300256' => 'stamer',
  '4008200256' => 'stamer',
  '1000800114' => 'ingram',
  '1000800240' => 'ingram',
  '4000600800' => 'mcgill',
  '4002000800' => 'mcgill',
  '1008200220' => 'levenson',
  '1000600220' => 'levenson',
  '1002800220' => 'levenson',
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

#
# submit event generation jobs
#

# run loop
for my $curr_runnu (keys %evg_gevgl_hash)  {

 # check whether to commit current run 
 if($runnu=~m/$curr_runnu/ || $runnu eq "all") {

    print "** submitting event generation run: $curr_runnu \n";

    #
    # get runnu-dependent info
    #
    $probe   = $evg_probepdg_hash       {$curr_runnu};
    $tgt     = $evg_tgtpdg_hash         {$curr_runnu};
    $ke      = $evg_kinetic_energy_hash {$curr_runnu};

    # submit subruns
    for($isubrun = 0; $isubrun < $nsubruns; $isubrun++) {

       # Run number key: PPTTTEEEEEMMxx
       $curr_subrunnu = 10000 * $curr_runnu + 100 * $model_enum + $isubrun;
       $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $fntemplate    = "$jobs_dir/inuke-$curr_subrunnu";
       $curr_seed     = $mcseed + $isubrun;
       $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $evgen_cmd     = "gtestIntranuke -n $nev_per_subrun -k $ke -p $probe -t $tgt -r $curr_subrunnu | $grep_pipe &> $fntemplate.evgen.log";
       $conv_cmd      = "gntpc -f ginuke -i gntp.$curr_subrunnu.ghep.root | $grep_pipe &> $fntemplate.conv.log";

       print "@@ exec: $evgen_cmd \n";

       #
       # submit
       #
  
       # PBS case
       if($batch_system eq 'PBS') {
          $batch_script  = "$fntemplate.pbs";
          open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
          print PBS "#!/bin/bash \n";
          print PBS "#PBS -N nucl-$curr_subrunnu \n";
          print PBS "#PBS -l cput=$time_limit \n";
          print PBS "#PBS -o $fntemplate.pbsout.log \n";
          print PBS "#PBS -e $fntemplate.pbserr.log \n";
          print PBS "source $genie_setup \n"; 
          print PBS "cd $jobs_dir \n";
          print PBS "export GSEED=$curr_seed \n";
          print PBS "$evgen_cmd \n";
          print PBS "$conv_cmd \n";
          close(PBS);
          `qsub -q $queue $batch_script`;
       } # PBS

       # LSF case
       if($batch_system eq 'LSF') {
          $batch_script  = "$fntemplate.sh";
          open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
          print LSF "#!/bin/bash \n";
          print LSF "#BSUB-j nucl-$curr_subrunnu \n";
          print LSF "#BSUB-q $queue \n";
          print LSF "#BSUB-c $time_limit \n";
          print LSF "#BSUB-o $fntemplate.lsfout.log \n";
          print LSF "#BSUB-e $fntemplate.lsferr.log \n";
          print LSF "source $genie_setup \n"; 
          print LSF "cd $jobs_dir \n";
          print LSF "export GSEED=$curr_seed \n";
          print LSF "$evgen_cmd \n";
          print LSF "$conv_cmd \n";
          close(LSF);
          `bsub < $batch_script`;
       } # LSF

       # no batch system, run jobs interactively
       if($batch_system eq 'none') {
          system("source $genie_setup; cd $jobs_dir; export GSEED=$curr_seed; $evgen_cmd; $conv_cmd");
       } # interactive mode

    } # loop over subruns

 } #checking whether to submit current run
} # loop over runs

