#!/usr/bin/perl

#---------------------------------------------------------------------------------------------------------------------
# Submit jobs to generate all data needed for validating GENIE's hadronization model.
# The generated data can be fed into GENIE's `gvld_hadronz_test' utility.
#
# Syntax:
#   perl submit_hadronization_validation_mc_jobs.pl <options>
#
# Options:
#    --version       : GENIE version number
#   [--nsubruns]     : number of subruns per run, default: 1
#   [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]   : production name, default: routine_validation
#   [--cycle]        : cycle in current production, default: 01
#   [--use-valgrind] : default: off
#   [--batch-system] : <PBS, LSF, slurm>, default: PBS
#   [--queue]        : default: prod
#   [--softw-topdir] : default: /opt/ppd/t2k/softw/GENIE
#
# Tested at the RAL/PPD Tier2 PBS batch farm.
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#---------------------------------------------------------------------------------------------------------------------
#
# EVENT SAMPLES:
#...................................................................................
# run number      |  init state      | energy   | event generator    | flux
#                 |                  | (GeV)    | list               |
#...................................................................................
#
# 1000xx          | numu    + n      | 0.5-80.  | HadronizationTest  | 1/E
# 1100xx          | numu    + p      | 0.5-80.  | HadronizationTest  | 1/E
# 1200xx          | numubar + n      | 0.5-80.  | HadronizationTest  | 1/E
# 1300xx          | numubar + p      | 0.5-80.  | HadronizationTest  | 1/E
#
# xx : Run ID, 01-99, 100k events each
#...................................................................................
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--nsubruns')      { $nsubruns      = $ARGV[$iarg+1]; }
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

$nsubruns       = 1                           unless defined $nsubruns;
$use_valgrind   = 0                           unless defined $use_valgrind;
$arch           = "SL5_64bit"                 unless defined $arch;
$production     = "routine_validation"        unless defined $production;
$cycle          = "01"                        unless defined $cycle;
$batch_system   = "PBS"                       unless defined $batch_system;
$queue          = "prod"                      unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE"  unless defined $softw_topdir;
$time_limit     = "60:00:00";
$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$jobs_dir       = "$softw_topdir/scratch/$genie_version-$production\_$cycle-hadronization";
$xspl_file      = "$softw_topdir/data/job_inputs/xspl/gxspl-vN-$genie_version.xml";
$mcseed         = 210921029;
$nev_per_subrun = 100000;

# inputs for event generation jobs
%evg_nupdg_hash = ( 
  '1000' =>   '14',
  '1100' =>   '14',
  '1200' =>  '-14',
  '1300' =>  '-14'
);
%evg_tgtpdg_hash = ( 
  '1000' =>  '1000000010',
  '1100' =>  '1000010010',
  '1200' =>  '1000000010',
  '1300' =>  '1000010010'
);
%evg_energy_hash = ( 
  '1000' =>  '0.5,80.0',
  '1100' =>  '0.5,80.0',
  '1200' =>  '0.5,80.0',
  '1300' =>  '0.5,80.0',
);
%evg_gevgl_hash = ( 
  '1000' =>  'HadronizationTest',
  '1100' =>  'HadronizationTest',
  '1200' =>  'HadronizationTest',
  '1300' =>  'HadronizationTest',
);
%evg_fluxopt_hash = ( 
  '1000' =>  '-f "1/x"',
  '1100' =>  '-f "1/x"',
  '1200' =>  '-f "1/x"',
  '1300' =>  '-f "1/x"',
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

# submit event generation jobs
#
for my $curr_runnu (keys %evg_gevgl_hash)  {

 # uncomment if you want to include a cmd line input to specify specific runs (and fill-in $runnu)
 # if($runnu=~m/$curr_runnu/ || $runnu eq "all") {

    print "** submitting event generation run: $curr_runnu \n";

    #
    # get runnu-dependent info
    #
    $nu      = $evg_nupdg_hash   {$curr_runnu};
    $tgt     = $evg_tgtpdg_hash  {$curr_runnu};
    $en      = $evg_energy_hash  {$curr_runnu};
    $gevgl   = $evg_gevgl_hash   {$curr_runnu};
    $fluxopt = $evg_fluxopt_hash {$curr_runnu};

    # submit subruns
    for($isubrun = 0; $isubrun < $nsubruns; $isubrun++) {

       $curr_subrunnu = 100 * $curr_runnu + $isubrun;
       $curr_seed     = $mcseed + $isubrun;
       $jntemplate    = "hadro-$curr_subrunnu";
       $fntemplate    = "$jobs_dir/$jntemplate";
       $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $evgen_opt     = "-n $nev_per_subrun -e $en -p $nu -t $tgt $fluxopt -r $curr_subrunnu --seed $curr_seed --cross-sections $xspl_file --event-generator-list $gevgl";
       $evgen_cmd     = "gevgen $evgen_opt | $grep_pipe &> $fntemplate.evgen.log";

       print "@@ exec: $evgen_cmd \n";

       #
       # submit
       #
  
       # PBS case
       if($batch_system eq 'PBS') {
           $batch_script = "$fntemplate.pbs";
           open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
           print PBS "#!/bin/bash \n";
           print PBS "#PBS -N $jntemplate \n";
           print PBS "#PBS -l cput=$time_limit \n";
           print PBS "#PBS -o $fntemplate.pbsout.log \n";
           print PBS "#PBS -e $fntemplate.pbserr.log \n";
           print PBS "source $genie_setup \n"; 
           print PBS "cd $jobs_dir \n";
           print PBS "$evgen_cmd \n";
           close(PBS);
           `qsub -q $queue $batch_script`;
       }#PBS

       # LSF case
       if($batch_system eq 'LSF') {
           $batch_script = "$fntemplate.sh";
           open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
           print LSF "#!/bin/bash \n";
           print PBS "#BSUB-j $jntemplate \n";
           print LSF "#BSUB-q $queue \n";
           print LSF "#BSUB-c $time_limit \n";
           print LSF "#BSUB-o $fntemplate.lsfout.log \n";
           print LSF "#BSUB-e $fntemplate.lsferr.log \n";
           print LSF "source $genie_setup \n"; 
           print LSF "cd $jobs_dir \n";
           print LSF "$evgen_cmd \n";
           close(LSF);
           `bsub < $batch_script`;
       }#LSF

       # slurm case
       if($batch_system eq 'slurm') {
           $batch_script = "$fntemplate.sh";
           open(SLURM, ">$batch_script") or die("Can not create the SLURM batch script");
           print SLURM "#!/bin/bash \n";
           print SLURM "#SBATCH-p $queue \n";
           print SLURM "#SBATCH-o $fntemplate.lsfout.log \n";
           print SLURM "#SBATCH-e $fntemplate.lsferr.log \n";
           print SLURM "source $genie_setup \n"; 
           print SLURM "cd $jobs_dir \n";
           print SLURM "$evgen_cmd \n";
           close(SLURM);
           `sbatch --job-name=$jntemplate $batch_script`;
       }#slurm

    } # loop over subruns
 # } #checking whether to submit current run
} # loop over runs

