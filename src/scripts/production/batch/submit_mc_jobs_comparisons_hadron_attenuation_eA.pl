#!/usr/bin/perl

#-----------------------------------------------------------------------------------------------------------
# Submit jobs to generate data needed for validating GENIE's hadronization + intranuclear transport model 
# (with emphasis  on the modeling of medium effects to the hadronization) in a set of comparisons against 
# the CLAS and HERMES nuclear attenuation data.
#
# Syntax:
#   perl submit_mc_jobs_comparisons_hadron_attenuation_eA.pl <options>
#
# Options:
#    --gen-version    : GENIE version number
#    --tune           : GENIE physics tune
#    --run            : runs to submit (eg --run 101102 / --run 101102,154002 / -run all)
#   [--model-enum]    : physics model enumeration, default: 0
#   [--nsubruns]      : number of subruns per run, default: 1
#   [--offset]        : subrun offset (for augmenting existing sample), default: 0
#   [--arch]          : <el7.x86_64, ...>, default: el7.x86_64
#   [--production]    : production name, default: <version>
#   [--cycle]         : cycle in current production, default: 01
#   [--use-valgrind]  : default: off
#   [--batch-system]  : <Slurm, PBS, LSF, none>, default: Slurm
#   [--queue]         : default: compute
#   [--softw-topdir]  : top level dir for softw installations, default: /user/costasa/projects/GENIE/softw/
#   [--jobs-topdir]   : top level dir for job files, default: /scratch/costasa/GENIE/
#
# EVENT SAMPLES:
#
# Run number key: ITTJJMxxx
# I   :  probe (1:e-, 2:e+)
# TT  :  nuclear target (01:H1, 02:D2, 03:He4, 06:C12, 07:N14, 08:O16, 10:Ne20, 26:Fe56, 36:Kr83, 54:Xe131)
# JJ  :  flux setting (01: 12.0 GeV [HERMES], 02: 27.6 GeV [HERMES])
# M   :  model enumeration
# xxx :  sub-run ID, 000-999, 50k events each
#
#...................................................................................
# run number     |  init state      | energy   | event gen.    | flux
#                |                  | (GeV)    | list          |
#...................................................................................
# 10202Mxxx      | e-    + D2       | 27.6     | EM            | monoenergetic
# 10302Mxxx      | e-    + He4      | 27.6     | EM            | monoenergetic
# 10702Mxxx      | e-    + N14      | 27.6     | EM            | monoenergetic
# 11002Mxxx      | e-    + Ne20     | 27.6     | EM            | monoenergetic
# 13602Mxxx      | e-    + Kr84     | 27.6     | EM            | monoenergetic
# 15402Mxxx      | e-    + Xe132    | 27.6     | EM            | monoenergetic
#
#
# Author:
#   Costas Andreopoulos <c.andreopoulos \st cern.ch>
#   University of Liverpool
#
# Copyright:
#   Copyright (c) 2003-2025, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#-----------------------------------------------------------------------------------------------------------
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--gen-version')   { $gen_version   = $ARGV[$iarg+1]; }
  if($_ eq '--tune')          { $tune          = $ARGV[$iarg+1]; }
  if($_ eq '--run')           { $runnu         = $ARGV[$iarg+1]; }
  if($_ eq '--model-enum')    { $model_enum    = $ARGV[$iarg+1]; }
  if($_ eq '--nsubruns')      { $nsubruns      = $ARGV[$iarg+1]; }
  if($_ eq '--offset')        { $offset        = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')  { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir  = $ARGV[$iarg+1]; }  
  if($_ eq '--jobs-topdir')   { $jobs_topdir   = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined GENIE Generator version. Use the --gen-version option]")
unless defined $gen_version;
die("** Aborting [Undefined GENIE physics tune. Use the --tune option]")
unless defined $tune;
die("** Aborting [You need to specify which runs to submit. Use the --run option]")
unless defined $runnu;

$model_enum       = "0"                                   unless defined $model_enum;
$nsubruns         = 1                                     unless defined $nsubruns;
$offset           = 0                                     unless defined $offset;
$use_valgrind     = 0                                     unless defined $use_valgrind;
$arch             = "el7.x86_64"                          unless defined $arch;
$production       = "prod"                                unless defined $production;
$cycle            = "01"                                  unless defined $cycle;
$batch_system     = "Slurm"                               unless defined $batch_system;
$queue            = "compute"                             unless defined $queue;
$softw_topdir     = "/user/costasa/projects/GENIE/softw/" unless defined $softw_topdir;
$jobs_topdir      = "/scratch/costasa/GENIE/"             unless defined $jobs_topdir;
$time_limit       = "10:00:00";
$gen_setup_script = "$softw_topdir/generator/builds/$arch/$gen_version-setup.sh";
$jobs_dir         = "$jobs_topdir/$gen_version-$tune-$production\_$cycle-hadroatten\_e";
$xspl_file        = "$softw_topdir/data/job_inputs/xspl/gxspl-eA-$genie_version.xml";
$mcseed           = 210921029;
$nev_per_subrun = 50000;

# inputs for event generation jobs
%evg_pdg_hash = ( 
  '10202' =>   '11',
  '10302' =>   '11',
  '10702' =>   '11',
  '11002' =>   '11',
  '13602' =>   '11',
  '15402' =>   '11'
);
%evg_tgtpdg_hash = ( 
  '10202' =>   '1000010020',
  '10302' =>   '1000020040',
  '10702' =>   '1000070140',
  '11002' =>   '1000100200',
  '13602' =>   '1000360840',
  '15402' =>   '1000541320'
);
%evg_energy_hash = ( 
  '10202' =>   '27.6',
  '10302' =>   '27.6',
  '10702' =>   '27.6',
  '11002' =>   '27.6',
  '13602' =>   '27.6',
  '15402' =>   '27.6'
);
%evg_gevgl_hash = ( 
  '10202' =>   'EM',
  '10302' =>   'EM',
  '10702' =>   'EM',
  '11002' =>   'EM',
  '13602' =>   'EM',
  '15402' =>   'EM'
);
%evg_fluxopt_hash = ( 
  '10202' =>   '',
  '10302' =>   '',
  '10702' =>   '',
  '11002' =>   '',
  '13602' =>   '',
  '15402' =>   ''
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

    # Get info that depends on the run number 
    # -------------------------------------------------------------

    $probe   = $evg_pdg_hash     {$curr_runnu};
    $tgt     = $evg_tgtpdg_hash  {$curr_runnu};
    $en      = $evg_energy_hash  {$curr_runnu};
    $gevgl   = $evg_gevgl_hash   {$curr_runnu};
    $fluxopt = $evg_fluxopt_hash {$curr_runnu};

    # loop over subruns
    for($isubrun = 0; $isubrun < $nsubruns; $isubrun++) {

       # Get info that depends on the run and subrun numbers
       # -------------------------------------------------------------
       # Run number key: ITTJJMxxx
       $curr_subrunnu     = 10000 * $curr_runnu + 1000 * $model_enum + $isubrun + $offset;
       $curr_seed         = $mcseed + $isubrun + $offset;
       $job_name          = "hadroatten\_e-$curr_subrunnu"; 
       $filename_basepath = "$jobs_dir/$job_name";

       # Form commands to run
       # -------------------------------------------------------------
#      $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $grep_pipe     = "grep -B 20 -A 30 -i fatal";
       $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $evgen_opt     = "-n $nev_per_subrun -e $en -p $probe -t $tgt $fluxopt --tune $tune -r $curr_subrunnu --seed $curr_seed --cross-sections $xspl_file --event-generator-list $gevgl";
       $evgen_cmd     = "gevgen $evgen_opt | $grep_pipe &> $fntemplate.evgen.log";
       $conv_cmd      = "gntpc -f gst -i gntp.$curr_subrunnu.ghep.root";

       print "@@ exec: $evgen_cmd \n";

       # Submit jobs
       # -------------------------------------------------------------

       # Slurm case
       if($batch_system eq 'slurm') {
           $batch_script  = "$filename_basepath.sh";
           open(SLURM, ">$batch_script") or die("Can not create the slurm batch script");
           print SLURM "#!/bin/bash \n";
           print SLURM "#SBATCH-p $queue \n";
           print SLURM "#SBATCH-J $job_name \n";
           print SLURM "#SBATCH-N 1 \n";
           print SLURM "#SBATCH-c 1 \n";
           print SLURM "#SBATCH-o $filename_basepath.slurmout.log \n";
           print SLURM "#SBATCH-e $filename_basepath.slurmerr.log \n";
           print SLURM "#SBATCH-t $time_limit \n";
           print SLURM "source $gen_setup_script \n";
           print SLURM "cd $jobs_dir \n";
           print SLURM "$evgen_cmd \n";
           print SLURM "$conv_cmd \n";
           close(SLURM);
           `sbatch $batch_script`;
       } #Slurm

       # PBS case
       if($batch_system eq 'PBS') {
           $batch_script  = "$filename_basepath.pbs";
           open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
           print PBS "#!/bin/bash \n";
           print PBS "#PBS -N $job_name \n";
           print PBS "#PBS -l cput=$time_limit \n";
           print PBS "#PBS -o $filename_basepath.pbsout.log \n";
           print PBS "#PBS -e $filename_basepath.pbserr.log \n";
           print PBS "source $gen_setup_script \n"; 
           print PBS "cd $jobs_dir \n";
           print PBS "$evgen_cmd \n";
           print PBS "$conv_cmd \n";
           close(PBS);
           `qsub -q $queue $batch_script`;
       } #PBS

       # LSF case
       if($batch_system eq 'LSF') {
           $batch_script  = "$filename_basepath.sh";
           open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
           print LSF "#!/bin/bash \n";
           print LSF "#BSUB-j $job_name \n";
           print LSF "#BSUB-q $queue \n";
           print LSF "#BSUB-c $time_limit \n";
           print LSF "#BSUB-o $filename_basepath.lsfout.log \n";
           print LSF "#BSUB-e $filename_basepath.lsferr.log \n";
           print LSF "source $gen_setup_script \n"; 
           print LSF "cd $jobs_dir \n";
           print LSF "$evgen_cmd \n";
           print LSF "$conv_cmd \n";
           close(LSF);
           `qsub < $batch_script`;
       } #LSF

       # no batch system, run jobs interactively
       if($batch_system eq 'none') {
          system("source $gen_setup_script; cd $jobs_dir; $evgen_cmd; $conv_cmd");
       } # interactive mode


    } # loop over subruns
  } #checking whether to submit current run
} # loop over runs

