#!/usr/bin/perl

#---------------------------------------------------------------------------------------------------------------------
# Submit jobs for validating GENIE against the MiniBooNE data releases.
# The outputs can be fed to appropriate apps of the GENIE/Comparisons product
# 
# Syntax:
#   shell% perl submit_miniboone_validation_mc_jobs.pl <options>
#
# Options:
#    --version        : GENIE version number
#    --run-type       : Comma separated list of run types
#   [--config-dir]    : Config directory, default is $GENIE/config 
#   [--nsubruns]      : number of subruns per run, default: 1
#   [--arch]          : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]    : default: routine_validation
#   [--cycle]         : default: 01
#   [--batch-system]  : <PBS, LyonPBS, LSF, slurm, HTCondor, HTCondor_PBS, none>, default: HTCondor_PBS
#   [--queue]         : default: prod
#   [--softw-topdir]  : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]   : top level dir for job files, default: /opt/ppd/t2k/scratch/GENIE/
#   [--spline-file]   : absoluyte path to xsec_spline_file, default: $softw_topdir/data/job_inputs/xspl/gxspl-vA-$genie_version.xml
#
#
# SAMPLES:
#........................................................................
#  run   | nev      |  init state             | energy     |   flux
#  type  | /subrun  |                         |            |
#........................................................................
#  10   | 100k     |  numu    + Hydrocarbon  | 0 - 10 GeV | (a)
#  20   | 100k     |  numubar + Hydrocarbon  | 0 - 10 GeV | (b)
#........................................................................
#
# Run numbers are constructed as TTXXX, where  TT is the run type and XXX is the subrun number
#
# Hydrocarbon: 85.7142% C12 + 14.2857% H1
# (a) april07_baseline_rgen610.6_flux_pospolarity_fluxes (flux_pos_pol_numu)
# (b) december2007_horn-174ka_rgen610.6_flux_negpolarity_fluxes (flux_neg_pol_numub)
#
# Author:
#   Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#---------------------------------------------------------------------------------------------------------------------
#
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')        { $genie_version   = $ARGV[$iarg+1]; }
  if($_ eq '--config-dir')     { $config_dir      = $ARGV[$iarg+1]; }
  if($_ eq '--run-type')       { $run_type        = $ARGV[$iarg+1]; }
  if($_ eq '--nsubruns')       { $nsubruns        = $ARGV[$iarg+1]; }
  if($_ eq '--arch')           { $arch            = $ARGV[$iarg+1]; }
  if($_ eq '--production')     { $production      = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')          { $cycle           = $ARGV[$iarg+1]; }
  if($_ eq '--ref-samples')    { $ref_sample_path = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')   { $use_valgrind    = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')   { $batch_system    = $ARGV[$iarg+1]; }
  if($_ eq '--queue')          { $queue           = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')   { $softw_topdir    = $ARGV[$iarg+1]; }  
  if($_ eq '--jobs-topdir')    { $jobs_topdir     = $ARGV[$iarg+1]; }
  if($_ eq '--spline-file')    { $xspl_file       = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined benchmark runs #. Use the --run-type option]")
unless defined $run_type;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$config_dir     = ""                            unless defined $config_dir;
$nsubruns        = 1                            unless defined $nsubruns;
$use_valgrind    = 0                            unless defined $use_valgrind;
$arch            = "SL6.x86_64"                 unless defined $arch;
$production      = "routine_validation"         unless defined $production;
$cycle           = "01"                         unless defined $cycle;
$batch_system    = "HTCondor_PBS"               unless defined $batch_system;
$queue           = "prod"                       unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE/"   unless defined $softw_topdir;
$jobs_topdir    = "/opt/ppd/t2k/scratch/GENIE/" unless defined $jobs_topdir;  
$ref_sample_path = 0                            unless defined $ref_sample_path;
$genie_setup     = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$xspl_file       = "$softw_topdir/data/job_inputs/xspl/gxspl-vA-$genie_version.xml" unless defined $xspl_file ;
$jobs_dir        = "$jobs_topdir/$genie_version-$production\_$cycle-miniboone";
$mcseed          = 210921029;

%nevents_hash = ( 
  '10' => '100000',
  '20' => '100000'
);

%nupdg_hash = ( 
  '10' =>   '14',
  '20' =>  '-14'
);

%tgtpdg_hash = ( 
  '10' => '1000060120[0.857142],1000010010[0.142857]',
  '20' => '1000060120[0.857142],1000010010[0.142857]'
);

%energy_hash = ( 
  '10' =>  '0,10',
  '20' =>  '0,10'
);

%flux_hash = ( 
  '10' =>  "$softw_topdir/comparisons/builds/$arch/vtrunk/data/fluxes/miniboone/miniboone_april07_baseline_rgen610.6_flux_pospolarity_fluxes.root,flux_pos_pol_numu",
  '20' =>  "$softw_topdir/comparisons/builds/$arch/vtrunk/data/fluxes/miniboone/miniboone_december2007_horn-174ka_rgen610.6_flux_negpolarity_fluxes.root,flux_neg_pol_numub"
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

print "Input runs: $run_type \n";

for my $curr_runtype (keys %nupdg_hash)  {
  print "Checking benchmark run: ...... $curr_runtype \n";

  if($run_type=~m/$curr_runtype/ || $run_type eq "all") {
    print "** matched -> submitting job \n";

    # get run type dependent info
    $nev   = $nevents_hash {$curr_runtype};
    $nu    = $nupdg_hash   {$curr_runtype};
    $tgt   = $tgtpdg_hash  {$curr_runtype};
    $en    = $energy_hash  {$curr_runtype};
    $flux  = $flux_hash    {$curr_runtype};

    # submit subruns
    for($subrunnu = 0; $subrunnu < $nsubruns; $subrunnu++) {

       $runnu   = 1000 * $curr_runtype + $subrunnu;
       $jobname = "miniboone-$runnu";
       $filename_template = "$jobs_dir/$jobname";
#      $grep_pipe         = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $grep_pipe         = "grep -B 20 -A 30 -i fatal";
       $valgrind_cmd      = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $curr_seed         = $mcseed + $runnu;
       $evgen_opt         = "-n $nev -p $nu -t $tgt -r $runnu --seed $curr_seed --cross-sections $xspl_file -e $en -f $flux";
       $evgen_cmd         = "gevgen $evgen_opt | $grep_pipe &> $filename_template.evgen.log";

       print "@@ exec: $evgen_cmd \n";

       #
       # submit
       #
  
       # PBS case
       if($batch_system eq 'PBS' || $batch_system eq 'HTCondor_PBS') {
          $batch_script  = "$filename_template.pbs";
          open(BATCH_SCRIPT, ">$batch_script") or die("Can not create the PBS batch script");
          print BATCH_SCRIPT "#!/bin/bash \n";
          print BATCH_SCRIPT "#PBS -N $jobname \n";
          print BATCH_SCRIPT "#PBS -o $filename_template.pbsout.log \n";
          print BATCH_SCRIPT "#PBS -e $filename_template.pbserr.log \n";
          print BATCH_SCRIPT "source $genie_setup $config_dir \n"; 
          print BATCH_SCRIPT "cd $jobs_dir \n";
          print BATCH_SCRIPT "$evgen_cmd \n";
          close(BATCH_SCRIPT);
          $job_submission_command = "qsub";
          if($batch_system eq 'HTCondor_PBS') {
            $job_submission_command = "condor_qsub";
          }
          `$job_submission_command -q $queue $batch_script`;
       } #PBS

       # LyonPBS case
       if($batch_system eq 'LyonPBS' ) {
         $batch_script = "$filename_template.pbs";
         open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
         print PBS "#!/bin/bash \n";
         print PBS "#\$ -P P_$ENV{'GROUP'} \n";
         print PBS "#\$ -N $jobname \n";
         print PBS "#\$ -o $filename_template.pbsout.log \n";
         print PBS "#\$ -e $filename_template.pbserr.log \n";
         print PBS "#\$ -l ct=6:00:00,sps=1 \n";
         print PBS "source $genie_setup $config_dir \n";
         print PBS "cd $jobs_dir \n";
         print PBS "$evgen_cmd \n";
         close(PBS);
         $job_submission_command = "qsub";
         `$job_submission_command  $batch_script`;
       } #LyonPBS 

       # LSF case
       if($batch_system eq 'LSF') {
          $batch_script  = "$filename_template.sh";
          open(BATCH_SCRIPT, ">$batch_script") or die("Can not create the LSF batch script");
          print BATCH_SCRIPT "#!/bin/bash \n";
          print BATCH_SCRIPT "#BSUB-j $jobname \n";
          print BATCH_SCRIPT "#BSUB-q $queue \n";
          print BATCH_SCRIPT "#BSUB-o $filename_template.lsfout.log \n";
          print BATCH_SCRIPT "#BSUB-e $filename_template.lsferr.log \n";
          print BATCH_SCRIPT "source $genie_setup $config_dir\n"; 
          print BATCH_SCRIPT "cd $jobs_dir \n";
          print BATCH_SCRIPT "$evgen_cmd \n";
          close(BATCH_SCRIPT);
          `bsub < $batch_script`;
       } #LSF
             
       # HTCondor
       if($batch_system eq 'HTCondor') {
          $batch_script = "$filename_template.htc";
          open(BATCH_SCRIPT, ">$batch_script") or die("Can not create the Condor submit description file: $batch_script");
          print BATCH_SCRIPT "Executable               = $softw_topdir/generator/builds/$arch/$genie_version/src/scripts/production/batch/htcondor_exec.sh \n";
          print BATCH_SCRIPT "Arguments                = $genie_setup $jobs_dir $gevgen_cmd \n";
          print BATCH_SCRIPT "Log                      = $filename_template.log \n";
          print BATCH_SCRIPT "Output                   = $filename_template.out \n";
          print BATCH_SCRIPT "Error                    = $filename_template.err \n";
          print BATCH_SCRIPT "Universe                 = vanilla \n";
          print BATCH_SCRIPT "Request_memory           = 2 GB \n";
          print BATCH_SCRIPT "Getenv                   = True \n";
          print BATCH_SCRIPT "should_transfer_files    = YES \n";
          print BATCH_SCRIPT "when_to_transfer_output  = ON_EXIT \n";
          print BATCH_SCRIPT "\nqueue\n";    
          print BATCH_SCRIPT "Queue \n";
          close(BATCH_SCRIPT);
          `condor_submit $batch_script`;
       } #HTCondor

       # slurm case
       if($batch_system eq 'slurm') {
          $batch_script  = "$filename_template.sh";
          open(BATCH_SCRIPT, ">$batch_script") or die("Can not create the SLURM batch script");
          print BATCH_SCRIPT "#!/bin/bash \n";
          print BATCH_SCRIPT "#SBATCH-p $queue \n";
          print BATCH_SCRIPT "#SBATCH-o $filename_template.lsfout.log \n";
          print BATCH_SCRIPT "#SBATCH-e $filename_template.lsferr.log \n";
          print BATCH_SCRIPT "source $genie_setup $config_dir\n"; 
          print BATCH_SCRIPT "cd $jobs_dir \n";
          print BATCH_SCRIPT "$evgen_cmd \n";
          close(BATCH_SCRIPT);
          `sbatch --job-name=$jobname $batch_script`;
       } #slurm
 
       # no batch system, run jobs interactively
       if($batch_system eq 'none') {
          system("source $genie_setup $config_dir; cd $jobs_dir; $evgen_cmd");
       } # interactive mode

    }#subrunnu

  }
}
