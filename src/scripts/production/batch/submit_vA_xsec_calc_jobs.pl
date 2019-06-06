#!/usr/bin/perl

#----------------------------------------------------------------------------------------------------------------
# Submit jobs for calculating GENIE x-section splines for all nuclear targets and at the energy range required 
# for generating the GENIE release validation samples. Note that other scripts are available for generating 
# x-section splines for the larger array of targets found in detector geometry descriptions of given expts.
# If needed, use the GENIE gspladd utility to merge the job outputs.
#
# Syntax:
#   shell% perl submit_vA_xsec_calc_jobs.pl <options>
#
# Options:
#    --version       : GENIE version number
#    --config-file   : Text file which specifies the list of initial states
#   [--arch]         : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]   : default: routine_validation
#   [--cycle]        : default: 01
#   [--use-valgrind] : default: off
#   [--batch-system]  : <PBS, LSF, slurm, HTCondor, HTCondor_PBS, none>, default: PBS
#   [--queue]        : default: prod
#   [--softw-topdir]  : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]   : top level dir for job files, default: /opt/ppd/t2k/scratch/GENIE/
#
# Author:
#   Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#----------------------------------------------------------------------------------------------------------------

use File::Path;
use File::Basename;

# inputs
#  
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')        { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--config-file')    { $config_file   = $ARGV[$iarg+1]; }
  if($_ eq '--arch')           { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')     { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')          { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')   { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')   { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')          { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')   { $softw_topdir  = $ARGV[$iarg+1]; }
  if($_ eq '--jobs-topdir')    { $jobs_topdir   = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;
die("** Aborting [Undefined configuration file. Use the --config-file option. See examples in ?]")
unless defined $config_file;

$use_valgrind   = 0                             unless defined $use_valgrind;
$arch           = "SL6.x86_64"                  unless defined $arch;
$production     = "routine_validation"          unless defined $production;
$cycle          = "01"                          unless defined $cycle;
$batch_system   = "PBS"                         unless defined $batch_system;
$queue          = "prod"                        unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE/"   unless defined $softw_topdir;
$jobs_topdir    = "/opt/ppd/t2k/scratch/GENIE/" unless defined $jobs_topdir;
$time_limit     = "60:00:00";

$genie_setup    = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$freenucsplines = "$softw_topdir/data/job_inputs/xspl/gxspl-vN-$genie_version.xml";
$config         = basename($config_file,".list");
$jobs_dir       = "$jobs_topdir/$genie_version-$production\_$cycle-xsec\_vA\_$config";

$nknots    =  0;
$emax      =  0;
@neutrinos = ();
@targets   = ();

open (CONFIG, $config_file);
$i=0;
while(<CONFIG>)
{
   chomp;
   # skip comment lines starting with # and empty lines that don't contain any number
   if( substr($_,0,1) ne '#' && m/[0-9]/ ) 
   { 
      s/ //g;  # rm empty spaces from $_
      #print "line = $_\n" ;
      $value = "";
      $idx = index($_,"#");
      if($idx == -1) { $value = $_; }
      else           { $value = substr($_,0,$idx); }

      if    ($i==0) { $emax   = $value; }
      elsif ($i==1) { $nknots = $value; }
      else {
          if( $value ==  12 || $value ==  14 || $value ==  16 ||
              $value == -12 || $value == -14 || $value == -16 ) 
          {
             push(@neutrinos, $value);
          }
          else 
          {
             push(@targets, $value);
          }
      }
      $i++;
   }
}
close(CONFIG);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

$nu_codes = "";
$i=0;
foreach (@neutrinos) 
{
    s/ //g;  # rm empty spaces from $_
    if($i == 0) { $nu_codes = $_; }
    else        { $nu_codes = $nu_codes . ",$_"; }
    $i++;
}

#
# loop over nuclear targets & submit jobs
#
foreach(@targets) 
{
    s/ //g;  # rm empty spaces from $_
    $tgt_code = $_;
    $jobname  = "vAxscalc-$tgt_code";
    $filename_template = "$jobs_dir/$jobname";
    $grep_pipe  = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
    $gmkspl_opt = "-p $nu_codes -t $_ -n $nknots -e $emax --input-cross-sections $freenucsplines --output-cross-sections gxspl_$tgt_code.xml";
    $gmkspl_cmd = "gmkspl $gmkspl_opt | $grep_pipe &> $filename_template.mkspl.log";
    print "@@ exec: $gmkspl_cmd \n";

    #
    # submit
    #
  
    # PBS case
    if($batch_system eq 'PBS' || $batch_system eq 'HTCondor_PBS') {
	$batch_script = "$filename_template.pbs";
	open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
	print PBS "#!/bin/bash \n";
        print PBS "#PBS -N $jobname \n";
        print PBS "#PBS -o $filename_template.pbsout.log \n";
        print PBS "#PBS -e $filename_template.pbserr.log \n";
	print PBS "source $genie_setup \n";
	print PBS "cd $jobs_dir \n";
	print PBS "$gmkspl_cmd \n";
        close(PBS);
        $job_submission_command = "qsub";
        if($batch_system eq 'HTCondor_PBS') {
           $job_submission_command = "condor_qsub";
        }
        `$job_submission_command -q $queue $batch_script`;
    } #PBS

    # LSF case
    if($batch_system eq 'LSF') {
	$batch_script = "$filename_template.sh";
	open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
	print LSF "#!/bin/bash \n";
        print LSF "#BSUB-j $jobname \n";
        print LSF "#BSUB-o $filename_template.lsfout.log \n";
        print LSF "#BSUB-e $filename_template.lsferr.log \n";
	print LSF "source $genie_setup \n";
	print LSF "cd $jobs_dir \n";
	print LSF "$gmkspl_cmd \n";
        close(LSF);
	`bsub < $batch_script`;
    } #LSF

    # HTCondor
    if($batch_system eq 'HTCondor') {
        $batch_script = "$filename_template.htc";
        open(HTC, ">$batch_script") or die("Can not create the Condor submit description file: $batch_script");
        print HTC "Universe               = vanilla \n";
        print HTC "Executable             = $softw_topdir/generator/builds/$arch/$genie_version/src/scripts/production/batch/htcondor_exec.sh \n";
        print HTC "Arguments              = $genie_setup $jobs_dir $gmkspl_cmd \n";
        print HTC "Log                    = $filename_template.log \n";
        print HTC "Output                 = $filename_template.out \n";
        print HTC "Error                  = $filename_template.err \n";
        print HTC "Request_memory         = 2 GB \n";
        print HTC "Queue \n";
        close(HTC);
        `condor_submit $batch_script`;
    } #HTCondor
    
    # slurm case
    if($batch_system eq 'slurm') {
        my $time_lim = `sinfo -h -p batch -o %l`;
        my ($days, $hours, $remainder) = $time_lim =~ /([0]+)-([0-9]+):(.*)/;
        my $newhours = $days * 24 + $hours;
        my $new_time_lim = "$newhours:$remainder";
        $time_limit = $new_time_lim lt $time_limit ? $new_time_lim : $time_limit;
	$batch_script = "$filename_template.sh";
	open(SLURM, ">$batch_script") or die("Can not create the SLURM batch script");
	print SLURM "#!/bin/bash \n";
        print SLURM "#SBATCH-p $queue \n";
        print SLURM "#SBATCH-o $filename_template.lsfout.log \n";
        print SLURM "#SBATCH-e $filename_template.lsferr.log \n";
        print SLURM "#SBATCH-t $time_limit \n";
	print SLURM "source $genie_setup \n";
	print SLURM "cd $jobs_dir \n";
	print SLURM "$gmkspl_cmd \n";
        close(SLURM);
	`sbatch --job-name=$jobname $batch_script`;
    } #slurm

    # run interactively
    if($batch_system eq 'none') {
        system("source $genie_setup; cd $jobs_dir; $gmkspl_cmd");
    }
}

