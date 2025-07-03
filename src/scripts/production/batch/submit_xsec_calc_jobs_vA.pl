#!/usr/bin/perl

#----------------------------------------------------------------------------------------------------------------
# Submit jobs for calculating the specified GENIE nuclear cross-section splines.
# The GENIE gspladd utility can be used to merge the job outputs.
#
# Syntax:
#   shell% perl submit_xsec_calc_jobs_vA.pl <options>
#
# Options:
#    --gen-version       : GENIE generator version number
#    --tune              : GENIE physics tune
#    --spline-list       : Text file listing the nuclear cross-section splines to be generated.
#   [--input-splines]    : Input free-nucleon splines that can speed-up the calculation of nuclear splines
#   [--arch]             : <el7.x86_64, ...>, default: el7.x86_64
#   [--production]       : default: routine_validation
#   [--cycle]            : default: 01
#   [--use-valgrind]     : default: off
#   [--batch-system]     : <Slurm, PBS, LSF, none>, default: Slurm
#   [--queue]            : default: compute
#   [--time-limit]       : default: 10:00:00
#   [--softw-topdir]     : top level dir for softw installations, default: /user/costasa/projects/GENIE/softw/
#   [--jobs-topdir]      : top level dir for job files, default: /scratch/costasa/GENIE/
#
# Examples:
#   % perl submit_xsec_calc_jobs_vA.pl --gen-version v3.06.00 \
#                  --tune G18_10a_02_11b --spline-list ./spline_lists/t2k.list 
#
# Author:
#   Costas Andreopoulos <c.andreopoulos \st cern.ch>
#   University of Liverpool
#
# Copyright:
#   Copyright (c) 2003-2025, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#----------------------------------------------------------------------------------------------------------------

use File::Path;
use File::Basename;

# inputs
#  
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--gen-version')    { $gen_version   = $ARGV[$iarg+1]; }
  if($_ eq '--tune')           { $tune          = $ARGV[$iarg+1]; }
  if($_ eq '--spline-list')    { $spline_list   = $ARGV[$iarg+1]; }
  if($_ eq '--input-splines')  { $input_splines = $ARGV[$iarg+1]; }
  if($_ eq '--arch')           { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')     { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')          { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')   { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')   { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')          { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--time-limit')     { $time_limit    = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')   { $softw_topdir  = $ARGV[$iarg+1]; }
  if($_ eq '--jobs-topdir')    { $jobs_topdir   = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined spline list. Use the --spline-list option. See examples in ?]")
unless defined $spline_list;
die("** Aborting [Undefined GENIE Generator version. Use the --gen-version option]")
unless defined $gen_version;
die("** Aborting [Undefined GENIE physics tune. Use the --tune option]")
unless defined $tune;

$use_valgrind      = 0                                       unless defined $use_valgrind;
$arch              = "el7.x86_64"                            unless defined $arch;
$production        = "routine_validation"                    unless defined $production;
$cycle             = "01"                                    unless defined $cycle;
$batch_system      = "Slurm"                                 unless defined $batch_system;
$queue             = "compute"                               unless defined $queue;
$time_limit        = "10:00:00"                              unless defined $time_limit;
$softw_topdir      = "/user/costasa/projects/GENIE/softw/"   unless defined $softw_topdir;
$jobs_topdir       = "/scratch/costasa/GENIE/"               unless defined $jobs_topdir;
$input_splines     = ""                                      unless defined $input_splines;
$gen_setup_script  = "$softw_topdir/generator/builds/$arch/$gen_version-setup.sh";
$splines_name      = basename($spline_list,".list");
$jobs_dir          = "$jobs_topdir/$gen_version-$tune-$production\_$cycle-xsec\_vA\_$splines_name";

$nknots    =  0;
$emax      =  0;
@neutrinos = ();
@targets   = ();

open (REQUIRED_SPLINES, $spline_list);
$i=0;
while(<REQUIRED_SPLINES>)
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
close(REQUIRED_SPLINES);

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
    $filename_basepath = "$jobs_dir/$jobname";

    $grep_pipe  = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
    $gmkspl_opt = "-p $nu_codes -t $_ -n $nknots -e $emax --tune $tune";
    $gmkspl_opt = "gmkspl_opt --output-cross-sections gxspl_$tgt_code.xml";
    $gmkspl_opt = "gmkspl_opt --input-cross-sections $input_splines" if -e $input_splines;
   #$gmkspl_cmd = "gmkspl $gmkspl_opt | $grep_pipe &> $filename_basepath.mkspl.log";
    $gmkspl_cmd = "gmkspl $gmkspl_opt &> $filename_basepath.mkspl.log";

    print "@@ exec: $gmkspl_cmd \n";

    # Submit jobs
    # -------------------------------------------------------------             

    # Slurm case
    if($batch_system eq 'Slurm') {
        ##my $time_lim = `sinfo -h -p batch -o %l`;
        ##my ($days, $hours, $remainder) = $time_lim =~ /([0]+)-([0-9]+):(.*)/;
        ##my $newhours = $days * 24 + $hours;
        ##my $new_time_lim = "$newhours:$remainder";
        ##$time_limit = $new_time_lim lt $time_limit ? $new_time_lim : $time_limit;
	$batch_script = "$filename_basepath.sh";
	open(SLURM, ">$batch_script") or die("Can not create the SLURM batch script");
	print SLURM "#!/bin/bash \n";
        print SLURM "#SBATCH-p $queue \n";
        print SLURM "#SBATCH-J $jobname \n";
        print SLURM "#SBATCH-N 1 \n";   
        print SLURM "#SBATCH-c 1 \n";   
        print SLURM "#SBATCH-o $filename_basepath.slurmout.log \n";
        print SLURM "#SBATCH-e $filename_basepath.slurmerr.log \n";
        print SLURM "#SBATCH-t $time_limit \n";
	print SLURM "source $gen_setup_script \n";
	print SLURM "cd $jobs_dir \n";
	print SLURM "$gmkspl_cmd \n";
        close(SLURM);
	`sbatch --job-name=$jobname $batch_script`;
    } #slurm

    # PBS case
    if($batch_system eq 'PBS') {
	$batch_script = "$filename_basepath.pbs";
	open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
	print PBS "#!/bin/bash \n";
        print PBS "#PBS -N $jobname \n";
        print PBS "#PBS -o $filename_basepath.pbsout.log \n";
        print PBS "#PBS -e $filename_basepath.pbserr.log \n";
	print PBS "source $gen_setup_script \n";
	print PBS "cd $jobs_dir \n";
	print PBS "$gmkspl_cmd \n";
        close(PBS);
        $job_submission_command = "qsub";
        `$job_submission_command -q $queue $batch_script`;
    } #PBS

    # LSF case
    if($batch_system eq 'LSF') {
	$batch_script = "$filename_basepath.sh";
	open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
	print LSF "#!/bin/bash \n";
        print LSF "#BSUB-j $jobname \n";
        print LSF "#BSUB-o $filename_basepath.lsfout.log \n";
        print LSF "#BSUB-e $filename_basepath.lsferr.log \n";
	print LSF "source $gen_setup_script \n";
	print LSF "cd $jobs_dir \n";
	print LSF "$gmkspl_cmd \n";
        close(LSF);
	`bsub < $batch_script`;
    } #LSF
    
    # run interactively
    if($batch_system eq 'none') {
        system("source $gen_setup_script; cd $jobs_dir; $gmkspl_cmd");
    }
}

