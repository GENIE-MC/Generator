#--------------------------------------------------------------------------------------------------
# A helper script to submit a single job with all commands specified directly on the command line,
# avoiding the need to manually create a separate batch script for one-off or ad hoc tasks.
#
# Syntax:
#    % perl submit_inline.pl <options>
#
# Options:
#    --cmd               : Command(s) to execute
#    --gen-version       : GENIE generator version number
#   [--job-name]         : default: tmp-(random_number)
#   [--arch]             : <el7.x86_64, ...>, default: el7.x86_64
#   [--production]       : default: prod
#   [--cycle]            : default: 01
#   [--batch-system]     : <Slurm, PBS, LSF>, default: Slurm
#   [--queue]            : default: compute
#   [--time-limit]       : default: 10:00:00
#   [--softw-topdir]     : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--job-topdir]       : top level dir for job files, default: /opt/ppd/t2k/softw/scratch/GENIE/
#
# Example:
#   % perl submit_inline.pl 
#     --cmd 'gevgen -p 14 -t 1000260560 -e 1 --cross-sections /some/path/xsec.xml \
#     --seed 1989298 --r 100; mv gntp.100.ghep.root /some/other/path/' \
#     --gen-version v3.6.0 \
#     --job-name numuFe100
#
# Author:
#   Costas Andreopoulos <c.andreopoulos \st cern.ch>
#   University of Liverpool
#
# Copyright:
#   Copyright (c) 2003-2025, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#--------------------------------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

$iarg=0;
foreach (@ARGV) {
  if($_ eq '--cmd')           { $cmd             = $ARGV[$iarg+1]; }
  if($_ eq '--gen-version')   { $gen_version     = $ARGV[$iarg+1]; }
  if($_ eq '--job-name')      { $job_name        = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch            = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production      = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle           = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')  { $batch_system    = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue           = $ARGV[$iarg+1]; }
  if($_ eq '--time-limit')    { $time_limit      = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir    = $ARGV[$iarg+1]; }
  if($_ eq '--job-topdir')    { $job_topdir      = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined GENIE command. Use the --cmd option]")
unless defined $cmd;
die("** Aborting [Undefined GENIE version. Use the --gen-version option]")
unless defined $gen_version;

$random_num = int (rand(999999999));

$arch              = "el7.x86_64"                            unless defined $arch;
$production        = "prod"                                  unless defined $production;
$cycle             = "01"                                    unless defined $cycle;
$batch_system      = "Slurm"                                 unless defined $batch_system;
$queue             = "compute"                               unless defined $queue;
$time_limit        = "10:00:00"                              unless defined $time_limit;
$softw_topdir      = "/user/costasa/projects/GENIE/softw/"   unless defined $softw_topdir;
$job_topdir        = "/scratch/costasa/GENIE/"               unless defined $job_topdir;
$job_name          = "tmp-$random_num"                       unless defined $job_name;
$job_dir           = "$job_topdir/$production\_$cycle-$gen_version-$job_name";
$gen_setup_script  = "$softw_topdir/generator/builds/$arch/$gen_version-setup.sh";
$filename_basepath = "$job_dir/$jobname";

# make the job directory
mkpath ($job_dir, {verbose => 1, mode=>0777});

# submit job
# -----------------------

# Slurm case
if($batch_system eq 'Slurm') {
	$batch_script = "$filename_basepath.sh";
        open(SLURM, ">$batch_script") or die("Can not create the SLURM batch script");
        print SLURM "#!/bin/bash \n";
        print SLURM "#SBATCH-p $queue \n";
        print SLURM "#SBATCH-J $job_name \n";
        print SLURM "#SBATCH-N 1 \n";
        print SLURM "#SBATCH-c 1 \n";
        print SLURM "#SBATCH-o $filename_basepath.slurmout.log \n";
        print SLURM "#SBATCH-e $filename_basepath.slurmerr.log \n";
        print SLURM "#SBATCH-t $time_limit \n";
        print SLURM "source $gen_setup_script \n";
        print SLURM "cd $job_dir \n";
        print SLURM "$gmkspl_cmd \n";
        close(SLURM);
        `sbatch $batch_script`;
} #slurm

# PBS case
if($batch_system eq 'PBS') {
        $batch_script = "$filename_basepath.pbs";
        open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
        print PBS "#!/bin/bash \n";
        print PBS "#PBS -N $job_name \n";
        print PBS "#PBS -o $filename_basepath.pbsout.log \n";
        print PBS "#PBS -e $filename_basepath.pbserr.log \n";
        print PBS "source $gen_setup_script \n";
        print PBS "cd $job_dir \n";
        print PBS "$cmd \n";
        close(PBS);
        `qsub -q $queue $batch_script`;
} #PBS

# LSF case
if($batch_system eq 'LSF') {
        $batch_script = "$filename_basepath.sh";
        open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
        print LSF "#!/bin/bash \n";
        print PBS "#BSUB-j $job_name \n";
        print LSF "#BSUB-q $queue \n";
        print LSF "#BSUB-o $filename_basepath.lsfout.log \n";
        print LSF "#BSUB-e $filename_basepath.lsferr.log \n";
        print LSF "source $gen_setup_script \n";
        print LSF "cd $job_dir \n";
        print LSF "$cmd \n";
        `bsub < $batch_script`;
} #LSF


