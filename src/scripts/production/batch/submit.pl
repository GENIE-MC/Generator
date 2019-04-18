#-------------------------------------------------------------------------------------------
# A script to submit a single job so as to avoid the manual creation of batch script files.
#
# Syntax:
#    % perl submit.pl <options>
#
# Options:
#    --cmd              : Command(s) to execute
#    --version          : GENIE version number
#   [--job-name]        : default: tmp-(random_number)
#   [--arch]            : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]      : default: routine_validation
#   [--cycle]           : default: 01
#   [--batch-system]    : <PBS, LSF, none>, default: PBS
#   [--queue]           : default: prod
#   [--softw-topdir]    : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]     : top level dir for job files, default: /opt/ppd/t2k/softw/scratch/GENIE/
#
# Example:
#   % perl submit.pl 
#     --cmd 'gevgen -p 14 -t 1000260560 -e 1 --cross-sections /some/path/xsec.xml \
#     --seed 1989298 --r 100; mv gntp.100.ghep.root /some/other/path/' \
#     --version v2.8.0 \
#     --job-name numuFe100
#
# Author:
#   Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#-------------------------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

$iarg=0;
foreach (@ARGV) {
  if($_ eq '--cmd')           { $cmd             = $ARGV[$iarg+1]; }
  if($_ eq '--version')       { $genie_version   = $ARGV[$iarg+1]; }
  if($_ eq '--job-name')      { $name            = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch            = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production      = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle           = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')  { $batch_system    = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue           = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir    = $ARGV[$iarg+1]; }
  if($_ eq '--jobs-topdir')   { $jobs_topdir     = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined GENIE command. Use the --cmd option]")
unless defined $cmd;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$random_num = int (rand(999999999));

$name           = "tmp-$random_num"             unless defined $name;
$arch           = "SL6.x86_64"                  unless defined $arch;
$production     = "routine_validation"          unless defined $production;
$cycle          = "01"                          unless defined $cycle;
$batch_system   = "PBS"                         unless defined $batch_system;
$queue          = "prod"                        unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE"    unless defined $softw_topdir;
$jobs_topdir    = "/opt/ppd/t2k/scratch/GENIE/" unless defined $jobs_topdir;
$genie_setup    = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$job_dir        = "$jobs_topdir/$genie_version-$production\_$cycle-$name";

# make the job directory
#
mkpath ($job_dir, {verbose => 1, mode=>0777});


# submit job
#

# PBS case
if($batch_system eq 'PBS') {
        $batch_script = "$job_dir/$name.pbs";
        open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
        print PBS "#!/bin/bash \n";
        print PBS "#PBS -N $name \n";
        print PBS "#PBS -o $job_dir/$name.pbsout.log \n";
        print PBS "#PBS -e $job_dir/$name.pbserr.log \n";
        print PBS "source $genie_setup \n";
        print PBS "cd $job_dir \n";
        print PBS "$cmd \n";
        close(PBS);
        `qsub -q $queue $batch_script`;
} #PBS

# LSF case
if($batch_system eq 'LSF') {
        $batch_script = "$job_dir/$name.sh";
        open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
        print LSF "#!/bin/bash \n";
        print PBS "#BSUB-j $name \n";
        print LSF "#BSUB-q $queue \n";
        print LSF "#BSUB-o $job_dir/$name.lsfout.log \n";
        print LSF "#BSUB-e $job_dir/$name.lsferr.log \n";
        print LSF "source $genie_setup \n";
        print LSF "cd $job_dir \n";
        print LSF "$cmd \n";
        `bsub < $batch_script`;
} #LSF


