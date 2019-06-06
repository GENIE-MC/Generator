#--------------------------------------------------------------------------------------------------
# Submit a GENIE/SK event generation job using a JNUBEAM ntuple-based neutrino flux description.
#
# Syntax:
#   shell$ perl submit_t2k_superk_fntp_mc_jobs.pl <options>
#
# Options:
#    --version           : GENIE version 
#    --flux-run          : Input flux run number
#   [--flux-version]     : JNUBEAM flux version, <07a, 10a, 10b, 10c, ...>, default: 10c
#   [--flux-config]      : JNUBEAM config, <nominal, yshift2mm,...>, default: nominal
#   [--flux-file-prefix] : JNUBEAM flux file prefix, default: nu.sk_horn250ka.
#   [--flux-file-suffix] : JNUBEAM flux file suffix, default: .root
#   [--arch]             : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]       : default: <version>
#   [--cycle]            : default: 01
#   [--use-valgrind]     : default: off
#   [--batch-system]     : <PBS, LSF>, default: PBS
#   [--queue]            : default: prod
#   [--softw-topdir]     : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE
#   [--jobs-topdir]      : top level dir for job files, default: /opt/ppd/t2k/scratch/GENIE/
#
# Example:
#   shell$ perl submit_t2k_superk_fntp_mc_jobs.pl --flux-run 180 \
#                          --version v2.4.0 --production mdc0 --cycle 01
# Author:
#   Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#--------------------------------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')           { $genie_version      = $ARGV[$iarg+1]; }
  if($_ eq '--flux-run')          { $flux_run           = $ARGV[$iarg+1]; }
  if($_ eq '--flux-version')      { $flux_version       = $ARGV[$iarg+1]; }
  if($_ eq '--flux-config')       { $flux_config        = $ARGV[$iarg+1]; }
  if($_ eq '--flux-file-prefix')  { $flux_file_prefix   = $ARGV[$iarg+1]; }
  if($_ eq '--flux-file-suffix')  { $flux_file_suffix   = $ARGV[$iarg+1]; }
  if($_ eq '--arch')              { $arch               = $ARGV[$iarg+1]; }
  if($_ eq '--production')        { $production         = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')             { $cycle              = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')      { $use_valgrind       = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')      { $batch_system       = $ARGV[$iarg+1]; }
  if($_ eq '--queue')             { $queue              = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')      { $softw_topdir       = $ARGV[$iarg+1]; }   
  if($_ eq '--jobs-topdir')       { $jobs_topdir        = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined run number. Use the --flux-run option]")
unless defined $flux_run;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$use_valgrind      = 0                             unless defined $use_valgrind;
$arch              = "SL6.x86_64"                  unless defined $arch;
$production        = "$genie_version"              unless defined $production;
$cycle             = "01"                          unless defined $cycle;
$batch_system      = "PBS"                         unless defined $batch_system;
$queue             = "prod"                        unless defined $queue;
$softw_topdir      = "/opt/ppd/t2k/softw/GENIE"    unless defined $softw_topdir;
$flux_version      = "10c"                         unless defined $flux_version;
$flux_config       = "nominal"                     unless defined $flux_config;
$flux_file_prefix  = "nu.sk_horn250ka."            unless defined $flux_file_prefix;
$flux_file_suffix  = ".root"                       unless defined $flux_file_suffix;
$jobs_topdir       = "/opt/ppd/t2k/scratch/GENIE/" unless defined $jobs_topdir;
$nevents           = "50000";   
$mcrun_base        = 10000000;
$mcseed_base       = 210921029;
$time_limit        = "25:00:00";
$inputs_dir        = "$softw_topdir/data/job_inputs";
$genie_setup       = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$geom_tgt_mix      = "1000080160[0.8879],1000010010[0.1121]";
$xspl_file         = "$inputs_dir/xspl/gxspl-t2k-$genie_version.xml";
$flux_dir          = "$inputs_dir/t2k_flux/$flux_version/sk/$flux_config";
$flux_file         = "$flux_dir/$flux_file_prefix$flux_run$flux_file_suffix";
$flux_det_loc      = "sk";
$job_dir           = "$jobs_topdir/skmc-$production\_$cycle";
$file_prefix       = "genie_sk";
$mcrun             = $mcrun_base  + $flux_run;
$mcseed            = $mcseed_base + $flux_run;

die("** Aborting [Can not find GENIE setup script: ..... $genie_setup]") 
unless -e $genie_setup;
die("** Aborting [Can not find flux file: .............. $flux_file]")   
unless -e $flux_file;
die("** Aborting [Can not find xsec file: .............. $xspl_file]")   
unless -e $xspl_file;

# make the jobs directory
# 
mkpath ($job_dir, {verbose => 1, mode=>0777}); 

print "@@@ Will submit job with MC run number = $mcrun (seed number = $mcseed)\n";

# form event generation and file conversion commands
#
$fntemplate    = "$job_dir/skjob-$mcrun";   
$ghep_file     = "$file_prefix.$production\_$cycle.$mcrun.ghep.root";
#$grep_pipe     = "grep -B 50 -A 50 -i \"warn\\|error\\|fatal\"";
$grep_pipe     = "grep -B 50 -A 50 -i fatal";
$valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
$evgen_opt     = "-g $geom_tgt_mix -f $flux_file,$flux_det_loc -r $mcrun --seed $mcseed -n $nevents --cross-sections $xspl_file";
$evgen_cmd     = "gevgen_t2k $evgen_opt | $grep_pipe &> $fntemplate.evgen.log";
$frenm_cmd     = "mv gntp.$mcrun.ghep.root $ghep_file";
$fconv_cmd     = "gntpc -f t2k_tracker -i $ghep_file --seed $mcseed &> $fntemplate.conv.log";

print "@@@ exec: $evgen_cmd \n";

#
# submit
#
  
# PBS case
if($batch_system eq 'PBS') {
  $batch_script  = "$fntemplate.pbs";
  open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
  print PBS "#!/bin/bash \n";
  print PBS "#PBS -N $mcrun\_sk-$production-$cycle \n";
  print PBS "#PBS -l cput=$time_limit \n";
  print PBS "#PBS -o $fntemplate.pbsout.log \n";
  print PBS "#PBS -e $fntemplate.pbserr.log \n";
  print PBS "source $genie_setup \n";
  print PBS "cd $job_dir \n";
  print PBS "$evgen_cmd \n";
  print PBS "$frenm_cmd \n";
  print PBS "$fconv_cmd \n";
  close(PBS);
  `qsub -q $queue $batch_script`;
}

# LSF case
if($batch_system eq 'LSF') {
  $batch_script  = "$fntemplate.sh";
  open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
  print LSF "#!/bin/bash \n";
  print LSF "#BSUB-j $mcrun\_sk-$production-$cycle \n";
  print LSF "#BSUB-q $queue \n";
  print LSF "#BSUB-c $time_limit \n";
  print LSF "#BSUB-o $fntemplate.lsfout.log \n";
  print LSF "#BSUB-e $fntemplate.lsferr.log \n";
  print LSF "source $genie_setup \n";
  print LSF "cd $job_dir \n";
  print LSF "$evgen_cmd \n";
  print LSF "$frenm_cmd \n";
  print LSF "$fconv_cmd \n";
  close(LSF);
  `bsub < $batch_script`;
}
