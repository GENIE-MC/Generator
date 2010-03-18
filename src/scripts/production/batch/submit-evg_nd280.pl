#-------------------------------------------------------------------------------------------------
# Submit GENIE/nd280 event generation job using the full geometry and a JNUBEAM ntuple-based flux 
# description.
#
# Syntax:
#   shell$ perl submit-evg_nd280.pl <options>
#
# Options:
#    --version             : GENIE version
#    --flux-run            : Input flux run number
#   [--flux-version]       : JNUBEAM flux version, <07a, 10>, default: 10
#   [--flux-file-prefix]   : JNUBEAM flux file prefix, default: nu.nd.
#   [--flux-file-suffix]   : JNUBEAM flux file suffix, default: .root
#   [--arch]               : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]         : default: <version>
#   [--cycle]              : default: 01
#   [--use-valgrind]       : default: off
#   [--batch-system]       : <PBS, >, default: PBS
#   [--queue]              : default: prod
#   [--softw-topdir]       : default: /opt/ppd/t2k/GENIE
#
# Example:
#   shell$ perl submit-evg_nd280.pl --flux-run 180 --version v2.4.0 --production mdc0 --cycle 01
#
# Tested at the RAL/PPD Tier2 PBS batch farm.
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#-------------------------------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')           { $genie_version    = $ARGV[$iarg+1]; }
  if($_ eq '--flux-run')          { $flux_run         = $ARGV[$iarg+1]; }
  if($_ eq '--flux-version')      { $flux_version     = $ARGV[$iarg+1]; }
  if($_ eq '--flux-file-prefix')  { $flux_file_prefix = $ARGV[$iarg+1]; }
  if($_ eq '--flux-file-suffix')  { $flux_file_suffix = $ARGV[$iarg+1]; }
  if($_ eq '--arch')              { $arch             = $ARGV[$iarg+1]; }
  if($_ eq '--production')        { $production       = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')             { $cycle            = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')      { $use_valgrind     = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')      { $batch_system     = $ARGV[$iarg+1]; }
  if($_ eq '--queue')             { $queue            = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')      { $softw_topdir     = $ARGV[$iarg+1]; } 
  $iarg++;
}
die("** Aborting [Undefined flux file run #. Use the --flux-run option]")
unless defined $flux_run;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$use_valgrind      = 0                         unless defined $use_valgrind;
$arch              = "SL5_64bit"               unless defined $arch;
$production        = "$genie_version"          unless defined $production;
$cycle             = "01"                      unless defined $cycle;
$batch_system      = "PBS"                     unless defined $batch_system;
$queue             = "prod"                    unless defined $queue;
$softw_topdir      = "/opt/ppd/t2k/GENIE"      unless defined $softw_topdir;
$flux_version      = "10"                      unless defined $flux_version;
$flux_file_prefix  = "nu.nd."                  unless defined $flux_file_prefix;
$flux_file_suffix  = ".root"                   unless defined $flux_file_suffix;
$job_pot           = "1E+18";
$mcrun_base        = 10000000;
$mcseed_base       = 210921029;
$time_limit        = "30:00:00";
$genie_setup       = "$softw_topdir/builds/$arch/$genie_version-setup";
$production_dir    = "$softw_topdir/scratch";
$inputs_dir        = "$softw_topdir/data/job_inputs";
$job_dir           = "$production_dir/nd280mc-$production\_$cycle";
$flux_dir          = "$inputs_dir/t2k_flux/$flux_version/nd";
$flux_file         = "$flux_dir/$flux_file_prefix$flux_run$flux_file_suffix";
$flux_det_loc      = "nd5";
$geom_file         = "$inputs_dir/t2k_geom/ND280.root";
$xspl_file         = "$inputs_dir/xspl/gxspl-t2k-$genie_version.xml";
$geom_dunits       = "clhep_def_density_unit";
$geom_lunits       = "mm";
$file_prefix       = "genie_nd280";
$mcrun             = $mcrun_base  + $flux_run; 
$mcseed            = $mcseed_base + $flux_run;

die("** Aborting [Can not find GENIE setup script: ....... $genie_setup]") 
unless -e $genie_setup;
die("** Aborting [Can not find flux file: ................ $flux_file]") 
unless -e $flux_file;
die("** Aborting [Can not find geometry file: ............ $geom_file]") 
unless -e $geom_file;
die("** Aborting [Can not find cross sections file: ...... $xspl_file]") 
unless -e $xspl_file;

print "@@@ Will submit job with MC run number = $mcrun (seed number = $mcseed)\n";

# make the jobs directory
#
mkpath ($job_dir, {verbose => 1, mode=>0777}); 
die("$job_dir doesn't exist") unless -d $job_dir;

# form event generation and file conversion commands
#
$job_file_base = "$job_dir/nd280job-$mcrun";
$ghep_file     = "$file_prefix.$production\_$cycle.$mcrun.ghep.root";
$grep_pipe     = "grep -B 50 -A 30 -i \"warn\\|error\\|fatal\" ";
$valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
$evgen_cmd     = "gT2Kevgen -g $geom_file -f $flux_file,$flux_det_loc -r $mcrun -L $geom_lunits -D $geom_dunits -E $job_pot | $grep_pipe &> $job_file_base.evgen.log";
$frenm_cmd     = "mv gntp.$mcrun.ghep.root $ghep_file";
$fconv_cmd     = "gntpc -f t2k_rootracker -i $ghep_file | $grep_pipe &> $job_file_base.conv.log";

print "@@@ exec: $evgen_cmd \n";

#
# submit
#

# PBS case
if($batch_system eq 'PBS') {
  $batch_script  = "$job_file_base.pbs";
  open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
  print PBS "#!/bin/bash \n";
  print PBS "#PBS -N $mcrun\_nd280-$production-$cycle \n";
  print PBS "#PBS -l cput=$time_limit \n";
  print PBS "#PBS -o $job_file_base.pbsout.log \n";
  print PBS "#PBS -e $job_file_base.pbserr.log \n";
  print PBS "source $genie_setup \n";
  print PBS "cd $job_dir \n";
  print PBS "export GSPLOAD=$xspl_file \n";
  print PBS "unset GEVGL \n";
  print PBS "export GSEED=$mcseed \n";
  print PBS "$evgen_cmd \n";
  print PBS "$frenm_cmd \n";
  print PBS "$fconv_cmd \n";
  close(PBS);
  `qsub -q $queue $batch_script`;
}

