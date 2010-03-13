#--------------------------------------------------------------------------------------------------
# Submit a GENIE/SK event generation job using a JNUBEAM ntuple-based neutrino flux description.
#
# For use at the RAL/PPD Tier2 PBS batch farm.
#
# Syntax:
#   shell$ perl submit-evg_sk_fntp.pl <options>
#
# Options:
#    --version             : GENIE version 
#    --flux-run            : Input flux run number
#   [--flux-version]       : JNUBEAM flux version, <07a, 10>, default: 10
#   [--flux-file-prefix]   : JNUBEAM flux file prefix, default: nu.sk_horn250ka.
#   [--flux-file-suffix]   : JNUBEAM flux file suffix, default: .root
#   [--arch]               : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]         : default: <version>
#   [--cycle]              : default: 01
#   [--use-valgrind]       : default: off
#   [--queue]              : default: prod
#   [--softw-topdir]       : default: /opt/ppd/t2k/GENIE
#
# Example:
#   shell$ perl submit-evg_sk_fntp.pl --flux-run 180 --version v2.4.0 --production mdc0 --cycle 01
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
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
  if($_ eq '--flux-file-prefix')  { $flux_file_prefix   = $ARGV[$iarg+1]; }
  if($_ eq '--flux-file-suffix')  { $flux_file_suffix   = $ARGV[$iarg+1]; }
  if($_ eq '--arch')              { $arch               = $ARGV[$iarg+1]; }
  if($_ eq '--production')        { $production         = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')             { $cycle              = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')      { $use_valgrind       = $ARGV[$iarg+1]; }
  if($_ eq '--queue')             { $queue              = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')      { $softw_topdir       = $ARGV[$iarg+1]; }   
  $iarg++;
}
die("** Aborting [Undefined run number. Use the --flux-run option]")
unless defined $flux_run;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$use_valgrind      = 0                          unless defined $use_valgrind;
$arch              = "SL5_64bit"                unless defined $arch;
$production        = "$genie_version"           unless defined $production;
$cycle             = "01"                       unless defined $cycle;
$queue             = "prod"                     unless defined $queue;
$softw_topdir      = "/opt/ppd/t2k/GENIE"       unless defined $softw_topdir;
$flux_version      = "10"                       unless defined $flux_version;
$flux_file_prefix  = "nu.sk_horn250ka."         unless defined $flux_file_prefix;
$flux_file_suffix  = ".root"                    unless defined $flux_file_suffix;
$nevents           = "50000";   
$mcrun_base        = 10000000;
$mcseed_base       = 210921029;
$time_limit        = "25:00:00";
$production_dir    = "$softw_topdir/scratch";
$inputs_dir        = "$softw_topdir/data/job_inputs";
$genie_setup       = "$softw_topdir/builds/$arch/$genie_version-setup";
$geom_tgt_mix      = "1000080160[0.8879],1000010010[0.1121]";
$xspl_file         = "$inputs_dir/xspl/gxspl-t2k-$genie_version.xml";
$flux_dir          = "$inputs_dir/t2k_flux/$fluxversion/sk";
$flux_file         = "$flux_dir/$flux_file_prfx$flux_run$flux_file_sufx";
$flux_det_loc      = "sk";
$job_dir           = "$production_dir/skmc-$production\_$cycle";
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

$batch_script  = "$job_dir/skjob-$mcrun.pbs";
$logfile_evgen = "$job_dir/skjob-$mcrun.evgen.log";   
$logfile_conv  = "$job_dir/skjob-$mcrun.conv.log";  
$logfile_pbse  = "$job_dir/skjob-$mcrun.pbs_e.log";
$logfile_pbso  = "$job_dir/skjob-$mcrun.pbs_o.log";
$ghep_file     = "$file_prefix.$production\_$cycle.$mcrun.ghep.root";
$grep_filt     = "grep -B 50 -A 50 -i \"warn\\|error\\|fatal\"";
$valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
$evgen_cmd     = "gT2Kevgen -g $geom_tgt_mix -f $flux_file,$flux_det_loc -r $mcrun -n $nevents | $grep_filt &> $logfile_evgen";
$frenm_cmd     = "mv gntp.$mcrun.ghep.root $ghep_file";
$fconv_cmd     = "gntpc -f t2k_tracker -i $ghep_file &> $logfile_conv";

# create the PBS script
#
open(PBS, ">$batch_script") or die("Can not create the PBS batch script");

print PBS "#!/bin/bash \n";
print PBS "#PBS -l cput=$time_limit \n";
print PBS "#PBS -N $mcrun\_sk-$production-$cycle \n";
print PBS "#PBS -o $logfile_pbso \n";
print PBS "#PBS -e $logfile_pbse \n";
print PBS "source $genie_setup \n";
print PBS "cd $job_dir \n";
print PBS "export GSPLOAD=$xspl_file \n";
print PBS "unset GEVGL \n";
print PBS "export GSEED=$mcseed \n";
print PBS "$evgen_cmd \n";
print PBS "$frenm_cmd \n";
print PBS "$fconv_cmd \n";

print "@@@ exec: $evgen_cmd \n";

# submit job
#
`qsub -q $queue $batch_script`;
