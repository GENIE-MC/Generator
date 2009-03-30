#--------------------------------------------------------------------------
# Submit a GENIE/SK event generation job using a JNUBEAM ntuple-based 
# neutrino flux description.
#
# For use at the RAL/PPD Tier2 PBS batch farm.
#
# Syntax:
#   shell$ perl submit-evg_sk_fntp.pl <options>
#
# Options:
#    --fluxrun       : Input flux run number
#    --version       : GENIE version 
#   [--production]   :
#   [--cycle]        :
#   [--use-valgrind] :
#
#
# Example:
#  shell$ perl submit-evg_sk_fntp.pl --fluxrun 180 --version v2.4.0
#                                    --production mdc0 --cycle 01
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#--------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--fluxrun')       { $fluxrun = $ARGV[$iarg+1]; }
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined run number. Use the --fluxrun option]")
unless defined $fluxrun;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$use_valgrind = 0         unless defined $use_valgrind;
$production   = "sktest"  unless defined $production;
$cycle        = "01"      unless defined $cycle;

$nevents        = "50000";   
$mcrun_base     = 10000000;
$mcseed_base    = 210921029;
$batch_queue    = "prod";
$time_limit     = "25:00:00";
$genie_inst_dir = "/opt/ppd/t2k/GENIE";
$production_dir = "/opt/ppd/t2k/GENIE/scratch";
$inputs_dir     = "/opt/ppd/t2k/GENIE/data/job_inputs";
$genie_setup    = "$genie_inst_dir/$genie_version-setup";
$geom_tgt_mix   = "1000080160[0.8879],1000010010[0.1121]";
$xspl_file      = "$inputs_dir/xspl/gxspl-t2k-$genie_version.xml";
$flux_dir       = "$inputs_dir/t2k_flux/07a/sk";
$flux_file_prfx = "nu.2.0deg.sk.";
$flux_file_sufx = ".root";
$flux_file      = "$flux_dir/$flux_file_prfx$fluxrun$flux_file_sufx";
$flux_det_loc   = "sk";
$job_dir        = "$production_dir/sk-$production\_$cycle";
$file_prefix    = "genie_sk";
$mcrun          = $mcrun_base  + $fluxrun;
$mcseed         = $mcseed_base + $fluxrun;

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
`qsub -q $batch_queue $batch_script`;
