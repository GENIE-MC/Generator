#----------------------------------------------------------------------------------------------------
# Submit a GENIE/SK event generation job using a histogram-based neutrino flux description
#
# Syntax:
#   shell% perl submit_t2k_superk_fhst_mc_jobs.pl <options>
#
# Options:
#    --version         : GENIE version
#    --run             : 1,2,3,...
#    --neutrino        : numu, numubar, nue, nuebar, nuesig
#   [--flux-version]   : JNUBEAM flux version, <07a, 10a, 10c, 11a...>, default: 11a
#   [--flux-config]    : JNUBEAM config, <nominal, yshift2mm,...>, default: nominal
#   [--flux-hist-file] : JNUBEAM flux histogram file, default: sk_flux_histograms.root
#   [--arch]           : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]     : default: <version>
#   [--cycle]          : default: 01
#   [--use-valgrind]   : default: off
#   [--batch-system]   : <PBS, LSF>, default: PBS
#   [--queue]          : default: prod
#   [--softw-topdir]   : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE
#   [--jobs-topdir]    : top level dir for job files, default: /opt/ppd/t2k/scratch/GENIE/
#
# Example:
#   shell& perl submit_t2k_superk_fhst_mc_jobs.pl --run 180 --neutrino numubar --version v2.5.1
#
# Author:
#   Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#----------------------------------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')        { $genie_version  = $ARGV[$iarg+1]; }
  if($_ eq '--run'     )       { $run            = $ARGV[$iarg+1]; }
  if($_ eq '--neutrino')       { $neutrino       = $ARGV[$iarg+1]; }
  if($_ eq '--flux-version')   { $flux_version   = $ARGV[$iarg+1]; }
  if($_ eq '--flux-config')    { $flux_config    = $ARGV[$iarg+1]; }
  if($_ eq '--flux-hist-file') { $flux_hist_file = $ARGV[$iarg+1]; }
  if($_ eq '--arch')           { $arch           = $ARGV[$iarg+1]; }
  if($_ eq '--production')     { $production     = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')          { $cycle          = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')   { $use_valgrind   = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')   { $batch_system   = $ARGV[$iarg+1]; }
  if($_ eq '--queue')          { $queue          = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')   { $softw_topdir   = $ARGV[$iarg+1]; }  
  if($_ eq '--jobs-topdir')    { $jobs_topdir   = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined run number. Use the --run option]")
unless defined $run;
die("** Aborting [Undefined neutrino type. Use the --neutrino option]")
unless defined $neutrino;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version; 

$use_valgrind   = 0                             unless defined $use_valgrind;
$arch           = "SL6.x86_64"                  unless defined $arch;
$production     = "$genie_version"              unless defined $production;
$cycle          = "01"                          unless defined $cycle;
$batch_system   = "PBS"                         unless defined $batch_system;
$queue          = "prod"                        unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE"    unless defined $softw_topdir;
$flux_version   = "11a"                         unless defined $flux_version;
$flux_config    = "nominal"                     unless defined $flux_config;
$flux_hist_file = "sk_flux_histograms.root"     unless defined $flux_hist_file;
$jobs_topdir    = "/opt/ppd/t2k/scratch/GENIE/" unless defined $jobs_topdir;
$nevents        = "2000";   
$time_limit     = "05:00:00";
$inputs_dir     = "$softw_topdir/data/job_inputs";
$genie_setup    = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$geom_tgt_mix   = "1000080160[0.8879],1000010010[0.1121]";
$xspl_file      = "$inputs_dir/xspl/gxspl-t2k-$genie_version.xml";
$flux_file      = "$inputs_dir/t2k_flux/$flux_version/sk/$flux_config/$flux_hist_file";
$job_dir        = "$jobs_topdir/skmc-$production\_$cycle-$neutrino";
$file_prefix    = "genie_sk";

%mcseed_base    = ( 'numu'    => '183221029',
                    'numubar' => '283221029',
                    'nue'     => '383221029',
                    'nuebar'  => '383221029',
                    'nuesig'  => '483221029' );
%mcrun_base     = ( 'numu'    => '10000000',
                    'numubar' => '10000000',
                    'nue'     => '10000000',
                    'nuebar'  => '10000000',
                    'nuesig'  => '10000000' );
%nu_pdg_code    = ( 'numu'    =>  '14',
                    'numubar' => '-14',
                    'nue'     =>  '12',
                    'nuebar'  => '-12',
                    'nuesig'  =>  '12' );
%flux_hist_name = ( 'numu'    => 'numu_flux',
                    'numubar' => 'numubar_flux',
                    'nue'     => 'nue_flux',
                    'nuebar'  => 'nuebar_flux',
                    'nuesig'  => 'numu_flux' );

die("** Aborting [Can not find GENIE setup script: ..... $genie_setup]") 
unless -e $genie_setup;
die("** Aborting [Can not find flux file: .............. $flux_file]")   
unless -e $flux_file;
die("** Aborting [Can not find xsec file: .............. $xspl_file]")   
unless -e $xspl_file;

# make the jobs directory
# 
mkpath ($job_dir, {verbose => 1, mode=>0777}); 

# form mcrun, mcseed numbers 
#
$mcrun  = $run + $mcrun_base  {$neutrino}; 
$mcseed = $run + $mcseed_base {$neutrino};

print "@@@ Will submit job with MC run number = $mcrun (seed number = $mcseed)\n";

# get neutrino code and flux histogram name
#
$nu  = $nu_pdg_code    {$neutrino};
$hst = $flux_hist_name {$neutrino};

# form event generation and file conversion commands
#
$fntemplate    = "$job_dir/skjob-$mcrun";
$ghep_file     = "$file_prefix.$production\_$cycle.$neutrino.$mcrun.ghep.root";
#$grep_pipe     = "grep -B 50 -A 50 -i \"warn\\|error\\|fatal\"";
$grep_pipe     = "grep -B 50 -A 50 -i fatal";
$valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
$evgen_opt     = "-g $geom_tgt_mix -f $flux_file,$nu\[$hst\] -r $mcrun --seed $mcseed -n $nevents --cross-sections $xspl_file";
$evgen_cmd     = "gevgen_t2k $evgen_opt | $grep_pipe &> $fntemplate.evgen.log";
$frenm_cmd     = "mv gntp.$mcrun.ghep.root $ghep_file";
$fconv_cmd     = "gntpc -f t2k_tracker -i $ghep_file --seed $mcseed";

print "@@@ exec: $evgen_cmd \n";

#
# submit
#

# PBS case
if($batch_system eq 'PBS') {
  $batch_script = "$fntemplate.pbs";
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
  $batch_script = "$fntemplate.sh";
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
