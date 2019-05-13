#!/usr/bin/perl

#-----------------------------------------------------------------------------------------------------------
# Script to submit jobs used for evaluating neutrino interaction systematics.
#
# Syntax:
#   perl submit_systematic_param_scan_1d.pl <options>
#
# Options:
#    --version         : GENIE version number.
#    --input-events    : Input event file.
#    --syst            : Systematics for which to submit weight calculation jobs, 
#                        eg `--syst MaCCQE', `--syst MaCCQE,MFP_pi', `--syst all' etc.
#                        Alternatively you can use `predefined' lists, eg  `--syst standard-t2k-analysis-list'
#   [--nevents]        : Event numbers to process. 
#                        To process the first 10000 events, type `--nevents 10000'.
#                        To process the 3000 events between 1000 and 4999, type `--nevents 1000,4999'.
#   [--output-weights] : Output weight file. Default: weights_<name of systematic>.root
#   [--arch]           : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]     : Production name. Default: <version>
#   [--cycle]          : Cycle in current production. Default: 01
#   [--use-valgrind]   : Use Valgrind? Default: off
#   [--batch-system]   : Batch system: <PBS, LSF, slurm, none>. Default: PBS
#   [--queue]          : Batch queue. Default: prod
#   [--softw-topdir]   : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]    : top level dir for job files, default: /opt/ppd/t2k/softw/scratch/GENIE/
#
# Author:
#   Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#-----------------------------------------------------------------------------------------------------------
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')        { $genie_version  = $ARGV[$iarg+1]; }
  if($_ eq '--syst')           { $syst           = $ARGV[$iarg+1]; }
  if($_ eq '--input-events')   { $input_events   = $ARGV[$iarg+1]; }
  if($_ eq '--nevents')        { $nevents        = $ARGV[$iarg+1]; }
  if($_ eq '--output-weights') { $output_weights = $ARGV[$iarg+1]; }
  if($_ eq '--arch')           { $arch           = $ARGV[$iarg+1]; }
  if($_ eq '--production')     { $production     = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')          { $cycle          = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')   { $use_valgrind   = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')   { $batch_system   = $ARGV[$iarg+1]; }
  if($_ eq '--queue')          { $queue          = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')   { $softw_topdir   = $ARGV[$iarg+1]; }  
  if($_ eq '--jobs-topdir')    { $jobs_topdir    = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;
die("** Aborting [You need to specify which systematics to study. Use the --syst option]")
unless defined $syst;
die("** Aborting [You need to specify an input event file. Use the --input-events option]")
unless defined $input-events;

$use_valgrind   = 0                             unless defined $use_valgrind;
$arch           = "SL6.x86_64"                  unless defined $arch;
$production     = "$genie_version"              unless defined $production;
$cycle          = "01"                          unless defined $cycle;
$batch_system   = "PBS"                         unless defined $batch_system;
$queue          = "prod"                        unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE"    unless defined $softw_topdir;
$jobs_topdir    = "/opt/ppd/t2k/scratch/GENIE/" unless defined $jobs_topdir;
$time_limit     = "60:00:00";
$genie_setup    = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$jobs_dir       = "$jobs_topdir/rwght1scan-$production\_$cycle";
$gexec          = "grwght1scan"; 

%def_ntwkdials = ( 
 'MaNCEL'              =>  '3',
 'EtaNCEL'             =>  '3',
 'NormCCQE'            =>  '3',
 'MaCCQE'              => '11',
 'MaCCQEshape'         => '11',
 'VecFFCCQEshape'      => '11',
 'NormCCRES'           =>  '3',
 'MaCCRESshape'        => '11',
 'MvCCRESshape'        => '11',
 'MaCCRES'             => '11',
 'MvCCRES'             => '11',
 'NormNCRES'           =>  '3',
 'MaNCRESshape'        => '11',
 'MvNCRESshape'        => '11',
 'MaNCRES'             => '11',
 'MvNCRES'             => '11',
 'MaCOHpi'             => '11',
 'R0COHpi'             => '11',
 'NonRESBGvpCC1pi'     =>  '3',
 'NonRESBGvpCC2pi'     =>  '3',
 'NonRESBGvpNC1pi'     =>  '3',
 'NonRESBGvpNC2pi'     =>  '3',
 'NonRESBGvnCC1pi'     =>  '3',
 'NonRESBGvnCC2pi'     =>  '3',
 'NonRESBGvnNC1pi'     =>  '3',
 'NonRESBGvnNC2pi'     =>  '3',
 'NonRESBGvbarpCC1pi'  =>  '3',
 'NonRESBGvbarpCC2pi'  =>  '3',
 'NonRESBGvbarpNC1pi'  =>  '3',
 'NonRESBGvbarpNC2pi'  =>  '3',
 'NonRESBGvbarnCC1pi'  =>  '3',
 'NonRESBGvbarnCC2pi'  =>  '3',
 'NonRESBGvbarnNC1pi'  =>  '3',
 'NonRESBGvbarnNC2pi'  =>  '3',
 'AhtBY'               => '11',
 'BhtBY'               => '11',
 'CV1uBY'              => '11',
 'CV2uBY'              => '11',
 'AhtBYshape'          => '11',
 'BhtBYshape'          => '11',
 'CV1uBYshape'         => '11',
 'CV2uBYshape'         => '11',
 'NormDISCC'           => '11',
 'RnubarnuCC'          => '11',
 'DISNuclMod'          => '11',
 'AGKYxF1pi'           => '11',
 'AGKYpT1pi'           => '11',
 'FormZone'            => '11',
 'MFP_pi'              => '11',
 'MFP_N'               => '11',
 'FrCEx_pi'            => '11',
 'FrElas_pi'           => '11',
 'FrInel_pi'           => '11',
 'FrAbs_pi'            => '11',
 'FrPiProd_pi'         => '11',
 'FrCEx_N'             => '11',
 'FrElas_N'            => '11',
 'FrInel_N'            => '11',
 'FrAbs_N'             => '11',
 'FrPiProd_N'          => '11',
 'CCQEPauliSupViaKF'   => '11',
 'CCQEMomDistroFGtoSF' => '11',
 'RDecBR1gamma'        =>  '3',
 'RDecBR1eta'          =>  '3',
 'Theta_Delta2Npi'     =>  '7'
);

%used_in_t2k_analysis = ( 
 'MaNCEL'              => '1',
 'EtaNCEL'             => '1',
 'NormCCQE'            => '0',
 'MaCCQE'              => '1',
 'MaCCQEshape'         => '0',
 'VecFFCCQEshape'      => '1',
 'NormCCRES'           => '0',
 'MaCCRESshape'        => '0',
 'MvCCRESshape'        => '0',
 'MaCCRES'             => '1',
 'MvCCRES'             => '1',
 'NormNCRES'           => '0',
 'MaNCRESshape'        => '0',
 'MvNCRESshape'        => '0',
 'MaNCRES'             => '1',
 'MvNCRES'             => '1',
 'MaCOHpi'             => '1',
 'R0COHpi'             => '1',
 'NonRESBGvpCC1pi'     => '1',
 'NonRESBGvpCC2pi'     => '0',
 'NonRESBGvpNC1pi'     => '1',
 'NonRESBGvpNC2pi'     => '0',
 'NonRESBGvnCC1pi'     => '1',
 'NonRESBGvnCC2pi'     => '0',
 'NonRESBGvnNC1pi'     => '1',
 'NonRESBGvnNC2pi'     => '0',
 'NonRESBGvbarpCC1pi'  => '0',
 'NonRESBGvbarpCC2pi'  => '0',
 'NonRESBGvbarpNC1pi'  => '0',
 'NonRESBGvbarpNC2pi'  => '0',
 'NonRESBGvbarnCC1pi'  => '0',
 'NonRESBGvbarnCC2pi'  => '0',
 'NonRESBGvbarnNC1pi'  => '0',
 'NonRESBGvbarnNC2pi'  => '0',
 'AhtBY'               => '1',
 'BhtBY'               => '1',
 'CV1uBY'              => '1',
 'CV2uBY'              => '1',
 'AhtBYshape'          => '0',
 'BhtBYshape'          => '0',
 'CV1uBYshape'         => '0',
 'CV2uBYshape'         => '0',
 'NormDISCC'           => '0',
 'RnubarnuCC'          => '0',
 'DISNuclMod'          => '0',
 'AGKYxF1pi'           => '1',
 'AGKYpT1pi'           => '1',
 'FormZone'            => '1',
 'MFP_pi'              => '1',
 'MFP_N'               => '1',
 'FrCEx_pi'            => '1',
 'FrElas_pi'           => '1',
 'FrInel_pi'           => '1',
 'FrAbs_pi'            => '1',
 'FrPiProd_pi'         => '1',
 'FrCEx_N'             => '1',
 'FrElas_N'            => '1',
 'FrInel_N'            => '1',
 'FrAbs_N'             => '1',
 'FrPiProd_N'          => '1',
 'CCQEPauliSupViaKF'   => '1',
 'CCQEMomDistroFGtoSF' => '1',
 'RDecBR1gamma'        => '1',
 'RDecBR1eta'          => '1',
 'Theta_Delta2Npi'     => '1'
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

#
# submit weight calculation jobs
#

# systematics loop
for my $curr_syst (keys %def_ntwkdials)  {

  # check whether to commit weight calculation run for current systematic
# print "checking whether to submit job for systematic: $curr_syst \n";

  $using_a_std_analysis_syst_list = 
         ($syst eq "standard-t2k-analysis-list");
  $syst_used_in_given_analysis = 
         ($syst eq "standard-t2k-analysis-list" && $used_in_t2k_analysis{$curr_syst} eq '1');

  $do_submit = 
       ((($syst=~m/$curr_syst/) || ($syst eq "all")) && ! $using_a_std_analysis_syst_list) ||
       $syst_used_in_given_analysis;

  if($do_submit) {
    print "** submitting weight calculation job \n";

    $fntemplate    = "$jobs_dir/syst-$curr_syst";
    $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
    $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
    $ntwkdials     = $def_ntwkdials{$curr_syst};
    $genie_cmd     = "$gexec -f $input_events -s $curr_syst -t $ntwkdials";
    $genie_cmd     = "$genie_cmd -n $nevents"        if defined $nevents;        # add optional argument
    $genie_cmd     = "$genie_cmd -o $output_weights" if defined $output_weights; # add optional argument
    $genie_cmd     = "$genie_cmd | $grep_pipe &> $fntemplate.log";

    print "@@ exec: $genie_cmd \n";

    #
    # submit
    #

    # PBS case
    if($batch_system eq 'PBS') {
       $batch_script  = "$fntemplate.pbs";
       open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
       print PBS "#!/bin/bash \n";
       print PBS "#PBS -N syst-$curr_syst \n";
       print PBS "#PBS -l cput=$time_limit \n";
       print PBS "#PBS -o $fntemplate.pbsout.log \n";
       print PBS "#PBS -e $fntemplate.pbserr.log \n";
       print PBS "source $genie_setup \n"; 
       print PBS "cd $jobs_dir \n";
       print PBS "$genie_cmd \n";
       close(PBS);
       `qsub -q $queue $batch_script`;
    } #PBS

    # LSF case
    if($batch_system eq 'LSF') {
       $batch_script  = "$fntemplate.sh";
       open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
       print LSF "#!/bin/bash \n";
       print LSF "#BSUB-j syst-$curr_syst \n";
       print LSF "#BSUB-q $queue \n";
       print LSF "#BSUB-c $time_limit \n";
       print LSF "#BSUB-o $fntemplate.lsfout.log \n";
       print LSF "#BSUB-e $fntemplate.lsferr.log \n";
       print LSF "source $genie_setup \n"; 
       print LSF "cd $jobs_dir \n";
       print LSF "$genie_cmd \n";
       close(LSF);
       `qsub < $batch_script`;
    } #LSF

    # slurm case
    if($batch_system eq 'slurm') {
       my $time_lim = `sinfo -h -p batch -o %l`;
       my ($days, $hours, $remainder) = $time_lim =~ /([0]+)-([0-9]+):(.*)/;
       my $newhours = $days * 24 + $hours;
       my $new_time_lim = "$newhours:$remainder";
       $time_limit = $new_time_lim lt $time_limit ? $new_time_lim : $time_limit;
       $batch_script  = "$fntemplate.sh";
       open(SLURM, ">$batch_script") or die("Can not create the SLURM batch script");
       print SLURM "#!/bin/bash \n";
       print SLURM "#SBATCH-p $queue \n";
       print SLURM "#SBATCH-o $fntemplate.lsfout.log \n";
       print SLURM "#SBATCH-e $fntemplate.lsferr.log \n";
       print SLURM "#SBATCH-t $time_limit \n";
       print SLURM "source $genie_setup \n"; 
       print SLURM "cd $jobs_dir \n";
       print SLURM "$genie_cmd \n";
       close(SLURM);
       `sbatch --job-name=syst-$curr_syst $batch_script`;
    } #slurm

    # no batch system, run jobs interactively
    if($batch_system eq 'none') {
       system("source $genie_setup; cd $jobs_dir; $genie_cmd");
    } # interactive mode

  } #checking whether to submit current systematic
} # loop over systematics

