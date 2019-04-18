#!/usr/bin/perl

#
#-----------------------------------------------------------------------------------------------------------
# Submit jobs to generate data needed for validating GENIE (e,e') differential cross-section modelling.
#
# Syntax:
#   perl submit_eA_xsec_validation_mc_jobs.pl <options>
#
# Options:
#    --version       : GENIE version number
#    --run           : runs to submit (eg --run 100601200200 / --run 100601200200,100801600680 / --run all / --run 12C)
#   [--model-enum]   : physics model enumeration, default: 0
#   [--nsubruns]     : number of subruns per run, default: 1
# @ [--Emin]         : minimum energy (eg `--run 56Fe --Emin 1.500' will only run the standard iron jobs above 1.5 GeV)
#   [--arch]         : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]   : production name, default: <model>_<version>
#   [--cycle]        : cycle in current production, default: 01
#   [--use-valgrind] : default: off
#   [--batch-system] : <PBS, LSF, slurm, none>, default: PBS
#   [--queue]        : default: prod
#   [--softw-topdir] : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]  : top level dir for job files, default: /opt/ppd/t2k/scratch/GENIE/
#
# Run number key: IZZZAAAEEEEEMxx
#
# I     : probe  (1:e-, 2:e+)
# ZZZ   : target Z (eg, 001:H1, 026:Fe56)
# AAA   : target A (eg: 056:Fe56, 238:U)
# EEEEE : energy used in MeV (eg 00680->0.68GeV, 02015->2.015GeV etc) 
# M     : physics model enumeration, 0-9
# xx    : sub-run ID, 00-99, 100k events each
#
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
  if($_ eq '--Emin')          { $Emin          = $ARGV[$iarg+1]; }
  if($_ eq '--nsubruns')      { $nsubruns      = $ARGV[$iarg+1]; }
  if($_ eq '--run')           { $run           = $ARGV[$iarg+1]; }
  if($_ eq '--model-enum')    { $model_enum    = $ARGV[$iarg+1]; }
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')  { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir  = $ARGV[$iarg+1]; }  
  if($_ eq '--jobs-topdir')   { $jobs_topdir   = $ARGV[$iarg+1]; }   
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;
die("** Aborting [You need to specify which runs to submit. Use the --run option]")
unless defined $run;

$model_enum     = "0"                                     unless defined $model_enum;
$nsubruns       = 1                                       unless defined $nsubruns;
$use_valgrind   = 0                                       unless defined $use_valgrind;
$arch           = "SL6.x86_64"                            unless defined $arch;
$production     = "$model_enum\_$genie_version"           unless defined $production;
$cycle          = "01"                                    unless defined $cycle;
$batch_system   = "PBS"                                   unless defined $batch_system;
$queue          = "prod"                                  unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE"              unless defined $softw_topdir;
$jobs_topdir    = "/opt/ppd/t2k/scratch/GENIE/"           unless defined $jobs_topdir;
$time_limit     = "60:00:00";
$xspl_file      = "$softw_topdir/data/job_inputs/xspl/gxspl-eA-$genie_version.xml";
$genie_setup    = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$jobs_dir       = "$jobs_topdir/eA-$production\_$cycle";
$mcseed         = 210921029;
$nev_per_subrun = 100000;

#
# MC runs required for generating predictions for all datasets described
# data/validation/eA/xsec/differential/qe/datasets.txt
#
@std_runnu = 
(
# 2D
#  ...
#  ... to be added
#  ...
# 3H
#  ...
#  ... to be added
#  ...
# 3He
#  ...
#  ... to be added
#  ...
# 4He
#  ...
#  ... to be added
#  ...
# 12C
  100601200160,  # (e-,e-')12C,  0.160 GeV
  100601200200,  # (e-,e-')12C,  0.200 GeV
  100601200240,  # (e-,e-')12C,  0.240 GeV
  100601200280,  # (e-,e-')12C,  0.280 GeV
  100601200320,  # (e-,e-')12C,  0.320 GeV
  100601200361,  # (e-,e-')12C,  0.361 GeV
  100601200400,  # (e-,e-')12C,  0.400 GeV
  100601200440,  # (e-,e-')12C,  0.440 GeV
  100601200480,  # (e-,e-')12C,  0.480 GeV
  100601200500,  # (e-,e-')12C,  0.500 GeV
  100601200519,  # (e-,e-')12C,  0.519 GeV
  100601200560,  # (e-,e-')12C,  0.560 GeV
  100601200620,  # (e-,e-')12C,  0.620 GeV
  100601200680,  # (e-,e-')12C,  0.680 GeV
  100601200730,  # (e-,e-')12C,  0.730 GeV
  100601200961,  # (e-,e-')12C,  0.961 GeV
  100601201108,  # (e-,e-')12C,  1.108 GeV
  100601201300,  # (e-,e-')12C,  1.300 GeV
  100601201500,  # (e-,e-')12C,  1.500 GeV
  100601201650,  # (e-,e-')12C,  1.650 GeV
  100601201930,  # (e-,e-')12C,  1.930 GeV
  100601202000,  # (e-,e-')12C,  2.000 GeV
  100601202015,  # (e-,e-')12C,  2.015 GeV
  100601202020,  # (e-,e-')12C,  2.020 GeV
  100601202130,  # (e-,e-')12C,  2.130 GeV
  100601202500,  # (e-,e-')12C,  2.500 GeV
  100601202700,  # (e-,e-')12C,  2.700 GeV
  100601203188,  # (e-,e-')12C,  3.188 GeV
  100601203595,  # (e-,e-')12C,  3.595 GeV
  100601203605,  # (e-,e-')12C,  3.605 GeV
  100601204045,  # (e-,e-')12C,  4.045 GeV
  100601204212,  # (e-,e-')12C,  4.212 GeV
  100601205120,  # (e-,e-')12C,  5.120 GeV
# 16O
  100801600700,  # (e-,e-')16O,  0.700 GeV
  100801600737,  # (e-,e-')16O,  0.737 GeV
  100801600880,  # (e-,e-')16O,  0.880 GeV
  100801601080,  # (e-,e-')16O,  1.080 GeV
  100801601200,  # (e-,e-')16O,  1.200 GeV
  100801601500,  # (e-,e-')16O,  1.500 GeV
#  ...
#  ... to be added
#  ...
# 27Al
#  ...
#  ... to be added
#  ...
# 40Ca
#  ...
#  ... to be added
#  ...
# 48Ca
#  ...
#  ... to be added
#  ...
# 56Fe
  102605600120,  # (e-,e-')56Fe,  0.120 GeV
  102605600155,  # (e-,e-')56Fe,  0.155 GeV
  102605600160,  # (e-,e-')56Fe,  0.160 GeV
  102605600174,  # (e-,e-')56Fe,  0.174 GeV
  102605600192,  # (e-,e-')56Fe,  0.192 GeV
  102605600200,  # (e-,e-')56Fe,  0.200 GeV
  102605600206,  # (e-,e-')56Fe,  0.206 GeV
  102605600221,  # (e-,e-')56Fe,  0.221 GeV
  102605600236,  # (e-,e-')56Fe,  0.236 GeV
  102605600240,  # (e-,e-')56Fe,  0.240 GeV
  102605600252,  # (e-,e-')56Fe,  0.252 GeV
  102605600271,  # (e-,e-')56Fe,  0.271 GeV
  102605600280,  # (e-,e-')56Fe,  0.280 GeV
  102605600293,  # (e-,e-')56Fe,  0.293 GeV
  102605600312,  # (e-,e-')56Fe,  0.312 GeV
  102605600320,  # (e-,e-')56Fe,  0.320 GeV
  102605600333,  # (e-,e-')56Fe,  0.333 GeV
  102605600340,  # (e-,e-')56Fe,  0.340 GeV
  102605600360,  # (e-,e-')56Fe,  0.360 GeV
  102605600400,  # (e-,e-')56Fe,  0.400 GeV
  102605600440,  # (e-,e-')56Fe,  0.440 GeV
  102605600480,  # (e-,e-')56Fe,  0.480 GeV
  102605600560,  # (e-,e-')56Fe,  0.560 GeV
  102605600620,  # (e-,e-')56Fe,  0.620 GeV
  102605600900,  # (e-,e-')56Fe,  0.900 GeV
  102605600961,  # (e-,e-')56Fe,  0.961 GeV
  102605601100,  # (e-,e-')56Fe,  1.100 GeV
  102605601108,  # (e-,e-')56Fe,  1.108 GeV
  102605601250,  # (e-,e-')56Fe,  1.250 GeV
  102605601299,  # (e-,e-')56Fe,  1.299 GeV
  102605601500,  # (e-,e-')56Fe,  1.500 GeV
  102605601650,  # (e-,e-')56Fe,  1.650 GeV
  102605602015,  # (e-,e-')56Fe,  2.015 GeV
  102605602020,  # (e-,e-')56Fe,  2.020 GeV
  102605602700,  # (e-,e-')56Fe,  2.700 GeV
  102605603188,  # (e-,e-')56Fe,  3.188 GeV
  102605603300,  # (e-,e-')56Fe,  3.300 GeV
  102605603595,  # (e-,e-')56Fe,  3.595 GeV
  102605603600,  # (e-,e-')56Fe,  3.600 GeV
  102605603605,  # (e-,e-')56Fe,  3.605 GeV
  102605603900,  # (e-,e-')56Fe,  3.900 GeV
  102605603995,  # (e-,e-')56Fe,  3.995 GeV
  102605604045,  # (e-,e-')56Fe,  4.045 GeV
  102605604212,  # (e-,e-')56Fe,  4.212 GeV
  102605604300,  # (e-,e-')56Fe,  4.300 GeV
  102605605120   # (e-,e-')56Fe,  5.120 GeV
#  ...
#  ... to be added
#  ...
# 197Au
#  ...
#  ... to be added
#  ...
# 208Pb
#  ...
#  ... to be added
#  ...
# 238U
#  ...
#  ... to be added
#  ...
);

@runnu = ();
if($run eq "all") {  
  push(@runnu, @std_runnu);  
}
else {  
  my @a = split(',', $run); 
  push(@runnu, @a);  
}

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

#
# submit event generation jobs
#

# run loop
foreach(@runnu) {
    my $curr_runnu = $_;
    print "** submitting event generation run: $curr_runnu \n";

    # match probe,Z,A,E from current run
    my ($m_probe, $m_Z, $m_A, $m_E) = $curr_runnu =~ /([0-9]{1})([0-9]{3})([0-9]{3})([0-9]{5})/;

    # convert to GENIE-expected probe and target PDG codes and energy
    $tgt_pdg = 1000000000 + $m_Z*10000 + $m_A*10;
    $E = $m_E/1000.;

    if    ($m_probe == 1) { $probe_pdg =  11; }
    elsif ($m_probe == 2) { $probe_pdg = -11; }
    else                  { die; }

    # submit subruns
    for($isubrun = 0; $isubrun < $nsubruns; $isubrun++) {

       # Run number key: IZZZAAAEEEEEMxx
       $curr_subrunnu = 1000 * $curr_runnu + 100 * $model_enum + $isubrun;

       $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $fntemplate    = "$jobs_dir/eA-$curr_subrunnu";
       $curr_seed     = $mcseed + $isubrun;
       $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $evgen_opt     = "-n $nev_per_subrun -e $E -p $probe_pdg -t $tgt_pdg -r $curr_subrunnu --seed $curr_seed --cross-sections $xspl_file --event-generator-list EM";
       $evgen_cmd     = "gevgen $evgen_opt | $grep_pipe &> $fntemplate.evgen.log";
       $conv_cmd      = "gntpc -f gst -i gntp.$curr_subrunnu.ghep.root | $grep_pipe &> $fntemplate.conv.log";

       print "@@ exec: $evgen_cmd \n";

       #
       # submit
       #
  
       # PBS case
       if($batch_system eq 'PBS') {
          $batch_script  = "$fntemplate.pbs";
          open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
          print PBS "#!/bin/bash \n";
          print PBS "#PBS -N nucl-$curr_subrunnu \n";
          print PBS "#PBS -l cput=$time_limit \n";
          print PBS "#PBS -o $fntemplate.pbsout.log \n";
          print PBS "#PBS -e $fntemplate.pbserr.log \n";
          print PBS "source $genie_setup \n"; 
          print PBS "cd $jobs_dir \n";
          print PBS "$evgen_cmd \n";
          print PBS "$conv_cmd \n";
          close(PBS);
          `qsub -q $queue $batch_script`;
       } # PBS

       # LSF case
       if($batch_system eq 'LSF') {
          $batch_script  = "$fntemplate.sh";
          open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
          print LSF "#!/bin/bash \n";
          print LSF "#BSUB-j nucl-$curr_subrunnu \n";
          print LSF "#BSUB-q $queue \n";
          print LSF "#BSUB-c $time_limit \n";
          print LSF "#BSUB-o $fntemplate.lsfout.log \n";
          print LSF "#BSUB-e $fntemplate.lsferr.log \n";
          print LSF "source $genie_setup \n"; 
          print LSF "cd $jobs_dir \n";
          print LSF "$evgen_cmd \n";
          print LSF "$conv_cmd \n";
          close(LSF);
          `bsub < $batch_script`;
       } # LSF

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
          print SLURM "$evgen_cmd \n";
          print SLURM "$conv_cmd \n";
          close(SLURM);
          `sbatch --job-name=nucl-$curr_subrunnu $batch_script`;
       } # slurm

       # no batch system, run jobs interactively
       if($batch_system eq 'none') {
          system("source $genie_setup; cd $jobs_dir; $evgen_cmd; $conv_cmd");
       } # interactive mode

    } # loop over subruns

} # loop over runs

