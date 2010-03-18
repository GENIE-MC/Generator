#!/usr/bin/perl

#-----------------------------------------------------------------------------------------------------------
# Submit jobs to generate data needed for validating GENIE's nuclear modelling
# giving empashis on comparisons with eletron QE scattering data.
# The generated data can be fed into GENIE's `gvld_e_qel_xsec' utility.
#
# Syntax:
#   perl submit-vld_nuclmod.pl <options>
#
# Options:
#    --version       : GENIE version number
#    --run           : runs to submit (eg --run 1060680 / --run 1060680,1062000 / -run all)
#   [--model-enum]   : physics model enumeration, default: 0
#   [--nsubruns]     : number of subruns per run, default: 1
#   [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]   : production name, default: <model>_<version>
#   [--cycle]        : cycle in current production, default: 01
#   [--use-valgrind] : default: off
#   [--batch-system] : <PBS, >, default: PBS
#   [--queue]        : default: prod
#   [--softw-topdir] : default: /opt/ppd/t2k/GENIE
#
# Tested at the RAL/PPD Tier2 PBS batch farm.
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#
# Nick Prouse <nicholas.prouse06 \at imperial.ac.uk>
# Imperial College London
#-----------------------------------------------------------------------------------------------------------
#
# EVENT SAMPLES:
#
# Run number key: ITTEEEEMxx
#
# I    : probe  (1:e-)
# TT   : target (01:H1, 02:D2, 06:C12, 08:O16, 26:Fe56)
# EEEE : energy used in MeV (eg 0680->0.68GeV, 2015->2.015GeV etc) 
# M    : physics model enumeration, 0-9
# xx   : sub-run ID, 00-99, 50k events each
#
#.................................................................
# run number   |  init state      | energy   | GEVGL             | 
#              |                  | (GeV)    | setting           |
#.................................................................
# 1060680Mxx   | e-    + C12      | 0.680    | EM                | 
# 1061501Mxx   | e-    + C12      | 1.501    | EM                | 
# 1062000Mxx   | e-    + C12      | 2.000    | EM                | 
#.................................................................
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--nsubruns')      { $nsubruns      = $ARGV[$iarg+1]; }
  if($_ eq '--run')           { $runnu         = $ARGV[$iarg+1]; }
  if($_ eq '--model-enum')    { $model_enum    = $ARGV[$iarg+1]; }
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')  { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir  = $ARGV[$iarg+1]; }  
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;
die("** Aborting [You need to specify which runs to submit. Use the --run option]")
unless defined $runnu;

$model_enum     = "0"                                     unless defined $model_enum;
$nsubruns       = 1                                       unless defined $nsubruns;
$use_valgrind   = 0                                       unless defined $use_valgrind;
$arch           = "SL5_64bit"                             unless defined $arch;
$production     = "$model_enum\_$genie_version"           unless defined $production;
$cycle          = "01"                                    unless defined $cycle;
$batch_system   = "PBS"                                   unless defined $batch_system;
$queue          = "prod"                                  unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/GENIE"                    unless defined $softw_topdir;
$time_limit     = "60:00:00";
$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$jobs_dir       = "$softw_topdir/scratch/vld\_nuclmod-$production\_$cycle";
$xspl_file      = "$softw_topdir/data/job_inputs/xspl/gxspl-emode-$genie_version.xml";
$mcseed         = 210921029;
$nev_per_subrun = 50000;

# inputs for event generation jobs
%evg_pdg_hash = ( 
  '1060680' =>   '11',
  '1061501' =>   '11',
  '1062000' =>   '11'
);
%evg_tgtpdg_hash = ( 
  '1060680' =>   '1000060120',
  '1061501' =>   '1000060120',
  '1062000' =>   '1000060120'
);
%evg_energy_hash = ( 
  '1060680' =>   '0.680',
  '1061501' =>   '1.501',
  '1062000' =>   '2.000'
);
%evg_gevgl_hash = ( 
  '1060680' =>   'EM',
  '1061501' =>   'EM',
  '1062000' =>   'EM'
);
%evg_fluxopt_hash = ( 
  '1060680' =>   '',
  '1061501' =>   '',
  '1062000' =>   ''
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

#
# submit event generation jobs
#

# run loop
for my $curr_runnu (keys %evg_gevgl_hash)  {

 # check whether to commit current run 
 if($runnu=~m/$curr_runnu/ || $runnu eq "all") {

    print "** submitting event generation run: $curr_runnu \n";

    #
    # get runnu-dependent info
    #
    $probe   = $evg_pdg_hash     {$curr_runnu};
    $tgt     = $evg_tgtpdg_hash  {$curr_runnu};
    $en      = $evg_energy_hash  {$curr_runnu};
    $gevgl   = $evg_gevgl_hash   {$curr_runnu};
    $fluxopt = $evg_fluxopt_hash {$curr_runnu};

    # submit subruns
    for($isubrun = 0; $isubrun < $nsubruns; $isubrun++) {

       # Run number key: ITTEEEEMxx
       $curr_subrunnu = 1000 * $curr_runnu + 100 * $model_enum + $isubrun;
       $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $job_file_base = "$jobs_dir/nucl-$curr_subrunnu";
       $curr_seed     = $mcseed + $isubrun;
       $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $evgen_cmd     = "gevgen -n $nev_per_subrun -s -e $en -p $probe -t $tgt -r $curr_subrunnu $fluxopt | $grep_pipe &> $job_file_base.evgen.log";
       $conv_cmd      = "gntpc -f gst -i gntp.$curr_subrunnu.ghep.root | $grep_pipe &> $job_file_base.conv.log";

       print "@@ exec: $evgen_cmd \n";

       #
       # submit
       #
  
       # PBS case
       if($batch_system eq 'PBS') {
          $batch_script  = "$job_file_base.pbs";
          open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
          print PBS "#!/bin/bash \n";
          print PBS "#PBS -l cput=$time_limit \n";
          print PBS "#PBS -o $job_file_base.pbsout.log \n";
          print PBS "#PBS -e $job_file_base.pbserr.log \n";
          print PBS "source $genie_setup \n"; 
          print PBS "cd $jobs_dir \n";
          print PBS "export GSPLOAD=$xspl_file \n";
          print PBS "export GEVGL=$gevgl \n";
          print PBS "export GSEED=$curr_seed \n";
          print PBS "$evgen_cmd \n";
          print PBS "$conv_cmd \n";
          close(PBS);
          `qsub -q $queue $batch_script`;
       } # PBS

    } # loop over subruns

 } #checking whether to submit current run
} # loop over runs

