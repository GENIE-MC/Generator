#!/usr/bin/perl

#-----------------------------------------------------------------------------------------------------------
# Submit jobs to generate the GENIE samples needed for tuning / validating the GENIE MEC generator.
#
# Syntax:
#   perl submit-vld_mec.pl <options>
#
# Options:
#    --version       : GENIE version number
#    --run           : runs to submit (eg --run 19090 / --run 19090,18080 / -run all)
#   [--model-enum]   : physics model enumeration, default: 0
#   [--nsubruns]     : number of subruns per run, default: 1
#   [--offset]       : subrun offset (for augmenting existing sample), default: 0
#   [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]   : production name, default: <version>
#   [--cycle]        : cycle in current production, default: 01
#   [--use-valgrind] : default: off
#   [--batch-system] : <PBS, LSF, slurm>, default: PBS
#   [--queue]        : default: prod
#   [--softw-topdir] : default: /opt/ppd/t2k/softw/GENIE
#
# Tested at the RAL/PPD Tier2 PBS batch farm.
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#-----------------------------------------------------------------------------------------------------------
#
# EVENT SAMPLES:
#
# Run number key: ITTJJMxxx
# I   :  probe (1: numu, 2: numubar)
# TT  :  nuclear target (06: Carbon, 08: Oxygen, 26: Iron, 80: Water, 90: MiniBooNE target mix)
# JJ  :  flux setting (01: 1 GeV, 80: T2K/SK numu flux, 81: T2K/SK numubar flux, 90: MiniBooNE numu flux
# M   :  model enumeration
# xxx :  sub-run ID, 000-999, 50k events each
#
#.............................................................................................
# run number     |  init state                   | event generator | flux
#                |                               | list            |
#.............................................................................................
#
# ...
# 10601Mxxx      | numu     + 12C                | CCMEC+CCQE      | 1 GeV monoenergetic
# 10801Mxxx      | numu     + 16O                | CCMEC+CCQE      | 1 GeV monoenergetic
# 12601Mxxx      | numu     + 56Fe               | CCMEC+CCQE      | 1 GeV monoenergetic
# 18080Mxxx      | numu     + H2O                | CCMEC+CCQE      | T2K/SK numu
# 28081Mxxx      | numubar  + H2O                | CCMEC+CCQE      | T2K/SK numubar
# 19090Mxxx      | numu     + MiniBooNE tgt mix  | CCMEC+CCQE      | MiniBooNE numu
# ...
#
#.............................................................................................
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--run')           { $runnu         = $ARGV[$iarg+1]; }
  if($_ eq '--model-enum')    { $model_enum    = $ARGV[$iarg+1]; }
  if($_ eq '--nsubruns')      { $nsubruns      = $ARGV[$iarg+1]; }
  if($_ eq '--offset')        { $offset        = $ARGV[$iarg+1]; }
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

$model_enum     = "0"                         unless defined $model_enum;
$nsubruns       = 1                           unless defined $nsubruns;
$offset         = 0                           unless defined $offset;
$use_valgrind   = 0                           unless defined $use_valgrind;
$arch           = "SL5_64bit"                 unless defined $arch;
$production     = "$genie_version"            unless defined $production;
$cycle          = "01"                        unless defined $cycle;
$batch_system   = "PBS"                       unless defined $batch_system;
$queue          = "prod"                      unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE"  unless defined $softw_topdir;
$time_limit     = "60:00:00";
$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$flux_dir       = "$softw_topdir/builds/$arch/$genie_version/data/flux/";
$jobs_dir       = "$softw_topdir/scratch/vld\_mec-$production\_$cycle";
$mcseed         = 210921029;
$nev_per_subrun = 50000;

# inputs for event generation jobs
%evg_pdg_hash = ( 
  '10601' =>   '14',
  '10801' =>   '14',
  '12601' =>   '14',
  '18080' =>   '14',
  '28081' =>   '-14',
  '19090' =>   '14'
);
%evg_tgtpdg_hash = ( 
  '10601' =>   '1000060120',
  '10801' =>   '1000080160',
  '12601' =>   '1000260560',
  '18080' =>   '1000080160[0.889],1000010010[0.111]',
  '28081' =>   '1000080160[0.889],1000010010[0.111]',
  '19090' =>   '1000060120[0.857],1000010010[0.143]'
);
%evg_energy_hash = ( 
  '10601' =>   '1',
  '10801' =>   '1',
  '12601' =>   '1',
  '18080' =>   '0,5',
  '28081' =>   '0,5',
  '19090' =>   '0,5'
);
%evg_gevgl_hash = ( 
  '10601' =>   'CCMEC+CCQE',
  '10801' =>   'CCMEC+CCQE',
  '12601' =>   'CCMEC+CCQE',
  '18080' =>   'CCMEC+CCQE',
  '28081' =>   'CCMEC+CCQE',
  '19090' =>   'CCMEC+CCQE'
);
%evg_fluxopt_hash = ( 
  '10601' =>   '',
  '10801' =>   '',
  '12601' =>   '',
  '18080' =>   '-f $flux_dir/t2ksk.root[numu_flux]',
  '28081' =>   '-f $flux_dir/t2ksk.root[numubar_flux]',
  '19090' =>   '-f $flux_dir/miniboone-numu.data'
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

       # Run number key: ITTJJMxxx
       $curr_subrunnu = 10000 * $curr_runnu + 1000 * $model_enum + $isubrun + $offset;
       $curr_seed     = $mcseed + $isubrun + $offset;
       $fntemplate    = "$jobs_dir/mec-$curr_subrunnu";
       $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $evgen_cmd     = "gevgen -n $nev_per_subrun -e $en -p $probe -t $tgt $fluxopt -r $curr_subrunnu --seed $curr_seed --event-generator-list $gevgl | $grep_pipe &> $fntemplate.evgen.log";
       $conv_cmd      = "gntpc -f gst -i gntp.$curr_subrunnu.ghep.root";

       print "@@ exec: $evgen_cmd \n";

       #
       # submit
       #

       # PBS case
       if($batch_system eq 'PBS') {
           $batch_script  = "$fntemplate.pbs";
           open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
           print PBS "#!/bin/bash \n";
           print PBS "#PBS -N mec-$curr_subrunnu \n";
           print PBS "#PBS -l cput=$time_limit \n";
           print PBS "#PBS -o $fntemplate.pbsout.log \n";
           print PBS "#PBS -e $fntemplate.pbserr.log \n";
           print PBS "source $genie_setup \n"; 
           print PBS "cd $jobs_dir \n";
           print PBS "$evgen_cmd \n";
           print PBS "$conv_cmd \n";
           close(PBS);
           `qsub -q $queue $batch_script`;
       } #PBS

       # LSF case
       if($batch_system eq 'LSF') {
           $batch_script  = "$fntemplate.sh";
           open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
           print LSF "#!/bin/bash \n";
           print LSF "#BSUB-j mec-$curr_subrunnu \n";
           print LSF "#BSUB-q $queue \n";
           print LSF "#BSUB-c $time_limit \n";
           print LSF "#BSUB-o $fntemplate.lsfout.log \n";
           print LSF "#BSUB-e $fntemplate.lsferr.log \n";
           print LSF "source $genie_setup \n"; 
           print LSF "cd $jobs_dir \n";
           print LSF "$evgen_cmd \n";
           print LSF "$conv_cmd \n";
           close(LSF);
           `qsub < $batch_script`;
       } #LSF

      # slurm case
       if($batch_system eq 'slurm') {
           $batch_script  = "$fntemplate.sh";
           open(SLURM, ">$batch_script") or die("Can not create the SLURM batch script");
           print SLURM "#!/bin/bash \n";
           print SLURM "#SBATCH-p $queue \n";
           print SLURM "#SBATCH-o $fntemplate.lsfout.log \n";
           print SLURM "#SBATCH-e $fntemplate.lsferr.log \n";
           print SLURM "source $genie_setup \n"; 
           print SLURM "cd $jobs_dir \n";
           print SLURM "$evgen_cmd \n";
           print SLURM "$conv_cmd \n";
           close(SLURM);
           `sbatch --job-name=mec-$curr_subrunnu $batch_script`;
       } #slurm

    } # loop over subruns
  } #checking whether to submit current run
} # loop over runs

