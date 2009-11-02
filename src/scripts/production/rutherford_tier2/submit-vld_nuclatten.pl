#!/usr/bin/perl

#-----------------------------------------------------------------------------------------------------------
# Submit jobs to generate data needed for validating GENIE's hadronization + intranuclear transport model 
# (with emphasis  on the modeling of medium effects to the hadronization) in a set of comparisons against 
# the CLAS and HERMES nuclear attenuation data.
#
# The generated data can be fed into GENIE's `gvld_medium_effects_hadronz_test' utility.
#
# Syntax:
#   perl submit-vld_nuclatten.pl <options>
#
# Options:
#  --version       : GENIE version number
# [--nsubruns]     : number of subruns per run, default: 1
# [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
# [--production]   : production name, default: hadnucvld_<version>
# [--cycle]        : cycle in current production, default: 01
# [--use-valgrind] : default: off
# [--queue]        : default: prod
# [--softw-topdir] : default: /opt/ppd/t2k/GENIE
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#-----------------------------------------------------------------------------------------------------------
#
# EVENT SAMPLES:
#
# Run number key: IZZZAAAJJxxx
# I   :  1->e-, 2->e+, 3->mu-, 4->mu+
# ZZZ :  nuclear target atomic number (eg 026 for Fe56)
# AAA :  nuclear target mass   number (eg 056 for Fe56)
# JJ  :  flux setting 
#            01 -> 12.0 GeV [HERMES]
#            02 -> 27.6 GeV [HERMES]
# xxx :  sub-run ID, 001-999, 20k events each
#
#...................................................................................
# run number      |  init state      | energy   | GEVGL              | flux
#                 |                  | (GeV)    | setting            |
#...................................................................................
#
# 100100202xxx    | e-    + D2       | 27.6     | EM                | -
# 100200402xxx    | e-    + He4      | 27.6     | EM                | -
# 100701402xxx    | e-    + N14      | 27.6     | EM                | -
# 101002002xxx    | e-    + Ne20     | 27.6     | EM                | -
# 103608302xxx    | e-    + Kr83     | 27.6     | EM                | -
# 105413102xxx    | e-    + Xe131    | 27.6     | EM                | -
#
#...................................................................................
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--nsubruns')      { $nsubruns      = $ARGV[$iarg+1]; }
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir  = $ARGV[$iarg+1]; }  
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$nsubruns       = 1                           unless defined $nsubruns;
$use_valgrind   = 0                           unless defined $use_valgrind;
$arch           = "SL5_64bit"                 unless defined $arch;
$production     = "hadnucvld\_$genie_version" unless defined $production;
$cycle          = "01"                        unless defined $cycle;
$queue          = "prod"                      unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/GENIE"        unless defined $softw_topdir;
$time_limit     = "60:00:00";
$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$jobs_dir       = "$softw_topdir/scratch/$production\_$cycle";
$xspl_file      = "$softw_topdir/data/job_inputs/xspl/gxspl-emode-$genie_version.xml";
$mcseed         = 210921029;
$nev_per_subrun = 20000;

# inputs for event generation jobs
%evg_pdg_hash = ( 
  '100100202' =>   '11',
  '100200402' =>   '11',
  '100701402' =>   '11',
  '101002002' =>   '11',
  '103608302' =>   '11',
  '105413102' =>   '11'
);
%evg_tgtpdg_hash = ( 
  '100100202' =>   '1000010020',
  '100200402' =>   '1000020040',
  '100701402' =>   '1000070140',
  '101002002' =>   '1000100200',
  '103608302' =>   '1000360830',
  '105413102' =>   '1000541310'
);
%evg_energy_hash = ( 
  '100100202' =>   '27.6',
  '100200402' =>   '27.6',
  '100701402' =>   '27.6',
  '101002002' =>   '27.6',
  '103608302' =>   '27.6',
  '105413102' =>   '27.6'
);
%evg_gevgl_hash = ( 
  '100100202' =>   'EM',
  '100200402' =>   'EM',
  '100701402' =>   'EM',
  '101002002' =>   'EM',
  '103608302' =>   'EM',
  '105413102' =>   'EM'
);
%evg_fluxopt_hash = ( 
  '100100202' =>   '',
  '100200402' =>   '',
  '100701402' =>   '',
  '101002002' =>   '',
  '103608302' =>   '',
  '105413102' =>   ''
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

# submit event generation jobs
#
for my $curr_runnu (keys %evg_gevgl_hash)  {

 # uncomment if you want to include a cmd line input to specify specific runs (and fill-in $runnu)
 # if($runnu=~m/$curr_runnu/ || $runnu eq "all") {

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

       $curr_subrunnu = 1000 * $curr_runnu + $isubrun;

       $batch_script  = "$jobs_dir/hdzvld-$curr_subrunnu.pbs";
       $logfile_evgen = "$jobs_dir/hdzvld-$curr_subrunnu.evgen.log";
       $logfile_conv  = "$jobs_dir/hdzvld-$curr_subrunnu.conv.log";
       $logfile_pbse  = "$jobs_dir/hdzvld-$curr_subrunnu.pbs_e.log";
       $logfile_pbso  = "$jobs_dir/hdzvld-$curr_subrunnu.pbs_o.log";

       $curr_seed     = $mcseed + $isubrun;
       $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $evgen_cmd     = "gevgen -n $nev_per_subrun -s -e $en -p $probe -t $tgt -r $curr_subrunnu $fluxopt | grep_pipe &> $logfile_evgen";

       # create the PBS script
       open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
       print PBS "#!/bin/bash \n";
       print PBS "#PBS -l cput=$time_limit \n";
       print PBS "#PBS -o $logfile_pbso \n";
       print PBS "#PBS -e $logfile_pbse \n";
       print PBS "source $genie_setup \n"; 
       print PBS "cd $jobs_dir \n";
       print PBS "export GSPLOAD=$xspl_file \n";
       print PBS "export GEVGL=$gevgl \n";
       print PBS "export GSEED=$curr_seed \n";
       print PBS "$evgen_cmd \n";
       close(PBS);

       print "EXEC: $evgen_cmd \n";

       #
       # submit job
       #
       `qsub -q $queue $batch_script`;

    } # loop over subruns
 # } #checking whether to submit current run
} # loop over runs

