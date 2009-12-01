#!/usr/bin/perl

#-----------------------------------------------------------------------------------------------------------
# Submit jobs to generate data needed for validating GENIE's nuclear modelling
# giving empashis on comparisons with eletron QE scattering data.
#
# The generated data can be fed into GENIE's `gvld_e_qel_xsec' utility.
#
# Syntax:
#   perl submit-vld_nuclmod.pl <options>
#
# Options:
#  --version       : GENIE version number
#  --run           : runs to submit (eg --run 10060120680 / --run 10060120680,10060122000 / -run all)
# [--model-enum]   : physics model enumeration, default: 00
# [--nsubruns]     : number of subruns per run, default: 1
# [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
# [--production]   : production name, default: hadnucvld_<version>
# [--cycle]        : cycle in current production, default: 01
# [--use-valgrind] : default: off
# [--system]       : <RAL_tier2, Imperial_batch, Interactive>, default: RAL_tier2
# [--queue]        : default: prod
# [--softw-topdir] : default: /opt/ppd/t2k/GENIE
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
# Run number key: IZZZAAAEEEEMMxxx
#
# I    : probe (1: e-)
# ZZZ  : nuclear target atomic number (eg 026 for Fe56)
# AAA  : nuclear target mass   number (eg 056 for Fe56)
# EEEE : energy used in MeV (eg 0680->0.68GeV, 2015->2.015GeV etc) 
# MM   : physics model enumeration, 00-99
# xxx  : sub-run ID, 000-999, 50k events each
#
#.......................................................................
# run number         |  init state      | energy   | GEVGL             | 
#                    |                  | (GeV)    | setting           |
#.......................................................................
# 10060120680MMxxx   | e-    + C12      | 0.680    | EM                | 
# 10060121501MMxxx   | e-    + C12      | 1.501    | EM                | 
# 10060122000MMxxx   | e-    + C12      | 2.000    | EM                | 
#.......................................................................
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--nsubruns')      { $nsubruns      = $ARGV[$iarg+1]; }
  if($_ eq '--runnu')         { $runnu         = $ARGV[$iarg+1]; }
  if($_ eq '--model-enum')    { $model_enum    = $ARGV[$iarg+1]; }
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--system')        { $system        = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir  = $ARGV[$iarg+1]; }  
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;
die("** Aborting [You need to specify which runs to submit. Use the --run option]")
unless defined $runnu;

$nsubruns       = 1                                       unless defined $nsubruns;
$use_valgrind   = 0                                       unless defined $use_valgrind;
$arch           = "SL5_64bit"                             unless defined $arch;
$production     = "nuclmod\_$model_enum\_$genie_version"  unless defined $production;
$cycle          = "01"                                    unless defined $cycle;
$queue          = "prod"                                  unless defined $queue;
$system         = "RAL_tier2"                             unless defined $system;
$softw_topdir   = "/opt/ppd/t2k/GENIE"                    unless defined $softw_topdir;
$time_limit     = "60:00:00";
$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$jobs_dir       = "$softw_topdir/scratch/$production\_$cycle";
$xspl_file      = "$softw_topdir/data/job_inputs/xspl/gxspl-emode-$genie_version.xml";
$mcseed         = 210921029;
$nev_per_subrun = 20000;

# inputs for event generation jobs
%evg_pdg_hash = ( 
  '10060120680' =>   '11',
  '10060121501' =>   '11',
  '10060122000' =>   '11'
);
%evg_tgtpdg_hash = ( 
  '10060120680' =>   '1000060120',
  '10060121501' =>   '1000060120',
  '10060122000' =>   '1000060120'
);
%evg_energy_hash = ( 
  '10060120680' =>   '0.680',
  '10060121501' =>   '1.501',
  '10060122000' =>   '2.000'
);
%evg_gevgl_hash = ( 
  '10060120680' =>   'EM',
  '10060121501' =>   'EM',
  '10060122000' =>   'EM'
);
%evg_fluxopt_hash = ( 
  '10060120680' =>   '',
  '10060121501' =>   '',
  '10060122000' =>   ''
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

       $curr_subrunnu = 100000 * $curr_runnu + 1000 * $model_enum + $isubrun;

       $logfile_evgen = "$jobs_dir/nucl-$curr_subrunnu.evgen.log";
       $logfile_conv  = "$jobs_dir/nucl-$curr_subrunnu.conv.log";

       $curr_seed     = $mcseed + $isubrun;
       $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $evgen_cmd     = "gevgen -n $nev_per_subrun -s -e $en -p $probe -t $tgt -r $curr_subrunnu $fluxopt | grep_pipe &> $logfile_evgen";
       $conv_cmd      = "gntpc -f gst -i gntp.$curr_subrunnu.ghep.root | grep_pipe &> $logfile_conv";

       #
       # specifics for the RAL tier2
       #
       if($system eq "RAL_tier2") {
          $batch_script  = "$jobs_dir/nucl-$curr_subrunnu.pbs";
          $logfile_pbse  = "$jobs_dir/nucl-$curr_subrunnu.pbs_e.log";
          $logfile_pbso  = "$jobs_dir/nucl-$curr_subrunnu.pbs_o.log";

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
          print PBS "$conv_cmd \n";
          close(PBS);

          print "EXEC: $evgen_cmd \n";

          # submit job
          `qsub -q $queue $batch_script`;
       } # RAL tier2

       #
       # specifics for the Imperial batch system
       #
       if($system eq "Imperial_batch") {

       }

       #
       # become unpopular:
       # run jobs at the interactive front-end
       #
       if($system eq "Interactive") {
          `source $genie_setup;
           cd $jobs_dir;
           export GSPLOAD=$xspl_file;
           export GEVGL=$gevgl;
           export GSEED=$curr_seed;
           $evgen_cmd;
           $conv_cmd`;
       }

    } # loop over subruns

 } #checking whether to submit current run
} # loop over runs

