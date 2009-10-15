#!/usr/bin/perl

#---------------------------------------------------------------------------------------------------------------------
# Submit jobs to generate all data needed for validating GENIE's cross section model 
#
# The generated data can be fed into GENIE's gvld_nuxsec_vs_world_data utility.
#
# Syntax:
#   perl submit-vld_xsec.pl <options>
#
# Options:
#  --version       : GENIE version number
# [--nsubruns]     : number of subruns per run, default: 1
# [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
# [--production]   : production name, default: xsecvld_<version>
# [--cycle]        : cycle in current production, default: 01
# [--use-valgrind] : default: off
# [--queue]        : default: prod
# [--softw-topdir] : default: /opt/ppd/t2k/GENIE
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#---------------------------------------------------------------------------------------------------------------------
#
# EVENT SAMPLES:
# The following samples will be submitted.
# The output GHEP event trees will be analyzed and converted to GST summary trees.
# The gvld_nuxsec_vs_world_data utility uses the samples to decompose the inclusive cross section
# to cross section for various exclusive channels
#......................................................................
# run number      |  init state      | energy   | processes
#                 |                  | (GeV)    | enabled
#......................................................................
#
# 1000xx          | numu    + n      | 0.1-120. | all
# 1100xx          | numu    + p      | 0.1-120. | all
# 1200xx          | numubar + n      | 0.1-120. | all
# 1300xx          | numubar + p      | 0.1-120. | all
#
# xx : Run ID, 01-99, 100k events each
#......................................................................
#
# CROSS SECTIONS:
# A job will be submitted to run GENIE's gspl2root utility, extract
# cross sections (from the pre-computed XML cross section spline file) 
# and save them into a single ROOT file (xsec.root).
# The file will contain the following folders (TDirectories):
# - nu_mu_n
# - nu_mu_H1
# - nu_mubar_n
# - nu_mubar_H1
# - nu_mu_Ne20
# - nu_mu_bar_Ne20
# - nu_mu_Al27
# - nu_mu_Si30     (Si30: proxy for freon (average A=30), used for coherent pi production plots only)
# - nu_mu_bar_Si30 ( -//- )
# Each folder will contain graphs for all cross sections for the given initial state.
# Expand as required if more published data are inluded in this GENIE physics benchmark test.
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

$nsubruns       = 1                         unless defined $nsubruns;
$use_valgrind   = 0                         unless defined $use_valgrind;
$arch           = "SL5_64bit"               unless defined $arch;
$production     = "xsecvld\_$genie_version" unless defined $production;
$cycle          = "01"                      unless defined $cycle;
$queue          = "prod"                    unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/GENIE"      unless defined $softw_topdir;
$time_limit     = "60:00:00";
$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$jobs_dir       = "$softw_topdir/scratch/$production\_$cycle";
$xspl_file      = "$softw_topdir/data/job_inputs/xspl/gxspl-t2k-$genie_version.xml";
$mcseed         = 210921029;
$nev_per_subrun = 100000;

# inputs for event generation jobs
%evg_nupdg_hash = ( 
  '1000' =>   '14',
  '1100' =>   '14',
  '1200' =>  '-14',
  '1300' =>  '-14'
);
%evg_tgtpdg_hash = ( 
  '1000' =>  '1000000010',
  '1100' =>  '1000010010',
  '1200' =>  '1000000010',
  '1300' =>  '1000010010'
);

%evg_energy_hash = ( 
  '1000' =>  '0.1,120.0',
  '1100' =>  '0.1,120.0',
  '1200' =>  '0.1,120.0',
  '1300' =>  '0.1,120.0',
);
%evg_gevgl_hash = ( 
  '1000' =>  'Default',
  '1100' =>  'Default',
  '1200' =>  'Default',
  '1300' =>  'Default',
);
%evg_fluxopt_hash = ( 
  '1000' =>  '-f "1/x"',
  '1100' =>  '-f "1/x"',
  '1200' =>  '-f "1/x"',
  '1300' =>  '-f "1/x"',
);

# inputs for xsec calculation jobs
%xsec_nupdg_hash = ( 
  '1' =>   '14',
  '2' =>   '14',
  '3' =>  '-14',
  '4' =>  '-14',
  '5' =>   '14',
  '6' =>  '-14',
  '7' =>   '14',
  '8' =>   '14',
  '9' =>  '-14'
);
%xsec_tgtpdg_hash = ( 
  '1' =>  '1000000010',
  '2' =>  '1000010010',
  '3' =>  '1000000010',
  '4' =>  '1000010010',
  '5' =>  '1000100200', # Ne20
  '6' =>  '1000100200', # Ne20
  '7' =>  '1000130270', # Al27
  '8' =>  '1000140300', # Si30, proxy for freon (average A=30), used for coherent pi production plots only
  '9' =>  '1000140300'  # Si30, proxy for freon (average A=30), used for coherent pi production plots only
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

#
# submit event generation jobs
#
for my $curr_runnu (keys %evg_gevgl_hash)  {

 # uncomment if you want to include a cmd line input to specify specific runs (and fill-in $runnu)
 # if($runnu=~m/$curr_runnu/ || $runnu eq "all") {

    print "** submitting event generation run: $curr_runnu \n";

    #
    # get runnu-dependent info
    #
    $nu      = $evg_nupdg_hash   {$curr_runnu};
    $tgt     = $evg_tgtpdg_hash  {$curr_runnu};
    $en      = $evg_energy_hash  {$curr_runnu};
    $gevgl   = $evg_gevgl_hash   {$curr_runnu};
    $fluxopt = $evg_fluxopt_hash {$curr_runnu};

    # submit subruns
    for($isubrun = 0; $isubrun < $nsubruns; $isubrun++) {

       $curr_subrunnu = 100 * $curr_runnu + $isubrun;

       $batch_script  = "$jobs_dir/xsecvld-$curr_subrunnu.pbs";
       $logfile_evgen = "$jobs_dir/xsecvld-$curr_subrunnu.evgen.log";
       $logfile_conv  = "$jobs_dir/xsecvld-$curr_subrunnu.conv.log";
       $logfile_pbse  = "$jobs_dir/xsecvld-$curr_subrunnu.pbs_e.log";
       $logfile_pbso  = "$jobs_dir/xsecvld-$curr_subrunnu.pbs_o.log";

       $curr_seed     = $mcseed + $isubrun;
       $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $evgen_cmd     = "gevgen -n $nev_per_subrun -s -e $en -p $nu -t $tgt -r $curr_subrunnu $fluxopt | grep_pipe &> $logfile_evgen";
       $conv_cmd      = "gntpc -f gst -i gntp.$curr_subrunnu.ghep.root | grep -B 100 -A 30 -i \"warn\\|error\\|fatal\" &> $logfile_conv";

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

       #
       # submit job
       #
       `qsub -q $queue $batch_script`;

    } # loop over subruns
 # } #checking whether to submit current run
} # loop over runs


#
# submit job to generate the cross section file
#
{
   $batch_script  = "$jobs_dir/xsecvld-xsec.pbs";
   $logfile_pbse  = "$jobs_dir/xsecvld-xsec.pbs_e.log";
   $logfile_pbso  = "$jobs_dir/xsecvld-xsec.pbs_o.log";

   # create the PBS script
   open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
   print PBS "#!/bin/bash \n";
   print PBS "#PBS -l cput=$time_limit \n";
   print PBS "#PBS -o $logfile_pbso \n";
   print PBS "#PBS -e $logfile_pbse \n";
   print PBS "source $genie_setup \n"; 
   print PBS "unset GEVGL \n"; 
   print PBS "unset GSPLOAD \n"; 
   print PBS "cd $jobs_dir \n";

   for my $curr_init_state (keys %xsec_nupdg_hash)  {
      $nu      = $xsec_nupdg_hash   {$curr_init_state};
      $tgt     = $xsec_tgtpdg_hash  {$curr_init_state};

      $gspl2root_cmd = "gspl2root -p $nu -t $tgt -f $xspl_file -o xsec.root";
      print PBS "$gspl2root_cmd \n";
      print "EXEC: $gspl2root_cmd \n";
   }
   close(PBS);

   #
   # submit job
   #
   `qsub -q $queue $batch_script`;
}
