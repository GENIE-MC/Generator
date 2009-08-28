#!/usr/bin/perl

#--------------------------------------------------------------------------------------
# Submit jobs for calculating GENIE free-nucleon cross section splines.
#
# For use at the RAL/PPD Tier2 PBS batch farm.
#
# Syntax:
#   shell% perl submit-xsec_freenuc.pl <options>
#
# Options:
#    --xsplset       : set of splines to generate
#    --version       : genie version number
#   [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]   : default: freenucxspl_<version>
#   [--cycle]        : default: 01
#   [--use-valgrind] : default: off
#   [--queue]        : default: prod
#   [--softw-topdir] : default: /opt/ppd/t2k/GENIE
#
# Examples:
#   shell% perl submit-xsection_jobs-freenuc.pl --version v2.6.0 --xsplset chm 
#   shell% perl submit-xsection_jobs-freenuc.pl --version v2.6.0 --xsplset chm,nue,qel
#   shell% perl submit-xsection_jobs-freenuc.pl --version v2.6.0 --xsplset all 
#
# Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#--------------------------------------------------------------------------------------

use File::Path;

# inputs
#  
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--xsplset')       { $xsplset       = $ARGV[$iarg+1]; }
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir  = $ARGV[$iarg+1]; }           
  $iarg++;
}
die("** Aborting [Undefined set of cross section splines #. Use the --xsplset option]")
unless defined $xsplset;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$use_valgrind   = 0                             unless defined $use_valgrind;
$arch           = "SL5_64bit"                   unless defined $arch;
$production     = "freenucxspl\_$genie_version" unless defined $production;
$cycle          = "01"                          unless defined $cycle;
$queue          = "prod"                        unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/GENIE"          unless defined $softw_topdir;
$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$jobs_dir       = "$softw_topdir/scratch/xsec-$production\_$cycle/";

$nkots = 500;
$emax  = 200;

%NUPDG =   ( 'chm'           => '12,-12,14,-14,16,-16',
             'nue'           => '12,-12,14,-14,16,-16',
             'qel'           => '12,-12,14,-14,16,-16',
             'dis_ebar_cc'   => '-12',
             'dis_ebar_nc'   => '-12',
             'dis_e_cc'      => '12',
             'dis_e_nc'      => '12',
             'dis_mubar_cc'  => '-14',
             'dis_mubar_nc'  => '-14',
             'dis_mu_cc'     => '14',
             'dis_mu_nc'     => '14',
             'dis_taubar_cc' => '-16',
             'dis_taubar_nc' => '-16',
             'dis_tau_cc'    => '16',
             'dis_tau_nc'    => '16',
             'res_ebar_cc'   => '-12',
             'res_ebar_nc'   => '-12',
             'res_e_cc'      => '12',
             'res_e_nc'      => '12',
             'res_mubar_cc'  => '-14',
             'res_mubar_nc'  => '-14',
             'res_mu_cc'     => '14',
             'res_mu_nc'     => '14',
             'res_taubar_cc' => '-16',
             'res_taubar_nc' => '-16',
             'res_tau_cc'    => '16',
             'res_tau_nc'    => '16' );

%TGTPDG =  ( 'chm'           => '1000010010,1000000010',
             'nue'           => '1000010010,1000000010',
             'qel'           => '1000010010,1000000010',
             'dis_ebar_cc'   => '1000010010,1000000010',
             'dis_ebar_nc'   => '1000010010,1000000010',
             'dis_e_cc'      => '1000010010,1000000010',
             'dis_e_nc'      => '1000010010,1000000010',
             'dis_mubar_cc'  => '1000010010,1000000010',
             'dis_mubar_nc'  => '1000010010,1000000010',
             'dis_mu_cc'     => '1000010010,1000000010',
             'dis_mu_nc'     => '1000010010,1000000010',
             'dis_taubar_cc' => '1000010010,1000000010',
             'dis_taubar_nc' => '1000010010,1000000010',
             'dis_tau_cc'    => '1000010010,1000000010',
             'dis_tau_nc'    => '1000010010,1000000010',
             'res_ebar_cc'   => '1000010010,1000000010',
             'res_ebar_nc'   => '1000010010,1000000010',
             'res_e_cc'      => '1000010010,1000000010',
             'res_e_nc'      => '1000010010,1000000010',
             'res_mubar_cc'  => '1000010010,1000000010',
             'res_mubar_nc'  => '1000010010,1000000010',
             'res_mu_cc'     => '1000010010,1000000010',
             'res_mu_nc'     => '1000010010,1000000010',
             'res_taubar_cc' => '1000010010,1000000010',
             'res_taubar_nc' => '1000010010,1000000010',
             'res_tau_cc'    => '1000010010,1000000010',
             'res_tau_nc'    => '1000010010,1000000010' );

%GEVGL =   ( 'chm'           => 'CHARM',
             'nue'           => 'NUE',
             'qel'           => 'QEL',
             'dis_ebar_cc'   => 'DIS-CC',
             'dis_ebar_nc'   => 'DIS-NC',
             'dis_e_cc'      => 'DIS-CC',
             'dis_e_nc'      => 'DIS-NC',
             'dis_mubar_cc'  => 'DIS-CC',
             'dis_mubar_nc'  => 'DIS-NC',
             'dis_mu_cc'     => 'DIS-CC',
             'dis_mu_nc'     => 'DIS-NC',
             'dis_taubar_cc' => 'DIS-CC',
             'dis_taubar_nc' => 'DIS-NC',
             'dis_tau_cc'    => 'DIS-CC',
             'dis_tau_nc'    => 'DIS-NC',
             'res_ebar_cc'   => 'RES-CC',
             'res_ebar_nc'   => 'RES-NC',
             'res_e_cc'      => 'RES-CC',
             'res_e_nc'      => 'RES-NC',
             'res_mubar_cc'  => 'RES-CC',
             'res_mubar_nc'  => 'RES-NC',
             'res_mu_cc'     => 'RES-CC',
             'res_mu_nc'     => 'RES-NC',
             'res_taubar_cc' => 'RES-CC',
             'res_taubar_nc' => 'RES-NC',
             'res_tau_cc'    => 'RES-CC',
             'res_tau_nc'    => 'RES-NC' );

%OUTXML =  ( 'chm'           => 'pgxspl-chm.xml',
             'nue'           => 'pgxspl-nue.xml',
             'qel'           => 'pgxspl-qel.xml',
             'dis_ebar_cc'   => 'pgxspl-dis_ebar_cc.xml',
             'dis_ebar_nc'   => 'pgxspl-dis_ebar_nc.xml',
             'dis_e_cc'      => 'pgxspl-dis_e_cc.xml',
             'dis_e_nc'      => 'pgxspl-dis_e_nc.xml',
             'dis_mubar_cc'  => 'pgxspl-dis_mubar_cc.xml',
             'dis_mubar_nc'  => 'pgxspl-dis_mubar_nc.xml',
             'dis_mu_cc'     => 'pgxspl-dis_mu_cc.xml',
             'dis_mu_nc'     => 'pgxspl-dis_mu_nc.xml',
             'dis_taubar_cc' => 'pgxspl-dis_taubar_cc.xml',
             'dis_taubar_nc' => 'pgxspl-dis_taubar_nc.xml',
             'dis_tau_cc'    => 'pgxspl-dis_tau_cc.xml',
             'dis_tau_nc'    => 'pgxspl-dis_tau_nc.xml',
             'res_ebar_cc'   => 'pgxspl-res_ebar_cc.xml',
             'res_ebar_nc'   => 'pgxspl-res_ebar_nc.xml',
             'res_e_cc'      => 'pgxspl-res_e_cc.xml',
             'res_e_nc'      => 'pgxspl-res_e_nc.xml',
             'res_mubar_cc'  => 'pgxspl-res_mubar_cc.xml',
             'res_mubar_nc'  => 'pgxspl-res_mubar_nc.xml',
             'res_mu_cc'     => 'pgxspl-res_mu_cc.xml',
             'res_mu_nc'     => 'pgxspl-res_mu_nc.xml',
             'res_taubar_cc' => 'pgxspl-res_taubar_cc.xml',
             'res_taubar_nc' => 'pgxspl-res_taubar_nc.xml',
             'res_tau_cc'    => 'pgxspl-res_tau_cc.xml',
             'res_tau_nc'    => 'pgxspl-res_tau_nc.xml' );

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

for my $curr_xsplset (keys %OUTXML)  {
  if($xsplset=~m/$curr_xsplset/ || $xsplset eq "all") {
    #
    # get runnu-dependent info 
    #
    $nu     = $NUPDG   {$curr_xsplset};
    $tgt    = $TGTPDG  {$curr_xsplset};
    $gevgl  = $GEVGL   {$curr_xsplset};
    $outxml = $OUTXML  {$curr_xsplset};

    #
    # create the PBS script 
    #
    $BATCH_SCRIPT = "$jobs_dir/job-$curr_xsplset.pbs";
    open(PBS, ">$BATCH_SCRIPT") or die("Can not create the PBS batch script");

    $logfile_pbse  = "$jobs_dir/job-xspl-$curr_xsplset.pbs_e.log";
    $logfile_pbso  = "$jobs_dir/job-xspl-$curr_xsplset.pbs_o.log";
    $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
    $cmd           = "gmkspl -p $nu -t $tgt -n $nkots -e $emax -o $outxml &> job-$curr_xsplset.log";

    print PBS "#!/bin/bash \n";
    print PBS "#PBS -o $logfile_pbso \n";
    print PBS "#PBS -e $logfile_pbse \n";
    print PBS "source $genie_setup \n";
    print PBS "cd $jobs_dir \n";
    print PBS "unset GSPLOAD \n";
    print PBS "export GEVGL=$gevgl \n";
    print PBS "$cmd \n";

    print "EXEC: $cmd \n";

    #
    # submit job 
    #
    `qsub -q $queue $BATCH_SCRIPT`;
  }
}
