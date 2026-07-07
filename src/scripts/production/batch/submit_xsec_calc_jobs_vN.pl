#!/usr/bin/perl

#----------------------------------------------------------------------------------------------------------------
# Submit jobs for calculating GENIE free-nucleon cross-sections used for constructing cross-section splines. 
# Generate data for all or specified processes, and split the task amongst several worker nodes.
# The GENIE gspladd utility can be used to merge the job outputs.
#
# Syntax:
#   shell% perl submit_xsec_calc_jobs_vN.pl <options>
#
# Options:
#    --gen-version     : GENIE generator version number
#    --tune            : GENIE physics tune
#    --xsec-spline-set : set of cross-section splines to generate
#   [--arch]           : <el7.x86_64, ...>, default: el7.x86_64
#   [--production]     : default: prod
#   [--cycle]          : default: 01
#   [--use-valgrind]   : default: off
#   [--batch-system]   : <Slurm, PBS, LSF, none>, default: Slurm
#   [--queue]          : default: compute
#   [--softw-topdir]   : top level dir for softw installations, default: /user/costasa/projects/GENIE/softw/
#   [--jobs-topdir]    : top level dir for job files, default: /scratch/costasa/GENIE/
#
# Examples:
#   % perl submit_xsec_calc_jobs_vN.pl --gen-version v3.06.00 --tune G18_10a_02_11b --xsec-spline-set chm 
#   % perl submit_xsec_calc_jobs_vN.pl --gen-version v3.06.00 --tune G18_10a_02_11b --xsec-spline-set chm,nue,qel
#   % perl submit_xsec_calc_jobs_vN.pl --gen-version v3.06.00 --tune G18_10a_02_11b --xsec-spline-set all 
#
# Author:
#   Costas Andreopoulos <c.andreopoulos \st cern.ch>
#   University of Liverpool
#
# Copyright:
#   Copyright (c) 2003-2025, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#-------------------------------------------------------------------------------------------------------------

use File::Path;

# inputs
#  
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--gen-version')     { $gen_version   = $ARGV[$iarg+1]; }
  if($_ eq '--tune')            { $tune          = $ARGV[$iarg+1]; }
  if($_ eq '--xsec-spline-set') { $xsplset       = $ARGV[$iarg+1]; }
  if($_ eq '--arch')            { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')      { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')           { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')    { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')    { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')           { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--time-limit')      { $time_limit    = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')    { $softw_topdir  = $ARGV[$iarg+1]; }           
  if($_ eq '--jobs-topdir')     { $jobs_topdir   = $ARGV[$iarg+1]; }           
  $iarg++;
}
die("** Aborting [Undefined set of cross section splines #. Use the --xsplset option]")
unless defined $xsplset;
die("** Aborting [Undefined GENIE Generator version. Use the --gen-version option]")
unless defined $gen_version;
die("** Aborting [Undefined GENIE physics tune. Use the --tune option]")
unless defined $tune;

$use_valgrind     = 0                                       unless defined $use_valgrind;
$arch             = "el7.x86_64"                            unless defined $arch;
$production       = "prod"                                  unless defined $production;
$cycle            = "01"                                    unless defined $cycle;
$batch_system     = "Slurm"                                 unless defined $batch_system;
$queue            = "compute"                               unless defined $queue;
$time_limit       = "10:00:00"                              unless defined $time_limit;  
$softw_topdir     = "/user/costasa/projects/GENIE/softw/"   unless defined $softw_topdir;
$jobs_topdir      = "/scratch/costasa/GENIE/"               unless defined $jobs_topdir;
$gen_setup_script = "$softw_topdir/generator/builds/$arch/$gen-version-setup.sh";
$jobs_dir         = "$jobs_topdir/$production\_$cycle-$gen_version-$tune-XSvN/";
$time_limit       = "10:00:00";

$nkots = 500;
$emax  = 500;

%NUPDG =   ( 'chm'           => '12,-12,14,-14,16,-16',
             'nue'           => '12,-12,14,-14,16,-16',
             'qel'           => '12,-12,14,-14,16,-16',
             'dfr'           => '12,-12,14,-14,16,-16',
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
             'dfr'           => '1000010010,1000000010',
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

%GEVGL =   ( 'chm'           => 'Charm',
             'nue'           => 'NuE',
             'qel'           => 'QE',
             'dfr'           => 'DFR',
             'dis_ebar_cc'   => 'CCDIS',
             'dis_ebar_nc'   => 'NCDIS',
             'dis_e_cc'      => 'CCDIS',
             'dis_e_nc'      => 'NCDIS',
             'dis_mubar_cc'  => 'CCDIS',
             'dis_mubar_nc'  => 'NCDIS',
             'dis_mu_cc'     => 'CCDIS',
             'dis_mu_nc'     => 'NCDIS',
             'dis_taubar_cc' => 'CCDIS',
             'dis_taubar_nc' => 'NCDIS',
             'dis_tau_cc'    => 'CCDIS',
             'dis_tau_nc'    => 'NCDIS',
             'res_ebar_cc'   => 'CCRES',
             'res_ebar_nc'   => 'NCRES',
             'res_e_cc'      => 'CCRES',
             'res_e_nc'      => 'NCRES',
             'res_mubar_cc'  => 'CCRES',
             'res_mubar_nc'  => 'NCRES',
             'res_mu_cc'     => 'CCRES',
             'res_mu_nc'     => 'NCRES',
             'res_taubar_cc' => 'CCRES',
             'res_taubar_nc' => 'NCRES',
             'res_tau_cc'    => 'CCRES',
             'res_tau_nc'    => 'NCRES' );

%OUTXML =  ( 'chm'           => 'partial-xspl-chm.xml',
             'nue'           => 'partial-xspl-nue.xml',
             'qel'           => 'partial-xspl-qel.xml',
             'dfr'           => 'partial-xspl-dfr.xml',
             'dis_ebar_cc'   => 'partial-xspl-dis_ebar_cc.xml',
             'dis_ebar_nc'   => 'partial-xspl-dis_ebar_nc.xml',
             'dis_e_cc'      => 'partial-xspl-dis_e_cc.xml',
             'dis_e_nc'      => 'partial-xspl-dis_e_nc.xml',
             'dis_mubar_cc'  => 'partial-xspl-dis_mubar_cc.xml',
             'dis_mubar_nc'  => 'partial-xspl-dis_mubar_nc.xml',
             'dis_mu_cc'     => 'partial-xspl-dis_mu_cc.xml',
             'dis_mu_nc'     => 'partial-xspl-dis_mu_nc.xml',
             'dis_taubar_cc' => 'partial-xspl-dis_taubar_cc.xml',
             'dis_taubar_nc' => 'partial xspl-dis_taubar_nc.xml',
             'dis_tau_cc'    => 'partial-xspl-dis_tau_cc.xml',
             'dis_tau_nc'    => 'partial-xspl-dis_tau_nc.xml',
             'res_ebar_cc'   => 'partial-xspl-res_ebar_cc.xml',
             'res_ebar_nc'   => 'partial-xspl-res_ebar_nc.xml',
             'res_e_cc'      => 'partial-xspl-res_e_cc.xml',
             'res_e_nc'      => 'partial-xspl-res_e_nc.xml',
             'res_mubar_cc'  => 'partial-xspl-res_mubar_cc.xml',
             'res_mubar_nc'  => 'partial-xspl-res_mubar_nc.xml',
             'res_mu_cc'     => 'partial-xspl-res_mu_cc.xml',
             'res_mu_nc'     => 'partial-xspl-res_mu_nc.xml',
             'res_taubar_cc' => 'partial-xspl-res_taubar_cc.xml',
             'res_taubar_nc' => 'partial-xspl-res_taubar_nc.xml',
             'res_tau_cc'    => 'partial-xspl-res_tau_cc.xml',
             'res_tau_nc'    => 'partial-xspl-res_tau_nc.xml' );

# make the jobs directory
#
print "@@ Creating job directory: $jobs_dir \n";
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

for my $curr_xsplset (keys %OUTXML)  {
  if($xsplset=~m/$curr_xsplset/ || $xsplset eq "all") {

    # Get runnu-dependent info 
    # -------------------------------------------------------------

    $nu     = $NUPDG   {$curr_xsplset};
    $tgt    = $TGTPDG  {$curr_xsplset};
    $gevgl  = $GEVGL   {$curr_xsplset};
    $outxml = $OUTXML  {$curr_xsplset};

    $job_name = "XSvN-$curr_xsplset"; 
    $filename_basepath = "$jobs_dir/$job_name"; 

    $grep_pipe     = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
    $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
    $gmkspl_opt    = "-p $nu -t $tgt -n $nkots -e $emax -o $outxml --event-generator-list $gevgl --tune $tune";
    $gmkspl_cmd    = "gmkspl $gmkspl_opt";

    print "@@ exec: $gmkspl_cmd \n";

    # Submit jobs
    # -------------------------------------------------------------
  
    # Slurm case
    if($batch_system eq 'Slurm') {
        ##my $time_lim = `sinfo -h -p batch -o %l`;
        ##my ($days, $hours, $remainder) = $time_lim =~ /([0]+)-([0-9]+):(.*)/;
        ##my $newhours = $days * 24 + $hours;
        ##my $new_time_lim = "$newhours:$remainder";
        ##$time_limit = $new_time_lim lt $time_limit ? $new_time_lim : $time_limit;
        $batch_script = "$filename_basepath.sh";
        open(SLURM, ">$batch_script") or die("Can not create the Slurm batch script");
        print SLURM "#!/bin/bash \n";
        print SLURM "#SBATCH-p $queue \n";
        print SLURM "#SBATCH-J $job_name \n";
        print SLURM "#SBATCH-N 1 \n";
        print SLURM "#SBATCH-c 1 \n";
        print SLURM "#SBATCH-o $filename_basepath.slurmout.log \n";
        print SLURM "#SBATCH-e $filename_basepath.slurmerr.log \n";
        print SLURM "#SBATCH-t $time_limit \n";
        print SLURM "source $gen_setup_script \n";
        print SLURM "cd $jobs_dir \n";
        print SLURM "$gmkspl_cmd | $grep_pipe &> $filename_basepath.mkspl.log \n";
        close(SLURM);
        `sbatch $batch_script`;
    } #Slurm

    # PBS case
    if($batch_system eq 'PBS') {
        $batch_script = "$filename_basepath.pbs";
        open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
        print PBS "#!/bin/bash \n";
        print PBS "#PBS -N $job_name \n";
        print PBS "#PBS -o $filename_basepath.pbsout.log \n";
        print PBS "#PBS -e $filename_basepath.pbserr.log \n";
        print PBS "source $gen_setup_script \n";
        print PBS "cd $jobs_dir \n";
        print PBS "$gmkspl_cmd | $grep_pipe &> $filename_basepath.mkspl.log \n";
        close(PBS);
        $job_submission_command = "qsub";
        `$job_submission_command -q $queue $batch_script`; 
    } #PBS 

    # LSF case
    if($batch_system eq 'LSF') {
        $batch_script = "$filename_basepath.sh";
        open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
        print LSF "#!/bin/bash \n";
        print PBS "#BSUB-j $job_name \n";
        print LSF "#BSUB-q $queue \n";
        print LSF "#BSUB-o $filename_basepath.lsfout.log \n";
        print LSF "#BSUB-e $filename_basepath.lsferr.log \n";
        print LSF "source $gen_setup_script \n";
        print LSF "cd $jobs_dir \n";
        print LSF "$gmkspl_cmd | $grep_pipe &> $filename_basepath.mkspl.log \n";
        close(LSF);
        `bsub < $batch_script`;
    } #LSF

    # no batch system, run jobs interactively
    if($batch_system eq 'none') {
       system("source $gen_setup_script $gen_version; cd $jobs_dir; $gmkspl_cmd");
    } # interactive mode


  }
}
