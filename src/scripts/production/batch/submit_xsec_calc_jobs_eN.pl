#!/usr/bin/perl

#----------------------------------------------------------------------------------------------------------------
# Submit jobs for calculating GENIE free-nucleon cross-sections used for constructing cross-section splines. 
# Generate data for all or specified processes, and split the task amongst several worker nodes.                    
# The GENIE gspladd utility can be used to merge the job outputs.
#
# Syntax:
#   shell% perl submit_xsec_calc_jobs_eN.pl <options>
#
# Options:
#    --gen-version     : GENIE generator version number
#    --tune            : GENIE physics tune
#    --xsec-spline-set : set of cross-section splines to generate
#   [--arch]           : <el7.x86_64, ...>, default: el7.x86_64
#   [--production]     : default: routine_validation
#   [--cycle]          : default: 01
#   [--use-valgrind]   : default: off
#   [--batch-system]   : <Slurm, PBS, LSF, none>, default: Slurm
#   [--queue]          : default: compute
#   [--softw-topdir]   : top level dir for softw installations, default: /user/costasa/projects/GENIE/softw/
#   [--jobs-topdir]    : top level dir for job files, default: /scratch/costasa/GENIE/
#
# Examples:
#   % perl submit_xsec_calc_jobs_eN.pl --gen-version v3.06.00 --tune G18_10a_02_11b --xsec-spline-set dis_eminusp,dis_eminusn
#   % perl submit_xsec_calc_jobs_eN.pl --gen-version v3.06.00 --tune G18_10a_02_11b --xsec-spline-set all 
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
$production       = "routine_validation"                    unless defined $production;
$cycle            = "01"                                    unless defined $cycle;
$batch_system     = "Slurm"                                 unless defined $batch_system;
$queue            = "compute"                               unless defined $queue;
$time_limit       = "10:00:00"                              unless defined $time_limit;  
$softw_topdir     = "/user/costasa/projects/GENIE/softw/"   unless defined $softw_topdir;
$jobs_topdir      = "/scratch/costasa/GENIE/"               unless defined $jobs_topdir;
$gen_setup_script = "$softw_topdir/generator/builds/$arch/$gen-version-setup.sh";
$jobs_dir         = "$jobs_topdir/$gen_version-$tune-$production\_$cycle-xsec\_eN/";
$time_limit       = "10:00:00";

$nkots = 100;
$emax  = 100;

%EPDG = ( 
    'qel_eminusp'    =>  '11' ,
    'qel_eminusn'    =>  '11' ,
    'qel_eplusp'     => '-11' ,
    'qel_eplusn'     => '-11' ,
    'res_eminusp'    =>  '11' ,
    'res_eminusn'    =>  '11' ,
    'res_eplusp'     => '-11' ,
    'res_eplusn'     => '-11' ,
    'dis_eminusp'    =>  '11' ,
    'dis_eminusn'    =>  '11' ,
    'dis_eplusp'     => '-11' ,
    'dis_eplusn'     => '-11'
);

%TGTPDG = (
    'qel_eminusp'    =>  '1000010010' ,
    'qel_eminusn'    =>  '1000000010' ,
    'qel_eplusp'     =>  '1000010010' ,
    'qel_eplusn'     =>  '1000000010' ,
    'res_eminusp'    =>  '1000010010' ,
    'res_eminusn'    =>  '1000000010' ,
    'res_eplusp'     =>  '1000010010' ,
    'res_eplusn'     =>  '1000000010' ,
    'dis_eminusp'    =>  '1000010010' ,
    'dis_eminusn'    =>  '1000000010' ,
    'dis_eplusp'     =>  '1000010010' ,
    'dis_eplusn'     =>  '1000000010'
);

%GEVGL = ( 
    'qel_eminusp'    =>  'EMQE'  ,
    'qel_eminusn'    =>  'EMQE'  ,
    'qel_eplusp'     =>  'EMQE'  ,
    'qel_eplusn'     =>  'EMQE'  ,
    'res_eminusp'    =>  'EMRES' ,
    'res_eminusn'    =>  'EMRES' ,
    'res_eplusp'     =>  'EMRES' ,
    'res_eplusn'     =>  'EMRES' ,
    'dis_eminusp'    =>  'EMDIS' ,
    'dis_eminusn'    =>  'EMDIS' ,
    'dis_eplusp'     =>  'EMDIS' ,
    'dis_eplusn'     =>  'EMDIS'
);

%OUTXML = ( 
    'qel_eminusp'    =>  'partial-xspl-e-p-qel.xml' ,
    'qel_eminusn'    =>  'partial-xspl-e-n-qel.xml' ,
    'qel_eplusp'     =>  'partial-xspl-e+p-qel.xml' ,
    'qel_eplusn'     =>  'partial-xspl-e+n-qel.xml' ,
    'res_eminusp'    =>  'partial-xspl-e-p-res.xml' ,
    'res_eminusn'    =>  'partial-xspl-e-n-res.xml' ,
    'res_eplusp'     =>  'partial-xspl-e+p-res.xml' ,
    'res_eplusn'     =>  'partial-xspl-e+n-res.xml' ,
    'dis_eminusp'    =>  'partial-xspl-e-p-dis.xml' ,
    'dis_eminusn'    =>  'partial-xspl-e-n-dis.xml' ,
    'dis_eplusp'     =>  'partial-xspl-e+p-dis.xml' ,
    'dis_eplusn'     =>  'partial-xspl-e+n-dis.xml'
);

# make the jobs directory
#
print "@@ Creating job directory: $jobs_dir \n";
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

for my $curr_xsplset (keys %OUTXML)  {
  if($xsplset=~m/$curr_xsplset/ || $xsplset eq "all") {

    # Get runnu-dependent info 
    # -------------------------------------------------------------

    $probe  = $EPDG    {$curr_xsplset};
    $tgt    = $TGTPDG  {$curr_xsplset};
    $gevgl  = $GEVGL   {$curr_xsplset};
    $outxml = $OUTXML  {$curr_xsplset};

    $jobname = "eNxscalc-$curr_xsplset"; 
    $filename_basepath = "$jobs_dir/$jobname"; 

    $grep_pipe     = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
    $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
    $gmkspl_opt    = "-p $probe -t $tgt -n $nkots -e $emax -o $outxml --event-generator-list $gevgl --tune $tune";
    $gmkspl_cmd    = "gmkspl $gmkspl_opt";

    print "@@ exec: $gmkspl_cmd \n";

    # Submit jobs
    # -------------------------------------------------------------
  
    # Slurm case
    if($batch_system eq 'Slurm') {
        $batch_script = "$filename_basepath.sh";
        open(SLURM, ">$batch_script") or die("Can not create the Slurm batch script");
        print SLURM "#!/bin/bash \n";
        print SLURM "#SBATCH-p $queue \n";
        print SLURM "#SBATCH-J $jobname \n";
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
        print PBS "#PBS -N $jobname \n";
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
        print PBS "#BSUB-j $jobname \n";
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
