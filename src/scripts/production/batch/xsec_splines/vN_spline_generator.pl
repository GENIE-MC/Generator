#!/usr/bin/perl

#----------------------------------------------------------------------------------------------------------------
# Submit jobs for calculating GENIE free-nucleon cross-section splines which can then be re-used for calculating 
# nuclear cross-sections. If needed, use the GENIE gspladd utility to merge the job outputs.
#
# Syntax:
#   shell% perl submit_vN_xsec_calc_jobs.pl <options>
#
# Options:
#    --version        : genie version number
#   [--config-dir]    : Config directory, default is $GENIE/config 
#   [--arch]          : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]    : default: routine_validation
#   [--cycle]         : default: 01
#   [--use-valgrind]  : default: off
#   [--batch-system]  : <PBS, LyonPBS, LSF, slurm, HTCondor, HTCondor_PBS, none>, default: PBS
#   [--queue]         : default: prod
#   [--softw-topdir]  : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]   : top level dir for job files, default: $PWD
#   [--gen-list]      : comma separated list of event generator list, default all
#
# Examples:
#   shell% perl submit_vN_xsec_calc_jobs.pl --version v2.6.0  
#   shell% perl submit_vN_xsec_calc_jobs.pl --version v2.6.0 
#   shell% perl submit_vN_xsec_calc_jobs.pl --version v2.6.0  
#
# Author:
#   Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2016, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#-------------------------------------------------------------------------------------------------------------

use File::Path;

# inputs
#  
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')        { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--config-dir')     { $config_dir    = $ARGV[$iarg+1]; }
  if($_ eq '--arch')           { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')     { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')          { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')   { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')   { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')          { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')   { $softw_topdir  = $ARGV[$iarg+1]; }           
  if($_ eq '--jobs-topdir')    { $jobs_topdir   = $ARGV[$iarg+1]; }           
  if($_ eq '--gen-list'   )    { $gen_list      = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$config_dir     = ""                            unless defined $config_dir;
$use_valgrind   = 0                             unless defined $use_valgrind;
$arch           = "SL6.x86_64"                  unless defined $arch;
$production     = "routine_validation"          unless defined $production;
$cycle          = "01"                          unless defined $cycle;
$batch_system   = "PBS"                         unless defined $batch_system;
$queue          = "prod"                        unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE/"   unless defined $softw_topdir;
$jobs_topdir    = $ENV{'PWD'}                   unless defined $jobs_topdir;
$genie_setup    = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$jobs_dir       = "$jobs_topdir/$genie_version-$production\_$cycle-xsec\_vN";


$nkots = 500;
$emax  = 500;

%nu_pdg = ( 've'      =>   12, 
            'vebar'   =>  -12, 
            'vmu'     =>   14, 
            'vmubar'  =>  -14, 
            'vtau'    =>   16, 
            'vtaubar' =>  -16 );

%tgt_pdg = ( 'n'  =>  1000000010, 
             'p'  =>  1000010010 );

if ( defined $gen_list ) {
###print $gen_list."\n";
@proc_v = split( ",", $gen_list ); 
}
else {
@proc_v = ( 'CCQE',     'NCEL', 
            'CCRES',    'NCRES', 
            'CCDIS',    'NCDIS', 
            'GLRES', 
            'CCDFR',    'NCDFR', 
            'CharmCCDIS', 'CharmCCQE', 
            'LambdaCCQE', 'SingleKaon', 
            'NuEElastic' );
}
 
# make the jobs directory
#
print "@@ Creating job directory: $jobs_dir \n";
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

foreach $nu ( keys %nu_pdg ) { 
  foreach $tgt ( keys %tgt_pdg ) { 
    foreach $proc ( @proc_v ) {

      $jobname = $nu."_on_".$tgt."_$proc"; 
      $filename_template = "$jobs_dir/$jobname"; 
      

      $grep_pipe     = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
      $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
      $gmkspl_opt    = "-p $nu_pdg{$nu} -t $tgt_pdg{$tgt} -n $nkots -e $emax -o $filename_template.xml --event-generator-list $proc";
      $gmkspl_cmd    = "gmkspl $gmkspl_opt";

      print "@@ exec: $gmkspl_cmd \n";

      # PBS case
      if($batch_system eq 'PBS' || $batch_system eq 'HTCondor_PBS') {
         $batch_script = "$filename_template.pbs";
         open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
         print PBS "#!/bin/bash \n";
         print PBS "#PBS -N $jobname \n";
         print PBS "#PBS -o $filename_template.pbsout.log \n";
         print PBS "#PBS -e $filename_template.pbserr.log \n";
         print PBS "source $genie_setup $config_dir \n";
         print PBS "cd $jobs_dir \n";
         print PBS "$gmkspl_cmd | $grep_pipe &> $filename_template.mkspl.log \n";
         close(PBS);
         $job_submission_command = "qsub";
         if($batch_system eq 'HTCondor_PBS') {
            $job_submission_command = "condor_qsub";
         }
         `$job_submission_command -q $queue $batch_script`;
       } #PBS / #HTCondor_PBS


      # LyonPBS case
      if($batch_system eq 'LyonPBS' ) {
         $batch_script = "$filename_template.pbs";
         open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
         print PBS "#!/bin/bash \n";
         print PBS "#\$ -P P_$ENV{'GROUP'} \n";
         print PBS "#\$ -N $jobname \n";
         print PBS "#\$ -o $filename_template.pbsout.log \n";
         print PBS "#\$ -e $filename_template.pbserr.log \n";
         print PBS "#\$ -l ct=6:00:00,sps=1 \n";
         print PBS "source $genie_setup $config_dir \n";
         print PBS "cd $jobs_dir \n";
         print PBS "$gmkspl_cmd \n";
         close(PBS);
         $job_submission_command = "qsub";
         `$job_submission_command  $batch_script`;
       } #LyonPBS 

       # LSF case
       if($batch_system eq 'LSF') {
    	 $batch_script = "$filename_template.sh";
  	 open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
 	 print LSF "#!/bin/bash \n";
 	 print PBS "#BSUB-j $jobname \n";
 	 print LSF "#BSUB-q $queue \n";
 	 print LSF "#BSUB-o $filename_template.lsfout.log \n";
 	 print LSF "#BSUB-e $filename_template.lsferr.log \n";
 	 print LSF "source $genie_setup $config_dir \n";
 	 print LSF "cd $jobs_dir \n";
 	 print LSF "$gmkspl_cmd | $grep_pipe &> $filename_template.mkspl.log \n";
 	 close(LSF);
 	 `bsub < $batch_script`;
      } #LSF

      # HTCondor
      if($batch_system eq 'HTCondor') {
	 $batch_script = "$filename_template.htc";
	 open(HTC, ">$batch_script") or die("Can not create the Condor submit description file: $batch_script");
	 print HTC "Universe               = vanilla \n";
	 print HTC "Executable             = $softw_topdir/generator/builds/$arch/$genie_version/src/scripts/production/batch/htcondor_exec.sh \n";
	 print HTC "Arguments              = $genie_setup $jobs_dir $gmkspl_cmd \n";
 	 print HTC "Log                    = $filename_template.log \n";
         print HTC "Output                 = $filename_template.out \n";
 	 print HTC "Error                  = $filename_template.err \n";
 	 print HTC "Request_memory         = 2 GB \n";
 	 print HTC "Queue \n";
 	 close(HTC);
 	 `condor_submit $batch_script`;
      } #HTCondor

      # slurm case
      if($batch_system eq 'slurm') {
 	 $batch_script = "$filename_template.sh";
 	 open(SLURM, ">$batch_script") or die("Can not create the slurm batch script");
 	 print SLURM "#!/bin/bash \n";
 	 print SLURM "#SBATCH-p $queue \n";
 	 print SLURM "#SBATCH-o $filename_template.slurmout.log \n";
 	 print SLURM "#SBATCH-e $filename_template.slurmerr.log \n";
         print SLURM "source $genie_setup $config_dir \n";
 	 print SLURM "cd $jobs_dir \n";
 	 print SLURM "$gmkspl_cmd | $grep_pipe &> $filename_template.mkspl.log \n";
 	 close(SLURM);
 	 `sbatch --job-name=$jobname $batch_script`;
      } #slurm

      # no batch system, run jobs interactively
      if($batch_system eq 'none') {
	system("source $genie_setup $config_dir; cd $jobs_dir; $gmkspl_cmd");
      } # interactive mode

    } 
  }
} 




    #
    # submit
    #





