#!/usr/bin/perl

#----------------------------------------------------------------------
# Submit jobs for calculating GENIE eA cross section splines to be used 
# with GENIE's validation programs comparing GENIE against electron 
# scattering data.
#
# Syntax:
#   shell% perl submit_eA_xsec_calc_jobs.pl <options>
#
# Options:
#    --version       : GENIE version number
#   [--arch]         : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]   : default: <version>
#   [--cycle]        : default: 01
#   [--use-valgrind] : default: off
#   [--batch-system] : <PBS, LSF, slurm, none>, default: PBS
#   [--queue]        : default: prod
#   [--softw-topdir] : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]  : top level dir for job files, default: /opt/ppd/t2k/scratch/GENIE/
#
# Author:
#   Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#----------------------------------------------------------------------

use File::Path;

# inputs
#  
$iarg=0;
foreach (@ARGV) {
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

$use_valgrind   = 0                             unless defined $use_valgrind;
$arch           = "SL6.x86_64"                  unless defined $arch;
$production     = "$genie_version"              unless defined $production;
$cycle          = "01"                          unless defined $cycle;
$batch_system   = "PBS"                         unless defined $batch_system;
$queue          = "prod"                        unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE"    unless defined $softw_topdir;
$jobs_topdir    = "/opt/ppd/t2k/scratch/GENIE/" unless defined $jobs_topdir;
$genie_setup    = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$jobs_dir       = "$jobs_topdir/xsec\_eA-$production\_$cycle/";
$time_limit     = "60:00:00";

$nkots     = 200;
$emax      =  35;
$probes    = "11";
%targets = (
	'H1'    =>  '1000010020',
	'H2'    =>  '1000010020',
	'H3'    =>  '1000010030',
	'He3'   =>  '1000020030',
	'He4'   =>  '1000020040',
	'C12'   =>  '1000060120',
	'N14'   =>  '1000070140', 
	'O16'   =>  '1000080160', 
	'Ne20'  =>  '1000100200', 
	'Al27'  =>  '1000130270', 
	'Ca40'  =>  '1000200400', 
	'Ca48'  =>  '1000200480', 
	'Fe56'  =>  '1000260560',
	'Kr83'  =>  '1000360830',
	'Xe131' =>  '1000541310',  
	'Au197' =>  '1000791970',  
	'Pb208' =>  '1000822080',  
	'U238'  =>  '1000922380'   );

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

#
# loop over nuclear targets & submit jobs
#
while( my ($tgt_name, $tgt_code) = each %targets ) {

    $fntemplate  = "$jobs_dir/job_$tgt_name";
    $grep_pipe   = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
    $gmkspl_opt  = "-p $probes -t $tgt_code -n $nkots -e $emax -o gxspl_emode_$tgt_name.xml --event-generator-list EM";
    $gmkspl_cmd  = "gmkspl $gmkspl_opt | $grep_pipe &> $fntemplate.mkspl.log";
    print "@@ exec: $gmkspl_cmd \n";

    # PBS case
    if($batch_system eq 'PBS') {
	$batch_script = "$fntemplate.pbs";
	open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
	print PBS "#!/bin/bash \n";
        print PBS "#PBS -N $tgt_name \n";
        print PBS "#PBS -o $fntemplate.pbsout.log \n";
        print PBS "#PBS -e $fntemplate.pbserr.log \n";
	print PBS "source $genie_setup \n";
	print PBS "cd $jobs_dir \n";
	print PBS "$gmkspl_cmd \n";
        close(PBS);
	`qsub -q $queue $batch_script`;
    } #PBS

    # LSF case
    if($batch_system eq 'LSF') {
	$batch_script = "$fntemplate.sh";
	open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
	print LSF "#!/bin/bash \n";
        print LSF "#BSUB-j $tgt_name \n";
        print LSF "#BSUB-q $queue \n";
        print LSF "#BSUB-o $fntemplate.lsfout.log \n";
        print LSF "#BSUB-e $fntemplate.lsferr.log \n";
	print LSF "source $genie_setup \n";
	print LSF "cd $jobs_dir \n";
	print LSF "$gmkspl_cmd \n";
        close(LSF);
	`bsub < $batch_script`;
    } #LSF

    # slurm case
    if($batch_system eq 'slurm') {
     my $time_lim = `sinfo -h -p batch -o %l`;
     my ($days, $hours, $remainder) = $time_lim =~ /([0]+)-([0-9]+):(.*)/;
     my $newhours = $days * 24 + $hours;
     my $new_time_lim = "$newhours:$remainder";
     $time_limit = $new_time_lim lt $time_limit ? $new_time_lim : $time_limit;
	$batch_script = "$fntemplate.sh";
	open(SLURM, ">$batch_script") or die("Can not create the SLURM batch script");
	print SLURM "#!/bin/bash \n";
        print SLURM "#SBATCH-p $queue \n";
        print SLURM "#SBATCH-o $fntemplate.lsfout.log \n";
        print SLURM "#SBATCH-e $fntemplate.lsferr.log \n";
        print SLURM "#SBATCH-t $time_limit \n";
	print SLURM "source $genie_setup \n";
	print SLURM "cd $jobs_dir \n";
	print SLURM "$gmkspl_cmd \n";
        close(SLURM);
	`sbatch --job-name=$tgt_name $batch_script`;
    } #slurm

    # no batch system, run jobs interactively
    if($batch_system eq 'none') {
        system("source $genie_setup; cd $jobs_dir; $gmkspl_cmd");
    } # interactive mode

}

