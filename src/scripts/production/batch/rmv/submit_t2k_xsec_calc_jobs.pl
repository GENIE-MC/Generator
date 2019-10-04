#!/usr/bin/perl

#----------------------------------------------------------------------------------------------------------------
# Submit jobs for calculating GENIE x-section splines for all the nuclear targets in T2K
#
# Syntax:
#   shell% perl submit_t2k_xsec_calc_jobs.pl <options>
#
# Options:
#    --version        : GENIE version number
#   [--arch]          : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]    : default: <version>
#   [--cycle]         : default: 01
#   [--use-valgrind]  : default: off
#   [--batch-system]  : <PBS, LSF, slurm, HTCondor, HTCondor_PBS, none>, default: PBS
#   [--queue]         : default: prod
#   [--softw-topdir]  : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/generator
#   [--jobs-topdir]   : top level dir for job files, default: /opt/ppd/t2k/softw/GENIE/scratch
#
# Author:
#   Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#----------------------------------------------------------------------------------------------------------------

use File::Path;

# inputs
#  
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')        { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--arch')           { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')     { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')          { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')   { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')   { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')          { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')   { $softw_topdir  = $ARGV[$iarg+1]; }
  if($_ eq '--jobs-topdir')    { $jobs_topdir   = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$use_valgrind   = 0                                    unless defined $use_valgrind;
$arch           = "SL6.x86_64"                         unless defined $arch;
$production     = "$genie_version"                     unless defined $production;
$cycle          = "01"                                 unless defined $cycle;
$batch_system   = "PBS"                                unless defined $batch_system;
$queue          = "prod"                               unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE/generator" unless defined $softw_topdir;
$jobs_topdir    = "/opt/ppd/t2k/softw/GENIE/scratch"   unless defined $jobs_topdir;  

$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$freenucsplines = "$softw_topdir/data/job_inputs/xspl/gxspl-vN-$genie_version.xml";
$jobs_dir       = "$jobs_topdir/xsec\_t2k-$production\_$cycle/";

$nkots     = 200;
$emax      =  35;
$neutrinos = "12,-12,14,-14";
%targets = (
	'H2'    =>  '1000010020',
	'B10'   =>  '1000050100',
	'B11'   =>  '1000050110',
	'C12'   =>  '1000060120',
	'C13'   =>  '1000060130', 
	'N14'   =>  '1000070140', 
	'N15'   =>  '1000070150', 
	'O16'   =>  '1000080160', 
	'O17'   =>  '1000080170', 
	'O18'   =>  '1000080180',
	'F19'   =>  '1000090190',
	'Na23'  =>  '1000110230',
	'Al27'  =>  '1000130270', 
	'Si28'  =>  '1000140280',
	'Si29'  =>  '1000140290',
	'Si30'  =>  '1000140300',
	'Cl35'  =>  '1000170350',
	'Cl37'  =>  '1000170370',
	'Ar36'  =>  '1000180360',
	'Ar38'  =>  '1000180380',
	'Ar40'  =>  '1000180400',
	'Ti46'  =>  '1000220460',
	'Ti47'  =>  '1000220470',
	'Ti48'  =>  '1000220480',
	'Ti49'  =>  '1000220490',
	'Ti50'  =>  '1000220500',
	'Fe54'  =>  '1000260540',
	'Fe56'  =>  '1000260560',
	'Fe57'  =>  '1000260570',
	'Fe58'  =>  '1000260580',
	'Co59'  =>  '1000270590',   
	'Cu63'  =>  '1000290630',
	'Cu65'  =>  '1000290650',
	'Zn64'  =>  '1000300640',
	'Zn66'  =>  '1000300660',
	'Zn67'  =>  '1000300670',
	'Zn68'  =>  '1000300680',
	'Zn70'  =>  '1000300700',
	'Pb204' =>  '1000822040',
	'Pb206' =>  '1000822060',
	'Pb207' =>  '1000822070',
	'Pb208' =>  '1000822080'  );


# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

#
# loop over nuclear targets & submit jobs
#
while( my ($tgt_name, $tgt_code) = each %targets ) {

    $filename_template = "$jobs_dir/job_$tgt_name";
    $grep_pipe     = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
    $cmd = "gmkspl -p $neutrinos -t $tgt_code -n $nkots -e $emax --input-cross-sections $freenucsplines --output-cross-sections gxspl_$tgt_name.xml | $grep_pipe  &> $filename_template.mkspl.log";
    print "@@ exec: $cmd \n";

    #
    # submit
    #

    # PBS case
    if($batch_system eq 'PBS' || $batch_system eq 'HTCondor_PBS') {
	$batch_script = "$filename_template.pbs";
	open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
	print PBS "#!/bin/bash \n";
        print PBS "#PBS -N $tgt_name \n";
        print PBS "#PBS -o $filename_template.pbsout.log \n";
        print PBS "#PBS -e $filename_template.pbserr.log \n";
	print PBS "source $genie_setup \n";
	print PBS "cd $jobs_dir \n";
	print PBS "$cmd \n";
        close(PBS);
        $job_submission_command = "qsub";
        if($batch_system eq 'HTCondor_PBS') {
           $job_submission_command = "condor_qsub"; 
        }
        `$job_submission_command -q $queue $batch_script`;
    } #PBS

    # LSF case
    if($batch_system eq 'LSF') {
	$batch_script = "$filename_template.sh";
	open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
	print LSF "#!/bin/bash \n";
        print LSF "#BSUB-j $tgt_name \n";
        print LSF "#BSUB-q $queue \n";
        print LSF "#BSUB-o $filename_template.pbsout.log \n";
        print LSF "#BSUB-e $filename_template.pbserr.log \n";
	print LSF "source $genie_setup \n";
	print LSF "cd $jobs_dir \n";
	print LSF "$cmd \n";
        close(LSF);
	`bsub < $batch_script`;
    } #LSF

}

