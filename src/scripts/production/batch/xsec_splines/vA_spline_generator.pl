#!/usr/bin/perl

#----------------------------------------------------------------------------------------------------------------
# Submit jobs for calculating GENIE x-section splines for all nuclear targets and at the energy range required
# for generating the GENIE release validation samples. Note that other scripts are available for generating
# x-section splines for the larger array of targets found in detector geometry descriptions of given expts.
# If needed, use the GENIE gspladd utility to merge the job outputs.
#
# Syntax:
#   shell% perl submit_vA_xsec_calc_jobs.pl <options>
#
# Options:
#    --version         : GENIE version number
#   [--config-dir]     : Config directory, default is $GENIE/config
#   [--tune]           : tune to be used for genie configuration
#   [--arch]           : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]     : default: routine_validation
#   [--cycle]          : default: 01
#   [--use-valgrind]   : default: off
#   [--batch-system]   : <PBS, LyonPBS, LSF, slurm, HTCondor, HTCondor_PBS, none>, default: PBS
#   [--queue]          : default: prod
#   [--softw-topdir]   : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]    : top level dir for job files, default: $PWD
#   [--freenucsplines] : Absolute path to free nucleon splines, default: $softw_topdir/data/job_inputs/xspl/gxspl-vN-$genie_version.xml
#   [--free-nuc-dir]   : Absolute path to free nuclear spline directory
#   [--gen-list]       : comma separated list of event generator list, default all
#   [--target-list]    : comma separated list of targets' PDGs, default De,He4,C12,O16,Ar40,Fe56,Pb207.
#                        Note that it needs the PDG code, not chemical name.
#   [--nu-list]        : comma separeted list of neutrino flavors. Both PDGs or names like vmubar,ve,vtau. Default all
#   [--e-max]         : maxmium energy of the splines in GeV. Default 200 GeV.
#   [--n-knots]       : number of knots per spline. Default 100.
#   [--with-priority]  : (boolean) set a prioirty to optimize bulk productions. Default false
#   [--run-one]        : (boolean) If called, one of the jobs is run as part of the script instead of submitted
#                        via the batch system. Default all the jobs are submitted
#
#
# Author:
#   Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#----------------------------------------------------------------------------------------------------------------

use File::Path;
use File::Basename;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')        { $genie_version  = $ARGV[$iarg+1]; }
  if($_ eq '--config-dir')     { $config_dir     = $ARGV[$iarg+1]; }
  if($_ eq '--tune')           { $tune           = $ARGV[$iarg+1]; }
  if($_ eq '--arch')           { $arch           = $ARGV[$iarg+1]; }
  if($_ eq '--production')     { $production     = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')          { $cycle          = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')   { $use_valgrind   = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')   { $batch_system   = $ARGV[$iarg+1]; }
  if($_ eq '--queue')          { $queue          = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')   { $softw_topdir   = $ARGV[$iarg+1]; }
  if($_ eq '--jobs-topdir')    { $jobs_topdir    = $ARGV[$iarg+1]; }
  if($_ eq '--freenucsplines') { $freenucsplines = $ARGV[$iarg+1]; }
  if($_ eq '--free-nuc-dir' )  { $free_nuc_dir   = $ARGV[$iarg+1]; }
  if($_ eq '--gen-list'   )    { $req_gen_list   = $ARGV[$iarg+1]; }
  if($_ eq '--target-list'   ) { $target_list	   = $ARGV[$iarg+1]; }
  if($_ eq '--nu-list'   )     { $req_nu_list    = $ARGV[$iarg+1]; }
  if($_ eq '--e-max' )         { $e_max	        = $ARGV[$iarg+1]; }
  if($_ eq '--n-knots' )       { $n_knots    	  = $ARGV[$iarg+1]; }
  if($_ eq '--with-priority' ) { $priority       = 1; }
  if($_ eq '--run-one' )       { $run_one        = 1; }
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
$freenucsplines = "$softw_topdir/data/job_inputs/xspl/gxspl-vN-$genie_version.xml" unless defined $freenucsplines;
$priority       = 0                             unless defined $priority ;
$e_max          = 200                           unless defined $e_max ;
$n_knots        = 100                           unless defined $n_knots ;


$genie_setup    = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$jobs_dir       = "$jobs_topdir/$genie_version-$production\_$cycle-xsec\_vA";

%nu_pdg_def = ( 've'      =>   12,
                'vebar'   =>  -12,
                'vmu'     =>   14,
                'vmubar'  =>  -14,
                'vtau'    =>   16,
                'vtaubar' =>  -16 );

%nu_name_def = ( 12 => 've'     ,
                -12 => 'vebar'  ,
                 14 => 'vmu'    ,
                -14 => 'vmubar' ,
                 16 => 'vtau'   ,
                -16 => 'vtaubar' );

##create the list of neutrino to be produced
if ( defined $req_nu_list ) {
  my @nu_temp_list = split( ",", $req_nu_list );
  @nu_list = ();
  foreach my $nu ( @nu_temp_list ) {
    if ( exists $nu_pdg_def{$nu} )   { push @nu_list, $nu ; }
    if ( exists $nu_name_def{$nu} )  { push @nu_list, $nu_name_def{$nu} ; }
  }
}
else {
  @nu_list = values %nu_name_def;
}
print "\n Neutrino List: @nu_list \n";


@nuclei_proc = ( 'none',
                 'WeakMEC',
                 'CCCOHPION',    'NCCOHPION',
                 'Fast',
                 'FastWithMEC'  ## this is supposed to be better with G16_01 comprehensive models
                 );


@nuclei_proc_def = ( 'none',
                     'WeakMEC',
                     'CCCOHPION',    'NCCOHPION',
                     'Fast'
	             );

# create the lsit of processes to be generated for composite nuclei
if ( defined $req_gen_list ) {
  my @temp_list = split( ",", $req_gen_list );
  @nuclei_proc_list = ();
  foreach my $proc ( @temp_list ) {
    if ( grep  {$_ eq $proc} @nuclei_proc ) { push @nuclei_proc_list, $proc ; }
  }
}
else {
  @nuclei_proc_list = @nuclei_proc_def ;
}
print "\n List of processes on composite nuclei: @nuclei_proc_list \n";

if ( defined $target_list ) {
@tgt_pdg = split( ",", $target_list );
}
else {
@tgt_pdg = ( 1000060120, #C12
             1000100200, #Ne20
             1000130270, #Al27
             1000140300, #Si30
	     1000180400, #Ar40
             1000260560, #Fe56
             1000822070, #Pb207
           );
}
print "\n Target List: @tgt_pdg \n";


#
# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

@batch_commands = ();
@direct_commands = () ;

foreach $nu ( @nu_list ) {
  foreach $tgt ( @tgt_pdg ) {
    foreach $proc ( @nuclei_proc_list ) {

      if ( $proc eq "none" ) {
        next ;
      }

      if ( $proc eq 'Fast' ) {
          $event_gen_list = 'FastOnNuclei' ;
      }
      else {
          $event_gen_list = $proc ;
      }


      $jobname  = $nu."_on_".$tgt."_$proc";
      $filename_template = "$jobs_dir/$jobname";
      $grep_pipe  = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
      if ( defined $free_nuc_dir ) {
        $in_splines=$free_nuc_dir."/"."total_xsec.xml";
      }
      else {
        $in_splines = $freenucsplines;
      }
      $gmkspl_opt = "-p $nu_pdg_def{$nu} -t $tgt -n $n_knots -e $e_max --input-cross-sections $in_splines --output-cross-sections $filename_template.xml --event-generator-list $event_gen_list --no-copy";
      $gmkspl_opt .= " --tune $tune " if ( defined $tune ) ;
      $gmkspl_cmd = "gmkspl $gmkspl_opt ";
      print "@@ exec: $gmkspl_cmd \n";

      #
      # submit
      #

      push( @direct_commands, "source $genie_setup $config_dir; cd $jobs_dir; $gmkspl_cmd" ) ;

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
	print PBS "$gmkspl_cmd \n";
        close(PBS);
        $job_submission_command = "qsub";
        if($batch_system eq 'HTCondor_PBS') {
           $job_submission_command = "condor_qsub";
        }

	push( @batch_commands, "$job_submission_command -q $queue $batch_script" ) ;

      } #PBS

      # LyonPBS case
      if($batch_system eq 'LyonPBS' ) {
         $batch_script = "$filename_template.pbs";
         open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
         print PBS "#!/bin/bash \n";
         print PBS "#\$ -P P_$ENV{'GROUP'} \n";
         print PBS "#\$ -N $jobname \n";
         print PBS "#\$ -o $filename_template.pbsout.log \n";
         print PBS "#\$ -e $filename_template.pbserr.log \n";
         print PBS "#\$ -l ct=30:00:00,sps=1,s_rss=4G \n";
         print PBS "#\$ -p -1 \n" if ( $priority ) ;
         print PBS "source $genie_setup $config_dir \n";
         print PBS "cd $jobs_dir \n";
         print PBS "$gmkspl_cmd \n";
         close(PBS);
         $job_submission_command = "qsub";

	 push( @batch_commands, "$job_submission_command  $batch_script " ) ;

       } #LyonPBS

      # LSF case
      if($batch_system eq 'LSF') {
	$batch_script = "$filename_template.sh";
	open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
	print LSF "#!/bin/bash \n";
        print LSF "#BSUB-j $jobname \n";
        print LSF "#BSUB-o $filename_template.lsfout.log \n";
        print LSF "#BSUB-e $filename_template.lsferr.log \n";
	print LSF "source $genie_setup $config_dir \n";
	print LSF "cd $jobs_dir \n";
	print LSF "$gmkspl_cmd \n";
        close(LSF);

	push( @batch_commands, "bsub < $batch_script" ) ;
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

	push( @batch_commands, "condor_submit $batch_script" ) ;

      } #HTCondor

      # slurm case
      if($batch_system eq 'slurm') {
	$batch_script = "$filename_template.sh";
	open(SLURM, ">$batch_script") or die("Can not create the SLURM batch script");
	print SLURM "#!/bin/bash \n";
        print SLURM "#SBATCH-p $queue \n";
        print SLURM "#SBATCH-o $filename_template.lsfout.log \n";
        print SLURM "#SBATCH-e $filename_template.lsferr.log \n";
	print SLURM "source $genie_setup $config_dir \n";
	print SLURM "cd $jobs_dir \n";
	print SLURM "$gmkspl_cmd \n";
        close(SLURM);

	push( @batch_commands, "sbatch --job-name=$jobname $batch_script" );
      } #slurm

      # run interactively
    }
  }
}


if ( $batch_system eq 'none' ) {
    ## run all of them interactively

    for my $run_cmd ( @direct_commands ) {
        `$run_cmd` ;
    }
}
else {
    ## submit all except the first
    foreach my $i ( 1 .. $#batch_commands ) {
        `$batch_commands[$i]` ;
    }

    # handle the first according to script options
    if ( defined $run_one ) {
        `$direct_commands[0]` ;
    }
    else {
        `$batch_commands[0]` ;
    }

}
