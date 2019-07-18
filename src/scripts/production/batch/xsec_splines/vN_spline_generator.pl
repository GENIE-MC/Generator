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
#   [--tune]          : Select a tune for Genie configuration
#   [--arch]          : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]    : default: routine_validation
#   [--cycle]         : default: 01
#   [--use-valgrind]  : default: off
#   [--batch-system]  : <PBS, LyonPBS, LSF, slurm, HTCondor, HTCondor_PBS, none>, default: PBS
#   [--queue]         : default: prod
#   [--softw-topdir]  : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]   : top level dir for job files, default: $PWD
#   [--gen-list]      : comma separated list of event generator list, default all
#   [--nu-list]       : comma separated list of neutrino. Both PDGs or names are ok, default all
#   [--e-max]         : maxmium energy of the splines in GeV. Default 200 GeV.
#   [--n-knots]       : number of knots per spline. Default 100.
#   [--with-priority] : (boolean) set a priority to optimize bulk productions. Default false
#   [--run-one]       : (boolean) If called, one of the jobs is run as part of the script instead of submitted
#                       via the batch system. Default all the jobs are submitted
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
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#-------------------------------------------------------------------------------------------------------------

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')        { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--config-dir')     { $config_dir    = $ARGV[$iarg+1]; }
  if($_ eq '--tune')           { $tune          = $ARGV[$iarg+1]; }
  if($_ eq '--arch')           { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')     { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')          { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')   { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')   { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')          { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')   { $softw_topdir  = $ARGV[$iarg+1]; }
  if($_ eq '--jobs-topdir')    { $jobs_topdir   = $ARGV[$iarg+1]; }
  if($_ eq '--gen-list'   )    { $req_gen_list  = $ARGV[$iarg+1]; }
  if($_ eq '--nu-list' )       { $req_nu_list	  = $ARGV[$iarg+1]; }
  if($_ eq '--e-max' )         { $e_max	        = $ARGV[$iarg+1]; }
  if($_ eq '--n-knots' )       { $n_knots    	  = $ARGV[$iarg+1]; }
  if($_ eq '--with-priority' ) { $priority	    = 1; }
  if($_ eq '--run-one' )       { $run_one	      = 1; }
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
$priority       = 0                             unless defined $priority ;
$e_max          = 200                           unless defined $e_max ;
$n_knots        = 100                           unless defined $n_knots ;


%nucleons_pdg = ( 'n'  =>  1000000010,
                  'p'  =>  1000010010 );

%nucleons_name = ( 1000000010 => 'n' ,
                   1000010010 => 'p' );

@nucleons_proc = ( 'none',
                   'CCRES', 'NCRES',
                   'CCDIS', 'NCDIS',
                   'CCDFR', 'NCDFR',
                   'Fast' );

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
    if ( exists $nu_pdg_def{$nu} ) { push @nu_list, $nu ; }
    if ( exists	$nu_name_def{$nu} ) { push @nu_list, $nu_name_def{$nu} ; }
  }
}
else {
  @nu_list = values %nu_name_def;
}
print "@nu_list \n";


if ( defined $req_gen_list ) {
  my @proc_temp_list = split( ",", $req_gen_list );
  @nucleons_proc_list = ();
  foreach my $proc ( @proc_temp_list ) {
    if ( grep  {$_ eq $proc} @nucleons_proc ) { push @nucleons_proc_list, $proc ; }
  }
}
else {
  @nucleons_proc_list = @nucleons_proc;
}
print "Process List: @nucleons_proc_list \n";

#
# make the jobs directory
#

@batch_commands = ();
@direct_commands = () ;

print "@@ Creating job directory: $jobs_dir \n";
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

foreach $nu ( @nu_list ) {
  foreach $tgt ( keys %nucleons_pdg ) {
    foreach $proc ( @nucleons_proc_list ) {

      if ( $proc eq "none" ) {
        next ;
      }

      if ( $proc eq 'CCQE' ) {
        if ( ($tgt eq 'n' ) && ( $nu_pdg_def{$nu} < 0 ) ) { next ; }
        if ( ($tgt eq 'p' ) && ( $nu_pdg_def{$nu} > 0 ) ) { next ; }
      }

      if ( ( $proc eq 'CCDFR' ) || ( $proc eq 'NCDFR' ) ) {
          if ( $tgt eq 'n' ) { next ; }
      }

      if ( $proc eq 'Fast' ) {
          $event_gen_list = 'FastOn' . ( uc $tgt );
      }
      else {
      	  $event_gen_list = $proc ;
      }

      $jobname = $nu."_on_".$tgt."_$proc";
      $filename_template = "$jobs_dir/$jobname";

      $grep_pipe     = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
      $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
      $gmkspl_opt    = "-p $nu_pdg_def{$nu} -t $nucleons_pdg{$tgt} -n $n_knots -e $e_max -o $filename_template.xml --event-generator-list $event_gen_list ";
      if ( defined $tune ) {
	  $gmkspl_opt.= " --tune $tune ";
      }
      $gmkspl_cmd    = "gmkspl $gmkspl_opt";

      print "@@ exec: $gmkspl_cmd \n";

      push( @direct_commands, "source $genie_setup $config_dir; cd $jobs_dir; $gmkspl_cmd" ) ;

      # PBS case
      if($batch_system eq 'PBS' || $batch_system eq 'HTCondor_PBS') {
         $batch_script = "$filename_template.pbs";
         open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
         print PBS "#!/bin/bash \n";
         print PBS "#PBS -N $jobname \n";
         print PBS "#PBS -o $filename_template.pbsout.log \n";
         print PBS "#PBS -e $filename_template.pbserr.log \n";
	 print PBS "#PBS -p -1 \n" if ( $priority ) ;
         print PBS "source $genie_setup $config_dir \n";
         print PBS "cd $jobs_dir \n";
         print PBS "$gmkspl_cmd \n";
         close(PBS);
         $job_submission_command = "qsub";
         if($batch_system eq 'HTCondor_PBS') {
            $job_submission_command = "condor_qsub";
         }

	 push( @batch_commands, "$job_submission_command -q $queue $batch_script" ) ;
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
         print PBS "#\$ -l ct=8:00:00,sps=1 \n";
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
 	 print PBS "#BSUB-j $jobname \n";
 	 print LSF "#BSUB-q $queue \n";
 	 print LSF "#BSUB-o $filename_template.lsfout.log \n";
 	 print LSF "#BSUB-e $filename_template.lsferr.log \n";
 	 print LSF "source $genie_setup $config_dir \n";
 	 print LSF "cd $jobs_dir \n";
 	 print LSF "$gmkspl_cmd | $grep_pipe &> $filename_template.mkspl.log \n";
 	 close(LSF);

	 push( @batch_commands, "bsub < $batch_script " ) ;

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
 	 push ( @batch_commands, "condor_submit $batch_script" ) ;
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

	 push( @batch_commands, "sbatch --job-name=$jobname $batch_script" ) ;

      } #slurm

    }
  }
}


#
# submit
#




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
