#!/usr/bin/perl

#----------------------------------------------------------------------------------------------------------------
# Submit jobs for calculating GENIE charged lepton - free-nucleon cross-section splines which can then be re-used for calculating
# nuclear cross-sections. If needed, use the GENIE gspladd utility to merge the job outputs.
#
# Syntax:
#   shell% perl submit_eN_xsec_calc_jobs.pl <options>
#
# Options:
#    --version        : genie version number
#   [--config-dir]    : Config directory, default is $GENIE/config
#   [--tune]          : Select a tune for Genie configuration
#   [--arch]          : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]    : default: routine_validation
#   [--cycle]         : default: 01
#   [--use-valgrind]  : default: off
#   [--batch-system]  : <PBS, LyonPBS, LSF, slurm, LyonSlurm, HTCondor, HTCondor_PBS, none>, default: PBS
#   [--queue]         : default: prod. LyonPBS default: P_gdrnu_genie
#   [--softw-topdir]  : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]   : top level dir for job files, default: $PWD
#   [--gen-list]      : comma separated list of event generator list, default all
#   [--chlepton-list] : comma separated list of charged leptons. Both PDGs or names are ok, default electron only
#   [--e-max]         : maxmium energy of the splines in GeV. Default 35 GeV.
#   [--n-knots]       : number of knots per spline. Default 100.
#   [--with-priority] : (boolean) set a priority to optimize bulk productions. Default false
#   [--run-one]       : (boolean) If called, one of the jobs is run as part of the script instead of submitted
#                       via the batch system. Default all the jobs are submitted
#   [--fnal-subjob]   : (boolean) Call if the job is a subjob of the main script to run at the FNAL grid
#   [--fnal-subfile]  : name of sumibmission file to be stored  
#
# Author:
#   Julia Tena-Vidal <j.tena-vidal \st liverpool.ac.uk>
#   University of Liverpool
#
# Copyright:
#   Copyright (c) 2022-2025, The GENIE Collaboration
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
  if($_ eq '--chlepton-list' ) { $req_chlepton_list = $ARGV[$iarg+1]; }
  if($_ eq '--e-max' )         { $e_max	        = $ARGV[$iarg+1]; }
  if($_ eq '--n-knots' )       { $n_knots    	= $ARGV[$iarg+1]; }
  if($_ eq '--with-priority' ) { $priority	= 1; }
  if($_ eq '--run-one' )       { $run_one	= 1; } 
  if($_ eq '--fnal-subjob' )   { $fnal_subjob	= 1; }
  if($_ eq '--fnal-subfile' )  { $fnal_subfile	= $ARGV[$iarg+1]; }
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
$queue_default  = "prod";
if ( $batch_system eq 'LyonPBS' ) {
    $queue_default  = "P_gdrnu_genie" ;
}
if ( $batch_system eq 'LyonSlurm' ) {
    $queue_default  = "htc" ;
}
if ( $batch_system eq 'FNAL' ) { 
    $queue_default = "genie" ; 
}
$queue          = $queue_default                unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE/"   unless defined $softw_topdir;
$jobs_topdir    = $ENV{'PWD'}                   unless defined $jobs_topdir;
$genie_setup    = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
if ( $batch_system eq 'FNAL' ){
    $genie_setup = "setup.sh";
}
$jobs_dir       = "$jobs_topdir/$genie_version-$production\_$cycle-xsec\_chleptonN";
$priority       = 0                             unless defined $priority ;
$e_max          = 35                            unless defined $e_max ;
$n_knots        = 100                           unless defined $n_knots ;
$fnal_subfile   = "fnal_chleptonN"              unless defined $fnal_subfile ;
$batch_script = "$fnal_subfile.fnal";
$xml_script = "$fnal_subfile.xml";


%nucleons_pdg = ( 'n'  =>  1000000010,
                  'p'  =>  1000010010 );

%nucleons_name = ( 1000000010 => 'n' ,
                   1000010010 => 'p' );

@nucleons_proc = ( 'none',
                   'EMRES',
                   'EMDIS',
                   'EMQE' );

%chlepton_pdg_def = ( 'e'      =>   11,
		      'ebar'   =>  -11 );
#		      'mu'     =>   13,
#		      'mubar'  =>  -13,
#		      'tau'    =>   15,
#		      'taubar' =>  -15 );

%chlepton_name_def = ( 11 => 'e'     ,
		       -11 => 'ebar'  );
#		       13 => 'mu'    ,
#		       -13 => 'mubar' ,
#		       15 => 'tau'   ,
#		       -15 => 'taubar' );

##create the list of neutrino to be produced
if ( defined $req_chlepton_list ) {
  my @chlepton_temp_list = split( ",", $req_chlepton_list );
  @chlepton_list = ();
  foreach my $chlepton ( @chlepton_temp_list ) {
    if ( exists $chlepton_pdg_def{$chlepton} ) { push @chlepton_list, $chlepton ; }
    if ( exists	$chlepton_name_def{$chlepton} ) { push @chlepton_list, $chlepton_name_def{$chlepton} ; }
  }
}
else {
    push @chlepton_list, 'e' ; 
}

print "@chlepton_list \n";

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

# Store information required to setup requirements at the FNAL grid
if($batch_system eq 'FNAL'){
    if( ! $fnal_subjob ) {
	unlink "$batch_script" if -e "$batch_script" ; 
	open(FNAL, ">", "$batch_script") or die("Can not create the slurm batch script");
	print FNAL "#!/bin/bash \n";
	print FNAL "source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups \n";
	print FNAL "setup fife_utils \n\n";
	$fnal_opt  = "-G $queue --memory=1GB --disk=20GB --expected-lifetime=8h -N 1 --role=Analysis ";
	$fnal_opt .= "--lines '+FERMIHTC_AutoRelease=True' ";
	$fnal_opt .= "--lines '+FERMIHTC_GraceMemory=4096' --lines '+FERMIHTC_GraceLifetime=6000' ";
	$fnal_opt .= "-f $jobs_topdir/$genie_setup ";
	$fnal_opt .= "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE ";
  
	print FNAL "jobsub_submit_dag $fnal_opt file://$xml_script\n"; 
	close(FNAL);

        # Open xml file
	unlink "$xml_script" if -e "$xml_script" ; 
	open(FNAL_XML, ">", "$xml_script") or die("Can not create the slurm batch script");
	print FNAL_XML "<parallel> \n";
	close(FNAL_XML);
    }
}

foreach $chlepton ( @chlepton_list ) {
  foreach $tgt ( keys %nucleons_pdg ) {
      foreach $proc ( @nucleons_proc_list ) {
	  
	  if ( $proc eq "none" ) {
	      next ;
	  }
	  else {
	      $event_gen_list = $proc ;
	  }

	  $jobname = $chlepton."_on_".$tgt."_$proc";
	  $filename_template = "$jobs_dir/$jobname";

	  $grep_pipe     = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
	  $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
	  $gmkspl_opt    = "-p $chlepton_pdg_def{$chlepton} -t $nucleons_pdg{$tgt} -n $n_knots -e $e_max --event-generator-list $event_gen_list " ;
	  if( $batch_system eq 'FNAL' ) {
	      $gmkspl_opt   .= "-o $jobname.xml ";
	  }
	  else {
	      $gmkspl_opt   .= "-o $filename_template.xml ";
	  }
	  $gmkspl_opt.= " --tune $tune " if ( defined $tune );
    	  $gmkspl_cmd    = "gmkspl $gmkspl_opt";

	  print "@@ exec: $gmkspl_cmd \n";
	  
	  # create sh file 
	  $shell_script = "$filename_template.sh";
	  unlink "$shell_script" if -e "$shell_script" ; 
	  open(COMMANDS, ">$shell_script") or die("Can not create the bash script");
	  print COMMANDS "#!/bin/bash \n";
	  if($batch_system ne 'FNAL'){
	      print COMMANDS "cd $jobs_dir \n";
	      print COMMANDS "source $genie_setup $config_dir \n";
	  } else {
	      print COMMANDS "cd \$CONDOR_DIR_INPUT\n";
	      print COMMANDS "source $genie_setup $config_dir \n";
	      print COMMANDS "cd \$CONDOR_DIR_INPUT\n";
	  }
	  print COMMANDS "$gmkspl_cmd \n";
	  print COMMANDS "ifdh cp -D $jobname.xml $jobs_dir \n" if( $batch_system == 'FNAL');
	  close(COMMANDS);
	  
	  # set executing privileges to the script 
	  `chmod ugo+x $filename_template.sh` ;
	  
	  push( @direct_commands, "bash $filename_template.sh" ) ;
	  
	  # PBS case
	  if($batch_system eq 'PBS' || $batch_system eq 'HTCondor_PBS') {
	      $batch_script = "$filename_template.pbs";
	      open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
	      print PBS "#!/bin/bash \n";
	      print PBS "#PBS -N $jobname \n";
	      print PBS "#PBS -o $filename_template.pbsout.log \n";
	      print PBS "#PBS -e $filename_template.pbserr.log \n";
	      print PBS "#PBS -p -1 \n" if ( $priority ) ;
	      print PBS "source $shell_script \n";
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
	      print PBS "#\$ -P $queue \n";
	      print PBS "#\$ -N $jobname \n";
	      print PBS "#\$ -o $filename_template.pbsout.log \n";
	      print PBS "#\$ -e $filename_template.pbserr.log \n";
	      print PBS "#\$ -l ct=8:00:00,sps=1 \n";
	      print PBS "#\$ -p -1 \n" if ( $priority ) ;
	      print PBS "source $shell_script \n";
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
	      print LSF "#BSUB-q $queue \n";
	      print LSF "#BSUB-o $filename_template.lsfout.log \n";
	      print LSF "#BSUB-e $filename_template.lsferr.log \n";
	      print LSF "source $shell_script \n";	 
	      close(LSF);
	      
	      push( @batch_commands, "bsub < $batch_script " ) ;
	      
	  } #LSF

	  # HTCondor
	  if($batch_system eq 'HTCondor') {
	      $batch_script = "$filename_template.htc";
	      open(HTC, ">$batch_script") or die("Can not create the Condor submit description file: $batch_script");
	      print HTC "Universe               = vanilla \n";
	      print HTC "Executable             = $shell_script \n";
	      print HTC "Log                    = $filename_template.log \n";
	      print HTC "Output                 = $filename_template.out \n";
	      print HTC "Error                  = $filename_template.err \n";
	      print HTC "Request_memory         = 2 GB \n";
	      print HTC "priority               = -1 \n" if ( $priority ) ;
	      print HTC "requirements           = (Opsys =?= \"LINUX\") && (AccessToData =?= True) && (OpSysAndVer =?= \"CentOS7\")  \n" ;
	      print HTC "Queue \n";
	      close(HTC);
	      push ( @batch_commands, "condor_submit $batch_script" ) ;
	  } #HTCondor
	  
	  # slurm case
	  if($batch_system eq 'slurm' || $batch_system eq 'LyonSlurm') {
	      $batch_script = "$filename_template.slr";
	      open(SLURM, ">$batch_script") or die("Can not create the slurm batch script");
	      print SLURM "#!/bin/bash \n";
	      print SLURM "#SBATCH -J $jobname \n";
	      print SLURM "#SBATCH -p $queue \n"; 
	      print SLURM "#SBATCH -o $filename_template.slurmout.log \n";
	      print SLURM "#SBATCH -e $filename_template.slurmerr.log \n";
	      print SLURM "#SBATCH -t 8:0:0 \n";
	      print SLURM "#SBATCH -L sps \n" if ($batch_system eq 'LyonSlurm');
	      print SLURM "#SBATCH --priority -1 \n" if ( $priority ) ; 
	      print SLURM "source $shell_script \n";
	      close(SLURM);
	      
	 push( @batch_commands, "sbatch $batch_script" ) ;
	      
	  } #slurm

	  # FNAL farm
	  if($batch_system eq 'FNAL'){
	      
	      open(FNAL_XML, ">>", "$xml_script") or die("Can not create the slurm batch script");
	      
	      $fnal_opt  = "-n --memory=1GB --disk=20GB --expected-lifetime=10h ";
	      $fnal_opt .= "-f $jobs_topdir/$genie_setup ";
	      $fnal_opt .= "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE ";  
	      $fnal_opt .= "--lines '+FERMIHTC_AutoRelease=True' ";
	      $fnal_opt .= "--lines '+FERMIHTC_GraceMemory=4096' --lines '+FERMIHTC_GraceLifetime=6000' ";
	      print FNAL_XML "jobsub_submit $fnal_opt file://$shell_script\n"; 
	      close FNAL_XML;
	  } #FNAL	  
      }
  }
}


#
# submit
#

if ( $batch_system eq 'none' ) {
    ## run all of them interactively
    for my $run_cmd ( @direct_commands ) {
	print "Executing: $run_cmd \n" ; 
	`$run_cmd` ;
    }
} elsif ($batch_system eq 'FNAL') {
    if( ! $fnal_subjob ){
	open(FNAL_XML, ">>", "$xml_script") or die("Can not create the slurm batch script");
	print FNAL_XML "</parallel> \n";
	close(FNAL_XML);
    }
} else {
    ## submit all except the first
    foreach my $i ( 1 .. $#batch_commands ) {
	`$batch_commands[$i]` ;
    }

    # handle the first according to script options
    if ( defined $run_one ) {
	print "Executing: $direct_commands[0] \n" ;
	`$direct_commands[0]` ;
    }
    else {
	`$batch_commands[0]` ;
    }

}