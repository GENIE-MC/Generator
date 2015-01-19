#!/usr/bin/perl

#---------------------------------------------------------------------------------------------------------------------
# Submit standard neutrino event generation jobs for GENIE release validation
# The outputs can be compared with outputs from past releases using GENIE's gvld_sample_comp utility.
# Sanity checks can be performed using GENIE's gvld_sample_scan utility.
#
# Syntax:
#   shell% perl submit_standard_neutrino_mc_test_jobs.pl <options>
#
# Options:
#    --version       : GENIE version number
#    --run           : Comma separated list of run numbers
#   [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]   : default: routine_validation
#   [--cycle]        : default: 01
#   [--ref-samples]  : Path for reference samples, default: no reference samples / no plots will be generated
#   [--use-valgrind] : default: off
#   [--batch-system] : <PBS, LSF, slurm, none>, default: PBS
#   [--queue]        : default: prod
#   [--softw-topdir] : default: /opt/ppd/t2k/softw/GENIE
#
# Examples:
#   shell% perl submit_standard_neutrino_mc_test_jobs.pl \
#               --production 2.5.1_prelease_tests --cycle 01 --version v2.5.1 --run 1001
#   shell% perl submit_standard_neutrino_mc_test_jobs.pl \
#                --production 2.5.1_prelease_tests --cycle 01 --version v2.5.1 --run 1000,1001,9203
#   shell% perl submit_standard_neutrino_mc_test_jobs.pl \
#                --production 2.5.1_prelease_tests --cycle 01 --version v2.5.1 --run all
#
# Tested at the RAL/PPD Tier2 PBS batch farm.
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#---------------------------------------------------------------------------------------------------------------------
#
# SAMPLES:
#......................................................................
#  run   | nev  |  init state      | energy  | processes
#  nu.   |      |                  | (GeV)   | enabled
#......................................................................
#  1000  | 100k |  numu    + n     |  0.5    | all       
#  1001  | 100k |  numu    + n     |    1    | all       
#  1002  | 100k |  numu    + n     |    5    | all       
#  1003  | 100k |  numu    + n     |   50    | all       
#  1101  | 100k |  numubar + p     |    1    | all       
#  1102  | 100k |  numubar + p     |    5    | all       
#  1103  | 100k |  numubar + p     |   50    | all       
#  2001  | 100k |  numu    + Fe56  |    1    | all       
#  2002  | 100k |  numu    + Fe56  |    5    | all       
#  2003  | 100k |  numu    + Fe56  |   50    | all       
#  2101  | 100k |  numubar + Fe56  |    1    | all      
#  2102  | 100k |  numubar + Fe56  |    5    | all       
#  2103  | 100k |  numubar + Fe56  |   50    | all       
#  9001  | 100k |  numu    + Fe56  |    5    | DIS charm 
#  9002  | 100k |  numu    + Fe56  |    5    | QEL charm 
#  9101  | 100k |  numu    + Fe56  |    2    | COH CC+NC 
#  9201  | 100k |  nue     + Fe56  |    1    | ve elastic
#  9202  | 100k |  numu    + Fe56  |    1    | ve elastic
#  9203  |  50k |  numu    + Fe56  |   20    | IMD
#  9204  |  50k |  nuebar  + Fe56  |   20    | IMD (annihilation)
#......................................................................
#

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--version')       { $genie_version   = $ARGV[$iarg+1]; }
  if($_ eq '--run')           { $runnu           = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch            = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production      = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle           = $ARGV[$iarg+1]; }
  if($_ eq '--ref-samples')   { $ref_sample_path = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')  { $use_valgrind    = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')  { $batch_system    = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue           = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir    = $ARGV[$iarg+1]; }  
  $iarg++;
}
die("** Aborting [Undefined benchmark runs #. Use the --run option]")
unless defined $runnu;
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$use_valgrind    = 0                          unless defined $use_valgrind;
$arch            = "SL5_64bit"                unless defined $arch;
$production      = "routine_validation"       unless defined $production;
$cycle           = "01"                       unless defined $cycle;
$batch_system    = "PBS"                      unless defined $batch_system;
$queue           = "prod"                     unless defined $queue;
$softw_topdir    = "/opt/ppd/t2k/softw/GENIE" unless defined $softw_topdir;
$ref_sample_path = 0                          unless defined $ref_sample_path;
$genie_setup     = "$softw_topdir/builds/$arch/$genie_version-setup";
$jobs_dir        = "$softw_topdir/scratch/$genie_version-$production\_$cycle-mctest";
$xspl_file       = "$softw_topdir/data/job_inputs/xspl/gxspl-vA-$genie_version.xml";
$mcseed          = 210921029;

%nevents_hash = ( 
  '1000' => '100000',
  '1001' => '100000',
  '1002' => '100000',
  '1003' => '100000',
  '1101' => '100000',
  '1102' => '100000',
  '1103' => '100000',
  '2001' => '100000',
  '2002' => '100000',
  '2003' => '100000',
  '2101' => '100000',
  '2102' => '100000',
  '2103' => '100000',
  '9001' => '100000',
  '9002' => '100000',
  '9101' => '100000',
  '9201' =>  '50000',
  '9202' =>  '50000',
  '9203' =>  '50000',  
  '9204' =>  '50000'  
);

%nupdg_hash = ( 
  '1000' =>  '14',
  '1001' =>  '14',
  '1002' =>  '14',
  '1003' =>  '14',
  '1101' => '-14',
  '1102' => '-14',
  '1103' => '-14',
  '2001' =>  '14',
  '2002' =>  '14',
  '2003' =>  '14',
  '2101' => '-14',
  '2102' => '-14',
  '2103' => '-14',
  '9001' =>  '14',
  '9002' =>  '14',
  '9101' =>  '14',
  '9201' =>  '12',
  '9202' =>  '14',
  '9203' =>  '14',  
  '9204' =>  '-12'  
);

%tgtpdg_hash = ( 
  '1000' => '1000000010',
  '1001' => '1000000010',
  '1002' => '1000000010',
  '1003' => '1000000010',
  '1101' => '1000010010',
  '1102' => '1000010010',
  '1103' => '1000010010',
  '2001' => '1000260560',
  '2002' => '1000260560',
  '2003' => '1000260560',
  '2101' => '1000260560',
  '2102' => '1000260560',
  '2103' => '1000260560',
  '9001' => '1000260560',
  '9002' => '1000260560',
  '9101' => '1000260560',
  '9201' => '1000260560',
  '9202' => '1000260560',
  '9203' => '1000260560', 
  '9204' => '1000260560' 
);

%energy_hash = ( 
  '1000' =>  '0.5',
  '1001' =>  '1.0',
  '1002' =>  '5.0',
  '1003' => '50.0',
  '1101' =>  '1.0',
  '1102' =>  '5.0',
  '1103' => '50.0',
  '2001' =>  '1.0',
  '2002' =>  '5.0',
  '2003' => '50.0',
  '2101' =>  '1.0',
  '2102' =>  '5.0',
  '2103' => '50.0',
  '9001' =>  '5.0',
  '9002' =>  '5.0',
  '9101' =>  '2.0',
  '9201' =>  '1.0',
  '9202' =>  '1.0',
  '9203' => '20.0',
  '9204' => '20.0'
);

%gevgl_hash = ( 
  '1000' => 'Default',
  '1001' => 'Default',
  '1002' => 'Default',
  '1003' => 'Default',
  '1101' => 'Default',
  '1102' => 'Default',
  '1103' => 'Default',
  '2001' => 'Default',
  '2002' => 'Default',
  '2003' => 'Default',
  '2101' => 'Default',
  '2102' => 'Default',
  '2103' => 'Default',
  '9001' => 'CharmCCDIS',
  '9002' => 'CharmCCQE',
  '9101' => 'COH',
  '9201' => 'NuEElastic',
  '9202' => 'NuEElastic',
  '9203' => 'IMD',
  '9204' => 'IMD'  
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

print "Input runs: $runnu \n";

for my $curr_runnu (keys %gevgl_hash)  {
  print "Checking benchmark run: ...... $curr_runnu \n";

  if($runnu=~m/$curr_runnu/ || $runnu eq "all") {
    print "** matched -> submitting job \n";

    #
    # get runnu-dependent info
    #
    $nev   = $nevents_hash {$curr_runnu};
    $nu    = $nupdg_hash   {$curr_runnu};
    $tgt   = $tgtpdg_hash  {$curr_runnu};
    $en    = $energy_hash  {$curr_runnu};
    $gevgl = $gevgl_hash   {$curr_runnu};

    $jntemplate    = "mctest-$curr_runnu";
    $fntemplate    = "$jobs_dir/$jntemplate";
    $grep_pipe     = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
    $valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
    $evgen_opt     = "-n $nev -e $en -p $nu -t $tgt -r $curr_runnu --seed $mcseed --cross-sections $xspl_file --event-generator-list $gevgl";
    $evgen_cmd     = "gevgen $evgen_opt | $grep_pipe &> $fntemplate.evgen.log";
    $conv_cmd      = "gntpc -f gst -i gntp.$curr_runnu.ghep.root | $grep_pipe &> $fntemplate.conv.log";
    $comp_cmd      = "gvld_sample_comp -f gntp.$curr_runnu.gst.root -r $ref_sample_path/gntp.$curr_runnu.gst.root | $grep_pipe &> $fntemplate.comp.log";

    print "@@ exec: $evgen_cmd \n";

    #
    # submit
    #
  
    # PBS case
    if($batch_system eq 'PBS') {
        $batch_script  = "$fntemplate.pbs";
        open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
        print PBS "#!/bin/bash \n";
        print PBS "#PBS -N $jntemplate \n";
        print PBS "#PBS -o $fntemplate.pbsout.log \n";
        print PBS "#PBS -e $fntemplate.pbserr.log \n";
        print PBS "source $genie_setup \n"; 
        print PBS "cd $jobs_dir \n";
        print PBS "$evgen_cmd \n";
        print PBS "$conv_cmd \n";
        if(-d $ref_sample_path) {
           print PBS "$comp_cmd \n";
        }
        close(PBS);
        `qsub -q $queue $batch_script`;
    } #PBS

    # LSF case
    if($batch_system eq 'LSF') {
        $batch_script  = "$fntemplate.sh";
        open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
        print LSF "#!/bin/bash \n";
        print PBS "#BSUB-j $jntemplate \n";
        print LSF "#BSUB-q $queue \n";
        print LSF "#BSUB-o $fntemplate.lsfout.log \n";
        print LSF "#BSUB-e $fntemplate.lsferr.log \n";
        print LSF "source $genie_setup \n"; 
        print LSF "cd $jobs_dir \n";
        print LSF "$evgen_cmd \n";
        print LSF "$conv_cmd \n";
        if(-d $ref_sample_path) {
           print LSF "$comp_cmd \n";
        }
        close(LSF);
        `bsub < $batch_script`;
    } #LSF

    # slurm case
    if($batch_system eq 'slurm') {
        $batch_script  = "$fntemplate.sh";
        open(SLURM, ">$batch_script") or die("Can not create the SLURM batch script");
        print SLURM "#!/bin/bash \n";
        print SLURM "#SBATCH-p $queue \n";
        print SLURM "#SBATCH-o $fntemplate.lsfout.log \n";
        print SLURM "#SBATCH-e $fntemplate.lsferr.log \n";
        print SLURM "source $genie_setup \n"; 
        print SLURM "cd $jobs_dir \n";
        print SLURM "$evgen_cmd \n";
        print SLURM "$conv_cmd \n";
        if(-d $ref_sample_path) {
           print SLURM "$comp_cmd \n";
        }
        close(SLURM);
        `sbatch --job-name=$jntemplate $batch_script`;
    } #slurm

    # no batch system, run jobs interactively
    if($batch_system eq 'none') {
        system("source $genie_setup; cd $jobs_dir; $evgen_cmd; $conv_cmd");
	if(-d $ref_sample_path){
	    system("$comp_cmd");
	}
    } # interactive mode

  }
}
