#---------------------------------------------------------------------------------------
# Automate running of all routine GENIE MC validation tasks.
#
# Use the cron daemon to run this script every few minutes or so. 
# It will monitor the progress, submit jobs as necessary, assemble the data outputs, 
# run the validation checks and summarize the generator status.
#
# Options:
#    --version          : genie version number
#   [--user]            : username, used for monitoring job progress, default: candreop
#   [--ref-data-topdir] : path to output files from a previous production used for reference
#   [--ref-data-label]  : 
#   [--arch]            : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]      : default: routine_validation
#   [--cycle]           : default: 01
#   [--batch-system]    : <PBS, LSF>, default: PBS (currently works only for PBS)
#   [--queue]           : default: prod
#   [--softw-topdir]  : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]   : top level dir for job files, default: /opt/ppd/t2k/softw/scratch/GENIE/
#
# Author:
#   Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#---------------------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

foreach (@ARGV) {
  if($_ eq '--version')         { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--user')            { $user          = $ARGV[$iarg+1]; }
  if($_ eq '--ref-data-topdir') { $ref_topdir    = $ARGV[$iarg+1]; }
  if($_ eq '--ref-data-label')  { $ref_label     = $ARGV[$iarg+1]; }
  if($_ eq '--arch')            { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')      { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')           { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')    { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')           { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')    { $softw_topdir  = $ARGV[$iarg+1]; }
  if($_ eq '--jobs-topdir')     { $jobs_topdir   = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$user           = "candreop"                    unless defined $user;
$ref_label      = "reference"                   unless defined $ref_label;
$arch           = "SL6.x86_64"                  unless defined $arch;
$production     = "routine_validation"          unless defined $production;
$cycle          = "01"                          unless defined $cycle;
$batch_system   = "PBS"                         unless defined $batch_system;
$queue          = "prod"                        unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE"    unless defined $softw_topdir;
$jobs_topdir    = "/opt/ppd/t2k/scratch/GENIE/" unless defined $jobs_topdir;

$out_data_dir   = "$softw_topdir/data/stage/$genie_version-$production\_$cycle";
$inp_data_dir   = "$softw_topdir/data/job_inputs/";
$genie_setup    = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$genie_topdir   = "$softw_topdir/generator/builds/$arch/$genie_version/";
$scripts_dir    = "$genie_topdir/src/scripts/production/batch/";
$scripts_dir2   = "$genie_topdir/src/scripts/production/misc/";
$status_file    = "$jobs_topdir/$genie_version-$production\_$cycle.status";

$std_args  = "--version $genie_version " .
             "--arch $arch " .
             "--softw-topdir $softw_topdir " .
             "--production $production " .
             "--cycle $cycle " .
             "--batch-system $batch_system " .
             "--queue $queue";

init_status_file($status_file);

$status = read_status_file($status_file);
print "Current status: $status \n";

#
# Create data directories
#
if($status eq "starting")
{
  print "Need to create directory structure for output data \n";
  mkpath ( $out_data_dir,                         {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/xsec/",                  {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/mctest/",                {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/mctest/ghep/",           {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/mctest/gst/",            {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/xsec_validation/",       {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/xsec_validation/ghep/",  {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/xsec_validation/gst/",   {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/hadronization/",         {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/hadronization/ghep/",    {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/intranuke/",             {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/intranuke/ghep/",        {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/reptest/",               {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/reptest/ghep/",          {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/flux/",                  {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/geom/",                  {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/reweight/",              {verbose => 1, mode=>0777});
  mkpath ("$out_data_dir/reports/",               {verbose => 1, mode=>0777});
  update_status_file($status_file,"done creating directory structure for output data");
  exit;
}

#
# ......................................................................................
# Calculate neutrino-nucleon cross-section splines for use in subsequent calculations
# and MC sample generation jobs. Compare the calculated neutrino-nucleon cross-sections 
# with reference calculations and generate report.
# ......................................................................................
#

#
# Run jobs needed
#
if($status eq "done creating directory structure for output data")
{
  print "Need to produce free-nucleon cross-section splines \n";
  $cmd = "perl $scripts_dir/submit_vN_xsec_calc_jobs.pl $std_args --xsplset all";
  execute_command($cmd);
  update_status_file($status_file,"calculating neutrino-nucleon cross-section splines");
  exit;
}

#
# Once all jobs are done, check for errors and merge all free-nucleon cross-section xml files
#
if($status eq "calculating neutrino-nucleon cross-section splines")
{
   $njobs = num_of_jobs_running($user,"vNxscalc");
   if($njobs > 0) {
      exit;
   }
   print "No neutrino-nucleon cross-section calculation job still running... \n";

   # Make sure there is one output XML file for each input PBS script
   print "Making sure all output files are present... \n";
   my $ret = check_number_of_xml_and_pbs_files("$jobs_topdir/$genie_version-$production\_$cycle-xsec\_vN/");
   if($ret != 0) {
     print "The number of output XML files does not match the number of input PBS files. \n";       
     print "Some of the batch jobs have failed. \n"; 
     print "Can not continue with validation runs unless this is sorted out.\n";       
     exit;
   }

   # Check log files for errors
   print "Checking for errors... \n";
   $nerr = `more $jobs_topdir/$genie_version-$production\_$cycle-xsec\_vN/*log | grep -i error | wc -l`
         + `more $jobs_topdir/$genie_version-$production\_$cycle-xsec\_vN/*log | grep -i fatal | wc -l`;
   if($nerr > 0)
   {
      print "I found $nerr errors in the neutrino-nucleon cross-section calculation job log-files. \n";       
      print "Can not continue with validation runs unless this is sorted out.\n";       
      exit;
   }

   # Merge all XML files to one and copy to standard locations
   print "Merging all free nucleon cross-section XML files... \n";
   $cmd = "source $genie_setup; " .
          "gspladd -d $jobs_topdir/$genie_version-$production\_$cycle-xsec\_vN/ -o gxspl-vN-$genie_version.xml; " .
          "cp gxspl-vN-$genie_version.xml $inp_data_dir/xspl/; " .
          "cp gxspl-vN-$genie_version.xml $out_data_dir/xsec/; " .
          "rm -f gxspl-vN-$genie_version.xml";
   execute_command($cmd);
   update_status_file($status_file,"done calculating neutrino-nucleon cross-section splines");
   exit;
}

#
# Convert calculated cross-sections from XML to ROOT format 
#
if($status eq "done calculating neutrino-nucleon cross-section splines")
{
   print "Converting all free-nucleon cross-section files from XML to ROOT format... \n";
   $batch_cmd = "source $genie_setup; " .
          "gspl2root -p 12,-12,14,-14,16,-16 -t 1000010010,1000000010 -o xsec.root -f $inp_data_dir/xspl/gxspl-vN-$genie_version.xml; " .
          "cp xsec.root $out_data_dir/xsec/";
   $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name xsconv $std_args";
   execute_command($cmd);
   update_status_file($status_file,"converting neutrino-nucleon cross-section data to ROOT format");
   exit;
}

#
# Compare calculated cross-sections with reference cross-sections and generate report
#
if($status eq "converting neutrino-nucleon cross-section data to ROOT format")
{
   $njobs = num_of_jobs_running($user,"xsconv");
   if($njobs > 0) {
      exit;
   }
   if (defined $ref_topdir)
   {
        print "Comparing free-nucleon cross-sections with reference calculation... \n";
        $batch_cmd = "source $genie_setup; " .
            "gvld_xsec_comp -f $out_data_dir/xsec/xsec.root,$genie_version -r $ref_topdir/xsec/xsec.root,$ref_label -o xsec.ps; " . 
            "ps2pdf14 xsec.ps; " .
            "cp xsec.pdf $out_data_dir/reports/";
        $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name xscomp $std_args";
        execute_command($cmd);
   }
   update_status_file($status_file,"comparing neutrino-nucleon cross-sections with reference calculations");
   exit;
}

#
# Check last report and make sure it is OK to continue
#
if($status eq "comparing neutrino-nucleon cross-sections with reference calculations")
{
   $njobs = num_of_jobs_running($user,"xscomp");
   if($njobs > 0) {
      exit;
   }

   # ...
   # ... check report
   # ...

   update_status_file($status_file,"done comparing neutrino-nucleon cross-sections with reference calculations");
   exit;
}

#
# ......................................................................................
# Calculate nuclear cross-sections needed for validation MC runs
# ......................................................................................
#

#
# Run jobs needed
#
if($status eq "done comparing neutrino-nucleon cross-sections with reference calculations")
{
  print "Need to produce nuclear cross-section splines \n";
  $cmd = "perl $scripts_dir/submit_vA_xsec_calc_jobs.pl $std_args --config-file $scripts_dir/xsec_splines/genie_test.list";
  execute_command($cmd);
  update_status_file($status_file,"calculating neutrino-nucleus cross-section splines");
  exit;
}

#
# Once all jobs are done, merge all nuclear cross-section xml files
#
if($status eq "calculating neutrino-nucleus cross-section splines")
{
   $njobs = num_of_jobs_running($user,"vAxscalc");
   if($njobs > 0) {
      exit;
   }
   print "No neutrino-nucleus cross-section calculation job still running... \n";

   # Make sure there is one output XML file for each input PBS script
   print "Making sure all output files are present... \n";
   my $ret = check_number_of_xml_and_pbs_files("$jobs_topdir/$genie_version-$production\_$cycle-xsec\_vA\_genie\_test/");
   if($ret != 0) {
     print "The number of output XML files does not match the number of input PBS files. \n";       
     print "Some of the batch jobs have failed. \n"; 
     print "Can not continue with validation runs unless this is sorted out.\n";       
     exit;
   }

   # Check log files for errors
   print "Checking for errors... \n";
   $nerr = `more $jobs_topdir/$genie_version-$production\_$cycle-xsec\_vA\_genie\_test/*log | grep -i error | wc -l`
         + `more $jobs_topdir/$genie_version-$production\_$cycle-xsec\_vA\_genie\_test/*log | grep -i fatal | wc -l`;
   if($nerr >= 0)
   {
      print "I found $nerr errors in the neutrino-nucleus cross-section calculation job log-files. \n";       
      print "Can not continue with validation runs unless this is sorted out.\n";       
      # ------------
      # Lots of errors are printed our by numerical integrators but things look ok
      # Disable exit for now since we will soon be removing all the old numerical algorithms
      # ------------
      # exit;
   }

   # Merge all XML files to one and copy to standard locations
   print "Merging all neutrino-nucleus cross-section XML files... \n";
   $cmd = "source $genie_setup; " .
          "gspladd -d $jobs_topdir/$genie_version-$production\_$cycle-xsec\_vA\_genie\_test/ -o gxspl-vA-$genie_version.xml; " .
          "cp gxspl-vA-$genie_version.xml $inp_data_dir/xspl/; " .
          "cp gxspl-vA-$genie_version.xml $out_data_dir/xsec/; " .
          "rm -f gxspl-vA-$genie_version.xml";
   execute_command($cmd);
   update_status_file($status_file,"done calculating neutrino-nucleus cross-section splines");
   exit;
}

#
# ......................................................................................
# Submit test MC runs, scan the log files, run sanity checks on the event samples
# and compare with reference samples.
# ......................................................................................
#

#
# Run jobs needed
#
if($status eq "done calculating neutrino-nucleus cross-section splines") 
{
  print "Generating standard neutrino MC jobs \n";
  $cmd = "perl $scripts_dir/submit_standard_neutrino_mc_test_jobs.pl $std_args --run all";
  execute_command($cmd);
  update_status_file($status_file,"running standard neutrino MC jobs");
  exit;
}

#
# Once test MC runs are completed, move all data files from the scratch area
#
if($status eq "running standard neutrino MC jobs") 
{
   $njobs = num_of_jobs_running($user,"mctest");
   if($njobs > 0) {
      exit;
   }

   #
   # Check log files for errors and make sure there is one output ROOT file for each input PBS script
   #
   print "Checking for errors... \n";
   $nerr = `more $jobs_topdir/$genie_version-$production\_$cycle-mctest/*log | grep -i error | wc -l`
         + `more $jobs_topdir/$genie_version-$production\_$cycle-mctest/*log | grep -i fatal | wc -l`;
   if($nerr > 0)
   {
      print "I found $nerr errors in the standard neutrino MC job log-files. \n";       
      print "Can not continue with validation runs unless this is sorted out.\n";       
#      exit;
   }
   print "Checking the number of output ROOT/GHEP files... \n";
   opendir my $dir, "$jobs_topdir/$genie_version-$production\_$cycle-mctest/" or die "Can not open directory: $!";
   my @evtfiles  = grep { /.ghep.root$/ } readdir $dir;
   rewinddir $dir;
   my @pbsfiles  = grep { /.pbs$/       } readdir $dir;
   closedir $dir;
   $num_evtfiles = @evtfiles;
   $num_pbsfiles = @pbsfiles;
   print "Found $num_evtfiles ROOT/GHEP event files and $num_pbsfiles PBS files. \n";
   if($num_evtfiles != $num_pbsfiles)   
   {
      print "The number of output event files doesn't match the number of input PBS files. \n";       
      print "Can not continue with validation runs unless this is sorted out.\n";       
      exit;
   }

   $cmd = "mv $jobs_topdir/$genie_version-$production\_$cycle-mctest/*ghep.root $out_data_dir/mctest/ghep/; " .
          "mv $jobs_topdir/$genie_version-$production\_$cycle-mctest/*gst.root  $out_data_dir/mctest/gst/";
   execute_command($cmd);
   update_status_file($status_file,"done running standard neutrino MC jobs");
   exit;
}

#
# Run sanity checks on the test MC runs and generate reports
#
if($status eq "done running standard neutrino MC jobs") 
{
  print "Running sanity checks on the outputs of the the standard neutrino MC jobs ... \n";
  opendir my $dir, "$out_data_dir/mctest/ghep/" or die "Can not open directory: $!";
  my @files = grep { !/^\./ } readdir $dir;
  closedir $dir;
  $ijob = 0;
  foreach(@files) {
     $batch_cmd = 
        "source $genie_setup; " .
        "gvld_sample_scan -f $out_data_dir/mctest/ghep/$_ -o $_.log " .
        "--add-event-printout-in-error-log --event-record-print-level 2 --max-num-of-errors-shown 10 " . 
        "--check-energy-momentum-conservation " .
        "--check-charge-conservation " .
        "--check-for-pseudoparticles-in-final-state " .
        "--check-for-off-mass-shell-particles-in-final-state " .
        "--check-for-num-of-final-state-nucleons-inconsistent-with-target " .
        "--check-vertex-distribution " .
        "--check-decayer-consistency; " .
        "mv $_.log $out_data_dir/reports/";
     $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name snchk-$ijob $std_args";
     execute_command($cmd);
     $ijob++;
  }
  update_status_file($status_file,"running sanity checks on standard neutrino MC jobs");
  exit;
}

#
# See the reports and make sure it is OK to continue
#
if($status eq "running sanity checks on standard neutrino MC jobs") 
{
   $njobs = num_of_jobs_running($user,"snchk");
   if($njobs > 0) {
      exit;
   }

   print "Checking results sanity checks on the outputs of the the standard neutrino MC jobs ... \n";

   # ...
   # ...
   # ...

   update_status_file($status_file,"done running sanity checks on standard neutrino MC jobs");
   exit;
}

#
# Compare the test MC samples with samples generated with a reference version of GENIE
#
if($status eq "done running sanity checks on standard neutrino MC jobs") 
{
  if (defined $ref_topdir)
  {
     opendir my $dir, "$out_data_dir/mctest/ghep/" or die "Can not open directory: $!";
     my @files = grep { !/^\./ } readdir $dir;
     closedir $dir;

     $ijob=0;
     foreach(@files) {
        $batch_cmd = 
           "source $genie_setup; " .
           "gvld_sample_comp -f $out_data_dir/mctest/ghep/$_ -r $ref_topdir/mctest/ghep/$_ -o x; " .
           "mv x $out_data_dir/reports/";
        $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name compmc-$ijob $std_args";
        system("$cmd");    
        $ijob++;
     }
  }
  update_status_file($status_file,"comparing standard neutrino MC jobs with reference samples");
  exit;
}

#
# ......................................................................................
# Repeatability test: Submit different jobs with same inputs and same seed and make sure 
# that the generated samples are identical.
# ......................................................................................
#

if($status eq "comparing standard neutrino MC jobs with reference samples") 
{ 
  # sample 1 
  {
   my $batch_cmd = 
      "source $genie_setup; " .
      "gevgen -p 14 -t 1000260560 -e 0.1,50 -f '1/x' --seed 123456 --run 100 --cross-sections $inp_data_dir/xspl/gxspl-vA-$genie_version.xml; " .
      "mv *ghep.root $out_data_dir/reptest/ghep/";
   my $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name repmc1 $std_args";
   execute_command($cmd);
  }

  # sample 2
  {
   my $batch_cmd = 
      "source $genie_setup; " .
      "gevgen -p 14 -t 1000260560 -e 0.1,50 -f '1/x' --seed 123456 --run 101 --cross-sections $inp_data_dir/xspl/gxspl-vA-$genie_version.xml; " .
      "mv *ghep.root $out_data_dir/reptest/ghep/";
   my $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name repmc2 $std_args";
   execute_command($cmd);
  }

  # sample 3 
  {
   my $batch_cmd = 
      "source $genie_setup; " .
      "gevgen -p 14 -t 1000260560 -e 0.1,50 -f '1/x' --seed 123456 --run 102 --cross-sections $inp_data_dir/xspl/gxspl-vA-$genie_version.xml; " .
      "mv *ghep.root $out_data_dir/reptest/ghep/";
   my $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name repmc2 $std_args";
   execute_command($cmd);
  }

  update_status_file($status_file,"running jobs for MC repeatability test");
  exit;
}

if($status eq "running jobs for MC repeatability test") 
{
   $njobs = num_of_jobs_running($user,"repmc");
   if($njobs > 0) {
     exit;
   }

   {
     my $batch_cmd = 
       "source $genie_setup; " .
       "gvld_repeatability_test --first-sample $out_data_dir/reptest/ghep/gntp.100.ghep.root " .
                               " --second-sample $out_data_dir/reptest/ghep/gntp.101.ghep.root " .
                               " --add-event-printout-in-error-log --max-num-of-errors-shown 10 -o reptest_runs100vs101.log; " .
       "mv reptest_runs100vs101.log $out_data_dir/reports/";
     my $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name reptest1 $std_args";
     execute_command($cmd);
   }

   {
     my $batch_cmd = 
       "source $genie_setup; " .
       "gvld_repeatability_test --first-sample $out_data_dir/reptest/ghep/gntp.100.ghep.root " .
                              " --second-sample $out_data_dir/reptest/ghep/gntp.102.ghep.root " .
                              " --add-event-printout-in-error-log --max-num-of-errors-shown 10 -o reptest_runs100vs102.log; " .
       "mv reptest_runs100vs102.log $out_data_dir/reports/";
     my $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name reptest2 $std_args";
     execute_command($cmd);
   }

   update_status_file($status_file,"running MC repeatability test");
   exit;
}

if($status eq "running MC repeatability test") 
{
   $njobs = num_of_jobs_running($user,"reptest");
   if($njobs > 0) {
      exit;
   }

   # ...
   # ... check log files and make sure it is OK to continue
   # ...

   update_status_file($status_file,"done running MC repeatability test");
   exit;
}

#
# ......................................................................................
# Compare the neutrino cross-section model with world data
# ......................................................................................
#

#
# Run jobs needed 
#
if($status eq "done running MC repeatability test") 
{
   $njobs = num_of_jobs_running($user,"compmc");
   if($njobs > 0) {
      exit;
   }
   print "Generating MC samples needed for comparing the cross-section model with data \n";
   $cmd = "perl $scripts_dir/submit_neutrino_xsec_validation_mc_jobs.pl $std_args --nsubruns 1";
   execute_command($cmd);
   update_status_file($status_file,"running MC jobs for cross-section model validation");
   exit;
}

#
# Once cross-section MC jobs are done, move all data files from the scratch area
#
if($status eq "running MC jobs for cross-section model validation") 
{
   $njobs = num_of_jobs_running($user,"xsecvld");
   if($njobs > 0) {
      exit;
   }
   $cmd = "mv $jobs_topdir/$genie_version-$production\_$cycle-xsec_validation/*ghep.root $out_data_dir/xsec_validation/ghep/; " .
          "mv $jobs_topdir/$genie_version-$production\_$cycle-xsec_validation/*gst.root  $out_data_dir/xsec_validation/gst/";
   execute_command($cmd);
   update_status_file($status_file,"done running MC jobs for cross-section model validation");
   exit;
}

#
# Run actual comparisons with cross-section data and compute error envelopes
#
if($status eq "done running MC jobs for cross-section model validation") 
{   
   # Create an XML file list containing the generated MC files and cross-section files needed for data/MC comparisons
   {
     my $cmd = "perl $scripts_dir2/make_genie_sim_file_list.pl " .
               "$out_data_dir/xsec_validation/ghep/:$out_data_dir/xsec/,$genie_version " .
               "> $out_data_dir/xsec_validation/file_list.xml";
     execute_command($cmd);
   }

   # Create an XML file list, as above, but include the reference samples too
   if (defined $ref_topdir)
   {
     my $cmd = "perl $scripts_dir2/make_genie_sim_file_list.pl " .
               "$out_data_dir/xsec_validation/ghep/:$out_data_dir/xsec/,$genie_version " .
               "$ref_topdir/xsec_validation/ghep/:$ref_topdir/xsec/,$ref_label " .
               "> $out_data_dir/xsec_validation/file_list_ref.xml";
     execute_command($cmd);
   }

   # Submit a single job to generate all GENIE/data comparisons for the current version
   {
     my $batch_cmd = 
        "source $genie_setup; " .
        "gvld_nu_xsec -g $out_data_dir/xsec_validation/file_list.xml -o genie_$genie_version-world_nu_xsec_data_comp-all; " .
        "ps2pdf14 genie_$genie_version-world_nu_xsec_data_comp-all.ps; " .
        "mv genie_$genie_version-world_nu_xsec_data_comp-all.pdf  $out_data_dir/reports/; " .
        "mv genie_$genie_version-world_nu_xsec_data_comp-all.root $out_data_dir/reports/; ";
     my $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name xseccmp1 $std_args";
     execute_command($cmd);
   }

   # Submit a single job to generate all GENIE/data comparisons for the current and reference versions
   if (defined $ref_topdir)
   {
     my $batch_cmd = 
         "source $genie_setup; " .
         "gvld_nu_xsec -g $out_data_dir/xsec_validation/file_list_ref.xml -o genie_$genie_version-world_nu_xsec_data_comp-all-withref; " .
         "ps2pdf14 genie_$genie_version-world_nu_xsec_data_comp-all-withref.ps; " .
         "mv genie_$genie_version-world_nu_xsec_data_comp-all-withref.pdf  $out_data_dir/reports/";
         "mv genie_$genie_version-world_nu_xsec_data_comp-all-withref.root $out_data_dir/reports/";
     my $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name xseccmp2 $std_args";
     execute_command($cmd);
   }

   #
   # Submit jobs to generate GENIE/data comparisons for the current version
   # and calculate error envelopes for the GENIE prediction. Because calculating error
   # envelopes is CPU-intensive submit one job per each GENIE/data comparison.
   #
   my @comparisons = ( 
     "numuCC_all", "numubarCC_all", 
     "numuCC_lowE", "numubarCC_lowE", "numuCC_highE", "numubarCC_highE", 
     "numuCC_minos", "numubarCC_minos", "numuCC_sciboone", "r_minos", 
     "numuCCQE_all", "numuCCQE_deuterium", "numuCCQE_heavy_target", 
     "numuCCQE_nomad_nucleon", "numuCCQE_nomad_nuclear", "numuCCQE_miniboone_nuclear", "numuCCQE_all_12C_nuclear", 
     "numubarCCQE_all", "numubarCCQE_deuterium", "numubarCCQE_heavy_target", 
     "numubarCCQE_nomad_nucleon", "numubarCCQE_nomad_nuclear", 
     "numuCCppip", "numuCCnpip", "numuCCppi0", "numuCCn2pip",  "numuCCppippi0",  "numuCCppippim", 
     "numuCCpi0_numuCCQE_k2k", 
     "numuNCcohpi0_Ne20", "numuCCcohpip_Ne20", "numubarCCcohpim_Ne20", "numuNCcohpi0_Al27", "numuNCcohpi0_Si30", "numuCCcohpip_Si30", "numubarCCcohpim_Si30", 
     "numuCC_dilepton_ratio_worldavg", "numubarCC_dilepton_ratio_worldavg", "numuCC_charm_ratio_worldavg", 
     "numuCC_dilepton_cdhs", "numuCC_dilepton_nomad", "numuCC_dilepton_e744_e770", "numuCC_dilepton_e744", "numuCC_dilepton_fnal15ft", "numuCC_dilepton_gargamelle"
   );

   foreach (@comparisons) {
      my $batch_cmd = 
         "source $genie_setup; " .
         "gvld_nu_xsec -g $out_data_dir/xsec_validation/file_list.xml -o genie_$genie_version-world_nu_xsec_data_comp-$_; " .
         "ps2pdf14 genie_$genie_version-world_nu_xsec_data_comp-$_.ps; " .
         "mv genie_$genie_version-world_nu_xsec_data_comp-$_.pdf  $out_data_dir/reports/";
         "mv genie_$genie_version-world_nu_xsec_data_comp-$_.root $out_data_dir/reports/";
      my $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name xseccmp-$_ $std_args";
      execute_command($cmd);
   }

   update_status_file($status_file,"running cross-section model comparisons with data");
   exit;
}

#
# ......................................................................................
# Compare the hadronization model with data
# ......................................................................................
#

#
# Run jobs needed for comparing the hadronization model with data
#
if($status eq "running cross-section model comparisons with data") 
{
   $njobs = num_of_jobs_running($user,"xseccmp");
   if($njobs > 0) {
      exit;
   }

   print "Generating MC samples needed for comparing the hadronization model with data \n";
   $cmd = "perl $scripts_dir/submit_hadronization_validation_mc_jobs.pl $std_args --nsubruns 1";
   execute_command($cmd);
   update_status_file($status_file,"running MC jobs for hadronization model validation");
   exit;
}

#
# Once hadronization MC jobs are done, move all data files from the scratch area
#
if($status eq "running MC jobs for hadronization model validation") 
{
   $njobs = num_of_jobs_running($user,"hadro");
   if($njobs > 0) {
      exit;
   }
   $cmd = "mv $jobs_topdir/$genie_version-$production\_$cycle-hadronization/*ghep.root $out_data_dir/hadronization/ghep/";
   execute_command($cmd);
   update_status_file($status_file,"done running MC jobs for hadronization model validation");
   exit;
}

#
# Run actual comparisons with hadronization data
#
if($status eq "done running MC jobs for hadronization model validation") 
{
   # Create an XML file list containing the generated MC files    
   {
     my $cmd = "perl $scripts_dir2/make_genie_sim_file_list.pl " .
              "$out_data_dir/hadronization/ghep/,$genie_version " .
              "> $out_data_dir/hadronization/file_list.xml";
     execute_command($cmd);
   }

   # Create an XML file list as above but include the reference MC files, if any.
   if (defined $ref_topdir)
   {
     my $cmd = "perl $scripts_dir2/make_genie_sim_file_list.pl " .
               "$out_data_dir/hadronization/ghep/,$genie_version " .
               "$ref_topdir/hadronization/ghep/,$ref_label " .
               "> $out_data_dir/hadronization/file_list_ref.xml";
     execute_command($cmd);
   }

   # Submit a single job to generate all GENIE/data comparisons for the current version
   {
     my $batch_cmd = 
       "source $genie_setup; " .
       "gvld_hadronz_test -g $out_data_dir/hadronization/file_list.xml -o genie_$genie_version-hadronization_test.ps; " .
       "ps2pdf14 genie_$genie_version-hadronization_test.ps; " .
       "mv genie_$genie_version-hadronization_test.pdf $out_data_dir/reports/";
     my $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name hadcmp1 $std_args";
     execute_command($cmd);
   }

   # Submit a single job to generate all GENIE/data comparisons for the current and reference versions
   if (defined $ref_topdir)
   {
     my $batch_cmd = 
         "source $genie_setup; " .
         "gvld_hadronz_test -g $out_data_dir/hadronization/file_list_ref.xml -o genie_$genie_version-hadronization_test-withref.ps; " .
         "ps2pdf14 genie_$genie_version-hadronization_test-withref.ps; " .
         "mv genie_$genie_version-hadronization_test-withref.pdf $out_data_dir/reports/";
     my $cmd = "perl $scripts_dir/submit.pl --cmd \'$batch_cmd\' --job-name hadcmp2 $std_args";
     execute_command($cmd);
   }

   update_status_file($status_file,"running hadronization model comparisons with data");
   exit;
}

#
# ......................................................................................
# Compare the intranuclear rescattering model with data
# ......................................................................................
#


#...
#...
#...


#
# ......................................................................................
# Compare the structure function model with data
# ......................................................................................
#

#...
#...
#...

#
# ......................................................................................
# Run event reweighting tests
# ......................................................................................
#

#...
#...
#...

#
# ......................................................................................
# Run flux driver tests
# ......................................................................................
#

#...
#...
#...

#
# ......................................................................................
# Run geometry navigation tests
# ......................................................................................
#

#...
#...
#...




#
########################################################################################
########################################################################################
#

sub num_of_jobs_running 
{
   $user   = shift;
   $jobtag = shift;
   @myjobs = `qstat -u $user`;
   $njobs = 0;
   foreach(@myjobs) {
     chomp($_);
     if( /$user/ && /$jobtag/ ) {
        $njobs++;
        print "Job still running: $_ \n";
     }
   }
   print "Number of jobs still running or on the queue: $njobs \n";
   return $njobs;
}

sub init_status_file
{
  $file = shift;
  if(! -e $file)
  {
    open(STF,">$file");
    print STF "starting";
    close STF;
  }
}

sub update_status_file
{
  $file = shift;
  $mesg = shift;
  open(STF,">$file");
  print STF $mesg;
  close STF;
}

sub read_status_file
{
  $file = shift;
  open(STF,"<$file");
  $status = <STF>; 
  chomp($status);
  close STF;
  return $status;
}

sub execute_command
{
  $cmd = shift;
  print "Executing: $cmd \n";
  my $ret = system("$cmd");
  if ($ret != 0) 
  { 
      die "** Failed executing: $cmd \n"; 
  }
}

sub check_number_of_xml_and_pbs_files
{
  $dirname = shift;
  print "Searching for XML and PBS files in directory: $dirname \n";
  opendir my $dir, $dirname or die "Can not open directory: $!";
  my @xmlfiles = grep { /.xml$/ } readdir $dir;
  rewinddir $dir;
  my @pbsfiles = grep { /.pbs$/ } readdir $dir;
  closedir $dir;
#  print "XML files found: @xmlfiles \n";
#  print "PBS files found: @pbsfiles \n";
  $nxmlfiles = @xmlfiles;
  $npbsfiles = @pbsfiles;
  print "Found $nxmlfiles XML files and $npbsfiles PBS files. \n";
  if($nxmlfiles != $npbsfiles)   
  {
    return 1;
  }  
  return 0;
}
