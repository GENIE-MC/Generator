#---------------------------------------------------------------------------------------
# Automate running of all routine GENIE MC validation tasks.
#
# Use cron to run this script every minute or so. It will monitor the progress, submit 
# jobs as necessary, assemble the data outputs, run the validation checks and summarize 
# the generator status.
#
# Options:
#    --version       : genie version number
#   [--user]          : username, used for monitoring job progress, default: candreop
#   [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]   : default: routine_validation
#   [--cycle]        : default: 01
#   [--batch-system] : <PBS, LSF>, default: PBS
#   [--queue]        : default: prod
#   [--softw-topdir] : default: /opt/ppd/t2k/softw/GENIE
#
# Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#---------------------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

foreach (@ARGV) {
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--user')          { $user          = $ARGV[$iarg+1]; }
  if($_ eq '--arch')          { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')    { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')         { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')  { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')         { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')  { $softw_topdir  = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;

$user           = "candreop"                    unless defined $user;
$arch           = "SL5_64bit"                   unless defined $arch;
$production     = "routine_validation"          unless defined $production;
$cycle          = "01"                          unless defined $cycle;
$batch_system   = "PBS"                         unless defined $batch_system;
$queue          = "prod"                        unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE"    unless defined $softw_topdir;

$out_data_dir   = "$softw_topdir/data/stage/$genie_version-$production\_$cycle";
$inp_data_dir   = "$softw_topdir/data/job_inputs/";
$job_dir        = "$softw_topdir/scratch/";
$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$genie_topdir   = "$softw_topdir/builds/$arch/$genie_version/";
$scripts_dir    = "$genie_topdir/src/scripts/production/batch/";
$status_file    = "$job_dir/$genie_version-$production\_$cycle.status";

$std_args  = "--version $genie_version --arch $arch --softw-topdir $softw_topdir " .
             "--production $production --cycle $cycle " .
             "--batch-system $batch_system --queue $queue";

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
  update_status_file($status_file,"done creating directory structure for output data");
  exit;
}

#
# Calculate neutrino-nucleon cross-section splines 
#
if($status eq "done creating directory structure for output data")
{
  print "Need to produce free-nucleon cross-section splines \n";
  $cmd = "perl $scripts_dir/submit_vN_xsec_calc_jobs.pl $std_args --xsplset all";
  print "Running $cmd \n";  
  system("$cmd $arg");
  update_status_file($status_file,"calculating neutrino-nucleon cross-section splines");
  exit;
}

#
# Once all jobs are done, merge all free-nucleon cross-section xml files
#
if($status eq "calculating neutrino-nucleon cross-section splines")
{
   $njobs = num_of_jobs_running($user,"vNxscalc");
   if($njobs > 0) {
      exit;
   }
   else {
      $cmd = "source $genie_setup; " .
             "gspladd -d $job_dir/$genie_version-$production\_$cycle-xsec\_vN/ -o gxspl-vN-$genie_version.xml; " .
             "cp gxspl-vN-$genie_version.xml $inp_data_dir/xspl/; " .
             "cp gxspl-vN-$genie_version.xml $out_data_dir/xsec/";
      system("$cmd");
      update_status_file($status_file,"done calculating neutrino-nucleon cross-section splines");
      exit;
   }
}

#
# Calculate nuclear cross-sections needed for validation MC runs
#
if($status eq "done calculating neutrino-nucleon cross-section splines")
{
  print "Need to produce nuclear cross-section splines \n";
  $cmd = "perl $scripts_dir/submit_vA_xsec_calc_jobs.pl $std_args --config-file $scripts_dir/xsec_splines/genie_test.list";
  print "Running $cmd \n";  
  system("$cmd $arg");
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
   else {
      $cmd = "source $genie_setup; " .
             "gspladd -d $job_dir/$genie_version-$production\_$cycle-xsec\_vA\_genie\_test/ -o gxspl-vA-$genie_version.xml; " .
             "cp gxspl-vA-$genie_version.xml $inp_data_dir/xspl/";
      system("$cmd");
      update_status_file($status_file,"done calculating neutrino-nucleus cross-section splines");
      exit;
   }
}

#
# Submit test MC runs
#
if($status eq "done calculating neutrino-nucleus cross-section splines") 
{
  print "Generating standard neutrino MC jobs \n";
  $cmd = "perl $scripts_dir/submit_standard_neutrino_mc_test_jobs.pl $std_args --run all";
  print "Running $cmd \n";  
  system("$cmd $arg");
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
   else {
      $cmd = "mv $job_dir/$genie_version-$production\_$cycle-mctest/*ghep.root $out_data_dir/mctest/ghep/; " .
             "mv $job_dir/$genie_version-$production\_$cycle-mctest/*gst.root $out_data_dir/mctest/gst/";
      system("$cmd");
      update_status_file($status_file,"done running standard neutrino MC jobs");
      exit;
   }
}


# ...............................................................
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
   print "There are $njobs jobs still running.\n";
   return $njobs;
}

# ...............................................................
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

# ...............................................................
sub update_status_file
{
  $file = shift;
  $mesg = shift;
  open(STF,">$file");
  print STF $mesg;
  close STF;
}

# ...............................................................
sub read_status_file
{
  $file = shift;
  open(STF,"<$file");
  $status = <STF>; 
  chomp($status);
  close STF;
  return $status;
}

# ...............................................................

