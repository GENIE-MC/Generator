#-----------------------------------------------------------------------------
# Automate running of all GENIE MC validation tasks
#
# Options:
#    --version       : genie version number
#   [--arch]         : <SL4_32bit, SL5_64bit>, default: SL5_64bit
#   [--production]   : default: <version>
#   [--cycle]        : default: 01
#   [--batch-system] : <PBS, LSF>, default: PBS
#   [--queue]        : default: prod
#   [--softw-topdir] : default: /opt/ppd/t2k/softw/GENIE
#
# Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#-----------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

foreach (@ARGV) {
  if($_ eq '--version')       { $genie_version = $ARGV[$iarg+1]; }
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

$user           = "candreop";
$arch           = "SL5_64bit"                   unless defined $arch;
$production     = "routine_validation"          unless defined $production;
$cycle          = "01"                          unless defined $cycle;
$batch_system   = "PBS"                         unless defined $batch_system;
$queue          = "prod"                        unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE"    unless defined $softw_topdir;

$genie_setup    = "$softw_topdir/builds/$arch/$genie_version-setup";
$status_file    = "$softw_topdir/scratch/$genie_version-$production\_$cycle.status";

$data_dir       = "$softw_topdir/data/job_inputs/";
$scripts_dir    = "$softw_topdir/builds/$arch/$genie_version/src/scripts/production/batch/";
$job_dir        = "$softw_topdir/scratch/";

$std_args  = "--version $genie_version --arch $arch --production $production --cycle $cycle --batch-system $batch_system --queue $queue --softw-topdir $softw_topdir";

init_status_file($status_file);

$status = read_status_file($status_file);
print "Current status: $status \n";

if($status eq "starting validation runs")
{
  print "Need to produce free-nucleon cross-section splines \n";
  $cmd = "perl $scripts_dir/submit_vN_xsec_calc_jobs.pl $std_args  --xsplset all";
  print "Running $cmd \n";  
  system("$cmd $arg");
  update_status_file($status_file,"calculating neutrino-nucleon cross-section splines");
  exit;
}

if($status eq "calculating neutrino-nucleon cross-section splines")
{
   $njobs = num_of_jobs_running($user,"vNxscalc");
   if($njobs > 0) {
      exit;
   }
   else {
      $cmd = "source $genie_setup; " .
             "gspladd -d $job_dir/$genie_version-$production\_$cycle-xsec\_vN/ -o gxspl-vN-$genie_version.xml; " .
             "cp gxspl-vN-$genie_version.xml $data_dir/xspl/";
      system("$cmd");
      update_status_file($status_file,"done calculating neutrino-nucleon cross-section splines");
      exit;
   }
}


if($status eq "done calculating neutrino-nucleon cross-section splines")
{
  print "Need to produce nuclear cross-section splines \n";
  $cmd = "perl $scripts_dir/submit_vA_xsec_calc_jobs.pl $std_args  --xsplset all";
  print "Running $cmd \n";  
  system("$cmd $arg");
  update_status_file($status_file,"calculating neutrino-nucleus cross-section splines");
  exit;
}


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

sub init_status_file
{
  $file = shift;
  if(! -e $file)
  {
    open(STF,">$file");
    print STF "starting validation runs";
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

