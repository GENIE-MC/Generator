#!/usr/bin/perl

#----------------------------------------------------------------------------------
# A perl script to download jnubeam hbook ntuples and convert them to ROOT format.
#
# Syntax:
#  perl fetch.pl 
#           --username <some_username> 
#           --passwd   <some_password> 
#           --url      <some_url> 
#           --prefix   <file_prefix> 
#           --suffix   <file_suffix>
#           --rmin     <min_run_nu> 
#           --rmax     <max_run_nu>
#          [--target   <local_disk_path (default: '.')>]
#
# Example:
#  You want to download/convert the 07a near detector flux files
#  (named as 'nu.nd280.<RUN_NUMBER>.hbk)
#  living at http://jnusrv01.kek.jp/internal/t2k/nubeam/flux/07a/nd/ 
#  You only want to download / convert the files corresponding to MC runs 20->80, 
#  and you want to copy the downloaded files in /some/local/datadir/
#  You need to type:
#
# perl fetch.pl 
#          --username t2k 
#          --passwd type_actual_passwd 
#          --url http://jnusrv01.kek.jp/internal/t2k/nubeam/flux/07a/nd/ 
#          --prefix nu.nd280. --suffix .hbk --rmin 20 --rmax 80
#          --target /some/local/datadir/
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#----------------------------------------------------------------------------------
#

$iarg=0;
foreach (@ARGV) {
   if($_ eq '--username') { $username   = $ARGV[$iarg+1]; }
   if($_ eq '--passwd')   { $passwd     = $ARGV[$iarg+1]; }
   if($_ eq '--url'   )   { $url        = $ARGV[$iarg+1]; }
   if($_ eq '--prefix')   { $file_prfx  = $ARGV[$iarg+1]; }
   if($_ eq '--suffix')   { $file_sufx  = $ARGV[$iarg+1]; }
   if($_ eq '--rmin')     { $run_min    = $ARGV[$iarg+1]; }
   if($_ eq '--rmax')     { $run_max    = $ARGV[$iarg+1]; }
   if($_ eq '--target')   { $target     = $ARGV[$iarg+1]; }
   $iarg++;
}
die("** Aborting [Undefined user name. Use the --username option]")
unless defined $username;
die("** Aborting [Undefined URL. Use the --url option]")
unless defined $url;
die("** Aborting [Undefined file prefix. Use the --prefix option]")
unless defined $file_prfx;
die("** Aborting [Undefined file suffix. Use the --suffix option]")
unless defined $file_sufx;
die("** Aborting [Undefined minimum run number. Use the --rmin option]")
unless defined $run_min;
die("** Aborting [Undefined maximum run number. Use the --rmax option]")
unless defined $run_max;

$passwd = "" unless defined $passwd;

$i=0;
for($i=$run_min; $i<=$run_max; $i++) {

# fetch hbook file from remote location
  $filename = $file_prfx.$i.$file_sufx;
  $fetch_cmd = "curl -O --user $username:$passwd $url/$filename";
  print "@@@ $fetch_cmd \n";
  system("$fetch_cmd");
  if( $? != 0 ) { 
    print "*** Couldn't download $url/$filename \n";
  }

# copy to local directory
  if( -defined $target ) {
     print "@@@ mv $filename  $target \n";
     `mv $filename  $target`;
  }
}

