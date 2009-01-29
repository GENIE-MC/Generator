#!/usr/bin/perl

#----------------------------------------------------------------------------------
# A perl script to download jnubeam hbook ntuples and convert them to ROOT format.
#
# Syntax:
#  perl fetch-t2k_flux.pl 
#           --passwd <usual_t2k_password> 
#           --url <web_url> 
#           --prefix <file_prefix> 
#          [--suffix <file_suffix (default:'.hbk')>] 
#           --rmin <min_run_nu> 
#           --rmax <max_run_nu>
#          [--target <local_disk_path (default: '.')>]
#
# Example:
#  You want to download/convert the 07a near detector flux files
#  (named as 'nu.nd280.<RUN_NUMBER>.hbk)
#  living at http://jnusrv01.kek.jp/internal/t2k/nubeam/flux/07a/nd/ 
#  You only want to download / convert the files corresponding to MC runs 20->80, 
#  and you want to copy all *hbk and *root files in /some/local/datadir/
#  You need to type:
#
# perl fetch-t2k_flux.pl 
#          --passwd type_actual_passwd 
#          --url http://jnusrv01.kek.jp/internal/t2k/nubeam/flux/07a/nd/ 
#          --prefix nu.nd280. --rmin 20 --rmax 80
#          --target /some/local/datadir/
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#----------------------------------------------------------------------------------
#

$user           = "t2k";
$passwd         = "";
$url            = "";
$file_prfx      = "";
$file_sufx      = ".hbk";
$run_min        = 1;
$run_max        = 999999;

$iarg=0;
foreach (@ARGV) {
   if($_ eq '--passwd')   { $passwd     = $ARGV[$iarg+1]; }
   if($_ eq '--url'   )   { $url        = $ARGV[$iarg+1]; }
   if($_ eq '--prefix')   { $file_prfx  = $ARGV[$iarg+1]; }
   if($_ eq '--suffix')   { $file_sufx  = $ARGV[$iarg+1]; }
   if($_ eq '--rmin')     { $run_min    = $ARGV[$iarg+1]; }
   if($_ eq '--rmax')     { $run_max    = $ARGV[$iarg+1]; }
   if($_ eq '--target')   { $target     = $ARGV[$iarg+1]; }
   $iarg++;
}

$i=0;
for($i=$run_min; $i<=$run_max; $i++) {

# fetch hbook file from remote location
  $filename_hbk = $file_prfx.$i.$file_sufx;
  $fetch_cmd = "curl -O --user $user:$passwd $url/$filename_hbk";
  print "@@@ $fetch_cmd \n";
  system("$fetch_cmd");
  if( $? != 0 ) { 
    print "*** Couldn't download $url/$filename_hbk \n";
  }

# convert to root
  $filename_root = $filename_hbk;
  $filename_root =~ s/.hbk/.root/;
  $convert_cmd = "h2root $filename_hbk";
  print "@@@ $convert_cmd \n";
  `$convert_cmd`;

# copy to local directory
  if( -defined $target ) {
     print "@@@ mv $filename_hbk  $target \n";
     `mv $filename_hbk  $target`;
     print "@@@ mv $filename_root $target \n";
     `mv $filename_root $target`;
  }
}

