#! /usr/bin/perl

#---------------------------------------------------------------------
# Given a mother directory and a daughter directory, the script tryies
# to copy (ln -s) all the files in the mother dir into the daughter
# If the file alredy exists, the link is not created.
# Mother dir is optional. Default is $GENIE/config
#
#
#   [--mother-dir]      : Source directory 
#   [--daughter-dir]    : Destination directory
#   [--spline-file]     : {no argument} Copies only files with format <vflavor>_on_<target>_<Process>.xml
#---------------------------------------------------------------------

$iarg=0;
foreach (@ARGV) {
  if($_ eq '--mother-dir')        { $mother_dir   = $ARGV[$iarg+1]; }
  if($_ eq '--daughter-dir')      { $daughter_dir = $ARGV[$iarg+1]; }
  if($_ eq '--spline-file')	  { $only_spline  = 1; }
  $iarg++;
}

$mother_dir=$ENV{'GENIE'}."/config"   unless defined $mother_dir;
$daughter_dir=$ENV{'PWD'}             unless defined $daughter_dir;

print "Selected source dir: $mother_dir \n";
print "Selected target dir: $daughter_dir \n";

opendir(DIR, $mother_dir) or die $!;

while (my $file = readdir(DIR) ) {

  # We only want files
  next unless (-f "$mother_dir/$file");

  # Use a regular expression to find files ending in .xml
  next unless ($file =~ m/\.xml$/);

  if ( defined $only_spline ) { 
    next unless ( index($file, "_on_") != -1 );

    my $under_counter = $file =~ tr/_//;
    next unless ( $under_counter == 3 );
  }

  if ( -f "$daughter_dir/$file" ) {
    print "$file already present \n";
  }
  else {
    #  `$move_command`;
      $link_command="ln -s $mother_dir/$file $daughter_dir";
      print "$link_command \n";
      `$link_command`;
  }

}

closedir(DIR);
