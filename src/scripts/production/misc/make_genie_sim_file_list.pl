#-------------------------------------------------------------------------------------------
# Create XML file list expected by GENIE validation programs
#
# Syntax:
# % perl make_genie_sim_file_list.pl <inputs>
#
# where inputs is a list of "file_paths,model_name" pairs.
# Multiple ':'-separated file paths can be specified.
#
# Example:
# %perl make_genie_sim_file_list.pl /some/path/,model_name /another/path1:/another/path2,another_model_name ...
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#-------------------------------------------------------------------------------------------

#!/usr/bin/perl

my %data = ();

foreach (@ARGV) {
  my ($file_paths, $model_name) = split(',', $_);
  $data{$model_name} .= $file_paths;
}

print "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?> \n\n";
print "<genie_simulation_outputs> \n";

while( my ($model_name, $file_paths) = each %data )
{
  #print "$model_name --> $file_paths \n";
  print "\n\n";
  print "\t <model name=\"$model_name\"> \n";

  @file_path_array = split(':',$file_paths);
  foreach(@file_path_array) {
    $dirname = $_;
    opendir my $dir, $dirname or die "Can not open directory: $!";
    my @files = grep { !/^\./ } readdir $dir;
    closedir $dir;
    foreach(@files) {
      $filename = $_;
      if($filename =~ m/.ghep.root$/) {
          print "\t\t <evt_file format=\"ghep\"> $dirname/$filename </evt_file> \n";
      }
      elsif($filename =~ m/.gst.root$/) {
          print "\t\t <evt_file format=\"gst\"> $dirname/$filename </evt_file> \n";
      } 
      elsif($filename =~ m/.root$/) {
          print "\t\t <xsec_file> $dirname/$filename </xsec_file> \n";
      }
    }
  }
  print "\t </model> \n\n";
}

print "</genie_simulation_outputs> \n";

