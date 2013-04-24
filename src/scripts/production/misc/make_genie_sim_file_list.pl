#-------------------------------------------------------------------------------------------
# Create XML file list expected by GENIE validation programs
#
# Syntax:
# % perl make_genie_sim_file_list.pl <inputs>
#
# where inputs is a list of file_path,model_name pairs
#
# Example:
# %perl make_genie_sim_file_list.pl /some/path/,model_name /another/path,another_model_name ...
#
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#-------------------------------------------------------------------------------------------

#!/usr/bin/perl

my %data = ();

foreach (@ARGV) {
  my ($file_path, $model_name) = split(',', $_);
  $data{$model_name} .= $file_path;
}

print "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?> \n\n";
print "<genie_simulation_outputs> \n";

while( my ($model_name, $file_path) = each %data )
{
  #print "$model_name --> $file_path \n";
  print "\n\n";
  print "\t <model name=\"$model_name\"> \n";
  opendir my $dir, $file_path or die "Can not open directory: $!";
  my @files = grep { !/^\./ } readdir $dir;
  closedir $dir;
  foreach(@files) {
      if(/.ghep.root$/) {
          print "\t\t <evt_file format=\"ghep\"> $_ </evt_file> \n";
      }
      elsif(/.gst.root$/) {
          print "\t\t <evt_file format=\"gst\"> $_ </evt_file> \n";
      } 
      else {
          print "\t\t <xsec_file> $_ </xsec_file> \n";
      }
  }

  print "\t </model> \n\n";
}

print "</genie_simulation_outputs> \n";

