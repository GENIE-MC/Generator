#!/usr/bin/perl
#---------------------------------------------------------------------
# Given a directory, the script looks for file in the form 
# *_on_*_*.xml 
# then it calls the appropriate gspladd to obtain *_on_*.xml comprehensive
# file. The file name will be total_xsec.xml
# The scpript also prodces a root output for the splines
# 
#  Options:
#  [--dir]             : Is the directory which contains all the xml files to be converted. Default: $PWD
#  [--tune]            : Tune option
#  
#---------------------------------------------------------------------

$iarg=0;
foreach (@ARGV) {
  if($_ eq '--dir')             { $dir            = $ARGV[$iarg+1]; }
  if($_ eq '--tune')            { $tune           = $ARGV[$iarg+1]; }
  $iarg++;
}

$dir=$ENV{'PWD'}      unless defined $dir;
$def_chlepton_list = "11"   unless defined $def_chlepton_list ;

opendir(DIR, $dir) or die $!;

%chleptons = ();
%tags = ();
%lists = ();

while (my $file = readdir(DIR) ) {

  # We only want files
  next unless (-f "$dir/$file");

  # Use a regular expression to find files ending in .xml
  next unless ($file =~ m/\.xml$/);

  next unless ( index($file, "_on_") != -1 );

  my $under_counter = $file =~ tr/_//; 
  next unless ( $under_counter == 3 );

  my $chlepton = substr($file, 0, index($file, '_') );
  if ( ! exists( $chleptons{$chlepton} ) ) { $chleptons{$chlepton}=1; } 

  my $tgt = substr($file, index($file, "on_")+3, rindex($file, '_')-index($file, "on_")-3 );
  if ( ! exists( $tgts{$tgt} ) ) { $tgts{$tgt}=1; } 

  my $list = substr($file, rindex($file, '_')+1, index($file,'.')-rindex($file, '_')-1 );
  if ( ! exists( $lists{$list} ) ) { $lists{$list}=1; } 

  print "$file \n";

}

closedir(DIR);

print " Charged leptons: \n";
foreach my $chlepton ( keys %chleptons ) {
    print "$chlepton ";
}

print "\n Targets: \n";
foreach my $tgt ( keys %tgts ) {
    print "$tgt ";
}

print "\n Processes: \n";
foreach my $proc ( keys %lists ) {
    print "$proc ";
}
print "\n";

# Charged chleptons
$chleptons{'e'}=11        if ( exists( $chleptons{'e'} ) );
$chleptons{'ebar'}=-11    if ( exists( $chleptons{'ebar'} ) );
$chleptons{'mu'}=13       if ( exists( $chleptons{'mu'} ) );
$chleptons{'mubar'}=-13   if ( exists( $chleptons{'mubar'} ) );
$chleptons{'tau'}=15      if ( exists( $chleptons{'tau'} ) );
$chleptons{'taubar'}=-15  if ( exists( $chleptons{'taubar'} ) );

$chlepton_list="";
foreach $chlepton ( keys %chleptons ) {
    $chlepton_list .= "," if ( $chlepton_list ne "" );
    $chlepton_list .= $chleptons{$chlepton};
}
##print $chlepton_list."\n";
$glob_file = $dir."/total_xsec.xml"; 

## Check if there is only one target and one lepton
$tgt_size = keys %tgts;
$tgt_size = keys %chleptons;

if ( $tgt_size == 1 && $chleptons ==1 ) {

    ###If there is only one target and one charged lepton, the global file has to be created manually
    my $temp_name = "$dir/".(keys %chleptons)[0]."_on".(keys %tgts)[0]."_EM.xml";
    
    my $cmd = "ln -s "; 
    
    $cmd .= "cp $temp_name $glob_file ";
    print "Executing: $cmd \n" ; 
    `$cmd` ;
    
} else {
    ## Combine all the splines into one
    $file_list = "";

    foreach my $tgt ( keys %tgts ) {
	
	$tmp_tgt;
	if ( length($tgt) == 1 ) {
	    $tmp_tgt = 1000010010 if ( $tgt eq 'p' );
	    $tmp_tgt = 1000000010 if ( $tgt eq 'n' );
	}   
	else { $tmp_tgt = $tgt; }

	$chlepton_file_list = "";	
	foreach my $chlepton ( keys %chleptons ) {
	    $temp_file = "$dir/".$chlepton."_on_".$tgt."_EM.xml";
	    if( -f $temp_file ) {
		$file_list .= "," if ( $file_list ne "" ) ; 
		$file_list .= "$temp_file";
	    }
	    else {
		print "Warning: missing file $temp_file";
	    }    
	}
    }
    
    $gspladd_opt    = " -o $glob_file -f $file_list ";
    $gspladd_cmd    = "gspladd $gspladd_opt";
    print "$gspladd_cmd \n \n";
    `$gspladd_cmd \n`;
}

if ( ! -f $glob_file ) {
    print "Error: $tmp_nu_file not created \n";
}



