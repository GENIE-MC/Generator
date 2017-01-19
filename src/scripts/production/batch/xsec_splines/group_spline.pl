#!/usr/bin/perl
#---------------------------------------------------------------------
# Given a directory, the script looks for file in the form 
# v*_on_*_*.xml 
# then it calls the appropriate gspladd to obtain v*_on_*.xml comprehensive
# file
#---------------------------------------------------------------------

$iarg=0;
foreach (@ARGV) {
  if($_ eq '--dir')        { $dir  = $ARGV[$iarg+1]; }
  $iarg++;
}

$dir=$ENV{'PWD'}   unless defined $dir;

opendir(DIR, $dir) or die $!;

%nus = ();
%tags = ();
%lists = ();

while (my $file = readdir(DIR) ) {

  # We only want files
  next unless (-f "$dir/$file");

  # Use a regular expression to find files ending in .txt
  next unless ($file =~ m/\.xml$/);

  next unless ( index($file, "_on_") != -1 );

  my $under_counter = $file =~ tr/_//; 
  next unless ( $under_counter == 3 );

  my $nu = substr($file, 0, index($file, '_') );
  if ( ! exists( $nus{$nu} ) ) { $nus{$nu}=1; } 

  my $tgt = substr($file, index($file, "on_")+3, rindex($file, '_')-index($file, "on_")-3 );
  if ( ! exists( $tgts{$tgt} ) ) { $tgts{$tgt}=1; } 

  my $list = substr($file, rindex($file, '_')+1, index($file,'.')-rindex($file, '_')-1 );
  if ( ! exists( $lists{$list} ) ) { $lists{$list}=1; } 

  print "$file \n";

}

closedir(DIR);

print " Neutrinos: \n";
foreach my $nu ( keys %nus ) {
    print "$nu ";
}

print "\n targets: \n";
foreach my $tgt ( keys %tgts ) {
    print "$tgt ";
}

print "\n Processes: \n";
foreach my $proc ( keys %lists ) {
    print "$proc ";
}
print "\n";

$nus{'ve'}=12        if ( exists( $nus{'ve'} ) );
$nus{'vebar'}=-12    if ( exists( $nus{'vebar'} ) );
$nus{'vmu'}=14       if ( exists( $nus{'vmu'} ) );
$nus{'vmubar'}=-14   if ( exists( $nus{'vmubar'} ) );
$nus{'vtau'}=16      if ( exists( $nus{'vtau'} ) );
$nus{'vtaubar'}=-16  if ( exists( $nus{'vtaubar'} ) );

$nu_list="";
foreach $nu ( keys %nus ) {
    $nu_list .= "," if ( $nu_list ne "" );
    $nu_list .= $nus{$nu};
}
##print $nu_list."\n";

$proc_list = "";
foreach $proc ( keys %lists ) {
  $proc_list .= "," if ($proc_list ne "" );
  $proc_list .= $proc ; 
} 

$tgt_file_list = "";
$tgt_list = "";

foreach my $tgt ( keys %tgts ) {

  $tmp_tgt;
  if ( length($tgt) == 1 ) {
      $tmp_tgt = 1000010010 if ( $tgt eq 'p' );
      $tmp_tgt = 1000000010 if ( $tgt eq 'n' );
  }   
  else { $tmp_tgt = $tgt; }

  $nu_file_list = "";

  foreach my $nu ( keys %nus ) {
  
    $proc_file_list = "";
    
    foreach my $proc ( keys %lists ) {

      $tmp_proc_file="$dir/".$nu."_on_".$tgt."_".$proc.".xml";
      if ( -f $tmp_proc_file ) {
	  $proc_file_list .= "," if ( $proc_file_list ne "" );
          $proc_file_list .= "$tmp_proc_file";
      }
      else {
	  print "Warning: Missing file $tmp_proc_file \n";
      }
    }  ##proc list

  ##  print "$list_file_list"."\n";
  ##  $jobname = "add_".$nu."_on_".$tgt; 
    $tmp_nu_file = "$dir/".$nu."_on_".$tgt.".xml"; 

    ##$grep_pipe     = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
    $gspladd_opt    = "-p $nus{$nu} -t $tmp_tgt  -o $tmp_nu_file -f $proc_file_list --event-generator-list $proc_list";
    $gspladd_cmd    = "gspladd $gspladd_opt";
    print "$gspladd_cmd \n \n";
    `$gspladd_cmd \n`;

    if ( -f $tmp_nu_file ) {
	$nu_file_list .= "," if ( $nu_file_list ne "" );
        $nu_file_list .= "$tmp_nu_file";
    }
    else { 
	print "Error: $tmp_nu_file not created \n";
    }
    
  } ## nu loop 

##  $jobname = "add_".$tgt; 
  $tmp_tgt_file = "$dir/".$tgt.".xml"; 

  ##$grep_pipe     = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
  $gspladd_opt    = "-p $nu_list -t $tmp_tgt -o $tmp_tgt_file -f $nu_file_list --event-generator-list $proc_list";
  $gspladd_cmd    = "gspladd $gspladd_opt";
  print "$gspladd_cmd \n \n";
  `$gspladd_cmd \n`;

  if ( -f $tmp_tgt_file ) {
    $tgt_file_list .= "," if ( $tgt_file_list ne "" );
    $tgt_file_list .= "$tmp_tgt_file";
    $tgt_list .= "," if ( $tgt_list ne "" );
    $tgt_list .= "$tmp_tgt";
  }
  else { 
    print "Error: $tmp_tgt_file not created \n";
  }
  
}  ##tgt loop

$glob_file = $dir."/total_xsec.xml"; 

$tgt_size = keys %tgts;
if ( $tgt_size > 1 ) {

    ##$grep_pipe     = "grep -B 100 -A 30 -i \"warn\\|error\\|fatal\"";
  $gspladd_opt    = "-p $nu_list -t $tgt_list  -o $glob_file -f $tgt_file_list --event-generator-list $proc_list";
  $gspladd_cmd    = "gspladd $gspladd_opt";
  print "$gspladd_cmd \n \n";
  `$gspladd_cmd \n`;
}
elsif ( $tgt_size == 1 ) {

  ###if there is only one target, the global file has to be created manually
  my $temp_name = "$dir/".(keys %tgts)[0].".xml";
  `ln -s $temp_name $glob_file` ;
}

##
## Creating a root file
##

if ( index($tgt_list, "1000010010") == -1 ) { 
  $tgt_list .= ",1000010010";
} 

if ( index($tgt_list, "1000000010") == -1 ) {
  $tgt_list .= ",1000000010";
}

my $cmd = "gspl2root ";
$cmd .= " -p $nu_list ";
$cmd .= " -t $tgt_list ";
$cmd .= " -f $glob_file ";
$cmd .= " -o $dir"."/total_xsec.root " ;
print "\n$cmd\n";
`$cmd`;

