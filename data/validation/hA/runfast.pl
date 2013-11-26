## FUNCTIONS:
## This script runs GENIE simulations with the gevgen_hadron command to produce gntp.inuke.***.ghep.root or gntp.***.ghep.root files.
## It then converts those files using the gntpc command if the user selected 'root' as the type.
## If the user selected 'xsec' as the type, it instead uses the ***.ghep.root file and the gtestINukeHadroXSec command to produce a text file with cross sections.



($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);  ## call time function; used to name files and generate initial run number
my @abbr = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
$day = sprintf("%02d", $mday % 100);
$yr = sprintf("%02d", $year % 100);

$iarg = 0;
foreach (@ARGV) {
    if ($_ eq '--n')       { $n         = $ARGV[$iarg+1]; } ## number of events per run
    if ($_ eq '--k')       { $k[0]      = $ARGV[$iarg+1]; } ## kinetic energy 1 (must be defined)
    if ($_ eq '--k1')      { $k[0]      = $ARGV[$iarg+1]; } ## kinetic energy 1 (must be defined)
    if ($_ eq '--k2')      { $k[1]      = $ARGV[$iarg+1]; } ## kinetic energy 2
    if ($_ eq '--k3')      { $k[2]      = $ARGV[$iarg+1]; } ## kinetic energy 3
    if ($_ eq '--k4')      { $k[3]      = $ARGV[$iarg+1]; } ## kinetic energy 4
    if ($_ eq '--k5')      { $k[4]      = $ARGV[$iarg+1]; } ## kinetic energy 5
    if ($_ eq '--k6')      { $k[5]      = $ARGV[$iarg+1]; } ## kinetic energy 6
    if ($_ eq '--min')     { $min_ke    = $ARGV[$iarg+1]; }
    if ($_ eq '--max')     { $max_ke    = $ARGV[$iarg+1]; }
    if ($_ eq '--s')       { $step_size = $ARGV[$iarg+1]; }
    if ($_ eq '--p')       { $prbpdg[0] = $ARGV[$iarg+1]; } ## probe 1 pdg code (must be defined)
    if ($_ eq '--p1')      { $prbpdg[0] = $ARGV[$iarg+1]; } ## probe 1 pdg code (must be defined)
    if ($_ eq '--p2')      { $prbpdg[1] = $ARGV[$iarg+1]; } ## probe 2 pdg code
    if ($_ eq '--r')       { $r         = $ARGV[$iarg+1]; } ## intial run number
    if ($_ eq '--t')       { $tgt[0]    = $ARGV[$iarg+1]; } ## target 1 (must be defined)
    if ($_ eq '--t1')      { $tgt[0]    = $ARGV[$iarg+1]; } ## target 1 (must be defined)
    if ($_ eq '--t2')      { $tgt[1]    = $ARGV[$iarg+1]; } ## target 2
    if ($_ eq '--t3')      { $tgt[2]    = $ARGV[$iarg+1]; } ## target 3
    if ($_ eq '--t4')      { $tgt[3]    = $ARGV[$iarg+1]; } ## target 4
    if ($_ eq '--t5')      { $tgt[4]    = $ARGV[$iarg+1]; } ## target 5
    if ($_ eq '--t6')      { $tgt[5]    = $ARGV[$iarg+1]; } ## target 5
    if ($_ eq '--msg')     { $msg       = $ARGV[$iarg+1]; } ## message thresholds
    if ($_ eq '--v')       { $version   = $ARGV[$iarg+1]; } ## GENIE version
    if ($_ eq '--m')       { $m[0]      = $ARGV[$iarg+1]; } ## GENIE model 1
    if ($_ eq '--m1')      { $m[0]      = $ARGV[$iarg+1]; } ## GENIE model 1
    if ($_ eq '--m2')      { $m[1]      = $ARGV[$iarg+1]; } ## GENIE model 2
    if ($_ eq '--a')       { $author    = $ARGV[$iarg+1]; } ## author for group of runs (will define all necessary parameters)
    if ($_ eq '--type')    { $type      = $ARGV[$iarg+1]; } ## choose to get a root file or a text file with cross sections
    if ($_ eq '--rm')      { $remove    = $ARGV[$iarg+1]; } ## choose to discard gntp files after they're used
    $iarg++;
};

if ($k[0]) {$k[0] = .001 * $k[0]};  ## instead of GeV, take command line energies in units of MeV; this way, manual input is same format for this as for the plotter
if ($k[1]) {$k[1] = .001 * $k[1]};
if ($k[2]) {$k[2] = .001 * $k[2]};
if ($k[3]) {$k[3] = .001 * $k[3]};
if ($k[4]) {$k[4] = .001 * $k[4]};
if ($k[5]) {$k[5] = .001 * $k[5]};
if ($min_ke) {$min_ke = .001 * $min_ke};
if ($max_ke) {$max_ke = .001 * $max_ke};
if ($step_size) {$step_size = .001 * $step_size};


%prb_input_hash = ( '2212' => '2212', '2112' => '2112', '211' => '211', '-211' => '-211', '111' => '111', '321' => '321', '-321' => '-321', '311' => '311', '-311' => '-311', '22' => '22', '13' => '13', '-13' => '13', 'p' => '2212', 
                    'n' => '2112', 'pip' => '211', 'pi+' => '211', 'pim' => '-211', 'pi-' => '-211', 'pi0' => '111', 'kp' => '321', 'k+' => '321', 'km' => '-321', 'k-' => '-321', 'k0' => '311', 'ak0' => '-311', 'gam' => '22', 'gamma' => '22',
                    'mup' => '-13', 'mu+' => '-13', 'mum' => '13', 'mu-' => '13');  ## hash to allow varied inputs for probes
if ($prbpdg[0]) {$prbpdg[0] = $prb_input_hash{$prbpdg[0]}}; 
if ($prbpdg[1]) {$prbpdg[1] = $prb_input_hash{$prbpdg[1]}};

%tgt_input_hash = (  ## a hash to allow different input formats for targets
        'h' => '1', '1' => '1',
	'd' => '1001', '1001' => '1001',  ## i'm just adding 1000 for less common isotopes
        'he' => '2', '2' => '2', 
        'li' => '3', '3' => '3',
	'li6' => '1003', '1003' => '1003',  ## same thing here
        'be' => '4', '4' => '4',
        'b' => '5', '5' => '5',
        'c' => '6', '6' => '6',
        'n' => '7', '7' => '7',
        'o' => '8', '8' => '8',
	'h2o' => '1008', '1008' => '1008',  ## adding 1000 to oxygen to make it mean water
        'al' => '13', '13' => '13',
        'ca' => '20', '20' => '20',
        'fe' => '26', '26' => '26',
        'co' => '27', '27' => '27',
        'ni' => '28', '28' => '28',
        'cu' => '29', '29' => '29',
        'zr' => '40', '40' => '40',
	'nb' => '41', '41' => '41',
	'sn' => '50', '50' => '50',
        'ta' => '73', '73' => '73',
        'pb' => '82', '82' => '82',
        'bi' => '83', '83' => '83',
);
if ($tgt[0]) {$tgt[0] = lc($tgt[0]); $tgt[0] = $tgt_input_hash{$tgt[0]}};
if ($tgt[1]) {$tgt[1] = lc($tgt[1]); $tgt[1] = $tgt_input_hash{$tgt[1]}};
if ($tgt[2]) {$tgt[2] = lc($tgt[2]); $tgt[2] = $tgt_input_hash{$tgt[2]}};
if ($tgt[3]) {$tgt[3] = lc($tgt[3]); $tgt[3] = $tgt_input_hash{$tgt[3]}};
if ($tgt[4]) {$tgt[4] = lc($tgt[4]); $tgt[4] = $tgt_input_hash{$tgt[4]}};
if ($tgt[5]) {$tgt[5] = lc($tgt[5]); $tgt[5] = $tgt_input_hash{$tgt[5]}};

$r = 100000 * $sec + 100 * $yday + $yr  unless defined $r;  ## default initial run number
$type = lc($type);
if ($type ne 'root' && $type ne 'xsec' && $type ne 'both') {error_exit("type")};
%xsec_hash = ('root' => 'no', 'xsec' => 'yes', 'both' => 'yes');
$xsec = $xsec_hash{$type};
%root_hash = ('root' => 'yes', 'xsec' => 'no', 'both' => 'yes');
$root = $root_hash{$type};

if ($author) {
    $author = lc ($author);
    $Author = ucfirst ($author);
};
%author_hash = (  ## acceptable author inputs are those listed at top of script
    'amian' => '1', 'baker' => '1', 'beck' => '1', 'bertrand' => '1', 'carman' => '1', 'chen' => '1', 'cochran' => '1',
    'franz' => '1', 'hautala' => '1', 'hayashi' => '1', 'ingram' => '1', 'iwamoto' => '1', 'kin' => '1', 'levenson' => '1',
    'mcgill' => '1', 'mckeown' => '1', 'meier' => '1', 'otsu' => '1', 'ouyang' => '1', 'roy' => '1', 'segel' => '1', 
    'slypen' => '1', 'stamer' => '1', 'tippawan' => '1', 'tyren' => '1', 'zumbro' => '1'
);
$valid_author = $author_hash {$author};
if ($valid_author ne '1' && $author ne '') { error_exit("author. The author you typed was not recognized") };
if ($msg) { $msg = lc ($msg) };
if ($m[0] eq '') { @m = qw( hA hN ) };  ## run both hA and hN models if user does not specify
%mM_hash = ('ha' => 'hA', 'Ha' => 'hA', 'HA' => 'hA', 'hA' => 'hA', 'hn' => 'hN', 'Hn' => 'hN', 'HN' => 'hN', 'hN' => 'hN');
$m[0] = $mM_hash{$m[0]};
if ($m[1]) {$m[1] = $mM_hash{$m[1]}}; 

if ($author eq '' && $type eq 'root') {
    if ($prbpdg[0] ne '2212' && $prbpdg[0] ne '2112' && $prbpdg[0] ne '211' && $prbpdg[0] ne '-211' && $prbpdg[0] ne '111' && $prbpdg[0] ne '311' &&  $prbpdg[0] ne '-311' && $prbpdg[0] ne '321' &&  $prbpdg[0] ne '-321' 
        &&  $prbpdg[0] ne '22' &&  $prbpdg[0] ne '13' &&  $prbpdg[0] ne '-13') {error_exit("probe")};
    error_exit("energy") unless defined $k[0];  ## exit if author undefined and beam energy undefined
    error_exit("target") unless defined $tgt[0];  ## exit if author undefined and target undefined
};

if ($type eq 'xsec' || $type eq 'both') {
    if ($prbpdg[0] ne '2212' && $prbpdg[0] ne '2112' && $prbpdg[0] ne '211' && $prbpdg[0] ne '-211' && $prbpdg[0] ne '111' && $prbpdg[0] ne '311' &&  $prbpdg[0] ne '-311' && $prbpdg[0] ne '321' &&  $prbpdg[0] ne '-321'
        &&  $prbpdg[0] ne '22' &&  $prbpdg[0] ne '13' &&  $prbpdg[0] ne '-13') {error_exit("probe")};
    error_exit("target") unless defined $tgt[0];
    error_exit("minimum energy") unless defined $min_ke;
    error_exit("maximum energy") unless defined $max_ke;
    error_exit("step size") unless defined $step_size;
};
    

%group_hash = (  ## link author input to all groups associated with that author
    'amian' => ['amian'],
    'baker' => ['baker_c', 'baker_ca'],
    'beck' => ['beck'],
    'bertrand' => ['bertrand'],
    'carman' => ['carman'],
    'chen' => ['chen'],
    'cochran' => ['cochran'],
    'franz' => ['franz'],
    'hautala' => ['hautala'],
    'hayashi' => ['hayashi'],
    'ingram' => ['ingram'],
    'iwamoto' => ['iwamoto_870', 'iwamoto_2100'],
    'kin' => ['kin'],
    'levenson' => ['levenson', 'levenson_c'],
    'mcgill' => ['mcgill'],    
    'mckeown' => ['mckeown'],
    'meier' => ['meier', 'meier_al'],
    'otsu' => ['otsu'],
    'ouyang' => ['ouyang'],
    'roy' => ['roy'],
    'segel' => ['segel'],
    'slypen' => ['slypen_c', 'slypen_fe'],
    'stamer' => ['stamer'],
    'tippawan' => ['tippawan'],
    'tyren' => ['tyren'],
    'zumbro' => ['zumbro'],
);

if ($author ne '') {  ## if author was specified, define parameters accordingly
foreach $group ( @{$group_hash {$author}} ) {

## AUTHOR PRESETS
%prbpdg1_hash = (                                 ## probe 1
    'amian' => '2212',
    'baker_c' => '2212',
    'baker_ca' => '2212',
    'beck' => '2212',
    'bertrand' => '2212',
    'carman' => '2212',
    'chen' => '2212',
    'cochran' => '2212',
    'franz' => '2112',
    'hautala' => '2212',
    'hayashi' => '2112',
    'ingram' => '211',
    'iwamoto_870' => '211',
    'iwamoto_2100' => '211',
    'kin' => '2212',
    'levenson' => '211',
    'levenson_c' => '211',
    'mcgill' => '2212',    
    'mckeown' => '211',
    'meier' => '2212',
    'meier_al' => '2212',
    'otsu' => '2212',
    'ouyang' => '-211',
    'roy' => '2212',
    'segel' => '2212',
    'slypen_c' => '2112',
    'slypen_fe' => '2112',
    'stamer' => '2212',
    'tippawan' => '2112',
    'tyren' => '2212',
    'zumbro' => '211',
);
$prbpdg[0] = $prbpdg1_hash {$group};
%prbpdg2_hash = (                                 ## probe 2
    'iwamoto_870' => '-211',
    'mckeown' => '-211',
);
if ($prbpdg2_hash{$group} ne '') {$prbpdg[1] = $prbpdg2_hash{$group}};
%target1_hash = (                                 ## target 1
    'amian' => '5',
    'baker_c' => '6',
    'baker_ca' => '20',
    'beck' => '26',
    'bertrand' => '26',
    'carman' => '6',
    'chen' => '82',
    'cochran' => '13',
    'franz' => '29',
    'hautala' => '6',
    'hayashi' => '6',
    'ingram' => '1008',
    'iwamoto_870' => '26',
    'iwamoto_2100' => '26',
    'kin' => '6',
    'levenson' => '2',
    'levenson_c' => '6',
    'mcgill' => '6',    
    'mckeown' => '13',
    'meier' => '82',
    'meier_al' => '13',
    'otsu' => '6',
    'ouyang' => '6',
    'roy' => '2',
    'segel' => '6',
    'slypen_c' => '6',
    'slypen_fe' => '26',
    'stamer' => '13',
    'tippawan' => '6',
    'tyren' => '6',
    'zumbro' => '6',
);
$tgt[0] = $target1_hash {$group};
%target2_hash = (                                 ## target 2
    'amian' => '4',
    'beck' => '82',
    'cochran' => '4',
    'hautala' => '20',
    'levenson' => '28',
    'mcgill' => '20',    
    'mckeown' => '4',
    'meier' => '6',
    'ouyang' => '83',
    'roy' => '28',
    'stamer' => '82',
);
if ($target2_hash{$group} ne '') {$tgt[1] = $target2_hash {$group}};
%target3_hash = (                                 ## target 3
    'amian' => '6',
    'cochran' => '6',
    'levenson' => '82',     
    'mckeown' => '6',
    'meier' => '26',
    'roy' => '73',
    'stamer' => '40',    
);
if ($target3_hash{$group} ne '') {$tgt[2] = $target3_hash{$group}};
%target4_hash = (                                 ## target 4
    'amian' => '8',
    'cochran' => '29',
    'mckeown' => '3',
    'meier' => '8',
);
if ($target4_hash{$group} ne '') {$tgt[3] = $target4_hash{$group}};
%target5_hash = (                                 ## target 5
    'amian' => '82',
    'cochran' => '82',    
    'mckeown' => '28',
);
if ($target5_hash{$group} ne '') {$tgt[4] = $target5_hash{$group}};
%target6_hash = (                                 ## target 6
    'cochran' => '1',    
    'mckeown' => '73',
);
if ($target6_hash{$group} ne '') {$tgt[5] = $target6_hash{$group}};
%k1_hash = (                                      ## beam energy 1
    'amian' => '.597',
    'baker_c' => '.318',
    'baker_ca' => '.320',
    'beck' => '.558',
    'bertrand' => '.065',
    'carman' => '.200',
    'chen' => '.290',
    'cochran' => '.730',
    'franz' => '.383',
    'hautala' => '.197',
    'hayashi' => '.147',
    'ingram' => '.114',
    'iwamoto_870' => '.870',
    'iwamoto_2100' => '2.100',
    'kin' => '.300',
    'levenson' => '.100',
    'levenson_c' => '.100',
    'mcgill' => '.800',    
    'mckeown' => '.100',
    'meier' => '.113',
    'meier_al' => '.256',
    'otsu' => '.392',
    'ouyang' => '.500',
    'roy' => '.500',
    'segel' => '.155',
    'slypen_c' => '.0265',
    'slypen_fe' => '.0255',
    'stamer' => '.256',
    'tippawan' => '.0956',
    'tyren' => '.185',
    'zumbro' => '.500',
);
$k[0] = $k1_hash {$group};
%k2_hash = (                                      ## beam energy 2
    'amian' => '.800',
    'franz' => '.425',
    'ingram' => '.163',
    'kin' => '.392',
    'levenson' => '.160',
    'levenson_c' => '.160',
    'mckeown' => '.160',
    'otsu' => '.400',
    'slypen_c' => '.050',
    'slypen_fe' => '.049',
    'stamer' => '.800',
    'iwamoto_870' => '',
);
if ($k2_hash{$group} ne '') {$k[1] = $k2_hash{$group}};
%k3_hash = (                                      ## beam energy 3
    'franz' => '.477',
    'ingram' => '.240',
    'levenson' => '.220',
    'levenson_c' => '.220',
    'mckeown' => '.220',
    'slypen_c' => '.0627',
    'slypen_fe' => '.0627',
    'iwamoto_870' => '',
);
if ($k3_hash{$group} ne '') {$k[2] = $k3_hash{$group}};
%k4_hash = (                                      ## beam energy 4
    'franz' => '.542',
    'levenson_c' => '.300', 
    'slypen_c' => '.0728',
    'iwamoto_870' => '',
);
if ($k4_hash{$group} ne '') {$k[3] = $k4_hash{$group}};
%k5_hash = (                                      ## beam energy 5
    'franz' => '.3174',
    'iwamoto_870' => '',
);
if ($k5_hash{$group} ne '') {$k[4] = $k5_hash{$group}};
%k6_hash = (                                      ## beam energy 6
    'franz' => '.3477',
    'iwamoto_870' => '',
);
if ($k6_hash{$group} ne '') {$k[5] = $k6_hash{$group}};
definitions();  ## call subroutine that defines all parameters that are not group-specific
execute();  ## call subroutine that executes gevgen and gntpc commands
clear_values();
};
} else {  ## if author was not defined, use the manual inputs
    definitions();  ## call subroutine that defines all paramaters that are not group-specific
    open_files();
    execute();  ## call subroutine that executes gevgen and gntpc commands (or gevgen and gtestINukeHadroXSec) 
};




### SUBROUTINES USED IN SCRIPT ###

## The definitions subroutine ##

sub definitions {

    $msg = 'laconic'                                 unless defined $msg;        ## default message thresholds
    $version = '280'                                 unless defined $version;    ## default GENIE version
    if ($version eq '266') {@m = qw(hA)};
    $n = 100000                                      unless defined $n;          ## default number of events per run      

    ## MESSAGE THRESHOLDS
    %msg_hash = (
        'laconic' => '--message-thresholds Messenger_laconic.xml',
        'normal' => '--message-thresholds Messenger.xml',
        'verbose' => '--message-thresholds Messenger_inuke_verbose.xml',
    );
    $msngr = $msg_hash {$msg};

    ## GENIE VERSION
    %version_hash = (
        '280' => '/usr/GENIE/setup_genie',
        '271' => '/usr/GENIEtrunk271/setup_genie',
        '266' => '/usr/GENIE_v266/setup_genie',
    );
    $genie_version = $version_hash {$version};

    %prefix_hash = (
        '280' => 'gntp.inuke',
        '271' => 'gntp.inuke',
        '266' => 'gntp',
    );
    $prefix = $prefix_hash {$version};

    ## PROBES
    %prb_hash = (
        '2212' => 'p',
        '2112' => 'n',
        '211' => 'pip',
        '-211' => 'pim',
        '111' => 'pi0',
        '22' => 'gam',
        '311' => 'k0',
        '-311' => 'ak0',
        '321' => 'kp',
        '-321' => 'km',
        '-13' => 'mup',
        '13' => 'mum',
    );
    $probe = $prb_hash{$probepdg};

    ## TARGETS
    %t_hash = (
        '1' => '1000010010',
	'1001' => '1000010020',
        '2' => '1000020040',
        '3' => '1000030070',
	'1003' => '1000030060',
        '4' => '1000040090',
        '5' => '1000050110',
        '6' => '1000060120',
        '7' => '1000070140',
        '8' => '1000080160',
        '1008' => '1000080160[.8881],1000010010[.1119]',
        '13' => '1000130270',
        '20' => '1000200400',
        '26' => '1000260560',
        '27' => '1000270590',
        '28' => '1000280580',
        '29' => '1000290630',
        '40' => '1000400900',
	'41' => '1000410930',
        '50' => '1000501200',
        '73' => '1000731810',
        '82' => '1000822080',
        '83' => '1000832090',
    );
    $tcode = $t_hash{$target};

    ## TARGETS
    %atom_hash = (
        '1' => 'H',
	'1001' => 'D',
        '2' => 'He',
        '3' => 'Li',
        '1003' => 'Li6',
        '4' => 'Be',
        '5' => 'B',
        '6' => 'C',
        '7' => 'N',
        '8' => 'O',
	'1008' => 'H2O',
        '13' => 'Al',
        '20' => 'Ca',
        '26' => 'Fe',
        '27' => 'Co',
        '28' => 'Ni',
        '29' => 'Cu',
        '40' => 'Zr',
	'41' => 'Nb',
	'50' => 'Sn',
        '73' => 'Ta',
        '82' => 'Pb',
        '83' => 'Bi',
    );
    $Atom = $atom_hash{$target};
    $atom = lc($Atom);

};


## The execute subroutine ##

sub execute {
    $np = @prbpdg; $nt = @tgt; $ne = @k; $nm = @m;
    foreach $probepdg (@prbpdg) {
	foreach $target (@tgt) {
	    definitions();
	    if ($type eq 'root') {
		foreach $nrg (@k) {
		    $nrgmev = $nrg * 1000;
		    foreach $mode (@m) {
			system ("source $genie_version; gevgen_hadron -p $probepdg -t $tcode -n $n -r $r -k $nrg -m $mode $msngr");
			system ("source $genie_version; gntpc -f ginuke -i $prefix.$r.ghep.root -o $abbr[$mon]\_$day\_$yr\_$probe\_$Atom\_$nrgmev\_v$version\_$mode.ginuke.root $msngr");
			if ($remove eq 'yes') {unlink ("$prefix.$r.ghep.root")};
			$r++;
		    };
		};
	    };
	    if ($xsec eq 'yes') {
		$nrg = $min_ke;
		$max_ke = $max_ke + .00000001;
		while ($nrg < $max_ke) {
		    print "NRG: $nrg\n";
		    $nrgmev = $nrg * 1000;
		    foreach $mode (@m) {
			if (-e gevgen_hadron_xsection.txt) {unlink ("gevgen_hadron_xsection.txt")};
			system ("source $genie_version; gevgen_hadron -p $probepdg -t $tcode -n $n -r $r -k $nrg -m $mode $msngr");
			if ($root eq 'yes') {system ("source $genie_version; gntpc -f ginuke -i $prefix.$r.ghep.root -o $abbr[$mon]\_$day\_$yr\_$probe\_$Atom\_$nrgmev\_v$version\_$mode.ginuke.root $msngr")};
			system ("source $genie_version; gtestINukeHadroXSec -f $prefix.$r.ghep.root -w");
			system ("gawk '!/#/ {print}' gevgen_hadron_xsection.txt >> $abbr[$mon]\_$day\_$yr\_$probe\_$Atom\_totxs_v$version\_$mode.txt");
			unlink ("gevgen_hadron_xsection.txt");
			if ($remove eq 'yes') {unlink ("$prefix.$r.ghep.root")};
			$r++;
		    };
		    $nrg = $nrg + $step_size;
		};
	    };
	};   
    };
};


## Subroutine to begin cross section files ##

sub open_files {
    if ($xsec eq 'yes') {
	foreach $probepdg (@prbpdg) {
	    foreach $target (@tgt) {
		foreach $mode (@m) {
		    definitions();
		    $do = 'maybe';
		    (-e "$abbr[$mon]\_$day\_$yr\_$probe\_$Atom\_totxs_v$version\_$mode.txt") ? ($do = 'no') : ($do = 'yes');
		    if ($do eq 'yes') {
			open (FILE, "> $abbr[$mon]\_$day\_$yr\_$probe\_$Atom\_totxs_v$version\_$mode.txt");
			print FILE "#KE     Undef   sig     CEx     sig     Elas    sig     Inelas  sig     Abs     sig     KO      sig     PiPro   sig     DCEx    sig     Reac    sig     Tot     sig     \n";
			close(FILE);
		    };
		};
	    };
	};
    };
};


## The incorrect usage subroutine ##

sub error_exit {
    print "\nThere was a problem with the command line arguments. (Invalid $_[0].) ";
    print "Would you like to get ****.ginuke.root files, text files with total cross sections, or both?\nEnter 'R' for root files, 'XS' for cross section files, or 'B' for both: ";
    $answer = <STDIN>; $answer = uc($answer); chomp ($answer);
    if ($answer ne 'R' && $answer ne 'XS' && $answer ne 'B') {$understood = 'no'};
    while ($understood eq 'no') {
	print "\nAnswer not recognized. Please enter 'R' for root files, 'XS' for cross section files, or 'B' for both: ";
	$answer = <STDIN>; chomp ($answer);
	if ($answer eq 'R' || $answer eq 'XS' || $answer eq 'B') {$understood = 'yes'};
    };
    %hash = ('R' => 'root files', 'XS' => 'text files with cross sections', 'B' => 'both root files and total xs text files');
    $choice = $hash{$answer};
    print "\nYou chose to get $choice. Here's how to do it:\n";
    if ($answer eq 'R') {
	print "\nTo use, type: perl runfast.pl --paramater your_input --paramater your_input --paramater your_input\n\n";
	print "Parameters:\n";
	print "**  --type  : specifies type of output file; enter 'root' to get ****.ginuke.root files\n";
	print "++  --a     : specifies author; see below for valid author inputs\n";
	print "    --v     : specifies GENIE version of root file; use no decimals; assumes 280 if not specified\n";
	print "    --m     : specifies first GENIE model; assumes hA if neither model is specified\n";
	print "    --m2    : specifies second GENIE model; assumes hN if neither model is specified\n";
	print "    --n     : specifies number of events; assumes 100000 if not specified\n";
	print "**  --p     : specifies probe; letter abbreviation or pdg code is acceptable\n";
	print "    --p2    : specifies an additional probe\n";
	print "**  --k     : specifies kinetic energy of probe in units of MeV\n";
	print "    --k2    : specifies an additional kinetic energy of probe\n";
	print "    --k3    : specifies a third kinetic energy\n";
	print "    --k4    : specifies a fourth kinetic energy\n";
	print "    --k5    : specifies a fifth kinetic energy\n";
	print "    --k6    : specifies a sixth kinetic energy\n";
	print "**  --t     : specifies a target; letter symbol or atomic number is acceptable\n";
	print "    --t2    : specifies an additional target\n";
	print "    --t3    : specifies a third target\n";
	print "    --t4    : specifies a fourth target\n";
	print "    --t5    : specifies a fifth target\n";
	print "    --t6    : specifies a sixth target\n";
	print "    --r     : specifies an initial run number; automatically generated if not specified\n";
	print "    --msg   : specifies message thresholds; choose laconic, normal, or verbose; assumes laconic if not specified\n";
	print "    --rm    : enter 'yes' to discard gntp files after they've been used\n";
	print "++ automatically defines all probes, energies, and targets\n";
	print "** necessary inputs\n\n";
	print "Valid Author Inputs:\n";
	print "amian, baker, beck, bertrand, carman, chen, cochran, franz, hautala, hayashi, ingram, iwamoto, kin,\n";
	print "levenson, mcgill, mckeown, meier, otsu, ouyang, roy, segel, slypen, stamer, tippawan, tyren, zumbro\n";
	die("\n");
    };
    if ($answer eq 'XS') {
	print "\nTo use, type: perl runfast.pl --paramater your_input --paramater your_input --paramater your_input\n\n";
	print "Parameters:\n";
	print "**  --type  : specifies type of output file; enter 'xsec' to get text files with total cross sections\n";
	print "    --v     : specifies GENIE version of root file; use no decimals; assumes 280 if not specified\n";
	print "    --m     : specifies first GENIE model; assumes hA if neither model is specified\n";
	print "    --m2    : specifies second GENIE model; assumes hN if neither model is specified\n";
	print "**  --p     : specifies probe; letter abbreviation or pdg code is acceptable\n";
	print "    --p2    : specifies an additional probe\n";
	print "**  --t     : specifies a target; letter symbol or atomic number is acceptable\n";
	print "    --t2    : specifies an additional target\n";
	print "    --t3    : specifies a third target\n";
	print "    --t4    : specifies a fourth target\n";
	print "    --t5    : specifies a fifth target\n";
	print "    --t6    : specifies a sixth target\n";
	print "**  --min   : specifies the minimum kinetic energy of probe in units of MeV\n";
	print "**  --max   : specifies the maximum kinetic energy of probe in units of MeV\n";
	print "**  --s     : specifies the step size of kinetic energy increments\n";
	print "    --n     : specifies number of events per step; assumes 100000 if not specified\n";
	print "    --r     : specifies an initial run number; automatically generated if not specified\n";
	print "    --msg   : specifies message thresholds; choose laconic, normal, or verbose; assumes laconic if not specified\n";
	print "    --rm    : enter 'yes' to discard gntp files after they've been used\n";
	print "** necessary inputs\n";
	die("\n");
    };
    if ($answer eq 'B') {
	print "\nTo use, type: perl runfast.pl --paramater your_input --paramater your_input --paramater your_input\n\n";
	print "Parameters:\n";
	print "**  --type  : specifies type of output file; enter 'both' to get root files and total xs text files\n";
	print "    --v     : specifies GENIE version of root file; use no decimals; assumes 280 if not specified\n";
	print "    --m     : specifies first GENIE model; assumes hA if neither model is specified\n";
	print "    --m2    : specifies second GENIE model; assumes hN if neither model is specified\n";
	print "**  --p     : specifies probe; letter abbreviation or pdg code is acceptable\n";
	print "    --p2    : specifies an additional probe\n";
	print "**  --t     : specifies a target; letter symbol or atomic number is acceptable\n";
	print "    --t2    : specifies an additional target\n";
	print "    --t3    : specifies a third target\n";
	print "    --t4    : specifies a fourth target\n";
	print "    --t5    : specifies a fifth target\n";
	print "    --t6    : specifies a sixth target\n";
	print "**  --min   : specifies the minimum kinetic energy of probe in units of MeV\n";
	print "**  --max   : specifies the maximum kinetic energy of probe in units of MeV\n";
	print "**  --s     : specifies the step size of kinetic energy increments\n";
	print "    --n     : specifies number of events per step; assumes 100000 if not specified\n";
	print "    --r     : specifies an initial run number; automatically generated if not specified\n";
	print "    --msg   : specifies message thresholds; choose laconic, normal, or verbose; assumes laconic if not specified\n";
	print "    --rm    : enter 'yes' to discard gntp files after they've been used\n";
	print "** necessary inputs\n";
	print "Note: Selecting 'both' will give you ****.ginuke.root files at each energy step\n";
	die("\n");
    };
};

## Subroutine to clear values before doing another group ##

sub clear_values {
    undef @prbpdg;
    undef @prb;
    undef @tgt;
    undef @atom;
    undef @t;
    undef @k;
};
