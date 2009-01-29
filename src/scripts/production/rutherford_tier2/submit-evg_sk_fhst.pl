#--------------------------------------------------------------------------
# Submit a GENIE/SK event generation job using a histogram-based neutrino 
# flux description
#
# For use at the RAL/PPD Tier2 PBS batch farm.
#
# Syntax:
#  perl submit-evg_sk_fhst.pl --run run_number --neutrino nu_type
#
# where: 
#   run_number: 1,2,3,...
#   nu_type:    numu, numubar, nue, nuesig
#
# Example:
#  perl submit-evg_sk_fhst.pl --run 180 --neutrino numubar
##
# Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
# STFC, Rutherford Appleton Lab
#--------------------------------------------------------------------------

#!/usr/bin/perl

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--run'     )  { $run      = $ARGV[$iarg+1]; }
  if($_ eq '--neutrino')  { $neutrino = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined run number. Use the --run option]")
unless defined $run;
die("** Aborting [Undefined neutrino type. Use the --neutrino option]")
unless defined $neutrino;

$production     = "test";
$cycle          = "01";
$nevents        = "2000";   
$batch_queue    = "prod";
$time_limit     = "05:00:00";
$genie_inst_dir = "/opt/ppd/t2k/GENIE";
$production_dir = "/opt/ppd/t2k/GENIE/scratch";
$inputs_top_dir = "/opt/ppd/t2k/GENIE/data/job_inputs";
$genie_setup    = "$genie_inst_dir/v2.4.0-setup";
$geom_tgt_mix   = "1000080160[0.8879],1000010010[0.1121]";
$xspl_file      = "$inputs_top_dir/xspl/gxspl-t2k-2.4.0.xml";
$flux_file      = "$inputs_top_dir/t2k_flux/07a/sk/hist.nu.sk.root";
$job_dir        = "$production_dir/sk-$production\_$cycle-$neutrino";
$file_prefix    = "genie_sk";

%mcseed_base    = ( 'numu'    => '183221029',
                    'numubar' => '283221029',
                    'nue'     => '383221029',
                    'nuesig'  => '483221029' );
%mcrun_base     = ( 'numu'    => '10000000',
                    'numubar' => '10000000',
                    'nue'     => '10000000',
                    'nuesig'  => '10000000' );
%nu_pdg_code    = ( 'numu'    =>  '14',
                    'numubar' => '-14',
                    'nue'     =>  '12',
                    'nuesig'  =>  '12' );
%flux_hist_name = ( 'numu'    => 'h100',
                    'numubar' => 'h200',
                    'nue'     => 'h300',
                    'nuesig'  => 'h100' );

die("** Aborting [Can not find GENIE setup script: ..... $genie_setup]") 
unless -e $genie_setup;
die("** Aborting [Can not find flux file: .............. $flux_file]")   
unless -e $flux_file;
die("** Aborting [Can not find xsec file: .............. $xspl_file]")   
unless -e $xspl_file;

# make the jobs directory
# 
mkpath ($job_dir, {verbose => 1, mode=>0777}); 

# form mcrun, mcseed numbers 
#
$mcrun  = $run + $mcrun_base  {$neutrino}; 
$mcseed = $run + $mcseed_base {$neutrino};

print "@@@ Will submit job with MC run number = $mcrun (seed number = $mcseed)\n";

# create the PBS script
#
$batch_script = "$job_dir/skjob-$mcrun.pbs";
open(PBS, ">$batch_script") or die("Can not create the PBS batch script");

$nu  = $nu_pdg_code    {$neutrino};
$hst = $flux_hist_name {$neutrino};

$ghep_file     = "$file_prefix.$production\_$cycle.$mcrun.ghep.root";
$grep_filt     = "grep -B 50 -A 50 -i \"warn\\|error\\|fatal\"";
$valgrind_cmd  = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
$evgen_cmd     = "gT2Kevgen -g $geom_tgt_mix -f $flux_file,$nu\[$hst\] -r $mcrun -n $nevents | $grep_filt &> skjob-$mcrun.log";
$frenm_cmd     = "mv gntp.$mcrun.ghep.root $ghep_file";
$fconv_cmd     = "gntpc -f 1 -i $ghep_file";

print PBS "#!/bin/bash \n";
print PBS "#PBS -l cput=$time_limit \n";
print PBS "#PBS -N $mcrun\_sk-$production-$cycle \n";
print PBS "#PBS -o $job_dir/skjob-$mcrun.pbs_o \n";
print PBS "#PBS -e $job_dir/skjob-$mcrun.pbs_e \n";
print PBS "source $genie_setup \n";
print PBS "cd $job_dir \n";
print PBS "export GSPLOAD=$xspl_file \n";
print PBS "unset GEVGL \n";
print PBS "export GSEED=$mcseed \n";
print PBS "$evgen_cmd \n";
print PBS "$frenm_cmd \n";
print PBS "$fconv_cmd \n";

print "@@@ exec: $evgen_cmd \n";

# submit job
#
`qsub -q $batch_queue $batch_script`;
