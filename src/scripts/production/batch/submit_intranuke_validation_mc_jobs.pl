#!/usr/bin/perl

#-----------------------------------------------------------------------------------------------------------
# Submit jobs to generate data needed for validating GENIE's hadron transport model
#
# Syntax:
#   perl submit_intranuke_validation_mc_jobs.pl <options>
#
# Options:
#    --version        : GENIE version number
#    --run            : Runs to submit (Can be a run number, or a comma separated list of run numbers.) 
#                       Use `--run all' to submit all jobs. 
#                       Can also specify runs used for comparisons with data from a specific author using 
#                       the author name, eg `--run iwamoto', or `--run iwamoto,ingram'.
#                       Can also specify runs by probe, eg `--run piplus', or `--run piplus,piminus,proton'.
#   [--inuke-model]   : Physics model, <hA, hN>, default: hA
#   [--model-enum]    : Physics model enumeration, default: 01
#   [--nsubruns]      : Number of subruns per run, default: 1
#   [--arch]          : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default: SL6.x86_64
#   [--production]    : Production name, default: <model>_<version>
#   [--cycle]         : Cycle in current production, default: 01
#   [--use-valgrind]  : Use valgrind? default: off
#   [--batch-system]  : <PBS, LSF, slurm, HTCondor, HTCondor_PBS, none>, default: PBS
#   [--queue]         : Batch queue, default: prod
#   [--softw-topdir]  : top level dir for softw installations, default: /opt/ppd/t2k/softw/GENIE/
#   [--jobs-topdir]   : top level dir for job files, default: /opt/ppd/t2k/scratch/GENIE/
#
# EVENT SAMPLES:
#
# Run number key: PPTTTEEEEEMMxx
#
# PP    : probe  (10:pi+, 11:pi0, 12:pi-, 20:K+, 21:K0, 22:K-, 30:gamma, 40:p, 41:n)
# TTT   : target (002: He4, 006:C12, 008:O16, 013:Al27, 020:Ca40, 026:Fe56: 028:Ni58, 029:Cu63, 082:Pb208, 083:Bi209)
# EEEEE : kinetic energy (MeV)
# MM    : physics model enumeration, default: 01
# xx    : sub-run ID, 00-99, 100k events each
#
#.......................................................................................
# run number       |  init state      | kin energy   | req.  | group of                | 
#                  |                  |   (GeV)      | stat  | runs                    |     
#                  |                  |              |(# evt)| runs                    |     
#.......................................................................................
# 1002600870MMxx   | pi+   + Fe56     |   0.870      | 0.5M  | iwamoto                 | 
# 1202600870MMxx   | pi-   + Fe56     |   0.870      | 1.0M  | iwamoto                 | 
# 1008200870MMxx   | pi+   + Pb208    |   0.870      | 1.0M  | iwamoto                 | 
# 1002602100MMxx   | pi+   + F56      |   2.100      | 1.0M  | iwamoto                 | 
# 1000601400MMxx   | pi+   + C12      |   1.400      | 1.0M  | shibata                 | 
# 1002901400MMxx   | pi+   + Cu63     |   1.400      | 1.0M  | shibata                 | 
# 1008201400MMxx   | pi+   + Pb208    |   1.400      | 1.0M  | shibata                 | 
# 1000600220MMxx   | pi+   + C12      |   0.220      | 1.0M  | mckeown,levenson        | 
# 1002800220MMxx   | pi+   + Ni58     |   0.220      | 1.0M  | mckeown,levenson        | 
# 1008200220MMxx   | pi+   + Pb208    |   0.220      | 1.0M  | mckeown,levenson        | 
# 1000600160MMxx   | pi+   + C12      |   0.160      | 1.0M  | mckeown,levenson        | 
# 1002800160MMxx   | pi+   + Ni58     |   0.160      | 1.0M  | mckeown,levenson        | 
# 1008200160MMxx   | pi+   + Pb208    |   0.160      | 1.0M  | mckeown,levenson        | 
# 1000600100MMxx   | pi+   + C12      |   0.100      | 1.0M  | mckeown,levenson        | 
# 1002800100MMxx   | pi+   + Ni58     |   0.100      | 1.0M  | mckeown,levenson        | 
# 1008200100MMxx   | pi+   + Pb208    |   0.100      | 1.0M  | mckeown,levenson        | 
# 1200600220MMxx   | pi-   + C12      |   0.220      | 1.0M  | mckeown                 | 
# 1202800220MMxx   | pi-   + Ni58     |   0.220      | 1.0M  | mckeown                 | 
# 1208200220MMxx   | pi-   + Pb208    |   0.220      | 1.0M  | mckeown                 | 
# 1200600160MMxx   | pi-   + C12      |   0.160      | 1.0M  | mckeown                 | 
# 1202800160MMxx   | pi-   + Ni58     |   0.160      | 1.0M  | mckeown                 | 
# 1208200160MMxx   | pi-   + Pb208    |   0.160      | 1.0M  | mckeown                 | 
# 1200600100MMxx   | pi-   + C12      |   0.100      | 1.0M  | mckeown                 | 
# 1202800100MMxx   | pi-   + Ni58     |   0.100      | 1.0M  | mckeown                 | 
# 1208200100MMxx   | pi-   + Pb208    |   0.100      | 1.0M  | mckeown                 |  
# 1000600500MMxx   | pi+   + C12      |   0.500      | 1.0M  | zumbro                  | 
# 1200600500MMxx   | pi-   + C12      |   0.500      | 1.0M  | ouyang                  | 
# 1208300500MMxx   | pi-   + Bi209    |   0.500      | 1.0M  | ouyang                  | 
# 1000600300MMxx   | pi+   + C12      |   0.300      | 1.0M  | levenson                | 
# 1000200300MMxx   | pi+   + He4      |   0.300      | 1.0M  | mckeown,levenson        | 
# 1000200220MMxx   | pi+   + He4      |   0.220      | 1.0M  | mckeown,levenson        | 
# 1000200160MMxx   | pi+   + He4      |   0.160      | 1.0M  | mckeown,levenson        | 
# 1000200100MMxx   | pi+   + He4      |   0.100      | 1.0M  | mckeown,levenson        | 
# 1000800114MMxx   | pi+   + O16      |   0.114      | 1.0M  | ingram                  | 
# 1000800163MMxx   | pi+   + O16      |   0.163      | 1.0M  | ingram                  | 
# 1000800240MMxx   | pi+   + O16      |   0.240      | 1.0M  | ingram                  | 
# 4000600800MMxx   | p     + C12      |   0.800      | 1.0M  | mcgill,amian            | 
# 4002901400MMxx   | p     + Cu63     |   1.400      | 1.0M  | shibata                 | 
# 4008201400MMxx   | p     + Pb208    |   1.400      | 1.0M  | shibata                 | 
# 4000601400MMxx   | p     + C12      |   1.400      | 1.0M  | shibata                 | 
# 4002000800MMxx   | p     + Ca40     |   0.800      | 1.0M  | mcgill                  | 
# 4002000800MMxx   | p     + Pb208    |   0.800      | 1.0M  | mcgill                  | 
# 4000600730MMxx   | p     + C12      |   0.730      | 3.0M  | cochran                 | 
# 4001300730MMxx   | p     + Al27     |   0.730      | 3.0M  | cochran                 | 
# 4002900730MMxx   | p     + Cu63     |   0.730      | 3.0M  | cochran                 | 
# 4008200730MMxx   | p     + Pb208    |   0.730      | 3.0M  | cochran                 | 
# 4000600800MMxx   | p     + C12      |   0.800      | 1.0M  | amian                   | 
# 4008200800MMxx   | p     + Pb208    |   0.800      | 1.0M  | amian                   | 
# 4000600597MMxx   | p     + C12      |   0.597      | 1.0M  | amian                   | 
# 4008200597MMxx   | p     + Pb208    |   0.597      | 1.0M  | amian                   | 
# 4002600558MMxx   | p     + Fe56     |   0.558      | 1.0M  | beck                    | 
# 4002600558MMxx   | p     + Pb208    |   0.558      | 1.0M  | beck                    | 
# 4000600300MMxx   | p     + C12      |   0.300      | 1.0M  | kin                     | 
# 4000600392MMxx   | p     + C12      |   0.392      | 1.0M  | kin                     | 
# 4002600065MMxx   | p     + Fe56     |   0.065      | 1.0M  | bertrand                | 
# 4001300256MMxx   | p     + Al27     |   0.256      | 1.0M  | stamer                  | 
# 4008200256MMxx   | p     + Pb208    |   0.256      | 1.0M  | stamer                  | 
# 4000600200MMxx   | p     + C12      |   0.200      | 0.5M  | carman                  | 
# 4000600197MMxx   | p     + C12      |   0.197      | 0.5M  | hautala                 | 
# 4000600113MMxx   | p     + C12      |   0.113      | 0.5M  | meier                   | 
# 4002600113MMxx   | p     + Fe56     |   0.113      | 0.5M  | meier                   | 
# 4008200113MMxx   | p     + Pb208    |   0.113      | 0.5M  | meier                   | 
# 1008200220MMxx   | pi+   + Pb208    |   0.220      | 1.0M  | levenson                | 
#.......................................................................................
#
# OUTPUTS:
#  - gntp.<inuke_mode>.<PPTTTEEEEEMMxx>.ghep.root   GHEP event file
#  - gntp.<inuke_mode>.<PPTTTEEEEEMMxx>.ginuke.root GINUKE summary/analysis ntuple
#
#
# Examples:
#
# (1) Submit (in an LSF farm) 10-subruns (100k events each) of run 4000600597, using GENIE v2.7.1:
#     % perl submit_intranuke_validation_mc_jobs.pl --version v2.7.1 \
#                     --nsubruns 10 --batch-system LSF --run 4000600597
#
# (2) Submit (in an LSF farm) 10-subruns (100k events each) of runs 4000600597,1002600870 and 1000800240,
#     using GENIE v2.7.1:
#     % perl submit_intranuke_validation_mc_jobs.pl --version v2.7.1 
#                --nsubruns 10 --batch-system LSF --run 4000600597,1002600870,1000800240
#
# (3) Submit (in an LSF farm) 10-subruns (100k events each) of *all* runs, using GENIE v2.7.1:
#     % perl submit_intranuke_validation_mc_jobs.pl --version v2.7.1 \
#                 --nsubruns 10 --batch-system LSF --run all
#
# (4) Submit (in an LSF farm) 10-subruns (100k events each) of all runs used for comparisons with the
#     `iwamoto' data, using GENIE v2.7.1:
#     % perl submit_intranuke_validation_mc_jobs.pl --version v2.7.1 \
#                 --nsubruns 10 --batch-system LSF --run iwamoto
#
# (5) Submit (in an LSF farm) 10-subruns (100k events each) of all runs with a pi+ or pi- probe,
#     using GENIE v2.7.1:
#     % perl submit_intranuke_validation_mc_jobs.pl --version v2.7.1 \
#                 --nsubruns 10 --batch-system LSF --run piplus,piminus
#
# Author:
#   Costas Andreopoulos <costas.andreopoulos \st stfc.ac.uk>
#   University of Liverpool & STFC Rutherford Appleton Laboratory
#
# Copyright:
#   Copyright (c) 2003-2019, The GENIE Collaboration
#   For the full text of the license visit http://copyright.genie-mc.org
#-----------------------------------------------------------------------------------------------------------

use File::Path;

# inputs
#
$iarg=0;
foreach (@ARGV) {
  if($_ eq '--nsubruns')       { $nsubruns      = $ARGV[$iarg+1]; }
  if($_ eq '--run')            { $runnu         = $ARGV[$iarg+1]; }
  if($_ eq '--inuke-model')    { $inuke_model   = $ARGV[$iarg+1]; }
  if($_ eq '--model-enum')     { $model_enum    = $ARGV[$iarg+1]; }
  if($_ eq '--version')        { $genie_version = $ARGV[$iarg+1]; }
  if($_ eq '--arch')           { $arch          = $ARGV[$iarg+1]; }
  if($_ eq '--production')     { $production    = $ARGV[$iarg+1]; }
  if($_ eq '--cycle')          { $cycle         = $ARGV[$iarg+1]; }
  if($_ eq '--use-valgrind')   { $use_valgrind  = $ARGV[$iarg+1]; }
  if($_ eq '--batch-system')   { $batch_system  = $ARGV[$iarg+1]; }
  if($_ eq '--queue')          { $queue         = $ARGV[$iarg+1]; }
  if($_ eq '--softw-topdir')   { $softw_topdir  = $ARGV[$iarg+1]; }  
  if($_ eq '--jobs-topdir')    { $jobs_topdir   = $ARGV[$iarg+1]; }
  $iarg++;
}
die("** Aborting [Undefined GENIE version. Use the --version option]")
unless defined $genie_version;
die("** Aborting [You need to specify which runs to submit. Use the --run option]")
unless defined $runnu;

$inuke_model    = "hA"                          unless defined $inuke_model;
$model_enum     = "01"                          unless defined $model_enum;
$nsubruns       = 1                             unless defined $nsubruns;
$use_valgrind   = 0                             unless defined $use_valgrind;
$arch           = "SL6.x86_64"                  unless defined $arch;
$production     = "$model_enum\_$genie_version" unless defined $production;
$cycle          = "01"                          unless defined $cycle;
$batch_system   = "PBS"                         unless defined $batch_system;
$queue          = "prod"                        unless defined $queue;
$softw_topdir   = "/opt/ppd/t2k/softw/GENIE/"   unless defined $softw_topdir;
$jobs_topdir    = "/opt/ppd/t2k/scratch/GENIE/" unless defined $jobs_topdir;  
$time_limit     = "60:00:00";
$genie_setup    = "$softw_topdir/generator/builds/$arch/$genie_version-setup";
$jobs_dir       = "$jobs_topdir/vld\_inuke-$production\_$cycle";
$mcseed         = 210921029;
$nev_per_subrun = 100000;

# inputs for event generation jobs
%evg_probepdg_hash = ( 
  '1002600870' => '211',
  '1202600870' => '-211',
  '1008200870' => '211',
  '1002602100' => '211',
  '1000601400' => '211',
  '1002901400' => '211',
  '1008201400' => '211',
  '1000600220' => '211',
  '1002800220' => '211',
  '1008200220' => '211',
  '1000600160' => '211',
  '1002800160' => '211',
  '1008200160' => '211',
  '1000600100' => '211',
  '1002800100' => '211',
  '1008200100' => '211',
  '1200600220' => '-211',
  '1202800220' => '-211',
  '1208200220' => '-211',
  '1200600160' => '-211',
  '1202800160' => '-211',
  '1208200160' => '-211',
  '1200600100' => '-211',
  '1202800100' => '-211',
  '1208200100' => '-211',
  '1000600500' => '211',
  '1200600500' => '211',
  '1208300500' => '211',
  '1000600300' => '211',
  '1000200300' => '211',
  '1000200220' => '211',
  '1000200160' => '211',
  '1000200100' => '211',
  '1000800114' => '211',
  '1000800163' => '211',
  '1000800240' => '211',
  '4000601400' => '211',
  '4002901400' => '211',
  '4008201400' => '211',
  '4000600800' => '2212',
  '4002000800' => '2212',
  '4000600730' => '2212',
  '4001300730' => '2212',
  '4002900730' => '2212',
  '4008200730' => '2212',
  '4000600597' => '2212',
  '4008200597' => '2212',
  '4000600800' => '2212',
  '4008200800' => '2212',
  '4002600558' => '2212',
  '4008200558' => '2212',
  '4000600300' => '2212',
  '4000600392' => '2212',
  '4002600065' => '2212',
  '4001300256' => '2212',
  '4008200256' => '2212',
  '4000600200' => '2212',
  '4000600197' => '2212',
  '4000600113' => '2212',
  '4002600113' => '2212',
  '4008200113' => '2212',
);
%evg_tgtpdg_hash = ( 
  '1002600870' => '1000260560',
  '1202600870' => '1000260560',
  '1008200870' => '1000822080',
  '1002602100' => '1000260560',
  '1000601400' => '1000060120',
  '1002901400' => '1000290630',
  '1008201400' => '1000822080',
  '1000600220' => '1000060120',
  '1002800220' => '1000280580',
  '1008200220' => '1000822080',
  '1000600160' => '1000060120',
  '1002800160' => '1000280580',
  '1008200160' => '1000822080',
  '1000600100' => '1000060120',
  '1002800100' => '1000280580',
  '1008200100' => '1000822080',
  '1200600220' => '1000060120',
  '1202800220' => '1000280580',
  '1208200220' => '1000822080',
  '1200600160' => '1000060120',
  '1202800160' => '1000280580',
  '1208200160' => '1000822080',
  '1200600100' => '1000060120',
  '1202800100' => '1000280580',
  '1208200100' => '1000822080',
  '1000600500' => '1000060120',
  '1200600500' => '1000060120',
  '1208300500' => '1000832090',
  '1000600300' => '1000060120',
  '1000200300' => '1000020040',
  '1000200220' => '1000020040',
  '1000200160' => '1000020040',
  '1000200100' => '1000020040',
  '1000800114' => '1000080160',
  '1000800163' => '1000080160',
  '1000800240' => '1000080160',
  '4000601400' => '1000060120',
  '4002901400' => '1000290630',
  '4008201400' => '1000822080',
  '4002000800' => '1000200400',
  '4000600800' => '1000060120',
  '4000600730' => '1000060120',
  '4001300730' => '1000130270',
  '4002900730' => '1000290630',
  '4008200730' => '1000822080',
  '4008200800' => '1000822080',
  '4000600597' => '1000060120',
  '4008200597' => '1000822080',
  '4002600558' => '1000260560',
  '4008200558' => '1000822080',
  '4000600300' => '1000060120',
  '4000600392' => '1000060120',
  '4002600065' => '1000260560',
  '4001300256' => '1000130270',
  '4008200256' => '1000822080',
  '4000600200' => '1000060120',
  '4000600197' => '1000060120',
  '4002600113' => '1000260560',
  '4008200113' => '1000822080',
);
%evg_kinetic_energy_hash = ( 
  '1002600870' => '0.870',
  '1202600870' => '0.870',
  '1008200870' => '0.870',
  '1002602100' => '2.100',
  '1000601400' => '1.400',
  '1002901400' => '1.400',
  '1008201400' => '1.400',
  '1000600220' => '0.220',
  '1002800220' => '0.220',
  '1008200220' => '0.220',
  '1000600160' => '0.160',
  '1002800160' => '0.160',
  '1008200160' => '0.160',
  '1000600100' => '0.100',
  '1002800100' => '0.100',
  '1008200100' => '0.100',
  '1202800220' => '0.220',
  '1208200220' => '0.220',
  '1200600160' => '0.160',
  '1202800160' => '0.160',
  '1208200160' => '0.160',
  '1200600100' => '0.100',
  '1202800100' => '0.100',
  '1208200100' => '0.100',
  '1000600500' => '0.500',
  '1200600500' => '0.500',
  '1208300500' => '0.500',
  '1000600300' => '0.300',
  '1000200300' => '0.300',
  '1000200220' => '0.220',
  '1000200160' => '0.160',
  '1000200100' => '0.100',
  '1000800114' => '0.114',
  '1000800163' => '0.163',
  '1000800240' => '0.240',
  '4000601400' => '1.400',
  '4002901400' => '1.400',
  '4008201400' => '1.400',
  '4000600800' => '0.800',
  '4002000800' => '0.800',
  '4000600730' => '0.730',
  '4001300730' => '0.730',
  '4002900730' => '0.730',
  '4008200730' => '0.730',
  '4008200800' => '0.800',
  '4000600597' => '0.597',
  '4008200597' => '0.597',
  '4002600558' => '0.558',
  '4008200558' => '0.558',
  '4000600300' => '0.300',
  '4000600392' => '0.392',
  '4002600065' => '0.065',
  '4001300256' => '0.256',
  '4008200256' => '0.256',
  '4000600800' => '0.200',
  '4000600197' => '0.197',
  '4000600113' => '0.113',
  '4002600113' => '0.113',
  '4008200113' => '0.113',
);
%vld_group_hash = ( 
  '1002600870' => 'iwamoto',
  '1202600870' => 'iwamoto',
  '1008200870' => 'iwamoto',
  '1002602100' => 'iwamoto',
  '1000601400' => 'shibata',
  '1002901400' => 'shibata',
  '1008201400' => 'shibata',
  '1000600220' => 'mckeown,levenson',
  '1002800220' => 'mckeown,levenson',
  '1008200220' => 'mckeown,levenson',
  '1000600160' => 'mckeown,levenson',
  '1002800160' => 'mckeown,levenson',
  '1008200160' => 'mckeown,levenson',
  '1000600100' => 'mckeown,levenson',
  '1002800100' => 'mckeown,levenson',
  '1008200100' => 'mckeown,levenson',
  '1000600220' => 'mckeown',
  '1002800220' => 'mckeown',
  '1008200220' => 'mckeown',
  '1000600160' => 'mckeown',
  '1002800160' => 'mckeown',
  '1008200160' => 'mckeown',
  '1000600100' => 'mckeown',
  '1002800100' => 'mckeown',
  '1008200100' => 'mckeown',
  '1000600500' => 'zumbro',
  '1200600500' => 'ouyang',
  '1200600500' => 'ouyang',
  '1000600300' => 'levenson',
  '1000200300' => 'mckeown,levenson',
  '1000200220' => 'mckeown,levenson',
  '1000200160' => 'mckeown,levenson',
  '1000200100' => 'mckeown,levenson',
  '1000800114' => 'ingram',
  '1000800163' => 'ingram',
  '1000800240' => 'ingram',
  '4000601400' => 'shibata',
  '4002901400' => 'shibata',
  '4008201400' => 'shibata',
  '4000600800' => 'mcgill,amian',
  '4002000800' => 'mcgill',
  '4000600730' => 'cochran',
  '4001300730' => 'cochran',
  '4002900730' => 'cochran',
  '4008200730' => 'cochran',
  '4008200800' => 'amian',
  '4000600597' => 'amian',
  '4008200597' => 'amian',
  '4002600558' => 'beck',
  '4008200558' => 'beck',
  '4000600300' => 'kin',
  '4000600392' => 'kin',
  '4002600065' => 'bertrand',
  '4001300256' => 'stamer',
  '4008200256' => 'stamer',
  '4000600200' => 'carman',
  '4000600197' => 'hautala',
  '4000600113' => 'meier',
  '4002600113' => 'meier',
  '4008200113' => 'meier',
);

# make the jobs directory
#
mkpath ($jobs_dir, {verbose => 1, mode=>0777});

#
# submit event generation jobs
#

# run loop
for my $curr_runnu (keys %evg_probepdg_hash)  {

 #
 # get runnu-dependent info
 #
 $probe   = $evg_probepdg_hash       {$curr_runnu};
 $tgt     = $evg_tgtpdg_hash         {$curr_runnu};
 $ke      = $evg_kinetic_energy_hash {$curr_runnu};
 $vldgrp  = $vld_group_hash          {$curr_runnu};

 # check whether to commit current run  
 print "checking whether to submit run: $curr_runnu \n";

 $do_submit = 
    ( $runnu=~m/$curr_runnu/                 ) || 
    ( $runnu eq "all"                        ) || 
    ( $vldgrp=~m/$runnu/                     ) ||
    ( $probe eq   '22' && $runnu=~m/gamma/   ) ||
    ( $probe eq  '211' && $runnu=~m/piplus/  ) ||
    ( $probe eq '-211' && $runnu=~m/piminus/ ) ||
    ( $probe eq  '111' && $runnu=~m/pi0/     ) ||
    ( $probe eq  '311' && $runnu=~m/Kplus/   ) ||
    ( $probe eq '-311' && $runnu=~m/Kminus/  ) ||
    ( $probe eq '2212' && $runnu=~m/proton/  ) ||
    ( $probe eq '2112' && $runnu=~m/neutron/ );

 if($do_submit) {
    print "** submitting event generation run: $curr_runnu \n";

    # submit subruns
    for($isubrun = 0; $isubrun < $nsubruns; $isubrun++) {

       # Run number key: PPTTTEEEEEMMxx
       $curr_subrunnu     = 10000 * $curr_runnu + 100 * $model_enum + $isubrun;
#      $grep_pipe         = "grep -B 20 -A 30 -i \"warn\\|error\\|fatal\"";
       $grep_pipe         = "grep -B 20 -A 30 -i fatal";
       $filename_template = "$jobs_dir/inuke-$inuke_model-$curr_subrunnu";
       $curr_seed         = $mcseed + $isubrun;
       $valgrind_cmd      = "valgrind --tool=memcheck --error-limit=no --leak-check=yes --show-reachable=yes";
       $gntp_prefix       = "gntp.$inuke_model";
       $evgen_cmd         = "gevgen_hadron -n $nev_per_subrun -m $inuke_model -k $ke -p $probe -t $tgt -r $curr_subrunnu --seed $curr_seed -o $gntp_prefix --message-thresholds Messenger_laconic.xml";
       $conv_cmd          = "gntpc -f ginuke -i $gntp_prefix.$curr_subrunnu.ghep.root --message-thresholds Messenger_laconic.xml";
       $evgen_cmd         = "$evgen_cmd | $grep_pipe &> $filename_template.evgen.log";
       $conv_cmd          = "$conv_cmd  | $grep_pipe &> $filename_template.conv.log ";

       print "@@ exec: $evgen_cmd \n";
       print "@@ exec: $conv_cmd  \n\n";

       #
       # submit
       #
  
       # PBS case
       if($batch_system eq 'PBS' || $batch_system eq 'HTCondor_PBS') {
          $batch_script  = "$filename_template.pbs";
          open(PBS, ">$batch_script") or die("Can not create the PBS batch script");
          print PBS "#!/bin/bash \n";
          print PBS "#PBS -N inuke-$curr_subrunnu \n";
          print PBS "#PBS -l cput=$time_limit \n";
          print PBS "#PBS -o $filename_template.pbsout.log \n";
          print PBS "#PBS -e $filename_template.pbserr.log \n";
          print PBS "source $genie_setup \n"; 
          print PBS "cd $jobs_dir \n";
          print PBS "$evgen_cmd \n";
          print PBS "$conv_cmd \n";
          close(PBS);
          $job_submission_command = "qsub";
          if($batch_system eq 'HTCondor_PBS') {
              $job_submission_command = "condor_qsub";
          }
          `$job_submission_command -q $queue $batch_script`;
       } # PBS

       # LSF case
       if($batch_system eq 'LSF') {
          $batch_script  = "$filename_template.sh";
          open(LSF, ">$batch_script") or die("Can not create the LSF batch script");
          print LSF "#!/bin/bash \n";
          print LSF "#BSUB-j inuke-$curr_subrunnu \n";
          print LSF "#BSUB-q $queue \n";
          print LSF "#BSUB-c $time_limit \n";
          print LSF "#BSUB-o $filename_template.lsfout.log \n";
          print LSF "#BSUB-e $filename_template.lsferr.log \n";
          print LSF "source $genie_setup \n"; 
          print LSF "cd $jobs_dir \n";
          print LSF "$evgen_cmd \n";
          print LSF "$conv_cmd \n";
          close(LSF);
          `bsub < $batch_script`;
       } # LSF

       # HTCondor
       if($batch_system eq 'HTCondor') {
           $batch_script = "$filename_template.htc";
           open(HTC, ">$batch_script") or die("Can not create the Condor submit description file: $batch_script");
           print HTC "Universe               = vanilla \n";
           print HTC "Executable             = $softw_topdir/generator/builds/$arch/$genie_version/src/scripts/production/batch/htcondor_exec.sh \n";
           print HTC "Arguments              = $genie_setup $jobs_dir $evgen_cmd $conv_cmd\n";
           print HTC "Log                    = $filename_template.log \n";
           print HTC "Output                 = $filename_template.out \n";
           print HTC "Error                  = $filename_template.err \n";
           print HTC "Request_memory         = 2 GB \n";
           print HTC "Queue \n";
           close(HTC);
          `condor_submit $batch_script`;
       } #HTCondor
    
       # slurm case
       if($batch_system eq 'slurm') {
          my $time_lim = `sinfo -h -p batch -o %l`;
          my ($days, $hours, $remainder) = $time_lim =~ /([0]+)-([0-9]+):(.*)/;
          my $newhours = $days * 24 + $hours;
          my $new_time_lim = "$newhours:$remainder";
          $time_limit = $new_time_lim lt $time_limit ? $new_time_lim : $time_limit;
          $batch_script  = "$filename_template.sh";
          open(SLURM, ">$batch_script") or die("Can not create the SLURM batch script");
          print SLURM "#!/bin/bash \n";
          print SLURM "#SBATCH-p $queue \n";
          print SLURM "#SBATCH-o $filename_template.lsfout.log \n";
          print SLURM "#SBATCH-e $filename_template.lsferr.log \n";
          print SLURM "#SBATCH-t $time_limit \n";
          print SLURM "source $genie_setup \n"; 
          print SLURM "cd $jobs_dir \n";
          print SLURM "$evgen_cmd \n";
          print SLURM "$conv_cmd \n";
          close(SLURM);
          `sbatch --job-name=inuke-$curr_subrunnu $batch_script`;
       } # slurm

       # no batch system, run jobs interactively
       if($batch_system eq 'none') {
          system("source $genie_setup; cd $jobs_dir; $evgen_cmd; $conv_cmd");
       } # interactive mode

    } # loop over subruns

 } #checking whether to submit current run
} # loop over runs

