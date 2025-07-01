#! /usr/bin/env python

"""\
This script is responsible to submit jobs to the grid
The script generates splines and events in parallel 

Author: 
      Julia Tena Vidal <jtenavidal \st tauex.tau.ac.il>
      Tel Aviv University
Copyright:
   Copyright (c) 2003-2025, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org
"""
import os, optparse, glob, tarfile

import sys
sys.path.insert(1, 'xsec_splines/')

import FNALGridUtils as FNAL

sys.path.insert(1, 'event_generation')
import hadronScatteringGenCommands as hadronA

op = optparse.OptionParser(usage=__doc__)
op.add_option("--version", dest="VERSION", default="master", help="Genie version. Default: %default")
op.add_option("--git-location", dest="GIT_LOCATION", default="https://github.com/GENIE-MC/Generator", help="Github location from where to get the GENIE Generator code. Defaulted to %default")
op.add_option("--git-branch", dest="BRANCH", default="master", help="Genie version branch name. Default: %default")
op.add_option("--cycle", dest="CYCLE", default="01", help="Cycle (default: %default)")
op.add_option("--arch", dest="ARCH", default='SL6.x86_64', help="arch number, default: %default")
op.add_option("--production", dest="PROD", default="routine_validation", help="Production (default: %default)")
op.add_option("--grid-system", dest="GRID", default="FNAL", help="Grid system is %default ")
op.add_option("--grid-group", dest="GROUP", default="genie", help="Grid group to run your jobs. I.E: FNAL groups are genie, dune, etc..")
op.add_option("--softw-topdir", dest="SOFTW", default=os.getenv('GENIE_MASTER_DIR'), help = "software top dir") 
op.add_option("--genie-topdir", dest="GENIE", default=os.getenv('GENIE'), help = "GENIE topdir: %default")
op.add_option("--jobs-topdir", dest="JOBSTD", default=os.getenv('PWD'), help="Top level dir for the job files (default: %default)")
op.add_option("--config-dir", dest="CONF", default='', help="Path to GENIE config dir")
op.add_option("--hadron-list", dest="PROBELIST", default='211', help = "Comma separated list of hadrons. Default: %default.") 
op.add_option("--tgt-list", dest="TGTLIST", default='all', help = "Comma separated list of Targets. Default: %default.") 
op.add_option("--ntotevents", dest="Events", type="int", default=100000, help="Number of total events, default: 100 k")
op.add_option("--nmaxevents",dest="NMax", type="int", default=400000,help="Max number of events to run per event generation, default %default")
op.add_option("--hadron-KE", dest="HadronKE", default="2", help="Comma separated list of hadron KineticEnergy. Default %default GeV")
op.add_option("--hadron-KE-min", dest="HadronKEMin", default="0", help="Hadron KineticEnergy minimum - used with Flux. Default %default GeV")
op.add_option("--hadron-KE-max", dest="HadronKEMax", default="0", help="Hadron KineticEnergy maximum - used with Flux. Default %default GeV")
op.add_option("--flux", dest="FLUX", default="none", help="Default none. To use an input root file, specify the location of the file and the TH1D name as follows: file.root,th1d_name. uniform is an option too.")
op.add_option("--minenergy-fluxrange", dest="MinEnergyFlux", default="0.1", help="Minimum electron energy. Default %default")
op.add_option("--maxenergy-fluxrange", dest="MaxEnergyFlux", default="100", help="Maximum electron energy. Default %default.")
op.add_option("--seed", dest="Seed", default=210921029, help="Set stargint point seed. Default %default")
op.add_option("--gen-runid", dest="RunID", default=0, help="Set Starting run id. Default %default")
op.add_option("--ginuke-output", dest="GINUKEOutput", default=False, action="store_true",help="Store ginuke root file.")
op.add_option("--no-ghep-output", dest="NoGHEPOutput", default=False, action="store_true",help="GHEP GENIE files is removed to reduce memory.")
op.add_option("--tune", dest="TUNE", default="G18_02a_02_11b", help="Tune to be compared against data (default: %default)")
op.add_option("--submit-jobs", dest="SUBMIT", default=False, action="store_true", help="Generate configuration and submit to grid" )
op.add_option("--job-lifetime", dest="JOBLIFE", default=15, help="Expected lifetime on the grid, default %default h")
op.add_option("--subjob-memory", dest="JOBMEM", default="2GB", help="Expected memory usage for each event generation job, default %default. Notice you must specify the units. ")
op.add_option("--subjob-disk", dest="JOBDISK", default="2GB", help="Expected disk space for each event generation job, default %default. Notice you must specify the units. ")
op.add_option("--mainjob-lifetime", dest="MAINLIFE", default=25, help="Expected lifetime on the grid for all the jobs to be finished, default %default h")
op.add_option("--mainjob-memory", dest="MAINMEM", default="20GB", help="Main job memory usage, default %default")
op.add_option("--mainjob-disk", dest="MAINDISK", default="20GB", help="Main job disk usate, default %default")
op.add_option("--store-comitinfo", dest="STORECOMMIT", default=False, action="store_true", help="Store command line in jobstopdir directory")
opts, args = op.parse_args()

if( opts.STORECOMMIT ) :
    output_file = opts.JOBSTD+"/input_options.txt"
    input_names = []
    input_variables = []
    for opt, value in opts.__dict__.items():
        input_variables.append( value )
        input_names.append(opt)

    if os.path.exists(output_file) :
        os.remove(output_file)

    with open(output_file,'w') as f:  
        f.write( "##################################################################################################")
        f.write( "# This document contains the input variables used to run the GENIE jobs stored in this directory #")
        f.write( "##################################################################################################")
        for i in range(len(input_names)) :
            f.write( str(input_names[i]) + " " + str(input_variables[i]) + "\n" )

# Print information
print ("Creating job substructure and submission scripts... \n")

# Check jobstopdir exists
if not os.path.exists(opts.JOBSTD) :
    print ( "Jobs top dir path "+opts.JOBSTD+" doesn't exist. Abort...\n")
    exit()

# Check jobstopdir exists
if not os.path.exists(opts.GENIE) :
    print ( "GENIE Path "+opts.GENIE+" doesn't exist. Abort...\n")
    exit()

if opts.GRID != 'FNAL':
    print ("For now only FNAL option is available. Abort... \n")
    exit()
    # In the future we might have other grids

# Check we run from pnfs persistent or scratch for FNAL:
if opts.GRID == 'FNAL':
    if 'pnfs' not in opts.JOBSTD : 
        print ("Not runing from pnfs:"+opts.JOBSTD+" . Jobs top dir must be in pnfs for the submission scrpits to work. Abort ...")
        exit()

message_thresholds = "$GALGCONF/Messenger.xml"
if opts.CONF : 
    if not os.path.exists(opts.CONF) : 
        print ( " GENIE Configuraion dir specified does not exist: " + opts.CONF + " . Abort ..." ) 
        exit()

    print( 'Using configuration files from ' + opts.CONF + ' ...' )
    # Check if we have a message thresholds file 
    if os.path.exists(opts.CONF+"/Messenger.xml"):
        message_thresholds = "$CONDOR_DIR_INPUT/conf/Messenger.xml"


#JobSub is made available through the UPS package jobsub_client
os.system("source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup" ) 

# Check version is not a path
temp_version = opts.VERSION.split('/') 
version = opts.VERSION

if len(temp_version ) : 
    version = temp_version[len(temp_version)-1]
    print ( ' Setting version to '+ version ) 

if opts.BRANCH: 
    print( ' Cloning GENIE Generator ' + opts.BRANCH ) 

# Define directories:
hNdir = opts.JOBSTD+'/'+version+'-'+opts.PROD+'_'+opts.CYCLE+'-hadroN/'

# configure setup 
if opts.GRID == 'FNAL' : 
    setup_file = opts.GENIE+'/src/scripts/production/python/setup_FNAL.sh'
    GENIE_setup_file = opts.GENIE+'/src/scripts/production/python/setup_GENIE.sh'
    pnfs_setup = opts.JOBSTD+"/setup_FNAL.sh"
    GENIE_pnfs_setup = opts.JOBSTD+"/setup_GENIE.sh"
    if os.path.exists(pnfs_setup) : 
        os.remove(pnfs_setup)
    if os.path.exists(GENIE_pnfs_setup) : 
        os.remove(GENIE_pnfs_setup)
    os.system('cp '+setup_file+' '+opts.JOBSTD )
    os.system('cp '+GENIE_setup_file+' '+opts.JOBSTD )
    grid_setup = opts.JOBSTD+"/setup_FNAL.sh"
    genie_setup = opts.JOBSTD+"/setup_GENIE.sh"
else : 
    genie_setup = opts.SOFTW+'/generator/builds/'+arch+'/'+version+'-setup'
    grid_setup = ""

# Check whether INCL/G4 have to be configured. Also keeping track of hA and hN to call the correct model in the production
configure_hA = False
configure_hN = False
configure_INCL = False
configure_G4 = False 
if "a" in opts.TUNE :
    configure_hA = True
elif "b" in opts.TUNE:
    configure_hN = True
elif "c" in opts.TUNE:
    configure_INCL = True
elif "d" in opts.TUNE:
    configure_G4 = True 

# Store commands with ID :
command_dict = {}

command_dict.update( hadronA.hadronScatteringGenCommands(opts.PROBELIST,opts.TGTLIST,opts.HadronKE,opts.Events,opts.HadronKEMin,opts.HadronKEMax,
                                                         opts.FLUX,opts.TUNE,opts.NMax,opts.Seed, opts.GINUKEOutput, opts.NoGHEPOutput,version,
                                                         opts.CONF, opts.ARCH, opts.PROD, opts.CYCLE,opts.GRID, opts.GROUP,opts.SOFTW,opts.GENIE,
                                                         opts.JOBSTD,grid_setup,genie_setup,message_thresholds,opts.JOBLIFE,opts.JOBMEM,opts.JOBDISK,opts.BRANCH,opts.GIT_LOCATION,
                                                         configure_hA, configure_hN, configure_INCL, configure_G4) )

if opts.GRID is 'FNAL':
    # Write xml file
    grid_name = FNAL.WriteXMLFile(command_dict, 0, 1, opts.JOBSTD)

    main_sub_name = FNAL.WriteMainSubmissionFile(opts.JOBSTD, opts.GENIE, opts.GROUP, grid_setup, genie_setup, grid_name, opts.MAINLIFE, opts.MAINMEM, opts.MAINDISK )

if opts.SUBMIT == True: 
    # SUBMIT JOB
    print ( "Submitting jobs to grid ... \n")
    os. system("source "+main_sub_name)
else : 
    print ( "In testing mode - jobs won't be submitted\n" )
