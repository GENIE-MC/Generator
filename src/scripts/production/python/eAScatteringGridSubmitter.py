#! /usr/bin/env python

"""\
This script is responsible to submit jobs to the grid
The script generates splines and events in parallel 

Author: 
      Julia Tena Vidal <jtenavidal \st tauex.tau.ac.il>
      Tel Aviv University
Copyright:
   Copyright (c) 2003-2022, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org
"""
import os, optparse, glob, tarfile

import sys
sys.path.insert(1, 'xsec_splines/')


import vNSplineCommands as vN
import vASplineCommands as vA
import GroupSplineCommands as group
import FNALGridUtils as FNAL

sys.path.insert(1, 'event_generation')
import eScatteringGenCommands as eA

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
op.add_option("--source-prod-dir", dest="MotherDir", default='', help="Jobs topdir used as a source for missing xsec splines.")
op.add_option("--total-xsec", dest="XSEC", default='', help="Total cross section file to use for event production")
op.add_option("--config-dir", dest="CONF", default='', help="Path to GENIE config dir")
op.add_option("--probe-list", dest="PROBELIST", default='11', help = "Comma separated list of lepton flavour (neutrino and electrons are handled). Default: %default.") 
op.add_option("--nu-tgt-list", dest="NUTGTLIST", default='all', help = "Comma separated list of Targets. Default: %default.") 
op.add_option("--e-tgt-list", dest="ETGTLIST", default='all', help = "Comma separated list of Targets. Default: %default.") 
op.add_option("--vN-gen-list", dest="vNList", default='all', help="Comma separated list of event generator list used for the free nucleon spline generation. Can be used to specify electron procecess as well")
op.add_option("--vA-gen-list", dest="vAList", default='all', help="Comma separated list of event generator list used for the nuclei spline generation.  Can be used to specify electron procecess as well")
op.add_option("--E-nu-maxsplines", dest="NuEMAXSPLINE", type="int", default=100, help="Maximum energy for the splines in GeV. Default: %default ")
op.add_option("--E-e-maxsplines", dest="EEMAXSPLINE", type="int", default=30, help="Maximum energy for the splines in GeV. Default: %default ")
op.add_option("--N-nu-knots", dest="NuKnots", type="int", default=100, help="Number of knots per neutrino spline. Default: %default")
op.add_option("--N-e-knots", dest="EKnots", type="int", default=30, help="Number of knots per electron spline. Default: %default")
op.add_option("--event-generator-list", dest="EvGenList", default='EM', help="Event generator list to be used for event generation. Default EM")
op.add_option("--nu-ntotevents", dest="NuEvents", type="int", default=10000, help="Number of total events, default: 100 k")
op.add_option("--e-ntotevents", dest="EEvents", type="int", default=100000, help="Number of total events, default: 100 k")
op.add_option("--nmaxevents",dest="NMax", type="int", default=400000,help="Max number of events to run per event generation, default %default")
op.add_option("--ebeam-energy", dest="EnergyBeam", default="2", help="Comma separated list of beam energy for electrons. Default %default GeV")
op.add_option("--gst-output", dest="GSTOutput", default=False, action="store_true",help="Store gst root file.")
op.add_option("--no-ghep-output", dest="NoGHEPOutput", default=False, action="store_true",help="GHEP GENIE files is removed to reduce memory.")
op.add_option("--starting-point", dest="start_ID", type="int", default=0, help="0 -> Free nucleon splines, 1 -> combine free nucl splines, 2 -> Compound nuclei splines, 3 -> Combine compound nuclei splines, 4 -> Event Production")
op.add_option("--stopping-point", dest="end_ID", type="int", default=9999, help="Numbers as above, Default: 9999") 
op.add_option("--tune", dest="TUNE", default="G18_02a_02_11b", help="Tune to be compared against data (default: %default)")
op.add_option("--submit-jobs", dest="SUBMIT", default=False, action="store_true", help="Generate configuration and submit to grid" )
op.add_option("--job-lifetime", dest="JOBLIFE", default=25, help="Expected lifetime on the grid for all the jobs to be finished, default %default h")
op.add_option("--job-lifetime-vN", dest="vNJOBLIFE", default=15, help="Expected lifetime on the grid for all the vN spline jobs to be finished, default %default h")
op.add_option("--job-lifetime-vA", dest="vAJOBLIFE", default=3, help="Expected lifetime on the grid for all the vA spline jobs to be finished, default %default h")
op.add_option("--job-lifetime-generation", dest="GENJOBLIFE", default=4, help="Expected lifetime on the grid for all the event generation jobs to be finished, default %default h")
op.add_option("--job-lifetime-group", dest="GROUPJOBLIFE", default=1, help="Expected lifetime on the grid for all the grouping jobs to be finished, default 1h")
op.add_option("--job-memory-generation", dest="GENJOBMEM", default=2, help="Expected memory usage for each event generation job, default %default GB")
op.add_option("--job-disk-generation", dest="GENJOBDISK", default=2000, help="Expected disk space for each event generation job, default %default MB")
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

if opts.CONF : 
    if not os.path.exists(opts.CONF) : 
        print ( " GENIE Configuraion dir specified does not exist: " + opts.CONF + " . Abort ..." ) 
        exit()

    print( 'Using configuration files from ' + opts.CONF + ' ...' )


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
vNdir = opts.JOBSTD+'/'+version+'-'+opts.PROD+'_'+opts.CYCLE+'-xsec_vN/'
vNsplines = vNdir+'total_xsec.xml'
vAdir = opts.JOBSTD+'/'+version+'-'+opts.PROD+'_'+opts.CYCLE+'-xsec_vA/'
vAsplines = vAdir+'total_xsec.xml'

if opts.XSEC : 
    if os.path.isfile(opts.XSEC) == False :
        print(" Input XSec file doesn't exist ")
        exit() 
    vAsplines = opts.XSEC 

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

# Store commands with ID :
command_dict = {}

# Get correct ID as requested by user:
loop_start = 0 
loop_end = 5
if opts.start_ID != 0 :
    loop_start = int(opts.start_ID)

if loop_end > int(opts.end_ID) : 
    loop_end= int(opts.end_ID)

total_time = 0 ; 
loop_i = loop_start
while loop_i < loop_end + 1: 
    # ID = 0 # vN splines
    if loop_i == 0 :
        command_dict.update( vN.vNSplineCommands(opts.PROBELIST,opts.vNList,opts.NuEMAXSPLINE,opts.EEMAXSPLINE,opts.NuKnots,opts.EKnots,opts.TUNE,version,opts.GRID,opts.GROUP,opts.CONF,opts.ARCH,opts.PROD,opts.CYCLE,opts.SOFTW,opts.GENIE,opts.JOBSTD,grid_setup,genie_setup,opts.vNJOBLIFE,opts.BRANCH,opts.GIT_LOCATION) )
        total_time += int(opts.vNJOBLIFE) 

    # ID = 1 # group vN splines
    if loop_i == 1 :
        vNMotherDir = ''
        if opts.MotherDir !='' : 
            vNMotherDir = opts.MotherDir+'/'+version+'-'+opts.PROD+'_'+opts.CYCLE+'-xsec_vN/'
        
        command_dict.update( group.GroupSplineCommands( True,vNdir,vNMotherDir,opts.TUNE,opts.vNList,version,opts.CONF,opts.GRID,opts.GROUP,opts.ARCH,opts.PROD,opts.CYCLE,opts.SOFTW,opts.GENIE,grid_setup,genie_setup,opts.JOBSTD,False, False, opts.GROUPJOBLIFE,opts.BRANCH,opts.GIT_LOCATION ) )
        total_time += int(opts.GROUPJOBLIFE)
 
    if loop_i == 2 : 
        # ID = 2 # vA splines
        command_dict.update( vA.vASplineCommands(opts.PROBELIST,opts.NUTGTLIST,opts.ETGTLIST,opts.vAList,opts.NuEMAXSPLINE,opts.EEMAXSPLINE,opts.NuKnots,opts.EKnots,opts.TUNE,vNsplines,version,opts.GRID,opts.GROUP,opts.CONF,opts.ARCH,opts.PROD,opts.CYCLE,opts.SOFTW,opts.GENIE,opts.JOBSTD,grid_setup,genie_setup,opts.vAJOBLIFE,opts.BRANCH,opts.GIT_LOCATION) )
        total_time += int(opts.vAJOBLIFE)

    if loop_i == 3 : 
        # ID = 3 # Group vA splines
        vAMotherDir = ''
        if opts.MotherDir !='' : 
            vAMotherDir = opts.MotherDir+'/'+version+'-'+opts.PROD+'_'+opts.CYCLE+'-xsec_vA/'
        command_dict.update( group.GroupSplineCommands( False,vAdir,vAMotherDir,opts.TUNE,opts.vAList,version,opts.CONF,opts.GRID,opts.GROUP,opts.ARCH,opts.PROD,opts.CYCLE,opts.SOFTW,opts.GENIE,grid_setup,genie_setup,opts.JOBSTD,False, False,opts.GROUPJOBLIFE,opts.BRANCH,opts.GIT_LOCATION ) )
        total_time += int(opts.GROUPJOBLIFE) 

    if loop_i == 4 : 
        # ID = 4 # Event generation commands
        command_dict.update( eA.eScatteringGenCommands(opts.PROBELIST,opts.ETGTLIST,opts.EnergyBeam,vAsplines,opts.EEvents,opts.TUNE, opts.EvGenList, opts.NMax, opts.GSTOutput, opts.NoGHEPOutput,version,opts.CONF, opts.ARCH, opts.PROD, opts.CYCLE,opts.GRID, opts.GROUP,opts.SOFTW,opts.GENIE,opts.JOBSTD,grid_setup,genie_setup,opts.GENJOBLIFE,opts.GENJOBMEM,opts.GENJOBDISK,opts.BRANCH,opts.GIT_LOCATION) )
        total_time += int(opts.GENJOBLIFE)
    
    loop_i += 1 

if total_time > int(opts.JOBLIFE) : 
    print ( "Total time of subjobs requested ("+str(total_time)+") is bigger than the job's expected time ("+str(opts.JOBLIFE)+") ... Abort ..." ) 
    exit() 

if opts.GRID is 'FNAL':
    if total_time > 96 or int(opts.JOBLIFE) > 96 : 
        print ( "Total time at the grid cannot exceed 96h ")
        exit() 

    # Write xml file
    grid_name = FNAL.WriteXMLFile(command_dict, loop_start, loop_end, opts.JOBSTD)

    main_sub_name = FNAL.WriteMainSubmissionFile(opts.JOBSTD, opts.GENIE, opts.GROUP, grid_setup, genie_setup, grid_name, opts.JOBLIFE)

if opts.SUBMIT == True: 
    # SUBMIT JOB
    print ( "Submitting jobs to grid ... \n")
    os. system("source "+main_sub_name)
else : 
    print ( "In testing mode - jobs won't be submitted\n" )
