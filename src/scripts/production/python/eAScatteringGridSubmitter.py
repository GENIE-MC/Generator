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
import os, optparse

import sys
sys.path.insert(1, 'xsec_splines/')

import vNSplineCommands as vN
import vASplineCommands as vA
import GroupSplineCommands as group
import GridUtils as utils

sys.path.insert(1, 'event_generation')
import eScatteringGenCommands as eA

op = optparse.OptionParser(usage=__doc__)
op.add_option("--version", dest="VERSION", default="master", help="Genie version number (default: %default)")
op.add_option("--cycle", dest="CYCLE", default="01", help="Cycle (default: %default)")
op.add_option("--arch", dest="ARCH", default='SL6.x86_64', help="arch number, default: %default")
op.add_option("--production", dest="PROD", default="routine_validation", help="Production (default: %default)")
op.add_option("--grid-system", dest="GRID", default="FNAL", help="Grid system is %default ")
op.add_option("--grid-group", dest="GROUP", default="genie", help="Grid group to run your jobs. I.E: FNAL groups are genie, dune, etc..")
op.add_option("--softw-topdir", dest="SOFTW", default=os.getenv('GENIE_MASTER_DIR'), help = "software top dir") 
op.add_option("--genie-topdir", dest="GENIE", default=os.getenv('GENIE'), help = "GENIE topdir: %default")
op.add_option("--jobs-topdir", dest="JOBSTD", default=os.getenv('PWD'), help="Top level dir for the job files (default: %default)")
op.add_option("--source-prod-dir", dest="MotherDir", default='', help="Jobs topdir used as a source for missing xsec splines.")
op.add_option("--config-dir", dest="CONF", default='', help="Path to GENIE config dir")
op.add_option("--nu-list", dest="NULIST", default='11', help = "Comma separated list of lepton flavour (neutrino and electrons are handled). Default: %default.") 
op.add_option("--tgt-list", dest="TGTLIST", default='all', help = "Comma separated list of Targets. Default: %default.") 
op.add_option("--gen-list", dest="GenList", default='all', help = "Comma separated list of event generator list to be used for splines, default all") 
op.add_option("--e-max", dest="EMAX", default=30, help="Maximum energy for the splines in GeV. Default: %default ")
op.add_option("--n-knots", dest="Knots", default=100, help="Number of knots per spline. Default: %default")
op.add_option("--event-generator-list", dest="EvGenList", default='all', help="Event generator list to be used for event generation. Default all")
op.add_option("--ntotevents", dest="NEvents", default=100000, help="Number of total events, default: 100 k")
op.add_option("--nmaxevents",dest="NMax", default=100000,help="Max number of events to run per event generation, default 100k")
op.add_option("--energy", dest="Energy", default="2", help="Comma separated list of beam energy for electrons. Default %default GeV")
op.add_option("--starting-point", dest="start_ID", default=0, help="0 -> Free nucleon splines, 1 -> combine free nucl splines, 2 -> Compound nuclei splines, 3 -> Combine compound nuclei splines, 4 -> Event Production")
op.add_option("--stopping-point", dest="end_ID", default=9999, help="Numbers as above, Default: 9999") 
op.add_option("--tune", dest="TUNE", default="G18_02a_02_11b", help="Tune to be compared against data (default: %default)")
op.add_option("--submit-jobs", dest="SUBMIT", default=False, action="store_true", help="Generate configuration and submit to grid" )
opts, args = op.parse_args()

op.add_option("--vN-gen-list", dest="vNList", default=opts.GenList, help="Comma separated list of event generator list used for the free nucleon spline generation. Default all or gen-list if defined")
op.add_option("--vA-gen-list", dest="vAList", default=opts.GenList, help="Comma separated list of event generator list used for the nuclei spline generation. Default all or gen-list if defined")
opts, args = op.parse_args()

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

# Define directories:
vNdir = opts.JOBSTD+'/'+opts.VERSION+'-'+opts.PROD+'_'+opts.CYCLE+'-xsec_vN/'
vNsplines = vNdir+'total.xsec'
vAdir = opts.JOBSTD+'/'+opts.VERSION+'-'+opts.PROD+'_'+opts.CYCLE+'-xsec_vA/'
vAsplines = vAdir+'total.xsec'

# configure setup 
if opts.GRID == 'FNAL' : 
    setup_file = opts.GENIE+'/src/scripts/production/python/setup_FNALGrid.sh'
    pnfs_setup = opts.JOBSTD+"/setup_FNALGrid.sh"
    if os.path.exists(pnfs_setup) : 
        os.remove(pnfs_setup)
    os.system('cp '+setup_file+' '+opts.JOBSTD )
    genie_setup = opts.JOBSTD+"/setup_FNALGrid.sh"
else : 
    genie_setup = opts.SOFTW+'/generator/builds/'+arch+'/'+version+'-setup'

# Store commands with ID :
command_dict = {}
# ID = 0 # vN splines
command_dict.update( vN.vNSplineCommands(opts.NULIST,opts.vNList,opts.EMAX,opts.Knots,opts.TUNE,opts.VERSION,opts.GRID,opts.GROUP,opts.CONF,opts.ARCH,opts.PROD,opts.CYCLE,opts.SOFTW,opts.GENIE,opts.JOBSTD,genie_setup) )

# ID = 1 # group vN splines
vNMotherDir = ''
if opts.MotherDir !='' : 
    vNMotherDir = opts.MotherDir+'/'+opts.VERSION+'-'+opts.PROD+'_'+opts.CYCLE+'-xsec_vN/'

command_dict.update( group.GroupSplineCommands( True,vNMotherDir,opts.TUNE,opts.VERSION,opts.CONF,opts.GRID,opts.GROUP,opts.ARCH,opts.PROD,opts.CYCLE,opts.SOFTW,opts.GENIE,genie_setup,opts.JOBSTD,vNdir,False, False ) )# THE LAST TWO TO BE CONFIGURED

# ID = 2 # vA splines
command_dict.update( vA.vASplineCommands(opts.NULIST,opts.TGTLIST,opts.vAList,opts.EMAX,opts.Knots,opts.TUNE,vNsplines,opts.VERSION,opts.GRID,opts.GROUP,opts.CONF,opts.ARCH,opts.PROD,opts.CYCLE,opts.SOFTW,opts.GENIE,opts.JOBSTD,genie_setup) )

# ID = 3 # Group vA splines
vAMotherDir = ''
if opts.MotherDir !='' : 
    vAMotherDir = opts.MotherDir+'/'+opts.VERSION+'-'+opts.PROD+'_'+opts.CYCLE+'-xsec_vA/'

command_dict.update( group.GroupSplineCommands( False,vAMotherDir,opts.TUNE,opts.VERSION,opts.CONF,opts.GRID,opts.GROUP,opts.ARCH,opts.PROD,opts.CYCLE,opts.SOFTW,opts.GENIE,genie_setup,opts.JOBSTD,vAdir,False, False ) )# THE LAST TWO TO BE CONFIGURED

# ID = 4 # Event generation commands
command_dict.update( eA.eScatteringGenCommands(opts.NULIST,opts.TGTLIST,opts.Energy,vAsplines,opts.NEvents,opts.TUNE, opts.EvGenList, opts.NMax, opts.VERSION, opts.CONF, opts.ARCH, opts.PROD, opts.CYCLE,opts.GRID, opts.GROUP,opts.SOFTW,opts.GENIE,opts.JOBSTD,genie_setup) )

# Get correct ID as requested by user:
loop_start = 0 
loop_end = command_dict.keys()[-1]
                    
if opts.start_ID != 0 :
    loop_start = opts.start_ID

if opts.end_ID < loop_end : 
    loop_end= opts.end_ID 

grid_name = utils.WriteXMLFile(command_dict, loop_start, loop_end, opts.JOBSTD)

main_sub_name = utils.WriteMainSubmissionFile(opts.JOBSTD, opts.GENIE, opts.GRID, opts.GROUP, genie_setup, grid_name )

if opts.SUBMIT == True: 
    # SUBMIT JOB
    print ( "Submitting jobs to grid ... \n")
    os. system("source "+main_sub_name)
else : 
    print ( "In testing mode - jobs won't be submitted\n" )
