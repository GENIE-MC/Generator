#! /usr/bin/env python

"""\
This script generates e-scattering events 

Author: 
      Julia Tena Vidal <jtenavidal \st tauex.tau.ac.il>
      Tel Aviv University
Copyright:
   Copyright (c) 2003-2022, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org
"""
import os 
import GridUtils

# Define Dictionaries

e_pdg_def = { 'e' : 11, 
              'ebar' : -11 }
e_name_def = { 11 : 'e', 
              -11: 'ebar' }

#Other required information
mcseed = 210921029 

# inputs for event generation jobs
evg_tgtpdg_hash = ['1000010020', '1000010030', '1000020030', '1000020040', '1000030060', '1000060120', '1000080160', '1000130270', 
                   '1000180400', '1000200400', '1000200480', '1000260560', '1000791970', '1000822080', '1000922380']

def eScatteringGenCommands( e_list = "11",tgt_list="1000060120", E_list="2", xspl_file="total_xsec.xml",ntotevents=1000000, tune='G18_02_02_11b',gen_list="EM", nmaxrun=100000,  
                            version='master', conf_dir='', arch='SL6.x86_64', production='routine_validation', cycle='01', grid_system='FNAL', group='genie', 
                            softw_topdir=os.getenv('GENIE_MASTER_DIR'), genie_topdir=os.getenv('GENIE'), jobs_topdir=os.getenv('PWD')) :

    jobs_dir = jobs_topdir+'/'+version+'-'+production+'_'+cycle+'-eScattering'
    # Make directory
    if not os.path.exists(jobs_dir) : 
        os.mkdir(jobs_dir)

    # configure setup 
    if grid_system == 'FNAL' : 
        genie_setup = genie_topdir+'src/scripts/production/python/setup_FNALGrid.sh' ## put correct path
    else : 
        genie_setup = softw_dopdir+'/generator/builds/'+arch+'/'+version+'-setup'

    # Electron list
    req_e_list = e_list.split(',')
    for i in range(len(req_e_list)):
        if req_e_list[i] not in e_pdg_def and e not in e_name_def : 
            print( "Configured lepton is not an electron:"+req_e_list[i])
            return 
        else if req_e_list[i] in e_pdg_def:
            req_e_list[i] = e_pdg_def[req_e_list[i]]

    req_tgt_list = tgt_list.split(',')
    req_E_list = E_list.split(',')
    req_gen_list = gen_list.split(',')
    
    nsubruns = ntotevents/nmaxrun
    if not nsubruns.is_integer() :
        nsubruns = 1+int(nsubruns)

    for e in req_e_list : 
        for tgt in req_tgt_list : 
            for E in req_E_list : 
                n_event_left = ntotevents; 
                nsubruns
                for isubrun in range(nsubruns) :
                    print "Events remaining ".$n_event_left."\n" ;
                    $nev = $n_event_left >= $nmax_per_run ? $nmax_per_run : $n_event_left ;
                    $n_event_left -= $nev ;
                    
                    # Set naming scheme: RREEEESS (RR: run number, EEEE: electron energy in MeV, SS: subrun id)
                    curr_subrune = 1E5 * ( ( 10 + int(tgt) * 10 + int(E) ) )+ isubrun; 
                    curr_seed         = mcseed + isubrun + int(tgt)
                    jobname           = "escattering-"+curr_subrune
                    filename_template = jobs_dir+"/"+jobname
                    
                    evgen_command = "gevgen -p "+e+"-n "+nev+" -e "+E+" -t "+tgt+" -r "+curr_subrune+" --seed "+curr_seed+" --cross-sections "+xspl_file+" --event-generator-list "+gen_list+" --tune "+tune

                    shell_file = GridUtils.CreateShellScript ( evgen_command , jobs_dir, filename_template, filename_template+".xml", grid_system, genie_setup, conf_dir, xspl_file ) 
                
                    if grid_system == 'FNAL' :
                        grid_command_options = GridUtils.FNALShellCommands(genie_setup)
                        command_list.append( "jobsub_submit "+grid_command_options+ " file://"+shell_file )


    ## Add command list to dictionary; Key is 4 => event production
    command_dict = {}
    command_dict[4] = command_list ; 

    return command_dict ; 
