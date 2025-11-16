#!/usr/bin/env python

"""\
This script generates hadron-scattering events 

Author: 
      Julia Tena Vidal <jtenavidal \st tauex.tau.ac.il>
      Tel Aviv University
Copyright:
   Copyright (c) 2003-2025, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org
"""
import os 
import FNALGridUtils as FNAL

# Define Dictionaries

hadron_pdg_def = { 'pip' : 211, 'pim' : -211, 'p' : 2212, 'n' : 2112 }

hadron_name_def = { 211: 'pip', -211: 'pim', 2212: 'p', 2112 : 'n' }

# inputs for event generation jobs
evg_tgtpdg_hash = ['1000010020', '1000010030', '1000020030', '1000020040', '1000030060', '1000060120', '1000080160', '1000130270', 
                   '1000180400', '1000200400', '1000200480', '1000260560', '1000791970', '1000822080', '1000922380']

def hadronScatteringGenCommands( hadron_list = "212",tgt_list="1000060120", KEBeam_list="2", ntotevents=1000000, KEBeamMin="2", KEBeamMax="3", Flux="none",
                                 tune='G18_02a_02_11b',nmaxrun=100000, mcseed=210921029, ginuke_output=False, no_ghep=False,version='master', 
                                 conf_dir='', arch='SL6.x86_64', production='routine_validation', cycle='01', grid_system='FNAL', group='genie', 
                                 softw_topdir=os.getenv('GENIE_MASTER_DIR'), genie_topdir=os.getenv('GENIE'), jobs_topdir=os.getenv('PWD'),
                                 grid_setup = os.getenv('GENIE')+'src/scripts/production/python/setup_FNAL.sh', 
                                 genie_setup= os.getenv('GENIE')+'src/scripts/production/python/setup_GENIE.sh',
                                 message_thresholds= os.getenv('GENIE')+'config/Messenger.xml',
                                 time='10', memory='1GB', disk='2GB',
                                 git_branch = "master", git_loc="https://github.com/GENIE-MC/Generator", 
                                 configure_hA= True, configure_hN = False, configure_INCL=False, configure_G4=False, GHEPMC3Output=False ) :

    jobs_dir = jobs_topdir+'/'+version+'-'+production+'_'+cycle+'-hadronScattering'
    # Make directory
    if not os.path.exists(jobs_dir) : 
        os.mkdir(jobs_dir)

    # Electron list
    final_hadron_list = []
    if hadron_list != 'all':
        req_hadron_list = hadron_list.split(',')
        for i in range(len(req_hadron_list)):
            if req_hadron_list[i] not in hadron_pdg_def and int(req_hadron_list[i]) not in hadron_name_def : 
                continue  
            if req_hadron_list[i] in hadron_pdg_def:
                final_hadron_list.append( hadron_pdg_def[req_hadron_list[i]])
            else :
                final_hadron_list.append( req_hadron_list[i] )
    else : 
        final_hadron_list.append(11) 

    if tgt_list != 'all':
        req_tgt_list = tgt_list.split(',')
    else : 
        req_tgt_list = evg_tgtpdg_hash 
    req_KEn_list = KEBeam_list.split(',')
    
    model = "hA"
    if configure_hA : 
        model = "hA2018"
    elif configure_hN :
        model = "hN2018"
    elif configure_INCL :
        model = "HINCL"
    elif configure_G4 : 
        model = "HG4BertCasc"

    true_nsubruns = ntotevents*1.0/nmaxrun
    nsubruns = int(round(ntotevents*1.0/nmaxrun))
    if( nsubruns < true_nsubruns ) : nsubruns += 1 
    if ntotevents <= nmaxrun : nsubruns = 1
    
    command_list = []
    for hadron in final_hadron_list : 
        for tgt in req_tgt_list : 
            for KE in req_KEn_list : 
                n_event_left = ntotevents
                
                for isubrun in range(nsubruns) :
                    if n_event_left >= nmaxrun : 
                        nev = nmaxrun
                    else:
                        nev = n_event_left
                    n_event_left -= nev
                    curr_subrune = str(abs(int(hadron)))+str(tgt)+str(isubrun); 
                    curr_seed         = int(mcseed) + isubrun + int(tgt)
                    jobname           = hadron_name_def[int(hadron)]+"_on_"+str(tgt)+"_"+str(int((float(KE)*1000)))+"MeV_"+str(isubrun)

                    evgen_command = "gevgen_hadron -p "+str(hadron)+" -n "+str(nev)+" -t "+str(tgt)+" -r "+curr_subrune+" --seed "+str(curr_seed)
                    if Flux is "none" : 
                        evgen_command += " -k "+KE
                    elif Flux is "uniform" :
                        evgen_command += " -k "+KEBeamMin+","+KEBeamMax
                    else :
                        evgen_command += " -k "+KEBeamMin+","+KEBeamMax+" -f "+Flux
                    evgen_command += " -m "+model + " -o "+jobname
                    evgen_command += " --message-thresholds "+message_thresholds 

                    out_files = [str(jobname+"*.ghep.root")]
                    if ginuke_output : 
                        evgen_command += " ; gntpc -i "+jobname+"*.ghep.root -o "+jobname+".ginuke.root -f ginuke --message-thresholds "+message_thresholds
                        out_files.append(str(jobname+".ginuke.root"))
                        if no_ghep :
                            out_files = [str(jobname+".ginuke.root")]

                    shell_file = ''                
                    if grid_system == 'FNAL' :
                        shell_file= FNAL.CreateShellScript ( evgen_command , jobs_dir, jobname, out_files, grid_setup, genie_setup, conf_dir, "", 
                                                             git_branch, git_loc, configure_INCL, configure_G4, GHEPMC3Output ) 
                        grid_command_options = FNAL.FNALShellCommands(grid_setup, genie_setup,time,memory,disk)
                        command_list.append( "jobsub_submit "+grid_command_options+ " file://"+shell_file )

    command_dict = {}
    command_dict[0] = command_list ; 

    return command_dict ; 
