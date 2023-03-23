#!/usr/bin/env python

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
import FNALGridUtils as FNAL

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

def eScatteringGenCommands( e_list = "11",tgt_list="1000060120", EBeam_list="2", xspl_file="total_xsec.xml",ntotevents=1000000, 
                            tune='G18_02_02_11b',gen_list="EM", nmaxrun=100000, gst_output=False, no_ghep=False,version='master', conf_dir='', arch='SL6.x86_64', 
                            production='routine_validation', cycle='01', grid_system='FNAL', group='genie', 
                            softw_topdir=os.getenv('GENIE_MASTER_DIR'), genie_topdir=os.getenv('GENIE'), jobs_topdir=os.getenv('PWD'),
                            grid_setup = os.getenv('GENIE')+'src/scripts/production/python/setup_FNAL.sh', 
                            genie_setup= os.getenv('GENIE')+'src/scripts/production/python/setup_GENIE.sh', time='10', memory='1', disk='2000',git_branch = "master", git_loc="https://github.com/GENIE-MC/Generator") :

    jobs_dir = jobs_topdir+'/'+version+'-'+production+'_'+cycle+'-eScattering'
    # Make directory
    if not os.path.exists(jobs_dir) : 
        os.mkdir(jobs_dir)

    # Electron list
    final_e_list = []
    if e_list != 'all':
        req_e_list = e_list.split(',')
        for i in range(len(req_e_list)):
            if req_e_list[i] not in e_pdg_def and int(req_e_list[i]) not in e_name_def : 
                continue  
            if req_e_list[i] in e_pdg_def:
                final_e_list.append( e_pdg_def[req_e_list[i]])
            else :
                final_e_list.append( req_e_list[i] )
    else : 
        final_e_list.append(11) 

    if tgt_list != 'all':
        req_tgt_list = tgt_list.split(',')
    else : 
        req_tgt_list = evg_tgtpdg_hash 
    req_En_list = EBeam_list.split(',')
    req_gen_list = gen_list.split(',')
        
    nsubruns = ntotevents/nmaxrun
    if ntotevents < nmaxrun : nsubruns = 1

    if not isinstance(nsubruns, int) :
        nsubruns = 1+int(nsubruns)

    xsec_filename = os.path.basename(xspl_file)

    if grid_system == 'FNAL' :
        input_xsec = "\$CONDOR_DIR_INPUT/"+xsec_filename
    else :
        input_xsec = free_nuc_dir+xsec_filename

    command_list = []
    for e in final_e_list : 
        for tgt in req_tgt_list : 
            for E in req_En_list : 
                n_event_left = ntotevents
                for isubrun in range(nsubruns) :
                    if n_event_left >= nmaxrun : 
                        nev = nmaxrun
                    else:
                        nev = n_event_left
                    n_event_left -= nev 
                    curr_subrune = "11"+str(tgt)+str(isubrun); 
                    curr_seed         = mcseed + isubrun + int(tgt)
                    jobname           = "e_on_"+str(tgt)+"_"+str(int((float(E)*1000)))+"MeV_"+str(isubrun)

                    evgen_command = "gevgen -p "+str(e)+" -n "+str(nev)+" -e "+E+" -t "+str(tgt)+" -r "+curr_subrune+" --seed "+str(curr_seed)
                    evgen_command += " --cross-sections "+input_xsec+" --event-generator-list "+gen_list+" --tune "+tune + " -o "+jobname+".ghep.root "
                    
                    out_files = [str(jobname+".ghep.root")]
                    if gst_output : 
                        evgen_command += " ; gntpc -i "+jobname+".ghep.root -o "+jobname+".gst.root -f gst "
                        out_files.append(str(jobname+".gst.root"))
                        if no_ghep :
                            out_files = [str(jobname+".gst.root")]

                    shell_file = ''                
                    if grid_system == 'FNAL' :
                        shell_file= FNAL.CreateShellScript ( evgen_command , jobs_dir, jobname, out_files, grid_setup, genie_setup, conf_dir, str(xspl_file), git_branch, git_loc ) 
                        grid_command_options = FNAL.FNALShellCommands(grid_setup, genie_setup,time,memory,disk)
                        command_list.append( "jobsub_submit "+grid_command_options+ " file://"+shell_file )

    ## Add command list to dictionary; Key is 4 => event production
    command_dict = {}
    command_dict[4] = command_list ; 

    return command_dict ; 
