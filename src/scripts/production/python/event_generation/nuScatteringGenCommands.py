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

nu_pdg_def = { 've'      :   12,
               'vebar'   :  -12,
               'vmu'     :   14,
               'vmubar'  :  -14,
               'vtau'    :   16,
               'vtaubar' :  -16 }

nu_name_def = { 12  : 've'     ,
                -12 : 'vebar'  ,
                 14 : 'vmu'    ,
                -14 : 'vmubar' ,
                 16 : 'vtau'   ,
                -16 : 'vtaubar' }

tgt_pdg = [1000010020, 1000010030, 1000020030, 1000020040, 1000060120, 1000080160, 1000130270, 1000200400, 1000200480, 1000260560, 1000791970, 1000822080, 1000922380 ]

#Other required information
mcseed = 210921029 

# inputs for event generation jobs
evg_tgtpdg_hash = ['1000010020', '1000010030', '1000020030', '1000020040', '1000030060', '1000060120', '1000080160', '1000130270', 
                   '1000180400', '1000200400', '1000200480', '1000260560', '1000791970', '1000822080', '1000922380']

def nuScatteringGenCommands( nu_list = "14",tgt_mix="", EFlux_min=0, EFlux_max=100, flux="\'1/x\'", xspl_file="total_xsec.xml",ntotevents=1000000, 
                            tune='G18_02_02_11b',gen_list="all", expname="general",nmaxrun=100000, gst_output=False, no_ghep=False,version='master', conf_dir='', arch='SL6.x86_64', 
                            production='routine_validation', cycle='01', grid_system='FNAL', group='genie', 
                            softw_topdir=os.getenv('GENIE_MASTER_DIR'), genie_topdir=os.getenv('GENIE'), jobs_topdir=os.getenv('PWD'),
                            grid_setup = os.getenv('GENIE')+'src/scripts/production/python/setup_FNAL.sh', 
                            genie_setup= os.getenv('GENIE')+'src/scripts/production/python/setup_GENIE.sh', 
                            time='10', git_branch = "master", git_loc="https://github.com/GENIE-MC/Generator") :

    jobs_dir = jobs_topdir+'/'+version+'-'+production+'_'+cycle+'-nuScattering'
    # Make directory
    if not os.path.exists(jobs_dir) : 
        os.mkdir(jobs_dir)

    if tgt_mix == 'all':
        tgt_mix = '1000060120'

    # Electron list
    final_nu_list = []
    if nu_list != 'all':
        req_nu_list = nu_list.split(',')
        for i in range(len(req_nu_list)):
            if req_nu_list[i] not in nu_pdg_def and int(req_nu_list[i]) not in nu_name_def : 
                continue  
            if req_nu_list[i] in nu_pdg_def:
                final_nu_list.append( nu_pdg_def[req_nu_list[i]])
            else :
                final_nu_list.append( req_nu_list[i] )
    else : 
        final_nu_list.append(14) 

    req_gen_list = gen_list.split(',')
        
    nsubruns = ntotevents/nmaxrun
    if ntotevents < nmaxrun : nsubruns = 1

    if not isinstance(nsubruns, int) :
        nsubruns = 1+int(nsubruns)

    if grid_system == 'FNAL' :
        input_xsec = "\$CONDOR_DIR_INPUT/total_xsec.xml"
    else :
        input_xsec = free_nuc_dir+"/total_xsec.xml"

    in_file_list = [str(xspl_file)]
    flux_path = ''
    if "root" in flux:
        flux_members = flux.split(",") # separate flux file location and histogram name
        if len(flux_members) is not 2 : 
            print( "Provide flux file path and histogram name" )
            exit()

        flux_path = flux_members[0] # file location 
        flux_histname = flux_members[1]
        flux_filename = os.path.basename( flux_members[0] )
        flux = "\$CONDOR_DIR_INPUT/"+flux_filename+","+flux_histname
        in_file_list.append(flux_path)

    command_list = []
    for nu in final_nu_list : 
        n_event_left = ntotevents
        for isubrun in range(nsubruns) :
            if n_event_left >= nmaxrun : 
                nev = nmaxrun
            else:
                nev = n_event_left
            n_event_left -= nev 
            curr_subrune = "14"+str(isubrun); 
            curr_seed         = mcseed + isubrun 
            jobname           = "nu_"+expname+"_"+str(isubrun)            
            evgen_command = "gevgen -p "+str(nu)+" -n "+str(nev)+" -e "+EFlux_min+","+EFlux_max+" -f " +flux+" -t "+str(tgt_mix)+" -r "+curr_subrune+" --seed "+str(curr_seed)
            evgen_command += " --cross-sections "+input_xsec+" --tune "+tune + " -o "+jobname+".ghep.root"

            if gen_list is not "all" : 
                evgen_command += " --event-generator-list "+gen_list+" "
                shell_file = ''                

                out_files = [str(jobname+".ghep.root")]
                if gst_output : 
                    evgen_command += " ; gntpc -i "+jobname+".ghep.root -o "+jobname+".gst.root -f gst "
                    out_files.append(str(jobname+".gst.root"))
                    if no_ghep :
                        out_files = [str(jobname+".gst.root")]
            if grid_system == 'FNAL' :
                shell_file= FNAL.CreateShellScript ( evgen_command , jobs_dir, jobname, out_files, grid_setup, genie_setup, conf_dir, in_file_list, git_branch, git_loc ) 
                grid_command_options = FNAL.FNALShellCommands(grid_setup, genie_setup,time)
                command_list.append( "jobsub_submit "+grid_command_options+ " file://"+shell_file )

    ## Add command list to dictionary; Key is 4 => event production
    command_dict = {}
    command_dict[4] = command_list ; 

    return command_dict ; 
