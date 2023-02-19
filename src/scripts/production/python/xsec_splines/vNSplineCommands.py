#! /usr/bin/env python

"""\
This script generates the scripts required to run vN splines. In addition it returns
the commands required to run the scripts in a specific grid. 

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
nucleons_pdg = {
    'n' : 1000000010,
    'p' : 1000010010 }
nucleons_name = {
    1000000010 : 'n',
    1000010010 : 'p' } 

nucleons_EW_proc = [ 'none', 'CCRES', 'NCRES', 'CCDIS', 'NCDIS', 'CCDFR', 'NCDFR', 'Fast' ]
nucleons_EM_proc = [ 'none', 'EMRES', 'EMDIS', 'EMQE' ]

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

e_pdg_def = { 'e' : 11, 
              'ebar' : -11 }

e_name_def = { 11 : 'e', 
              -11: 'ebar' }

def vNSplineCommands( probe_list='all', gen_list='all', nu_E_max=200, e_E_max=30, nu_n_knots=100, e_n_knots=100, 
                      tune='G18_02_02_11b', version='master', 
                      grid_system='FNAL', group='genie', conf_dir='', arch='SL6.x86_64', production='routine_validation', cycle='01',  
                      softw_topdir=os.getenv('GENIE_MASTER_DIR'), genie_topdir=os.getenv('GENIE'), jobs_topdir=os.getenv('PWD'), 
                      grid_setup= os.getenv('GENIE')+'src/scripts/production/python/setup_FNAL.sh',
                      genie_setup= os.getenv('GENIE')+'src/scripts/production/python/setup_GENIE.sh', time=15, git_branch = 'master', git_loc="https://github.com/GENIE-MC/Generator") :

    jobs_dir = jobs_topdir+'/'+version+'-'+production+'_'+cycle+'-xsec_vN'

    req_nu_list = []
    req_e_list = []
    in_files=[]
    if( probe_list != 'all' ) :
        req_particle_list = probe_list.split(',')
        for particle in req_particle_list:
            if particle in nu_pdg_def : 
                req_nu_list.append(nu_pdg_def[particle])
            if particle in e_pdg_def : 
                req_e_list.append(e_pdg_def[particle])
            if int(particle) in nu_name_def : 
                req_nu_list.append(int(particle))
            if int(particle) in e_name_def : 
                req_e_list.append(int(particle))
    else : 
        for key in nu_name_def : 
            req_nu_list.append(key)
        for key in e_name_def : 
            req_e_list.append(key)

    req_EW_list = []
    req_EM_list = []
    if ( gen_list != 'all' ) :
        req_gen_list = gen_list.split(',')
        for process in req_gen_list :
            if process in nucleons_EW_proc :
                req_EW_list.append(process)
            if process in nucleons_EM_proc : 
                req_EM_list.append(process)
    else : 
        req_EW_list = nucleons_EW_proc 
        req_EM_list = nucleons_EM_proc

    # Make directory
    if not os.path.exists(jobs_dir) : 
        os.mkdir(jobs_dir)

    # grid options 
    command_list = []

    grid_command_options = ''
    if grid_system == 'FNAL' :
        grid_command_options = FNAL.FNALShellCommands(grid_setup, genie_setup,time)
    else :
        print( "Only FNAL grid is available" )
        return 

    # Create neutrino spline commands:
    grid_sub_cmd = []     
    shell_file_list = []
    for nu in req_nu_list : 
        for target in nucleons_pdg : 
            for process in req_EW_list :
                if process == 'none' : continue 
                if process == 'CCQE' : 
                    if target == 'n' and nu_pdg_def[nu] < 0 : continue  
                    elif target == 'p' and nu_pdg_def[nu] > 0 : continue 
                if process == 'CCDFR' or process == 'NCDFR' :
                    if target == 'n' : continue 
                if process == 'Fast' :
                    if target == 'n' : event_gen_list = 'FastOnN'
                    elif target == 'p' : event_gen_list = 'FastOnP'
                else :
                    event_gen_list = process 
                
                # Job definition
                job_name = nu_name_def[nu]+'_on_'+target+'_'+process
                filename_template = ''
                if grid_system == 'FNAL' : 
                    filename_template = job_name
                else :
                    filename_template = jobs_dir + '/' + job_name

                if os.path.exists( jobs_dir + '/' + filename_template+".xml" ) : 
                    # Need to remove xml files before re-generating them
                    os.remove( jobs_dir + '/' + filename_template+".xml" )

                gmkspl_cmd = "gmkspl -p "+str(nu)+ " -t "+ str(nucleons_pdg[target]) + " -n "+ str(nu_n_knots) + " -e "+ str(nu_E_max) + " --tune " + tune 
                gmkspl_cmd += " -o "+ filename_template+".xml --event-generator-list " + event_gen_list   
                
                shell_file = ''
                if grid_system == 'FNAL' :
                    shell_file = FNAL.CreateShellScript ( gmkspl_cmd , jobs_dir, filename_template, filename_template+".xml", grid_setup, genie_setup, conf_dir, in_files, git_branch, git_loc ) 
                    command_list.append( "jobsub_submit "+grid_command_options+ " file://"+shell_file )

    # Create electron spline commands:
    
    for e in req_e_list : 
        for target in nucleons_pdg : 
            for process in req_EM_list :
                if process == 'none' : continue 
                event_gen_list = process 
                
                # Job definition
                job_name = e_name_def[e]+'_on_'+target+'_'+process
                if grid_system == 'FNAL' : 
                    output_spline = job_name
                else :
                    output_spline = jobs_dir + '/' + job_name
                
                if os.path.exists( jobs_dir + '/' + output_spline+".xml" ) : 
                    # Need to remove xml files before re-generating them
                    os.remove( jobs_dir + '/' + output_spline+".xml" )

                gmkspl_cmd = "gmkspl -p "+str(e)+ " -t "+ str(nucleons_pdg[target]) + " -n "+ str(e_n_knots) + " -e "+ str(e_E_max) + " --tune " + tune 
                gmkspl_cmd += " -o "+ output_spline+".xml --event-generator-list " + event_gen_list   
                
                shell_file = ''
                if grid_system == 'FNAL' :
                    shell_file = FNAL.CreateShellScript ( gmkspl_cmd , jobs_dir, output_spline, output_spline+".xml", grid_setup, genie_setup, conf_dir, in_files, git_branch, git_loc ) 
                    command_list.append( "jobsub_submit "+grid_command_options+ " file://"+shell_file )


    ## Add command list to dictionary; Key is 0 => free nucleon spline calculation
    command_dict = {}
    command_dict[0] = command_list ; 
    return command_dict ; 
