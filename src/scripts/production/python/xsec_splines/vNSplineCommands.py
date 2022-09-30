#! /usr/bin/env python

"""\
This script generates the scripts required to run vN splines. In addition it returns
the commands required to run the scripts in a specific grid. 
The input variables required by the vNSPlineCommands function are:
- version : genie version number
- config_dir : configuration directory
- tune : GENIE configuration tune
- arch : <SL4.x86_32, SL5.x86_64, SL6.x86_64, ...>, default SL6.x86_64
- production : Default routine_validation
- cycle : Default 01
- use_valgrind : Default false 
- grid_system : FNAL, (more to be added in the future). Default FNAL
- group : group at the grid
- softw_topdir : top level dir for softw installations. Default pwd($GENIE_BASEDIR)
- jobs_topdir : top level dir for job files, default: $PWD
- gen_list : comma separated list of event generator list, default all. For electron mode, all corresponds to all EM modes.
- nu_list : comma separated list of neutrino types. Default all neutrinos. It is possible to use electrons instead.
- e_max : maxmium energy of the splines in GeV. Default 200 GeV.
- n_knots : number of knots per spline. Default 100.  
- with_priority : (boolean) set a priority to optimize bulk productions. Default false 
- run_one : (boolean) If called, one of the jobs is run as part of the script instead of submitted via the batch system. Default all the jobs are submitted                   

"""
import os 

# Define Dictionaries
nucleons_pdg = {
    'n' : 1000000010,
    'p' : 1000010010 }
nucleons_name = {
    1000000010 : 'n',
    1000010010 : 'p' } 

nucleons_EW_proc = [ 'none', 'CCRES', 'NCRES', 'CCDIS', 'NCDIS', 'CCDFR', 'NCDFR', 'Fast' ]
nucleons_EM_proc = [ 'none', 'EMRES', 'EMDIS', 'EMQEL' ]

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

def CreateShellScript ( command , jobs_dir, file_name, grid_system, genie_setup, conf_dir ) :
    shell_file = jobs_dir+"/"+file_name+".sh"
    if os.path.exists(shell_file):
        os.remove(shell_file)

    script = open( shell_file, 'w' ) 
    script.write("#!/bin/bash \n")
    if grid_system == 'FNAL':
        script.write("cd $CONDOR_DIR_INPUT \n")
    else :
        script.write("cd "+jobs_dir)
    script.write("source "+genie_setup+" "+conf_dir+" \n")
    if grid_system == 'FNAL':
        script.write("cd $CONDOR_DIR_INPUT \n")
    script.write(command+"\n")
    if grid_system == 'FNAL':
        script.write("ifdh cp -D "+file_name+".xml "+jobs_dir+" \n")

    script.close()
    return shell_file 

def vNSplineCommands( version='master', conf_dir='', tune='G18_02_02_11b', arch='SL6.x86_64', production='routine_validation', cycle='01', grid_system='FNAL', group='genie', 
                      softw_topdir=os.getenv('GENIE_MASTER_DIR'), genie_topdir=os.getenv('GENIE'), jobs_topdir=os.getenv('PWD'), gen_list='all', nu_list='all', 
                      e_max=200, n_knots=100, with_priority=False, run_one=False ) :

    jobs_dir = jobs_topdir+'/'+version+'-'+production+'_'+cycle+'-xsec_vN'

    # configure setup 
    if grid_system == 'FNAL' : 
        genie_setup = genie_topdir+'src/scripts/production/python/setup_FNALGrid.sh' ## put correct path
    else : 
        genie_setup = softw_dopdir+'/generator/builds/'+arch+'/'+version+'-setup'

    req_nu_list = []
    req_e_list = []
    if( nu_list != 'all' ) :
        req_particle_list = nu_list.split(',')
        for particle in req_particle_list:
            if particle in nu_pdg_def : 
                req_nu_list.append(particle)
            if particle in e_pdg_def : 
                req_e_list.append(particle)
            if particle in nu_name_def : 
                req_nu_list.append(nu_name_def[particle])
            if particle in e_name_def : 
                req_e_list.append(e_name_def[particle])
    else : 
        for key in nu_pdg_def : 
            req_nu_list.append(key)
        for key in e_pdg_def : 
            req_e_list.append(key)

    if len(req_nu_list) : 
        print ( 'Neutrino list = ')
        for i in range(len(req_nu_list)): print(req_nu_list[i]); 
 
    if len(req_e_list) : 
        print ( 'Electron list = ') 
        for i in range(len(req_e_list)): print(req_e_list[i]); 
        
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

    if len(req_EW_list) : 
        print( 'EW processes =' )
        for i in range(len(req_EW_list)): print(req_EW_list[i]); 
 
    if len(req_EM_list) : 
        print( 'EM processes =' )
        for i in range(len(req_EM_list)): print(req_EM_list[i]); 

    # Make directory
    if not os.path.exists(jobs_dir) : 
        os.mkdir(jobs_dir)

    # grid options 
    command_list = []
    grid_command_options = ''
    if grid_system == 'FNAL' :
        grid_command_options = "-n --memory=1GB --disk=20GB --expected-lifetime=10h -f "+genie_setup 
        grid_command_options += " --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --lines '+FERMIHTC_AutoRelease=True' "
        grid_command_options += "--lines '+FERMIHTC_GraceMemory=4096' --lines '+FERMIHTC_GraceLifetime=6000' "
    
    # Create neutrino spline commands:
    grid_sub_cmd = []     
    shell_file_list = []
    for nu in req_nu_list : 
        for target in nucleons_pdg : 
            for process in req_EW_list :
                if process == 'none' : continue 
                if process == 'CCQE' : 
                    if target == 'n' and nu_pdg_def[nu] < 0 : continue  
                    if target == 'p' and nu_pdg_def[nu] > 0 : continue 
                if process == 'CCDFR' or process == 'NCDFR' :
                    if target == 'n' : continue 
                if process == 'FastOn' : 
                    event_gen_list = process+target 
                else :
                    event_gen_list = process 
                
                # Job definition
                job_name = nu+'_on_'+target+'_'+process
                if grid_system == 'FNAL' : 
                    filename_template = job_name
                else :
                    filename_template = jobs_dir + '/' + job_name

                gmkspl_cmd = "gmkspl -p "+str(nu_pdg_def[nu])+ " -t "+ str(nucleons_pdg[target]) + " -n "+ str(n_knots) + " -e "+ str(e_max) + " --tune " + tune 
                gmkspl_cmd += " -o "+ filename_template+".xml --event-generator-list " + event_gen_list   
                
                shell_file = CreateShellScript ( gmkspl_cmd , jobs_dir, filename_template, grid_system, genie_setup, conf_dir ) 
                if grid_system == 'FNAL' :
                    command_list.append( "jobsub_submit "+grid_command_options+ " file://"+shell_file )


    # Create electron spline commands:
    
    for e in req_e_list : 
        for target in nucleons_pdg : 
            for process in req_EM_list :
                if process == 'none' : continue 
                event_gen_list = process 
                
                # Job definition
                job_name = e+'_on_'+target+'_'+process
                if grid_system == 'FNAL' : 
                    output_spline = job_name
                else :
                    output_spline = jobs_dir + '/' + job_name
                
                gmkspl_cmd = "gmkspl -p "+str(e_pdg_def[e])+ " -t "+ str(nucleons_pdg[target]) + " -n "+ str(n_knots) + " -e "+ str(e_max) + " --tune " + tune 
                gmkspl_cmd += " -o "+ filename_template+".xml --event-generator-list " + event_gen_list   
                
                shell_file = CreateShellScript ( gmkspl_cmd , jobs_dir, filename_template, grid_system, genie_setup, conf_dir ) 
                if grid_system == 'FNAL' :
                    command_list.append( "jobsub_submit "+grid_command_options+ " file://"+shell_file )


    ## Add command list to dictionary; Key is 0 => free nucleon spline calculation
    command_dict = {}
    command_dict[0] = command_list ; 
    return command_dict ; 
