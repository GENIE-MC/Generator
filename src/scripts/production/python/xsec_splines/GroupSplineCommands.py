#! /usr/bin/env python

"""\
 Given a directory, the script looks for file in the form v*_on_*_*.xml 
 then it calls the appropriate gspladd to obtain v*_on_*.xml comprehensive
 file. If a spline is missing it copies the files from a specified directory, otherwise it aborts. 
 The file name will be total_xsec.xml
 The scpript also prodces a root output for the splines

Author: 
      Julia Tena Vidal <jtenavidal \st tauex.tau.ac.il>
      Tel Aviv University
Copyright:
   Copyright (c) 2003-2022, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org

"""
import os, glob, re
import FNALGridUtils as FNAL

# Define Dictionaries
nucleons_pdg = {
    'n' : 1000000010,
    'p' : 1000010010 }
nucleons_name = {
    1000000010 : 'n',
    1000010010 : 'p' } 
tgt_pdg = [1000010020, 1000010030, 1000020030, 1000020040, 1000060120, 1000080160, 1000130270, 1000200400, 1000200480, 1000260560, 1000791970, 1000822080, 1000922380 ]

nucleons_EW_proc = [ 'none', 'CCRES', 'NCRES', 'CCDIS', 'NCDIS', 'CCDFR', 'NCDFR', 'Fast' ]
nucleons_EM_proc = [ 'none', 'EMRES', 'EMDIS', 'EMQEL' ]
nuclei_EW_proc = [ 'none', 'WeakMEC', 'CCCOHPION', 'NCCOHPION', 'Fast' ]
nuclei_EM_proc = [ 'none', 'EMMEC', 'EMQE', 'EMRES', 'EMDIS' ]

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


def GroupSplineCommands( group_vN=False, xml_dir=os.getenv('PWD'), mother_dir='', tune='G18_02_02_11b', gen_list='all',version='master', conf_dir='', grid_system='FNAL', group='genie', 
                         arch='SL6.x86_64', production='routine_validation', cycle='01', softw_topdir=os.getenv('GENIE_MASTER_DIR'),
                         genie_topdir=os.getenv('GENIE'), grid_setup = os.getenv('GENIE')+'src/scripts/production/python/setup_FNAL.sh',
                         genie_setup = os.getenv('GENIE')+'src/scripts/production/python/setup_GENIE.sh', 
                         jobs_topdir=os.getenv('PWD'), add_list=False, add_nucleons = False, time=2, git_branch="master", git_loc="https://github.com/GENIE-MC/Generator" ) :
    
    # Store root output only for vA spilnes:
    root_output = False 
    if group_vN == False: 
        root_output = True

    if not os.path.exists(xml_dir) :
        print ( xml_dir+" doesn't exist")
        return 

    store_total_xsec = False
    if gen_list == 'none' :
        store_total_xsec = True 


    if group_vN == True : 
        process_name = "group_vN"
        job_ID = 1 
    else : 
        process_name = "group_vA"
        job_ID = 3 

    if mother_dir != '' : 
        if os.path.exists(mother_dir) :
            xml_files_motherdir = glob.glob(mother_dir+"/*.xml")
        else :
            print ( mother_dir+" doesn't exist")
            return  

        #Given a mother directory and a daughter directory, the script tryies
        # to copy (ln -s) all the files in the mother dir into the daughter
        # If the file alredy exists, the link is not created.
        for xml_file in xml_files_motherdir :          
            # Check if exist in xml_dir 
            xml_file_name = os.path.basename(xml_file)
            if os.path.exists(xml_dir+"/"+xml_file_name[:-4]+".sh") : continue 
            
            if xml_file_name[:-4] == 'total_xsec' : 
                if store_total_xsec == True : 
                    os.link(xml_file,xml_dir+"/"+xml_file_name) # link xml files
                    continue 
            os.link(xml_file,xml_dir+"/"+xml_file_name) # link xml files
            os.link(xml_file[:-4]+".sh",xml_dir+"/"+xml_file_name[:-4]+".sh") # link sh files
        if store_total_xsec == True : 
            temp_command_dict = {}
            temp_command_dict[job_ID] = []
            return temp_command_dict 

    # Get names of sh files: these determine the name of the future xml files
    xml_files_dir = glob.glob(xml_dir+"/*.sh")
    
    # Store nu, tgt and process that have a corresponding xml file
    dir_nu_list = []
    dir_e_list = []
    dir_nu_tgt_list = []
    dir_e_tgt_list = []
    dir_EW_process_list = []
    dir_EM_process_list = []
    ## sore for later
    lepton_list = []
    tgt_list = []
    in_xml_files = []
    for xml_file in xml_files_dir : 
        xml_file = os.path.basename(xml_file)[:-3]
        in_xml_files.append(xml_dir+xml_file+".xml")
        xml_content = xml_file.split("_")
        if len(xml_content) < 4 : continue

        if xml_content[0] in nu_pdg_def : 
            if xml_content[0] not in dir_nu_list : 
                dir_nu_list.append(xml_content[0])
            if xml_content[2] not in dir_nu_tgt_list: 
                dir_nu_tgt_list.append(xml_content[2])
            if xml_content[3] not in dir_EW_process_list:
                dir_EW_process_list.append(xml_content[3])
        elif xml_content[0] in e_pdg_def : 
            if xml_content[0] not in dir_e_list : 
                dir_e_list.append(xml_content[0])
            if xml_content[2] not in dir_e_tgt_list : 
                dir_e_tgt_list.append(xml_content[2])
            if xml_content[3] not in dir_EM_process_list:
                dir_EM_process_list.append(xml_content[3])

        ## For root output 
        if root_output : 
            if xml_content[0] not in lepton_list : 
                lepton_list.append(xml_content[0]) 
            if xml_content[2] not in tgt_list : 
                tgt_list.append(xml_content[2]) 

    dict_target = {}
    for target in dir_nu_tgt_list : 
        dict_nu = {}
        for nu in dir_nu_list : 
            dict_nu[nu] = []
            for process in dir_EW_process_list : 
                if process == 'CCDFR' or process == 'NCDFR' :
                    if target == 'n' or target == '1000000010': continue
                dict_nu[nu].append(nu+"_on_"+target+"_"+process+".xml")
        # Add all files to merge here:
        dict_target[target] = dict_nu 
    
    for target in dir_e_tgt_list : 
        dict_e={}
        if target in dict_target : 
            dict_e = dict_target[target]
        for e in dir_e_list : 
            dict_e[e] = []
            for process in dir_EM_process_list : 
                dict_e[e].append(e+"_on_"+target+"_"+process+".xml")

        # Add all files to merge here:
        dict_target[target] = dict_e 

    if grid_system == 'FNAL' : 
        path = "$CONDOR_DIR_INPUT/"
    else : 
        path = xml_dir

    commands = []
    com_total = "gspladd -o "+path+"total_xsec.xml -f "        
    for tgt in dict_target : 
        com_nu = "gspladd -o "+path+tgt+".xml -f "
        for nu in dict_target[tgt]:
            com_proc = "gspladd -o "+path+nu+"_on_"+tgt+".xml -f "
            for file_proc in dict_target[tgt][nu] : 
                com_proc += path+file_proc + ","
            com_proc = com_proc[:-1]
            commands.append(com_proc) 
            com_nu += path+nu+"_on_"+tgt+".xml,"
        com_nu = com_nu[:-1]
        if len(dict_target[tgt]) == 1 : 
            commands.append("ifdh cp "+path+nu+"_on_"+tgt+".xml "+path+tgt+".xml")
        else :
            commands.append(com_nu) 
        com_total += tgt+".xml,"
    com_total = com_total[:-1]

    ## if only one target simply rename
    if len(dict_target) == 1 : 
        commands.append("ifdh cp "+path+tgt+".xml "+path+"total_xsec.xml")
    else :
        commands.append(com_total) 

    out_files = [ "total_xsec.xml" ] 
    if os.path.exists( xml_dir + '/total_xsec.xml' ) :
        if store_total_xsec == False :
            # Need to remove xml files before re-generating them                                                                                                                                                     
            os.remove( xml_dir + '/total_xsec.xml' )
        
    if root_output :
        # Check if file exists - and remove
        if os.path.exists( xml_dir + '/total_xsec.root' ) :
            # Need to remove xml files before re-generating them                                                                                                                                              
            os.remove( xml_dir + '/total_xsec.root' )

        str_probe_list = ''
        str_tgt_list = ''
        for tgt in tgt_list:
            str_tgt_list += tgt+","
        for lepton in lepton_list : 
            if lepton in nu_pdg_def : 
                str_probe_list += str(nu_pdg_def[lepton])+","
            elif lepton in e_pdg_def : 
                str_probe_list += str(e_pdg_def[lepton])+","
        str_probe_list = str_probe_list[:-1]
        str_tgt_list = str_tgt_list[:-1]

        ## Create an output file with all the splines in root format
        commands.append( "gspl2root -p "+str_probe_list+" -t "+str_tgt_list+" -f "+path+"total_xsec.xml -o "+path+"total_xsec.root --tune "+tune )
        out_files.append("total_xsec.root")

    # Call Commands
    shell_file = ''
    command_list = []
    if grid_system == 'FNAL' :
        shell_file=FNAL.CreateShellScript ( commands , xml_dir, process_name, out_files, grid_setup, genie_setup, conf_dir, in_xml_files, git_branch, git_loc ) 
        grid_command_options = FNAL.FNALShellCommands(grid_setup, genie_setup, time)
        command_list.append( "jobsub_submit "+grid_command_options+ " file://"+shell_file )

    ## Add command list to dictionary; 
    command_dict = {}
    command_dict[job_ID] = command_list 
    return command_dict 
