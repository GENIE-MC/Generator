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
import GridUtils 

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


def GroupSplineCommands( version='master', conf_dir='', tune='G18_02_02_11b', arch='SL6.x86_64', production='routine_validation', cycle='01', grid_system='FNAL', group='genie', 
                         softw_topdir=os.getenv('GENIE_MASTER_DIR'), genie_topdir=os.getenv('GENIE'), jobs_topdir=os.getenv('PWD'), xml_dir=os.getenv('PWD'), mother_dir='', 
                         group_vN=False,add_list=False, root_output = False, add_nucleons = False ) :

    if mother_dir != '' : 
        if os.path.exists(mother_dir) :
            xml_files_motherdir = glob.glob(mother_dir+"/*.xml")
        else :
            print ( mother_dir+"doesn't exist")
            return  

        #Given a mother directory and a daughter directory, the script tryies
        # to copy (ln -s) all the files in the mother dir into the daughter
        # If the file alredy exists, the link is not created.
        for xml_file in xml_files_motherdir :          
            # Check if exist in xml_dir 
            xml_file_name = os.path.basename(xml_file)
            if os.path.exists(xml_dir+"/"+xml_file_name) : continue 
            os.link(xml_file,xml_dir+"/"+xml_file_name)
            
    if os.path.exists(xml_dir) :
        xml_files_dir = glob.glob(xml_dir+"/*.xml")
        if len(xml_files_dir) ==0 : 
            return 
    else : 
        print ( xml_dir+"doesn't exist")
        return 

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
    for xml_file in xml_files_dir : 
        xml_file = os.path.basename(xml_file)[:-4]
        xml_content = xml_file.split("_")
        if len(xml_content) < 4 : continue

        if xml_content[0] in nu_pdg_def : 
            dir_nu_list.append(xml_content[0])
            dir_nu_tgt_list.append(xml_content[2])
            dir_EW_process_list.append(xml_content[3])
        elif xml_content[0] in e_pdg_def : 
            dir_e_list.append(xml_content[0])
            dir_e_tgt_list.append(xml_content[2])
            dir_EM_process_list.append(xml_content[3])

        ## For root output 
        if root_output : 
            if xml_content[0] not in lepton_list : 
                lepton_list.append(xml_content[0]) 
            if xml_content[2] not in tgt_list : 
                tgt_list.append(xml_content[2]) 

    dict_target = {}
    for target in dir_nu_tgt_list : 
        if len(target) == 1 : 
            if target == 'p' : 
                target = '1000010010'
            if target == 'n' :
                target = '1000000010'

        dict_nu = {}
        for nu in dir_nu_list : 
            dict_nu[nu] = []
            for process in dir_EW_process_list : 
                if os.path.isfile(xml_dir+"/"+nu+"_on_"+target+"_"+process+".xml") :
                    dict_nu[nu].append(nu+"_on_"+target+"_"+process+".xml")
        
        for e in dir_e_list : 
            dict_nu[e] = []
            for process in dir_EM_process_list : 
                if os.path.isfile(xml_dir+"/"+e+"_on_"+target+"_"+process+".xml") :
                    dict_nu[e].append(e+"_on_"+target+"_"+process+".xml")

        # Add all files to merge here:
        dict_target[target] = dict_nu 

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
            commands.append("mv "+path+nu+"_on_"+tgt+".xml "+tgt+".xml")
        else :
            commands.append(com_nu) 
        com_total += tgt+".xml,"
    com_total = com_total[:-1]

    ## if only one target simply rename
    if len(dict_target) == 1 : 
        commands.append("mv "+path+tgt+".xml "+path+"total_xsec.xml")
    else :
        commands.append(com_total) 

    out_files = [ "total_xsec.xml" ] 
    if root_output :
        str_nu_list = ''
        str_tgt_list = ''
        for tgt in tgt_list:
            str_tgt_list += tgt+","
        for lepton in lepton_list : 
            str_nu_list += lepton+","
        str_nu_list = str_nu_list[:-1]
        str_tgt_list = str_tgt_list[:-1]

        ## Create an output file with all the splines in root format
        commands.append( "gspl2root -p "+str_nu_list+" -t "+str_tgt_list+" -f "+path+"total_xsec.xml -o "+path+"total_xsec.root --tune "+tune )
        out_files.append("total_xsec.root")

    # configure setup 
    if grid_system == 'FNAL' : 
        genie_setup = genie_topdir+'src/scripts/production/python/setup_FNALGrid.sh' ## put correct path
    else : 
        genie_setup = softw_dopdir+'/generator/builds/'+arch+'/'+version+'-setup'

    if groupvN == True : 
        process_name = "group_vN"
        job_ID = 1 
    else : 
        process_name = "group_vA"
        job_ID = 3 

    # Call Commands
    shell_file = GridUtils.CreateShellScript ( commands , jobs_topdir, process_name, out_files, grid_system, genie_setup, conf_dir, xml_files_dir ) 

    if grid_system == 'FNAL' :
        command_list.append( "jobsub_submit "+grid_command_options+ " file://"+shell_file )

    ## Add command list to dictionary; 
    command_dict = {}
    command_dict[job_ID] = command_list ; 
    return command_dict ; 
