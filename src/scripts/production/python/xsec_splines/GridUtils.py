#! /usr/bin/env python
import os 

def CreateShellScript ( commands , jobs_dir, shell_name, out_files, grid_system, genie_setup, conf_dir, in_files ) :
    shell_file = jobs_dir+"/"+shell_name+".sh"
    if os.path.exists(shell_file):
        os.remove(shell_file)

    script = open( shell_file, 'w' ) 
    script.write("#!/bin/bash \n")
    if grid_system == 'FNAL':
        script.write("cd $CONDOR_DIR_INPUT \n")
        if isinstance(in_files, list) :
            for i in range(len(in_files)):
                script.write("ifdh cp -D "+in_files[i]+"  $CONDOR_DIR_INPUT \n")
        else : 
            script.write("ifdh cp -D "+in_files+"  $CONDOR_DIR_INPUT \n")

    else :
        script.write("cd "+jobs_dir)
    script.write("source "+genie_setup+" "+conf_dir+" \n")
    if grid_system == 'FNAL':
        script.write("cd $CONDOR_DIR_INPUT \n")

    if isinstance(commands,list):
        for command in commands : 
            script.write(command+"\n")
    else : 
        script.write(commands+"\n")

    if grid_system == 'FNAL':
        if isinstance(out_files, list) :
            for file_name in out_files : 
                script.write("ifdh cp -D "+file_name+" "+jobs_dir+" \n")
        else :
            script.write("ifdh cp -D "+out_files+" "+jobs_dir+" \n")

    script.close()
    return shell_file 

def FNALShellCommands(genie_setup):
    grid_command_options = "-n --memory=1GB --disk=20GB --expected-lifetime=10h -f "+genie_setup 
    grid_command_options += " --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --lines '+FERMIHTC_AutoRelease=True' "
    grid_command_options += "--lines '+FERMIHTC_GraceMemory=4096' --lines '+FERMIHTC_GraceLifetime=6000' "

    return grid_command_options 
