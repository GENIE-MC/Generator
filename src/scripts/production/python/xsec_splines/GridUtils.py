#! /usr/bin/env python
import os 

def CreateShellScript ( command , jobs_dir, file_name, grid_system, genie_setup, conf_dir, in_files ) :
    shell_file = jobs_dir+"/"+file_name+".sh"
    if os.path.exists(shell_file):
        os.remove(shell_file)

    script = open( shell_file, 'w' ) 
    script.write("#!/bin/bash \n")
    if grid_system == 'FNAL':
        script.write("cd $CONDOR_DIR_INPUT \n")
        if len(in_files) : 
            for i in range(len(in_files)):
                script.write("ifdh cp -D "+in_files[i]+"  $CONDOR_DIR_INPUT \n")
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

def FNALShellCommands(genie_setup):
    grid_command_options = "-n --memory=1GB --disk=20GB --expected-lifetime=10h -f "+genie_setup 
    grid_command_options += " --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE --lines '+FERMIHTC_AutoRelease=True' "
    grid_command_options += "--lines '+FERMIHTC_GraceMemory=4096' --lines '+FERMIHTC_GraceLifetime=6000' "

    return grid_command_options 
