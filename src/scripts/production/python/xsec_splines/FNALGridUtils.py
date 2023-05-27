#! /usr/bin/env python
import os, glob

def CreateShellScript ( commands , jobs_dir, shell_name, out_files, grid_setup, genie_setup, conf_dir, in_files, git_branch, git_loc ) :
    shell_file = jobs_dir+"/"+shell_name+".sh"

    if os.path.exists(shell_file):
        os.remove(shell_file)
        
    script = open( shell_file, 'w' ) 
    script.write("#!/bin/bash \n")
    script.write("cd $CONDOR_DIR_INPUT ;\n")
    script.write("source "+os.path.basename(grid_setup)+" ; \n")
    if conf_dir is not '' : 
        conf_files = glob.glob(conf_dir+"/*.xml")
        for conf_i in conf_files : 
            script.write("mkdir $CONDOR_DIR_INPUT/conf ;\n")
            script.write("ifdh cp -D "+conf_i+"  $CONDOR_DIR_INPUT/conf ;\n")
        conf_dir = "$CONDOR_DIR_INPUT/conf"
    script.write("source "+os.path.basename(genie_setup)+" "+git_loc+" "+git_branch+" "+conf_dir+" ;\n")
    script.write("cd $CONDOR_DIR_INPUT ;\n")

    if isinstance(in_files, list) :
        for i in range(len(in_files)):
            script.write("ifdh cp -D "+in_files[i]+"  $CONDOR_DIR_INPUT ;\n")
    else : 
        script.write("ifdh cp -D "+in_files+"  $CONDOR_DIR_INPUT ;\n")

    if isinstance(commands,list):
        for command in commands : 
            script.write(command+" ;\n")
    else : 
        script.write(commands+" ;\n")

    if isinstance(out_files, list) :
        for file_name in out_files : 
            script.write("ifdh cp -D "+file_name+" "+jobs_dir+" ;\n")
    else :
        script.write("ifdh cp -D "+out_files+" "+jobs_dir+" ;\n")

    script.close()
    return shell_file 

def FNALShellCommands(grid_setup, genie_setup, hours = 10, memory=1, disk=500, GraceMemory=4096, GraceLifeTime=6000):
    grid_command_options = " -n --memory="+str(memory)+"GB --disk="+str(disk)+"MB --expected-lifetime="+str(hours)+"h " 
    grid_command_options += " --OS=SL7 --lines '+FERMIHTC_AutoRelease=True' -f "+grid_setup+" -f "+genie_setup 
    grid_command_options += " --lines '+FERMIHTC_GraceMemory="+str(GraceMemory)+"' --lines '+FERMIHTC_GraceLifetime="+str(GraceLifeTime)+"' --mail_on_error "

    return grid_command_options 

def WriteXMLFile(commands_dict, start, end, jobs_dir, file_name='grid_submission.xml'):
    grid_file = jobs_dir+"/"+file_name
    if os.path.exists(grid_file):
        os.remove(grid_file)

    script = open( grid_file, 'w' )
    in_serial = False 

    for id in range(start,end) :
        command_list = commands_dict[id]
        if id < end - 1 : 
            command_list_next = commands_dict[id+1]
        else : 
            command_list_next = command_list 

        if len(command_list) == 0 : 
            continue
        if len(command_list) == 1 : # serial
            if in_serial == False: 
                script.write("<serial>\n")
                in_serial = True
 
            script.write(command_list[0]+"\n")

            if ( len(command_list_next) != 1 ) or ( id == end - 1 and in_serial == True ) : 
                script.write("</serial>\n")
                in_serial = False

        else : 
            script.write("<parallel>\n")
            for i in range(len(command_list)) : 
                script.write(command_list[i]+"\n")
            script.write("</parallel>\n")
                
    script.close()
    return grid_file

def WriteMainSubmissionFile(jobs_dir, genie_topdir, group, grid_setup='/src/scripts/production/python/setup_FNAL.sh', genie_setup='/src/scripts/production/python/setup_GENIE.sh', in_file_name='grid_submission.xml', expectedlife=60,  memory=1, disk=500, out_file_name='fnal_dag_submit.fnal', jobs=1, role="Analysis"):

    fnal_file = jobs_dir+"/"+out_file_name
    if os.path.exists(fnal_file):
        os.remove(fnal_file)

    script = open( fnal_file, 'w' ) 
    script.write("#!/bin/bash\n")
    script.write("source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups ;\n")
    script.write("setup fife_utils ;\n")
    script.write("jobsub_submit -G "+group+" --OS=SL7 --memory="+str(memory)+"GB --disk="+str(disk)+"MB --expected-lifetime="+str(expectedlife)+"h -N "+str(jobs)+" --role="+role+" --dag file://"+in_file_name+";" )

    return fnal_file
