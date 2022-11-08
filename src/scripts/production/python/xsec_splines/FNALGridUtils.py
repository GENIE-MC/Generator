#! /usr/bin/env python
import os 

def CreateShellScript ( commands , jobs_dir, shell_name, out_files, genie_setup, conf_dir, in_files ) :
    shell_file = jobs_dir+"/"+shell_name+".sh"

    if os.path.exists(shell_file):
        os.remove(shell_file)

    script = open( shell_file, 'w' ) 
    script.write("#!/bin/bash \n")
    script.write("cd $CONDOR_DIR_INPUT ;\n")
    script.write("source "+os.path.basename(genie_setup)+" "+conf_dir+" ;\n")
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

def FNALShellCommands(genie_setup, hours = 10, memory=1, disk=20, GraceMemory=4096,GraceLifeTime=6000):
    grid_command_options = "-n --memory="+str(memory)+"GB --disk="+str(disk)+"GB --expected-lifetime="+str(hours)+"h -f "+genie_setup 
    grid_command_options += " --OS=SL7 --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC --lines '+FERMIHTC_AutoRelease=True' "
    grid_command_options += "--lines '+FERMIHTC_GraceMemory="+str(GraceMemory)+"' --lines '+FERMIHTC_GraceLifetime="+str(GraceLifeTime)+"' "

    return grid_command_options 

def WriteXMLFile(commands_dict, start, end, jobs_dir, file_name='grid_submission.xml'):
    grid_file = jobs_dir+"/"+file_name
    if os.path.exists(grid_file):
        os.remove(grid_file)

    script = open( grid_file, 'w' ) 

    for id in range(start,end+1) :
        command_list = commands_dict[id]
        
        if len(command_list) == 1 : # serial
            script.write("<serial>\n")
            script.write(command_list[0]+"\n")
            script.write("</serial>\n")
        else : 
            script.write("<parallel>\n")
            for i in range(len(command_list)) : 
                script.write(command_list[i]+"\n")
            script.write("</parallel>\n")
                
    script.close()
    return grid_file

def WriteMainSubmissionFile(jobs_dir, genie_topdir, group, genie_setup='/src/scripts/production/python/setup_FNALGrid.sh', in_file_name='grid_submission.xml', expectedlife=60, out_file_name='fnal_dag_submit.fnal', memory=1, disk=20, jobs=1, role="Analysis"):

    fnal_file = jobs_dir+"/"+out_file_name
    if os.path.exists(fnal_file):
        os.remove(fnal_file)

    script = open( fnal_file, 'w' ) 
    script.write("#!/bin/bash\n")
    script.write("source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups ;\n")
    script.write("setup fife_utils ;\n")
    script.write("jobsub_submit_dag -g --group "+group+" --OS=SL7 --memory="+str(memory)+"GB --disk="+str(disk)+"GB --expected-lifetime="+str(expectedlife)+"h -N "+str(jobs)+" --role="+role+" -f "+genie_setup+" --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC  file://"+in_file_name+";" )
    
    return fnal_file
