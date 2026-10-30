#! /usr/bin/env python

"""\
Utilities to run the GENIE production scripts at the CC-IN2P3 computing centre (Lyon).
CC-IN2P3 uses the Slurm batch system and the job top directory is expected to be in a
shared area mounted on the worker nodes (i.e. /sps/<group>/...). Contrary to FNAL, GENIE
is not built in each job: the jobs use the compiled GENIE given by --genie-topdir. Each job runs in its
local scratch area, copies the input files there, and copies the output files back to
the jobs directory, mirroring what is done with ifdh at FNAL.

The FNAL DAG (serial/parallel xml) is replaced by a bash script that submits each stage
with sbatch and makes it depend on the successful completion of the previous stage.

Author:
      Julia Tena Vidal <jtenavidal \st tauex.tau.ac.il>
      Tel Aviv University
Copyright:
   Copyright (c) 2003-2026, The GENIE Collaboration
   For the full text of the license visit http://copyright.genie-mc.org
"""
import os, glob, re

# Working directory of each job on the worker node. It is set by the job script itself
JOB_DIR = "$GENIE_JOB_DIR"

def CreateShellScript ( commands , jobs_dir, shell_name, out_files, grid_setup, genie_setup, conf_dir, in_files, genie_topdir, configure_INCL, configure_G4, GHEPMC3Output=False ) :
    # GENIE is never built on the worker nodes: genie_topdir is an existing GENIE build accessible from them.
    # It must have been compiled with the options needed by the jobs (Geant4, HepMC3, ...)
    if configure_INCL :
        print( "INCL is not available at LYON. Abort ..." )
        exit()

    shell_file = jobs_dir+"/"+shell_name+".sh"

    if os.path.exists(shell_file):
        os.remove(shell_file)

    script = open( shell_file, 'w' )
    script.write("#!/bin/bash \n")
    script.write("export GENIE_JOB_DIR=${TMPDIR:-/tmp}/genie_${SLURM_JOB_ID:-$$} ;\n")
    script.write("mkdir -p $GENIE_JOB_DIR ;\n")
    script.write("trap 'rm -rf $GENIE_JOB_DIR' EXIT ;\n")
    script.write("cd $GENIE_JOB_DIR ;\n")
    script.write("source "+grid_setup+" ; \n")
    if conf_dir != '' :
        conf_files = glob.glob(conf_dir+"/*.xml")
        script.write("mkdir $GENIE_JOB_DIR/conf ;\n")
        for conf_i in conf_files :
            script.write("cp "+conf_i+"  $GENIE_JOB_DIR/conf ;\n")
        conf_dir = "$GENIE_JOB_DIR/conf"
    script.write("source "+genie_setup+" "+genie_topdir+" \""+conf_dir+"\" || exit 1 ;\n")
    script.write("cd $GENIE_JOB_DIR ;\n")
    # From here on any failure must make the job fail, so that the dependent jobs are not run
    script.write("set -e ;\n")

    if isinstance(in_files, list) :
        for i in range(len(in_files)):
            script.write("cp "+in_files[i]+"  $GENIE_JOB_DIR ;\n")
    else :
        if len(in_files) != 0 :
            script.write("cp "+in_files+"  $GENIE_JOB_DIR ;\n")

    if isinstance(commands,list):
        for command in commands :
            script.write(command+" ;\n")
    else :
        script.write(commands+" ;\n")

    if isinstance(out_files, list) :
        for file_name in out_files :
            script.write("cp "+file_name+" "+jobs_dir+" ;\n")
    else :
        script.write("cp "+out_files+" "+jobs_dir+" ;\n")

    script.close()
    os.chmod(shell_file, 0o755)
    return shell_file

def SlurmMemory( memory ) :
    # FNAL style units (1GB, 500MB) to Slurm units (1G, 500M)
    memory = str(memory).strip().upper()
    if memory.isdigit() :
        return memory+"M"
    return re.sub(r'B$', '', memory)

def LyonShellCommands(group, jobs_dir, job_name, hours = 10, memory="1GB", disk="500MB", partition="htc"):
    # disk is not requested: the CC-IN2P3 Slurm configuration rejects --tmp
    grid_command_options = " --job-name="+job_name+" --partition="+partition+" --ntasks=1 --time="+str(int(float(hours)))+":00:00 "
    grid_command_options += " --mem="+SlurmMemory(memory)
    grid_command_options += " --output="+jobs_dir+"/"+job_name+".%j.log "
    if group :
        grid_command_options += " --account="+group
    # Jobs reading from or writing to /sps must declare it
    if jobs_dir.startswith('/sps') :
        grid_command_options += " --licenses=sps"

    return grid_command_options

def SubmitCommand(group, jobs_dir, shell_file, hours = 10, memory="1GB", disk="500MB"):
    job_name = os.path.basename(shell_file)[:-3]
    return "sbatch "+LyonShellCommands(group, jobs_dir, job_name, hours, memory, disk)+" "+shell_file

def WriteMainSubmissionFile(commands_dict, start, end, jobs_dir, max_parallel=0, out_file_name='lyon_slurm_submit.sh'):
    # Equivalent of the FNAL dag: every stage only starts if all the jobs of the previous stage succeeded.
    # At most max_parallel jobs of each stage run at the same time (0: no limit): job k of a stage starts once
    # job k-max_parallel of the same stage has ended
    lyon_file = jobs_dir+"/"+out_file_name
    if os.path.exists(lyon_file):
        os.remove(lyon_file)

    script = open( lyon_file, 'w' )
    # The script can be run with bash or sourced, so it must not exit the shell (no set -e).
    # If a submission fails, the jobs already submitted are cancelled and nothing else is submitted
    script.write("#!/bin/bash\n")
    script.write("# Submits the GENIE jobs to the CC-IN2P3 Slurm batch system\n")
    script.write("# At most MAX_PARALLEL_JOBS jobs of each stage run at the same time (0: no limit).\n")
    script.write("# It can be changed when submitting: MAX_PARALLEL_JOBS=100 bash "+out_file_name+"\n")
    script.write("MAX_PARALLEL_JOBS=${MAX_PARALLEL_JOBS:-"+str(int(max_parallel))+"}\n")
    script.write("submit() {\n")
    script.write("  local jid deps=\"\"\n")
    script.write("  [ -z \"$SUBMIT_FAILED\" ] || return 1\n")
    script.write("  [ -z \"$PREV_IDS\" ] || deps=\"afterok$PREV_IDS\"\n")
    script.write("  if [ \"$MAX_PARALLEL_JOBS\" -gt 0 ] && [ ${#STAGE_IDS[@]} -ge \"$MAX_PARALLEL_JOBS\" ] ; then\n")
    script.write("    deps=\"${deps:+$deps,}afterany:${STAGE_IDS[${#STAGE_IDS[@]}-MAX_PARALLEL_JOBS]}\"\n")
    script.write("  fi\n")
    script.write("  if [ -n \"$deps\" ] ; then\n")
    script.write("    jid=$(sbatch --parsable --kill-on-invalid-dep=yes --dependency=$deps \"$@\")\n")
    script.write("  else\n")
    script.write("    jid=$(sbatch --parsable \"$@\")\n")
    script.write("  fi\n")
    script.write("  if [ $? -ne 0 ] || [ -z \"$jid\" ] ; then\n")
    script.write("    SUBMIT_FAILED=1\n")
    script.write("    return 1\n")
    script.write("  fi\n")
    script.write("  jid=${jid%%;*}\n")
    script.write("  CURR_IDS=\"$CURR_IDS:$jid\"\n")
    script.write("  SUBMITTED_IDS=\"$SUBMITTED_IDS $jid\"\n")
    script.write("  STAGE_IDS+=($jid)\n")
    script.write("}\n")
    script.write("cancel_submitted() {\n")
    script.write("  echo \"Job submission failed. Cancelling the jobs already submitted:$SUBMITTED_IDS\"\n")
    script.write("  [ -z \"$SUBMITTED_IDS\" ] || scancel $SUBMITTED_IDS\n")
    script.write("}\n")
    script.write("PREV_IDS=\"\"\n")
    script.write("SUBMITTED_IDS=\"\"\n")
    script.write("SUBMIT_FAILED=\"\"\n")

    for id in range(start,end+1) :
        if id not in commands_dict :
            continue
        command_list = commands_dict[id]
        if len(command_list) == 0 :
            continue
        script.write("\n# Stage "+str(id)+": "+str(len(command_list))+" job(s)\n")
        script.write("CURR_IDS=\"\"\n")
        script.write("STAGE_IDS=()\n")
        for command in command_list :
            script.write(re.sub(r'^sbatch ', 'submit ', command)+"\n")
        script.write("if [ -n \"$SUBMIT_FAILED\" ] ; then cancel_submitted ; return 1 2> /dev/null || exit 1 ; fi\n")
        script.write("echo \"Stage "+str(id)+" submitted: ${CURR_IDS#:}\"\n")
        script.write("PREV_IDS=$CURR_IDS\n")

    script.close()
    os.chmod(lyon_file, 0o755)
    return lyon_file
