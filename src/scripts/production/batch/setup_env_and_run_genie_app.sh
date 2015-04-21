#
# Simple shell script used as 'Executable' with HTCondor
#
# Syntax:
#   ./setup_env_and_run_genie_app.sh [/path/genie_setup_script] [/job/dir] [executable --argument1  --argument2 ...]
#
# C.Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

#!/bin/bash

echo "Running setup_env_and_run_genie_app.sh: "

genie_setup=$1
jobs_dir=$2
executable=$3

shift
shift
shift

args=$@

echo " - GENIE setup script: $genie_setup" 
echo " - Job directory: $jobs_dir"
echo " - Executable: $executable"
echo " - Arguments: $args"

source $genie_setup
cd $jobs_dir
$executable $@

