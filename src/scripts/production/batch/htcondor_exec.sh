#
# Simple shell script used as 'Executable' with HTCondor batch submission scripts.
# Configures GENIE, cd's to the job directory and executes the input GENIE commands.
#
# Syntax:
#   ./htcondor_exec.sh [/path/genie_setup_script] [/job/dir] [executable --argument1  --argument2 ...] [another_executable --argument1  --argument2 ...] ...
#
# C.Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
#

#!/bin/bash

#echo "Number of arguments: ${#@}"

#
# first two arguments are the setup script and the job directory
#

genie_setup=$1
jobs_dir=$2

echo "htcondor_exec.sh / GENIE setup script: $genie_setup" 
echo "htcondor_exec.sh / Job directory: $jobs_dir"

source $genie_setup
cd $jobs_dir

shift
shift

#
# remaining arguments are genie executables and their arguments
#

# Known commands to look
known_exec=(
  'gevcomp' 
  'gevdump' 
  'gevgen' 
  'gevgen_atmo' 
  'gevgen_hadron' 
  'gevgen_ndcy' 
  'gevgen_t2k' 
  'gevpick' 
  'gevscan' 
  'gmkspl' 
  'gmxpl' 
  'gntpc' 
  'gspl2root' 
  'gspladd' 
  'gxscomp'
 )

#echo "Number of known commands: ${#known_exec[@]}"

genie_command=""
genie_args=""

do_run=0

# loop over arguments
for x in $@
do
	new_command=0

	# loop over known GENIE commands
	for y in ${known_exec[@]}
	do
		# current argument is a known GENIE command
    		if [ "$x" ==  "$y" ]
		then
		   #		   
		   if [ "$do_run" -eq 1 ]
		   then
			echo "htcondor_exec.sh / Executing:  ${genie_command} ${genie_args}"
			${genie_command} ${genie_args}

		   fi
		   new_command=1	   
                   do_run=1
		   genie_command="$x"	
                   genie_args=""
		   break
		fi
	done

	if [ "$new_command" -eq 1 ] 
	then
		continue
	fi

        genie_args="${genie_args} $x";
done

echo "htcondor_exec.sh / Executing:  ${genie_command} ${genie_args}"
${genie_command} ${genie_args}

