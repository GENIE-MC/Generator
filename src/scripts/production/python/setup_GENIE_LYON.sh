#!/bin/bash 

# GENIE setup for the CC-IN2P3 (Lyon) worker nodes. Requires setup_LYON.sh to be sourced first.
# GENIE is not built here: the jobs use an existing GENIE build, compiled with the setup_LYON.sh environment
# Arguments: GENIE top dir (compiled), config dir (optional)

export GENIE=$1
GENIE_CONFIG_DIR=$2

if [ ! -x "$GENIE/bin/gevgen" ] ; then
    echo "No compiled GENIE found in '$GENIE'"
    return 1
fi
echo Using compiled GENIE from $GENIE

export GENIE_BIN=$GENIE/bin
export GENIE_LIB=$GENIE/lib
export PATH=$GENIE/bin:$PATH 
export LD_LIBRARY_PATH=$GENIE_LIB:$LD_LIBRARY_PATH

if [ ! -z "$GENIE_CONFIG_DIR" ] ; then
    echo Requested Genie config $GENIE_CONFIG_DIR 
    if [ -d "$GENIE_CONFIG_DIR" ]; then 
        export GALGCONF=$GENIE_CONFIG_DIR
    elif [ -d "$GENIE/$GENIE_CONFIG_DIR" ]; then
        export GALGCONF=$GENIE/$GENIE_CONFIG_DIR
    fi
    echo Genie config directory $GALGCONF
else
    # GENIE default configuration. Needed as the jobs use $GALGCONF/Messenger.xml
    export GALGCONF=$GENIE/config
fi
