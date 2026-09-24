#!/usr/bin/bash

# Robert Hatcher's colours
if [ -z "$PS1" ]; then
  # if $- contains "i" then interactive session
  export ESCCHAR="\x1B" # or \033 # Mac OS X bash doesn't support \e as esc?
  export OUTBLACK="${ESCCHAR}[0;30m"
  export OUTBLUE="${ESCCHAR}[0;34m"
  export OUTGREEN="${ESCCHAR}[0;32m"
  export OUTCYAN="${ESCCHAR}[0;36m"
  export OUTRED="${ESCCHAR}[0;31m"
  export OUTPURPLE="${ESCCHAR}[0;35m"
  export OUTORANGE="${ESCCHAR}[0;33m" # orange, more brownish?
  export OUTLTGRAY="${ESCCHAR}[0;37m"
  export OUTDKGRAY="${ESCCHAR}[1;30m"
  # labelled "light but appear in some cases to show as "bold"
  export OUTLTBLUE="${ESCCHAR}[1;34m"
  export OUTLTGREEN="${ESCCHAR}[1;32m"
  export OUTLTCYAN="${ESCCHAR}[1;36m"
  export OUTLTRED="${ESCCHAR}[1;31m"
  export OUTLTPURPLE="${ESCCHAR}[1;35m"
  export OUTYELLOW="${ESCCHAR}[1;33m"
  export OUTWHITE="${ESCCHAR}[1;37m"
  export OUTNOCOL="${ESCCHAR}[0m" # No Color
fi
# use as:   echo -e "${OUTRED} this is red ${OUTNOCOL}"

# Take in arguments and construct a DAG

export b0=$(basename $0)
export MAXCALLS=100
export DOREWRITE=0

usage() {
    echo -e "${OUTCYAN}"
    cat >&2 <<EOF
Submit a series of gmkspl calls to the grid, each run being one nucleon throw.

There are two primary steps. 
The first initialises a working area and writes preliminary config files.
These are then edited before jobs are submitted to be run.
The second step consists of running the produced DAG file.

This script handles the first step. Arguments:

 ${b0} 
       -T | --top         Top level directory to write output files to.
       -v | --version     GENIE version on UPS
       -q | --qualifier   GENIE qualifier on UPS
            --setup       Which UPS to setup with
	    --tune        GENIE tune
	    --genlist     Generator list (see EventGeneratorListAssembler.xml in \$GENIE for enum)
	    --rewrite     Whether to rewrite top level directory
       -e | --emax        Maximum energy in GeV
       -n | --nknots      N knots for gmkspl call
       -r | --nthrows     Number of times to call gmkspl. Max ${MAXCALLS} per job. 
                          Higher numbers will get distributed among separate jobs.
            --config      The config file containing points to sample.
            		  This will be copied to the 'cfg' directory in top-level output.

Use "${b0} -h" to show this help message.
EOF
    echo -e "${OUTNOCOL}"
    exit 0
}

process_args() {

    PRINTUSAGE=0
    
    TEMP=$(getopt -n $0 -s bash -a \
     --longoptions="help top: version: qualifier: setup: tune: genlist: xmlpath: \
     rewrite emax: nknots: nthrows: config:" \
     -o hT:v:q:p:t:e:n:r: -- "$@") || exit 1

    eval set -- "${TEMP}"
    unset TEMP

    let iarg=0
    set -u
    while [ $# -gt 0 ]; do
	let iarg=${iarg}+1
	case "$1" in
	    "--"             ) shift                      ; break ;;
	    -h | --help      ) PRINTUSAGE=1                       ;;
	    -T | --top       ) export OUTPUTTOP="$2"      ; shift ;;
	    -v | --version   ) export GENIEVERSION="$2"   ; shift ;;
	    -q | --qualifier ) export GENIEQUALIFIER="$2" ; shift ;;
	    --setup          ) export INITSETUPSTR="$2"   ; shift ;;
	    --tune           ) export GENIETUNE="$2"      ; shift ;;
	    --genlist        ) export GENLIST="$2"        ; shift ;;
	    --rewrite        ) export DOREWRITE=1                 ;;
	    -e | --emax      ) export EMAX="$2"           ; shift ;;
	    -n | --nknots    ) export NKNOTS="$2"         ; shift ;;
	    -r | --nthrows   ) export NTHROWSTOTAL="$2"   ; shift ;;
	    --config         ) export CONFIGFILE="$2"     ; shift ;;
	    -*               ) echo "unknown flag $opt ($1)" ; PRINTUSAGE=1 ;;
	esac
	shift # eat up the arg we just used
    done
    set +u

    if [[ -z ${GENIETUNE} ]] ; then
	echo -e "${OUTRED}ERROR:${OUTNOCOL} Please specify a GENIE tune using --tune."
	PRINTUSAGE=1
    fi

    if [[ ${PRINTUSAGE} -eq 1 ]] ; then
	usage
    fi

    echo -e "${OUTYELLOW}Using tune: ${GENIETUNE} with version, qualifier = ${GENIEVERSION}, ${GENIEQUALIFIER}${OUTNOCOL}"
    if [[ ! -z ${GENLIST} ]] && [[ ${GENLIST} != "Default" ]] ; then
	echo -e "${OUTLTPURPLE} Using generator list: ${GENLIST}${OUTNOCOL}"
    fi
}

process_args "$@"

# Check how many throws to do
export MULTIPLIER=1 # how many jobs for each point
export NTHROWS=${NTHROWSTOTAL}
if [[ ${NTHROWSTOTAL} -gt ${MAXCALLS} ]] ; then
    export MULTIPLIER=$(echo "(${NTHROWSTOTAL}-1) / ${MAXCALLS} + 1" | bc)
    export NTHROWS=${MAXCALLS}
    echo -e "${OUTCYAN}Passed ${NTHROWSTOTAL} repetitions of gmkspl per job, splitting into ${MULTIPLIER} copies of ${MAXCALLS}${OUTNOCOL}"
fi

if [ ! -d ${OUTPUTTOP} ] || [ ${DOREWRITE} -eq 1 ] ; then
    echo -e "${OUTGREEN}Making directory ${OUTPUTTOP} -- this is where your job and config files will live${OUTNOCOL}"
    if [[ -d ${OUTPUTTOP} ]] ; then
	rm -rf ${OUTPUTTOP}
    fi
    mkdir -p ${OUTPUTTOP}
    mkdir -p ${OUTPUTTOP}/cfg
    mkdir -p ${OUTPUTTOP}/bin
    mkdir -p ${OUTPUTTOP}/work-products
elif [ -d ${OUTPUTTOP} ] && [ ${DOREWRITE} -eq 0 ] ; then
    echo -e "${OUTRED}Directory ${OUTPUTTOP} exists, not overwriting. Pass --rewrite if you want this to be rewritten.${OUTNOCOL}"
    exit 1
fi

# Find the number of points to iterate over.
export NLINES=$(awk '!/^[[:space:]]*($|#)/' ${CONFIGFILE} | wc -l)
export LASTPROCESS=$(echo "${NLINES} * ${MULTIPLIER} - 1" | bc) # minus one cause zero indexed
# Copy the config file over
cp ${CONFIGFILE} ${OUTPUTTOP}/cfg

# which CVMFS?
CVMFS_SETUP=
case ${INITSETUPSTR} in
    "uboone") CVMFS_SETUP=/cvmfs/uboone.opensciencegrid.org/products/setup_uboone_mcc9.sh   ;;
    "sbnd")   CVMFS_SETUP=/cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh       ;;
    "icarus") CVMFS_SETUP=/cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh ;;
    "dune")   CVMFS_SETUP=/cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh       ;;
    *) echo "Unknown setup script, exiting" ; exit 1 ;;
esac

# write the binary file
export EXECFILE=${OUTPUTTOP}/bin/make_splines_ZExp.sh
cat > ${EXECFILE} <<EOF
define_cfg()
{
  # This configuration is for:
  # GENIEVERSION   = "${GENIEVERSION}"
  # GENIEQUALIFIER = "${GENIEQUALIFIER}"
  # setup from ${CVMFS_SETUP}

  source ${CVMFS_SETUP}
  setup genie ${GENIEVERSION} -q ${GENIEQUALIFIER}
  setup ifdhc v2_8_0 -q e26:p3915:prof
  setup ifdhc_config v2_8_0
  ups active

  # details of the actual jobs
  export JOBOUTPUTTOP=${OUTPUTTOP}
  export JOBGENIETUNE=${GENIETUNE}
  export JOBGENLIST=${GENLIST}
  export JOBXMLPATH=\${WORKDIR}/cfg
  export JOBEMAX=${EMAX}
  export JOBNKNOTS=${NKNOTS}
  export JOBNTHROWS=${NTHROWS}
  export JOBCONFIG=$(basename ${CONFIGFILE})

  export MULTIPLIER=\$1
  export SUBPROCESS=\$2

  # this would be in a loop
  # I submitted N jobs total, so SUBPROCESS goes from 0 to N-1	
  # These jobs are grouped into M groups of K jobs per group, that are at the same point.
  # MULTIPLIER = K

  # Within the group, to avoid duplication, we need to know which seeds to use.
  # SUBPROCESS // MULTIPLIER goes from 0 to M-1
  # SUBPROCESS %  MULTIPLIER goes from 0 to K-1
  # So SUBPROCESS % MULTIPLIER indexes where along the same group we are and defines the seed

  export JOBGROUP=\$(echo "\${SUBPROCESS} / \${MULTIPLIER}" | bc)
  export JOBINDEX=\$(echo "\${SUBPROCESS} % \${MULTIPLIER}" | bc)
  export FIRSTSEED=\$(echo "\${JOBNTHROWS} * \${JOBINDEX}" | bc)
  export LASTSEED=\$(echo "\${FIRSTSEED} + \${JOBNTHROWS} - 1" | bc)

  # Fetch config files
  mkdir -p \${JOBXMLPATH}
  ifdh cp -r -D ${OUTPUTTOP}/cfg/ \${JOBXMLPATH}/

  # Read the config file and stop at the correct line (which JOBGROUP ?)
  # skip any comment lines (allowing for whitespace) and empty lines
  CFGLINE=\$(awk -v target=\$(echo "\${JOBGROUP}+1" | bc) '!/^[[:space:]]*($|#)/ { if (++n == target) { print; exit } }' \${JOBXMLPATH}/\${JOBCONFIG})
  echo Selected line \${CFGLINE}

  # Read the neutrino and target PDGs from the config line
  # config should be nupdg tgtpdg
  export JOBNUPDG=\$(echo \${CFGLINE} | awk -F " " '{print \$1}')
  export JOBTGTPDG=\$(echo \${CFGLINE} | awk -F " " '{print \$2}')
  export JOBPRETTYNU=\$(echo \${CFGLINE} | awk -F " " '{print \$3}')
  export JOBPRETTYTGT=\$(echo \${CFGLINE} | awk -F " " '{print \$4}')
}

setup_python()
{
  # Setup a virtual environment
  python3 -m venv .
  source bin/activate
  # Pull in pip, numpy, pandas, fastparquet
  python3 -m pip install numpy
  python3 -m pip install pandas==2.2.3 # newer than that and meson seems to complain
  python3 -m pip install fastparquet "pandas<2.3"
  python3 -m pip install  tqdm colorama
}

export BASEDIR=\$(pwd)
export WORKDIR=\$(mktemp -d -p \${BASEDIR})
mkdir -p \$WORKDIR && cd \$WORKDIR

define_cfg \$1 \$2

# Prepend to GXMLPATH
export GXMLPATH=\${JOBXMLPATH}:\${GXMLPATH}
# copy in NewQELXSec and specify 1 nucleon
awk '{sub(/NumNucleonThrows"> *[0-9]+/, "NumNucleonThrows\"> 1")}1' \${GENIE}/config/NewQELXSec.xml > \${JOBXMLPATH}/NewQELXSec.xml
# (temporarily) copy in average_splines.py from GPVMs... This should live in tagged contrib
ifdh cp -D /exp/sbnd/app/users/kplows/gen_workshop/BuildEventGenerators/production/splines-interpolation/test_configs/average_splines.py \${JOBXMLPATH}/

# build the command
GENIECMD="gmkspl"
GENIEARGS=" -p \${JOBNUPDG} -t \${JOBTGTPDG} -e \${JOBEMAX} -n \${JOBNKNOTS} --tune \${JOBGENIETUNE}"
if [[ \${JOBGENLIST} != "Default" ]] && [[ ! -z \${JOBGENLIST} ]] ; then
  GENIEARGS=\${GENIEARGS}" --event-generator-list \${JOBGENLIST}"
fi

mkdir -p \${WORKDIR}/xmls
for i in \$(seq \${FIRSTSEED} \${LASTSEED}) ; do
    fullcmd=\${GENIECMD}\${GENIEARGS}" --seed \${i} -o \${WORKDIR}/xmls/seed-\${i}.xml"
    echo \$fullcmd
    eval \$fullcmd 2>&1 | tee -a run-\${SUBPROCESS}.log
done

# Set up python and run the average_splines script
setup_python
export OUTNAME=throws_\${JOBPRETTYNU}_\${JOBPRETTYTGT}_idx\${JOBINDEX}
python3 \${JOBXMLPATH}/average_splines.py -b \${WORKDIR}/xmls -o \${WORKDIR}/\${OUTNAME}
#python3 \${GENIE}/src/contrib/splines/average_splines.py -b \${WORKDIR}/xmls -o \${WORKDIR}/\${OUTNAME}

# ifdh cp the resulting files
ifdh cp -D \$(find \${WORKDIR} -iname "\${OUTNAME}.parquet") ${OUTPUTTOP}/work-products/
ifdh cp -D \$(find \${WORKDIR} -iname "\${OUTNAME}.xml") ${OUTPUTTOP}/work-products/
#ifdh cp -D \$(find \${WORKDIR} -iname "run-\${SUBPROCESS}.log") ${OUTPUTTOP}/work-products/

# cleanup
cd \${BASEDIR} && rm -rf \${WORKDIR}/
EOF
chmod u+x ${EXECFILE}

# Now write a DAG file.
JOBSUBMAIN="jobsub_submit -n -G $(id -ng)  --resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE  --append_condor_requirements='(TARGET.HAS_CVMFS_sbn_opensciencegrid_org==true)' --singularity-image=/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest"
RESOURCES="--expected-lifetime 8h --disk 8GB --memory 2GB"
JOBSUBFULL=${JOBSUBMAIN}" "${RESOURCES}
echo -e "<parallel>" > genie-splines.dag
for i in $(seq 0 ${LASTPROCESS}) ; do
    echo -e "${JOBSUBFULL} file://${EXECFILE} ${MULTIPLIER} ${i}" >> genie-splines.dag
done
echo -e "</parallel>" >> genie-splines.dag

mv genie-splines.dag ${OUTPUTTOP}/cfg/

echo -e "${OUTCYAN}Exit status 0${OUTNOCOL}"
