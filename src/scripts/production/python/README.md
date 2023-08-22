# Grid Submitters 
In this directory you can find a serie of python scripts that facilitate to submit GENIE jobs in Grids (such as the FNAL Grid). The Sumitter scripts are responsible to launch the correct commands to the nodes. In order to do so, it calls a number of functions in the expected order and writes a `fnal_dag_submit.fnal` file. The submitters are also responsible to create the directory substructure where your final files will be. It also creates the scripts that will be launched by `fnal_dag_submit.fnal`. In order to use this script you need to have a 'GENIE' enviromental variable. To run jobs from the FNAL grid, you must run from the pnfs area. If these two conditions are not satisfied the code will exit.

Two submitters are availalbe: 
- `XMLGridSubmitter.py`: submits the jobs to obtain splines for any lepton
- `eAScatteringGridSubmitter.py`: submits the jobs to create splines and generate e-A events
other scripts can be easily added. Does not support neutrino scattering event generation
- `vAScatteringGridSubmitter.py`: submits the jobs to create splines and generate nu-A events
other scripts can be easily added. Does not support electron scattering event generation

NOTE: avoid using the scripts as a black box. Job lifetime and memory should be optimized for each job. Runing with the default might bring your jobs to hold. 

## Requirements
- You must have the $GENIE enviromental variable to your GENIE-MC/Generator directory. This is required only to run the python scripts, no previous compilation of the GENIE-MC/Generator code is required.
- Your jobs-topdir must be in the pnfs area

## Possible options
An example of a Submitter code is `eAScatteringGridSubmitter.py`. This code will launch the splines (for free nucleon and free nuclei) and electron scattering events. The options are: 

- version: Genie version number. Default `master`.
- git-location: link to the github repository from which to clone GENIE-Generator. By default https://github.com/GENIE-MC/Generator is used. 
- git-branch: This is the name of the git branch you want to pull from the GENIE Generator repository. Default is master
- cycle : defult 01
- production: default routine_validation
- grid-system: for now only FNAL available
- grid-group: you can specify your grid group here. I.E for FNAL: genie, dune,...
- softw-topdir: base directory for all your GENIE code. 
- genie-topdir: path to your GENIE/Generator code
- jobs-topdir: directory where you want to have your output files and submission commands. For FNAL it must be in pnfs.
- source-prod-dir: Jobs topdir used as a source for missing cross section splines. It must have the same directory structure (see below). 
- config-dir: configuration directory containing alternative xml files to use for model configuration
- probe-list: list of neutrinos **and electrons** that you want to use. Separate them by a comma. In eA Submitter code it is defaulted to electron.
- nu-tgt-list: comma separated list of target pdg codes for the generation of neutrino splines and events
- e-tgt-list: comma separated list of target pdg codes for the generation of electron splines and events
- vN and vA generator list: Comma separated list of event generator list to be used for all splines for vN and vA splines
- event-generator-list: to be used for event generation.
- nu-ntotevents: total number of events to run for neutrinos
- e-ntotevents: total number of events to run for electrons
- nmaxevents: maximum number of events to run per job
- ebeam-energy : comma separated list of beam energy for electrons. Monoenergetic beam
- starting-point: 0) Free nucleon splines, 1) combine free nucl splines, 2) Compound nuclei splines, 3) Combine compound nuclei splines, 4) Event Production
- stopping point: Same as above
- tune: tune to be used for spline and event generation
- submit-jobs: this is defaulted to false so that the user can inspect and double check the submission files. However if you want to just launch them, use this option to automatically run the whole chain.

## Jobs top dir final structure
The directory generated with the default parameters is going to have the following format:

- fnal_dag_submit.fnal: it contains the main jobsubmit command responsible to launch all the jobs. You can source it automatically with the `--submit-jobs` option. Otherwise, you can source it manually. 
- grid_submission.xml : it contains the list of commands to be submitted by the fnal_dag_submit script. It deals with parallel jobs. 
- master-routine_validation_01-xsec_vN/: this directory contains the scripts to run neutrino(electron)-nucleon splines on parallel. 
- master-routine_validation_01-xsec_vA/: this directory contains the scripts to run neutrino(electron)-nuclei splines from the neutrino(electron)-nucleon splines.
- master-routine_validation_01-eScattering/ : this directory contains the scripts to run events with GENIE and it will contain the output files with your events. 
- setup_FNAL.sh : this setup will be sourced in each node. It sets the requirements needed to run the jobs at the FNAL grid. 
- setup_GENIE.sh : this setup script deals with the GENIE configuration. It will clone the code from `https://github.com/GENIE-MC/Generator` and build the code in every node. It also deals with the branch version and configuration directory setup, in case you want a different configuration for your tunes. 

## How to run the scripts
To submit jobs to run electrons on carbon and oxigen simply do:
`python eAScatteringGridSubmitter.py --probe-list 11 --tgt-list 1000060120,1000080160 --energy 1.6,2.2,3 --submit-jobs --jobs-topdir <jobs-topdir>`
The flag --submit-jobs will indicate that we want to submit the jobs directly. It is recommended to not use it and inspect the files carefuly before submission. To submit without this option, simply do 'source <jobs-topdir>/fnal_dag_submit.fnal', where <jobs-topdir> is the specified jobs directory in the /pnfs/ area.  

If you only want to generate events (you already have the splines), do:
`python eAScatteringGridSubmitter.py --probe-list 11 --tgt-list 1000060120,1000080160 --energy 1.6,2.2,3 --submit-jobs --starting-point 4 --source-prod-dir MySplines/ --jobs-topdir <jobs-topdir>` 

and MySplines put your spline (named total_xsec.xml) in MySplines/master-routine_validation_01-xsec_vA/total_xsec.xml

If you have a specific github version to use, specify the git branch name:
`python eAScatteringGridSubmitter.py --probe-list 11 --tgt-list 1000060120,1000080160 --energy 1.6,2.2,3 --submit-jobs --starting-point 4 --source-prod-dir MySplines/ --git-branch my_branch_name --jobs-topdir <jobs-topdir>`


In some cases, the user might desire to change the model configuration. For instance, we might want to understand the effect of FSI in the model. One could opt to turn FSI effects off. For the G18_10a_02_11b CMC follow the following steps:
   - In your jobs-topdir create a conf/ directory
   - Copy the $GENIE/config/G18_10a/ModelConfiguration.xml in your <jobs-topdir>/conf directory
   - Turn FSI off by changing the <param type="bool" name="HadronTransp-Enable"> true </param> flag to false. 
   - Add the option --conf-dir <jobs-topdir>/conf in your main submission script option

For instance: 
`python eAScatteringGridSubmitter.py --probe-list 11 --tgt-list 1000060120,1000080160 --energy 1.6,2.2,3 --submit-jobs --starting-point 4 --source-prod-dir MySplines/ --git-branch my_branch_name --jobs-topdir <jobs-topdir> --tune G18_10a_02_11b --conf-dir <jobs-topdir>/conf`

