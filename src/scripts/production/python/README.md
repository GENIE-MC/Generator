# Grid Submitters 
In this directory you can find a serie of python scripts that facilitate to submit GENIE jobs in Grids (such as the FNAL Grid). The Sumitter scripts are responsible to launch the correct commands to the nodes. In order to do so, it calls a number of functions in the expected order and writes a `fnal_dag_submit.fnal` file, which is sourced automatically by the code in order to submit your jobs to the file. The submitters are also responsible to create the directory substructure where your final files will be. It also creates the scripts that will be launched by `fnal_dag_submit.fnal`. In order to use this script you need to have a 'GENIE' enviromental variable. To run jobs from the FNAL grid, you must run from the pnfs area. If these two conditions are not satisfied the code will exit.

## Possible options
An example of a Submitter code is `eAScatteringGridSubmitter.py`. This code will launch the splines (for free nucleon and free nuclei) and electron scattering events. The options are: 
- Version: Genie version number. Default `master`
- Cycle : defult 01
- Production: default routine_validation
- Grid System: for now only FNAL available
- Grid group: you can specify your grid group here. I.E for FNAL: genie, dune,...
- Software top dir: base directory for all your GENIE code. 
- GENIE topdir: path to your GENIE/Generator code
- Jobs top dir: directory where you want to have your output files and submission commands. For FNAL it must be in pnfs.
- Source top dir: Jobs topdir used as a source for missing cross section splines. It must have the same directory structure (see below). 
- Config dir: GENIE configuration directory. Useful when your configuration files are different from the ones in github. 
- Nu list: list of neutrinos **and electrons** that you want to use. Separate them by a comma. In eA Submitter code it is defaulted to electron.
- Tgt list: comma separated list of target pdg codes
- Gen list: Comma separated list of event generator list to be used for all splines
- vN and vA generator list: same as above but different for vN and vA splines
- event-generator-list: to be used for event generation.
- ntotevents: total number of events to run
- nmaxevents: maximum number of events to run per job
- energy : comma separated list of beam energy for electrons. Monoenergetic beam
- starting-point: 0) Free nucleon splines, 1) combine free nucl splines, 2) Compound nuclei splines, 3) Combine compound nuclei splines, 4) Event Production
- stopping point: Same as above
- tune: tune to be used for spline and event generation
- submit-jobs: this is defaulted to false so that the user can inspect and double check the submission files. However if you want to just launch them, use this option to automatically run the whole chain.

## Jobs top dir final structure
The directory generated with the default parameters is going to have the following format:
fnal_dag_submit.fnal  master-routine_validation_01-eScattering/
grid_submission.xml   master-routine_validation_01-xsec_vA/
group_vA.sh           master-routine_validation_01-xsec_vN/
group_vN.sh           setup_FNALGrid.sh

## How to run the scripts
To submit jobs to run electrons on carbon and oxigen simply do:
`python eAScatteringGridSubmitter.py --nu-list 11 --tgt-list 1000060120,1000080160 --energy 1.6,2.2,3 --submit-jobs` 

If you only want to generate events (you already have the splines), do:
`python eAScatteringGridSubmitter.py --nu-list 11 --tgt-list 1000060120,1000080160 --energy 1.6,2.2,3 --submit-jobs --starting-point 4 --source-prod-dir MySplines/` 

and MySplines put your spline (named total_xsec.xml) in MySplines/master-routine_validation_01-xsec_vA/total_xsec.xml

