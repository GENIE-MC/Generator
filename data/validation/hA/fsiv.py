#!/usr/bin/env python
"""
    Author: Tomasz Golan

    For help run one from the following commands:
	'./fsiv.py --help' or './fsiv.py -h'
	'python fsiv.py --help' or 'python fsiv.py -h'
"""

import os, re, sys, getopt, time
from subprocess import call

#GLOBAL VARIABLES

DAT_DIR = os.environ.get('FSIV_DAT_DIR') or 'data_files'	#set a directory with data files (use 'pwd/data' if $FSIV_DAT_DIR is not defined)
SIM_DIR = os.environ.get('FSIV_SIM_DIR') or 'sim_files'		#set a directory for root/xsec files (use 'pwd/sim_files' if $FSIV_SIM_DIR is not defined)
PNG_DIR = os.environ.get('FSIV_PNG_DIR') or 'png_files'		#set a directory for plots (use 'pwd/png_files' if $FSIV_PNG_DIR is not defined)
LOG_DIR = os.environ.get('FSIV_LOG_DIR') or 'log_files'		#set a direcotry for log files (use 'pwd/log_files' if $FSIV_LOG_DIR is not defined)   
BSH_DIR = os.environ.get('FSIV_BSH_DIR') or 'bash_scripts'	#set a direcotry for bash scripts (use 'pwd/bash_files' if $FSIV_BSH_DIR is not defined)   

authors = []		#storage for the list of available authors
tot_proc_list = []	#storage for the list of available processes for total cross section

date_ = time.strftime("%b") + "_" + time.strftime("%d") + "_" + time.strftime("%y")

#CLASSES

class Author:
    #store available processes for given author
    def __init__ (self, name): 
	self.name = name								#author's name
	self.pion = self.nucl = self.kaon = self.tot = self.nrg = self.ang = False	#author's switchers for simluations/plots to be done
	self.tot_proc = []								#storage for processes for total cross section data ([particle, tk_min, tk_max, target, fate])
    def __cmp__(self, name): return cmp(self.name, name)				#compare authors by name
    def __repr__(self): return repr(self.name)
    def print_ (self, dest = None):							#print switchers for given author
	print >> dest, '-'*22
	print >> dest, self.name,':'
	print >> dest, '-'*22
	print >> dest, 'Pion mode:\t',self.pion
	print >> dest, 'Nucleon mode:\t',self.nucl
	print >> dest, 'Kaon mode:\t',self.kaon
	print >> dest, 'Xsec data:\t',self.tot
	print >> dest, 'Energy data:\t',self.nrg
	print >> dest, 'Angle data:\t',self.ang
    def check (self, particle):								#return particle switcher for given particle
	if todo_authors and self.name not in todo_authors: return False
	if (particle == 'pion'): return self.pion
	elif (particle == 'nucleon'): return self.nucl
	elif (particle == 'kaon'): return self.kaon

class TotalProcess:
    #store information about particle, target and energy range for total cross section
    def __init__ (self, proc):					#create proccess from proc = [particle, [tk list], xsec_max, target, fate]
	self.particle = proc[0]
	if self.particle in ['p','n']: self.particle_kind = 'nucleon'
	elif self.particle in ['kp', 'km']: self.particle_kind = 'kaon'
	else: self.particle_kind = 'pion'
	self.tk_data = set(proc[1])
	self.tk_max = [max(proc[1])]
	self.xsec_max = [proc[2]]
	self.target = proc[3]
	self.fate = [proc[4]]
	self.name = self.particle + '_' + self.target
    def __cmp__(self, x): return cmp(self.name, x)		#compare processes by name = particle_target
    def __repr__(self): return repr(self.name)
    def update (self,x):					#update process information
	self.tk_data.update(x.tk_data)
	if x.fate[0] not in self.fate:				#add new fate (and new xsec max) if necessary
	    self.fate.append(x.fate[0])
	    self.xsec_max.append(x.xsec_max[0])
	    self.tk_max.append(x.tk_max[0])
	else:							#if fate exists -> update xsec max and energy max [for plot range]
	    index = self.fate.index(x.fate[0])
	    if x.xsec_max[0] > self.xsec_max[index]: self.xsec_max[index] = x.xsec_max[0]
	    if x.tk_max[0] > self.tk_max[index]: self.tk_max[index] = x.tk_max[0]
    def print_ (self, dest = None):
	print >> dest, ' * ', self.particle, '+', self.target
	print >> dest, '\tTk found in data = ', str(sorted(self.tk_data))
	print >> dest, '\tTk prepared for simulation = ', str(sorted(self.tk_sim))
	print >> dest, '\tFates:'
	for f in self.fate:
	    index = self.fate.index(f)
	    print >> dest, '\t\t', f, '\t-> ', 'Tk max = ', self.tk_max[index], ', xsec max = ', self.xsec_max[index]

#FUNCTIONS
	
def help():
    print '\n','*'*110
    print """
fsiv.py runs 'runfast.pl' and 'intranukeplotter.pl' in a loop (assuming they are in the same directory)

the data are taken from the directory defined in $FSIV_DAT_DIR environmental variable
if $FSIV_DAT_DIR is not defined 'pwd/data_files' directory is taken instead

the output root/xsec files are saved into the directory defined in $FSIV_SIM_DIR environmental variable
if $FSIV_SIM_DIR is not defined 'pwd/sim_files' is used

the output png files are saved into the directory defined in $FSIV_PNG_DIR environmental variable
if $FSIV_PNG_DIR is not defined 'pwd/png_files' is used

the bash scripts are saved into the directory defined in $FSIV_BSH_DIR environmental variable
if $FSIV_BSH_DIR is not defined 'pwd/bash_scripts' is used

the log files are saved into the directory defined in $FSIV_LOG_DIR environmental variable
if $FSIV_LOG_DIR is not defined 'pwd/log_files' is used
    
by default all available data sets are done unless otherwise is specified

available options are:

    -p, --pion    \t -> turn on pion validation
    -n, --nucleon \t -> turn on nucleon validation
    -k, --kaon    \t -> turn on kaon validation

    *** if none from aboves is chosen they all are set on ***

    -t, --total   \t -> turn on total cross section validation
    -e, --energy  \t -> turn on dsigma / denergy validation
    -a, --angle   \t -> turn on dsigma / dangle validation

    *** if none from aboves is chosen they all are set on ***

    -s, --simulation \t -> turn on 'runfast.pl' to generate root files
    -f, --plot       \t -> turn on 'intrenukeplotter.pl' to generate plots

    *** if none from aboves is chosen they all are set on ***

    --nof_events_diff=   -> number of events for differential cross section (default=1,000,000)
    --nof_events_total=  -> number of events for total cross section (default=100,000)
    
    --command=  \t -> command
    --FNAL	\t -> set command for FNAL grid (i.e. jobsub -e GENIE -e LD_LIBRARY_PATH -e PATH)

    --authors=  \t -> list of authors to be done (only for differential xsec)
    --isospins= \t -> list of particles to be done (only for total xsec): pip,pim,p,n,kp,km

Examples:

    1. run simulations and create plots for all available pion data sets for total cross section

	./fsiv.py --pion --total

    2. create only plots for author1 and author2

	./fsiv.py -n --plots --authors='author1 author2'

    3. run simulations and create plots for all particles for dsigma/denergy data sets with custom command

	./fsiv.py --energy --command='my_command --option'

    4. run only simulations for total cross section for pip and kp with 1000 events

	./fsiv.py -t --isospin='pip kp' --nof_events_total=1000

"""
    print '*'*110,'\n'

def init ():
    #set global switchers respect to the command line arguments (see help for details)
    global do_pion, do_nucl, do_kaon, do_tot, do_nrg, do_ang, do_sims, do_plot, nof_events_diff, nof_events_total, command, log_file, todo_authors, todo_isospins
    command = None				#command to run bash scripts (e.g. jobsub -q script_name.sh)

    nof_events_diff = '1000000'			#number of events for runfast.pl for differential cross section
    nof_events_total = '100000'			#number of events for runfast.pl for total cross section

    do_pion = do_nucl = do_kaon = False		#switchers for particles
    do_tot = do_nrg = do_ang = False		#switchers for types of plots
    do_sims = do_plot = False			#switchers for simulations and plots

    todo_authors  = None			#the list of authors to be done for diff xsec (if not defined -> all are done)
    todo_isospins = None			#the list of targets to be done (if not defined -> all are done)

    try:
	opts, args = getopt.getopt(sys.argv[1:], "hpnkteasf", ["help", "pion", "nucleon","kaon","total","energy","angle","simulation","plot","nof_events_diff=","nof_events_total=","command=","authors=","isospins=","FNAL"])
    except getopt.GetoptError as err:
	help()
	print 'INPUT ERROR: ', str(err), '\n'
	sys.exit(2)
    
    for o, a in opts:
	if o in ("-h", "--help"):
	    help()
	    sys.exit(0)
	elif o in ("-p", "--pion"): do_pion = True
	elif o in ("-n", "--nucleon"): do_nucl = True
	elif o in ("-k", "--kaon"): do_kaon = True
	elif o in ("-t", "--total"): do_tot = True
	elif o in ("-e", "--energy"): do_nrg = True
	elif o in ("-a", "--angle"): do_ang = True
	elif o in ("-s", "--simulation"): do_sims = True
	elif o in ("-f", "--plot"): do_plot = True
	elif o in ("--nof_events_diff"): nof_events_diff = a
	elif o in ("--nof_events_total"): nof_events_total = a
	elif o in ("--command"): command = a
	elif o in ("--authors"): todo_authors = a.split(' ')
	elif o in ("--isospins"): todo_isospins = a.split(' ')
	elif o in ("--FNAL"): command = "jobsub -e GENIE -e LD_LIBRARY_PATH -e PATH -e HOME"
	
    if not do_pion + do_nucl + do_kaon: do_pion = do_nucl = do_kaon = True	#if no particle is chosen all are set on
    if not do_tot + do_nrg + do_ang: do_tot = do_nrg = do_ang = True		#if no xsec type is chosen all are set on
    if not do_sims + do_plot: do_sims = do_plot = True				#if no task is selected do both: run simulations and prepare plots

    print_options()
    if not opts: print "\nYou did not choose any option. Default settings will be used. Run ./fsiv.py -h to see available options."
    while True:
	c = raw_input("\nDo you want to proceed ['Y' or 'N']? ")
	if c.lower() == 'n': sys.exit()
	if c.lower() == 'y': break;

    call(["mkdir", "-p", LOG_DIR])
    log_file = open(LOG_DIR + "/fsiv_" + date_ + ".log", "w")

def print_options (dest=None):
    print >> dest, '\nThe following options were chosen:\n'
    print >> dest, 'Pion mode:\t',do_pion
    print >> dest, 'Nucleon mode:\t',do_nucl
    print >> dest, 'Kaon mode:\t',do_kaon
    print >> dest, '\nXsec data:\t',do_tot
    print >> dest, 'Energy data:\t',do_nrg
    print >> dest, 'Angle data:\t', do_ang
    print >> dest, '\nSimulations:\t', do_sims
    print >> dest, 'Plots:\t\t', do_plot
    print >> dest, '\nNof events:\t', nof_events_diff, '(differential), ', nof_events_total, '(total)'
    print >> dest, 'Command:\t', command
    print >> dest, 'Author cut:\t', todo_authors
    print >> dest, 'Isospin cut:\t', todo_isospins
    print >> dest, '\nData directory:\t', DAT_DIR
    print >> dest, 'ROOT/xsec dir:\t', SIM_DIR
    print >> dest, 'PNG directory:\t', PNG_DIR
    print >> dest, 'Logs directory:\t', LOG_DIR
    print >> dest, 'Scripts dir:\t', BSH_DIR
    
def read_authors(dir_):
    #read authors list from given directory
    global wrongfiles
    wrongfiles=[]
    fates = ['total', 'tot', 'reac', 'elas', 'inel', 'abs', 'cex']
    pions = ['pip', 'pim']
    nucls = ['p', 'n', 'pp', 'np']
    kaons = ['kp', 'km']
    for root, dirs, files in os.walk(dir_):
	if '.svn' in dirs:
	    dirs.remove('.svn')							#skip .svn directory
	for file in files:
	    if file.endswith(".dat"):						#look only at files with 'dat' extension
		filename  = re.split('[- .]',file)				#extract process info assuimng the following convention: author-particle-target-fate.dat (total) or author-particle-*.dat (differential)
		if len(filename) < 5: 						#skip file with less then 4 parts in the filename
		    wrongfiles.append(file)
		    continue
		author   = filename[0]
		particle = filename[1]
		if author not in authors: authors.append(Author(author))	#add author to the list (if not exist)
		index = authors.index(author)					#set index for given author
		if particle in pions: authors[index].pion = True		#add particle mode corresponding to filename
		elif particle in nucls: authors[index].nucl = True
		elif particle in kaons: authors[index].kaon = True
		else:
		    wrongfiles.append(file)
		    continue
		if filename[3] in fates:					#is total cross section data
		    authors[index].tot = True					#add total xsec distribution for this author
		    target = filename[2]
		    fate = filename[3]
		    tk = map(float,list(line.split()[0] for line in open(root+'/'+file, 'r') \
			    if line.strip() and not line.strip().startswith('#')))		#take first column of any non-empty and non-commented line (kinetic energy)
		    xs = map(float,list(line.split()[1] for line in open(root+'/'+file, 'r') \
			    if line.strip() and not line.strip().startswith('#')))		#take second column of any non-empty and non-commented line (xsec)
		    authors[index].tot_proc.append( \
			    [particle, \
			    tk,
			    max(xs), \
			    target, \
			    fate])						#add proccess information in the following format: [particle, [tk list], xsec_max, target, fate]
		elif file.endswith("angdist.dat"): authors[index].ang = True	#add dsigma/dangle distribution for this author if the file ends with 'angdist.dat'
		else: authors[index].nrg = True;				#add disgma/denergy distribution if not 'total' nor 'angdist'

def auto_set_tk ():
    #remove tk points where to dense and add where to rare
    for proc in tot_proc_list:					#for each available process for total xsec
	res = []						#storage for result
	temp = set([int(x) for x in proc.tk_data])		#take data tk points (and round them to int)
	temp = filter(lambda x: x > 10, temp)			#remove all points less than 10
	if len(temp) < 3: res = temp				#do nothing if there are less than 3 points
	else:
	    temp = sorted(temp, reverse=True)			#sort points (so pop() takes the lowest value)
	    res.append(temp.pop())				#save first tk point
	    first_tk = sum_tk = temp.pop()			#take the lowest
	    counter = 1						#and set counter to 1
	    while True:						#start loop
		next_tk = temp.pop()				#take next one
		if not temp:					#if it is the last one
		    res.append(sum_tk/counter)			#save the current average tk
		    res.append(next_tk)				#save the last tk value
		    break					#and stop the loop
		if float(next_tk)/float(first_tk) < 1.25:	#if difference between first_tk and next_tk is smaller than 25%
		    sum_tk += next_tk				#add next_tk to the sum
		    counter += 1				#and increase counter
		else:						#if not
		    res.append(sum_tk/counter)			#save current average tk
		    counter = 1					#set counter to 1
		    first_tk = sum_tk = next_tk			#next_tk is new first_tk

	res.sort()
	todo = True
	while todo:						#start loop
	    for x in res:					#go through all tk points
		index = res.index(x)
		if not index == len(res)-1:			#if it is not last point
		    if  float(res[index+1]) / float(x) > 2.0:	#if the difference between two adjacent points is greater than 200%
			res.append(int((x+res[index+1])/2))	#add the average of them
			res.sort()
			break					#and start for loop from begining
		else:						#stop main loop if reach last point
		    todo = False
		    break
	    
	res.sort()
	proc.tk_sim = set(res)

def prepare_total():
    #create a list of particle-nuclei processes to be simulated
    for author in filter(lambda a: a.tot, authors):				#for all authors with total xsec data
	for tp in author.tot_proc:						#take each type of process
	    temp = TotalProcess(tp)						#add or update if already exists
	    if temp not in tot_proc_list: tot_proc_list.append(temp)
	    else: tot_proc_list[tot_proc_list.index(temp)].update(temp)
    auto_set_tk()

def log(text = None):
    #print 'text' into a log file
    if text: print >> log_file, text
    #if 'text' is not given, print the list of chosen options, found authors and list of processes to be simulated for total xsec
    else:
	print_options(log_file)
	print >> log_file, '\nThe following non-data files were found in ' + DAT_DIR + ': ', wrongfiles
	print >> log_file, '\nThe following authors were found in ' + DAT_DIR + ':\n'
	for author in authors: author.print_(log_file)
	print >> log_file, '-'*22, '\n'
	print >> log_file, 'The following processes/fates are available for total cross section:\n'
	for x in tot_proc_list: x.print_(log_file)
    log_file.flush()

def make_dirs(xsec,particle):
    #create folders for logs, bash scripts, root/xsec files and plots
    call(["mkdir", "-p", LOG_DIR+"/runfast/"+xsec+"/"+particle])
    call(["mkdir", "-p", LOG_DIR+"/intranukeplotter/"+xsec+"/"+particle])
    call(["mkdir", "-p", PNG_DIR+"/"+xsec+"/"+particle])
    call(["mkdir", "-p", BSH_DIR+"/"+xsec+"/"+particle])
    call(["mkdir", "-p", SIM_DIR+"/"+xsec+"/"+particle])

def create_dirs():
    #call make_dirs for turn on processes
    if do_tot:	
	if do_pion: make_dirs("total", "pion")
	if do_nucl: make_dirs("total", "nucleon")
	if do_kaon: make_dirs("total", "kaon")
    if do_nrg :
	if do_pion: make_dirs("dbldiff_dist", "pion")
	if do_nucl: make_dirs("dbldiff_dist", "nucleon")
	if do_kaon: make_dirs("dbldiff_dist", "kaon")
    if do_ang:
	if do_pion: make_dirs("ang_dist", "pion")
	if do_nucl: make_dirs("ang_dist", "nucleon")
	if do_kaon: make_dirs("ang_dist", "kaon")

def run_bs (bs):
    #run bash script
    print "\n\trunning " + bs + "\n"
    log("\trunning " + bs)
    if command: call(command + " 'bash " + bs + "'", shell=True)	#if custom command is defined add it in front of 'bash bash_script.sh'
    else: call("bash " + bs, shell=True)

def run_total (particle):
    log("\nTotal xsec validation for " + particle + ":\n")
    for proc in filter(lambda p: p.particle_kind == particle,tot_proc_list):									#go through all available processes for given particle
	if not todo_isospins or proc.particle in todo_isospins:											#check isospin filter
	    bash_scrt = BSH_DIR + "/total/" + particle + "/" + proc.particle + "_" + proc.target + "_" + date_ + ".sh"				#create bash script file
	    bash_file = open(bash_scrt, "w")
	    print >> bash_file, "#!/bin/bash"
	    if do_sims:																#if simulations are on
		runfast = "perl runfast.pl --type totxs --rm yes" \
		    + " --p " + proc.particle \
		    + " --t " + proc.target \
		    + " --el '" + str(sorted(proc.tk_sim))[1:-1] + "'" \
		    + " --n " + nof_events_total \
		    + " --rootdir " + SIM_DIR + "/total/" + particle \
		    + " 1>/dev/null 2>" + LOG_DIR + "/runfast/total/" + particle + "/" + proc.particle + "_" + proc.target + "_" + date_ + ".log"		#prepare runfast.pl command
		print >> bash_file, runfast													#put the command to bash sctipt
	    if do_plot:																#if plots are on
		print >> bash_file, "shopt -s nocaseglob"											#turn off case sensivity
		temp = "f=$(ls " + SIM_DIR + "/total/" + particle + "/*_" + proc.particle + "_" + proc.target + "_*.txt -Art | tail -n 1)"	#bash script will find the recent *_particle_target_*.txt file
		print >> bash_file, temp
		print >> bash_file, "m=$(echo ${f##*/} | cut -d'_' -f 1)"									#and extract its date for intranukeplotter
		print >> bash_file, "d=$(echo ${f##*/} | cut -d'_' -f 2)"
		print >> bash_file, "y=$(echo ${f##*/} | cut -d'_' -f 3)"
		for fate in proc.fate:														#for each available fate for given process
		    index = proc.fate.index(fate)
		    plotter = "perl intranukeplotter.pl --type totxs --rm yes" \
			+ " --stype " + fate \
			+ " --p " + proc.particle \
			+ " --t " + proc.target \
			+ " --vmax " + str(1.3*proc.xsec_max[index]) \
			+ " --hmax " + str(1.1*proc.tk_max[index]) \
			+ " --dorf $m'_'$d'_'$y" \
			+ " --rootdir " + SIM_DIR + "/total/" + particle \
			+ " --datadir " + DAT_DIR + "/total/" + proc.particle_kind \
			+ " --pngdir " + PNG_DIR + "/total/" + proc.particle_kind \
			+ " > " + LOG_DIR + "/intranukeplotter/total/" + particle + "/" + proc.particle + "_" + proc.target + "_" + date_ + ".log"	#prepare intranukeplotter.pl command
		    print >> bash_file, plotter			#put the command to bash script
	    bash_file.close()
	    run_bs(bash_scrt)					#run bash script
	    

def run_diff (particle,ang=False):
    if ang: log("\ndsigma/dtheta validation for " + particle + ":\n")
    else: log("\ndsigma/dE validation for " + particle + ":\n")
    dir_ = "dbldiff_dist"
    type_ = "nrg"
    if  ang:
	dir_ = "ang_dist"
	type_ = "ang"
    for author in filter(lambda a: a.check(particle),authors): 									#go through all available authors for given particle
	if (not ang and author.nrg) or (ang and author.ang):
	    bash_scrt = BSH_DIR + "/" + dir_ + "/" + particle + "/" + author.name + "_" + date_ + ".sh"				#create bash script file
	    bash_file = open(bash_scrt, "w")
	    print >> bash_file, "#!/bin/bash"
	    if do_sims:														#if simulations are on
		runfast = "perl runfast.pl --type root --rm yes --name yes" \
		    + " --n " + nof_events_diff \
		    + " --a " + author.name \
		    + " --rootdir " + SIM_DIR + "/" + dir_ + "/" + particle \
		    + " 1>/dev/null 2>" + LOG_DIR + "/runfast/" + dir_ + "/" + particle + "/" + author.name + "_" + date_ + ".log"		#prepare runfast.pl command
		print >> bash_file, runfast											#put the command to bash script
	    if do_plot:														#if plots are on
		print >> bash_file, "shopt -s nocaseglob"									#turn off case sensivity
		temp = "f=$(ls " + SIM_DIR + "/" + dir_ + "/" + particle + "/" + author.name + "_*.root -Art | tail -n 1)"	#bash script will find the recent author_*.root file
		print >> bash_file, temp
		print >> bash_file, "m=$(echo ${f##*/} | cut -d'_' -f 2)"							#and extract its date for intranukeplotter
		print >> bash_file, "d=$(echo ${f##*/} | cut -d'_' -f 3)"
		print >> bash_file, "y=$(echo ${f##*/} | cut -d'_' -f 4)"
		plotter = "perl intranukeplotter.pl --type " + type_ + " --name yes --rm yes" \
			+ " --rootdir " + SIM_DIR + "/" + dir_ + "/" + particle \
			+ " --datadir " + DAT_DIR + "/" + dir_ + "/" + particle \
			+ " --pngdir " + PNG_DIR + "/" + dir_ + "/" + particle \
			+ " --dorf $m'_'$d'_'$y"  \
			+ " --a " + author.name \
			+ " >" + LOG_DIR + "/intranukeplotter/" + dir_ + "/" + particle + "/" + author.name + "_" + date_ + ".log"	#prepare intrnukeplotter.pl command
		print >> bash_file, plotter		#put the command to bash script
	    bash_file.close()
	    run_bs(bash_scrt)				#run bash script

#MAIN PROGRAM
    
if __name__ == "__main__":

    init()
    read_authors(DAT_DIR)
    prepare_total()
    log()
    create_dirs()

    if do_tot:
	if do_pion: run_total ('pion')
	if do_nucl: run_total ('nucleon')
	if do_kaon: run_total ('kaon')
	
    if do_nrg:
	if do_pion: run_diff ('pion')
	if do_nucl: run_diff ('nucleon')
	if do_kaon: run_diff ('kaon')

    if do_ang:
	if do_pion: run_diff ('pion',True)
	if do_nucl: run_diff ('nucleon',True)
	if do_kaon: run_diff ('kaon',True)

    log('\nDONE')
