#!/usr/bin/python
import subprocess
import threading
import os
import sys

mode_prefix = sys.argv[1]
nucZs = [6]
nucAs = [12]
if mode_prefix == "EFF":
	nucZs = [1, 2, 6,  10, 13, 18, 26, 82 ]
	nucAs = [2, 4, 12, 20, 27, 40, 56, 208]


energy = "3.0"

num_ev = "50000"
neut_pdgs = ["14", "-14"]
if mode_prefix == "EFF":
	neut_pdgs = ["14"]
spline_dir = "../../../../splines/"
spline_suffix = "QE_splines.xml"
start_run = 500000
end_run   = 500010
data_dir = "./data/" + mode_prefix + "/"

spline = spline_dir + mode_prefix + "_" + spline_suffix

cmds = []


subprocess.call(["mkdir", data_dir]) 

for neut_pdg in neut_pdgs:
	data_dir = "./data/" + mode_prefix + "/" + "nu_" + neut_pdg + "/"
	subprocess.call(["mkdir", data_dir]) 
	for runnum in range(start_run, end_run + 1):
		for i in range(len(nucZs)):
			nucZ = nucZs[i]
			nucA = nucAs[i]
			folder = data_dir + str(nucZ) + "_" + str(nucA)
			cmd = ["mkdir", folder]
			subprocess.call(cmd)
			target = str(1000000000+nucA*10+nucZ*10000)
			
			cmd1 = ["gevgen", "-r", str(runnum), "-n", num_ev, "-p", neut_pdg, "-t", target,"-e", energy,
		          "--cross-sections", spline, "--seed", str(runnum), "--event-generator-list", "CCQE"]
			cmd2 = ["gntpc", "-i", "gntp." + str(runnum) + ".ghep.root", "-f", "gst"]
			cmds.append([cmd1, cmd2, folder])

num_threads = 32

mutex = threading.Lock()

def generate_genie_files():
	while True:
		mutex.acquire()
		if len(cmds) == 0:
			mutex.release()
			return
		cmd = cmds.pop()
		mutex.release()
		subprocess.call(cmd[0], cwd = cmd[2])
		subprocess.call(cmd[1], cwd = cmd[2])

threads = []
for i in range(num_threads):
	t = threading.Thread(target=generate_genie_files)
	t.daemon=True
	t.start()
	threads.append(t)
for t in threads:
	t.join()
