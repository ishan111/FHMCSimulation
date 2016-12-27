"""@docstring
@brief This is an example script for creating binary systems for ALTSEP
@author Nathan A. Mahynski
@date 12/27/2016
@filename altsep.py
"""

def altsep (settings):
	"""
	Create settings for binary FS-LJ system inside a cylindrical pore.

	Parameters
	----------
	settings : tuple
		Settings

	"""

	
	
	return

if __name__ == "__main__":

	import os, sys, shutil
	FHMCLIB = "/home/nam4/"
	sys.path.append(FHMCLIB)
	import FHMCAnalysis
	import FHMCAnalysis.moments.win_patch.windows as win
	
	HELPLIB = "/home/nam4/FHMCSimulation/helper/"
	sys.path.append(HELPLIB)
	import window_helper as hP

	# Overwrite existing inputs
	overwrite = True

	# Establish bounds for windows
	ntot_max = 600
	final_window_width = 20
	num_windows = 20
	num_overlap = 6

	# Need installation information
	install_dir = "/home/nam4/FHMCSimulation/"
	binary = install_dir+"/bin/fhmc_tmmc"
	git_head = install_dir+"/.git/logs/HEAD"
	jobs_per = 12
	scratch_dir = "/scratch/nam4/"
	q = "medium"
	tag = "fs_binary_example"
	hours = 72

	# Window settings
	input_name = "input.json"
	prefix = "./"
	bounds = win.ntot_window_scaling (ntot_max, final_window_width, num_windows, num_overlap)

	# System settings
	beta = 1.0

	for w in range(num_windows):
		dname = prefix+"/"+str(w+1)
		if ((str(w+1) in os.listdir(prefix)) and overwrite):
			shutil.rmtree(dname)
		os.makedirs(dname)

		hP.make_input (dname+"/"+input_name, (bounds[w], beta), altsep)
		hP.make_sleeper (dname+"/sleeper.sh")

	hP.raritan_sbatch (num_windows, binary, git_head, tag, prefix, input_name, jobs_per, q, hours, scratch_dir)
