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
	settings : dict
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

	# System settings
	eps_11 = 1.0
	sig_11 = 1.0
	
        sett = {}
	sett["beta"] = 1.0*eps_11
        sett["sig"] = [sig_11, 1.2*sig_11] 
        sett["eps"] = [eps_11, 1.0*eps_11]
	sett["D_cyl"] = 9.0*sig_11
	sett["lam_w"] = []
        sett["eps_w"] = []

	ff_range = max(np.array(sett["lam"])*np.array(sett["sig"]))
	fw_range = max(np.array(sett["lam_w"])*np.array(sett["sig"]))
	Lxy = D_cyl + 2*max(ff_range, fw_range) 
	ntot_max = 600
	Lz = 
	sett["box"] = [Lxy, Lxy, Lz]

	# Establish bounds for windows
	final_window_width = 20
	num_windows = 24
	num_overlap = 6

	# Need installation information
	install_dir = "/home/nam4/FHMCSimulation/"
	binary = install_dir+"/bin/fhmc_tmmc"
	git_head = install_dir+"/.git/logs/HEAD"
	jobs_per = 12
	scratch_dir = "/scratch/nam4/"
	q = "mml"
	hours = 72
	tag = "fs_binary_example"

	# Window settings
	input_name = "input.json"
	
	for dMu2 in [-2.9, -1.4, 0.0, 1.4, 2.9]:	
		prefix = "./dMu2_"+str(dMu2)
		bounds = win.ntot_window_scaling (ntot_max, final_window_width, num_windows, num_overlap)

		if (!os.path.isdir(prefix)):
			os.makedirs(prefix)

		sett["dMu2"] = dMu2
		for w in range(num_windows):
			dname = prefix+"/"+str(w+1)
			if ((str(w+1) in os.listdir(prefix)) and overwrite):
				shutil.rmtree(dname)
			os.makedirs(dname)
	
			sett["bounds"] = bounds[w]
			hP.make_input (dname+"/"+input_name, sett, altsep)
			hP.make_sleeper (dname+"/sleeper.sh")

		hP.raritan_sbatch (num_windows, binary, git_head, tag, prefix, input_name, jobs_per, q, hours, scratch_dir)
