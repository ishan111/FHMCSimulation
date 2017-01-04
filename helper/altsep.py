"""@docstring
@brief This is an example script for creating binary systems for ALTSEP
@author Nathan A. Mahynski
@date 12/27/2016
@filename altsep.py
"""

import numpy as np

def random_word (length):
   return ''.join(random.choice(string.lowercase) for i in range(length))

def binary_fslj_pore (settings):
	"""
	Create settings for binary FS-LJ system inside a cylindrical pore.
	Basis of sigma_11 = epsilon_11 = 1.0 (Fluid-Fluid).

	Parameters
	----------
	settings : dict
		Settings continaining {"T_eps11", "sig22_sig11", "eps22_eps11", "lam1w", "lam2w", "eps1w_eps11", "eps2w_eps11", "bounds", "r_pore", "mu"}

	"""

	# Basis
	sig11 = 1.0
	eps11 = 1.0

	# Mixing rules
	eta_e = 1.0
	eta_s = 1.0

	# Pore length
	Lz = 22.0

	# Estimated packing efficiency at max filling
	eta_p = 0.63

	info = {}

	# Simulation information
	info["num_species"] = 2
	info["beta"] = 1.0/settings["T_eps11"]*eps11
	info["mu"] = [settings["mu"][0], settings["mu"][1]]
	info["seed"] = -10
	info["min_N"] = [int(0), int(0)]
	info["window"] = [settings["bounds"][0], settings["bounds"][1]]
	info["restart_file"] = ""
	info["num_expanded_states"] = int(1)
	info["tmmc_sweep_size"] = int(200)
	info["total_tmmc_sweeps"] = int(1e3)
	info["wala_sweep_size"] = int(1e6)
	info["num_crossover_visits"] = int(1000)
	info["lnF_start"] = 1.0
	info["lnF_end"] = 1.0e-8
	info["wala_g"] = 0.5
	info["wala_s"] = 0.8

	# Barrier(s)
	info["barriers"] = {}
	info["barriers"]["cylindrical_pore1"] = {}
	info["barriers"]["cylindrical_pore1"]["type"] = "cylinder_z"
	info["barriers"]["cylindrical_pore1"]["species"] = 1
	info["barriers"]["cylindrical_pore1"]["radius"] = r_pore
	info["barriers"]["cylindrical_pore1"]["sigma"] = sig11
	info["barriers"]["cylindrical_pore1"]["epsilon"] = settings["eps1w_eps11"]*eps11
	info["barriers"]["cylindrical_pore1"]["width"] = settings["lam1w"]*info["barriers"]["cylindrical_pore1"]["sigma"]

	info["barriers"]["cylindrical_pore2"] = {}
	info["barriers"]["cylindrical_pore2"]["type"] = "cylinder_z"
	info["barriers"]["cylindrical_pore2"]["species"] = 2
	info["barriers"]["cylindrical_pore2"]["radius"] = r_pore
	info["barriers"]["cylindrical_pore2"]["sigma"] = settings["sig22_sig11"]*sig11
	info["barriers"]["cylindrical_pore2"]["epsilon"] = settings["eps2w_eps11"]*eps11
	info["barriers"]["cylindrical_pore2"]["width"] = settings["lam2w"]*info["barriers"]["cylindrical_pore2"]["sigma"]

	# Determine global max particle bounds based on max packing efficiency stipulated
	maxN1 = int(np.ceil(eta_p*np.pi*((r_pore - info["barriers"]["cylindrical_pore1"]["sigma"]/2.0)**2)*Lz))
	maxN2 = int(np.ceil(eta_p*np.pi*((r_pore - info["barriers"]["cylindrical_pore2"]["sigma"]/2.0)**2)*Lz))
	info["max_N"] = [maxN1, maxN2]

	# Monte Carlo moves
	info["moves"] = {}
	info["moves"]["ins_del_1"] = 0.6
	info["moves"]["translate_1"] = 0.4
	info["moves"]["max_translation_1"] = 0.2
	info["moves"]["ins_del_2"] = 0.6
	info["moves"]["translate_2"] = 0.4
	info["moves"]["max_translation_2"] = 0.2
	info["moves"]["swap_1_2"] = 0.2

	# Pair potentials
	info["ppot_1_1"] = "fs_lennard_jones"
	info["ppot_1_1_params"] = {}
	info["ppot_1_1_params"]["sigma"] = sig11
	info["ppot_1_1_params"]["r_cut"] = 3.0*info["ppot_1_1_params"]["sigma"]
	info["ppot_1_1_params"]["epsilon"] = eps11
	info["ppot_1_1_params"]["cell_list"] = True

	info["ppot_2_2"] = "fs_lennard_jones"
	info["ppot_2_2_params"] = {}
	info["ppot_2_2_params"]["sigma"] = settings["sig22_sig11"]*sig11
	info["ppot_2_2_params"]["r_cut"] = 3.0*info["ppot_2_2_params"]["sigma"]
	info["ppot_2_2_params"]["epsilon"] = settings["eps22_eps11"]*eps11
	info["ppot_2_2_params"]["cell_list"] = True

	info["ppot_1_2"] = "fs_lennard_jones"
	info["ppot_1_2_params"] = {}
	info["ppot_1_2_params"]["sigma"] = eta_s*(info["ppot_1_1_params"]["sigma"] + info["ppot_2_2_params"]["sigma"])/2.0
	info["ppot_1_2_params"]["r_cut"] = 3.0*info["ppot_1_2_params"]["sigma"]
	info["ppot_1_2_params"]["epsilon"] = eta_e*np.sqrt(info["ppot_1_1_params"]["epsilon"]*info["ppot_2_2_params"]["epsilon"])
	info["ppot_1_2_params"]["cell_list"] = True

	# Compute box size that prevents any periodic interactions (technically cylinder doesn't create that effect, but do anyway)
	ff_range = np.max([info["ppot_1_1_params"]["r_cut"], info["ppot_1_2_params"]["r_cut"], info["ppot_2_2_params"]["r_cut"]])
	fw_range = np.max([info["barriers"]["cylindrical_pore1"]["width"], info["barriers"]["cylindrical_pore2"]["width"]])
	Lxy = (2*r_pore + np.max([ff_range, fw_range]))*1.05 # 5% fudge factor

	# Ensure that cell list will always be created (need 3 cells in each direction)
	Lxy = np.max([Lxy, 3.0*ff_range*1.05]) # 5% fudge factor

	info["box"] = [Lxy, Lxy, Lz]
	info["barriers"]["cylindrical_pore1"]["x"] = Lxy/2.0
	info["barriers"]["cylindrical_pore1"]["y"] = Lxy/2.0
	info["barriers"]["cylindrical_pore2"]["x"] = Lxy/2.0
	info["barriers"]["cylindrical_pore2"]["y"] = Lxy/2.0

	return info

if __name__ == "__main__":
	print "altsep.py"

	"""

	* Tutorial: Below is an example of a script to generate input for FHMCSimulation for the ALTSEP project

	import os, sys, shutil
	import random, string
	FHMCLIB = "/home/nam4/"
	sys.path.append(FHMCLIB)
	import FHMCAnalysis
	import FHMCAnalysis.moments.win_patch.windows as win
	import FHMCSimulation.helper.window_helper as hP
	import FHMCSimulation.helper.altsep as altsep

	# System settings - user defined
	settings = {
	"T_eps11" : 1.0,
	"sig22_sig11" : 1.0,
	"eps22_eps11" : 1.0,
	"lam1w" : 1.5,
	"lam2w" : 1.5,
	"eps1w_eps11" : 5.0,
	"eps2w_eps11" : 5.0,
	"r_pore" : 4.5
	}

	# Overwrite existing inputs
	overwrite = True

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
	tag = "altsep-"+altsep.random_word(6)

	# Window settings
	input_name = "input.json"

	for dMu2 in [-2.94, -1.1, 0.0, 1.1, 2.94]:
		prefix = "./dMu2_"+str(dMu2)
		bounds = win.ntot_window_scaling (ntot_max, final_window_width, num_windows, num_overlap)

		if (not os.path.isdir(prefix)):
			os.makedirs(prefix)

		settings["mu"] = [0.0, dMu2]
		for w in range(num_windows):
			dname = prefix+"/"+str(w+1)
			if ((str(w+1) in os.listdir(prefix)) and overwrite):
				shutil.rmtree(dname)
			os.makedirs(dname)

			settings["bounds"] = bounds[w]
			hP.make_input (dname+"/"+input_name, settings, altsep.binary_fslj_pore)
			hP.make_sleeper (dname+"/sleeper.sh")

		hP.raritan_sbatch (num_windows, binary, git_head, tag, prefix, input_name, jobs_per, q, hours, scratch_dir)
	"""
