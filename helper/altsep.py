"""@docstring
@brief This is an example script for creating binary systems for ALTSEP
@author Nathan A. Mahynski
@date 12/27/2016
@filename altsep.py
"""

import numpy as np
import random, string

def random_word (length):
   return ''.join(random.choice(string.lowercase) for i in range(length))

def binary_fslj (settings):
	"""
	Create settings for binary FS-LJ system in the bulk.
	Basis of sigma_11 = epsilon_11 = 1.0 (Fluid-Fluid).

	Parameters
	----------
	settings : dict
		Settings continaining {"T_eps11", "sig22_sig11", "eps22_eps11", "bounds", "mu"}

	"""

	# Basis
	sig11 = 1.0
	eps11 = 1.0

	# Mixing rules
	eta_e = 1.0
	eta_s = 1.0

	# Estimated packing efficiency at max filling
	eta_p = 0.50 # 0.53 = volume of sphere / volume of cube

	info = {}

	# Simulation information
	info["box"] = [10.0, 10.0, 10.0]
	info["num_species"] = 2
	info["beta"] = 1.0/(settings["T_eps11"]*eps11)
	info["mu"] = [settings["mu"][0], settings["mu"][1]]
	info["seed"] = -10
	info["min_N"] = [int(0), int(0)]
	info["window"] = [settings["bounds"][0], settings["bounds"][1]]
	info["restart_file"] = ""
	info["num_expanded_states"] = int(1)
	info["tmmc_sweep_size"] = int(100)
	info["total_tmmc_sweeps"] = int(1e8)  # essentially infinite (will patch based on checkpoints)
	info["wala_sweep_size"] = int(1e6)
	info["num_crossover_visits"] = int(100) # int(500) # int(1e3)
	info["lnF_start"] = 1.0
	info["lnF_end"] = 4.0e-6 # 1.0e-8
	info["wala_g"] = 0.5
	info["wala_s"] = 0.8
	info["delta_u_hist"] = 5.0
	info["max_order"] = 3
	info["use_ke"] = False

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
	if (np.min(info["box"])/3.0 < info["ppot_1_1_params"]["r_cut"]):
		info["ppot_1_1_params"]["cell_list"] = False
	else:
		info["ppot_1_1_params"]["cell_list"] = True

	info["ppot_2_2"] = "fs_lennard_jones"
	info["ppot_2_2_params"] = {}
	info["ppot_2_2_params"]["sigma"] = settings["sig22_sig11"]*sig11
	info["ppot_2_2_params"]["r_cut"] = 3.0*info["ppot_2_2_params"]["sigma"]
	info["ppot_2_2_params"]["epsilon"] = settings["eps22_eps11"]*eps11
	if (np.min(info["box"])/3.0 < info["ppot_2_2_params"]["r_cut"]):
		info["ppot_2_2_params"]["cell_list"] = False
	else:
		info["ppot_2_2_params"]["cell_list"] = True

	info["ppot_1_2"] = "fs_lennard_jones"
	info["ppot_1_2_params"] = {}
	info["ppot_1_2_params"]["sigma"] = eta_s*(info["ppot_1_1_params"]["sigma"] + info["ppot_2_2_params"]["sigma"])/2.0
	info["ppot_1_2_params"]["r_cut"] = 3.0*info["ppot_1_2_params"]["sigma"]
	info["ppot_1_2_params"]["epsilon"] = eta_e*np.sqrt(info["ppot_1_1_params"]["epsilon"]*info["ppot_2_2_params"]["epsilon"])
	if (np.min(info["box"])/3.0 < info["ppot_1_2_params"]["r_cut"]):
		info["ppot_1_2_params"]["cell_list"] = False
	else:
		info["ppot_1_2_params"]["cell_list"] = True

	# Determine global max particle bounds based on max packing efficiency stipulated
	maxN1 = int(np.ceil(eta_p*info["box"][0]*info["box"][1]*info["box"][2]/(4./3.*np.pi*(info["ppot_1_1_params"]["sigma"]/2.0)**3)))
	maxN2 = int(np.ceil(eta_p*info["box"][0]*info["box"][1]*info["box"][2]/(4./3.*np.pi*(info["ppot_2_2_params"]["sigma"]/2.0)**3)))
	info["max_N"] = [maxN1+maxN2, maxN1+maxN2] # Produce an overestimate of upper bound so that bounds on Ntot do not exceed either of these numbers
	info["__maxN1__"] = maxN1
	info["__maxN2__"] = maxN2

	ff_range = np.max([info["ppot_1_1_params"]["r_cut"], info["ppot_1_2_params"]["r_cut"], info["ppot_2_2_params"]["r_cut"]])
	if (np.min(info["box"])/2.0 <= ff_range): raise Exception ("Must increase box size, range > L/2")

	return info

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
	eta_p = 0.50 # 0.53 = volume of sphere / volume of cube

	info = {}

	# Simulation information
	info["num_species"] = 2
	info["beta"] = 1.0/(settings["T_eps11"]*eps11)
	info["mu"] = [settings["mu"][0], settings["mu"][1]]
	info["seed"] = -10
	info["min_N"] = [int(0), int(0)]
	info["window"] = [settings["bounds"][0], settings["bounds"][1]]
	info["restart_file"] = ""
	info["num_expanded_states"] = int(1)
	info["tmmc_sweep_size"] = int(100)
	info["total_tmmc_sweeps"] = int(1e8) # essentially infinite (will patch based on checkpoints)
	info["wala_sweep_size"] = int(1e4) # int(1e6)
	info["num_crossover_visits"] = int(5*info["tmmc_sweep_size"]) # int(100)
	info["lnF_start"] = 1.0
	info["lnF_end"] = 4.0e-6 # 1.0e-8
	info["wala_g"] = 0.5
	info["wala_s"] = 0.8
	info["delta_u_hist"] = 5.0
	info["max_order"] = 3
	info["use_ke"] = False

	# Barrier(s)
	info["barriers"] = {}
	info["barriers"]["cylindrical_pore1"] = {}
	info["barriers"]["cylindrical_pore1"]["type"] = "cylinder_z"
	info["barriers"]["cylindrical_pore1"]["species"] = 1
	info["barriers"]["cylindrical_pore1"]["radius"] = settings["r_pore"]
	info["barriers"]["cylindrical_pore1"]["sigma"] = sig11
	info["barriers"]["cylindrical_pore1"]["epsilon"] = settings["eps1w_eps11"]*eps11
	info["barriers"]["cylindrical_pore1"]["width"] = settings["lam1w_sig11"]*info["barriers"]["cylindrical_pore1"]["sigma"]

	info["barriers"]["cylindrical_pore2"] = {}
	info["barriers"]["cylindrical_pore2"]["type"] = "cylinder_z"
	info["barriers"]["cylindrical_pore2"]["species"] = 2
	info["barriers"]["cylindrical_pore2"]["radius"] = settings["r_pore"]
	info["barriers"]["cylindrical_pore2"]["sigma"] = settings["sig22_sig11"]*sig11
	info["barriers"]["cylindrical_pore2"]["epsilon"] = settings["eps2w_eps11"]*eps11
	info["barriers"]["cylindrical_pore2"]["width"] = settings["lam2w_sig11"]*info["barriers"]["cylindrical_pore2"]["sigma"]

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

    # Determine global max particle bounds based on max packing efficiency stipulated
	maxN1 = int(np.ceil(eta_p*((settings["r_pore"] - info["barriers"]["cylindrical_pore1"]["sigma"]/2.0)**2)*Lz*6.0/info["ppot_1_1_params"]["sigma"]**3))
	maxN2 = int(np.ceil(eta_p*((settings["r_pore"] - info["barriers"]["cylindrical_pore2"]["sigma"]/2.0)**2)*Lz*6.0/info["ppot_2_2_params"]["sigma"]**3))
	info["max_N"] = [maxN1+maxN2, maxN1+maxN2] # Produce an overestimate of upper bound so that bounds on Ntot do not exceed either of these numbers
	info["__maxN1__"] = maxN1
	info["__maxN2__"] = maxN2

	# Compute box size that prevents any periodic interactions (technically cylinder doesn't create that effect, but do anyway)
	ff_range = np.max([info["ppot_1_1_params"]["r_cut"], info["ppot_1_2_params"]["r_cut"], info["ppot_2_2_params"]["r_cut"]])
	fw_range = np.max([info["barriers"]["cylindrical_pore1"]["width"], info["barriers"]["cylindrical_pore2"]["width"]])
	Lxy = (2*settings["r_pore"] + np.max([ff_range, fw_range]))*1.05 # 5% fudge factor

	# Ensure that cell list will always be created (need 3 cells in each direction)
	Lxy = np.max([Lxy, 3.0*ff_range*1.05]) # 5% fudge factor

	info["box"] = [Lxy, Lxy, Lz]
	info["barriers"]["cylindrical_pore1"]["x"] = Lxy/2.0
	info["barriers"]["cylindrical_pore1"]["y"] = Lxy/2.0
	info["barriers"]["cylindrical_pore2"]["x"] = Lxy/2.0
	info["barriers"]["cylindrical_pore2"]["y"] = Lxy/2.0

	if (np.min(info["box"])/2.0 <= ff_range): raise Exception ("Must increase box size, range > L/2")

	return info

if __name__ == "__main__":
	print "altsep.py"

	"""

	* Tutorial: Below is an example of a script to generate input for FHMCSimulation for the ALTSEP project

	import os, sys, shutil, math
	FHMCLIB = "/home/nam4/"
	sys.path.append(FHMCLIB)
	import FHMCAnalysis
	import FHMCAnalysis.moments.win_patch.windows as win
	import FHMCSimulation.helper.window_helper as hP
	import FHMCSimulation.helper.altsep as altsep

	# System settings - user defined
	settings = {
	"T_eps11" : 1.1,
	"sig22_sig11" : 1.5,
	"eps22_eps11" : 1.0,
	"lam1w" : 1.5,
	"lam2w" : 1.5,
	"eps1w_eps11" : 5.0,
	"eps2w_eps11" : 5.0,
	"r_pore" : 4.5
	}

	# Overwrite existing inputs
	overwrite = True

	# Need installation information
	install_dir = "/home/nam4/FHMCSimulation/"
	binary = install_dir+"/bin/fhmc_tmmc"
	git_head = install_dir+"/.git/logs/HEAD"
	jobs_per = 8
	scratch_dir = "/scratch/nam4/"
	q = "mml"
	hours = 72

	# Window settings
	input_name = "input.json"

	remove = ["__maxN1__", "__maxN2__"]

	for betaDMu2 in [-2.94, -1.1, 0.0, 1.1, 2.94]:
        dMu2 = betaDMu2*settings["T_eps11"] # eps11 = 1.0
		prefix = "./dMu2=%2.4f"%dMu2
		settings["mu"] = [0.0, dMu2]
		settings["bounds"] = [0, 0] # Dummy setting for now

		info = altsep.binary_fslj_pore(settings)
        ntot_max = max(info["__maxN1__"],info["__maxN2__"])

        # Establish bounds for windows
        final_window_width = 0.03*ntot_max # 3% heuristic
        num_windows = 24
        num_overlap = 6 # Usually a robust number that always works

		bounds = win.ntot_window_scaling (ntot_max, final_window_width, num_windows, num_overlap)

		if (not os.path.isdir(prefix)):
			os.makedirs(prefix)

		for w in range(num_windows):
			dname = prefix+"/"+str(w+1)
			if ((str(w+1) in os.listdir(prefix)) and overwrite):
				shutil.rmtree(dname)
			os.makedirs(dname)

			settings["bounds"] = bounds[w]
			hP.make_input (dname+"/"+input_name, settings, altsep.binary_fslj_pore, remove)
			hP.make_sleeper (dname+"/sleeper.sh")

		tag = "altsep-"+altsep.random_word(6)
		hP.raritan_sbatch (num_windows, binary, git_head, tag, prefix, input_name, jobs_per, q, hours, scratch_dir)
	"""
