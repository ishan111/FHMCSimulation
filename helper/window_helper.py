"""@docstring
@brief Sample generation of windows and respective input files using FHMCAnalysis library
@author Nathan A. Mahynski
@date 12/23/2016
@filename window_helper.py
"""

import re, json
import math as m

def pure_settings (settings):
	"""
	Example of settings for the prototypical pure "lambda 1.5" square-well fluid

	Parameters
	----------
	settings : dict
		Dictionary containing (beta, bounds) where, beta = 1/T* for the simulation, bounds = (low, high) for Ntot

	Returns
	-------
	dict
		Information for input file to be read by FHMCSimulation binary

	"""

	bounds = settings["bounds"]
	beta = settings["beta"]

	info = {}

	info["num_species"] = 1
	info["beta"] = beta
	info["box"] = [8.0, 8.0, 8.0]
	info["mu"] = [0.0*i for i in xrange(info["num_species"])]
	info["seed"] = -10
	info["max_N"] = [int(1200)]
	info["min_N"] = [int(0)]
	info["window"] = [bounds[0], bounds[1]]
	info["restart_file"] = ""
	info["num_expanded_states"] = int(1)
	info["tmmc_sweep_size"] = int(100)
	info["total_tmmc_sweeps"] = int(1e3)
	info["wala_sweep_size"] = int(1e6)
	info["num_crossover_visits"] = int(100)
	info["lnF_start"] = 1.0
	info["lnF_end"] = 1.0e-8
	info["wala_g"] = 0.5
	info["wala_s"] = 0.8
	info["prob_pr_ins_del_1"] = 0.7
	info["prob_pr_displace_1"] = 0.3
	info["prob_eq_ins_del_1"] = 0.7
	info["prob_eq_displace_1"] = 0.3
	info["max_pr_displacement_1"] = 0.2
	info["max_eq_displacement_1"] = 0.2
	info["ppot_1_1"] = "square_well"
	info["ppot_1_1_params"] = [1.0, 0.5, 1.0] # [sig, width, eps]
	info["ppot_1_1_use_cell_list"] = True

	return info

def make_input (filename, settings, generator):
	"""
	Example production of input file for FHMCSimulation.

	Parameters
	----------
	filename : str
		Name of file to produce
	settings : dict
		Dictionary of settings, user defined
	generator : function
		Takes settings tuple as only argument and returns json input as dictionaryi, e.g., pure_settings()

	"""

	info = generator(settings)
	f = open(filename, 'w')
	json.dump(info, f, sort_keys=True, indent=4)
	f.close()

def make_sleeper (filename):
	"""
	Example of rsync "sleeper" that syncs remote directory to local one every 15 minutes.

	Parameters
	----------
	filename : str
		Name of file to produce

	"""

	sync_every = 15 # minutes

	f = open(filename, 'w')
	f.write('for ((i = 0; i < 100000; ++i)); do\n\tfor ((j = 0; j < 100000; ++j)); do\n\t\tsleep $((1*'+str(int(sync_every))+'*60));\n\t\trsync -a $1 $2;\n\tdone;\ndone;')
	f.close()

def gibbs_qsub (num_windows, binary, git_head, tag, prefix, input_name="input.json", jobs_per_node=12, q="medium", scratch_dir="/tmp/nam4/"):
	"""
	Example of submission script for PBS system, in this case, for gibbs.nist.gov.
	This produces qsub_X.pbs files, for as many X are necessary.

	Parameters
	----------
	num_windows : int
		Number of total windows
	binary : str
		Path to binary (FHMCSimulation) to execute
	git_head : str
		Path to git head log
	tag : str
		Name of this job
	prefix : str
		Directory to place qsub files
	input_name : str
		Name of local input file inside window directory (default=input.json)
	jobs_per_node : int
		Number of jobs per node to allow (default=12)
	q : str
		Name of the queue to submit to (default="medium")
	scratch_dir : str
		Absolute path to scratch space for user
	
	
	"""

	base_string = "#!/bin/bash\n#PBS -l nodes=1:ppn=__PPNVIS__\n#PBS -q __QUEUE__\n#PBS -V\n#PBS -N __TAGNAME__\n#PBS -M nathan.mahynski@nist.gov\n#PBS -m a\n\n# move to dir pbs was launched from\ncd $PBS_O_WORKDIR;\n\n# report\necho \"Running on $(hostname)\";\necho \"Time: $(date)\";\necho \"Starting directory: $PWD\";\nheaddir=$PWD;\npids=\"\";\nfor i in {__MINWIN__..__MAXWIN__}; do\n\t# create and move to tmpdir\n\ttmpdir=__SCRATCHDIR__/$PBS_JOBID/$i;\n\thomedir=$headdir/$i;\n\tmkdir -p $tmpdir;\n\techo \"Moving to temporary directory: $tmpdir\";\n\tcp -r $homedir/* $tmpdir/;\n\tcd $tmpdir;\n\n\t# get info about binary\n\ttail -1 __GITHEAD__ > binary.info;\n\n\t# run sleeper and binary\n\tsh sleeper.sh ./ $homedir &\n\t__BINARY__ __INPUTNAME__ 2>> err >> log &\n\tpids=\"$pids $!\"\ndone\n\n# wait for all process ids\nfor p in $pids; do\n\twait $p;\ndone\n\n# final sync and clean up\ncd $headdir;\nsids=\"\";\nfor i in {__MINWIN__..__MAXWIN__}; do\n\ttmpdir=__SCRATCHDIR__/$PBS_JOBID/$i;\n\thomedir=$headdir/$i;\n\techo \"Final sync from $tmpdir to $homedir\";\n\trsync -a $tmpdir/ $homedir/;\n\techo \"Removing $tmpdir\";\n\trm -r $tmpdir;\n\tsids=\"$sids $!\"\ndone\n\n# wait for all syncs to finish\nfor s in $sids; do\n\twait $s;\ndone\n\necho \"Finished on $(date)\";"

	new_string = re.sub('__QUEUE__', str(q), base_string)
	new_string = re.sub('__SCRATCHDIR__', str(scratch_dir), new_string)
	new_string = re.sub('__GITHEAD__', str(git_head), new_string)
	new_string = re.sub('__BINARY__', str(binary), new_string)
	new_string = re.sub('__INPUTNAME__', str(input_name), new_string)

	jobs_remaining = num_windows
	for idx in xrange(0, int(m.ceil(num_windows/float(jobs_per)))):
		start = idx*jobs_per+1
		end = idx*jobs_per+min([jobs_remaining,jobs_per])
		i_string = re.sub('__MINWIN__', str(start), new_string)
		i_string = re.sub('__MAXWIN__', str(end), i_string)
		i_string = re.sub('__TAGNAME__', str(tag+"_"+str(idx+1)), i_string)
		i_string = re.sub('__PPNVIS__', str(2*(end-start+1)), i_string) # double to deal with hyperthreading
		jobs_remaining -= jobs_per

		f = open(prefix+"/qsub_"+str(idx+1)+".pbs", 'w')
		f.write(i_string)
		f.close()

	return

def raritan_sbatch (num_windows, binary, git_head, tag, prefix, input_name="input.json", jobs_per_node=12, q="mml", hours=72, scratch_dir="/scratch/nam4/"):
	"""
	Example of submission script for PBS system, in this case, for raritan.nist.gov.
	This produces sbatch_X.sb files, for as many X are necessary.

	Parameters
	----------
	num_windows : int
		Number of total windows
	binary : str
		Path to binary (FHMCSimulation) to execute
	git_head : str
		Path to git head log
	tag : str
		Name of this job
	prefix : str
		Directory to place sbatch files	
	input_name : str
		Name of local input file inside window directory (default=input.json)
	jobs_per_node : int
		Number of jobs per node to allow (default=12)
	q : str
		Name of the queue to submit to (default="mml")
	hours : int
		Number of hours to submit the job for (default=72)
	scratch_dir : str
		Absolute path to scratch space for user (default="/scratch/nam4/")

	"""

	sdays = int(m.floor(float(hours)/24.))
	shours = int(m.floor(float(hours-24*sdays)))
	sminutes = int(m.floor(float(hours-24*sdays-shours)*60))

	if (sdays > 14 or (sdays == 14 and (shours > 0 or sminutes > 0))):
		raise Exception ("Cannot exceed 14 days on raritan queues")

	base_string = "#!/bin/bash\n#\n#SBATCH -n __PPNVIS__\n#SBATCH -N 1\n#SBATCH -p __QUEUE__\n#SBATCH -t __SDAYS__-__SHOURS__:__SMINUTES__\n#SBATCH -m __TOTMEM__\n#SBATCH -o hostname_%j.out\n#SBATCH -e hostname_%j.err\n#SBATCH --mail-type=FAIL\n#SBATCH --mail-user=nathan.mahynski@nist.gov\n#SBATCH -J __TAGNAME__\n#SBATCH --export=all\n\n# move to dir pbs was launched from\ncd $SLURM_SUBMIT_DIR;\n\n# report\necho \"Running on $(hostname)\";\necho \"Time: $(date)\";\necho \"Starting directory: $PWD\";\nheaddir=$PWD;\npids=\"\";\nfor i in {__MINWIN__..__MAXWIN__}; do\n\t# create and move to tmpdir\n\ttmpdir=__SCRATCHDIR__/$SLURM_JOBID/$i;\n\thomedir=$headdir/$i;\n\tmkdir -p $tmpdir;\n\techo \"Moving to temporary directory: $tmpdir\";\n\tcp -r $homedir/* $tmpdir/;\n\tcd $tmpdir;\n\n\t# get info about binary\n\ttail -1 __GITHEAD__ > binary.info;\n\n\t# run sleeper and binary\n\tsh sleeper.sh ./ $homedir &\n\t__BINARY__ __INPUTNAME__ 2>> err >> log &\n\tpids=\"$pids $!\"\ndone\n\n# wait for all process ids\nfor p in $pids; do\n\twait $p;\ndone\n\n# final sync and clean up\ncd $headdir;\nsids=\"\";\nfor i in {__MINWIN__..__MAXWIN__}; do\n\ttmpdir=__SCRATCHDIR__/$SLURM_JOBID/$i;\n\thomedir=$headdir/$i;\n\techo \"Final sync from $tmpdir to $homedir\";\n\trsync -a $tmpdir/ $homedir/;\n\techo \"Removing $tmpdir\";\n\trm -r $tmpdir;\n\tsids=\"$sids $!\"\ndone\n\n# wait for all syncs to finish\nfor s in $sids; do\n\twait $s;\ndone\n\necho \"Finished on $(date)\";"

	new_string = re.sub('__QUEUE__', str(q), base_string)
	new_string = re.sub('__SCRATCHDIR__', str(scratch_dir), new_string)
	new_string = re.sub('__GITHEAD__', str(git_head), new_string)
	new_string = re.sub('__BINARY__', str(binary), new_string)
	new_string = re.sub('__INPUTNAME__', str(input_name), new_string)
	new_string = re.sub('__SDAYS__', str(sdays), new_string)
	new_string = re.sub('__SHOURS__', str(shours), new_string)
	new_string = re.sub('__SMINUTES__', str(sminutes), new_string)

	jobs_remaining = num_windows
	for idx in xrange(0, int(m.ceil(num_windows/float(jobs_per)))):
		start = idx*jobs_per+1
		end = idx*jobs_per+min([jobs_remaining,jobs_per])
		totmem = 100*(end-start+1) # MB/job, assuming 100MB per job
		totmem = max(totmem, 4000) # Impose a 4GB min reservation on raritan
		i_string = re.sub('__MINWIN__', str(start), new_string)
		i_string = re.sub('__MAXWIN__', str(end), i_string)
		i_string = re.sub('__TAGNAME__', str(tag+"_"+str(idx+1)), i_string)
		i_string = re.sub('__PPNVIS__', str(2*(end-start+1)), i_string) # double to deal with hyperthreading
		i_string = re.sub('__TOTMEM__', str(totmem), i_string)
		jobs_remaining -= jobs_per

		f = open(prefix+"/sbatch_"+str(idx+1)+".sb", 'w')
		f.write(i_string)
		f.close()

	return

if __name__ == "__main__":
	print "window_helper.py"

	"""
	Tutorial:
	
	* Below is an example of a script to use these functions to produce windows

	import os, sys, shutil
	FHMCLIB = "/home/nam4/Desktop/sandbox/"
	sys.path.append(FHMCLIB)
	import FHMCAnalysis
	import FHMCAnalysis.moments.win_patch.windows as win
	
	HELPLIB = "./"
	sys.path.append(HELPLIB)
	import window_helper as hP

	# Overwrite existing inputs
	overwrite = True

	# Simulation settings
	beta = 1.0

	# Establish bounds for windows
	ntot_max = 500
	final_window_width = 20
	num_windows = 30
	num_overlap = 6

	# Need installation information
	install_dir = "/home/nam4/FHMCSimulation/"
	binary = install_dir+"/bin/fhmc_tmmc"
	git_head = install_dir+"/.git/logs/HEAD"
	jobs_per = 12
	scratch_dir = "/tmp/nam4/"
	q = "medium"
	tag = "example"

	# Window settings
	input_name = "input.json"
	prefix = "./"
	bounds = win.ntot_window_scaling (ntot_max, final_window_width, num_windows, num_overlap)
	if (!os.path.isdir(prefix)):
		os.makedirs(prefix)

	sett = {}
	sett["beta"] = beta
	for w in range(num_windows):
		dname = prefix+"/"+str(w+1)
		if ((str(w+1) in os.listdir(prefix)) and overwrite):
			shutil.rmtree(dname, pure_settings)
		os.makedirs(dname)

		sett["bounds"] = bounds[w]
		hP.make_input (dname+"/"+input_name, sett, hP.pure_settings)
		hP.make_sleeper (dname+"/sleeper.sh")

	hP.gibbs_qsub (num_windows, binary, git_head, tag, './', input_name, jobs_per, q, scratch_dir)
	"""

