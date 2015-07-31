import numpy as np
import os

if __name__ == "__main__":
	T = np.linspace(1.27, 1.40, 4)
	MU = np.linspace(-3.6, -2.9, 4)
	SEED = 10
	dU = 1
	dN = 1
	L = 6
	Nmax = 2*L*L*L
	Umin = Nmax*12*(-1.0)
	
	index = 80
	suffix = "a"
	tracker = open("tracker.dat", "a")
	for t in T:
		for mu in MU:
			os.makedirs ("./"+str(index)+str(suffix)+"/")
			os.chdir("./"+str(index)+str(suffix)+"/")
			file = open("run.sh", "w")
			file.write("../../../../bin/gcmc ==N:1 ==L:"+str(L)+" ==T:"+str(t)+" ==max1:"+str(Nmax)+" ==sweepSize:100 ==prodSweep:100000 ==equilSweep:3000 ==seed:"+str(SEED)+" ==mu1:"+str(mu)+"\n")
			file.close()
			file2 = open("analyze.py", "w")
			file2.write("import sys\nimport os\nsys.path.append(\"/Users/nathanmahynski/Documents/synced_folder/projects/active/star_gcmc/code/engine/tests/\")\nfrom mc_analysis_tools import * \n")
			file2.write("files = os.listdir(\"./\")\ninFile = \"none\"\nfor file in files:\n\tif (\"thermo\" in file):\n\t\tinFile = file\n\t\tbreak\n")
			file2.write("create_histogram (input=inFile, output=\"hist\", dim=2, min=[0, "+str(Umin)+"], bin_widths=["+str(dN)+","+str(dU)+"], L="+str(L)+", mu=["+str(mu)+"], T="+str(t)+", name=[\"number\", \"energy\"])\nconvert_hist_nam2azp (input=\"hist\", output=\"azpHist\")\n")
			file2.write("mc_stats (input=inFile, output=\"stats.dat\", args=2, blocks=3)\n")
			file2.close()
			os.chdir("../")
			tracker.write(str(index)+str(suffix)+"\t"+str(mu)+"\t"+str(t)+"\n")
			index += 1
	tracker.close()