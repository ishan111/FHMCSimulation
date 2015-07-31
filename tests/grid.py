import numpy as np
import os

if __name__ == "__main__":
	T = np.linspace(1.27, 1.40, 4)
	MU = np.linspace(-10, 5, 20)
	SEED = 10
	dU = 0.5
	dN = 1
	L = 10
	Nmax = 2000
	Umin = Nmax*12*(-1.0)
	
	index = 0
	for t in T:
		for mu in MU:
			#os.makedirs ("./"+str(index)+"/")
			os.chdir("./"+str(index)+"/")
			#file = open("run.sh", "w")
			#file.write("../../../bin/gcmc ==N:1 ==L:"+str(L)+" ==T:"+str(t)+" ==max1:"+str(Nmax)+" ==sweepSize:100 ==prodSweep:100000 ==equilSweep:3000 ==seed:"+str(SEED)+" ==mu1:"+str(mu)+"\n")
			#file.close()
			file2 = open("analyze.py", "w")
			file2.write("import sys\nimport os\nsys.path.append(\"/Users/nathanmahynski/Documents/synced_folder/projects/active/star_gcmc/code/engine/tests/\")\nfrom mc_analysis_tools import * \n")
			file2.write("files = os.listdir(\"./\")\ninFile = \"none\"\nfor file in files:\n\tif (\"thermo\" in file):\n\t\tinFile = file\n\t\tbreak\n")
			file2.write("create_histogram (input=inFile, output=\"hist\", dim=2, min=[0, "+str(Umin)+"], bin_widths=["+str(dN)+","+str(dU)+"], L="+str(L)+", mu=["+str(mu)+"], T="+str(t)+", name=[\"number\", \"energy\"])\nconvert_hist_nam2azp (input=\"hist\", output=\"azpHist\")")
			file2.close()
			os.chdir("../")
			index += 1