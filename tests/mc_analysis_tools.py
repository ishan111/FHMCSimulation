## @package mc_analysis_tools
#	Standard toolkit for post-processing monte carlo simulations including histograms and statistics
#

import sys
import numpy as np
import os

## Reads in a file with columns of instantaneous measurements and produces a histogram file in a "C"-like array format
#	NOTE: this does NOT check that there aren't gaps inside the histogram which could cause problems with reweighting
#		This can be fixed by either coarsening the bins, or by measuring the results at EVERY step in the MC program
#		This program operates in its local directory.  This should be viewed as an "emergency" backup scheme since it is not very efficient.
#
#	@param **kwargs includes:
#\verbatim
#	input			::	name distinguishing the input file(s), there may be multiple and will be sorted
#	output			::	name of the output file
#	dim			::	number of columns (dimensions) in file
#	bin_widths		::	widths of each dimensions bins
#	L		::	box size
#	mu		::	chemical potential(s)
#	T		::	temperature
#	min			:: (optional) set the start of each dimension's beginning at min[dim], must be specified after dim.  Use this if you need to create consistent histograms between different runs
#\endverbatim
def create_histogram (**kwargs):
	iset = 0
	for key in kwargs:
		if key == "input":
			filetag = kwargs[key]
		elif key == "output":
			outfile = kwargs[key]
		elif key == "L":
			L = kwargs[key]
		elif key == "mu":
			mu = kwargs[key]
		elif key == "T":
			temperature = kwargs[key]
		elif key == "name":
			colname = kwargs[key]
		elif key == "dim":
			dim = kwargs[key]
			if (iset != 0):
				print "Must specify min after dim.\n"
				sys.exit()
			else:
				mymin = []
				for i in range(0, dim):
					mymin.append(1.0e9)
		elif key == "bin_widths":
			bw = kwargs[key]
		elif key == "min":
			mymin = kwargs[key]	
			iset = 1
		else:
			print "Unknown keyword "+str(key)
			sys.exit()
	
	# notify
	print "\n=========================================================\n"
	print "Obtaining Histogram for *"+filetag+"*"
	print "\n=========================================================\n\n"
	
	# check the output file is valid
	icheck = 0
	try:
		out = open(outfile, 'r')
	except IOError:
		icheck = 1

	if (icheck == 0):
		print "Output file "+outfile+" already exists.\n"
		sys.exit()

	# find all the input files
	a = os.listdir(".")
	
	# list and sort the files
	ifiles = []
	for name in a:
		if filetag in name:
			ifiles.append(name)
	ifiles.sort()
	
	# allocate
	imax = []
	imin = []
	nbins = []
	for i in range(0, dim):
		imax.append(-1.0e9)
		imin.append(1.0e9)
		nbins.append(0)

	# scan through all files to find limits and allocate appropriately
	for filename in ifiles:
		# open file
		try:
			file = open(filename, 'r')
		except IOError:
			print "Could not open "+filename
			sys.exit()
		
		# check
		if (len(bw) != dim):
			print "Error: number of widths specified != dim in "+filename
			sys.exit()
	
		# read the contents
		raw_data = file.read().strip().split()
		file.close()
	
		# check the data
		if (len(raw_data)%dim != 0):
			print "File appears incomplete.\n"
			sys.exit()
		
		# allocate, do this for every file and do NOT store since this may overload python with data
		data = []
		for i in range(0, dim):
			data.append([])
		
		# parse the contents
		index = 0
		for i in range(0, len(raw_data)/dim):
			for j in range(0, dim):
				data[j].append(float(raw_data[index]))
				index += 1

		# find limits and number of bins to specify
		for i in range(0, dim):
			if (np.max(data[i]) > imax[i]):
				imax[i] = np.max(data[i])
			if (np.min(data[i]) < imin[i]):
				imin[i] = np.min(data[i])
			if (iset == 1):
				if (imin[i] < mymin[i]):
					print "Cannot use min value for dimension "+str(i+1)+"; smaller value found, "+str(imin[i])+", in "+filename
					sys.exit()
				else:
					imin[i] = mymin[i]
	
	# calculate bins and total dimension
	tot_dim = 1
	for i in range(0, dim):	
		# for rounding purposes, we add one (i.e. if min-max is integer multiple of bin width, max value will be out of bin range)
		nbins[i] = int(np.ceil((float(imax[i])-float(imin[i]))/bw[i]))+1
		tot_dim = tot_dim*nbins[i]

	# allocate histogram
	hist = []
	for i in range(0, tot_dim):
		hist.append(0)

	# go through each file
	for filename in ifiles:
		# get data from the file
		file = open(filename, 'r')
		raw_data = file.read().strip().split()
		file.close()

		# allocate
		data = []
		for i in range(0, dim):
			data.append([])
		
		# parse the contents
		index = 0
		for i in range(0, len(raw_data)/dim):
			for j in range(0, dim):
				data[j].append(float(raw_data[index]))
				index += 1
	
		# sort the data into a histogram in "C" indexed format, that is indexing by 1st dimension, 2nd dimension, etc.
		for i in range(0, len(raw_data)/dim):
			# find the correct index
			index = 0
			factor = 1
			bin = []

			for j in range(0, dim):
				#bin.append(int(np.floor(float(data[j][i]-imin[j])/float(bw[j]))))
				bin.append(int(np.floor(float(data[j][i]-(imin[j]-0.5*bw[j]))/float(bw[j]))))
				if (bin[len(bin)-1] >= nbins[j]):
					bin[len(bin)-1] = nbins[j]-1
				index += bin[j]*factor
				factor = factor*nbins[j]
			hist[index] += 1
		
	# output the data 
	try:
		out = open(outfile, 'w')
	except IOError:
		print "Could not open "+outfile
		sys.exit()
	
	out.write("Sweep: "+str(len(raw_data)/dim)+"\n")
	out.write("Box: "+str(L)+"\t"+str(L)+"\t"+str(L)+"\n")
	out.write("Mu: ")
	for m in mu:
		out.write(str(m)+"\t")
	out.write("\n")
	out.write("Temp: "+str(temperature)+"\n")
	out.write("Dim\tName\tMax\tMin\tBins\n")
	for i in range(0, dim):
		# DO NOT REPORT IMAX, for other routines to work, the max has to report the bin centers of the histogram not the actual value inside
		# this is because bin_widths not nbins was taken as input.  for programs which set nbins and calculate the width this may be different
		# the above routines took as a lower bound reference imin[i] so this is the center of the lower bin
		# the upper bin center is calculated relative to that
		# now we can take a raw set of data, convert it to "my" histogram format, then convert to azp with the code later correctly
		# this will now work for the case where we do not "fix the grid" as would happen if histogram coming directly from code or if mymin option used above
		out.write(str(i+1)+"\t"+str(colname[i])+"\t"+str(imin[i]+(nbins[i]-1)*bw[i])+"\t"+str(imin[i])+"\t"+str(nbins[i])+"\n")
	out.write("\n")
	for i in range(0, tot_dim):
		out.write(str(hist[i])+"\t")
	
	out.close()
	
## Statistics for a MC run
#
#	@param **kwargs includes:
#\verbatim
#	input		::	characteristic name of input file(s) (there may be multiple which will be sorted in order)
#	output		::	output filename
#	args		::	number of columns in the input file
#	blocks		::	number of statistical blocks
#\endverbatim
def mc_stats (**kwargs):	
	for key in kwargs:
		if key == "input":
			filetag = kwargs[key]
		elif key == "output":
			outfile = kwargs[key]
		elif key == "args":
			args = kwargs[key]
			if (args < 1):
				print "args < 1.\n"
				sys.exit()
		elif key == "blocks":
			blocks = kwargs[key]
			if (blocks < 1):
				print "blocks < 2.\n"
				sys.exit()
		else:
			print "Unknown keyword "+str(key)
			sys.exit()
		
	# notify
	print "\n=========================================================\n"
	print "Obtaining Monte Carlo Statistics for *"+filetag+"*"
	print "\n=========================================================\n\n"

	# check the output file is valid
	icheck = 0
	try:
		out = open(outfile, 'r')
	except IOError:
		icheck = 1

	if (icheck == 0):
		print "Output file "+outfile+" already exists.\n"
		sys.exit()

	# find all the input files
	a = os.listdir(".")
	
	# list and sort the files
	ifiles = []
	for name in a:
		if filetag in name:
			ifiles.append(name)
	ifiles.sort()	
		
	# loop through files to find the total number of entries (store as double because may overflow integer storage capacity)
	tot_entries = 1.0
	for filename in ifiles:
		# open file
		try:
			file = open(filename, 'r')
		except IOError:
			print "Could not open "+filename
			sys.exit()
	
		print "Checked "+filename
		# read the contents
		raw_data = file.read().strip().split()
		file.close()
	
		# check the data
		ilen = len(raw_data)
		if (ilen%args != 0):
			print "File appears incomplete.\n"
			sys.exit()
		
		tot_entries += ilen/args
	
	# find size of each block (is a double)
	nper = np.floor(tot_entries/blocks)

	# allocate
	bave = []
	nave = []
	stdev = []
	for i in range(0, args):
		bave.append([])
		nave.append(0.0)
		stdev.append(0.0)
		for j in range(0, blocks):
			bave[i].append(0.0)
			
	# loop through files to find the average
	tot_iter = 0.0
	iter = 0.0
	iblock = 0
	iquit = False
	for filename in ifiles:
		# open file and read the contents
		print "Analyzing "+filename
		file = open(filename, 'r')
		raw_data = file.read().strip().split()
		file.close()	
		
		# parse the contents
		index = 0
		for i in range(0, len(raw_data)/args):
			for j in range(0, args):
				bave[j][iblock] += float(raw_data[index])
				index += 1
			
			# check the block
			if (iter+1.0 >= nper):
				iter = 0.0
				iblock += 1
			else:
				iter += 1.0
				
			# check if we have reached the end of where we want to record
			if (tot_iter+1.0 >= nper*blocks):
				iquit = True
				break
			tot_iter += 1.0
			
		if (iquit):
			break
			
	# average
	for i in range(0, args):
		for j in range(0, blocks):
			bave[i][j] = bave[i][j]/nper
			nave[i] += bave[i][j]/float(blocks)
	
	# get standard error
	for i in range(0, args):
		for j in range(0, blocks):
			stdev[i] += (nave[i]-bave[i][j])**2
		stdev[i] = 1.96/(blocks)**(1.0/2.0)*(stdev[i]/float(blocks-1))**(1.0/2.0)	
		
	# report
	try:
		out = open(outfile, 'w')
	except IOError:
		print "Could not open "+outfile
		sys.exit()
		
	out.write("Block:\t")
	for i in range(0, args):
		out.write("Prop "+str(i+1)+"\t+/-\t")
	out.write("\n")
	for i in range(0, blocks):
		out.write(str(i+1)+"\t")
		for j in range(0, args):
			out.write(str(bave[j][i])+"\t\t")
		out.write("\n")	
	for i in range(0, args):	
		out.write("----------------\t----------------\t----------------\t")
	out.write("\n\t")
	for i in range(0, args):
		out.write(str(nave[i])+"\t"+str(stdev[i])+"\t")
	
	out.close()
	
## Takes columns of data and produces a histogram in "AZP" format for binary system ONLY, so that it can be read with "entropy" series
#	This DOES assume there is a third column floating around that is from the energy (spac program's output)
#	This should be viewed as an "emergency backup" method since this routine is not very efficient
#
# @param **kwargs includes:
#\verbatim
#	input			::	tag for input file(s) (may be multiple which will be sorted)
#	output			:: 	output filename
#	bin_widths		::	bin widths for histogram for each dimension (2)
#	mu				::	chemical potentials (2)
#	box				::	box size (3)
#\endverbatim
def azp_hist (**kwargs):
	for key in kwargs:
		if key == "input":
			filetag = kwargs[key]
		elif key == "output":
			outfile = kwargs[key]
		elif key == "bin_widths":
			bw = kwargs[key]
		elif key == "mu":
			mu = kwargs[key]
		elif key == "box":
			box = kwargs[key]	
		else:
			print "Unknown keyword "+str(key)
			sys.exit()
	
	# notify
	print "\n=========================================================\n"
	print "Obtaining AZP-style histogram for *"+filetag+"*"
	print "\n=========================================================\n\n"
	
	# check
	if (len(bw) != 2):
		print "Error: number of widths specified != dim."
		sys.exit()
	if (len(mu) != 2):
		print "Error: not enough chemical potentials specified."
		sys.exit()
	if (len(box) != 3):
		print "Error: not enough box dimensions specified."
		sys.exit()
			
	# check the output file is valid
	icheck = 0
	try:
		out = open(outfile, 'r')
	except IOError:
		icheck = 1

	if (icheck == 0):
		print "Output file "+outfile+" already exists.\n"
		sys.exit()

	# find all the input files
	a = os.listdir(".")
	
	# list and sort the files
	ifiles = []
	for name in a:
		if filetag in name:
			ifiles.append(name)
	ifiles.sort()
	
	# allocate
	imax = []
	imin = []
	nbins = []
	for i in range(0, 2):
		imax.append(-1.0e9)
		imin.append(1.0e9)
		nbins.append(0)
	
	# loop through files to get limits
	for filename in ifiles:
		# open file
		try:
			file = open(filename, 'r')
		except IOError:
			print "Could not open "+filename
			sys.exit()
	
		print "Checked "+filename
		# read the contents
		raw_data = file.read().strip().split()
		file.close()
	
		# check the data
		if (len(raw_data)%(2+1) != 0):
			print "File appears incomplete.\n"
			sys.exit()
	
		# parse the contents
		index = 0
		data = []
		for i in range(0, 2):
			data.append([])
		for i in range(0, len(raw_data)/(2+1)):
			for j in range(0, 2):
				data[j].append(float(raw_data[index]))
				index += 1
			index += 1
		
		# convert to np array for speed
		data = np.array(data)
	
		# find limits and number of bins to specify
		for	i in range(0, 2):
			if (imax[i] < np.max(data[i])):
				imax[i] = np.max(data[i])
			if (imin[i] > np.min(data[i])):
				imin[i] = np.min(data[i])
	
	# get bins
	for	i in range(0, 2):	
		# for rounding purposes, we add one (i.e. if min-max is integer multiple of bin width, max value will be out of bin range)
		nbins[i] = int(np.ceil(float(imax[i]-imin[i])/bw[i]))+1
	
	# allocate and initialize histogram
	hist = []
	for i in range(0, nbins[0]):
		hist.append([])
		for	j in range(0, nbins[1]):
			hist[i].append(0)
	
	# loop through files to get limits
	for filename in ifiles:
		# open file and read the contents
		file = open(filename, 'r')
		raw_data = file.read().strip().split()
		file.close()
	
		print "Analyzing "+filename
		# parse the contents
		index = 0
		data = []
		for i in range(0, 2):
			data.append([])
		for i in range(0, len(raw_data)/(2+1)):
			for j in range(0, 2):
				data[j].append(float(raw_data[index]))
				index += 1
			index += 1
		
		# convert to np array for speed
		data = np.array(data)
	
		# sort the data into a histogram into AZP format (2D grid)
		index = 0
		for i in range(0, len(raw_data)/(2+1)):
			hist[int(np.floor(float(data[0][index]-imin[0])/bw[0]))][int(np.floor(float(data[1][index]-imin[1])/bw[1]))] += 1
			index += 1
	
	# output the data
	try:
		out = open(outfile, 'w')
	except IOError:
		print "Could not open "+outfile
		sys.exit()
	
	out.write("mu1   mu2     x- y- zdim  ! note that N->N1, U->N2\n")
	for i in range(0, 2):
		out.write(str(mu[i])+"\t")
	out.write(str(bw[1])+"\t")
	for i in range(0, 3):
		out.write(str(box[i])+"\t")
	out.write("\n")
	
	for i in range(0, nbins[0]):
		# find range
		istart = 0
		iend = nbins[1]-1
		if1 = True
		if2 = True
		for j in range(0, nbins[1]):
			if (hist[i][j] == 0 and if1):
				istart += 1
			else:
				if1 = False
			if (hist[i][nbins[1]-1-j] == 0 and if2):
				iend -= 1
			else:
				if2 = False
			
		# check range
		if (istart > iend):
			print "Error in finding non-zero center N2 = ("+str(int(imin[1]+istart*bw[1]))+","+str(int(imin[1]+iend*bw[1]))+") : N1 = "+str(int(imin[0]+i*bw[0]))+".\n"
			sys.exit()
			
		# check no zeros in between
		for j in range(istart, iend+1):
			if (hist[i][j] == 0):
				print "Warning: zeros between bounds N2 = ("+str(int(imin[1]+istart*bw[1]))+","+str(int(imin[1]+iend*bw[1]))+") : N1 = "+str(int(imin[0]+i*bw[0]))+".\n"
				#sys.exit()
		
		# write
		out.write(str(int(imin[0]+i*bw[0]))+"\t"+str(iend-istart+1)+"\t"+str(int(imin[1]+istart*bw[1]))+"\n")
		for j in range(istart, iend+1):
			out.write(str(hist[i][j])+"\t")
		out.write("\n")
		
	out.close()	
	
## Takes C-style histogram from spac's output and produces an "azp"-style histogram for reweighting (based on first two dimensions it finds in C-file)
#
# @param **kwargs includes:
#\verbatim
#	input			::	C-style histogram name
#	output			:: 	output filename
#\endverbatim
def convert_hist_nam2azp (**kwargs):
	for key in kwargs:
		if key == "input":
			filename = kwargs[key]
		elif key == "output":
			outfile = kwargs[key]	
		else:
			print "Unknown keyword "+str(key)
			sys.exit()

	# open file
	try:
		file = open(filename, 'r')
	except IOError:
		print "Could not open "+filename
		sys.exit()
	
	# read header and get info
	file.readline()
	a = file.readline().strip().split()
	box = []
	for i in range(0, 3):
		box.append(float(a[1+i]))
	a = file.readline().strip().split()
	mu = []
	for i in range(0, len(a)-1):
		mu.append(float(a[1+i]))	
	a = file.readline().strip().split()
	temperature = float(a[len(a)-1])
	file.readline()

	igo = True
	ndims = 0
	dim_min = []
	dim_max = []
	dim_bins = []
	name = []
	bw = []
	tot_dim = 1
	while (igo):
		a = file.readline().strip().split()
		if (a == []):
			igo = False
		else:
			ndims += 1
			dim_max.append(float(a[2]))
			dim_min.append(float(a[3]))
			dim_bins.append(int(a[4]))
			if (dim_bins[ndims-1] == 1):
				bw.append(1.0)	# default bin width to unity if only 1 bin of it
			else:
				bw.append((dim_max[ndims-1]-dim_min[ndims-1])/float(dim_bins[ndims-1]-1))
			name.append(a[1])
			tot_dim = tot_dim*dim_bins[ndims-1]
	
	# allocate and initialize hist
	hist = []
	for i in range(0, dim_bins[0]):
		hist.append([])
		for j in range(0, dim_bins[1]):
			hist[i].append(0)
			
	print "Using "+name[0]+" and "+name[1]+" for histogram."
	
	hist_data = file.readline().strip().split()
	tot_index = 0
	while (hist_data != []):
		# analyze current line
		index = 0
		for i in range(0, len(hist_data)):
			pos = get_array_loc(tot_index, ndims, tot_dim, dim_bins)
			# regardless of other "dimensions", all the data is collapsed using the first two dimensions only
			hist[pos[0]][pos[1]] += int(hist_data[index])
			index += 1
			tot_index += 1
			
		# read the next line
		hist_data = file.readline().strip().split()

	# print to file
	try:
		out = open(outfile, 'w')
	except IOError:
		print "Could not open "+outfile
		sys.exit()
	
	out.write("mu1(or Temp)   mu2(or mu1)     x- y- zdim  ! note that N->N1, U->N2 (T specified first then)\n")
	# azp needs T specified first
	if (len(mu) == 1):
		out.write(str(temperature)+"\t")
	for i in range(0, len(mu)):
		out.write(str(mu[i])+"\t")
	out.write(str(bw[1])+"\t")
	for i in range(0, 3):
		out.write(str(box[i])+"\t")
	out.write("\n")
	
	ican = 1
	init = 0
	for i in range(0, dim_bins[0]):
		# find range
		istart = 0
		iend = dim_bins[1]-1
		if1 = True
		if2 = True
		for j in range(0, dim_bins[1]):
			if (hist[i][j] == 0 and if1):
				istart += 1
			else:
				if1 = False
			if (hist[i][dim_bins[1]-1-j] == 0 and if2):
				iend -= 1
			else:
				if2 = False
			
		# check range
		if (istart > iend and (istart != dim_bins[1] or iend != -1)):
			print "Error in finding non-zero center "+name[1]+" = ("+str(int(dim_min[1]+istart*bw[1]))+","+str(int(dim_min[1]+iend*bw[1]))+") : "+name[0]+" = "+str(int(dim_min[0]+i*bw[0]))
			sys.exit()
		
		"""
		# check no zeros in between if found a non-zero row (ie iend >= istart)
		for j in range(istart, iend+1):
			if (hist[i][j] == 0):
				print "Warning: zeros between bounds N2 = ("+str(int(dim_min[1]+istart*bw[1]))+","+str(int(dim_min[1]+iend*bw[1]))+") : N1 = "+str(int(dim_min[0]+i*bw[0]))
				#sys.exit()
		"""
		# write if have a non-zero row
		if (iend >= istart):
			if (ican == 0 and init == 1):
				print "Warning, skipping row of "+name[0]+" between "+str(i+1)+" and "+str(i+2)
			#out.write(str(int(dim_min[0]+i*bw[0]))+"\t"+str(iend-istart+1)+"\t"+str(int(dim_min[1]+istart*bw[1]))+"\n")
			out.write(str(int(dim_min[0]+i*bw[0]))+"\t"+str(iend-istart+1)+"\t"+str(dim_min[1]+istart*bw[1])+"\n")
			for j in range(istart, iend+1):
				out.write(str(hist[i][j])+"\t")
			out.write("\n")
			ican = 1
			init = 1
		else:
			ican = 0
	
## Find the coordinates of a multidimensional array, given the equivalent index in a 1-D array
#
#	@param index known index
#	@param ndims number of dimensions
#	@param tot_dim length of 1D array
#	@param nbins length of each dimension (remember to include 0 if it is within range)
def get_array_loc (index, ndims, tot_dim, nbins):
	ans = []
	for i in range(0, ndims):
		ans.append(0)
	
	tot_dim = tot_dim/nbins[ndims-1]
	for i in range(0, ndims):
		ans[ndims-1-i] = int(np.floor(index/(float(tot_dim))))
		index -= tot_dim*ans[ndims-1-i]
		if (i < ndims-1):
			tot_dim = tot_dim/nbins[(ndims-1)-(i+1)]
			
	if (index != 0):
		print "Error in finding dimensions.\n"+str(index)+str(ans)
		sys.exit()
	else:
		return ans
		
