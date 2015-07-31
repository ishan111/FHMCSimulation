# organize the grid for gnuplot visualization
import os
import sys
import operator

try:
	file = open("tracker.dat", 'r')
except IOError:
	print "Error: couldn't open tracker.dat for organization."
	sys.exit()
	
# read tracker.dat and organize into a library based on unique tag
data = file.read().strip().split()
file.close()

# info[tag] = (mu, temp) 
info = {}
mu_sort = {}
for i in range(0, len(data)/3):
	info[data[3*i]] = [float(data[3*i+2]), float(data[3*i+1])]
	mu_sort[data[3*i]] = float(data[3*i+2])
	
# sort info based on colloid chemical potential (will maintain different entries for repeated muc values)
sorted_mu = sorted(mu_sort.iteritems(), key=operator.itemgetter(1))

# loop over all unique temps, create a new dictionary to sort for each unique temp, based on mu_gnp
TEMP = []

# get number of unique mu_gnp's
for i in range(0, len(sorted_mu)):
	if (i == 0):
		current_mu = sorted_mu[i][1]
		TEMP.append(current_mu)
	else:
		# check if we need to move to next mu
		if (sorted_mu[i][1] != current_mu):
			current_mu = sorted_mu[i][1]
			TEMP.append(current_mu)

# create a dictionary for each unique mu_gnp
t_sort = []
sorted_t = []
for i in range(0, len(TEMP)):
	t_sort.append({})
	sorted_t.append([])

# fill those dictionaries
index = 0
for i in range(0, len(sorted_mu)):
	if (i == 0):
		current_mu = TEMP[index]
	else:
		if (sorted_mu[i][1] != current_mu):
			index += 1
			current_mu = TEMP[index]
	
	t_sort[index][sorted_mu[i][0]] = info[sorted_mu[i][0]][1]

# sort each of the dictionaries
for i in range(0, len(TEMP)):
	sorted_t[i] = sorted(t_sort[i].iteritems(), key=operator.itemgetter(1))

# now sorted in order of ascending mu and temp (for each mu) 
for i in range(0, len(TEMP)):
	# create a new file for each mu
	try:
		ifile = open("temp_"+str(TEMP[i])+"_summary.dat", 'w')
	except IOError:
		print "Error: couldn't open summary file for T = "+str(TEMP[i])
		sys.exit()
		
	# write header
	ifile.write("# T = "+str(TEMP[i])+"\n")
	ifile.write("# mu\tN\t+/-\t<E>\t+/-\n")
	
	# get rest of sorted data
	for j in range(0, len(sorted_t[i])):
		try:
			file = open(sorted_t[i][j][0]+"/stats.dat", 'r')
		except:
			print sorted_t[i][j][0]+" has no stats yet"
			continue
		
		newdat = file.read().strip().split()
		file.close()
		
		ifile.write(str(sorted_t[i][j][1])+"\t"+newdat[len(newdat)-4]+"\t"+newdat[len(newdat)-3]+"\t"+newdat[len(newdat)-2]+"\t"+newdat[len(newdat)-1]+"\n")		
		
	ifile.close()	

# provide a gnuplot script to plot this data
try:
	file = open("plot_summary.dat", 'w')
except IOError:
	print "Error: couldn't open plot_summary.dat."
	sys.exit()
	
file.write("set term post enh color solid \"Times-Roman\" 16\n")
file.write("set output \"results.eps\"\n\n")
file.write("set xlabel \"{/Symbol m}\"\n")
file.write("set ylabel \"N\"\n")
file.write("plot ")
for i in range(0, len(TEMP)):
	if (i != 0):
		file.write(",\t")
	file.write("\"temp_"+str(TEMP[i])+"_summary.dat\" u 1:2:3 with yerrorbars title\"T = "+str(TEMP[i])+"\"")
file.close()

