#include "bias.h"

/*!
 * Initialize the tmmc object.
 * 
 * \param [in] nSpec Number of species in the simulation.
 * \param [in] Nmax Vector of upper bound for number of particles of each species.
 * \param [in] Nmin Vector of lower bound for number of particles of each species.
 */
tmmc::tmmc (const int nSpec, const std::vector <int> &Nmax, const std::vector <int> &Nmin) {
	if (nSpec > 0) {
		nSpec_ = nSpec;
	} else {
		throw customException ("Number of species in simulation must be > 0, cannot initialize tmmc");
	}
	
	if (Nmax.size() != Nmin.size()) {
		throw customException ("Nmin and Nmax vectors have unequal sizes, cannot initialize tmmc");
	}
	
	if (Nmax.size() != nSpec_) {
		throw customException ("Nmin and Nmax must specify bounds for all species, cannot initialize tmmc");
	}
	
	W_.resize(nSpec_, 0);
	
	__BIAS_INT_TYPE__ size = 1;
	for (unsigned int i = 0; i < nSpec_; ++i) {
		if (Nmin[i] > Nmax[i]) {
			throw customException ("Nmin > Nmax for species "+sstr(i+1));
		}
		W_[i] = (Nmax[i] - Nmin[i] +1);
		size *= W_[i];
	}
	
	Nmin_ = Nmin;
	Nmax_ = Nmax;
	
	// attempt to allocate memory for collection matrix and initializes it all to 0
	try {
		C_.resize((1+2*nSpec_)*size, 0);
	} catch (const std::bad_alloc &ce) {
		throw customException ("Out of memory, cannot allocate space for collection matrix in tmmc");
	}

	// attempt to allocate memory for probability matrix and initializes it all to 0
	try {
		P_.resize((1+2*nSpec_)*size, 0);
	} catch (const std::bad_alloc &ce) {
		throw customException ("Out of memory, cannot allocate space for probability matrix in tmmc");
	}
	
	// attempt to allocate memory for lnPI matrix and initializes it all to 0
	try {
		lnPI_.resize(size, 0.0);
	} catch (const std::bad_alloc &ce) {
		throw customException ("Out of memory, cannot allocate space for macrostate distribution matrix in tmmc");
	}
}

/*!
 * Get the difference offset between two states of the simulation. 
 * This assumes ONLY one species has been modified and by a value of 1.  
 * Multiple changes, or those > 1, WILL result in unexpected behavior.
 * Default resets specId and addOrSubtract to 0 when called.
 * 
 * \param [in] Nstart Vector of the initial number of atoms of each type (in order) in the system
 * \param [in] Nend Vector of the final number of atoms of each type (in order) in the system
 * \param [out] specId Index of the first species which has changed between Nstart and Nend
 * \param [out] addOrSubtract +1 if specId has increased, -1 if it has decreased, and 0 if it has not changed
 */
void tmmc::difference_ (const std::vector <int> &Nstart, const std::vector <int> &Nend, int &specId, int &addOrSubtract) {
	specId = 0;
	addOrSubtract = 0;
	for (unsigned int i = 0; i < Nend.size(); ++i) {
		addOrSubtract = Nend[i] - Nstart[i];
		specId = i;
		if (addOrSubtract != 0) {
			// break once the (first) difference has been found
			break;
		}
	}
	return;
}

/*!
 * For a given multidimensional array which has been cast into 1D, find the address that refers to a given transition.
 * 
 * \param [in] Nstart Vector of the number of each species (in order) initially (before MC move)
 * \param [in] Nend Vector of the number of each species (in order) after the MC move
 */
const __BIAS_INT_TYPE__ tmmc::getTransitionAddress (const std::vector <int> &Nstart, const std::vector <int> &Nend) {
	// Layout of y = [0, +1, -1, +1, -1, ... ] where each pair of +1/-1 is ordered from specId = 0 to N-1
	int specId = 0, addOrSubtract = 0;
	difference_(Nstart, Nend, specId, addOrSubtract);
			
	int y = 0;
	if (addOrSubtract != 0) {
		y += 2*(specId+1);
		if (addOrSubtract == 1) {
			y -= 1;
		} else {
			throw customException ("Bad addOrSubtract value: "+sstr(addOrSubtract));
		}
	}
	
	__BIAS_INT_TYPE__ x = (Nstart[nSpec_-1] - Nmin_[nSpec_-1]);
	for (int i = nSpec_-2; i >= 0; --i) {
		x *= W_[i];
		x += (Nstart[i] - Nmin_[i]);
	}
	
	return x*(1+2*nSpec_) + y;
}

/*!
 * For a given multidimensional array which has been cast into 1D, find the address that refers to a given transition. 
 * This is overloaded to allow for the calculation of the differences between initial and final states previously without repeating the calculation.
 * 
 * \param [in] Nstart Vector of the number of each species (in order) initially (before MC move)
 * \param [in] specId Index of species that was changed
 * \param [in] addOrSubtract Indicates if specId was increased (+1), decreased (-1), or no change (0)
 */
const __BIAS_INT_TYPE__ tmmc::getTransitionAddress (const std::vector <int> &Nstart, const int specId, const int addOrSubtract) {		
	int y = 0;
	if (addOrSubtract != 0) {
		y += 2*(specId+1);
		if (addOrSubtract == 1) {
			y -= 1;
		} else {
			throw customException ("Bad addOrSubtract value: "+sstr(addOrSubtract));
		}
	}
	
	__BIAS_INT_TYPE__ x = (Nstart[nSpec_-1] - Nmin_[nSpec_-1]);
	for (int i = nSpec_-2; i >= 0; --i) {
		x *= W_[i];
		x += (Nstart[i] - Nmin_[i]);
	}
	
	return x*(1+2*nSpec_) + y;
}

/*!
 * Get the address in lnPI that corresponds to a given macrostate.
 * 
 * \param [in] Nval Vector of the number of each species (in order)
 */
const __BIAS_INT_TYPE__ tmmc::getAddress (const std::vector <int> &Nval) {		
	__BIAS_INT_TYPE__ x = (Nval[nSpec_-1] - Nmin_[nSpec_-1]);
	for (int i = nSpec_-2; i >= 0; --i) {
		x *= W_[i];
		x += (Nval[i] - Nmin_[i]);
	}
	return x;
}

/*!
 * Update the collection matrix.
 * 
 * \param [in] Nstart Vector of the number of each species (in order) initially (before MC move)
 * \param [in] Nend Vector of the number of each species (in order) after the MC move
 * \param [in] pa Unbiased Metropolis criterion for making a MC move (i.e. pa = min(1, exp(...)))
 */
void tmmc::updateC (const std::vector <int> &Nstart, const std::vector <int> &Nend, const double pa) {
	int specId = 0, addOrSubtract = 0;
	difference_ (Nstart, Nend, specId, addOrSubtract);
	const __BIAS_INT_TYPE__ i = getTransitionAddress (Nstart, specId, addOrSubtract);
	__BIAS_INT_TYPE__ j = i;
	if (addOrSubtract != 0) {
		j = getTransitionAddress (Nstart, specId, 0);
	} 
	C_[i] += pa;
	C_[j] += (1-pa);
}

/*!
 * Calculate the (natural logarithm of the) macrostate density matrix via the probability matrix.
 */
void tmmc::calculatePI () {
	for (__BIAS_INT_TYPE__ i = 0; i < C_.size(); i += (1+2*nSpec_)) {
		double sum = 0.0;
		for (unsigned int j = 0; j < 1+2*nSpec_; ++j) {
			sum += C_[i+j];
		}
		if (sum > 0) {
			for (unsigned int j = 0; j < 1+2*nSpec_; ++j) {
				P_[i+j] = C_[i+j] / sum;
			}
		} else {
			for (unsigned int j = 0; j < 1+2*nSpec_; ++j) {
				P_[i+j] = 0;
			}
		}
	}
	
	// Reset first value to zero just to start fresh. Since only ratios matter this is perfectly fair.
	lnPI_[0] = 0.0;
	if (nSpec_ == 1) {
		std::vector <int> N1 (1), N2 (1);
		__BIAS_INT_TYPE__ address1, address2;
		for (__BIAS_INT_TYPE__ i = 0; i < lnPI_.size()-1; ++i) {
			N1[0] = Nmin_[0] + i;
			address1 = getTransitionAddress (N1, 0, +1);
			N2[0] = N1[0] + 1;
			address2 = getTransitionAddress (N2, 0, -1);
			lnPI_[i+1] = lnPI_[i] + log(P_[address1]/P_[address2]);
		}
	} else {
		// not built yet
		exit(-1);
	}
}

/*!
 * Print the UN-NORMALIZED biasing function (lnPI) and possibly collection matrix to files. 
 * Will overwrite the files if another with that name exists. 
 * Prints in netCDF format if enabled.
 * 
 * \param [in] fileName Name of the file to print to.  Will append with "_lnPI" and "_C" for biasing function and collection matrix, respectively.
 * \param [in] printC Defaults to false, but if true will also print the collection matrix.
 */
void tmmc::print (const std::string fileName, bool printC) {
#ifdef NETCDF_CAPABLE
	// Print collection matrix
	if (printC) {
		const std::string name = fileName + "_C.nc"
		NcFile outFile(name.c_str(), NcFile::replace);
		NcDim probDim = outFile.addDim("vectorized_position", C_.size());
		NcVar probVar = outFile.addVar("C", ncDouble, probDim);
		const std::string dummyName = "number_species:";
		probVar.putAtt(dummyName.c_str(), sstr(nSpec_).c_str());
		for (unsigned int i = 0; i < nSpec_; ++i) {
	  	    const std::string attName = "species_"+sstr(i+1)+"_upper_bound:";
	    	probVar.putAtt(attName.c_str(), sstr(Nmax_[i]).c_str());
		}
		for (unsigned int i = 0; i < nSpec_; ++i) {
		    const std::string attName = "species_"+sstr(i+1)+"_lower_bound:";
	     	probVar.putAtt(attName.c_str(), sstr(Nmin_[i]).c_str());
		}
		probVar.putVar(&C_[0]);
	}
	
	// Print lnPI (bias) matrix
	const std::string name = fileName + "_lnPI.nc"
	NcFile outFile(name.c_str(), NcFile::replace);
	NcDim probDim = outFile.addDim("vectorized_position", lnPI_.size());
	NcVar probVar = outFile.addVar("lnPI", ncDouble, probDim);
	const std::string dummyName = "number_species:";
	probVar.putAtt(dummyName.c_str(), sstr(nSpec_).c_str());
	for (unsigned int i = 0; i < nSpec_; ++i) {
  	    const std::string attName = "species_"+sstr(i+1)+"_upper_bound:";
    	probVar.putAtt(attName.c_str(), sstr(Nmax_[i]).c_str());
	}
	for (unsigned int i = 0; i < nSpec_; ++i) {
	    const std::string attName = "species_"+sstr(i+1)+"_lower_bound:";
     	probVar.putAtt(attName.c_str(), sstr(Nmin_[i]).c_str());
	}
	probVar.putVar(&lnPI_[0]);
#else
	// Print collection matrix
	if (printC) {
		std::ofstream of;
		of.open(fileName+"_C.dat", 'w');
		of << "# Collection matrix in single row (vectorized) notation." << std::endl;
		of << "# Number of species:" << nSpec_ << std::endl;
		for (unsigned int i = 0; i < nSpec_; ++i) {
			of << "# species_"+sstr(i+1)+"_upper_bound: " << Nmax_[i] << std::endl;
		}
		for (unsigned int i = 0; i < nSpec_; ++i) {
			of << "# species_"+sstr(i+1)+"_lower_bound: " << Nmin_[i] << std::endl;
		}
		for (long long int i = 0; i < C_.size(); ++i) {
			of << C_[i] << std::endl;
		}
		of.close();
	}
	
	// Print lnPI (bias) matrix
	std::ofstream of;
	of.open(fileName+"_lnPI.dat", 'w');
	of << "# lnPI (bias) matrix in single row (vectorized) notation." << std::endl;
	of << "# Number of species:" << nSpec_ << std::endl;
	for (unsigned int i = 0; i < nSpec_; ++i) {
		of << "# species_"+sstr(i+1)+"_upper_bound: " << Nmax_[i] << std::endl;
	}
	for (unsigned int i = 0; i < nSpec_; ++i) {
		of << "# species_"+sstr(i+1)+"_lower_bound: " << Nmin_[i] << std::endl;
	}
	for (long long int i = 0; i < lnPI_.size(); ++i) {
		of << lnPI_[i] << std::endl;
	}
	of.close();
#endif
}

/*!
 * Read the collection matrix from a file.  This assumes the user has already guaranteed that the bounds are consistent,
 * e.g. Nmin and Nmax, as it will not check this automatically.  Also assumes file was generated by this code.  "Hand made"
 * ones might have formatting issues since parsing is done based on tokens.
 * 
 * \param [in] fileName Name of file containing the collection matrix.  Must include file extension.
 */
void tmmc::readC (const std::string fileName) {
#ifdef NETCDF_CAPABLE
	NcFile dataFile (fileName.c_str(), NcFile::read);
	NcVar C_data = dataFile.getVar("C");
	if (C_data.isNull()) throw customException("Collection matrix was empty, cannot read");
	C_data.getVar(&C_[0]);
#else
	std::string line;
	std::ifstream inF (fileName);
	
	// Skip file header
	bool header = true;
	while (header) {
		std::getline (inF, line);
		if (line.compare(0,1,"#",0,1) != 0) {
			header = false;
		}
	}
	
	// Read line by line, parsing based on token
	C_[0] = atof(line.c_str());
	long long int index = 1;
	while (inF >> C_[index]) {
		index++;
	}
#endif
}

/*!
 * Read the macrostate distribution (biasing function) from a file.  This assumes the user has already guaranteed that the bounds are consistent,
 * e.g. Nmin and Nmax, as it will not check this automatically.  Also assumes file was generated by this code.  "Hand made"
 * ones might have formatting issues since parsing is done based on tokens.
 * 
 * \param [in] fileName Name of file containing lnPI.  Must include file extension.
 */
void tmmc::readlnPI (const std::string fileName) {
#ifdef NETCDF_CAPABLE
	NcFile dataFile (fileName.c_str(), NcFile::read);
	NcVar lnPI_data = dataFile.getVar("lnPI");
	if (lnPI_data.isNull()) throw customException("Macrostate distribution matrix (biasing function) was empty, cannot read");
	lnPI_data.getVar(&lnPI_[0]);
#else
	std::string line;
	std::ifstream inF (fileName);
	
	// Skip file header
	bool header = true;
	while (header) {
		std::getline (inF, line);
		if (line.compare(0,1,"#",0,1) != 0) {
			header = false;
		}
	}
	
	// Read line by line, parsing based on token
	lnPI_[0] = atof(line.c_str());
	long long int index = 1;
	while (inF >> lnPI_[index]) {
		index++;
	}
#endif
}

/*!
 * Wang-Landau biasing constructor.
 * 
 * \param [in] lnF Factor by which the estimate of the density of states in updated each time it is visited.
 * \param [in] g Factor by which lnF is reduced (multiplied) once "flatness" has been achieved.
 * \param [in] s Factor by which the min(H) must be within the mean of H to be considered "flat", e.g. 0.8 --> min is within 20% error of mean
 * \param [in] nSpec Number of species in the simulation.
 * \param [in] Nmax Vector of upper bound for number of particles of each species.
 * \param [in] Nmin Vector of lower bound for number of particles of each species. 
 */
wala::wala (const double lnF, const double g, const double s, const int nSpec, const std::vector <int> &Nmax, const std::vector <int> &Nmin) {
	if (lnF < 0) {
		throw customException ("lnF in Wang-Landau cannot be < 0");
	}
	lnF_ = lnF;
	
	if (g <= 0 || g >= 1) {
		throw customException ("In Wang-Landau 0 < g < 1");
	}
	g_ = g;
	
	if (s <= 0 || s >= 1) {
		throw customException ("In Wang-Landau 0 < s < 1");
	}
	s_ = s;
	
	if (nSpec > 0) {
		nSpec_ = nSpec;
	} else {
		throw customException ("Number of species in simulation must be > 0, cannot initialize wala");
	}
		
	if (Nmax.size() != Nmin.size()) {
		throw customException ("Nmin and Nmax vectors have unequal sizes, cannot initialize wala");
	}
		
	if (Nmax.size() != nSpec_) {
		throw customException ("Nmin and Nmax must specify bounds for all species, cannot initialize wala");
	}
		
	W_.resize(nSpec_, 0);
		
	__BIAS_INT_TYPE__ size = 1;
	for (unsigned int i = 0; i < nSpec_; ++i) {
		if (Nmin[i] > Nmax[i]) {
			throw customException ("Nmin > Nmax for species "+sstr(i+1));
		}
		W_[i] = (Nmax[i] - Nmin[i] +1);
		size *= W_[i];
	}
		
	Nmin_ = Nmin;
	Nmax_ = Nmax;
	
	// attempt to allocate memory for macrostate distribution matrix and initializes it all to 0
	try {
		lnPI_.resize(size, 0.0);
	} catch (const std::bad_alloc &ce) {
		throw customException ("Out of memory, cannot allocate space for macrostate distribution matrix in wala");
	}
	
	// initialize the visited-states histogram
	try {
		H_.resize(size, 0.0);
	} catch (const std::bad_alloc &ce) {
		throw customException ("Out of memory, cannot allocate space for visited-states histogram in wala");
	}
}

/*!
 * For multidimensional Wang-Landau biasing, get the 1D coordinate of the bias for multidimensional data.
 * 
 * \param [in] Nval Vector of the number of each species (in order) 
 */
const __BIAS_INT_TYPE__ wala::getAddress (const std::vector <int> &Nval) {
	__BIAS_INT_TYPE__ x = (Nval[nSpec_-1] - Nmin_[nSpec_-1]);
	for (int i = nSpec_-2; i >= 0; --i) {
		x *= W_[i];
		x += (Nval[i] - Nmin_[i]);
	}
	return x;
}

/*!
 * Update the estimate of the macrostate distribution.
 * 
 * \param [in] Nval Vector of current number of particles in the system
 */
void wala::update (const std::vector <int> &Nval) {
	__BIAS_INT_TYPE__ address = getAddress (Nval);
	lnPI_[address] += lnF_;
	H_[address] += 1.0;
}

/*!
 * Evaluate if the visited states histogram is approxiamtely "flat"
 * 
 * \return Returns whether the histogram is flat or not.
 */
bool wala::evaluateFlatness () {
	double min = H_[0], lnMean = -DBL_MAX;
	for (unsigned int i = 0; i < H_.size(); ++i) {
		if (H_[i] < min) {
			min = H_[i];
		}
		
		// summing so many doubles may overrun DBL_MAX, so instead track the lnMean
		lnMean = specExp(lnMean, log(H_[i]));
	}
	lnMean -= log(H_.size());
	
	if (log(min) - lnMean > log(s_)) {
		return true;
	}
	return false;
}

/*!
 * This should only be called when the "flatness" criterion is met. This then resets the visited-states histogram H, and decrements lnF.
 */
void wala::iterateForward () {
	lnF_ = lnF_*g_;
	std::fill(H_.begin(), H_.end(), 0);
}

/*!
 * Print the UN-NORMALIZED biasing function (lnPI) and possible the visted-states histogram to files. 
 * Will overwrite the files if another with that name exists. 
 * Prints in netCDF format if enabled.
 * 
 * \param [in] fileName Name of the file to print to.  Will append with "_lnPI" and "_H" for the macrostate distribution and visited-states histogram, respectively.
 * \param [in] printH Defaults to false, but if true will also print the visited-states histogram.
 */
void wala::print (const std::string fileName, bool printH) {
#ifdef NETCDF_CAPABLE
	// Print visited-states histogram
	if (printH) {
		const std::string name = fileName + "_H.nc"
		NcFile outFile(name.c_str(), NcFile::replace);
		NcDim probDim = outFile.addDim("vectorized_position", H_.size());
		NcVar probVar = outFile.addVar("H", ncDouble, probDim);
		const std::string dummyName = "number_species:";
		probVar.putAtt(dummyName.c_str(), sstr(nSpec_).c_str());
		for (unsigned int i = 0; i < nSpec_; ++i) {
	  	    const std::string attName = "species_"+sstr(i+1)+"_upper_bound:";
	    	probVar.putAtt(attName.c_str(), sstr(Nmax_[i]).c_str());
		}
		for (unsigned int i = 0; i < nSpec_; ++i) {
		    const std::string attName = "species_"+sstr(i+1)+"_lower_bound:";
	     	probVar.putAtt(attName.c_str(), sstr(Nmin_[i]).c_str());
		}
		probVar.putVar(&H_[0]);
	}
	
	// Print lnPI (bias) matrix
	const std::string name = fileName + "_lnPI.nc"
	NcFile outFile(name.c_str(), NcFile::replace);
	NcDim probDim = outFile.addDim("vectorized_position", lnPI_.size());
	NcVar probVar = outFile.addVar("lnPI", ncDouble, probDim);
	const std::string dummyName = "number_species:";
	probVar.putAtt(dummyName.c_str(), sstr(nSpec_).c_str());
	for (unsigned int i = 0; i < nSpec_; ++i) {
  	    const std::string attName = "species_"+sstr(i+1)+"_upper_bound:";
    	probVar.putAtt(attName.c_str(), sstr(Nmax_[i]).c_str());
	}
	for (unsigned int i = 0; i < nSpec_; ++i) {
	    const std::string attName = "species_"+sstr(i+1)+"_lower_bound:";
     	probVar.putAtt(attName.c_str(), sstr(Nmin_[i]).c_str());
	}
	probVar.putVar(&lnPI_[0]);
#else
	// Print visited-states histogram
	if (printH) {
		std::ofstream of;
		of.open(fileName+"_H.dat", 'w');
		of << "# Visited-states histogram in single row (vectorized) notation." << std::endl;
		of << "# Number of species:" << nSpec_ << std::endl;
		for (unsigned int i = 0; i < nSpec_; ++i) {
			of << "# species_"+sstr(i+1)+"_upper_bound: " << Nmax_[i] << std::endl;
		}
		for (unsigned int i = 0; i < nSpec_; ++i) {
			of << "# species_"+sstr(i+1)+"_lower_bound: " << Nmin_[i] << std::endl;
		}
		for (long long int i = 0; i < H_.size(); ++i) {
			of << H_[i] << std::endl;
		}
		of.close();
	}
	
	// Print lnPI (bias) matrix
	std::ofstream of;
	of.open(fileName+"_lnPI.dat", 'w');
	of << "# lnPI (bias) matrix in single row (vectorized) notation." << std::endl;
	of << "# Number of species:" << nSpec_ << std::endl;
	for (unsigned int i = 0; i < nSpec_; ++i) {
		of << "# species_"+sstr(i+1)+"_upper_bound: " << Nmax_[i] << std::endl;
	}
	for (unsigned int i = 0; i < nSpec_; ++i) {
		of << "# species_"+sstr(i+1)+"_lower_bound: " << Nmin_[i] << std::endl;
	}
	for (long long int i = 0; i < lnPI_.size(); ++i) {
		of << lnPI_[i] << std::endl;
	}
	of.close();
#endif
}

/*!
 * Read the macrostate distribution (biasing function) from a file.  This assumes the user has already guaranteed that the bounds are consistent,
 * e.g. Nmin and Nmax, as it will not check this automatically.  Also assumes file was generated by this code.  "Hand made"
 * ones might have formatting issues since parsing is done based on tokens.
 * 
 * \param [in] fileName Name of file containing lnPI.  Must include file extension.
 */
void wala::readlnPI (const std::string fileName) {
#ifdef NETCDF_CAPABLE
	NcFile dataFile (fileName.c_str(), NcFile::read);
	NcVar lnPI_data = dataFile.getVar("lnPI");
	if (lnPI_data.isNull()) throw customException("Macrostate distribution matrix (biasing function) was empty, cannot read");
	lnPI_data.getVar(&lnPI_[0]);
#else
	std::string line;
	std::ifstream inF (fileName);
	
	// Skip file header
	bool header = true;
	while (header) {
		std::getline (inF, line);
		if (line.compare(0,1,"#",0,1) != 0) {
			header = false;
		}
	}
	
	// Read line by line, parsing based on token
	lnPI_[0] = atof(line.c_str());
	long long int index = 1;
	while (inF >> lnPI_[index]) {
		index++;
	}
#endif
}