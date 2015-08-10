#include "bias.h"

/*!
 * Initialize the tmmc object.
 * 
 * \param [in] nSpec Number of species in the simulation.
 * \param [in] Nmax Upper bound for total number of particles.
 * \param [in] Nmin Lower bound for total number of particles.
 */
//tmmc::tmmc (const int nSpec, const int Nmax, const int Nmin) {
tmmc::tmmc (const int Nmax, const int Nmin) {	
	if (Nmin > Nmax) {
		throw customException ("Nmin > Nmax in TMMC bias");
	}
	__BIAS_INT_TYPE__ size = (Nmax - Nmin + 1);
	Nmin_ = Nmin;
	Nmax_ = Nmax;
	
	// attempt to allocate memory for collection matrix and initializes it all to 0
	try {
		C_.resize(3*size, 0);
	} catch (const std::bad_alloc &ce) {
		throw customException ("Out of memory, cannot allocate space for collection matrix in tmmc");
	}

	// attempt to allocate memory for probability matrix and initializes it all to 0
	try {
		P_.resize(3*size, 0);
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
 * Check if the the collection matrix has been filled (i.e. contains no zeros, except when crossing imposed bounds).
 * Technically, if in an ideal gas state, dU ~ 0 so pacc = 1.0, therefore probability of remaining in same state is 0.
 * However, this matrix should be considered filled, if transitions to N+1 and N-1 for each N are found (except at bounds).
 */
bool tmmc::checkFullyVisited () {
	for (__BIAS_INT_TYPE__ i = 0; i < C_.size(); i += 3) {
		if (i == 0) {
			// lower bound, so only +1 move must be sampled
			if (!(C_[i+1] > 0)) {
				return false;
			}			
		} else if (i == C_.size()-3) {
			// upper bound, so only -1 move must be sampled
			if (!(C_[i+2] > 0)) {
				return false;
			}
		} else {
			// midpoints, both +1 and -1 moves must be sampled
			if (!(C_[i+1] > 0) || !(C_[i+2] > 0)) {
				return false;
			}
		}
	}
	return true;
}

/*!
 * For a given multidimensional array which has been cast into 1D, find the address that refers to a given transition.
 * 
 * \param [in] Nstart Number of total species initially (before MC move)
 * \param [in] Nend Number of total species (in order) after the MC move
 */
const __BIAS_INT_TYPE__ tmmc::getTransitionAddress (const int Nstart, const int Nend) {
	// Layout of y = [0, +1, -1] 
	int addOrSubtract = (Nend - Nstart), y = 0;
	if (addOrSubtract == 0) {
		y = 0;
	} else if (addOrSubtract == 1) {
		y = 1;
	} else if (addOrSubtract == -1) {
		y = 2;
	} else {
		throw customException ("Illegal addOrSubtract value");
	}
	__BIAS_INT_TYPE__ x = Nstart - Nmin_;
	return x*3 + y;
}

/*!
 * Get the address in lnPI that corresponds to a given macrostate.
 * 
 * \param [in] Nval Number of total atoms
 */
const __BIAS_INT_TYPE__ tmmc::getAddress (const int Nval) {
	__BIAS_INT_TYPE__ x = Nval - Nmin_;
	return x;
}

/*!
 * Update the collection matrix.
 * 
 * \param [in] Nstart Total number of atoms initially (before MC move)
 * \param [in] Nend Total number of atoms after the MC move
 * \param [in] pa Unbiased Metropolis criterion for making a MC move (i.e. pa = min(1, exp(...)))
 */
void tmmc::updateC (const int Nstart, const int Nend, const double pa) {
	const int i = getTransitionAddress(Nstart, Nend), j = getTransitionAddress(Nstart, Nstart);
	C_[i] += pa;
	C_[j] += (1-pa);
}

/*!
 * Calculate the (natural logarithm of the) macrostate density matrix via the probability matrix.
 */
void tmmc::calculatePI () {
	for (__BIAS_INT_TYPE__ i = 0; i < C_.size(); i += 3) {
		double sum = 0.0;
		for (unsigned int j = 0; j < 3; ++j) {
			sum += C_[i+j];
		}
		if (sum > 0) {
			for (unsigned int j = 0; j < 3; ++j) {
				P_[i+j] = C_[i+j] / sum;
			}
		} else {
			// This state has not been visited at all if sum = 0.  However, at high system densities this could be the "correct"
			// result so this error may need to be discarded later in favor of:
			// P_[i+j] = 0;
			// However, having this throw an exception is also a good way to find that upper bound where the system is completely packed
			// so I will keep it this way for now.
			throw customException ("Cannot compute TMMC macrostate distribution because probability matrix contains zeros");
		}
	}
	
	// Reset first value to zero just to start fresh. Since only ratios matter this is perfectly fair.
	lnPI_[0] = 0.0;
	__BIAS_INT_TYPE__ address1, address2;
	for (__BIAS_INT_TYPE__ i = 0; i < lnPI_.size()-1; ++i) {
		address1 = getTransitionAddress(Nmin_+i, Nmin_+i+1);
		address2 = getTransitionAddress(Nmin_+i+1, Nmin_+i);
		if (!(P_[address1] > 0) || !(P_[address2] > 0)) {
			throw customException ("Cannot compute TMMC macrostate distribution because probability matrix contains zeros at address: P["+sstr(address1)+"] = "+sstr(P_[address1])+", P["+sstr(address2)+"] = "+sstr(P_[address2]));
		}
		lnPI_[i+1] = lnPI_[i] + log(P_[address1]/P_[address2]); // this is why P_ cannot be zero
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
		const std::string attName = "species_total_upper_bound:";
		probVar.putAtt(attName.c_str(), sstr(Nmax_).c_str());
		const std::string attName = "species_total_lower_bound:";
		probVar.putAtt(attName.c_str(), sstr(Nmin_).c_str());
		probVar.putVar(&C_[0]);
	}
	
	// Print lnPI (bias) matrix
	const std::string name = fileName + "_lnPI.nc"
	NcFile outFile(name.c_str(), NcFile::replace);
	NcDim probDim = outFile.addDim("vectorized_position", lnPI_.size());
	NcVar probVar = outFile.addVar("lnPI", ncDouble, probDim);
	const std::string attName = "species_total_upper_bound:";
	probVar.putAtt(attName.c_str(), sstr(Nmax_).c_str());
	const std::string attName = "species_total_lower_bound:";
	probVar.putAtt(attName.c_str(), sstr(Nmin_).c_str());
	probVar.putVar(&lnPI_[0]);
#else
	// Print collection matrix
	if (printC) {
		std::ofstream of;
		of.open(fileName+"_C.dat", std::ofstream::out);
		of << "# Collection matrix in single row (vectorized) notation." << std::endl;
		of << "# species_total_upper_bound: " << Nmax_ << std::endl;
		of << "# species_total_lower_bound: " << Nmin_ << std::endl;
		for (long long int i = 0; i < C_.size(); ++i) {
			of << C_[i] << std::endl;
		}
		of.close();
	}
	
	// Print lnPI (bias) matrix
	std::ofstream of;
	of.open(fileName+"_lnPI.dat", std::ofstream::out);
	of << "# lnPI (bias) matrix in single row (vectorized) notation." << std::endl;
	of << "# species_total_upper_bound: " << Nmax_ << std::endl;
	of << "# species_total_lower_bound: " << Nmin_ << std::endl;
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
	__BIAS_INT_TYPE__ index = 1;
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
 * \param [in] Nmax Upper bound for total number of particles.
 * \param [in] Nmin Lower bound for total number of particles. 
 */
//wala::wala (const double lnF, const double g, const double s, const int nSpec, const int Nmax, const int Nmin) {
wala::wala (const double lnF, const double g, const double s, const int Nmax, const int Nmin) {
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
		
	if (Nmin > Nmax) {
		throw customException ("Nmin > Nmax in Wang-Landau object");
	}
	
	__BIAS_INT_TYPE__ size = (Nmax - Nmin + 1);
	
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
 * For multidimensional Wang-Landau biasing, get the 1D coordinate of the macrostate distribution estimate (bias) for multidimensional data.
 * 
 * \param [in] Nval Total number of atoms in the system
 */
const __BIAS_INT_TYPE__ wala::getAddress (const int Nval) {
	__BIAS_INT_TYPE__ x = Nval - Nmin_;
	return x;
}

/*!
 * Update the estimate of the macrostate distribution.
 * 
 * \param [in] Nval Total current number of atoms in the system
 */
void wala::update (const int Nval) {
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
		const std::string attName = "species_total_upper_bound:";
		probVar.putAtt(attName.c_str(), sstr(Nmax_).c_str());
		const std::string attName = "species_upper_lower_bound:";
		probVar.putAtt(attName.c_str(), sstr(Nmin_).c_str());
		probVar.putVar(&H_[0]);
	}
	
	// Print lnPI (bias) matrix
	const std::string name = fileName + "_lnPI.nc"
	NcFile outFile(name.c_str(), NcFile::replace);
	NcDim probDim = outFile.addDim("vectorized_position", lnPI_.size());
	NcVar probVar = outFile.addVar("lnPI", ncDouble, probDim);
	const std::string dummyName = "number_species:";
	probVar.putAtt(dummyName.c_str(), sstr(nSpec_).c_str());
	const std::string attName = "species_total_upper_bound:";
	probVar.putAtt(attName.c_str(), sstr(Nmax_).c_str());
	const std::string attName = "species_total_lower_bound:";
	probVar.putAtt(attName.c_str(), sstr(Nmin_).c_str());
	probVar.putVar(&lnPI_[0]);
#else
	// Print visited-states histogram
	if (printH) {
		std::ofstream of;
		of.open(fileName+"_H.dat", std::ofstream::out);
		of << "# Visited-states histogram in single row (vectorized) notation." << std::endl;
		of << "# species_total_upper_bound:" << Nmax_ << std::endl;
		of << "# species_total_lower_bound:" << Nmin_ << std::endl;
		for (long long int i = 0; i < H_.size(); ++i) {
			of << H_[i] << std::endl;
		}
		of.close();
	}
	
	// Print lnPI (bias) matrix
	std::ofstream of;
	of.open(fileName+"_lnPI.dat", std::ofstream::out);
	of << "# lnPI (bias) matrix in single row (vectorized) notation." << std::endl;
	of << "# species_total_upper_bound:" << Nmax_ << std::endl;
	of << "# species_total_lower_bound:" << Nmin_ << std::endl;
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
	__BIAS_INT_TYPE__ index = 1;
	while (inF >> lnPI_[index]) {
		index++;
	}
#endif
}