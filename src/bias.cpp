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
	
	// attempt to allocate memory for lnPI matrix and initializes it all to 1
	try {
		lnPI_.resize(size, 1.0);
	} catch (const std::bad_alloc &ce) {
		throw customException ("Out of memory, cannot allocate space for macrostate distribution matrix in tmmc");
	}
}

/*!
 * Get the difference offset between two states of the simulation. 
 * This assumes ONLY one species has been modified and by a value of 1.  
 * Multiple changes, or those > 1, WILL result in unexpected behavior.
 * Default resets specId and addOrSubtract to 0 when called.
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
	
	// Reset first value to unity just to start fresh. Since only ratios matter this is perfectly fair.
	lnPI_[0] = 1.0;
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
	
	// attempt to allocate memory for macrostate distribution matrix and initializes it all to 1
	try {
		lnPI_.resize(size, 1.0);
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
 * This should only be called when the "flatness" criterion is met. This then reset the visited-states histogram H, and decrements lnF.
 */
void wala::iterateForward () {
	lnF_ = lnF_*g_;
	std::fill(H_.begin(), H_.end(), 0);
}