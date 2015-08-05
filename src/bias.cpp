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
	
	int size = 1;
	for (unsigned int i = 0; i < Nmin.size(); ++i) {
		if (Nmin[i] > Nmax[i]) {
			throw customException ("Nmin > Nmax for species "+sstr(i+1));
		}
		W_[i] = (Nmax[i] - Nmin[i] +1);
		size *= W[i];
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
const int tmmc::getAddress (const std::vector <int> &Nstart, const std::vector <int> &Nend) {
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
	
	long long int x = (Nstart[nSpec_-1] - Nmin[nSpec_-1]);
	for (unsigned int i = nSpec_-2; i >= 0; --i) {
		x *= W_[i];
		x += (Nstart[i] - Nmin[i]);
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
const int tmmc::getAddress (const std::vector <int> &Nstart, const int specId, const int addOrSubtract) {		
	int y = 0;
	if (addOrSubtract != 0) {
		y += 2*(specId+1);
		if (addOrSubtract == 1) {
			y -= 1;
		} else {
			throw customException ("Bad addOrSubtract value: "+sstr(addOrSubtract));
		}
	}
	
	long long int x = (Nstart[nSpec_-1] - Nmin[nSpec_-1]);
	for (unsigned int i = nSpec_-2; i >= 0; --i) {
		x *= W_[i];
		x += (Nstart[i] - Nmin[i]);
	}
	
	return x*(1+2*nSpec_) + y;
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
	const int i = getAddress (Nstart, specId, addOrSubtract);
	if (addOrSubtract == 0) {
		C_[i] += (1-pa);
	} else {
		C_[i] += pa;
	}
}

/*!
 * Calculate the probability matrix.
 */
void tmmc::calculateP () {
	
}