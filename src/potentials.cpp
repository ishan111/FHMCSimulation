#include "potentials.h"

/*!
 * Function for saving arbitrary potential into ASCII file.
 * 
 * \param [in] filename Name of ASCII file to save potential to 
 * \param [in] start r value to start printing U(r) from 
 * \param [in] dr Increment to move in r between prints
 */
void pairPotential::savePotential(std::string filename, double start, double dr)
{
	if (dr <= 0.0) {
		throw customException("The value for dr must be positive");
	}
	
	double r = start;
	std::ofstream outData(filename.c_str());
		
	atom a1, a2;
	std::vector <double> p1 (3, 0), box(3, 2.1*rcut()); // box must be atleast 2 rcut, +0.1 for good measure
	a1.pos = p1;
	a2.pos = p1;
	
	while (r < rcut()) {
		a2.pos[2] = r;
		outData << r << "\t" << energy(&a1, &a2, box) << std::endl;		
		r += dr;
	}
	
	outData.close();
}

/*!
 * Set the parameters in the Lennard-Jones equation.
 * 
 * \param [in] params Vector of inputs: {epsilon, sigma, r_cut, u_shift}
 */
void lennardJones::setParameters (const std::vector < double > params) {
	if (params.size() != 4) {
		throw customException ("For lennardJones must specify 4 parameters: epsilon, sigma, r_cut, u_shift");
	} else {
		if (params[0] < 0) {
			throw customException ("For lennardJones, epsilon > 0");
		}
		if (params[1] < 0) {
			throw customException ("For lennardJones, sigma > 0");
		}
		if (params[2] < 0) {
			throw customException ("For lennardJones, r_cut > 0");
		}
		
		paramsAreSet_ = true;
		params_ = params;
		
		useTailCorrection = true;
	}
}

/*!
 * Return the energy of two particles.
 * \f[ U(r) = 4 \epsilon \left( \left \frac{ \sigma }{ r } \right)^{12} - \left( \frac{ sigma }{ r } \right)^6 \right) + U_{shift} \quad r < r_{cut}
	\f]
 *
 * \param [in] a1 Atom 1
 * \param [in] a2 Atom 2
 * \param [in] box Simulation box dimensions
 * 
 * \return U(r)
 */
double lennardJones::energy (const atom* a1, const atom* a2, const std::vector < double > &box) {
	if (!paramsAreSet_) {
		throw customException ("For lennardJones parameters not set");
	}
	
	if (a1->mState != 1 || a2->mState != 1) {
		throw customException ("Have not implemented expanded ensemble for lennardJones potential yet");
	} 
	
	const double r = sqrt(pbc_dist2(a1->pos, a2->pos, box));
	
	double r1 = (params_[1]/r), r3 = r1*r1*r1, r6 = r3*r3, r12 = r6*r6;
	if (r < params_[2]) {
		return 4.0*params_[0]*(r12 - r6) + params_[3];
	} else {
		return 0.0;
	}
}

/*!
 * Calculate the tail correction with the approximation g(r) = 1 for r_{cut} > 1
 * as explained in Frenkel & Smit in eq. (3.2.5)
 * 
 * \param [in] rhoBath NUmber density of the surrounding fluid
 * 
 * \return U_tail
 */
double lennardJones::tailCorrection(const double rhoBath) {
	const double r3 = (params_[1]*params_[1]*params_[1])/(params_[2]*params_[2]*params_[2]);
	const double r9 = r3*r3*r3;
	
	return (8.0/3.0*PI*rhoBath*params_[0]*params_[1]*params_[1]*params_[1]*(r9/3.0 - r3));
}

/*!
 * Return the value of r_{cut} for this potential.
 * 
 * \return rcut
 */
double lennardJones::rcut () {
	if (!paramsAreSet_) {
		throw customException ("For lennardJones parameters not set");
	} else {
		return params_[2];
	}
}

/*!
 * Set the parameters in the tabulated potential
 * 
 * \param [in] params Vector of inputs: {r_cut, r_shift, u_shift, u_infinity}
 */
void tabulated::setParameters (const std::vector < double > params) {
	if (params.size() != 4) {
		throw customException ("For tabulated must specify 4 parameters: r_cut, r_shift, u_shift, u_infinity");
	} else {
		if (params[0] < 0) {
			throw customException ("For tabulated, r_cut > 0");
		}
		
		paramsAreSet_ = true;
		params_ = params;
		
		useTailCorrection = false;
	}
}

/*!
 * Load a tabulated potential from an external file and store it internally
 * 
 * \param [in] filename Name of ASCII file to read (r, U(r)) from 
 */
void tabulated::loadPotential(std::string filename)
{
	std::cout<<"Loading pair potential from file: "<<filename<<std::endl;
	// first check, if file exists
	if (fileExists(filename))
	{
		std::cout<<"File found, processing."<<std::endl;
		table.clear();
		double r, pot;
		unsigned int lineCounter = 0;
		
		std::ifstream inputData(filename.c_str());
		while (inputData>>r>>pot) {
			if (lineCounter == 0) {
				start = r;
			}
			else if (lineCounter == 1) {
				dr = r - start;
			}
			
			table.push_back(pot);
			lineCounter++;
		}
		inputData.close();
		
		if (lineCounter < 2) {
			paramsAreSet_ = false;
			std::cerr<<"Tabulated potential "<<filename<<" needs at least 2 entries, cannot setup potential."<<std::endl;
			return;
		}
		
		// if parameters are not set, set default parameters
		if (!paramsAreSet_) {
			params_.assign(4, 0.0);
			params_[0] = start + (table.size()-1)*dr;
			params_[1] = start;
			params_[2] = 0.0;
			params_[3] = 0.0;
			paramsAreSet_ = true;
		}
	}
	else {
		std::cerr<<"File "<<filename<<" not found, cannot setup potential."<<std::endl;
		paramsAreSet_ = false;
	}
}

/*!
 * Return the energy of two particles.
 * Use linear interpolation to calculate energy from tabulated values
 * 
 * \param [in] a1 Atom 1
 * \param [in] a2 Atom 2
 * \param [in] box Simulation box dimensions
 * 
 * \return U(r)
 */
double tabulated::energy (const atom* a1, const atom* a2, const std::vector < double > &box) {
	if (!paramsAreSet_) {
		throw customException ("For tabulated parameters not set");
	}
	
	if (a1->mState != 1 || a2->mState != 1) {
		throw customException ("Have not implemented expanded ensemble for tabulated potential yet");
	} 
	
	const double r = sqrt(pbc_dist2(a1->pos, a2->pos, box));
	
	if (r < params_[1]) {
		std::cerr<<"distance r too small in energy calculation in tabulated potential. Returning value at r="<<start<<std::endl;
	} else if (r > params_[0]) {
		return params_[3];
	} else {
		const unsigned int lowerIndex = floor((r-params_[1])/dr);
		const unsigned int upperIndex = ceil((r-params_[1])/dr);
		
		const double upperFraction = (r-params_[1])/dr-lowerIndex;
		const double lowerFraction = 1.0-upperFraction;
		
		return (lowerFraction*table[lowerIndex] + upperFraction*table[upperIndex] + params_[2]);
	}
}

/*!
 * Tail correction for a tabulated potential always returns 0 since no information about what the potential is after its cutoff radius.
 * 
 * \param [in] Number density of the surrounding fluid
 * 
 * \return U_tail
 */
double tabulated::tailCorrection(const double rhoBath) {
	return 0.0;
}

/*!
 * Return the value of r_{cut} for this potential.
 * 
 * \return rcut
 */
double tabulated::rcut () {
	if (!paramsAreSet_) {
		throw customException ("For tabulated parameters not set");
	} else {
		return params_[0];
	}
}

/*!
 * Set the parameters in the square-well equation.
 * 
 * \param [in] params Vector of inputs: {sigma, wellwidth, welldepth (magnitude)}
 */
void squareWell::setParameters (const std::vector < double > params) {
	if (params.size() != 3) {
		throw customException ("For squareWell must specify 3 parameters: sigma, wellwidth, welldepth");
	} else {
		if (params[0] < 0) {
			throw customException ("For squareWell, sigma > 0");
		}
		if (params[1] < 0) {
			throw customException ("For squareWell, wellwidth > 0");
		}
		if (params[2] < 0) {
			throw customException ("For squareWell, welldepth (magnitude) > 0");
		}
		
		// save parameters as sigma, (sigma+wellWidth), -wellDepth to speed up energy calculation
		params_ = params;
		params_[1] = params[0] + params[1];
		params_[2] = -params[2];
		paramsAreSet_ = true;
		
		useTailCorrection = false;
	}
}

/*!
 * Return the energy of two particles.
 * 
 * \param [in] a1 Atom 1
 * \param [in] a2 Atom 2
 * \param [in] box Simulation box dimensions
 * 
 * \return U(r)
 */
double squareWell::energy (const atom* a1, const atom* a2, const std::vector < double > &box) {
	if (!paramsAreSet_) {
		throw customException ("For squareWell parameters not set");
	}
	
	if (a1->mState != 1 || a2->mState != 1) {
		throw customException ("Have not implemented expanded ensemble for squareWell potential yet");
	} 
	
	const double r = sqrt(pbc_dist2(a1->pos, a2->pos, box));
	
	if (r < params_[0]) {
		return NUM_INFINITY;
	} else if (r < params_[1]) {
		return params_[2];
	} else {
		return 0.0;
	}
}

/*!
 * Tail correction for a square well potential always returns 0.
 * 
 * \param [in] Number density of the surrounding fluid
 * 
 * \return U_tail
 */
double squareWell::tailCorrection(const double rhoBath) {
	return 0.0;
}

/*!
 * Return the value of r_{cut} for this potential.
 * 
 * \return r_cut
 */
double squareWell::rcut () {
	if (!paramsAreSet_) {
		throw customException ("For squareWell parameters not set");
	} else {
		return (params_[1]);
	}
}

/*!
 * Set the parameters in the hard-core potential.
 * 
 * \param [in] params Vector of inputs: {sigma, infinity}
 */
void hardCore::setParameters (const std::vector < double > params) {
	if (params.size() != 1) {
		throw customException ("For hardCore must specify 1 parameter: sigma");
	} else {
		if (params[0] < 0) {
			throw customException ("For hardCore, sigma > 0");
		}
		
		params_ = params;
		paramsAreSet_ = true;
		
		useTailCorrection = false;
	}
}

/*!
 * Return the energy of two particles.
 * 
 * \param [in] a1 Atom 1
 * \param [in] a2 Atom 2
 * \param [in] box Simulation box dimensions
 * 
 * \return U(r)
 */
double hardCore::energy (const atom* a1, const atom* a2, const std::vector < double > &box) {
	if (!paramsAreSet_) {
		throw customException ("For hardCore parameters not set");
	}
	
	if (a1->mState != 1 || a2->mState != 1) {
		throw customException ("Have not implemented expanded ensemble for hardCore potential yet");
	} 
	
	const double r = sqrt(pbc_dist2(a1->pos, a2->pos, box));
	
	if (r < params_[0]) {
		return NUM_INFINITY;
	}
	else {
		return 0.0;
	}
}

/*!
 * Tail correction for a hard core potential always returns 0 radius.
 * 
 * \param [in] Number density of the surrounding fluid
 * 
 * \return U_tail
 */
double hardCore::tailCorrection(const double rhoBath) {
	return 0.0;
}

/*!
 * Return the value of r_{cut} for this potential.
 * 
 * \return r_cut
 */
double hardCore::rcut () {
	if (!paramsAreSet_) {
		throw customException ("For hardCore parameters not set");
	} else {
		return (params_[0]);
	}
}
