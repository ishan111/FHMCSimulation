#include "potentials.h"

/*!
 * Function for saving arbitrary potential into textfile
*/
void pairPotential::savePotential(std::string filename, double start, double dr)
{
	if (dr <= 0.0) {
		throw customException("The value for dr must be positive");
	}
	
	double r = start;
	std::ofstream outData(filename.c_str());
		
	while (r < rcut()) {
		outData<<r<<"\t"<<energy(r)<<std::endl;		
		r += dr;
	}
	
	outData.close();
}

/*!
 * Set the parameters in the Lennard-Jones equation.
 * \param [in] params Vector of inputs: {epsilon, sigma, rcut, ushift}
 */
void lennardJones::setParameters (const std::vector < double > params) {
	if (params.size() != 4) {
		throw customException ("For lennardJones must specify 4 parameters: epsilon, sigma, rcut, ushift");
	} else {
		if (params[0] < 0) {
			throw customException ("For lennardJones, epsilon > 0");
		}
		if (params[1] < 0) {
			throw customException ("For lennardJones, sigma > 0");
		}
		if (params[2] < 0) {
			throw customException ("For lennardJones, rcut > 0");
		}
		
		paramsAreSet_ = true;
		params_ = params;
		
		useTailCorrection = true;
	}
}

/*!
 * Return the energy of two particles separate by a distance r.
 * \f[ U(r) = 4 \epsilon \left( \left \frac{ \sigma }{ r } \right)^{12} - \left( \frac{ sigma }{ r } \right)^6 \right) + U_{shift} \quad r < r_{cut}
	\f]
 * \param [in] r Scalar separation, needs to be the minimum image
 */
double lennardJones::energy (const double r) {
	if (!paramsAreSet_) {
		throw customException ("For lennardJones parameters not set");
	}
	double r1 = (params_[1]/r), r3 = r1*r1*r1, r6 = r3*r3, r12 = r6*r6;
	if (r < params_[2]) {
		return 4.0*params_[0]*(r12 - r6) + params_[3];
	} else {
		return 0.0;
	}
}

/*!
 * Calculate the tail correction with the approximation g(r) =1 for r_{cut} > 1
 * as explained in Frenkel&Smit in eq. (3.2.5)
*/
double lennardJones::tailCorrection(const double rhoBath) {
	const double r3 = (params_[1]*params_[1]*params_[1])/(params_[2]*params_[2]*params_[2]);
	const double r9 = r3*r3*r3;
	
	return (8.0/3.0*PI*rhoBath*params_[0]*params_[1]*params_[1]*params_[1]*(r9/3.0 - r3));
}

/*!
 * Return the value of r_{cut} for this potential.
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
 * \param [in] params Vector of inputs: {rcut, rshift, ushift, uinfinity}
 */
void tabulated::setParameters (const std::vector < double > params) {
	if (params.size() != 4) {
		throw customException ("For tabulated must specify 4 parameters: rcut, rshift, ushift, uinfinity");
	} else {
		if (params[0] < 0) {
			throw customException ("For tabulated, rcut > 0");
		}
		
		paramsAreSet_ = true;
		params_ = params;
		
		useTailCorrection = false;
	}
}

/*!
 * Load a tabulated potential from an external file and store it internally
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
 * Return the energy of two particles separate by a distance r.
 * Use linear interpolation to calculate energy from tabulated values
 * \param [in] r Scalar separation, needs to be the minimum image
 */
double tabulated::energy (const double r) {
	if (!paramsAreSet_) {
		throw customException ("For tabulated parameters not set");
	}
	
	if (r < params_[1]) {
		std::cerr<<"distance r too small in energy calculation in tabulated potential. Returning value at r="<<start<<std::endl;
	}
	else if (r > params_[0]) {
		return params_[3];
	}
	else {
		const unsigned int lowerIndex = floor((r-params_[1])/dr);
		const unsigned int upperIndex = ceil((r-params_[1])/dr);
		
		const double upperFraction = (r-params_[1])/dr-lowerIndex;
		const double lowerFraction = 1.0-upperFraction;
		
		return (lowerFraction*table[lowerIndex] + upperFraction*table[upperIndex] + params_[2]);
	}
}

double tabulated::tailCorrection(const double rhoBath) {
	return 0.0;
}

/*!
 * Return the value of r_{cut} for this potential.
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
 * \param [in] params Vector of inputs: {sigma, wellwidth, welldepth}
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
		
		// save parameters as sigma, (sigma+wellWidth), -wellDepth to speed up energy calculation
		params_ = params;
		params_[1] = params[0] + params[1];
		params_[2] = -params[2];
		paramsAreSet_ = true;
		
		useTailCorrection = false;
	}
}

/*!
 * Return the energy of two particles separate by a distance r.
 * \param [in] r Scalar separation, needs to be the minimum image
 */
double squareWell::energy (const double r) {
	if (!paramsAreSet_) {
		throw customException ("For squareWell parameters not set");
	}
	if (r < params_[0]) {
		return 1.0E12;
	} else if (r < params_[1]) {
		return params_[2];
	}
	else {
		return 0.0;
	}
}

double squareWell::tailCorrection(const double rhoBath) {
	return 0.0;
}

/*!
 * Return the value of r_{cut} for this potential.
 * \return rcut
 */
double squareWell::rcut () {
	if (!paramsAreSet_) {
		throw customException ("For squareWell parameters not set");
	} else {
		return (params_[1]);
	}
}


/*!
 * Set the parameters in the hard-core potential 
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
 * Return the energy of two particles separate by a distance r.
 * \param [in] r Scalar separation, needs to be the minimum image
 */
double hardCore::energy (const double r) {
	if (!paramsAreSet_) {
		throw customException ("For hardCore parameters not set");
	}
	if (r < params_[0]) {
		return 1.0E12;
	}
	else {
		return 0.0;
	}
}

double hardCore::tailCorrection(const double rhoBath) {
	return 0.0;
}

/*!
 * Return the value of r_{cut} for this potential.
 * \return rcut
 */
double hardCore::rcut () {
	if (!paramsAreSet_) {
		throw customException ("For hardCore parameters not set");
	} else {
		return (params_[0]);
	}
}
