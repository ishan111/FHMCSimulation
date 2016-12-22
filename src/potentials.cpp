#include "potentials.h"

/*!
 * Function for saving arbitrary potential into ASCII file.
 *
 * \param [in] filename Name of ASCII file to save potential to
 * \param [in] start r value to start printing U(r) from
 * \param [in] dr Increment to move in r between prints
 */
void pairPotential::savePotential(std::string filename, double start, double dr) {
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
 * Set the parameters in the Force-Shift Lennard-Jones equation.
 *
 * \param [in] params Vector of inputs: {epsilon, sigma, r_cut, Mtot}
 */
void fsLennardJones::setParameters (const std::vector < double > params) {
	if (params.size() != 4) {
		throw customException ("For fsLennardJones must specify 5 parameters: epsilon, sigma, r_cut, Mtot");
	} else {
		if (params[0] < 0) {
			throw customException ("For fsLennardJones, epsilon > 0");
		}
		if (params[1] < 0) {
			throw customException ("For fsLennardJones, sigma > 0");
		}
		if (params[2] < 0) {
			throw customException ("For fsLennardJones, r_cut > 0");
		}
		if (int(params[3]) < 1) {
			throw customException ("For fsLennardJones, total expanded ensemble states, Mtot >= 1");
		}

		paramsAreSet_ = true;
		params_ = params;

		useTailCorrection = false;

		// use a "constant volume" scheme to distribute the stages
		sigmaM_.resize(int(params[3]), 0);
		for (int i = 0; i < sigmaM_.size(); ++i) {
			if (i == 0) {
				// fully inserted
				sigmaM_[i] = params[1];
			} else {
				// use volume scaling so each stage is separated from its neighbors by the same dV
				double lastSigma = 0;
				if (i == 1) {
					lastSigma = 0;
				} else {
					lastSigma = sigmaM_[i-1];
				}
				sigmaM_[i] = pow(params[1]*params[1]*params[1]/(8.0*int(params[3])) + lastSigma*lastSigma*lastSigma, 1./3.);
			}
		}

		// scale energy linearly across the stages
		epsM_.resize(int(params[3]), 0);
		for (int i = 0; i < epsM_.size(); ++i) {
			if (i == 0) {
				// fully inserted
				epsM_[i] = params[0];
			} else {
				epsM_[i] = i*(params[0]/int(params[3]));
			}
		}
	}
}

/*!
 * Return the energy of two particles.
 * \f[ U(r) = U_{LJ}(r) - \frac{\del U_{LJ}(r_{cut})}{\del r}(r - r_{cut}) - U_{LJ}(r_{cut}) \quad r < r_{cut}
	\f]
 *
 * \param [in] a1 Atom 1
 * \param [in] a2 Atom 2
 * \param [in] box Simulation box dimensions
 *
 * \return U(r)
 */
double fsLennardJones::energy (const atom* a1, const atom* a2, const std::vector < double > &box) {
	if (!paramsAreSet_) {
		throw customException ("For fsLennardJones parameters not set");
	}

	const double r_sq = pbcDist2(a1->pos, a2->pos, box);

	// only one of these atoms (at most) should be "partially" inserted
	int mState = 0;
	if (a1->mState != 0) {
		mState = a1->mState;
	}
	if (a2->mState != 0) {
		mState = a2->mState;
	}

	const double r2 = (sigmaM_[mState]*sigmaM_[mState]/r_sq), r6 = r2*r2*r2, r12 = r6*r6;
	const double r2c = (sigmaM_[mState]*sigmaM_[mState]/(params_[2]*params_[2])), r6c = r2c*r2c*r2c, r12c = r6c*r6c;
	if (r_sq < params_[2]*params_[2]) {
		return 4.0*epsM_[mState]*(r12 - r6) - 4.0*epsM_[mState]*(r12c - r6c) - (sqrt(r_sq) - params_[2])*(-24.0*epsM_[mState]/params_[2]*(2.0*r12c - r6c));
	} else {
		return 0.0;
	}
}

/*!
 * Tail correction for Force-Shifted Lennard-Jones is, by definition, zero always.
 *
 * \param [in] rhoBath Number density of the surrounding fluid
 *
 * \return U_tail
 */
double fsLennardJones::tailCorrection(const double rhoBath) {
	return 0.0;
}

/*!
 * Return the value of r_{cut} for this potential.
 *
 * \return rcut
 */
double fsLennardJones::rcut () {
	if (!paramsAreSet_) {
		throw customException ("For fsLennardJones parameters not set");
	} else {
		return params_[2];
	}
}

/*!
 * Set the parameters in the Lennard-Jones equation.
 *
 * \param [in] params Vector of inputs: {epsilon, sigma, r_cut, u_shift, Mtot}
 */
void lennardJones::setParameters (const std::vector < double > params) {
	if (params.size() != 5) {
		throw customException ("For lennardJones must specify 5 parameters: epsilon, sigma, r_cut, u_shift, Mtot");
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
		if (int(params[4]) < 1) {
			throw customException ("For lennardJones, total expanded ensemble states, Mtot >= 1");
		}

		paramsAreSet_ = true;
		params_ = params;

		useTailCorrection = true;

		// use a "constant volume" scheme to distribute the stages
		sigmaM_.resize(int(params[4]), 0);
		for (int i = 0; i < sigmaM_.size(); ++i) {
			if (i == 0) {
				// fully inserted
				sigmaM_[i] = params[1];
			} else {
				// use volume scaling so each stage is separated from its neighbors by the same dV
				double lastSigma = 0;
				if (i == 1) {
					lastSigma = 0;
				} else {
					lastSigma = sigmaM_[i-1];
				}
				sigmaM_[i] = pow(params[1]*params[1]*params[1]/(8.0*int(params[4])) + lastSigma*lastSigma*lastSigma, 1./3.);
			}
		}

		// scale energy linearly across the stages
		epsM_.resize(int(params[4]), 0);
		for (int i = 0; i < epsM_.size(); ++i) {
			if (i == 0) {
				// fully inserted
				epsM_[i] = params[0];
			} else {
				epsM_[i] = i*(params[0]/int(params[4]));
			}
		}

		// scale energy linearly across the stages
		uShiftM_.resize(int(params[4]), 0);
		for (int i = 0; i < epsM_.size(); ++i) {
			if (i == 0) {
				// fully inserted
				uShiftM_[i] = params[3];
			} else {
				uShiftM_[i] = i*(params[3]/int(params[4]));
			}
		}
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

	const double r_sq = pbcDist2(a1->pos, a2->pos, box);

	// only one of these atoms (at most) should be "partially" inserted
	int mState = 0;
	if (a1->mState != 0) {
		mState = a1->mState;
	}
	if (a2->mState != 0) {
		mState = a2->mState;
	}

	double r2 = (sigmaM_[mState]*sigmaM_[mState]/r_sq), r6 = r2*r2*r2, r12 = r6*r6;
	if (r_sq < params_[2]*params_[2]) {
		return 4.0*epsM_[mState]*(r12 - r6) + uShiftM_[mState];
	} else {
		return 0.0;
	}
}

/*!
 * Calculate the tail correction with the approximation g(r) = 1 for r_{cut} > 1
 * as explained in Frenkel & Smit in eq. (3.2.5).  Tail corrections only account for number of fully inserted particles
 * so I have chosen not to scale this part of the energy with expanded ensemble stage.
 *
 * \param [in] rhoBath Number density of the surrounding fluid
 *
 * \return U_tail
 */
double lennardJones::tailCorrection(const double rhoBath) {
	if (rhoBath < 0) {
		return 0;
	}
	const double r3 = (params_[1]*params_[1]*params_[1])/(params_[2]*params_[2]*params_[2]);
	const double r9 = r3*r3*r3;

	return 2.0*(8.0/3.0*PI*rhoBath*params_[0]*params_[1]*params_[1]*params_[1]*(r9/3.0 - r3));
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
 * \param [in] params Vector of inputs: {r_cut, r_shift, u_shift, u_infinity, Mtot}
 */
void tabulated::setParameters (const std::vector < double > params) {
	if (params.size() != 5) {
		throw customException ("For tabulated must specify 5 parameters: r_cut, r_shift, u_shift, u_infinity, Mtot");
	} else {
		if (params[0] < 0) {
			throw customException ("For tabulated, r_cut > 0");
		}
		if (int(params[4]) < 1) {
			throw customException ("For tabulated, total expanded ensemble states, Mtot >= 1");
		}

		paramsAreSet_ = true;
		params_ = params;

		useTailCorrection = false;

		// scale energy by a constant factor
		mScale.resize(int(params[4]), 0);
		for (int i = 0; i < mScale.size(); ++i) {
			if (i == 0) {
				mScale[i] = 1.0;
			} else {
				mScale[i] = 1.0/int(params[4])*i;
			}
		}
	}
}

/*!
 * Load a tabulated potential from an external file and store it internally
 *
 * \param [in] filename Name of ASCII file to read (r, U(r)) from
 */
void tabulated::loadPotential(std::string filename) {
	std::cout << "Loading pair potential from file: " << filename<<std::endl;
	// first check, if file exists
	if (fileExists(filename)) {
		std::cout << "File found, processing." << std::endl;
		table.clear();
		double r, pot;
		unsigned int lineCounter = 0;

		std::ifstream inputData(filename.c_str());
		while (inputData >> r >> pot) {
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
			std::cerr << "Tabulated potential " << filename << " needs at least 2 entries, cannot setup potential." << std::endl;
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
	} else {
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

	const double r = sqrt(pbcDist2(a1->pos, a2->pos, box));
	double en = 0.0;

	// only one of these atoms (at most) should be "partially" inserted
	int mState = 0;
	if (a1->mState != 0) {
		mState = a1->mState;
	}
	if (a2->mState != 0) {
		mState = a2->mState;
	}

	if (r < params_[1]) {
		std::cerr << "Distance r too small in energy calculation in tabulated potential. Returning value at r = " << start << std::endl;
		en = table[0];
	} else if (r > params_[0]) {
		en = params_[3];
	} else {
		const unsigned int lowerIndex = floor((r-params_[1])/dr);
		const unsigned int upperIndex = ceil((r-params_[1])/dr);
		const double upperFraction = (r-params_[1])/dr-lowerIndex;
		const double lowerFraction = 1.0-upperFraction;
		en = (lowerFraction*table[lowerIndex] + upperFraction*table[upperIndex] + params_[2])*mScale[mState];
	}

	return en;
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
 * \param [in] params Vector of inputs: {sigma, wellwidth, welldepth (magnitude), Mtot}
 */
void squareWell::setParameters (const std::vector < double > params) {
	if (params.size() != 4) {
		throw customException ("For squareWell must specify 4 parameters: sigma, wellwidth, welldepth, Mtot");
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
		if (int(params[3]) < 1) {
			throw customException ("For squareWell, total expanded ensemble states, Mtot >= 1");
		}

		useTailCorrection = false;

		// use a "constant volume" scheme to distribute the stages
		sigmaM_.resize(int(params[3]), 0);
		rangeM_.resize(int(params[3]), 0);
		for (int i = 0; i < sigmaM_.size(); ++i) {
			if (i == 0) {
				// fully inserted
				sigmaM_[i] = params[0];
				rangeM_[i] = params[0] + params[1];
			} else {
				// use volume scaling so each stage is separated from its neighbors by the same dV
				double lastSigma = 0;
				if (i == 1) {
					lastSigma = 0;
				} else {
					lastSigma = sigmaM_[i-1];
				}
				sigmaM_[i] = pow(params[0]*params[0]*params[0]/(8.0*int(params[3])) + lastSigma*lastSigma*lastSigma, 1./3.);
				rangeM_[i] = sigmaM_[i] + params[1];
			}
		}

		// scale energy linearly across the stages
		epsM_.resize(int(params[3]), 0);
		for (int i = 0; i < epsM_.size(); ++i) {
			if (i == 0) {
				// fully inserted
				epsM_[i] = -params[2];
			} else {
				epsM_[i] = -i*(params[2]/int(params[3]));
			}
		}

		// save parameters as sigma, (sigma+wellWidth), -wellDepth to speed up energy calculation
		params_ = params;
		params_[1] = params[0] + params[1]; // max rcut
		params_[2] = -params[2];
		paramsAreSet_ = true;
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

	const double r = sqrt(pbcDist2(a1->pos, a2->pos, box));

	int mState = 0;
	if (a1->mState != 0) {
		mState = a1->mState;
	}
	if (a2->mState != 0) {
		mState = a2->mState;
	}

	if (r < sigmaM_[mState]) {
		return NUM_INFINITY;
	} else if (r < rangeM_[mState]) {
		return epsM_[mState];
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
 * \param [in] params Vector of inputs: {sigma, Mtot}
 */
void hardCore::setParameters (const std::vector < double > params) {
	if (params.size() != 2) {
		throw customException ("For hardCore must specify 2 parameter: {sigma, Mtot}");
	} else {
		if (params[0] < 0) {
			throw customException ("For hardCore, sigma >= 0");
		}
		if (int(params[1]) < 1) {
			throw customException ("For hardCore, total expanded ensemble state, Mtot >= 1");
		}

		params_ = params;
		paramsAreSet_ = true;

		useTailCorrection = false;

		// use a "constant volume" scheme to distribute the stages
		sigmaM_.resize(int(params[1]), 0);
		for (int i = 0; i < sigmaM_.size(); ++i) {
			if (i == 0) {
				// fully inserted
				sigmaM_[i] = params[0];
			} else {
				// use volume scaling so each stage is separated from its neighbors by the same dV
				double lastSigma = 0;
				if (i == 1) {
					lastSigma = 0;
				} else {
					lastSigma = sigmaM_[i-1];
				}
				sigmaM_[i] = pow(params[0]*params[0]*params[0]/(8.0*int(params[1])) + lastSigma*lastSigma*lastSigma, 1./3.);
			}
		}
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

	int mState = 0;
	if (a1->mState != 0) {
		mState = a1->mState;
	}
	if (a2->mState != 0) {
		mState = a2->mState;
	}

	const double r = sqrt(pbcDist2(a1->pos, a2->pos, box));

	if (r < sigmaM_[mState]) {
		return NUM_INFINITY;
	} else {
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
		if (fabs(params_[0]) < 1.0e-12) {
			// in case sigma = 0 (used for ideal gas case) just return finite value for cell lists
			return 1.0;
		} else {
			return (params_[0]);
		}
	}
}
