#ifndef POTENTIALS_H_
#define POTENTIALS_H_

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include "global.h"
#include "utilities.h"
#include "atom.h"

#define NUM_INFINITY std::numeric_limits<double>::max()

/*!
 * Abstract base class which defines all pair potentials.
 * Note that tail corrections should be "double" the "standard tail corrections" which try to account for
 * double-counting ahead of time.  The system only iterates over unique pairs, so do NOT do this.
 * Also, when doing MC insertions/deletions this double counting is never a problem so easier to do it this way.
 */
class pairPotential {
public:
	pairPotential () { paramsAreSet_ = false; }
	virtual ~pairPotential () {;}

	bool useTailCorrection;
	virtual double energy (const atom* a1, const atom* a2, const std::vector < double > &box) = 0;
	virtual double tailCorrection (const double rhoBath) = 0;
	virtual void setParameters (const std::vector < double > params) = 0;
	virtual double rcut () = 0;	//!< All potentials should be able to return their r_{cut} values so neighbor lists, etc. can use them

	void savePotential(std::string filename, double start, double dr);
	std::vector < double > params_; //!< Parameters (constants) that are needed to calculate U(r)
	bool paramsAreSet_; //!< Logical check if the paramters for this potential have been specified by the user
};

/*!
 * Lennard-Jones Potential
 * Parameters should be specified in the following order: { epsilon, sigma, r_cut, u_shift, Mtot }
 */
class lennardJones : public pairPotential {
public:
	~lennardJones () {;}
	void setParameters (const std::vector < double > params);
	double energy (const atom* a1, const atom* a2, const std::vector < double > &box);
	double tailCorrection (const double rhoBath);
	double rcut ();

private:
	std::vector < double > epsM_; //!< Epsilon as a function of the expanded ensemble state
	std::vector < double > sigmaM_; //!< Sigma as a function of the expanded ensemble state
	std::vector < double > uShiftM_; //!< Energy shift as a function of the expanded ensemble state
};

/*!
 * Force-Shifted Lennard-Jones Potential
 * Parameters should be specified in the following order: { epsilon, sigma, r_cut, Mtot }
 */
class fsLennardJones : public pairPotential {
public:
	~fsLennardJones () {;}
	void setParameters (const std::vector < double > params);
	double energy (const atom* a1, const atom* a2, const std::vector < double > &box);
	double tailCorrection (const double rhoBath);
	double rcut ();

private:
	std::vector < double > epsM_; //!< Epsilon as a function of the expanded ensemble state
	std::vector < double > sigmaM_; //!< Sigma as a function of the expanded ensemble state
};

/*!
 * Tabulated Potential
 * Parameters should be specified in the following order: { r_cut, r_shift, u_shift, u_infinity, Mtot }
 */
class tabulated : public pairPotential {
public:
	~tabulated () {;}
	void setParameters (const std::vector < double > params);
	void loadPotential(std::string filename);
	double energy (const atom* a1, const atom* a2, const std::vector < double > &box);
	double tailCorrection (const double rhoBath);
	double rcut ();

private:
	double start; //!< r To start from
	double dr; //!< Increment for r
	std::vector <double> table, mScale;
};

/*!
 * Square-well potential
 * Parameters should be specified in the following order: { sigma, wellwidth, welldepth, Mtot }
 */
class squareWell : public pairPotential {
public:
	~squareWell () {;}
	void setParameters (const std::vector < double > params);
	double energy (const atom* a1, const atom* a2, const std::vector < double > &box);
	double tailCorrection (const double rhoBath);
	double rcut ();

private:
	std::vector < double > sigmaM_; //!< Hard sphere overlap as a function of expanded ensemble stage
	std::vector < double > rangeM_;	//!< Range of interaction as a function of expanded ensemble stage
	std::vector < double > epsM_; //!< Epsilon as a function of the expanded ensemble state
};

/*!
 * Hard-core potential
 * Parameters should be specified in the following order: { sigma, Mtot }
 */
class hardCore : public pairPotential {
public:
	~hardCore () {;}
	void setParameters (const std::vector <double> params);
	double energy (const atom* a1, const atom* a2, const std::vector < double > &box);
	double tailCorrection (const double rhoBath);
	double rcut ();

private:
	std::vector < double > sigmaM_; //!< Hard sphere overlap as a function of expanded ensemble stage
};

#endif
