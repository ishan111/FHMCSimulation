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

#define NUM_INFINITY std::numeric_limits<double>::max()

/*!
 * Abstract base class which defines all pair potentials.
 */
class pairPotential {
public:
	pairPotential () { paramsAreSet_ = false; }
	virtual ~pairPotential () {;}
	
	bool useTailCorrection;
	virtual double energy (const double r) = 0;	
	virtual double tailCorrection (const double rhoBath) = 0;
	virtual void setParameters (const std::vector < double > params) = 0;
	virtual double rcut () = 0;	//!< All potentials should be able to return their r_{cut} values so neighbor lists, etc. can use them
	
	inline double energy (const std::vector < double > &p1, const std::vector < double > &p2, const std::vector < double > &box) { return energy(sqrt(pbc_dist2(p1, p2, box))); } //!< Calculates the minimum image distance and invokes the energy call 
	void savePotential(std::string filename, double start, double dr);
//protected:
	std::vector < double > params_; //!< Parameters (constants) that are needed to calculate U(r)
	bool paramsAreSet_; //!< Logical check if the paramters for this potential have been specified by the user
};

/*!
 * Lennard-Jones Potential
 * Parameters should be specified in the following order: { epsilon, sigma, r_cut, u_shift }
 */
class lennardJones : public pairPotential {
public:
	~lennardJones () {;}
	void setParameters (const std::vector < double > params);
	double energy (const double r);
	double tailCorrection (const double rhoBath);
	double rcut ();
};

/*!
 * Tabulated Potential
 * Parameters should be specified in the following order: { r_cut, r_shift, u_shift, u_infinity }
 */
class tabulated : public pairPotential {
private:
	std::vector <double> table;
	double start; //!< r To start from
	double dr; //!< Increment for r
public:
	~tabulated () {;}
	void setParameters (const std::vector < double > params);
	void loadPotential(std::string filename);
	double energy (const double r);
	double tailCorrection (const double rhoBath);
	double rcut ();
};

/*!
 * Square-well potential
 * Parameters should be specified in the following order: { sigma, wellwidth, welldepth }
 */
class squareWell : public pairPotential {
public:
	~squareWell () {;}
	void setParameters (const std::vector < double > params);
	double energy (const double r);
	double tailCorrection (const double rhoBath);
	double rcut ();
};

/*!
 * Hard-core potential
 * Parameters should be specified in the following order: { sigma }
 */
class hardCore : public pairPotential {
public:
	~hardCore () {;}
	void setParameters (const std::vector <double> params);
	double energy (const double r);
	double tailCorrection (const double rhoBath);
	double rcut ();
};

#endif
