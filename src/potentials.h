#ifndef POTENTIALS_H_
#define POTENTIALS_H_

#include "global.h"
#include "utilities.h"
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream> 

/*!
 * Abstract base class which defines all pair potentials.
 */
class pairPotential {
public:
	pairPotential () { paramsAreSet_ = false; }
	bool useTailCorrection;
	virtual double energy (const double r) = 0;	
	virtual double tailCorrection (const double rhoBath) = 0;
	virtual void setParameters (const std::vector < double > params) = 0;
	virtual double rcut () = 0;	//!< All potentials should be able to return their r_{cut} values so neighbor lists, etc. can use them
	
	inline double energy (const std::vector < double > &p1, const std::vector < double > &p2, const std::vector < double > &box) { return energy(sqrt(pbc_dist2(p1, p2, box))); } //!< Calculates the minimum image distance and invokes the energy call 
	void savePotential(std::string filename, double start, double dr);
//protected:
	std::vector < double > params_;
	bool paramsAreSet_;
};

/*!
 * Lennard-Jones Potential
 * Parameters should be specified in the following order: { epsilon, sigma, rcut, ushift }
 */
class lennardJones : public pairPotential {
public:
	void setParameters (const std::vector < double > params);
	double energy (const double r);
	double tailCorrection (const double rhoBath);
	double rcut ();
};

/*!
 * Tabulated Potential
 * Parameters should be specified in the following order: { rcut, rshift, ushift, uinfinity }
 */
class tabulated : public pairPotential {
private:
	std::vector<double> table;
	double start;
	double dr;
public:
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
	void setParameters (const std::vector <double> params);
	double energy (const double r);
	double tailCorrection (const double rhoBath);
	double rcut ();
};

#endif
