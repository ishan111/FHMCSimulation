#ifndef BIAS_H_
#define BIAS_H_

#include <vector>
#include <string>
#include <climits>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "global.h"
#include "utilities.h"

// netCDF if enabled
#ifdef NETCDF_CAPABLE
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;
#endif

//! For biasing, need to address things with large indices 
#define __BIAS_INT_TYPE__ long long int

/*! Jacobi method for evaluating natural logarithm of a sum of exponentials
 * E.g. ln( exp(a) + exp(b) ) = specExp(a, b)
 */
#define specExp(a, b) std::max(a, b) + log(1.0 + exp(-std::fabs(a-b)))

//! Transition Matrix Monte Carlo biasing class
class tmmc {
public:
	tmmc () {};
	tmmc (const int Nmax, const int Nmin, const int Mtot, const long long int tmmcSweepSize, const std::vector <double> box);
	
	void updateC (const int Nstart, const int Nend, const int Mstart, const int Mend, const double pa);
	void calculatePI ();
	void print (const std::string fileName, bool printC = false);
	void readC (const std::string fileName);
	void setlnPI (const std::vector < double > &lnPIguess) { if (lnPIguess.size() == lnPI_.size()) lnPI_ = lnPIguess; } //!< Blindly assign a guess of the macrostate distribution
	void iterateForward ();
	void dumpVisited (const std::string fileName);
	bool checkFullyVisited ();
	const __BIAS_INT_TYPE__ getTransitionAddress (const int Nstart, const int Nend, const int Mstart, const int Mend);
	const __BIAS_INT_TYPE__ getAddress (const int Nval, const int Mval);	
	const long long int numSweeps () { return nSweeps_; }
	const double getBias (const int address) { return -lnPI_[address]; }
	const std::vector < double > getC () { return C_; } //!< Return the collection matrix as it is
	const std::vector < double > getHC () { return HC_; }
	
private:
	int Nmax_; //!< Maximum number of total species in a simulation
	int Nmin_; //!< Minimum number of total species in a simulation
	int Mtot_; //!< Total number of expanded ensemble states in a simulation	
	long long int nSweeps_; //!< Number of sweeps that have been completed
	long long int tmmcSweepSize_; //!< Number of times the each point in the collection matrix must be visited for a "sweep" to be considered finished
	std::vector < double > C_; //!< Collection matrix
	std::vector < double > P_; //!< Probability matrix
	std::vector < double > lnPI_; //!< Estimated (natural logarithm of) macrostate density used as bias
	std::vector < double > box_; //!< Size of the simulation box this is originating from
	std::vector < double > HC_; //!< Count the number of time each state and transition has been sampled 
};

//! Wang-Landau biasing class
class wala {
public:
	wala () {};
	wala (const double lnF, const double g, const double s, const int Nmax, const int Nmin, const int Mtot, const std::vector <double> box);
	
	void update (const int Nval, const int Mval);
	void iterateForward ();
	void print (const std::string fileName, const bool printH = false);
	void readlnPI (const std::string fileName);
	void setlnPI (const std::vector < double > &lnPIguess) { if (lnPIguess.size() == lnPI_.size()) lnPI_ = lnPIguess; } //!< Blindly assign a guess of the macrostate distribution
	bool evaluateFlatness ();
	const __BIAS_INT_TYPE__ getAddress (const int Nval, const int Mval);
	const double lnF () { return lnF_; }
	const double getBias (const int address) { return -lnPI_[address]; }
	const std::vector <double> getlnPI () { return lnPI_; } //!< Return the current estimate of the macrostate distribution
	const std::vector <double> getH () { return H_; } //!< Return the visited-states histogram
	
private:
	int Mtot_; //!< Total number of expanded ensemble states in a simulation	
	double lnF_; //!< Factor to add to the macrostate density matrix each time a state is visited
	double g_; //!< Factor by which lnF_ is occasionally reduced by as the simulation progresses
	double s_; //!< Ratio between minimum of the visited states histogram and its mean which determines "flatness"
	long long int Nmax_; //!< Maximum total number of atoms in a simulation
	long long int Nmin_; //!< Minimum total number of atoms in a simulation
	std::vector <double> H_; //!< Histogram of visited states
	std::vector <double> lnPI_; //!< Estimated (natural logarithm of) macrostate density used as bias
	std::vector <double> box_; //!< Size of the simulation box this is originating from
};

#endif
