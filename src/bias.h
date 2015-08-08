#ifndef BIAS_H_
#define BIAS_H_

#include <vector>
#include <string>
#include <climits>
#include <cfloat>
#include <cmath>
#include <fstream>
#include "global.h"
#include "utilities.h"

//! For biasing, need to address things with large indices 
#define __BIAS_INT_TYPE__ long long int

/*! Jacobi method for evaluating natural logarithm of a sum of exponentials
 * E.g. ln( exp(a) + exp(b) ) = specExp(a, b)
 */
#define specExp(a, b) std::max(a, b) + log(1.0 + exp(-std::abs(a-b)))

//! Transition Matrix Monte Carlo biasing class
class tmmc {
public:
	tmmc () {};
	tmmc (const int nSpec, const std::vector <int> &Nmax, const std::vector <int> &Nmin);

	void updateC (const std::vector <int> &Nstart, const std::vector <int> &Nend, const double pa);
	void calculatePI ();
	void print (const std::string fileName, bool printC = false);
	void readC (const std::string fileName);
	void readlnPI (const std::string fileName);
	const __BIAS_INT_TYPE__ getTransitionAddress (const std::vector <int> &Nstart, const std::vector <int> &Nend);
	const __BIAS_INT_TYPE__ getTransitionAddress (const std::vector <int> &Nstart, const int specId, const int addOrSubtract);
	const __BIAS_INT_TYPE__ getAddress (const std::vector <int> &Nval);	
	const double getBias (const int address) { return -lnPI_[address]; }
	
private:
	void difference_ (const std::vector <int> &Nstart, const std::vector <int> &Nend, int &specId, int &addOrSubtract);
	int nSpec_; //!< Number of species in simulation
	std::vector <int> Nmax_; //!< Maximum number of each species in a simulation
	std::vector <int> Nmin_; //!< Minimum number of each species in a simulation
	std::vector <int> W_; //!< Width of each dimension (Nmax - Nmin + 1)
	std::vector <double> C_; //!< Collection matrix
	std::vector <double> P_; //!< Probability matrix
	std::vector <double> lnPI_; //!< Estimated (natural logarithm of) macrostate density used as bias
};

//! Wang-Landau biasing class
class wala {
public:
	wala () {};
	wala (const double lnF, const double g, const double s, const int nSpec, const std::vector <int> &Nmax, const std::vector <int> &Nmin);
	
	void update (const std::vector <int> &Nval);
	void iterateForward ();
	void print (const std::string fileName, const bool printH = false);
	void readlnPI (const std::string fileName);
	bool evaluateFlatness ();
	const __BIAS_INT_TYPE__ getAddress (const std::vector <int> &Nval);
	const double lnF () { return lnF_; }
	const double getBias (const int address) { return -lnPI_[address]; }
	
private:
	int nSpec_; //!< Number of species in simulation
	double lnF_; //!< Factor to add to the macrostate density matrix each time a state is visited
	double g_; //!< Factor by which lnF_ is occasionally reduced by as the simulation progresses
	double s_; //!< Ratio between minimum of the visited states histogram and its mean which determines "flatness"
	std::vector <int> Nmax_; //!< Maximum number of each species in a simulation
	std::vector <int> Nmin_; //!< Minimum number of each species in a simulation
	std::vector <int> W_; //!< Width of each dimension (Nmax - Nmin + 1)
	std::vector <double> H_; //!< Histogram of visited states
	std::vector <double> lnPI_; //!< Estimated (natural logarithm of) macrostate density used as bias
};

#endif
