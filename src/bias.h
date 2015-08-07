#ifndef BIAS_H_
#define BIAS_H_

#include <vector>
#include <string>
#include <climits>
#include <cfloat>
#include <cmath>
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
	
	const __BIAS_INT_TYPE__ getTransitionAddress (const std::vector <int> &Nstart, const std::vector <int> &Nend);
	const __BIAS_INT_TYPE__ getTransitionAddress (const std::vector <int> &Nstart, const int specId, const int addOrSubtract);
	const __BIAS_INT_TYPE__ getAddress (const std::vector <int> &Nval);
	
	void updateC (const std::vector <int> &Nstart, const std::vector <int> &Nend, const double pa);
	void calculatePI ();
	const double getBias (const int address) { return -lnPI_[address]; }
	
private:
	void difference_ (const std::vector <int> &Nstart, const std::vector <int> &Nend, int &specId, int &addOrSubtract);
	int nSpec_; // number of species in simulation
	std::vector <int> Nmax_, Nmin_, W_;
	std::vector <double> C_, P_, lnPI_;
};

//! Wang-Landau biasing class
class wala {
public:
	wala () {};
	wala (const double lnF, const double g, const double s, const int nSpec, const std::vector <int> &Nmax, const std::vector <int> &Nmin);
	
	const __BIAS_INT_TYPE__ getAddress (const std::vector <int> &Nval);
	const double lnF () { return lnF_; }
	const double getBias (const int address) { return -lnPI_[address]; }
	bool evaluateFlatness ();
	void update (const std::vector <int> &Nval);
	void iterateForward ();
	
private:
	int nSpec_;
	double lnF_, g_, s_;
	std::vector <int> Nmax_, Nmin_, W_;
	std::vector <double> H_, lnPI_;
};

#endif
