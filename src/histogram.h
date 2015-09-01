#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include "global.h"
#include "utilities.h"

class histogram {
public:
	~histogram () {;}
    	histogram (const std::vector <double> lbound, const std::vector <double> ubound, const std::vector <long long unsigned int> nbins);
    
    	void print (const std::string fileName);
    	void increment (const long long unsigned int address, const double val);
    	void increment (const std::vector <double> &coords, const double val);
    	const long long unsigned int getAddress (const std::vector <double> &coords);
    	const std::vector <double> getCoords (long long unsigned int address);
    	std::vector <double> getHistogram () { return h_; } //!< Return the current histogram
    
private:
    	int dim_; //!< Dimensionality of the histogram
    	long long unsigned int size_; //!< Total size of histogram
    	std::vector <long long unsigned int> nbins_; //!< Number of bins along each dimension
    	std::vector <long long unsigned int> widths_; //!< Total widths (number of bins) of each iteratively projected dimension
    	std::vector <double> h_; //!< Histogram
    	std::vector <double> ubound_; //!< Upper bound for each dimension
    	std::vector <double> lbound_; //!< Lower bound for each dimension
    	std::vector <double> delta_; //!< Magnitude of the width of a bin in each dimension
};

#endif
