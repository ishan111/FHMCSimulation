#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include <vector>
#include <deque>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include "global.h"
#include "utilities.h"

class histogram {
public:
	~histogram () {;}
		histogram () {;}
    	histogram (const std::vector <double> lbound, const std::vector <double> ubound, const std::vector <long long unsigned int> nbins);

    	void print (const std::string fileName);
    	void increment (const long long unsigned int address, const double val);
    	void increment (const std::vector <double> &coords, const double val);
    	const long long unsigned int getAddress (const std::vector <double> &coords);
    	const std::vector <double> getCoords (long long unsigned int address);
    	std::vector <double> getRawHistogram () { return h_; } //!< Return the current histogram
    	std::vector <double> getCounter () { return counter_; } //!< Return the current histogram counter

private:
    	int dim_; //!< Dimensionality of the histogram
    	long long unsigned int size_; //!< Total size of histogram
    	std::vector <long long unsigned int> nbins_; //!< Number of bins along each dimension
    	std::vector <long long unsigned int> widths_; //!< Total widths (number of bins) of each iteratively projected dimension
    	std::vector <double> h_; //!< Histogram
    	std::vector <double> counter_; //!< Histogram counter
    	std::vector <double> ubound_; //!< Upper bound for each dimension
    	std::vector <double> lbound_; //!< Lower bound for each dimension
    	std::vector <double> delta_; //!< Magnitude of the width of a bin in each dimension
};

class dynamic_one_dim_histogram {
	public:
		~dynamic_one_dim_histogram () {;}
		dynamic_one_dim_histogram () {;}
		dynamic_one_dim_histogram (const double lb, const double ub, const double delta);

		void reinitialize (const double lb, const double ub, const double delta);
		void trim_edges ();
		void prepend_bins (const unsigned int nbins);
		void append_bins (const unsigned int nbins);
		void record (const double value);
		const double get_delta () { return delta_; }
		const double get_lb () { return lb_; }
		const double get_ub () { return ub_; }
		std::deque < double > get_hist () { return h_; }

	private:
		int nbins_; //!< Number of bins
		double lb_; //!< Lower bound
		double ub_; //!< Upper bound
		double delta_; //!< Width of a bin
		std::deque < double > h_; //!< Histogram
		void initialize_ (const double lb, const double ub, const double delta);
};

#endif
