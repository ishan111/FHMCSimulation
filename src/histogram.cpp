#include "histogram.h"
#include <iostream>

/*!
 * Instantiate a multidimensional histogram.  Bounds and widths must be specified for each dimension.
 * A bin is considered "centered" on its value.
 *
 * \param [in] lbound Vector of lower bounds for each dimension
 * \param [in] ubound Vector of upper bounds for each dimension
 * \param [in] nbins Number of bins to use along each dimension
 */
histogram::histogram (const std::vector <double> lbound, const std::vector <double> ubound, const std::vector <long long unsigned int> nbins) {
    if (lbound.size() != ubound.size()) {
        throw customException ("Upper and lower bounds for histogram do have the same size");
    }
    if (nbins.size() != lbound.size()) {
        throw customException ("Number of bins for histogram's dimensions does not have the same size as its bounds");
    }

    dim_ = nbins.size();
    widths_.resize(dim_, 0);
    delta_.resize(dim_, 0);
 
    size_ = 1;
    for (unsigned int i = 0; i < dim_; ++i) {
        if (lbound[i] >= ubound[i]) {
            throw customException ("Lower bound >= upper bound illegal for a histogram");
        }
        if (nbins[i] <= 1) {
            throw customException ("Must > 1 bins for each dimensions in the histogram");
        }
        size_ *= nbins[i];
        delta_[i] = (ubound[i] - lbound[i])/(nbins[i]-1);
         
        // build projected widths
        if (i == 0) {
            widths_[i] = 1;
        } else {
            widths_[i] = widths_[i-1]*nbins[i-1];
        }
    }
        
    lbound_ = lbound;
    ubound_ = ubound;
    nbins_ = nbins;
    
    // initialize the histogram to 0
    try {
        h_.resize(size_, 0);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory for histogram");
    }
}

/*!
 * Increment the histogram by a given value at a given address.
 *
 * \param [in] address Address of the histogram to increment
 * \param [in] val Value to add to the histogram at address
 */
void histogram::increment (const long long unsigned int address, const double val) {
    if (address < size_) {
        h_[address] += val;
    } else {
        throw customException ("Histogram address out of bounds");
    }
}

/*!
 * Increment the histogram by a given value at a given coordinate.
 *
 * \param [in] coords Vector of coordinates correponding to a location in the histogram to increment
 * \param [in] val Value to add to the histogram at address
 */
void histogram::increment (const std::vector <double> &coords, const double val) {
    long long unsigned int address = 0;
    try {
        address = getAddress (coords);
    } catch (customException &ce) {
        throw customException ("Histogram address out of bounds");
    }
    h_[address] += val;
}

/*!
 * Get the linear address of the multidimensional coordinate.
 *
 * \param [in] coords Coordinates
 */
const long long unsigned int histogram::getAddress (const std::vector <double> &coords) {
    if (coords.size() != dim_) {
        throw customException ("Illegal number of coordinate dimensions, cannot locate histogram address");
    }
    long long unsigned int address = 0;
    for (unsigned int i = 0; i < dim_; ++i) {
        address += round((coords[i] - lbound_[i])/delta_[i])*widths_[i]; // will work safely for integers too
    }
    return address;
}

/*!
 * Given an address, return the (center of the) coordinate this refers to.
 *
 * \param [in] address Address to check
 */
const std::vector <double> histogram::getCoords (long long unsigned int address) {
    std::vector <double> coords (dim_, 0);
    if (address >= size_) {
        throw customException ("Histogram address out of bounds");
    }
    
    for (unsigned int i = dim_-1; i > 0; --i) {
        long long int diff = floor(address/widths_[i]);
        coords[i] = diff*delta_[i] + lbound_[i];
        address -= diff*widths_[i];
    }
    coords[0] = address*delta_[0] + lbound_[0];
    
    return coords;
}

/*!
 * Print a histogram to file
 */
void histogram::print (const std::string fileName) {
	// Print histogram
	std::ofstream of;
	of.open(fileName.c_str(), std::ofstream::out);
	if (!of.is_open()) {
		throw customException ("Unable to write histogram to "+fileName);
	}
	of << "# Histogram in single row (vectorized) notation." << std::endl;
	for (unsigned int i = 0; i < dim_; ++i) {
		of << "# dim_"+sstr(i+1)+"_upper_bound:" << ubound_[i] << std::endl;
		of << "# dim_"+sstr(i+1)+"_lower_bound:" << lbound_[i] << std::endl;
		of << "# dim_"+sstr(i+1)+"_number_of_bins:" << nbins_[i] << std::endl;
	}
	
	for (unsigned long long int i = 0; i < h_.size(); ++i) {
		of << h_[i] << std::endl;
	}
	of.close();
}
 