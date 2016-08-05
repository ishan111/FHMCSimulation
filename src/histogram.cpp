#include "histogram.h"
#include <iostream>

/*!
 * Trim the size of the histogram to remove leading and trailing zeros.
 */
void dynamic_one_dim_histogram::trim_edges () {
	long unsigned int leading = 0, trailing = 0;
	for (std::deque < double >::iterator it = h_.begin(); it != h_.end(); ++it) {
		if (*it <= 0) {
			leading++;
		} else {
			break;
		}
	}
	for (std::deque < double >::reverse_iterator rit = h_.rbegin(); rit != h_.rend(); ++rit) {
		if (*rit <= 0) {
			trailing++;
		} else {
			break;
		}
	}

	if (leading + trailing >= h_.size()) {
		throw customException ("Cannot trim dynamic_one_dim_histogram because it is empty");
	}
	
	nbins_ -= (leading + trailing);
	lb_ += leading*delta_;
	ub_ -= trailing*delta_;
	
	for (unsigned int i = 0; i < leading; ++i) {
		h_.pop_front();
	}
	for (unsigned int i = 0; i < trailing; ++i) {
		h_.pop_back();
	}
}

/*!
 * Initialize histogram and its bounds. 
 */
void dynamic_one_dim_histogram::initialize_ (const double lb, const double ub, const double delta) {
	if (lb > ub) {
        throw customException ("Lower bound >= upper bound illegal for a dynamic_one_dim_histogram");
    }
    if (delta <= 0) {
        throw customException ("Bin width must be > 0 for a dynamic_one_dim_histogram");
    }
	delta_ = delta; 
 	lb_ = lb;
 	ub_ = ub;
 	nbins_ = ceil((ub - lb)/delta); 
 	if (fabs(round((ub - lb)/delta) - ((ub - lb)/delta)) < 1.0e-6) {
 		nbins_++; // include endpoint
 	}
    
    // initialize the histogram to 0
    h_.resize(0);
    try {
        h_.resize(nbins_, 0);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory for dynamic_one_dim_histogram");
    }
}

/*!
 * Re-initialize histogram and its bounds.  All entries are zeroed.
 */
void dynamic_one_dim_histogram::reinitialize (const double lb, const double ub, const double delta) {
	try {
		initialize_ (lb, ub, delta);
	} catch (customException &ce) {
		throw customException ("Unable to re-initialize dynamic_one_dim_histogram");
	}
}

/*!
 * Instantiate a 1D histogram that grow as needed to record values.  
 * A bin is considered "centered" on its value.
 *
 * \param [in] lb Lower bound
 * \param [in] ub Upper bound
 * \param [in] delta Bin width
 */
dynamic_one_dim_histogram::dynamic_one_dim_histogram (const double lb, const double ub, const double delta) {
	try {
		initialize_ (lb, ub, delta);
	} catch (customException &ce) {
		throw customException ("Unable to initialize dynamic_one_dim_histogram");
	}
}

/*!
 * Add a number of bins to the beginning of the histogram.
 *
 * \param [in] nbins Number of bins to add to the beginning
 */
void dynamic_one_dim_histogram::prepend_bins (const unsigned int nbins) {
	for (unsigned int i = 0; i < nbins; i++) {
		try {
			h_.push_front (0.0);
			nbins_++;
			lb_ -= delta_;
		} catch (std::bad_alloc &ba) {
        	throw customException ("Out of memory for dynamic_one_dim_histogram");
   		}
	}
}

/*!
 * Add a number of bins to the end of the histogram.
 *
 * \param [in] nbins Number of bins to add to the end
 */
void dynamic_one_dim_histogram::append_bins (const unsigned int nbins) {
	for (unsigned int i = 0; i < nbins; i++) {
		try {
			h_.push_back (0.0);
			nbins_++;
			ub_ += delta_;
		} catch (std::bad_alloc &ba) {
        	throw customException ("Out of memory for dynamic_one_dim_histogram");
   		}
	}
}

/*!
 * Record an entry in the histogram at the bin corresponding to where value falls.
 *
 * \param [in] value Raw value, bin this corresponds to is internally calculated
 */
void dynamic_one_dim_histogram::record (const double value) {
	int bin = round((value - lb_)/delta_); // this "centers" the bin
	if (bin < 0) {
		// prepend and fill
		try {
			prepend_bins(-bin);
		} catch (customException &ce) {
			std::string a = "Unable to prepend dynamic_one_dim_histogram: ", b = ce.what();
			throw customException (a+b);
		}
		bin = round((value - lb_)/delta_); 
	} else if (bin >= nbins_) {
		// append and fill
		try {
			append_bins (bin - nbins_ + 1);
		} catch (customException &ce) {
			std::string a = "Unable to append dynamic_one_dim_histogram: ", b = ce.what();
			throw customException (a+b);
		}
		bin = round((value - lb_)/delta_); 
	}
	h_[bin] += 1.0;
}

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
    	try {
        	counter_.resize(size_, 0);
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
        	counter_[address] += 1.0;
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
    	counter_[address] += 1.0;
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
 
