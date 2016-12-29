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
 * Initialize histogram and its bounds. Aligns against the lower bound, and uses this as the bin's "center."  The upper bound is internally re-calculated from the lower bound and the value of delta.
 * Bins include values between [lo, hi) where a bin's value = (low+hi)/2.
 *
 * \param [in] lb Lowest value that must be covered by histogram. This is used to align the histogram when ub and lb not integer delta's apart.
 * \param [in] ub Largest value that must be covered by histogram. Less relevant than lb, which is used for alignment.  This value is only used to determined the number of bins necessary and will be changed internally.
 * \param [in] delta Bin width.
 */
void dynamic_one_dim_histogram::initialize_ (const double lb, const double ub, const double delta) {
	if (lb > ub) {
        throw customException ("Lower bound >= upper bound illegal for a dynamic_one_dim_histogram");
    }
    if (delta <= 0) {
        throw customException ("Bin width must be > 0 for a dynamic_one_dim_histogram");
    }
	tol_ = std::numeric_limits < double >::epsilon();

	delta_ = delta;
 	lb_ = lb;

	const double x = (ub - (lb_- delta_/2.0))/delta_;
	if (fabs(round(x) - x) < tol_) {
		// ub is on the "edge" of the high end of a bin, round for numerical stability and then include next bin, because bins are [lo, hi)
		nbins_ = round(x)+1;
	} else {
		nbins_ = ceil(x);
	}
	ub_ = (nbins_-1)*delta_+lb_;

    // initialize the histogram to 0
	try {
        h_.assign(nbins_, 0);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory for dynamic_one_dim_histogram");
    }
}

/*!
 * Set the histogram.  Intended to be used to restart from a checkpoint.
 *
 * \param [in] h Histogram to set to.
 */
void dynamic_one_dim_histogram::set_hist (const std::deque < double > h) {
	if (h.size() != h_.size()) {
		throw customException ("Histogram using to set is not the same as inherent, aborting");
	} else {
		h_ = (std::deque < double >)h;
	}
}

/*!
 * Re-initialize histogram and its bounds.  All entries are zeroed. Aligns against the lower bound, and uses this as the bin's "center."
 * The upper bound is internally re-calculated from the lower bound and the value of delta. Bins include values between [lo, hi) where a bin's value = (low+hi)/2.
 *
 * \param [in] lb Lowest value that must be covered by histogram. This is used to align the histogram when ub and lb not integer delta's apart.
 * \param [in] ub Largest value that must be covered by histogram. Less relevant than lb, which is used for alignment.  This value is only used to determined the number of bins necessary.
 * \param [in] delta Bin width.
 */
void dynamic_one_dim_histogram::reinitialize (const double lb, const double ub, const double delta) {
	try {
		initialize_ (lb, ub, delta);
	} catch (customException &ce) {
		throw customException ("Unable to re-initialize dynamic_one_dim_histogram");
	}
}

/*!
 * Instantiate a 1D histogram that grow as needed to record values. Aligns against the lower bound, and uses this as the bin's "center."
 * The upper bound is internally re-calculated from the lower bound and the value of delta. Bins include values between [lo, hi) where a bin's value = (low+hi)/2.
 *
 * \param [in] lb Lowest value that must be covered by histogram. This is used to align the histogram when ub and lb not integer delta's apart.
 * \param [in] ub Largest value that must be covered by histogram. Less relevant than lb, which is used for alignment.  This value is only used to determined the number of bins necessary.
 * \param [in] delta Bin width.
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
	if (std::abs(bin) < tol_) {
		// prevent -0 case and set to 0
		bin = 0;
	}

	if (bin < 0) {
		// prepend and fill
		try {
			prepend_bins(-bin);
		} catch (customException &ce) {
			std::string a = "Unable to prepend dynamic_one_dim_histogram: ", b = ce.what();
			throw customException (a+b);
		}
	} else if (bin >= nbins_) {
		// append and fill
		try {
			append_bins (bin - nbins_ + 1);
		} catch (customException &ce) {
			std::string a = "Unable to append dynamic_one_dim_histogram: ", b = ce.what();
			throw customException (a+b);
		}
	}

	// re-calculate after lb potentially adjusted
	bin = round((value - lb_)/delta_);
	if (std::abs(bin) < tol_) {
			// prevent -0 case and set to 0
			bin = 0;
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
histogram::histogram (const std::vector < double > lbound, const std::vector < double > ubound, const std::vector < long long unsigned int > nbins) {
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
    	if (lbound[i] > ubound[i]) {
        	throw customException ("Lower bound > upper bound illegal for a histogram");
    	}

		if (lbound[i] < ubound[i]) {
	        if (nbins[i] <= 1) {
        		throw customException ("Must > 1 bins for each dimensions in the histogram");
        	}
	        size_ *= nbins[i];
        	delta_[i] = (ubound[i] - lbound[i])/(nbins[i]-1);
		} else {
			// special case when upper and lower bound are the same (nbins = 1)
			if (nbins[i] != 1) {
				throw customException ("nbins must be 1 if upper and lower bounds are equal in histogram");
			}
			size_ *= nbins[i];
			delta_[i] = 1.0; // arbitrary
		}

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
 * Print RAW (UN-NORMALIZED) histogram to file.
 *
 * \param [in] fileName Name of file to print to
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
		of << "# dim_"+std::to_string(i+1)+"_upper_bound:" << ubound_[i] << std::endl;
		of << "# dim_"+std::to_string(i+1)+"_lower_bound:" << lbound_[i] << std::endl;
		of << "# dim_"+std::to_string(i+1)+"_number_of_bins:" << nbins_[i] << std::endl;
	}
	for (unsigned long long int i = 0; i < h_.size(); ++i) {
		of << h_[i] << std::endl;
	}
	of.close();
}

/*!
 * Assign the histogram and its corresponding counter.
 *
 * \param [in] h histogram
 * \param [in] ctr Counter
 */
void histogram::set (const std::vector <double> &h, const std::vector <double> &ctr) {
	if (h.size() != ctr.size()) {
		throw customException ("Cannot set the histogram since counter and histogram have different lengths");
	}
	if (h.size() != h_.size()) {
		throw customException ("Cannot set the histogram since new histogram has different length compared to current one");
	}
	if (ctr.size() != counter_.size()) {
		throw customException ("Cannot set the histogram since new counter has different length compared to current one");
	}
	h_ = h;
	counter_ = ctr;
}
