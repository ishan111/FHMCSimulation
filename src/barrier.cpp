#include "barrier.h"

/*!
 * Instantiate a hard wall with boundaries in the +/- z direction.
 *
 * \param [in] lb z-Position of the lower wall
 * \param [in] ub z-Position of the upper wall
 * \param [in] sigma Hard-sphere diameter the species this wall interacts with can approach within
 * \param [in] M Total number of expanded ensemble states possible for this atom type (defaults to 1)
 */
hardWallZ::hardWallZ (const double lb, const double ub, const double sigma, const int M) {
    if (lb >= ub) {
        throw customException ("hardWallZ must have lower bound < upper bound");
    }
    if (sigma < 0) {
        throw customException ("hardWallZ must have sigma >= 0");
    }
    if (M < 1) {
        throw customException ("hardWallZ must have M >= 1");
    }
    
    sigma_ = sigma;
    ub_ = ub;
    lb_ = lb;
    M_ = M;
}

/*!
 * Return whether or not a point falls between the walls (subject to a hard-sphere exclusion radius). Sigma is scaled linearly with expanded ensemble state.
 *
 * \param [in] a1 Pointer to atom with position to test - this does NOT need to be in the simulation box a priori
 * \param [in] box Simulation box
 */
bool hardWallZ::inside (const atom *a1, const std::vector < double > &box) {
    std::vector < double > p = a1->pos;
    pbc (p, box);
    
    double sig = sigma_;
    if (a1->mState > 0) {
        sig = (sigma_/M_)*a1->mState;
    }
    if (a1->mState < 0 || a1->mState > M_-1) {
        throw customException ("mState out of bounds for hardWallZ");
    }
    
    if (p[2] >= ub_ - sig/2.0 || p[2] <= lb_ + sig/2.0) {
        return false;
    } else {
        return true;
    }
}

/*!
 * Interaction energy with the wall. Sigma is scaled linearly with expanded ensemble state.
 *
 * \param [in] a1 Pointer to atom with position to test - this does NOT need to be in the simulation box a priori
 * \param [in] box Simulation box
 */
double hardWallZ::energy (const atom *a1, const std::vector < double > &box) {
    std::vector < double > p = a1->pos;
    pbc (p, box);
    
    double sig = sigma_;
    if (a1->mState > 0) {
        sig = (sigma_/M_)*a1->mState;
    }
    if (a1->mState < 0 || a1->mState > M_-1) {
        throw customException ("mState out of bounds for hardWallZ");
    }
    
    if (p[2] >= ub_ - sig/2.0 || p[2] <= lb_ + sig/2.0) {
        return NUM_INFINITY;
    } else {
        return 0.0;
    }
}

/*!
 * Instantiate a square well wall with boundaries in the +/- z direction.
 *
 * \param [in] lb z-Position of the lower wall
 * \param [in] ub z-Position of the upper wall
 * \param [in] sigma Hard-sphere diameter the species this wall interacts with can approach within
 * \param [in] range Distance normal to the wall's surface where there is an interaction
 * \param [in] eps Magnitude of the wall interaction (U = -eps)
 * \param [in] M Total number of expanded ensemble states possible for this atom type (defaults to 1)
 */
squareWellWallZ::squareWellWallZ (const double lb, const double ub, const double sigma, const double range, const double eps, const int M) {
    if (lb >= ub) {
        throw customException ("squareWellWallZ must have lower bound < upper bound");
    }
    if (sigma < 0) {
        throw customException ("squareWellWallZ must have sigma >= 0");
    }
    if (range < 0) {
        throw customException ("squareWellWallZ must have range >= 0");
    }
    if (eps < 0) {
        throw customException ("squareWellWallZ must have eps >= 0");
    }
    if (sigma/2.0 >= range) {
        throw customException ("squareWellWallZ must have sigma/2 < range to have a finite range of interaction");
    }
    if (M < 1) {
        throw customException ("squareWellWallZ must have M >= 1");
    }
    
    sigma_ = sigma;
    range_ = range;
    ub_ = ub;
    lb_ = lb;
    M_ = M;
}

/*!
 * Return whether or not a point falls between the walls (subject to a hard-sphere exclusion radius).  Sigma is scaled linearly with expanded ensemble state.
 *
 * \param [in] a1 Pointer to atom with position to test - this does NOT need to be in the simulation box a priori
 * \param [in] box Simulation box
 */
bool squareWellWallZ::inside (const atom *a1, const std::vector < double > &box) {
    std::vector < double > p = a1->pos;
    pbc (p, box);
    
    double sig = sigma_;
    if (a1->mState > 0) {
        sig = (sigma_/M_)*a1->mState;
    }
    if (a1->mState < 0 || a1->mState > M_-1) {
        throw customException ("mState out of bounds for squareWellWallZ");
    }
    
    if (p[2] >= ub_ - sig/2.0 || p[2] <= lb_ + sig/2.0) {
        return false;
    } else {
        return true;
    }
}

/*!
 * Interaction energy with the wall. Sigma and epsilon are scaled linearly with expanded ensemble state.
 *
 * \param [in] a1 Pointer to atom with position to test - this does NOT need to be in the simulation box a priori
 * \param [in] box Simulation box
 */
double squareWellWallZ::energy (const atom *a1, const std::vector < double > &box) {
    std::vector < double > p = a1->pos;
    pbc (p, box);
    double U = 0.0;
    
    double sig = sigma_, eps = eps_;
    if (a1->mState > 0) {
        sig = (sigma_/M_)*a1->mState;
        eps = (eps_/M_)*a1->mState;
    }
    if (a1->mState < 0 || a1->mState > M_-1) {
        throw customException ("mState out of bounds for squareWellWallZ");
    }
    
    // return infinity if out of bounds
    if (p[2] >= ub_ - sig/2.0 || p[2] <= lb_ + sig/2.0) {
        return NUM_INFINITY;
    }
    
    // interaction with top wall
    if (p[2] > ub_ - range_) {
        U += -eps;
    }
    
    // interaction with lower wall
    if (p[2] < lb_ + range_) {
        U += -eps;
    }
    
    return U;
}

/*!
 * Add a hard wall to interact with.
 *
 * \param [in] lb z-Position of the lower wall
 * \param [in] ub z-Position of the upper wall
 * \param [in] sigma Hard-sphere diameter the species this wall interacts with can approach within
 * \param [in] M Total number of expanded ensemble states possible for this atom type (defaults to 1)
 */
void compositeBarrier::addhardWallZ (const double lb, const double ub, const double sigma, const int M) {
    if (sysBarriers_.begin() == sysBarriers_.end()) {
        try {
            sysBarriers_.resize(1);
        } catch (std::bad_alloc &ba) {
            throw customException ("Unable to allocate space for a new barrier");
        }
    } else {
        try {
            sysBarriers_.resize(sysBarriers_.size()+1);
        } catch (std::bad_alloc &ba) {
            throw customException ("Unable to allocate space for a new barrier");
        }
    }
    try {
        sysBarriers_[sysBarriers_.size()-1] = new hardWallZ (lb, ub, sigma, M);
    } catch (customException &ce) {
        throw customException ("Cannot add hardWallZ to composite barrier: "+sstr(ce.what()));
    }
}

/*!
 * Add a square well wall to interact with.
 *
 * \param [in] lb z-Position of the lower wall
 * \param [in] ub z-Position of the upper wall
 * \param [in] sigma Hard-sphere diameter the species this wall interacts with can approach within
 * \param [in] range Distance normal to the wall's surface where there is an interaction
 * \param [in] eps Magnitude of the wall interaction (U = -eps)
 * \param [in] M Total number of expanded ensemble states possible for this atom type (defaults to 1)
 */
void compositeBarrier::addSquareWellWallZ (const double lb, const double ub, const double sigma, const double range, const double eps, const int M) {
    if (sysBarriers_.begin() == sysBarriers_.end()) {
        try {
            sysBarriers_.resize(1);
        } catch (std::bad_alloc &ba) {
            throw customException ("Unable to allocate space for a new barrier");
        }
    } else {
        try {
            sysBarriers_.resize(sysBarriers_.size()+1);
        } catch (std::bad_alloc &ba) {
            throw customException ("Unable to allocate space for a new barrier");
        }
    }
    try {
        sysBarriers_[sysBarriers_.size()-1] = new squareWellWallZ (lb, ub, sigma, range, eps, M);
    } catch (customException &ce) {
        throw customException ("Cannot add squareWellWallZ to composite barrier: "+sstr(ce.what()));
    }
}

/*!
 * Deallocate any system barriers present.
 */
compositeBarrier::~compositeBarrier () {
    if (sysBarriers_.begin() != sysBarriers_.end()) {
        for (unsigned int i = 0; i < sysBarriers_.size(); ++i) {
            delete sysBarriers_[i];
        }
        sysBarriers_.clear();
    }
}

/*!
 * Test if inside ALL the barriers.  Returns false if outside any single one.
 *
 * \param [in] a1 Pointer to atom with position to test - this does NOT need to be in the simulation box a priori
 * \param [in] box Simulation box
 */
bool compositeBarrier::inside (const atom *a1, const std::vector < double > &box) {
    for (std::vector < barrier* >::iterator it = sysBarriers_.begin(); it != sysBarriers_.end(); ++it) {
        if (!(*it)->inside (a1, box)) {
            return false;
        }
    }
    return true;
}

/*!
 * Find the total energy of interaction from ALL the barriers.
 *
 * \param [in] a1 Pointer to atom with position to test - this does NOT need to be in the simulation box a priori
 * \param [in] box Simulation box
 */
double compositeBarrier::energy (const atom *a1, const std::vector < double > &box) {
    double U = 0.0;
    for (std::vector < barrier* >::iterator it = sysBarriers_.begin(); it != sysBarriers_.end(); ++it) {
        double dU = (*it)->energy (a1, box);
        if (dU == NUM_INFINITY) {
            return NUM_INFINITY;
        } else {
            U += dU;
        }
    }
    return U;
}