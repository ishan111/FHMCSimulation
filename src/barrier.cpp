#include "barrier.h"

/*!
 * Instantiate a hard wall with boundaries in the +/- z direction.
 *
 * \param [in] lb z-Position of the lower wall
 * \param [in] ub z-Position of the upper wall
 * \param [in] sigma Hard-sphere diameter the species this wall interacts with can approach within
 */
hardWallZ::hardWallZ (const double lb, const double ub, const double sigma) {
    if (lb >= ub) {
        throw customException ("hardWallZ must have lower bound < upper bound");
    }
    if (sigma < 0) {
        throw customException ("hardWallZ must have sigma >= 0");
    }
    
    sigma_ = sigma;
    ub_ = ub;
    lb_ = lb;
}

/*!
 * Return whether or not a point falls between the walls (subject to a hard-sphere exclusion radius).
 *
 * \param [in] point (x, y, z) Point to test - this does NOT need to be in the simulation box a priori
 * \param [in] box Simulation box
 */
bool hardWallZ::inside (const std::vector < double > &point, const std::vector < double > &box) {
    std::vector < double > p = point;
    pbc (p, box);
    if (p[2] >= ub_ - sigma_/2.0 || p[2] <= lb_ + sigma_/2.0) {
        return false;
    } else {
        return true;
    }
}

/*!
 * Interaction energy with the wall.
 *
 * \param [in] point (x, y, z) Point to test - this does NOT need to be in the simulation box a priori
 * \param [in] box Simulation box
 */
double hardWallZ::energy (const std::vector < double > &point, const std::vector < double > &box) {
    std::vector < double > p = point;
    pbc (p, box);
    if (p[2] >= ub_ - sigma_/2.0 || p[2] <= lb_ + sigma_/2.0) {
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
 */
squareWellWallZ::squareWellWallZ (const double lb, const double ub, const double sigma, const double range, const double eps) {
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
    
    sigma_ = sigma;
    range_ = range;
    ub_ = ub;
    lb_ = lb;
}

/*!
 * Return whether or not a point falls between the walls (subject to a hard-sphere exclusion radius).
 *
 * \param [in] point (x, y, z) Point to test - this does NOT need to be in the simulation box a priori
 * \param [in] box Simulation box
 */
bool squareWellWallZ::inside (const std::vector < double > &point, const std::vector < double > &box) {
    std::vector < double > p = point;
    pbc (p, box);
    if (p[2] >= ub_ - sigma_/2.0 || p[2] <= lb_ + sigma_/2.0) {
        return false;
    } else {
        return true;
    }
}

/*!
 * Interaction energy with the wall.
 *
 * \param [in] point (x, y, z) Point to test - this does NOT need to be in the simulation box a priori
 * \param [in] box Simulation box
 */
double squareWellWallZ::energy (const std::vector < double > &point, const std::vector < double > &box) {
    std::vector < double > p = point;
    pbc (p, box);
    double U = 0.0;
    
    // return infinity if out of bounds
    if (p[2] >= ub_ - sigma_/2.0 || p[2] <= lb_ + sigma_/2.0) {
        return NUM_INFINITY;
    }
    
    // interaction with top wall
    if (p[2] > ub_ - range_) {
        U += -eps_;
    }
    
    // interaction with lower wall
    if (p[2] < lb_ + range_) {
        U += -eps_;
    }
    
    return U;
}