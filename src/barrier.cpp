#include "barrier.h"

/*!
 * Instantiate a right triangle feature with its base in the xy-plane.  It is raised in the z-direction and extends in the x-direction.
 *
 * \param [in] width Width of triangle's feature
 * \param [in] theta Elevation angle of the feature in radians (0, PI)
 * \param [in] lamW Attractive range ratio relative to hard sphere in contact with the feature (akin to square well), must be >= 1
 * \param [in] eps Attraction strength to feature
 * \param [in] sigma Hard sphere diameter of interaction with the feature
 * \param [in] sep Distance between features
 * \param [in] offset Offset from x = 0 position of the first feature
 * \param [in] box System box size to check the feature (as specified) is periodic in the box
 * \param [in] top If true, feature is on the "top", else is on the bottom (default)
 * \param [in] Number of expanded ensemble states to recognize (default = 1).
 */
rightTriangleZ::rightTriangleZ (const double width, const double theta, const double lamW, const double eps, const double sigma, const double sep, const double offset, const std::vector < double > &box, bool top, const int M) {
    if (sep < 0.0) {
        throw customException ("rightTriangle sep is out of bounds");
        exit(SYS_FAILURE);
    }
    // Check the feature fits periodically in the box
    if ( fmod(box[0]/(sep+width), 1) != 0.0 ) {
        throw customException ("rightTriangle width+separation is not commensurate with the box size");
        exit(SYS_FAILURE);
    }
    if (theta <= 0.0 || theta >= PI) {
        throw customException ("rightTriangle elevation angle is out of bounds");
        exit(SYS_FAILURE);
    }
    if (sigma <= 0.0) {
        throw customException ("rightTriangle sigma is out of bounds");
        exit(SYS_FAILURE);
    }
    if (lamW < 1.0) {
        throw customException ("rightTriangle lamW is out of bounds");
        exit(SYS_FAILURE);
    }
    if (eps < 1.0) {
        throw customException ("rightTriangle eps is out of bounds");
        exit(SYS_FAILURE);
    }
    if (M < 1) {
        throw customException ("rightTriangle M value is out of bounds");
        exit(SYS_FAILURE);
    }
    
    // store variables internally and do conversion to conventions in which derivation was done
    lam_ = width;
    theta_ = theta;
    sigma_ = sigma;
    lamWall_ = lamW*sigma - sigma;
    eps_ = eps;
    lamP_ = sep;
    xOffset_ = offset;
    top_ = top;
    box_.resize(3);
    box_ = box;
    M_ = M;
    
    // precompute points and trigonometry
    cosTheta_ = cos(theta_);
    sinTheta_ = sin(theta_);
    a_ = lam_*cos(theta_)*sin(theta_);
    b_ = lam_*sin(theta_)*cos(PI/2.0-theta_);
    c_ = lam_*cos(theta_)*cos(theta_);
    
    m_ = c_ - (lamWall_ + sigma_)*sinTheta_;
    n_ = a_ + (lamWall_ + sigma_)*cosTheta_;
    u_ = c_ + sigma_/2.0*(cosTheta_ - sinTheta_);
    v_ = a_ + sigma_/2.0*(cosTheta_ + sinTheta_);
    p_ = c_ + (lamWall_ + sigma_)*cosTheta_;
    q_ = a_ + (lamWall_ + sigma_)*sinTheta_;
    j_ = c_ + b_ + sigma_/2.0*cosTheta_;
    k_ = sigma_/2.0*sinTheta_;
    
    s_ = (q_ - n_)/(p_ - m_);
    
    lbounds_.resize(7);
    ubounds_.resize(7);
    
    lbounds_[0] = -(lamWall_+sigma_)*sinTheta_;
    lbounds_[1] = -sigma_/2.0*sinTheta_;
    lbounds_[2] = 0.0;
    lbounds_[3] = m_;
    lbounds_[4] = u_;
    lbounds_[5] = p_;
    lbounds_[6] = c_ + b_;
    lbounds_[7] = j_;
    
    ubounds_[0] = -sigma_/2.0*sinTheta_;
    ubounds_[1] = 0.0;
    ubounds_[2] = m_;
    ubounds_[3] = u_;
    ubounds_[4] = p_;
    ubounds_[5] = b_ + c_;
    ubounds_[6] = j_;
    ubounds_[7] = b_ + c_ + (lamWall_ + sigma_)*cosTheta_;
}

/*!
 * Return the x-position of an atom that is within the first period of this feature.
 *
 * \param [in] atom Pointer to atom to examine
 */
double rightTriangleZ::getRange_ (double x) {
    double L = lam_ + lamP_;
    while (x > lam_ + (lamWall_ + sigma_)*cosTheta_ + lamP_) {
        x -= L;
    }
    while (x < -(lamWall_ + sigma_)*sinTheta_) {
        x += L;
    }
    return x;
}

/*!
 * Get the energy of a position assuming it is located within window 0.
 *
 * \param [in] atom Pointer to atom to examine
 */
double rightTriangleZ::getU0_ (const double x, const double z) {
    double ub = (a_/c_)*(x + (lamWall_ + sigma_)*sinTheta_) + (lamWall_ + sigma_)*cosTheta_;
    double lb = (-a_/b_)*(x + (lamWall_ + sigma_)*sinTheta_) + (lamWall_ + sigma_)*cosTheta_;
    if (z < lb) {
        return 0.0;
    } else if (z <= ub) {
        return -eps_;
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 1.
 *
 * \param [in] atom Pointer to atom to examine
 */
double rightTriangleZ::getU1_ (const double x, const double z) {
    double ub = (a_/c_)*(x + (lamWall_ + sigma_)*sinTheta_) + (lamWall_ + sigma_)*cosTheta_;
    double lb1 = (a_/c_)*(x + sigma_/2.0*sinTheta_) + sigma_/2.0*cosTheta_;
    double lb2 = (-a_/b_)*(x + (lamWall_ + sigma_)*sinTheta_) + (lamWall_ + sigma_)*cosTheta_;
    if (z < lb2) {
        return 0.0;
    } else if (z < lb1) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_;
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 2.
 *
 * \param [in] atom Pointer to atom to examine
 */
double rightTriangleZ::getU2_ (const double x, const double z) {
    double ub = (a_/c_)*(x + (lamWall_ + sigma_)*sinTheta_) + (lamWall_ + sigma_)*cosTheta_;
    double lb = (a_/c_)*(x + sigma_/2.0*sinTheta_) + sigma_/2.0*cosTheta_;
    if (z < lb) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_;
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 3.
 *
 * \param [in] atom Pointer to atom to examine
 */
double rightTriangleZ::getU3_ (const double x, const double z) {
    double ub = s_*(x - m_) + n_;
    double lb = (a_/c_)*(x + sigma_/2.0*sinTheta_) + sigma_/2.0*cosTheta_;
    if (z < lb) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_;
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 4.
 *
 * \param [in] atom Pointer to atom to examine
 */
double rightTriangleZ::getU4_ (const double x, const double z) {
    double ub = s_*(x - m_) + n_;
    double lb = (-a_/b_)*(x - u_) + v_;
    if (z < lb) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_;
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 5.
 *
 * \param [in] atom Pointer to atom to examine
 */
double rightTriangleZ::getU5_ (const double x, const double z) {
    double ub = (-a_/b_)*(x - p_) + q_;
    double lb = (-a_/b_)*(x - u_) + v_;
    if (z < lb) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_;
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 6.
 *
 * \param [in] atom Pointer to atom to examine
 */
double rightTriangleZ::getU6_ (const double x, const double z) {
    double ub = (-a_/b_)*(x - p_) + q_;
    double lb1 = (-a_/b_)*(x - u_) + v_;
    double lb2 = (a_/c_)*(x - j_) + k_;
    if (z < lb2) {
        return 0.0;
    } else if (z < lb1) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_;
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 7.
 *
 * \param [in] atom Pointer to atom to examine
 */
double rightTriangleZ::getU7_ (const double x, const double z) {
    double ub = (-a_/b_)*(x - p_) + q_;
    double lb = (a_/c_)*(x - j_) + k_;
    if (z < lb) {
        return 0.0;
    } else if (z <= ub) {
        return -eps_;
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 8.
 *
 * \param [in] atom Pointer to atom to examine
 * \param [in] box System box size. Not actually used, but will be checked that it is identical to value at class instantiation.
 */
double rightTriangleZ::energy (const atom *a1, const std::vector < double > &box) {
    for (unsigned int i = 0; i < box_.size(); ++i) {
        if (box[i] != box_[i]) {
            throw customException ("System box size has changed from when rightTriangleZ was instantiated");
            exit(SYS_FAILURE);
        }
    }
    
    double z = a1->pos[2], x = a1->pos[0];
    if (top_) {
        z = box_[2] - a1->pos[2];
    }
    
    // early check
    if (z > std::max(n_, q_)) { // n > q if theta < PI/4
        return 0.0;
    }
    
    x = getRange_ (x);
    double U = 0.0;
    
    if (x < ubounds_[0] and x >= lbounds_[0]) {
        U += getU0_ (x, z);
    } else if (x < ubounds_[1] and x >= lbounds_[1]) {
        U += getU1_ (x, z);
    } else if (x < ubounds_[2] and x >= lbounds_[2]) {
        U += getU2_ (x, z);
    } else if (x < ubounds_[3] and x >= lbounds_[3]) {
        U += getU3_ (x, z);
    } else if (x < ubounds_[4] and x >= lbounds_[4]) {
        U += getU4_ (x, z);
    } else if (x < ubounds_[5] and x >= lbounds_[5]) {
        U += getU5_ (x, z);
    } else if (x < ubounds_[6] and x >= lbounds_[6]) {
        U += getU6_ (x, z);
    } else if (x < ubounds_[7] and x >= lbounds_[7]) {
        U += getU7_ (x, z);
    } else {
        U += 0.0;
    }
    
    // Check interactions with periodic features
    if (U < NUM_INFINITY) {
        double xn = 0.0;
        bool periodic = false;
        if (x < 0+(lamWall_ + sigma_)*cosTheta_ - lamP_) {
            xn = x+ (lamP_ + lam_);
            periodic = true;
        } else if (x > lam_ - (lamWall_ + sigma_)*sinTheta_ + lamP_) {
            xn = x - (lamP_ + lam_);
            periodic = true;
        }
            
        if (periodic) {
            if (xn < ubounds_[0] and xn >= lbounds_[0]) {
                U += getU0_ (xn, z);
            } else if (xn < ubounds_[1] and xn >= lbounds_[1]) {
                U += getU1_ (xn, z);
            } else if (xn < ubounds_[2] and xn >= lbounds_[2]) {
                U += getU2_ (xn, z);
            } else if (xn < ubounds_[3] and xn >= lbounds_[3]) {
                U += getU3_ (xn, z);
            } else if (xn < ubounds_[4] and xn >= lbounds_[4]) {
                U += getU4_ (xn, z);
            } else if (xn < ubounds_[5] and xn >= lbounds_[5]) {
                U += getU5_ (xn, z);
            } else if (xn < ubounds_[6] and xn >= lbounds_[6]) {
                U += getU6_ (xn, z);
            } else if (xn < ubounds_[7] and xn >= lbounds_[7]) {
                U += getU7_ (xn, z);
            } else {
                U += 0.0;
            }
        }
    }
    
    return U;
}

/*!
 * Check if an atom is overlapping the feature.  The term "inside" can be a bit confusing here.
 * This function returns true if an atom does NOT overlap the feature (have infinite interaction energy),
 * which might be considered as being "outside" the feature depending on how you look at it.
 * However, to be consistent with the expected behavior of this virtual function, this is how
 * this function must behave.
 *
 * \param [in] a1 Atom whose position to test.
 * \param [in] box System box size. Not actually used, but will be checked that it is identical to value at class instantiation.
 */
bool rightTriangleZ::inside (const atom *a1, const std::vector < double > &box) {
    double U = NUM_INFINITY;
    try {
        U = energy (a1, box);
    } catch (customException &ce) {
        throw customException ("Unable to test if inside rightTriangleZ : "+sstr(ce.what()));
        exit (SYS_FAILURE);
    }
    
    if (U < NUM_INFINITY) {
        return true;
    } else {
        return false;
    }
}

// also add to compositeBarrier as well
// fix rightTriangleZ to account for possibility that an atom has different M (must change sigma and therefore init must be adjusted to use a separate function to initialize variables - or else variables must be looked up on the fly based on M)


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
    
    eps_ = eps;
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
void compositeBarrier::addHardWallZ (const double lb, const double ub, const double sigma, const int M) {
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
 * Test if inside ALL the barriers.  Returns false if outside any single one, but defaults to true (infinitely far away walls/barriers).
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