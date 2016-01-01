#include "barrier.h"

/*!
 * Instantiate a right triangle feature with its base in the xy-plane.  It is raised in the z-direction and extends in the x-direction. Expanded ensembles scale the range of interaction as well as the magnitude.
 *
 * \param [in] width Width of triangle's feature
 * \param [in] theta Elevation angle of the feature in radians (0, PI)
 * \param [in] lamW Attractive range ratio relative to hard sphere in contact with the feature (akin to a square well potential's lambda), must be >= 1
 * \param [in] eps Attraction strength to feature
 * \param [in] sigma Hard sphere diameter of interaction with the feature
 * \param [in] sep Distance between features
 * \param [in] offset Offset from x = 0 position of the first feature.  A positive offset shifts the feature in the +x direction.
 * \param [in] box System box size to check the feature (as specified) is periodic in the box
 * \param [in] zbase Z-coordinate of XY plane that defines the base of the feature.  To avoid periodicity effects be sure it is > 0 and less than Lz, but this depends on other interactions so it cannot be checked automatically here.
 * \param [in] top If true, feature is on the "top", else is on the bottom (default)
 * \param [in] Number of expanded ensemble states to recognize (default = 1)
 */
rightTriangleZ::rightTriangleZ (const double width, const double theta, const double lamW, const double eps, const double sigma, const double sep, const double offset, const std::vector < double > &box, const double zbase, bool top, const int M) {
    if (sep < 0.0) {
        throw customException ("rightTriangle sep is out of bounds");
    }
    // Check the feature fits periodically in the box
    if ( fmod(box[0]/(sep+width), 1) != 0.0 ) {
        throw customException ("rightTriangle width+separation is not commensurate with the box size");
    }
    if ((int)round(box[0]/(sep+width)) == 1) {
        throw customException ("Only 1 rightTriangle feature in box - there could exist (self) periodicity problems");
    }
    if (theta <= 0.0 || theta >= PI/2.0) {
        throw customException ("rightTriangle elevation angle is out of bounds");
    }
    if (sigma <= 0.0) {
        throw customException ("rightTriangle sigma is out of bounds");
    }
    if (lamW < 1.0) {
        throw customException ("rightTriangle lamW is out of bounds");
    }
    if (eps < 0.0) {
        throw customException ("rightTriangle eps is out of bounds");
    }
    if (M < 1) {
        throw customException ("rightTriangle M value is out of bounds");
    }
    if (zbase < 0.0 || zbase > box[2]) {
        throw customException ("rightTriangle zbase value is out of bounds");
    }
    
    // store variables internally and do conversion to conventions in which derivation was done
    zbase_ = zbase;
    lam_ = width;
    lamP_ = sep;
    theta_ = theta;
    xOffset_ = offset;
    top_ = top;
    box_.resize(3);
    box_ = box;
    M_ = M;
    cosTheta_ = cos(theta_);
    sinTheta_ = sin(theta_);
    
    // precompute points and trigonometry
    a_ = lam_*cos(theta_)*sin(theta_);
    b_ = lam_*sin(theta_)*cos(PI/2.0-theta_);
    c_ = lam_*cos(theta_)*cos(theta_);
    
    try {
        sigma_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        lamWall_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        eps_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        m_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        n_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        u_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        v_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        p_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        q_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        j_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        k_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        s_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        lbounds_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    for (unsigned int i = 0; i < M_; ++i) {
        try {
            lbounds_[i].resize(8);
        } catch (std::bad_alloc &ba) {
            throw customException ("Out of memory");
        }
    }
    try {
        ubounds_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    for (unsigned int i = 0; i < M_; ++i) {
        try {
            ubounds_[i].resize(8);
        } catch (std::bad_alloc &ba) {
            throw customException ("Out of memory");
        }
    }
    
    for (unsigned int i = 0; i < M_; ++i) {
        if (i == 0) {
            sigma_[i] = sigma;
            eps_[i] = eps;
        } else {
            sigma_[i] = sigma_[0]/M_*i;
            eps_[i] = eps_[0]/M_*i;
        }
        
        lamWall_[i] = lamW*sigma_[i] - sigma_[i];
        m_[i] = c_ - (lamWall_[i] + sigma_[i])*sinTheta_;
        n_[i] = a_ + (lamWall_[i] + sigma_[i])*cosTheta_;
        u_[i] = c_ + sigma_[i]/2.0*(cosTheta_ - sinTheta_);
        v_[i] = a_ + sigma_[i]/2.0*(cosTheta_ + sinTheta_);
        p_[i] = c_ + (lamWall_[i] + sigma_[i])*cosTheta_;
        q_[i] = a_ + (lamWall_[i] + sigma_[i])*sinTheta_;
        j_[i] = c_ + b_ + sigma_[i]/2.0*cosTheta_;
        k_[i] = sigma_[i]/2.0*sinTheta_;
        s_[i] = (q_[i] - n_[i])/(p_[i] - m_[i]);
        
        lbounds_[i][0] = -(lamWall_[i]+sigma_[i])*sinTheta_;
        lbounds_[i][1] = -sigma_[i]/2.0*sinTheta_;
        lbounds_[i][2] = 0.0;
        lbounds_[i][3] = m_[i];
        lbounds_[i][4] = u_[i];
        lbounds_[i][5] = p_[i];
        lbounds_[i][6] = c_ + b_;
        lbounds_[i][7] = j_[i];
        
        ubounds_[i][0] = -sigma_[i]/2.0*sinTheta_;
        ubounds_[i][1] = 0.0;
        ubounds_[i][2] = m_[i];
        ubounds_[i][3] = u_[i];
        ubounds_[i][4] = p_[i];
        ubounds_[i][5] = b_ + c_;
        ubounds_[i][6] = j_[i];
        ubounds_[i][7] = b_ + c_ + (lamWall_[i] + sigma_[i])*cosTheta_;
    }
}

/*!
 * Return the x-position of an atom that is within the first period of this feature. This function takes into account the x-offset of the feature automatically for all other routines.
 *
 * \param [in] x X position of an atom accounting for box pbc, but not the feature's periodicity
 * \param [in] mv Expanded ensemble state of the atom
 */
double rightTriangleZ::getRange_ (double x, const int mv) {
    double L = lam_ + lamP_;
    
    // x is currently in the "box frame" so account for xOffset
    while (x > lam_ + (lamWall_[mv] + sigma_[mv])*cosTheta_ + lamP_ + xOffset_) {
        x -= L;
    }
    while (x < -(lamWall_[mv] + sigma_[mv])*sinTheta_ + xOffset_) {
        x += L;
    }
    
    // shift back to "origin" of feature
    x -= xOffset_;
    
    // in the case of having only one feature, these loops could place the particle out of the box
    // however, this is not a problem since that convention is not assumed during this calculation
    
    return x;
}

/*!
 * Get the energy of a position assuming it is located within window 0.
 *
 * \param [in] x X position, wrapped to periodicity of feature
 * \param [in] z Z position, shifted to account for feature's zbase
 * \param [in] mv Expanded ensemble state of the atom
 */
double rightTriangleZ::getU0_ (const double x, const double z, const int mv) {
    double ub = (a_/c_)*(x + (lamWall_[mv] + sigma_[mv])*sinTheta_) + (lamWall_[mv] + sigma_[mv])*cosTheta_;
    double lb = (-a_/b_)*(x + (lamWall_[mv] + sigma_[mv])*sinTheta_) + (lamWall_[mv] + sigma_[mv])*cosTheta_;
    if (z < lb) {
        return 0.0;
    } else if (z <= ub) {
        return -eps_[mv];
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 1.
 *
 * \param [in] x X position, wrapped to periodicity of feature
 * \param [in] z Z position, shifted to account for feature's zbase
 * \param [in] mv Expanded ensemble state of the atom
 */
double rightTriangleZ::getU1_ (const double x, const double z, const int mv) {
    double ub = (a_/c_)*(x + (lamWall_[mv] + sigma_[mv])*sinTheta_) + (lamWall_[mv] + sigma_[mv])*cosTheta_;
    double lb1 = (a_/c_)*(x + sigma_[mv]/2.0*sinTheta_) + sigma_[mv]/2.0*cosTheta_;
    double lb2 = (-a_/b_)*(x + (lamWall_[mv] + sigma_[mv])*sinTheta_) + (lamWall_[mv] + sigma_[mv])*cosTheta_;
    if (z < lb2) {
        return 0.0;
    } else if (z < lb1) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_[mv];
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 2.
 *
 * \param [in] x X position, wrapped to periodicity of feature
 * \param [in] z Z position, shifted to account for feature's zbase
 * \param [in] mv Expanded ensemble state of the atom
 */
double rightTriangleZ::getU2_ (const double x, const double z, const int mv) {
    double ub = (a_/c_)*(x + (lamWall_[mv] + sigma_[mv])*sinTheta_) + (lamWall_[mv] + sigma_[mv])*cosTheta_;
    double lb = (a_/c_)*(x + sigma_[mv]/2.0*sinTheta_) + sigma_[mv]/2.0*cosTheta_;
    if (z < lb) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_[mv];
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 3.
 *
 * \param [in] x X position, wrapped to periodicity of feature
 * \param [in] z Z position, shifted to account for feature's zbase
 * \param [in] mv Expanded ensemble state of the atom
 */
double rightTriangleZ::getU3_ (const double x, const double z, const int mv) {
    double ub = s_[mv]*(x - m_[mv]) + n_[mv];
    double lb = (a_/c_)*(x + sigma_[mv]/2.0*sinTheta_) + sigma_[mv]/2.0*cosTheta_;
    if (z < lb) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_[mv];
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 4.
 *
 * \param [in] x X position, wrapped to periodicity of feature
 * \param [in] z Z position, shifted to account for feature's zbase
 * \param [in] mv Expanded ensemble state of the atom
 */
double rightTriangleZ::getU4_ (const double x, const double z, const int mv) {
    double ub = s_[mv]*(x - m_[mv]) + n_[mv];
    double lb = (-a_/b_)*(x - u_[mv]) + v_[mv];
    if (z < lb) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_[mv];
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 5.
 *
 * \param [in] x X position, wrapped to periodicity of feature
 * \param [in] z Z position, shifted to account for feature's zbase
 * \param [in] mv Expanded ensemble state of the atom
 */
double rightTriangleZ::getU5_ (const double x, const double z, const int mv) {
    double ub = (-a_/b_)*(x - p_[mv]) + q_[mv];
    double lb = (-a_/b_)*(x - u_[mv]) + v_[mv];
    if (z < lb) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_[mv];
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 6.
 *
 * \param [in] x X position, wrapped to periodicity of feature
 * \param [in] z Z position, shifted to account for feature's zbase
 * \param [in] mv Expanded ensemble state of the atom
 */
double rightTriangleZ::getU6_ (const double x, const double z, const int mv) {
    double ub = (-a_/b_)*(x - p_[mv]) + q_[mv];
    double lb1 = (-a_/b_)*(x - u_[mv]) + v_[mv];
    double lb2 = (a_/c_)*(x - j_[mv]) + k_[mv];
    if (z < lb2) {
        return 0.0;
    } else if (z < lb1) {
        return NUM_INFINITY;
    } else if (z <= ub) {
        return -eps_[mv];
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 7.
 *
 * \param [in] x X position, wrapped to periodicity of feature
 * \param [in] z Z position, shifted to account for feature's zbase
 * \param [in] mv Expanded ensemble state of the atom
 */
double rightTriangleZ::getU7_ (const double x, const double z, const int mv) {
    double ub = (-a_/b_)*(x - p_[mv]) + q_[mv];
    double lb = (a_/c_)*(x - j_[mv]) + k_[mv];
    if (z < lb) {
        return 0.0;
    } else if (z <= ub) {
        return -eps_[mv];
    } else {
        return 0.0;
    }
}

/*!
 * Get the energy of a position assuming it is located within window 8.
 *
 * \param [in] atom Pointer to atom to examine
 * \param [in] box System box size. Will be checked that it is identical to value at class instantiation.
 */
double rightTriangleZ::energy (const atom *a1, const std::vector < double > &box) {
    for (unsigned int i = 0; i < box_.size(); ++i) {
        if (box[i] != box_[i]) {
            throw customException ("System box size has changed from when rightTriangleZ was instantiated");
        }
    }

    const int mv = a1->mState;
    std::vector < double > p = a1->pos;
    pbc (p, box);
    
    if (mv < 0 || mv > M_-1) {
        throw customException ("mState out of bounds for rightTriangleZ");
    }
    
    // Shift z and x into origin frame of feature
    double z = p[2], x = p[0];
    if (top_) {
        z = zbase_ - p[2];
    } else {
        z = p[2] - zbase_;
    }
    x = getRange_ (x, mv);
    
    // early check
    if (z > std::max(n_[mv], q_[mv])) { // n > q if theta < PI/4
        return 0.0;
    }

    double U = 0.0;
    
    if (x < ubounds_[mv][0] and x >= lbounds_[mv][0]) {
        U = getU0_ (x, z, mv);
    } else if (x < ubounds_[mv][1] and x >= lbounds_[mv][1]) {
        U = getU1_ (x, z, mv);
    } else if (x < ubounds_[mv][2] and x >= lbounds_[mv][2]) {
        U = getU2_ (x, z, mv);
    } else if (x < ubounds_[mv][3] and x >= lbounds_[mv][3]) {
        U = getU3_ (x, z, mv);
    } else if (x < ubounds_[mv][4] and x >= lbounds_[mv][4]) {
        U = getU4_ (x, z, mv);
    } else if (x < ubounds_[mv][5] and x >= lbounds_[mv][5]) {
        U = getU5_ (x, z, mv);
    } else if (x < ubounds_[mv][6] and x >= lbounds_[mv][6]) {
        U = getU6_ (x, z, mv);
    } else if (x < ubounds_[mv][7] and x >= lbounds_[mv][7]) {
        U = getU7_ (x, z, mv);
    } else {
        U = 0.0;
    }
    
    // Check interactions with periodic features
    if (U < NUM_INFINITY) {
        double xn = 0.0;
        bool periodic = false;
        if (x < 0 + (lamWall_[mv] + sigma_[mv])*cosTheta_ - lamP_) {
            xn = x + (lamP_ + lam_);
            periodic = true;
        } else if (x > lam_ - (lamWall_[mv] + sigma_[mv])*sinTheta_ + lamP_) {
            xn = x - (lamP_ + lam_);
            periodic = true;
        }
            
        if (periodic) {
            double dU = 0.0;
            if (xn < ubounds_[mv][0] and xn >= lbounds_[mv][0]) {
                dU = getU0_ (xn, z, mv);
            } else if (xn < ubounds_[mv][1] and xn >= lbounds_[mv][1]) {
                dU = getU1_ (xn, z, mv);
            } else if (xn < ubounds_[mv][2] and xn >= lbounds_[mv][2]) {
                dU = getU2_ (xn, z, mv);
            } else if (xn < ubounds_[mv][3] and xn >= lbounds_[mv][3]) {
                dU = getU3_ (xn, z, mv);
            } else if (xn < ubounds_[mv][4] and xn >= lbounds_[mv][4]) {
                dU = getU4_ (xn, z, mv);
            } else if (xn < ubounds_[mv][5] and xn >= lbounds_[mv][5]) {
                dU = getU5_ (xn, z, mv);
            } else if (xn < ubounds_[mv][6] and xn >= lbounds_[mv][6]) {
                dU = getU6_ (xn, z, mv);
            } else if (xn < ubounds_[mv][7] and xn >= lbounds_[mv][7]) {
                dU = getU7_ (xn, z, mv);
            } else {
                dU = 0.0;
            }
            
            if (dU < NUM_INFINITY) {
                U += dU;
            } else {
                return NUM_INFINITY;
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
        U = energy (a1, box); // takes care of pbc internally
    } catch (customException &ce) {
        throw customException ("Unable to test if inside rightTriangleZ : "+sstr(ce.what()));
        //exit (SYS_FAILURE);
    }
    
    if (U < NUM_INFINITY) {
        return true;
    } else {
        return false;
    }
}

/*!
 * Instantiate a hard wall with boundaries in the +/- z direction. Expanded ensembles scale the range of interaction via the sigma parameter.
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
 * Instantiate a square well wall with boundaries in the +/- z direction. Expanded ensembles primarily scale the the magnitude of interaction.  The repulsive boundary scales with sigma at the boundary, but the attractive cutoff remains fixed relative to the boundary.
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
 * Add a rightTriangleZ feature to interact with.
 *
 * \param [in] width Width of triangle's feature
 * \param [in] theta Elevation angle of the feature in radians (0, PI)
 * \param [in] lamW Attractive range ratio relative to hard sphere in contact with the feature (akin to square well), must be >= 1
 * \param [in] eps Attraction strength to feature
 * \param [in] sigma Hard sphere diameter of interaction with the feature
 * \param [in] sep Distance between features
 * \param [in] offset Offset from x = 0 position of the first feature
 * \param [in] box System box size to check the feature (as specified) is periodic in the box
 * \param [in] zbase Z-coordinate of XY plane that defines the base of the feature.  To avoid periodicity effects be sure it is > 0 and less than Lz, but this depends on other interactions so it cannot be checked automatically here.
 * \param [in] top If true, feature is on the "top", else is on the bottom (default)
 * \param [in] Number of expanded ensemble states to recognize (default = 1)
 */
void compositeBarrier::addRightTriangleZ (const double width, const double theta, const double lamW, const double eps, const double sigma, const double sep, const double offset, const std::vector < double > &box, const double zbase, bool top, const int M) {
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
        sysBarriers_[sysBarriers_.size()-1] = new rightTriangleZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        throw customException ("Cannot add rightTriangleZ to composite barrier: "+sstr(ce.what()));
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
        if (dU < NUM_INFINITY) {
            U += dU;
        } else {
            return NUM_INFINITY;
        }
    }
    return U;
}