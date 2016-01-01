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
 * \param [in] M Number of expanded ensemble states to recognize (default = 1)
 */
rightTriangleXZ::rightTriangleXZ (const double width, const double theta, const double lamW, const double eps, const double sigma, const double sep, const double offset, const std::vector < double > &box, const double zbase, bool top, const int M) {
    if (sep < 0.0) {
        throw customException ("rightTriangle sep is out of bounds");
    }
    // Check the feature fits periodically in the box
    if ( fmod(box[0]/(sep+width), 1) != 0.0 ) {
        throw customException ("rightTriangle width+separation is not commensurate with the box size");
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
    if (offset >= box[0] || offset < 0) {
        throw customException ("rightTriangle offset value is out of bounds");
    }
    
    // store variables internally and do conversion to conventions in which derivation was done
    zbase_ = zbase;
    width_ = width;
    sep_ = sep;
    theta_ = theta;
    xOffset_ = offset;
    top_ = top;
    M_ = M;
    x_max_image_ = (int) round(box[0]/(sep+width)) - 1;
    
    // precompute points and trigonometry
    a_ = width_*cos(theta_)*sin(theta_);
    b_ = width_*sin(theta_)*cos(PI/2.0-theta_);
    c_ = width_*cos(theta_)*cos(theta_);

    try {
        sigma_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        eps_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    try {
        box_.resize(3);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    box_ = box;
    try {
        ub_check_.resize(M_, -NUM_INFINITY);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    
    try {
        ub_seg_x_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    for (unsigned int i = 0; i < M_; ++i) {
        try {
            ub_seg_x_[i].resize(4);
        } catch (std::bad_alloc &ba) {
            throw customException ("Out of memory");
        }
    }
    try {
        lb_seg_x_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    for (unsigned int i = 0; i < M_; ++i) {
        try {
            lb_seg_x_[i].resize(4);
        } catch (std::bad_alloc &ba) {
            throw customException ("Out of memory");
        }
    }
    try {
        ub_seg_z_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    for (unsigned int i = 0; i < M_; ++i) {
        try {
            ub_seg_z_[i].resize(4);
        } catch (std::bad_alloc &ba) {
            throw customException ("Out of memory");
        }
    }
    try {
        lb_seg_z_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    for (unsigned int i = 0; i < M_; ++i) {
        try {
            lb_seg_z_[i].resize(4);
        } catch (std::bad_alloc &ba) {
            throw customException ("Out of memory");
        }
    }
    try {
        ub_slope_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    for (unsigned int i = 0; i < M_; ++i) {
        try {
            ub_slope_[i].resize(5);
        } catch (std::bad_alloc &ba) {
            throw customException ("Out of memory");
        }
    }
    try {
        lb_slope_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    for (unsigned int i = 0; i < M_; ++i) {
        try {
            lb_slope_[i].resize(5);
        } catch (std::bad_alloc &ba) {
            throw customException ("Out of memory");
        }
    }
    try {
        ub_int_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    for (unsigned int i = 0; i < M_; ++i) {
        try {
            ub_int_[i].resize(5);
        } catch (std::bad_alloc &ba) {
            throw customException ("Out of memory");
        }
    }
    try {
        lb_int_.resize(M_);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory");
    }
    for (unsigned int i = 0; i < M_; ++i) {
        try {
            lb_int_[i].resize(5);
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
        
        ub_seg_x_[i][0] = -lamW_*sigma_[i]*sin(theta_);
        ub_seg_x_[i][1] = c_ - lamW_*sigma_[i]*sin(theta_);
        ub_seg_x_[i][2] = c_ + lamW_*sigma_[i]*cos(theta_);
        ub_seg_x_[i][3] = width_ + lamW_*sigma_[i]*cos(theta_);
        
        lb_seg_x_[i][0] = -sigma_[i]/2.0*sin(theta_);
        lb_seg_x_[i][1] = c_ - sigma_[i]/2.0*sin(theta_);
        lb_seg_x_[i][2] = c_ + sigma_[i]/2.0*cos(theta_);
        lb_seg_x_[i][3] = width_ + sigma_[i]/2.0*cos(theta_);
        
        ub_seg_z_[i][0] = lamW_*sigma_[i]*cos(theta_);
        ub_seg_z_[i][1] = a_ + lamW_*sigma_[i]*cos(theta_);
        ub_seg_z_[i][2] = a_ + lamW_*sigma_[i]*sin(theta_);
        ub_seg_z_[i][3] = lamW_*sigma_[i]*sin(theta_);
        
        lb_seg_z_[i][0] = sigma_[i]/2.0*cos(theta_);
        lb_seg_z_[i][1] = a_ + sigma_[i]/2.0*cos(theta_);
        lb_seg_z_[i][2] = a_ + sigma_[i]/2.0*sin(theta_);
        lb_seg_z_[i][3] = sigma_[i]/2.0*sin(theta_);
        
        ub_check_[i] = *std::max_element(ub_seg_z_[i].begin(), ub_seg_z_[i].end());
        
        ub_slope_[i][0] = 0.0;
        ub_slope_[i][1] = (ub_seg_z_[i][1] - ub_seg_z_[i][0])/(ub_seg_x_[i][1] - ub_seg_x_[i][0]);
        ub_slope_[i][2] = (ub_seg_z_[i][2] - ub_seg_z_[i][1])/(ub_seg_x_[i][2] - ub_seg_x_[i][1]);
        ub_slope_[i][3] = (ub_seg_z_[i][3] - ub_seg_z_[i][2])/(ub_seg_x_[i][3] - ub_seg_x_[i][2]);
        ub_slope_[i][4] = 0.0;
        
        ub_int_[i][0] = 0.0;
        ub_int_[i][1] = ub_seg_z_[i][0];
        ub_int_[i][2] = ub_seg_z_[i][1];
        ub_int_[i][3] = ub_seg_z_[i][2];
        ub_int_[i][4] = 0.0;
        
        lb_slope_[i][0] = 0.0;
        lb_slope_[i][1] = (lb_seg_z_[i][1] - lb_seg_z_[i][0])/(lb_seg_x_[i][1] - lb_seg_x_[i][0]);
        lb_slope_[i][2] = (lb_seg_z_[i][2] - lb_seg_z_[i][1])/(lb_seg_x_[i][2] - lb_seg_x_[i][1]);
        lb_slope_[i][3] = (lb_seg_z_[i][3] - lb_seg_z_[i][2])/(lb_seg_x_[i][3] - lb_seg_x_[i][2]);
        lb_slope_[i][4] = 0.0;
        
        lb_int_[i][0] = 0.0;
        lb_int_[i][1] = lb_seg_z_[i][0];
        lb_int_[i][2] = lb_seg_z_[i][1];
        lb_int_[i][3] = lb_seg_z_[i][2];
        lb_int_[i][4] = 0.0;
    }
}

/*!
 * Get the interaction energy with a feature at a given period.
 *
 * \param [in] dx Difference between x-coordinate and offset, when x has already been replaced in the box
 * \param [in] dz Z-distance to the feature's surface
 * \param [in] x_shift How much to shift imaginary feature (in period units) to obtain coordinates of one interested in
 * \param [in] m Expanded ensemble state of the system
 */
double rightTriangleXZ::featureInteraction_ (const double dx, const double dz, const int x_shift, const int m) {
    double ub = 0, lb = 0, U = 0.0;
    
    if (dx < ub_seg_x_[m][0]+x_shift) {
        ub = 0.0;
    } else if (dx < ub_seg_x_[m][1]+x_shift) {
        ub = ub_slope_[m][1]*(dx-ub_seg_x_[m][0]-x_shift) + ub_int_[m][1];
    } else if (dx < ub_seg_x_[m][2]+x_shift) {
        ub = ub_slope_[m][2]*(dx-ub_seg_x_[m][1]-x_shift) + ub_int_[m][2];
    } else if (dx < ub_seg_x_[m][3]+x_shift) {
        ub = ub_slope_[m][3]*(dx-ub_seg_x_[m][2]-x_shift) + ub_int_[m][3];
    } else {
        ub = 0.0;
    }
  
    if (dx < lb_seg_x_[m][0]+x_shift) {
        lb = 0.0;
    } else if (dx < lb_seg_x_[m][1]+x_shift) {
        lb = lb_slope_[m][1]*(dx-lb_seg_x_[m][0]-x_shift) + lb_int_[m][1];
    } else if (dx < lb_seg_x_[m][2]+x_shift) {
        lb = lb_slope_[m][2]*(dx-lb_seg_x_[m][1]-x_shift) + lb_int_[m][2];
    } else if (dx < lb_seg_x_[m][3]+x_shift) {
        lb = lb_slope_[m][3]*(dx-lb_seg_x_[m][2]-x_shift) + lb_int_[m][3];
    } else {
        lb = 0.0;
    }
                    
    if (dz >= ub) {
        U = 0.0;
    } else if (dz >= lb) {
        U = -eps_[m];
    } else {
        U = NUM_INFINITY;
    }
    
    return U;
}

/*!
 * Get the energy of a position assuming it is located within window 8.
 *
 * \param [in] atom Pointer to atom to examine
 * \param [in] box System box size. Will be checked that it is identical to value at class instantiation.
 */
double rightTriangleXZ::energy (const atom *a1, const std::vector < double > &box) {
    for (unsigned int i = 0; i < box_.size(); ++i) {
        if (box[i] != box_[i]) {
            throw customException ("System box size has changed from when rightTriangleXZ was instantiated");
        }
    }

    const int mv = a1->mState;
    std::vector < double > p = a1->pos;
    
    
    if (mv < 0 || mv > M_-1) {
        throw customException ("mState out of bounds for rightTriangleZ");
    }
    
    // First find nearest feature (the one right below)
    pbc (p, box);
    double dx = p[0] - xOffset_, dz = 0.0;
    int x_image = int(floor(dx/(width_+sep_)));
    double x_shift = x_image*(width_+sep_);
    
    if (top_) {
        dz = zbase_ - p[2];
    } else {
        dz = p[2] - zbase_;
    }
    
    if (dz > ub_check_[mv]) {
        return 0.0;
    }
    
    double U = featureInteraction_ (dx, dz, x_shift, mv);
    if (U == NUM_INFINITY) {
        return U;
    }
    
    // Must check all neighboring images, including images beyond edge of each box for periodicity effects
    int x_i = x_image+1;
    if (x_i > x_max_image_+1) {
        x_i = -1;
    }
    while (x_i != x_image) { // Stop once one complete cycle is finished
        x_shift = x_i*(width_+sep_);
        double dU = featureInteraction_ (dx, dz, x_shift, mv);
        if (dU == NUM_INFINITY) {
            return NUM_INFINITY;
        } else {
            U += dU;
        }
        x_i += 1;
        if (x_i > x_max_image_+1) {
            x_i = -1;
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
bool rightTriangleXZ::inside (const atom *a1, const std::vector < double > &box) {
    double U = NUM_INFINITY;
    try {
        U = energy (a1, box); // takes care of pbc internally
    } catch (customException &ce) {
        throw customException ("Unable to test if inside rightTriangleXZ : "+sstr(ce.what()));
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
void compositeBarrier::addRightTriangleXZ (const double width, const double theta, const double lamW, const double eps, const double sigma, const double sep, const double offset, const std::vector < double > &box, const double zbase, bool top, const int M) {
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
        sysBarriers_[sysBarriers_.size()-1] = new rightTriangleXZ (width, theta, lamW, eps, sigma, sep, offset, box, zbase, top, M);
    } catch (customException &ce) {
        throw customException ("Cannot add rightTriangleXZ to composite barrier: "+sstr(ce.what()));
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