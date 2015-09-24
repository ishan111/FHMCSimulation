#include "barrier"

/*!
 * Instantiate a hard wall with boundaries in the +/- z direction.
 *
 * \param [in] lb z-Position of the lower wall
 * \param [in] ub z-Position of the upper wall
 */
hardWallZ::hardWallZ (const double lb, const double ub) {
    if (lb >= ub) {
        throw customException ("hardWallZ must have lower bound < upper bound");
    }
    
    ub_ = ub;
    lb_ = lb;
}

/*!
 * Return whether or not a point falls between the walls (subject to a hard-sphere exclusion radius).
 *
 * \param [in] point (x, y, z) Point to test - this does NOT need to be in the simulation box a priori
 * \param [in] params Vector containing one entry: the hard-sphere DIAMETER of the particle (how close to allow the sphere to approach the wall's surface)
 * \param [in] box Simulation box
 */
bool hardWallZ::inside (const std::vector < double > &point, const std::vector < double > params, const std::vector < double > &box) {
    std::vector < double > p = point;
    pbc (p, box);
    if (p.z => ub_ - params[0]/2.0 || pz <= lb_ + params[0]/2.0) {
        return false;
    } else {
        return true;
    }
}