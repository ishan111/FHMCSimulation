#ifndef BARRIER_H_
#define BARRIER_H_

#include <cmath>
#include <vector>
#include "global.h"
#include "utilities.h"
#include "potentials.h"
#include "atom.h"

/*!
 * Virtual base class for barriers.  It is intended that there should be a separate barrier class defined for each species type.
 */
class barrier {
public:
    virtual ~barrier () {;}
    virtual bool inside (const atom *a1, const std::vector < double > &box) = 0;
    virtual double energy (const atom *a1, const std::vector < double > &box) = 0;
    
protected:
    int M_;
};

/*!
 * Parallel hard walls in the z-direction.
 */
class hardWallZ : public barrier {
public:
    ~hardWallZ () {};
    hardWallZ (const double lb, const double ub, const double sigma, const int M = 1);
    bool inside (const atom *a1, const std::vector < double > &box);
    double energy (const atom *a1, const std::vector < double > &box);
    
private:
    double lb_; //!< Lower bound for wall in the z-direction
    double ub_; //!< Uppber bound for all in the z-direction
    double sigma_; //!< Hard-sphere diameter the species this wall interacts with can approach within
};

/*!
 * Parallel square-well walls in the z-direction.
 */
class squareWellWallZ : public barrier {
public:
    ~squareWellWallZ () {};
    squareWellWallZ (const double lb, const double ub, const double sigma, const double range, const double eps,  const int M = 1);
    bool inside (const atom *a1, const std::vector < double > &box);
    double energy (const atom *a1, const std::vector < double > &box);
    
private:
    double lb_; //!< Lower bound for wall in the z-direction
    double ub_; //!< Uppber bound for all in the z-direction
    double range_; //!< Distance normal to the wall's surface where there is an interaction
    double eps_; //!< Magnitude of the interaction
    double sigma_; //!< Hard-sphere diameter the species this wall interacts with can approach within
};

class rightTriangleXZ : public barrier {
public:
    ~rightTriangleXZ () {};
    rightTriangleXZ (const double width, const double theta, const double lamW, const double eps, const double sigma, const double sep, const double offset, const std::vector < double > &box, const double zbase, bool top = false, const int M = 1);
    bool inside (const atom *a1, const std::vector < double > &box);
    double energy (const atom *a1, const std::vector < double > &box);
    
private:
    double featureInteraction_ (const double dx, const double dz, const int x_shift, const int m);
    
    bool top_; //!< If true, feature is on the "top", else is on the bottom (default)
    int x_max_image_; //!< Max image index of the feature along the x-direction
    double width_; //!< Width of triangle's feature
    double theta_; //!< Elevation angle of the feature
    double sep_; //!< Distance between features
    double xOffset_; //!< Offset from x = 0 position of the first feature
    double zbase_; //!< Z-coordinate of XY plane that defines the base of the feature
    double lamW_; //!< Range of interaction with wall in terms of sigma
    double a_, b_, c_;

    std::vector < double > ub_check_; //!< Max range of interaction as a function of M
    std::vector < double > eps_; //!< Attraction strength to feature as a function of M
    std::vector < double > sigma_; //!< Hard sphere diameter of interaction with the feature as a function of M
    std::vector < double > box_; //!< Box size the feature was initialized with
    std::vector < std::vector < double > > ub_seg_x_, lb_seg_x_, ub_seg_z_, lb_seg_z_, ub_slope_, ub_int_, lb_slope_, lb_int_;
};

/*!
 * Class which tracks all barriers (superimposed) which interact with a given species.
 */
class compositeBarrier {
public:
    compositeBarrier () {};
    ~compositeBarrier ();
    
    void addHardWallZ (const double lb, const double ub, const double sigma, const int M = 1);
    void addSquareWellWallZ (const double lb, const double ub, const double sigma, const double range, const double eps, const int M = 1);
    void addRightTriangleXZ (const double width, const double theta, const double lamW, const double eps, const double sigma, const double sep, const double offset, const std::vector < double > &box, const double zbase, bool top = false, const int M = 1);
    
    bool inside (const atom *a1, const std::vector < double > &box);
    double energy (const atom *a1, const std::vector < double > &box);
    
private:
    std::vector < barrier* > sysBarriers_;
};

#endif
