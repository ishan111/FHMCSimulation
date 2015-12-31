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

class rightTriangleZ : public barrier {
public:
    ~rightTriangleZ () {};
    rightTriangleZ (const double width, const double theta, const double lamW, const double eps, const double sigma, const double sep, const double offset, const std::vector < double > &box, bool top = false, const int M = 1);
    bool inside (const atom *a1, const std::vector < double > &box);
    double energy (const atom *a1, const std::vector < double > &box);
    
private:
    double getRange_ (double x);
    double getU0_ (const double x, const double z);
    double getU1_ (const double x, const double z);
    double getU2_ (const double x, const double z);
    double getU3_ (const double x, const double z);
    double getU4_ (const double x, const double z);
    double getU5_ (const double x, const double z);
    double getU6_ (const double x, const double z);
    double getU7_ (const double x, const double z);
    
    double lam_; //!< Width of triangle's feature
    double theta_; //!< Elevation angle of the feature
    double lamWall_; //!< Attractive range beyond the hard sphere in contact with the feature
    double eps_; //!< Attraction strength to feature
    double sigma_; //!< Hard sphere diameter of interaction with the feature
    double lamP_; //!< Distance between features
    double xOffset_; //!< Offset from x = 0 position of the first feature
    bool top_; //!< If true, feature is on the "top", else is on the bottom (default)
    double cosTheta_, sinTheta_, a_, b_, c_, m_, n_, p_, q_, u_, v_, j_, k_, s_;
    std::vector < double > lbounds_, ubounds_, box_;
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
    
    bool inside (const atom *a1, const std::vector < double > &box);
    double energy (const atom *a1, const std::vector < double > &box);
    
private:
    std::vector < barrier* > sysBarriers_;
};

#endif
