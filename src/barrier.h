#ifndef BARRIER_H_
#define BARRIER_H_

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
