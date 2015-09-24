#ifndef BARRIER_H_
#define BARRIER_H_

#include <vector>
#include "global.h"
#include "utilities.h"

/*!
 * Virtual base class for barriers.
 */
class barrier {
public:
	barrier () {};
	virtual ~barrier () = 0;
    bool inside (const std::vector < double > &point, const std::vector < double > params, const std::vector < double > &box) = 0;
};

/*!
 * Parallel hard walls in the z-direction.
 */
class hardWallZ : public barrier {
public:
    ~hardWallZ () {};
    hardWallZ (const double lb, const double ub);
    bool inside (const std::vector < double > &point, const std::vector < double > params, const std::vector < double > &box);
    
private:
    double lb_; //!< Lower bound for wall in the z-direction
    double ub_; //!< Uppber bound for all in the z-direction
};
#endif
