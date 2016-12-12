#include "quaternion.h"

/*!
 * Set the quaternion's value.
 *
 * \param [in] q 4 component vector, [x, x, y, z]
 */
void quaternion::set (std::vector < double > &q) {
    if (q.size() != 4) {
        throw customException ("Quaternion must have 4 elements");
    }
    for (unsigned int i = 0; i < q.size(); ++i) {
        q_[i] = q;
    }
}
