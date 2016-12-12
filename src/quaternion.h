#ifndef QUATERNION_H_
#define QUATERNION_H_

#include <vector>
#include "global.h"

class quaternion {
public:
    quaternion () { setRandomRot(); q_.resize(4, 0.0); }
    ~quaternion () {};

    void set (std::vector < double > &q); // Set the quaternion's values
    void setRandomRot (); // Choose a quaternion corresponding to a random 3D rotation

    quaternion translate (); // Translate the quaternion's "imaginary" part in 3D space
    quaternion conjugate (); // Conjugate the quaternion
    quaternion normalize (); // Normalize the quaternion
    quaternion inverse (); // Compute the inverse of the quaternion

    quaternion operator+ (const quaternion &p); // Add two quaternions
    quaternion operator- (const quaternion &p); // Subtract two quaternions
    quaternion operator* (const quaternion &p); // Multiply two quaternions

private:
    std::vector < double >  q_; //!< Quanternion vector
};

#endif
