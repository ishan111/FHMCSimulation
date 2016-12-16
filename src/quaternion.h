#ifndef QUATERNION_H_
#define QUATERNION_H_

#include <vector>
#include <cmath>
#include "global.h"
#include "utilities.h"

class quaternion {
public:
    quaternion () { q_.resize(4, 0.0); std::vector < double > dummy (3, 0); R_.resize(3, dummy); setRandomRot (); }
    ~quaternion () {};

    void set (const std::vector < double > &q); // Set the quaternion's values
    void setAxisAngle (const std::vector < double > &u, const double angle); // Set quaternion based on angle-axis representations
    void setRandomRot (); // Choose a quaternion corresponding to a random 3D rotation

    void conjugate (); // Conjugate the quaternion
    void normalize (); // Normalize the quaternion
    void inverse (); // Compute the inverse of the quaternion
    void translate (const std::vector < double > &t); // Translate the quaternion's "imaginary" part in 3D space

    double getNorm (); // Return the norm of the quaternion
    std::vector < double > get () { return q_; } // Return the quaternion's values
    std::vector < double > rotateVec (const std::vector < double > &vec); // Return the quaternion's values

    /*!
     * Add a quaternion to this one.
     *
     * \param [in] p Quaternion to add to this one.
     */
    quaternion operator+ (const quaternion &p) {
        quaternion ans;
        std::vector < double > sum (4, 0);
        for (unsigned int i = 0; i < 4; ++i) {
            sum[i] = this->q_[i] + p.q_[i];
        }
        ans.set(sum);
        return ans;
    }

    /*!
     * Subtract a quaternion from this one.
     *
     * \param [in] p Quaternion to subtract from this one.
     */
    quaternion operator- (const quaternion &p) {
        quaternion ans;
        std::vector < double > sum (4, 0);
        for (unsigned int i = 0; i < 4; ++i) {
            sum[i] = this->q_[i] - p.q_[i];
        }
        ans.set(sum);
        return ans;
    }

    /*!
     * Multiply this quaternion by another.  Does this in order of self*other.
     *
     * \param [in] other Quaternion to multiply this one by.
     *
     * \returns self*other
     */
    quaternion operator* (const quaternion &other) {
        quaternion ans;
        std::vector < double > prod (4, 0);
        prod[0] = this->q_[0]*other.q_[0] - this->q_[1]*other.q_[1] - this->q_[2]*other.q_[2] - this->q_[3]*other.q_[3];
        prod[1] = this->q_[0]*other.q_[1] + this->q_[1]*other.q_[0] + this->q_[2]*other.q_[3] - this->q_[3]*other.q_[2];
        prod[2] = this->q_[0]*other.q_[2] - this->q_[1]*other.q_[3] + this->q_[2]*other.q_[0] + this->q_[3]*other.q_[1];
        prod[3] = this->q_[0]*other.q_[3] + this->q_[1]*other.q_[2] - this->q_[2]*other.q_[1] + this->q_[3]*other.q_[0];
        ans.set(prod);
        return ans;
    }

    /*!
     * Check if two quaternions are equal to each other to within a tolerance of 1.0e-12.
     *
     * \param [in] other Other quaternion to compare with.
     */
    bool operator== (const quaternion &other) {
        for (unsigned int i = 0; i < 4; ++i) {
            if (fabs(this->q_[i] - other.q_[i]) > 1.0e-12) {
                return false;
            }
        }
        return true;
    }

    /*!
     * Check if two quaternions are not equal to each other to within a tolerance of 1.0e-12.
     *
     * \param [in] other Other quaternion to compare with.
     */
    bool operator!= (const quaternion &other) {
        return !(*this == other);
    }

private:
    std::vector < double >  q_; //!< Quanternion vector
    std::vector < std::vector < double > > R_; //!< Rotation matrix corresponding to q_
    void assignq_ (const std::vector < double > &qv); // Assign components of quaternion
    void assignr_ (const std::vector < double > &qv); // Assign rotation matrix of quaternion
};

#endif
