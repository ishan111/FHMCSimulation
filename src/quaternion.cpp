#include "quaternion.h"

/*!
 * Set the quaternion's value.
 *
 * \param [in] qv 4 component vector, [x, x, y, z]
 */
void quaternion::set (const std::vector < double > &qv) {
    assignq_ (qv);
}

/*!
 * Assign a quaternion's rotation matrix based on (unnormalized) q.
 *
 * \param [in] qv 4 component vector, [x, x, y, z]
 */
void quaternion::assignr_ (const std::vector < double > &qv) {
    quaternion qnorm;
    qnorm.assignq_(qv);
    qnorm.normalize();
    std::vector < double > qn = qnorm.get();

    R_[0][0] = 1.0 - 2.0*qn[2]*qn[2] - 2.0*qn[3]*qn[3];
    R_[0][1] = 2.0*(qn[1]*qn[2] - qn[3]*qn[0]);
    R_[0][2] = 2.0*(qn[1]*qn[3] + qn[2]*qn[0]);

    R_[1][0] = 2.0*(qn[1]*qn[2] + qn[3]*qn[0]);
    R_[1][1] = 1.0 - 2.0*qn[1]*qn[1] - 2.0*qn[3]*qn[3];
    R_[1][2] = 2.0*(qn[2]*qn[3] - qn[1]*qn[0]);

    R_[2][0] = 2.0*(qn[1]*qn[3] - qn[2]*qn[0]);
    R_[2][1] = 2.0*(qn[2]*qn[3] + qn[1]*qn[0]);
    R_[2][2] = 1.0 - 2.0*qn[1]*qn[1] - 2.0*qn[2]*qn[2];
}

/*!
 * Assign a quaternion's components.
 *
 * \param [in] qv 4 component vector, [x, x, y, z]
 */
void quaternion::assignq_ (const std::vector < double > &qv) {
    if (qv.size() != 4) {
        throw customException ("Quaternion must have 4 elements");
    }
    for (unsigned int i = 0; i < qv.size(); ++i) {
        q_[i] = qv[i];
    }
}

/*!
 * Transform quaternion into its conjugate.
 */
void quaternion::conjugate () {
    for (unsigned int i = 1; i < 4; ++i) {
        q_[i] = -q_[i];
    }
}

/*!
 * Normalize a quaternion.
 *
 * \param [in] t Translation vector.
 */
void quaternion::normalize () {
    const double norm = getNorm();
    for (unsigned int i = 0; i < 4; ++i) {
        q_[i] /= norm;
    }
}

/*!
 * Transform a quanternion into its inverse.
 */
void quaternion::inverse () {
    conjugate();
    double n2 = 0.0;
    for (unsigned int i = 0; i < 4; ++i) {
        n2 += q_[i]*q_[i];
    }
    for (unsigned int i = 0; i < 4; ++i) {
        q_[i] /= n2;
    }
}

/*!
 * Translate a quaternion in 3D space.
 *
 * \param [in] t Translation vector.
 */
void quaternion::translate (const std::vector < double > &t) {
    if (t.size() != 3) {
        throw customException ("Translation vector for Quaternion must have 3 elements");
    }
    for (unsigned int i = 1; i < 4; ++i) {
        q_[i] += t[i-1];
    }
}

/*!
 * Return the norm of the quaternion.
 */
double quaternion::getNorm () {
    double n2 = 0.0;
    for (unsigned int i = 0; i < 4; ++i) {
        n2 += q_[i]*q_[i];
    }
    return std::sqrt(n2);
}

/*!
 * Use this quaternion to operate on a vector, that causes its rotation.
 *
 * \param [in] vec Vector to rotate.
 */
std::vector < double > quaternion::rotateVec (const std::vector < double > &vec) {
    std::vector < double > ans (3, 0);
    assignr_ (q_);

    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            ans[i] += R_[i][j]*vec[j];
        }
    }

    return ans;
}

/*!
 * Pick a quaternion that corresponds to a random 3D rotation.
 */
 void quaternion::setRandomRot () {
     // spherical coordinates to get random unit vector
     const double theta = std::acos(2.0*rng(&RNG_SEED) - 1.0);
     const double phi = 2.0*PI*rng(&RNG_SEED), sin_phi = std::sin(phi);

     // randomly pick angle
     const double angle = 2.0*PI*rng(&RNG_SEED), s = std::sin(angle/2.0);
     q_[0] = std::cos(angle/2.0);
     q_[1] = s*1.0*sin_phi*std::cos(theta);
     q_[2] = s*1.0*sin_phi*std::sin(theta);
     q_[3] = s*1.0*std::cos(phi);
 }

 /*!
  * Set the quaternion based on angle, axis desired.  Does not need to be normalized.
  *
  * \param [in] angle Angle of rotation desired relative to axis (in radians)
  * \param [in] aixs 3D axis to rotate about (in right-handed coordinates)
  */
void quaternion::setAxisAngle (const std::vector < double > &u, const double angle) {
    if (u.size() != 3) {
        throw customException ("Axis for Quaternion must have 3 elements");
    }
    q_[0] = std::cos(angle/2.0);
    const double s = std::sin(angle/2.0);
    q_[1] = s*u[0];
    q_[2] = s*u[1];
    q_[3] = s*u[2];
}
