#include "atom.h"
#include "utilities.h"


/*!
 * Instantiate an atom with centers oriented relative to the center.
 *
 * \param [in] ncenters Number of centers
 * /param [in] relOr Vectors pointing to rigid centers from atom's pos
 */
atom::atom (int ncenters, std::vector < std::vector < double > > relOr) {
	try {
		std::vector < double > dummy (3, 0);
		vecToCenters.resize(ncenters, dummy);
	} catch (std::bad_alloc &ba) {
		throw customException ("Out of memory for atom centers");
	}

	if (relOr.size() != ncenters) {
		throw customException ("Must specify exactly one 3D orientation for each center");
	} else {
		for (unsigned int i = 0; i < ncenters; ++i) {
			if (relOr[i].size() != 3) {
				throw customException ("Error - 3D atom center orientations must be in 3D");
			}
			for (unsigned int j = 0; j < 3; ++j) {
				vecToCenters[i][j] = relOr[i][j];
			}
		}
	}
}

/*!
 * Rotate centers with Euler angles.
 *
 * \param [in] alpha Radians to rotate centers by around x-axis
 * \param [in] beta Radians to rotate centers by around y-axis
 * \param [in] gamma Radians to rotate centers by around z-axis
 */
/*void atom::rotateCenters (double alpha, double beta, double gamma) {
	// assert valid ranges - https://en.wikipedia.org/wiki/Euler_angles
	if (alpha < -PI or alpha >= PI) {
		throw customException ("Invalid range for alpha");
	}
	if (gamma < -PI or gamma >= PI) {
		throw customException ("Invalid range for gamma");
	}
	if (beta < 0 or beta > PI) {
		throw customException ("Invalid range for beta");
	}

	std::vector < std::vector < double > > R = rotationMatrix(alpha, beta, gamma);

	for (unsigned int i = 0; i < vecToCenters.size(); ++i) {
		const double mag = sqrt(vecToCenters[i][0]*vecToCenters[i][0] + vecToCenters[i][1]*vecToCenters[i][1] + vecToCenters[i][2]*vecToCenters[i][2]);
		std::vector < double > tmpVec (3, 0), ans (3, 0);

		// scale to unit vector
		for (unsigned int j = 0; j < 3; ++j) {
			tmpVec[j] = vecToCenters[i][j]/mag;
		}

		// rotate
		for (unsigned int j = 0; j < 3; ++j) {
			for (unsigned int k = 0; k < 3; ++k) {
				ans[j] += R[j][k]*tmpVec[k];
			}
		}

		// scale back
		for (unsigned int j = 0; j < 3; ++j) {
			vec_to_centers[i][j] = mag*ans[j];
		}
	}
}*/
