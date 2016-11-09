#include "atom.h"
#include "utilities.h"


/*!
 * Instantiate an atom with centers oriented relative to the center.
 *
 * \param [in] ncenters Number of centers
 * /param [in] rel_or Vectors pointing to rigid centers from atom's pos
 */
atom::atom (int ncenters, std::vector < std::vector < double > > rel_or) {
	try {
		std::vector < double > dummy (3, 0);
		vec_to_centers.resize(ncenters, dummy);
	} catch (std::bad_alloc &ba) {
		throw customException ("Out of memory for atom centers");
	}

	if (rel_or.size() != ncenters) {
		throw customException ("Must specify exactly one 3D orientation for each center");
	} else {
		for (unsigned int i = 0; i < ncenters; ++i) {
			if (rel_or[i].size() != 3) {
				throw customException ("Error - 3D atom center orientations must be in 3D");
			}
			for (unsigned int j = 0; j < 3; ++j) {
				vec_to_centers[i][j] = rel_or[i][j];
			}
		}
	}
}

/*!
 * Rotate centers
 *
 * \param [in] alpha Radians to rotate centers by around x-axis
 * \param [in] beta Radians to rotate centers by around y-axis
 * \param [in] gamma Radians to rotate centers by around z-axis
 */
void atom::rotate_centers (double alpha, double beta, double gamma) {
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

	for (unsigned int i = 0; i < vec_to_centers.size(); ++i) {
		const double mag = sqrt(vec_to_centers[i][0]*vec_to_centers[i][0] + vec_to_centers[i][1]*vec_to_centers[i][1] + vec_to_centers[i][2]*vec_to_centers[i][2]);
		std::vector < double > tmpVec (3, 0), ans (3, 0);

		// scale to unit vector
		for (unsigned int j = 0; j < 3; ++j) {
			tmpVec[j] = vec_to_centers[i][j]/mag;
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
}
