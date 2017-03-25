#include "atom.h"
#include "utilities.h"

/*!
 * Rotate all atom centers using a quaternion.
 *
 * \param [in] q Quaternion to rotate the centers with
 */
void atom::rotateCenters (quaternion &q) {
	for (std::vector < std::vector < double > >::iterator it = vecToCenters.begin(); it != vecToCenters.end(); ++it) {
		*it = q.rotateVec(*it);
	}
}

/*!
 * Instantiate an atom with centers oriented relative to the center.
 *
 * \param [in] nCenters Number of centers
 * /param [in] relOr Vectors pointing to rigid centers from atom's pos
 */
atom::atom (unsigned int nCenters, std::vector < std::vector < double > > relOr) {
	pos.resize(3, 0);
	mState = 0;

	try {
		std::vector < double > dummy (3, 0);
		vecToCenters.resize(nCenters, dummy);
	} catch (std::bad_alloc &ba) {
		throw customException ("Out of memory for atom centers");
	}

	if (relOr.size() != nCenters) {
		throw customException ("Must specify exactly one 3D orientation for each center");
	} else {
		for (unsigned int i = 0; i < nCenters; ++i) {
			if (relOr[i].size() != 3) {
				throw customException ("Error - 3D atom center orientations must be in 3D");
			}
			for (unsigned int j = 0; j < 3; ++j) {
				vecToCenters[i][j] = relOr[i][j];
			}
		}
	}
}
