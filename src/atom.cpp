#include "atom.h"
#include "utilities.h"


/*!
 * Instantiate an atom with centers oriented relative to the center.
 *
 * \param [in] ncenters Number of centers
 * /param [in] relOr Vectors pointing to rigid centers from atom's pos
 */
atom::atom (unsigned int ncenters, std::vector < std::vector < double > > relOr) {
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
