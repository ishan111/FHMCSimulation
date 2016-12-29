#ifndef ATOM_H_
#define ATOM_H_

#include <vector>
#include "global.h"
#include "quaternion.h"

/*!
 * Atom class with rigid internal degrees of freedom besides its center of mass.
 */

class atom {
public:
	atom () { pos.resize(3, 0); mState = 0; } // Assumes no centers
	atom (unsigned int ncenters, std::vector < std::vector < double > > rel_or); // Instantiate with rigid centers oriented according to these vectors
	~atom () {};

	void rotateCenters (const quaternion &q); // Rotate the atom's center(s) using quaternions

	int mState; //!< State of fraction insertion of the atom in the expanded ensemble, 0 = fully inserted
	std::vector < double > pos; //!< 3D position
	std::vector < std::vector < double > > vecToCenters; //!< Vectors pointing from pos to rigid centers
};

#endif
