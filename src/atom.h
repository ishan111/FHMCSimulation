#ifndef ATOM_H_
#define ATOM_H_

#include <vector>
#include "global.h"

/*!
 * Atom class with rigid internal degrees of freedom besides its center of mass.
 */

class atom {
public:
	atom () { pos.resize(3, 0); mState = 0; } // assumes no centers
	atom (int ncenters, std::vector < std::vector < double > > rel_or); // Instantiate with rigid centers oriented according to these vectors
	~atom () {};

	int mState; //!< State of fraction insertion of the atom in the expanded ensemble, 0 = fully inserted
	std::vector < double > pos; //!< 3D position
	std::vector < std::vector < double > > centers; //!< Rigid centers

	void rotateCenters (double alpha, double beta, double gamma); //!< Rotate centers
};

#endif
