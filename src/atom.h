#ifndef ATOM_H_
#define ATOM_H_

#include <vector>

/*!
 * Atom class. 
 */
class atom {
public:
	atom () { pos.resize(3, 0); mState = 0; }
	~atom () {};
	
	int mState; //!< State of fraction insertion of the atom in the expanded ensemble, 0 = fully inserted
	std::vector < double > pos; //!< 3D position
};

#endif
