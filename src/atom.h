#ifndef ATOM_H_
#define ATOM_H_

#include <vector>

/*!
 * Atom class. 
 */
class atom {
public:
	atom () { pos.resize(3, 0); }
	~atom () {};
	
	std::vector < double > pos; //!< 3D position
};

#endif
