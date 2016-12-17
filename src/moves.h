#ifndef MOVES_H_
#define MOVES_H_

#define MOVE_SUCCESS 1 //!< Returned by mcMove if move was successful
#define MOVE_FAILURE 0 //!< Returned by mcMove if move was not successful

#include <vector>
#include <string>
#include "system.h"
#include "global.h"

/*!
 * Virtual base class for all Monte Carlo moves.
 */
class mcMove {
public:
	mcMove () {};
	mcMove (const int typeIndex, const std::string tag) { typeIndex_ = typeIndex; name_ = tag+sstr(typeIndex); }
   	virtual ~mcMove () = 0;
    virtual int make (simSystem &sys) = 0; //!< Make a MC move, return MOVE_SUCCESS or MOVE_FAILURE
	const bool changeN () { return changeN_; }	//!< Returns whether or not the move has the ability to change the net number of particles in the system
	const int whatType () { return typeIndex_; } //!< Returns the index referring to the atom type this move operates on
	const std::string myName () { return name_; }	//!< Return the name of this move

protected:
	bool changeN_; 	//!< Does this move have the capacity to create a net change in the total number of particles?
	int typeIndex_;   //!< Species index this move will operate on
	std::string name_;	//!< Move name
};

#endif
