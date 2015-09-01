#ifndef MOVES_H_
#define MOVES_H_

#define MOVE_SUCCESS 1 //!< Returned by mcMove if move was successful
#define MOVE_FAILURE 0 //!< Returned by mcMove if move was not successful

#include <vector>
#include <string>
#include "system.h"
#include "utilities.h"
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

/*!
 * Class that tracks and decides which moves whould be made.  However, it does NOT store the moves themselves so they should be fixed in memory elsewhere.
 */
class moves {
public:
    moves (const int M = 1);
    ~moves ();
    
    void makeMove (simSystem &sys);	
    void addMove (mcMove *newMove, const double probability);
    std::vector < std::vector < double > > reportMoveStatistics ();
    std::vector < double > reportProbabilities () { return normProbabilities_; } //!< Echo the normalized probabilities of each move in the object
	const std::vector < mcMove* > includedMoves () { return moves_; } //!< Returns a vector of pointers to move objects currently being used
	
private:
    std::vector < double > normProbabilities_; //!< Sum of un-normalized probability of each move included
    std::vector < double > rawProbabilities_; //!< Un-normalized probabilty of each move
    std::vector < std::vector < double > > succeeded_; //!< Number of times each move was successful
    std::vector < std::vector < double > > attempted_; //!< Number of times each move was attempted
    std::vector < mcMove* > moves_; //!< Vector of pointers to all moves used
    int M_; //!< Number of stages for insert/delete moves
};

#endif
