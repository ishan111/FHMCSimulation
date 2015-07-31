#ifndef MOVES_H_
#define MOVES_H_

#define MOVE_SUCCESS 1 //!< Returned by mcMove if move was successful
#define MOVE_FAILURE 0 //!< Returned by mcMove if move was not successful

#include "system.h"
#include <vector>
#include <string>
#include <boost/lexical_cast.hpp>

/*!
 * Virtual base class for all Monte Carlo moves.
 */
class mcMove {
public:
	mcMove () {};
	mcMove (const int typeIndex, const std::string tag) { typeIndex_ = typeIndex; name_ = tag + boost::lexical_cast<std::string>(typeIndex); }  
    virtual ~mcMove () = 0;
    virtual int make (simSystem &sys) = 0; //!< Make a MC move, return MOVE_SUCCESS or MOVE_FAILURE
	const int whatType () { return typeIndex_; }
	const std::string myName () { return name_; }	//!< Return the name of this move
	
protected:
	int typeIndex_;   //!< Species index this move will operate on
	std::string name_;	//!< Move name
};

/*!
 * Class that tracks and decides which moves whould be made.
 */
class moves {
public:
    moves ();
    ~moves ();
    
    void makeMove (simSystem &sys);	
    void addMove (mcMove *newMove, const double probability);
    std::vector < double > reportMoveStatistics ();
    std::vector < double > reportProbabilities () { return normProbabilities_; } //!< Echo the normalized probabilities of each move in the object
	const std::vector < mcMove* > includedMoves () { return moves_; }
	
private:
    std::vector < double > normProbabilities_, rawProbabilities_, succeeded_, attempted_;
    std::vector < mcMove* > moves_;
};

#endif