#ifndef MOVER_H_
#define MOVER_H_

#include <vector>
#include <string>
#include "system.h"
#include "insert.h"
#include "delete.h"
#include "swap.h"
#include "translate.h"
#include "moves.h"
#include "utilities.h"
#include "global.h"

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
	void print (const std::string filename);

	void addInsert (const int index, const double prob);
	void addDelete (const int index, const double prob);
	void addSwap(const int index1, const int index2, const double prob);
	void addTranslate (const int index, const double prob);

private:
    std::vector < double > normProbabilities_; //!< Sum of un-normalized probability of each move included
    std::vector < double > rawProbabilities_; //!< Un-normalized probabilty of each move
    std::vector < std::vector < double > > succeeded_; //!< Number of times each move was successful
    std::vector < std::vector < double > > attempted_; //!< Number of times each move was attempted
    std::vector < mcMove* > moves_; //!< Vector of pointers to all moves used
    int M_; //!< Number of stages for insert/delete moves

	std::vector < insertParticle > insert_;
	std::vector < deleteParticle > delete_;
	std::vector < translateParticle > translate_;
	std::vector < swapParticles > swap_;

};

#endif
