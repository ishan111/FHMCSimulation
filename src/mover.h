#ifndef MOVER_H_
#define MOVER_H_

#include <memory>
#include <vector>
#include <string>
#include <fstream>
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
    std::vector < std::vector < double > > reportMoveStatistics ();
    std::vector < double > reportProbabilities () { return normProbabilities_; } //!< Echo the normalized probabilities of each move in the object
	void print (const std::string filename);
	void addInsert (const int index, const double prob);
	void addDelete (const int index, const double prob);
	void addSwap(const int index1, const int index2, const double prob);
	void addTranslate (const int index, const double prob, const double maxD, const std::vector < double > &box);
    int getM () { return M_; }
    void setM (const int M);

    //void addMove (mcMove *newMove, const double probability); //DEPR
    //const std::vector < mcMove* > includedMoves () { return moves_; } //!< Returns a vector of pointers to move objects currently being used //DEPR
    
private:
    std::vector < double > normProbabilities_; //!< Sum of un-normalized probability of each move included
    std::vector < double > rawProbabilities_; //!< Un-normalized probabilty of each move
    std::vector < std::vector < double > > succeeded_; //!< Number of times each move was successful
    std::vector < std::vector < double > > attempted_; //!< Number of times each move was attempted
    int M_; //!< Number of stages for insert/delete moves

    /*std::vector < mcMove* > moves_; //!< Vector of pointers to all moves used //DEPR
	std::vector < insertParticle > insert_; //DEPR
	std::vector < deleteParticle > delete_; //DEPR
	std::vector < translateParticle > translate_; //DEPR
	std::vector < swapParticles > swap_; //DEPR*/

    std::vector < std::shared_ptr < mcMove > > ownedMoves_;
    void addOn_ (bool changeN, const double probability);
};

#endif
