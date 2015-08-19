#include "moves.h"

mcMove::~mcMove () {
    ;
}

moves::moves () {
	;
}

moves::~moves () {
	;
}

/*!
 * Add a new move to the object.
 *
 * \param [in] newMove Pointer to a newly instantiated move.  This is stored as a pointer, so the move cannot be moved in memory later.
 * \param [in] probability Unnormalized probability of making this move.
 */
void moves::addMove (mcMove *newMove, const double probability) {
	// add new move to the class
	moves_.push_back(newMove);
	rawProbabilities_.push_back(probability);
	normProbabilities_.resize(rawProbabilities_.size());
	succeeded_.resize(rawProbabilities_.size(), 0.0);
    attempted_.resize(rawProbabilities_.size(), 0.0);
    
	// update move probabilities
	double sum = 0.0;
	for (unsigned int i = 0; i < rawProbabilities_.size(); ++i) {
		sum += rawProbabilities_[i];
	}
	
	normProbabilities_[0] = rawProbabilities_[0]/sum;
	for (unsigned int i = 1; i < rawProbabilities_.size(); ++i) {
		normProbabilities_[i] = rawProbabilities_[i]/sum + normProbabilities_[i-1];
	}
	
	// for exactness, specify the upper bound
	normProbabilities_[normProbabilities_.size()-1] = 1.0;
}

/*!
 * Choose a move to make. If in an expanded ensemble, will restrict moves which change the number of particles to the atom type 
 * that is currently on partially in the system.
 *  
 * \param [in] sys simSystem object to make a move in.
 */
void moves::makeMove (simSystem &sys) {
	int moveChosen = -1, succ = 0;
	bool done = false;
	while (!done) {
		const double ran = rng (&RNG_SEED);
		for (unsigned int i = 0; i < normProbabilities_.size(); ++i) {
			if (ran < normProbabilities_[i]) {
				if (sys.getTotalM() > 1) {
					// expanded ensemble has to check the moves because have to only work on the partially inserted atom 
					if ((moves_[i]->changeN() == true) && (moves_[i]->whatType() != sys.getFractionalAtomType()) && (sys.getCurrentM() > 0)) {
						// reject this choice because we must only insert/delete the type that is already partially inserted IFF we are *already* in a partially inserted state
						// choose a new move
						
						done = false;
						break;
					} else {
						try {
							succ = moves_[i]->make(sys);
						} catch (customException &ce) {
							std::string a = "Failed to make a move properly";
							std::string b = ce.what();
							throw customException(a+b);
						}
						done = true;
			    		moveChosen = i;
			    		break;
					}
				} else {
					// without expanded ensemble, inserts/deletes can proceed unchecked
					try {
						succ = moves_[i]->make(sys);
					} catch (customException &ce) {
						std::string a = "Failed to make a move properly";
						std::string b = ce.what();
						throw customException(a+b);
					}
					done = true;
		    			moveChosen = i;
		    			break;
				}
			}
		}
	}
	
	if (moveChosen < 0) {
		throw customException("Failed to choose a move properly");
	}
	
    attempted_[moveChosen] += 1.0;
    succeeded_[moveChosen] += succ;
}

/*!
 * Report the statistics on the success/failure of each move made so far.
 *
 * \return ans Number of Success / Total Attempts for each move
 */
std::vector < double > moves::reportMoveStatistics () {
    std::vector < double > ans = succeeded_;
    if (attempted_.begin() == attempted_.end()) {
        throw customException ("No moves added to system");
    }
    for (unsigned int i = 0; i < attempted_.size(); ++i) {
        ans[i] /= attempted_[i];
    }
    return ans;
}
