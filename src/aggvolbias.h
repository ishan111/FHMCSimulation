#ifndef AGGVOLBIAS_H_
#define AGGVOLBIAS_H_

#include <cmath>
#include "system.h"
#include "global.h"
#include "moves.h"
#include "atom.h"
#include "utilities.h"

class aggVolBias : public mcMove {
public:
    aggVolBias () { changeN_ = false; pBias_ = 0; }
    aggVolBias (const int typeIndex, const int typeIndex2, const double pBias, const double rc1, const double rc2, const std::string tag);
    int make (simSystem &sys);

private:
    double getTempEnergy_ (simSystem &sys, const std::vector < double > &box, const double V, const int chosenAtomType, atom* chosenAtom);
	int typeIndex2_; //!< Atom type for pk K
	double pBias_; //!< Probability to choose a particle "in" J over others to move
	double rc1_; //!< Bonding radius around pk J
	double rc2_; //!< Bonding radius around pk K
};

#endif
