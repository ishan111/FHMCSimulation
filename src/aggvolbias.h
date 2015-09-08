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
    aggVolBias (const int typeIndex, const int typeIndex2, const double pBias, const double rc1, const double rc2 const std::string tag);
    int make (simSystem &sys);

private:
    double getTempEnergy_ (const simSystem &sys, const std::vector < double > &box, const int chosenAtomType, const atom* chosenAtom);
	int typeIndex2_;
	double pBias_;
};

#endif
