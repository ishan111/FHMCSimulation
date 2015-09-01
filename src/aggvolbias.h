#ifndef AGGVOLBIAS_H_
#define AGGVOLBIAS_H_

#include <cmath>
#include "system.h"
#include "global.h"
#include "moves.h"
#include "atom.h"

class aggVolBias : public mcMove {
public:
        aggVolBias () { changeN_ = false; pBias_ = 0; }
        aggVolBias (const int typeIndex, const int typeIndex2, const double pBias; const std::string tag);
        int make (simSystem &sys);

private:
	int typeIndex2_;
	double pBias_;
};

#endif
