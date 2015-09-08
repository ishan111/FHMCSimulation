#ifndef AGGVOLBIAS_H_
#define AGGVOLBIAS_H_

#include <cmath>
#include "system.h"
#include "global.h"
#include "moves.h"
#include "atom.h"
#include "utilities.h"

class aggVolBias3 : public mcMove {
public:
    	aggVolBias3 () { changeN_ = false; pBias_ = 0; }
    	aggVolBias3 (const int typeIndex, const int typeIndex2, const double pBias, const std::vector < double > rc1, const std::vector < double > rc2, const std::string tag);
    	int make (simSystem &sys);

private:
    	double getTempEnergy_ (simSystem &sys, const std::vector < double > &box, const double V, const int chosenAtomType, atom* chosenAtom);
	int typeIndex2_; //!< Atom type for pk K
	double pBias_; //!< Probability to choose a particle "in" J over others to move
	double rc1min_; //!< Minimum bonding radius around pk J
	double rc2min_; //!< Minimum bonding radius around pk K
	double rc1max_; //!< Maximum bonding radius around pk J
	double rc2max_; //!< Maximum bonding radius around pk K
};

class aggVolBiasInsert : public mcMove {
public:
	aggVolBiasInsert () { changeN_ = true; pBias_ = 0; }
	aggVolBiasInsert (const int typeIndex, const double pBias, const std::vector < double > rc, const std::string tag);
	int make (simSystem &sys);

private:
	double pBias_;
	double rcmin_; //!< Minimum bonding radius around pk J
        double rcmax_; //!< Maximum bonding radius around pk J
};

class aggVolBiasDelete : public mcMove {
public:
        aggVolBiasDelete () { changeN_ = true; pBias_ = 0; }
        aggVolBiasDelete (const int typeIndex, const double pBias, const std::vector < double > rc, const std::string tag);
        int make (simSystem &sys);

private:
        double pBias_;
        double rcmin_; //!< Minimum bonding radius around pk J
        double rcmax_; //!< Maximum bonding radius around pk J
};

#endif
