#ifndef SWAP_H_
#define SWAP_H_

#include <cmath>
#include <sstream>
#include <iostream>
#include "system.h"
#include "global.h"
#include "moves.h"
#include "atom.h"

class swapParticles : public mcMove {
public:
	swapParticles () {};
	swapParticles (const int typeIndex1, const int typeIndex2, const std::string tag) { typeIndex_ = typeIndex1; typeIndex2_ = typeIndex2; name_ = tag + boost::lexical_cast<std::string>(typeIndex1) + "_" + boost::lexical_cast<std::string>(typeIndex2); } //!< Instantiate a new move, also give a name which is the combination of a user-defined tag + the particle indices it operates on
	int make (simSystem &sys);

private:
	int typeIndex2_; //!< Index of the particle type to swap with typeIndex_
};

#endif