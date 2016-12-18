#ifndef SWAP_H_
#define SWAP_H_

#include <cmath>
#include <sstream>
#include <string>
#include <iostream>
#include "system.h"
#include "global.h"
#include "moves.h"
#include "atom.h"
#include "barrier.h"

class swapParticles : public mcMove {
public:
	swapParticles () { changeN_ = false; }
	swapParticles (const int typeIndex1, const int typeIndex2, const std::string tag) { typeIndex_ = typeIndex1; typeIndex2_ = typeIndex2; name_ = tag + std::to_string(typeIndex1) + "_" + std::to_string(typeIndex2); changeN_ = false; } //!< Instantiate a new move, also give a name which is the combination of a user-defined tag + the particle indices it operates on
	int make (simSystem &sys);

private:
	int typeIndex2_; //!< Index of the particle type to swap with typeIndex_
};

#endif
