#ifndef INSERT_H_
#define INSERT_H_

#include <cmath>
#include <sstream>
#include <iostream>
#include "system.h"
#include "global.h"
#include "moves.h"
#include "atom.h"
#include "barrier.h"

class insertParticle : public mcMove {
public:
	insertParticle () { changeN_ = true; }
	insertParticle (const int typeIndex, const std::string tag) { typeIndex_ = typeIndex; name_ = tag + sstr(typeIndex); changeN_ = true; } //!< Instantiate a new move, also give a name which is the combination of auser-defined tag + the particle index it operates on
	int make (simSystem &sys);
};

#endif
