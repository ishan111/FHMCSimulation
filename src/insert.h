#ifndef INSERT_H_
#define INSERT_H_

#include <cmath>
#include <sstream>
#include <iostream>
#include "system.h"
#include "global.h"
#include "moves.h"
#include "atom.h"

class insertParticle : public mcMove {
public:
	insertParticle () {};
    insertParticle (const int typeIndex, const std::string tag) { typeIndex_ = typeIndex; name_ = tag + boost::lexical_cast<std::string>(typeIndex); } //!< Instantiate a new move, also give a name which is the combination of auser-defined tag + the particle index it operates on
	int make (simSystem &sys);
};

#endif