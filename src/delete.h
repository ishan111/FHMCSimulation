#ifndef DELETE_H_
#define DELETE_H_

#include <iostream>
#include <sstream>
#include <string>
#include "system.h"
#include "global.h"
#include "moves.h"
#include "atom.h"

class deleteParticle : public mcMove {
public:
	deleteParticle () { changeN_ = true; }
    deleteParticle (const int typeIndex, const std::string tag) { typeIndex_ = typeIndex; name_ = tag + boost::lexical_cast<std::string>(typeIndex); } //!< Instantiate a new move, also give a name which is the combination of auser-defined tag + the particle index it operates on
	int make (simSystem &sys);
};

#endif
