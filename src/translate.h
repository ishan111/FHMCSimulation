#ifndef TRANSLATE_H_
#define TRANSLATE_H_

#include "system.h"
#include "global.h"
#include "moves.h"
#include <string>

class translateParticle : public mcMove {
public:
	translateParticle () {};
    translateParticle (const int typeIndex, const std::string tag) { typeIndex_ = typeIndex; name_ = tag + boost::lexical_cast<std::string>(typeIndex); } //!< Instantiate a new move, also give a name which is the combination of auser-defined tag + the particle index it operates on
	int make (simSystem &sys);
};

#endif
