#ifndef TRANSLATE_H_
#define TRANSLATE_H_

#include <string>
#include <iostream>
#include <sstream>
#include <string>
#include "system.h"
#include "global.h"
#include "moves.h"
#include "atom.h"

class translateParticle : public mcMove {
public:
	translateParticle () { changeN_ = false; }
	translateParticle (const int typeIndex, const std::string tag) { typeIndex_ = typeIndex; name_ = tag + std::to_string(typeIndex); maxD_ = 0.1; changeN_ = false; } //!< Instantiate a new move, also give a name which is the combination of auser-defined tag + the particle index it operates on
	int make (simSystem &sys);
	void setMaxTranslation (const double maxD, const std::vector < double > &box);
	const double getMaxTranslation () { return maxD_; } //!< Return the max translation allowed in a single move

private:
	double maxD_; //!< Maximum translation allowed in a given move, defaults to 0.1
};

#endif
