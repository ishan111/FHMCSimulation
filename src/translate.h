#ifndef TRANSLATE_H_
#define TRANSLATE_H_

#include <string>
#include <iostream>
#include <sstream>
#include "system.h"
#include "global.h"
#include "moves.h"
#include "atom.h"

class translateParticle : public mcMove {
public:
	translateParticle () { changeN_ = false; }
	translateParticle (const int typeIndex, const std::string tag) { typeIndex_ = typeIndex; name_ = tag + sstr(typeIndex); maxD_ = 0.1; changeN_ = false; } //!< Instantiate a new move, also give a name which is the combination of auser-defined tag + the particle index it operates on
	int make (simSystem &sys);
	void setMaxDisplacement (const double maxD, const std::vector < double > &box);
	const double getMaxDisplacement () { return maxD_; } //!< Return the max displacement allowed in a single move
	
private:
	double maxD_; //!< Maximum displacement allowed in a given move, defaults to 0.1
};

#endif
