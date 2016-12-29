#ifndef INPUT_H_
#define INPUT_H_

#include <string>
#include "mover.h"
#include "system.h"
#include "insert.h"
#include "delete.h"
#include "translate.h"
#include "swap.h"
#include "utilities.h"

void checkBounds (simSystem &sys);
void initializeSystemBarriers (simSystem &sys, const rapidjson::Document &doc);
void setPairPotentials (simSystem &sys, const rapidjson::Document &doc);
void setup (simSystem &sys, const std::string filename);
simSystem initialize (const std::string filename, moves* eqMoves, moves* prMoves);

#endif
