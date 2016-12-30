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

simSystem initialize (const std::string filename, moves* eqMoves, moves* prMoves);
void setConfig (simSystem &sys, const std::string filename);

void setMoves (simSystem &sys, const rapidjson::Document &doc, moves* usedMovesEq, moves* usedMovesPr);
void setBarriers (simSystem &sys, const rapidjson::Document &doc);
void setPairPotentials (simSystem &sys, const rapidjson::Document &doc);
void checkBounds (simSystem &sys);

#endif
