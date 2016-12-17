#ifndef WALA_H_
#define WALA_H_

#include "system.h"
#include "restart.h"
#include "utilities.h"
#include "moves.h"

void performWALA (simSystem &sys, restartInfo &res, const move *usedMovesEq);

#endif
