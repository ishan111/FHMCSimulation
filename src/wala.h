#ifndef WALA_H_
#define WALA_H_

#include "system.h"
#include "checkpoint.h"
#include "utilities.h"
#include "mover.h"

void performWALA (simSystem &sys, checkpoint &res, moves *usedMovesEq);

#endif
