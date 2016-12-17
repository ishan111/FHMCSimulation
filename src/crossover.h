#ifndef CROSSOVER_H_
#define CROSSOVER_H_

#include "system.h"
#include "restart.h"
#include "utilities.h"
#include "mover.h"

void performCrossover (simSystem &sys, restartInfo &res, moves *usedMovesEq);

#endif
