#ifndef CROSSOVER_H_
#define CROSSOVER_H_

#include <string>
#include <memory>
#include "system.h"
#include "checkpoint.h"
#include "utilities.h"
#include "mover.h"
#include "sanity.h"

void performCrossover (simSystem &sys, checkpoint &res, moves *usedMovesEq);

#endif
