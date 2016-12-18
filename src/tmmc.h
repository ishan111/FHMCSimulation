#ifndef TMMC_H_
#define TMMC_H_

#include "system.h"
#include "checkpoint.h"
#include "utilities.h"
#include "mover.h"

void performTMMC (simSystem &sys, checkpoint &res, moves *usedMovesPr);

#endif
