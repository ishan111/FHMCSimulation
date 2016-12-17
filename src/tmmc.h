#ifndef TMMC_H_
#define TMMC_H_

#include "system.h"
#include "restart.h"
#include "utilities.h"
#include "mover.h"

void performTMMC (simSystem &sys, restartInfo &res, moves *usedMovesPr);

#endif
