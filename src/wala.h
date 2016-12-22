#ifndef WALA_H_
#define WALA_H_

#include <memory>
#include "system.h"
#include "checkpoint.h"
#include "utilities.h"
#include "mover.h"
#include "sanity.h"

void performWALA (simSystem &sys, checkpoint &res, moves *usedMovesEq);

#endif
