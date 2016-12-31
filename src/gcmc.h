#ifndef GCMC_H_
#define GCMC_H_

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <exception>
#include <cmath>
#include <iomanip>
#include <memory>
#include <algorithm>

#include "system.h"
#include "utilities.h"
#include "mover.h"
#include "sanity.h"
#include "global.h"

void performGCMC (simSystem &sys, moves *usedMovesEq, moves *usedMovesPr);

#endif
