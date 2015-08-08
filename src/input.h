#ifndef INPUT_H_
#define INPUT_H_

#include <vector>
#include <string>
#include <time.h>
#include <utility>
#include <sstream>     
#include <fstream>
#include <stdio.h> 
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include "input.h"
#include "global.h"
#include "utilities.h"

std::vector < std::string > parseInput (int argc, char * const argv[], time_t rawtime);

#endif