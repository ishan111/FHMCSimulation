#ifndef RESTART_H_
#define RESTART_H_

#include <string>
#include <direct.h>
#include "system.h"
#include "utilities.h"

/*!
 * Information to restart the simulation
 */
class restartInfo {
public:
	restartInfo () {};
	restartInfo (const std::string filename, const int frequency);
    ~restartInfo () {};

	void dump (const simSystem &sys);
	void load (const std::string filename);

    bool tmmcDone, crossoverDone, walaDone;
	int freq;
	std::string fname;
};

#endif
