#ifndef RESTART_H_
#define RESTART_H_

#include <string>
#include "system.h"

/*!
 * Information to restart the simulation
 */
class restartInfo {
public:
	restartInfo () {};
	restartInfo (const std::string filename);
    ~restartInfo () {};

	void saveState (const simSystem &sys);
	void readState (const std::string filename);

    bool tmmcDone, crossoverDone, walaDone;
};

#endif
