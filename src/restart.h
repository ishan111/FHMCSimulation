#ifndef RESTART_H_
#define RESTART_H_

#include <string>
#include <vector>
#include "global.h"
#include "system.h"
#include "utilities.h"

// JSON interface from local distro of rapidjson
#include "rapidjson/include/rapidjson/document.h"
#include "rapidjson/include/rapidjson/writer.h"
#include "rapidjson/include/rapidjson/stringbuffer.h"
#include "rapidjson/include/rapidjson/filereadstream.h"
#include "rapidjson/include/rapidjson/prettywriter.h"

/*!
 * Information to restart the simulation
 */
class restartInfo {
public:
	restartInfo () { tmmcDone = false; crossoverDone = false; walaDone = false; hasCheckpoint = false; restartFromTMMC = false; restartFromWALA = false; };
	restartInfo (const std::string filename, const int frequency);
    ~restartInfo () {};

	void dump (const simSystem &sys);
	void load (const std::string filename);

	bool hasCheckpoint; //!< At least one checkpoint has been made that the system can restart from
    bool tmmcDone, crossoverDone, walaDone; //!< Progress of each stage
	long int freq; //!< Frequency (in seconds) that the system should print a new instantaneous snapshot of itself
	std::string fname; //!< Name of the restart snaphot file containing the information to reinitialize the system (json)

	// Manual TMMC restarts
	bool restartFromTMMC; //!< Flag to signal manual restart from beginning of TMMC
	std::string restartFromTMMCFile; //!< Filename to manual restart from beginning of TMMC with a given C matrix

	// Manual Wang-Landau restarts
	bool restartFromWALA; //!< Flag to signal manual restart from beginning of WALA
	std::string restartFromWALAFile; //!< Filename to manual restart from beginning of WALA with a given lnPI matrix
};

#endif
