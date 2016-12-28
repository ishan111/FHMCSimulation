#ifndef CHECKPOINT_H_
#define CHECKPOINT_H_

#include <string>
#include <vector>
#include <ctime>
#include "global.h"
#include "system.h"
#include "utilities.h"

// JSON interface from local distro of rapidjson
#include "rapidjson/include/rapidjson/document.h"
#include "rapidjson/include/rapidjson/writer.h"
#include "rapidjson/include/rapidjson/stringbuffer.h"
#include "rapidjson/include/rapidjson/filewritestream.h"
#include "rapidjson/include/rapidjson/filereadstream.h"
#include "rapidjson/include/rapidjson/prettywriter.h"

/*!
 * Information to restart/checkpoint the simulation
 */
class checkpoint {
public:
	checkpoint () { tmmcDone = false; crossoverDone = false; walaDone = false; hasCheckpoint = false; resFromTMMC = false; resFromWALA = false; resFromCross = false; takeSnaps = false; freq = -1; dir = "checkpt"; }
	checkpoint (const std::string directory, const long int frequency, simSystem &sys, const bool snaps=false, const bool override=false);
    ~checkpoint () {};

	void load (simSystem &sys, const bool override);
	void dump (simSystem &sys, const long long int moveCounter=0, const long long int sweepCounter=0, const bool refine=true);
	bool check (simSystem &sys, const long long int moveCounter=0, const long long int sweepCounter=0, const bool refine=true);

	bool hasCheckpoint; //!< At least one checkpoint has been made that the system can restart from
	bool takeSnaps; //!< Save snapshot of the system each time a record is made
	bool tmmcDone, crossoverDone, walaDone; //!< Progress of each stage, regardless of where the checkpoint indicated to start from
	bool resFromWALA, resFromCross, resFromTMMC; //!< Flags corresponding to which stage the checkpoint indicated to restart from

	long int freq; //!< Frequency (in seconds) that the system should print a new instantaneous snapshot of itself, does not load from checkpoints but is assigned when instantiated (is dumped though)
	long long int moveCounter; //!< Tracks the number of moves in a given sweep that have executed
	long long int sweepCounter; //!< Tracks the number of sweeps that have executed
	double wala_lnF; //!< Current value of lnF from WALA

	std::string dir; //!< Name of the checkpoint directory containing the information to reinitialize the system (json)
	std::string chkptName; //!< Name of checkpoint file

	std::vector < double > elb, eub; //!< Upper and lower energy bounds for energy histogram

private:
	time_t lastCheckPt_, now_; //!< Time last checkpoint was taken and current time
};

#endif
