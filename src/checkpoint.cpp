#include "checkpoint.h"

/*!
 * Read system state from a file.
 *
 * \param [in] dir Directory where system state was saved
 * \param [in] frequency Frquency to take snapshots/checkpoints of the system (< 0 disables)
 */
checkpoint::checkpoint (const std::string directory, const int frequency) {
    tmmcDone = false;
    crossoverDone = false;
    walaDone = false;
    hasCheckpoint = false;
    restartFromTMMC = false;
    restartFromWALA = false;

    dir = directory;
    freq = frequency;

    chkptName = dir+"/state.json";
    if (fileExists(chkptName)) {
        load (chkptName);
    } else {
        std::string command = "mkdir -p "+dir+" && touch "+chkptName;
        system(command.c_str());
    }
}

/*!
 * Read state of a system from a json file.
 *
 * \param [in] filename Filename of where system state was saved
 */
void checkpoint::load (const std::string filename) {
    hasCheckpoint = true;
}

/*!
 * Save the state of a system to a json file.
 *
 * \param [in] filename Filename of where system state is to be saved
 */
void checkpoint::dump (const simSystem &sys) {
    hasCheckpoint = true;
}

/*!
 * Check how long it has been since last checkpoint, and write new one if has exceeded frequency.
 */
 void checkpoint::check (const simSystem &sys) {
     
     dump(sys);
 }
