#include "restart.h"

/*!
 * Read system state from a file.
 *
 * \param [in] filename Filename of where system state was saved
 */
restartInfo::restartInfo (const std::string filename, const int frequency) {
    tmmcDone = false;
    crossoverDone = false;
    walaDone = false;
    hasCheckpoint = false;
    restartFromTMMC = false;
    restartFromWALA = false;

    fname = filename;
    if (frequency < 1) {
        throw customException ("Invalid restart save frequency");
    }
    freq = frequency;

    if (fileExists(fname)) {
        load (fname);
    } else {
        std::vector < std::string > xplode = splitstr(fname, '/');
        std::string xyz = "";
        for (unsigned int i = 0; i < xplode.size()-1; ++i) {
            xyz += xplode[i]+"/";
        }
        std::string command = "mkdir -p "+xyz+" && touch "+fname;
        system(command.c_str());
    }
}

/*!
 * Read state of a system from a json file.
 *
 * \param [in] filename Filename of where system state was saved
 */
void restartInfo::load (const std::string filename) {
    hasCheckpoint = true;
}

/*!
 * Save the state of a system to a json file.
 *
 * \param [in] filename Filename of where system state is to be saved
 */
void restartInfo::dump (const simSystem &sys) {
    hasCheckpoint = true;
}
