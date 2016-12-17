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

    fname = filename;
    assert (frequency >= 1);
    freq = frequency;

    if (fileExists(fname)) {
        load (fname);
    } else {
        mkdir(fname);
    }
}

/*!
 * Read state of a system.
 *
 * \param [in] filename Filename of where system state was saved
 */
void restartInfo::load (const std::string filename) {
    hasCheckpoint = true;
}

/*!
 * Save the state of a system.
 *
 * \param [in] filename Filename of where system state is to be saved
 */
void restartInfo::dump (const simSystem &sys) {
    hasCheckpoint = true;
}
