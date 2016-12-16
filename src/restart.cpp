#include "restart.h"

/*!
 * Read system state from a file.
 *
 * \param [in] filename Filename of where system state was saved
 */
restartInfo::restartInfo (const std::string filename) {
    tmmcDone = false;
    crossoverDone = false;
    walaDone = false;

    readState (filename);
}

/*!
 * Read state of a system.
 *
 * \param [in] filename Filename of where system state was saved
 */
void restartInfo::readState (const std::string filename) {
    ;
}
