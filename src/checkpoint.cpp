#include "checkpoint.h"

/*!
 * Read system state from a file.
 *
 * \param [in] dir Directory where system state was saved
 * \param [in] frequency Frquency to take snapshots/checkpoints of the system (< 0 disables)
 * \param [in] snaps Take snapshots each time a record is made? (default = false)
 */
checkpoint::checkpoint (const std::string directory, const int frequency, bool snaps) {
    tmmcDone = false;
    crossoverDone = false;
    walaDone = false;
    hasCheckpoint = false;
    restartFromTMMC = false;
    restartFromWALA = false;
    takeSnaps = snaps;

    dir = directory;
    freq = frequency;

    chkptName = dir+"/state.json";
    if (fileExists(chkptName)) {
        load (chkptName);
    } else {
        std::string command = "mkdir -p "+dir+" && touch "+chkptName;
        system(command.c_str());
    }

    time(&lastCheckPt_); // take time when object was instantiated as initial time
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
 * \param [in] sys System to checkpoint
 */
void checkpoint::dump (simSystem &sys) {
    FILE* fp = fopen(chkptName.c_str(), "r");
    char writeBuffer[65536];
	rapidjson::FileWriteStream os(fp, writeBuffer, sizeof(writeBuffer));

	rapidjson::Document doc;
    rapidjson::Value _tmmcDone(tmmcDone), _crossoverDone(crossoverDone), _walaDone(walaDone);
    rapidjson::Value _restartFromTMMC(restartFromTMMC), _restartFromWALA(restartFromWALA);
    rapidjson::Value _hasCheckpoint(hasCheckpoint), _takeSnaps(takeSnaps);
    rapidjson::Value _freq(freq), _restartFromTMMCFile, _restartFromWALAFile, _dir;

    char bufferT[300], bufferW[300], bufferD[300];
    int len = sprintf(bufferT, "%s", restartFromTMMCFile.c_str());
    _restartFromTMMCFile.SetString(bufferT, len, doc.GetAllocator());
    len = sprintf(bufferW, "%s", restartFromWALAFile.c_str());
    _restartFromWALAFile.SetString(bufferW, len, doc.GetAllocator());
    len = sprintf(bufferD, "%s", dir.c_str());
    _dir.SetString(bufferD, len, doc.GetAllocator());

    // also need info to restore extMom, etc.

    rapidjson::Writer < rapidjson::FileWriteStream > writer(os);
    doc.Accept(writer);
    fclose(fp);

    if (takeSnaps) {
        sys.printSnapshot(dir+"/snaps.xyz", getTimeStamp(), false);
    }

    time(&lastCheckPt_);
    hasCheckpoint = true;
}

/*!
 * Check how long it has been since last checkpoint, and write new one if has exceeded frequency.
 *
 * \param [in] sys System to checkpoint
 */
void checkpoint::check (simSystem &sys) {
    if (freq > 0) {
        if (std::abs(difftime(time(&now), lastCheckPt_)) >= freq) {
            dump(sys);
        }
    }
}
