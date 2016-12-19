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
        const int succ = system(command.c_str());
        if (succ != 0) {
            throw customException("Unable to remove previous checkpoint file");
        }
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
    rapidjson::StringBuffer s;
    rapidjson::PrettyWriter < rapidjson::StringBuffer > writer(s);

    // Write restart/checkpoint options
    writer.StartObject();
    writer.String("tmmcDone");
    writer.Bool(tmmcDone);

    writer.String("crossoverDone");
    writer.Bool(crossoverDone);

    writer.String("walaDone");
    writer.Bool(walaDone);

    writer.String("restartFromTMMC");
    writer.Bool(restartFromTMMC);

    writer.String("restartFromWALA");
    writer.Bool(restartFromWALA);

    writer.String("hasCheckpoint");
    writer.Bool(hasCheckpoint);

    writer.String("takeSnaps");
    writer.Bool(takeSnaps);

    writer.String("freq");
    writer.Int(freq);

    writer.String("dir");
    writer.String(dir.c_str());

    if (walaDone && crossoverDone) {
        // in final TMMC stage or just finished the TMMC (end of simulation)
        sys.getTMMCBias()->print(dir+"/tmmc", true);
        sys.printEnergyHistogram(dir+"/eHist", false); // Un-normalized Energy histogram
        sys.printPkHistogram(dir+"/pkHist", false); // Un-normalized Particle histogram
        sys.printExtMoments(dir+"/extMom", false); // Un-normalized Extensive moments, plus counter (number of times each recorded)
        writer.String("extMomCounter");
        std::vector < double > ctr = sys.extMomCounter();
        writer.StartArray();
        for (std::vector < double >::iterator it = ctr.begin(); it < ctr.end(); ++it) {
            writer.Double(*it);
        }
        writer.EndArray();
    } else if (walaDone && !crossoverDone && !tmmcDone) {
        // in crossover stage
        sys.getTMMCBias()->print(dir+"/tmmc", true);
        sys.getWALABias()->print(dir+"/wala", true);
        // energy upper and lower bounds for histogram
        std::vector < double > elb = sys.getELB(), eub = sys.getEUB();
        writer.String("energyHistogramLB");
        writer.StartArray();
        for (std::vector < double >::iterator it = elb.begin(); it < elb.end(); ++it) {
            writer.Double(*it);
        }
        writer.EndArray();
        writer.String("energyHistogramUB");
        for (std::vector < double >::iterator it = eub.begin(); it < eub.end(); ++it) {
            writer.Double(*it);
        }
        writer.EndArray();
    } else if (!walaDone && !crossoverDone && !tmmcDone) {
        // in WALA stage
        sys.getWALABias()->print(dir+"/wala", true);
        // energy upper and lower bounds for histogram
        std::vector < double > elb = sys.getELB(), eub = sys.getEUB();
        writer.String("energyHistogramLB");
        writer.StartArray();
        for (std::vector < double >::iterator it = elb.begin(); it < elb.end(); ++it) {
            writer.Double(*it);
        }
        writer.EndArray();
        writer.String("energyHistogramUB");
        for (std::vector < double >::iterator it = eub.begin(); it < eub.end(); ++it) {
            writer.Double(*it);
        }
        writer.EndArray();
    } else {
        throw customException ("Uncertain which stage simulation is in, so cannot checkpoint");
    }

    writer.EndObject();
    std::ofstream outData(chkptName.c_str());
    outData << s.GetString() << std::endl;

    if (takeSnaps) {
        // this only prints M = 0 atoms (fully inserted)
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
