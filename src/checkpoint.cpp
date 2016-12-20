#include "checkpoint.h"

/*!
 * Read system state from a file.  If checkpoint directory is found, data is loaded from it.
 *
 * \param [in] dir Directory where system state was saved
 * \param [in] frequency Frquency to take snapshots/checkpoints of the system (< 0 disables)
 * \param [in] snaps Take snapshots each time a record is made to make a movie? (default = false)
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
 * \param [in] sys System to checkpoint
 */
void checkpoint::load (simSystem &sys) {
    rapidjson::Document doc;
    try {
        FILE* fp = fopen(chkptName.c_str(), "r");
    	char readBuffer[65536];
    	rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
    	doc.ParseStream(is);
    	fclose(fp);

        tmmcDone = doc["tmmcDone"].GetBool();
        crossoverDone = doc["crossoverDone"].GetBool();
        walaDone = doc["walaDone"].GetBool();
        restartFromTMMC = doc["restartFromTMMC"].GetBool();
        restartFromWALA = doc["restartFromWALA"].GetBool();
        hasCheckpoint = doc["hasCheckpoint"].GetBool();
        takeSnaps = doc["takeSnaps"].GetBool();
        freq = doc["freq"].GetInt();
        dir = doc["dir"].GetString();

        if (walaDone && crossoverDone) {
            // in final TMMC stage or just finished the TMMC (end of simulation)
            sys.startTMMC(sys.tmmcSweepSize, sys.getTotalM());
            sys.getTMMCBias()->readC(dir+"/tmmc_C.dat");
            sys.getTMMCBias()->readHC(dir+"/tmmc_HC.dat");
            sys.getTMMCBias()->calculatePI();

            std::vector < double > ctr (doc["extMomCounter"].Size(), 0);
            for (unsigned int i = 0; i < doc["extMomCounter"].Size(); ++i) {
                ctr[i] = doc["extMomCounter"][i].GetDouble();
            }
            sys.restartEnergyHistogram(dir+"/eHist");
            sys.restartPkHistogram(dir+"/pkHist");
            sys.restartExtMoments(dir+"/extMom", ctr);
        } else if (walaDone && !crossoverDone && !tmmcDone) {
            // in crossover stage
            sys.startTMMC(sys.tmmcSweepSize, sys.getTotalM());
            currentLnF = doc["wala_lnF"].GetDouble()
            sys.startWALA (currentLnF, sys.wala_g, sys.wala_s, sys.getTotalM());

            sys.getTMMCBias()->readC(dir+"/tmmc_C.dat");
            sys.getTMMCBias()->readHC(dir+"/tmmc_HC.dat");
            sys.getWALABias()->readlnPI(dir+"/wala_lnPI.dat");
            sys.getWALABias()->readH(dir+"/wala_H.dat");

            // energy upper and lower bounds for histogram
            elb.resize(doc["energyHistogramLB"].Size(), 0);
            for (unsigned int i = 0; i < doc["energyHistogramLB"].Size(); ++i) {
                elb[i] = doc["energyHistogramLB"][i].GetDouble();
            }
            sys.setELB(elb);

            eub.resize(doc["energyHistogramUB"].Size(), 0);
            for (unsigned int i = 0; i < doc["energyHistogramUB"].Size(); ++i) {
                eub[i] = doc["energyHistogramUB"][i].GetDouble();
            }
            sys.setEUB(eub);
        } else if (!walaDone && !crossoverDone && !tmmcDone) {
            // in WALA stage
            currentLnF = doc["wala_lnF"].GetDouble()
            sys.startWALA (currentLnF, sys.wala_g, sys.wala_s, sys.getTotalM());

            sys.getWALABias()->readlnPI(dir+"/wala_lnPI.dat");
            sys.getWALABias()->readH(dir+"/wala_H.dat");

            // energy upper and lower bounds for histogram
            elb.resize(doc["energyHistogramLB"].Size(), 0);
            for (unsigned int i = 0; i < doc["energyHistogramLB"].Size(); ++i) {
                elb[i] = doc["energyHistogramLB"][i].GetDouble();
            }
            sys.setELB(elb);

            eub.resize(doc["energyHistogramUB"].Size(), 0);
            for (unsigned int i = 0; i < doc["energyHistogramUB"].Size(); ++i) {
                eub[i] = doc["energyHistogramUB"][i].GetDouble();
            }
            sys.setEUB(eub);
        } else {
            throw customException ("Uncertain which stage simulation is in, so cannot checkpoint");
        }

        sys.readRestart(dir+"/snap.xyz");
    } catch () {
        throw customException ("Unable to load checkppoint");
    } else {
        hasCheckpoint = true;
    }
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
        sys.getTMMCBias()->print(dir+"/tmmc", true, true);
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
        sys.getTMMCBias()->print(dir+"/tmmc", true, true);
        sys.getWALABias()->print(dir+"/wala", true);

        writer.String("wala_lnF");
        writer.Double(sys.getWALABias()->lnF());

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

    sys.printSnapshot(dir+"/snap.xyz", getTimeStamp(), true); // instantaneous snapshot
    if (takeSnaps) {
        // this only prints M = 0 atoms (fully inserted) to create a movie
        sys.printSnapshot(dir+"/movie.xyz", getTimeStamp(), false);
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
