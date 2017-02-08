#include "checkpoint.h"

/*!
 * Read system state from a file.  If checkpoint directory is found and json file is valid, data is loaded from it.
 * Note that although this stores the frequency, it does not use the value it dumps to file.
 * This class always uses the frequency given when the class is instantiated.
 *
 * \param [in] dir Directory where system state was saved
 * \param [in] frequency Frquency to take snapshots/checkpoints of the system (< 0 disables)
 * \param [in] sys System to checkpoint
 * \param [in] snaps Take snapshots each time a record is made to make a movie? (default = false)
 * \param [in] override Manually override exceptions, use with extreme caution (default=false)
 */
checkpoint::checkpoint (const std::string directory, const long int frequency, simSystem &sys, const bool snaps, const bool override) {
    tmmcDone = false;
    crossoverDone = false;
    walaDone = false;
    hasCheckpoint = false;
    takeSnaps = snaps;
    dir = directory;
    freq = frequency;
    moveCounter = 0;
    sweepCounter = 0;
    resFromWALA = false;
    resFromCross = false;
    resFromTMMC = false;

    chkptName = dir+"/state.json";
    if (fileExists(chkptName)) {
        try {
            load(sys, override);
        } catch (std::exception &ex) {
            std::string a = "Unable to load checkpoint "+chkptName+" : ";
            std::string b = ex.what();
            throw customException (a+b);
        }
    } else {
        // Forcible skip to TMMC stage if want to manually start TMMC
        if (sys.restartFromTMMC){
            walaDone = true;
            crossoverDone = true;
        }
    }

    time(&lastCheckPt_); // Take time when object was instantiated as initial time so that check() has a point of reference
}

/*!
 * Read state of a system from a json file.
 *
 * \param [in] sys System to checkpoint
 * \param [in] override Manually override exceptions, use with extreme caution (default=false)
 */
void checkpoint::load (simSystem &sys, const bool override) {
    if (!fileExists(chkptName) && !override) {
        throw customException ("No checkpoint by the name: "+chkptName);
    }

    rapidjson::Document doc;
    try {
        parseJson (chkptName, doc);

        tmmcDone = doc["tmmcDone"].GetBool();
        crossoverDone = doc["crossoverDone"].GetBool();
        walaDone = doc["walaDone"].GetBool();
        hasCheckpoint = doc["hasCheckpoint"].GetBool();
        takeSnaps = doc["takeSnaps"].GetBool();
        dir = doc["dir"].GetString();
        moveCounter = (long long int)doc["moveCounter"].GetDouble();
        sweepCounter = (long long int)doc["sweepCounter"].GetDouble();
        sys.walaTotalStepCounter = (long double)doc["walaTotalStepCounter"].GetDouble();
    	sys.crossoverTotalStepCounter = (long double)doc["crossoverTotalStepCounter"].GetDouble();
    	sys.tmmcTotalStepCounter = (long double)doc["tmmcTotalStepCounter"].GetDouble();

        if (walaDone && crossoverDone) { // In final TMMC stage or just finished the TMMC (end of simulation)
            resFromTMMC = true;
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
        } else if (walaDone && !crossoverDone && !tmmcDone) { // In crossover stage
            resFromCross = true;
            sys.startTMMC(sys.tmmcSweepSize, sys.getTotalM());
			for (unsigned int i = 0; i < sweepCounter; ++i) sys.getTMMCBias()->iterateForward (); // Have to iterate forward to "catch up" the internal counter

            wala_lnF = doc["wala_lnF"].GetDouble();
            sys.startWALA (wala_lnF, sys.wala_g, sys.wala_s, sys.getTotalM()); // Starting WALA from wala_lnF is already "caught up" so no need to iterate forward

            sys.getTMMCBias()->readC(dir+"/tmmc_C.dat");
            sys.getTMMCBias()->readHC(dir+"/tmmc_HC.dat");
            sys.getWALABias()->readlnPI(dir+"/wala_lnPI.dat");
            sys.getWALABias()->readH(dir+"/wala_H.dat");

            // Energy upper and lower bounds for histogram
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
        } else if (!walaDone && !crossoverDone && !tmmcDone) { // In WALA stage and printed 1st checkpoint already
            resFromWALA = true;
            wala_lnF = doc["wala_lnF"].GetDouble();
            sys.startWALA (wala_lnF, sys.wala_g, sys.wala_s, sys.getTotalM());

            sys.getWALABias()->readlnPI(dir+"/wala_lnPI.dat");
            sys.getWALABias()->readH(dir+"/wala_H.dat");

            // Energy upper and lower bounds for histogram
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
            if (!override) {
                sendErr("Uncertain which stage simulation is in, so cannot checkpoint");
                exit(SYS_FAILURE);
            }
        }

        sys.readConfig(dir+"/snap.xyz");
        hasCheckpoint = true;
    } catch (std::exception &ex) {
        std::string msg = ex.what();
        if (!override) {
            hasCheckpoint = false;
            sendErr("Unable to load checkpoint "+msg);
            exit(SYS_FAILURE);
        } else {
            sendErr("Overriding the following errors to load checkpoint "+msg);
        }
    }

    sendMsg("Checkpoint loaded from "+chkptName);
}

/*!
 * Save the state of a system to a json file. Creates the checkpoint directory if it doesn't exist.
 *
 * \param [in] sys System to checkpoint
 * \param [in] moveCounter Number of moves out of a given sweep that have executed
 * \param [in] sweepCounter Number of loops/sweeps that have executed
 * \param [in] refine Refine the histogram boundaries before printing any? (default=false)
 */
void checkpoint::dump (simSystem &sys, const long long int moveCounter, const long long int sweepCounter, const bool refine) {
    if (!fileExists(chkptName)) {
        std::string command = "mkdir -p "+dir+" && touch "+chkptName;
        const int succ = system(command.c_str());
        if (succ != 0) {
            throw customException("Unable to initialize checkpoint");
        }
    }

    rapidjson::StringBuffer s;
    rapidjson::PrettyWriter < rapidjson::StringBuffer > writer(s);
    hasCheckpoint = true;

    // Write restart/checkpoint options
    writer.StartObject();
    writer.String("tmmcDone");
    writer.Bool(tmmcDone);

    writer.String("crossoverDone");
    writer.Bool(crossoverDone);

    writer.String("walaDone");
    writer.Bool(walaDone);

    writer.String("hasCheckpoint");
    writer.Bool(hasCheckpoint);

    writer.String("takeSnaps");
    writer.Bool(takeSnaps);

    writer.String("freq");
    writer.Int64(freq);

    writer.String("dir");
    writer.String(dir.c_str());

    writer.String("moveCounter");
    writer.Double(moveCounter);

    writer.String("sweepCounter");
    writer.Double(sweepCounter);

    writer.String("walaTotalStepCounter");
    writer.Double(sys.walaTotalStepCounter);

    writer.String("crossoverTotalStepCounter");
    writer.Double(sys.crossoverTotalStepCounter);

    writer.String("tmmcTotalStepCounter");
    writer.Double(sys.tmmcTotalStepCounter);

    if (walaDone && crossoverDone) { // In final TMMC stage or just finished the TMMC (end of simulation)
        sys.getTMMCBias()->print(dir+"/tmmc", true, true);
        if (refine) {
            sys.refineEnergyHistogramBounds();
        }
        sys.printEnergyHistogram(dir+"/eHist", false); // Un-normalized Energy histogram
        if (refine) {
            sys.refinePkHistogramBounds();
        }
        sys.printPkHistogram(dir+"/pkHist", false); // Un-normalized Particle histogram
        sys.printExtMoments(dir+"/extMom", false); // Un-normalized Extensive moments, plus counter (number of times each recorded)
        writer.String("extMomCounter");
        std::vector < double > ctr = sys.extMomCounter();
        writer.StartArray();
        for (std::vector < double >::iterator it = ctr.begin(); it < ctr.end(); ++it) {
            writer.Double(*it);
        }
        writer.EndArray();
    } else if (walaDone && !crossoverDone && !tmmcDone) { // In crossover stage
        sys.getTMMCBias()->print(dir+"/tmmc", true, true);
        sys.getWALABias()->print(dir+"/wala", true);

        writer.String("wala_lnF");
        writer.Double(sys.getWALABias()->lnF());

        // Energy upper and lower bounds for histogram
        std::vector < double > elb = sys.getELB(), eub = sys.getEUB();
        writer.String("energyHistogramLB");
        writer.StartArray();
        for (std::vector < double >::iterator it = elb.begin(); it < elb.end(); ++it) {
            writer.Double(*it);
        }
        writer.EndArray();
        writer.String("energyHistogramUB");
        writer.StartArray();
        for (std::vector < double >::iterator it = eub.begin(); it < eub.end(); ++it) {
            writer.Double(*it);
        }
        writer.EndArray();
    } else if (!walaDone && !crossoverDone && !tmmcDone) { // In WALA stage
        sys.getWALABias()->print(dir+"/wala", true);

        writer.String("wala_lnF");
        writer.Double(sys.getWALABias()->lnF());

        // Energy upper and lower bounds for histogram
        std::vector < double > elb = sys.getELB(), eub = sys.getEUB();
        writer.String("energyHistogramLB");
        writer.StartArray();
        for (std::vector < double >::iterator it = elb.begin(); it < elb.end(); ++it) {
            writer.Double(*it);
        }
        writer.EndArray();
        writer.String("energyHistogramUB");
        writer.StartArray();
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

    sys.printSnapshot(dir+"/snap.xyz", getTimeStamp(), true); // Instantaneous snapshot
    if (takeSnaps) { // This only prints M = 0 atoms (fully inserted) to create a movie
        sys.printSnapshot(dir+"/movie.xyz", getTimeStamp(), false);
    }

    time(&lastCheckPt_);
}

/*!
 * Check how long it has been since last checkpoint, and write new one if has exceeded frequency.
 *
 * \param [in] sys System to checkpoint
 * \param [in] moveCounter Number of moves out of a given sweep that have executed
 * \param [in] sweepCounter Number of loops/sweeps that have executed
 * \param [in] refine Refine the histogram boundaries before printing any? (default=false)
 *
 * \returns bool Is a checkpoint being generated or not
 */
bool checkpoint::check (simSystem &sys, const long long int moveCounter, const long long int sweepCounter, const bool refine) {
    if (freq > 0) {
        if (std::abs(difftime(time(&now_), lastCheckPt_)) >= freq) {
            dump(sys, moveCounter, sweepCounter, refine);
            return true;
        }
    }
    return false;
}
