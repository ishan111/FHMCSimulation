#include "tmmc.h"

/*!
 * Perform TMMC stage of simulation
 *
 * \param [in] sys System to simulate
 * \param [in] res Restart/checkpoint information
 * \param [in] usedMovesPr Move class to use
 */
void performTMMC (simSystem &sys, checkpoint &res, moves *usedMovesPr) {
    if (res.tmmcDone) {
        throw customException ("Checkpoint indicates TMMC already finished");
    }
    sendMsg("Beginning TMMC");

    // Specifically for printing progress every 1% of the simulation
    const long long int numSweepSnaps = 100;
	unsigned long long int sweepPrint = sys.totalTMMCSweeps, printCounter = 0;
	if (sys.totalTMMCSweeps >= numSweepSnaps) {
		sweepPrint /= numSweepSnaps;
	}

    res.tmmcDone = false;
    unsigned long long int sweep = 0;
    unsigned long long int counter = 0;
    if (!res.resFromTMMC) {
        if (sys.useWALA) {
            throw customException ("WALA not deactivated, cannot proceeed with TMMC");
        }
        if (sys.restartFromTMMC) {
            if (!sys.useTMMC) {
                // In case where not restarting from checkpoint, but starting fresh from TMMC, crossover stage was skipped so TMMC has not been activated yet
                sys.startTMMC (sys.tmmcSweepSize, sys.getTotalM());
            }
    		try {
    			sys.getTMMCBias()->readC(sys.restartFromTMMCFile); // read collection matrix
    			sys.getTMMCBias()->calculatePI();
                sendMsg("Restarted TMMC from collection matrix from "+sys.restartFromTMMCFile);
    		} catch (customException &ce) {
    			sys.stopTMMC(); // Deallocate
                std::string msg = ce.what();
                sendErr("Failed to initialize from TMMC collection matrix : "+msg);
    			exit(SYS_FAILURE);
    		}
    	}
    } else {
        printCounter = res.moveCounter;
        sweep = res.sweepCounter;
    }

    sendMsg("Starting progress stage from "+numToStr(printCounter)+"/"+numToStr(std::min(numSweepSnaps, sys.totalTMMCSweeps)));
    sendMsg("Starting from "+numToStr(sweep)+"/"+numToStr(sys.totalTMMCSweeps)+" total TMMC sweeps");

    unsigned long long int checkPoint = sys.tmmcSweepSize*(sys.totNMax() - sys.totNMin() + 1)*3; // how often to check full traversal of collection matrix
	while (sweep < sys.totalTMMCSweeps) {
		bool done = false;

		// Perform a sweep
		while (!done) {
			try {
				usedMovesPr->makeMove(sys);
                sys.tmmcTotalStepCounter += 1.0;
			} catch (customException &ce) {
                sendErr(ce.what());
				exit(SYS_FAILURE);
			}

			// Only record properties of the system when it is NOT in an intermediate state
			if (sys.getCurrentM() == 0) {
				sys.recordEnergyHistogram();
				sys.recordPkHistogram();
				sys.recordExtMoments();
			}

			// Check if sweep is done
			if (counter%checkPoint == 0) {
				done = sys.getTMMCBias()->checkFullyVisited();
				counter = 0;
			}

			counter++;
            res.check(sys, printCounter, sweep, false);
		}

		sys.getTMMCBias()->iterateForward(); // Reset the counting matrix and increment total sweep number
		sweep++;

        sendMsg("Finished "+numToStr(sweep)+"/"+numToStr(sys.totalTMMCSweeps)+" total TMMC sweeps");

		sys.getTMMCBias()->calculatePI(); // Update biasing function from collection matrix

		// Periodically write out checkpoints to monitor convergence properties later - all are used in FHMCAnalysis at this point (12/22/16)
		if (sweep%sweepPrint == 0) {
            try {
                printCounter++;
    			sys.getTMMCBias()->print("tmmc-Checkpoint-"+numToStr(printCounter), false, false); // true, false);
                sys.refineEnergyHistogramBounds();
    			sys.printEnergyHistogram("eHist-Checkpoint-"+numToStr(printCounter));
                sys.refinePkHistogramBounds();
                sys.printPkHistogram("pkHist-Checkpoint-"+numToStr(printCounter));
                sys.printExtMoments("extMom-Checkpoint-"+numToStr(printCounter));
                usedMovesPr->print("tmmc.stats");
            } catch (std::exception &ex) {
                const std::string msg = ex.what();
                throw customException ("Unable to print checkpoint : "+msg);
            }
		}
	}

    // Print final results
    try {
        sys.getTMMCBias()->print("final", false, false);
        sys.refineEnergyHistogramBounds();
        sys.printEnergyHistogram("final_eHist");
        sys.refinePkHistogramBounds();
        sys.printPkHistogram("final_pkHist");
        sys.printExtMoments("final_extMom");
        sys.printSnapshot("final.xyz", "last configuration");
        usedMovesPr->print("tmmc.stats");
        sanityChecks(sys);
    } catch (std::exception &ex) {
        const std::string msg = ex.what();
        throw customException ("Unable to print final TMMC results : "+msg);
    }

    res.tmmcDone = true; // Do not need to dump a checkpoint
    sendMsg("Completed "+numToStr(sys.tmmcTotalStepCounter)+" total MC steps as part of TMMC stage");
    sendMsg("Total MC steps taken in simulation: "+numToStr(sys.walaTotalStepCounter+sys.crossoverTotalStepCounter+sys.tmmcTotalStepCounter));
}
