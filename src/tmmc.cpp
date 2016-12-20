#include "tmmc.h"

/*!
 * Perform TMMC stage of simulation
 *
 * \param [in] sys System to simulate
 * \param [in] res Restart/checkpoint information
 * \param [in] usedMovesPr Move class to use
 */
void performTMMC (simSystem &sys, checkpoint &res, moves *usedMovesPr) {
    if (sys.useWALA) {
        throw customException ("WALA not deactivated, cannot proceeed with TMMC");
    }

    std::cout << "Beginning TMMC at " << getTimeStamp() << std::endl;
    res.tmmcDone = false;

	if (res.restartFromTMMC) {
        if (!sys.useTMMC) {
            // if restating, may not have done crossover stage so TMMC not activated yet
            sys.startTMMC (sys.tmmcSweepSize, sys.getTotalM());
        }
		try {
			sys.getTMMCBias()->readC(res.restartFromTMMCFile); // read collection matrix
			sys.getTMMCBias()->calculatePI();
	        std::cout << "Restarted TMMC from collection matrix from " << res.restartFromTMMCFile << std::endl;
		} catch (customException &ce) {
			sys.stopTMMC(); // deallocate
			std::cerr << "Failed to initialize from TMMC collection matrix: " << ce.what() << std::endl;
			exit(SYS_FAILURE);
		}
	}

    const int numSweepSnaps = 100;
	long long int sweep = 0;
	long long int sweepPrint = sys.totalTMMCSweeps, printCounter = 0;
	if (sys.totalTMMCSweeps > numSweepSnaps) {
		sweepPrint /= numSweepSnaps;
	}
	while (sweep < sys.totalTMMCSweeps) {
		bool done = false;
		unsigned long long int counter = 0;
		unsigned long long int checkPoint = sys.tmmcSweepSize*(sys.totNMax() - sys.totNMin() + 1)*3; // how often to check full traversal of collection matrix
		// perform a sweep
		while (!done) {
			try {
				usedMovesPr->makeMove(sys);
			} catch (customException &ce) {
				std::cerr << ce.what() << std::endl;
				exit(SYS_FAILURE);
			}

			// only record properties of the system when it is NOT in an intermediate state
			if (sys.getCurrentM() == 0) {
				sys.recordEnergyHistogram();
				sys.recordPkHistogram();
				sys.recordExtMoments();
			}

			// check if sweep is done
			if (counter%checkPoint == 0) {
				done = sys.getTMMCBias()->checkFullyVisited();
				counter = 0;
			}

            res.check();
			counter++;
		}

		sys.getTMMCBias()->iterateForward (); // reset the counting matrix and increment total sweep number
		sweep++;

		std::cout << "Finished " << sweep << "/" << sys.totalTMMCSweeps << " total TMMC sweeps at " << getTimeStamp() << std::endl;

		// Update biasing function from collection matrix
		sys.getTMMCBias()->calculatePI();

		// Periodically write out checkpoints and report statistics
		if (sweep%sweepPrint == 0) {
			printCounter++;
			sys.getTMMCBias()->print("tmmc-Checkpoint-"+std::to_string(printCounter), true);
			sys.refineEnergyHistogramBounds();
			sys.printEnergyHistogram("eHist-Checkpoint-"+std::to_string(printCounter));
            sys.refinePkHistogramBounds();
            sys.printPkHistogram("pkHist-Checkpoint-"+std::to_string(printCounter));
            sys.printExtMoments("extMom-Checkpoint-"+std::to_string(printCounter));
            usedMovesPr->print("tmmc.stats");
		}
	}

    res.tmmcDone = true;
    usedMovesPr->print("tmmc.stats");
}
