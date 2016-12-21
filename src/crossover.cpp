#include "crossover.h"

/*!
 * Perform crossover from Wang-Landau stage of simulation to TMMC
 *
 * \param [in] sys System to simulate
 * \param [in] res Restart/checkpoint information
 * \param [in] usedMovesEq Move class to use
 */
void performCrossover (simSystem &sys, checkpoint &res, moves *usedMovesEq) {
    if (res.crossoverDone) {
        throw customException ("Checkpoint indicates crossover already finished");
    }
    std::cout << "Crossing over to build TMMC matrix at " << getTimeStamp() << std::endl;

    res.crossoverDone = false;
    long long int timesFullyVisited = 0, moveStart = 0;
    if (!res.resFromCross) {
        if (!sys.useWALA) {
            throw customException ("WALA not configured, cannot proceeed with crossover");
        }
        sys.startTMMC (sys.tmmcSweepSize, sys.getTotalM());
    } else {
        timesFullyVisited = res.sweepCounter;
        moveStart = res.moveCounter;
    }

    while (timesFullyVisited < sys.nCrossoverVisits) {
        for (unsigned long long int move = moveStart; move < sys.wlSweepSize; ++move) {
            try {
                usedMovesEq->makeMove(sys);
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
                exit(SYS_FAILURE);
            }
            if (sys.getCurrentM() == 0) {
                sys.checkEnergyHistogramBounds ();
            }
            res.check(sys, move, timesFullyVisited);
        }

        if (sys.getTMMCBias()->checkFullyVisited()) {
            try {
                sys.getTMMCBias()->calculatePI();
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
                sys.getTMMCBias()->print("tmmc-crossover-fail", true);
                sys.getTMMCBias()->dumpVisited("tmmc-crossover-fail-visited");
                exit(SYS_FAILURE);
            }
            sys.getTMMCBias()->iterateForward (); // reset the counting matrix and increment total sweep number
            timesFullyVisited = sys.getTMMCBias()->numSweeps();
            std::cout << "Times C fully visited = " << timesFullyVisited << " at " << getTimeStamp() << std::endl;
        }

        // Check if bias has flattened out, just for continuous improvement
        bool flat = sys.getWALABias()->evaluateFlatness();
        if (flat) {
            // If flat, need to reset H and reduce lnF
            sys.getWALABias()->iterateForward();
            std::cout << "lnF = " << sys.getWALABias()->lnF() << " at " << getTimeStamp() << std::endl;
        }
    }

    // Switch over to TMMC completely
    sys.stopWALA();

    std::cout << "Switching over to TMMC completely, ending Wang-Landau" << std::endl;
    try {
        sys.getTMMCBias()->calculatePI();
        //sys.getTMMCBias()->print("tmmc-beginning-Checkpoint", true);
    } catch (customException &ce) {
        std::cerr << ce.what() << std::endl;
        sys.getTMMCBias()->print("tmmc-beginning-fail", true);
        sys.getTMMCBias()->dumpVisited("tmmc-beginning-fail-visited");
        exit(SYS_FAILURE);
    }

    // if doing initial WL "equilibration" re-initialize the histogram using bounds
    sys.reInitializeEnergyHistogram();

    res.crossoverDone = true;
    usedMovesEq->print("crossover.stats");
}
