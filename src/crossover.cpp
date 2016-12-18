#include "crossover.h"

/*!
 * Perform crossover from Wang-Landau stage of simulation to TMMC
 *
 * \param [in] sys System to simulate
 * \param [in] res Restart/checkpoint information
 * \param [in] usedMovesEq Move class to use
 */
void performCrossover (simSystem &sys, checkpoint &res, moves *usedMovesEq) {
    std::cout << "Crossing over to build TMMC matrix at " << getTimeStamp() << std::endl;
    res.crossoverDone = false;

    // After a while, combine to initialize TMMC collection matrix
    sys.startTMMC (sys.tmmcSweepSize, sys.getTotalM());

    // actually this should run until all elements of the collection matrix have been populated
    int timesFullyVisited = 0;
    while (timesFullyVisited < sys.nCrossoverVisits) {
        for (unsigned int move = 0; move < sys.wlSweepSize; ++move) {
            try {
                usedMovesEq->makeMove(sys);
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
                exit(SYS_FAILURE);
            }
            if (sys.getCurrentM() == 0) {
                sys.checkEnergyHistogramBounds ();
            }
        }

        // Check if collection matrix is ready to take over, not necessarily at points where WL is flat
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
            sys.getWALABias()->print("wl-crossover-Checkpoint-"+std::to_string(timesFullyVisited), true);
            sys.getTMMCBias()->print("tmmc-crossover-Checkpoint-"+std::to_string(timesFullyVisited), true);

            std::cout << "Times C fully visited = " << timesFullyVisited << " at " << getTimeStamp() << std::endl;
        }

        // Check if bias has flattened out
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
        sys.getTMMCBias()->print("tmmc-beginning-Checkpoint", true);
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
