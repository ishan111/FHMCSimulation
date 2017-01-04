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
    sendMsg("Crossing over to build TMMC matrix");

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

    sendMsg("Starting from lnF = "+numToStr(sys.getWALABias()->lnF()));
    sendMsg("Starting from "+numToStr(moveStart)+" moves in current sweep");
    sendMsg("Starting from "+numToStr(timesFullyVisited)+" out of "+numToStr(sys.nCrossoverVisits)+" sweeps");

    while (timesFullyVisited < sys.nCrossoverVisits) {
        for (long long int move = moveStart; move < sys.wlSweepSize; ++move) {
            try {
                usedMovesEq->makeMove(sys);
                sys.crossoverTotalStepCounter += 1.0;
            } catch (customException &ce) {
                sendErr(ce.what());
                exit(SYS_FAILURE);
            }
            if (sys.getCurrentM() == 0) {
                sys.checkEnergyHistogramBounds ();
            }
            res.check(sys, move, timesFullyVisited, false);
        }

        if (sys.getTMMCBias()->checkFullyVisited()) {
            try {
                sys.getTMMCBias()->calculatePI();
            } catch (customException &ce) {
                sendErr(ce.what());
                sys.getTMMCBias()->print("tmmc-crossover-fail", true);
                sys.getTMMCBias()->dumpVisited("tmmc-crossover-fail-visited");
                exit(SYS_FAILURE);
            }
            sys.getTMMCBias()->iterateForward (); // Reset the counting matrix and increment total sweep number
            timesFullyVisited = sys.getTMMCBias()->numSweeps();
            sendMsg("Times C fully visited = "+numToStr(timesFullyVisited));
            usedMovesEq->print("crossover.stats");
        }

        // Check if bias has flattened out, just for continuous improvement
        bool flat = sys.getWALABias()->evaluateFlatness();
        if (flat) {
            sys.getWALABias()->iterateForward(); // If flat, need to reset H and reduce lnF
            sendMsg("Wang-Landau is now flat, new lnF = "+numToStr(sys.getWALABias()->lnF()));
        }
    }

    // Switch over to TMMC completely
    sendMsg("Switching over to TMMC completely, ending Wang-Landau");
    sys.stopWALA();
    try {
        sys.getTMMCBias()->calculatePI();
    } catch (customException &ce) {
        sendErr(ce.what());
        sys.getTMMCBias()->print("tmmc-beginning-fail", true);
        sys.getTMMCBias()->dumpVisited("tmmc-beginning-fail-visited");
        exit(SYS_FAILURE);
    }

    sys.reInitializeEnergyHistogram(); // If doing initial WL "equilibration" re-initialize the histogram using bounds
    sanityChecks(sys);
    res.crossoverDone = true; // Do not need to dump a checkpoint
    sendMsg("Completed "+numToStr(sys.crossoverTotalStepCounter)+" total MC steps as part of crossover stage");
    sendMsg("Total MC steps taken in simulation: "+numToStr(sys.walaTotalStepCounter+sys.crossoverTotalStepCounter));
}
