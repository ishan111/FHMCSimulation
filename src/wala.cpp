#include "wala.h"

/*!
 * Perform Wang-Landau stage of simulation
 *
 * \param [in] sys System to simulate
 * \param [in] res Restart/checkpoint information
 * \param [in] usedMovesEq Move class to use
 */
void performWALA (simSystem &sys, checkpoint &res, moves *usedMovesEq) {
    if (res.walaDone) {
        throw customException ("Checkpoint indicates WALA already finished");
    }
    sendMsg("Beginning Wang-Landau");

    bool flat = false;
    double lnF = sys.lnF_start;
    long long int moveStart = 0;

    if (!res.resFromWALA) {
        if (sys.useTMMC or sys.useWALA) {
            throw customException ("WALA or TMMC already active, cannot proceeed with WALA");
        }

        sys.startWALA (lnF, sys.wala_g, sys.wala_s, sys.getTotalM());
        if (sys.restartFromWALA) {
            // Specified to start WALA from a lnPI guess, and this is not a restart from a checkpoint
            try {
                sys.getWALABias()->readlnPI(sys.restartFromWALAFile);
            } catch (customException &ce) {
                sendErr(ce.what());
                exit(SYS_FAILURE);
            }
            sendMsg("Read initial lnPI for Wang-Landau from "+sys.restartFromWALAFile);
        }
    } else {
        lnF = sys.getWALABias()->lnF(); // Checkpoint re-initialized to starting value
        moveStart = res.moveCounter; // Loop in current lnF stage
    }

    sendMsg("Initial lnF = "+numToStr(sys.getWALABias()->lnF()));
    sendMsg("Starting from "+numToStr(moveStart)+" moves");

    while (lnF > sys.lnF_end) {
        for (unsigned long long int move = moveStart; move < sys.wlSweepSize; ++move) {
            try {
                usedMovesEq->makeMove(sys);
                sys.walaTotalStepCounter += 1.0;
            } catch (customException &ce) {
                sendErr(ce.what());
                exit(SYS_FAILURE);
            }
            if (sys.getCurrentM() == 0){
                sys.checkEnergyHistogramBounds ();
            }
            res.check(sys, move, 0, false);
        }

        // Check if bias has flattened out
        flat = sys.getWALABias()->evaluateFlatness();
        if (flat) {
            sys.getWALABias()->iterateForward(); // If flat, need to reset H and reduce lnF
            lnF = sys.getWALABias()->lnF();
            flat = false;
            sendMsg("Wang-Landau is now flat, new lnF = "+numToStr(lnF));
            usedMovesEq->print("wala.stats");
        }
    }

    sanityChecks(sys);

    res.walaDone = true;
    res.dump(sys, sys.wlSweepSize, 0, false); // Also dump checkpoint

    sendMsg("Completed "+numToStr(sys.walaTotalStepCounter)+" total MC steps as part of Wang-Landau stage");
    sendMsg("Total MC steps taken in simulation: "+numToStr(sys.walaTotalStepCounter));
}
