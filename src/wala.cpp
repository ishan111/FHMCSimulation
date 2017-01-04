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
    long double total_step_counter = 0;
    long long int moveStart = 0;

    if (!res.resFromWALA) {
        if (sys.useTMMC or sys.useWALA) {
            throw customException ("WALA or TMMC already active, cannot proceeed with WALA");
        }

        sys.startWALA (lnF, sys.wala_g, sys.wala_s, sys.getTotalM());
        if (sys.restartFromWALA) {
            // specified to start WALA from a lnPI guess, and this is not a restart from a checkpoint
            try {
                sys.getWALABias()->readlnPI(sys.restartFromWALAFile);
            } catch (customException &ce) {
                sendErr(ce.what());
                exit(SYS_FAILURE);
            }
            sendMsg("Read initial lnPI for Wang-Landau from "+sys.restartFromWALAFile);
        }
    } else {
        lnF = sys.getWALABias()->lnF(); // checkpoint re-initialized to starting value
        moveStart = res.moveCounter;
    }

    sendMsg("Initial lnF = "+numToStr(sys.getWALABias()->lnF()));
    sendMsg("Starting from "+numToStr(moveStart)+" moves");

    while (lnF > sys.lnF_end) {
        for (unsigned long long int move = moveStart; move < sys.wlSweepSize; ++move) {
            try {
                usedMovesEq->makeMove(sys);
            } catch (customException &ce) {
                sendErr(ce.what());
                exit(SYS_FAILURE);
            }
            if (sys.getCurrentM() == 0){
                sys.checkEnergyHistogramBounds ();
            }
            res.check(sys, move);
        }

        // Check if bias has flattened out
        flat = sys.getWALABias()->evaluateFlatness();
        if (flat) {
            sys.getWALABias()->iterateForward(); // if flat, need to reset H and reduce lnF
            lnF = sys.getWALABias()->lnF();
            flat = false;
            sendMsg("Wang-Landau is now flat, new lnF = "+numToStr(lnF));
            usedMovesEq->print("wala.stats");
        }
    }

    sanityChecks(sys);
    res.walaDone = true;
}
