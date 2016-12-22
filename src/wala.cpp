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
    std::cout << "Beginning Wang-Landau portion at " << getTimeStamp() << std::endl;

    bool flat = false;
    double lnF = sys.lnF_start;
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
                std::cerr << ce.what() << std::endl;
                exit(SYS_FAILURE);
            }
            std::cout << "Read initial lnPI for Wang-Landau from " << sys.restartFromWALAFile << std::endl;
        }
    } else {
        lnF = sys.getWALABias()->lnF(); // checkpoint re-initialized to starting value
        moveStart = res.moveCounter;
    }

    std::cout << "Initial lnF = " << sys.getWALABias()->lnF() << " at " << getTimeStamp() << std::endl;

    while (lnF > sys.lnF_end) {
        for (unsigned long long int move = moveStart; move < sys.wlSweepSize; ++move) {
            try {
                usedMovesEq->makeMove(sys);
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
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
            std::cout << "lnF = " << lnF << " at " << getTimeStamp() << std::endl;
        }
    }

    sanityChecks(sys);
    
    res.walaDone = true;
    usedMovesEq->print("wala.stats");
}
