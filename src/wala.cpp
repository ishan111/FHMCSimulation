#include "wala.h"

/*!
 * Perform Wang-Landau stage of simulation
 *
 * \param [in] sys System to simulate
 * \param [in] res Restart information
 * \param [in] usedMovesEq Move class to use
 */
void performWALA (simSystem &sys, restartInfo &res, moves *usedMovesEq) {
    std::cout << "Beginning Wang-Landau portion at " << getTimeStamp() << std::endl;
    res.walaDone = false;

    //bool highSnap = false, lowSnap = false;

    // Initially do a WL simulation
    bool flat = false;
    double lnF = sys.lnF_start;
    sys.startWALA (lnF, sys.wala_g, sys.wala_s, sys.getTotalM());

    std::cout << "Initial lnF = " << sys.lnF_start << " at " << getTimeStamp() << std::endl;

    if (res.restartFromWALA) {
        try {
            sys.getWALABias()->readlnPI(res.restartFromWALAFile);
        } catch (customException &ce) {
            std::cerr << ce.what() << std::endl;
            exit(SYS_FAILURE);
        }
        std::cout << "Read initial lnPI for Wang-Landau from " << res.restartFromWALAFile << std::endl;
    }

    long long int counter = 0;
    while (lnF > sys.lnF_end) {
        for (unsigned int move = 0; move < sys.wlSweepSize; ++move) {
            try {
                usedMovesEq->makeMove(sys);
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
                exit(SYS_FAILURE);
            }

            if (sys.getCurrentM() == 0){
                sys.checkEnergyHistogramBounds ();
            }
        }

        // Check if bias has flattened out
        flat = sys.getWALABias()->evaluateFlatness();
        if (flat) {
            counter++;

            // Periodically write out checkpoints - before iterateForward() which destroys H matrix
            sys.getWALABias()->print("wl-Checkpoint-"+std::to_string(counter), true);
            sys.getWALABias()->iterateForward(); // if flat, need to reset H and reduce lnF
            lnF = sys.getWALABias()->lnF();
            flat = false;

            std::cout << "lnF = " << lnF << " at " << getTimeStamp() << std::endl;
        }

        // also check to print out snapshots with 10% of bounds to be used for other restarts
        /*if (!highSnap) {
            if (sys.getTotN() > sys.totNMax() - (sys.totNMax()-sys.totNMin())*0.1) {
                sys.printSnapshot("high.xyz", "snapshot near upper bound");
                highSnap = true;
            }
        }
        if (!lowSnap) {
            if (sys.getTotN() < sys.totNMin() + (sys.totNMax()-sys.totNMin())*0.1 && sys.getTotN() > 0) {
                sys.printSnapshot("low.xyz", "snapshot near lower bound");
                lowSnap = true;
            }
        }*/
    }

    res.walaDone = true;
    usedMovesEq->print("wala.stats");
}
