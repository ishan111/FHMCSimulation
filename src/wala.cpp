#include "wl.h"

/*!
 * Perform Wang-Landau stage of simulation
 *
 * \param [in] sys System to simulate
 * \param [in] res Restart information
 */
void performWALA (simSystem &sys, restartInfo &res) {
    std::cout << "Beginning Wang-Landau portion" << std::endl;

    char *timestamp;
	time_t rawtime;
	struct tm * timeinfo;

    bool highSnap = false, lowSnap = false;

    // Initially do a WL simulation
    bool flat = false;
    double lnF = lnF_start;
    sys.startWALA (lnF, g, s, sys.getTotalM()); //!< Using Shen and Errington method this syntax is same for single and multicomponent

    time_t rawtime_t;
    time (&rawtime_t);
    struct tm * timeinfo_t;
    timeinfo_t = localtime (&rawtime_t);
    char dummy_t [80];
    strftime (dummy_t,80,"%d/%m/%Y %H:%M:%S",timeinfo_t);
    std::cout << "Initial lnF = " << lnF_start << " at " << dummy_t << std::endl;

    if (restartFromWALA) {
        try {
            sys.getWALABias()->readlnPI(restartFromWALAFile);
        } catch (customException &ce) {
            std::cerr << ce.what() << std::endl;
            exit(SYS_FAILURE);
        }
        std::cout << "Read initial lnPI for Wang-Landau from " << restartFromWALAFile << std::endl;
    }

    long long int counter = 0;
    while (lnF > lnF_end) {
        for (unsigned int move = 0; move < wlSweepSize; ++move) {
            try {
                usedMovesEq.makeMove(sys);
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
                exit(SYS_FAILURE);
            }
            if (sys.getCurrentM() == 0){
                sys.check_energy_histogram_bounds ();
            }
        }

        // Check if bias has flattened out
        flat = sys.getWALABias()->evaluateFlatness();
        if (flat) {
            counter++;

            // Periodically write out checkpoints - before iterateForward() which destroys H matrix
            sys.getWALABias()->print("wl-Checkpoint-"+sstr(counter), true);

            // if flat, need to reset H and reduce lnF
            sys.getWALABias()->iterateForward();
            lnF = sys.getWALABias()->lnF();
            flat = false;

            time_t rawtime_tmp;
            time (&rawtime_tmp);
            struct tm * timeinfo_tmp;
            timeinfo_tmp = localtime (&rawtime_tmp);
            char dummy_tmp [80];
            strftime (dummy_tmp,80,"%d/%m/%Y %H:%M:%S",timeinfo_tmp);
            std::cout << "lnF = " << lnF << " at " << dummy_tmp << std::endl;
        }

        // also check to print out snapshots with 10% of bounds to be used for other restarts
        if (!highSnap) {
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
        }
    }
}
