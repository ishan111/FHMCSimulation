#include "crossover.h"

/*!
 * Perform crossover from Wang-Landau stage of simulation to TMMC
 *
 * \param [in] sys System to simulate
 * \param [in] res Restart information
 */
void performCrossover (simSystem &sys, const restartInfo &res) {
    std::cout << "Crossing over to build TMMC matrix" << std::endl;

    char *timestamp;
	time_t rawtime;
	struct tm * timeinfo;

    // After a while, combine to initialize TMMC collection matrix
    sys.startTMMC (tmmcSweepSize, sys.getTotalM());

    // actually this should run until all elements of the collection matrix have been populated
    int timesFullyVisited = 0;
    while (timesFullyVisited < nCrossoverVisits) {
        for (unsigned int move = 0; move < wlSweepSize; ++move) {
            try {
                usedMovesEq.makeMove(sys);
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
                exit(SYS_FAILURE);
            }
            if (sys.getCurrentM() == 0) {
                sys.check_energy_histogram_bounds ();
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
            sys.getWALABias()->print("wl-crossover-Checkpoint-"+sstr(timesFullyVisited), true);
            sys.getTMMCBias()->print("tmmc-crossover-Checkpoint-"+sstr(timesFullyVisited), true);

            time_t rawtime_tmp;
            time (&rawtime_tmp);
            struct tm * timeinfo_tmp;
            timeinfo_tmp = localtime (&rawtime_tmp);
            char dummy_tmp [80];
            strftime (dummy_tmp,80,"%d/%m/%Y %H:%M:%S",timeinfo_tmp);
            std::cout << "Times C fully visited = " << timesFullyVisited << " at " << dummy_tmp << std::endl;
        }

        // Check if bias has flattened out
        flat = sys.getWALABias()->evaluateFlatness();
        if (flat) {
            // If flat, need to reset H and reduce lnF
            sys.getWALABias()->iterateForward();

            time_t rawtime_tmp;
            time (&rawtime_tmp);
            struct tm * timeinfo_tmp;
            timeinfo_tmp = localtime (&rawtime_tmp);
            char dummy_tmp [80];
            strftime (dummy_tmp,80,"%d/%m/%Y %H:%M:%S",timeinfo_tmp);
            std::cout << "lnF = " << sys.getWALABias()->lnF() << " at " << dummy_tmp << std::endl;
        }
    }

    // Report move statistics for equilibration stage
    char eq_statName [80];
    strftime (eq_statName,80,"%Y_%m_%d_%H_%M_%S-equilibration-stats.log",timeinfo);
    std::ofstream eq_statFile (eq_statName);
    std::vector < std::vector < double > > eq_stats = usedMovesEq.reportMoveStatistics();
    eq_statFile << " ---------- Move Statistics --------- " << std::endl << " Move\t\% Success" << std::endl;
    for (unsigned int i = 0; i < eq_stats.size(); ++i) {
        double prod = 1.0;
        for (unsigned int j = 0; j < eq_stats[i].size(); ++j) {
            prod *= eq_stats[i][j];
            eq_statFile << usedMovesEq.includedMoves()[i]->myName() << " (from M = " << j << ")\t" << eq_stats[i][j]*100.0 << std::endl;
        }
        if (eq_stats[i].size() > 1) {
            eq_statFile << "-------------------------------------\nProduct of percentages (%) = " << prod*100 << "\n-------------------------------------" << std::endl;
        }
    }
    eq_statFile << std::endl;
    eq_statFile.close();

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
}
