#include "tmmc.h"

/*!
 * Perform TMMC stage of simulation
 *
 * \param [in] sys System to simulate
 * \param [in] res Restart information
 */
void performTMMC (simSystem &sys, const restartInfo &res) {
    std::cout << "Beginning TMMC" << std::endl;
	if (restartFromTMMC) {
		sys.startTMMC (tmmcSweepSize, sys.getTotalM()); // this was otherwise started during the crossover phase if WL was used
		try {
			sys.getTMMCBias()->readC(restartFromTMMCFile); // read collection matrix
			sys.getTMMCBias()->calculatePI();
	        std::cout << "Restarted TMMC from collection matrix from " << restartFromTMMCFile << std::endl;
		} catch (customException &ce) {
			sys.stopTMMC(); // deallocate
			std::cerr << "Failed to initialize from TMMC collection matrix: " << ce.what() << std::endl;
			exit(SYS_FAILURE);
		}
	}

	long long int sweep = 0;
	int sweepPrint = totalTMMCSweeps, numSweepSnaps = 100, printCounter = 0;
	if (totalTMMCSweeps > numSweepSnaps) {
		sweepPrint /= numSweepSnaps;
	}
	while (sweep < totalTMMCSweeps) {
		bool done = false;
		unsigned long long int counter = 0;
		unsigned long long int checkPoint = tmmcSweepSize*(sys.totNMax() - sys.totNMin() + 1)*3; // how often to check full traversal of collection matrix
		// perform a sweep
		while (!done) {
			try {
				usedMovesPr.makeMove(sys);
			} catch (customException &ce) {
				std::cerr << ce.what() << std::endl;
				exit(SYS_FAILURE);
			}

			// only record properties of the system when it is NOT in an intermediate state
			if (sys.getCurrentM() == 0) {
				// record energy histogram at each Ntot
				sys.recordEnergyHistogram();

				// record N1, N2, etc. histogram at each Ntot
				sys.recordPkHistogram();

				// record extensive moments
				sys.recordExtMoments();
			}

			// check if sweep is done
			if (counter%checkPoint == 0) {
				done = sys.getTMMCBias()->checkFullyVisited();
				counter = 0;
			}

			counter++;
		}

		sys.getTMMCBias()->iterateForward (); // reset the counting matrix and increment total sweep number
		sweep++;

		time_t rawtime_tmp;
		time (&rawtime_tmp);
		struct tm * timeinfo_tmp;
		timeinfo_tmp = localtime (&rawtime_tmp);
		char dummy_tmp [80];
		strftime (dummy_tmp,80,"%d/%m/%Y %H:%M:%S",timeinfo_tmp);
		std::cout << "Finished " << sweep << "/" << totalTMMCSweeps << " total TMMC sweeps at " << dummy_tmp << std::endl;

		// Update biasing function from collection matrix
		sys.getTMMCBias()->calculatePI();

		// Periodically write out checkpoints and report statistics
		if (sweep%sweepPrint == 0) {
			printCounter++;
			sys.getTMMCBias()->print("tmmc-Checkpoint-"+sstr(printCounter), true);
			sys.refine_energy_histogram_bounds();
			sys.printEnergyHistogram("eHist-Checkpoint-"+sstr(printCounter));
            sys.refine_pk_histogram_bounds();
            sys.printPkHistogram("pkHist-Checkpoint-"+sstr(printCounter));
            sys.printExtMoments("extMom-Checkpoint-"+sstr(printCounter));

            char statName [80];
            strftime (statName,80,"%Y_%m_%d_%H_%M_%S-stats.log",timeinfo);
            std::ofstream statFile (statName);
            std::vector < std::vector < double > > stats = usedMovesPr.reportMoveStatistics();
            statFile << " ---------- Move Statistics --------- " << std::endl << " Move\t\% Success" << std::endl;
            for (unsigned int i = 0; i < stats.size(); ++i) {
                double prod = 1.0;
                for (unsigned int j = 0; j < stats[i].size(); ++j) {
                    prod *= stats[i][j];
                    statFile << usedMovesPr.includedMoves()[i]->myName() << " (from M = " << j << ")\t" << stats[i][j]*100.0 << std::endl;
                }
                if (stats[i].size() > 1) {
                    statFile << "-------------------------------------\nProduct of percentages (%) = " << prod*100 << "\n-------------------------------------" << std::endl;
                }
            }
            statFile << std::endl;
            statFile.close();
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
