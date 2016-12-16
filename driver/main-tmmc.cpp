/*!
 * Perform (un/)biased multicomponent GCMC.
 *
 * \author Nathan A. Mahynski
 * \date 08/07/15
 *
 * \mainpage
 *
 * \section Dependencies
 *
 * \section Compiling
 * Run tests and run main
 *
 * Complile with -DFLUID_PHASE_SIMULATIONS if simulations are purely in the fluid phase.
 * This will allow tail corrections to be enabled which are only valid assuming a converging g(r) at large r.
 *
 * Compile with -DNETCDF_CAPABLE if netCDF libraries are installed and can be compiled against.  Certain data will be output to these arrays
 * instead of ASCII files if so.
 *
 * Tests compile with cmake and by default, do not link to netcdf libraries since they do not use them.  Everything tested in ASCII.
 *
 * \section Input
 *
 */

#include <iostream>
#include <time.h>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "../src/system.h"
#include "../src/utilities.h"
#include "../src/global.h"
#include "../src/insert.h"
#include "../src/delete.h"
#include "../src/translate.h"
#include "../src/swap.h"
#include "../src/moves.h"
#include "../src/barrier.h"
#include "../src/input.h"
#include "../src/restart.h"

// JSON interface from local distro of rapidjson
/*#include "../src/rapidjson/include/rapidjson/document.h"
#include "../src/rapidjson/include/rapidjson/writer.h"
#include "../src/rapidjson/include/rapidjson/stringbuffer.h"
#include "../src/rapidjson/include/rapidjson/filereadstream.h"
#include "../src/rapidjson/include/rapidjson/prettywriter.h"*/

/*!
 * Usage: ./binary_name input.json
 */
int main (int argc, char * const argv[]) {
	std::cout << "Beginning simulation at " << getTimeStamp() << std::endl;

	simSystem sys = initialize (argv[1]); // read json file and create system class

	// if restart/ exists, default to use this information to restart the simulation
	restartInfo *res;
    if (fileExists("restart/sys.state")) {
        res = new restartInfo ("restart/sys.state");
    } else {
        // otherwise generate initial configuarion
        setup (sys, argv[1]);
    }

	if (!res->walaDone) {
		// perform Wang-Landau simulation
		performWALA (sys, res);
		performCrossover (sys, res);
		performTMMC (sys, res);
	} else if (!res->crossoverDone) {
		// crossover to TMMC
		performCrossover (sys, res);
		performTMMC (sys, res);
	} else if (!res->tmmcDone) {
		// perform tmmc portion of the simulation
		performTMMC (sys, res);
	}








    // Report move statistics for final TMMC ("production") stage
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

 	// print out restart file (xyz)
    sys.printSnapshot("final.xyz", "last configuration");

    // Print out energy histogram at each Ntot
    sys.refine_energy_histogram_bounds();
    sys.printEnergyHistogram("final_eHist");

    // Print out pk number histogram at Ntot
    sys.refine_pk_histogram_bounds();
    sys.printPkHistogram("final_pkHist");

    // Print out extensive moments
    sys.printExtMoments("final_extMom");

    // Print out final macrostate distribution
    sys.getTMMCBias()->print("final", false);


    // Finished
    time (&rawtime);
 	timeinfo = localtime (&rawtime);
    strftime (timestamp,80,"%d/%m/%Y %H:%M:%S",timeinfo);
    std::cout << "Finished simulation at " << timestamp << std::endl;
	return SAFE_EXIT;

	if (fileExists("restart/sys.state")) {
        delete res;
    }
}
