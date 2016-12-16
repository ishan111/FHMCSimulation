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

/*!
 * Usage: ./binary_name input.json
 */
int main (int argc, char * const argv[]) {
	std::cout << "Beginning simulation at " << getTimeStamp() << std::endl;

	moves *usedMovesEq, *usedMovesPr;
	simSystem sys = initialize (argv[1], usedMovesEq, usedMovesPr); // read json file and create system class

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

	usedMovesPr.print("tmmc_stats.log");

    sys.printSnapshot("final.xyz", "last configuration");
    sys.refine_energy_histogram_bounds();
    sys.printEnergyHistogram("final_eHist");
    sys.refine_pk_histogram_bounds();
    sys.printPkHistogram("final_pkHist");
    sys.printExtMoments("final_extMom");
    sys.getTMMCBias()->print("final", false);

	if (fileExists("restart/sys.state")) {
        delete res;
    }
	delete usedMovesEq;
	delete usedMovesPr;

    time (&rawtime);
 	timeinfo = localtime (&rawtime);
    strftime (timestamp,80,"%d/%m/%Y %H:%M:%S",timeinfo);
    std::cout << "Finished simulation at " << timestamp << std::endl;
	return SAFE_EXIT;
}
