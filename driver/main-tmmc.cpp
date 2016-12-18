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


#include "../src/fhmc.h"

/*!
 * Usage: ./binary_name input.json
 */
int main (int argc, char * const argv[]) {
	const std::string restartFile = "restart/state.json";
	std::cout << "Beginning simulation at " << getTimeStamp() << std::endl;

	moves usedMovesEq, usedMovesPr;
	simSystem sys = initialize (argv[1], &usedMovesEq, &usedMovesPr); // read json file and create system class

	// usedMovesEq, usedMovesPr need to be global variables so they are static in memory (or atleast the vectors of moves need to be)
	// pair potentials also need to be static somewhow - compare against gcmc (old) driver
	// maybe used a reserve statement somewhow
	// shared pointers? -->  must move to heap!! example: http://stackoverflow.com/questions/23060406/does-stdmove-invalidate-pointers

	// in mover, make moves_ --> shared pointer to the vectors of each move type?
	// in system, ppot needs to be share_ptr? to potential classes, which should be made_shared when created by addPotential()
	// wl, tmmc biases in system? consider replacing new/delete with just simple asigniments so tmmc* tmmcBias; --> tmmc tmmcBias; e.g.
	//http://stackoverflow.com/questions/25405034/stdshared-ptr-of-abstract-class-to-instantiate-derived-class

	// cellLists.... -> rely on atoms not moving in system, preallocated to max_atoms, cellLists are reserved so ok

	// 1. start with making ppot's shared_ptr so they can be stored locally in system

	// if restart/ exists, default to use this information to restart the simulation
	restartInfo res (restartFile, 900);
    if (!res.hasCheckpoint) {
        setup (sys, argv[1]);
    }

	// TODO: finish adding restarting information so system knows how to pick up where it left off

	if (!res.walaDone) {
		// perform Wang-Landau simulation
		performWALA (sys, res, &usedMovesEq);
		performCrossover (sys, res, &usedMovesEq);
		performTMMC (sys, res, &usedMovesPr);
	} else if (!res.crossoverDone) {
		// crossover to TMMC
		performCrossover (sys, res, &usedMovesEq);
		performTMMC (sys, res, &usedMovesPr);
	} else if (!res.tmmcDone) {
		// perform tmmc portion of the simulation
		performTMMC (sys, res, &usedMovesPr);
	}

    sys.printSnapshot("final.xyz", "last configuration");
    sys.refineEnergyHistogramBounds();
    sys.printEnergyHistogram("final_eHist");
    sys.refinePkHistogramBounds();
    sys.printPkHistogram("final_pkHist");
    sys.printExtMoments("final_extMom");
    sys.getTMMCBias()->print("final", false);

    std::cout << "Finished simulation at " << getTimeStamp() << std::endl;
	return SAFE_EXIT;
}
