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
	std::cout << "Beginning simulation at " << getTimeStamp() << std::endl;

	moves usedMovesEq, usedMovesPr;
	simSystem sys = initialize (argv[1], &usedMovesEq, &usedMovesPr); // read json file and create system class

	// if restart/ exists, default to use this information to restart the simulation
	checkpoint cpt ("checkpt", 900);
    if (!cpt.hasCheckpoint) {
        setup (sys, argv[1]);
    }

	// 1.
	// get checlpoint to write vals it has to json
	// 2.
	// from input parsing, if restart from TMMC or another one, set previous stages to completion so it goes right to it?
	// store this into cpt object
	// move restartFromWALA reading, etc. in input sfrom inistialize to setup and then store in cpt
	// remove double reading from setup (like box) and just set values read from initialize
	// 3.
	// make reader
	// 4.
	// restore extMom and eHist and pkHist also from snapshot?
	// perhaps just add read() member in addition to print() ones for histogram classes - no need to add to checkpoint class
	// just have checkpoint class use native print() function then

	if (!cpt.walaDone) {
		// perform Wang-Landau simulation
		performWALA (sys, cpt, &usedMovesEq);
		performCrossover (sys, cpt, &usedMovesEq);
		performTMMC (sys, cpt, &usedMovesPr);
	} else if (!cpt.crossoverDone) {
		// crossover to TMMC
		performCrossover (sys, cpt, &usedMovesEq);
		performTMMC (sys, cpt, &usedMovesPr);
	} else if (!cpt.tmmcDone) {
		// perform tmmc portion of the simulation
		performTMMC (sys, cpt, &usedMovesPr);
	}

	// Perform sanity checks
	sanityChecks(sys);

	// Print final results
    sys.printSnapshot("final.xyz", "last configuration");
    sys.refineEnergyHistogramBounds();
    sys.printEnergyHistogram("final_eHist");
    sys.refinePkHistogramBounds();
    sys.printPkHistogram("final_pkHist");
    sys.printExtMoments("final_extMom");
    sys.getTMMCBias()->print("final", false);

	// Dump a final checkpoint
	cpt.dump(sys);

    std::cout << "Finished simulation at " << getTimeStamp() << std::endl;
	return SAFE_EXIT;
}
