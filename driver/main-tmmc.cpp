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
	simSystem sys = initialize (argv[1], &usedMovesEq, &usedMovesPr);

	// If checkpoint exists, default to use this information to restart the simulation
	checkpoint cpt ("checkpt", 1, sys);
    if (!cpt.hasCheckpoint) {
        setup (sys, argv[1]);
    }

	// Check not falsely restarting the simulation
	if (cpt.tmmcDone) {
		std::cerr << "TMMC stage already finished, terminating" << std::endl;
		return SAFE_EXIT;
	}

	// Choose stage based on what is completed, not where restart is from in case not restarting from checkpoint
	if (!cpt.walaDone) {
		// Perform Wang-Landau simulation
		performWALA (sys, cpt, &usedMovesEq);
		performCrossover (sys, cpt, &usedMovesEq);
		performTMMC (sys, cpt, &usedMovesPr);
	} else if (!cpt.crossoverDone) {
		// Crossover to TMMC simulation
		performCrossover (sys, cpt, &usedMovesEq);
		performTMMC (sys, cpt, &usedMovesPr);
	} else if (!cpt.tmmcDone) {
		// Perform TMMC portion of the simulation
		performTMMC (sys, cpt, &usedMovesPr);
	} else {
		std::cerr << "Error in establishing which stage to begin from" << std::endl;
		return SYS_FAILURE;
	}

	// Dump a final checkpoint to indicate the simulation has finished, so do not restart it
	cpt.dump(sys);

    std::cout << "Finished simulation at " << getTimeStamp() << std::endl;
	return SAFE_EXIT;
}
