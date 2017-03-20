#include "../src/fhmc.h"

/*!
 * Usage: ./binary_name input.json
 */
int main (int argc, char * const argv[]) {
	sendMsg ("Beginning simulation");

	moves usedMovesEq, usedMovesPr;
	simSystem sys = initialize (argv[1], &usedMovesEq, &usedMovesPr);

	// If checkpoint exists, default to use this information to restart the simulation
	checkpoint cpt ("checkpt", 600, sys);
    if (!cpt.hasCheckpoint) {
        setConfig (sys, argv[1]);
    }

	// Check not falsely restarting the simulation
	if (cpt.tmmcDone) {
		sendMsg ("TMMC stage already finished");
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
		sendErr ("Error in establishing which stage to begin from");
		return SYS_FAILURE;
	}

	// Dump a final checkpoint to indicate the simulation has finished, so do not restart it
	cpt.dump(sys);

	sendMsg ("Finished simulation");
	return SAFE_EXIT;
}
