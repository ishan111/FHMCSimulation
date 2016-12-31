#include "../src/fhmc.h"

/*!
 * Usage: ./binary_name input.json
 */
int main (int argc, char * const argv[]) {
	sendMsg ("Beginning simulation");

	moves usedMovesEq, usedMovesPr;
	simSystem sys = initialize (argv[1], &usedMovesEq, &usedMovesPr);
	setConfig (sys, argv[1]);

	try {
		performGCMC (sys, &usedMovesEq, &usedMovesPr);
	} catch (std::exception &ex) {
		sendErr (ex.what());
		exit(SYS_FAILURE);
	}
	
	sendMsg ("Finished simulation");
	return SAFE_EXIT;
}
