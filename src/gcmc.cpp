#include "gcmc.h"

/*!
 * Perform an unbiased GCMC simulation
 *
 * \param [in] sys System to simulate
 * \param [in] usedMovesEq Move class to use for equilibration phase
 * \param [in] usedMovesPr Move class to use for production phase
 */
void performGCMC (simSystem &sys, moves *usedMovesEq, moves *usedMovesPr) {
    // Equilibration Phase
    std::ofstream ofs;
    ofs.open ("thermoEq.log", std::ofstream::out);
    if (!ofs.is_open()) throw customException ("Unable to open log file to record equilibration GCMC thermodynamics");
    ofs << "# Step\tU\tN_tot\t";
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
        ofs << "N_"+numToStr(i+1)+"\t";
	}
    ofs << std::endl;

    double ctr = 0, thermo_ctr = 0; // Count number of times system was in M = 0 state
    for (double move = 0; move < sys.gcmcEqSteps; move += 1.0) {
		try {
			usedMovesEq->makeMove(sys);
		} catch (customException &ce) {
            sendErr(ce.what());
            exit(SYS_FAILURE);
		}

        // Only record properties of the system when it is NOT in an intermediate state
        if (thermo_ctr < sys.gcmcThermoFreq && sys.getCurrentM() == 0) {
            thermo_ctr += 1.0;
        } else if (thermo_ctr >= sys.gcmcThermoFreq && sys.getCurrentM() == 0) {
            double sys_ntot = 0;
            for (unsigned int i = 0; i < sys.nSpecies(); ++i) { sys_ntot += sys.numSpecies[i]; }
            ofs << move << "\t" << sys.energy() << "\t" << sys_ntot << "\t";
            for (unsigned int i = 0; i < sys.nSpecies(); ++i) { ofs << sys.numSpecies[i]/sys_ntot << "\t"; }
            ofs << std::endl;
            thermo_ctr = 0.0;
        }

        if (ctr < sys.gcmcSnapFreq) {
            ctr += 1.0;
        } else {
            ctr = 0.0;
            usedMovesEq->print("equilibration.stats");
        }
	}
    ofs.close();
    usedMovesEq->print("equilibration.stats");
    sanityChecks(sys);

    // Production phase
    ofs.open ("thermoPr.log", std::ofstream::out);
    if (!ofs.is_open()) throw customException ("Unable to open log file to record production GCMC thermodynamics");
    ofs << "# Step\tU\tN_tot\t";
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) { ofs << "N_"+numToStr(i+1)+"\t"; }
    ofs << std::endl;

    ctr = 0;
    thermo_ctr = 0;
    for (double move = 0; move < sys.gcmcPrSteps; move += 1.0) {
		try {
			usedMovesPr->makeMove(sys);
		} catch (customException &ce) {
            sendErr(ce.what());
            exit(SYS_FAILURE);
		}

        // Only record properties of the system when it is NOT in an intermediate state
        if (thermo_ctr < sys.gcmcThermoFreq && sys.getCurrentM() == 0) {
            thermo_ctr += 1.0;
        } else if (thermo_ctr >= sys.gcmcThermoFreq && sys.getCurrentM() == 0) {
            double sys_ntot = 0;
            for (unsigned int i = 0; i < sys.nSpecies(); ++i) { sys_ntot += sys.numSpecies[i]; }
            ofs << move << "\t" << sys.energy() << "\t" << sys_ntot << "\t";
            for (unsigned int i = 0; i < sys.nSpecies(); ++i) { ofs << sys.numSpecies[i]/sys_ntot << "\t"; }
            ofs << std::endl;
            thermo_ctr = 0.0;
        }

        // Can print snapshot regardless of M state, since partial atoms are neglected from printing routine
        if (ctr < sys.gcmcSnapFreq) {
            ctr += 1.0;
        } else {
            ctr = 0.0;
            usedMovesPr->print("production.stats");
            sys.printSnapshot("movie.xyz", numToStr(move), false);
            sendMsg ("Completed "+numToStr(move)+"/"+numToStr(sys.gcmcPrSteps)+" production steps");
        }
	}
    ofs.close();
    usedMovesPr->print("production.stats");
    sys.printSnapshot("movie.xyz", numToStr(sys.gcmcPrSteps), false);
    sanityChecks(sys);
}
