/*!
 * A general purpose grand canonical Monte Carlo code
 * 
 * \author Nathan A. Mahynski
 * \date 08/01/14
 */

#include <iostream>
#include "system.h"
#include "utilities.h"
#include "global.h"
#include <boost/lexical_cast.hpp>
#include "insert.h"
#include "delete.h"
#include "translate.h"
#include "moves.h"
#include "input.h"
#include <time.h>
#include <fstream>
#include <cmath>

/*! Only uncomment this if simulations are purely in the fluid phase.  
 *This will enable tail corrections which are only valid assuming a converging g(r) at large r.
 */
//#define FLUID_PHASE_SIMULATIONS

int main (int argc, char * const argv[]) {
	// get time stamp
	time_t rawtime;
	time (&rawtime);
	struct tm * timeinfo;
	timeinfo = localtime (&rawtime);
	
    // user input parsing can be added later
	RNG_SEED = -10;
	const int sysN = 1;
	const double sysBeta = 1.0;
	const std::vector < double > sysBox (3, 6);
	std::vector < double > sysMu (sysN, 2.5);
	std::vector < int > sysMax (sysN, 300);
	simSystem sys (sysN, sysBeta, sysBox, sysMu, sysMax);
	
	const int tmmcSweepSize = 1e6, totalTMMCSweep = 1000, wlSweepSize = 1e6;
	
	// specify pair potentials and add to system
    const double eps = 1.0, lambda = 0.3, sigma = 1.0;
    std::vector < double > params (3);
    params[2] = -eps;
    params[1] = sigma+lambda;
    params[0] = sigma;
    squareWell sqW;
	sqW.setParameters(params);
    sqW.savePotential("sqW_potential.dat", 0.01, 0.01);
    sys.addPotential (0, 0, &sqW, true);

    // specify moves to use for the system
    moves usedMovesEq, usedMovesPr;
	std::vector < double > moveProb (sysN, 0.2);	// add these to input parser in the future
	std::vector < insertParticle > insertions (sysN);
	std::vector < deleteParticle > deletions (sysN);
	std::vector < translateParticle > translations (sysN);
	for (unsigned int i = 0; i < sysN; ++i) {
		insertParticle newIns (i, "insert");
		deleteParticle newDel (i, "delete");	
		translateParticle newTranslate(i, "translate");
		insertions[i] = newIns;
		deletions[i] = newDel;
		translations[i] = newTranslate;
		usedMovesEq.addMove (&insertions[i], moveProb[i]);
		usedMovesEq.addMove (&deletions[i], moveProb[i]); 
		usedMovesEq.addMove (&translations[i], 3*moveProb[i]); // probs are symmetric in each direction
		usedMovesPr.addMove (&insertions[i], moveProb[i]);
		usedMovesPr.addMove (&deletions[i], moveProb[i]); 
		usedMovesPr.addMove (&translations[i], 3*moveProb[i]); // probs are symmetric in each direction
	}

	// check all pair potentials have been set and all r_cut < L/2
	double minL = sys.box()[0];
	for (unsigned int i = 1; i < 2; ++i) {
		minL = std::min(minL, sys.box()[i]);
	}
	for (unsigned int i = 0; i < sysN; ++i) {
		for (unsigned int j = 0; j < sysN; ++j) {
			if (!sys.potentialIsSet(i, j)) {
				std::cerr << "Not all pair potentials are set" << std::endl;
				exit(SYS_FAILURE);
			}
			if (!(sys.ppot[i][j]->rcut() < minL/2.0)) {
				std::cerr << "Pair potential r_cut for species " << i << ", " << j << " is > L/2" << std::endl;
				exit(SYS_FAILURE);
			}
		}
	}
	
	// Initially do a WL simulation
	const double g = 0.5, s = 0.8;
	const std::vector <int> maxN (1, sys.maxSpecies(0)), minN (1, 0);
	double lnF = 1;
	bool flat = false;
	sys.startWALA (lnF, g, s, sysN, maxN, minN);
	
	while (lnF > 2.0e-18) {
		for (unsigned int move = 0; move < wlSweepSize; ++move) {
			try {
				usedMovesEq.makeMove(sys);
			} catch (customException &ce) {
				std::cerr << ce.what() << std::endl;
				exit(SYS_FAILURE);
			}
		}
			
		// Check if bias has flattened out
		flat = sys.getWALABias()->evaluateFlatness();
		if (flat) {
			// if flat, need to reset H and reduce lnF
			sys.getWALABias()->iterateForward();
		}
		lnF = sys.getWALABias()->lnF();
		
		// Periodically write out checkpoints
	}
	
	// After a while, combine to initialize TMMC collection matrix
// have to add this to use still
	sys.startTMMC (sysN, maxN, minN);
	int count = 0;
	while (count < 2) {
		for (unsigned int move = 0; move < wlSweepSize; ++move) {
			try {
				usedMovesEq.makeMove(sys);
			} catch (customException &ce) {
				std::cerr << ce.what() << std::endl;
				exit(SYS_FAILURE);
			}
		}
				
		// Check if bias has flattened out
		flat = sys.getWALABias()->evaluateFlatness();
		if (flat) {
			// If flat, need to reset H and reduce lnF
			sys.getWALABias()->iterateForward();
			count++;
		}
		
		// Periodically write out checkpoints
	}
	
	// Print checkpoint here
	
	// Switch over to TMMC completely
	sys.stopWALA();
	for (unsigned int sweep = 0; sweep < totalTMMCSweep; ++sweep) {
		for (unsigned int move = 0; move < tmmcSweepSize; ++move) {
			try {
				usedMovesPr.makeMove(sys);
			} catch (customException &ce) {
				std::cerr << ce.what() << std::endl;
				exit(SYS_FAILURE);
			}
		}
					
		// Update biasing function
		
		// Periodically write out checkpoints
	}
		
	// Sanity checks
	if (sys.nSpecies() != sys.atoms.size()) {
        std::cerr << "Error: Number of components changed throughout simulation" << std::endl;
        exit(SYS_FAILURE);
    }
    const double tol = 1.0e-9;
    const double scratchEnergy = sys.scratchEnergy(), incrEnergy = sys.energy();
    if (fabs(scratchEnergy - incrEnergy) > tol) {
        std::cerr << "Error: scratch energy calculation = " << scratchEnergy << ", but incremental = " << incrEnergy << ", |diff| = " << fabs(scratchEnergy - incrEnergy) << std::endl;
        exit(SYS_FAILURE);
    }
    
    // Report move statistics for final TMMC ("production") stage
	char statName [80];
	strftime (statName,80,"%Y_%m_%d_%H_%M_%S-stats.log",timeinfo);
	std::ofstream statFile (statName);
    std::vector < double > stats = usedMovesPr.reportMoveStatistics();
    statFile << " ----- Move Statistics ----- " << std::endl << " Move\t\% Success" << std::endl;
    for (unsigned int i = 0; i < stats.size(); ++i) {
        statFile << usedMovesPr.includedMoves()[i]->myName() << "\t" << stats[i]*100.0 << std::endl;
    }
    statFile << std::endl;
    statFile.close();
	
    // print out restart file (xyz)
    sys.printSnapshot("final.xyz", "last configuration");
    
	return SAFE_EXIT;
}
