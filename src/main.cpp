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

//! Only uncomment this if simulations are purely in the fluid phase.  This will enable tail corrections which are only valid assuming a converging g(r) at large r.
//#define FLUID_PHASE_SIMULATIONS

int main (int argc, char * const argv[]) {
	// get time stamp
	time_t rawtime;
	time (&rawtime);
	struct tm * timeinfo;
	timeinfo = localtime (&rawtime);
	
    // user input parsing
    std::vector < std::string > args;
    try {
        args = parseInput (argc, argv, rawtime);
    } catch (customException &ce) {
		std::cerr << ce.what() << std::endl;
        exit(SYS_FAILURE);
    }
	
    // define system & simulation parameters
    const int sysN = boost::lexical_cast<int>(args[0]);
    const double sysBeta = boost::lexical_cast<double>(args[2]);
    const std::vector < double > sysBox (3, boost::lexical_cast<double>(args[1]));
    std::vector < double > sysMu (sysN);
    for (unsigned int i = 0; i < sysMu.size(); ++i) {
        sysMu[i] = boost::lexical_cast<double>(args[3+i]);
    }
    std::vector < int > sysMax (sysN);
    for (unsigned int i = 0; i < sysMax.size(); ++i) {
        sysMax[i] = boost::lexical_cast<int>(args[3+sysMu.size()+i]);
    }
    simSystem sys (sysN, sysBeta, sysBox, sysMu, sysMax);
	const unsigned int sweepSize = boost::lexical_cast<int>(args[3+sysMu.size()+sysMu.size()]);
    const unsigned int prodSweep = boost::lexical_cast<int>(args[3+sysMu.size()+sysMu.size()+1]);
	const unsigned int equilSweep = boost::lexical_cast<int>(args[3+sysMu.size()+sysMu.size()+2]);
	RNG_SEED = boost::lexical_cast<int>(args[3+sysMu.size()+sysMu.size()+3]);
	std::string restartFileName = args[3+sysMu.size()+sysMu.size()+4];
	if (restartFileName != "none") {
		sys.readRestart(restartFileName);
	}
	
	/* -------------------------------------- */
	// specify pair potentials and add to system - for now this is manual
	// Lennard Jones Potential
	/*const double eps = 1.0, sig = 1.0, rcut = 3.5, ushift = 0;
	std::vector < double > paramsAA (4);
	paramsAA[0] = eps;
	paramsAA[1] = sig;
	paramsAA[2] = rcut;
	paramsAA[3] = ushift;
    lennardJones ppAA;
	ppAA.setParameters(paramsAA);
	ppAA.savePotential("potentialAA.dat", 0.01, 0.01);

	std::vector < double > paramsBB (4);
	paramsBB[0] = 0.5*eps;
	paramsBB[1] = sig;
	paramsBB[2] = rcut;
	paramsBB[3] = ushift;
    lennardJones ppBB;
	ppBB.setParameters(paramsBB);
	ppBB.savePotential("potentialBB.dat", 0.01, 0.01);

	std::vector < double > paramsAB (4);
	paramsAB[0] = sqrt(eps*0.5*eps);
	paramsAB[1] = sig;
	paramsAB[2] = rcut;
	paramsAB[3] = ushift;
    lennardJones ppAB;
	ppAB.setParameters(paramsAB);
	ppAB.savePotential("potentialAB.dat", 0.01, 0.01);*/
	
	// star-star, colloid-colloid, and star-colloid interactions
	tabulated ppAA;
	ppAA.loadPotential("potStarStar.dat");
	ppAA.savePotential("potentialAA.dat", 0.00, 0.01);

	std::vector < double > paramsBB (1);
	paramsBB[0] = 8.0;
    hardCore ppBB;
    ppBB.setParameters(paramsBB);
	ppBB.useTailCorrection = false;
	ppBB.savePotential("potentialBB.dat", 0.00, 0.01);
	
	tabulated ppAB;
	ppAB.loadPotential("potStarColloid.dat");
	ppAB.savePotential("potentialAB.dat", 0.00, 0.01);
	
	for (unsigned int i = 0; i < sysN; ++i) {
		for (unsigned int j = 0; j < sysN; ++j) {
			if (i == j) {
				if (i == 0) {
					sys.addPotential (i, j, &ppAA, true);
				}
				else if (i == 1) {
					sys.addPotential (i, j, &ppBB, true);
				}
			}
			else {
				sys.addPotential (i, j, &ppAB, true);
			}
		}
	}
	/* -------------------------------------- */
	

    // specify moves to use for a binary system
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
		usedMovesEq.addMove (&translations[i], moveProb[i]); // probs are symmetric in each direction
		usedMovesPr.addMove (&insertions[i], moveProb[i]);
		usedMovesPr.addMove (&deletions[i], moveProb[i]); 
		usedMovesPr.addMove (&translations[i], moveProb[i]); // probs are symmetric in each direction
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
	
    // equilibration
	for (unsigned int sweep = 0; sweep < equilSweep; ++sweep) {
        for (unsigned int move = 0; move < sweepSize; ++move) {
            try {
                usedMovesEq.makeMove(sys);              
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
                exit(SYS_FAILURE);
            }
        }
    }
	
	// production
	char prodLog [80];
	strftime (prodLog,80,"%Y_%m_%d_%H_%M_%S-thermo.log",timeinfo);
	std::ofstream prodFile (prodLog);
	for (unsigned int sweep = 0; sweep < prodSweep; ++sweep) {
        for (unsigned int move = 0; move < sweepSize; ++move) {
            try {
                usedMovesPr.makeMove(sys);
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
                exit(SYS_FAILURE);
            }
        }
        for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
            prodFile << sys.numSpecies[i] << "\t";
        }
        prodFile << sys.energy() << std::endl;
    }
	prodFile.close();
	
    // check total energy and number of particles in the system at the end
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
    
    // report move statistics
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
