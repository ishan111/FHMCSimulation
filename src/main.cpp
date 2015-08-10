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
 * \section Input
 * 
 */

#include <iostream>
#include <time.h>
#include <fstream>
#include <cmath>
#include <boost/lexical_cast.hpp>
#include "system.h"
#include "utilities.h"
#include "global.h"
#include "insert.h"
#include "delete.h"
#include "translate.h"
#include "swap.h"
#include "moves.h"

// JSON interface from local distro of rapidjson
#include "rapidjson/include/rapidjson/document.h"
#include "rapidjson/include/rapidjson/writer.h"
#include "rapidjson/include/rapidjson/stringbuffer.h"
#include "rapidjson/include/rapidjson/filereadstream.h"
#include "rapidjson/include/rapidjson/prettywriter.h"

/*! 
 * Only uncomment this if simulations are purely in the fluid phase.  
 * This will allow tail corrections to be enabled which are only valid assuming a converging g(r) at large r.
 */
//#define FLUID_PHASE_SIMULATIONS

/*!
 * Uncomment this if netCDF libraries are installed and can be compiled against.  Data will be output to these arrays
 * instead of ASCII files if so.
 */
//#define NETCDF_CAPABLE

/*!
 * Usage: ./binary_name inputFile.json
 */
int main (int argc, char * const argv[]) {
	// Get time stamp
	time_t rawtime;
	time (&rawtime);
	struct tm * timeinfo;
	timeinfo = localtime (&rawtime);
	char timestamp [80];
	strftime (timestamp,80,"%d/%m/%Y %H:%M:%S",timeinfo);
	std::cout << "Beginning simulation at " << timestamp << std::endl;
	
	/* -------------------- BEGIN INPUT -------------------- */
	
	// Parse input JSON file
	FILE* fp = fopen(argv[1], "r");
	char readBuffer[65536];
	rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));
	rapidjson::Document doc;
	doc.ParseStream(is);	
	fclose(fp);
		
	// Assert that this is a JSON document
	assert(doc.IsObject());
	
	// Check each member exists and is in the correct format
	assert(doc.HasMember("num_species"));
	assert(doc["num_species"].IsInt());
	assert(doc.HasMember("beta"));
	assert(doc["beta"].IsNumber());
	
	assert(doc.HasMember("box"));
	assert(doc["box"].IsArray());
	assert(doc["box"].Size() == 3);
	std::vector < double > sysBox (3, 0);
	for (rapidjson::SizeType i = 0; i < doc["box"].Size(); ++i) {
		assert(doc["box"][i].IsNumber());
		sysBox[i] = doc["box"][i].GetDouble();
	}

	assert(doc.HasMember("mu"));
	assert(doc["mu"].IsArray());
	assert(doc["mu"].Size() == doc["num_species"].GetInt());
	std::vector < double > sysMu (doc["mu"].Size(), 0);
	for (rapidjson::SizeType i = 0; i < doc["mu"].Size(); ++i) {
		assert(doc["mu"][i].IsNumber());
		sysMu[i] = doc["mu"][i].GetDouble();
	}

	assert(doc.HasMember("seed"));
	assert(doc["seed"].IsInt());
	RNG_SEED = doc["seed"].GetInt();
	
	assert(doc.HasMember("max_N"));
	assert(doc["max_N"].IsArray());
	assert(doc["max_N"].Size() == doc["num_species"].GetInt());
	std::vector < int > sysMax (doc["max_N"].Size(), 0);
	for (rapidjson::SizeType i = 0; i < doc["max_N"].Size(); ++i) {
		assert(doc["max_N"][i].IsInt());
		sysMax[i] = doc["max_N"][i].GetInt();
	}

	assert(doc.HasMember("min_N"));
	assert(doc["min_N"].IsArray());
	assert(doc["min_N"].Size() == doc["num_species"].GetInt());
	std::vector < int > sysMin (doc["min_N"].Size(), 0);
	for (rapidjson::SizeType i = 0; i < doc["min_N"].Size(); ++i) {
		assert(doc["min_N"][i].IsInt());
		sysMin[i] = doc["min_N"][i].GetInt();
	}	

	simSystem sys (doc["num_species"].GetInt(), doc["beta"].GetDouble(), sysBox, sysMu, sysMax, sysMin);
	
	std::vector < int > sysWindow;
	if (doc.HasMember("window")) {
		assert(doc["window"].IsArray());
		assert(doc["window"].Size() == 2);
		sysWindow.resize(2, 0);
		sysWindow[0] = doc["window"][0].GetInt();
		sysWindow[1] = doc["window"][1].GetInt();
	}

	if (sysWindow.begin() != sysWindow.end()) {
		sys.setTotNBounds(sysWindow);
	}
	
	assert(doc.HasMember("restart_file"));
	assert(doc["restart_file"].IsString());
	const std::string restart_file = doc["restart_file"].GetString();

	assert(doc.HasMember("tmmc_sweep_size"));
	assert(doc["tmmc_sweep_size"].IsNumber());
	double tmpT = doc["tmmc_sweep_size"].GetDouble(); // possibly in scientific notation
	const long long int tmmcSweepSize = tmpT; // convert

	assert(doc.HasMember("total_tmmc_sweeps"));
	assert(doc["total_tmmc_sweeps"].IsNumber());
	double tmpS = doc["total_tmmc_sweeps"].GetDouble(); // possibly in scientific notation
	const long long int totalTMMCSweeps = tmpS; // convert

	assert(doc.HasMember("wala_sweep_size"));
	assert(doc["wala_sweep_size"].IsNumber());
	double tmpW = doc["wala_sweep_size"].GetDouble(); // possibly in scientific notation
	const long long int wlSweepSize = tmpW; // convert

	assert(doc.HasMember("wala_g"));
	assert(doc["wala_g"].IsNumber());
	const double g = doc["wala_g"].GetDouble();

	assert(doc.HasMember("wala_s"));
	assert(doc["wala_s"].IsNumber());
	const double s = doc["wala_s"].GetDouble();
	
	std::vector < double > ref (sys.nSpecies(), 0);
	std::vector < std::vector < double > > probEqSwap (sys.nSpecies(), ref), probPrSwap (sys.nSpecies(), ref);
	std::vector < double > probPrInsDel (sys.nSpecies(), 0), probPrDisp (sys.nSpecies(), 0);
	std::vector < double > probEqInsDel (sys.nSpecies(), 0), probEqDisp (sys.nSpecies(), 0);
	std::vector < double > maxPrD (sys.nSpecies(), 0), maxEqD (sys.nSpecies(), 0);
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		std::string dummy = "prob_pr_ins_del_" + sstr(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		probPrInsDel[i] = doc[dummy.c_str()].GetDouble();
	}
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		std::string dummy = "prob_pr_displace_" + sstr(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		probPrDisp[i] = doc[dummy.c_str()].GetDouble();
		dummy = "max_pr_displacement_" + sstr(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		maxPrD[i] = doc[dummy.c_str()].GetDouble();
	}
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		std::string dummy = "prob_eq_ins_del_" + sstr(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		probEqInsDel[i] = doc[dummy.c_str()].GetDouble();
	}
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		std::string dummy = "prob_eq_displace_" + sstr(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		probEqDisp[i] = doc[dummy.c_str()].GetDouble();
		dummy = "max_eq_displacement_" + sstr(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		maxEqD[i] = doc[dummy.c_str()].GetDouble();
	}
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		for (unsigned int j = i+1; j < sys.nSpecies(); ++j) {
			std::string name1 = "prob_pr_swap_"+sstr(i+1)+"_"+sstr(j+1);
			std::string name2 = "prob_pr_swap_"+sstr(j+1)+"_"+sstr(i+1);
			std::string moveName = "";
			bool foundIJ = false;
			if (doc.HasMember(name1.c_str())) {
				moveName = name1;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && !foundIJ) {
				moveName = name2;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && foundIJ) {
				std::cerr << "Input file doubly specifies production swap move probability for species pair ("+sstr(i+1)+", "+sstr(j+1)+")" << std::endl;
				exit(SYS_FAILURE);
			} else {
				std::cerr << "Input file does not specify production swap move probability for species pair ("+sstr(i+1)+", "+sstr(j+1)+")" << std::endl;
				exit(SYS_FAILURE);
			}
			assert(doc[moveName.c_str()].IsNumber());
			probPrSwap[i][j] = doc[moveName.c_str()].GetDouble();
			probPrSwap[j][i] = doc[moveName.c_str()].GetDouble();
		}
	}
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		for (unsigned int j = i+1; j < sys.nSpecies(); ++j) {
			std::string name1 = "prob_eq_swap_"+sstr(i+1)+"_"+sstr(j+1);
			std::string name2 = "prob_eq_swap_"+sstr(j+1)+"_"+sstr(i+1);
			std::string moveName = "";
			bool foundIJ = false;
			if (doc.HasMember(name1.c_str())) {
				moveName = name1;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && !foundIJ) {
				moveName = name2;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && foundIJ) {
				std::cerr << "Input file doubly specifies equilibration swap move probability for species pair ("+sstr(i+1)+", "+sstr(j+1)+")" << std::endl;
				exit(SYS_FAILURE);
			} else {
				std::cerr << "Input file does not specify equilibration swap move probability for species pair ("+sstr(i+1)+", "+sstr(j+1)+")" << std::endl;
				exit(SYS_FAILURE);
			}
			assert(doc[moveName.c_str()].IsNumber());
			probEqSwap[i][j] = doc[moveName.c_str()].GetDouble();
			probEqSwap[j][i] = doc[moveName.c_str()].GetDouble();
		}
	}	
	std::vector < pairPotential* > ppotArray (sys.nSpecies()*(sys.nSpecies()-1)/2);
	std::vector < std::string > ppotType (sys.nSpecies()*(sys.nSpecies()-1)/2 + sys.nSpecies());
	int ppotIndex = 0, ppotTypeIndex = 0;
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		for (unsigned int j = i; j < sys.nSpecies(); ++j) {
			std::string name1 = "ppot_"+sstr(i+1)+"_"+sstr(j+1);
			std::string name2 = "ppot_"+sstr(j+1)+"_"+sstr(i+1);
			std::string ppotName = "", dummy = "";
			bool foundIJ = false;
			if (doc.HasMember(name1.c_str())) {
				ppotName = name1;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && !foundIJ) {
				ppotName = name2;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && foundIJ) {
				std::cerr << "Input file doubly specifies pair potential for species pair ("+sstr(i+1)+", "+sstr(j+1)+")" << std::endl;
				exit(SYS_FAILURE);
			} else {
				std::cerr << "Input file does not specify pair potential for species pair ("+sstr(i+1)+", "+sstr(j+1)+")" << std::endl;
				exit(SYS_FAILURE);
			} 
			assert(doc[ppotName.c_str()].IsString());
			ppotType[ppotTypeIndex] = doc[ppotName.c_str()].GetString();
			dummy = ppotName+"_params";
			assert(doc.HasMember(dummy.c_str()));
			assert(doc[dummy.c_str()].IsArray());
			std::vector < double > params (doc[dummy.c_str()].Size(), 0);
			for (unsigned int k = 0; k < params.size(); ++k) {
				assert(doc[dummy.c_str()][k].IsNumber());
				params[k] = doc[dummy.c_str()][k].GetDouble();
			}
			bool useCellList = false; // default
			dummy = ppotName+"_use_cell_list";
			if (doc.HasMember(dummy.c_str())) {
				assert(doc[dummy.c_str()].IsBool());
				useCellList = doc[dummy.c_str()].GetBool();
			}
			if (ppotType[ppotTypeIndex] == "square_well") {
				try {
					ppotArray[ppotIndex] = new squareWell;
					ppotArray[ppotIndex]->setParameters(params);
				} catch (customException &ce) {
					std::cerr << ce.what() << std::endl;
					exit(SYS_FAILURE);
				}
				ppotArray[ppotIndex]->savePotential(ppotName+".dat", 0.01, 0.01);
				sys.addPotential (i, j, ppotArray[ppotIndex], useCellList);
			} else if (ppotType[ppotTypeIndex] == "lennard_jones") {
				try {
					ppotArray[ppotIndex] = new lennardJones;
					ppotArray[ppotIndex]->setParameters(params);
				} catch (customException &ce) {
					std::cerr << ce.what() << std::endl;
					exit(SYS_FAILURE);
				}
				ppotArray[ppotIndex]->savePotential(ppotName+".dat", 0.01, 0.01);
				sys.addPotential (i, j, ppotArray[ppotIndex], useCellList);				
			} else if (ppotType[ppotTypeIndex] == "hard_sphere") {
				try {
					ppotArray[ppotIndex] = new hardCore;
					ppotArray[ppotIndex]->setParameters(params);
				} catch (customException &ce) {
					std::cerr << ce.what() << std::endl;
					exit(SYS_FAILURE);
				}
				ppotArray[ppotIndex]->savePotential(ppotName+".dat", 0.01, 0.01);
				sys.addPotential (i, j, ppotArray[ppotIndex], useCellList);
			} else if (ppotType[ppotTypeIndex] == "tabulated") {
				try {
					ppotArray[ppotIndex] = new tabulated;
					ppotArray[ppotIndex]->setParameters(params);
				} catch (customException &ce) {
					std::cerr << ce.what() << std::endl;
					exit(SYS_FAILURE);
				}
				ppotArray[ppotIndex]->savePotential(ppotName+".dat", 0.01, 0.01);
				sys.addPotential (i, j, ppotArray[ppotIndex], useCellList);
			} else {
				std::cerr << "Unrecognized pair potential name for species "<< ppotTypeIndex << std::endl;
				exit(SYS_FAILURE);
			}
			ppotTypeIndex++;
			ppotIndex++;
		}
	}

	// check all pair potentials have been set and all r_cut < L/2
	double minL = sys.box()[0];
	for (unsigned int i = 1; i < 2; ++i) {
		minL = std::min(minL, sys.box()[i]);
	}
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		for (unsigned int j = 0; j < sys.nSpecies(); ++j) {
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
	
	// specify moves to use for the system
    moves usedMovesEq, usedMovesPr;
	std::vector < insertParticle > eqInsertions (sys.nSpecies()), prInsertions (sys.nSpecies());
	std::vector < deleteParticle > eqDeletions (sys.nSpecies()), prDeletions (sys.nSpecies());
	std::vector < translateParticle > eqTranslations (sys.nSpecies()), prTranslations (sys.nSpecies());
	std::vector < swapParticles > eqSwaps (sys.nSpecies()*(sys.nSpecies()-1)/2), prSwaps (sys.nSpecies()*(sys.nSpecies()-1)/2);
	
	int swapCounter = 0;
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		insertParticle newIns (i, "insert");
		eqInsertions[i] = newIns;
		usedMovesEq.addMove (&eqInsertions[i], probEqInsDel[i]);
		
		deleteParticle newDel (i, "delete");
		eqDeletions[i] = newDel;
		usedMovesEq.addMove (&eqDeletions[i], probEqInsDel[i]);
		
		translateParticle newTranslate (i, "translate");
		eqTranslations[i] = newTranslate;
		usedMovesEq.addMove (&eqTranslations[i], probEqDisp[i]);
		
		insertParticle newIns2 (i, "insert");
		prInsertions[i] = newIns2;
		usedMovesPr.addMove (&prInsertions[i], probPrInsDel[i]);
		
		deleteParticle newDel2 (i, "delete");
		prDeletions[i] = newDel2;
		usedMovesPr.addMove (&prDeletions[i], probPrInsDel[i]);		
		
		translateParticle newTranslate2 (i, "translate");
		prTranslations[i] = newTranslate2;
		usedMovesPr.addMove (&prTranslations[i], probPrDisp[i]);
		
		for (unsigned int j = i+1; j < sys.nSpecies(); ++j) {
			swapParticles newSwap (i, j, "swap");
			eqSwaps[swapCounter] = newSwap;
			usedMovesEq.addMove (&eqSwaps[swapCounter], probEqSwap[i][j]);
			
			swapParticles newSwap2 (i, j, "swap");
			prSwaps[swapCounter] = newSwap2;
			usedMovesPr.addMove (&prSwaps[swapCounter], probPrSwap[i][j]);
			
			swapCounter++;
		}
	}
	
	/* -------------------- END INPUT -------------------- */
	
	std::cout << "Beginning Wang-Landau portion" << std::endl;
	
	// Initially do a WL simulation
	double lnF = 1;
	bool flat = false;
	sys.startWALA (lnF, g, s, sys.totNMax(), sys.totNMin()); //!< Using Shen and Errington method this syntax is same for single and multicomponent
	while (lnF > 2.0e-18) {
		for (unsigned int move = 0; move < wlSweepSize; ++move) {
			try {
				usedMovesEq.makeMove(sys);
			} catch (customException &ce) {
				std::cerr << ce.what() << std::endl;
				exit(SYS_FAILURE);
			}
			
			// record U
			sys.recordU();
		}
			
		// Check if bias has flattened out
		flat = sys.getWALABias()->evaluateFlatness();
		if (flat) {
			// if flat, need to reset H and reduce lnF
			sys.getWALABias()->iterateForward();
			lnF = sys.getWALABias()->lnF();
			flat = false;
			
			time_t rawtime_tmp;
			time (&rawtime_tmp);
			struct tm * timeinfo_tmp;
			timeinfo_tmp = localtime (&rawtime_tmp);
			char dummy_tmp [80];
			strftime (dummy_tmp,80,"%d/%m/%Y %H:%M:%S",timeinfo);
			std::cout << "lnF = " << lnF << " at " << dummy_tmp << std::endl;
			
			// Periodically write out checkpoints
			sys.getWALABias()->print("wl-Checkpoint", true);
		}
	}
	
	std::cout << "Crossing over to build TMMC matrix" << std::endl;
	
	// After a while, combine to initialize TMMC collection matrix
	sys.startTMMC (sys.totNMax(), sys.totNMin());
	
	std::cout << "Assigning initial macrostate density guess from Wang-Landau portion" << std::endl;
	
	// Initial guess from Wang-Landau density of states
	sys.getTMMCBias()->setLnPI(sys.getWALABias()->getlnPI());
	
	// actually this should run until all elements of the collection matrix have been populated
	bool fullyVisited = false;
	while (!fullyVisited) {
		for (unsigned int move = 0; move < wlSweepSize; ++move) {
			try {
				usedMovesEq.makeMove(sys);
			} catch (customException &ce) {
				std::cerr << ce.what() << std::endl;
				exit(SYS_FAILURE);
			}
			
			// record U
			sys.recordU();
		}
				
		// Check if bias has flattened out
		flat = sys.getWALABias()->evaluateFlatness();
		if (flat) {
			// If flat, need to reset H and reduce lnF
			sys.getWALABias()->iterateForward();
			
			time_t rawtime_tmp;
			time (&rawtime_tmp);
			struct tm * timeinfo_tmp;
			timeinfo_tmp = localtime (&rawtime_tmp);
			char dummy_tmp [80];
			strftime (dummy_tmp,80,"%d/%m/%Y %H:%M:%S",timeinfo_tmp);
			std::cout << "lnF = " << sys.getWALABias()->lnF() << " at " << dummy_tmp << std::endl;	
			
			// Periodically write out checkpoints
			sys.getWALABias()->print("wl-crossover-Checkpoint", true);
			sys.getTMMCBias()->print("tmmc-crossover-Checkpoint", true);
		}

		// Check if collection matrix is ready to take over, not necessarily at points where WL is flat
		fullyVisited = sys.getTMMCBias()->checkFullyVisited();
	}

	std::cout << "Switching over to TMMC completely, ending Wang-Landau" << std::endl;
	sys.getTMMCBias()->print("tmmc-beginning-Checkpoint", true);
	
	// Switch over to TMMC completely
	sys.stopWALA();
	
	std::cout << "Beginning TMMC" << std::endl;
	
	for (unsigned int sweep = 0; sweep < totalTMMCSweeps; ++sweep) {
		for (unsigned int move = 0; move < tmmcSweepSize; ++move) {
			try {
				usedMovesPr.makeMove(sys);
			} catch (customException &ce) {
				std::cerr << ce.what() << std::endl;
				exit(SYS_FAILURE);
			}
			
			// record U
			sys.recordU();
		}
		
		if (sweep%(totalTMMCSweeps/100) == 0) {
			time_t rawtime_tmp;
			time (&rawtime_tmp);
			struct tm * timeinfo_tmp;
			timeinfo_tmp = localtime (&rawtime_tmp);
			char dummy_tmp [80];
			strftime (dummy_tmp,80,"%d/%m/%Y %H:%M:%S",timeinfo_tmp);
			std::cout << "Finished " << sweep << "/" << totalTMMCSweeps << " total TMMC sweeps at " << dummy_tmp << std::endl;
		}
		
		// Update biasing function from collection matrix
		sys.getTMMCBias()->calculatePI();
		
		// Periodically write out checkpoints
		sys.getTMMCBias()->print("tmmc-Checkpoint", true);
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
    sys.printSnapshot(restart_file+".xyz", "last configuration");
    
    // Print out energy histogram
    sys.printU("energyHistogram");
    
    // Print out final macrostate distribution
    sys.getTMMCBias()->print("final", false);
	
	// Free pair potential pointers
	for (unsigned int i = 0; i < ppotArray.size(); ++i) {
		delete ppotArray[i];
	}
    ppotArray.clear();
    
    // Finished
    time (&rawtime);
    timeinfo = localtime (&rawtime);
    strftime (timestamp,80,"%d/%m/%Y %H:%M:%S",timeinfo);
    std::cout << "Finished simulation at " << timestamp << std::endl;
	return SAFE_EXIT;
}
