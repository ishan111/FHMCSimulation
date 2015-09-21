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
 * Compile with -DNETCDF_CAPABLE if netCDF libraries are installed and can be compiled against.  Data will be output to these arrays
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
	std::cout << "Parsed " << argv[1] << std::endl;
	
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
	
	int Mtot = 1;
	if (doc.HasMember("num_expanded_states")) {
		assert(doc["num_expanded_states"].IsInt());
		Mtot = doc["num_expanded_states"].GetInt();
	}	
	
	simSystem sys (doc["num_species"].GetInt(), doc["beta"].GetDouble(), sysBox, sysMu, sysMax, sysMin, Mtot);

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
	
	std::string restart_file = "";
	if (doc.HasMember("restart_file")) {
		assert(doc["restart_file"].IsString());
		restart_file = doc["restart_file"].GetString();
	}

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
	
	double lnF_start = 1.0; // default for lnF_start
	if (doc.HasMember("lnF_start")) {
		assert(doc["lnF_start"].IsNumber());
		lnF_start = doc["lnF_start"].GetDouble(); // bounds are checked later
	}
	
	double lnF_end = 2.0e-18; // default for lnF_end
	if (doc.HasMember("lnF_end")) {
		assert(doc["lnF_end"].IsNumber());
		lnF_end = doc["lnF_end"].GetDouble();
		if (lnF_end >= 1.0) {
			std::cerr << "Terminal lnF factor for Wang-Landau must be < 1" << std::endl;
			exit(SYS_FAILURE);
		}
	}
	if (lnF_end >= lnF_start) {
		std::cerr << "lnF_end must be < lnF_start for Wang-Landau to proceed forward" << std::endl;
		exit(SYS_FAILURE);
	}
	
	bool restartFromWALA = false;
	std::string restartFromWALAFile = "";
	if (doc.HasMember("restart_from_wala_lnPI")) {
		assert(doc["restart_from_wala_lnPI"].IsString());
		restartFromWALAFile = doc["restart_from_wala_lnPI"].GetString();
		if (restartFromWALAFile != "") {
			restartFromWALA = true;
		}
	}
	
	// restarting from TMMC overrides WL by skipping that portion altogether
	bool restartFromTMMC = false;
	std::string restartFromTMMCFile = "";
	if (doc.HasMember("restart_from_tmmc_C")) {
		assert(doc["restart_from_tmmc_C"].IsString());
		restartFromTMMCFile = doc["restart_from_tmmc_C"].GetString();
		if (restartFromTMMCFile != "") {
			restartFromTMMC = true;
		}
	}
	
	// number of times the TMMC C matrix has to be traversed during the WALA --> TMMC crossover
	int nCrossoverVisits = 5; // default
	if (doc.HasMember("num_crossover_visits")) {
		assert(doc["num_crossover_visits"].IsInt());
		nCrossoverVisits = doc["num_crossover_visits"].GetInt();
		if (nCrossoverVisits < 1) {
			std::cerr << "Must allow the collection matrix to be traversed at least once in the crossover from Wang-Landau to TMMC" << std::cerr;
			exit(SYS_FAILURE);
		}
	}	

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

	std::vector < pairPotential* > ppotArray (sys.nSpecies()*(sys.nSpecies()-1)/2 + sys.nSpecies());
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
			std::vector < double > params (doc[dummy.c_str()].Size()+1, 0);
			for (unsigned int k = 0; k < params.size()-1; ++k) {
				assert(doc[dummy.c_str()][k].IsNumber());
				params[k] = doc[dummy.c_str()][k].GetDouble();
			}	
			params[params.size()-1] = Mtot;

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
	moves usedMovesEq (sys.getTotalM()), usedMovesPr (sys.getTotalM());
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
		newTranslate.setMaxDisplacement (maxEqD[i], sys.box());
		eqTranslations[i] = newTranslate;
		usedMovesEq.addMove (&eqTranslations[i], probEqDisp[i]);
		
		insertParticle newIns2 (i, "insert");
		prInsertions[i] = newIns2;
		usedMovesPr.addMove (&prInsertions[i], probPrInsDel[i]);
		
		deleteParticle newDel2 (i, "delete");
		prDeletions[i] = newDel2;
		usedMovesPr.addMove (&prDeletions[i], probPrInsDel[i]);		
		
		translateParticle newTranslate2 (i, "translate");
		newTranslate2.setMaxDisplacement (maxPrD[i], sys.box());
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
	
	// Read from restart file if specified
	if (restart_file != "") {
		std::cout << "Reading initial configuration from " << restart_file << std::endl;
		try {
			sys.readRestart(restart_file);
		} catch (customException &ce) {
			std::cerr << ce. what() << std::endl;
			for (unsigned int i = 0; i < ppotArray.size(); ++i) {
				delete ppotArray[i];
			}
			ppotArray.clear();
			exit(SYS_FAILURE);
		}
	} else if (restart_file == "" && sys.totNMin() > 0) {
        	std::cout << "Automatically generating the initial configuration" << std::endl;
		
		// have to generate initial configuration manually - start with mu = INF
        	std::vector < double > initMu (doc["num_species"].GetInt(), 1.0e9);
        	simSystem initSys (doc["num_species"].GetInt(), 1/50., sysBox, initMu, sysMax, sysMin, Mtot); // beta =  1/T, so low beta to have high T

        	// add the same potentials
        	int initInd = 0;
        	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
            		for (unsigned int j = i; j < sys.nSpecies(); ++j) {
                		initSys.addPotential (i, j, ppotArray[initInd], true); // default use of cell list even if otherwise not in main simulation
                		initInd++;
            		}
        	}

		// iteratively add each individual species, assume we want an equimolar mixture to start from
		//for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		for (unsigned int i = 0; i < 1; ++i) {
			std::cout << "Initializing species " << i << " configurations" << std::endl;
			
			// insert this species i
			moves initMove (initSys.getTotalM());
			insertParticle initIns (i, "insert");
            		initMove.addMove (&initIns, 1.0);
			std::vector < translateParticle > initTrans (i+1);

			// also add displacment moves for all species present
			for (unsigned int j = 0; j <= i; ++j) {
				std::cout << "Added translation moves for initialization of species " << j << std::endl;
				translateParticle newTrans (j, "translate");
                		newTrans.setMaxDisplacement (1.0, initSys.box()); // allow large displacements if necessary
				initTrans[j] = newTrans;
          			initMove.addMove (&initTrans[j], 2.0); // move more than insert so this relaxes better (qualitative observation)
			}

			// now do simuation until within proper range
			int targetNum = sys.totNMin();/*/sys.nSpecies();
			if (i == sys.nSpecies() - 1) {
				// to account for integer rounding
				targetNum = sys.totNMin() - (sys.nSpecies()-1)*(sys.totNMin()/sys.nSpecies());
			}*/
			std::cout << "Target number = " << targetNum << std::endl;
			int tmpCounter = 0, statusPrint = 10e6;
			while (initSys.numSpecies[i] < targetNum) {
				// this creates an equimolar mixture
				try {
                    			initMove.makeMove(initSys);
                		} catch (customException &ce) {
                    			for (unsigned int i = 0; i < ppotArray.size(); ++i) {
                        			delete ppotArray[i];
                    			}
                    			ppotArray.clear();
                    			std::cerr << "Failed to create an initial configuration: " << ce.what() << std::endl;
                    			exit(SYS_FAILURE);
                		}
				tmpCounter++;
				if (tmpCounter%statusPrint == 0) {
					tmpCounter = 0;
					std::cout << "Grew " << initSys.numSpecies[i] << " atoms of type " << i << " so far" << std::endl;
				}
			}
		}

		// print snapshot from initSys
		initSys.printSnapshot("auto-init.xyz", "auto-generated initial configuration");

		// read into sys
		try {
			sys.readRestart("auto-init.xyz");
       		} catch (customException &ce) {
            		std::cerr << "Failed to read auto-generated initialization file: " << ce. what() << std::endl;
            		for (unsigned int i = 0; i < ppotArray.size(); ++i) {
                		delete ppotArray[i];
            		}
            		ppotArray.clear();
            		exit(SYS_FAILURE);
        	}
   	}


	bool highSnap = false, lowSnap = false;
					
	if (!restartFromTMMC) {
		std::cout << "Beginning Wang-Landau portion" << std::endl;
		
		// Initially do a WL simulation
		bool flat = false;
		double lnF = lnF_start;
		sys.startWALA (lnF, g, s, sys.getTotalM()); //!< Using Shen and Errington method this syntax is same for single and multicomponent
		
		time_t rawtime_t;
		time (&rawtime_t);
		struct tm * timeinfo_t;
		timeinfo_t = localtime (&rawtime_t);
		char dummy_t [80];
		strftime (dummy_t,80,"%d/%m/%Y %H:%M:%S",timeinfo_t);
		std::cout << "Initial lnF = " << lnF_start << " at " << dummy_t << std::endl;
			
		if (restartFromWALA) {
			try {
				sys.getWALABias()->readlnPI(restartFromWALAFile);
			} catch (customException &ce) {
				for (unsigned int i = 0; i < ppotArray.size(); ++i) {
					delete ppotArray[i];
				}
				ppotArray.clear();
				std::cerr << ce.what() << std::endl;
				exit(SYS_FAILURE);
			}
			std::cout << "Read initial lnPI for Wang-Lnadau from " << restartFromWALAFile << std::endl;
		}
		
		long long int counter = 0;
		while (lnF > lnF_end) {
			for (unsigned int move = 0; move < wlSweepSize; ++move) {
				try {
					usedMovesEq.makeMove(sys);
				} catch (customException &ce) {
					for (unsigned int i = 0; i < ppotArray.size(); ++i) {
                                                delete ppotArray[i];
                                        }
                                        ppotArray.clear();
					std::cerr << ce.what() << std::endl;
					exit(SYS_FAILURE);
				}	
			}

			// Check if bias has flattened out
			flat = sys.getWALABias()->evaluateFlatness();
			if (flat) {
				counter++;
				
				// Periodically write out checkpoints - before iterateForward() which destroys H matrix
				sys.getWALABias()->print("wl-Checkpoint-"+sstr(counter), true);
						
				// if flat, need to reset H and reduce lnF
				sys.getWALABias()->iterateForward();
				lnF = sys.getWALABias()->lnF();
				flat = false;
				
				time_t rawtime_tmp;
				time (&rawtime_tmp);
				struct tm * timeinfo_tmp;
				timeinfo_tmp = localtime (&rawtime_tmp);
				char dummy_tmp [80];
				strftime (dummy_tmp,80,"%d/%m/%Y %H:%M:%S",timeinfo_tmp);
				std::cout << "lnF = " << lnF << " at " << dummy_tmp << std::endl;
			}
		
			// also check to print out snapshots with 10% of bounds to be used for other restarts
			if (!highSnap) {
				if (sys.getTotN() > sys.totNMax() - (sys.totNMax()-sys.totNMin())*0.1) {
					sys.printSnapshot("high.xyz", "snapshot near upper bound");
					highSnap = true;
				}
			}
			if (!lowSnap) {
				if (sys.getTotN() < sys.totNMin() + (sys.totNMax()-sys.totNMin())*0.1 && sys.getTotN() > 0) {
					sys.printSnapshot("low.xyz", "snapshot near lower bound");
					lowSnap = true;
				}
			}
		}
	
		std::cout << "Crossing over to build TMMC matrix" << std::endl;
	
		// After a while, combine to initialize TMMC collection matrix
		sys.startTMMC (tmmcSweepSize, sys.getTotalM());
	
		// actually this should run until all elements of the collection matrix have been populated
		int timesFullyVisited = 0;
		while (timesFullyVisited < nCrossoverVisits) { 
			for (unsigned int move = 0; move < wlSweepSize; ++move) {
				try {
					usedMovesEq.makeMove(sys);
				} catch (customException &ce) {
					std::cerr << ce.what() << std::endl;
					for (unsigned int i = 0; i < ppotArray.size(); ++i) {
						delete ppotArray[i];
					}
					ppotArray.clear();
					exit(SYS_FAILURE);
				}
			}
			
			// Check if collection matrix is ready to take over, not necessarily at points where WL is flat
			if (sys.getTMMCBias()->checkFullyVisited()) {
				sys.getTMMCBias()->iterateForward (); // reset the counting matrix and increment total sweep number
				timesFullyVisited = sys.getTMMCBias()->numSweeps(); 	
				sys.getWALABias()->print("wl-crossover-Checkpoint-"+sstr(timesFullyVisited), true);
				sys.getTMMCBias()->print("tmmc-crossover-Checkpoint-"+sstr(timesFullyVisited), true);

				time_t rawtime_tmp;
				time (&rawtime_tmp);
				struct tm * timeinfo_tmp;
				timeinfo_tmp = localtime (&rawtime_tmp);
				char dummy_tmp [80];
				strftime (dummy_tmp,80,"%d/%m/%Y %H:%M:%S",timeinfo_tmp);
				std::cout << "Times C fully visited = " << timesFullyVisited << " at " << dummy_tmp << std::endl;
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
			}
		}

		// Switch over to TMMC completely
		sys.stopWALA();
		
		std::cout << "Switching over to TMMC completely, ending Wang-Landau" << std::endl;
		sys.getTMMCBias()->calculatePI();
		sys.getTMMCBias()->print("tmmc-beginning-Checkpoint", true);
	}
	
	std::cout << "Beginning TMMC" << std::endl;
	if (restartFromTMMC) {
		sys.startTMMC (tmmcSweepSize, sys.getTotalM()); // this was otherwise started during the crossover phase if WL was used
		try {
			sys.getTMMCBias()->readC(restartFromTMMCFile); // read collection matrix
		} catch (customException& ce) {
			sys.stopTMMC(); // deallocate
			std::cerr << "Failed to initialize from TMMC collection matrix: " << ce.what() << std::endl;
			for (unsigned int i = 0; i < ppotArray.size(); ++i) {
				delete ppotArray[i];
			}
			ppotArray.clear();
			exit(SYS_FAILURE);
		}
		sys.getTMMCBias()->calculatePI();
		std::cout << "Restarted TMMC from collection matrix from " << restartFromTMMCFile << std::endl;
	}
	
	long long int sweep = 0;
	int sweepPrint = totalTMMCSweeps, numSweepSnaps = 100, printCounter = 0;
	if (totalTMMCSweeps > numSweepSnaps) {
		sweepPrint /= numSweepSnaps;
	}
	while (sweep < totalTMMCSweeps) {
		bool done = false;
		unsigned long long int counter = 0;
		unsigned long long int checkPoint = tmmcSweepSize*(sys.totNMax() - sys.totNMin() + 1)*3; // how often to check full traversal of collection matrix
		// perform a sweep
		while (!done) {
			try {
				usedMovesPr.makeMove(sys);
			} catch (customException &ce) {
				std::cerr << ce.what() << std::endl;
				for (unsigned int i = 0; i < ppotArray.size(); ++i) {
					delete ppotArray[i];
				}
				ppotArray.clear();
				exit(SYS_FAILURE);
			}
			
			// only record properties of the system when it is NOT in an intermediate state
			if (sys.getCurrentM() == 0) {
				// record U
				sys.recordU();
			
				// record composition
				sys.recordComposition();
			}
	
			// check if sweep is done
			if (counter%checkPoint == 0) {
				done = sys.getTMMCBias()->checkFullyVisited();
				counter = 0;
			}

			counter++;
		}
		
		sys.getTMMCBias()->iterateForward (); // reset the counting matrix and increment total sweep number
		sweep++;
				
		time_t rawtime_tmp;
		time (&rawtime_tmp);
		struct tm * timeinfo_tmp;
		timeinfo_tmp = localtime (&rawtime_tmp);
		char dummy_tmp [80];
		strftime (dummy_tmp,80,"%d/%m/%Y %H:%M:%S",timeinfo_tmp);
		std::cout << "Finished " << sweep << "/" << totalTMMCSweeps << " total TMMC sweeps at " << dummy_tmp << std::endl;
		
		// Update biasing function from collection matrix
		sys.getTMMCBias()->calculatePI();
		
		// Periodically write out checkpoints
		if (sweep%sweepPrint == 0) {
			printCounter++;
			sys.getTMMCBias()->print("tmmc-Checkpoint-"+sstr(printCounter), true);
			sys.printU("energy-Checkpoint-"+sstr(printCounter));
			sys.printComposition("composition-Checkpoint-"+sstr(printCounter));
		}

		// also check to print out snapshots with 10% of bounds to be used for other restarts
		if (!highSnap) {
			if (sys.getTotN() > sys.totNMax() - (sys.totNMax()-sys.totNMin())*0.1) {
				sys.printSnapshot("high.xyz", "snapshot near upper bound");
				highSnap = true;
			}
		}
		if (!lowSnap) {
			if (sys.getTotN() < sys.totNMin() + (sys.totNMax()-sys.totNMin())*0.1 && sys.getTotN() > 0) {
				sys.printSnapshot("low.xyz", "snapshot near lower bound");
				lowSnap = true;
			}
		}
	}
		
	// Sanity checks
	if (sys.nSpecies() != sys.atoms.size()) {
        	std::cerr << "Error: Number of components changed throughout simulation" << std::endl;
		for (unsigned int i = 0; i < ppotArray.size(); ++i) {
			delete ppotArray[i];
		}
		ppotArray.clear();
        	exit(SYS_FAILURE);
    	}
	if (sys.getTotalM() > 1) {
		if (sys.getFractionalAtom()->mState != sys.getCurrentM()) {
			std::cerr << "Expanded ensemble state deviates between atom ("+sstr(sys.getFractionalAtom()->mState)+") and system log ("+sstr(sys.getCurrentM())+")" << std::endl;
			exit(SYS_FAILURE);
			for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
				int end = sys.numSpecies[i];
				if (i == sys.getFractionalAtomType()) {
					end++;			
				}
				for (unsigned int j = 0; j < end; ++j) {
					if (&sys.atoms[i][j] != sys.getFractionalAtom()) {
						if (sys.atoms[i][j].mState != 0) {
							std::cerr << "Atom ("+sstr(i)+", "+sstr(j)+") has non-zero expanded ensemble state ("+sstr(sys.atoms[i][j].mState)+")" << std::endl;	
							exit(SYS_FAILURE);				
						}				
					} else {
						if (sys.atoms[i][j].mState != sys.getCurrentM()) {
							std::cerr << "Fractional atom ("+sstr(i)+", "+sstr(j)+")'s expanded ensemble state ("+sstr(sys.atoms[i][j].mState)+") does not match system's ("+sstr(sys.getCurrentM())+")" << std::endl;	
							exit(SYS_FAILURE);				
						}
					}								
				}
			}
		}
	} else {
		for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
			for (unsigned int j = 0; j < sys.numSpecies[i]; ++j) {
				if (sys.atoms[i][j].mState != 0) {
					std::cerr << "Atom ("+sstr(i)+", "+sstr(j)+") has non-zero expanded ensemble state ("+sstr(sys.atoms[i][j].mState)+")" << std::endl;	
					exit(SYS_FAILURE);				
				}												
			}
		}	
	}

    	// Report move statistics for final TMMC ("production") stage
	char statName [80];
	strftime (statName,80,"%Y_%m_%d_%H_%M_%S-stats.log",timeinfo);
	std::ofstream statFile (statName);
    	std::vector < std::vector < double > > stats = usedMovesPr.reportMoveStatistics();
    	statFile << " ---------- Move Statistics --------- " << std::endl << " Move\t\% Success" << std::endl;
    	for (unsigned int i = 0; i < stats.size(); ++i) {
        	double prod = 1.0;
        	for (unsigned int j = 0; j < stats[i].size(); ++j) {
            		prod *= stats[i][j];
            		statFile << usedMovesPr.includedMoves()[i]->myName() << " (from M = " << j << ")\t" << stats[i][j]*100.0 << std::endl;
        	}
        	if (stats[i].size() > 1) {
            		statFile << "-------------------------------------\nProduct of percentages (%) = " << prod*100 << "\n-------------------------------------" << std::endl;
        	}
    	}
    	statFile << std::endl;
    	statFile.close();
	
 	// print out restart file (xyz)
    	sys.printSnapshot("final.xyz", "last configuration");
    
    	// Print out energy histogram
    	sys.printU("energyHistogram");
    
   	// Print out composition histogram
    	sys.printComposition("compositionHistogram");
    
    	// Print out final macrostate distribution
    	sys.getTMMCBias()->print("final", false);
	
	// Still allow for printing of all data, even if there is an error, in order to interrogate the results anyway
	const double tol = 1.0e-6;
	const double scratchEnergy = sys.scratchEnergy(), incrEnergy = sys.energy();
    	if (fabs(scratchEnergy - incrEnergy) > tol) {
        	std::cerr << "Error: scratch energy calculation = " << std::setprecision(20) << scratchEnergy << ", but incremental = " << std::setprecision(20) << incrEnergy << ", |diff| = " << std::setprecision(20) << fabs(scratchEnergy - incrEnergy) << std::endl;
        	for (unsigned int i = 0; i < ppotArray.size(); ++i) {
            		delete ppotArray[i];
        	}
        	ppotArray.clear();
        	exit(SYS_FAILURE);
    	} else {
        	std::cout << "Passed: Final scratch energy - incremental = " << std::setprecision(20) << scratchEnergy - incrEnergy << std::endl;
    	}

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
