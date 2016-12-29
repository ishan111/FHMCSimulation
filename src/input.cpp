#include "input.h"

/*!
 * Check the following bounds on system.  If any fail, an emergency exit is performed.
 * 1. Check all pair potentials are set.
 * 2. Check that rcut < L/2 for all potentials
 *
 * \param [in] sys System to check bounds on
 */
void checkBounds (simSystem &sys) {
    // check all pair potentials have been set and all r_cut < L/2
    double minL = sys.box()[0];
    for (unsigned int i = 1; i < 3; ++i) {
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
}

/*!
 * Parse a json input file and initialize system object accordingly.
 *
 * \params [in] filename Input JSON document's filename
 * \params [in] usedMovesEq Pointer to move object that will be used during "equilibration" (WL)
 * \params [in] usedMovesPr Pointer to move object that will be used during "production" (TMMC)
 */
simSystem initialize (const std::string filename, moves* usedMovesEq, moves* usedMovesPr) {
	rapidjson::Document doc;
    parseJson (filename, doc);

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

	double duh = 10.0;
	if (doc.HasMember("delta_u_hist")) {
		assert(doc["delta_u_hist"].IsNumber());
		duh = doc["delta_u_hist"].GetDouble();
	}

	int max_order = 2;
	if (doc.HasMember("max_order")) {
		assert(doc["max_order"].IsNumber());
		max_order = doc["max_order"].GetInt();
	}

	bool use_ke = false;
	if (doc.HasMember("use_ke")) {
		assert(doc["use_ke"].IsBool());
		use_ke = doc["use_ke"].GetBool();
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

	simSystem sys (doc["num_species"].GetInt(), doc["beta"].GetDouble(), sysBox, sysMu, sysMax, sysMin, Mtot, duh, max_order);
	if (use_ke) {
		sys.toggleKE();
		if (sys.addKECorrection() == false) {
			throw customException ("Unable to set KE flag");
		}
	}

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

	assert(doc.HasMember("tmmc_sweep_size"));
	assert(doc["tmmc_sweep_size"].IsNumber());
	double tmpT = doc["tmmc_sweep_size"].GetDouble(); // possibly in scientific notation
	sys.tmmcSweepSize = tmpT; // convert

	assert(doc.HasMember("total_tmmc_sweeps"));
	assert(doc["total_tmmc_sweeps"].IsNumber());
	double tmpS = doc["total_tmmc_sweeps"].GetDouble(); // possibly in scientific notation
	sys.totalTMMCSweeps = tmpS; // convert

	assert(doc.HasMember("wala_sweep_size"));
	assert(doc["wala_sweep_size"].IsNumber());
	double tmpW = doc["wala_sweep_size"].GetDouble(); // possibly in scientific notation
	sys.wlSweepSize = tmpW; // convert

	assert(doc.HasMember("wala_g"));
	assert(doc["wala_g"].IsNumber());
	sys.wala_g = doc["wala_g"].GetDouble();

	assert(doc.HasMember("wala_s"));
	assert(doc["wala_s"].IsNumber());
	sys.wala_s = doc["wala_s"].GetDouble();

	if (doc.HasMember("lnF_start")) {
		assert(doc["lnF_start"].IsNumber());
		sys.lnF_start = doc["lnF_start"].GetDouble(); // bounds are checked later
	}

	if (doc.HasMember("lnF_end")) {
		assert(doc["lnF_end"].IsNumber());
		sys.lnF_end = doc["lnF_end"].GetDouble();
		if (sys.lnF_end >= 1.0) {
			std::cerr << "Terminal lnF factor for Wang-Landau must be < 1" << std::endl;
			exit(SYS_FAILURE);
		}
	}
	if (sys.lnF_end >= sys.lnF_start) {
		std::cerr << "lnF_end must be < lnF_start for Wang-Landau to proceed forward" << std::endl;
		exit(SYS_FAILURE);
	}

	sys.restartFromWALA = false;
	sys.restartFromWALAFile = "";
	if (doc.HasMember("restart_from_wala_lnPI")) {
		assert(doc["restart_from_wala_lnPI"].IsString());
		sys.restartFromWALAFile = doc["restart_from_wala_lnPI"].GetString();
		if (sys.restartFromWALAFile != "") {
			sys.restartFromWALA = true;
		}
	}

	// restarting from TMMC overrides WL by skipping that portion altogether
	sys.restartFromTMMC = false;
	sys.restartFromTMMCFile = "";
	if (doc.HasMember("restart_from_tmmc_C")) {
		assert(doc["restart_from_tmmc_C"].IsString());
		sys.restartFromTMMCFile = doc["restart_from_tmmc_C"].GetString();
		if (sys.restartFromTMMCFile != "") {
			sys.restartFromTMMC = true;
		}
	}

	// number of times the TMMC C matrix has to be traversed during the WALA --> TMMC crossover
	if (doc.HasMember("num_crossover_visits")) {
		assert(doc["num_crossover_visits"].IsNumber());
		sys.nCrossoverVisits = doc["num_crossover_visits"].GetDouble(); // convert
		if (sys.nCrossoverVisits < 1) {
			std::cerr << "Must allow the collection matrix to be traversed at least once in the crossover from Wang-Landau to TMMC" << std::endl;
			exit(SYS_FAILURE);
		}
	}

    setMoves (sys, doc, usedMovesEq, usedMovesPr);
    setPairPotentials (sys, doc);

    checkBounds (sys);
    std::cout << "System from " << filename << " passed bounds checks at " << getTimeStamp() << std::endl;

    setSystemBarriers (sys, doc);
    std::cout << "Initialized barriers from " << filename << " at " << getTimeStamp() << std::endl;

    std::cout << "Successfully read valid parameters from " << filename << " at " << getTimeStamp() << std::endl;
    return sys;
}

/*!
 * Assign the Monte Carlo moves based on the JSON input file.
 *
 * \param [in] sys Simulation system that has been initialized
 * \param [in] doc JSON document corresponding to input file
 * \params [in] usedMovesEq Pointer to move object that will be used during "equilibration" (WL)
 * \params [in] usedMovesPr Pointer to move object that will be used during "production" (TMMC)
 */
void setMoves (simSystem &sys, const rapidjson::Document &doc, moves* usedMovesEq, moves* usedMovesPr) {
    std::vector < double > ref (sys.nSpecies(), 0);
	std::vector < std::vector < double > > probEqSwap (sys.nSpecies(), ref), probPrSwap (sys.nSpecies(), ref);
	std::vector < double > probPrInsDel (sys.nSpecies(), 0), probPrDisp (sys.nSpecies(), 0);
	std::vector < double > probEqInsDel (sys.nSpecies(), 0), probEqDisp (sys.nSpecies(), 0);
	std::vector < double > maxPrD (sys.nSpecies(), 0), maxEqD (sys.nSpecies(), 0);
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		std::string dummy = "prob_pr_ins_del_" + std::to_string(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		probPrInsDel[i] = doc[dummy.c_str()].GetDouble();
	}

	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		std::string dummy = "prob_pr_displace_" + std::to_string(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		probPrDisp[i] = doc[dummy.c_str()].GetDouble();
		dummy = "max_pr_displacement_" + std::to_string(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		maxPrD[i] = doc[dummy.c_str()].GetDouble();
	}

	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		std::string dummy = "prob_eq_ins_del_" + std::to_string(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		probEqInsDel[i] = doc[dummy.c_str()].GetDouble();
	}

	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		std::string dummy = "prob_eq_displace_" + std::to_string(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		probEqDisp[i] = doc[dummy.c_str()].GetDouble();
		dummy = "max_eq_displacement_" + std::to_string(i+1);
		assert(doc.HasMember(dummy.c_str()));
		assert(doc[dummy.c_str()].IsNumber());
		maxEqD[i] = doc[dummy.c_str()].GetDouble();
	}

	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		for (unsigned int j = i+1; j < sys.nSpecies(); ++j) {
			std::string name1 = "prob_pr_swap_"+std::to_string(i+1)+"_"+std::to_string(j+1);
			std::string name2 = "prob_pr_swap_"+std::to_string(j+1)+"_"+std::to_string(i+1);
			std::string moveName = "";
			bool foundIJ = false;
			if (doc.HasMember(name1.c_str())) {
				moveName = name1;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && !foundIJ) {
				moveName = name2;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && foundIJ) {
				std::cerr << "Input file doubly specifies production swap move probability for species pair ("+std::to_string(i+1)+", "+std::to_string(j+1)+")" << std::endl;
				exit(SYS_FAILURE);
			} else {
				std::cerr << "Input file does not specify production swap move probability for species pair ("+std::to_string(i+1)+", "+std::to_string(j+1)+")" << std::endl;
				exit(SYS_FAILURE);
			}
			assert(doc[moveName.c_str()].IsNumber());
			probPrSwap[i][j] = doc[moveName.c_str()].GetDouble();
			probPrSwap[j][i] = doc[moveName.c_str()].GetDouble();
		}
	}

	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		for (unsigned int j = i+1; j < sys.nSpecies(); ++j) {
			std::string name1 = "prob_eq_swap_"+std::to_string(i+1)+"_"+std::to_string(j+1);
			std::string name2 = "prob_eq_swap_"+std::to_string(j+1)+"_"+std::to_string(i+1);
			std::string moveName = "";
			bool foundIJ = false;
			if (doc.HasMember(name1.c_str())) {
				moveName = name1;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && !foundIJ) {
				moveName = name2;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && foundIJ) {
				std::cerr << "Input file doubly specifies equilibration swap move probability for species pair ("+std::to_string(i+1)+", "+std::to_string(j+1)+")" << std::endl;
				exit(SYS_FAILURE);
			} else {
				std::cerr << "Input file does not specify equilibration swap move probability for species pair ("+std::to_string(i+1)+", "+std::to_string(j+1)+")" << std::endl;
				exit(SYS_FAILURE);
			}
			assert(doc[moveName.c_str()].IsNumber());
			probEqSwap[i][j] = doc[moveName.c_str()].GetDouble();
			probEqSwap[j][i] = doc[moveName.c_str()].GetDouble();
		}
	}

    usedMovesEq->setM(sys.getTotalM());
    usedMovesPr->setM(sys.getTotalM());
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
        usedMovesEq->addInsert(i, probEqInsDel[i]);
        usedMovesPr->addInsert(i, probPrInsDel[i]);

        usedMovesEq->addDelete(i, probEqInsDel[i]);
        usedMovesPr->addDelete(i, probPrInsDel[i]);

        usedMovesEq->addTranslate(i, probEqDisp[i], maxEqD[i], sys.box());
        usedMovesPr->addTranslate(i, probPrDisp[i], maxPrD[i], sys.box());

        for (unsigned int j = i+1; j < sys.nSpecies(); ++j) {
            usedMovesEq->addSwap(i, j, probEqSwap[i][j]);
            usedMovesPr->addSwap(i, j, probPrSwap[i][j]);
        }
    }
}

/*!
 * Assign the pair potentials based on the JSON input file.
 *
 * \param [in] sys Simulation system that has been initialized
 * \param [in] doc JSON document corresponding to input file
 */
void setPairPotentials (simSystem &sys, const rapidjson::Document &doc) {
    int Mtot = 1;
    if (doc.HasMember("num_expanded_states")) {
        Mtot = doc["num_expanded_states"].GetInt();
    }

    //std::vector < pairPotential* > ppotArray (sys.nSpecies()*(sys.nSpecies()-1)/2 + sys.nSpecies());
	std::vector < std::string > ppotType (sys.nSpecies()*(sys.nSpecies()-1)/2 + sys.nSpecies());
	int ppotTypeIndex = 0;
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		for (unsigned int j = i; j < sys.nSpecies(); ++j) {
			std::string name1 = "ppot_"+std::to_string(i+1)+"_"+std::to_string(j+1), name2 = "ppot_"+std::to_string(j+1)+"_"+std::to_string(i+1);
			std::string ppotName = "", dummy = "";
			bool foundIJ = false;
			if (doc.HasMember(name1.c_str())) {
				ppotName = name1;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && !foundIJ) {
				ppotName = name2;
				foundIJ = true;
			} else if (doc.HasMember(name2.c_str()) && foundIJ) {
				std::cerr << "Input file doubly specifies pair potential for species pair ("+std::to_string(i+1)+", "+std::to_string(j+1)+")" << std::endl;
				exit(SYS_FAILURE);
			} else {
				std::cerr << "Input file does not specify pair potential for species pair ("+std::to_string(i+1)+", "+std::to_string(j+1)+")" << std::endl;
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

            sys.addPotential(i, j, ppotType[ppotTypeIndex], params, useCellList);
            sys.ppot[i][j]->savePotential(ppotName+".dat", 0.01, 0.01);

			ppotTypeIndex++;
		}
	}
}

/*!
 * Setup a system's initial configuration as necessary.
 * 1. If "restart_file" in input json file, read initial config from there.
 * 2. In not, randomly generate initial configuration.
 *
 * \param [in] sys System to initialize
 * \param [in] filename Input JSON filename
 */
void setup (simSystem &sys, const std::string filename) {
    rapidjson::Document doc;
    parseJson (filename, doc);

    std::string restart_file = "";
	if (doc.HasMember("restart_file")) {
		assert(doc["restart_file"].IsString());
		restart_file = doc["restart_file"].GetString();
	}

	std::vector < double > sysBox = sys.box();

    double duh = 10.0;
	if (doc.HasMember("delta_u_hist")) {
		assert(doc["delta_u_hist"].IsNumber());
		duh = doc["delta_u_hist"].GetDouble();
	}

	int max_order = 2;
	if (doc.HasMember("max_order")) {
		assert(doc["max_order"].IsNumber());
		max_order = doc["max_order"].GetInt();
	}

	bool use_ke = sys.addKECorrection();
    int Mtot = sys.getTotalM();

	std::vector < double > sysMu (doc["mu"].Size(), 0);
	for (rapidjson::SizeType i = 0; i < doc["mu"].Size(); ++i) {
		sysMu[i] = doc["mu"][i].GetDouble();
	}
	std::vector < int > sysMax (doc["max_N"].Size(), 0);
	for (rapidjson::SizeType i = 0; i < doc["max_N"].Size(); ++i) {
		sysMax[i] = doc["max_N"][i].GetInt();
	}
	std::vector < int > sysMin (doc["min_N"].Size(), 0);
	for (rapidjson::SizeType i = 0; i < doc["min_N"].Size(); ++i) {
		sysMin[i] = doc["min_N"][i].GetInt();
	}

    // Read from restart file if specified
	if (restart_file != "") {
		try {
			sys.readConfig(restart_file);
		} catch (customException &ce) {
			std::cerr << ce. what() << std::endl;
		}
	} else if (restart_file == "" && sys.totNMin() > 0) {
		std::cout << "Automatically generating the initial configuration" << std::endl;

		// have to generate initial configuration manually - start with mu = INF
        std::vector < double > initMu (doc["num_species"].GetInt(), 1.0e2);
		simSystem initSys (doc["num_species"].GetInt(), 1/10., sysBox, initMu, sysMax, sysMin, Mtot, duh, max_order); // beta =  1/T, so low beta to have high T
		if (use_ke) {
			initSys.toggleKE();
			if (initSys.addKECorrection() == false) {
				throw customException ("Unable to set KE flag");
			}
		}

        // add the same potentials
        setPairPotentials (initSys, doc);
		setSystemBarriers (initSys, doc);

        std::vector < int > initialization_order (sys.nSpecies(), 0), check_init (sys.nSpecies(), 0);
    	std::vector < double > init_frac (sys.nSpecies(), 1.0);
    	double sum = 0.0;
    	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
    		initialization_order[i] = i;
    		if (i > 0) init_frac[i] = 0.0;
    		sum += init_frac[i];
    	}
    	if (doc.HasMember("init_order")) {
    		assert(doc["init_order"].IsArray());
    		assert(doc["init_order"].Size() == doc["num_species"].GetInt());

    		for (rapidjson::SizeType i = 0; i < doc["init_order"].Size(); ++i) {
    			assert(doc["init_order"][i].IsInt());
    			initialization_order[i] = doc["init_order"][i].GetInt();
    			if (initialization_order[i] < 0 || initialization_order[i] >= sys.nSpecies()) {
    				std::cerr << "Order of initialization goes out of bounds, should include 0 <= i < nSpec" << std::endl;
    				exit(SYS_FAILURE);
    			}
    			if (check_init[initialization_order[i]] != 0) {
    				std::cerr << "Order of initialization repeats itself" << std::endl;
    				exit(SYS_FAILURE);
    			} else {
    				check_init[initialization_order[i]] = 1;
    			}
    		}
    	}
        if (doc.HasMember("init_frac")) {
    		assert(doc["init_frac"].IsArray());
    		assert(doc["init_frac"].Size() == doc["num_species"].GetInt());
    		sum = 0.0;
    		for (rapidjson::SizeType i = 0; i < doc["init_frac"].Size(); ++i) {
    			assert(doc["init_frac"][i].IsNumber());
    			init_frac[i] = doc["init_frac"][i].GetDouble();
    			if (init_frac[i] < 0 || init_frac[i] >= 1.0) {
    				std::cerr << "Initialization fraction out of bounds" << std::endl;
    				exit(SYS_FAILURE);
    			}
    			sum += init_frac[i];
    		}
    	}
    	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
    		init_frac[i] /= sum;
    	}

		// iteratively add each individual species, assume we want an equimolar mixture to start from
		int added = 0;
		for (unsigned int idx = 0; idx < sys.nSpecies(); ++idx) {
			unsigned int i = initialization_order[idx];
			std::cout << "Initializing species " << i << " configurations" << std::endl;

			// insert this species i
			moves initMove (initSys.getTotalM());
            initMove.addInsert(i, 1.0);

			// also add displacment moves for all species present
			for (unsigned int j = 0; j <= idx; ++j) {
				std::cout << "Added translation moves for initialization of species " << initialization_order[j] << std::endl;
                initMove.addTranslate(initialization_order[j], 2.0, 1.0, initSys.box());
			}

			// now do simuation until within proper range
			int targetNum = sys.totNMin()*init_frac[idx];
			if (idx == sys.nSpecies() - 1) {
				// to account for integer rounding
				targetNum = sys.totNMin() - added;
			}
			added += targetNum;

			std::cout << "Target number = " << targetNum << " for species " << i+1 << std::endl;
			int tmpCounter = 0, statusPrint = 10e6;
			while (initSys.numSpecies[i] < targetNum) {
				try {
                    initMove.makeMove(initSys);
                } catch (customException &ce) {
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

		// print snapshot from Reading initial configuration
		initSys.printSnapshot("auto-init.xyz", "auto-generated initial configuration");

		// read into sys
		try {
			sys.readConfig("auto-init.xyz");
		} catch (customException &ce) {
			std::cerr << "Failed to read auto-generated initialization file: " << ce. what() << std::endl;
        }
   	}
}

/*!
 * Initialize the barriers in a system by parsing the input document.  This function is defined separately since it must be done several times.
 *
 * \params [in, out] sys System to initialize with barriers
 * \params [in] doc Input JSON document
 */
void setSystemBarriers (simSystem &sys, const rapidjson::Document &doc) {
	// get Mtot, first from doc, otherwise try sys, but they should be the same
	int Mtot = 1;
	if (doc.HasMember("num_expanded_states")) {
		assert(doc["num_expanded_states"].IsInt());
		Mtot = doc["num_expanded_states"].GetInt();
	} else {
		Mtot = sys.getTotalM();
	}

	// Hard wall (expect parameters: {lb, ub, sigma})
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		bool convention0 = false;
		std::string dummy = "hardWallZ_" + std::to_string(i+1);
		std::vector < double > wallParams (3, 0);
		if (doc.HasMember(dummy.c_str())) {
			assert(doc[dummy.c_str()].IsArray());
			assert(doc[dummy.c_str()].Size() == 3);
			for (unsigned int j = 0; j < 3; ++j) {
				wallParams[j] = doc[dummy.c_str()][j].GetDouble();
			}
			try {
				sys.speciesBarriers[i].addHardWallZ (wallParams[0], wallParams[1], wallParams[2], Mtot);
			} catch (customException &ce) {
				std::cerr << ce.what() << std::endl;
				exit(SYS_FAILURE);
			}
			convention0 = true;
		}
		for (unsigned int j = 1; j <= MAX_BARRIERS_PER_SPECIES; ++j) {
			// alternatively allow multiple walls to specified with a suffix up to a max
			std::string dummy = "hardWallZ_" + std::to_string(i+1) + "_" + std::to_string(j);
            if (doc.HasMember(dummy.c_str())) {
				if (convention0) {
					std::cerr << "Error: multiple barrier naming conventions used for the same species" << std::endl;
					exit(SYS_FAILURE);
                }
				if (doc.HasMember(dummy.c_str())) {
					assert(doc[dummy.c_str()].IsArray());
					assert(doc[dummy.c_str()].Size() == 3);
					for (unsigned int j = 0; j < 3; ++j) {
						wallParams[j] = doc[dummy.c_str()][j].GetDouble();
					}
                    try {
						sys.speciesBarriers[i].addHardWallZ (wallParams[0], wallParams[1], wallParams[2], Mtot);
					} catch (customException &ce) {
						std::cerr << ce.what() << std::endl;
						exit(SYS_FAILURE);
                    }
                }
            }
        }
	}

    // Square well wall (expect parameters: {lb, ub, sigma, range, eps})
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		bool convention0 = false;
		std::string dummy = "squareWellWallZ_" + std::to_string(i+1);
		std::vector < double > wallParams (5, 0);
		if (doc.HasMember(dummy.c_str())) {
			assert(doc[dummy.c_str()].IsArray());
			assert(doc[dummy.c_str()].Size() == 5);
			for (unsigned int j = 0; j < 5; ++j) {
				wallParams[j] = doc[dummy.c_str()][j].GetDouble();
			}
			try {
				sys.speciesBarriers[i].addSquareWellWallZ (wallParams[0], wallParams[1], wallParams[2], wallParams[3], wallParams[4], Mtot);
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
                exit(SYS_FAILURE);
            }
            convention0 = true;
		}
		for (unsigned int j = 1; j <= MAX_BARRIERS_PER_SPECIES; ++j) {
			// alternatively allow multiple walls to specified with a suffix up to a max
			std::string dummy = "squareWellWallZ_" + std::to_string(i+1) + "_" + std::to_string(j);
			if (doc.HasMember(dummy.c_str())) {
				if (convention0) {
					std::cerr << "Error: multiple barrier naming conventions used for the same species" << std::endl;
					exit(SYS_FAILURE);
				}
				if (doc.HasMember(dummy.c_str())) {
					assert(doc[dummy.c_str()].IsArray());
					assert(doc[dummy.c_str()].Size() == 5);
					for (unsigned int j = 0; j < 5; ++j) {
						wallParams[j] = doc[dummy.c_str()][j].GetDouble();
					}
					try {
						sys.speciesBarriers[i].addSquareWellWallZ (wallParams[0], wallParams[1], wallParams[2], wallParams[3], wallParams[4], Mtot);
					} catch (customException &ce) {
						std::cerr << ce.what() << std::endl;
						exit(SYS_FAILURE);
					}
                }
            }
        }
	}

	// cylinderZ (expect parameters: {x, y, radius, width, sigma, eps})
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		bool convention0 = false;
		std::string dummy = "cylinderZ_" + std::to_string(i+1);
		std::vector < double > wallParams (6, 0);
		if (doc.HasMember(dummy.c_str())) {
			assert(doc[dummy.c_str()].IsArray());
			assert(doc[dummy.c_str()].Size() == 6);
			for (unsigned int j = 0; j < 6; ++j) {
				wallParams[j] = doc[dummy.c_str()][j].GetDouble();
			}
			try {
				sys.speciesBarriers[i].addCylinderZ (wallParams[0], wallParams[1], wallParams[2], wallParams[3], wallParams[4], wallParams[5], Mtot);
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
                exit(SYS_FAILURE);
            }
            convention0 = true;
		}
		for (unsigned int j = 1; j <= MAX_BARRIERS_PER_SPECIES; ++j) {
			// alternatively allow multiple walls to specified with a suffix up to a max
			std::string dummy = "cylinderZ_" + std::to_string(i+1) + "_" + std::to_string(j);
			if (doc.HasMember(dummy.c_str())) {
				if (convention0) {
					std::cerr << "Error: multiple barrier naming conventions used for the same species" << std::endl;
					exit(SYS_FAILURE);
				}
				if (doc.HasMember(dummy.c_str())) {
					assert(doc[dummy.c_str()].IsArray());
					assert(doc[dummy.c_str()].Size() == 6);
					for (unsigned int j = 0; j < 6; ++j) {
						wallParams[j] = doc[dummy.c_str()][j].GetDouble();
					}
					try {
						sys.speciesBarriers[i].addCylinderZ (wallParams[0], wallParams[1], wallParams[2], wallParams[3], wallParams[4], wallParams[5], Mtot);
					} catch (customException &ce) {
						std::cerr << ce.what() << std::endl;
						exit(SYS_FAILURE);
					}
                }
            }
        }
	}

    // rightTriangleXZ (expect parameters: {width, theta, lamW, eps, sigma, sep, offset, zbase, top})
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
        bool convention0 = false;
        std::string dummy = "rightTriangleXZ_" + std::to_string(i+1);
        std::vector < double > wallParams (8, 0);
        bool top = false;
        assert(doc.HasMember("box"));
        assert(doc["box"].IsArray());
        assert(doc["box"].Size() == 3);
        std::vector < double > sysBox (3, 0);
        for (rapidjson::SizeType j = 0; j < doc["box"].Size(); ++j) {
            assert(doc["box"][j].IsNumber());
            sysBox[j] = doc["box"][j].GetDouble();
        }
        if (doc.HasMember(dummy.c_str())) {
            assert(doc[dummy.c_str()].IsArray());
            assert(doc[dummy.c_str()].Size() == 9);
            for (unsigned int j = 0; j < 8; ++j) {
                assert (doc[dummy.c_str()][j].IsDouble());
                wallParams[j] = doc[dummy.c_str()][j].GetDouble();
            }
            assert (doc[dummy.c_str()][8].IsBool());
            top = doc[dummy.c_str()][8].GetBool();
            try {
                sys.speciesBarriers[i].addRightTriangleXZ (wallParams[0], wallParams[1], wallParams[2], wallParams[3], wallParams[4], wallParams[5], wallParams[6], sysBox, wallParams[7], top, Mtot);
            } catch (customException &ce) {
                std::cerr << ce.what() << std::endl;
                exit(SYS_FAILURE);
            }
            convention0 = true;
        }
        for (unsigned int j = 1; j <= MAX_BARRIERS_PER_SPECIES; ++j) {
            // alternatively allow multiple walls to specified with a suffix up to a max
            std::string dummy = "rightTriangleXZ_" + std::to_string(i+1) + "_" + std::to_string(j);
            if (doc.HasMember(dummy.c_str())) {
                if (convention0) {
                    std::cerr << "Error: multiple barrier naming conventions used for the same species" << std::endl;
                    exit(SYS_FAILURE);
                }
                if (doc.HasMember(dummy.c_str())) {
                    assert(doc[dummy.c_str()].IsArray());
                    assert(doc[dummy.c_str()].Size() == 9);
                    for (unsigned int k = 0; k < 8; ++k) {
                        assert (doc[dummy.c_str()][k].IsNumber());
                        wallParams[k] = doc[dummy.c_str()][k].GetDouble();
                    }
                    assert (doc[dummy.c_str()][8].IsBool());
                    top = doc[dummy.c_str()][8].GetBool();
                    try {
                        sys.speciesBarriers[i].addRightTriangleXZ (wallParams[0], wallParams[1], wallParams[2], wallParams[3], wallParams[4], wallParams[5], wallParams[6], sysBox, wallParams[7], top, Mtot);
                    } catch (customException &ce) {
                        std::cerr << ce.what() << std::endl;
                        exit(SYS_FAILURE);
                    }
                }
            }
        }
    }
}
