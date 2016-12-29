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
                sendErr("Not all pair potentials are set");
                exit(SYS_FAILURE);
            }
            if (!(sys.ppot[i][j]->rcut() < minL/2.0)) {
                sendErr("Pair potential r_cut for species "+numToStr(i)+", "+numToStr(j)+" is > L/2");
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
            sendErr("Terminal lnF factor for Wang-Landau must be < 1");
			exit(SYS_FAILURE);
		}
	}
	if (sys.lnF_end >= sys.lnF_start) {
        sendErr("lnF_end must be < lnF_start for Wang-Landau to proceed forward");
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
            sendErr("Must allow the collection matrix to be traversed at least once in the crossover from Wang-Landau to TMMC");
			exit(SYS_FAILURE);
		}
	}
    setMoves (sys, doc, usedMovesEq, usedMovesPr);
    setPairPotentials (sys, doc);

    checkBounds (sys);
    sendMsg("System from "+filename+" passed bounds checks");

    setSystemBarriers (sys, doc);
    sendMsg("Initialized barriers from "+filename);

    sendMsg("Successfully read valid parameters from "+filename);
    return sys;
}

/*!
 * Assign the Monte Carlo moves based on the JSON input file.  Uses same information to specify "production" and "equilibration" phases.
 *
 * \param [in] sys Simulation system that has been initialized
 * \param [in] doc JSON document corresponding to input file
 * \params [in] usedMovesEq Pointer to move object that will be used during "equilibration" (WL + Crossover)
 * \params [in] usedMovesPr Pointer to move object that will be used during "production" (TMMC)
 */
void setMoves (simSystem &sys, const rapidjson::Document &doc, moves* usedMovesEq, moves* usedMovesPr) {
    std::vector < double > ref (sys.nSpecies(), 0);
	std::vector < std::vector < double > > probPrSwap (sys.nSpecies(), ref);
	std::vector < double > probPrInsDel (sys.nSpecies(), 0), probPrDisp (sys.nSpecies(), 0), maxPrD (sys.nSpecies(), 0);

    if (!doc.HasMember("moves")) throw customException("Input file does not have Monte Carlo moves specified");
    if (!doc["moves"].IsObject()) throw customException("Input file does not have Monte Carlo moves specified as correct JSON document");

    // Insert/Delete moves
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		std::string dummy = "ins_del_" + numToStr(i+1);
        if (!doc["moves"].HasMember(dummy.c_str())) throw customException("Input file does not have insert/delete move specified for species "+numToStr(i+1));
        if (!doc["moves"][dummy.c_str()].IsNumber()) throw customException("Input file does not correctly specify insert/delete move probability for species "+numToStr(i+1));
		probPrInsDel[i] = doc["moves"][dummy.c_str()].GetDouble();
	}

    // Translation moves
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
        std::string dummy = "displace_" + numToStr(i+1);
        if (!doc["moves"].HasMember(dummy.c_str())) throw customException("Input file does not have displacement move specified for species "+numToStr(i+1));
        if (!doc["moves"][dummy.c_str()].IsNumber()) throw customException("Input file does not correctly specify displacement move probability for species "+numToStr(i+1));
		probPrDisp[i] = doc["moves"][dummy.c_str()].GetDouble();

        dummy = "max_displacement_" + numToStr(i+1);
        if (!doc["moves"].HasMember(dummy.c_str())) throw customException("Input file does not have displacement magnitude specified for species "+numToStr(i+1));
        if (!doc["moves"][dummy.c_str()].IsNumber()) throw customException("Input file does not correctly specify displacement move magnitude for species "+numToStr(i+1));
		maxPrD[i] = doc["moves"][dummy.c_str()].GetDouble();
	}

    // Swap moves
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		for (unsigned int j = i+1; j < sys.nSpecies(); ++j) {
			std::string name1 = "swap_"+numToStr(i+1)+"_"+numToStr(j+1);
			std::string name2 = "swap_"+numToStr(j+1)+"_"+numToStr(i+1);
			std::string moveName = "";
			bool foundIJ = false;

			if (doc["moves"].HasMember(name1.c_str())) {
				moveName = name1;
				foundIJ = true;
			} else if (doc["moves"].HasMember(name2.c_str()) && !foundIJ) {
				moveName = name2;
				foundIJ = true;
			} else if (doc["moves"].HasMember(name2.c_str()) && foundIJ) {
                sendErr("Input file doubly specifies production swap move probability for species pair ("+numToStr(i+1)+", "+numToStr(j+1)+")");
				exit(SYS_FAILURE);
			} else {
                sendErr("Input file does not specify production swap move probability for species pair ("+numToStr(i+1)+", "+numToStr(j+1)+")");
				exit(SYS_FAILURE);
			}

            if (!doc["moves"][moveName.c_str()].IsNumber()) throw customException("Input file does not correctly specify swap move probability for species pair ("+numToStr(i+1)+", "+numToStr(j+1)+")");
			probPrSwap[i][j] = doc["moves"][moveName.c_str()].GetDouble();
			probPrSwap[j][i] = doc["moves"][moveName.c_str()].GetDouble();
		}
	}

    usedMovesEq->setM(sys.getTotalM());
    usedMovesPr->setM(sys.getTotalM());
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
        usedMovesEq->addInsert(i, probPrInsDel[i]);
        usedMovesPr->addInsert(i, probPrInsDel[i]);

        usedMovesEq->addDelete(i, probPrInsDel[i]);
        usedMovesPr->addDelete(i, probPrInsDel[i]);

        usedMovesEq->addTranslate(i, probPrDisp[i], maxPrD[i], sys.box());
        usedMovesPr->addTranslate(i, probPrDisp[i], maxPrD[i], sys.box());

        for (unsigned int j = i+1; j < sys.nSpecies(); ++j) {
            usedMovesEq->addSwap(i, j, probPrSwap[i][j]);
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
                sendErr("Input file doubly specifies pair potential for species pair ("+numToStr(i+1)+", "+numToStr(j+1)+")");
				exit(SYS_FAILURE);
			} else {
                sendErr("Input file does not specify pair potential for species pair ("+numToStr(i+1)+", "+numToStr(j+1)+")");
				exit(SYS_FAILURE);
			}

            if (!doc[ppotName.c_str()].IsString()) throw customException ("Pair potential is not a name for ("+numToStr(i+1)+","+numToStr(j+1)+")");
			ppotType[ppotTypeIndex] = doc[ppotName.c_str()].GetString();
            dummy = ppotName+"_params";
            if (!doc.HasMember(dummy.c_str())) throw customException ("Input file missing pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+")");
            if (!doc[dummy.c_str()].IsObject()) throw customException ("Pair potential's parameters are not valid json document");

            std::vector < double > params;

            bool useCellList = false; // default
            if (doc[dummy.c_str()].HasMember("cell_list")) {
                useCellList = doc[dummy.c_str()]["cell_list"].GetBool();
            }

            if (ppotType[ppotTypeIndex] == "square_well") {
                // Expects sigma, width, depth, cell_list
                if (!doc[dummy.c_str()].HasMember("sigma")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"sigma\"");
                if (!doc[dummy.c_str()].HasMember("width")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"width\"");
                if (!doc[dummy.c_str()].HasMember("depth")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"depth\"");

                if (!doc[dummy.c_str()]["sigma"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"sigma\" is not a number");
                if (!doc[dummy.c_str()]["width"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"width\" is not a number");
                if (!doc[dummy.c_str()]["depth"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"depth\" is not a number");

                params.push_back(doc[dummy.c_str()]["sigma"].GetDouble());
                params.push_back(doc[dummy.c_str()]["width"].GetDouble());
                params.push_back(doc[dummy.c_str()]["depth"].GetDouble());
            } else if (ppotType[ppotTypeIndex] == "lennard_jones") {
                // Expects epsilon, sigma, r_cut, u_shift
                if (!doc[dummy.c_str()].HasMember("epsilon")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"epsilon\"");
                if (!doc[dummy.c_str()].HasMember("sigma")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"sigma\"");
                if (!doc[dummy.c_str()].HasMember("r_cut")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"r_cut\"");
                if (!doc[dummy.c_str()].HasMember("u_shift")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"u_shift\"");

                if (!doc[dummy.c_str()]["epsilon"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"epsilon\" is not a number");
                if (!doc[dummy.c_str()]["sigma"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"sigma\" is not a number");
                if (!doc[dummy.c_str()]["r_cut"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"r_cut\" is not a number");
                if (!doc[dummy.c_str()]["u_shift"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"u_shift\" is not a number");

                params.push_back(doc[dummy.c_str()]["epsilon"].GetDouble());
                params.push_back(doc[dummy.c_str()]["sigma"].GetDouble());
                params.push_back(doc[dummy.c_str()]["r_cut"].GetDouble());
                params.push_back(doc[dummy.c_str()]["u_shift"].GetDouble());
            } else if (ppotType[ppotTypeIndex] == "fs_lennard_jones") {
                // Expects epsilon, sigma, r_cut
                if (!doc[dummy.c_str()].HasMember("epsilon")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"epsilon\"");
                if (!doc[dummy.c_str()].HasMember("sigma")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"sigma\"");
                if (!doc[dummy.c_str()].HasMember("r_cut")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"r_cut\"");

                if (!doc[dummy.c_str()]["epsilon"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"epsilon\" is not a number");
                if (!doc[dummy.c_str()]["sigma"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"sigma\" is not a number");
                if (!doc[dummy.c_str()]["r_cut"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"r_cut\" is not a number");

                params.push_back(doc[dummy.c_str()]["epsilon"].GetDouble());
                params.push_back(doc[dummy.c_str()]["sigma"].GetDouble());
                params.push_back(doc[dummy.c_str()]["r_cut"].GetDouble());
            } else if (ppotType[ppotTypeIndex] == "hard_sphere") {
                // Expects sigma
                if (!doc[dummy.c_str()].HasMember("sigma")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"sigma\"");

                if (!doc[dummy.c_str()]["sigma"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"sigma\" is not a number");

                params.push_back(doc[dummy.c_str()]["sigma"].GetDouble());
            } else if (ppotType[ppotTypeIndex] == "tabulated") {
                // Expects r_cut, r_shift, u_shift, u_infinity
                if (!doc[dummy.c_str()].HasMember("r_cut")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"r_cut\"");
                if (!doc[dummy.c_str()].HasMember("r_shift")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"r_shift\"");
                if (!doc[dummy.c_str()].HasMember("u_shift")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"u_shift\"");
                if (!doc[dummy.c_str()].HasMember("u_infinity")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"u_infinity\"");

                if (!doc[dummy.c_str()]["r_cut"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"r_cut\" is not a number");
                if (!doc[dummy.c_str()]["r_shift"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"r_shift\" is not a number");
                if (!doc[dummy.c_str()]["u_shift"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"u_shift\" is not a number");
                if (!doc[dummy.c_str()]["u_infinity"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"u_infinity\" is not a number");

                params.push_back(doc[dummy.c_str()]["r_cut"].GetDouble());
                params.push_back(doc[dummy.c_str()]["r_shift"].GetDouble());
                params.push_back(doc[dummy.c_str()]["u_shift"].GetDouble());
                params.push_back(doc[dummy.c_str()]["u_infinity"].GetDouble());
            } else {
                throw customException ("Unrecognized pair potential "+ppotType[ppotTypeIndex]);
            }

            params.push_back(Mtot);

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
void setConfig (simSystem &sys, const std::string filename) {
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
            sendErr(ce.what());
		}
	} else if (restart_file == "" && sys.totNMin() > 0) {
        sendMsg("Automatically generating the initial configuration");

		// have to generate initial configuration manually - start with mu = INF
        std::vector < double > initMu (doc["num_species"].GetInt(), 1.0e2);

		simSystem initSys (doc["num_species"].GetInt(), doc["beta"].GetDouble()/100.0, sysBox, initMu, sysMax, sysMin, Mtot, duh, max_order); // beta =  1/T, so low beta to have high T
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
                    sendErr("Order of initialization goes out of bounds, should include 0 <= i < nSpec");
    				exit(SYS_FAILURE);
    			}
    			if (check_init[initialization_order[i]] != 0) {
                    sendErr("Order of initialization repeats itself");
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
                    sendErr("Initialization fraction out of bounds");
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
            sendMsg("Initializing species "+numToStr(i)+" configurations");

			// insert this species i
			moves initMove (initSys.getTotalM());
            initMove.addInsert(i, 1.0);

			// also add displacment moves for all species present
			for (unsigned int j = 0; j <= idx; ++j) {
                sendMsg("Added translation moves for initialization of species "+numToStr(initialization_order[j]));
                initMove.addTranslate(initialization_order[j], 2.0, 1.0, initSys.box());
			}

			// now do simuation until within proper range
			int targetNum = sys.totNMin()*init_frac[idx];
			if (idx == sys.nSpecies() - 1) {
				// to account for integer rounding
				targetNum = sys.totNMin() - added;
			}
			added += targetNum;

            sendMsg("Target number = "+numToStr(targetNum)+" for species "+numToStr(i+1));
			int tmpCounter = 0, statusPrint = 10e6;
			while (initSys.numSpecies[i] < targetNum) {
				try {
                    initMove.makeMove(initSys);
                } catch (customException &ce) {
                    std::string msg = ce.what();
                    sendErr("Failed to create an initial configuration because "+msg);
                    exit(SYS_FAILURE);
                }
				tmpCounter++;
				if (tmpCounter%statusPrint == 0) {
					tmpCounter = 0;
                    sendMsg("Grew "+numToStr(initSys.numSpecies[i])+" atoms of type "+numToStr(i)+" so far");
				}
			}
		}

		// print snapshot from Reading initial configuration
		initSys.printSnapshot("auto-init.xyz", "auto-generated initial configuration");

		// read into sys
		try {
			sys.readConfig("auto-init.xyz");
		} catch (customException &ce) {
            std::string msg = ce.what();
            sendErr("Failed to read auto-generated initialization file because "+msg);
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
                sendErr(ce.what());
				exit(SYS_FAILURE);
			}
			convention0 = true;
		}
		for (unsigned int j = 1; j <= MAX_BARRIERS_PER_SPECIES; ++j) {
			// alternatively allow multiple walls to specified with a suffix up to a max
			std::string dummy = "hardWallZ_" + std::to_string(i+1) + "_" + std::to_string(j);
            if (doc.HasMember(dummy.c_str())) {
				if (convention0) {
                    sendErr("Error, multiple barrier naming conventions used for the same species");
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
						sendErr(ce.what());
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
                sendErr(ce.what());
                exit(SYS_FAILURE);
            }
            convention0 = true;
		}
		for (unsigned int j = 1; j <= MAX_BARRIERS_PER_SPECIES; ++j) {
			// alternatively allow multiple walls to specified with a suffix up to a max
			std::string dummy = "squareWellWallZ_" + std::to_string(i+1) + "_" + std::to_string(j);
			if (doc.HasMember(dummy.c_str())) {
				if (convention0) {
                    sendErr("Error, multiple barrier naming conventions used for the same species");
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
						sendErr(ce.what());
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
                sendErr(ce.what());
                exit(SYS_FAILURE);
            }
            convention0 = true;
		}
		for (unsigned int j = 1; j <= MAX_BARRIERS_PER_SPECIES; ++j) {
			// alternatively allow multiple walls to specified with a suffix up to a max
			std::string dummy = "cylinderZ_" + std::to_string(i+1) + "_" + std::to_string(j);
			if (doc.HasMember(dummy.c_str())) {
				if (convention0) {
                    sendErr("Error, multiple barrier naming conventions used for the same species");
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
						sendErr(ce.what());
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
                sendErr(ce.what());
                exit(SYS_FAILURE);
            }
            convention0 = true;
        }
        for (unsigned int j = 1; j <= MAX_BARRIERS_PER_SPECIES; ++j) {
            // alternatively allow multiple walls to specified with a suffix up to a max
            std::string dummy = "rightTriangleXZ_" + std::to_string(i+1) + "_" + std::to_string(j);
            if (doc.HasMember(dummy.c_str())) {
                if (convention0) {
                    sendErr("Error, multiple barrier naming conventions used for the same species");
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
                        sendErr(ce.what());
                        exit(SYS_FAILURE);
                    }
                }
            }
        }
    }
}
