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
    if (!doc.HasMember("num_species")) throw customException("\"num_species\" is not specified in "+filename);
    if (!doc["num_species"].IsInt()) throw customException("\"num_species\" is not an integer in "+filename);

    if (!doc.HasMember("beta")) throw customException("\"beta\" is not specified in "+filename);
    if (!doc["beta"].IsNumber()) throw customException("\"beta\" is not a number in "+filename);

    if (!doc.HasMember("box")) throw customException("\"box\" is not specified in "+filename);
    if (!doc["box"].IsArray()) throw customException("\"box\" is not an array in "+filename);
    if (doc["box"].Size() != 3) throw customException("\"box\" is not a length 3 array in "+filename);
	std::vector < double > sysBox (3, 0);
	for (rapidjson::SizeType i = 0; i < doc["box"].Size(); ++i) {
        if (!doc["box"][i].IsNumber()) throw customException("box index "+numToStr(i)+" is not a number in "+filename);
		sysBox[i] = doc["box"][i].GetDouble();
	}

	double duh = 10.0;
	if (doc.HasMember("delta_u_hist")) {
        if (!doc["delta_u_hist"].IsNumber()) throw customException("\"delta_u_hist\" is not a number in "+filename);
		duh = doc["delta_u_hist"].GetDouble();
	}

	int maxOrder = 2;
	if (doc.HasMember("max_order")) {
        if (!doc["max_order"].IsInt()) throw customException("\"max_order\" is not an integer in "+filename);
		maxOrder = doc["max_order"].GetInt();
	}

	bool useKe = false;
	if (doc.HasMember("use_ke")) {
        if (!doc["use_ke"].IsBool()) throw customException("\"use_ke\" is not a boolean in "+filename);
		useKe = doc["use_ke"].GetBool();
	}

    if (!doc.HasMember("mu")) throw customException("\"mu\" is not specified in "+filename);
    if (!doc["mu"].IsArray()) throw customException("\"mu\" is not an array in "+filename);
    if (doc["mu"].Size() != doc["num_species"].GetInt()) throw customException("\"mu\" is not specified for each species in "+filename);
	std::vector < double > sysMu (doc["mu"].Size(), 0);
	for (rapidjson::SizeType i = 0; i < doc["mu"].Size(); ++i) {
        if (!doc["mu"][i].IsNumber()) throw customException("\"mu\" for species "+numToStr(i+1)+" is not a number in "+filename);
		sysMu[i] = doc["mu"][i].GetDouble();
	}

    if (!doc.HasMember("seed")) throw customException("\"seed\" is not specified in "+filename);
    if (!doc["seed"].IsInt()) throw customException("\"seed\" is not an integer in "+filename);
	RNG_SEED = doc["seed"].GetInt();

    if (!doc.HasMember("max_N")) throw customException("\"max_N\" is not specified in "+filename);
    if (!doc["max_N"].IsArray()) throw customException("\"max_N\" is not an array in "+filename);
	if (doc["max_N"].Size() != doc["num_species"].GetInt()) throw customException("\"max_N\" is not specified for each species in "+filename);
	std::vector < int > sysMax (doc["max_N"].Size(), 0);
	for (rapidjson::SizeType i = 0; i < doc["max_N"].Size(); ++i) {
        if (!doc["max_N"][i].IsInt()) throw customException("\"max_N\" of species "+numToStr(i+1)+" is not an integer in "+filename);
		sysMax[i] = doc["max_N"][i].GetInt();
	}

    if (!doc.HasMember("min_N")) throw customException("\"min_N\" is not specified in "+filename);
    if (!doc["min_N"].IsArray()) throw customException("\"min_N\" is not an array in "+filename);
	if (doc["min_N"].Size() != doc["num_species"].GetInt()) throw customException("\"min_N\" is not specified for each species in "+filename);
	std::vector < int > sysMin (doc["min_N"].Size(), 0);
	for (rapidjson::SizeType i = 0; i < doc["min_N"].Size(); ++i) {
        if (!doc["min_N"][i].IsInt()) throw customException("\"min_N\" of species "+numToStr(i+1)+" is not an integer in "+filename);
		sysMin[i] = doc["min_N"][i].GetInt();
	}

	int Mtot = 1;
	if (doc.HasMember("num_expanded_states")) {
        if (!doc["num_expanded_states"].IsInt()) throw customException("\"num_expanded_states\" is not an integer in "+filename);
		Mtot = doc["num_expanded_states"].GetInt();
	}

	simSystem sys (doc["num_species"].GetInt(), doc["beta"].GetDouble(), sysBox, sysMu, sysMax, sysMin, Mtot, duh, maxOrder);
	if (useKe) {
		sys.toggleKE();
		if (sys.addKECorrection() == false) {
			throw customException ("Unable to set KE flag");
		}
	}

	std::vector < int > sysWindow;
	if (doc.HasMember("window")) {
        if (!doc["window"].IsArray()) throw customException("\"window\" is not an array in "+filename);
        if (doc["window"].Size() != 2) throw customException("\"window\" should have 2 entries (min,max) in "+filename);
		sysWindow.resize(2, 0);
        if (!doc["window"][0].IsInt()) throw customException("\"window\" min is not an integer in "+filename);
        if (!doc["window"][1].IsInt()) throw customException("\"window\" max is not an integer in "+filename);
		sysWindow[0] = doc["window"][0].GetInt();
		sysWindow[1] = doc["window"][1].GetInt();
	}

	if (sysWindow.begin() != sysWindow.end()) {
		sys.setTotNBounds(sysWindow);
	}

    if (!doc.HasMember("tmmc_sweep_size")) throw customException("\"tmmc_sweep_size\" is not specified in "+filename);
    if (!doc["tmmc_sweep_size"].IsNumber()) throw customException("\"tmmc_sweep_size\" is not a number in "+filename);
	double tmpT = doc["tmmc_sweep_size"].GetDouble(); // Possibly in scientific notation
	sys.tmmcSweepSize = tmpT; // Convert

    if (!doc.HasMember("total_tmmc_sweeps")) throw customException("\"total_tmmc_sweeps\" is not specified in "+filename);
    if (!doc["total_tmmc_sweeps"].IsNumber()) throw customException("\"total_tmmc_sweeps\" is not a number in "+filename);
	double tmpS = doc["total_tmmc_sweeps"].GetDouble(); // Possibly in scientific notation
	sys.totalTMMCSweeps = tmpS; // Convert

    if (!doc.HasMember("wala_sweep_size")) throw customException("\"wala_sweep_size\" is not specified in "+filename);
    if (!doc["wala_sweep_size"].IsNumber()) throw customException("\"wala_sweep_size\" is not a number in "+filename);
	double tmpW = doc["wala_sweep_size"].GetDouble(); // Possibly in scientific notation
	sys.wlSweepSize = tmpW; // Convert

    if (!doc.HasMember("wala_g")) throw customException("\"wala_g\" is not specified in "+filename);
    if (!doc["wala_g"].IsNumber()) throw customException("\"wala_g\" is not a number in "+filename);
	sys.wala_g = doc["wala_g"].GetDouble();

    if (!doc.HasMember("wala_s")) throw customException("\"wala_s\" is not specified in "+filename);
    if (!doc["wala_s"].IsNumber()) throw customException("\"wala_s\" is not a number in "+filename);
	sys.wala_s = doc["wala_s"].GetDouble();

	if (doc.HasMember("lnF_start")) {
        if (!doc["lnF_start"].IsNumber()) throw customException("\"lnF_start\" is not a number in "+filename);
		sys.lnF_start = doc["lnF_start"].GetDouble(); // Bounds are checked later
	}

	if (doc.HasMember("lnF_end")) {
        if (!doc["lnF_end"].IsNumber()) throw customException("\"lnF_end\" is not a number in "+filename);
		sys.lnF_end = doc["lnF_end"].GetDouble();
		if (sys.lnF_end >= 1.0) {
            sendErr("Terminal lnF factor for Wang-Landau must be < 1");
			exit(SYS_FAILURE);
		}
        if (sys.lnF_end <= 0.0) {
            sendErr("Terminal lnF factor for Wang-Landau must be a positive number");
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
        if (!doc["restart_from_wala_lnPI"].IsString()) throw customException("\"restart_from_wala_lnPI\" filename is not a string in "+filename);
		sys.restartFromWALAFile = doc["restart_from_wala_lnPI"].GetString();
		if (sys.restartFromWALAFile != "") {
			sys.restartFromWALA = true;
		}
	}

	// Restarting from TMMC overrides WL by skipping that portion altogether
	sys.restartFromTMMC = false;
	sys.restartFromTMMCFile = "";
	if (doc.HasMember("restart_from_tmmc_C")) {
        if (!doc["restart_from_tmmc_C"].IsString()) throw customException("\"restart_from_tmmc_C\" filename is not a string in "+filename);
		sys.restartFromTMMCFile = doc["restart_from_tmmc_C"].GetString();
		if (sys.restartFromTMMCFile != "") {
			sys.restartFromTMMC = true;
		}
	}

	// Number of times the TMMC C matrix has to be traversed during the WALA --> TMMC crossover
	if (doc.HasMember("num_crossover_visits")) {
        if (!doc["num_crossover_visits"].IsNumber()) throw customException("\"num_crossover_visits\" is not a number in "+filename);
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

    sendMsg("Successfully read valid parameters from "+filename);
    return sys;
}

/*!
 * Assign the Monte Carlo moves based on the JSON input file.  Uses same information to specify "production" and "equilibration" phases.
 * Clears any existing information and will overwrite with information from doc.
 *
 * \param [in] sys Simulation system that has been initialized
 * \param [in] doc JSON document corresponding to input file
 * \params [in] usedMovesEq Pointer to move object that will be used during "equilibration" (WL + Crossover)
 * \params [in] usedMovesPr Pointer to move object that will be used during "production" (TMMC)
 */
void setMoves (simSystem &sys, const rapidjson::Document &doc, moves* usedMovesEq, moves* usedMovesPr) {
    usedMovesEq->clearAll();
    usedMovesPr->clearAll();

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
 * Overwrites any existing pair potential information with new settings from doc.
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
			std::string ppotName = "", dummy = "", tabFile = "";
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
                // Expects sigma, width, epsilon, cell_list
                if (!doc[dummy.c_str()].HasMember("sigma")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"sigma\"");
                if (!doc[dummy.c_str()].HasMember("width")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"width\"");
                if (!doc[dummy.c_str()].HasMember("epsilon")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"epsilon\"");

                if (!doc[dummy.c_str()]["sigma"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"sigma\" is not a number");
                if (!doc[dummy.c_str()]["width"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"width\" is not a number");
                if (!doc[dummy.c_str()]["epsilon"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"epsilon\" is not a number");

                params.push_back(doc[dummy.c_str()]["sigma"].GetDouble());
                params.push_back(doc[dummy.c_str()]["width"].GetDouble());
                params.push_back(doc[dummy.c_str()]["epsilon"].GetDouble());
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
                // Also must specify file to load potential from
                if (!doc[dummy.c_str()].HasMember("r_cut")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"r_cut\"");
                if (!doc[dummy.c_str()].HasMember("r_shift")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"r_shift\"");
                if (!doc[dummy.c_str()].HasMember("u_shift")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"u_shift\"");
                if (!doc[dummy.c_str()].HasMember("u_infinity")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"u_infinity\"");
                if (!doc[dummy.c_str()].HasMember("filename")) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") is missing \"filename\"");

                if (!doc[dummy.c_str()]["r_cut"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"r_cut\" is not a number");
                if (!doc[dummy.c_str()]["r_shift"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"r_shift\" is not a number");
                if (!doc[dummy.c_str()]["u_shift"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"u_shift\" is not a number");
                if (!doc[dummy.c_str()]["u_infinity"].IsNumber()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"u_infinity\" is not a number");
                if (!doc[dummy.c_str()]["filename"].IsString()) throw customException ("Pair potential parameters for ("+numToStr(i+1)+","+numToStr(j+1)+") parameter \"filename\" is not a string");

                params.push_back(doc[dummy.c_str()]["r_cut"].GetDouble());
                params.push_back(doc[dummy.c_str()]["r_shift"].GetDouble());
                params.push_back(doc[dummy.c_str()]["u_shift"].GetDouble());
                params.push_back(doc[dummy.c_str()]["u_infinity"].GetDouble());
                tabFile = doc[dummy.c_str()]["filename"].GetString();
            } else {
                throw customException ("Unrecognized pair potential "+ppotType[ppotTypeIndex]);
            }

            params.push_back(Mtot);
            sys.addPotential(i, j, ppotType[ppotTypeIndex], params, useCellList, tabFile);
            sys.ppot[i][j]->savePotential(ppotName+".dat", 0.01, 0.01);

			ppotTypeIndex++;
		}
	}
}

/*!
 * Setup a system's initial configuration as necessary.
 * Will empty a system if there are currently any particles present and overwrite with new information.
 * 1. If "restart_file" in input json file, read initial config from there.
 * 2. In not, randomly generate initial configuration.
 *
 * \param [in] sys System to initialize
 * \param [in] filename Input JSON filename
 */
void setConfig (simSystem &sys, const std::string filename) {
    rapidjson::Document doc;
    parseJson (filename, doc);

    // Get a few things from file not easily accessible from system object
    std::string restart_file = "";
	if (doc.HasMember("restart_file")) {
        restart_file = doc["restart_file"].GetString();
    }

    std::vector < int > sysMax (doc["max_N"].Size(), 0);
	for (rapidjson::SizeType i = 0; i < doc["max_N"].Size(); ++i) {
		sysMax[i] = doc["max_N"][i].GetInt();
	}
	std::vector < int > sysMin (doc["min_N"].Size(), 0);
	for (rapidjson::SizeType i = 0; i < doc["min_N"].Size(); ++i) {
		sysMin[i] = doc["min_N"][i].GetInt();
	}

    // Rest from existing system
    int Mtot = sys.getTotalM();
    int maxOrder = sys.getMaxOrder();
    bool useKe = sys.addKECorrection();
    double duh = 10.0;
    std::vector < double > sysBox = sys.box();

    // Read from restart file if specified
	if (restart_file != "") {
		try {
			sys.readConfig(restart_file);
		} catch (customException &ce) {
            sendErr(ce.what());
		}
	} else if (restart_file == "" && sys.totNMin() > 0) {
        sendMsg("Automatically generating the initial configuration");

		// Have to generate initial configuration manually - start with mu = INF
        std::vector < double > initMu (doc["num_species"].GetInt(), 1.0e2);

		simSystem initSys (doc["num_species"].GetInt(), doc["beta"].GetDouble()/100.0, sysBox, initMu, sysMax, sysMin, Mtot, duh, maxOrder); // beta =  1/T, so low beta to have high T
		if (useKe) {
			initSys.toggleKE();
			if (initSys.addKECorrection() == false) {
				throw customException ("Unable to set KE flag");
			}
		}

        // Add the same potentials
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
            if (!doc["init_order"].IsArray()) throw customException("\"init_order\" is not an array in "+filename);
            if (doc["init_order"].Size() != doc["num_species"].GetInt()) throw customException("\"init_order\" not specified for each species in "+filename);

    		for (rapidjson::SizeType i = 0; i < doc["init_order"].Size(); ++i) {
                if (!doc["init_order"][i].IsInt()) throw customException("\"init_order\" is not an integer for species "+numToStr(i+1)+" in "+filename);
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
            if (!doc["init_frac"].IsArray()) throw customException("\"init_frac\" is not an array in "+filename);
            if (doc["init_frac"].Size() != doc["num_species"].GetInt()) throw customException("\"init_frac\" not specified for each species in "+filename);

    		sum = 0.0;
    		for (rapidjson::SizeType i = 0; i < doc["init_frac"].Size(); ++i) {
                if (!doc["init_frac"][i].IsNumber()) throw customException("\"init_frac\" is not a number for species "+numToStr(i+1)+" in "+filename);
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

		// Iteratively add each individual species, assume we want an equimolar mixture to start from
		int added = 0;
		for (unsigned int idx = 0; idx < sys.nSpecies(); ++idx) {
			unsigned int i = initialization_order[idx];
            sendMsg("Initializing species "+numToStr(i)+" configurations");

			// Insert this species i
			moves initMove (initSys.getTotalM());
            initMove.addInsert(i, 1.0);

			// Also add displacment moves for all species present
			for (unsigned int j = 0; j <= idx; ++j) {
                sendMsg("Added translation moves for initialization of species "+numToStr(initialization_order[j]));
                initMove.addTranslate(initialization_order[j], 2.0, 1.0, initSys.box());
			}

			// Now do simuation until within proper range
			int targetNum = sys.totNMin()*init_frac[idx];
			if (idx == sys.nSpecies() - 1) {
				// To account for integer rounding
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

		// Print snapshot from Reading initial configuration
		initSys.printSnapshot("auto-init.xyz", "auto-generated initial configuration");

		// Read into sys
		try {
			sys.readConfig("auto-init.xyz");
		} catch (customException &ce) {
            std::string msg = ce.what();
            sendErr("Failed to read auto-generated initialization file because "+msg);
        }
   	}
}

/*!
 * Initialize the barriers in a system by parsing the input document.
 * Clears any existing information and will overwrite with information from doc.
 *
 * \params [in, out] sys System to initialize with barriers
 * \params [in] doc Input JSON document
 */
void setSystemBarriers (simSystem &sys, const rapidjson::Document &doc) {
	int Mtot = sys.getTotalM();

    if (doc.HasMember("barriers")) {
        // Clear any existing barriers
        for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
            sys.speciesBarriers[i].clearAll();
        }

        // Iterate over all barriers specified for this species
        for (rapidjson::Value::ConstMemberIterator itr = doc["barriers"].MemberBegin(); itr != doc["barriers"].MemberEnd(); ++itr) {
            // Get barrier type and name
            std::string barrName = itr->name.GetString();
            if (!itr->value.IsObject()) throw customException ("Barrier "+barrName+" is not in a valid json document");
            if (!itr->value.HasMember("type")) throw customException ("Barrier "+barrName+" does not specify a type");
            if (!itr->value["type"].IsString()) throw customException ("Barrier "+barrName+" type is not a string");
            std::string barrType = itr->value["type"].GetString();

            // Get the species this barrier interacts with
            if (!itr->value.HasMember("species")) throw customException ("Barrier "+barrName+" does not specify a species to interact with");
            if (!itr->value["species"].IsInt()) throw customException ("Barrier "+barrName+" species is not an integer");
            const int species = itr->value["species"].GetInt();
            if (species < 1 || species > sys.nSpecies()) throw customException ("Barrier "+barrName+" species is not valid for this system");

            // Depending on barrier type, read parameters and initialize
            if (barrType == "hard_wall_z") {
                // Expects lb, ub, sigma
                if (!itr->value.HasMember("lb")) throw customException (barrName+" does not contain \"lb\" parameter");
                if (!itr->value.HasMember("ub")) throw customException (barrName+" does not contain \"ub\" parameter");
                if (!itr->value.HasMember("sigma")) throw customException (barrName+" does not contain \"sigma\" parameter");

                if (!itr->value["lb"].IsNumber()) throw customException ("\"lb\" for "+barrName+" is not a number");
                if (!itr->value["ub"].IsNumber()) throw customException ("\"ub\" for "+barrName+" is not a number");
                if (!itr->value["sigma"].IsNumber()) throw customException ("\"sigma\" for "+barrName+" is not a number");

                const double lbBarr = itr->value["lb"].GetDouble();
                const double ubBarr = itr->value["ub"].GetDouble();
                const double sigmaBarr = itr->value["sigma"].GetDouble();

    			try {
    				sys.speciesBarriers[species-1].addHardWallZ (lbBarr, ubBarr, sigmaBarr, Mtot);
    			} catch (customException &ce) {
                    sendErr(ce.what());
    				exit(SYS_FAILURE);
    			}
            } else if (barrType == "square_well_wall_z") {
                // Expect lb, ub, sigma, range, epsilon
                if (!itr->value.HasMember("lb")) throw customException (barrName+" does not contain \"lb\" parameter");
                if (!itr->value.HasMember("ub")) throw customException (barrName+" does not contain \"ub\" parameter");
                if (!itr->value.HasMember("sigma")) throw customException (barrName+" does not contain \"sigma\" parameter");
                if (!itr->value.HasMember("range")) throw customException (barrName+" does not contain \"range\" parameter");
                if (!itr->value.HasMember("epsilon")) throw customException (barrName+" does not contain \"epsilon\" parameter");

                if (!itr->value["lb"].IsNumber()) throw customException ("\"lb\" for "+barrName+" is not a number");
                if (!itr->value["ub"].IsNumber()) throw customException ("\"ub\" for "+barrName+" is not a number");
                if (!itr->value["sigma"].IsNumber()) throw customException ("\"sigma\" for "+barrName+" is not a number");
                if (!itr->value["range"].IsNumber()) throw customException ("\"range\" for "+barrName+" is not a number");
                if (!itr->value["epsilon"].IsNumber()) throw customException ("\"epsilon\" for "+barrName+" is not a number");

                const double lbBarr = itr->value["lb"].GetDouble();
                const double ubBarr = itr->value["ub"].GetDouble();
                const double sigmaBarr = itr->value["sigma"].GetDouble();
                const double rangeBarr = itr->value["range"].GetDouble();
                const double epsBarr = itr->value["epsilon"].GetDouble();

    			try {
    				sys.speciesBarriers[species-1].addSquareWellWallZ (lbBarr, ubBarr, sigmaBarr, rangeBarr, epsBarr, Mtot);
    			} catch (customException &ce) {
                    sendErr(ce.what());
    				exit(SYS_FAILURE);
    			}
            } else if (barrType == "cylinder_z") {
                // Expect x, y, radius, width, sigma, epsilon
                if (!itr->value.HasMember("x")) throw customException (barrName+" does not contain \"x\" parameter");
                if (!itr->value.HasMember("y")) throw customException (barrName+" does not contain \"y\" parameter");
                if (!itr->value.HasMember("radius")) throw customException (barrName+" does not contain \"radius\" parameter");
                if (!itr->value.HasMember("width")) throw customException (barrName+" does not contain \"width\" parameter");
                if (!itr->value.HasMember("sigma")) throw customException (barrName+" does not contain \"sigma\" parameter");
                if (!itr->value.HasMember("epsilon")) throw customException (barrName+" does not contain \"epsilon\" parameter");

                if (!itr->value["x"].IsNumber()) throw customException ("\"x\" for "+barrName+" is not a number");
                if (!itr->value["y"].IsNumber()) throw customException ("\"y\" for "+barrName+" is not a number");
                if (!itr->value["radius"].IsNumber()) throw customException ("\"radius\" for "+barrName+" is not a number");
                if (!itr->value["width"].IsNumber()) throw customException ("\"width\" for "+barrName+" is not a number");
                if (!itr->value["sigma"].IsNumber()) throw customException ("\"sigma\" for "+barrName+" is not a number");
                if (!itr->value["epsilon"].IsNumber()) throw customException ("\"epsilon\" for "+barrName+" is not a number");

                const double xBarr = itr->value["x"].GetDouble();
                const double yBarr = itr->value["y"].GetDouble();
                const double radiusBarr = itr->value["radius"].GetDouble();
                const double widthBarr = itr->value["width"].GetDouble();
                const double sigmaBarr = itr->value["sigma"].GetDouble();
                const double epsBarr = itr->value["epsilon"].GetDouble();

    			try {
    				sys.speciesBarriers[species-1].addCylinderZ (xBarr, yBarr, radiusBarr, widthBarr, sigmaBarr, epsBarr, Mtot);
    			} catch (customException &ce) {
                    sendErr(ce.what());
    				exit(SYS_FAILURE);
    			}
            } else if (barrType == "right_triangle_xz") {
                // Expect parameters width, theta, lamW, epsilon, sigma, sep, offset, zbase, top
                if (!itr->value.HasMember("width")) throw customException (barrName+" does not contain \"width\" parameter");
                if (!itr->value.HasMember("theta")) throw customException (barrName+" does not contain \"theta\" parameter");
                if (!itr->value.HasMember("lamW")) throw customException (barrName+" does not contain \"lamW\" parameter");
                if (!itr->value.HasMember("epsilon")) throw customException (barrName+" does not contain \"epsilon\" parameter");
                if (!itr->value.HasMember("sigma")) throw customException (barrName+" does not contain \"sigma\" parameter");
                if (!itr->value.HasMember("sep")) throw customException (barrName+" does not contain \"sep\" parameter");
                if (!itr->value.HasMember("offset")) throw customException (barrName+" does not contain \"offset\" parameter");
                if (!itr->value.HasMember("zbase")) throw customException (barrName+" does not contain \"zbase\" parameter");
                if (!itr->value.HasMember("top")) throw customException (barrName+" does not contain \"top\" parameter");

                if (!itr->value["width"].IsNumber()) throw customException ("\"width\" for "+barrName+" is not a number");
                if (!itr->value["theta"].IsNumber()) throw customException ("\"theta\" for "+barrName+" is not a number");
                if (!itr->value["lamW"].IsNumber()) throw customException ("\"lamW\" for "+barrName+" is not a number");
                if (!itr->value["epsilon"].IsNumber()) throw customException ("\"epsilon\" for "+barrName+" is not a number");
                if (!itr->value["sigma"].IsNumber()) throw customException ("\"sigma\" for "+barrName+" is not a number");
                if (!itr->value["sep"].IsNumber()) throw customException ("\"sep\" for "+barrName+" is not a number");
                if (!itr->value["offset"].IsNumber()) throw customException ("\"offset\" for "+barrName+" is not a number");
                if (!itr->value["zbase"].IsNumber()) throw customException ("\"zbase\" for "+barrName+" is not a number");
                if (!itr->value["top"].IsBool()) throw customException ("\"top\" for "+barrName+" is not a boolean");

                const double widthBarr = itr->value["width"].GetDouble();
                const double thetaBarr = itr->value["theta"].GetDouble();
                const double lamwBarr = itr->value["lamW"].GetDouble();
                const double epsBarr = itr->value["epsilon"].GetDouble();
                const double sigmaBarr = itr->value["sigma"].GetDouble();
                const double sepBarr = itr->value["sep"].GetDouble();
                const double offsetBarr = itr->value["offset"].GetDouble();
                const double zbaseBarr = itr->value["zbase"].GetDouble();
                const double topBarr = itr->value["top"].GetBool();

    			try {
    				sys.speciesBarriers[species-1].addRightTriangleXZ (widthBarr, thetaBarr, lamwBarr, epsBarr, sigmaBarr, sepBarr, offsetBarr, sys.box(), zbaseBarr, topBarr, Mtot);
    			} catch (customException &ce) {
                    sendErr(ce.what());
    				exit(SYS_FAILURE);
    			}
            } else {
                throw customException ("Unrecognized barrier type "+barrType+" from barrier "+barrName);
            }
        }
        sendMsg("Initialized barriers");
    } else {
        sendMsg("No barriers to initialize");
    }
}
