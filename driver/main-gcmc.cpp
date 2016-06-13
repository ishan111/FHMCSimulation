/*!
 * Perform (un/)biased multicomponent GCMC.
 * 
 * \author Nathan A. Mahynski
 * \date 02/24/16
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
#include "../src/system.h"
#include "../src/utilities.h"
#include "../src/global.h"
#include "../src/insert.h"
#include "../src/delete.h"
#include "../src/translate.h"
#include "../src/swap.h"
#include "../src/moves.h"
#include "../src/barrier.h"

// JSON interface from local distro of rapidjson
#include "../src/rapidjson/include/rapidjson/document.h"
#include "../src/rapidjson/include/rapidjson/writer.h"
#include "../src/rapidjson/include/rapidjson/stringbuffer.h"
#include "../src/rapidjson/include/rapidjson/filereadstream.h"
#include "../src/rapidjson/include/rapidjson/prettywriter.h"

/*!
 * Initialize the barriers in a system by parsing the input document.  This function is defined separately since it must be done several times.
 *
 * \params [in, out] sys System to initialize with barriers
 * \params [in] doc Input JSON document
 */
void initializeSystemBarriers (simSystem &sys, const rapidjson::Document &doc) {
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
        	std::string dummy = "hardWallZ_" + sstr(i+1);
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
                std::string dummy = "hardWallZ_" + sstr(i+1) + "_" + sstr(j);
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
        	std::string dummy = "squareWellWallZ_" + sstr(i+1);
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
                std::string dummy = "squareWellWallZ_" + sstr(i+1) + "_" + sstr(j);
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
    
    // rightTriangleXZ (expect parameters: {width, theta, lamW, eps, sigma, sep, offset, zbase, top})
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
        bool convention0 = false;
        std::string dummy = "rightTriangleXZ_" + sstr(i+1);
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
            std::string dummy = "rightTriangleXZ_" + sstr(i+1) + "_" + sstr(j);
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

    	// in the future, can be multiple instances of the same barrier, but for the above, assume only 1
}

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

	assert(doc.HasMember("eq_steps"));
	assert(doc["eq_steps"].IsNumber());
	double eqs = doc["eq_steps"].GetDouble(); // possibly in scientific notation
	const long long int eqSteps = eqs; // convert

	assert(doc.HasMember("pr_steps"));
	assert(doc["pr_steps"].IsNumber());
	double prs = doc["pr_steps"].GetDouble(); // possibly in scientific notation
	const long long int prSteps = prs; // convert

	assert(doc.HasMember("snap_freq"));
	assert(doc["snap_freq"].IsNumber());
	double sfreq = doc["snap_freq"].GetDouble(); // possibly in scientific notation
	const long long int movePrint = sfreq; // convert
	
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
    
    	// Add barriers for each species
	initializeSystemBarriers (sys, doc);

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
        	std::vector < double > initMu (doc["num_species"].GetInt(), 1.0e2);
        	simSystem initSys (doc["num_species"].GetInt(), 1/10., sysBox, initMu, sysMax, sysMin, Mtot); // beta =  1/T, so low beta to have high T

        	// add the same potentials
        	int initInd = 0;
        	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
            		for (unsigned int j = i; j < sys.nSpecies(); ++j) {
                		initSys.addPotential (i, j, ppotArray[initInd], true); // default use of cell list even if otherwise not in main simulation
                		initInd++;
            		}
        	}

		initializeSystemBarriers (initSys, doc);

		// iteratively add each individual species, assume we want an equimolar mixture to start from
		//for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		for (unsigned int i = 1; i < 2; ++i) {
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
			int targetNum = sys.totNMin(); //sys.nSpecies();
			/*if (i == sys.nSpecies() - 1) {
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
        
	for (unsigned long long int move = 0; move < eqSteps; ++move) {
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

        // Report move statistics for equilibration stage
        char eq_statName [80];
        strftime (eq_statName,80,"%Y_%m_%d_%H_%M_%S-equilibration-stats.log",timeinfo);
        std::ofstream eq_statFile (eq_statName);
        std::vector < std::vector < double > > eq_stats = usedMovesEq.reportMoveStatistics();
        eq_statFile << " ---------- Move Statistics --------- " << std::endl << " Move\t\% Success" << std::endl;
        for (unsigned int i = 0; i < eq_stats.size(); ++i) {
            double prod = 1.0;
            for (unsigned int j = 0; j < eq_stats[i].size(); ++j) {
                prod *= eq_stats[i][j];
                eq_statFile << usedMovesEq.includedMoves()[i]->myName() << " (from M = " << j << ")\t" << eq_stats[i][j]*100.0 << std::endl;
            }
            if (eq_stats[i].size() > 1) {
                eq_statFile << "-------------------------------------\nProduct of percentages (%) = " << prod*100 << "\n-------------------------------------" << std::endl;
            }
        }
        eq_statFile << std::endl;
        eq_statFile.close();
	
	long double sys_U = 0.0, measure_count = 0.0, sys_ntot = 0.0;
	std::vector < double > sys_comp (sys.nSpecies(), 0.0);
	long long int printCounter = 0;
	for (unsigned long long int move = 0; move < prSteps; ++move) {
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
			sys_U += sys.energy();
	
			// record composition
			for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
				sys_comp[i] += sys.numSpecies[i];
				sys_ntot += sys.numSpecies[i];
			}

			measure_count += 1.0;
		}

		// Periodically write out checkpoints and report statistics
		if (move%movePrint == 0) {
			sys.printSnapshot("snap_"+sstr(printCounter)+".xyz", "");
			printCounter++;
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
	
    	// Report move statistics for final "production" stage
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
    
	std::cout << "<U> : " << sys_U/measure_count << std::endl;
	std::cout << "<N_tot> : " << sys_ntot/measure_count << std::endl;
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		std::cout << "x_"+sstr(i+1)+" = <N_"+sstr(i+1)+">/<N_tot> : " << sys_comp[i]/sys_ntot << std::endl;
	}
	
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
