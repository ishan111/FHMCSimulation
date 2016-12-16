#include "input.h"

/*!
 * Parse a json input file and initialize system object accordingly.
 *
 * \params [in] filename Input JSON document's filename
 * \params [in, out] sys System to initialize with barriers
 */
void parse_json (const std::string filename, simSystem &sys) {
    ;
}

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

	// cylinderZ (expect parameters: {x, y, radius, width, sigma, eps})
    for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		bool convention0 = false;
		std::string dummy = "cylinderZ_" + sstr(i+1);
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
			std::string dummy = "cylinderZ_" + sstr(i+1) + "_" + sstr(j);
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
}
