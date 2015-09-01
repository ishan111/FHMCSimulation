#include "aggvolbias.h"

/*!
 * Instantiate an aggregation-volume bias move to swap particles.
 *
 * \param [in] typeIndex Type to consider in swapping in the vicinity of
 * \param [in] typeIndex2 Second type to consider in swapping in the vicinity of
 * \param [in] pBias Biasing probability
 * \param [in] rc1 Radius around particle of typeIndex to consider as being the "in" region
 * \param [in] rc2 Radius around particle of typeIndex2 to consider as being the "in" region
 * \param [in] tag Name modifier to identify move to user
 */
aggVolBias::aggVolBias (const int typeIndex, const int typeIndex2, const double pBias, const double rc1, const double rc2, const std::string tag) {
	typeIndex_ = typeIndex; 
	typeIndex2_ = typeIndex2; 
	name_ = tag + sstr(typeIndex); 
	if (pBias < 1 && pBias > 0) {
		pBias_ = pBias;
	} else {
		throw customException ("Bias probability must be > 0 and < 1 for aggVolBias");
	}
    if (!(rc1 > 0.0)) {
        throw customException ("Neighborhood radius for species 1 in aggVolBias must be > 0");
    } else {
        rc1_ = rc1;
    }
    if (!(rc2 > 0.0)) {
        throw customException ("Neighborhood radius for species 2 in aggVolBias must be > 0");
    } else {
        rc2_ = rc2;
    }
    name_ = tag + sstr(typeIndex) + "_" + sstr(typeIndex2);
	changeN_ = false; 
}

/*!
 * Use aggregation-volume bias to swap two particles in the system.  All other information is stored in the simSystem object.
 * 
 * \param [in] sys System object to attempt to swap particles in.
 * 
 * \return MOVE_SUCCESS if particles are swapped, otherwise MOVE_FAILURE if not.  Will throw exceptions if there was an error.
 */
int aggVolBias::make (simSystem &sys) {
	// choose a particle of typeIndex_
	int nAvail1 = sys.numSpecies[typeIndex_];
	if (sys.getCurrentM() > 0 && sys.getFractionalAtomType() == typeIndex_) {
		nAvail1++;
	}

	// choose a particle of typeIndex2_
	int nAvail2 = sys.numSpecies[typeIndex2_];
	if (sys.getCurrentM() > 0 && sys.getFractionalAtomType() == typeIndex2_) {
                nAvail2++;
        }

	// update biases even upon failure to try move (also reject if types are same and only one overall particle)
	if (nAvail1 < 1 || nAvail2 < 1 || (typeIndex_ == typeIndex2_ && nAvail1 == 1)) {
		if (sys.useWALA) {
                        sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
                }
                if (sys.useTMMC) {
                        sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), 0.0);
                }
		return MOVE_FAILURE;
	}

	// choose particles, j and k - at present we are guaranteed at least 2 unique particles in the system
	int pkJ = (int) floor(rng (&RNG_SEED) * nAvail1), pkK = 0;
	bool jkOverlap = true;
	while (jkOverlap) {
		pkK = (int) floor(rng (&RNG_SEED) * nAvail2);
		
		// establish if j and k overlap their "in" regions
		const double rc = (rc1_ + rc2_)/2.0;
		
		// if pkK and pkJ are withing rcut of each other, they overlap
		if (pbc_dist2(sys.atoms[typeIndex_][pkJ]->pos, sys.atoms[typeIndex2_][pkK]->pos, sys.box()) > rc*rc) {
			jkOverlap = false;
		}

		// However, we also technically need to ensure that if typeIndex_ == typeIndex2_, j and k are not the same particle.  
		// Happily, this is already guaranteed by the above distance check (a particle is a distance 0 from itself)
	}
    
    const std::vector < double > box = sys.box();
    double minL = box[0];
    double Vtot = 1.0;
    for (unsigned int i = 0; i < box.size(); ++i) {
        Vtot *= box[i];
        minL = std::min(minL, box[i]);
    }
    
    // sanity check for rc's
    if (!(rc1_ < minL) {
        throw customException ("Neighborhood radius for species 1 in aggVolBias must be < box/2");
    }
    if (!(rc2_ < minL) {
        throw customException ("Neighborhood radius for species 2 in aggVolBias must be < box/2");
    }
    
    // based on the choices made below, the unbiased acceptance probability will be set in each case
    double p_u = 1.0, bias = 1.0, dU = 0.0;

    if (rng(&RNG_SEED) < pBias_) {
        // choose a particle "in" k or "out" j with equal probability
        // note that "out" j includes the "in" k region as well
        if (rng (&RNG_SEED) < 0.5) {
            // choose particle "in" k
            
        } else {
            // choose particle "out" of j
            
        }
    } else {
		// choose a particle "in" j, this chosen particle can be of any type
        std::vector < atom* > neighborAtoms;
        neighborAtoms.reserve(100); // 100 is just an arbitrary number to help accelerate
        for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
            int totJatoms = sys.numSpecies[i];
            
            // account for if in expanded ensemble and have an additional partially inserted particle floating around
            if (sys.getCurrentM() > 0 && sys.getFactionalAtomType() == i) {
                totJatoms++;
            }
            
            for (unsigned int j = 0; j < totJatoms; ++j) {
                if (pbc_dist2(sys.atoms[i][j]->pos, sys.atoms[typeIndex1_][pkJ]->pos, sys.box()) < rc1_*rc1_) {
                    if (sys.atoms[typeIndex1_][pkJ] != sys.atoms[i][j]) {
                        neighborAtoms.push_back(&sys.atoms[i][j]);
                    }
                }
            }
        }
        
        // reject move if no neighbors
        if (neighborAtoms.begin() == neighborAtoms.end()) {
            if (sys.useWALA) {
                sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
            }
            if (sys.useTMMC) {
                sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), 0.0);
            }
            return MOVE_FAILURE;
        }
        
        // otherwise choose an atom
        const int atomIndex = (int) floor(rng (&RNG_SEED) * neighborAtoms.size());
        atom* chosenAtomInJ = neighborAtoms[atomIndex];
        
        // what to do with chosenAtomInJ?
        if (rng (&RNG_SEED) < 0.5) {
            // move this chosenAtomInJ "out" of pkJ
            std::vector < double > outJpos (3, 0.0);
            bool inJ = true;
            while (inJ) {
                for (unsigned int i = 0; i < box.size(); ++i) {
                    outJpos[i] = rng (&RNG_SEED) * box[i];
                }
                
                // check that this position is "out" of J
                if (pbc_dist2(outJpos, sys.atoms[typeIndex_][pkJ]->pos, sys.box()) > rc1_*rc1_) {
                    inJ = false;
                }
            }
            
            // calculate energy of chosenAtomInJ in current configuration
            
            // "move" chosenAtomInJ to new location and get energy (amounts to same algorithm as a translation)
            
            // assign p_u
            
        } else {
            // move this chosenAtomInJ "in" pkK
            
            // randomly choose a radius from [0, rc2) <-- Could be a bit smarter with this (adjust V_in later) if HS, but for now just leave as in
            
            // then choose a point randomly on the surface of that sphere to place chosenAtomInJ
        }
	}
	
	// remember to upated TMMC and WL regardless of what happens


	return MOVE_FAILURE;
}
