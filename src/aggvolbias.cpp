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
    if (nAvail1 < 1 || nAvail2 < 1 || (typeIndex_ == typeIndex2_ && nAvail1 <= 1)) {
        if (sys.useWALA) {
            sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
        }
        if (sys.useTMMC) {
            sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), 0.0);
        }
        return MOVE_FAILURE;
    }

    // choose particles, j and k - at present we are guaranteed at least 2 unique particles in the system
    int pkJ = (int) floor(rng (&RNG_SEED) * nAvail1), pkK = 0, iterMax = 25, iter = 0;
    bool jkOverlap = true;
    while (jkOverlap && iter < iterMax) {
        // if there are only a few particles in the system, they may all overlap so reject after a number of trials
        iter++;
        
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
    
    // check for failure to find pkJ and pkK
    if (iter >= iterMax) {
        if (sys.useWALA) {
            sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
        }
        if (sys.useTMMC) {
            sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), 0.0);
        }
        return MOVE_FAILURE;
    }
    
    const std::vector < double > box = sys.box();
    double minL = box[0];
    double Vtot = 1.0;
    for (unsigned int i = 0; i < box.size(); ++i) {
        Vtot *= box[i];
        minL = std::min(minL, box[i]);
    }
    
    // sanity check for rc's
    if (!(rc1_ < minL/2.0) {
        throw customException ("Neighborhood radius for species 1 in aggVolBias must be < box/2");
    }
    if (!(rc2_ < minL/2.0) {
        throw customException ("Neighborhood radius for species 2 in aggVolBias must be < box/2");
    }
    
    // based on the choices made below, the unbiased acceptance probability will be set in each case
    double p_u = 1.0, bias = 1.0, dU = 0.0;
    int chosenAtomType = 0;
    atom* chosenAtom;
    atom tmpNewAtom; // stores position the chosenAtom is being proposed to be moved to
        
    if (rng(&RNG_SEED) < pBias_) {
        // choose a particle either "in" k or "out" j with equal probability
        // note that "out" j includes the "in" k region as well
        if (rng (&RNG_SEED) < 0.5) {
            // choose particle "in" k, this chosen particle can be of any type
            std::vector < atom* > neighborAtoms;
            neighborAtoms.reserve(100); // 100 is just an arbitrary number to help accelerate
            std::vector < int > nEachType (sys.nSpecies(), 0);
            for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
                int totKatoms = sys.numSpecies[i];
                
                // account for if in expanded ensemble and have an additional partially inserted particle floating around
                if (sys.getCurrentM() > 0 && sys.getFactionalAtomType() == i) {
                    totKatoms++;
                }
                
                for (unsigned int j = 0; j < totKatoms; ++j) {
                    if (pbc_dist2(sys.atoms[i][j]->pos, sys.atoms[typeIndex2_][pkK]->pos, sys.box()) < rc2_*rc2_) {
                        if (sys.atoms[typeIndex2_][pkK] != sys.atoms[i][j]) { // since J and K do not overlap do not have to check if this is pkJ
                            neighborAtoms.push_back(&sys.atoms[i][j]);
                            nEachType[i] += 1;
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
            chosenAtom = neighborAtoms[atomIndex];
            int tot = 0;
            while (tot < atomIndex) {
                tot += nEachType[chosenAtomType];
                chosenAtomType++;
                if (chosenAtomType > sys.nSpecies()) { // > not >= because chosenAtomType will be decremented next
                    throw customException ("Error, could not properly identify the type of the atom chosen to be moved in aggVolBias");
                }
            }
            chosenAtomType--;
            if (chosenAtomType < 0) { // in case atomIndex = 0 and above loop did not execute at all
                chosenAtomType = 0;
            }
        } else {
            // choose particle "out" of j, this chosen particle can be of any type
            
            // just pick atoms (of any type) at random and see if it is "out" of j
            // this should be faster than establishing, a priori, all atoms which are "out"
            // and then picking from that list since *most* atoms should be "out"
            // unlike when we have to pick from "in" ones which are more rare
            
            int iter = 0, iterMax = 25;
            bool inJ = true;
            while (inJ && iter < iterMax) {
                iter++;
                const int ranSpec = (int) floor(rng (&RNG_SEED) * sys.nSpecies());
                int availAtoms = sys.numSpecies[ranSpec];
                // account for if in expanded ensemble and have an additional partially inserted particle floating around
                if (sys.getCurrentM() > 0 && sys.getFactionalAtomType() == ranSpec) {
                    availAtoms++;
                }
                int ranIndex = (int) floor(rng (&RNG_SEED) * availAtoms);
                
                // check this atom is neither pkJ nor pkK
                if (sys.atoms[ranSpec][ranIndex] != sys.atoms[typeIndex_][pkJ] && sys.atoms[ranSpec][ranIndex] != sys.atoms[typeIndex2_][pkK]) {
                    if (!(pbc_dist2(sys.atoms[ranSpec][ranIndex]->pos, sys.atoms[typeIndex1_][pkJ]->pos, sys.box()) < rc1_*rc1_)) {
                        chosenAtom = sys.atoms[ranSpec][ranIndex];
                        chosenAtomType = ranSpec;
                        inJ = false;
                    }
                }
            }
            
            // check if we could not locate a particle "out" of j
            if (iter >= iterMax) {
                if (sys.useWALA) {
                    sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
                }
                if (sys.useTMMC) {
                    sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), 0.0);
                }
                return MOVE_FAILURE;
            }
        }
        
        // move the chosenAtom "in" j
        
        // (1) get energy of chosenAtom in current state
        double oldEnergy = 0.0;
        try {
            oldEnergy = getTempEnergy_ (sys, box, chosenAtomType, chosenAtom);
        } catch (customException &ce) {
            throw customException (ce.what());
        }
        
        // (2) "move" chosenAtom - randomly choose a radius from [0, rc1) around pkJ <-- Could be a bit smarter with this (adjust V_in later) if HS, but for now just leave as in
        const double magnitude = rng (&RNG_SEED)*rc1_;
        
        // then choose a point randomly on the surface of that sphere to place chosenAtom
        std::vector < double > surfaceVec = random3DSurfaceVector (const double magnitude), origPos = chosenAtom->pos;
        for (unsigned int i = 0; i < newPos.size(); ++i) {
            chosenAtom->pos[i] = sys.atoms[typeIndex_][pkJ].pos[i] + surfaceVec[i];
            
            // apply periodic boundary conditions
            if (chosenAtom->pos[i] >= box[i]) {
                chosenAtom->pos[i] -= box[i];
            } else if (chosenAtom->pos[i] < 0) {
                chosenAtom->pos[i] += box[i];
            }
        }
        
        // store the newly proposed position for later
        tmpNewAtom.pos = chosenAtom->pos;
        
        double newEnergy = 0.0;
        try {
            newEnergy = getTempEnergy_ (sys, box, chosenAtomType, chosenAtom);
        } catch (customException &ce) {
            throw customException (ce.what());
        }
        
        // restore the atom's original state before continuing
        chosenAtom->pos = origPos;
        
        // assign p_u
        dU = newEnergy - oldEnergy;
        p_u =
        
        
        
    } else {
		// choose a particle "in" j, this chosen particle can be of any type
        std::vector < atom* > neighborAtoms;
        neighborAtoms.reserve(100); // 100 is just an arbitrary number to help accelerate
        std::vector < int > nEachType (sys.nSpecies(), 0);
        for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
            int totJatoms = sys.numSpecies[i];
            
            // account for if in expanded ensemble and have an additional partially inserted particle floating around
            if (sys.getCurrentM() > 0 && sys.getFactionalAtomType() == i) {
                totJatoms++;
            }
            
            for (unsigned int j = 0; j < totJatoms; ++j) {
                if (pbc_dist2(sys.atoms[i][j]->pos, sys.atoms[typeIndex1_][pkJ]->pos, sys.box()) < rc1_*rc1_) {
                    if (sys.atoms[typeIndex1_][pkJ] != sys.atoms[i][j]) { // since J and K do not overlap do not have to check this is pkK
                        neighborAtoms.push_back(&sys.atoms[i][j]);
                        nEachType[i] += 1;
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
        chosenAtom = neighborAtoms[atomIndex];
        int tot = 0;
        while (tot < atomIndex) {
            tot += nEachType[chosenAtomType];
            chosenAtomType++;
            if (chosenAtomType > sys.nSpecies()) { // > not >= because chosenAtomType will be decremented next
                throw customException ("Error, could not properly identify the type of the atom chosen to be moved in aggVolBias");
            }
        }
        chosenAtomType--;
        if (chosenAtomType < 0) { // in case atomIndex = 0 and above loop did not execute at all
            chosenAtomType = 0;
        }
        
        // what to do with chosenAtom?
        if (rng (&RNG_SEED) < 0.5) {
            // move this chosenAtom "out" of pkJ
            bool inJ = true;
            while (inJ) {
                for (unsigned int i = 0; i < box.size(); ++i) {
                    tmpNewAtom.pos[i] = rng (&RNG_SEED) * box[i];
                }
                
                // check that this position is "out" of J
                if (pbc_dist2(tmpNewAtom.pos, sys.atoms[typeIndex_][pkJ]->pos, sys.box()) > rc1_*rc1_) {
                    inJ = false;
                }
            }
            
            // calculate energy of chosenAtom in current configuration
            double oldEnergy = 0.0;
            try {
                oldEnergy = getTempEnergy_ (sys, box, chosenAtomType, chosenAtom);
            } catch (customException &ce) {
                throw customException (ce.what());
            }
            
            // "move" chosenAtomInJ to new location and get energy (amounts to same algorithm as a translation)
            std::vector < double > origPos = chosenAtom->pos;
            chosenAtom->pos = tmpNewAtom.pos;
            
            double newEnergy = 0.0;
            try {
                newEnergy = getTempEnergy_ (sys, box, chosenAtomType, chosenAtom);
            } catch (customException &ce) {
                throw customException (ce.what());
            }
            
            // restore atom to original position
            chosenAtom->pos = origPos;
            
            // assign p_u
            dU = newEnergy - oldEnergy;
            p_u =
            
            
            
            
        } else {
            // move this chosenAtom "in" pkK
            
            // calculate the energy of its old location
            double oldEnergy = 0.0;
            try {
                oldEnergy = getTempEnergy_ (sys, box, chosenAtomType, chosenAtom);
            } catch (customException &ce) {
                throw customException (ce.what());
            }
            
            // randomly choose a radius from [0, rc2) around pkK <-- Could be a bit smarter with this (adjust V_in later) if HS, but for now just leave as in
            const double magnitude = rng (&RNG_SEED)*rc2_;
            
            // then choose a point randomly on the surface of that sphere to place chosenAtom
            std::vector < double > surfaceVec = random3DSurfaceVector (const double magnitude), origPos = chosenAtom->pos;
            for (unsigned int i = 0; i < newPos.size(); ++i) {
                chosenAtom->pos[i] = sys.atoms[typeIndex2_][pkK].pos[i] + surfaceVec[i];
                
                // apply periodic boundary conditions
                if (chosenAtom->pos[i] >= box[i]) {
                    chosenAtom->pos[i] -= box[i];
                } else if (chosenAtom->pos[i] < 0) {
                    chosenAtom->pos[i] += box[i];
                }
            }
            
            // store for later
            tmpNewAtom.pos = chosenAtom->pos;
            
            // move, recalculate energy
            double newEnergy = 0.0;
            try {
                newEnergy = getTempEnergy_ (sys, box, chosenAtomType, chosenAtom);
            } catch (customException &ce) {
                throw customException (ce.what());
            }
            
            // restore the chosenAtoms position
            chosenAtom->pos = origPos;
            
            // assign p_u
            dU = newEnergy - oldEnergy;
            p_u =
            
            
            
        }
	}
	
    // biasing
    double bias = calculateBias(sys, sys.getTotN(), sys.getCurrentM()); // N_tot doesn't change throughout this move
        
    // tmmc gets updated the same way, regardless of whether the move gets accepted
    if (sys.useTMMC) {
        sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), std::min(1.0, p_u)); // since the total number of atoms isn't changing, can use getTotN() as both initial and final states
    }
        
    if (rng (&RNG_SEED) < p_u*bias) {
        try {
            // atoms were not modified overall so must do that here
            std::vector < double > origPos = chosenAtom->pos;
            chosenAtom->pos = tmpNewAtom.pos; // move the atom
            sys.translateAtom(chosenAtomType, chosenAtom, origPos);
        } catch (customException &ce) {
            std::string a = "Failed to move atom in aggVolBias: ", b = ce.what();
            throw customException (a+b);
        }
        sys.incrementEnergy(dU);
            
        // update Wang-Landau bias, if used
        if (sys.useWALA) {
            sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
        }
        
        return MOVE_SUCCESS;
    }
        
    // update Wang-Landau bias (even if moved failed), if used
    if (sys.useWALA) {
        sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
    }
        
    return MOVE_FAILURE;
}
        
/*!
 * Helper function for aggregation volume bias calculation which just calculates system's interaction energy of a particle before or after a move.
 *
 * \param [in] sys System undergoing the proposed change
 * \param [in] box System box size
 * \param [in] chosenAtomType Type index of the atom chosen to be moved
 * \param [in] chosenAtom Pointer to the atom chosen to move
 *
 * \return tmpEnergy Energy of the system in its current configuration
 */
double aggVolBias::getTempEnergy_ (const simSystem &sys, const std::vector < double > &box, const int chosenAtomType, const atom* chosenAtom) {
    double tmpEnergy = 0.0;
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
        // get positions of neighboring atoms around chosenAtom
        std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, chosenAtomType, chosenAtom);
        for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
            try {
                tmpEnergy += sys.ppot[spec][chosenAtomType]->energy(neighborAtoms[i], chosenAtom, box);
            } catch (customException& ce) {
                std::string a = "Cannot perform aggVolBias move because of energy error: ", b = ce.what();
                throw customException (a+b);
            }
        }
        
        // add tail correction to potential energy
#ifdef FLUID_PHASE_SIMULATIONS
        if (sys.ppot[spec][chosenAtomType]->useTailCorrection) {
            if (!(sys.getCurrentM() > 0 && sys.getFractionalAtom() == chosenAtom)) {
                // then chosenAtom is not a partially inserted particle and tail interactions must be included
                if (spec == chosenAtomType) {
                    if (sys.numSpecies[spec]-1 > 0) {
                        tmpEnergy += sys.ppot[spec][chosenAtomType]->tailCorrection((sys.numSpecies[spec]-1)/V);
                    }
                } else {
                    if (sys.numSpecies[spec] > 0) {
                        tmpEnergy += sys.ppot[spec][chosenAtomType]->tailCorrection(sys.numSpecies[spec]/V);
                    }
                }
            }
        }
#endif
    }
    
    return tmpEnergy;
}