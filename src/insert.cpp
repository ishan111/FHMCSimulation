#include "insert.h"

/*!
 * Insert a particle into the system.  All other information is stored in the simSystem object.
 * 
 * \param [in] sys System object to attempt to insert a particle into.
 * 
 * \return MOVE_SUCCESS if inserted a particle, otherwise MOVE_FAILURE if did not.  Will throw exceptions if there was an error.
 */
int insertParticle::make (simSystem &sys) {
    // check if at upper bound for this specific species type
    if (sys.numSpecies[typeIndex_] >= sys.maxSpecies(typeIndex_)) {
        return MOVE_FAILURE;
    }

    // also check if at upper bound for total number of atoms
    if (sys.getTotN() >= sys.totNMax()) {
    	return MOVE_FAILURE;
    }
    
	const std::vector < double > box = sys.box();
	double V = 1.0;
	for (unsigned int i = 0; i < box.size(); ++i) {
		V *= box[i];
	}
	
	double insEnergy = 0.0;
	int origState = 0;
	bool createdAtom = false;

	// reference N state
	long long int nHigh = sys.numSpecies[typeIndex_]+1;
	
    atom* newAtom;
    if (sys.getCurrentM() == 0) {
    	// attempt to insert a brand new one
    	newAtom = new atom;
    	createdAtom = true;
    	origState = 0;
    	if (sys.getTotalM() > 1) {
    		newAtom->mState = 1; // incremental insertion state if doing expanded ensemble, otherwise leave as 0
    	}
    	for (unsigned int i = 0; i < box.size(); ++i) {
    		newAtom->pos[i] = rng (&RNG_SEED) * box[i];
    	}
    } else {
    	// continue to try to insert the partially inserted one
    	// which always stored at the "end" of the list
    	// don't increment the state yet
    	newAtom = sys.getFractionalAtom(); // mcMove object guarantees we are only making this move if the fractional atom if type typeIndex_ 
    	origState = newAtom->mState;

    	// if doing expanded ensemble and one is already partially inserted, have to get baseline, else this baseline is 0
    	for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
        	// get positions of neighboring atoms around newAtom
        	std::vector < atom* > neighborAtoms = sys.getNeighborAtoms (spec, typeIndex_, newAtom);
            for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
    			try {
    				insEnergy -= sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], newAtom, box);	
    			} catch (customException& ce) {
    				std::string a = "Cannot insert because of energy error: ", b = ce.what();
    				if (createdAtom) {
    					delete newAtom;
    				}
    				throw customException (a+b);
    			}
            }
            
            // neglect all tail corrections for partially inserted particles
        }
        
        // now increment the expanded ensemble state after baseline has been calculated
        newAtom->mState += 1;
        if (newAtom->mState == sys.getTotalM()) {
        	newAtom->mState = 0;
        }
    }
    
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    	// get positions of neighboring atoms around newAtom
    	std::vector < atom* > neighborAtoms = sys.getNeighborAtoms (spec, typeIndex_, newAtom);
        for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
			try {
				insEnergy += sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], newAtom, box);	
			} catch (customException& ce) {
				std::string a = "Cannot insert because of energy error: ", b = ce.what();
				if (createdAtom) {
					delete newAtom;
				}
				throw customException (a+b);
			}
        }
        // add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
        	if (newAtom->mState == 0) { // the mState was updated above to be what the atom will be if move is accepted
        		// if current atom becomes a full atom, include tail corrections
        		insEnergy += sys.ppot[spec][typeIndex_]->tailCorrection(sys.numSpecies[spec]/V);
        	} 
        }
#endif
    }
    
    // restore the original mState of the newAtom
    newAtom->mState = origState;
    
    // biasing
    double dN = 1.0/sys.getTotalM();
    const double p_u = pow(V/nHigh, dN)*exp(sys.beta()*(sys.mu(typeIndex_)*dN - insEnergy));
    int nTotFinal = sys.getTotN(), mFinal = sys.getCurrentM() + 1;
    if (sys.getCurrentM() == sys.getTotalM()-1) {
    	nTotFinal++;
    	mFinal = 0;
    }
    double bias = calculateBias(sys, nTotFinal, mFinal);
   
    // tmmc gets updated the same way, regardless of whether the move gets accepted
    if (sys.useTMMC) {
    	sys.tmmcBias->updateC (sys.getTotN(), nTotFinal, sys.getCurrentM(), mFinal, std::min(1.0, p_u)); 
    }
    
	// metropolis criterion
	if (rng (&RNG_SEED) < p_u*bias) {
        try {
            sys.insertAtom(typeIndex_, newAtom);
        } catch (customException &ce) {
            std::string a = "Failed to insert atom: ", b = ce.what();
			if (createdAtom) {
				delete newAtom;
			}
            throw customException (a+b);
        }
		sys.incrementEnergy(insEnergy);
		
		// update Wang-Landau bias, if used
		if (sys.useWALA) {
			sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
		}
		
		if (createdAtom) {
			delete newAtom;
		}
		
        return MOVE_SUCCESS;
    }
    
	// update Wang-Landau bias (even if moved failed), if used
	if (sys.useWALA) {
		sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
	}
	
	if (createdAtom) {
		delete newAtom;
	}
	
	return MOVE_FAILURE;
}
