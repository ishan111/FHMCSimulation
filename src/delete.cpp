#include "delete.h"

/*!
 * Delete a particle from the system.  All other information is stored in the simSystem object.
 * 
 * \param [in] sys System object to attempt to remove a particle from.
 * 
 * \return MOVE_SUCCESS if deleted a particle, otherwise MOVE_FAILURE if did not.  Will throw exceptions if there was an error.
 */
int deleteParticle::make (simSystem &sys) {
	// check if any can be deleted from this species
    if (sys.numSpecies[typeIndex_] < sys.minSpecies(typeIndex_)) {
        return MOVE_FAILURE;
    }
    if (sys.numSpecies[typeIndex_] == sys.minSpecies(typeIndex_) && sys.getCurrentM() == 0) {
        return MOVE_FAILURE;
    }

    // also check if at global bound on total number of particles
    if (sys.getTotN() < sys.totNMin()) {
    	return MOVE_FAILURE;
    }
    if (sys.getTotN() == sys.totNMin() && sys.getCurrentM() == 0) { // move class guarantees only operating on the correct species already
        return MOVE_FAILURE;
    }
    
	const std::vector < double > box = sys.box();
	double V = 1.0;
	for (unsigned int i = 0; i < box.size(); ++i) {
		V *= box[i];
	}
	
	double delEnergy = 0.0;

	// initial guess at the N state we are coming from
        long long int nHigh = sys.numSpecies[typeIndex_];

	atom* chosenAtom; 
    if (sys.getCurrentM() == 0) {
    	// pick a brand new one to delete
    	chosenAtom = &sys.atoms[typeIndex_][(int) floor(rng (&RNG_SEED) * sys.numSpecies[typeIndex_])];
	nHigh = sys.numSpecies[typeIndex_];
    } else {
    	// continue to try to delete the partially deleted one
    	chosenAtom = sys.getFractionalAtom(); // mcMove guarantees this move only being made if fractional atom of type typeIndex_
    	nHigh = sys.numSpecies[typeIndex_] + 1; // again mcMove guarantees this species is the fractional one, reference has to be at the next fully inserted level
	}
    
    // get baseline as the particle currently is
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
        // get positions of neighboring atoms around chosenAtom
        std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, chosenAtom);
        for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
            try {
				delEnergy -= sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], chosenAtom, box);
			}
			catch (customException& ce) {
				std::string a = "Cannot delete because of energy error: ", b = ce.what();
				throw customException (a+b);
			}
        }
        // add tail correction to potential energy -- only enable for fluid phase simulations       
#ifdef FLUID_PHASE_SIMULATIONS
        if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
        	if (chosenAtom->mState == 0) { 
        		// if current atom is a full atom right now, include tail corrections
        		if (spec == typeIndex_) {
        			delEnergy -= sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec]-1)/V);
        		} else {
        			delEnergy -= sys.ppot[spec][typeIndex_]->tailCorrection(sys.numSpecies[spec]/V);
        		}
        	} 
        }
#endif
    }
    
    // if the particle is about to be completely removed, no further calculation is required
    if (chosenAtom->mState != 1 && sys.getTotalM() > 1) { // if 1, it is just completely removed - otherwise have to do calculation since in expanded ensemble if M > 1
    	// temporarily decrement the expanded ensemble state on the atom
    	int orig_state = chosenAtom->mState; 	
    	chosenAtom->mState -= 1;
    	if (chosenAtom->mState < 0) {
    		chosenAtom->mState = sys.getTotalM() - 1;
    	}
    	
        for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
            // get positions of neighboring atoms around chosenAtom
            std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, chosenAtom);
            for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
                try {
    				delEnergy += sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], chosenAtom, box);
    			} catch (customException& ce) {
    				std::string a = "Cannot delete because of energy error: ", b = ce.what();
    				throw customException (a+b);
    			}
            }
            // no tail corrections for partially inserted particles
        }
        
        // restore the expanded ensemble state
        chosenAtom->mState = orig_state;
    }
    
    // biasing
    double dN = 1.0/sys.getTotalM();
    const double p_u = pow(nHigh/V, dN)*exp(sys.beta()*(-sys.mu(typeIndex_)*dN - delEnergy));
    int nTotFinal = sys.getTotN(), mFinal = sys.getCurrentM() - 1;
    if (sys.getCurrentM() == 0) {
    	nTotFinal--;
    	mFinal = sys.getTotalM() - 1;
    }
    double bias = calculateBias(sys, nTotFinal, mFinal); 
    
    // tmmc gets updated the same way, regardless of whether the move gets accepted
    if (sys.useTMMC) {
    	sys.tmmcBias->updateC (sys.getTotN(), nTotFinal, sys.getCurrentM(), mFinal, std::min(1.0, p_u)); // also has to be function of N and M now
    }
    
	// metropolis criterion
	if (rng (&RNG_SEED) < p_u*bias) {
		int counter = 0;
		for (std::vector<atom>::iterator it = sys.atoms[typeIndex_].begin(); it != sys.atoms[typeIndex_].end(); ++it) {
			if (&(*it) == chosenAtom) {
				break;
			} else {
				counter++;
			}
		}
	    try {
            sys.deleteAtom(typeIndex_, counter);
        } catch (customException &ce) {
            std::string a = "Failed to delete atom: ", b = ce.what();
            throw customException (a+b);
        }
		sys.incrementEnergy(delEnergy);	
		
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
