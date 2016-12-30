#include "delete.h"

/*!
 * Delete a particle from the system.  All other information is stored in the simSystem object.
 *
 * \param [in] sys System object to attempt to remove a particle from.
 *
 * \return MOVE_SUCCESS if deleted a particle, otherwise MOVE_FAILURE if did not.  Will throw exceptions if there was an error.
 */
int deleteParticle::make (simSystem &sys) {
	bool earlyReject = false;

	// Check if any can be deleted from this species
	if (sys.numSpecies[typeIndex_] < sys.minSpecies(typeIndex_)) {
		earlyReject = true;
    }
    if (sys.numSpecies[typeIndex_] == sys.minSpecies(typeIndex_) && sys.getCurrentM() == 0) {
    	earlyReject = true;
    }

    // Also check if at global bound on total number of particles
    if (sys.getTotN() < sys.totNMin()) {
    	earlyReject = true;
    }
    if (sys.getTotN() == sys.totNMin() && sys.getCurrentM() == 0) { // Move class guarantees only operating on the correct species already
    	earlyReject = true;
    }

	// Updates to biasing functions must be done even if at bounds
    if (earlyReject) {
    	if (sys.useWALA) {
        	sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
    	}
        if (sys.useTMMC) {
        	int nTotFinal = sys.getTotN(), mFinal = sys.getCurrentM() - 1;
        	if (sys.getCurrentM() == 0) {
            	nTotFinal--;
        		mFinal = sys.getTotalM() - 1;
            }
        	sys.tmmcBias->updateC (sys.getTotN(), nTotFinal, sys.getCurrentM(), mFinal, 0.0);
    	}
        return MOVE_FAILURE;
    }

	const std::vector < double > box = sys.box();
	double V = 1.0;
	for (unsigned int i = 0; i < box.size(); ++i) {
		V *= box[i];
	}

	double delEnergy = 0.0;

	// Initial guess at the N state we are coming from
    long long int nHigh = sys.numSpecies[typeIndex_];

	atom* chosenAtom;
    if (sys.getCurrentM() == 0) {
    	// Pick a brand new one to delete
    	chosenAtom = &sys.atoms[typeIndex_][(int) floor(rng (&RNG_SEED) * sys.numSpecies[typeIndex_])];
        nHigh = sys.numSpecies[typeIndex_];
    } else {
    	// Continue to try to delete the partially deleted one
    	chosenAtom = sys.getFractionalAtom(); // mcMove guarantees this move only being made if fractional atom of type typeIndex_
    	nHigh = sys.numSpecies[typeIndex_] + 1; // Again mcMove guarantees this species is the fractional one, reference has to be at the next fully inserted level
	}

    // Get baseline as the particle currently is
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    	// Get positions of neighboring atoms around chosenAtom
    	std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, chosenAtom);
    	for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
    		try {
				delEnergy -= sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], chosenAtom, box);
			} catch (customException& ce) {
				std::string a = "Cannot delete because of energy error: ", b = ce.what();
				throw customException (a+b);
			}
        }

        // Add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
        	if (chosenAtom->mState == 0) {
        		// If current atom is a full atom right now, include tail corrections
        		if (spec == typeIndex_) {
        			if (sys.numSpecies[spec]-1 > 0) {
                    	delEnergy -= sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec]-1)/V);
                	}
        		} else {
                    if (sys.numSpecies[spec] > 0) {
                        delEnergy -= sys.ppot[spec][typeIndex_]->tailCorrection(sys.numSpecies[spec]/V);
                    }
                }
        	}
        }
#endif
    }

    // Also account for any wall or barrier interaction
    delEnergy -= sys.speciesBarriers[typeIndex_].energy(chosenAtom, box);

    // If the particle is about to be completely removed, no further calculation is required
    if (chosenAtom->mState != 1 && sys.getTotalM() > 1) { // If 1, it is just completely removed - otherwise have to do calculation since in expanded ensemble if M > 1
    	// Temporarily decrement the expanded ensemble state on the atom
    	int orig_state = chosenAtom->mState;
    	chosenAtom->mState -= 1;
    	if (chosenAtom->mState < 0) {
    		chosenAtom->mState = sys.getTotalM() - 1;
    	}

        for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
			// Get positions of neighboring atoms around chosenAtom
            std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, chosenAtom);
        	for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
                try {
    				delEnergy += sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], chosenAtom, box);
    			} catch (customException& ce) {
    				std::string a = "Cannot delete because of energy error: ", b = ce.what();
    				throw customException (a+b);
    			}
            }
        	// No tail corrections for partially inserted particles
        }

        // Also account for any wall or barrier interaction
        delEnergy += sys.speciesBarriers[typeIndex_].energy(chosenAtom, box);

        // Restore the expanded ensemble state
        chosenAtom->mState = orig_state;
    }

    // Biasing
    double dN = 1.0/sys.getTotalM();
    double p_u = 1.0;
    if (sys.addKECorrection()) {
        const double Lambda3 = pow(H2_2PI*sys.beta()/sys.mass(typeIndex_), 1.5);
        p_u = pow(Lambda3*nHigh/V, dN)*exp(sys.beta()*(-sys.mu(typeIndex_)*dN - delEnergy));
    } else {
        p_u = pow(nHigh/V, dN)*exp(sys.beta()*(-sys.mu(typeIndex_)*dN - delEnergy));
    }

    int nTotFinal = sys.getTotN(), mFinal = sys.getCurrentM() - 1;
    if (sys.getCurrentM() == 0) {
    	nTotFinal--;
    	mFinal = sys.getTotalM() - 1;
    	if (sys.addKECorrection()) {
    		delEnergy -= 1.5/sys.beta();
    	}
    }
    double bias = calculateBias(sys, nTotFinal, mFinal);

    // TMMC gets updated the same way, regardless of whether the move gets accepted
    if (sys.useTMMC) {
    	sys.tmmcBias->updateC (sys.getTotN(), nTotFinal, sys.getCurrentM(), mFinal, std::min(1.0, p_u)); // Also has to be function of N and M now
    }

	// Metropolis criterion
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

		// Update Wang-Landau bias, if used
		if (sys.useWALA) {
			sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
		}

        return MOVE_SUCCESS;
    }

	// Update Wang-Landau bias (even if moved failed), if used
	if (sys.useWALA) {
		sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
	}
	return MOVE_FAILURE;
}
