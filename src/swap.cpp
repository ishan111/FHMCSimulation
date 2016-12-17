#include "swap.h"

/*!
 * Swap two particles of different types in the system.  All other information is stored in the simSystem object.
 *
 * \param [in] sys System object to attempt to swap particles in.
 *
 * \return MOVE_SUCCESS if particles are swapped, otherwise MOVE_FAILURE if not.  Will throw exceptions if there was an error.
 */
int swapParticles::make (simSystem &sys) {
	// Choose an atom of each type to try to exchange
	int n1Avail = sys.numSpecies[typeIndex_];
	if (sys.getCurrentM() > 0 && sys.getFractionalAtomType() == typeIndex_) {
		n1Avail++;
	}

	int n2Avail = sys.numSpecies[typeIndex2_];
	if (sys.getCurrentM() > 0 && sys.getFractionalAtomType() == typeIndex2_) {
		n2Avail++;
	}

	if (n1Avail < 1 || n2Avail < 1) {
		// updates to biasing functions must be done even if at bounds
        if (sys.useWALA) {
            sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
        }
        if (sys.useTMMC) {
            sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), 0.0);
        }
        return MOVE_FAILURE;
	}

	// because the locations are effectively being swapped, it is fair to select a partially inserted atom to be involved
	const int a1 = (int) floor(rng (&RNG_SEED) * n1Avail);
	const int a2 = (int) floor(rng (&RNG_SEED) * n2Avail);
	atom a1_orig = sys.atoms[typeIndex_][a1], a1_new = a1_orig;
	atom a2_orig = sys.atoms[typeIndex2_][a2], a2_new = a2_orig;

	// positions will be exchanged, but no other property should change
	a1_new.pos = a2_orig.pos;
	a2_new.pos = a1_orig.pos;

	const std::vector < double > box = sys.box();
	double V = 1.0;
	for (unsigned int i = 0; i < box.size(); ++i) {
        	V *= box[i];
    	}

    	double delEnergy = 0.0;
    	for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
        	// get positions of neighboring atoms around a1
        	std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, &sys.atoms[typeIndex_][a1]);
        	for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
        		if (neighborAtoms[i] == &sys.atoms[typeIndex2_][a2]) {
        			// skip their interaction as they were already deleted - this is fine because their net pairwise interaction doesn't change over the
        			// course of the simulation - this means we don't have to actually re-insert the other atom / get its interation somehow later
        			continue;
        		} else {
				try {
					delEnergy += sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], &sys.atoms[typeIndex_][a1], box);
				} catch (customException& ce) {
					std::string a = "Cannot delete because of energy error: ", b = ce.what();
					throw customException (a+b);
				}
        		}
        	}
        	// add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        	if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
        		if (!(sys.getCurrentM() > 0 && sys.getFractionalAtom () == &sys.atoms[typeIndex_][a1])) {
        			// then a1 is not the partially inserted particle and tail interactions must be included
        			if (spec == typeIndex_) {
                    			if (sys.numSpecies[spec]-1 > 0) {
                        			delEnergy += sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec]-1)/V);
                    			}
        			} else {
                    			if (sys.numSpecies[spec] > 0) {
                        			delEnergy += sys.ppot[spec][typeIndex_]->tailCorrection(sys.numSpecies[spec]/V);
                    			}
        			}
        		}
        	}
#endif
	}

    // account for wall interaction energy
    delEnergy += sys.speciesBarriers[typeIndex_].energy(&sys.atoms[typeIndex_][a1], box);

    	for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
        	// get positions of neighboring atoms around a2
        	std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex2_, &sys.atoms[typeIndex2_][a2]);
        	for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
        		if (neighborAtoms[i] == &sys.atoms[typeIndex_][a1]) {
        			// skip their interaction as they were already deleted - this is fine because their net pairwise interaction doesn't change over the
        			// course of the simulation - this means we don't have to actually re-insert the other atom / get its interation somehow later
        			continue;
        		} else {
				try {
					delEnergy += sys.ppot[spec][typeIndex2_]->energy(neighborAtoms[i], &sys.atoms[typeIndex2_][a2], box);
				} catch (customException& ce) {
					std::string a = "Cannot delete because of energy error: ", b = ce.what();
					throw customException (a+b);
				}
        		}
        	}
        	// add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        	if (sys.ppot[spec][typeIndex2_]->useTailCorrection) {
        		if (!(sys.getCurrentM() > 0 && sys.getFractionalAtom () == &sys.atoms[typeIndex2_][a2])) {
        			// then a2 is not the partially inserted particle and tail interactions must be included
        			if (spec == typeIndex2_) {
                    			if (sys.numSpecies[spec]-1 > 0) {
                        			delEnergy += sys.ppot[spec][typeIndex2_]->tailCorrection((sys.numSpecies[spec]-1)/V);
                    			}
        			} else {
                    			if (sys.numSpecies[spec] > 0) {
                        			delEnergy += sys.ppot[spec][typeIndex2_]->tailCorrection(sys.numSpecies[spec]/V);
                    			}
        			}
        		}
        	}
#endif
    	}

        // account for wall interaction energy
        delEnergy += sys.speciesBarriers[typeIndex2_].energy(&sys.atoms[typeIndex2_][a2], box);

    	double insEnergy = 0.0;
    	for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    		// get positions of neighboring atoms around a1's (a2's) new (old) location
    		std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, &a1_new);
        	for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
        		// with these new "copy atoms" getNeighborAtoms can't guarantee it doesn't point to old self, so must check
        		if ((neighborAtoms[i] == &sys.atoms[typeIndex2_][a2]) || (neighborAtoms[i] == &sys.atoms[typeIndex_][a1])) {
        			continue;
        		} else {
				try {
					insEnergy += sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], &a1_new, box);
				} catch (customException& ce) {
					std::string a = "Cannot insert because of energy error: ", b = ce.what();
					throw customException (a+b);
				}
        		}
        	}
        	// add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        	if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
        		if (!(sys.getCurrentM() > 0 && sys.getFractionalAtom () == &sys.atoms[typeIndex_][a1])) {
        			// then a1 is not the partially inserted particle and tail interactions must be included
        			if (spec == typeIndex_) {
                    			if (sys.numSpecies[spec]-1 > 0) {
                        			insEnergy += sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec]-1)/V);
                   			}
        			} else {
                    			if (sys.numSpecies[spec] > 0) {
                        			insEnergy += sys.ppot[spec][typeIndex_]->tailCorrection(sys.numSpecies[spec]/V);
                    			}
        			}
        		}
        	}
#endif
    	}
        // account for wall interaction energy
        insEnergy += sys.speciesBarriers[typeIndex_].energy(&a1_new, box);

    	for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    		// get positions of neighboring atoms around a2's (a1's) new (old) location
    		std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex2_, &a2_new);
        	for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
        		if ((neighborAtoms[i] == &sys.atoms[typeIndex_][a1]) || (neighborAtoms[i] == &sys.atoms[typeIndex2_][a2])) {
        			continue;
        		} else {
				try {
					insEnergy += sys.ppot[spec][typeIndex2_]->energy(neighborAtoms[i], &a2_new, box);
				} catch (customException& ce) {
					std::string a = "Cannot insert because of energy error: ", b = ce.what();
					throw customException (a+b);
				}
        		}
        	}
        	// add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        	if (sys.ppot[spec][typeIndex2_]->useTailCorrection) {
        		if (!(sys.getCurrentM() > 0 && sys.getFractionalAtom () == &sys.atoms[typeIndex2_][a2])) {
        			// then a2 is not the partially inserted particle and tail interactions must be included
        			if (spec == typeIndex2_) {
                    			if (sys.numSpecies[spec]-1 > 0) {
                        			insEnergy += sys.ppot[spec][typeIndex2_]->tailCorrection((sys.numSpecies[spec]-1)/V);
                    			}
        			} else {
                    			if (sys.numSpecies[spec] > 0) {
                        			insEnergy += sys.ppot[spec][typeIndex2_]->tailCorrection(sys.numSpecies[spec]/V);
                    			}
        			}
        		}
        	}
#endif
    	}
        // account for wall interaction energy
        insEnergy += sys.speciesBarriers[typeIndex2_].energy(&a2_new, box);

	// Biasing
    	const double p_u = exp(-sys.beta()*(insEnergy - delEnergy));
    	double bias = calculateBias(sys, sys.getTotN(), sys.getCurrentM());

    	// tmmc gets updated the same way, regardless of whether the move gets accepted
    	if (sys.useTMMC) {
    		sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), std::min(1.0, p_u));
    	}

	if (rng (&RNG_SEED) < p_u*bias) {
		sys.incrementEnergy(insEnergy - delEnergy);

		// swap the particles by deleting/reinserting

        // could probably make this faster by using translateAtom instead of full insert/delete?

		// -a1 completely
	    try {
	    	sys.deleteAtom(typeIndex_, a1, true);
	    } catch (customException &ce) {
	    	std::string a = "Failed to delete atom during swapping: ", b = ce.what();
	    	throw customException (a+b);
	    }

	    // -a2 completely
	    try {
	    	sys.deleteAtom(typeIndex2_, a2, true);
	    } catch (customException &ce) {
	    	std::string a = "Failed to delete atom during swapping: ", b = ce.what();
	    	throw customException (a+b);
	    }

	    // +a1_new completely
	    try {
	    	sys.insertAtom(typeIndex_, &a1_new, true);
	    } catch (customException &ce) {
	    	std::string a = "Failed to insert atom during swapping: ", b = ce.what();
	    	throw customException (a+b);
	    }

	    // +a2_new completely
	    try {
	    	sys.insertAtom(typeIndex2_, &a2_new, true);
	    } catch (customException &ce) {
	    	std::string a = "Failed to insert atom during swapping: ", b = ce.what();
	    	throw customException (a+b);
	    }

		// update Wang-Landau bias, if used
		if (sys.useWALA) {
			sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
		}

		// double check
		/*if (n1_orig != sys.numSpecies[typeIndex_]) {
			throw customException ("Number of species 1 atoms do not match before and after swap move");
		}
		if (n2_orig != sys.numSpecies[typeIndex2_]) {
			throw customException ("Number of species 2 atoms do not match before and after swap move");
		}
		if (fracType != sys.getFractionalAtomType()) {
			throw customException ("Fractional type has changed during course of swap move");
		}
		if (m_orig != sys.getCurrentM()) {
			throw customException ("Expanded ensemble state of system has changed during course of swap move");
		}*/

	return MOVE_SUCCESS;
    	}

	// update Wang-Landau bias (even if moved failed), if used
	if (sys.useWALA) {
		sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
	}

	return MOVE_FAILURE;
}
