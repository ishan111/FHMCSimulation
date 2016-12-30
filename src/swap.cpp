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
		// Updates to biasing functions must be done even if at bounds
        if (sys.useWALA) {
            sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
        }
        if (sys.useTMMC) {
            sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), 0.0);
        }
        return MOVE_FAILURE;
	}

	// Because the locations are effectively being swapped, it is fair to select a partially inserted atom to be involved
	const int a1 = (int) floor(rng (&RNG_SEED) * n1Avail);
	const int a2 = (int) floor(rng (&RNG_SEED) * n2Avail);
	atom a1_orig = sys.atoms[typeIndex_][a1], a1_new = a1_orig;
	atom a2_orig = sys.atoms[typeIndex2_][a2], a2_new = a2_orig;

	// Positions will be exchanged, but no other property should change
	a1_new.pos = a2_orig.pos;
	a2_new.pos = a1_orig.pos;

	const std::vector < double > box = sys.box();
	double V = 1.0;
	for (unsigned int i = 0; i < box.size(); ++i) {
        V *= box[i];
    }

	double dU = 0.0;

    double delEnergy = 0.0;
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    	// Get positions of neighboring atoms around a1
        std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, &sys.atoms[typeIndex_][a1]);
        for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
        	if (neighborAtoms[i] == &sys.atoms[typeIndex2_][a2]) {
        		// Skip their interaction as they were already deleted - this is fine because their net pairwise interaction doesn't change over the course of the simulation - this means we don't have to actually re-insert the other atom / get its interation somehow later
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
        // Add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
        	if (!(sys.getCurrentM() > 0 && sys.getFractionalAtom () == &sys.atoms[typeIndex_][a1])) {
        		// Then a1 is not the partially inserted particle and tail interactions must be included
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

    // Account for wall interaction energy
    delEnergy += sys.speciesBarriers[typeIndex_].energy(&sys.atoms[typeIndex_][a1], box);

    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    	// Get positions of neighboring atoms around a2
        std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex2_, &sys.atoms[typeIndex2_][a2]);
        for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
        	if (neighborAtoms[i] == &sys.atoms[typeIndex_][a1]) {
        		// Skip their interaction as they were already deleted - this is fine because their net pairwise interaction doesn't change over the course of the simulation - this means we don't have to actually re-insert the other atom / get its interation somehow later
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
        // Add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        if (sys.ppot[spec][typeIndex2_]->useTailCorrection) {
        	if (!(sys.getCurrentM() > 0 && sys.getFractionalAtom () == &sys.atoms[typeIndex2_][a2])) {
        		// Then a2 is not the partially inserted particle and tail interactions must be included
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

    // Account for wall interaction energy
    delEnergy += sys.speciesBarriers[typeIndex2_].energy(&sys.atoms[typeIndex2_][a2], box);

    double insEnergy = 0.0;

	// Account for wall interaction energy first to be more efficient
	dU = sys.speciesBarriers[typeIndex_].energy(&a1_new, box);
	if (dU < NUM_INFINITY) {
		insEnergy += dU;
	} else {
		insEnergy = NUM_INFINITY;
	}

	if (insEnergy < NUM_INFINITY) {
		// Account for wall interaction energy
	    dU = sys.speciesBarriers[typeIndex2_].energy(&a2_new, box);
		if (dU < NUM_INFINITY) {
			insEnergy += dU;
		} else {
			insEnergy = NUM_INFINITY;
		}
	}

	if (insEnergy < NUM_INFINITY) {
		for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
	    	// Get positions of neighboring atoms around a1's (a2's) new (old) location
	    	std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, &a1_new);
	    	for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
	    		// With these new "copy atoms" getNeighborAtoms can't guarantee it doesn't point to old self, so must check
	    		if ((neighborAtoms[i] == &sys.atoms[typeIndex2_][a2]) || (neighborAtoms[i] == &sys.atoms[typeIndex_][a1])) {
	    			continue;
	    		} else {
					try {
						dU = sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], &a1_new, box);
					} catch (customException& ce) {
						std::string a = "Cannot insert because of energy error: ", b = ce.what();
						throw customException (a+b);
					}
					if (dU < NUM_INFINITY) {
						insEnergy += dU;
					} else {
						insEnergy = NUM_INFINITY;
						break;
					}
	    		}
	    	}
			if (insEnergy == NUM_INFINITY) break; // Don't add anything if "infinite" already

	    	// Add tail correction to potential energy -- only enable for fluid phase simulations
	#ifdef FLUID_PHASE_SIMULATIONS
	    	if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
	        	if (!(sys.getCurrentM() > 0 && sys.getFractionalAtom () == &sys.atoms[typeIndex_][a1])) {
	        		// Then a1 is not the partially inserted particle and tail interactions must be included
	        		if (spec == typeIndex_) {
	                	if (sys.numSpecies[spec]-1 > 0) {
	                    	insEnergy += sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec]-1)/V); // Never infinite
						}
	        		} else {
	                	if (sys.numSpecies[spec] > 0) {
	                		insEnergy += sys.ppot[spec][typeIndex_]->tailCorrection(sys.numSpecies[spec]/V); // Never infinite
	        			}
	    			}
				}
	    	}
	#endif
	    }
	}

	if (insEnergy < NUM_INFINITY) {
		for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
	    	// Get positions of neighboring atoms around a2's (a1's) new (old) location
	    	std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex2_, &a2_new);
	        for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
	        	if ((neighborAtoms[i] == &sys.atoms[typeIndex_][a1]) || (neighborAtoms[i] == &sys.atoms[typeIndex2_][a2])) {
	        		continue;
	        	} else {
					try {
						dU = sys.ppot[spec][typeIndex2_]->energy(neighborAtoms[i], &a2_new, box);
					} catch (customException& ce) {
						std::string a = "Cannot insert because of energy error: ", b = ce.what();
						throw customException (a+b);
					}
					if (dU < NUM_INFINITY) {
						insEnergy += dU;
					} else {
						insEnergy = NUM_INFINITY;
						break;
					}
	        	}
	        }
			if (insEnergy == NUM_INFINITY) break; // Don't add anything if "infinite" already

	        // Add tail correction to potential energy -- only enable for fluid phase simulations
	#ifdef FLUID_PHASE_SIMULATIONS
	        if (sys.ppot[spec][typeIndex2_]->useTailCorrection) {
	        	if (!(sys.getCurrentM() > 0 && sys.getFractionalAtom () == &sys.atoms[typeIndex2_][a2])) {
	        		// Then a2 is not the partially inserted particle and tail interactions must be included
	        		if (spec == typeIndex2_) {
	                	if (sys.numSpecies[spec]-1 > 0) {
	                        insEnergy += sys.ppot[spec][typeIndex2_]->tailCorrection((sys.numSpecies[spec]-1)/V); // Never infinite
	                	}
	        		} else {
	                    if (sys.numSpecies[spec] > 0) {
	                        insEnergy += sys.ppot[spec][typeIndex2_]->tailCorrection(sys.numSpecies[spec]/V); // Never infinite
	                    }
	        		}
	        	}
	        }
	#endif
	    }
	}

	// Biasing
    double p_u = 0.0;
	if (insEnergy < NUM_INFINITY) {
		p_u = exp(-sys.beta()*(insEnergy - delEnergy));
	}
    double bias = calculateBias(sys, sys.getTotN(), sys.getCurrentM());

    // TMMC gets updated the same way, regardless of whether the move gets accepted
    if (sys.useTMMC) {
    	sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), std::min(1.0, p_u));
    }

	if (rng (&RNG_SEED) < p_u*bias) { // Swap the particles by deleting/reinserting
		sys.incrementEnergy(insEnergy - delEnergy);

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
