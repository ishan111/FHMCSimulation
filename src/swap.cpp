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
	if (sys.numSpecies[typeIndex_] < 1 || sys.numSpecies[typeIndex2_] < 1) {
		return MOVE_FAILURE;
	}
	
	// because the locations are effectively being swapped, it is fair to select a partially inserted atom to be involved
	const int a1 = (int) floor(rng (&RNG_SEED) * sys.numSpecies[typeIndex_]);
	const int a2 = (int) floor(rng (&RNG_SEED) * sys.numSpecies[typeIndex2_]);
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
	    
	// Swapping = COMPLETELY delete both, one at a time, then add both at opposite positions, one at a time
	// Delete a1 - take it from some, potentially partially inserted state, and entirely remove it
    double delEnergy = 0.0;
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
        // get positions of neighboring atoms around a1
        std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, &sys.atoms[typeIndex_][a1]);
        for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
            try {
				delEnergy += sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], &sys.atoms[typeIndex_][a1], box);
			}
			catch (customException& ce) {
				std::string a = "Cannot delete because of energy error: ", b = ce.what();
				throw customException (a+b);
			}
        }
        // add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
        	if (spec == typeIndex_) {
				delEnergy += sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec]-1)/V);
			}
			else {
				delEnergy += sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec])/V);
			}
		}
#endif
    }
    
    // Remove a1 from cell lists, etc.
    try {
    	sys.deleteAtom(typeIndex_, a1, true);
    } catch (customException &ce) {
    	std::string a = "Failed to delete atom during swapping: ", b = ce.what();
    	throw customException (a+b);
    }
    
	// Delete a2
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
        // get positions of neighboring atoms around a2
        std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex2_, &sys.atoms[typeIndex2_][a2]);
        for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
            try {
				delEnergy += sys.ppot[spec][typeIndex2_]->energy(neighborAtoms[i], &sys.atoms[typeIndex2_][a2], box);
			}
			catch (customException& ce) {
				std::string a = "Cannot delete because of energy error: ", b = ce.what();
				throw customException (a+b);
			}
        }
        // add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        if (sys.ppot[spec][typeIndex2_]->useTailCorrection) {
        	if (spec == typeIndex2_) {
				delEnergy += sys.ppot[spec][typeIndex2_]->tailCorrection((sys.numSpecies[spec]-1)/V);
			}
			else {
				delEnergy += sys.ppot[spec][typeIndex2_]->tailCorrection((sys.numSpecies[spec])/V);
			}
		}
#endif
    }

    // Remove a2 from cell lists, etc.
    try {
    	sys.deleteAtom(typeIndex2_, a2, true);
    } catch (customException &ce) {
    	std::string a = "Failed to delete atom during swapping: ", b = ce.what();
    	throw customException (a+b);
    }

	// Insert new_a1
    double insEnergy = 0.0;
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    	// get positions of neighboring atoms around a1's (a2's) new (old) location
    	std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, &a1_new);
        for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
			try {
				insEnergy += sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], &a1_new, box);	
			} catch (customException& ce) {
				std::string a = "Cannot insert because of energy error: ", b = ce.what();
				throw customException (a+b);
			}
        }
        // add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        if (sys.ppot[spec][typeIndex_]->useTailCorrection){
			insEnergy += sys.ppot[spec][typeIndex_]->tailCorrection(sys.numSpecies[spec]/V);
		}
#endif
    }
    
    // Add a1 into the system in a2's original place
    try {
    	sys.insertAtom(typeIndex_, &a1_new, true); 
    } catch (customException &ce) {
    	std::string a = "Failed to insert atom during swapping: ", b = ce.what();
    	throw customException (a+b);
    }

    // Inster new_a2
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    	std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex2_, &a2_new);
        for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
			try {
				insEnergy += sys.ppot[spec][typeIndex2_]->energy(neighborAtoms[i], &a2_new, box);	
			} catch (customException& ce) {
				std::string a = "Cannot insert because of energy error: ", b = ce.what();
				throw customException (a+b);
			}
        }
        // add tail correction to potential energy -- only enable for fluid phase simulations
#ifdef FLUID_PHASE_SIMULATIONS
        if (sys.ppot[spec][typeIndex2_]->useTailCorrection){
			insEnergy += sys.ppot[spec][typeIndex2_]->tailCorrection(sys.numSpecies[spec]/V);
		}
#endif
    }

    // Add a2 to the system
    try {
    	sys.insertAtom(typeIndex2_, &a2_new, true);
    } catch (customException &ce) {
    	std::string a = "Failed to insert atom during swapping: ", b = ce.what();
    	throw customException (a+b);
    }
      
    // Biasing
    const double p_u = exp(-sys.beta()*(insEnergy - delEnergy));
    double bias = calculateBias(sys, sys.getTotN()); // sys.numSpecies already contains the currently proposed modifications
    
    // tmmc gets updated the same way, regardless of whether the move gets accepted
    if (sys.useTMMC) {
    	sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), std::min(1.0, p_u)); // since the total number of atoms isn't changing, can use getTotN() as both initial and final states
    }
    
	if (rng (&RNG_SEED) < p_u*bias) {
	   sys.incrementEnergy(insEnergy - delEnergy);	
		
		// update Wang-Landau bias, if used
		if (sys.useWALA) {
			sys.getWALABias()->update(sys.getTotN());
		}
			
        return MOVE_SUCCESS;
    }
	
	// Undo to the swap if move was rejected
	// The "new" atoms are now at the ends of each of the vector for each type
	try {
    	sys.deleteAtom(typeIndex_, sys.numSpecies[typeIndex_]-1);
    } catch (customException &ce) {
    	std::string a = "Failed to delete atom during swapping: ", b = ce.what();
    	throw customException (a+b);
    }
    try {
        sys.deleteAtom(typeIndex2_, sys.numSpecies[typeIndex2_]-1);
    } catch (customException &ce) {
       	std::string a = "Failed to delete atom during swapping: ", b = ce.what();
        throw customException (a+b);
    }
    try {
    	sys.insertAtom(typeIndex_, &a1_orig);
    } catch (customException &ce) {
        std::string a = "Failed to insert atom during swapping: ", b = ce.what();
        throw customException (a+b);
    }
    try {
    	sys.insertAtom(typeIndex2_, &a2_orig);
    } catch (customException &ce) {
        std::string a = "Failed to insert atom during swapping: ", b = ce.what();
        throw customException (a+b);
    }
        
	// update Wang-Landau bias (even if moved failed), if used
	if (sys.useWALA) {
		sys.getWALABias()->update(sys.getTotN());
	}
	
	return MOVE_FAILURE;
}
