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
	const int a1 = (int) floor(rng (&RNG_SEED) * sys.numSpecies[typeIndex_]);
	const int a2 = (int) floor(rng (&RNG_SEED) * sys.numSpecies[typeIndex2_]);
	atom a1_orig = sys.atoms[typeIndex_][a1];
	atom a2_orig = sys.atoms[typeIndex2_][a2];

	const std::vector < double > box = sys.box();
	double V = 1.0;
	for (unsigned int i = 0; i < box.size(); ++i) {
        V *= box[i];
    }
	    
	// Swapping = delete both, one at a time, then add both at opposite positions, one at a time
	// Delete a1
    double delEnergy = 0.0;
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
        // get positions of neighboring atoms around a1
        std::vector < std::vector< double > > neighborPositions = sys.getNeighborPositions(spec, typeIndex_, &sys.atoms[typeIndex_][a1]);
        for (unsigned int i = 0; i < neighborPositions.size(); ++i) {
            try {
				delEnergy += sys.ppot[spec][typeIndex_]->energy(neighborPositions[i], sys.atoms[typeIndex_][a1].pos, box);
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
        std::vector < std::vector< double > > neighborPositions = sys.getNeighborPositions(spec, typeIndex2_, &sys.atoms[typeIndex2_][a2]);
        for (unsigned int i = 0; i < neighborPositions.size(); ++i) {
            try {
				delEnergy += sys.ppot[spec][typeIndex2_]->energy(neighborPositions[i], sys.atoms[typeIndex2_][a2].pos, box);
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

	// Insert a1 at a2's original location
    double insEnergy = 0.0;
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    	// get positions of neighboring atoms around a1's (a2's) new (old) location
    	std::vector < std::vector < double > > neighborPositions = sys.getNeighborPositions(spec, typeIndex_, &a2_orig);
        for (unsigned int i = 0; i < neighborPositions.size(); ++i) {
			try {
				insEnergy += sys.ppot[spec][typeIndex_]->energy(neighborPositions[i], a2_orig.pos, box);	
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
    
    // Add a1 to the system
    try {
    	sys.insertAtom(typeIndex_, &a2_orig);
    } catch (customException &ce) {
    	std::string a = "Failed to insert atom during swapping: ", b = ce.what();
    	throw customException (a+b);
    }

	// Insert a2 at a1's original location
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    	// get positions of neighboring atoms around a2's (a1's) new (old) location
    	std::vector < std::vector < double > > neighborPositions = sys.getNeighborPositions(spec, typeIndex2_, &a1_orig);
        for (unsigned int i = 0; i < neighborPositions.size(); ++i) {
			try {
				insEnergy += sys.ppot[spec][typeIndex2_]->energy(neighborPositions[i], a1_orig.pos, box);	
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
    	sys.insertAtom(typeIndex2_, &a1_orig);
    } catch (customException &ce) {
    	std::string a = "Failed to insert atom during swapping: ", b = ce.what();
    	throw customException (a+b);
    }
       
    // Biasing
    const double p_u = exp(-sys.beta()*(insEnergy - delEnergy));
    double bias = calculateBias(sys, sys.getTotN(), p_u); // sys.numSpecies already contains the currently proposed modifications
    
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
