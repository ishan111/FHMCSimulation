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
    if (sys.numSpecies[typeIndex_] <= sys.minSpecies(typeIndex_)) {
        return MOVE_FAILURE;
    }
    // also check if at global bound on total number of particles
    if (sys.getTotN() <= sys.totNMin()) {
       	return MOVE_FAILURE;
    }
    
	// choose a random particle (index) of that type
	const int chosenAtom = (int) floor(rng (&RNG_SEED) * sys.numSpecies[typeIndex_]);

	// attempt to delete that one
	const std::vector < double > box = sys.box();
    double V = 1.0;
    for (unsigned int i = 0; i < box.size(); ++i) {
        V *= box[i];
    }
        
    double delEnergy = 0.0;
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
        // get positions of neighboring atoms around chosenAtom
        std::vector < std::vector< double > > neighborPositions = sys.getNeighborPositions(spec, typeIndex_, &sys.atoms[typeIndex_][chosenAtom]);
        for (unsigned int i = 0; i < neighborPositions.size(); ++i) {
            try {
				delEnergy -= sys.ppot[spec][typeIndex_]->energy(neighborPositions[i], sys.atoms[typeIndex_][chosenAtom].pos, box);
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
				delEnergy -= sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec]-1)/V);
			}
			else {
				delEnergy -= sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec])/V);
			}
		}
#endif
    }
    
    // biasing
    const double p_u = sys.numSpecies[typeIndex_]/V*exp(sys.beta()*(-sys.mu(typeIndex_) - delEnergy));
    int nTotFinal = sys.getTotN() - 1;
    double bias = calculateBias(sys, nTotFinal, p_u);
    
	// metropolis criterion
	if (rng (&RNG_SEED) < p_u*bias) {
	    try {
            sys.deleteAtom(typeIndex_, chosenAtom);
        } catch (customException &ce) {
            std::string a = "Failed to delete atom: ", b = ce.what();
            throw customException (a+b);
        }
		sys.incrementEnergy(delEnergy);	
		
		// update Wang-Landau bias, if used
		if (sys.useWALA) {
			sys.getWALABias()->update(sys.getTotN());
		}
				
        return MOVE_SUCCESS;
    }
    
	// update Wang-Landau bias (even if moved failed), if used
	if (sys.useWALA) {
		sys.getWALABias()->update(sys.getTotN());
	}
		
	return MOVE_FAILURE;
}
