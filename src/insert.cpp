#include "insert.h"

/*!
 * Insert a particle into the system.  All other information is stored in the simSystem object.
 * 
 * \param [in] sys System object to attempt to insert a particle into.
 * 
 * \return MOVE_SUCCESS if inserted a particle, otherwise MOVE_FAILURE if did not.  Will throw exceptions if there was an error.
 */
int insertParticle::make (simSystem &sys) {
    // check if at upper bound
    if (sys.numSpecies[typeIndex_] >= sys.maxSpecies(typeIndex_)) {
        return MOVE_FAILURE;
    }
    
	// attempt to insert a new one
    atom newAtom;
    const std::vector < double > box = sys.box();
    double V = 1.0;
    for (unsigned int i = 0; i < box.size(); ++i) {
        newAtom.pos[i] = rng (&RNG_SEED) * box[i];
        V *= box[i];
    }
    
    // compute energy to insert
    double insEnergy = 0.0;
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    	// get positions of neighboring atoms around newAtom
    	std::vector < std::vector < double > > neighborPositions = sys.getNeighborPositions(spec, typeIndex_, &newAtom);
        for (unsigned int i = 0; i < neighborPositions.size(); ++i) {
			try {
				insEnergy += sys.ppot[spec][typeIndex_]->energy(neighborPositions[i], newAtom.pos, box);	
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
    
    // biasing
    const double p_u = V/(sys.numSpecies[typeIndex_]+1.0)*exp(sys.beta()*(sys.mu(typeIndex_) - insEnergy));
    std::vector <int> Nend = sys.numSpecies;
    Nend[typeIndex_] += 1;
    double bias = calculateBias(sys, Nend, p_u);
   
	// metropolis criterion
	if (rng (&RNG_SEED) < p_u*bias) {
        try {
            sys.insertAtom(typeIndex_, &newAtom);
        } catch (customException &ce) {
            std::string a = "Failed to insert atom: ", b = ce.what();
            throw customException (a+b);
        }
		sys.incrementEnergy(insEnergy);
		
		// update Wang-Landau bias, if used
		if (sys.useWALA) {
			sys.getWALABias()->update(sys.numSpecies);
		}
		
        return MOVE_SUCCESS;
    }
    
	// update Wang-Landau bias (even if moved failed), if used
	if (sys.useWALA) {
		sys.getWALABias()->update(sys.numSpecies);
	}
	
	return MOVE_FAILURE;
}
