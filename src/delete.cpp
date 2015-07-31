#include "system.h"
#include "delete.h"
#include "global.h"
#include "moves.h"
#include "atom.h"
#include <iostream>
#include <sstream>
#include <string>

/*!
 * Delete a particle from the system.  All other information is stored in the simSystem object.
 * \param [in] sys System object to attempt to remove a particle from.
 * \return MOVE_SUCCESS if deleted a particle, otherwise MOVE_FAILURE if did not.  Will throw exceptions if there was an error.
 */
int deleteParticle::make (simSystem &sys) {
	// check if any exist to be deleted
    if (sys.numSpecies[typeIndex_] < 1) {
        //std::string itype = static_cast<std::ostringstream*>( &(std::ostringstream() << typeIndex_) )->str();
        //std::cerr << "Hit lower bound for species type "+itype << std::endl;
        return MOVE_FAILURE;
    }
    
	// choose a random particle of that type
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
      std::vector< std::vector<double> > neighborPositions = sys.getNeighborPositions(spec, typeIndex_, &sys.atoms[typeIndex_][chosenAtom]);
    	 for (unsigned int i = 0; i < neighborPositions.size(); ++i) {
			try {
				delEnergy -= sys.ppot[spec][typeIndex_]->energy(neighborPositions[i], sys.atoms[typeIndex_][chosenAtom].pos, box);
			}
			catch (customException& ce) {
				std::string a = "Cannot delete because of energy error: ", b = ce.what();
				throw customException (a+b);
			}
        }
        // add tail correction to potential energy
        if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
        	if (spec == typeIndex_) {
				delEnergy -= sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec]-1)/V);
			}
			else {
				delEnergy -= sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec])/V);
			}
		}
    }
    
	// metropolis criterion
	if (rng (&RNG_SEED) < sys.numSpecies[typeIndex_]/V*exp(sys.beta()*(-sys.mu(typeIndex_) - delEnergy))) {
	    try {
            sys.deleteAtom(typeIndex_, chosenAtom);
        } catch (customException &ce) {
            std::string a = "Failed to delete atom: ", b = ce.what();
            throw customException (a+b);
        }
		sys.incrementEnergy(delEnergy);	
        return MOVE_SUCCESS;
    }
    
	return MOVE_FAILURE;
}
