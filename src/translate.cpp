#include "system.h"
#include "translate.h"
#include "global.h"
#include "moves.h"
#include "atom.h"
#include <iostream>
#include <sstream>
#include <string>

/*!
 * translate a particle in the system.  All other information is stored in the simSystem object.
 * \param [in] sys System object to attempt to remove a particle from.
 * \return MOVE_SUCCESS if translated a particle, otherwise MOVE_FAILURE if did not.  Will throw exceptions if there was an error.
 */
int translateParticle::make (simSystem &sys) {
	// check if any exist to be translated
    if (sys.numSpecies[typeIndex_] < 1) {
        return MOVE_FAILURE;
    }
    
	// choose a random particle of that type
	const int chosenAtom = (int) floor(rng (&RNG_SEED) * sys.numSpecies[typeIndex_]);
 
	// attempt to translate that one
	const std::vector < double > box = sys.box();
    double V = 1.0;
    for (unsigned int i = 0; i < box.size(); ++i) {
        V *= box[i];
    }
        
    double oldEnergy = 0.0;
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
        // get positions of neighboring atoms around chosenAtom
        std::vector< std::vector<double> > neighborPositions = sys.getNeighborPositions(spec, typeIndex_, &sys.atoms[typeIndex_][chosenAtom]);
    	for (unsigned int i = 0; i < neighborPositions.size(); ++i) {
			try {
				oldEnergy += sys.ppot[spec][typeIndex_]->energy(neighborPositions[i], sys.atoms[typeIndex_][chosenAtom].pos, box);
			}
			catch (customException& ce) {
				std::string a = "Cannot translate because of energy error: ", b = ce.what();
				throw customException (a+b);
			}
        }
        // add tail correction to potential energy
        if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
			oldEnergy += sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec])/V);
		}
    }
    
    // store old position and move particle along random direction in interval [-0.1:0.1]
    //std::cout<<"Moving: "<<chosenAtom<<": "<<std::endl;
    std::vector<double> oldPos = sys.atoms[typeIndex_][chosenAtom].pos;
    for (unsigned int i=0; i<sys.atoms[typeIndex_][chosenAtom].pos.size(); i++) {
    	sys.atoms[typeIndex_][chosenAtom].pos[i] += 0.2*(0.5-rng (&RNG_SEED));
    	
    	// apply periodic boundary conditions
    	if (sys.atoms[typeIndex_][chosenAtom].pos[i] >= box[i]) {
    		sys.atoms[typeIndex_][chosenAtom].pos[i] -= box[i];
    	}
    	else if (sys.atoms[typeIndex_][chosenAtom].pos[i] < 0) {
    		sys.atoms[typeIndex_][chosenAtom].pos[i] += box[i];
    	}
    }
    
    // calculate energy at new position
    double newEnergy = 0.0;
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
        // get positions of neighboring atoms around chosenAtom
        std::vector< std::vector<double> > neighborPositions = sys.getNeighborPositions(spec, typeIndex_, &sys.atoms[typeIndex_][chosenAtom]);
    	for (unsigned int i = 0; i < neighborPositions.size(); ++i) {
			try {
				newEnergy += sys.ppot[spec][typeIndex_]->energy(neighborPositions[i], sys.atoms[typeIndex_][chosenAtom].pos, box);
			}
			catch (customException& ce) {
				std::string a = "Cannot delete because of energy error: ", b = ce.what();
				throw customException (a+b);
			}
        }
        // add tail correction to potential energy
        if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
			newEnergy += sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec])/V);
		}
    }
    
	// metropolis criterion
	double dEnergy = newEnergy - oldEnergy;
	
	if (rng (&RNG_SEED) < exp(-sys.beta()*dEnergy)) {
	    try {
            sys.translateAtom(typeIndex_, chosenAtom, oldPos);
        } catch (customException &ce) {
            std::string a = "Failed to translate atom: ", b = ce.what();
            throw customException (a+b);
        }
		sys.incrementEnergy(dEnergy);	
        return MOVE_SUCCESS;
    }
    
    // if move failed, reset position
    for (unsigned int i=0; i<sys.atoms[typeIndex_][chosenAtom].pos.size(); i++) {
    	sys.atoms[typeIndex_][chosenAtom].pos[i] = oldPos[i];
    }
    
	return MOVE_FAILURE;
}
