#include "translate.h"

/*!
 * Translate a particle in the system.  All other information is stored in the simSystem object.
 *
 * \param [in] sys System object to attempt to remove a particle from.
 *
 * \return MOVE_SUCCESS if translated a particle, otherwise MOVE_FAILURE if did not.  Will throw exceptions if there was an error.
 */
int translateParticle::make (simSystem &sys) {
	bool earlyReject = false;

	// Check if any exist to be translated
	if (sys.getFractionalAtomType() == typeIndex_) {
		if (sys.numSpecies[typeIndex_] == 0 && sys.getCurrentM() == 0) {
        	earlyReject = true;
    	}
	} else {
		if (sys.numSpecies[typeIndex_] == 0) {
			earlyReject = true;
		}
	}

	// Updates to biasing functions must be done even if at bounds
	if (earlyReject) {
        if (sys.useWALA) {
        	sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
        }
        if (sys.useTMMC) {
        	sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), 0.0);
        }
        return MOVE_FAILURE;
    }

	// Choose a random particle of that type
	int chosenAtom = 0;
	if (sys.getCurrentM() > 0 && sys.getFractionalAtomType() == typeIndex_) {
		// We are moving a species that has a partially inserted atom, so account for that in the choice
		chosenAtom = (int) floor(rng (&RNG_SEED) * (sys.numSpecies[typeIndex_]+1));
	} else {
		// All atoms of this type are fully inserted
		chosenAtom = (int) floor(rng (&RNG_SEED) * sys.numSpecies[typeIndex_]);
	}

	// Attempt to translate that one
	const std::vector < double > box = sys.box();
    double V = 1.0;
    for (unsigned int i = 0; i < box.size(); ++i) {
    	V *= box[i];
    }

    double oldEnergy = 0.0;
    for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
    	// Get positions of neighboring atoms around chosenAtom
        std::vector < atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, &sys.atoms[typeIndex_][chosenAtom]);
    	for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
			try {
				oldEnergy += sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], &sys.atoms[typeIndex_][chosenAtom], box);
			} catch (customException& ce) {
				std::string a = "Cannot translate because of energy error: ", b = ce.what();
				throw customException (a+b);
			}
        }

        // Add tail correction to potential energy
#ifdef FLUID_PHASE_SIMULATIONS
        if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
        	if (!(sys.getCurrentM() > 0 && sys.getFractionalAtom () == &sys.atoms[typeIndex_][chosenAtom])) {
        		// Then chosenAtom is not a partially inserted particle and tail interactions must be included
        		if (spec == typeIndex_) {
                    if (sys.numSpecies[spec]-1 > 0) {
                    	oldEnergy += sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec]-1)/V);
                	}
        		} else {
                    if (sys.numSpecies[spec] > 0) {
                    	oldEnergy += sys.ppot[spec][typeIndex_]->tailCorrection(sys.numSpecies[spec]/V);
                    }
        		}
        	}
        }
#endif
    }

    // Account for wall or barrier interactions
    oldEnergy += sys.speciesBarriers[typeIndex_].energy(&sys.atoms[typeIndex_][chosenAtom], box);

    // Store old position and move particle along random direction in interval [-maxD_:maxD_]
    std::vector<double> oldPos = sys.atoms[typeIndex_][chosenAtom].pos;
    for (unsigned int i = 0; i< sys.atoms[typeIndex_][chosenAtom].pos.size(); ++i) {
    	sys.atoms[typeIndex_][chosenAtom].pos[i] += 2.0*maxD_*(0.5-rng (&RNG_SEED));

    	// Apply periodic boundary conditions
    	if (sys.atoms[typeIndex_][chosenAtom].pos[i] >= box[i]) {
    		sys.atoms[typeIndex_][chosenAtom].pos[i] -= box[i];
	    } else if (sys.atoms[typeIndex_][chosenAtom].pos[i] < 0) {
    		sys.atoms[typeIndex_][chosenAtom].pos[i] += box[i];
    	}
    }

	// Calculate energy at new position
    double newEnergy = 0.0;
	double dU = 0.0;

	// Account for wall or barrier interactions first to be more efficient
	dU = sys.speciesBarriers[typeIndex_].energy(&sys.atoms[typeIndex_][chosenAtom], box);
	if (dU < NUM_INFINITY) {
		newEnergy += dU;
	} else {
		newEnergy = NUM_INFINITY;
	}

	if (newEnergy < NUM_INFINITY) {
		for (unsigned int spec = 0; spec < sys.nSpecies(); ++spec) {
	    	// Get positions of neighboring atoms around chosenAtom
	    	std::vector< atom* > neighborAtoms = sys.getNeighborAtoms(spec, typeIndex_, &sys.atoms[typeIndex_][chosenAtom]);
	    	for (unsigned int i = 0; i < neighborAtoms.size(); ++i) {
				try {
					dU = sys.ppot[spec][typeIndex_]->energy(neighborAtoms[i], &sys.atoms[typeIndex_][chosenAtom], box);
				} catch (customException& ce) {
					std::string a = "Cannot delete because of energy error: ", b = ce.what();
					throw customException (a+b);
				}
				if (dU < NUM_INFINITY) {
					newEnergy += dU;
				} else {
					newEnergy = NUM_INFINITY;
					break;
				}
	        }
			if (newEnergy == NUM_INFINITY) break; // Don't add anything if "infinite" already

	        // Add tail correction to potential energy
	#ifdef FLUID_PHASE_SIMULATIONS
		    if (sys.ppot[spec][typeIndex_]->useTailCorrection) {
	        	if (!(sys.getCurrentM() > 0 && sys.getFractionalAtom () == &sys.atoms[typeIndex_][chosenAtom])) {
	        		// Then chosenAtom is not a partially inserted particle and tail interactions must be included
	        		if (spec == typeIndex_) {
	                    if (sys.numSpecies[spec]-1 > 0) {
		                    newEnergy += sys.ppot[spec][typeIndex_]->tailCorrection((sys.numSpecies[spec]-1)/V); // Never infinite
	        	        }
	        		} else {
	                	if (sys.numSpecies[spec] > 0) {
	                        newEnergy += sys.ppot[spec][typeIndex_]->tailCorrection(sys.numSpecies[spec]/V); // Never infinite
		                }
	        		}
	        	}
	        }
	#endif
	    }
	}

	// Biasing
	double p_u = 0.0;
	if (newEnergy < NUM_INFINITY) {
		p_u = exp(-sys.beta()*(newEnergy - oldEnergy));
	}
	double bias = calculateBias(sys, sys.getTotN(), sys.getCurrentM()); // N_tot doesn't change throughout this move

	// TMMC gets updated the same way, regardless of whether the move gets accepted
    if (sys.useTMMC) {
    	sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), std::min(1.0, p_u)); // Since the total number of atoms isn't changing, can use getTotN() as both initial and final states
    }

	if (rng (&RNG_SEED) < p_u*bias) {
	    try {
        	sys.translateAtom(typeIndex_, chosenAtom, oldPos);
	    } catch (customException &ce) {
        	std::string a = "Failed to translate atom: ", b = ce.what();
            throw customException (a+b);
        }
		sys.incrementEnergy(newEnergy - oldEnergy);

		// Update Wang-Landau bias, if used
		if (sys.useWALA) {
			sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
		}
	    return MOVE_SUCCESS;
    }

   	// If move failed, reset position
    for (unsigned int i = 0; i < sys.atoms[typeIndex_][chosenAtom].pos.size(); ++i) {
    	sys.atoms[typeIndex_][chosenAtom].pos[i] = oldPos[i];
    }

    // Update Wang-Landau bias (even if moved failed), if used
   	if (sys.useWALA) {
   		sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
    }
	return MOVE_FAILURE;
}

/*!
 * Set the maximum translation in any single move. Should be postive number lss than half the box size.
 *
 * \param [in] maxD Maximium translation
 * \param [in] box Box dimensions
 */
void translateParticle::setMaxTranslation (const double maxD, const std::vector < double > &box) {
	for (unsigned int i = 0; i < box.size(); ++i) {
		if (maxD >= box[i]/2.) {
			throw customException ("Max translation too large");
		}
	}
	if (maxD > 0) {
		maxD_ = maxD;
	} else {
		throw customException ("Max translation must be positive");
	}
}
