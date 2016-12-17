#include "system.h"

/*!
 * Increase the expanded ensemble state of the system by 1.  Accounts for the periodicity of [0, M)
 */
void simSystem::incrementMState () {
	Mcurrent_++;
	if (Mcurrent_ == Mtot_) {
		Mcurrent_ = 0;
	}
}

/*!
 * Decrease the expanded ensemble state of the system by 1.  Accounts for the periodicity of [0, M)
 */
void simSystem::decrementMState () {
	Mcurrent_--;
	if (Mcurrent_ < 0) {
		Mcurrent_ += Mtot_;
	}
}

/*!
 * Set the bounds on the total number of particles in a system.  If not set manually, this defaults to the sum of the bounds given for
 * each individual species in the system.  Therefore, for single component simulations, this is identical to [minSpecies(0), maxSpecies(0)]
 * unless otherwise set.  These bounds are intended to be used to create "windows" so that specific simulations can sample subregions
 * of [minSpecies(0), maxSpecies(0)] and be stitched together with histogram reweighting later.
 *
 * However, this routine will ALSO cause the system to reevaluate its bounds.  If these total bounds are outside any individual
 * bound for each atom type, nothing will change.  However, if the upper bound for total atoms is less than an upper bound for
 * a specific species, that species will have its bounds changed to match the total maximum.  As a result sys.atoms can change so this
 * routine should be called at the beginning of a simulation, never during. The total minimum will also be checked.
 * That is, if the sum of the minimum for all species is still higher than this, an exception will be throw since the system will never
 * reach such a low density anyway.  Most likely the user has made a mistake.
 *
 * Be sure to initialize other objects, such as biases, AFTER this routine has been called since it will adjust the allowable number of
 * particles in the system.
 *
 * \param [in] bounds Vector of [min, max]
 */
void simSystem::setTotNBounds (const std::vector < int > &bounds) {
	if (bounds.size() != 2) {
		throw customException ("Bounds on total N must supplied as vector of <minN, maxN>");
	}
	if (bounds[0] < 0) {
		throw customException ("Lower bound on total particles must be > 0");
	}
	if (bounds[0] > bounds[1]) {
		throw customException ("Upper bound must be greater than lower bound for total number of particles in the system");
	}
	totNBounds_ = bounds;

	int totMin = 0;
	for (unsigned int i = 0; i < nSpecies_; ++i) {
		if (maxSpecies_[i] > totNBounds_[1]) {
			maxSpecies_[i] = totNBounds_[1];
		}
		totMin += minSpecies_[i];
	}
	if (totMin > totNBounds_[0]) {
		// this isn't the end of the world, but for now, alert the user in case something is wrong
		throw customException ("Lower total N bound is lower than the sum of all individual lower bounds, region cannot be completely sampled");
	}

	// recheck bounds and possibly resize
	int tmpTot = 0;
    	for (unsigned int i = 0; i < nSpecies_; ++i) {
        	if (maxSpecies_[i] < minSpecies_[i]) {
            		throw customException ("Max species < Min species");
        	}
		try {
			atoms[i].resize(maxSpecies_[i]);
		} catch (std::exception &e) {
			throw customException (e.what());
		}
		if (numSpecies[i] > atoms[i].size()) {
			numSpecies[i] = atoms.size();
		}
		tmpTot += numSpecies[i];
	}
    	totN_ = tmpTot;

    	// Allocate space for energy matrix - this will only be recorded when the system is within the specific window we are looking for
    	// Because of implementation of Shen and Errington method, this syntax is the same for single and multicomponent systems
    	long long int size = totNBounds_[1] - totNBounds_[0] + 1;
    	/*try {
    		AverageU_.resize(size, 0);
    	} catch (std::bad_alloc &ba) {
    		throw customException ("Out of memory for energy record");
    	}
	try {
    		averageU2_.resize(size, 0);
    	} catch (std::bad_alloc &ba) {
    		throw customException ("Out of memory for energy^2 record");
    	}
    	try {
    		numAverageU_.resize(size, 0);
    	} catch (std::bad_alloc &ba) {
        	throw customException ("Out of memory for energy record");
    	}
	try {
    		numAverageFlucts_.resize(size, 0);
    	} catch (std::bad_alloc &ba) {
        	throw customException ("Out of memory for fluctuations record");
    	}
	for (unsigned int i = 0; i < averageN_.size(); ++i) {
		try {
        		averageN_[i].resize(size, 0);
    		} catch (std::bad_alloc &ba) {
        		throw customException ("Out of memory for composition histogram");
    		}
	}
	for (unsigned int i = 0; i < averageUNi_.size(); ++i) {
		try {
        		averageUNi_[i].resize(size, 0);
    		} catch (std::bad_alloc &ba) {
        		throw customException ("Out of memory for <UN_i> histogram");
    		}
	}
	for (unsigned int i = 0; i < averageNiNj_.size(); ++i) {
		try {
        		averageNiNj_[i].resize(size, 0);
    		} catch (std::bad_alloc &ba) {
        		throw customException ("Out of memory for <N_iN_j> histogram");
    		}
	}
	try {
        	numAverageN_.resize(size, 0);
    	} catch (std::bad_alloc &ba) {
        	throw customException ("Out of memory for composition histogram");
    	}*/

    energyHistogram_.resize(0);
    energyHistogram_lb_.resize(size, -5.0);
    energyHistogram_ub_.resize(size, 5.0);

    for (unsigned int i = 0; i < size; ++i) {
    	try {
    		dynamic_one_dim_histogram dummyHist (energyHistogram_lb_[i], energyHistogram_ub_[i], energyHistDelta_);
    		energyHistogram_.resize(i+1, dummyHist);
    	} catch (std::bad_alloc &ba) {
    		throw customException ("Out of memory for energy histogram for each Ntot");
    	}
    }

    pkHistogram_.resize(0);
    dynamic_one_dim_histogram dummyPkHist (0.0, totNBounds_[1], 1.0);
    try {
    	std::vector < dynamic_one_dim_histogram > tmp (totNBounds_[1]-totNBounds_[0]+1, dummyPkHist);
		pkHistogram_.resize(nSpecies_, tmp);
	 } catch (std::bad_alloc &ba) {
		throw customException ("Out of memory for particle histogram for each Ntot");
   	}

   	// initialize moments
	std::vector < double > lbn (6,0), ubn(6,0);
	std::vector < long long unsigned int > nbn (6,0);
	ubn[0] = nSpecies_-1;
	ubn[1] = max_order_;
	ubn[2] = nSpecies_-1;
	ubn[3] = max_order_;
	ubn[4] = max_order_;
	ubn[5] = totNBounds_[1]-totNBounds_[0];

	nbn[0] = nSpecies_;
	nbn[1] = max_order_+1;
	nbn[2] = nSpecies_;
	nbn[3] = max_order_+1;
	nbn[4] = max_order_+1;
	nbn[5] = size;

	histogram hnn (lbn, ubn, nbn);
	extensive_moments_ = hnn;
}

/*!
 * Insert an atom into the system. Does all the bookkeepping behind the scenes.
 *
 * \param [in] typeIndex What type the atom is (>= 0)
 * \param [in] newAtom Pointer to new atom.  A copy is stored in the system so the original may be destroyed.
 * \param [in] override Override command that prevents the expanded ensemble state from being changed.  Used during swap moves where "insertions" are temporary.
 */
void simSystem::insertAtom (const int typeIndex, atom *newAtom, bool override) {
	if (typeIndex < nSpecies_ && typeIndex >= 0) {
		if (numSpecies[typeIndex] < maxSpecies_[typeIndex]) {
			if (Mtot_ > 1 && !override) {
        			// expanded ensemble behavior, "normal" insertion and deletion
        			if (Mcurrent_ > 0) { // further inserting an atom that already partially exists in the system
        				// ensure the system pointer is correct if currently a partially inserted atom
        				if (fractionalAtom_ != newAtom || typeIndex != fractionalAtomType_) {
        					throw customException ("Fractional atom pointer does not point to atom believed to be inserted");
        				}

	        			// increment expanded state
        				fractionalAtom_->mState++;
        				Mcurrent_++;

 	       				// check if now fully inserted
        				if (fractionalAtom_->mState == Mtot_) {
        					fractionalAtom_->mState = 0;
	        				Mcurrent_ = 0;
        					totN_++;
        					numSpecies[typeIndex]++;
        				}
        			} else {
	        			// inserting a new atom for the first time
        				atoms[typeIndex][numSpecies[typeIndex]] = (*newAtom);

        				// assign fractional atom
        				fractionalAtom_ = &atoms[typeIndex][numSpecies[typeIndex]];
	        			fractionalAtomType_ = typeIndex;

        				// increment expanded state
        				fractionalAtom_->mState = 1;
        				Mcurrent_ = 1;

					// add particle into appropriate cell lists
					for (unsigned int i = 0; i < nSpecies_; ++i) {
						if (useCellList_[typeIndex][i]) {
							cellList* cl = cellListsByPairType_[typeIndex][i];
							cl->insertParticle(&atoms[typeIndex][numSpecies[typeIndex]]); // numSpecies[typeIndex] is the number of fully inserted ones, this partially inserted one comes after that
						}
					}
        			}
        		} else if (Mtot_ > 1 && override) {
	        		// expanded ensemble behavior, but now amidst a "swap move" rather than an actual insertion or deletion.
        			// for this, insertions involve just putting the atom "back" into the system / cellLists after being artificially completely removed

        			// ensure we insert at the proper "end"
        			int end = numSpecies[typeIndex];
        			if (Mcurrent_ > 0 && typeIndex == fractionalAtomType_ && newAtom->mState == 0) {
        				end++; // insert after the partially inserted one since newAtom is NOT the partial one
	        		}
        			atoms[typeIndex][end] = (*newAtom);

        			// if we just added a partially inserted/deleted particle back to the system, need to update the pointer
	        		if (atoms[typeIndex][end].mState != 0) {
        				fractionalAtom_ = &atoms[typeIndex][end];
        				fractionalAtomType_ = typeIndex;

        				// set the system's mState back to that of the atom just inserted, iff it was the partial one
        				Mcurrent_ = atoms[typeIndex][end].mState;
        			} else {
        				totN_++; // we just added a "full" atom
        				numSpecies[typeIndex]++; // we just added a "full" atom
        			}

	        		// put newAtom into the cell lists whatever its state
               			for (unsigned int i = 0; i < nSpecies_; ++i) {
        	        		if (useCellList_[typeIndex][i]) {
                 				cellList* cl = cellListsByPairType_[typeIndex][i];
                 				cl->insertParticle(&atoms[typeIndex][end]);
                 			}
                 		}
        		} else {
	        		// direct insertion (no expanded ensemble)
        	       		atoms[typeIndex][numSpecies[typeIndex]] = (*newAtom);
               			numSpecies[typeIndex]++;
                		totN_++;

	               		// add particle into appropriate cell lists
        	       		for (unsigned int i = 0; i < nSpecies_; ++i) {
               				if (useCellList_[typeIndex][i]) {
                				cellList* cl = cellListsByPairType_[typeIndex][i];
                				cl->insertParticle(&atoms[typeIndex][numSpecies[typeIndex] - 1]);
                			}
                		}
        		}
        	} else {
            		std::string index = static_cast<std::ostringstream*>( &(std::ostringstream() << typeIndex) )->str();
            		throw customException ("Reached upper bound, cannot insert an atom of type index "+index);
        	}
	} else {
        	throw customException ("That species index does not exist, cannot insert an atom");
    	}
}

/*!
 * Delete an atom from the system. Does all the bookkeepping behind the scenes.
 *
 * \param [in] typeIndex What type the atom is (>= 0)
 * \param [in] atomIndex Which atom *index* of type typeIndex to destroy (>= 0)
 * \param [in] Optional override command which allows the system to delete a particle even it goes below the minimum allowed. E.g. during a swap move.
 */
void simSystem::deleteAtom (const int typeIndex, const int atomIndex, bool override) {
	if (typeIndex < nSpecies_ && typeIndex >= 0) {
       		if ((numSpecies[typeIndex] > minSpecies_[typeIndex]) || ((numSpecies[typeIndex] == minSpecies_[typeIndex]) && (Mcurrent_ > 0)) || override) {
        		if (override) {
        			// doing a swap move
        			if (Mtot_ > 1) {
        				// expanded ensemble and not necessarily deleting the partial atom

					int end = numSpecies[typeIndex] - 1;
					if (fractionalAtomType_ == typeIndex && Mcurrent_ > 0) {
						// we are deleting a particle which has to watch out for the partial atom
						end++;
					}

					if (atoms[typeIndex][atomIndex].mState == 0) {
						// if we are removing a "full" particle, have to decrement Ntot, else not
						numSpecies[typeIndex]--;
						totN_--;
					} else {
						// but if removing the partial particle, M is affected
						Mcurrent_ = 0; // regardless of how M was originally, the partial particle is now "entirely" gone
					}

					bool replace = false;
					if (&atoms[typeIndex][end] == fractionalAtom_) {
						// then the fractional atom is about to be used to replace a "full" one
						replace = true;
					}

					// have to entirely remove the particle
					for (unsigned int i = 0; i < nSpecies_; ++i) {
						if (useCellList_[typeIndex][i]) {
							cellList* cl = cellListsByPairType_[typeIndex][i];
							cl->swapAndDeleteParticle(&atoms[typeIndex][atomIndex], &atoms[typeIndex][end]);
						}
					}

					atoms[typeIndex][atomIndex] = atoms[typeIndex][end];    // "replacement" operation

					if (replace) {
						fractionalAtom_ = &atoms[typeIndex][atomIndex];	// update the pointer if necessary
					}
        			} else {
        				// no expanded ensemble, just delete particle from appropriate cell list
                    			for (unsigned int i = 0; i < nSpecies_; ++i) {
                    				if (useCellList_[typeIndex][i]) {
                    					cellList* cl = cellListsByPairType_[typeIndex][i];
                    					cl->swapAndDeleteParticle(&atoms[typeIndex][atomIndex], &atoms[typeIndex][numSpecies[typeIndex] - 1]);
                    				}
                   			}

                   			atoms[typeIndex][atomIndex] = atoms[typeIndex][numSpecies[typeIndex] - 1];    // "replacement" operation
                   			numSpecies[typeIndex]--;
                   			totN_--;
        			}
        		} else {
        			// not doing a swap move, just a "regular" deletion
        			if (Mtot_ > 1) {
        				// expanded ensemble
               				if (Mcurrent_ == 1) {
               					// when we delete this atom, it is entirely gone

               					// first ensure the system pointer is correct if currently a partially inserted atom
            					if (fractionalAtom_ != &atoms[typeIndex][atomIndex] || typeIndex != fractionalAtomType_) {
            						throw customException ("Fractional atom pointer does not point to atom belived to be inserted");
            					}

            					// decrement expanded state
            					fractionalAtom_->mState = 0;
            					Mcurrent_ = 0;

            					// since deleting partial particle, do not update Ntot, etc.
            					// however, do have to remove from cellLists
            					int end = numSpecies[typeIndex]; // includes space for the partially inserted one currently in cellList
                        			for (unsigned int i = 0; i < nSpecies_; ++i) {
                        				if (useCellList_[typeIndex][i]) {
                        					cellList* cl = cellListsByPairType_[typeIndex][i];
                        					cl->swapAndDeleteParticle(&atoms[typeIndex][atomIndex], &atoms[typeIndex][end]);
                        				}
                        			}

						atoms[typeIndex][atomIndex] = atoms[typeIndex][end];    // "replacement" operation

               				} else if (Mcurrent_ == 0) {
               					// have to decrement Ntot, but keep in cell lists
               					numSpecies[typeIndex]--;
               					totN_--;

            					// this is a new fractional atom
            					fractionalAtom_ = &atoms[typeIndex][atomIndex];
            					fractionalAtomType_ = typeIndex;

            					// decrement expanded state
            					fractionalAtom_->mState = Mtot_-1;
            					Mcurrent_ = Mtot_-1;
               				} else {
               					// further deleting an atom that already partially exists in the system, but remains in cell lists

               					// first ensure the system pointer is correct if currently a partially inserted atom
            					if (fractionalAtom_ != &atoms[typeIndex][atomIndex] || typeIndex != fractionalAtomType_) {
            						throw customException ("Fractional atom pointer does not point to atom belived to be inserted");
            					}

               					// decrement expanded state
               					fractionalAtom_->mState -= 1;
               					Mcurrent_ -= 1;
                			}
	        		} else {
        				// no expanded ensemble, just delete particle from appropriate cell list
                    			for (unsigned int i = 0; i < nSpecies_; ++i) {
                    				if (useCellList_[typeIndex][i]) {
                    					cellList* cl = cellListsByPairType_[typeIndex][i];
                    					cl->swapAndDeleteParticle(&atoms[typeIndex][atomIndex], &atoms[typeIndex][numSpecies[typeIndex] - 1]);
                    				}
                    			}

                    			atoms[typeIndex][atomIndex] = atoms[typeIndex][numSpecies[typeIndex] - 1];    // "replacement" operation
                   			numSpecies[typeIndex]--;
                   			totN_--;
        			}
        		}
       		} else {
            		std::string index = static_cast<std::ostringstream*>( &(std::ostringstream() << typeIndex) )->str();
           		throw customException ("System going below minimum allowable number of atoms, cannot delete an atom of type index "+index);
       		}
    	} else {
        	throw customException ("That species index does not exist, cannot delete an atom");
    	}
}

/*!
 * Translate an atom in the system. Does all the bookkeeping behind the scenes.
 * Do nothing if there is no cell list defined for the type
 *
 * \param [in] typeIndex What type the atom is (>= 0)
 * \param [in] atomIndex Which atom *index* of type typeIndex to translate (>= 0)
 * \param [in] oldPos Old position of the atom.  The current/new position should already be stored in the atom at sys.atoms[typeIndex][atomIndex]
 */
void simSystem::translateAtom (const int typeIndex, const int atomIndex, std::vector<double> oldPos) {
	if (typeIndex < nSpecies_ && typeIndex >= 0) {
        	if (atomIndex >= 0) {
        		// delete particle from appropriate cell list, move to new one
            		for (unsigned int i=0; i<nSpecies_; i++) {
            			if (useCellList_[typeIndex][i]) {
	            			cellList* cl = cellListsByPairType_[typeIndex][i];
	        	    		cl->translateParticle(&atoms[typeIndex][atomIndex], oldPos);
        	    		}
            		}
        	} else {
            		std::string index = static_cast<std::ostringstream*>( &(std::ostringstream() << typeIndex) )->str();
            		throw customException ("Number of those atoms in system is out of bounds, cannot translate an atom of type index "+index);
        	}
    	} else {
       		throw customException ("That species index does not exist, cannot translate the atom");
    	}
}

/*!
 * Toggle KE adjustment to energy setting
 */
void simSystem::toggleKE() {
	if (toggleKE_ == false) {
		toggleKE_ = true;
	} else {
		toggleKE_ = false;
	}
}

/*
 * Destructor frees biases if used and not turned off already.
 */
simSystem::~simSystem () {
	if (useTMMC) {
		delete tmmcBias;
	}
	if (useWALA) {
		delete wlBias;
	}
	for (unsigned int i = 0; i < ppot.size(); ++i) {
		for (unsigned int j = 0; j < ppot[i].size(); ++j) {
			delete ppot[i][j];
		}
		ppot[i].clear();
	}
	ppot.clear();
}

/*!
 * Initialize the system. Sets the use of both WL and TMMC biasing to false.
 *
 * \param [in] nSpecies Number of unqiue species types to allow in the system
 * \param [in] beta Inverse temperature (1/kT)
 * \param [in] box Box dimensions [x, y, z]
 * \param [in] mu Chemical potential of each species
 * \param [in] maxSpecies Maximum number of each species to allow in the system
 * \param [in] Mtot Total number of expanded ensemble states
 * \param [in] energyHistDelta Bin width of energy histogram at each Ntot (optional, default = 10.0)
 * \param [in] max_order Maximum order to record correlations to (default = 2)
 */
simSystem::simSystem (const unsigned int nSpecies, const double beta, const std::vector < double > box, const std::vector < double > mu, const std::vector < int > maxSpecies, const std::vector < int > minSpecies, const int Mtot, const double energyHistDelta, const int max_order) {
	if ((box.size() != 3) || (nSpecies != mu.size()) || (maxSpecies.size() != nSpecies)) {
		throw customException ("Invalid system initialization parameters");
		exit(SYS_FAILURE);
	} else {
		nSpecies_ = nSpecies;
		maxSpecies_ = maxSpecies;
		minSpecies_ = minSpecies;
		box_ = box;
		mu_ = mu;
		beta_ = beta;
	}

	lnF_start = 1.0; // default for lnF_start
	lnF_end = 2.0e-18; // default for lnF_end
	toggleKE_ = false; //default, do NOT adjust energy by kinetic contribution of 3/2kT per atom (just record PE)
	totalTMMCSweeps = 0;
	wlSweepSize = 0;
	wala_g = 0.5;
	wala_s = 0.8;
	nCrossoverVisits = 5;

	if (max_order < 1){
		throw customException ("max_order must be >= 1");
	}
	max_order_ = max_order;

	if (energyHistDelta <= 0) {
		throw customException ("energyHistDelta must be > 0");
	}
	energyHistDelta_ = energyHistDelta;

	for (unsigned int i = 0; i < 3; ++i) {
		if (box_[i] <= 0) {
			throw customException ("Box dimensions must be > 0");
		}
	}

	if (Mtot < 1) {
		throw customException ("Total fractional states for expanded ensemble must be >= 1");
	}
	Mtot_ = Mtot;
	Mcurrent_ = 0; // always start from fully inserted state

	try {
		ppot.resize(nSpecies);
	} catch (std::exception &e) {
		throw customException (e.what());
	}
	for (unsigned int i = 0; i < nSpecies; ++i) {
		try {
			ppot[i].resize(nSpecies);
		} catch (std::exception &e) {
			throw customException (e.what());
		}
	}

	try {
		mass_.resize(nSpecies, 1.0);
	} catch (std::exception &e) {
		throw customException (e.what());
	}

	try {
		ppotSet_.resize(nSpecies);
	} catch (std::exception &e) {
		throw customException (e.what());
	}
	for (unsigned int i = 0; i < nSpecies; ++i) {
		try {
			ppotSet_[i].resize(nSpecies, false);
		} catch (std::exception &e) {
			throw customException (e.what());
		}
	}

	// Wall potentials for each species, if there are any?
	try {
		speciesBarriers.resize(nSpecies);
	} catch (std::exception &e) {
		throw customException (e.what());
	}

	// Prepare vectors and matrices for cell lists.
	// It is crucial to reserve the correct number of cellLists in advance
	// since cellListsByPairType uses the addresses of cellLists. Otherwise,
	// if dynamic memory reallocation takes place, the pointers do not
	// correspond to initial values anymore, causing the simulation to crash.
	cellLists_.reserve(nSpecies_*nSpecies_);

	try {
		useCellList_.resize(nSpecies);
		cellListsByPairType_.resize(nSpecies);
	} catch (std::exception &e) {
		throw customException (e.what());
	}
	for (unsigned int i = 0; i < nSpecies; ++i) {
		try {
			useCellList_[i].resize(nSpecies);
			cellListsByPairType_[i].assign(nSpecies, NULL);
		} catch (std::exception &e) {
			throw customException (e.what());
		}
	}

	totN_ = 0;
    	try {
		numSpecies.resize(nSpecies, 0);
	} catch (std::exception &e) {
		throw customException (e.what());
	}

    	try {
        	atoms.resize(nSpecies);
   	} catch (std::exception &e) {
        	throw customException (e.what());
    	}
    	for (unsigned int i = 0; i < nSpecies; ++i) {
        	if (minSpecies_[i] < 0) {
            		throw customException ("Min species < 0");
        	}
        	if (maxSpecies_[i] < minSpecies_[i]) {
            		throw customException ("Max species < Min species");
        	}
		try {
			atoms[i].resize(maxSpecies_[i]);
		} catch (std::exception &e) {
			throw customException (e.what());
		}
	}

    	energy_ = 0.0;

    	useTMMC = false;
    	useWALA = false;

    	totNBounds_.resize(2, 0);
    	for (unsigned int i = 0; i < nSpecies_; ++i) {
    		totNBounds_[0] += minSpecies_[i];
    		totNBounds_[1] += maxSpecies_[i];
    	}

    	// allocate space for average U storage matrix - Shen and Errington method implies this size is always the same for
    	// both single and multicomponent mixtures
    	long long int size = totNBounds_[1] - totNBounds_[0] + 1;
    	energyHistogram_lb_.resize(size, -5.0);
    	energyHistogram_ub_.resize(size, 5.0);
    	for (unsigned int i = 0; i < size; ++i) {
    		energyHistogram_lb_[i] = -5.0;
    		energyHistogram_ub_[i] = 5.0;
    		try {
    			dynamic_one_dim_histogram dummyHist (energyHistogram_lb_[i], energyHistogram_ub_[i], energyHistDelta_);
    			energyHistogram_.resize(i+1, dummyHist);
    		} catch (std::bad_alloc &ba) {
    			throw customException ("Out of memory for energy histogram for each Ntot");
    		}
    	}
    	pkHistogram_.resize(0);
    	dynamic_one_dim_histogram dummyPkHist (0.0, totNBounds_[1], 1.0);
		try {
			std::vector < dynamic_one_dim_histogram > tmp (totNBounds_[1]-totNBounds_[0]+1, dummyPkHist);
			pkHistogram_.resize(nSpecies_, tmp);
		} catch (std::bad_alloc &ba) {
    		throw customException ("Out of memory for particle histogram for each Ntot");
    	}

    	/*try {
        	numAverageU_.resize(size, 0);
    	} catch (std::bad_alloc &ba) {
     		throw customException ("Out of memory for energy record");
    	}
	try {
        	numAverageFlucts_.resize(size, 0);
    	} catch (std::bad_alloc &ba) {
     		throw customException ("Out of memory for fluctuation record");
    	}

    	try {
       		AverageU_.resize(size, 0);
    	} catch (std::bad_alloc &ba) {
        	throw customException ("Out of memory for energy record");
    	}
	try {
       		averageU2_.resize(size, 0);
    	} catch (std::bad_alloc &ba) {
        	throw customException ("Out of memory for energy^2 record");
    	}

	try {
                numAverageN_.resize(size, 0);
        } catch (std::bad_alloc &ba) {
                throw customException ("Out of memory for particle histogram");
        }
	try {
		averageN_.resize(nSpecies_);
	} catch (std::bad_alloc &ba) {
		throw customException ("Out of memory for particle histogram");
	}
	try {
		averageUNi_.resize(nSpecies_);
	} catch (std::bad_alloc &ba) {
		throw customException ("Out of memory for particle histogram");
	}
	try {
		averageNiNj_.resize(nSpecies_*nSpecies_);
	} catch (std::bad_alloc &ba) {
		throw customException ("Out of memory for <N_iN_j> histogram");
	}

	for (unsigned int i = 0; i < nSpecies_; ++i) {
		try {
			averageN_[i].resize(size, 0);
		} catch (std::bad_alloc &ba) {
			throw customException ("Out of memory for particle histogram");
		}
		try {
			averageUNi_[i].resize(size, 0);
		} catch (std::bad_alloc &ba) {
			throw customException ("Out of memory for <UN_i> histogram");
		}
	}
	for (unsigned int i = 0; i < nSpecies_*nSpecies_; ++i) {
		try {
			averageNiNj_[i].resize(size, 0);
		} catch (std::bad_alloc &ba) {
			throw customException ("Out of memory for <N_iN_j> histogram");
		}
	}*/

	// initialize moments
	std::vector < double > lbn (6,0), ubn(6,0);
	std::vector < long long unsigned int > nbn (6,0);
	ubn[0] = nSpecies_-1;
	ubn[1] = max_order_;
	ubn[2] = nSpecies_-1;
	ubn[3] = max_order_;
	ubn[4] = max_order_;
	ubn[5] = totNBounds_[1]-totNBounds_[0];

	nbn[0] = nSpecies_;
	nbn[1] = max_order_+1;
	nbn[2] = nSpecies_;
	nbn[3] = max_order_+1;
	nbn[4] = max_order_+1;
	nbn[5] = size;

	histogram hnn (lbn, ubn, nbn);
	extensive_moments_ = hnn;
}

/*!
 * Save the instantaneous fluctuation properties, e.g X_i*Y_j, for various properties, X, Y, such as energy, etc.
 * Only records values when N_tot in range of [min, max].
 */
/*void simSystem::recordFluctuation () {
	// only record if in range (removes equilibration stage to get in this range, if there was any)
        if (totN_ >= totNBounds_[0] && totN_ <= totNBounds_[1]) {
                const long int address = totN_-totNBounds_[0];
		numAverageFlucts_[address] += 1.0;
		averageU2_[address] += energy_*energy_;
		for (unsigned int i = 0; i < nSpecies_; ++i) {
                	averageUNi_[i][address] += numSpecies[i]*energy_;
			for (unsigned int j = 0; j < nSpecies_; ++j) {
				averageNiNj_[j+i*nSpecies_][address] += numSpecies[i]*numSpecies[j];
			}
        	}
	}
}*/

/*!
 * Print the average properties for fluctuation calculations to file for every N_tot within range recorded. Will overwrite the file if another with that name exists. Prints in netCDF format if enabled.
 *
 * \param [in] fileName Name of the file to print to
 */
/*void simSystem::printFluctuation (const std::string fileName) {
	std::vector < std::vector < long double > > ave_ninj, ave_uni;
	std::vector < long double > ave_uu;
	try {
		ave_ninj.resize(averageNiNj_.size());
	} catch (std::bad_alloc &ba) {
		throw customException ("Out of memory for fluctuation printing");
	}
	for (unsigned int i = 0; i < ave_ninj.size(); ++i) {
		try {
                	ave_ninj[i].resize(averageNiNj_[i].size(), 0);
        	} catch (std::bad_alloc &ba) {
                	throw customException ("Out of memory for fluctuation printing");
        	}
	}
	for (unsigned int i = 0; i < averageNiNj_.size(); ++i) {
		for (unsigned int j = 0; j < averageNiNj_[i].size(); ++j) {
			ave_ninj[i][j] = averageNiNj_[i][j]/numAverageFlucts_[j];
		}
	}

	try {
		ave_uni.resize(averageUNi_.size());
	} catch (std::bad_alloc &ba) {
		throw customException ("Out of memory for fluctuation printing");
	}
	for (unsigned int i = 0; i < ave_uni.size(); ++i) {
		try {
                	ave_uni[i].resize(averageUNi_[i].size(), 0);
        	} catch (std::bad_alloc &ba) {
                	throw customException ("Out of memory for fluctuation printing");
        	}
	}
	for (unsigned int i = 0; i < averageUNi_.size(); ++i) {
		for (unsigned int j = 0; j < averageUNi_[i].size(); ++j) {
			ave_uni[i][j] = averageUNi_[i][j]/numAverageFlucts_[j];
		}
	}

	try {
		ave_uu.resize(averageU2_.size());
	} catch (std::bad_alloc &ba) {
		throw customException ("Out of memory for fluctuation printing");
	}
	for (unsigned int i = 0; i < averageU2_.size(); ++i) {
		ave_uu[i] = averageU2_[i]/numAverageFlucts_[i];
	}

#ifdef NETCDF_CAPABLE
	// If netCDF libs are enabled, write to this format
    	const std::string name = fileName + ".nc";
        NcFile outFile(name.c_str(), NcFile::replace);
        NcDim probDim = outFile.addDim("vectorized_position", averageU2_[0].size());
        std::vector < NcVar > probVar (1+nSpecies_+nSpecies_*nSpecies_); // U^2, UN_i ..., NiNj's
	std::string vName = "average_UU";
	int idx = 0;
	probVar[idx] = outFile.addVar(vName.c_str(), ncdouble, probDim);
	idx++;
	for (unsigned int i = 0; i < nSpecies_; ++i) {
		vName = "average_UN_"+sstr(i+1);
		probVar[idx] = outFile.addVar(vName.c_str(), ncdouble, probDim);
		idx++;
	}
	for (unsigned int i = 0; i < nSpecies_; ++i) {
		for (unsigned int j = 0; j < nSpecies_; ++j) {
			vName = "average_N_"+sstr(i)+"_N_"+sstr(j);
			probVar[idx+j+i*nSpecies_] = outFile.addVar(vName.c_str(), ncdouble, probDim);
		}
	}

	std::string attName = "species_total_upper_bound";
        probVar[0].putAtt(attName.c_str(), sstr(totNBounds_[1]).c_str());
        attName = "species_total_lower_bound";
        probVar[0].putAtt(attName.c_str(), sstr(totNBounds_[0]).c_str());
        const std::string dummyName = "number_species";
        probVar[0].putAtt(dummyName.c_str(), sstr(nSpecies_).c_str());
	attName = "volume";
        double V = box_[0]*box_[1]*box_[2];
        probVar[0].putAtt(attName.c_str(), sstr(V).c_str());

	idx = 0;
	probVar[idx].putVar(&ave_uu);
	idx++;
        for (unsigned int i = 0; i < nSpecies_; ++i) {
		probVar[idx].putVar(&ave_uni[i]);
		idx++;
	}
	for (unsigned int i = 0; i < nSpecies_; ++i) {
		for (unsigned int j = 0; j < nSpecies_; ++j) {
			probVar[idx+j+i*nSpecies_].putVar(&ave_ninj[j+i*nSpecies_]);
		}
	}

#else
	// Without netCDF capabilities, just print to ASCII file
        std::ofstream of;
        std::string name = fileName+".dat";
        of.open(name.c_str(), std::ofstream::out);
        of << "# <X_i*Y_i> as a function of N_tot." << std::endl;
        of << "# Number of species: " << nSpecies_ << std::endl;
        of << "# species_total_upper_bound: " << totNBounds_[1] << std::endl;
        of << "# species_total_lower_bound: " << totNBounds_[0] << std::endl;
        double V = box_[0]*box_[1]*box_[2];
        of << "# volume: " << std::setprecision(15) << V << std::endl;
	of << "# <U*U>\t{<U*N_1>\t<U*N_2>\t...}\t{<N_1*N_1>\t<N_1*N_2>\t...<N_1*N_n>\t<N_2*N_1>...\t<N_2*N_n>\t...\t<N_n*N_1>\t...\t<N_n*N_n>" << std::endl;
	for (unsigned long long int i = 0; i < ave_uu.size(); ++i) {
		of << ave_uu[i] << "\t";
		for (unsigned long long int j = 0; j < nSpecies_; ++j) {
			of << ave_uni[j][i] << "\t";
		}
		for (unsigned long long int j = 0; j < nSpecies_; ++j) {
			for (unsigned long long int k = 0; k < nSpecies_; ++k) {
				of << ave_ninj[k+j*nSpecies_][i] << "\t";
			}
		}
		of << std::endl;
	}
        of.close();
#endif
}*/

/*!
 * Save the instantaneous number of each component as a function of the total number of particles in the system.
 * Only records values when N_tot in range of [min, max].
 */
/*void simSystem::recordComposition () {
	// only record if in range (removes equilibration stage to get in this range, if there was any)
        if (totN_ >= totNBounds_[0] && totN_ <= totNBounds_[1]) {
                const long int address = totN_-totNBounds_[0];
		numAverageN_[address] += 1.0;
		for (unsigned int i = 0; i < nSpecies_; ++i) {
                	averageN_[i][address] += numSpecies[i];
        	}
	}
}*/

/*!
 * Print the average composition to file, <N_i> for every N_tot within range recorded. Will overwrite the file if another with that name exists. Prints in netCDF format if enabled.
 *
 * \param [in] fileName Name of the file to print to
 */
/*void simSystem::printComposition (const std::string fileName) {
	std::vector < std::vector < double > > aveX;
	try {
		aveX.resize(averageN_.size());
	} catch (std::bad_alloc &ba) {
		throw customException ("Out of memory for composition printing");
	}
	for (unsigned int i = 0; i < aveX.size(); ++i) {
		try {
                	aveX[i].resize(averageN_[i].size(), 0);
        	} catch (std::bad_alloc &ba) {
                	throw customException ("Out of memory for composition printing");
        	}
	}
	for (unsigned int i = 0; i < averageN_.size(); ++i) {
		for (unsigned int j = 0; j < averageN_[i].size(); ++j) {
			aveX[i][j] = averageN_[i][j]/numAverageN_[j];
		}
	}

#ifdef NETCDF_CAPABLE
	// If netCDF libs are enabled, write to this format
    	const std::string name = fileName + ".nc";
        NcFile outFile(name.c_str(), NcFile::replace);
        NcDim probDim = outFile.addDim("vectorized_position", averageN_[0].size());
        std::vector < NcVar > probVar (nSpecies_);
	for (unsigned int i = 0; i < nSpecies_; ++i) {
		std::string vName = "average_N_"+sstr(i+1);
		probVar[i] = outFile.addVar(vName.c_str(), ncdouble, probDim);
	}
	std::string attName = "species_total_upper_bound";
        probVar[0].putAtt(attName.c_str(), sstr(totNBounds_[1]).c_str());
        attName = "species_total_lower_bound";
        probVar[0].putAtt(attName.c_str(), sstr(totNBounds_[0]).c_str());
        const std::string dummyName = "number_species";
        probVar[0].putAtt(dummyName.c_str(), sstr(nSpecies_).c_str());
	attName = "volume";
        double V = box_[0]*box_[1]*box_[2];
        probVar[0].putAtt(attName.c_str(), sstr(V).c_str());
        for (unsigned int i = 0; i < nSpecies_; ++i) {
		probVar[i].putVar(&aveX[i]);
	}
#else
	// Without netCDF capabilities, just print to ASCII file
        std::ofstream of;
        std::string name = fileName+".dat";
        of.open(name.c_str(), std::ofstream::out);
        of << "# <N_i> as a function of N_tot." << std::endl;
        of << "# Number of species: " << nSpecies_ << std::endl;
        of << "# species_total_upper_bound: " << totNBounds_[1] << std::endl;
        of << "# species_total_lower_bound: " << totNBounds_[0] << std::endl;
        double V = box_[0]*box_[1]*box_[2];
        of << "# volume: " << std::setprecision(15) << V << std::endl;
	of << "# <N_1>\t<N_2>\t...\t<N_n>" << std::endl;
        for (unsigned long long int i = 0; i < aveX[0].size(); ++i) {
		for (unsigned long long int j = 0; j < nSpecies_; ++j) {
			of << aveX[j][i] << "\t";
		}
                of << std::endl;
        }
        of.close();
#endif
}*/

/*!
 * Record the extensive moment at a given Ntot.
 */
void simSystem::recordExtMoments () {
	// only record if in range (removes equilibration stage to get in this range, if there was any)
	if (totN_ >= totNBounds_[0] && totN_ <= totNBounds_[1]) {
		double val = 0.0;
		std::vector < double > coords (6,0);
		coords[5] = totN_-totNBounds_[0];
		for (unsigned int i = 0; i < nSpecies_; ++i) {
			coords[0] = i;
			for (unsigned int j = 0; j <= max_order_; ++j) {
				coords[1] = j;
				for (unsigned int k = 0; k < nSpecies_; ++k) {
					coords[2] = k;
					for (unsigned int m = 0; m <= max_order_; ++m) {
						coords[3] = m;
						for (unsigned int p = 0; p <= max_order_; ++p) {
							coords[4] = p;
							val = pow(numSpecies[i], j)*pow(numSpecies[k], m)*pow(energy_, p);
							extensive_moments_.increment (coords, val);
						}
					}
				}
			}
		}
	}
}

/*!
 * Print the (normalized) extensive energy histogram for each Ntot. netCDF4 not enabled
 *
 * \param [in] fileName Name of the file to print to
 */
void simSystem::printExtMoments (const std::string fileName) {
#ifdef NETCDF_CAPABLE
    throw customException ("Cannot record the extensive moments for each Ntot in netCDF4 format.");
#else
	// Without netCDF capabilities, just print to ASCII file
	std::ofstream of;
	std::string name = fileName+".dat";
	of.open(name.c_str(), std::ofstream::out);
	of << "# <N_i^j*N_k^m*U^p> as a function of N_tot." << std::endl;
	of << "# number_of_species: " << nSpecies_ << std::endl;
	of << "# max_order: " << max_order_ << std::endl;
	of << "# species_total_upper_bound: " << totNBounds_[1] << std::endl;
	of << "# species_total_lower_bound: " << totNBounds_[0] << std::endl;
	double V = box_[0]*box_[1]*box_[2];
	of << "# volume: " << std::setprecision(15) << V << std::endl;
	of << "#\tN_tot\t";
	for (unsigned int i = 0; i < nSpecies_; ++i) {
		for (unsigned int j = 0; j <= max_order_; ++j) {
			for (unsigned int k = 0; k < nSpecies_; ++k) {
				for (unsigned int m = 0; m <= max_order_; ++m) {
					for (unsigned int p = 0; p <= max_order_; ++p) {
						of << "N_"+sstr(i+1)+"^"+sstr(j)+"*N_"+sstr(k+1)+"^"+sstr(m)+"*U^"+sstr(p)+"\t";
					}
				}
			}
		}
	}
	of << std::endl;
	std::vector <double> h = extensive_moments_.getRawHistogram ();
	std::vector <double> ctr = extensive_moments_.getCounter ();
	std::vector <double> coords (6,0);
	long unsigned int idx = 0;
	for (unsigned int n = 0; n < totNBounds_[1]-totNBounds_[0]+1; ++n) {
		of << n+totNBounds_[0] << "\t";
		coords[5] = n;
		for (unsigned int i = 0; i < nSpecies_; ++i) {
			coords[0] = i;
			for (unsigned int j = 0; j <= max_order_; ++j) {
				coords[1] = j;
				for (unsigned int k = 0; k < nSpecies_; ++k) {
					coords[2] = k;
					for (unsigned int m = 0; m <= max_order_; ++m) {
						coords[3] = m;
						for (unsigned int p = 0; p <= max_order_; ++p) {
							coords[4] = p;
							idx = extensive_moments_.getAddress(coords);
							of << h[idx]/ctr[idx] << "\t";
						}
					}
				}
			}
		}
		of << std::endl;
	}
	of.close();
#endif
}

/*!
 * Record the energy histogram for the system at a given Ntot.
 * Only records values when N_tot in range of [min, max].
 */
void simSystem::recordEnergyHistogram () {
	// only record if in range (removes equilibration stage to get in this range, if there was any)
	if (totN_ >= totNBounds_[0] && totN_ <= totNBounds_[1]) {
		const int address = totN_-totNBounds_[0];
		energyHistogram_[address].record(energy_);
	}
}

/*!
 * Monitor the energy histogram bounds at each Ntot
 */
void simSystem::checkEnergyHistogramBounds () {
	if (totN_ >= totNBounds_[0] && totN_ <= totNBounds_[1]) {
		const int address = totN_-totNBounds_[0];
		energyHistogram_lb_[address] = std::min(energyHistogram_lb_[address], energy_);
		energyHistogram_ub_[address] = std::max(energyHistogram_ub_[address], energy_);
	}
}

/*!
 * Check the histogram entries and trim off zero-valued entries and bounds
 */
void simSystem::refineEnergyHistogramBounds () {
	for (std::vector < dynamic_one_dim_histogram >::iterator it = energyHistogram_.begin(); it != energyHistogram_.end(); ++it) {
		try {
			it->trim_edges();
		} catch (customException &ce) {
			std::string a = "Unable to trim edges in energyHistogram at each Ntot: ", b = ce.what();
			throw customException (a+b);
		}
	}
}

/*!
 * Re-initialize the energy histogram with internal estimates of bounds.
 */
void simSystem::reInitializeEnergyHistogram () {
	double lb = 0.0, ub = 0.0;
	if (energyHistogram_lb_.size() != energyHistogram_ub_.size()) {
			throw customException ("Bad energy histogram bound sizes");
		}
		if (energyHistogram_lb_.size() != totNMax() - totNMin() + 1) {
			throw customException ("Bad energy histogram bound sizes");
		}
		for (unsigned int i = 0; i < totNMax() - totNMin() + 1; ++i) {
			if (energyHistogram_lb_[i] > energyHistogram_ub_[i]) {
				throw customException ("Bad energy histogram bound sizes");
			}
			// "standardize" the bounds against U = 0 for to "align" the bins, already done for pkHistogram
			// this way the energy is reported as the limit of the edge of the aligned bins
			if (energyHistogram_lb_[i] < 0) {
				lb = floor((energyHistogram_lb_[i] - 0.0)/energyHistDelta_);
			} else {
				lb = ceil((energyHistogram_lb_[i] - 0.0)/energyHistDelta_);
			}
			if (energyHistogram_ub_[i] < 0) {
				ub = floor((energyHistogram_ub_[i] - 0.0)/energyHistDelta_);
			} else {
				ub = ceil((energyHistogram_ub_[i] - 0.0)/energyHistDelta_);
			}

			try {
				energyHistogram_[i].reinitialize(lb,ub,energyHistDelta_);
			} catch (customException &ce) {
				throw customException ("Unable to reinitialize the energyHistogram");
			}
		}
}

/*!
 * Print the (normalized) energy histogram for each Ntot. netCDF4 not enabled
 *
 * \param [in] fileName Name of the file to print to
 */
void simSystem::printEnergyHistogram (const std::string fileName) {
#ifdef NETCDF_CAPABLE
    throw customException ("Cannot record the energyHistogram for each Ntot in netCDF4 format.");
#else
	// Without netCDF capabilities, just print to ASCII file
	std::ofstream of;
	std::string name = fileName+".dat";
	of.open(name.c_str(), std::ofstream::out);
	of << "# <P(U)> as a function of N_tot." << std::endl;
	of << "# number_of_species: " << nSpecies_ << std::endl;
	of << "# species_total_upper_bound: " << totNBounds_[1] << std::endl;
	of << "# species_total_lower_bound: " << totNBounds_[0] << std::endl;
	double V = box_[0]*box_[1]*box_[2];
	of << "# volume: " << std::setprecision(15) << V << std::endl;
	of << "# Bin widths for each" << std::endl;
	for (std::vector < dynamic_one_dim_histogram >::iterator it = energyHistogram_.begin(); it != energyHistogram_.end(); ++it) {
		of << it->get_delta() << "\t";
	}
	of << std::endl;
	of << "# Bin lower bound for each" << std::endl;
	for (std::vector < dynamic_one_dim_histogram >::iterator it = energyHistogram_.begin(); it != energyHistogram_.end(); ++it) {
		of << it->get_lb() << "\t";
	}
	of << std::endl;
	of << "# Bin upper bound for each" << std::endl;
	for (std::vector < dynamic_one_dim_histogram >::iterator it = energyHistogram_.begin(); it != energyHistogram_.end(); ++it) {
		of << it->get_ub() << "\t";
	}
	of << std::endl;
	of << "# Normalized histogram for each" << std::endl;
	for (std::vector < dynamic_one_dim_histogram >::iterator it = energyHistogram_.begin(); it != energyHistogram_.end(); ++it) {
		std::deque <double> h = it->get_hist();
		double sum = 0.0;
		for (std::deque <double>::iterator it2 = h.begin(); it2 != h.end(); ++it2) {
			sum += *it2;
		}
		for (std::deque <double>::iterator it2 = h.begin(); it2 != h.end(); ++it2) {
			of << *it2/sum << "\t";
		}
		of << std::endl;
	}
	of.close();
#endif
}

/*!
 * Record the particle number histogram for the system at a given Ntot.
 * Only records values when N_tot in range of [min, max].
 */
void simSystem::recordPkHistogram () {
	// only record if in range (removes equilibration stage to get in this range, if there was any)
	if (totN_ >= totNBounds_[0] && totN_ <= totNBounds_[1]) {
		const int address = totN_-totNBounds_[0];
		for (unsigned int i = 0; i < nSpecies_; ++i) {
			pkHistogram_[i][address].record(numSpecies[i]);
		}
	}
}

/*!
 * Check the histogram entries and trim off zero-valued entries and bounds
 */
void simSystem::refinePkHistogramBounds () {
	for (std::vector < std::vector < dynamic_one_dim_histogram > >::iterator it = pkHistogram_.begin(); it != pkHistogram_.end(); ++it) {
		for (std::vector < dynamic_one_dim_histogram >::iterator it2 = it->begin(); it2 != it->end(); ++it2) {
			try {
				it2->trim_edges();
			} catch (customException &ce) {
				throw customException ("Unable to trim edges in pkHistogram at each Ntot");
			}
		}
	}
}

/*!
 * Print the (normalized) particle number histogram for each Ntot. netCDF4 not enabled
 *
 * \param [in] fileName Name of the file to print to
 */
void simSystem::printPkHistogram (const std::string fileName) {
#ifdef NETCDF_CAPABLE
    throw customException ("Cannot record the pkHistogram for each Ntot in netCDF4 format.");
#else
	// Without netCDF capabilities, just print to ASCII file
	for (unsigned int i = 0; i < nSpecies_; ++i) {
		std::ofstream of;
		std::string name = fileName+"_"+sstr(i+1)+".dat";
		of.open(name.c_str(), std::ofstream::out);
		of << "# <P(N_" << i+1 << ")> as a function of N_tot." << std::endl;
		of << "# number_of_species: " << nSpecies_ << std::endl;
		of << "# species_total_upper_bound: " << totNBounds_[1] << std::endl;
		of << "# species_total_lower_bound: " << totNBounds_[0] << std::endl;
		double V = box_[0]*box_[1]*box_[2];
		of << "# volume: " << std::setprecision(15) << V << std::endl;
		of << "# Bin widths for each species index " << std::endl;
		for (std::vector < dynamic_one_dim_histogram >::iterator it = pkHistogram_[i].begin(); it != pkHistogram_[i].end(); ++it) {
			of << it->get_delta() << "\t";
		}
		of << std::endl;
		of << "# Bin lower bound for each species index " << std::endl;
		for (std::vector < dynamic_one_dim_histogram >::iterator it = pkHistogram_[i].begin(); it != pkHistogram_[i].end(); ++it) {
			of << it->get_lb() << "\t";
		}
		of << std::endl;
		of << "# Bin upper bound for each species index " << std::endl;
		for (std::vector < dynamic_one_dim_histogram >::iterator it = pkHistogram_[i].begin(); it != pkHistogram_[i].end(); ++it) {
			of << it->get_ub() << "\t";
		}
		of << std::endl;
		of << "# Normalized histogram for each species index " << std::endl;
		for (std::vector < dynamic_one_dim_histogram >::iterator it = pkHistogram_[i].begin(); it != pkHistogram_[i].end(); ++it) {
			std::deque <double> h = it->get_hist();
			double sum = 0.0;
			for (std::deque <double>::iterator it2 = h.begin(); it2 != h.end(); ++it2) {
				sum += *it2;
			}
			for (std::deque <double>::iterator it2 = h.begin(); it2 != h.end(); ++it2) {
				of << *it2/sum << "\t";
			}
			of << std::endl;
		}
		of.close();
	}
#endif
}

/*!
 * Save the instantaneous energy of the system as a function of the number of particles in the system.
 * Only records values when N_tot in range of [min, max].
 */
/*void simSystem::recordU () {
	// only record if in range (removes equilibration stage to get in this range, if there was any)
	if (totN_ >= totNBounds_[0] && totN_ <= totNBounds_[1]) {
		const int address = totN_-totNBounds_[0];
		AverageU_[address] += energy_;
		numAverageU_[address] += 1.0;
	}
}*/

/*!
 * Print the average energy to file, <U> for every N_tot within range is recorded.  Will overwrite the file if another with that name exists. Prints in netCDF format if enabled.
 *
 * \param [in] fileName Name of the file to print to
 */
/*void simSystem::printU (const std::string fileName) {
	std::vector < double > aveU (AverageU_.size(), 0);
	for (long long int i = 0; i < AverageU_.size(); ++i) {
		aveU[i] = AverageU_[i]/numAverageU_[i];
	}

#ifdef NETCDF_CAPABLE
    	// If netCDF libs are enabled, write to this format
    	const std::string name = fileName + ".nc";
  	NcFile outFile(name.c_str(), NcFile::replace);
	NcDim probDim = outFile.addDim("vectorized_position", aveU.size());
	NcVar probVar = outFile.addVar("averageU", ncDouble, probDim);
	const std::string dummyName = "number_species";
	probVar.putAtt(dummyName.c_str(), sstr(nSpecies_).c_str());
	std::string attName = "species_total_upper_bound";
	probVar.putAtt(attName.c_str(), sstr(totNBounds_[1]).c_str());
	attName = "species_total_lower_bound";
	probVar.putAtt(attName.c_str(), sstr(totNBounds_[0]).c_str());
	attName = "volume";
	double V = box_[0]*box_[1]*box_[2];
	probVar.putAtt(attName.c_str(), sstr(V).c_str());
	probVar.putVar(&aveU[0]);
#else
	// Without netCDF capabilities, just print to ASCII file
	std::ofstream of;
	std::string name = fileName+".dat";
	of.open(name.c_str(), std::ofstream::out);
	of << "# <U> as a function of N_tot." << std::endl;
	of << "# Number of species: " << nSpecies_ << std::endl;
	of << "# species_total_upper_bound: " << totNBounds_[1] << std::endl;
	of << "# species_total_lower_bound: " << totNBounds_[0] << std::endl;
	double V = box_[0]*box_[1]*box_[2];
	of << "# volume: " << std::setprecision(15) << V << std::endl;
	for (long long int i = 0; i < aveU.size(); ++i) {
		of << aveU[i] << std::endl;
	}
	of.close();
#endif
}*/

/*!
 * Add a pair potential to the system which governs the pair (spec1, spec2). However, it only stores the pointer so the object must be
 * fixed in memory somewhere else throughout the simulation.
 *
 * \param [in] spec1 Species index 1 (>= 0)
 * \param [in] spec2 Species index 2 (>= 0)
 * \param [in] ppot_name Name of pair potential
 * \param [in] params Vector of parameters which define pair potential
 * \param [in] bool Optional argument of whether or not to build and maintain a cell list for this pair (spec1, spec2)
 */
 //}, pairPotential *pp, bool useCellList) {
void simSystem::addPotential (const int spec1, const int spec2, const std::string ppot_name, const std::vector < double > &params, const bool useCellList) {
	if (spec1 >= nSpecies_) {
		throw customException ("Trying to define pair potential for species (1) that does not exist yet");
		return;
	}
	if (spec2 >= nSpecies_) {
		throw customException ("Trying to define pair potential for species (2) that does not exist yet");
		return;
	}

	if (ppot_name == "square_well") {
		try {
			ppot[spec1][spec2] = new squareWell;
			ppot[spec1][spec2]->setParameters(params);
			ppot[spec2][spec1] = new squareWell;
			ppot[spec2][spec1]->setParameters(params);
		} catch (customException &ce) {
			std::cerr << ce.what() << std::endl;
			exit(SYS_FAILURE);
		}
	} else if (ppot_name == "lennard_jones") {
		try {
			ppot[spec1][spec2] = new lennardJones;
			ppot[spec1][spec2]->setParameters(params);
			ppot[spec2][spec1] = new lennardJones;
			ppot[spec2][spec1]->setParameters(params);
		} catch (customException &ce) {
			std::cerr << ce.what() << std::endl;
			exit(SYS_FAILURE);
		}
	} else if (ppot_name == "fs_lennard_jones") {
		try {
			ppot[spec1][spec2] = new fsLennardJones;
			ppot[spec1][spec2]->setParameters(params);
			ppot[spec2][spec1] = new fsLennardJones;
			ppot[spec2][spec1]->setParameters(params);
		} catch (customException &ce) {
			std::cerr << ce.what() << std::endl;
			exit(SYS_FAILURE);
		}
	} else if (ppot_name == "hard_sphere") {
		try {
			ppot[spec1][spec2] = new hardCore;
			ppot[spec1][spec2]->setParameters(params);
			ppot[spec2][spec1] = new hardCore;
			ppot[spec2][spec1]->setParameters(params);
		} catch (customException &ce) {
			std::cerr << ce.what() << std::endl;
			exit(SYS_FAILURE);
		}
	} else if (ppot_name == "tabulated") {
		try {
			ppot[spec1][spec2] = new tabulated;
			ppot[spec1][spec2]->setParameters(params);
			ppot[spec2][spec1] = new tabulated;
			ppot[spec2][spec1]->setParameters(params);
		} catch (customException &ce) {
			std::cerr << ce.what() << std::endl;
			exit(SYS_FAILURE);
		}
	} else {
		std::cerr << "Unrecognized pair potential name for species " << spec1 << ", " << spec2 << std::endl;
		exit(SYS_FAILURE);
	}

	ppotSet_[spec1][spec2] = true;
	ppotSet_[spec2][spec1] = true;

	if (useCellList) {
		std::cout << "Setting up cell list for interactions between type " << spec1 << " and " << spec2 << std::endl;
		// add creation of cell lists
		if ((ppot[spec1][spec2]->rcut() > box_[0]/3.0) || (ppot[spec1][spec2]->rcut() > box_[1]/3.0) || (ppot[spec1][spec2] ->rcut() > box_[2]/3.0)) {
			std::cerr << "Cutoff (" << ppot[spec1][spec2]->rcut() << ") larger than 1.0/3.0 boxsize, disabling cell lists for this interaction." << std::endl;
			useCellList_[spec1][spec2] = false;
			useCellList_[spec2][spec1] = false;
		} else {
			std::cout << "Creating Cell list with rcut = " << ppot[spec1][spec2]->rcut() << std::endl;
			useCellList_[spec1][spec2] = true;
			useCellList_[spec2][spec1] = true;

			std::vector <atom*> dummyList(0);

			if (cellListsByPairType_[spec1][spec2] == NULL) {
				cellLists_.push_back(cellList(box_, ppot[spec1][spec2]->rcut(), dummyList));
				cellListsByPairType_[spec1][spec2] = &cellLists_[cellLists_.size()-1];
			}
			if (cellListsByPairType_[spec2][spec1] == NULL) {
				cellLists_.push_back(cellList(box_, ppot[spec2][spec1]->rcut(), dummyList));
				cellListsByPairType_[spec2][spec1] = &cellLists_[cellLists_.size()-1];
			}
		}
	} else {
		useCellList_[spec1][spec2] = false;
		useCellList_[spec2][spec1] = false;
	}
}

/*!
 * Print an XYZ file of the instantaneous system configuration.  This can be read in at a later time via readRestart() function.
 *
 * \param [in] filename File to store XYZ coordinates to
 * \param [in] comment Comment line for the file
 */
void simSystem::printSnapshot (std::string filename, std::string comment) {
	std::ofstream outfile (filename.c_str());

    	int tot = 0;
    	for (unsigned int j = 0; j < nSpecies_; ++j) {
        	tot += numSpecies[j]; // only count fully inserted species
    	}

    	outfile << tot << std::endl;
    	outfile << comment << std::endl;

    	for (unsigned int j = 0; j < nSpecies_; ++j) {
        	long long int num = numSpecies[j];
		if (Mcurrent_ > 1 && fractionalAtomType_ == j) {
        		num += 1; // account for partially inserted atom
		}
		for (unsigned int i = 0; i < num; ++i) {
			if (atoms[j][i].mState == 0) { // only print fully inserted atoms
            			outfile << j << "\t" <<  std::setprecision(15) << atoms[j][i].pos[0] << "\t" << std::setprecision(15) << atoms[j][i].pos[1] << "\t" << std::setprecision(15) << atoms[j][i].pos[2] << std::endl;
        		}
		}
    	}

    	outfile.close();
}

/*!
 * Read an XYZ file as the system's initial configuration.  Note that the number of species, etc. must already be specified in the constructor.
 * Will also reset and calculate the energy from scratch so these potentials should be set before reading in a restart file.
 *
 * \param [in] filename File to read XYZ coordinates from
 */
void simSystem::readRestart (std::string filename) {
	std::ifstream infile (filename.c_str());
	std::string line;
	std::vector < atom > sysatoms;
	std::vector < int > index;
	int natoms = 0;
	int lineIndex = 0;
	while(std::getline(infile,line)) {
		std::stringstream lineStream(line);
		if (lineIndex == 0) {
			lineStream >> natoms;
			index.resize(natoms);
			sysatoms.resize(natoms);
		} else if (lineIndex > 1) {
			lineStream >> index[lineIndex-2] >> sysatoms[lineIndex-2].pos[0] >> sysatoms[lineIndex-2].pos[1] >> sysatoms[lineIndex-2].pos[2];
		}
		lineIndex++;
	}
	infile.close();

	// check if within global bounds
	if (sysatoms.size() > totNBounds_[1] || sysatoms.size() < totNBounds_[0]) {
		throw customException ("Number of particles ("+sstr(sysatoms.size())+") in the restart file out of target range ["+sstr(totNBounds_[0])+", "+sstr(totNBounds_[1])+"]");
	}

	// sort by type
	std::map < int, int > types;
	for (unsigned int j = 0; j < natoms; ++j) {
		if (types.find(index[j]) != types.end()) {
			types[index[j]] += 1;
		} else {
			types[index[j]] = 1;
		}
	}
	int maxType = -1;
	for (std::map<int,int>::iterator it = types.begin(); it != types.end(); ++it) {
		maxType = std::max(maxType, it->first);
		if (it->first < 0 || it->first >= nSpecies_) {
			throw customException ("Restart file corrupted, types out of range");
		}
	}

	// check that pair potentials exist so energy can be calculated
	for (unsigned int i = 0; i < nSpecies_; ++i) {
		for (unsigned int j = 0; j < nSpecies_; ++j) {
			if (!potentialIsSet(i, j)) {
				throw customException("Not all pair potentials are set, so cannot initial from file");
			}
		}
	}

	energy_ = 0.0;

	for (unsigned int j = 0; j < sysatoms.size(); ++j) {
		try {
			// "partially" insert each atom so it goes through all the stages
			insertAtom (index[j], &sysatoms[j]);
			for (unsigned int k = 1; k < Mtot_; ++k) {
				insertAtom (index[j], fractionalAtom_); // this will check that within each species own max and min, global bounds handled above
			}
		} catch (customException &ce) {
			std::string a = "Could not initialize system from restart file, ", b = ce.what();
			throw customException (a+b);
		}
	}

	// recalculate system's initial energy
	energy_ = scratchEnergy();
}

/*!
 * Return the list of neighbors of type A, around a particle of type B which is passed
 *
 * \param [in] typeIndexA Index of first atom type
 * \param [in] typeIndexB Index of second atom type
 * \param [in] atom Pointer to atom to find neighbors around
 *
 * \return neighbor_list
 */
std::vector < atom* > simSystem::getNeighborAtoms(const unsigned int typeIndexA, const unsigned int typeIndexB, atom* _atom) {
	std::vector < atom* > neighbors;

	int end = numSpecies[typeIndexA];
	if (Mcurrent_ > 0 && typeIndexA == fractionalAtomType_) {
		// account for partial atom too
		end++;
	}
	neighbors.reserve(end);

	// if no cell lists are defined for this interaction, return all particles
	if (!useCellList_[typeIndexA][typeIndexB]) {
		for (unsigned int i = 0; i < end; ++i) {
			if (_atom != &atoms[typeIndexA][i]) { // watch out for self in case typeA = typeB
				neighbors.push_back(&atoms[typeIndexA][i]);
			}
		}
	} else if (useCellList_[typeIndexA][typeIndexB]) {
		cellList* cl = cellListsByPairType_[typeIndexA][typeIndexB];
		const unsigned int cellIndex = cl->calcIndex(_atom->pos[0], _atom->pos[1], _atom->pos[2]);

		// loop over own cell
		for (unsigned int i = 0; i < cl->cells[cellIndex].size(); ++i) {
			if (_atom != cl->cells[cellIndex][i]) {
				neighbors.push_back(cl->cells[cellIndex][i]);
			}
		}

		// loop over neighboring cells
		for (unsigned int i = 0; i < cl->neighbours[cellIndex].size(); ++i) {
			const unsigned int neighborCellIndex = cl->neighbours[cellIndex][i];
			for (unsigned int j = 0; j < cl->cells[neighborCellIndex].size(); ++j) {
				if (_atom != cl->cells[neighborCellIndex][j]) {
					neighbors.push_back(cl->cells[neighborCellIndex][j]);
				}
			}
		}
	}

	return neighbors;
}

/*!
 * Recalculate the energy of the system from scratch.
 *
 * \returns totU Total energy of the system
 */
const double simSystem::scratchEnergy () {
   	double totU = 0.0;
   	double V = 1.0;

	for (unsigned int i = 0; i < box_.size(); ++i) {
    		V *= box_[i];
    	}

    	for (unsigned int spec1 = 0; spec1 < nSpecies_; ++spec1) {
        	int num1 = 0, adj1 = 0;
        	try {
            		num1 = numSpecies[spec1];
        	} catch (customException &ce) {
            		std::string a = "Cannot recalculate energy from scratch: ", b = ce.what();
            		throw customException (a+b);
        	}

		// possibly have fractionally inserted atom
		if (fractionalAtomType_ == spec1 && Mcurrent_ > 0) {
			adj1 = 1;
		}

		// wall/barrier interactions
		for (unsigned int j = 0; j < num1+adj1; ++j) {
			double dU = 0.0;
			try {
                        	dU = speciesBarriers[spec1].energy(&atoms[spec1][j], box_);
                        } catch (customException &ce) {
                                std::string a = "Cannot recalculate energy from scratch: ", b = ce.what();
                                throw customException (a+b);
                        }
			if (dU == NUM_INFINITY) {
				return NUM_INFINITY;
			} else {
				totU += dU;
			}
                }

        	// interactions with same type
        	for (unsigned int j = 0; j < num1+adj1; ++j) {
	        	for (unsigned int k = j+1; k < num1+adj1; ++k) {
       				try {
                    			totU += ppot[spec1][spec1]->energy(&atoms[spec1][j], &atoms[spec1][k], box_);
                		} catch (customException &ce) {
                    			std::string a = "Cannot recalculate energy from scratch: ", b = ce.what();
                    			throw customException (a+b);
                		}
            		}
        	}

        	// add tail correction to potential energy but only for atoms fully inserted
#ifdef FLUID_PHASE_SIMULATIONS
        	if ((ppot[spec1][spec1]->useTailCorrection) && (num1 > 1)) {
        		totU += (num1)*0.5*ppot[spec1][spec1]->tailCorrection((num1-1)/V);
        	}
#endif

        	// interactions with other unique types
        	for (unsigned int spec2 = spec1+1; spec2 < nSpecies_; ++spec2) {
            		int num2 = 0, adj2 = 0;
            		try {
                		num2 = numSpecies[spec2];
            		} catch (customException &ce) {
                		std::string a = "Cannot recalculate energy from scratch: ", b = ce.what();
                		throw customException (a+b);
            		}

			if (fractionalAtomType_ == spec2 && Mcurrent_ > 0) {
				adj2 = 1;
			}

            		for (unsigned int j = 0; j < num1+adj1; ++j) {
                		for (unsigned int k = 0; k < num2+adj2; ++k) {
                    			try {
                        			totU += ppot[spec1][spec2]->energy(&atoms[spec1][j], &atoms[spec2][k], box_);
					} catch (customException &ce) {
                        			std::string a = "Cannot recalculate energy from scratch: ", b = ce.what();
                        			throw customException (a+b);
                    			}
                		}
            		}

            		// add tail correction to potential energy but only bewteen fully inserted species
#ifdef FLUID_PHASE_SIMULATIONS
            		if ((ppot[spec1][spec2]->useTailCorrection) && (num2 > 0) && (num1 > 0)) {
                		totU += (num1)*ppot[spec1][spec2]->tailCorrection(num2/V);
            		}
#endif
        	}
    	}

    	if (toggleKE_ == true) {
    		double ns = 0.0;
    		for (unsigned int i = 0; i < nSpecies_; ++i) {
    			ns += numSpecies[i];
    		}
    		totU += 1.5/beta_*ns; // only adjust for FULLY-INSERTED ATOMS
    	}

    	return totU;
}

/*!
 * Returns the absolute maximum number of a given species type allowed in the system.
 *
 * \param [in] index  Species index to query
 *
 * \return maxSpecies Maximum number of them allowed
 */
const int simSystem::maxSpecies (const int index) {
	if (maxSpecies_.begin() == maxSpecies_.end()) {
        	throw customException ("No species in the system, cannot report a maximum");
    	}
    	if (maxSpecies_.size() <= index) {
        	throw customException ("System does not contain that species, cannot report a maximum");
    	} else  {
        	return maxSpecies_[index];
    	}
}

/*!
 * Returns the absolute minimum number of a given species type allowed in the system.
 *
 * \param [in] index  Species index to query
 *
 * \return minSpecies Minimum number of them allowed
 */
const int simSystem::minSpecies (const int index) {
	if (minSpecies_.begin() == minSpecies_.end()) {
        	throw customException ("No species in the system, cannot report a minimum");
    	}
    	if (minSpecies_.size() <= index) {
        	throw customException ("System does not contain that species, cannot report a minimum");
    	} else  {
        	return minSpecies_[index];
    	}
}

/*!
 * Return a pointer to the TMMC biasing object, if using TMMC, else throws an exception.
 *
 * \return tmmc Pointer to TMMC biasing object being used.
 */
tmmc* simSystem::getTMMCBias () {
	if (useTMMC == true) {
		return tmmcBias;
	} else {
		throw customException ("Not using TMMC");
	}
}

/*!
 * Return a pointer to the TMMC biasing object, if using TMMC, else throws an exception.
 *
 * \return wala Pointer to WALA biasing object being used.
 */
wala* simSystem::getWALABias () {
	if (useWALA == true) {
		return wlBias;
	} else {
		throw customException ("Not using WALA");
	}
}

/*
 * Start using a Wang-Landau bias in the simulation. Throws an exception if input values are illegal or there is another problem (e.g. memory).
 *
 * \param [in] lnF Factor by which the estimate of the density of states in updated each time it is visited.
 * \param [in] g Factor by which lnF is reduced (multiplied) once "flatness" has been achieved.
 * \param [in] s Factor by which the min(H) must be within the mean of H to be considered "flat", e.g. 0.8 --> min is within 20% error of mean
 * \param [in] Mtot Total number of expanded ensemble state allowed within the system
 */
void simSystem::startWALA (const double lnF, const double g, const double s, const int Mtot) {
	// initialize the wala object
	try {
		wlBias = new wala (lnF, g, s, totNBounds_[1], totNBounds_[0], Mtot, box_);
	} catch (customException& ce) {
		throw customException ("Cannot start Wang-Landau biasing in system: "+sstr(ce.what()));
	}

	useWALA = true;
}

/*!
 * Start using a transition-matrix in the simulation. Throws an exception if input values are illegal or there is another problem (e.g. memory).
 *
 * \param [in] tmmcSweepSize Number of times each transition in the collection matrix must be visited for a "sweep" to be completed
 * \param [in] Mtot Total number of expanded ensemble state allowed within the system
 */
void simSystem::startTMMC (const long long int tmmcSweepSize, const int Mtot) {
	// initialize the tmmc object
	try {
		tmmcBias = new tmmc (totNBounds_[1], totNBounds_[0], Mtot, tmmcSweepSize, box_);
	} catch (customException& ce) {
		throw customException ("Cannot start TMMC biasing in system: "+sstr(ce.what()));
	}

	useTMMC = true;
}

/*!
 * Calculate the bias based on a systems current state and the next state being proposed.
 *
 * \param [in] sys System object containing the current state of the system.
 * \param [in] nTotFinal Total atoms in the proposed final state.
 * \param [in] mFinal Final value of the expanded ensemble state of the system.
 * \param [in] p_u Ratio of the system's partition in the final and initial state (e.g. unbiased p_acc = min(1, p_u)).
 *
 * \return rel_bias The value of the relative bias to apply in the metropolis criteria during sampling
 */
const double calculateBias (simSystem &sys, const int nTotFinal, const int mFinal) {
	double rel_bias = 1.0;

	if (sys.useTMMC && !sys.useWALA) {
		// TMMC biasing
		const __BIAS_INT_TYPE__ address1 = sys.tmmcBias->getAddress(sys.getTotN(), sys.getCurrentM()), address2 = sys.tmmcBias->getAddress(nTotFinal, mFinal);
		const double b1 = sys.tmmcBias->getBias (address1), b2 = sys.tmmcBias->getBias (address2);
		rel_bias = exp(b2-b1);
    	} else if (!sys.useTMMC && sys.useWALA) {
    		// Wang-Landau Biasing
    		const __BIAS_INT_TYPE__ address1 = sys.wlBias->getAddress(sys.getTotN(), sys.getCurrentM()), address2 = sys.wlBias->getAddress(nTotFinal, mFinal);
    		const double b1 = sys.wlBias->getBias (address1), b2 = sys.wlBias->getBias (address2);
    		rel_bias = exp(b2-b1);
    	} else if (sys.useTMMC && sys.useWALA) {
    		// Crossover phase where we use WL but update TMMC collection matrix
	    	const int address1 = sys.wlBias->getAddress(sys.getTotN(), sys.getCurrentM()), address2 = sys.wlBias->getAddress(nTotFinal, mFinal);
    		const double b1 = sys.wlBias->getBias (address1), b2 = sys.wlBias->getBias (address2);
    		rel_bias = exp(b2-b1);
    	} else {
    		// No biasing
    		rel_bias = 1.0;
    	}

	return rel_bias;
}
