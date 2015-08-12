#include "system.h"

// netCDF if enable
#ifdef NETCDF_CAPABLE
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;
#endif

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
    try {
    	AverageU_.resize(size, 0); 
    } catch (std::bad_alloc &ba) {
    	throw customException ("Out of memory for energy record");
    }
    try {
    	numAverageU_.resize(size, 0);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory for energy record");
    }
}

/*!
 * Insert an atom into the system. Does all the bookkeepping behind the scenes.
 * 
 * \param [in] typeIndex What type the atom is (>= 0)
 * \param [in] newAtom Pointer to new atom.  A copy is stored in the system so the original may be destroyed.
 */
void simSystem::insertAtom (const int typeIndex, atom *newAtom) {
    if (typeIndex < nSpecies_ && typeIndex >= 0) {
        if (numSpecies[typeIndex] < maxSpecies_[typeIndex]) {
            atoms[typeIndex][numSpecies[typeIndex]] = (*newAtom);
            numSpecies[typeIndex]++;
            totN_++;
           // add particle into appropriate cell list
           for (unsigned int i=0; i<nSpecies_; i++)
           {
           		if (useCellList_[typeIndex][i])
            	{
            		cellList* cl = cellListsByPairType_[typeIndex][i];
            		cl->insertParticle(&atoms[typeIndex][numSpecies[typeIndex]- 1]);
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
        if ((numSpecies[typeIndex] > minSpecies_[typeIndex]) || override) {
        	// delete particle from appropriate cell list
            for (unsigned int i=0; i<nSpecies_; i++)
            {
            	if (useCellList_[typeIndex][i])
            	{
            		cellList* cl = cellListsByPairType_[typeIndex][i];
            		cl->swapAndDeleteParticle(&atoms[typeIndex][atomIndex], &atoms[typeIndex][numSpecies[typeIndex] - 1]);
            	}
            }
        
            atoms[typeIndex][atomIndex] = atoms[typeIndex][numSpecies[typeIndex] - 1];    // "replacement" operation
            numSpecies[typeIndex]--;
            totN_--;
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
            for (unsigned int i=0; i<nSpecies_; i++)
            {
            	if (useCellList_[typeIndex][i])
            	{
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
}

/*!
 * Initialize the system. Sets the use of both WL and TMMC biasing to false.
 * 
 * \param [in] nSpecies Number of unqiue species types to allow in the system
 * \param [in] beta Inverse temperature (1/kT)
 * \param [in] box Box dimensions [x, y, z]
 * \param [in] mu Chemical potential of each species
 * \param [in] maxSpecies Maximum number of each species to allow in the system
 */
simSystem::simSystem (const unsigned int nSpecies, const double beta, const std::vector < double > box, const std::vector < double > mu, const std::vector < int > maxSpecies, const std::vector < int > minSpecies) {
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
    try {
        numAverageU_.resize(size, 0);
    } catch (std::bad_alloc &ba) {
     	throw customException ("Out of memory for energy record");
    }
    try {
        AverageU_.resize(size, 0);
    } catch (std::bad_alloc &ba) {
        throw customException ("Out of memory for energy record");
    }
}

/*!
 * Save the instantaneous energy of the system as a function of the number of particles in the system.
 * Only records values when N_tot in range of [min, max].
 */
void simSystem::recordU () {
	// only record if in range (removes equilibration stage to get in this range, if there was any)
	if (totN_ >= totNBounds_[0] && totN_ <= totNBounds_[1]) {
		const int address = totN_-totNBounds_[0];
		AverageU_[address] += energy_;
		numAverageU_[address] += 1.0;
	}
}

/*!
 * Print the average energy to file, <U> for every N_tot within range is recorded.  Will overwrite the file if another with that name exists. Prints in netCDF format if enabled.
 * 
 * \param [in] fileName Name of the file to print to
 */
void simSystem::printU (const std::string fileName) {
	std::vector < double > aveU (AverageU_.size(), 0);
	for (long long int i = 0; i < AverageU_.size(); ++i) {
		aveU[i] = AverageU_[i]/numAverageU_[i]; 
	}
	
#ifdef NETCDF_CAPABLE
    // If netCDF libs are enabled, write to this format
    const std::string name = fileName + ".nc";
  	NcFile outFile(name.c_str(), NcFile::replace);
	NcDim probDim = outFile.addDim("vectorized_position", aveU.size());
	NcVar probVar = outFile.addVar("<U>", ncDouble, probDim);
	const std::string dummyName = "number_species:";
	probVar.putAtt(dummyName.c_str(), sstr(nSpecies_).c_str());
	const std::string attName = "species_total_upper_bound:";
	probVar.putAtt(attName.c_str(), sstr(totNBounds_[1]).c_str());
	const std::string attName = "species_total_lower_bound:";
	probVar.putAtt(attName.c_str(), sstr(totNBounds_[0]).c_str());
	probVar.putVar(&aveU[0]);
#else
	// Without netCDF capabilities, just print to ASCII file
	std::ofstream of;
	of.open(fileName+".dat", std::ofstream::out);
	of << "# <U> as a function of N_tot." << std::endl;
	of << "# Number of species:" << nSpecies_ << std::endl;
	of << "# species_total_upper_bound:" << totNBounds_[1] << std::endl;
	of << "# species_total_lower_bound:" << totNBounds_[0] << std::endl;
	for (long long int i = 0; i < aveU.size(); ++i) {
		of << aveU[i] << std::endl;
	}
	of.close();
#endif
}

/*!
 * Add a pair potential to the system which governs the pair (spec1, spec2). However, it only stores the pointer so the object must be 
 * fixed in memory somewhere else throughout the simulation.
 * 
 * \param [in] spec1 Species index 1 (>= 0)
 * \param [in] spec2 Species index 2 (>= 0)
 * \param [in] pp Pointer to the pair potential to add.  Will be stored as a pointer, so original cannot be destroyed in memory.
 * \param [in] bool Optional argument of whether or not to build and maintain a cell list for this pair (spec1, spec2)
 */
void simSystem::addPotential (const int spec1, const int spec2, pairPotential *pp, bool useCellList) {
	if (spec1 >= nSpecies_) {
		throw customException ("Trying to define pair potential for species (1) that does not exist yet");
		return;
	}
	if (spec2 >= nSpecies_) {
		throw customException ("Trying to define pair potential for species (2) that does not exist yet");
		return;
	}
	ppot[spec1][spec2] = pp;
	ppot[spec2][spec1] = pp;
	ppotSet_[spec1][spec2] = true;
	ppotSet_[spec2][spec1] = true;
	
	if (useCellList)
	{
		std::cout<<"Setting up cell list for interactions between type "<<spec1<<" and "<<spec2<<std::endl;
		// add creation of cell lists
		if ((pp->rcut() > box_[0]/3.0) || (pp->rcut() > box_[1]/3.0) || (pp->rcut() > box_[2]/3.0))
		{
			std::cerr<<"Cutoff ("<<pp->rcut()<<") larger than 1.0/3.0 boxsize, disabling cell lists for this interaction."<<std::endl;
			useCellList_[spec1][spec2] = false;
			useCellList_[spec2][spec1] = false;
		}
		else
		{
			std::cout<<"Creating Cell list with rcut="<<pp->rcut()<<std::endl;
			useCellList_[spec1][spec2] = true;
			useCellList_[spec2][spec1] = true;
			
			std::vector<atom*> dummyList(0);
			
			if (cellListsByPairType_[spec1][spec2] == NULL)
			{
				cellLists_.push_back(cellList(box_, pp->rcut(), dummyList));
				cellListsByPairType_[spec1][spec2] = &cellLists_[cellLists_.size()-1];
			}
			if (cellListsByPairType_[spec2][spec1] == NULL)
			{
				cellLists_.push_back(cellList(box_, pp->rcut(), dummyList));
				cellListsByPairType_[spec2][spec1] = &cellLists_[cellLists_.size()-1];
			}
		}
	}
	else
	{
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
        tot += numSpecies[j];
    }
    
    outfile << tot << std::endl;
    outfile << comment << std::endl;
    
    for (unsigned int j = 0; j < nSpecies_; ++j) {
        const int num = numSpecies[j];
        for (unsigned int i = 0; i < num; ++i) {
            outfile << j << "\t" << atoms[j][i].pos[0] << "\t" << atoms[j][i].pos[1] << "\t" << atoms[j][i].pos[2] << std::endl;
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
		throw customException ("Number of particles in the restart file out of target range");
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
			insertAtom (index[j], &sysatoms[j]); // this will check that within each species own max and min, global bounds handled above
		}
		catch (customException &ce) {
			std::string a = "Could not initialize system from restart file, ", b = ce.what();
			throw customException (a+b);
		}
	}
	
	// recalculate system's initial energy
	energy_ = scratchEnergy();
}

/*!
 * Return the list of neighbors of type A, for a particle of type B at position pos
 * 
 * \param [in] typeIndexA Index of first atom type
 * \param [in] typeIndexB Index of second atom type
 * \param [in] atom Pointer to atom to find neighbors around
 * 
 * \return neighbor_list
 */
std::vector< std::vector<double> > simSystem::getNeighborPositions(const unsigned int typeIndexA, const unsigned int typeIndexB, atom* _atom)
{
	std::vector< std::vector<double> > neighbors;
	neighbors.reserve(numSpecies[typeIndexA]);
	
	// if no cell lists are defined for this interaction, return all particles
	if (!useCellList_[typeIndexA][typeIndexB])
	{
		for (unsigned int i=0; i<numSpecies[typeIndexA]; i++)
		{
			if (_atom != &atoms[typeIndexA][i])
			{
				neighbors.push_back(atoms[typeIndexA][i].pos);
			}
		}
	}
	else if (useCellList_[typeIndexA][typeIndexB])
	{
		cellList* cl = cellListsByPairType_[typeIndexA][typeIndexB];
		const unsigned int cellIndex = cl->calcIndex(_atom->pos[0], _atom->pos[1], _atom->pos[2]);

		// loop over own cell
		for (unsigned int i=0; i<cl->cells[cellIndex].size(); i++)
		{
			if (_atom != cl->cells[cellIndex][i])
			{
				neighbors.push_back(cl->cells[cellIndex][i]->pos);
			}
		}
		// loop over neighboring cells
		for (unsigned int i=0; i<cl->neighbours[cellIndex].size(); i++)
		{
			const unsigned int neighborCellIndex = cl->neighbours[cellIndex][i];
			
			for (unsigned int j=0; j<cl->cells[neighborCellIndex].size(); j++)
			{
				if (_atom != cl->cells[neighborCellIndex][j])
				{
					neighbors.push_back(cl->cells[neighborCellIndex][j]->pos);
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
        int num1;
        try {
            num1 = numSpecies[spec1];
        } catch (customException &ce) {
            std::string a = "Cannot recalculate energy from scratch: ", b = ce.what();
            throw customException (a+b);
        }
		
        // interactions with same type
        for (unsigned int j = 0; j < num1; ++j) {
            for (unsigned int k = j+1; k < num1; ++k) {
                try {
                    totU += ppot[spec1][spec1]->energy(atoms[spec1][j].pos, atoms[spec1][k].pos, box_); 
                } catch (customException &ce) {
                    std::string a = "Cannot recalculate energy from scratch: ", b = ce.what();
                    throw customException (a+b);
                }
            }
        }
        
        // add tail correction to potential energy
#ifdef FLUID_PHASE_SIMULATIONS
        if ((ppot[spec1][spec1]->useTailCorrection) && (num1 > 1)) {
        	totU += (num1)*0.5*ppot[spec1][spec1]->tailCorrection((num1-1)/V);
        }
#endif        
        // interactions with other types
        for (unsigned int spec2 = 0; spec2 < nSpecies_; ++spec2) {
            if (spec2 > spec1) { // only compute unique interactions
                int num2;
                try {
                    num2 = numSpecies[spec2];
                } catch (customException &ce) {
                    std::string a = "Cannot recalculate energy from scratch: ", b = ce.what();
                    throw customException (a+b);
                }
                
                for (unsigned int j = 0; j < num1; ++j) {
                    for (unsigned int k = 0; k < num2; ++k) {
                        try {
                            totU += ppot[spec1][spec2]->energy(atoms[spec1][j].pos, atoms[spec2][k].pos, box_);
                        } catch (customException &ce) {
                            std::string a = "Cannot recalculate energy from scratch: ", b = ce.what();
                            throw customException (a+b);
                        }
                    }
                }
                // add tail correction to potential energy
#ifdef FLUID_PHASE_SIMULATIONS
                if ((ppot[spec1][spec2]->useTailCorrection) && (num2 > 0) && (num1 > 0)) {
                	totU += (num1)*ppot[spec1][spec2]->tailCorrection(num2/V);
        		}
#endif
            }
        }
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
 * \param [in] Nmax Upper bound for total number of particles in the system.
 * \param [in] Nmin Lower bound for total number of particles in the system. 
 */
void simSystem::startWALA (const double lnF, const double g, const double s, const int Nmax, const int Nmin) { 
	// initialize the wala object
	try {
		wlBias = new wala (lnF, g, s, Nmax, Nmin);
	} catch (customException& ce) {
		throw customException ("Cannot start Wang-Landau biasing in system: "+sstr(ce.what()));
	}
	
	useWALA = true;
}

/*!
 * Start using a transition-matrix in the simulation. Throws an exception if input values are illegal or there is another problem (e.g. memory).
 * 
 * \param [in] Nmax Upper bound for total number of particles in the system.
 * \param [in] Nmin Lower bound for total number of particles in the system.
 */
void simSystem::startTMMC (const int Nmax, const int Nmin) { 
	// initialize the tmmc object
	try {
		tmmcBias = new tmmc (Nmax, Nmin);
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
 * \param [in] p_u Ratio of the system's partition in the final and initial state (e.g. unbiased p_acc = min(1, p_u)).
 * 
 * \return rel_bias The value of the relative bias to apply in the metropolis criteria during sampling
 */
const double calculateBias (simSystem &sys, const int nTotFinal) { //, const double p_u) {
	double rel_bias = 1.0;
	
	if (sys.useTMMC && !sys.useWALA) {
		// TMMC biasing
		const __BIAS_INT_TYPE__ address1 = sys.tmmcBias->getAddress(sys.getTotN()), address2 = sys.tmmcBias->getAddress(nTotFinal);
		const double b1 = sys.tmmcBias->getBias (address1), b2 = sys.tmmcBias->getBias (address2);
		rel_bias = exp(b2-b1);
    } else if (!sys.useTMMC && sys.useWALA) {
    	// Wang-Landau Biasing
    	const __BIAS_INT_TYPE__ address1 = sys.wlBias->getAddress(sys.getTotN()), address2 = sys.wlBias->getAddress(nTotFinal); 
    	const double b1 = sys.wlBias->getBias (address1), b2 = sys.wlBias->getBias (address2);
    	rel_bias = exp(b2-b1);
    } else if (sys.useTMMC && sys.useWALA) {
    	// Crossover phase where we use WL but update TMMC collection matrix
    	const int address1 = sys.wlBias->getAddress(sys.getTotN()), address2 = sys.wlBias->getAddress(nTotFinal);
    	const double b1 = sys.wlBias->getBias (address1), b2 = sys.wlBias->getBias (address2);
    	rel_bias = exp(b2-b1);
    } else {
    	// No biasing
    	rel_bias = 1.0;
    }
	
	return rel_bias;
}