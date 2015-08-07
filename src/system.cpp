#include "system.h"
#include "global.h"
#include <stdlib.h>
#include <exception>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>

/*!
 * Insert an atom into the system. Does all the bookkeepping behind the scenes.
 * \param [in] typeIndex What type the atom is (>= 0)
 * \param [in] newAtom Pointer to new atom.  A copy is stored in the system so the original may be destroyed.
 */
void simSystem::insertAtom (const int typeIndex, atom *newAtom) {
    if (typeIndex < nSpecies_) {
        if (numSpecies[typeIndex] < maxSpecies_[typeIndex]) {
            atoms[typeIndex][numSpecies[typeIndex]] = (*newAtom);
            numSpecies[typeIndex]++;
            
           // add particle into appropriate cell list
           for (unsigned int i=0; i<nSpecies_; i++)
           {
           		if (useCellList_[typeIndex][i])
            	{
            		cellList* cl = cellListsByPairType_[typeIndex][i];
            		cl->insertParticle(&atoms[typeIndex][numSpecies[typeIndex]-1]);
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
 * \param [in] typeIndex What type the atom is (>= 0)
 * \param [in] atomIndex Which atom index of type typeIndex to destroy (>= 0)
 */
void simSystem::deleteAtom (const int typeIndex, const int atomIndex) {
    if (typeIndex < nSpecies_) {
        if (numSpecies[typeIndex] > 0) {
        
        	// delete particle from appropriate cell list
            for (unsigned int i=0; i<nSpecies_; i++)
            {
            	if (useCellList_[typeIndex][i])
            	{
            		cellList* cl = cellListsByPairType_[typeIndex][i];
            		cl->swapAndDeleteParticle(&atoms[typeIndex][atomIndex], &atoms[typeIndex][numSpecies[typeIndex]-1]);
            	}
            }
        
            atoms[typeIndex][atomIndex] = atoms[typeIndex][numSpecies[typeIndex]-1];    // "replacement" operation
            numSpecies[typeIndex]--;
        } else {
            std::string index = static_cast<std::ostringstream*>( &(std::ostringstream() << typeIndex) )->str();
            throw customException ("No atoms left in system, cannot delete an atom of type index "+index);
        }
    } else {
        throw customException ("That species index does not exist, cannot delete an atom");
    }
}

/*!
 * Translate an atom in the system. Does all the bookkeepping behind the scenes.
 * Do nothing if there is no cell list defined for the type
 * \param [in] typeIndex What type the atom is (>= 0)
 * \param [in] atomIndex Which atom index of type typeIndex to translate (>= 0)
 */
void simSystem::translateAtom (const int typeIndex, const int atomIndex, std::vector<double> oldPos) {
    if (typeIndex < nSpecies_) {
        if (numSpecies[typeIndex] > 0) {
        
        	// delete particle from appropriate cell list
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
            throw customException ("No atoms left in system, cannot translate an atom of type index "+index);
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
		delete [] tmmcBias;
	}
	if (useWALA) {
		delete [] wlBias;
	}
}

/*!
 * Initialize the system. Sets the use of both WL and TMMC biasing to false.
 * \param [in] nSpecies Number of unqiue species types to allow in the system
 * \param [in] beta Inverse temperature (1/kT)
 * \param [in] box Box dimensions [x, y, z]
 * \param [in] mu Chemical potential of each species
 * \param [in] maxSpecies Maximum number of each species to allow in the system
 */
simSystem::simSystem (const unsigned int nSpecies, const double beta, const std::vector < double > box, const std::vector < double > mu, const std::vector < int > maxSpecies) {
	if ((box.size() != 3) || (nSpecies != mu.size()) || (maxSpecies.size() != nSpecies)) {
		throw customException ("Invalid system initialization parameters");
		exit(SYS_FAILURE);
	} else {
		nSpecies_ = nSpecies;
        maxSpecies_ = maxSpecies;
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
        if (maxSpecies_[i] <= 0) {
            throw customException ("Max species <= 0");
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
}

/*!
 * Add a pair potential to the system which governs the pair (spec1, spec2).
 * \param [in] spec1 Species index 1 (>= 0)
 * \param [in] spec2 Species index 2 (>= 0)
 * \param [in] pp Pointer to the pair potential to add.  Will be stored as a pointer, so original cannot be destroyed in memory.
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
 * Print an XYZ file of the instantaneous system configuration.
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
            outfile << "A" << j+1 << "\t" << atoms[j][i].pos[0] << "\t" << atoms[j][i].pos[1] << "\t" << atoms[j][i].pos[2] << std::endl;
        }
    }
    
    outfile.close();
}

/*!
 * Read an XYZ file as the system's initial configuration.  Note that the number of species, etc. must already be specified in the constructor.
 * \param [in] filename File to read XYZ coordinates from
 */
void simSystem::readRestart (std::string filename) {
	std::ifstream infile (filename.c_str());
	std::string line;
	int natoms = 0;
	infile >> natoms;
	std::getline(infile, line); // comment line 
	
	std::vector < atom > sysatoms (natoms);
	std::vector < int > index (natoms);
	std::map < int, int > types;
	for (unsigned int j = 0; j < natoms; ++j) {
        infile >> index[j] >> sysatoms[j].pos[0] >> sysatoms[j].pos[1] >> sysatoms[j].pos[2];
    }
	infile.close();
	
	// sort by type
	for (unsigned int j = 0; j < natoms; ++j) {
		if (types.find(index[j]) != types.end()) {
			types[index[j]] += 1;
		} else {
			types[index[j]] = 1;
		}
	}
	int maxType = -1;
	for (std::map<int,int>::iterator it=types.begin(); it != types.end(); ++it) {
		maxType = std::max(maxType, it->first);
		if (it->first < 0 || it->first >= nSpecies_) {
			throw customException ("Restart file corrupted, types out of range");
		}
	}
	
	// ensure system is empty
	for (unsigned int j = 0; j < nSpecies_; ++j) {
		atoms[j].resize(0);
	}
	
	for (unsigned int j = 0; j < sysatoms.size(); ++j) {
		try {
			insertAtom (index[j], &sysatoms[j]);
		}
		catch (customException &ce) {
			std::string a = "Could not initialize system from restart file, ", b = ce.what();
			throw customException (a+b);
		}
	}
	
	// double check
	for (unsigned int j = 0; j < nSpecies_; ++j) {
		if (atoms[j].begin() != atoms[j].end()) {
			if (atoms[j].size() != types[j]) {
				throw customException ("Failed to properly insert old atoms into system");
			}
		} else {
			if (0 != types[j]) {
				throw customException ("Failed to properly insert old atoms into system");
			}
		}
	}
}

/*!
 * Return the list of neighbors of type A, for a particle of type B at position pos
 * \returns neighbors list
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
        if ((ppot[spec1][spec1]->useTailCorrection) && (num1 > 1)) {
        	totU += (num1)*0.5*ppot[spec1][spec1]->tailCorrection((num1-1)/V);
        }
        
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
                if ((ppot[spec1][spec2]->useTailCorrection) && (num2 > 0) && (num1 > 0)) {
                	totU += (num1)*ppot[spec1][spec2]->tailCorrection(num2/V);
        		}

            }
        }
    }
    
    return totU;
}

/*!
 * Returns the absolute maximum number of a given species type allowed in the system.
 * \param [in] index  Species index to query
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
 * Return a pointer to the TMMC biasing object, if using TMMC, else throws an exception.
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
 * \param [in] Nmax Vector of upper bound for number of particles of each species.
 * \param [in] Nmin Vector of lower bound for number of particles of each species. 
 */
void simSystem::startWALA (const double lnF, const double g, const double s, const std::vector <int> &Nmax, const std::vector <int> &Nmin) { 
	// initialize the wala object
	try {
		wlBias = new wala (lnF, g, s, nSpecies_, Nmax, Nmin);
	} catch (customException& ce) {
		throw customException ("Cannot start Wang-Landau biasing in system: "+sstr(ce.what()));
	}
	
	useWALA = true;
}

/*!
 * Start using a transition-matrix in the simulation. Throws an exception if input values are illegal or there is another problem (e.g. memory).
 * 
 * \param [in] Nmax Vector of upper bound for number of particles of each species.
 * \param [in] Nmin Vector of lower bound for number of particles of each species.
 */
void simSystem::startTMMC (const std::vector <int> &Nmax, const std::vector <int> &Nmin) { 
	// initialize the tmmc object
	try {
		tmmcBias = new tmmc (nSpecies_, Nmax, Nmin);
	} catch (customException& ce) {
		throw customException ("Cannot start TMMC biasing in system: "+sstr(ce.what()));
	}
		
	useTMMC = true; 
}
