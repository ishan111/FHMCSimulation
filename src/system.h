#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <exception>
#include <map>
#include <cmath>
#include "potentials.h"
#include "atom.h"
#include "cellList.h"
#include "bias.h"
#include "global.h"

// netCDF if enabled
#ifdef NETCDF_CAPABLE
#include <netcdf>
using namespace netCDF;
using namespace netCDF::exceptions;
#endif

/*!
 * System information for the simulation.
 */
class simSystem {
public:
	simSystem (const unsigned int nSpecies, const double beta, const std::vector < double > box, const std::vector < double > mu, const std::vector < int > maxSpecies, const std::vector < int > minSpecies);
	~simSystem ();
    
	void incrementEnergy (const double dU) { energy_ += dU; } //!< Increment the system's energy
	void addPotential (const int spec1, const int spec2, pairPotential *pp, bool useCellList=false);
	void printSnapshot (std::string filename, std::string comment);
	void insertAtom (const int typeIndex, atom *newAtom);
    void deleteAtom (const int typeIndex, const int atomIndex, bool override=false);
    void translateAtom (const int typeIndex, const int atomIndex, std::vector<double> oldPos);
    void readRestart (std::string filename);
    void recordU ();
    void printU (const std::string fileName);
    void startWALA (const double lnF, const double g, const double s); //!< Start using Wang-Landau and instantiate the bias object
   	void stopWALA () { useWALA = false; delete wlBias; } //!< Stop using Wang-Landau and free the bias object
   	void startTMMC (); //!< Start using TMMC and instantiate the bias object
   	void stopTMMC () { useTMMC = false; delete tmmcBias; } //!< Stop using TMMC and free the bias object
    void setTotNBounds (const std::vector < int > &bounds);
   	bool potentialIsSet (const int spec1, const int spec2) { return ppotSet_[spec1][spec2]; }	//!< Boolean which returns whether or not a pair has had its potential specified by the user yet
    const int nSpecies () { return nSpecies_; } //!< Return the number of different species in the system
    const int maxSpecies (const int index);
    const int minSpecies (const int index);
    const int totNMax () { return totNBounds_[1]; } //!< Return upper bound on the total number of atoms in the system
    const int totNMin () { return totNBounds_[0]; } //!< Return lower bound on the total number of atoms in the system
    const int getTotN () { return totN_; } //!< Return a sum of the total number of atoms currently in the system
    const double energy () { return energy_; } //!< Return the system's instantaneous energy
    const double scratchEnergy ();
    const double beta () { return beta_; } //!< Return 1/kT
    const double mu (const int index) { return mu_[index]; } //!< Return the chemical potential for a given species' index
    const std::vector < double > box () { return box_; } //!< Return the system box dimensions
    std::vector< std::vector <double> > getNeighborPositions(const unsigned int typeIndexA, const unsigned int typeIndexB, atom* _atom);
    tmmc* getTMMCBias (); //!< Return pointer to the TMMC bias
    wala* getWALABias (); //!< Return pointer to the Wang-Landau bias	
    
    bool useTMMC; //!< Logical stating whether or not to use TMMC biasing
    bool useWALA; //!< Logical stating whether or not to use Wang-Landau biasing
    tmmc* tmmcBias; //!< TMMC biasing function
    wala* wlBias; //!< WL biasing function
    std::vector < int > numSpecies;		//!< Total number of each type of atom the system contains
    std::vector < std::vector < atom > > atoms;	//!< Atoms in a matrix by type, and particle index, respectively that a system CAN hold but not all are actually "in" the system
    std::vector < std::vector < pairPotential* > > ppot;	//!< Matrix of pair potentials for atom types i, j

private:
    int nSpecies_; //!< Number of species types allowed in the simulation (single component = 1, multicomponent > 1)
    int totN_; //!< Sum total of all atoms in the system
    double beta_; //!< Inverse temperature, really 1/kT
    double energy_; //!< Instantaneous energy of the system
    std::vector < int > totNBounds_; //!< For multicomponent mixtures, biases use Shen and Errington method which uses bounds on the total number of particles in the system
    std::vector < int > maxSpecies_; //!< Maximum number of each species allowed in the simulation
    std::vector < int > minSpecies_; //!< Minimum number of each species allowed in the simulation
	std::vector < double > box_; //!< System box size
	std::vector < double > mu_; //!< Chemical potential of each species
	std::vector < long double > AverageU_; //!< Accumulator for U for each N_tot macrostate observed between the [min, max] ranges specified
	std::vector < long double > numAverageU_; //!< Number of times U has been recorded at each macrostate
	std::vector < std::vector < bool > > ppotSet_; //!< Matrix of pair potentials between type i and j
    std::vector < std::vector < bool > > useCellList_;  //!< Matrix of whether or not to use cell lists to calculate potentials for pair type (i,j)
    std::vector <cellList> cellLists_; // this vector stores the actual cell lists for the inserted potentials
    std::vector < std::vector <cellList*> > cellListsByPairType_; // this matrix stores pointers to the actual cell lists for all pair types
};

const double calculateBias (simSystem &sys, const int nTotFinal);

#endif