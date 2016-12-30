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
#include <iomanip>
#include <memory>
#include "potentials.h"
#include "atom.h"
#include "cellList.h"
#include "bias.h"
#include "global.h"
#include "histogram.h"
#include "barrier.h"

/*!
 * System information for the simulation.
 */
class simSystem {
public:
	simSystem () {;}
	simSystem (const unsigned int nSpecies, const double beta, const std::vector < double > box, const std::vector < double > mu, const std::vector < int > maxSpecies, const std::vector < int > minSpecies, const int Mtot, const double energyHistDelta = 10.0, const int max_order = 2);
	~simSystem ();

    bool addKECorrection () { return toggleKE_; }
    void toggleKE ();
	void incrementEnergy (const double dU) { energy_ += dU; } //!< Increment the system's energy
	void addPotential (const int spec1, const int spec2, const std::string ppot_name, const std::vector < double > &params, const bool useCellList=false, const std::string tabFile=""); //(const int spec1, const int spec2, pairPotential *pp, bool useCellList=false);
	void printSnapshot (std::string filename, std::string comment, bool overwrite=true);
	void insertAtom (const int typeIndex, atom *newAtom, bool override=false);
	void deleteAtom (const int typeIndex, const int atomIndex, bool override=false);
	void translateAtom (const int typeIndex, const int atomIndex, std::vector<double> oldPos);
   	void readConfig (std::string filename);

	void setEUB (const std::vector < double > eub) { energyHistogram_ub_ = eub; } //!< Assign energy upper bound from restart
	void setELB (const std::vector < double > elb) { energyHistogram_lb_ = elb; } //!< Assign energy lower bound from restart
	void checkEnergyHistogramBounds ();
   	void refineEnergyHistogramBounds ();
	void recordEnergyHistogram ();
	void reInitializeEnergyHistogram ();
	void printEnergyHistogram (const std::string fileName, const bool normalize=true);
	void restartEnergyHistogram (const std::string prefix);

	void refinePkHistogramBounds ();
	void recordPkHistogram ();
	void printPkHistogram (const std::string fileName, const bool normalize=true);
	void restartPkHistogram (const std::string prefix);

	void recordExtMoments ();
	void printExtMoments (const std::string fileName, const bool normalize=true);
	void restartExtMoments (const std::string prefix, const std::vector < double > &ctr);

	void startWALA (const double lnF, const double g, const double s, const int Mtot); //!< Start using Wang-Landau and instantiate the bias object
   	void stopWALA () { useWALA = false; delete wlBias; } //!< Stop using Wang-Landau and free the bias object
   	void startTMMC (const long long int tmmcSweepSize, const int Mtot); //!< Start using TMMC and instantiate the bias object
   	void stopTMMC () { useTMMC = false; delete tmmcBias; } //!< Stop using TMMC and free the bias object

	void setTotNBounds (const std::vector < int > &bounds);

	void incrementMState ();
   	void decrementMState ();

	bool potentialIsSet (const int spec1, const int spec2) { return ppotSet_[spec1][spec2]; }	//!< Boolean which returns whether or not a pair has had its potential specified by the user yet

	const int nSpecies () { return nSpecies_; } //!< Return the number of different species in the system
	const int maxSpecies (const int index);
	const int minSpecies (const int index);
	const int totNMax () { return totNBounds_[1]; } //!< Return upper bound on the total number of atoms in the system
	const int totNMin () { return totNBounds_[0]; } //!< Return lower bound on the total number of atoms in the system
	const int getTotN () { return totN_; } //!< Return a sum of the total number of atoms currently in the system
	const int getMaxOrder () { return max_order_; } //!< Return the max order the extensive moments are being stored out to

	const int getCurrentM () { return Mcurrent_; } //!< Return the system's current expanded ensemble fractional state
	const int getTotalM () { return Mtot_; } //!< Return the total number of fractional states available to species in the expanded ensemble
	const int getFractionalAtomType () { return fractionalAtomType_; } //!< Return the atom type of the fractional atom

	const double energy () { return energy_; } //!< Return the system's instantaneous energy
	const double scratchEnergy ();

	const double beta () { return beta_; } //!< Return 1/kT
	const double mu (const int index) { return mu_[index]; } //!< Return the chemical potential for a given species' index
    const double mass (const int index) { return mass_[index]; } //!< Return the mass of a given species' index

	bool useTMMC; //!< Logical stating whether or not to use TMMC biasing
	bool useWALA; //!< Logical stating whether or not to use Wang-Landau biasing
	bool restartFromWALA, restartFromTMMC; //!< Flags to restart from WALA or TMMC initially

	long long int tmmcSweepSize; //!< Size of a sweep in TMMC
	long long int totalTMMCSweeps; //!< Total number of sweeps to perform during TMMC
	long long int wlSweepSize; //!< Size of Wang-Landau sweep
	long long int nCrossoverVisits; //!< Number of crossovers that must occur before switching from WALA to TMMC

	double lnF_start, lnF_end; //!< Starting and ending lnF for Wang-Landau
	double wala_g, wala_s; //!< Wang-Landau g and s factors

	const std::vector < double > extMomCounter () { return extensive_moments_.getCounter(); } //!< Get counter for extensive moments needed for restarting system from a checkpoint
	const std::vector < double > box () { return box_; } //!< Return the system box dimensions
	std::vector < double > getELB () { return energyHistogram_lb_; } //!< Returns current tally of energy min at each Ntot for checkpointing
	std::vector < double > getEUB () { return energyHistogram_ub_; } //!< Returns current tally of energy max at each Ntot for checkpointing
	std::vector < atom* > getNeighborAtoms (const unsigned int typeIndexA, const unsigned int typeIndexB, atom* _atom);
	tmmc* getTMMCBias (); //!< Return pointer to the TMMC bias
	wala* getWALABias (); //!< Return pointer to the Wang-Landau bias
	atom* getFractionalAtom () { return fractionalAtom_; } //!< Returns a pointer the atom in the system that is currently only fractionally inserted/deleted

	tmmc* tmmcBias; //!< TMMC biasing function
	wala* wlBias; //!< WL biasing function

	std::string restartFromWALAFile, restartFromTMMCFile; //!< Files to restart from WALA or TMMC initially

	std::vector < int > numSpecies;	//!< Total number of each type of atom the system contains
	std::vector < std::vector < atom > > atoms;	//!< Atoms in a matrix by type, and particle index, respectively that a system CAN hold but not all are actually "in" the system
	std::vector < std::vector < std::shared_ptr < pairPotential > > > ppot; //!< Matrix of pair potentials for atom types i, j
	std::vector < compositeBarrier > speciesBarriers; //!< Barriers, if any, for each species

private:
	atom* fractionalAtom_; //!< Pointer to the atom in the system that is currently only fractionally inserted/deleted

	bool toggleKE_; //!< Flag for the use of kinetic energy corrections

	int fractionalAtomType_; //!< Type of atom that is currently fractionally inserted
	int nSpecies_; //!< Number of species types allowed in the simulation (single component = 1, multicomponent > 1)
   	int Mcurrent_; //!< Fractional level of insertion of the current atom in an "expanded" state, all species have the same Mtot_
	int Mtot_; //!< Number of fractional states available to each atom of each species in the expanded ensemble, all species are identical
	int totN_; //!< Sum total of all atoms in the system
	int max_order_; //!< Maximum order for correlations

	double beta_; //!< Inverse temperature, really 1/kT
	double energy_; //!< Instantaneous energy of the system
	double energyHistDelta_; //!< Bin width for energy histograms at each Ntot

	std::vector < int > totNBounds_; //!< For multicomponent mixtures, biases use Shen and Errington method which uses bounds on the total number of particles in the system
	std::vector < int > maxSpecies_; //!< Maximum number of each species allowed in the simulation
	std::vector < int > minSpecies_; //!< Minimum number of each species allowed in the simulation

	std::vector < double > box_; //!< System box size
	std::vector < double > mass_; //!< Mass of each species.  Only used for deBroglie calculation if KE correction used
	std::vector < double > mu_; //!< Chemical potential of each species
	std::vector < double > energyHistogram_lb_; //!< Lowest energy observed for each Ntot, is adjusted to become lb for energyHistogram_
	std::vector < double > energyHistogram_ub_; //!< Highest energy observed for each Ntot, is adjusted in energyHistogram_

	std::vector < std::vector < bool > > ppotSet_; //!< Matrix of pair potentials between type i and j
	std::vector < std::vector < bool > > useCellList_;  //!< Matrix of whether or not to use cell lists to calculate potentials for pair type (i,j)

	std::vector < cellList > cellLists_; // this vector stores the actual cell lists for the inserted potentials
	std::vector < std::vector < cellList* > > cellListsByPairType_; // this matrix stores pointers to the actual cell lists for all pair types

	std::vector < dynamic_one_dim_histogram > energyHistogram_; //!< Histogram of energy at each Ntot
	std::vector < std::vector < dynamic_one_dim_histogram > > pkHistogram_; //!< Histogram of particle numbers at each Ntot [species_idx][Ntot]
	histogram extensive_moments_; //<! N_i^jN_k^mU^p[Ntot] matrix
};

const double calculateBias (simSystem &sys, const int nTotFinal, const int mFinal);

#endif
