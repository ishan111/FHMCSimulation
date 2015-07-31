#ifndef SYSTEM_H_
#define SYSTEM_H_

#include "potentials.h"
#include "atom.h"
#include "cellList.h"
#include <string>
#include <vector>

/*!
 * System information for the simulation.
 */
class simSystem {
public:
	simSystem (const unsigned int nSpecies, const double beta, const std::vector < double > box, const std::vector < double > mu, const std::vector < int > maxSpecies);
	~simSystem () {};
    
    const int nSpecies () { return nSpecies_; }
    const int maxSpecies (const int index);
    const double energy () { return energy_; }
    void incrementEnergy (const double dU) { energy_ += dU; }
    const double scratchEnergy ();
    const double beta () { return beta_; }
    const double mu (const int index) { return mu_[index]; }
	const std::vector < double > box () { return box_; }
	void addPotential (const int spec1, const int spec2, pairPotential *pp, bool useCellList=false);
    void printSnapshot (std::string filename, std::string comment);
    void insertAtom (const int typeIndex, atom *newAtom);
    void deleteAtom (const int typeIndex, const int atomIndex);
    void translateAtom (const int typeIndex, const int atomIndex, std::vector<double> oldPos);
    void readRestart (std::string filename);
	
	std::vector< std::vector<double> > getNeighborPositions(const unsigned int typeIndexA, const unsigned int typeIndexB, atom* _atom);
		
    // these are updated very often so leave open to user
    std::vector < std::vector < atom > > atoms;	//!< Atoms in a matrix by type, and particle index, respectively that a system CAN hold but not all are actually "in" the system
    std::vector < int > numSpecies;		//!< Total number of each type of atom the system contains
    std::vector < std::vector < pairPotential* > > ppot;	//!< Matrix of pair potentials for atom types i, j
	bool potentialIsSet (const int spec1, const int spec2) { return ppotSet_[spec1][spec2]; }	//!< Boolean which returns whether or not a pair has had its potential specified by the user yet
    
private:
    int nSpecies_;
    std::vector < int > maxSpecies_;
    double beta_, energy_;
	std::vector < double > box_;
	std::vector < double > mu_;
    std::vector < std::vector < bool > > ppotSet_;
    std::vector < std::vector < bool > > useCellList_;
    std::vector <cellList> cellLists_; // this vector stores the actual cell lists for the inserted potentials
    std::vector < std::vector <cellList*> > cellListsByPairType_; // this matrix stores pointers to the actual cell lists for all pair types
};

#endif
