#include "aggvolbias.h"

/*!
 * Instantiate an aggregation-volume bias move to swap particles.
 *
 * \param [in] typeIndex Type to consider in swapping

 */
aggVolBias::aggVolBias (const int typeIndex, const int typeIndex2, const double pBias; const std::string tag) { 
	typeIndex_ = typeIndex; 
	typeIndex2_ = typeIndex2; 
	name_ = tag + sstr(typeIndex); 
	if (pBias < 1 && pBias > 0) {
		pBias_ = pBias;
	} else {
		throw customException ("Bias probability must be > 0 and < 1 for aggVolBias");
	}
	changeN_ = false; 
}

/*!
 * Use aggregation-volume bias to swap two particles in the system.  All other information is stored in the simSystem object.
 * 
 * \param [in] sys System object to attempt to swap particles in.
 * 
 * \return MOVE_SUCCESS if particles are swapped, otherwise MOVE_FAILURE if not.  Will throw exceptions if there was an error.
 */
int aggVolBias::make (simSystem &sys) {
	// choose a particle of typeIndex_
	int nAvail1 = sys.numSpecies[typeIndex_];
	if (sys.getCurrentM() > 0 && sys.getFractionalAtomType() == typeIndex_) {
		nAvail1++;
	}

	// choose a particle of typeIndex2_
	int nAvail2 = sys.numSpecies[typeIndex2_];
	if (sys.getCurrentM() > 0 && sys.getFractionalAtomType() == typeIndex2_) {
                nAvail2++;
        }

	// update biases even upon failure to try move (also reject if types are same and only one overall particle)
	if (nAvail1 < 1 || nAvail2 < 1 || (typeIndex_ == typeIndex2_ && nAvail1 == 1)) {
		if (sys.useWALA) {
                        sys.getWALABias()->update(sys.getTotN(), sys.getCurrentM());
                }
                if (sys.useTMMC) {
                        sys.tmmcBias->updateC (sys.getTotN(), sys.getTotN(), sys.getCurrentM(), sys.getCurrentM(), 0.0);
                }
		return MOVE_FAILURE;
	}

	// choose particles, j and k - at present we are guaranteed at least 2 unique particles in the system
	int pkJ = (int) floor(rng (&RNG_SEED) * nAvail1), pkK = 0;
	bool jkOverlap = true;
	while (jkOverlap) {
		pkK = (int) floor(rng (&RNG_SEED) * nAvail2);
		
		// establish if j and k overlap their "in" regions
		const double rc = sys->ppot[typeIndex_][typeIndex2_]->rcut();
		
		// if pkK and pkJ are withing rcut of each other, they overlap
		if (pbc_dist2(sys.atoms[typeIndex_][pkJ]->pos, sys.atoms[typeIndex2_][pkK]->pos, sys.box()) > rc*rc) {
			jkOverlap = false;
		}

		// However, we also technically need to ensure that if typeIndex_ == typeIndex2_, j and k are not the same particle.  
		// Happily, this is already guaranteed by the above distance check (a particle is a distance 0 from itself)
	}

	if (rng(&RNG_SEED) < pBias_) {
		// choose a particle "in" j
		
		
	} else {
		// choose a particle "in" k or "out" j with equal probability
		// note that "out" j includes the "in" k region as well
	
		if (rng (&RNG_SEED) < 0.5) {
			// choose particle "in" k

		} else {
			// choose particle "out" of j
		
		}
	}
	
	// remember to upated TMMC and WL regardless of what happens


	return MOVE_FAILURE;
}
