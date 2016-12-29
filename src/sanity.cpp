#include "sanity.h"

/*!
 * Perform sanity checks after simulation
 *
 * \param [in] sys System that was simulated
 */
void sanityChecks (simSystem &sys) {
	if (sys.nSpecies() != sys.atoms.size()) {
		sendErr("Error, number of components changed throughout simulation");
        exit(SYS_FAILURE);
	} else {
		sendMsg("Passed number of species present check");
	}

	long long int ns = 0;
	for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
		ns += sys.numSpecies[i];
	}
	if (ns != sys.getTotN()) {
		sendErr("Sum of fully inserted atoms deviates from total counter");
		exit(SYS_FAILURE);
	} else {
		sendMsg("Passed sum of atoms consistency check with total counter");
	}

	if (sys.getTotalM() > 1) {
		if (sys.getFractionalAtom()->mState != sys.getCurrentM()) {
			sendErr("Expanded ensemble state deviates between atom ("+numToStr(sys.getFractionalAtom()->mState)+") and system ("+numToStr(sys.getCurrentM())+")");
			exit(SYS_FAILURE);
		}
		for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
			int end = sys.numSpecies[i];
			if (i == sys.getFractionalAtomType()) {
				end++;
			}
			for (unsigned int j = 0; j < end; ++j) {
				if (&sys.atoms[i][j] != sys.getFractionalAtom()) {
					if (sys.atoms[i][j].mState != 0) {
						sendErr("Atom ("+numToStr(i)+", "+numToStr(j)+") has non-zero expanded ensemble state ("+numToStr(sys.atoms[i][j].mState)+")");
						exit(SYS_FAILURE);
					}
				} else {
					if (sys.atoms[i][j].mState != sys.getCurrentM()) {
						sendErr("Fractional atom ("+numToStr(i)+", "+numToStr(j)+")'s expanded ensemble state ("+numToStr(sys.atoms[i][j].mState)+") does not match system's ("+numToStr(sys.getCurrentM())+")");
						exit(SYS_FAILURE);
					}
				}
			}
		}
	} else {
		for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
			for (unsigned int j = 0; j < sys.numSpecies[i]; ++j) {
				if (sys.atoms[i][j].mState != 0) {
					sendErr("Atom ("+numToStr(i)+", "+numToStr(j)+") has non-zero expanded ensemble state ("+numToStr(sys.atoms[i][j].mState)+")");
					exit(SYS_FAILURE);
				}
			}
		}
	}
	sendMsg("Passed expanded ensemble state check for all atoms");

	const double tol = 1.0e-6;
	const double scratchEnergy = sys.scratchEnergy(), incrEnergy = sys.energy();
    if (fabs(scratchEnergy - incrEnergy) > tol) {
		sendErr("Error, scratch energy calculation = "+numToStr(scratchEnergy)+", but incremental = "+numToStr(incrEnergy)+", |diff| = "+numToStr(fabs(scratchEnergy - incrEnergy)));
        exit(SYS_FAILURE);
    } else {
		sendMsg("Passed, final scratch energy - incremental = "+numToStr(scratchEnergy - incrEnergy));
    }
}
