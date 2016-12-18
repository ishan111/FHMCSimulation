#include "sanity.h"

/*!
 * Perform sanity checks after simulation
 *
 * \param [in] sys System that was simulated
 */
void sanityChecks (simSystem &sys) {
    // Sanity checks
	if (sys.nSpecies() != sys.atoms.size()) {
        std::cerr << "Error: Number of components changed throughout simulation" << std::endl;
        exit(SYS_FAILURE);
	}
	if (sys.getTotalM() > 1) {
		if (sys.getFractionalAtom()->mState != sys.getCurrentM()) {
			std::cerr << "Expanded ensemble state deviates between atom ("+std::to_string(sys.getFractionalAtom()->mState)+") and system log ("+std::to_string(sys.getCurrentM())+")" << std::endl;
			exit(SYS_FAILURE);
			for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
				int end = sys.numSpecies[i];
				if (i == sys.getFractionalAtomType()) {
					end++;
				}
				for (unsigned int j = 0; j < end; ++j) {
					if (&sys.atoms[i][j] != sys.getFractionalAtom()) {
						if (sys.atoms[i][j].mState != 0) {
							std::cerr << "Atom ("+std::to_string(i)+", "+std::to_string(j)+") has non-zero expanded ensemble state ("+std::to_string(sys.atoms[i][j].mState)+")" << std::endl;
							exit(SYS_FAILURE);
						}
					} else {
						if (sys.atoms[i][j].mState != sys.getCurrentM()) {
							std::cerr << "Fractional atom ("+std::to_string(i)+", "+std::to_string(j)+")'s expanded ensemble state ("+std::to_string(sys.atoms[i][j].mState)+") does not match system's ("+std::to_string(sys.getCurrentM())+")" << std::endl;
							exit(SYS_FAILURE);
						}
					}
				}
			}
		}
	} else {
		for (unsigned int i = 0; i < sys.nSpecies(); ++i) {
			for (unsigned int j = 0; j < sys.numSpecies[i]; ++j) {
				if (sys.atoms[i][j].mState != 0) {
					std::cerr << "Atom ("+std::to_string(i)+", "+std::to_string(j)+") has non-zero expanded ensemble state ("+std::to_string(sys.atoms[i][j].mState)+")" << std::endl;
					exit(SYS_FAILURE);
				}
			}
		}
	}

    // Still allow for printing of all data, even if there is an error, in order to interrogate the results anyway
	const double tol = 1.0e-6;
	const double scratchEnergy = sys.scratchEnergy(), incrEnergy = sys.energy();
    if (fabs(scratchEnergy - incrEnergy) > tol) {
        std::cerr << "Error: scratch energy calculation = " << std::setprecision(20) << scratchEnergy << ", but incremental = " << std::setprecision(20) << incrEnergy << ", |diff| = " << std::setprecision(20) << fabs(scratchEnergy - incrEnergy) << std::endl;
        exit(SYS_FAILURE);
    } else {
        std::cout << "Passed: Final scratch energy - incremental = " << std::setprecision(20) << scratchEnergy - incrEnergy << std::endl;
    }
}
