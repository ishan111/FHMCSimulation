#include "utilities.h"

int RNG_SEED = -1024;	//!< Default RNG seed

// for rng
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/*!
 * Random number generator (from Numerical Recipes)
 * 
 * \param [in] idum seed
 * 
 * \return temp Pseudo-random number between [0, 1)
 */
double rng (int *idum) {
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		} iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

/*!
 * Replace a position inside a box assuming periodic boundary conditions.
 * 
 * \param [in, out] pos Position to be placed in box
 * \param [in] box Box dimensions
 */
void pbc (std::vector < double > &pos, const std::vector < double > &box) {
	// generally while loops are faster than round statements
	for (unsigned int i = 0; i < pos.size(); ++i) {
		while (pos[i] < 0.0) {
			pos[i] += box[i];
		}
		while (pos[i] >= box[i]) {
			pos[i] -= box[i];
		}
	}
}

/*!
 * Calculate the minimum image distance squared between two coordinates assuming periodic boundary conditions.  Coordinates do not have to be in the box to begin with.
 *
 * \param [in] \p1 Position 1
 * \param [in] \p1 Position 2
 * \param [in] \box Box size
 *
 * \return d2 (distance squared)
 */
double pbc_dist2 (const std::vector < double > &p1, const std::vector < double > &p2, const std::vector < double > &box) {
    double d2 = 0.0;
    for (unsigned int i = 0; i < p2.size(); ++i) {
    	double dr = p2[i] - p1[i];
		while (dr < -box[i]/2.0) {
			dr += box[i];
		}
		while (dr > box[i]/2.0) {
			dr -= box[i];
		}
		d2 += dr*dr;
	}
    
    return d2;
}

/*!
 * Function to check whether a given file exists or not.
 * 
 * \param [in] fileName Name of file to check
 * 
 * \return If file exists
 */
bool fileExists(std::string fileName)
{
	struct stat stFileInfo;
	
	if (stat(fileName.c_str(),&stFileInfo) == 0)
		return true;
	else 
		return false;
}
