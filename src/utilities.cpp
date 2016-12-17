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
 * C++ explode a string.
 */
std::vector < std::string > splitstr (const std::string &s, char delim) {
    std::stringstream ss(s);
    std::string item;
    std::vector <std::string> tokens;
    while (std::getline(ss, item, delim)) {
        tokens.push_back(item);
    }
    return tokens;
}

/*!
 * Obtain instantaneous timestamp.
 *
 * \return char* timestamp
 */
std::string getTimeStamp () {
    time_t rawtime;
    time (&rawtime);
    struct tm * timeinfo;
    timeinfo = localtime (&rawtime);
    char timestamp [80];
    strftime (timestamp,80,"%d/%m/%Y %H:%M:%S",timeinfo);
    std::string ans (timestamp);
    return ans;
}

/*!
 * Compute product of 2 matrices, W = UV
 *
 * \param [in] mat1 Matrix U
 * \param [in] mat2 Matrix V
 * \return W
 */
std::vector < std::vector < double > > matrixProduct (std::vector < std::vector < double > > &mat1, std::vector < std::vector < double > > &mat2) {
    std::vector < double > dummy (3, 0);
    std::vector < std::vector < double > > W (3, dummy);

    for (unsigned int i = 0; i < 3; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            W[i][j] = mat1[i][j]*mat2[j][i];
        }
    }

    return W;
}

/*!
 * Create 3D rotation matrix from angles (x, then y, then z)
 *
 * \param [in] alpha Radians to rotate centers by around x-axis
 * \param [in] beta Radians to rotate centers by around y-axis
 * \param [in] gamma Radians to rotate centers by around z-axis
 */
std::vector < std::vector < double > > rotationMatrix (const double alpha, const double beta, const double gamma) {

    std::vector < double > dummy (3, 0);
	std::vector < std::vector < double > > Rx (3, dummy), Ry(3, dummy), Rz(3, dummy), Ryx, Rzyx;

	// https://en.wikipedia.org/wiki/Rotation_matrix#General_rotations
	Rx[0][0] = 1.0;
	Rx[1][1] = cos(alpha);
	Rx[1][2] = -sin(alpha);
	Rx[2][1] = sin(alpha);
	Rx[2][2] = cos(alpha);

	Ry[0][0] = cos(beta);
	Ry[0][2] = sin(beta);
	Ry[1][1] = 1.0;
	Ry[2][0] = -sin(beta);
	Ry[2][2] = cos(beta);

	Rz[0][0] = cos(gamma);
	Rz[0][1] = -sin(gamma);
	Rz[1][0] = sin(gamma);
	Rz[1][1] = cos(gamma);
	Rz[2][2] = 1.0;

	Ryx = matrixProduct(Ry, Rx);
	Rzyx = matrixProduct(Rz, Ryx);

    return Rzyx;
}

/*!
 * Chose a random vector sampled from a random distribution on the surface of a sphere.
 *
 * \param [in] magnitude Mangitude of the vector (radius of sphere)
 *
 * \return ans double3 containing coordinates of vector
 */
std::vector < double > random3DSurfaceVector (const double magnitude) {
    int success = 0;
    std::vector < double > ans (3, 0);
    while (success == 0) {
        double r1 = rng(&RNG_SEED), r2 = rng(&RNG_SEED), x1 = 1.0-2.0*r1, x2 = 1.0-2.0*r2;
        double sum2 = x1*x1+x2*x2;
        if (sum2 < 1.0) {
            ans[0] = 2*x1*sqrt(1-sum2)*magnitude;
            ans[1] = 2*x2*sqrt(1-sum2)*magnitude;
            ans[2] = (1-2.0*sum2)*magnitude;
            success = 1;
        }
    }
    return ans;
}

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
double pbcDist2 (const std::vector < double > &p1, const std::vector < double > &p2, const std::vector < double > &box) {
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
bool fileExists(std::string fileName) {
	struct stat stFileInfo;

	if (stat(fileName.c_str(),&stFileInfo) == 0)
		return true;
	else
		return false;
}
