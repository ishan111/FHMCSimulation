#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <string>
#include <sys/stat.h>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <time.h>

//! For converting to string
#define sstr( x ) dynamic_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()

void pbc (std::vector < double > &pos, const std::vector < double > &box);
bool fileExists(std::string);
std::string getTimeStamp ();
double rng (int *idum);
double pbc_dist2 (const std::vector < double > &p1, const std::vector < double > &p2, const std::vector < double > &box);
std::vector < double > random3DSurfaceVector (const double magnitude);
std::vector < std::vector < double > > rotationMatrix (const double alpha, const double beta, const double gamma);
std::vector < std::vector < double > > matrixProduct (std::vector < std::vector < double > > &mat1, std::vector < std::vector < double > > &mat2);

#endif
