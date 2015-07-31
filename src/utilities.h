#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <string>
#include <sys/stat.h> 
#include <cstdio>

double pbc_dist2 (const std::vector < double > &p1, const std::vector < double > &p2, const std::vector < double > &box);
void pbc (std::vector < double > &pos, const std::vector < double > &box);
double rng (int *idum);
bool fileExists(std::string);

#endif
