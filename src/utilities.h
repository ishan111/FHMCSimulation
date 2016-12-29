#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <vector>
#include <string>
#include <sys/stat.h>
#include <cstdio>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <time.h>
#include <iostream>
#include "global.h"

// JSON interface from local distro of rapidjson
#include "rapidjson/include/rapidjson/document.h"
#include "rapidjson/include/rapidjson/writer.h"
#include "rapidjson/include/rapidjson/stringbuffer.h"
#include "rapidjson/include/rapidjson/filereadstream.h"
#include "rapidjson/include/rapidjson/prettywriter.h"

void parseJson (const std::string filename, rapidjson::Document &doc);
void pauseCode (long int dur);
void pbc (std::vector < double > &pos, const std::vector < double > &box);
bool fileExists(std::string fileName);
std::string getTimeStamp ();
double rng (int *idum);
double pbcDist2 (const std::vector < double > &p1, const std::vector < double > &p2, const std::vector < double > &box);
std::vector < double > random3DSurfaceVector (const double magnitude);
std::vector < std::vector < double > > rotationMatrix (const double alpha, const double beta, const double gamma);
std::vector < std::vector < double > > matrixProduct (std::vector < std::vector < double > > &mat1, std::vector < std::vector < double > > &mat2);
std::vector < std::string > splitstr (const std::string &s, char delim);

#endif
