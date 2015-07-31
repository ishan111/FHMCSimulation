#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include "input.h"
#include "global.h"
#include <utility>
#include <sstream>
#include <stdio.h>      
#include <time.h>
#include <fstream>
#include "utilities.h"

using namespace boost;

/*!
 * Parse the input string from the command line.
 * \param [in] argc Command line
 * \param [in] argv Command line
 * \param [in] rawtime Time stamp from the time the binary was invoked
 * \return ans Vector of strings with parsed information, see ==help option from command line
 */
std::vector < std::string > parseInput (int argc, char * const argv[], time_t rawtime) {
	// report in time-stamped log file 
	struct tm * timeinfo;
	char buffer [80];
	timeinfo = localtime (&rawtime);
	strftime (buffer,80,"%Y_%m_%d_%H_%M_%S-input.log",timeinfo);
	std::ofstream logfile (buffer);
	
	// reflect input to use
	std::string input_args, restartFileName = "none";
	std::vector < std::string > parsed;
	logfile << "# Input: " << std::endl;
	for (unsigned int i = 0; i < argc; ++i) {
		std::string dummy(argv[i]);
		logfile << argv[i] << " ";
		input_args.append(dummy);
	}
	logfile << std::endl;
	logfile.close();
	
	// variables needed
	int N = 0;
	int seed = RNG_SEED;	// take default from utilities.cpp
	double L = -1, beta = -1;
	std::vector < double > mu;
    std::vector < int > maxAtom;
    std::vector < std::pair < int, double > > tmpMu;
    std::vector < std::pair < int, int > > tmpMax;
	long long int prodSweep = -1, equilSweep = -1, sweepSize = -1;
	
	split (parsed, input_args, is_any_of("=,:, "), token_compress_on);
	for (unsigned int i = 1; i < parsed.size(); ++i) {
        std::string dummy = parsed[i].substr(0, parsed[i].size()-1);
		int dummyIndex = atoi(parsed[i].substr(parsed[i].size()-1, parsed[i].size()).c_str()) - 1; // -1 to begin counting from 1 for user
	
		if (parsed[i] == "N") {
			N = (int) atoi(parsed[i+1].c_str());
			if (N < 1) {
				throw customException("Illegal number of components");
			}
			i++;
		} else if (parsed[i] == "L") {
			L = (double) atof(parsed[i+1].c_str());
			if (L < 0.0) {
				throw customException("Illegal box size");
			}
			i++;
		} else if (parsed[i] == "T") {
			beta = 1.0/((double) atof(parsed[i+1].c_str()));
			if (beta < 0.0) {
				throw customException("Illegal beta");
			}
			i++;
		} else if (dummy == "mu") {
			std::pair < int, double > newPair = std::make_pair (dummyIndex, (double) atof(parsed[i+1].c_str()));
			tmpMu.push_back(newPair);
			i++;
		} else if (dummy == "max") {
			std::pair < int, int > newPair = std::make_pair (dummyIndex, (int) atoi(parsed[i+1].c_str()));
			tmpMax.push_back(newPair);
			i++;
		} else if (parsed[i] == "sweepSize") {
			sweepSize = atoi(parsed[i+1].c_str());
			if (sweepSize < 1) {
				throw customException("Illegal sweepSize");
			}
			i++;
		} else if (parsed[i] == "prodSweep") {
			prodSweep = atoi(parsed[i+1].c_str());
			if (prodSweep < 1) {
				throw customException("Illegal prodSweep");
			}
			i++;
		} else if (parsed[i] == "equilSweep") {
			equilSweep = atoi(parsed[i+1].c_str());
			if (equilSweep < 1) {
				throw customException("Illegal equilSweep");
			}
			i++;
		} else if (parsed[i] == "seed") {
			seed = atoi(parsed[i+1].c_str());
			i++;
		} else if (parsed[i] == "restart") {
			restartFileName = parsed[i+1];
			i++;
		} else if (parsed[i] == "help") {
			std::cerr << "./omcs ==N: ==L: ==T: ==muX: ==maxX: ==sweepSize: ==prodSweep: ==equilSweep: ==seed: ==restart:" << std::endl;
			std::cerr << "\tParameters:" << std::endl;
			std::cerr << "\t\t==N:\t Number of components" << std::endl;
			std::cerr << "\t\t==L:\t Box dimensions" << std::endl;
			std::cerr << "\t\t==T:\t Temperature (kT)" << std::endl;
			std::cerr << "\t\t==muX:\t Chemical potential of component X (X >= 1)" << std::endl;
            std::cerr << "\t\t==maxX:\t Max number of atoms to allow of component X (X >= 1), defaults to 0 if not specified" << std::endl;
			std::cerr << "\t\t==sweepSize:\t Number of steps in a MC sweep" << std::endl;
			std::cerr << "\t\t==prodSweep:\t Number of sweeps in production phase" << std::endl;
			std::cerr << "\t\t==equilSweep:\t Number of steps in equilibration phase" << std::endl;
			std::cerr << "\t\t==seed:\t (optional) RNG seed, default is -1024" << std::endl;
			std::cerr << "\t\t==restart:\t (optional) XYZ file to read restart configuration from" << std::endl;
			throw customException ("Failed to start");
		} else if (parsed[i] != "./omcs") {
			std::string a = "Unrecognized input " + parsed[i] + ", try ./omcs --help";
			throw customException ("Failed to start: "+a);
		}
	}
	
	// check ranges
	if (N < 1 || L < 0 || beta < 0 || tmpMu.size() != N || tmpMax.size() != N) {
		throw customException ("Failed to specify valid parameters");
	}
	mu.resize(N);
	std::vector < int > muFound (N, 0);
	for (unsigned int i = 0; i < tmpMu.size(); ++i) {
		mu[tmpMu[i].first] = tmpMu[i].second;
		muFound[tmpMu[i].first] = 1;
	}
	for (unsigned int i = 0; i < muFound.size(); ++i) {
		if (muFound[tmpMu[i].first] != 1) {
			std::string number = static_cast<std::ostringstream*>( &(std::ostringstream() << tmpMu[i].first+1) )->str();
			throw customException ("Failed to specify mu for species number "+number);
		}
	}
    maxAtom.resize(N);
    std::vector < int > maxFound (N, 0);	// default to zero to protect the program from user error
	for (unsigned int i = 0; i < tmpMax.size(); ++i) {
		maxAtom[tmpMax[i].first] = tmpMax[i].second;
        if (tmpMax[i].second <= 0) {
            std::string number = static_cast<std::ostringstream*>( &(std::ostringstream() << tmpMax[i].first+1) )->str();
			throw customException ("Illegal max number of atoms for species number "+number);
        }
		maxFound[tmpMax[i].first] = 1;
	}
	for (unsigned int i = 0; i < maxFound.size(); ++i) {
		if (maxFound[tmpMax[i].first] != 1) {
			std::string number = static_cast<std::ostringstream*>( &(std::ostringstream() << tmpMax[i].first+1) )->str();
			throw customException ("Failed to specify max number of atoms for species number "+number);
		}
	}
    
	// return strings to user
	std::vector < std::string > ans (3+N+N+5);
	std::string strN = static_cast<std::ostringstream*>( &(std::ostringstream() << N) )->str();
	ans[0] = strN;
	std::string strL = static_cast<std::ostringstream*>( &(std::ostringstream() << L) )->str();
	ans[1] = strL;
	std::string strBeta = static_cast<std::ostringstream*>( &(std::ostringstream() << beta) )->str();
	ans[2] = strBeta;
	for (unsigned int i = 0; i < N; ++i) {
		std::string strMu = static_cast<std::ostringstream*>( &(std::ostringstream() << mu[i]) )->str();
		ans[3+i] = strMu;
	}
    for (unsigned int i = 0; i < N; ++i) {
		std::string strMax = static_cast<std::ostringstream*>( &(std::ostringstream() << maxAtom[i]) )->str();
		ans[3+N+i] = strMax;
	}
	std::string strSweepSize = static_cast<std::ostringstream*>( &(std::ostringstream() << sweepSize) )->str();
	ans[3+N+N] = strSweepSize;
	std::string strProdSweep = static_cast<std::ostringstream*>( &(std::ostringstream() << prodSweep) )->str();
	ans[3+N+N+1] = strProdSweep;
	std::string strEquilSweep = static_cast<std::ostringstream*>( &(std::ostringstream() << equilSweep) )->str();
	ans[3+N+N+2] = strEquilSweep;
	std::string strRNGseed = static_cast<std::ostringstream*>( &(std::ostringstream() << seed) )->str();
	ans[3+N+N+3] = strRNGseed;
	ans[3+N+N+4] = restartFileName;
	return ans;
}