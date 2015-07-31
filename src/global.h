#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <string>
#include <exception>

#define SYS_FAILURE -1	//!< Returned from main the simualtion encountered a major error
#define SAFE_EXIT 0		//!< Returned from main the simulation ran without error
#define PI 3.14159265359

extern int RNG_SEED;	// defined elsewhere, but this allows other routines to access it

/*!
 * Error handling via custom exceptions.
 */
class customException : public std::exception {
public:
	const char* what() const throw() {return msg_.c_str();}	//!< Return the user's message
	customException (std::string m="custom exception occurred"):msg_(m){} //!< Instantiate an exeption with a user-defined error message
	~customException () throw() {};	//!< Throw the message to the user
	
private:
	std::string msg_;
};

#endif
