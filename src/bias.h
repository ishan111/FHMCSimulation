#ifndef BIAS_H_
#define  BIAS_H_

#include <vector>
#include <string>

class tmmc {
public:
	tmmc (const int nSpec, const std::vector <int> &Nmax, const std::vector <int> &Nmin);
	~tmmc () {};
	
	const int getAddress (const std::vector <int> &Nstart, const std::vector <int> &Nend);
	void updateC (const std::vector <int> &Nstart, const std::vector <int> &Nend, const double pa);
	void calculateP ();
	
private:
	void difference_ (const std::vector <int> &Nstart, const std::vector <int> &Nend, int &specId, int &addOrSubtract);
	const int tmmc::getAddress (const std::vector <int> &Nstart, const std::vector <int> &Nend);
	int nSpec_; // number of species in simulation
	std::vector <int> Nmax_, Nmin_, W_, C_, P_;
};

#endif
