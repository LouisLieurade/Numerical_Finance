#ifndef VANDERCORPUT_H
#define VANDERCORPUT_H

#include <vector>
#include "../RandomGenerator/PseudoGenerator.h"


class VanDerCorput : public PseudoGenerator
{
public:
	VanDerCorput();
	~VanDerCorput();
	virtual double Generate();
protected:
	int m_d;
	std::vector<int> m_p;
	int m_n;
	int m_i; // index of p
	int m_n_rec;
};


#endif 