#include "stdafx.h"
#include "VanDerCorput.h"


VanDerCorput::VanDerCorput() : PseudoGenerator()
{
	m_n = 100;
	m_d = 10;
	m_p = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47 };
	m_i = 0;
}


VanDerCorput::~VanDerCorput()
{
}

double VanDerCorput::Generate()
{
	int p = m_p[m_i];
	double phi = 0;
	int ak;
	int puissance = 0;
	m_n_rec = m_n;

	while (m_n_rec != 0)
	{
		ak = m_n_rec % p;
		m_n_rec = (m_n_rec - ak) / p;
		phi += ak / std::pow(p, ++puissance);
	}

	if (++m_i >= m_d)
	{
		m_i = 0;
		m_n++;
	}

	return phi;
}