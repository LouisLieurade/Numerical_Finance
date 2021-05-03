#include "stdafx.h"
#include "MonteCarlo.h"
#include <iostream>

MonteCarlo::MonteCarlo(Matrix<double> S0, Matrix<double> weight, Matrix<double> sigma, double rate, int n_sim, bool pseudo_random_method)
	:m_S0(S0), m_weight(weight), m_sigma(sigma), m_rate(rate), m_n_sim(n_sim), m_pseudo_random_method(pseudo_random_method)
{
	UniformGenerator* Unif;
	if (m_pseudo_random_method)
		Unif = new VanDerCorput();
	else
		Unif = new EcuyerCombined();

	m_NormBox = new NormalBoxMuller(0., 1., Unif);
	m_assets = S0.get_rows();
	m_variance = m_sigma.get_diagonal();
}


MonteCarlo::~MonteCarlo()
{
}
