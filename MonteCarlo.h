#pragma once
#include "../RandomGenerator/Normal.h"
#include "../RandomGenerator/EcuyerCombined.h"
#include "../RandomGenerator/VanDerCorput.h"
#include "../Matrix/Matrix.h"

class MonteCarlo
{
public:
	MonteCarlo(Matrix<double> S0, Matrix<double> weight, Matrix<double> sigma, double rate, int n_sim, bool pseudo_random_method);
	~MonteCarlo();
	//void diffusion() virtual;

protected:
	Matrix<double> m_S0; 
	Matrix<double> m_weight;
	Matrix<double> m_sigma;
	Matrix<double> m_variance;
	double m_rate;
	int m_assets;
	int m_n_sim;
	NormalBoxMuller* m_NormBox;
	bool m_pseudo_random_method;
};

