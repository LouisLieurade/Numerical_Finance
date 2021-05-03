// NumericalFinance.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include <iostream>

#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/EcuyerCombined.h"
#include "../RandomGenerator/FiniteSet.h"
#include "../RandomGenerator/Exponential.h"
#include "../RandomGenerator/Normal.h"
#include "../RandomGenerator/Poisson.h"
#include "../SDE/BSEuler1D.h"
#include "../PDE/PDEGrid2D.h"
#include "../SDE/MonteCarlo.h"
#include "../SDE/MC_European.h"
#include "../Matrix/Matrix.h"
//#include "VanDerCorput.h"


void test_MonteCarlo(bool antithetic_method = false, bool pseudo_variate_control = false, bool pseudo_random_method = false)
{

	int n_assets = 3;
	int n_sim = 1000;
	double S00[3] = { 100., 100., 100. };
	double correl[9] = { 1., 0.4, 0.6, 0.4, 1., 0.8, 0.6, 0.8, 1. };
	//double correl[9] = { 1., 0.99, 0.99, 0.99, 1., 0.99, 0.99, 0.99, 1. };
	//double correl[9] = { 1., 0, 0, 0, 1., 0, 0, 0, 1. };
	double vol[3] = { 0.3, 0.4, 0.6 };
	Matrix<double> Correlation(n_assets, n_assets, correl);
	Matrix<double> Vol(n_assets, 1, vol);
	Matrix<double> sigma(n_assets);
	for (int i = 0; i < n_assets; i++)
		for (int j = 0; j < n_assets; j++)
			sigma.set_elem_at(i, j, Correlation.elem_at(i, j) * Vol.elem_at(i, 0) * Vol.elem_at(j, 0));


	double weight0[3] = { 0.3, 0.3, 0.4 };
	Matrix<double> S0(n_assets, 1, S00);
	Matrix<double> weight(n_assets, 1, weight0);
	double rate = 0.01;
	//std::vector<double> exercise_date = { 0.1, 0.2, 0.3, 0.4, 0.5 };
	std::vector<double> exercise_date = { 0.5 };
	double K = 100.;
	MC_European asset(S0, weight, sigma, rate, n_sim, exercise_date, antithetic_method, pseudo_variate_control, pseudo_random_method);


	asset.multiple_diffusions();
	asset.compute_price(K);


	double prix = asset.get_mean_price(K);
	double var = asset.get_variance();

	//double prix = asset.get_pseudo_variate_control_price(K);
	//double var = asset.get_pseudo_variate_control_variance();

	std::cout << "Bravo votre projet est termine, le prix est de : " << prix << std::endl;
	std::cout << "Bravo à vous deux, la variance est de : " << var << std::endl;
	std::cout << "L'intervalle de confiance pour le prix est : [" << prix << " - " << 1.96 * std::pow(var / n_sim, 0.5) << "; " << prix << " + " << 1.96 * std::pow(var / n_sim, 0.5)<< "] a 95%." << std::endl;
}




int main()
{
	bool antithetic_method = true;
	bool pseudo_variate_control = true;
	bool pseudo_random_method = true;
	test_MonteCarlo(antithetic_method, pseudo_variate_control, pseudo_random_method);

	
}

