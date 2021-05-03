#pragma once
#include "MonteCarlo.h"
#include <vector>
#include "../Matrix/Matrix.h"
#include "../RandomGenerator/Normal.h"
#include "../RandomGenerator/EcuyerCombined.h"
#include "../RandomGenerator/VanDerCorput.h"
#include "../BlackScholes/CallBasket.h"

class MC_European : public MonteCarlo
{
public:
	MC_European(Matrix<double> S0, Matrix<double> weight, Matrix<double> sigma, double rate, int n_sim, std::vector<double> exercise_date, bool antithetic_method, bool pseudo_variate_control, bool pseudo_random_method);
	~MC_European();
	
	void multiple_diffusions();
	void compute_price(double K);
	double get_mean_price(double K);
	double get_variance();


protected:
	//std::vector<double> m_S;
	std::vector<double> m_exercise_date;
	Matrix<double> m_S;
	Matrix<double> m_S_weighted;
	Matrix<double> m_Price;
	Matrix<double> m_Price_prim;
	Matrix<double> m_S_exp_weighted;
	Matrix<double> m_cholesky;
	int m_T;
	bool m_antithetic_method;
	bool m_pseudo_variate_control;
	
	

private:
	Matrix<double> get_weighted_average();
	Matrix<double> get_exp_weighted_average();
	void diffusion(bool antithetic);

	double payoff(const Matrix<double> &S, int i, int t, double K);
	Matrix<double> payoff(const Matrix<double> &S, int t, double K);
	Matrix<double> LongstaffSchwartz(int L, int k, double K, Matrix<double> &S);
	double get_pseudo_variate_control_price(double K);
	double get_pseudo_variate_control_variance();
	
};

