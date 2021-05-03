#include "stdafx.h"
#include "MC_European.h"
#include <iostream>
#include <algorithm>



MC_European::MC_European(Matrix<double> S0, Matrix<double> weight, Matrix<double> sigma, double rate, int n_sim, std::vector<double> exercise_date, bool antithetic_method, bool pseudo_variate_control, bool pseudo_random_method)
	:MonteCarlo(S0, weight, sigma, rate, n_sim, pseudo_random_method), m_exercise_date(exercise_date), m_antithetic_method(antithetic_method), m_pseudo_variate_control(pseudo_variate_control)
{
	m_exercise_date.insert(m_exercise_date.begin(), 0.);
	m_T = int(m_exercise_date.size());
	m_S = Matrix<double>(m_assets, m_T);
	m_S.fill_column(0, m_S0);
	m_S_weighted = Matrix<double>(m_n_sim, m_T);
	m_S_exp_weighted = Matrix<double>(m_n_sim, m_T);
	m_Price = Matrix<double>(m_n_sim, 1);
	m_Price_prim = Matrix<double>(m_n_sim, 1);
	m_cholesky = m_sigma.Cholesky();

}


MC_European::~MC_European()
{
}

void MC_European::diffusion(bool antithetic)
{
	Matrix<double> bm(m_assets, 1);
	double time = 0;
	
	for (int i = 1; i < m_T; i++)
	{
		double delta_t = m_exercise_date[i] - m_exercise_date[i-1];
		// antithetic variance reduction
		bm = m_NormBox->Generate(m_assets) * std::pow(-1., m_antithetic_method);
		// correlate the brownian motian
		bm = m_cholesky.dot(bm);
		m_S.fill_column(i, m_S.column(i-1) * (((m_variance * 0.5) * (-1.) + m_rate) * delta_t + bm * std::pow(delta_t, 0.5)).exp());
	}
}


void MC_European::multiple_diffusions()
{
	for (int i = 0; i < int(m_n_sim / 2); i++)
	{
		diffusion(m_antithetic_method);
		m_S_weighted.fill_row(i, get_weighted_average());
		m_S_exp_weighted.fill_row(i, get_exp_weighted_average());
	}

	if (m_antithetic_method)
	{
		m_antithetic_method = 1 - m_antithetic_method;
		delete m_NormBox;
		UniformGenerator* Unif;
		if (m_pseudo_random_method)
			Unif = new VanDerCorput();
		else
			Unif = new EcuyerCombined();
		m_NormBox = new NormalBoxMuller(0., 1., Unif);
	}


	for (int i = int(m_n_sim / 2); i < m_n_sim; i++)
	{
		diffusion(m_antithetic_method);
		m_S_weighted.fill_row(i, get_weighted_average());
		m_S_exp_weighted.fill_row(i, get_exp_weighted_average());
	}
}



double MC_European::payoff(const Matrix<double> &S, int i, int t, double K)
{
	double Spot = S.elem_at(i, t);
	double res = std::max(Spot - K, 0.);
	return res;
}

Matrix<double> MC_European::payoff(const Matrix<double> &S, int t, double K)
{
	Matrix<double> res(m_n_sim, 1);
	for (int i = 0; i < m_n_sim; i++)
	{
		res.set_elem_at(i, 0, payoff(S, i, t, K));
	}
	return res;
}



void MC_European::compute_price(double K)
{
	
	int L = 2; // degree of the polynom
	// case for k = m_T - 1, the price is the payoff

	m_Price = payoff(m_S_weighted, m_T - 1, K);
	m_Price_prim = payoff(m_S_exp_weighted, m_T - 1, K) * std::exp(-m_rate * m_exercise_date.back());

	std::cout << "moyenne payoff = " << m_Price.mean() << std::endl;
	// loop for all the k, backwardation 
	for (int k = m_T - 2; k >= 1; k--)
	{
		Matrix<double> alpha = LongstaffSchwartz(L, k, K, m_S_weighted);
		//Matrix<double> alpha_prim = LongstaffSchwartz(L, k, K, m_S_exp_weighted);

		// The price at step k is the max between, the payoff at k and the expectation at k+1 discounted
		Matrix<double> P(m_n_sim, L);
		for (int i = 0; i < m_n_sim; i++)
			for (int l = 0; l < L; l++)
				P.set_elem_at(i, l, std::pow(m_S_weighted.elem_at(i, k), l+1));

		/*
		Matrix<double> P_prim(m_n_sim, L);
		for (int i = 0; i < m_n_sim; i++)
			for (int l = 0; l < L; l++)
				P_prim.set_elem_at(i, l, std::pow(m_S_exp_weighted.elem_at(i, k), l + 1));
		*/

		m_Price = maximum(payoff(m_S_weighted, k, K), P.dot(alpha));
		//m_Price_prim = P_prim.dot(alpha_prim);
	}
	// when k is equal to 0, we are back to a european call
	// the price is then the expectation discounted. We already take the expectation in the get_mean_price function. 
	// So here, we only need to discount m_Price
	m_Price = m_Price * std::exp(-m_rate * (m_exercise_date[1] - m_exercise_date[0]));


}

double MC_European::get_mean_price(double K)
{
	double prix;
	if (m_pseudo_variate_control)
		prix = get_pseudo_variate_control_price(K);
	else
		prix = m_Price.mean();
	return prix;
}

double MC_European::get_pseudo_variate_control_price(double K)
{
	double BS_price = Call_Basket(m_S0, m_sigma, m_exercise_date.back(), m_rate, K, m_weight).price();
	std::cout << "BlackScholes BasketCall = " << BS_price << std::endl;
	std::cout << "Berdumean MonteCarlo Price = " << m_Price.mean() << std::endl;
	std::cout << "BasketCall MonteCarlo Price = " <<  m_Price_prim.mean() << std::endl;
	double prix = (m_Price - m_Price_prim).mean() + BS_price;
	return prix;
}

double MC_European::get_variance()
{
	double var;
	if (m_pseudo_variate_control)
		var = get_pseudo_variate_control_variance();
	else
		var = m_Price.variance();
	return var;
}

double MC_European::get_pseudo_variate_control_variance()
{
	double var = (m_Price - m_Price_prim).variance();
	return var;
}


Matrix<double> MC_European::get_weighted_average()
{
	Matrix<double> S_weighted(m_exercise_date.size(), 1);
	for (int t = 0; t < int(m_exercise_date.size()); t++)
	{
		double weighted_average = m_S.column(t).transpose().dot(m_weight).elem_at(0, 0);
		S_weighted.set_elem_at(t, 0, weighted_average);
	}
	return S_weighted;
}

Matrix<double> MC_European::get_exp_weighted_average()
{
	Matrix<double> S_exp_weighted(m_exercise_date.size(), 1);
	for (int t = 0; t < int(m_exercise_date.size()); t++)
	{
		double exp_weighted_average = 0; 
		for (int i = 0; i < m_assets; i++)
		{
			exp_weighted_average += m_weight.elem_at(i) * std::log(m_S.elem_at(i, t));
		}
		exp_weighted_average = std::exp(exp_weighted_average);
		S_exp_weighted.set_elem_at(t, 0, exp_weighted_average);
	}
	return S_exp_weighted;
}

Matrix<double> MC_European::LongstaffSchwartz(int L, int k, double K, Matrix<double> &S)
{
	// P is the matrix of underlying power t
	Matrix<double> P(m_n_sim, L);
	// E is the matrix of discounted payoff / price
	Matrix<double> E(m_n_sim, 1);

	for (int i = 0; i < m_n_sim; i++)
		for (int j = 1; j <= L; j++)
			P.set_elem_at(i, j-1, std::pow(S.elem_at(i, k), j));
	
	for (int j = 0; j < m_n_sim; j++)
	{
		double delta_t = m_exercise_date[k + 1] - m_exercise_date[k];
		E.set_elem_at(j, 0, std::exp(-m_rate * delta_t) * m_Price.elem_at(j, 0));
	}

	// in order to find alpha, we compute the right term
	Matrix<double> Optimal_P(L, 1);
	for (int j = 0; j < m_n_sim; j++)
	{
		Optimal_P +=  P.transpose().column(j) * E.elem_at(j, 0);
	}
	
	// H computation (the left term)

	Matrix<double> H(L, L);
	for (int i0 = 0; i0 < L; i0++)
	{
		for (int i = 0; i < L; i++)
		{
			double coef = 0;
			for (int j = 0; j < m_n_sim; j++)
			{
				coef += P.elem_at(j, i) * P.elem_at(j, i0);
			}
			H.set_elem_at(i0, i, coef);
		}
	}

	//H.display();
	//H.inverse().display();
	//Optimal_P.display();

	Matrix<double> res = H.inverse().dot(Optimal_P);
	//res.display();
	return res;
		
}
