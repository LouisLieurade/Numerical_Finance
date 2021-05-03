#include "stdafx.h"
#include "../BlackScholes/CallBasket.h"
#include <cmath>
#include <iostream>

double N(double x) // Phi(-∞, x) aka N(x)
{
	return std::erfc(-x / std::sqrt(2)) / 2;
}

Call_Basket::Call_Basket(Matrix<double> S, Matrix<double> sigma, double T, double rate, double K, Matrix<double> weight)
	:m_S(S), m_sigma(sigma), m_T(T), m_rate(rate), m_K(K), m_weight(weight)
{
	m_assets = m_S.get_rows();
}


Call_Basket::~Call_Basket()
{
}

double Call_Basket::get_spot()
{
	double res = 1.;
	for (int i = 0; i < m_assets; i++)
	{
		res *= std::pow(m_S.elem_at(i, 0), m_weight.elem_at(i, 0));
	}
	return res;
}

double Call_Basket::get_rate()
{
	double sum = 0.;
	double sigma_i;
	for (int i = 0; i < m_assets; i++)
	{
		sigma_i = 0.;
		for (int j = 0; j <= i; j++)
		{
			sigma_i += m_sigma.elem_at(i, j) * m_sigma.elem_at(i, j);
		}
		sum += m_weight.elem_at(i, 0) * sigma_i;
	}
	double res = (m_weight.transpose().dot(m_sigma).dot(m_sigma.transpose()).dot(m_weight)).elem_at(0, 0) * 0.5 + m_rate - 0.5 * sum;
	return res;
}

double Call_Basket::get_vol()
{
	Matrix<double> cholesky = m_sigma.Cholesky();
	double res = std::pow((m_weight.transpose().dot(cholesky).dot(cholesky.transpose()).dot(m_weight)).elem_at(0, 0), 0.5);
	return res;
}

double Call_Basket::price()
{
	double S = get_spot();
	double r = get_rate();
	double sigma = get_vol();

	double d1 = (std::log(S / m_K) + (r + sigma * sigma / 2.) * m_T) / (sigma * std::pow(m_T, 0.5));
	double d2 = d1 - sigma * std::pow(m_T, 0.5);

	double res = S * N(d1) - m_K * std::exp(-r * m_T) * N(d2);
	return res;
}

