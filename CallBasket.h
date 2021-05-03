#pragma once
#include "../Matrix/Matrix.h"

class Call_Basket
{
public:
	Call_Basket(Matrix<double> S, Matrix<double> sigma, double T, double rate, double K, Matrix<double> weight);
	~Call_Basket();
	double get_spot();
	double get_rate();
	double get_vol();
	double price();

protected:
	Matrix<double> m_S;
	Matrix<double> m_sigma; 
	double m_T;
	double m_rate;
	double m_K;
	Matrix<double> m_weight;
	int m_assets;
};

