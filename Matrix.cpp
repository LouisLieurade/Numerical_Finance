#include "stdafx.h"

#include "Matrix.h"
#include <iostream>
#include <algorithm>

template <typename T>
Matrix<T>::Matrix(){
    //ctor
}

template <typename T>
Matrix<T>::Matrix(int M, int N): m_rows(M), m_cols(N), m_data(new T[M * N]) {
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++)
            m_data[i * m_cols + j] =  T();
}

template <typename T>
Matrix<T>::Matrix(int N): m_rows(N), m_cols(N), m_data(new T[N * N]) {
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++)
            m_data[i * m_cols + j] =  T(i==j);
}

template <typename T>
Matrix<T>::Matrix(int M, int N, T* dat): m_rows(M), m_cols(N), m_data(new T[M * N]) {
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++)
            m_data[i * m_cols + j] =  dat[i * m_cols + j];
}

/*
template <typename T>
Matrix<T>::Matrix(const Matrix &B)
{
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++)
			m_data[i * m_cols + j] = B.m_data[i * m_cols + j];
}
*/

template <typename T>
std::string Matrix<T>::to_string() const{
    std::string res="";
    for(int i = 0; i < m_rows; i++){
        res+="[";
        for(int j = 0; j < m_cols; j++){
            res+=std::to_string(m_data[i * m_cols + j]);
            res+=", ";
        }
        res.erase(res.end()-2, res.end());
        res+="]\n";
    }
    return res;
}

template <typename T>
void Matrix<T>::display() const
{
	std::string temp = this->to_string(); 
	std::cout << temp << std::endl;
}

template <typename T>
void Matrix<T>::display_end() const
{
	std::string res = "";
	for (int i = 0; i < m_rows; i++) {
		res += "[";
		for (int j = m_cols - 15; j < m_cols; j++) {
			res += std::to_string(m_data[i * m_cols + j]);
			res += ", ";
		}
		res.erase(res.end() - 2, res.end());
		res += "]\n";
	}
	std::cout << res << std::endl;
}

template <typename T>
int Matrix<T>::get_rows() const
{
	return m_rows;
}
template <typename T>
int Matrix<T>::get_columns() const
{
	return m_cols;
}


template <typename T>
T Matrix<T>::elem_at(int i, int j) const{
    return m_data[i * m_cols + j];
}

template <typename T>
T Matrix<T>::elem_at(int i) const {
	return elem_at(i, 0);
}


template <typename T>
Matrix<T>* Matrix<T>::set_elem_at(int i, int j, T val){
    m_data[i * m_cols + j]=val;
    return this;
}

template <typename T>
Matrix<T> Matrix<T>::transpose() const{
    Matrix<T>* res = new Matrix<T>(m_cols,m_rows);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++){
            res->set_elem_at(j,i,elem_at(i,j));
        }
    return *res;
}

template <typename T>
void Matrix<T>::switch_row(int row1, int row2)
{
	for (int j = 0; j < m_cols; j++)
	{
		double temp = elem_at(row1, j);
		set_elem_at(row1, j, elem_at(row2, j));
		set_elem_at(row2, j, temp);
	}
}


template <typename T>
// in order the pivot is equal to 1 
void Matrix<T>::reduction(int row, double alpha)
{
	for (int j = 0; j < m_cols; j++)
	{
		set_elem_at(row, j, elem_at(row, j) * alpha);
	}
}

template <typename T>
void Matrix<T>::transvection(int row1, int row2, double alpha)
{
	for (int j = 0; j < m_cols; j++)
	{
		set_elem_at(row1, j, elem_at(row1, j) + elem_at(row2, j) * alpha);
	}
}

template <typename T>
// We create 0 abose and below the pivot
void Matrix<T>::zeros(int row, int col)
{
	for (int k = 0; k < m_rows; k++)
	{
		double c = elem_at(k, col);
		if (k != row)
		{
			for (int l = 0; l < m_cols; l++)
			{
				set_elem_at(k, l, elem_at(k, l) - c * elem_at(row, l));
			}
		}
	}
}

template <typename T>
void Matrix<T>::inverse_row(Matrix<T>* M, int row, int col)
{
	for (int k = 0; k < m_rows; k++)
	{
		double c = elem_at(k, col);
		if (k != row)
		{
			for (int l = 0; l < m_cols; l++)
			{
				M->set_elem_at(k, l, M->elem_at(k, l) - c * M->elem_at(row, l));
			}
		}
	}
}

template <typename T>
// we make sure the column is not full of 0, otherwise, it means that the matrix is not inversible
bool Matrix<T>::is_still_inversible(int row, int col)
{
	for (int l = row; l < m_rows; l++)
	{
		if (elem_at(l, col) != 0)
			return true;
	}
	return false;
}


template <typename T>
Matrix<T> Matrix<T>::inverse()
{
	Matrix<T>* Id = new Matrix<T>(m_cols);
	Matrix<T>* Rhs = new Matrix<T>(m_rows, m_cols, m_data);

	int l = 0;
	int i = 0;
	for (int j = 0; j < m_cols; j++)
	{
		//Id->display();
		
		if (!Rhs->is_still_inversible(i, j))
		{
			throw std::exception("The coef matrix is not inversible");
		}
		
		
		l = i;
		while (Rhs->elem_at(l, j) == 0)
		{
			l++;
		}
		
		//Id->display();
		Id->switch_row(l, i);
		Rhs->switch_row(l, i);
		
		Id->reduction(i, 1 / Rhs->elem_at(i, j));
		Rhs->reduction(i, 1 / Rhs->elem_at(i, j));
		
		Rhs->inverse_row(Id, i, j);
		Rhs->zeros(i, j);
		i++;
		//Id->display();
	}
	delete Rhs; 
	return *Id;
}


template <typename T>
Matrix<T> Matrix<T>::dot(const Matrix &B) const
{
    Matrix<T>* res = new Matrix<T>(m_rows,B.m_cols);

	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < B.m_cols; j++)
		{
			T tot = 0.;
			for (int k = 0; k < m_cols; k++)
			{
				tot += (elem_at(i, k)*B.elem_at(k, j));
			}
			res->set_elem_at(i, j, tot);
		}
	}
    return *res;
}

template <typename T>
Matrix<T> Matrix<T>::column(int j)
{
	Matrix<T>* res = new Matrix<T>(m_rows, 1);
	for (int i = 0; i < m_rows; i++)
		res->set_elem_at(i, 0, elem_at(i, j));
	return *res;
}

template <typename T>
void Matrix<T>::fill_column(int j, const Matrix<T> data)
{
	for (int i = 0; i < m_rows; i++)
		set_elem_at(i, j, data.elem_at(i, 0));
}

template <typename T>
void Matrix<T>::fill_row(int i, const Matrix<T> data)
{
	for (int j = 0; j < m_cols; j++)
		set_elem_at(i, j, data.elem_at(j, 0));
}

//Hadamard product
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix &B) const{
    Matrix<T>* res = new Matrix<T>(m_rows,m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++){
            res->set_elem_at(i,j,B.elem_at(i,j) * elem_at(i,j));
        }
    return *res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T &scal) const{
    Matrix<T>* res = new Matrix<T>(m_rows,m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++){
            res->set_elem_at(i,j,scal * elem_at(i,j));
        }
    return *res;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix &B) const{
    Matrix<T>* res = new Matrix<T>(m_rows,m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++){
            res->set_elem_at(i,j,B.elem_at(i,j) + elem_at(i,j));
        }
    return *res;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const T &scal) const{
    Matrix<T>* res = new Matrix<T>(m_rows,m_cols);
    for(int i = 0; i < m_rows; i++)
        for(int j = 0; j < m_cols; j++){
            res->set_elem_at(i,j,scal + elem_at(i,j));
        }
    return *res;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix &B) const {
	Matrix<T>* res = new Matrix<T>(m_rows, m_cols);
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++) {
			res->set_elem_at(i, j, elem_at(i, j) - B.elem_at(i, j));
		}
	return *res;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const T &scal) const {
	Matrix<T>* res = new Matrix<T>(m_rows, m_cols);
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++) {
			res->set_elem_at(i, j, elem_at(i, j) - scal);
		}
	return *res;
}


template <typename T>
Matrix<T> Matrix<T>::operator+=(Matrix &B)
{
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++) {
			set_elem_at(i, j, B.elem_at(i, j) + elem_at(i, j));
		}
	return *this;
}

/*
template <typename T>
Matrix<T> Matrix<T>::operator=(const Matrix &B)
{
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++)
			m_data[i * m_cols + j] = B.m_data[i * m_cols + j];
	return *this;
}
*/

template <typename T>
bool Matrix<T>::operator==(const Matrix &B) const
{
	if (m_rows != B.m_rows)
		return false;
	if (m_cols != B.m_cols)
		return false;
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++)
			if (elem_at(i, j) != B.elem_at(i, j))
				return false;
	return true;
}

template <typename T>
bool Matrix<T>::isnull() const
{
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++)
			if (elem_at(i, j) != 0)
				return false;
	return true;
}

template <typename T>
T Matrix<T>::operator()(int i, int j) const
{
	//return this->elem_at(i, j);
	return m_data[i * m_cols + j];
}

template <typename T>
T Matrix<T>::operator()(int i) const
{
	//return this->elem_at(i, j);
	return m_data[i * m_cols];
}



template <typename T>
double Matrix<T>::mean()
{
	double res = 0;
	for (int i = 0; i < m_rows; i++)
	{
		for (int j = 0; j < m_cols; j++)
		{
			res += elem_at(i, j);
		}
	}
	return res / (m_cols*m_rows);
}	

template <typename T>
double Matrix<T>::variance()
{
	double res = 0;
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++)
			res += std::pow(elem_at(i, j), 2);
	res = (res / (m_cols*m_rows)) - std::pow(mean(), 2) ;
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::get_diagonal()
{
	Matrix<T> res(m_rows, 1);
	for (int i = 0; i < m_rows; i++)
		res.set_elem_at(i, 0, elem_at(i, i));
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::exp()
{
	for (int i = 0; i < m_rows; i++)
		for (int j = 0; j < m_cols; j++)
			set_elem_at(i, j, std::exp(elem_at(i, j)));
	return *this;
}

template <typename T>
Matrix<T> Matrix<T>::Cholesky()
{
	if (m_rows != m_cols)
	{
		throw std::exception("The matrix is not square");
	}
	// Initialisation
	Matrix<T> res(m_rows, m_cols);
	// first element filled
	res.set_elem_at(0, 0, std::pow(elem_at(0, 0), 0.5));
	// first column filled
	for (int j = 1; j < m_rows; j++)
	{
		res.set_elem_at(j, 0, elem_at(0, j) / res.elem_at(0, 0));
	}
	// recursive scheme to compute column j+1 knowing column j 
	for (int i = 1; i < m_rows; i++)
	{
		// compute first the diagonal term
		double lik_sum = 0; 
		for (int k = 0; k <= i - 1; k++)
		{
			lik_sum += std::pow(res.elem_at(i, k), 2);
		}
		double nombre = std::pow(elem_at(i, i) - lik_sum, 0.5);
		res.set_elem_at(i, i, std::pow(elem_at(i, i) - lik_sum, 0.5));
		//  compute the full column under the diagonal term
		// loop for j the column
		for (int j = i + 1; j < m_rows; j++)
		{
			lik_sum = 0;
			for (int k = 0; k <= i - 1; k++)
			{
				lik_sum += res.elem_at(i, k) * res.elem_at(j, k);
			}
			res.set_elem_at(j, i, (elem_at(i, j) - lik_sum) / res.elem_at(i, i));
		}
	}

	return res;
}



template <typename T>
Matrix<T>::~Matrix(){
	m_data = nullptr;
}

Matrix<double> maximum(const Matrix<double> &A, const Matrix<double> &B)
{
	Matrix<double> res(A.get_rows(), A.get_columns());
	double n;
	for (int i = 0; i < A.get_rows(); i++)
	{
		for (int j = 0; j < A.get_columns(); j++)
		{
			n = std::max(A.elem_at(i, j), B.elem_at(i, j));
			res.set_elem_at(i, j, n);
		}
	}
	return res;
}

template class Matrix<int>;
template class Matrix<long>;
template class Matrix<float>;
template class Matrix<double>;


