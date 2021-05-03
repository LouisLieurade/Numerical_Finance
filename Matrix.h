#ifndef MATRIX_H
#define MATRIX_H
#include <string>


template <typename T=double>
class Matrix
{
    public:
        Matrix();
        explicit Matrix(int N);
        explicit Matrix(int M, int N);
        explicit Matrix(int M, int N, T* dat);
		//Matrix(const Matrix& B);
		//Matrix operator=(const Matrix &B);


        std::string to_string() const;
		void display() const;
		void display_end() const;
		int get_rows() const;
		int get_columns() const;

        T elem_at(int i, int j) const;
		T elem_at(int i) const;
        Matrix* set_elem_at(int i, int j, T val);

		Matrix inverse();
		void switch_row(int row1, int row2);
		void reduction(int row, double alpha);
		void transvection(int row1, int row2, double alpha);  
		void zeros(int row, int col);
		void inverse_row(Matrix<T>* M, int row, int col);
		bool is_still_inversible(int row, int col);

        Matrix transpose() const;
        Matrix dot(const Matrix &B) const;
		Matrix column(int j);
		void fill_column(int j, const Matrix data);
		void fill_row(int i, const Matrix data);

        Matrix operator*(const Matrix &B) const;
        Matrix operator*(const T &scal) const;
        Matrix operator+(const Matrix &B) const;
        Matrix operator+(const T &scal) const;
		Matrix operator-(const Matrix &B) const;
		Matrix operator-(const T &scal) const;
		Matrix operator+=(Matrix &B);
		bool operator==(const Matrix &B) const;
		bool isnull() const;
		T operator()(int i, int j) const; 
		T operator()(int i) const;
	

		double mean();
		double variance();
		Matrix get_diagonal();
		Matrix exp();
		Matrix Cholesky();
		//void test(); 
        virtual ~Matrix();

    protected:

    private:
        int m_rows;
        int m_cols;
        T* m_data;

};

Matrix<double> maximum(const Matrix<double> &A, const Matrix<double> &B);

#endif // MATRIX_H
