#pragma once
#include<mpi.h>
#include<vector>
#include <cmath>

class Matrix
{
private:
	std::vector< std::vector<double> > matrix;
	int N;
public:
	Matrix(int N);

	int size() {
		return N;
	}
	const std::vector<double>& operator[](int index) const;
	std::vector<double>& operator[](int index);
};
//
//class Value_vector
//{
//private:
//	std::vector<double> value_vector;
//	int N;
//public:
//	void initialize_vector();
//};

class Solving_Linear_Equations
{
private:
	static const double ti_plus;
	static const double ti_minus;

	int N;
	Matrix A;
	std::vector<double> b;
	std::vector<double> x;
	double norm_denominator;
	//Value_vector b;
	//Value_vector x;
private:
	double multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column) const;
	double find_norm(const std::vector<double>& row) const;
	void proximity_function();
	bool accuracy_check(double epsilon) const;
public:
	Solving_Linear_Equations(Matrix A, std::vector<double> x, std::vector<double> b);
	
	std::vector<double> execute(double epsilon);
};
