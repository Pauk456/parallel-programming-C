#pragma once
#include<mpi.h>
#include<vector>
#include<cmath>
#include"Matrix.h"

class Solving_Linear_Equations_interface
{
private:
	virtual double multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column) const = 0;
	virtual double find_norm(const std::vector<double>& row) const = 0;
	virtual void proximity_function() = 0;
	virtual bool accuracy_check(double epsilon) const = 0;
public:
	static const double ti_plus;
	static const double ti_minus;

	virtual std::vector<double> execute(double epsilon) = 0;
	virtual ~Solving_Linear_Equations_interface() = default;
};


class Solving_Linear_Equations_usual : Solving_Linear_Equations_interface
{
private:
	int N;
	Matrix A;
	std::vector<double> b;
	std::vector<double> x;
	double norm_denominator;
private:
	double multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column) const override;
	double find_norm(const std::vector<double>& row) const override;
	void proximity_function() override;
	bool accuracy_check(double epsilon) const override;
public:
	Solving_Linear_Equations_usual(Matrix A, std::vector<double> x, std::vector<double> b);
	virtual std::vector<double> execute(double epsilon) override;
};
