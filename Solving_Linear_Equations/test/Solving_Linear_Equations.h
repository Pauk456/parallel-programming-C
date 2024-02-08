#pragma once
#include<mpi.h>
#include<vector>
#include<cmath>
#include <chrono>
#include"Matrix.h"

class Solving_Linear_Equations_virtual
{
protected:
	double ti = ti_plus;
	static const double ti_plus;
	static const double ti_minus;
	int N;
	Matrix A;
	std::vector<double> b;
	std::vector<double> x;
	double norm_denominator;

	virtual void proximity_function() = 0;
	virtual bool accuracy_check(double epsilon) const = 0;
	double multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column) const;
	double find_norm(const std::vector<double>& row) const;
public:
	virtual std::vector<double> execute(double epsilon);
	Solving_Linear_Equations_virtual(Matrix A, std::vector<double> x, std::vector<double> b);
	virtual ~Solving_Linear_Equations_virtual() = default;
};

class Solving_Linear_Equations_usual : public Solving_Linear_Equations_virtual
{
private:
	void proximity_function() override;
	bool accuracy_check(double epsilon) const override;
public:
	Solving_Linear_Equations_usual(Matrix A, std::vector<double> x, std::vector<double> b) 
		: Solving_Linear_Equations_virtual(A, x, b) {};
};

class Solving_Linear_Equations_parallel : public Solving_Linear_Equations_virtual
{
private:
	void proximity_function() override;
	bool accuracy_check(double epsilon) const override;
public:
	Solving_Linear_Equations_parallel(Matrix A, std::vector<double> x, std::vector<double> b)
		: Solving_Linear_Equations_virtual(A, x, b) {};
};
