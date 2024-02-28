#pragma once
#include<mpi.h>
#include<iostream>
#include<vector>
#include<cmath>
#include <chrono>
#include"Matrix.h"

extern std::chrono::microseconds durationGLOBAL;

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
	double norm_denominator = 0;

	virtual void proximity_function() = 0;
	virtual bool accuracy_check(double epsilon) = 0;
	double multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column) const;
	double find_norm(const std::vector<double>& row) const;
public:
	Solving_Linear_Equations_virtual(Matrix A) : A(A), N(A.size()) {};
	Solving_Linear_Equations_virtual(Matrix A, std::vector<double> x, std::vector<double> b);
	virtual ~Solving_Linear_Equations_virtual() = default;

	virtual void print_result() = 0;
	virtual void execute(double epsilon);
};

class Solving_Linear_Equations_usual : public Solving_Linear_Equations_virtual
{
private:
	void proximity_function() override;
	bool accuracy_check(double epsilon) override;
public:
	Solving_Linear_Equations_usual(Matrix A, std::vector<double> x, std::vector<double> b) 
		: Solving_Linear_Equations_virtual(A, x, b) {};

	void print_result() override;
};

class Solving_Linear_Equations_parallel_first : public Solving_Linear_Equations_virtual
{
private:
	std::vector<double> result_part;
	std::vector<double> result;

	int size, rank, ibeg, iend, count_for_process;
	void proximity_function() override;
	bool accuracy_check(double epsilon) override;
public:
	Solving_Linear_Equations_parallel_first(Matrix A, std::vector<double> x, std::vector<double> b);

	void print_result() override;
};

class Solving_Linear_Equations_parallel_second : public Solving_Linear_Equations_virtual
{
private:
	std::vector<double> x_process;
	std::vector<double> result;

	int size, rank, ibeg, iend, count_for_process, destination, sender;

	int find_norm_b();
	void proximity_function() override;
	bool accuracy_check(double epsilon) override;
	double multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column, int offset, int count) const;
public:
	Solving_Linear_Equations_parallel_second(Matrix A, std::vector<double> x, std::vector<double> b);

	void print_result() override;
};

