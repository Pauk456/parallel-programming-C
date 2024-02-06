#include "Solving_Linear_Equations.h"

const double Solving_Linear_Equations_virtual::ti_plus = 0.01;
const double Solving_Linear_Equations_virtual::ti_minus = -0.01;

Solving_Linear_Equations_virtual::Solving_Linear_Equations_virtual(Matrix A, std::vector<double> x, std::vector<double> b)
	: A(A), N(A.size()), x(x), b(b)
{
	norm_denominator = find_norm(b);
}

double Solving_Linear_Equations_virtual::multiply_row_by_column(const std::vector<double>& row, const std::vector<double>& column) const
{
	double result = 0.0;
	for (int i = 0; i < N; ++i) {
		result += row[i] * column[i];
	}
	return result;
}

double Solving_Linear_Equations_virtual::find_norm(const std::vector<double>& row) const
{
	double result = 0.0;
	for (int i = 0; i < row.size(); i++)
	{
		result += row[i] * row[i];
	}
	result = std::pow(result, 0.5);
	return result;
}
