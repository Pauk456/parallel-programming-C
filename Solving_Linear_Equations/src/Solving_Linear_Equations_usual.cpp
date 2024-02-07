#include "Solving_Linear_Equations.h"

void Solving_Linear_Equations_usual::proximity_function()
{
	std::vector<double> new_x(N);
	for (int i = 0; i < N; i++)
	{
		new_x[i] = multiply_row_by_column(A[i], x);
	}

	for (int i = 0; i < N; i++)
	{
		new_x[i] = x[i] - ti_plus * (new_x[i] - b[i]);
	}
	x = new_x;
}

bool Solving_Linear_Equations_usual::accuracy_check(double epsilon) const
{
	std::vector<double> result(N);
	for (int i = 0; i < N; i++)
	{
		result[i] = multiply_row_by_column(A[i], x);
	}

	for (int i = 0; i < N; i++)
	{
		result[i] = (result[i] - b[i]);
	}

	double norm_numerator = find_norm(result);
	return norm_numerator / norm_denominator < epsilon ? true : false;
}

